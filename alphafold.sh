#!/bin/bash
set -e

alphafold_version='2.3.2_conda'
db_version='2020-05-14'

echo "\
---------------------------------------------------------
AlphaFold2-Multimer pipeline for a SLURM compute cluster.
---------------------------------------------------------
Author: Eugene Valkov, NCI/NIH.
AlphaFold2 version: $alphafold_version                        
Template release date cutoff: $db_version
---------------------------------------------------------
"

# Check if the 'sbatch' command is available and get its version
check_slurm=$(sbatch -V 2>/dev/null)

if [ -z "$check_slurm" ]; then
    echo -e "This pipeline requires access to a SLURM-enabled compute environment."
    exit 1
else
    echo -e "SLURM is available. Version: $check_slurm"
fi


function show_help() {
	echo ""
  	echo "Usage: $0 [options]"
  	echo ""
  	echo "Options:"
  	echo "  -d DIRECTORY   Specify the directory path. (Required)"
  	echo "  -f FILE        Specify the sequence file path. (Required)"
	echo "  -m MODE        Specify 'quick' or 'thorough'. ([1] or 5 predictions/model)."
  	echo "  -h             Show help information."
}

directory=""
seqfile=""
mode="quick"
num_pred_per_model=1

while getopts "hd:f:m:" opt; do
  case $opt in
    h) show_help
       exit 0
    ;;
    d) directory="$OPTARG"
    ;;
    f) seqfile="$OPTARG"
    ;;
    m) if [[ "$OPTARG" == "quick" ]]; then
	 mode="$OPTARG"
         num_pred_per_model=1
       elif [[ "$OPTARG" == "thorough" ]]; then
         mode="$OPTARG"
         num_pred_per_model=5
       else
         echo "Invalid mode. Choose 'quick' or 'thorough'."
         exit 1
       fi
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
       exit 1
    ;;
  esac
done

# Check if both -d and -f options were provided
if [ -z "$directory" ] || [ -z "$seqfile" ]; then
  echo -e "Both -d (directory) and -f (file) options must be provided."
  show_help
  exit 1
fi

# Check if the directory exists, is a directory, and is writable
if [ ! -d "$directory" ] || [ ! -w "$directory" ]; then
	echo "The specified directory does not exist, is not a directory, or is not writable: $directory"
	exit 1
else
	mkdir -p ""$directory"$USER/"
        procdir=""$directory"$USER/"
fi

# Check if the file exists and is a file
if [ ! -f "$seqfile" ]; then
  echo "The specified file does not exist or is not a file: $seqfile"
  exit 1
fi

echo "Directory for output: $procdir"
echo "File provided is: $seqfile"
echo "$num_pred_per_model predictions/model will be generated."

# Checks if you are logged into the head node for the cluster
submithost=`echo $HOSTNAME`
if [ ! "$submithost" = "fsitgl-head01p.ncifcrf.gov" ]; then
	echo -e "\nYou must be logged in to FRCE cluster to use this script.\n"
	exit
fi


if [ "$seqfile" ]; then
	af2dir=`echo $seqfile | sed 's|.*/||' | sed 's/\.[^.]*$//'`
	cp $seqfile "$procdir""$af2dir".fa
	num_seqs=`grep '^>' $seqfile | wc -l`
	chain_names=`grep '^>' $seqfile`
	if [ ! "$chain_names" ]; then
        	echo -e "$seqfile is not in fasta format:"
		cat $seqfile
                exit
        else
		echo -e "\
Found $num_seqs protein chains:
$chain_names"
	fi
else
	echo -e "Fasta sequence not provided!\n"
        exit
fi

len=`sed '/^>/d' $seqfile | tr -d '\n' | wc -c`

# Checks if the sequence contains only single-letter amino acid codes

valid_codes="GPVALIMCFYWHKRNQEDST"
found_error=0

if [ ! "$seqfile" = "" ]; then
awk -v valid="$valid_codes" '/^>/ {next} {
  gsub(/[^[:alnum:]]/, "")
  for (i = 1; i <= length; i++) {
    char = substr($0, i, 1)
    if (index(valid, char) == 0) {
      print "Non-amino acid letter in line:", NR, "Character:", char
      found_error=1
      exit 1
    }
  }
}

END {
  if (found_error == 0) {
    print "The file contains only valid amino acid codes."
  }
}' $seqfile
fi


echo "\
#!/bin/bash
#SBATCH --job-name=$af2dir
#SBATCH --output=$af2dir.out
#SBATCH --partition=gpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=15G" > $procdir"$af2dir"_af2.sh

if (( "$len" <= 249 )); then
  echo -e "Found $len residues, setting 6h time limit."
  echo "\
#SBATCH --gres=gpu:2
#SBATCH --time=06:00:00
" >> $procdir"$af2dir"_af2.sh
elif (( "$len" >= 250 && "$len" <= 999 )); then
  echo -e "Found $len residues, setting 24h time limit."
  echo "\
#SBATCH --gres=gpu:v100:2
#SBATCH --time=24:00:00
" >> $procdir"$af2dir"_af2.sh
elif (( "$len" >= 1000 && "$len" <= 1299 )); then
  echo -e "Found $len residues, setting 36h time limit."
  echo "\
#SBATCH --gres=gpu:v100:2
#SBATCH --time=48:00:00
" >> $procdir"$af2dir"_af2.sh
elif (( "$len" >= 1300 && "$len" <= 2499 )); then
  echo -e "Found $len residues, setting 48h time limit."
  echo "\
#SBATCH --gres=gpu:v100:2
#SBATCH --time=2-00:00:00
" >> $procdir"$af2dir"_af2.sh
elif (( "$len" >= 2500 )); then
  num_pred_per_model=1
  echo -e "Large prediction so will run in 'quick' mode only."
  echo -e "Found $len residues, setting 72h time limit."
  echo -e "The GPU is set to offload memory."
  echo -e "Only $num_pred_per_model prediction/model will be generated." 
  echo "\
#SBATCH --gres=gpu:v100:2
#SBATCH --time=3-00:00:00

export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=\"4.0\"
" >> $procdir"$af2dir"_af2.sh 
fi

echo "\
module load alphafold/$alphafold_version

run --fasta_paths=$procdir$af2dir.fa \\
--output_dir=$procdir \\
--db_preset=full_dbs \\
--max_template_date=$db_version \\
--models_to_relax=best \\
--model_preset=multimer \\" >> $procdir"$af2dir"_af2.sh

echo "--num_multimer_predictions_per_model=$num_pred_per_model
" >> $procdir"$af2dir"_af2.sh

echo "\
if [ ! -e "$procdir$af2dir/ranked_0.pdb" ]; then
	# workaround with libcrypto not accessing the correct openssl 
	export LD_LIBRARY_PATH=/lib64:\$LD_LIBRARY_PATH
	tail -50 $af2dir.out | mutt -s \"$af2dir\" -e 'my_hdr From:AlphaFold2 (AlphaFold2)' -- $USER@nih.gov
	exit
fi
" >> $procdir"$af2dir"_af2.sh


echo "\
set bgcolor white
open $procdir$af2dir/ranked_0.pdb
cartoon style protein modeh tube rad 2 sides 24
cartoon style width 2 thick 0.2
rainbow chain palette RdYlBu-5
lighting simple shadows false intensity 0.5
show #1 cartoons
view all
hide atoms
movie record size 1500,1500 format png supersample 4 directory $procdir$af2dir
show #1 models
turn y 6 360 models #1
wait 60
stop
movie stop
" > $procdir"$af2dir"_chimera_movie.cxc

echo "\
module load ChimeraX/1.5
ChimeraX --offscreen --script $procdir"$af2dir"_chimera_movie.cxc --exit
convert -dispose previous -delay 10 -loop 0 -dither None -colors 256 -layers Optimize -resize 500x500 -filter Lanczos -coalesce $procdir$af2dir/chimovie*.png $procdir$af2dir/animated.gif
rm $procdir$af2dir/chimovie*.png
" >> $procdir"$af2dir"_af2.sh


echo "\
set bgcolor white
open $procdir$af2dir/ranked*.pdb
cartoon style protein modeh tube rad 2 sides 24
cartoon style width 2 thick 0.2
rainbow chain palette RdYlBu-5
lighting simple shadows false intensity 0.5
view all
hide atoms
show cartoons
alphafold pae #1 file $procdir$af2dir/pae_model_1_multimer_v3_pred_0.json palette paegreen plot false
hide all models
show #1 models
" > $procdir"$af2dir"_chimera_align.cxc


echo "\
# Initialize variables to track the chain with the most atoms
max_atoms=0
max_chain=\"\"

# Process the PDB file and store the result in a variable
readarray -t lines < <(grep '^ATOM' $procdir$af2dir/ranked_0.pdb | cut -c 22 | sort | uniq -c | sort -nr)

# Iterate over the lines
for line in \"\${lines[@]}\"; do
    read -r atoms chain <<< \"\$line\"
    if (( atoms > max_atoms )); then
        max_atoms=\$atoms
        max_chain=\$chain
    fi
done

echo \"\
matchmaker all to #1/\$max_chain pairing bs
save $procdir"$af2dir"_chimera_align.cxs
\" >> $procdir"$af2dir"_chimera_align.cxc
ChimeraX --offscreen --script $procdir"$af2dir"_chimera_align.cxc --exit
" >> $procdir"$af2dir"_af2.sh


echo "\
set bgcolor white
open $procdir$af2dir/ranked_0.pdb
cartoon style protein modeh tube rad 2 sides 24
cartoon style width 2 thick 0.2
rainbow chain palette RdYlBu-5
lighting simple shadows false intensity 0.5
show #1 cartoons
view all
hide atoms
hide #1 models
" > $procdir"$af2dir"_chimera_pae.cxc


echo "\
# Initialize an empty array
chains=()
# Extract chain identifiers and store in an array
while IFS= read -r line; do
	chains+=(\"\$line\")
done < <(grep '^ATOM' $procdir$af2dir/ranked_0.pdb | cut -c 22 | sort -u)

# Print the array in the desired format
count=1
# Loop through the list
for (( i=0; i<\${#chains[@]}; i++ )); do
	for (( j=i+1; j<\${#chains[@]}; j++ )); do
		count=\$((count + 1))
		echo \"\
combine #1 modelId #\$count name 'interface chains \${chains[i]} to \${chains[j]}'
sel #\$count/\${chains[i]},\${chains[j]}
select ~sel & ##selected
alphafold pae #\$count file $procdir$af2dir/pae_model_1_multimer_v3_pred_0.json plot false
alphafold contacts #\$count/\${chains[i]} to #\$count/\${chains[j]} distance 8 palette paecontacts range 0,30 radius 0.05 dashes 1
del sel
save $procdir"$af2dir"_chimera_pae.cxs
\" >> $procdir"$af2dir"_chimera_pae.cxc
	done
done
ChimeraX --offscreen --script $procdir"$af2dir"_chimera_pae.cxc --exit
" >> $procdir"$af2dir"_af2.sh


echo "\
mv $procdir"$af2dir"*.cxs $procdir"$af2dir"/
mv $procdir"$af2dir"_af2.sh $procdir"$af2dir"/
mv $procdir"$af2dir".fa $procdir"$af2dir"/

export_dir=\""$af2dir"_$(date +"%Y%m%d_%H%M%S")\"
mkdir $procdir\$export_dir
cp $procdir$af2dir/ranked*.pdb $procdir\$export_dir
cp $procdir$af2dir/pae_model_1_multimer_v3_pred_0.json $procdir\$export_dir
" >> $procdir"$af2dir"_af2.sh


echo "\
sed -e 's|$procdir$af2dir/||g' \\
    -e '/save /d' \\
    -e 's/ plot false//' \\
    $procdir"$af2dir"_chimera_align.cxc > $procdir\$export_dir/"$af2dir"_chimera_align.cxc

sed -e 's|$procdir$af2dir/||g' \\
    -e '/save /d' \\
    -e 's/ plot false//' \\
    $procdir"$af2dir"_chimera_pae.cxc > $procdir\$export_dir/"$af2dir"_chimera_pae.cxc

tar cjf $procdir"$af2dir".tar.bz2 -C $procdir \$export_dir
mv $procdir"$af2dir"*.cxc $procdir"$af2dir"/
" >> $procdir"$af2dir"_af2.sh

echo "\
maxsize=\$((7*1024*1024)) # 7 MB in bytes
size1=\$(stat -c%s \"$procdir"$af2dir".tar.bz2\")
size2=\$(stat -c%s \"$procdir$af2dir/animated.gif\")
totalsize=\$((size1 + size2))

# workaround with libcrypto not accessing the correct openssl 
export LD_LIBRARY_PATH=/lib64:\$LD_LIBRARY_PATH

if [ \$totalsize -lt \$maxsize ]; then
	echo -e \"<img src=\"cid:animated.gif\" />\" | mutt -e 'set content_type=text/html' -s \"$af2dir\" -a $procdir$af2dir/animated.gif -a $procdir"$af2dir".tar.bz2 -e 'my_hdr From:AlphaFold2 (AlphaFold2)' -- $USER@nih.gov
else
	echo -e \"<img src=\"cid:animated.gif\" />\" | mutt -e 'set content_type=text/html' -s \"$af2dir\" -a $procdir$af2dir/animated.gif -e 'my_hdr From:AlphaFold2 (AlphaFold2)' -- $USER@nih.gov
fi
rm $procdir"$af2dir".tar.bz2
rm -rf $procdir\$export_dir
" >> $procdir"$af2dir"_af2.sh

if [ -e $HOME/.netrc ]; then
	echo "\
# Start lftp session
lftp ftp.box.com << EOF
cd Alphafold
mirror -R --no-symlinks $procdir$af2dir
bye
EOF
" >> $procdir"$af2dir"_af2.sh
	else
        echo "\
rsync -vagu $procdir$af2dir $storage_dir
" >> $procdir"$af2dir"_af2.sh
fi

sbatch $procdir"$af2dir"_af2.sh
