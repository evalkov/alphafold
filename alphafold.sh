#!/bin/bash

set -e

alphafold_version='2.3.2_conda'
db_version='2020-05-14'

################### DO NOT EDIT BELOW ################

mkdir -p "/scratch/cluster_scratch/$USER/alphafold"
procdir="/scratch/cluster_scratch/$USER/alphafold/"

if [[ -r /mnt/projects/RNABL-GRID-SPE/active/Valkov/alphafold/ ]]; then
	storage_dir="/mnt/projects/RNABL-GRID-SPE/active/Valkov/alphafold/"
else
	storage_dir="$procdir"
fi

# Checks if you are logged into the head node for the cluster
submithost=`echo $HOSTNAME`
if [ ! "$submithost" = "fsitgl-head01p.ncifcrf.gov" ]; then
	echo -e "\nYou must be logged in to FRCE cluster to use this script.\n"
	exit
fi

echo -e "\

Script to submit AlphaFold2 jobs to the FRCE cluster.
(C) Eugene Valkov, National Cancer Institute, U.S.A.
AlphaFold version: $alphafold_version                        
Template release date cutoff: $db_version

Usage: alphafold fastafile.fa (add -quick or -thorough flags to generate 5 or 25 predictions)\n"

# Initialize flag variables
quick_flag=false
thorough_flag=false

# Process command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -quick)
            quick_flag=true
            shift
            ;;
        -thorough)
            thorough_flag=true
            shift
            ;;
        *)
            # Treat any other argument as a file to be read
            if [ -f "$1" ]; then
                echo "Reading file: $1"
                seqfile="$1"
            else
                echo -e "Fasta sequence not provided!\n"
                exit
            fi
            shift
            ;;
    esac
done


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
}' "$seqfile"
fi

echo "\
#!/bin/bash
#SBATCH --job-name=$af2dir
#SBATCH --output="$af2dir".out
#SBATCH --partition=gpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user="$USER"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=15G" > "$procdir"/"$af2dir"_af2.sh

if (( "$len" <= 249 )); then
  echo -e "Found $len residues, setting 6h time limit."
  echo "\
#SBATCH --gres=gpu:2
#SBATCH --time=06:00:00
" >> "$procdir"/"$af2dir"_af2.sh
elif (( "$len" >= 250 && "$len" <= 999 )); then
  echo -e "Found $len residues, setting 24h time limit."
  echo "\
#SBATCH --gres=gpu:v100:2
#SBATCH --time=24:00:00
" >> "$procdir"/"$af2dir"_af2.sh
elif (( "$len" >= 1000 && "$len" <= 1299 )); then
  echo -e "Found $len residues, setting 36h time limit."
  echo "\
#SBATCH --gres=gpu:v100:2
#SBATCH --time=48:00:00
" >> "$procdir"/"$af2dir"_af2.sh
elif (( "$len" >= 1300 && "$len" <= 2499 )); then
  echo -e "Found $len residues, setting 48h time limit."
  echo "\
#SBATCH --gres=gpu:v100:2
#SBATCH --time=2-00:00:00
" >> "$procdir"/"$af2dir"_af2.sh
elif (( "$len" >= 2500 )); then
  echo -e "Found $len residues, setting 72h time limit."
  echo -e "The GPU is set to offload memory."
  echo -e "Only 5 models will be predicted." 
  echo "\
#SBATCH --gres=gpu:v100:2
#SBATCH --time=3-00:00:00

export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=\"4.0\"
" >> "$procdir"/"$af2dir"_af2.sh 
fi

echo "\
module load alphafold/"$alphafold_version"
module load pymol/2.6.0

run --fasta_paths="$procdir"/"$af2dir".fa \\
--output_dir="$procdir" \\
--db_preset=full_dbs \\
--max_template_date="$db_version" \\
--models_to_relax=best \\
--model_preset=multimer \\" >> "$procdir"/"$af2dir"_af2.sh


# Checking flags to set number of predictions per model
if [ "$quick_flag" = true ]; then
        echo -e "The -quick flag is recognized."
        echo -e "Setting multimer predictions per model to 1."
        num_pred_per_model=1
elif [ "$thorough_flag" = true ]; then
        echo -e "The -thorough flag is recognized."
        echo -e "Setting multimer predictions per model to 5."
        num_pred_per_model=5
else
	if (( "$len" <= 2499 )); then
		num_pred_per_model=5
	elif (( "$len" >= 2500 )); then
		num_pred_per_model=1
	fi
fi

echo  "--num_multimer_predictions_per_model=$num_pred_per_model" >> "$procdir"/"$af2dir"_af2.sh

echo "\
if [ ! -e ""$procdir"/"$af2dir"/ranked_0.pdb" ]; then
	tail -50 "$af2dir".out | mutt -s \"$af2dir\" -e 'my_hdr From:AlphaFold2 (AlphaFold2)' -b valkove2@nih.gov -- "$USER"@nih.gov
	exit
fi
" >> "$procdir"/"$af2dir"_af2.sh

echo "\
ln -s "$procdir"/"$af2dir"/result_model_1_*_0.pkl "$procdir"/"$af2dir"/result_model_1.pkl
ln -s "$procdir"/"$af2dir"/result_model_2_*_0.pkl "$procdir"/"$af2dir"/result_model_2.pkl
ln -s "$procdir"/"$af2dir"/result_model_3_*_0.pkl "$procdir"/"$af2dir"/result_model_3.pkl
ln -s "$procdir"/"$af2dir"/result_model_4_*_0.pkl "$procdir"/"$af2dir"/result_model_4.pkl
ln -s "$procdir"/"$af2dir"/result_model_5_*_0.pkl "$procdir"/"$af2dir"/result_model_5.pkl
mv "$procdir"/"$af2dir".fa "$procdir"/"$af2dir"/
cp "$af2dir".out "$procdir"/"$af2dir"/
mv "$procdir"/"$af2dir"_af2.sh "$procdir"/"$af2dir"/
" >> "$procdir"/"$af2dir"_af2.sh

echo "\
load "$procdir"/"$af2dir"/ranked_0.pdb
show cartoon
set cartoon_cylindrical_helices, 1
spectrum chain
set antialias, 1
orient
zoom complete=1
clip slab, 500
set ray_shadows, 0
viewport 1000,1000
set ray_trace_mode, 0
mset 1 x60
movie.nutate(1,60,angle=120)
mpng "$procdir"/"$af2dir"/test.png
" > "$procdir"/"$af2dir"_pymol.pml

echo "\
cp "$procdir"/"$af2dir"_pymol.pml "$procdir"/"$af2dir"/
pymol -qc "$procdir"/"$af2dir"/"$af2dir"_pymol.pml
convert -dispose previous -delay 10 -loop 0 "$procdir"/"$af2dir"/test*.png -coalesce -scale 800x800 "$procdir"/"$af2dir"/animated.gif
rm "$procdir"/"$af2dir"/test*.png
" >> "$procdir"/"$af2dir"_af2.sh

echo "\
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
import pickle

def get_pae_plddt(model_names):
    out = {}
    for i, name in enumerate(model_names):
        d = pickle.load(open(name, 'rb'))
        out[f'model_{i+1}'] = {'plddt': d['plddt'], 'pae':d['predicted_aligned_error']}
    return out

def generate_output_images(feature_dict, out_dir, name, pae_plddt_per_model):
    msa = feature_dict['msa']
    seqid = (np.array(msa[0] == msa).mean(-1))
    seqid_sort = seqid.argsort()
    non_gaps = (msa != 21).astype(float)
    non_gaps[non_gaps == 0] = np.nan
    final = non_gaps[seqid_sort] * seqid[seqid_sort, None]

    ##################################################################
    plt.figure(figsize=(28, 8))
    ##################################################################
    plt.subplot(1, 2, 1)
    plt.title(\"Sequence coverage\")
    plt.imshow(final,
               interpolation='nearest', aspect='auto',
               cmap=\"rainbow_r\", vmin=0, vmax=1, origin='lower')
    plt.plot((msa != 21).sum(0), color='black')
    plt.xlim(-0.5, msa.shape[1] - 0.5)
    plt.ylim(-0.5, msa.shape[0] - 0.5)
    plt.colorbar(label=\"Sequence identity to query\")
    plt.xlabel(\"Positions\")
    plt.ylabel(\"Sequences\")

    ##################################################################
    plt.subplot(1, 2, 2)
    plt.title(\"Predicted LDDT per position\")
    for model_name, value in pae_plddt_per_model.items():
        plt.plot(value[\"plddt\"], label=model_name)
    plt.legend()
    plt.ylim(0, 100)
    plt.ylabel(\"Predicted LDDT\")
    plt.xlabel(\"Positions\")
    svg_filename = f\"{out_dir}/{name+('_' if name else '')}coverage_LDDT.svg\"
    plt.savefig(svg_filename, dpi=600)  # Save as SVG
    png_filename = f\"{out_dir}/{name+('_' if name else '')}coverage_LDDT.png\"
    plt.savefig(png_filename, dpi=100)  # Save as PNG
    ##################################################################

    ##################################################################
    num_models = 5
    plt.figure(figsize=(6 * num_models, 4))
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        plt.subplot(1, num_models, n + 1)
        plt.title(model_name)
        plt.imshow(value[\"pae\"], label=model_name, cmap=\"coolwarm\", vmin=0, vmax=30)
        plt.colorbar()
    svg_filename = f\"{out_dir}/{name+('_' if name else '')}PAE.svg\"
    plt.savefig(svg_filename, dpi=600)  # Save as SVG
    png_filename = f\"{out_dir}/{name+('_' if name else '')}PAE.png\"
    plt.savefig(png_filename, dpi=100)  # Save as PNG
    ##################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir', dest='input_dir', required=True)
parser.add_argument('--name', dest='name')
parser.set_defaults(name='')
parser.add_argument('--output_dir', dest='output_dir')
parser.set_defaults(output_dir='')
args = parser.parse_args()

feature_dict = pickle.load(open(f'{args.input_dir}/features.pkl', 'rb'))
is_multimer = ('result_model_1_multimer.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
is_ptm = ('result_model_1_ptm.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
model_names = [f'{args.input_dir}/result_model_{f}{\"_multimer\" if is_multimer else \"_ptm\" if is_ptm else \"\"}.pkl' for f in range(1, 6)]

pae_plddt_per_model = get_pae_plddt(model_names)
generate_output_images(feature_dict, args.output_dir if args.output_dir else args.input_dir, args.name, pae_plddt_per_model)
" > "$procdir"/"$af2dir"_vis.py

echo "\
cp "$procdir"/"$af2dir"_vis.py "$procdir"/"$af2dir"/
python3 "$procdir"/"$af2dir"/"$af2dir"_vis.py --input_dir "$procdir"/"$af2dir"/ --name "$af2dir"
tar -C "$procdir"/"$af2dir"/ -cvjf "$procdir"/"$af2dir"/"$af2dir"_top_ranked.tar.bz2 ranked_0.pdb ranked_1.pdb ranked_2.pdb ranked_3.pdb ranked_4.pdb 
echo -e \"<img src=\"cid:animated.gif\" />\" | mutt -e 'set content_type=text/html' -s \"$af2dir\" -a "$procdir"/"$af2dir"/"$af2dir"_top_ranked.tar.bz2 -a "$procdir"/"$af2dir"/*.png -a "$procdir""$af2dir"/animated.gif -e 'my_hdr From:AlphaFold2 (AlphaFold2)' -b eugene.valkov@gmail.com -- "$USER"@nih.gov
rm "$procdir"/"$af2dir"/"$af2dir"_top_ranked.tar.bz2
rm "$procdir"/"$af2dir"_vis.py
rm "$procdir"/"$af2dir"_pymol.pml
rm "$procdir"/"$af2dir"/result_model_1.pkl
rm "$procdir"/"$af2dir"/result_model_2.pkl
rm "$procdir"/"$af2dir"/result_model_3.pkl
rm "$procdir"/"$af2dir"/result_model_4.pkl
rm "$procdir"/"$af2dir"/result_model_5.pkl
" >> "$procdir"/"$af2dir"_af2.sh

echo "\

Top-scoring predictions will be emailed to $USER@nih.gov."


if [ -e "$HOME"/.boxpassword ]; then
	username=`echo "$USER"@nih.gov`
	password=`awk '{print $1}' $HOME/.boxpassword`
	echo "\
set ftp:ssl-force true; 
set mirror:parallel-directories true;
connect ftp://ftp.box.com;
user '$username' '$password';
cd Alphafold;
mirror -R --no-symlinks "$procdir"/"$af2dir";
bye" > $HOME/."$af2dir".lftpconfig
	chmod 600 $HOME/."$af2dir".lftpconfig
        echo "\
lftp -f $HOME/."$af2dir".lftpconfig
rm $HOME/."$af2dir".lftpconfig
" >> "$procdir"/"$af2dir"_af2.sh
echo "\
All files will be copied to box.com under $USER@nih.gov account."
	else
        echo "\
rsync -vagu "$procdir"/"$af2dir" "$storage_dir"
" >> "$procdir"/"$af2dir"_af2.sh
        echo "\
All files will be copied to:
$storage_dir$af2dir"
fi

sbatch "$procdir"/"$af2dir"_af2.sh
