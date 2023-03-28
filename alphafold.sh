#!/bin/bash

alphafold_version='2.3.1_conda'
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

Usage: alphafold fastafile.fa\n"

if [ ! "$1" = "" ]; then
	seqfile="$1"
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
elif [ "$1" = "" ]; then
	echo -e "Fasta sequence not provided!\n"
	exit
fi

len=`sed '/^>/d' $1 | tr -d '\n' | wc -c`


if (( "$len" <= 250 )); then
  echo -e "Found $len residues, setting 3h time limit."
  echo "\
#!/bin/bash
#SBATCH --job-name=$af2dir
#SBATCH --output="$af2dir".out
#SBATCH --partition=gpu
#SBATCH --gres=gpu:p100:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=03:00:00" > "$procdir"/"$af2dir"_af2.sh
elif (( "$len" >= 250 && "$len" <= 1500 )); then
  echo -e "Found $len residues, setting 16h time limit."
  echo "\
#!/bin/bash
#SBATCH --job-name=$af2dir
#SBATCH --output="$af2dir".out
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G
#SBATCH --time=16:00:00" > "$procdir"/"$af2dir"_af2.sh
elif (( "$len" >= 1500 )); then
  echo -e "Found $len residues, setting 36h time limit."
  echo "\
#!/bin/bash
#SBATCH --job-name=$af2dir
#SBATCH --output="$af2dir".out
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100:2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=15G
#SBATCH --time=36:00:00" > "$procdir"/"$af2dir"_af2.sh
fi

echo "\
#SBATCH --mail-type=ALL
#SBATCH --mail-user="$USER"

module load alphafold/"$alphafold_version"
module load pymol/2.3.0

run --fasta_paths="$procdir"/"$af2dir".fa \
	--output_dir="$procdir" \
	--db_preset=full_dbs \
	--num_multimer_predictions_per_model=1 \
	--max_template_date="$db_version" \
	--model_preset=multimer

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
from matplotlib import pyplot as plt
import argparse
import pickle

def get_pae_plddt(model_names):
    out = {}
    for i,name in enumerate(model_names):
        d = pickle.load(open(name,'rb'))
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
    plt.figure(figsize=(14, 4), dpi=100)
    ##################################################################
    plt.subplot(1, 2, 1)
    plt.title(\"Sequence coverage\")
    plt.imshow(final,
               interpolation='nearest', aspect='auto',
               cmap=\"rainbow_r\", vmin=0, vmax=1, origin='lower')
    plt.plot((msa != 21).sum(0), color='black')
    plt.xlim(-0.5, msa.shape[1] - 0.5)
    plt.ylim(-0.5, msa.shape[0] - 0.5)
    plt.colorbar(label=\"Sequence identity to query\", )
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
    plt.savefig(f\"{out_dir}/{name+('_' if name else '')}coverage_LDDT.png\")
    ##################################################################

    ##################################################################
    num_models = 5
    plt.figure(figsize=(3 * num_models, 2), dpi=100)
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        plt.subplot(1, num_models, n + 1)
        plt.title(model_name)
        plt.imshow(value[\"pae\"], label=model_name, cmap=\"bwr\", vmin=0, vmax=30)
        plt.colorbar()
    plt.savefig(f\"{out_dir}/{name+('_' if name else '')}PAE.png\")
    ##################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir',dest='input_dir',required=True)
parser.add_argument('--name',dest='name')
parser.set_defaults(name='')
parser.add_argument('--output_dir',dest='output_dir')
parser.set_defaults(output_dir='')
args = parser.parse_args()

# print(os.listdir(args.input_dir))
feature_dict = pickle.load(open(f'{args.input_dir}/features.pkl','rb'))
is_multimer = ('result_model_1_multimer.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
is_ptm = ('result_model_1_ptm.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
model_names = [f'{args.input_dir}/result_model_{f}{\"_multimer\" if is_multimer else \"_ptm\" if is_ptm else \"\"}.pkl' for f in range(1,6)]

pae_plddt_per_model = get_pae_plddt(model_names)
generate_output_images(feature_dict, args.output_dir if args.output_dir else args.input_dir, args.name, pae_plddt_per_model)
" > "$procdir"/"$af2dir"_vis.py

echo "\
<img src=\"cid:animated.gif\" />
" > "$procdir"/"$af2dir"_mail.htm

echo "\
cp "$procdir"/"$af2dir"_vis.py "$procdir"/"$af2dir"/
cp "$procdir"/"$af2dir"_mail.htm "$procdir"/"$af2dir"/
python3 "$procdir"/"$af2dir"/"$af2dir"_vis.py --input_dir "$procdir"/"$af2dir"/ --name "$af2dir"
tar -C "$procdir"/"$af2dir"/ -czvf "$procdir"/"$af2dir"/"$af2dir"_top_ranked.tar.gz ranked_0.pdb ranked_1.pdb ranked_2.pdb ranked_3.pdb ranked_4.pdb
mutt -e 'set content_type=text/html' -s \"$af2dir\" -a "$procdir"/"$af2dir"/"$af2dir"_top_ranked.tar.gz -a "$procdir"/"$af2dir"/*.png -a "$procdir""$af2dir"/animated.gif -e 'my_hdr From:AlphaFold2 (AlphaFold2)' -b valkove2@nih.gov -- "$USER"@nih.gov < "$procdir"/"$af2dir"/"$af2dir"_mail.htm
rm "$procdir"/"$af2dir"/"$af2dir"_top_ranked.tar.gz
rm "$procdir"/"$af2dir"_mail.htm
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
connect ftp://ftp.box.com;
user '$username' '$password';
cd Alphafold;
mirror -R "$procdir"/"$af2dir";
bye" > $HOME/.lftpconfig
	chmod 600 $HOME/.lftpconfig
        echo "\
lftp -f $HOME/.lftpconfig
rm $HOME/.lftpconfig
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
