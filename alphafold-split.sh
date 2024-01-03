#!/bin/bash
#SBATCH --job-name=alphafold
#SBATCH --cpus-per-task=16
#SBATCH --time=8:00:00
#SBATCH --mem=40g

SEQ_FILE="/home/valkove2/fastafiles/hexTTP_NOT1.fa"
OUT_DIR="/mnt/beegfs/valkov/alignments"
AF2_DIR="/mnt/alphafold/2.3.2"

module load alphafold/2.3.2_conda

AF2_DIR="/mnt/alphafold/2.3.2"

python `which run_alphafold_msa.py` \
--run_relax=false \
--use_gpu_relax=false \
--model_preset=multimer \
--db_preset=full_dbs \
--max_template_date=2020-05-14 \
--fasta_paths=$SEQ_FILE \
--output_dir=$OUT_DIR \
--data_dir=$AF2_DIR \
--pdb_seqres_database_path=$AF2_DIR/pdb_seqres/pdb_seqres.txt \
--uniref30_database_path=$AF2_DIR/uniref30/UniRef30_2021_06/UniRef30_2021_06 \
--uniprot_database_path=$AF2_DIR/uniprot/uniprot.fasta \
--uniref90_database_path=$AF2_DIR/uniref90/uniref90.fasta \
--mgnify_database_path=$AF2_DIR/mgnify/mgy_clusters_2022_05.fa \
--template_mmcif_dir=$AF2_DIR/pdb_mmcif/mmcif_files \
--obsolete_pdbs_path=$AF2_DIR/pdb_mmcif/obsolete.dat \
--bfd_database_path=$AF2_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt

INF_RUN="alphafold_inference.sh"

# Paste the content of this script with the common settings into the inference script
awk '/^#!\/bin\/bash$/,/--bfd_database_path=/' $0 > $INF_RUN

# Extend the run time to 24hr
sed -i '/#SBATCH --time=8:00:00/c\#SBATCH --time=24:00:00' $INF_RUN

# Add lines for using the gpu nodes
sed -i '/#SBATCH --mem=40g/a #SBATCH --partition=gpu\n#SBATCH --gres=gpu:v100:2' $INF_RUN

# Change the Python execution script
sed -i 's/python `which run_alphafold_msa.py`/python `which run_alphafold_predict.py`/' $INF_RUN

# Modify the relax and gpu_relax settings to 'true'
sed -i 's/false/true/g' $INF_RUN

# Insert the lines with # of predictions/model and pre-computed msa
sed -i '/run_alphafold_predict.py/a --num_multimer_predictions_per_model=1 \\\n--use_precomputed_msas=true \\' $INF_RUN

# Submit the inference script
sbatch --dependency=afterok:$SLURM_JOB_ID $INF_RUN

# Send an email notification when the job is complete
NAME=$(echo "$SEQ_FILE" | awk -F'/' '{print $NF}')
mutt -s "AlphaFold for $NAME is complete." -- $USER@nih.gov
