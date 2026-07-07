#!/bin/bash
#SBATCH -J motif_enum
#SBATCH -o motif_enum_%j.out
#SBATCH -e motif_enum_%j.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ashley.scruse@morehouse.edu

module load python3

source ~/envs/subnetwork-motifs/bin/activate

cd $SLURM_SUBMIT_DIR/..

python -m src.significance.enumerate \
    --graph data/processed/yeast_graph.pt \
    --pi results/predicted_pi.csv \
    --max-k 4 \
    --output results/motif_catalog.csv
