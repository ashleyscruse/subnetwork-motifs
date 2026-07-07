#!/bin/bash
#SBATCH -J gat_train
#SBATCH -o gat_train_%j.out
#SBATCH -e gat_train_%j.err
#SBATCH -p rtx
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ashley.scruse@morehouse.edu

module load python3
module load cuda

source ~/envs/subnetwork-motifs/bin/activate

cd $SLURM_SUBMIT_DIR/..

python -m src.training.train \
    --data data/processed/yeast_graph.pt \
    --epochs 200 \
    --lr 0.001 \
    --folds 5 \
    --output results/
