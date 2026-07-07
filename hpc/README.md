# HPC Setup (TACC Frontera)

## Access

The Morehouse Supercomputing Facility (MSCF) has an allocation on TACC. Use your TACC credentials to log in.

```bash
ssh username@frontera.tacc.utexas.edu
```

## Environment Setup

```bash
# Load modules
module load python3
module load cuda

# Create virtual environment
python3 -m venv ~/envs/subnetwork-motifs
source ~/envs/subnetwork-motifs/bin/activate

# Install dependencies
pip install torch torch-geometric networkx pandas numpy scipy scikit-learn matplotlib optuna
```

## Job Submission

TACC uses SLURM for job scheduling.

```bash
# Submit a job
sbatch train_gat.slurm

# Check job status
squeue -u $USER

# Cancel a job
scancel <job_id>
```

## Container Notes

TACC uses Apptainer (not Docker). If a containerized pipeline is needed:

```bash
apptainer pull docker://pytorch/pytorch:latest
apptainer exec pytorch_latest.sif python train.py
```

## GPU Nodes

GAT training requires GPU. Request GPU nodes with:
```
#SBATCH -p rtx        # RTX queue (Frontera)
#SBATCH -N 1
#SBATCH -n 1
```

Motif enumeration is CPU-bound. Use the normal compute queue.
