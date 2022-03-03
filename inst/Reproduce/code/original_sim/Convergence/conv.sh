#!/bin/bash
#SBATCH -J ConvOT #A single job name for the array
#SBATCH -e SlurmOutput/conv_%A_%a_out.txt #Standard error
#SBATCH -o SlurmOutput/conv_%A_%a_out.txt #Standard output
#SBATCH -p shared #Partition serial_requeue, shared
#SBATCH -t 0-20:00:00 #Running time of 0 day(s), 20 hours
#SBATCH --mem-per-cpu 5GB #Memory request #12000
#SBATCH -n 1 #Number of cores
#SBATCH --mail-type=ALL #mail when start and finish
#SBATCH --mail-user=edunipace@mail.harvard.edu #mail to harvard email

module load gcc/9.3.0-fasrc01 R/4.1.0-fasrc01
module load Anaconda3/2020.11
module load libpng/1.6.25-fasrc01

mkdir -p TextOutput/conv/${SLURM_ARRAY_JOB_ID}
mkdir -p TextOutput/pen/${SLURM_ARRAY_JOB_ID}

source activate COT

R CMD BATCH --no-save convergence_sim.R TextOutput/pen/${SLURM_ARRAY_JOB_ID}/conv_pen_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

conda deactivate
