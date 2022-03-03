#!/bin/bash
#SBATCH -J CausalOT #A single job name for the array
#SBATCH -e SlurmOutput/ot_%A_%a_out.txt #Standard error
#SBATCH -o SlurmOutput/ot_%A_%a_out.txt #Standard output
#SBATCH -p shared #Partition short, medium long
#SBATCH -t 0-0:45:00 #Running time of 0 day(s), 0 hour(s), 45 minutes
#SBATCH --mem-per-cpu 2GB #Memory request #2000
#SBATCH -n 1 #Number of cores
#SBATCH --mail-type=ALL #mail when start and finish
#SBATCH --mail-user=edunipace@g.harvard.edu #mail to harvard email

module load gcc/9.3.0-fasrc01 R/4.1.0-fasrc01
module load Anaconda3/2020.11
module load libpng/1.6.25-fasrc01

export NOBS=512

mkdir -p TextOutput/SINGLE_ARRAY/${NOBS}/${SLURM_ARRAY_JOB_ID}

source activate COT # Conda environment

echo $PYTHON_PATH

export PYTHON_PATH="~/.conda/envs/COT/bin/python"

echo $PYTHON_PATH

echo "SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"

export ARRAY_SET=1
R CMD BATCH --no-save hainmueller_sim.R TextOutput/SINGLE_ARRAY/${NOBS}/${SLURM_ARRAY_JOB_ID}/OT_SINGLE_ARRAY_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${ARRAY_SET}.txt

#done
conda deactivate

