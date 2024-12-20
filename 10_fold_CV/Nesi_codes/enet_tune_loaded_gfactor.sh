#!/bin/bash -e
#SBATCH --account=uoo03493
#SBATCH --job-name    enet
#SBATCH --cpus-per-task 8	
#SBATCH --mem 51200MB
#SBATCH --time        02:00:00
#SBATCH --output=/nesi/nobackup/uoo03493/Yue/stacking_gfactor_10foldcv/slurm_out/enet.%J.out # Include the job ID in the names of
#SBATCH --error=/nesi/nobackup/uoo03493/Yue/stacking_gfactor_10foldcv/slurm_error/enet.%J.err # the output and error files

module load R/4.3.2-foss-2023a

DATA=$1

# Help R to flush errors and show overall job progress by printing
# "executing" and "finished" statements.
echo "Executing R ..."
srun Rscript enet_tune_loaded_gfactor.R $DATA
echo "R finished."