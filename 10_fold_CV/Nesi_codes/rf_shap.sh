#!/bin/bash -e
#SBATCH --account=uoo03493
#SBATCH --job-name    random_forest_shap
#SBATCH --partition=milan
#SBATCH --cpus-per-task 20	
#SBATCH --mem 153600MB
#SBATCH --time        15:00:00
#SBATCH --output=/nesi/nobackup/uoo03493/Yue/stacking_gfactor_modelling/slurm_out/rf_shap.%J.out # Include the job ID in the names of
#SBATCH --error=/nesi/nobackup/uoo03493/Yue/stacking_gfactor_modelling/slurm_error/rf_shap.%J.err # the output and error files

module load R/4.3.2-foss-2023a

DATA=$1

# Help R to flush errors and show overall job progress by printing
# "executing" and "finished" statements.
echo "Executing R ..."
srun Rscript rf_shap.R $DATA
echo "R finished."