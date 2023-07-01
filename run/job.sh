#!/bin/bash
#SBATCH --job-name=Hypopt_0deg_128x64x64_15000
#SBATCH --mem=18G                    # Memory requested in megabytes. If omitted, the default is 1024 MB.
#SBATCH --time=4-0:0:0      # How long will your job run for? If omitted, the default is 3 hours.
#SBATCH --account=hpc5182
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aidan.sheedy@queensu.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=42
#SBATCH --cpus-per-task=1
 
# commands for your job go here
mpiexec -np $SLURM_NTASKS ./hypopt -maxItr 15000