#!/bin/sh

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=63
#SBATCH --mem=240G
#SBATCH -o log/%x_%j.out
#SBATCH --mail-user=frederic.ouimet.23@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=message/%x_%j-result.txt
#SBATCH --error=message/%x_%j-error.txt


## Load the R environment
module load StdEnv/2023 r/4.3.1

## Run the R script directly
Rscript Wishart_density_estimation.R
