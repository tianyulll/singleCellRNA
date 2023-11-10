#!/bin/bash
#SBATCH --error=log/int_%j.err.txt
#SBATCH --output=log/int_%j.txt
#SBATCH --mem=320G
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16


module load R

Rscript /home/lut2/project/singleCellRNA/integrat_test.r