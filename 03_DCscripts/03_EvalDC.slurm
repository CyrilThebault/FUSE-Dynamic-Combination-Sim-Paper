#!/bin/bash
#SBATCH --job-name=EvalDC
#SBATCH --partition=cpu2022,cpu2021,cpu2019
#SBATCH --output=/home/cyril.thebault/Out/Bassin-%A-%a.out
#SBATCH --error=/home/cyril.thebault/Err/Bassin-%A-%a.error
#SBATCH --mem=5G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-01:00:00
#SBATCH --array=1-559


module load R
module load osgeo/gdal
module load lib/netcdf

## Instruction BASH permettant de récupérer la ligne correspondante au numéro d’ARRAY correspondant
input=${SLURM_ARRAY_TASK_ID}
## On lance le programme R avec en  argument le chemin vers le fichier correspondant au bassin SLURM_ARRAY_TASK_ID
Rscript 03_EvalDC ${input}
