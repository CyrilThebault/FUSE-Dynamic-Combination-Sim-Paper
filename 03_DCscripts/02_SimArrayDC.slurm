#!/bin/bash
#SBATCH --job-name=SimArrayDC
#SBATCH --partition=cpu2022,cpu2021,cpu2019
#SBATCH --output=/home/cyril.thebault/Out/Bassin-%A-%a.out
#SBATCH --error=/home/cyril.thebault/Err/Bassin-%A-%a.error
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --array=1


module load R/4.3.1
module load osgeo/gdal
module load lib/netcdf


Rscript 02_SimArrayDC.R
