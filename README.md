# FUSE-Dynamic-Combination-Sim-Paper

This repository allows for the reproduction of results from the paper:

**"Varying the Combination of Hydrological Models in Time and Space: Towards a More Accurate Representation of Streamflow Across Large Domains"**  
*Thébault, C., Knoben, W. J. M., Addor, N., Newman, A.J., Clark, M. P. (2025)*

## Repository Structure

The repository is organized as follows:

- **`00_DATA/`**: Contains the results from individual FUSE models (`-1`, `1` folders), and dynamic combination approach (`HQLQ_XX_XX` folders). The lists of FUSE decisions, parameters for the dynamic combination and catchments are also provided here. Observed streamflow and results from sampling uncertainty are included. The data can be download directly from: http://www.hydroshare.org/resource/cd5cd3116bc544488c223122b4af5516

- **`01_FUSEscripts/`**: Includes scripts to run individual FUSE models. These scripts help set up, execute and evaluate the 78 hydrological models used in the study.

- **`02_DCscripts/`**: Contains scripts to run and evaluate the dynamic combination approach.

- **`03_Visualisation/`**: Includes scripts to generate the figures used in the manuscript, leveraging the data stored in `00_DATA/` (outputs of the previous stages). Include also a folder named `UtilitaryScripts` (script to calculate sampling uncertainty and to analyse temporal evolution of model selection)

- **`99_Figures/`**: Include the figures generated with `04_Visualisation/`.

- **`Shp/`**: Include USA boundaries shapefile (from “North American Atlas - Political Boundaries” (Commission for Environmental Cooperation, 2022)) and gauges locations (from CAMELS dataset (Addor et al, 2017)) used in `04_Visualisation/`. The folder can be download directly from: http://www.hydroshare.org/resource/cd5cd3116bc544488c223122b4af5516

Metrics.R file include various functions to calculate metrics for streamflow evaluation. 

InstallRPackages.R installs the various R packages used in this work.

## Reproducing the Results

To reproduce the results presented in the paper, follow these steps:

1. **Run FUSE models**  
   Execute the scripts in `01_FUSEscripts/` to generate individual model outputs.

2. **Apply dynamic combination approach**  
   The scripts in `02_DCscripts/` implement the time-varying model combination methodology.

3. **Generate some analysis**  
   Run the scripts in `04_Visualisation/UtilitaryScripts/` to generate specific files needed for the visualisation.

4. **Generate visualizations**  
   Run the scripts in `04_Visualisation/` to reproduce the figures presented in the manuscript.

## Requirements

To run the scripts, ensure you have the following dependencies installed:

- R (version 4.3.1 used here)  
- Required R packages: run InstallRPackages.R script
- On a HPC: load the different module needed to run FUSE and to install R packages (e.g. gcc, R, netcdf, proj, gdal, geos, hdf5, openblas -- HPC dependant)

## Citation

If you use this repository or its outputs in your research, please cite:

> Thébault, C, Knoben W. J. M., Addor, N., Newman, A.J., Clark, M. P. (2025). *Varying the Combination of Hydrological Models in Time and Space: Towards a More Accurate Representation of Streamflow Across Large Domains*.
