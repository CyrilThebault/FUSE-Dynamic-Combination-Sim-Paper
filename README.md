# FUSE-Dynamic-Combination-Sim-Paper

This repository allows for the reproduction of results from the paper:

**"A Time-Varying Combination of Hydrological Models: Towards a Better Representation of Streamflow Across Large Domains"**  
*Thébault, C, Knoben W. J. M., Addor, N., Clark, M. P. (2025)*

## Repository Structure

The repository is organized as follows:

- **`00_DATA/`**: Contains the results from individual FUSE models (-1, 1, Comp folders), multi-model mosaic (SamplingUncertainty folder) and dynamic combination approach (WA folder). The lists of FUSE decisions, parameters for the dynamic combination and catchments are also provided here. The data can be download directly from: http://www.hydroshare.org/resource/86bcf167b8ff433baaf3be892f74394a

- **`01_FUSEscripts/`**: Includes scripts to run individual FUSE models. These scripts help set up, execute and evaluate the 78 hydrological models used in the study.

- **`02_MMMscripts/`**: Provides scripts to reproduce the multi-model mosaic, which was derived from the paper "Technical note: How many models do we need to simulate hydrologic processes across large geographical domains?" (Knoben et al., 2025).

- **`03_DCscripts/`**: Contains scripts to run and evaluate the dynamic combination approach.

- **`04_Visualisation/`**: Includes scripts to generate the figures used in the manuscript, leveraging the data stored in `00_DATA/` (outputs of the previous stages).

- **`99_Figures/`**: Include the figures generated with `04_Visualisation/`.

- **`Shp/`**: Include USA boundaries shapefile (from “North American Atlas - Political Boundaries” (Commission for Environmental Cooperation, 2022)) and gauges locations (from CAMELS dataset (Addor et al, 2017)) used in `04_Visualisation/`. The folder can be download directly from: http://www.hydroshare.org/resource/86bcf167b8ff433baaf3be892f74394a

Metrics.R file include various functions to calculate metrics for streamflow evaluation. 

InstallRPackages.R installs the various R packages used in this work.

## Reproducing the Results

To reproduce the results presented in the paper, follow these steps:

1. **Run FUSE models**  
   Execute the scripts in `01_FUSEscripts/` to generate individual model outputs.

2. **Perform sampling uncertainty analysis to design the multi-model mosaic**  
   Use the scripts in `02_MMMscripts/` to analyze sampling uncertainty and design the multi-model mosaic.

3. **Apply dynamic combination approach**  
   The scripts in `03_DCscripts/` implement the time-varying model combination methodology.

4. **Generate visualizations**  
   Run the scripts in `04_Visualisation/` to reproduce the figures presented in the manuscript.

## Requirements

To run the scripts, ensure you have the following dependencies installed:

- R (version 4.3.1 used here)  
- Required R packages: run InstallRPackages.R script
- On a HPC: load the different module needed to run FUSE and to install R packages (e.g. gcc, R, netcdf, proj, gdal, geos, hdf5, openblas -- HPC dependant)

## Citation

If you use this repository or its outputs in your research, please cite:

> Thébault, C, Knoben W. J. M., Addor, N., Clark, M. P. (2025). *A Time-Varying Combination of Hydrological Models: Towards a Better Representation of Streamflow Across Large Domains*.
