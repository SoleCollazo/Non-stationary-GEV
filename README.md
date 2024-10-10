# Non-stationary-GEV in R

### Overview
This repository contains an R script, NSGEV_rio.R, which performs an analysis using the Non-Stationary Generalized Extreme Value (NSGEV) model. The script is designed for environmental or climate-related data, where extremes (such as rainfall or temperature) are modeled and analyzed.

### Features
* Non-Stationary Generalized Extreme Value (NSGEV) model implementation.
* Works with environmental data such as climate or hydrological data.
* Visualizes the results and trends over time.

### Prerequisites
Before running the script, make sure you have the following dependencies installed in your R environment:

install.packages(c("ROOPSD", "lubridate", "dplyr", "ggplot2", "trend", "readxl", "scales", "pals", "ismev"))

### Usage
1) Clone the repository:

git clone https://github.com/SoleCollazo/Non-stationary-GEV.git

cd Non-stationary-GEV

2) Prepare your data:

Ensure your data is formatted correctly before using the script. The script expects time series data for environmental variables. Modify the script as needed to load your dataset.

3) Run the script:

Open the script in RStudio or run it in your R console:

source("NSGEV_rio.R")

### Output
* Plots: Visualizations of the extreme value trends over time.
* Model Summary: Coefficients and statistics of the fitted NSGEV model.
