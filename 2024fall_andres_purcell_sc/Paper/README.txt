# README

## Overview
This repository contains all the necessary scripts and data to reproduce the results presented in the manuscript and its supplementary information (SI). The organization of files and instructions are outlined below to ensure smooth usage and reproducibility.
---

## File Structure
### 1. **Coding Scripts**
- **Decay_rates/**  
  Contains scripts to calculate decay rates. Adjust the parameters as needed and submit jobs to generate the desired data.

- **Energy_relaxation/**  
  Contains scripts for simulating energy relaxation processes. Adjust parameters as needed before submitting the jobs.

- **Spectrum/**  
  Contains scripts to reproduce spectral data. Customize parameters as required and submit the jobs.
---

### 2. **Plotting Scripts**
- **Plotting_scripts/**  
  - Includes all plotting scripts to reproduce figures and subplots from the manuscript and SI.  
  - The manuscript and SI data are organized by their corresponding figure and subplot labels.  
  - Ensure the file paths in the scripts are updated to match the directory structure on your local machine.
---

## Usage Instructions
1. **Reproducing Data:**
   - Navigate to the relevant directory (**Decay_rates**, **Energy_relaxation**, or **Spectrum**).
   - Modify the parameters in the scripts according to your requirements.
   - Submit the jobs using your preferred job submission system (e.g., SLURM, PBS, or another scheduler).

2. **Generating Plots:**
   - Navigate to the **Plotting_scripts/** directory.
   - Update the file paths in the scripts to correspond with your local file system.
   - Run the scripts to generate the figures and subplots as presented in the manuscript.
---

## Notes
- Ensure all required dependencies and libraries are installed before running the scripts.
- Consult the manuscript or SI for details about the data and plots if needed.
- For any issues or questions, refer to the contact information provided in the manuscript.
---

**Happy coding!**

