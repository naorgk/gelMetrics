# Instructions for Use
## *GelMetrics: An algorithm for analyzing the dynamics of gel-like phase separated condensates*  

### This repository contains all the code and data required to recreate the figures from the paper.
>**Tested on Matlab R2025b**

### Prerequisites
- Curve Fitting Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- In addition, ImageJ FIJI (specifically, the Mosaic plugin) was used to read microscopy image files, detect bright spots, and convert to the formats required for the scripts.
### Modes
- To run on simulated data, run the script: `SimulationRunWithAnalysis.m`  
  You can switch between simulation types using the `simType` variable defined within the script.

- To run on experimental data, run the script: `ExperimentalAnalysis.m`  
  You can switch between different data sets by un/commenting the relevant signals and backgrounds variables which load the signal files

- New data pipeline:
  1. Obtain microscopy timelapse images. From our experience, 8-12 images (each containing 360 timepoints) is mostly sufficient.
  2. In ImageJ FIJI, use the Particle Tracker 2D/3D (included in Mosaic plugin) to generate CSV files of spatiotemporal granule coordinates. We provide a macro file (`BatchParticleTracker.ijm`) written to be run through ImageJ FIJI: `Plugins` &#8594; `Macros` &#8594; `Run...`  
      > For a detailed explanation about this part of the analysis, see our previous work in: https://github.com/naorgk/slncRNA_Analysis/blob/main/Full%20pipeline%20walkthrough.pdf

  3. Run script `read_tables.m` on CSV files generated in step II. The script should be in the same folder as CSV files. The script expects the specific structure of the Mosaic Particle Tracker.
  4. Run script `Collect_intensity_data.m` on microscopy images and `.mat` files generated in step III. The script assumes the same alphabetical order for both `.tiff` and `.mat` files. The script will generate two `.mat` files corresponding to data & background.
  5. Run script `ExperimentalAnalysis.m` on the `.mat` files generated in step IV.
