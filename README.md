# BMPintegral
code for the BMP signaling histories project
Most analysis and figure were produced with MATLAB version 2019b or newer

Code for preprocessing, segmenting, and quantifying live and fixed-cell image data as described in the manuscript is in the folder "image_analysis"
Code for linking live to fixed-cell data and consolidating into a single data structure is in the subfolder "live_fixed_analysis"
Code for automated single-cell tracking and manual validation and correction of automated tracking results is in the "tracking" subfolder
"external" contains code from outside sources, e.g., publicly available GitHub repositories and MATLAB central File Exchange
"helpers" contains commonly used functions primarily for making, formatting, and saving figures
Code and associated data to make most of the figures in the manuscript is contained in the "figures" folder
To add the code in this repository and its subfolders to your MATLAB path, run the 'setup.m' script from the top directory of the repository.

Within the figures folder, subfolders contain code and data for specific figures and figure panels.
Files in "bulkRNAseq" produce panels fig 5EFG and SI fig 5EF
Files in "level_duration_integral_thresholds" reproduce fig 4 A-C,E-H and SI fig 4 D
Files in "radial_history_analysis" reproduce fig 1E-N, SI fig 1E,G,P-S
Files in single_cell_history_analysis reproduce fig 3 and SI fig 3 panels
Files in SOX2_SMAD4_model reproduce fig 5 HIJ, SI fig 5GH
Files in "time_series_IF" reproduce fig 5 ABC, SI fig 5ACD