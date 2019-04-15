# MATLAB_code
Collection of MATLAB scripts written during PhD work

## Included scripts
* HYPR_recon: Implementation of HYPR MRI recon. Allows for imaging of vascular dynamics with high SNR and good temporal sensitivity.
* characterize_gating2: Gathers rudimentary statistics on gating files from MRI scans and generates histogram of the cardiac gating. Allows for user-guided cleaning of corrupted gating files.
* flow_analysis_2D: Tool allowing user to draw ROI around vessel and receive time-resolved flow measurements through the ROI.
* hypr_phantom: Shepp-Logan phantom with dynamic changes to contrast in selected regions to emulate changes in signal over the cardiac cycle.
* ke_mimics: Imports ventricle mask from MIMICS, applies the mask to 4D flow derived velocity images, and calculates kinetic energy.
* velocity_histogram_TR: Generates time-resolved velocity histograms in region of imported mask to profile flow distributions.
* vorticity: Calculates vorticity in imported mask from corresponding 4D flow images.
