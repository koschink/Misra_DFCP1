# DFCP1 omegasome image analysis

Image-analysis and plotting scripts used for the analysis of omegasome dynamics in Misra et al., *Structural basis for the ATP-dependence of omegasome biogenesis by DFCP1*.

## Analysis workflow

1. **`MeasureTracks_Gallery_SaveROIs_Overlay.py` (Fiji/ImageJ, Jython).** Starting from a time-lapse image and manually annotated point ROIs, extracts a centered image gallery for each tracked omegasome. Saves galleries (TIFF), ROI sets and overlays, and a CSV containing fluorescence measurements and track coordinates.
2. **`MeasureOmegasomes_new.py` (Fiji/ImageJ, Jython).** Segments the manually centered omegasome in each gallery frame using thresholding, selects the object in the central region, and exports per-frame area, Feret diameters, aspect ratio and fluorescence intensity in two channels as CSV files.
3. **`Analysis_new_V7.py` (Python 3).** Uses the measurement CSVs to analyze omegasome diameter trajectories, maximum size, lifetime and constriction kinetics. Produces plots and tabular exports.
4. **`IntensityAnalysis_V3.py` (Python 3).** Analyzes time-dependent fluorescence intensities, with per-track baseline subtraction or min–max normalization, and plots intensity trajectories stratified by omegasome size or lifetime.

The objects are **manually tracked and centered** before automated segmentation. The central-object selection is intentional and is not a whole-field object-detection procedure.

## Usage

Run the first two scripts in **Fiji's Jython script environment** and the final two scripts with **Python 3** (using `numpy`, `pandas`, `matplotlib` and `seaborn`). Edit the directory paths, experimental groups, channel assignments, pixel size and frame interval in each script's configuration section to match the dataset.

The scripts are provided as used for the analysis; they contain experiment-specific settings and have not been refactored into a general-purpose package. Original image data and measurement files are not included here.
