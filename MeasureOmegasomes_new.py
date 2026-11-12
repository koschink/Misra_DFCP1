from ij import IJ
from ij import *
from ij import ImagePlus
from ij import ImageStack
from ij.measure import *
from ij.plugin import *
from ij.process import *
from ij.plugin import ChannelSplitter
from ij.measure import Measurements
from ij.measure import ResultsTable
from ij.plugin import ImageCalculator
from ij.plugin import RGBStackConverter
from ij.plugin import StackCombiner
from ij.plugin import Duplicator
from ij.plugin.frame import RoiManager
from ij.process import ImageProcessor
from ij.process import ImageStatistics
from ij.process import EllipseFitter
from ij.gui import Roi
from ij.gui import PointRoi
from ij.gui import OvalRoi
from ij.gui import Overlay
import math
import time
from java.awt import *
from java.awt import Font
import itertools
import glob
import os

# -------------------------------------------------------------------
# USER CONFIG
# -------------------------------------------------------------------

# Root path with subfolders <directory>/<cell_line>/Galleries/*.tif
datapath = r"C:\DFCP1\Analysis"

# Example directories and strains – edit as needed
dirnames = ["251016", "251024", "251028"]
cell_lines = ["243", "244", "245", "246"]

# ROI size (currently unused in this script, but kept for compatibility)
roisize = 10

# Pixel size and frame rate (not used in this script, but often useful)
Pixelsize = 0.080
framerate = 0.33

# DEFAULT logical channel mapping (ImageJ is 1-based)
# "primary" = channel used for segmentation and treated as Channel 1 in CSV
# "secondary" = other marker, treated as Channel 2 in CSV
DEFAULT_PRIMARY_CHANNEL = 2
DEFAULT_SECONDARY_CHANNEL = 1

# Per-directory RAW channel mapping:
#   - directory_primary_channel_map[directory]   = raw channel index of primary channel
#   - directory_secondary_channel_map[directory] = raw channel index of secondary channel
#
# Adjust these if acquisition order differs between dates.
directory_primary_channel_map = {
    "251016": 2,  # primary marker in channel 2 for 251016
    "251024": 1,  # primary marker in channel 1 for 251024
    "251028": 1,  # primary marker in channel 2 for 251028

    
}

directory_secondary_channel_map = {
    "251016": 1,  # secondary marker in channel 1 for 251016
    "251024": 2,  # secondary marker in channel 2 for 251024
    "251028": 2,  # secondary marker in channel 2 for 251028

}

# Default marker names if a directory is not listed
DEFAULT_PRIMARY_NAME = "DFCP1"
DEFAULT_SECONDARY_NAME = "OtherMarker"

# Per-directory marker names for primary / secondary channels
# Edit to reflect your constructs
directory_marker_name_map = {
    "251016": {"primary": "mNG-DFCP1", "secondary": "LC3"},
    "251024": {"primary": "mNG-DFCP1", "secondary": "LC3"},
    "251028": {"primary": "mNG-DFCP1", "secondary": "LC3"},

}

# -------------------------------------------------------------------
# HELPER FUNCTIONS
# -------------------------------------------------------------------

def normList2(L, normalizeTo=100):
    """normalize values of a list to make its max = normalizeTo"""
    vMax = max(L)
    return [x / (vMax * 1.0) * normalizeTo for x in L]

def CalcVelocity(distancelist, framerate):
    return [x * framerate for x in distancelist]

def GenerateRGB(imp):
    channels1 = imp.getNChannels()
    slices1 = slices2 = imp.getNSlices()
    frames1 = frames2 = imp.getNFrames()
    c1 = imp.getChannel()
    z1 = imp.getSlice()
    t2 = imp.getFrame()
    imp2 = imp.createHyperStack("title2", 1, slices2, frames2, 24)
    RGBStackConverter().convertHyperstack(imp, imp2)
    return imp2

def CombineStacks(imp1, imp2):
    imp1stack = imp1.getStack()
    imp2stack = imp2.getStack()
    imp3stack = StackCombiner().combineVertically(imp1stack, imp2stack)
    imp3 = ImagePlus("Combined Stacks", imp3stack)
    print imp3
    return imp3

def ExtractChannel(imp, channel):
    """Return a duplicate of the selected channel from a hyperstack."""
    slicenumber = imp.getNSlices()
    timepoints = imp.getNFrames()
    extracted = Duplicator().run(imp, channel, channel,
                                 1, slicenumber,
                                 1, timepoints)
    return extracted

def dist(p0, p1):
    """Calculates the distance between two xy coordinates."""
    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)

def Generate_segmented_manual(imp, channel, threshold_method, filtersize):
    """Manual thresholding helper (unused in automatic pipeline)."""
    imp_Threshold_1 = ExtractChannel(imp, channel)
    imp_Threshold_1.setTitle("ThresholdChannel" + str(channel))
    imp_Threshold_1.show()
    IJ.run(imp_Threshold_1, "Median...", "radius=" + str(filtersize))
    IJ.run(imp_Threshold_1, "Subtract Background...", "rolling=50")
    IJ.run("Threshold...")
    WaitForUserDialog("Threshold", "set a threshold").show()
    thres_min = imp_Threshold_1.getProcessor().getMinThreshold()
    thres_max = imp_Threshold_1.getProcessor().getMaxThreshold()
    IJ.setThreshold(imp_Threshold_1, thres_min, thres_max)
    IJ.run(imp_Threshold_1, "Convert to Mask", "")
    return imp_Threshold_1

def ThresholdImage(imp, threshold_method, filtersize):
    IJ.run(imp, "Median...", "sigma=" + str(filtersize))
    IJ.setAutoThreshold(imp, threshold_method + " dark stack")
    IJ.run(imp, "Convert to Mask",
           "method=" + threshold_method + " background=Dark black")
    return imp

# -------------------------------------------------------------------
# MAIN LOOP
# -------------------------------------------------------------------

for dirnum, directory in enumerate(dirnames):
    print(directory)
    datapath1 = datapath + "/" + directory
    print(datapath1)

    # Per-directory primary/secondary mapping
    primary_channel = directory_primary_channel_map.get(
        directory, DEFAULT_PRIMARY_CHANNEL
    )
    secondary_channel = directory_secondary_channel_map.get(
        directory, DEFAULT_SECONDARY_CHANNEL
    )

    marker_names = directory_marker_name_map.get(directory, {})
    primary_name = marker_names.get("primary", DEFAULT_PRIMARY_NAME)
    secondary_name = marker_names.get("secondary", DEFAULT_SECONDARY_NAME)

    print("  Using PRIMARY channel (logical Channel 1):", primary_channel)
    print("  Using SECONDARY channel (logical Channel 2):", secondary_channel)
    print("  PRIMARY marker name:", primary_name)
    print("  SECONDARY marker name:", secondary_name)

    for cellline_num, cell_line in enumerate(cell_lines):
        savepath = datapath + "/" + directory + "/" + cell_line + "/Measurements_2/"
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        readpath = datapath + "/" + directory + "/" + cell_line + "/Galleries/"
        readpath = readpath + "*.tif"
        print(readpath)
        files1 = glob.glob(readpath)

        print(files1)  # these files will be processed

        for imagefile in files1:
            filename = os.path.normpath(imagefile)
            imp = IJ.openImage(filename)
            imp.show()
            imp = IJ.getImage()

            # Use PRIMARY channel for segmentation
            imp_seg = ExtractChannel(imp, primary_channel)
            imp_seg.show()
            ThresholdImage(imp_seg, "Otsu", 1)

            IJ.run(imp_seg, "Dilate", "stack")
            # IJ.run(imp_seg, "Erode", "stack")

            # Add particles to ROI Manager
            IJ.run(imp_seg, "Analyze Particles...",
                   "size=3-Infinity pixel exclude add stack")

            rm = RoiManager.getInstance()
            if rm is None:
                rm = RoiManager()

            # Filter ROIs to central region
            roicount = rm.getCount()
            added_rois = 0
            for roi_index in range(roicount):
                rm.select(imp_seg, roi_index)
                roi1 = imp_seg.getRoi()
                bounds = roi1.getBounds()
                XCoord = bounds.x + (bounds.width / 2)
                YCoord = bounds.y + (bounds.height / 2)

                # Hard-coded central window from 20–30 px in x and y
                if XCoord < 20:
                    continue
                if XCoord > 30:
                    continue
                if YCoord < 20:
                    continue
                if YCoord > 30:
                    continue

                IJ.run(imp_seg, "Add Selection...", "")
                added_rois += 1

            rm.runCommand(imp, "Deselect")
            rm.runCommand("Delete")

            if added_rois > 0:
                # transfer filtered ROIs to ROI Manager linked to imp
                IJ.run("To ROI Manager", "")

                movietime = []
                MeanPrimary = []
                MeanSecondary = []
                Area = []
                max_diameter = []
                min_diameter = []

                rois_to_measure = rm.getCount()

                for measure_roi in range(rois_to_measure):
                    rm.select(imp, measure_roi)
                    roi2 = imp.getRoi()

                    movietime.append(imp.getT())

                    # Measure geometry + PRIMARY intensity
                    imp.setC(primary_channel)
                    imp.setRoi(roi2)
                    feret = roi2.getFeretValues()
                    stats = imp.getStatistics(Measurements.MEAN |
                                              Measurements.AREA |
                                              Measurements.FERET |
                                              Measurements.CENTROID)
                    max_diameter.append(feret[0])
                    min_diameter.append(feret[2])
                    Area.append(stats.area)
                    MeanPrimary.append(stats.mean)

                    # Measure SECONDARY channel intensity
                    imp.setC(secondary_channel)
                    stats2 = imp.getStatistics(Measurements.MEAN |
                                               Measurements.AREA |
                                               Measurements.FERET |
                                               Measurements.CENTROID)
                    MeanSecondary.append(stats2.mean)

                # Normalize (logical Channel 1 is PRIMARY)
                NormPrimary = normList2(MeanPrimary)
                NormSecondary = normList2(MeanSecondary)

                ort = ResultsTable()
                ort.setPrecision(3)

                count = len(MeanPrimary)
                dataname = imp.getShortTitle() + ".csv"

                for i in range(count):
                    ort.incrementCounter()
                    ort.addValue("Frame", i)
                    ort.addValue("Area", Area[i])
                    ort.addValue("Max Diameter", max_diameter[i])
                    ort.addValue("Min Diameter", min_diameter[i])
                    ort.addValue("Aspect ratio",
                                 (max_diameter[i] / min_diameter[i]))

                    # Logical channels:
                    # Channel 1 / NormCh 1 = PRIMARY
                    # Channel 2 / NormCh 2 = SECONDARY
                    ort.addValue("Channel 1", MeanPrimary[i])
                    ort.addValue("Channel 2", MeanSecondary[i])
                    ort.addValue("NormCh 1", NormPrimary[i])
                    ort.addValue("NormCh 2", NormSecondary[i])

                    # Marker identity and raw channel indices
                    ort.addValue("Chan1_marker", primary_name)
                    ort.addValue("Chan2_marker", secondary_name)
                    ort.addValue("Chan1_raw", primary_channel)
                    ort.addValue("Chan2_raw", secondary_channel)

                    ort.addValue("Timepoint", str(movietime[i]))
                    ort.addValue("Strain", str(cell_line))
                    ort.addValue("Filename", str(filename))
                    ort.addValue("Dataname", str(dataname))

                ort.show("Measured intensities")
                print dataname
                ort.saveAs(savepath + "/" + dataname)

                rm.runCommand(imp, "Deselect")
                rm.runCommand("Delete")

            imp.changes = False
            imp_seg.changes = False
            IJ.run("Close All", "")