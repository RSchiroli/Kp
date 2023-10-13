# Kp

Code for automated computation of partition coefficient in coacervates and analysis.

# GENERAL

Julia 1.8.5 or later version is needed to run the code.
The code can be ran via Jupyter Notebook, however it's not suggested to do that since the large amount of RAM required.
Run it directly from terminal allow quicker analysis.


General comment on acquisition process and data needed for the analysis:
 -  the code works analyzing a choosen number 'nPlane' around the maximum intensity one. So it's in general needed a Zstack with nPlane/2 planes above and below the maximum intensity.
 -  nPlane>1 leads in some cases to errors, it is suggested to work with only one plane, it means that the intensities of dilute and dense phase will be extracted for any plane but the labelling of single droplets and Kp computed only in the maximum intensity one.
 -  a mask is required to run the analysis. To acquire the mask Fiji's plugin LabKit is needed. Alternatively the code could be modified to deal with manual intensity threshold (see below).

# MASK ACQUISITION

To construct the mask for the image open Fiji, then from the menu plugin->LabKit->open image with LabKit. It's important to train the classifier with the right order in the label sequence, since it determines the values of pixels in that class.
- label 1 (background) = dilute phase
- label 2 (forecast) = dense phase
- label 3 = blur
  
Then Segmentation->train classifier, Segmentation->Save segmentation results as TIF. Segmentation->Save classifier if you want to export the classifier. You can then apply it to other images with similar features without performing the training. To apply it you just need to open a new image in LabKit, then Segmentation->open classifier.

# COMPUTE KP

Before start you need to include the code. Open the terminal in the code directory. Then open julia from the terminal and include the files:

- include("Kp.jl") is always fine.
- include("BockStats.jl") if only the analysis is needed.

Function _RunKpAnalysis(imgName,maskName;Gain::Int64=100,threshold::Int64=0,nPlane::Int64=0,returnAll::Bool=false)_ from Kp.jl

  -  imgName::String = "full_path_to_your_image_in_string_format/image_name.tif"
  -  maskName::String = "full_path_to_your_mask_in_string_format/mask_name.tif"

optional parameters (they are called with VarName=VarValue):
  -  Gain::Int64 = gain used in the aquisition (it's relevant only if you are interested in restore the original intensity on your image, for Kp computation it's irrelevant)
  -  threshold::Int64 = threshold number of pixel in a droplet. Droplets below the threshold won't be analysed
  -  nPlane::Int64 = number of planes around the maximum intensity plane. Better to use the default values (nPlane=0)
  -  returnAll::Bool = if equal to True the function returns avgKp, stdKp, allKp, Sizes, label. If false it only save them in a .csv file.

If the code is called with imgName = "path/name.tif", results will be saved in "path/name.csv".

# ANALYSIS

function _BlockStats(CSVpathName)_ from BockStats.jl

Function to extract Kp from a single ZStack. Blocks are constructed dividing droplets by size. It takes as an input the full path to CSV file produced by _RunKpAnalysis_ and returns two arrays "results, stats" where:
results = avgKpBlock, StdKpBlock, AvgSizeBlock, StdSizeBlock
stats = intensityDilutePhaseInThePlane, NumberOfDropletsInThePlane, NumberOfPixels in the dilute phase


function _MultipleImagesBlockStats(CSVpathName,Title)_ from BockStats.jl

Function to extract Kp from multiple ZStack (of the same system). Blocks are constructed building an array with all the droplets from different ZStacks and dividing them by size. It takes as an input the full path to CSV file produced by _RunKpAnalysis_ and title should be used in saved plots. It returns AvgKpBlocks, StdKpBlocks, AvgSizeBlocks, StdSizeBlocks containing statistics for any block. In addition it creat a folder "path/Title" where "path" is the path to the first CSV passed as an argument and save in it the following list of plots (y vs x):
  - avgKp vs avgSize in the block;
  - avgKp vs avgSize in the block dividing results for image;
  - avgSize in the block n vs n;
  - number of droplets vs avgSize dropets;
  - dilute Intensity vs number of pixel in the dilute intensity;
  - dilute intensity vs avgKp.


function _MultipleImagesStats(CSVpathName,Title)_ from BockStats.jl

Function to extract Kp from a single ZStack. Blocks are constructed creating an array with all the droplets in different ZStacks and dividing them by size. It takes as an input the full path to CSV file produced by _RunKpAnalysis_ and saves results in a CSV file located in the path to the first file passed as an argument.
In any line of the file the are AvgKp, StdKp, AvgSize, StdSize of any image analysed, a line of 0.0 and then results averaged over all the images.

# MANUAL INTENSITY THRESHOLD

To use manual threshold on the intensity the easiest way is to pass as an argument of _RunKpAnalysis_ the path ro raw data both for imgName and maskName.
Then it is sufficient to modify the function _Segment!_ from ComputeIntensity.jl to classify dilute, dense and blur pixels based on the manualt threshold. A better version can be implemented avoiding call twice the same image.
