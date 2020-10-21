# Spine imaging
This repository is accompanying the paper:  
***Functional clustering of dendritic activity during decision-making***   
Kerlin, A., Mohar, B., Flickinger, D., MacLennan, B.J., Davis, C., Spruston, N., and Svoboda, K. (2018). [BioRxiv 440396](https://www.biorxiv.org/content/early/2018/10/10/440396)

## Visualization (SpineVis)
The accompanying website [SpineVis](spinevis.janelia.org) has all the sessions included in the paper.  
The paper's supplement shows how to find the underlying session and data for each figure.
See the documentation for the website's API [here](https://github.com/boazmohar/spinevis/blob/master/docs/API.md)


## Analysis code
The python package we used is called `prep` within it are several modules:
1. Session.py - Is the main class that handles reading data from [ScanImage](https://vidriotechnologies.com/scanimage/) 
2. SpineSession.py - a derived class from Session.py that knows how to handle [mROI](http://scanimage.vidriotechnologies.com/display/SI2018/Multiple+Region+of+Interest+%28MROI%29+Imaging) data from ScanImage alongside a `Sp.mat` file that defines the trajectory and other metadata of a spine imaging session.
3. Registration.py - Has registration related functions
4. Timecourses.py - Has code for extraction of time courses after registration and tracing. Includes functions to measure distances among spines and dendrites.
5. Embedding.py - Has code to embed  small 2d fields back to 3d space for visualization and tracing.
6. IO.py - Has read and write related functions
7. Utils.py - Has misc. helper functions
8. Step.py, Steps.py, Pipeline.py - Preprocessing related functions

## Code requirements

### Spark
To run this code, you would need a [Spark cluster](https://spark.apache.org/) running to preform parallel computing (tested with version 2.1 and 2.2, should work with a newer versions as well). It will be possible to run it on a single computer [(See here)](https://spark.apache.org/docs/latest/spark-standalone.html), but it would take a very long time to process	 a single session. We routinely use 100-300 cores to run registration, time course extraction and other analysis.

### Python

This code was tested with both python 2.7 and 3.6. See the dependencies in the `requirements` folder with the `Conda_*.txt` flies for using Anaconda or a `requirements.txt` for pip.
### Folder structure
The code was run with some assumption as to the location of the raw / processed data.

#### Raw data
We assume the following structure of the raw data: ***[base_folder]/[date]_[animal]/[Run#]***
For example: ***/Users/username/Data/151123_BMWR30/Run1***
Inside would be:

 - `*.tif` files each representing a behavioral trial.
Each file would be 3d array of images that is:
	1. (z) #planes per volume x  #volumes per trial
	2. x pixels
	3. y pixels
 - `Sp.mat` file has information from our trajectory design GUI about the location and timing of each imaging field (see the folder ***Matlab\ Segmentation***).
 - `behavior.mat` file with information regarding the type of the behavioral trial, when and where the animal licked etc.

#### Procced data
After registration, a few intimidate files are saved.
1. Pre-processed version of the raw data (CleanBinary)
2. Registered version of the data (RegBInary)
3. Metadata about the session and parameters used (Step_x.p)
4. Registration intermediates  in folders: `Run1Intermediates`, `Run1View`

At this point these have defaults that are right for our cluster folder structure and could be changed.

After tracing of the registered session (see the ***Matlab\workflow.m*** and ***Matlab\DendSegmentation\***) is complete a few more files are made in a different location based on the following  structure:
***[base_folder]/Database/[Animal]/[FOV]/[Date_Run]***
For example: ***/Users/username/Data/Database/BMWR30/FOV2/151123_Run1***
Inside there would be the following:
1. The volume used to trace spines and dendrites: `expended_new.tif`
2. An swc file for the centerline of the dendrites: `dendrite.swc`
3. A mat file with information about the masks: `prepareMasksAutoSWC.mat`
4. A image of the imaged dendrites as they appear in the original volume of that field of view (used to register the session to the field of view): `Session.tif`

## Example session for analysis
If you choose to run the code on a single machine  locally you will be able to see the results of most of the steps but only on a subset of the data due to high demand of memory our analysis requires.

###  Example data
Please download the `data.zip` file from [figshare](https://figshare.com/s/3d6d65a09a3b3bd7af1e).
Please change the paths in the notebook according to the location of the files after you extract them.

### Notebooks
An example notebook `151123_BMWR30_Example.ipynb` shows how to run locally.
If you have sufficient resources, you can run the same notebook but with loading all of the data and with default parameters for the registration steps (number of clusters etc.).
