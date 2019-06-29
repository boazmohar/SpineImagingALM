# Spine imaging
This repository is accompanying the paper:
***Functional clustering of dendritic activity during decision-making*** 
Kerlin, A., Mohar, B., Flickinger, D., MacLennan, B.J., Davis, C., Spruston, N., and Svoboda, K. (2018). [BioRxiv 440396](https://www.biorxiv.org/content/early/2018/10/10/440396)
## Requierments

### Spark
To run this code you would need a [Spark cluster](https://spark.apache.org/) running to preform parallel computing (tested with version 2.1 and 2.2, sholud work with a newer versions as well). It will be possible to run it on a single computer [(See here)](https://spark.apache.org/docs/latest/spark-standalone.html), but it would take a very long time to procces a single session. We rotinly used 100-300 cores to run registration,  time course extraction and other analysis.

### Python

This code was tested with both python 2.7 and 3.6. See the dependencies in the `requiermnts` folder with the `Conda_*.txt` flies for using Anaconda or a `requiermnts.txt` for pip.
### Folder structure
The code was run with some assumption as to the location of the raw / proccesed data.

#### Raw data
We assume the following structure of the raw data: ***[base_folder]/[date]_[animal]/[Run#]***
For example: ***/Users/username/Data/151123_BMWR30/Run1***
Inside would be:

 - `*.tif` files each representing a behavioral trial.
Each file would be 3d array of images that is:
	1. (z) #planes per volume x  #volumes per trial
	2. x pixels
	3. y pixels
 - `Sp.mat` file has information from our [trajectory design GUI](www.github.com) about the location and timing of each imaging field.
 - `behavior.mat` file with information regarding the type of the behavioral trial, when and where the animal licked etc.

#### Procced data
After registration, a few intermidate files are saved.
1. Pre-procssed version of the raw data (CleanBinary)
2. Registered version of the data (RegBInary)
3. Metadata about the session and parameters used (Step_x.p)
4. Registration intermidiates in folders: `Run1Intermediates`, `Run1View`

At this point these have defualts that are right for our cluster folder structure and could be changed.

After tracing of the registered session is complete a few more files are made in a different location based on the followind structure:
***[base_folder]/Database/[Animal]/[FOV]/[Date_Run]***
For example: ***/Users/username/Data/Database/BMWR30/FOV2/151123_Run1***
Inside there would be the following:
1. The volume used to trace spines and dendrties: `expended_new.tif`
2. An swc file for the centerline of the dendrties: `dendrite.swc`
3. A mat file with information about the masks: `prepareMasksAutoSWC.mat`
4. A image of the imaged dendrites as they appere in the original volume of that field of view (used to register the session to the field of view): `Session.tif`

## Example session for analysis
If you choose to run the code on a single mechine localy you will be able to see the results of most of the steps but only on a subset of the data due to high demand of memory our analysis requires.

### Notebooks
An example notebook `151123_BMWR30_Example.ipynb` shows how to run localy.
If yout have sufficent resorces you can run the same notebook but with loading all of the data and with deafult parameters for the registration steps (number of clusters etc.).
 
###  Example data
Please download the `data.zip` file from [figshare](https://figshare.com/s/3d6d65a09a3b3bd7af1e).
Please change the paths in the notebook according to the location of the files after you extract them.
 

