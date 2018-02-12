# SFF_Algorithm

## installation
Requires OpenCV, MathGL (cuda support not needed)
To compile, use :
<code>
cmake <srcfolder>
make
</code>
Creates Cmakefiles with "cmake ./<srcfolder>" and compiles source files with "make"
.
Try "./ProjSFF help all" if you need help and detailed explanation on parameters. In case you should use a .txt file to store all your parameters, syntax is :

     1/0; -f; -fsi; -fs; ext; -ffi; -fdi; -fli \\n
     1/0; -filenames vector;;; \\n
     1/0; -gt; gta; gtb; \\n
     1/0; -o; \\n
     1/0; -focusdepths(a vector of imagefocus); \\n
     1/0; -sc; -gw; -n1; -n2; -nga; -ngs\\n
     1/0; -psharp; \\n
     1/0; -pdepth; oversampling \\n
     1/0; -pnd energy data; -pnr energyregul; \\n
     1/0; -poptimization; connexity ; maxiteration \\n
     1/0; -lambda_r vector;;; \\n
     1/0; -lambda_d vector;;; \\n

The 1/0 are set to 1 if the line is active, 0 if it is inactive and that the program has to use default values.
It is then possible to run SFF with "./ProjSFF -D [path_of_your_.txt_file]". 

* -optf is the only parameter that must be called as an independant parameter. It represents a vector of boolean, each of them being linked to a special parameter, for testing purposes. -optf 0 may be the good default value.
     - bit 0 adds points on the borders of the interpolating interval (bad behaviour), 
     - bit 1 enable "showInterpolation" (see bugs), 
     - bit 2 enables to write the noisy images in a new folder.


## HowTo General introduction

The executable is *ProjSFF* running the main program.

A second executable, named *SFFDataMgmt*, is usefull for modifying all parameters of every targeted samples at once, following a template from a "datafile.txt". It is designed by argument -M, while the path of the targeted datafile is the input from argument -f.
"./SFFDataMgmt -M data1.param.txt -f data2.txt"

It's currently the way our algorithm is called. For convenience, it is possible to run that in a script with gnu_parallel. For the moment, the only thing to care about is that the data.txt files storing informations about the data also store parameters values, and so, the set of parameters of each set of images must be treated sequentially. About the text output concatenated in the files logout.txt, we warn that there is no warranty of atomicity of the read/write operation which also could be a problem. However, it is not difficult to recover these information from subfolders.




## About

### Bugs to fix

#### Various
- ShowInterpolationAt is partly broken, due to a conflict when using showImage from openCV which would ruin the code. Temporary patch forbids the call to showImage in IOWizard, and the coordinates of the pixels is hardcoded. 

#### Concerning parameters 
- Loading the dataset : in the "*.param.txt" file, current versions of the program worked with the vector of filepaths and focus values (Rows 2 and 5). It had first been implemented with Row 1 (with data about file syntax) but it is now deprecated.
- Loading the groundtruth : to compute properly our algorithm, Row 3 may be mandatory. If not needed, any image of the dataset with gta=0 may work. 
- Oversamping (integer) over the depth values is not working (only acceptable value for oversampling is 1). The idea for the moment is that for higher values, somehow, only some components of the shapness profiles get oversampled eg the raw sharpness only for instance, making the results totally uncomparable to the groundtruth. It has almost never been tested. By the way, other objects like for instance the gtmat_label, must be computed differently, and ideas/methods/attributes like "mean step between focus values" (evalClass), "DepthToRank vector" (depthEstimator) may have to be partly redefined.
- nga, ngs - both these parameters deal with gaussian noise. Actually, nga must be set to 1 in order to get the expected values for ngs (sigma) so this feature should be redefined.
- connexity can be set to 4 or 8
- maxiteration is currently not in use.
- there are some limitations in the energy dynamics. For instance, with energy data being nL1Rw1 and lambda d being 10^6 on some samples, the graph cut algorithm reaches its maximum bound value and is therefore not computed, which crashes the program (todo). We suggest to sort values of lambda in ascending order, when possible. Additionally, due to some constraints about the regularization term being a metric and float rounding, we suggest that lambda r  are to be typed in as integers.
- The manual call of parameters (./ProjSFF -ps 1) is not supported anymore, but may work for some.

### ToDoList
- Fixing all the above features-limitations.
- Expanding the inputs to -M so that *SFFDataMgmt* can be replaced totally by *ProjSFF -M model -f dataset*, which would allow parallel processing of multiple parameters over the same image set.
- When calling "./ProjSFF help all", the program bugs instead of bypassing the SFF algorithm, and it's not beautifull.
- The exit value of ProjSFF is not used. It could however be used as an input for gnu_parallel
- The program crashes when graph cut optimisation returns an error.



### History
#### [old]HowTo-Scripts
For convenience, some scripts can run such executables alternatively.
Those scripts are currently written in *relative paths*.
*Metascript.sh* detects all *.param.txt files ($1) and calls :
    *ChangeParameters.sh* finds .txt in (LIST_FOLD) and changes them to ($1).
    *script_autofolders.sh* does the same in the appropriate subfolder and compute.


