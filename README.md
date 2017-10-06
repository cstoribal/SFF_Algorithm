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
Try "./ProjSFF help all" if you need help. In case you should use a .txt file to store all your parameters, syntax is :

     1/0; -f; -fsi; -fs; ext; -ffi; -fdi; -fli \\n
     1/0; -filenames vector;;; \\n
     1/0; -gt; gta; gtb; \\n
     1/0; -o; \\n
     1/0; -focusdepths(a vector of imagefocus); \\n
     1/0; -sc; -gw; -n1; -n2; -nga; -ngs\\n
     1/0; -psharp; \\n
     1/0; -pdepth; \\n
     1/0; -pnd energy data; -pnr energyregul; \\n
     1/0; -poptimization; \\n
     1/0; -lambda_r vector;;; \\n
     1/0; -lambda_d vector;;; \\n



## HowTo - General introduction

The executables are *ProjSFF* running the main program.

A second executable, named *SFFDataMgmt*, is usefull for modifying all parameters of every targeted samples at once, following a template from a "datafile.txt". It is designed by argument -M, while the path of the targeted datafile is the input from argument -f.
"./SFFDataMgmt -M data1.param.txt -f data2.txt"

For convenience, some scripts can run such executables alternatively.
Those scripts are currently written in *relative paths*.
*Metascript.sh* detects all *.param.txt files ($1) and calls :
    *ChangeParameters.sh* finds .txt in (LIST_FOLD) and changes them to ($1).
    *script_autofolders.sh* does the same (1 directory ahead) and compute.



