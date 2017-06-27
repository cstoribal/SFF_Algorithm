# SFF_Algorithm

## installation
Requires OpenCV, MathGL (cuda support not needed)
To compile, use :
<code>
cmake <srcfolder>
make
</code>
And try "./ProjSFF help all" if you need help. In case you should use a .txt file to store all your parameters, syntax is :

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



