/****************************
*	miscdef.h
*	cstoribal
*	12-04-17
****************************/

/*
Defines - stores generic function, typedef & structures
*/

#ifndef MISCDEF_H_
#define MISCDEF_H_


#include <iostream>
#include <string>

#include <opencv2/core.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fstream>
#include <vector>
#include <math.h>

#define CPING(f) (cout << f << endl)
#define CPING2(x,y) (cout << x << "  " << y << endl)

#define COUT(f) (cout << f << endl)
#define COUT2(x,y) (cout << x << "  " << y << endl)

#define SFF_PRECISION double
#if SFF_PRECISION == double
    typedef double      fType;
    typedef cv::Mat1d   Mat1T;
    typedef cv::Mat3d   Mat3T;
    //typedef std::Vec1d     Vec1T;
    //typedef std::Vec3d     Vec3T;
    #define CV_TF       CV_64F   
    #define CV_TFC1     CV_64FC1 
    #define CV_TFC3     CV_64FC3 
#else
#if SFF_PRECISION == float
    typedef float       fType;
    typedef cv::Mat1f   Mat1T;
    typedef cv::Mat3f   Mat3T;
    //typedef std::Vec1f     Vec1T;
    //typedef std::Vec3f     Vec3T;
    #define CV_TF       CV_32F   
    #define CV_TFC1     CV_32FC1 
    #define CV_TFC3     CV_32FC3 
#endif
#endif
//#include "energy2.h"


namespace std {
    template<typename T>
    std::string to_string2(const T &n) {
        std::ostringstream s;
        s << n;
        string str=s.str();
        /* 
        int offset{1}; 
        if (str.find_last_not_of('0') == str.find('.'))
            { offset = 0; } 
        str.erase(str.find_last_not_of('0') + offset, string::npos);
        */
        return str; 
    }
}






struct struct_input{
    // holds data from parsing file. Can hold almost all parameters ? Must have something for "set", eg, 1st row of each vector, 1 or 0.
    int file1_set; //0
    std::string file1_path;
    int         file1_sizeind;
    std::string file1_sep;
    std::string file1_ext;
    int         file1_firsti;
    int         file1_deltai;
    int         file1_lasti;
    

    int file2_set; //1
    std::vector<std::string> file2;

    int groundt_set;//2
    std::string gtpath;
    double gta,gtb;
    
    int outputf_set;//3
    std::string outputfolder;
    
    int focus_set;//4
    std::vector<double> focus;
 
    int preproc_set; //5
    double scale;
    int gauss; //windowsize
    double noise1,noise2,noise3;

    int sharp_set; //6
    std::string sharp;

    int depth_set;//7
    std::string depth;

    int nrj_set;     //8
    std::string nrj_d,nrj_r;
    
    int opti_set;  //9
    std::string opti;

    int optir_set;   //10
    std::vector<double> optirange;

};
typedef struct struct_input tdf_input;


struct tagged_img{
    std::vector<cv::Mat1d> ivmat;
    double focus;
    double dpth;
    int dim; //dimension of ivmat vector
    std::string name; //name of file
    int rank; // file number
}; 
typedef std::vector<struct tagged_img> typedef_imgset;





struct paramdepth{      // type d'interpolation, paramètres. Mal nommé.
    std::string type;
    std::vector<cv::Mat1d> vmat;
    int degree;
};
typedef struct paramdepth tdfp_depth;




struct energy_sidestruct{
    std::string typeD, typeR; //types de mesure d'énergie
    std::string info;
    int dim1, dim2;
    int autoroi; // set 1 if algorithm must seek for neighbors modified pixels (resulting in a double roi augmentation)
    cv::Mat1d dmat; // Matrice de profondeur initiale, lue à partir des données locales uniquement
    cv::Mat1d rmat; // Reliability matrix. Calculée à partir des données de sharpness
    double scale; // tmp
}; 
typedef struct energy_sidestruct tdfp_energy ;



union interpol_func
{
    bool (*fptr)(void*,cv::Point);
    void* context;
};









#endif // MISCDEF_H_














