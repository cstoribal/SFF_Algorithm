/****************************
*	ProjSFF
*	sharpnessOperator.h
*	cstoribal
*	09-04-17
****************************/

/*
Computing sharpness matrix
*/

#ifndef SHARPNESSOPERATOR_H_
#define SHARPNESSOPERATOR_H_


#include <iostream>
#include <string>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fstream>
#include <vector>

#include "../misc/miscdef.h"

using namespace cv;
using namespace std;


class Sharpness_Operator{


public:
    Sharpness_Operator();
    ~Sharpness_Operator();
    
    string type;
    
    bool optypeSelector(string type);
    bool computeOp(string optype, const typedef_imgset & iset, typedef_imgset& sset);
    bool compute(typedef_imgset iset, typedef_imgset& sset);
    
    
    bool compute_SMLAP(vector<Mat1d> ivmat, vector<Mat1d>& smat);
    bool compute_DLAP(vector<Mat1d> ivmat, vector<Mat1d>& smat);


    bool combinRule(vector<Mat1d> in_svmat, vector<Mat1d>& out_svmat);




};


#endif  // SHARPNESSOPERATOR_H_
