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

#include "../io/logs.h"
#include "../misc/miscdef.h"

using namespace cv;
using namespace std;


class Sharpness_Operator{


public:
    Sharpness_Operator();
    ~Sharpness_Operator();
    
    string type;
    bool setlogs(MyLog* myLog);
    bool optypeSelector(string type);
    bool computeOp(string optype, const tdf_imgset & iset, tdf_imgset& sset);
    bool compute(tdf_imgset iset, tdf_imgset& sset);
    
    
    bool compute_SMLAP(vector<Mat1T> ivmat, vector<Mat1T>& smat);
    bool compute_DLAP(vector<Mat1T> ivmat, vector<Mat1T>& smat);
    bool compute_3DLAP(tdf_imgset & ioset);
    bool compute_STA2(vector<Mat1T> ivmat, vector<Mat1T>& smat);


    bool combinRule(vector<Mat1T> in_svmat, vector<Mat1T>& out_svmat);

private:
    MyLog* myLog;

    int height;
    int width;





};


#endif  // SHARPNESSOPERATOR_H_
