/****************************
*	testClass.h
*	cstoribal
*	09-04-17
****************************/

/*
Aimed at holding test functions
*/

#ifndef TESTCLASS_H_
#define TESTCLASS_H_


#include <iostream>
#include <iomanip>
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
#include <cmath>




#include "miscdef.h"

using namespace cv;
using namespace std;


class TestClass{


public:
    TestClass();
    ~TestClass();

    
    bool coutMatType(Mat matrix);
    // pas de multiplication de matrices
    bool matrixOp1(Mat& matrix); // multiplication scalaire ok
    bool matrixOp2(Mat& matrix); // seuillage OK
    bool matrixOp3(Mat& M);
    bool rgbToGray(Mat3f M, Mat1f& N);
    bool cmatToArraymat(Mat3f M, vector<Mat1f>& Vm);
    bool fillSharpPoly(typedef_imgset & sharpSet);
    bool fillAMatrix(Mat1d & imat);
    int polyfit();
    

};





#endif // TESTCLASS_H_
