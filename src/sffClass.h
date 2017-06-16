/****************************
*	sffClass.h
*	cstoribal
*	10-04-17
****************************/

/*
Aimed at holding test functions
*/

#ifndef SFFCLASS_H_
#define SFFCLASS_H_


#include <iostream>
#include <string>

#include <opencv2/core.hpp>

#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fstream>
#include <vector>
#include <cstdio>
#include <fstream>
#include <cassert>
#include <utility>

#include "gco/GCoptimization.h"

#include "IOWizard.h"
#include "sharpnessOperator.h"
#include "depthEstimator.h"
#include "energy2.h"
#include "testClass.h"
#include "optimization.h"
#include "evaluation.h"

using namespace cv;
using namespace std;



class MySFF{
public:
    MySFF();~MySFF();
    TestClass testClass;
    IOWizard ioWizard;
    Sharpness_Operator sharpOP;
    DepthClass depthEst;
    EnergyClass energyClass;
    OptiClass optiClass;
    EvalClass evalClass;


    tdf_input input_prts;

    typedef_imgset imageSet;
    typedef_imgset sharpSet;
    tdfp_depth   depth_parameters;
    tdfp_opti     opti_prts;
    
    Mat1T gt_dmat;    // ground truth matrix
    Mat1T dmat;       // local depth matrix
    Mat1T dmat_rank;  // how to compute multifocus
    Mat1T dmat_score; // 
    //Mat3f imat;
    //vector<Mat1f> ivmat;
    Mat1T gmat;
    Mat3T image_MF; // Multifocus	
    // Mat1T smat;
    Mat1T rmat; // Matrice de l'image régularisée 

    fType rmse;
    fType q_rmse;
    fType psnr;
    Point A_tst;  




// Fonctions
    bool loadProblem(int argc, char** argv);
    
    bool preTreat(void);
    
    bool doSharpness(void);

    bool doDepth(void);
   
    bool showInterpolationAt(void);   
    bool clickInterpolation(Mat image, int timer);
    bool showInterpolation(Point A);

    bool setMultifocus(void);

    bool setMultifocusRmat(void);

    bool testEnergy(void);

    bool optimize(void);

    bool evaluate(void);

    //bool testPrintDepths(void);
    
    
private:
    int dim1, dim2;
    int nb_labels;
    

};


bool forwarder(void* context, Point A);

#endif // SFFCLASS_H_
