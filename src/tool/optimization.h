/****************************
*	ProjSFF
*	optimization.h
*	cstoribal
*	03-05-17
****************************/

/*

*/

#ifndef OPTI_H_
#define OPTI_H_

#include <iostream>
#include <string>
#include <sstream>

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
#include "../gco/GCoptimization.h"
#include "../tool/energy2.h"
#include "../io/logs.h"

using namespace std;
using namespace cv;









class OptiClass{
public:
    OptiClass(); ~OptiClass();
    string name_opti; // sets what kind of optimization method is used 

    bool setlogs(MyLog* mylog);

    bool set_param(tdfp_opti popti);
    bool set_optimization(void);
    bool compute_optimization(void);
    bool writebackmatrix(Mat1T & do_mat);
    

private:
// common
    bool set;
    MyLog* myLog;
    EnergyClass *energyClass;

    int nb_pixels;
    int nb_labels;
    int width;
    int height;
    vector<fType> labels;

    // ce sont des labels, donc des entiers ?
    vector<eType> data_in; //poids data
    vector<int> data_out;  //labels output


    //nap-shield
    bool build_rank4xy(void);
    Mat1i          getrank;
    vector<Point>  getxy;

    bool convert_mat2labvec(const vector<Mat1E> & vmat, vector<eType> & vect);
    bool convert_vec2mat(const vector<int> & vect, Mat1T & vmat);

    
// used by graph cuts
    GCoptimization *gco;
    string gcotype;

    vector<eType> smoothvect;  //label*label
    int nb_neighbors;
    Mat1b neighbor_mat; //other way of describing neighbors?
    vector<vector<Point> >  neighborhoods; //wow
    vector<vector<eType> > weights;       //if needed


    bool set_optimization_gco_grid(void);
    bool compute_gco_grid(void);

    bool set_optimization_gco_gen(void);
    bool compute_gco_gen(void);

// others
    



    
};


#endif //OPTI_H_
