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

#include "miscdef.h"
#include "gco/GCoptimization.h"
#include "energy2.h"

using namespace std;
using namespace cv;



struct set_param_opti{
    string      type;
    EnergyClass *energyclass;
    int nb_pixels, nb_labels, width, height;
    vector<double> labels; // allows conversion label -> ranklabel
    
    
    
};
typedef struct set_param_opti tdfp_opti;












class OptiClass{
public:
    OptiClass(); ~OptiClass();
    string name_opti; // sets what kind of optimization method is used 

    bool set_param(tdfp_opti popti);
    bool set_optimization(void);
    bool compute_optimization(void);
    bool writebackmatrix(Mat1d & do_mat);
    

private:
// common
    bool set;
    EnergyClass *energyClass;

    int nb_pixels;
    int nb_labels;
    int width;
    int height;
    vector<double> labels;


    vector<double> data_in;
    vector<double> data_out;


    //nap-shield
    bool build_rank4xy(void);
    Mat1i          getrank;
    vector<Point>  getxy;

    bool convert_mat2labvec(const vector<Mat1d> & vmat, vector<double> & vect);
    bool convert_vec2mat(const vector<double> & vect, Mat1d & vmat);

    
// used by graph cuts
    GCoptimization *gco;
    string gcotype;

    vector<double> smoothvect;  //label*label
    int nb_neighbors;
    Mat1b neighbor_mat; //other way of describing neighbors?
    vector<vector<Point> >  neighborhoods; //wow
    vector<vector<double> > weights;       //if needed


    bool set_optimization_gco_grid(void);
    bool compute_gco_grid(void);

    bool set_optimization_gco_gen(void);
    bool compute_gco_gen(void);

// others
    



    
};


#endif //OPTI_H_
