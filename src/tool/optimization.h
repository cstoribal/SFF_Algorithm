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
#include <cmath>

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
    bool do_optimization(void);
    //bool set_optimization(void);
    //bool compute_optimization(void);
    bool writebackmatrix(Mat1T & do_mat);
    bool set_allneighbors(void);
    bool reset(void);
    

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
    int connexity;

    // ce sont des labels, donc des entiers ?
    vector<eType> data_in; //poids data
    vector<int> data_out;  //labels output


    //nap-shield
    bool build_rank4xy(void);
    Mat1i          getrank;
    vector<Point>  getxy;

    bool convert_mat2labvec(const vector<Mat1E> & vmat, vector<eType> & vect);
    bool convert_vec2mat(const vector<int> & vect, Mat1T & vmat);
    bool convert_vec2mat(const vector<int> & vect, Mat1i & vmat);

    
// used by graph cuts
    GCoptimization *gco;
    string gcotype;

    vector<eType> smoothvect;  //label*label
    int nb_neighbors;
    Mat1E neighbor_mat; //other way of describing neighbors? weights
    //vector<vector<Point> >  neighborhoods; //wow
    //vector<vector<eType> > weights;       //if needed
    //Neighborsystem nbs_
    bool nbs_set;
    int* nbs_nb;
    int* nbs_n1D;
    int** nbs_n;
    eType* nbs_w1D;
    eType** nbs_w;
    eType* nbs_wk1D;
    eType** nbs_wk;


    bool set_optimization_gco_grid(void);
    bool compute_gco_grid(void);

    bool set_optimization_gco_gen(void);
    bool compute_gco_gen(void);

    bool compute_opt_binary(void);
    bool compute_opt_multiscale(void);

// others
    



    
};


#endif //OPTI_H_
