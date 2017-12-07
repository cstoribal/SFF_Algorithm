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
#include "../tool/depthEstimator.h"
#include "../tool/evaluation.h"
#include "../io/logs.h"
#include "../io/IOWizard.h"
#include "../tool/optiplan.h"

extern "C"{
    #include <unistd.h>
}

using namespace std;
using namespace cv;



class OptiClass{
public:
    OptiClass(); ~OptiClass();
    string name_opti; // sets what kind of optimization method is used 
    string selected_typename; // sets what kind of optimization method is used 

    bool setlogs(IOWizard* _ioW, MyLog* mylog);
    bool set_param(tdfp_opti popti, const Mat1T & _gt_dmat);
    bool set_blindestimation(const Mat1T & _blindmat);
    bool set_optiplan(OptiPlan* _p_OptiPlanner);
    bool do_optimization(void); //deprecated


    //bool select_optimization_method(int method);	//TODO later
    bool select_optimization_method(std::string _type);	//TODO later
    bool do_all_optimizations(void);
    //bool set_optimization(void);
    //bool compute_optimization(void);
    bool writebackmatrix(Mat1T & do_mat);
    bool set_allneighbors(void);
    bool reset(fType l_d=-1, fType l_r=0, bool _new_EdataMatrix=true);
    

private:
// common
    bool set;
    MyLog* myLog;
    IOWizard* ioW;
    EnergyClass *energyClass;
    DepthClass	*depthClass;
    EvalClass   *evalClass;

    int nb_pixels;
    int nb_labels;
    int width;
    int height;
    std::vector<fType> labels;
    int connexity;
    size_t actual_idx_method;
    
    fType lambda; // pour nommer fichiers //labmda = l_r/l_d
    bool newfolders,new_EdataMatrix;
    // ce sont des labels, donc des entiers ?
    std::vector<eType> data_in; //poids data
    std::vector<int> data_out;  //labels output

    std::vector<size_t> v_selected_method;
    std::vector<std::string> v_types;
    std::vector<vector<fType> > vv_rmse;	//TODO size + shows

    //nap-shield
    bool build_rank4xy(void);
    cv::Mat1i          getrank;
    std::vector<Point>  getxy;

    bool convert_mat2labvec(const vector<Mat1E> & vmat, vector<eType> & vect);
    bool convert_vec2mat(const vector<int> & vect, Mat1T & vmat);
    bool convert_vec2mat(const vector<int> & vect, Mat1i & vmat);

    
// used by graph cuts
    GCoptimization *gco;
    string gcotype;

    std::vector<eType> smoothvect;  //label*label
    std::vector<Mat1E> Edata;
    //vector<unsigned int> histogram;
    int nb_neighbors;
    Mat1E neighbor_mat; //other way of describing neighbors? weights
    cv::Mat1i regularized_labelmat;
    Mat1T     regularized_depthmat;
    Mat1T     gt_dmat;

    bool nbs_set;
    int* nbs_nb;
    int* nbs_n1D, *nbs_nk1D;
    eType* nbs_w1D, *nbs_wk1D;
    int* nbs_nbk;     // working
    int** nbs_N;      // working
    eType** nbs_W;   // working

    // get opti_plan
    OptiPlan* p_OptiPlanner; // renvoie les seuils & centroides.
    
    // Functions
    bool set_optimization_gco_grid(void);
    bool compute_gco_grid(void);
    bool compute_opt_custom(void);

// Visualisation
    Mat1T blindmat;
    bool show_rmse(void);	//TODO ioW show... and mkdir.
    bool show_all_rmse(void);
    bool gnuplot_vect(FILE* gnuplot, vector<fType> vect);
    bool compute_write_cross_RMSE(void);

// others
    bool error(std::string text=""); 
};


#endif //OPTI_H_
