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


/*class OptiIterate{
public:
    OptiIterate(); ~OptiIterate();
    
    bool setup(MyLog* mylog, int* sortedlabel_in, int** adaptindex_in, int maxseek_in, int nblabels_in, int maxiteration_in, int max_sortedrank); // builds the vectors
    bool update(int i); // iterate one more time. start's -1+1, set active & flags

    bool active; //set to false when all instances are at maxseek or labelmax
    bool flag;   // something needs to be checked ? init to 0

    int getlabelp(int state); // pixel call for label number
    int getlabeln(int state); // pixel call for label number

    int getoutlabel(int state);      // résultat = dernier résultat à 1
    int getnext_outlabel(int state); // variante, suppose le résultat = prochain label à 1
                   // => SI il y a un seek supplémentaire à state, le prendre
                   // => SI il n'y a pas de seek supplémentaire, prendre la sortie 0 ou 1


private:
    // at setup
    MyLog* myLog;
    bool set;
    int maxseek;
    int maxiteration; // to get sizes
    int* sorted_label_img;
    int** adapt_index;
    int nblabels;
    int max_sortedrank;

    int getprevious(int current, int iter);
    bool duplicate_to(int prev_iter, int prev_state, int iter, int state);
    bool display(void); // just to show-off the work. And see what's wrong here
    
    // other
    int iteration;    // init to zero
    vector<vector<int> > path; // init to zero
    vector<vector<int> > lmin; // init to zero           -  those values are forbidden
    vector<vector<int> > lmax; // init to maximum label
    vector<vector<int> > seek; // get local depth of iteration
    vector<vector<int> > labelp; // get label corresponding to binary process (01001101)
    vector<vector<int> > labeln; // get label corresponding to binary process (01001101)
    vector<vector<bool> > enabled;// set to false when over maxseek
    vector<vector<int> > labelout;

    int lastiteration;
    
    
};

*/


class OptiClass{
public:
    OptiClass(); ~OptiClass();
    string name_opti; // sets what kind of optimization method is used 
    string selected_typename; // sets what kind of optimization method is used 

    bool setlogs(IOWizard* _ioW, MyLog* mylog);
    bool set_param(tdfp_opti popti, const Mat1T & _gt_dmat);
    bool set_optiplan(OptiPlan* _p_OptiPlanner);
    bool do_optimization(void); //deprecated


    //bool select_optimization_method(int method);	//TODO later
    bool select_optimization_method(std::string _type);	//TODO later
    bool do_all_optimizations(void);
    //bool set_optimization(void);
    //bool compute_optimization(void);
    bool writebackmatrix(Mat1T & do_mat);
    bool set_allneighbors(void);
    bool reset(fType l_d=-1, fType l_r=0);
    

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
    //vector<unsigned int> histogram;
    int nb_neighbors;
    Mat1E neighbor_mat; //other way of describing neighbors? weights
    cv::Mat1i regularized_labelmat;
    Mat1E     regularized_depthmat;
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
    bool show_rmse(void);	//TODO ioW show... and mkdir.
    bool show_all_rmse(void);
    bool gnuplot_vect(FILE* gnuplot, vector<fType> vect);

// others
    bool error(std::string text=""); 
};


#endif //OPTI_H_
