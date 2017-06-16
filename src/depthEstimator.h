/****************************
*	ProjSFF
*	depthEstimator.h
*	cstoribal
*	18-04-17
****************************/



#ifndef DEPTH_H_
#define DEPTH_H_

#include <string>
#include <vector>

#include "miscdef.h"

using namespace std;

class DepthClass{
public:
    DepthClass(); ~DepthClass();
    
    std::string type; // 8poly, gaussian, ...
    bool set_param(string settype);
    bool buildEstimation(void);
    bool buildEstimation(const typedef_imgset & sharpSet, tdfp_depth & pdmat);
    bool buildDepthmat(const tdfp_depth & dparam, cv::Mat1d & dmat, cv::Mat1d & dmat_rank, cv::Mat1d & dmat_score);
    int getRankFromDepth(fType input);
    int getNbLabels(void);
    vector<fType> getLabels(void);
    
    bool showInterpolationAt(const vector<cv::Point> & vP, const typedef_imgset & SharpSet, const tdfp_depth & dparam, const cv::Mat1d & dmat, const string & folder);

    
    
private:
    bool set;
    fType focus_min;
    fType focus_max;
    vector<fType> focus;
    int set_dim[3];  // Image(w,h), Number of images
    int oversampling;
    vector<vector<fType> > DepthToRank; // holds interpolated depth [0] and corresponding picture [1]

////Polynomial interpolation
    bool f_poly(const typedef_imgset & sharpSet, tdfp_depth & pdmat);
    bool d_poly(const tdfp_depth & dparam, cv::Mat1d & dmat, cv::Mat1d & dmat_rank, cv::Mat1d & dmat_score);
    int degree; // degré du polynôme
    
    bool interpolate(const vector<fType> & x, const vector<fType> & y, int n, int N, vector<fType> X, vector<fType> & z);
    bool s_poly_ij(const vector<cv::Point> & vP, const typedef_imgset & sharpSet, const tdfp_depth & dparam, const cv::Mat1d & dmat, const string & folder);

    
////Gaussian interpolation




};





#endif //DEPTH_H_
