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

#include "../io/logs.h"
#include "../misc/miscdef.h"

using namespace std;

class DepthClass{
public:
    DepthClass(); ~DepthClass();
    
    std::string type; // 8poly, gaussian, ...
    bool setlogs(MyLog* mylog);
    bool set_param(string settype);
    bool buildEstimation(void);
    bool buildEstimation(const tdf_imgset & sharpSet, tdfp_depth & pdmat);
    bool buildDepthmat(const tdfp_depth & dparam, Mat1T & dmat, Mat1T & dmat_rank, Mat1T & dmat_score);
    bool          getVmatSharpI(vector<Mat1T> & vmat);
    int           getRankFromDepth(fType input);
    bool getMRankFromMDepth(const Mat1T& minput, cv::Mat1i & moutput);
    bool getMLabelFromMDepth(const Mat1T& minput, cv::Mat1i & moutput);
    int           getNbLabels(void);
    fType         getOversampling(void);
    vector<fType> getLabels(void);
    vector<fType> getMeanFocusStep(void);
    vector<fType> getMeanLogFocusStep(void);
    
    
    bool showInterpolationAt(const vector<cv::Point> & vP, const tdf_imgset & SharpSet, const tdfp_depth & dparam, const Mat1T & dmat, const string & folder);

    
    
private:
    MyLog* myLog;
    bool set;
    bool vmat_sharp_i_set;
    vector<Mat1T> vmat_sharp_i;
    string steptype; //linear, logarithmic
    fType focus_min;
    fType focus_max;
    vector<fType> focus;
    int set_dim[3];  // Image(w,h), Number of images
    int oversampling;
    vector<vector<fType> > DepthToRank; // holds interpolated depth [0] and corresponding picture [1]

////Polynomial interpolation
    bool f_poly(const tdf_imgset & sharpSet, tdfp_depth & pdmat);
    bool d_poly(const tdfp_depth & dparam, Mat1T & dmat, Mat1T & dmat_rank, Mat1T & dmat_score);
    int degree; // degré du polynôme
    
    bool interpolate(const vector<fType> & x, const vector<fType> & y, int n, int N, vector<fType> X, vector<fType> & z);
    bool s_poly_ij(const vector<cv::Point> & vP, const tdf_imgset & sharpSet, const tdfp_depth & dparam, const Mat1T & dmat, const string & folder);

    
////Gaussian interpolation




};





#endif //DEPTH_H_
