/***************************
	evaluation.h
	cstoribal
	24-05-17
***************************/

/*
This class holds :
Public
 - f_ rmse (data) 
 - f_ Q
 - f_ Qr
 - Bruit de quantification

Private
 - A shared pointer to the groundtruth (a ~ cv::mat)
 - width, height
 - Liste des labels
 - 

*/


#ifndef EVALUATION_H_
#define EVALUATION_H_

#include "../io/logs.h"
#include "../misc/miscdef.h"
#include "../tool/depthEstimator.h"

using namespace std;
using namespace cv;



class EvalClass{
public:
    EvalClass();~EvalClass();
    bool setup(MyLog* mylog, DepthClass* depthClass);
    bool set_parameters(const Mat1T & gt_dmat, const vector<fType> & labels);
    
    bool compute_RMSE_label(const Mat1T & dmat, fType & rmse, fType & q_rmse );
    bool compute_RMSE_label(const Mat1T & dmat, fType & rmse);
    bool compute_RMSE(const Mat1T & dmat, fType & rmse, fType & q_rmse );
    bool compute_PSNR(const Mat1T & dmat, fType & psnr);

    Mat1T check_diff(const Mat1T & dmat);

    bool compute_quantization_noise(void);
    
private:
    MyLog* myLog;
    DepthClass* p_depth;

    bool set;
    Mat1T gt_dmat;
    int height, width;
    vector<fType> labels;
    fType step;

    Rect roi;
    
    
    
    
    
};


#endif // EVALUATION_H_
