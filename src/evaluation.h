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

#include "miscdef.h"

using namespace std;
using namespace cv;



class EvalClass{
public:
    EvalClass();~EvalClass();

    bool set_parameters(const Mat1T & gt_dmat, const vector<fType> & labels);
    
    bool compute_RMSE(const Mat1T & dmat, fType & rmse, fType & q_rmse );
    bool compute_PSNR(const Mat1T & dmat, fType & psnr);

    bool compute_quantization_noise(void);
    
private:
    bool set;
    Mat1T gt_dmat;
    int height, width;
    vector<fType> labels;
    
    
    
    
    
};


#endif // EVALUATION_H_
