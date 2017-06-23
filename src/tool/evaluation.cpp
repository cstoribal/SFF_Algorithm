/***************************
	evaluation.cpp
	cstoribal
	24-05-17
***************************/

/*
Given the groundtruth, this class is able to
compute as many quality measures as needed.
RMSE, PSNR, UQI
*/

#include "evaluation.h"


EvalClass::EvalClass(){ this->set = false; }
EvalClass::~EvalClass(){}

bool EvalClass::setlogs(MyLog* mylog){
    this->myLog = mylog;
    return true;
}

bool EvalClass::set_parameters(const Mat1T & gt_dmat, const vector<fType> & labels){
    this->gt_dmat = gt_dmat;
    height = gt_dmat.rows;
    width  = gt_dmat.cols;
    this->labels = labels;

    this->set = true;
    return true;
}


bool EvalClass::compute_RMSE(const Mat1T & dmat, fType & rmse, fType & q_rmse ){
    Mat1T tmpmat;
    tmpmat = dmat-gt_dmat;
    pow(tmpmat,2,tmpmat);
    rmse = sqrt(cv::sum(tmpmat)[0]/(height*width));
    q_rmse = 1/rmse;
    cout << "RMSE is : " << rmse << " and Q is : " << q_rmse << endl; 
    return true;
}


bool EvalClass::compute_PSNR(const Mat1T & dmat, fType & psnr ){


    return true;
}
