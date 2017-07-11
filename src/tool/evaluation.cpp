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

bool EvalClass::setup(MyLog* mylog, DepthClass* depthClass){
    this->myLog = mylog;
    this->p_depth = depthClass;
    return true;
}

bool EvalClass::set_parameters(const Mat1T & gt_dmat, const vector<fType> & labels){
    this->gt_dmat = gt_dmat;
    height = gt_dmat.rows;
    width  = gt_dmat.cols;
    this->labels = labels;
    this->step = p_depth->getMeanFocusStep()[0]; //vnl

    this->set = true;
    return true;
}


bool EvalClass::compute_RMSE_label(const Mat1T & dmat, fType & rmse, fType & q_rmse ){
    Mat1T tmpmat;
    tmpmat = dmat/this->step-gt_dmat/this->step;
    pow(tmpmat,2,tmpmat);
    rmse = sqrt(cv::sum(tmpmat)[0]/(height*width));
    q_rmse = 1/rmse;

    string logout = " * RMSE_label is : " + to_string2(rmse) + " and Q is : " + to_string2(q_rmse) + "\n";
    myLog->a(logout);
    return true;
}


bool EvalClass::compute_RMSE(const Mat1T & dmat, fType & rmse, fType & q_rmse ){
    Mat1T tmpmat;
    tmpmat = dmat-gt_dmat;
    pow(tmpmat,2,tmpmat);
    rmse = sqrt(cv::sum(tmpmat)[0]/(height*width));
    q_rmse = 1/rmse;

    string logout = " * RMSE is : " + to_string2(rmse) + " and Q is : " + to_string2(q_rmse) + "\n";
    myLog->a(logout);
    return true;
}


bool EvalClass::compute_PSNR(const Mat1T & dmat, fType & psnr ){


    return true;
}
