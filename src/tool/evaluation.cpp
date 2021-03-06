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

    if(height>160 & width>160){
        this->roi = Rect(30,30,width-60,height-60);
        roi_surface=(width-60)*(height-60);
    } else {
        myLog->a("Image is too small, no roi applied to evaluation \n");
        this->roi = Rect(0,0,width,height);
        roi_surface=width*height;
    }
    
    
    if(1){
        this->roi=Rect(0,0,width,height); 
        roi_surface=width*height;//TODO disabled roi
    }
      

    this->set = true;
    return true;
}


bool EvalClass::compute_RMSE_label(const Mat1T & dmat, fType & rmse){
    Mat1T tmpmat;
    tmpmat = dmat/this->step-gt_dmat/this->step;
    pow(tmpmat,2,tmpmat);
    tmpmat=Mat(tmpmat,roi);
    rmse = sqrt(cv::sum(tmpmat)[0]/(roi_surface))/(fType)this->labels.size();

    string logout = " - RMSE is : " + to_string2(rmse) + "  ";
    CPING(logout);
    myLog->a(logout);
    return true;
}


bool EvalClass::compute_PSNR(const Mat1T & dmat, fType & psnr ){
    Mat1T tmpmat;
    tmpmat = dmat-gt_dmat;
    pow(tmpmat,2,tmpmat);
    tmpmat=Mat(tmpmat,roi);
    psnr = cv::sum(tmpmat)[0]/(roi_surface);
    psnr = 10*log10( (this->labels[this->labels.size()-1]-this->labels[0])*
		(this->labels[this->labels.size()-1]-this->labels[0]) /psnr);
    string logout = " - PSNR is : " + to_string2(psnr) + "\n";
    myLog->a(logout);
    
    return true;
}


bool EvalClass::compute_RMSE(const Mat1T & dmat, fType & rmse, fType & q_rmse ){
    //Deprecated
    assert(0);
    Mat1T tmpmat;
    tmpmat = dmat-gt_dmat;
    pow(tmpmat,2,tmpmat);
    tmpmat=Mat(tmpmat,roi);
    rmse = sqrt(cv::sum(tmpmat)[0]/(roi_surface));
    q_rmse = 1/rmse;

    string logout = " * RMSE is : " + to_string2(rmse) + " and Q is : " + to_string2(q_rmse) + "\n";
    myLog->a(logout);
    return true;
}

Mat1T EvalClass::check_diff(const Mat1T & dmat){
    Mat1T outmat = abs(dmat - gt_dmat);
    return outmat;
}

bool EvalClass::compute_RMSE_label(const Mat1T & dmat, fType & rmse, fType & q_rmse ){
    //deprecated
    CPING("warning, bad version of rmse");
    assert(0);
    Mat1T tmpmat;
    tmpmat = dmat/this->step-gt_dmat/this->step;
    pow(tmpmat,2,tmpmat);
    tmpmat=Mat(tmpmat,roi);
    rmse = sqrt(cv::sum(tmpmat)[0]/(roi_surface));
    q_rmse = 1/rmse;

    string logout = " * RMSE_label is : " + to_string2(rmse) + " and Q is : " + to_string2(q_rmse) + "\n";
    myLog->a(logout);
    return true;
}
