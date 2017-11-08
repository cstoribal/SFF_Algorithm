/****************************
*	utils.cpp
*	cstoribal
*	07-11-17
****************************/

/*
Various functions that are usefull for quickly getting information on images
*/

#include "utils.h"




bool Utils::build_int_histogram(const cv::Mat1i & mat_in,
		int nb_labels, std::vector<size_t> & hist){
    //faire un cvminmaxloc pour vérifier erreur éventuelle
    hist.resize(nb_labels);
    for(int i=0; i<nb_labels; i++){ hist[i]=0;}
    for(int i=0; i<mat_in.rows; i++) for(int j=0; j<mat_in.cols; j++){
        hist[mat_in.at<int>(i,j)]++;
    }
    
    return true;
}


float Utils::compute_rmse_label(const cv::Mat1i & mat_in1,
			const cv::Mat1i & mat_in2, size_t nb_labels){
    if(mat_in1.cols!=mat_in2.cols || mat_in1.rows!=mat_in2.rows || mat_in1.cols==0 ){
        return Utils::errormsg("RMSE failure, dimensions matrix must agree");}
    cv::Mat1f tmpmat1, tmpmat2;
    mat_in1.convertTo(tmpmat1,CV_32F,1.0);
    mat_in2.convertTo(tmpmat2,CV_32F,1.0);
    tmpmat1=tmpmat1-tmpmat2;
    pow(tmpmat1,2,tmpmat1);
    /*
    Rect roi;
    if(mat_in1.cols>160 & mat_in1.rows>160) roi = Rect(30,30,mat_in1.cols-60,mat_in1.rows-60);
    else {
        roi = Rect(0,0,mat_in1.cols,mat_in1.rows);
    }
    tmpmat1=cv::Mat1f(tmpmat1,roi);
    */
    float inv_nbpixels = 1.0f/(float)(mat_in1.cols*mat_in1.rows);
    float rmse = sqrt(cv::sum(tmpmat1)[0]*inv_nbpixels)/(float)nb_labels;
    return rmse;
}



bool Utils::errormsg(std::string text){
    CPING(std::string("error")+text);
    return false;
}
