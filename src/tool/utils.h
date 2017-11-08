/****************************
*	utils.h
*	cstoribal
*	07-11-17
****************************/

/*
blablabla
*/

#ifndef _IMGPROCESS_H_
#define _IMGPROCESS_H_

#include "../misc/miscdef.h"

namespace Utils{
    
    bool build_int_histogram(const cv::Mat1i & mat_in, 
		int nb_labels, std::vector<size_t> & hist);

    float compute_rmse_label(const cv::Mat1i & mat_in1,
			const cv::Mat1i & mat_in2, size_t nb_labels);

    bool errormsg(std::string text);

}






#endif
