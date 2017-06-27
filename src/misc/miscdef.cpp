/****************************
*	miscdef.cpp
*	cstoribal
*	27-05-17
****************************/

/*
Defines - stores generic function, typedef & structures
*/




#include "miscdef.h"



namespace cv {
    bool minMaxLoc(cv::Mat matrix, float* fmin, float* fmax)
    {
        double dmin, dmax;
        float vmin, vmax;
        cv::minMaxLoc(matrix, &dmin, &dmax);
        *fmin = static_cast<float>(dmin);
        *fmax = static_cast<float>(dmax);
        return true; 
    }
}
