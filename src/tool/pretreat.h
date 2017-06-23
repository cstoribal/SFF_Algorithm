/****************************
*	pretreat.h
*	cstoribal
*	15-06-17
****************************/

/*
Outputs at the same time a precise log of each simulation encountered
and the results under the format of a csv
*/

#ifndef PRETREAT_H_
#define PRETREAT_H_

#include <cstdlib>
#include <ctime>
#include <vector>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>


#include "../misc/miscdef.h"
#include "../io/logs.h"

using namespace cv;
using namespace std;

/*
fType randftype(fType a, fType b)
{
    fType r = ( ((fType) rand()) / (fType) RAND_MAX ) * (a - b) + b;
    return r;
}
*/ //TODO fix problem  of multiple definition here, if needed.
    


class PreTreatment{
public:
    PreTreatment(); ~PreTreatment();
    
    bool setlogs(MyLog* mylog);
    bool set_param(const tdf_input & prts);
    
    bool compute_noises(Mat_<fType> & image);

private:
    //Setup related
    MyLog* myLog;
    fType noise_a, noise_b, noise_ca, noise_cs;
       // multiplicatif normal gaussien (amplitude&sigma)
    //Defined on demand
    int dimx, dimy, dimc;


    bool compute_nmult(Mat_<fType> & image);
    bool compute_nunif(Mat_<fType> & image);
    bool compute_ngauss(Mat_<fType> & image);

};



#endif // PRETREAT_H_
