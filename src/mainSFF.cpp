/**************************
*	ProjSFF
*	mainSFF.cpp
*	cstoribal
*	07-04-17
**************************/

#include <iostream>
#include <string>

#include <ctime> 
#include <stdio.h>
//#include <mgl2/qt.h>
#include <mgl2/mgl.h>





#include "sffClass.h"


using namespace cv;
using namespace std;
vector<bool> MiscClass::optional_features(16);
MySFF mySFF; // global eurk variable


int main( int argc, char** argv )
{
    
    //TestClass myTest;
    //myTest.gethalfindex(58);
    //return 0;
    

    mySFF.setlogs();
    MyLog* myLog = &mySFF.myLog; //alias
    //std::time_t t0 = std::time(0);
    //std::time_t t1 = std::time(0);
    myLog->time_i();
    mySFF.loadProblem(argc, argv);
    double t0 = myLog->time_r(0);
    //myLog->a(to_string2(t0)+" seconds for building dataset\n");
    mySFF.preTreat();
    myLog->time_r(1);
    //myLog->a(to_string2(t0)+" seconds for pretreatment\n");
    mySFF.doSharpness();
    myLog->time_r(2);
    //myLog->a(to_string2(t0)+" seconds for computing sharpness\n");

    mySFF.doDepth();
    myLog->time_r(3);
    //myLog->a(to_string2(t0)+" seconds for computing interpolation\n");
    
    //mySFF.testEnergy();
    mySFF.ioWizard.img_setscale(1);
    //mySFF.showInterpolation(Point(140,10));
    if(MiscClass::optional_features[1]){mySFF.ioWizard.showImage("scale","Idist.png",Mat::zeros(1,1,CV_TF),1);}//mySFF.dmat,1);}
    mySFF.showInterpolationAt();
    //mySFF.showInterpolation(Point(10,11));
    mySFF.ioWizard.img_unsetscale();
    Mat1T tmpmat;
    log(mySFF.dmat_score,tmpmat);



    myLog->time_i();
    mySFF.setMultifocus(); //old method. 
    mySFF.ioWizard.writeImage("I_Multifocus.png",mySFF.image_MF);
    myLog->time_r(4);

    mySFF.ioWizard.img_setscale(1);
    //myLog->a(to_string2(t0)+" seconds for misc (multifocus & dmatscore)\n");

    mySFF.prepare_optimization_plan();
    
    mySFF.optimize2();

    //mySFF.setMultifocusRmat();
    mySFF.ioWizard.img_unsetscale();  
    //mySFF.ioWizard.writeImage("I_Multifocus_regu.png",mySFF.image_MF);
    
    //Mat1T reliability_matrix;
    //mySFF.energyClass.eParams.rmat.copyTo(reliability_matrix);
        
    return 0;

}




