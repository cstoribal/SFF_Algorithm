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

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>


#include "testClass.h"
#include "sffClass.h"
// #include "depthEstimator.h"


using namespace cv;
using namespace std;

MySFF mySFF; // global eurk variable


int main( int argc, char** argv )
{
    

    std::time_t t0 = std::time(0);
    std::time_t t1 = std::time(0);
    
    TestClass testClass;    
    //MySFF mySFF;
    string tmp;
    //vector<vector<string> > testvec;
    //Mat3b imat;
    //Mat1b imatgrey;
    //mySFF.ioWizard.parsefile2vect("data.txt",testvec);
    
    mySFF.loadProblem(argc, argv);

    
    t1 = std::time(0);
    printf(" %li seconds for buildingdataset\n", t1-t0);
    mySFF.preTreat(); //TODO add noise !

    t0 = std::time(0);
    printf(" %li seconds for pretreatment\n", t0-t1);

    mySFF.doSharpness();
    t1 = std::time(0);
    printf(" %li seconds for computing sharpness\n", t1-t0);

    mySFF.doDepth();
    t0 = std::time(0);
    printf(" %li seconds for computing interpolation\n", t0-t1);
    
    //mySFF.testEnergy();

    //mySFF.ioWizard.showImage("scale",mySFF.dmat,1000);
    mySFF.ioWizard.writeImage("Idist.png",mySFF.dmat);

    Mat1d tmpmat;
    log(mySFF.dmat_score,tmpmat);
    //mySFF.ioWizard.showImage("scale",tmpmat,5000);
    mySFF.ioWizard.writeImage("Iscore.png",tmpmat);
 
    mySFF.setMultifocus();
    //mySFF.ioWizard.showImage("scale",mySFF.image_MF,10);
    mySFF.ioWizard.writeImage("I_Multifocus.png",mySFF.image_MF);

    //mySFF.clickInterpolation(mySFF.image_MF,0);
    mySFF.ioWizard.showImage("scale",mySFF.dmat,100);
    mySFF.showInterpolationAt();
    t0 = std::time(0);
    printf(" %li seconds for misc\n", t1-t0);

    mySFF.optimize();
    t1 = std::time(0);
    printf(" %li seconds for computing optimisation\n", t0-t1);
    //mySFF.ioWizard.showImage("scale",mySFF.rmat,1000);
    mySFF.ioWizard.writeImage("I_dregu.png",mySFF.rmat);

    mySFF.setMultifocusRmat();
    //mySFF.ioWizard.showImage("scale",mySFF.image_MF,10);
    mySFF.ioWizard.writeImage("I_Multifocus_regu.png",mySFF.image_MF);

    //mySFF.evaluate();
    
    //mySFF.ioWizard.showImage("scale", mySFF.smat, 4000);
    //Mat imgtmp;
    //cvtColor(mySFF.smat,imgtmp,CV_GRAY2BGR);
    //mySFF.ioWizard.showImage("scale", mySFF.imat/2+imgtmp/2, 10000);
    //mySFF.ioWizard.showImage("scale", mySFF.imat, 4000);

    //cout<<mySFF.p_sharpOP<<endl;
    
    
    
    
    
    return 0;

}




