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


#include "misc/testClass.h"
#include "sffClass.h"


using namespace cv;
using namespace std;

MySFF mySFF; // global eurk variable


int main( int argc, char** argv )
{

    TestClass myTest;
    myTest.seekmedian(51);
    return 0;

    mySFF.setlogs();
    MyLog* myLog = &mySFF.myLog; //alias
    std::time_t t0 = std::time(0);
    std::time_t t1 = std::time(0);
    
    mySFF.loadProblem(argc, argv);

    
    t1 = std::time(0);
    myLog->a(to_string2(t1-t0)+" seconds for building dataset\n");
    mySFF.preTreat();
/*
    for(int i=0; i<mySFF.imageSet.size(); i++)
    {
        Mat tmpmat2;
        merge(mySFF.imageSet[i].ivmat,tmpmat2);
        mySFF.ioWizard.showImage("scxle",tmpmat2,10);
    }
*/
    t0 = std::time(0);
    myLog->a(to_string2(t0-t1)+" seconds for pretreatment\n");

    mySFF.doSharpness();
    t1 = std::time(0);
    myLog->a(to_string2(t1-t0)+" seconds for computing sharpness\n");

    mySFF.doDepth();
    t0 = std::time(0);
    myLog->a(to_string2(t0-t1)+" seconds for computing interpolation\n");
    
    //mySFF.testEnergy();

    mySFF.ioWizard.writeImage("Idist.png",mySFF.dmat);
    mySFF.ioWizard.img_unsetscale();
    Mat1T tmpmat;
    log(mySFF.dmat_score,tmpmat);
    mySFF.ioWizard.showImage("scale","score",tmpmat,4000);
    mySFF.ioWizard.writeImage("Iscore.png",tmpmat);
 
    mySFF.setMultifocus();
    //mySFF.ioWizard.showImage("scale",mySFF.image_MF,10);
    mySFF.ioWizard.writeImage("I_Multifocus.png",mySFF.image_MF);

    mySFF.ioWizard.img_setscale(1);
    //mySFF.clickInterpolation(mySFF.image_MF,0);
    mySFF.ioWizard.showImage("scale","dmat",mySFF.dmat,1000);
    mySFF.showInterpolationAt();
    t1 = std::time(0);
    myLog->a(to_string2(t1-t0)+" seconds for misc\n");

    mySFF.optimize();
    t0 = std::time(0);
    myLog->a(to_string2(t0-t1)+" seconds for computing optimisation\n");
    //mySFF.ioWizard.showImage("scale","rmat",mySFF.rmat,1000);
    mySFF.ioWizard.writeImage("I_dregu.png",mySFF.rmat);

    mySFF.setMultifocusRmat();
    mySFF.ioWizard.img_unsetscale();
    //mySFF.ioWizard.showImage("scale",mySFF.image_MF,10);
    mySFF.ioWizard.writeImage("I_Multifocus_regu.png",mySFF.image_MF);

    //mySFF.evaluate();
    
    //mySFF.ioWizard.showImage("scale", mySFF.smat, 4000);
    //Mat imgtmp;
    //cvtColor(mySFF.smat,imgtmp,CV_GRAY2BGR);
    //mySFF.ioWizard.showImage("scale", mySFF.imat/2+imgtmp/2, 10000);
    //mySFF.ioWizard.showImage("scale", mySFF.imat, 4000);

    //cout<<mySFF.p_sharpOP<<endl;
    
    
    
    //mySFF.myLog.write();

        
    return 0;

}




