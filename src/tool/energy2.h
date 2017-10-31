/****************************
*	ProjSFF
*	energy.h
*	cstoribal
*	26-04-17
****************************/

/*
Computes both energy associated with data terms
and to regularisation 
*/

#ifndef ENERGY2_H_
#define ENERGY2_H_

#define E_DATA 1
#define E_REGU 2
#define E_BOTH 3

#include <iostream>
#include <string>
#include <sstream>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>

#include <opencv2/imgproc.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fstream>
#include <vector>
#include <cmath>

#include "../tool/depthEstimator.h"
#include "../io/logs.h"
#include "../misc/miscdef.h"
#include "../io/IOWizard.h"

using namespace std;
using namespace cv;

class EnergyClass{
public:
    EnergyClass(); ~EnergyClass();

    bool setlogs(MyLog* mylog);
    bool set_class(DepthClass* p_depthclass, IOWizard* p_io);
    bool set_parameters(const string & typeD, const string & typeR, const tdf_imgset & sharpStruct, const Mat1T & dmat, fType scale_d, fType scale_r);

    //Mat1T  ed_mat ;
    //Mat1T  er_mat ;
    //fType  eData  ;
    //fType  eRegu  ;

/*
// Ajout d'un paramètre stype : stype :: 1=data, 2=Regularisation, 3=both
    // set roi with vector<Point>, prendre vector<Point>() pour roi = img
    bool computeMatEnergy(int stype, const Mat1T & drmat, 
		const vector<Point> & P);;
    bool updateEnergy(int stype);
    // update mat just stores exxxxij into .at(i,j) and registers var_e

    // ici stype = 1 Data ou 2 Regul
    bool get_P_Qneighbors(int stype, const vector<Point> & P, vector<Point> & Q);
    bool get_P_ROI(int stype, const vector<Point> & P, Rect & mROI);
    bool get_ROI_from_ROI(int stype, const Rect & roiin, Rect & roiout);
*/
    
    bool getCrossLabelMatrix(const vector<fType> & lvect, Mat1E & lmat);
    bool getDataEnergy_3DMatrix(const vector<fType> & lvect, vector<Mat1E> & emat);
    bool get_pointer_dmat( Mat1T* & matrix); 

/*
    Rect  roi_D,      roi_R;
    Mat1b roi_Dmat,   roi_Rmat;
    Mat1b flag_Dmat,  flag_Rmat; // use ROI to locate points in the flagMAT with the good scale
    Mat1T eDataij,    eReguij; 
    Mat1T delta_eDij, delta_eRij; //nouvelles matrices d'énergie
    Mat1T zeromatrix;
    fType delta_eD,  delta_eR;

    Mat1T accumat;
*/
private:
    //typedef_structurespécifique
    tdfp_energy eParams;
    
    MyLog* myLog;
    DepthClass* p_depthClass;
    IOWizard* p_ioW;

    fType scale_to_etype;


    
    
/////////////////////////////////////////////////////////////////////
//	DataNRJ
/////////////////////////////////////////////////////////////////////
/*
    bool d_absdiff(const Mat1T & drmat, const Rect & roi);
    bool d_normeL2(const Mat1T & drmat, const Rect & roi);



/////////////////////////////////////////////////////////////////////
//	DataREG
/////////////////////////////////////////////////////////////////////
    bool r_absdiff(const Mat1T & drmat, const Rect & roi);
*/

 
/////////////////////////////////////////////////////////////////////
//	DNRJ - Lab2Mat wrning - totally ignore previous system.
/////////////////////////////////////////////////////////////////////
    bool e_normeL1(const vector<fType> & lvect, vector<Mat1T> & emat);
    bool e_normeL2(const vector<fType> & lvect, vector<Mat1T> & emat);
    bool e_Nsharp( const vector<fType> & lvect, vector<Mat1T> & emat);
    bool e_nLx_Rw1(const vector<fType> & lvect, vector<Mat1T> & emat);
    
/////////////////////////////////////////////////////////////////////
//	REG - Lab2Lab
/////////////////////////////////////////////////////////////////////
    bool l_checkmetric(Mat1E & lmat);
    bool l_absdiff(const vector<fType> & lvect, Mat1E & lmat);
    bool l_normeL2(const vector<fType> & lvect, Mat1E & lmat);
    
    
    
///////////////////
//      Utility
///////////////////
    eType round_etype(fType input);
    Mat1E round_etype(Mat1T input);
    vector<Mat1E> round_etype(vector<Mat1T> input);
    
    
};




#endif // ENERGY2_H_
