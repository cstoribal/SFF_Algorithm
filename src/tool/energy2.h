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

#include "../misc/miscdef.h"

using namespace std;
using namespace cv;

class EnergyClass{
public:
    EnergyClass(); ~EnergyClass();

    bool set_parameters(const string & typeD, const string & typeR, const tdf_imgset & sharpStruct, const Mat1d dmat, int nb_labels, double scale);

    Mat1d  ed_mat;
    Mat1d  er_mat;
    double eData ;
    double eRegu ;


// Ajout d'un paramètre stype : 1=data, 2=Regularisation, 3=both
    //bool computeMatEnergy(int stype, const Mat1d & drmat);
    bool computeMatEnergy(int stype, const Mat1d & drmat, 
		const vector<Point> & P);
    //bool computeMatEnergy(int stype, const Mat1d & drmat);
    bool updateEnergy(int stype);//, const Mat1d & drmat, 
				    // const vector<Point> & P);
    //bool updateSumEnergy(int stype);//, const vector<Point> & P);
    bool getCrossLabelMatrix(const vector<double> & lvect, Mat1d & lmat);
    // update mat just stores exxxxij into .at(i,j) and registers var_e
    // update sum adds var_exij to exxxx 
    // if updateSumEnergy with vector NULL, sum over the whole matrix

    //bool sumEnergyData(int stype); //over the whole matrix

    // ici stype = 1 Data ou 2 Regul
    bool get_P_Qneighbors(int stype, const vector<Point> & P, vector<Point> & Q);
    bool get_P_ROI(int stype, const vector<Point> & P, Rect & mROI);
    bool get_ROI_from_ROI(int stype, const Rect & roiin, Rect & roiout);

    //Ajouter un gestionnaire pour un changement de vector<Point> points
    // nb remplacer exxxxij par vector<double> exxxxij ?
    
    // Penser surcharge d'opérateurs. Non, mettre un if, avec un vecteur Nul
    // (argument optionnel, si nul, alors ROI= toute la matrice)
    // tester la taille pour "optimiser un peu" le calcul.
    // egalement penser en terme de ROI
    // get var_exij ?


    Rect  roi_D,      roi_R;
    Mat1b roi_Dmat,   roi_Rmat;
    Mat1b flag_Dmat,  flag_Rmat; // use ROI to locate points in the flagMAT with the good scale
    Mat1d eDataij,    eReguij; 
    Mat1d delta_eDij, delta_eRij; //nouvelles matrices d'énergie
    Mat1d zeromatrix;
    double delta_eD,  delta_eR;

    Mat1d accumat;

private:
    //typedef_structurespécifique
    tdfp_energy eParams;
    


    
    
/////////////////////////////////////////////////////////////////////
//	DataNRJ
/////////////////////////////////////////////////////////////////////
    bool d_absdiff(const Mat1d & drmat, const Rect & roi);
    bool d_normeL2(const Mat1d & drmat, const Rect & roi);











/////////////////////////////////////////////////////////////////////
//	DataREG
/////////////////////////////////////////////////////////////////////
    bool r_absdiff(const Mat1d & drmat, const Rect & roi);
    
    
    
    
    
    
    
/////////////////////////////////////////////////////////////////////
//	REG - Lab2Lab
/////////////////////////////////////////////////////////////////////
    bool l_absdiff(const vector<double> & lvect, Mat1d & lmat);
    bool l_normeL2(const vector<double> & lvect, Mat1d & lmat);
    
    
    
    
    
    
};




#endif // ENERGY2_H_
