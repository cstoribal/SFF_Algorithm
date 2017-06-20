/****************************
*	ProjSFF
*	sharpClass.cpp
*	cstoribal
*	10-04-17
****************************/

/*
Computing sharpness matrix
*/



#include "sharpnessOperator.h"

using namespace std;
using namespace cv;



Sharpness_Operator::Sharpness_Operator() {

}

Sharpness_Operator::~Sharpness_Operator(){

}

bool Sharpness_Operator::optypeSelector(string type){
    this->type = type;
    return true;
}

bool Sharpness_Operator::computeOp(string optype, const tdf_imgset & iset, tdf_imgset& sset){
    tdf_imgset tmpset;
    for(int i=0;i<iset.size();i++)
    {   
        struct tagged_img tmpimg;
        tmpset.push_back(tmpimg);
        tmpset[i].name = iset[i].name;
        tmpset[i].rank = iset[i].rank;
        tmpset[i].focus= iset[i].focus;
        tmpset[i].dpth = iset[i].dpth;
        if(optype=="SMLAP") compute_SMLAP(iset[i].ivmat, tmpset[i].ivmat);
    }

    sset = tmpset;
    return true;
}

bool Sharpness_Operator::compute(tdf_imgset iset, tdf_imgset& sset){
    computeOp(this->type, iset, sset);
    return true;
}


bool Sharpness_Operator::compute_SMLAP(vector<Mat1d> ivmat, vector<Mat1d>& smat){
    // Build temporary modlaplacian
    // Sum Laplacians
    // over !

    Mat lap, lapb, lapc;
    vector<Mat1d> svmat;
    Mat Lx = (Mat_<double>(3,1) << -1, 2, -1);
    Mat Ly = (Mat_<double>(1,3) << -1, 2, -1);
    
    for(int i=0; i<ivmat.size(); i++)
    {
        filter2D(ivmat[i],lap,-1,Lx, Point(-1,-1), 0, BORDER_REPLICATE);
        filter2D(ivmat[i],lapb,-1,Ly, Point(-1,-1), 0, BORDER_REPLICATE);
    
        lap = abs(lap)+abs(lapb);

        // (lap*0).copyTo(lap,lap<seuil);
        
        lapb = Mat::ones(6, 6, CV_64F); //TODO define window dimensions
        filter2D(lap,lapc,-1,lapb,Point(-1,-1),0,BORDER_REPLICATE);
        
        svmat.push_back(lapc);
    }
    
    combinRule(svmat,smat);	//TODO Select a rule for
					//combining different Laplacians

    double seuil = 0;//1;			//TODO fixer seuil
    smat[0] = smat[0].setTo(0,smat[0]<seuil);
    return true;
}


bool Sharpness_Operator::combinRule(vector<Mat1d> in_svmat, vector<Mat1d>& out_svmat){
    vector<Mat1d> smattmp;
    smattmp.push_back(Mat::zeros(in_svmat[0].rows,in_svmat[0].cols,CV_64F));
    for(int i=0; i<in_svmat.size(); i++)
    {
        smattmp[0] += in_svmat[i];
    }
    out_svmat = smattmp;
    return true;
}
