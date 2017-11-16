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

bool Sharpness_Operator::setlogs(MyLog* mylog){
    this->myLog = mylog;
    return true;
}

bool Sharpness_Operator::optypeSelector(string type){
    this->type = type;
    return true;
}

bool Sharpness_Operator::computeOp(string optype, const tdf_imgset & iset, tdf_imgset& sset){
    this->height = iset[0].ivmat[0].rows;
    this->width  = iset[0].ivmat[0].cols;
    sset.resize(iset.size());
    for(int i=0;i<iset.size();i++)
    {   
        sset[i].name = iset[i].name;
        sset[i].rank = iset[i].rank;
        sset[i].focus= iset[i].focus;
        sset[i].dpth = iset[i].dpth;
        if(optype=="SMLAP"
          |optype=="3DLAP") compute_SMLAP(iset[i].ivmat, sset[i].ivmat);
        if(optype=="STA2")  compute_STA2( iset[i].ivmat, sset[i].ivmat);
    }
    
    if(optype=="3DLAP") compute_3DLAP(sset);
    
    
    
    
    return true;
}

bool Sharpness_Operator::compute(tdf_imgset iset, tdf_imgset& sset){
    computeOp(this->type, iset, sset);
    return true;
}


bool Sharpness_Operator::compute_SMLAP(vector<Mat1T> ivmat, vector<Mat1T>& smat){
    // Build temporary modlaplacian
    // Sum Laplacians
    // over !

    Mat lap, lapb, lapc;
    vector<Mat1T> svmat;
    Mat Lx = (Mat_<fType>(3,1) << -1, 2, -1);
    Mat Ly = (Mat_<fType>(1,3) << -1, 2, -1);
    
    for(int i=0; i<ivmat.size(); i++)
    {
        filter2D(ivmat[i],lap,-1,Lx, Point(-1,-1), 0, BORDER_REPLICATE);
        filter2D(ivmat[i],lapb,-1,Ly, Point(-1,-1), 0, BORDER_REPLICATE);
    
        lap = abs(lap)+abs(lapb);

        // (lap*0).copyTo(lap,lap<seuil);
        
        lapb = Mat::ones(7, 7, CV_TF); //TODO define window dimensions
        filter2D(lap,lapc,-1,lapb,Point(-1,-1),0,BORDER_REPLICATE);
        
        svmat.push_back(lapc);
    }
    
    combinRule(svmat,smat);	//TODO Select a rule for
				//combining different Laplacians

    fType seuil = 0;//1;			//TODO fixer seuil
    smat[0] = smat[0].setTo(0,smat[0]<seuil);
    return true;
}

bool Sharpness_Operator::compute_3DLAP(tdf_imgset & ioset){
    // Ayea
    int indsize = ioset.size();
    tdf_imgset tmpset;
    tmpset=ioset;
    tmpset[0].ivmat[0]=(2*ioset[0].ivmat[0]+ioset[1].ivmat[0])/3;
    tmpset[indsize-1].ivmat[0]=(2*ioset[indsize-1].ivmat[0]+ioset[indsize-2].ivmat[0])/1;

    for(int i=1; i<indsize-1; i++)
    {
        tmpset[i].ivmat[0] = (ioset[i-1].ivmat[0] + ioset[i].ivmat[0]
						 + ioset[i+1].ivmat[0])/1;
    }
    ioset = tmpset;
    return true;
}

bool Sharpness_Operator::compute_STA2(vector<Mat1T> ivmat, vector<Mat1T>& smat){
    // 
    int winsize = 10;
    for(int i=0; i<this->height;i++) for(int j=0; j<this->width; j++)
    {
        Rect roi = Rect( Point((0>j-10)?0:(j-10),(0>i-10)?0:(i-10)),
	Point(  (this->width-1 <=j+10)?this->width-1 :(j+10),
		(this->height-1<=i+10)?this->height-1:(i+10)) );
        
        
    }
    
    
    
    
    
    return true;
}


bool Sharpness_Operator::combinRule(vector<Mat1T> in_svmat, vector<Mat1T>& out_svmat){
    vector<Mat1T> smattmp;
    smattmp.push_back(Mat::zeros(in_svmat[0].rows,in_svmat[0].cols,CV_TF));
    for(int i=0; i<in_svmat.size(); i++)
    {
        smattmp[0] += in_svmat[i];
    }
    out_svmat = smattmp;
    return true;
}
