/****************************
*	ProjSFF
*	energy.cpp
*	cstoribal
*	26-04-17
****************************/

#include "energy2.h"

EnergyClass::EnergyClass(){ }
EnergyClass::~EnergyClass(){}



bool EnergyClass::set_parameters(const string & typeD, const string & typeR, const typedef_imgset & sharpStruct, const Mat1d dmat, int nb_labels, double scale){
    // also initialises Matrix flag_Mat (& ROI as part), exxxxij, var_exij
    this->eParams.typeD = typeD;
    this->eParams.typeR = typeR;
    this->eParams.autoroi = 1;
    this->eParams.dmat  = dmat;
    this->eParams.dim1  = dmat.rows;
    this->eParams.dim2  = dmat.cols;
    this->eParams.rmat = this->eParams.dmat.clone(); //TODO
    this->eParams.scale = scale;
    
    this->flag_Dmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_8U);
    this->flag_Rmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_8U);
    this->eDataij   = Mat::zeros(eParams.dim1,eParams.dim2, CV_64F);
    this->eReguij   = Mat::zeros(eParams.dim1,eParams.dim2, CV_64F);
    this->delta_eDij= Mat::zeros(eParams.dim1,eParams.dim2, CV_64F);
    this->delta_eRij= Mat::zeros(eParams.dim1,eParams.dim2, CV_64F);
    this->ed_mat    = Mat::zeros(eParams.dim1,eParams.dim2, CV_64F);
    this->er_mat    = Mat::zeros(eParams.dim1,eParams.dim2, CV_64F);
    this->zeromatrix= Mat::zeros(eParams.dim1,eParams.dim2, CV_64F);
    this->accumat   = Mat::zeros(eParams.dim1,eParams.dim2, CV_64F);
    
    this->eData=0;
    this->eRegu=0;
    return true;
}





bool EnergyClass::computeMatEnergy(int stype, const Mat1d & drmat, 
			const vector<Point> & P){

    if(stype==E_DATA or stype==E_BOTH){
        if(P.size()==0){
            roi_D=Rect(0,0,eParams.dim2,eParams.dim1);
            roi_Dmat=flag_Dmat;
            roi_Dmat=roi_Dmat*0+255;
        }
    //if P.size < seuil de flemme
    // get augmented vector<Point> P to Q 
    // set vector<Point>Q at 1 on flag_Dmat
    // get ROI from vector<Point>Q

    //if P.size > seuil de flemme
    // get ROI from vector<Point>P
    // set ROI at 1
    // get augmented ROI from ROI
    
    //call d_ function with ROI
        else 
        if(eParams.autoroi==1){ // TODO if seuil <> /TODO pour l'instant on ne passe pas par ici
            get_P_ROI(1,P,roi_D);
            roi_Dmat = Mat( flag_Dmat, roi_D );
            roi_Dmat = roi_Dmat*0+255;
            get_ROI_from_ROI(1,roi_D,roi_D);
        }

    //Computeselection
        if(eParams.typeD=="absdiff")d_absdiff(drmat, roi_D);
        if(eParams.typeD=="normeL2")d_normeL2(drmat, roi_D);
    }

    //tmp r_compute.
    if(stype==E_REGU or stype==E_BOTH){
        if(P.size()==0){
            roi_R=Rect(0,0,eParams.dim2,eParams.dim1);
            roi_Rmat=flag_Rmat;
            roi_Rmat=roi_Rmat*0+255;
        } 
        // TODO else compute ROI - not done yet actually
        if(eParams.typeR == "absdiff") r_absdiff(drmat, roi_R);
    }

    return true;

}

bool EnergyClass::updateEnergy(int stype){
    //CPING2(eParams.dmat.rows,eParams.dmat.cols);
    if(stype==E_DATA or stype==E_BOTH){
        ed_mat = ed_mat + delta_eDij;
        //cout << zeromatrix.type() << endl;
        //cout << delta_eDij.type() << endl;
        //cout << flag_Dmat.type() << endl;

        zeromatrix.copyTo(delta_eDij,flag_Dmat);
        zeromatrix.copyTo(eDataij, flag_Dmat);
        flag_Dmat = flag_Dmat<flag_Dmat;
        eData = eData + delta_eD;
        delta_eD = 0;
    }

    if(stype==E_REGU or stype==E_BOTH){
        er_mat = er_mat + delta_eRij;
        zeromatrix.copyTo(delta_eRij,flag_Rmat);
        zeromatrix.copyTo(eReguij, flag_Rmat);
        flag_Rmat = flag_Rmat<flag_Rmat;

        eRegu = eRegu + delta_eR;
        delta_eR = 0;
    }





    return true;
}



bool EnergyClass::getCrossLabelMatrix(const vector<double> & lvect, Mat1d & lmat){
    if(eParams.typeR == "absdiff") l_absdiff(lvect,lmat);
    if(eParams.typeR == "normeL2") l_normeL2(lvect,lmat);
    return true;
}



/////////////////////////////////////////////////////////////////////
///	ROIs
/////////////////////////////////////////////////////////////////////

bool EnergyClass::get_P_ROI(int stype, const vector<Point> & P, Rect & mROI){
    //tmp
    int xmax,xmin,ymax,ymin;
    xmax = P[0].x;xmin=xmax;ymax=P[0].y;ymin=ymax;
    for(int i=0;i<P.size();i++){
        xmax = (xmax<P[i].x)? P[i].x:xmax;
        xmin = (xmin>P[i].x)? P[i].x:xmin;
        ymax = (ymax<P[i].y)? P[i].y:ymax;
        ymin = (ymin>P[i].y)? P[i].y:ymin;
    }
    mROI = Rect(xmin,ymin, xmax-xmin, ymax-ymin);
    
    return true;
}

bool EnergyClass::get_ROI_from_ROI(int stype, const Rect & roiin, Rect & roiout){
    int w=5;
    int h=5;

    if(stype == 1){
        
        
    }
    else if(stype==2)
    {
        
        
    }
    roiout = roiin - Point(w/2,h/2);
    roiout = roiout + Size(w,h);
    cout << "ROI is" << endl;
    cout << roiout.x << endl;
    cout << roiout.y << endl;
    cout << roiout.height << endl;
    cout << roiout.width << endl;
    // Limit ROI to bounds
    return true;
}
    

/////////////////////////////////////////////////////////////////////
//	DataNRJ
/////////////////////////////////////////////////////////////////////


bool EnergyClass::d_absdiff(const Mat1d & drmat, const Rect & roi){
    Mat1d eDataij_roi = Mat(eDataij,roi);
    Mat1d tmpmat= Mat(drmat,roi) - Mat(eParams.dmat,roi);
    tmpmat = cv::abs(tmpmat);
    tmpmat.copyTo(eDataij_roi,flag_Dmat);
    tmpmat = tmpmat - ed_mat;
    tmpmat.copyTo(delta_eDij,flag_Dmat);
    delta_eD = sum(delta_eDij)[0];
    //done
    
    
    return true;
}

bool EnergyClass::d_normeL2(const Mat1d & drmat, const Rect & roi){
    Mat1d eDataij_roi = Mat(eDataij,roi);
    Mat1d tmpmat= Mat(drmat,roi) - Mat(eParams.dmat,roi);
    tmpmat = tmpmat.mul(tmpmat);
    tmpmat.copyTo(eDataij_roi,flag_Dmat);
    tmpmat = tmpmat - ed_mat;
    tmpmat.copyTo(delta_eDij,flag_Dmat);
    delta_eD = sum(delta_eDij)[0];
    //done
    
    
    return true;
}





/////////////////////////////////////////////////////////////////////
//	ReguNRJ
/////////////////////////////////////////////////////////////////////


bool EnergyClass::r_absdiff(const Mat1d & drmat, const Rect & roi){
    Mat1d eReguij_roi = Mat(eReguij,roi);
    Mat1d tmpmat;
    zeromatrix.copyTo(accumat);
    Mat mL = (Mat_<double>(3,1) << -1, 1, 0);
    Mat mR = (Mat_<double>(3,1) <<  0, 1,-1);
    Mat mU = (Mat_<double>(1,3) << -1, 1, 0);
    Mat mD = (Mat_<double>(1,3) <<  0, 1,-1);
    
    filter2D(drmat,tmpmat,-1,mL, Point(-1,-1), 0, BORDER_REPLICATE);
    accumat = abs(tmpmat);
    filter2D(drmat,tmpmat,-1,mR, Point(-1,-1), 0, BORDER_REPLICATE);
    accumat += abs(tmpmat);
    filter2D(drmat,tmpmat,-1,mU, Point(-1,-1), 0, BORDER_REPLICATE);
    accumat += abs(tmpmat);
    filter2D(drmat,tmpmat,-1,mD, Point(-1,-1), 0, BORDER_REPLICATE);
    accumat += abs(tmpmat);
    accumat -= er_mat;
    accumat.copyTo(delta_eRij,flag_Rmat);
    delta_eR = sum(delta_eRij)[0];
    //
    
    return true;    
    
}



/////////////////////////////////////////////////////////////////////
//	REG -Lab2Lab
/////////////////////////////////////////////////////////////////////


bool EnergyClass::l_absdiff(const vector<double> & lvect, Mat1d & lmat){
    //
    lmat = Mat::zeros(lvect.size(),lvect.size(), CV_64F);
    for(int i=0;i<lvect.size();i++){
    for(int j=0;j<lvect.size();j++){
        lmat.at<double>(i,j) = this->eParams.scale*abs( lvect[i] - lvect[j] );
    }}

    return true;
}


bool EnergyClass::l_normeL2(const vector<double> & lvect, Mat1d & lmat){
    lmat = Mat::zeros(lvect.size(),lvect.size(), CV_64F);
    for(int i=0;i<lvect.size();i++){
    for(int j=0;j<lvect.size();j++){
        lmat.at<double>(i,j) = this->eParams.scale*pow(lvect[i] - lvect[j],2);
    }}

    return true;
}

