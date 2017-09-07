/****************************
*	ProjSFF
*	energy.cpp
*	cstoribal
*	26-04-17
****************************/

#include "energy2.h"

EnergyClass::EnergyClass(){ }
EnergyClass::~EnergyClass(){}

bool EnergyClass::setlogs(MyLog* mylog){
    this->myLog = mylog;
    return true;
}

bool EnergyClass::set_class(DepthClass* p_depthclass, IOWizard* p_io){
    this->p_depthClass = p_depthclass;
    this->p_ioW = p_io;
    return true;
}

bool EnergyClass::set_parameters(const string & typeD, const string & typeR, const tdf_imgset & sharpStruct, const Mat1T & dmat, fType scale_d, fType scale_r){
    // TODO MODIFIER avec un energy_prts. ici c'est nul.

    this->eParams.typeD = typeD;
    this->eParams.typeR = typeR;
    this->eParams.autoroi = 1;
    this->eParams.dmat  = dmat;      //  needed
    this->eParams.dim1  = dmat.rows; // ?
    this->eParams.dim2  = dmat.cols; 

    this->eParams.scale_d = scale_d; //ok
    this->eParams.scale_r = scale_r; //ok

    this->eParams.rmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    this->eParams.rmatset = "NULL"; //TODO manage intelligent update.
                                    //TODO kill all useless functions.
                                    //TODO add a storage vor vectmat data
                                    //TODO so that it's computed once4all
                                    //      excepted scale_d factor
                                    //      IDEM for scale_r
    /*
    this->flag_Dmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_8U);
    this->flag_Rmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_8U);
    this->eDataij   = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    this->eReguij   = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    this->delta_eDij= Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    this->delta_eRij= Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    this->ed_mat    = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    this->er_mat    = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    this->zeromatrix= Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    this->accumat   = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    
    this->eData=0;
    this->eRegu=0;
    */
    
    if(is_same<eType,int>::value){
        this->scale_to_etype = 10/(this->p_depthClass->getMeanFocusStep()[0]);
        myLog->a("Setting scale to ints : ");
        myLog->a(to_string2(this->scale_to_etype)+"\n"); 
    }
    else{
        this->scale_to_etype = 1;
    }
    
    return true;
}




/*
bool EnergyClass::computeMatEnergy(int stype, const Mat1T & drmat, 
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
*/


bool EnergyClass::getCrossLabelMatrix(const vector<fType> & lvect, Mat1E & lmat){
    if(eParams.typeR == "absdiff") l_absdiff(lvect,lmat);
    if(eParams.typeR == "normeL2") l_normeL2(lvect,lmat);

    if(!l_checkmetric(lmat)){
        myLog->a("regularisation term is not a metric\n");
        COUT("regularisation term is not a metric");
        return false;
    }
    return true;
}

bool EnergyClass::getDataEnergy_3DMatrix(const vector<fType> & lvect, vector<Mat1E> & emat){
    //TODO check initialisation
    //TODO add feature checking if energy>>overflow
    
    
    eParams.typeD = "nL1_Rw1"; //TODO TMP
    vector<Mat1T> e1mat;
    if(true) //NotSet.
    {
    if(eParams.typeD == "absdiff"
     ||eParams.typeD == "normeL1") e_normeL1(lvect,e1mat);
    if(eParams.typeD == "normeL2") e_normeL2(lvect,e1mat);
    if(eParams.typeD == "Nsharp")  e_Nsharp( lvect,e1mat);
    if(eParams.typeD == "nL1_Rw1"
     ||eParams.typeD == "nL2_Rw1") e_nLx_Rw1(lvect,e1mat);
    }
    
    for(int k=0; k<emat.size();k++){
        emat[k] = round_etype(e1mat[k]*this->eParams.scale_d);
    }
    return true;
}

bool EnergyClass::get_pointer_dmat( Mat1T* & matrix){
    matrix = &(this->eParams.dmat);
    return true;
}

/*
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


bool EnergyClass::d_absdiff(const Mat1T & drmat, const Rect & roi){
    Mat1T eDataij_roi = Mat(eDataij,roi);
    Mat1T tmpmat= Mat(drmat,roi) - Mat(eParams.dmat,roi);
    tmpmat = this->eParams.scale_d*(cv::abs(tmpmat));
    tmpmat.copyTo(eDataij_roi,flag_Dmat);
    tmpmat = tmpmat - ed_mat;
    tmpmat.copyTo(delta_eDij,flag_Dmat);
    delta_eD = sum(delta_eDij)[0];
    //done
    
    
    return true;
}

bool EnergyClass::d_normeL2(const Mat1T & drmat, const Rect & roi){
    Mat1T eDataij_roi = Mat(eDataij,roi);
    Mat1T tmpmat= Mat(drmat,roi) - Mat(eParams.dmat,roi);
    tmpmat = this->eParams.scale_d*tmpmat.mul(tmpmat);
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


bool EnergyClass::r_absdiff(const Mat1T & drmat, const Rect & roi){
    Mat1T eReguij_roi = Mat(eReguij,roi);
    Mat1T tmpmat;
    zeromatrix.copyTo(accumat);
    Mat mL = (Mat_<fType>(3,1) << -1, 1, 0);
    Mat mR = (Mat_<fType>(3,1) <<  0, 1,-1);
    Mat mU = (Mat_<fType>(1,3) << -1, 1, 0);
    Mat mD = (Mat_<fType>(1,3) <<  0, 1,-1);
    
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

*/





/////////////////////////////////////////////////////////////////////
//	DNRJ - Lab2Mat wrning - totally ignore previous system.
/////////////////////////////////////////////////////////////////////


bool EnergyClass::e_normeL1(const vector<fType> & lvect, vector<Mat1T> & emat){
    emat = vector<Mat1T>( lvect.size() );
    Mat1T labelmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    Mat1T tmpmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    for(int k=0; k<lvect.size(); k++)
    {
        labelmat = labelmat*0+lvect[k];
        tmpmat = this->scale_to_etype*cv::abs(labelmat-this->eParams.dmat);
        tmpmat.copyTo(emat[k]);
        /*
        Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
        for(int i=0; i<this->eParams.dim1; i++) 
        for(int j=0; j<this->eParams.dim2; j++)
        {
            emat[k].at<eType>(i,j) = round_etype(tmpmat.at<fType>(i,j));
        }
        */
    }

    return true;
}

bool EnergyClass::e_normeL2(const vector<fType> & lvect, vector<Mat1T> & emat){
    emat = vector<Mat1T>( lvect.size() );
    Mat1T labelmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    Mat1T tmpmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    for(int k=0; k<lvect.size(); k++)
    {
        labelmat = labelmat*0+lvect[k];
        tmpmat = this->scale_to_etype
                 *cv::abs(labelmat-this->eParams.dmat);
        tmpmat = tmpmat.mul(tmpmat);
        tmpmat.copyTo(emat[k]);
        /*
        Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
        for(int i=0; i<this->eParams.dim1; i++) 
        for(int j=0; j<this->eParams.dim2; j++)
        {
            emat[k].at<eType>(i,j) = round_etype(tmpmat.at<fType>(i,j));
        }
        */
    }
    return true;
}

bool EnergyClass::e_Nsharp( const vector<fType> & lvect, vector<Mat1T> & emat){
    // Prendre les sharpmatrix (à priori m^ dimensions)
    // trouver maxmax
    // maxmax-sharpmat = okay
    fType maxmax, fmin, fmax;
    vector<Mat1T> sharp3Dvmat;
    Mat1T tmpmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
    emat.resize(lvect.size());
    p_depthClass->getVmatSharpI(sharp3Dvmat);
    minMaxLoc(sharp3Dvmat[0],&fmin,&maxmax);
    for(int k=0; k<lvect.size(); k++)
        {
            minMaxLoc(sharp3Dvmat[k],&fmin,&fmax);
            maxmax= (fmax>maxmax)?fmax:maxmax;           
        }
    for(int k=0; k<lvect.size(); k++)
        {
            tmpmat = (1-sharp3Dvmat[k]/maxmax)*
			this->scale_to_etype*
			lvect.size();
            emat[k] = tmpmat.clone();
        }
    return true;
}


bool EnergyClass::e_nLx_Rw1( const vector<fType> & lvect, vector<Mat1T> & emat){
    // to be scaled, labeldynamic nblabel*(1/scale_to_etype)
    //               sharpnessdynamic (max, min?)
    // L1 norm reweighted with reliability 1 == (Max-Min)/RArea/(MaxMax)
    //						*labels (scaling factor)
    // unit is labels
    if(eParams.rmatset != "rw1") //1step setting reliability matrix
    {
        fType fmin, fmax, maxmax,minmin,auc; //auc area under curve
        vector<Mat1T> sharp3Dvmat;
        p_depthClass->getVmatSharpI(sharp3Dvmat);
        int nblabel=sharp3Dvmat.size();
        
        Mat1T accumulator, maxmat, minmat;
        accumulator = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
        maxmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
        minmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_TF);
        Mat1b flagmat = Mat::zeros(eParams.dim1,eParams.dim2, CV_8U);
        //init min & max with a lower/upper bound
        minMaxLoc(sharp3Dvmat[0],&minmin,&maxmax);
        sharp3Dvmat[0].copyTo(maxmat);
        sharp3Dvmat[0].copyTo(minmat);
        for(int k=0; k<nblabel; k++)
        {
            minMaxLoc(sharp3Dvmat[k],&fmin,&fmax);
            flagmat= (sharp3Dvmat[k]>maxmat);//?sharp3Dvmat[k]:maxmat;
            sharp3Dvmat[k].copyTo(maxmat,flagmat);
            flagmat = (sharp3Dvmat[k]<minmat);//?sharp3Dvmat[k]:minmat;
            sharp3Dvmat[k].copyTo(minmat,flagmat);
            maxmax= (fmax>maxmax)?fmax:maxmax;
            minmin= (fmin<minmin)?fmin:minmin;            
        }

        for(int k=0; k<nblabel; k++)
        {
            accumulator += sharp3Dvmat[k]-minmat;            
        }
        
        Mat1T matrix;
        divide((maxmat-minmat).mul(maxmat-minmat)*nblabel, (accumulator*maxmax), matrix);
        //minMaxLoc(matrix,&minmin,&maxmax);
        //CPING("rmat");CPING2(minmin,maxmax);
        this->eParams.rmat = matrix/mean(matrix)[0];
        this->eParams.rmatset= "rw1";
    }


    if(eParams.typeD == "nL1_Rw1") this->e_normeL1(lvect, emat);
    if(eParams.typeD == "nL2_Rw1") this->e_normeL2(lvect, emat);
    
    for(int k=0; k<emat.size(); k++)
    {
        emat[k] = emat[k].mul(this->eParams.rmat);
    }
    //fType maxmax,minmin;
    //minMaxLoc(this->eParams.rmat,&minmin,&maxmax);
    //CPING("rmat");CPING2(minmin,maxmax);
    return true;
}



////
bool EnergyClass::l_checkmetric(Mat1E & lmat){
    for(int i=0;i<lmat.rows;i++){
        for(int j=0;j<lmat.cols;j++){
            if(i==j && lmat.at<eType>(i,j)!=0){
                
                COUT("Separation1");
                return false;
            }
            if(i!=j && lmat.at<eType>(i,j)==0){
                
                COUT("Separation2");
                return false;
            }
            if(lmat.at<eType>(i,j)<0){
                COUT("positivité");
                return false;
            }
            if( lmat.at<eType>(i,j)!=lmat.at<eType>(i,j) ){
                COUT("symmétrie");
                return false;
            }
            for(int k=0; k<lmat.rows; k++){
                if( lmat.at<eType>(i,j)>
		(lmat.at<eType>(i,k)+lmat.at<eType>(k,j)) ){
                COUT("Inégalité triangulaire");
                printf("delta       %.40f \n",this->eParams.scale_r);
                printf("i,j  %2i,%2i  %.40f \n",i,j,lmat.at<eType>(i,j) );
                printf("i,k  %2i,%2i  %.40f \n",i,k,lmat.at<eType>(i,k) );
                printf("k,j  %2i,%2i  %.40f \n",k,j,lmat.at<eType>(k,j) );
                printf("delta       %.40f \n",lmat.at<eType>(i,j)-lmat.at<eType>(i,k)-lmat.at<eType>(k,j) );
                return false;
                }
            }
        }
    }
    return true;    
}


/////////////////////////////////////////////////////////////////////
//	REG -Lab2Lab
/////////////////////////////////////////////////////////////////////


bool EnergyClass::l_absdiff(const vector<fType> & lvect, Mat1E & lmat){
    //
    lmat = Mat::zeros(lvect.size(),lvect.size(), CV_TE);
    //myLog->a("\n regularisation term matrix \n");
    for(int i=0;i<lvect.size();i++)
    {
    for(int j=0;j<lvect.size();j++)
    {
        lmat.at<eType>(i,j) = this->eParams.scale_r*abs(i-j);
        //lmat.at<eType>(i,j) = round_etype(this->scale_to_etype*this->eParams.scale_r*abs( lvect[i] - lvect[j] ));
        //myLog->a(to_string2(lmat.at<eType>(i,j)) +"; ");
    }
    //myLog->a("\n");
    }
    
    //myLog->a("\n");
    //myLog->write();


    return true;
}


bool EnergyClass::l_normeL2(const vector<fType> & lvect, Mat1E & lmat){
    lmat = Mat::zeros(lvect.size(),lvect.size(), CV_TE);
    for(int i=0;i<lvect.size();i++)    for(int j=0;j<lvect.size();j++){
        lmat.at<eType>(i,j) = round_etype(this->eParams.scale_r*pow(this->scale_to_etype*(lvect[i] - lvect[j]),2));
    }

    return true;
}



/////////////////////////////
//    Utility
/////////////////////////////

eType EnergyClass::round_etype(fType input){
    if(is_same<eType,fType>::value)
        return input;
    if(is_same<eType,int>::value)
        return round(input);
    else
        return static_cast<eType>(input);
}

Mat1E EnergyClass::round_etype(Mat1T input){
    if(is_same<eType,fType>::value)
        return input;

    Mat1E output = Mat::zeros(eParams.dim1,eParams.dim2, CV_TE);

    if(is_same<eType,int>::value)
    {
        for(int i=0; i<input.rows; i++) 
        for(int j=0; j<input.cols; j++)
        {
            output.at<eType>(i,j) = round(input.at<fType>(i,j));
        }
    }
    else
    {
        for(int i=0; i<input.rows; i++) 
        for(int j=0; j<input.cols; j++)
        {
            output.at<eType>(i,j) = static_cast<eType>(input.at<fType>(i,j));
        }
    }
    return output;
}

vector<Mat1E> EnergyClass::round_etype(vector<Mat1T> input){
    vector<Mat1E> output;
    output.resize(input.size());

    if(is_same<eType,fType>::value)
    {
        for(int k=0; k<input.size(); k++){
            output[k]=input[k];
        }
    }
    else
    {
        for(int k=0; k<input.size(); k++){
            output[k]=round_etype(input[k]);
        }
    }
    return output;
}
