/****************************
*	ProjSFF
*	depthEstimator.cpp
*	cstoribal
*	18-04-17
****************************

/*
this class computes the local approximation of depth.
We first have to read sharpness measures from an imageset,
then various interpolation methods can be used,
in order to determine the estimated depth of each pixel separately.
Various optional features as reliability estimation can be performed.
nb: different indices of reliability may coexist in the SFF algorithm.
*/


/*
HowTo : implement various interpolation methods
- we need a structure to store polynomial coefficients if used
- we need to store gaussian parameters (scale, variance, mean...)
- we should then have a structure in which one parameter ('type')
  gives us the kind of parameters used to store data esimated.

Prerequisite : polymorphism ?

Data type input, data type output ?
- seems like input will simply be the image sharpness
- intermediate output is likely to be a matrix of vectors -> new typedef, matrice de paramètres.
- we'll output a float matrix with estimated depth ?
- existence of a reliability matrix is not excluded.
*/




#include "depthEstimator.h"

using namespace std;




////////////////////////////////////////////////////
// PUBLIC //
////////////////////////////////////////////////////




DepthClass::DepthClass():set(false){
    this->vmat_sharp_i_set = false;
}

DepthClass::~DepthClass(){}

bool DepthClass::setlogs(MyLog* mylog, IOWizard* _myIo){
    this->myIo  = _myIo;
    this->myLog = mylog;
    return true;
}

bool DepthClass::set_param(string settype){
    this->type = settype;
    this->degree = 8;
    this->oversampling = 1;
    this->set = true;
    myLog->av("depth detection is "+to_string2(settype)+"\n");
    myLog->av("  degree set to "+to_string2(degree)+"\n");
    return true;
}

bool DepthClass::buildEstimation(void){
    return true;
}

bool DepthClass::buildEstimation(const tdf_imgset & sharpSet, tdfp_depth & pmat){
    this->set_dim[0] = sharpSet[0].ivmat[0].rows;
    this->set_dim[1] = sharpSet[0].ivmat[0].cols;
    this->set_dim[2] = sharpSet.size();

    this->sharpSetStored = sharpSet;

    for(int k=0;k<set_dim[2];k++)
    {
        this->focus.push_back(sharpSet[k].focus);
    }
    if(type=="polynome"||type=="polymod")
        f_poly(sharpSet, pmat);
    if(type=="argmax")
        f_argmax(sharpSet, pmat);
    if(type=="gauss")
        f_gauss(sharpSet, pmat);
    return true;
}

bool DepthClass::buildDepthmat(const tdfp_depth & dparam, Mat1T & dmat, Mat1T & dmat_rank, Mat1T & dmat_score){
    cv::Mat1i dmat_label;

    if(type=="polynome"||type=="polymod") 
			d_poly(dparam, dmat, dmat_rank, dmat_score, dmat_label);
    if(type=="argmax")  d_argmax(dparam, dmat, dmat_rank, dmat_score);
    if(type=="gauss")   d_gauss(dparam, dmat, dmat_rank, dmat_score, dmat_label);



    if(type=="polymod") d_polymod(dmat, dmat_rank, dmat_score, dmat_label);
    cout << "depthmap built" << endl;
    return true;
}


bool DepthClass::getVmatSharpI(vector<Mat1T> & vmat){
    if(!vmat_sharp_i_set){
        myLog->a("getVmatSharpI failed, depthmap not set");
        return false;
    }
    
    vmat = this->vmat_sharp_i;
    return true;
}

int DepthClass::getRankFromDepth(fType input){
    int rank = -1;
    int i=0;
    while(rank == -1 and i<DepthToRank.size() ){
        if(HalfDepthVect[i]>=input)
            rank = (int)round(DepthToRank[i][1]);
        i++;
    }
    if(rank == -1)
    {
        //cout << "could not find rank, value set to last image" << endl;
        rank = (int)round(DepthToRank[DepthToRank.size()-1][1]);
    }
    return rank;
}

bool DepthClass::getMRankFromMDepth(const Mat1T& minput, cv::Mat1i & moutput){
    // renvoie la matrice des rangs (parmi les images d'origine)
    moutput = cv::Mat::zeros(minput.rows,minput.cols,CV_32S);
    cv::Mat1i m_zero = cv::Mat::zeros(minput.rows,minput.cols,CV_32S);
    cv::Mat1i m_hold = cv::Mat::zeros(minput.rows,minput.cols,CV_32S);
    cv::Mat1i m_mask = cv::Mat::zeros(minput.rows,minput.cols,CV_32S);
    fType tempvalA,tempvalB;
    int   temprank;
    for(int i=1; i<this->DepthToRank.size();i++){
        m_mask = m_mask*0;
        //TODO halfdepth
        temprank=this->DepthToRank[i][1];
        tempvalA=this->HalfDepthVect[i-1];//[0];
        tempvalB=this->HalfDepthVect[i];//[0];
        m_hold=m_zero+temprank;
        m_hold.copyTo(m_mask, (minput<=tempvalB) );
        m_zero.copyTo(m_mask, (minput<=tempvalA) );
        m_mask.copyTo(moutput,(m_mask>0));
    }
    return true;
}

bool DepthClass::getMLabelFromMDepth(const Mat1T& minput, cv::Mat1i & moutput){
    // renvoie la matrice des labels (parmi les labels interpolés)
    moutput = cv::Mat::zeros(minput.rows,minput.cols,CV_32S);
    cv::Mat1i m_zero = cv::Mat::zeros(minput.rows,minput.cols,CV_32S);
    cv::Mat1i m_hold = cv::Mat::zeros(minput.rows,minput.cols,CV_32S);
    cv::Mat1i m_mask = cv::Mat::zeros(minput.rows,minput.cols,CV_32S);
    fType tempvalA,tempvalB;
    //int   temprank;
    for(int i=1; i<this->DepthToRank.size();i++){
        m_mask = m_mask*0;
        //temprank=this->DepthToRank[i][1];
        // TODO Halfdepth
        tempvalA=this->HalfDepthVect[i-1];//[0];
        tempvalB=this->HalfDepthVect[i];//[0];
        m_hold=m_zero+i;
        m_hold.copyTo(m_mask, (minput<=tempvalB) );
        m_zero.copyTo(m_mask, (minput<=tempvalA) );
        m_mask.copyTo(moutput,(m_mask>0));
    }
    return true;
}

int DepthClass::getNbLabels(void){
    return DepthToRank.size();
}

vector<fType> DepthClass::getLabels(void){
    vector<fType> A(this->DepthToRank.size());
    for(int i=0; i<this->DepthToRank.size(); i++){
        A[i] = this->DepthToRank[i][0];
    }
    return A;
}

fType DepthClass::getOversampling(void){
    return this->oversampling;
}

vector<fType> DepthClass::getMeanFocusStep(void){ // how to check DTR is built
    vector<fType> mean(1);
    fType accumul = 0;
    int k=0;
    for(int i = 1; i<DepthToRank.size(); i++){
        accumul += DepthToRank[i][0]-DepthToRank[i-1][0];
        k++;
    }
    mean[0] = accumul/k;
    return mean;
}

vector<fType> DepthClass::getMeanLogFocusStep(void){
    //fit  Bexp(Ax)
    vector<fType> ab; ab.resize(2);
    int k=0;
    fType accumul=0;
    for(int i = 1; i<DepthToRank.size(); i++){
        accumul += log(DepthToRank[i][0])-log(DepthToRank[i-1][0]);
        k++;
    }
    ab[0]=accumul/k;
    accumul=0;
    for(int i = 0; i<DepthToRank.size(); i++){
        accumul += log(DepthToRank[i][0])-i*ab[0];
    }
    ab[1] = accumul/DepthToRank.size();
    return ab;
}


bool DepthClass::showInterpolationAt(const vector<cv::Point> & vP, const tdf_imgset & SharpSet, const tdfp_depth & dparam, const Mat1T & dmat, const string & folder){
    if(type=="polynome"||type=="polymod") s_poly_ij(vP, SharpSet, dparam, dmat, folder);
    if(type=="argmax")   s_argmax_ij(vP, SharpSet, dparam, dmat, folder);
    if(type=="gauss")   show_interpol_generic(vP, SharpSet, dparam, dmat, folder);
    
    return true ;
}





//////////////////////////
// PRIVATE //
//////////////////////////

bool DepthClass::build_HalfDepthToRank(void){ //>TODO check indice
    HalfDepthVect.resize(DepthToRank.size());
    for(int i=0; i<DepthToRank.size()-1; i++){
        HalfDepthVect[i] = (DepthToRank[i][0]+DepthToRank[i+1][0])/2.0f ;
    }
    HalfDepthVect[DepthToRank.size()-1]=DepthToRank[DepthToRank.size()-1][0];
    return true;
}



bool DepthClass::f_poly(const tdf_imgset & sharpSet, tdfp_depth & dparam){
    // gives the parameters of each i*j polynom (coefficients) from 
    if(sharpSet.size() == 0) return false;
    int N = set_dim[2];
    int dim[] = {set_dim[0],set_dim[1],degree+1};
    //vector<Mat1d> tmpmat(degree+1);
    vector<fType> x(N); 
    vector<fType> y(N);
    vector<fType> z(degree+1);
    dparam.vmat = vector<Mat1T>(degree+1);
    for(int k=0;k<degree+1;k++){
        dparam.vmat[k] = cv::Mat::zeros(set_dim[0],set_dim[1], CV_TF);}
    for(int k=0;k<N; k++)
        {x[k]=sharpSet[k].focus;}
    

    vector<fType> X(2*degree+1);
    for(int i=0;i<2*degree+1;i++)
    {
        X[i]=0;
        for(int j=0;j<N;j++)
            X[i]=X[i]+pow(x[j],i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    }
    
    
    for(int i=0; i<set_dim[0]; i++){
        for(int j=0; j<set_dim[1]; j++){

            for(int k=0; k<set_dim[2]; k++){
                y[k]=sharpSet[k].ivmat[0].at<fType>(i,j);
            }
            this->interpolate(x,y,degree,set_dim[2],X,z);
            
            for(int k=0; k<degree+1; k++){
                dparam.vmat[k].at<fType>(i,j) = z[k];
            }
        }
    }
    dparam.type   = this->type;
    dparam.degree = this->degree;
    
    return true;
}

bool DepthClass::f_argmax(const tdf_imgset & sharpSet, tdfp_depth & dparam){
    if(sharpSet.size()==0) return false;
    int N = set_dim[2];
    int dim[] = {set_dim[0],set_dim[1],set_dim[2]};

    dparam.vmat = vector<Mat1T>(dim[2]);
    for(int k=0;k<dim[2];k++){
        dparam.vmat[k] = cv::Mat::zeros(set_dim[0],set_dim[1], CV_TF);
    }

    for(int i=0; i<set_dim[0]; i++){
        for(int j=0; j<set_dim[1]; j++){
            for(int k=0; k<dim[2]; k++){
                dparam.vmat[k].at<fType>(i,j) = sharpSet[k].ivmat[0].at<fType>(i,j);
            }
        }
    }
    dparam.type   = this->type;
    dparam.degree = this->set_dim[2];
    return true;
}

bool DepthClass::f_gauss(const tdf_imgset & sharpSet, tdfp_depth & dparam){
    if(sharpSet.size()==0) return false;
    f_argmax(sharpSet,dparam);
    return true;
}



bool DepthClass::build_DR(void){
    fType inv_oversampling = 1.0f/(fType)oversampling;
    vector<vector<fType> > DR( (set_dim[2]-1)*oversampling+1, vector<fType>(2));
    
    for(int k = 0; k<set_dim[2]-1; k++){
    for(int j = 0; j<oversampling; j++){
        DR[oversampling*k+j][0] = (fType)(j*focus[k+1]+(oversampling-j)*focus[k])*inv_oversampling;
        if(j<=oversampling/2) DR[oversampling*k+j][1] = k;
        else DR[oversampling*k+j][1] = k+1;
    }}
    DR[oversampling*(set_dim[2]-1)+0][0] = (fType)focus[ (set_dim[2]-1) ];
    DR[oversampling*(set_dim[2]-1)+0][1] = (fType) (set_dim[2]-1) ;

    this->DepthToRank = DR;
    build_HalfDepthToRank();
    CPING2("DR SET : ", this->DepthToRank.size());
    return true;
}


bool DepthClass::get_dmats_from_sharpset(Mat1T & dmat, Mat1T & dmat_rank, Mat1T & dmat_score, cv::Mat1i & dmat_label){
    if(!vmat_sharp_i_set){return false;}

    dmat = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    dmat = dmat - 20.0;
    dmat_score = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    dmat_rank = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    dmat_label = cv::Mat::zeros(set_dim[0],set_dim[1],CV_32S);
    
    Mat1T matarg = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    Mat1T matrnk = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    Mat1T mat1T_ones = cv::Mat::ones(set_dim[0],set_dim[1],CV_TF);
    cv::Mat1b mat_mask = cv::Mat::zeros(set_dim[0],set_dim[1],CV_8U);
    cv::Mat1i mat_ones = cv::Mat::ones(set_dim[0],set_dim[1],CV_32S);
    cv::Mat1i mat_labels = cv::Mat::ones(set_dim[0],set_dim[1],CV_32S);

    vmat_sharp_i[0].copyTo(dmat_score);
    dmat_label = dmat_label*0.0f;
    dmat       = dmat*0.0f+DepthToRank[0][0];
    for(int k=0; k<set_dim[2]-1; k++){
        for(int f=1; f<=oversampling;f++){
            mat_labels = mat_ones*(k*oversampling+f);

            mat_mask = ( vmat_sharp_i[k*oversampling+f]>dmat_score);
            vmat_sharp_i[k*oversampling+f].copyTo(dmat_score,mat_mask);
            mat_labels.copyTo(dmat_label,mat_mask);
            
            matarg = mat1T_ones * DepthToRank[k*oversampling+f][0];
            matarg.copyTo(dmat,mat_mask);
            matrnk = mat1T_ones * DepthToRank[k*oversampling+f][1];
            matrnk.copyTo(dmat_rank, mat_mask);
        }
    }
    CPING("echo world");
    return true;
}




bool DepthClass::d_gauss(const tdfp_depth & dparam, Mat1T & dmat, Mat1T & dmat_rank, Mat1T & dmat_score, cv::Mat1i & dmat_label){
    
    //dmat_label = cv::Mat::zeros(set_dim[0],set_dim[1],CV_32S);

    int lbar = oversampling*(set_dim[2]-1)+1;
    vector<Mat1T> vectmat(lbar);
    for(int i=0; i<lbar; i++)
    {
        vectmat[i] = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    }

    //SETTING DR
    build_DR();
    /// Vectmat

    
    vector<fType> v_h(lbar);
    fType _fi = 0.0f;
    std::generate(v_h.begin(), v_h.end(), [&_fi] {return _fi++;});
    Mat1T mat_h   = Mat1T(1,       lbar, v_h.data()).clone();
    Mat1T mat_v   = Mat1T(set_dim[2], 1, v_h.data()).clone();
    _fi = 1.0f;
    fType _fincrement = 0.2f;//12.0f/(fType)set_dim[2];
    std::generate(v_h.begin(), v_h.end(), [&_fi,_fincrement]
        {_fi+=_fincrement; return _fi;});
    Mat1T mat_sig = Mat1T(set_dim[2], 1, v_h.data()).clone();
    Mat1T mat_H, mat_V;  mat_h.copyTo(mat_H);  mat_v.copyTo(mat_V);
    Mat1T mat_SIG;       mat_sig.copyTo(mat_SIG);
    Mat1T mat_ones = cv::Mat::ones(set_dim[2],lbar,CV_TF);
    for(int i=0; i<set_dim[2]-1; i++){
        vconcat(mat_H, mat_h, mat_H);
    }
    for(int i=0; i<lbar-1; i++){
        hconcat(mat_V,   mat_v,   mat_V);
        hconcat(mat_SIG, mat_sig, mat_SIG);
    }
    Mat1T mat_HV = cv::Mat::zeros(set_dim[2],lbar,CV_TF);
    mat_HV  = cv::abs(mat_H-mat_V);
    myIo->img_unsetscale();
    fType inv_sqrt2 = 1.0/(fType)sqrt(2);
    mat_SIG = mat_SIG*inv_sqrt2;
    mat_HV  = mat_HV/mat_SIG;
    mat_HV  = -mat_HV.mul(mat_HV);
    cv::exp(mat_HV,mat_HV);
    mat_SIG = (mat_ones / mat_SIG) * (1.0f/(M_PI*2.0f)) ;
    mat_HV  = mat_HV.mul(mat_SIG);
    //mod mat_HV
    fType accu=0.0f;
    for(int j=0; j<lbar; j++){
        accu += mat_HV.at<fType>(set_dim[2]-1,j);
        mat_HV.at<fType>(set_dim[2]-1,j) = accu;
    }
    
    vector<fType> v_coeffs(set_dim[2]);
    vector<fType> v_sum(lbar);
    Mat1T         m_coeffs;
    Mat1T         M_coeffs;

    //for(int i=0; i<lbar; i++)for(int j=0; j<set_dim[1]; j++){
    //    CPING2("mat_HV",mat_HV.at<fType>(i,j));
    //}
    
    for(int i=0; i<set_dim[0]; i++) for(int j=0; j<set_dim[1]; j++){
        for(int k=0; k<set_dim[2]; k++){
            v_coeffs[k]=dparam.vmat[k].at<fType>(i,j);
        }
        m_coeffs = Mat1T(set_dim[2], 1, v_coeffs.data()).clone();
        m_coeffs.copyTo(M_coeffs);
        for(int i=0; i<lbar-1; i++){
            hconcat(M_coeffs, m_coeffs, M_coeffs);
        }
        M_coeffs = M_coeffs.mul(mat_HV);
        cv::reduce(M_coeffs, v_sum, 0, CV_REDUCE_SUM, CV_TF);
        for(int k=0; k<lbar; k++){
            vectmat[k].at<fType>(i,j)=v_sum[k];
        }
    }
    
    if(!vmat_sharp_i_set)
    {
        this->vmat_sharp_i = vectmat;
        vmat_sharp_i_set=true;
    }
    get_dmats_from_sharpset(dmat, dmat_rank, dmat_score, dmat_label);
    
    return true;
}

bool DepthClass::d_argmax(const tdfp_depth & dparam, Mat1T & dmat, Mat1T & dmat_rank, Mat1T & dmat_score){
    // Créer un vecteur des points en lesquels évaluer chaque profondeur
    // D abscisses des profondeurs
    // taille degree*oversampling*set_dim[2]

    //int dim[] = {set_dim[0],set_dim[1]};

    cv::Mat1i dmat_label = cv::Mat::zeros(set_dim[0],set_dim[1],CV_32S);
    vector<Mat1T> vectmat;
    fType inv_oversampling = 1.0f/(fType)oversampling;
    vectmat.resize(oversampling*(set_dim[2]-1)+1);
    for(int i=0; i<oversampling*(set_dim[2]-1)+1; i++)
    {
        vectmat[i] = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    }
    
    dmat = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    dmat = dmat + 20.0;
    dmat_score = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    dmat_rank = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);

    //tmpmat = cv::Mat1d::Mat(set_dim[0],set_dim[1],CV_TF);
    CPING2(dmat.rows,dmat.cols);
   // vector<fType> X( (set_dim[2]-1)*oversampling*degree);
    vector<vector<fType> > DR( (set_dim[2]-1)*oversampling+1, vector<fType>(2));
    
    for(int k = 0; k<set_dim[2]-1; k++){
    for(int j = 0; j<oversampling; j++){
        DR[oversampling*k+j][0] = (fType)(j*focus[k+1]+(oversampling-j)*focus[k])*inv_oversampling;
        if(j<=oversampling/2) DR[oversampling*k+j][1] = k;
        else DR[oversampling*k+j][1] = k+1;
    }}
    DR[oversampling*(set_dim[2]-1)+0][0] = (fType)focus[ (set_dim[2]-1) ];
    DR[oversampling*(set_dim[2]-1)+0][1] = (fType) (set_dim[2]-1) ;
    this->DepthToRank = DR;
    build_HalfDepthToRank();
    CPING2("DR SET : ", this->DepthToRank.size());
    
    
    Mat1T matarg = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    Mat1T matrnk = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    Mat1T mat1T_ones = cv::Mat::ones(set_dim[0],set_dim[1],CV_TF);
    cv::Mat1b mat_mask = cv::Mat::zeros(set_dim[0],set_dim[1],CV_8U);
    cv::Mat1i mat_ones = cv::Mat::ones(set_dim[0],set_dim[1],CV_32S);
    cv::Mat1i mat_labels = cv::Mat::ones(set_dim[0],set_dim[1],CV_32S);

    dparam.vmat[0].copyTo(vectmat[0]);

    for(int k=0; k<set_dim[2]-1; k++){
        for(int f=1; f<=oversampling;f++){
            vectmat[k*oversampling+f] = 
		(dparam.vmat[k]*(oversampling-f) + dparam.vmat[k+1]*f)
		*inv_oversampling;
        }
    }
    vectmat[0].copyTo(dmat_score);
    dmat_label = dmat_label*0.0f;
    dmat       = dmat*0.0f;
    for(int k=0; k<set_dim[2]-1; k++){
        for(int f=1; f<=oversampling;f++){
            mat_labels = mat_ones*(k*oversampling+f);

            mat_mask = ( vectmat[k*oversampling+f]>dmat_score);
            vectmat[k*oversampling+f].copyTo(dmat_score,mat_mask);
            mat_labels.copyTo(dmat_label,mat_mask);
            
            matarg = mat1T_ones * DR[k*oversampling+f][0];
            matarg.copyTo(dmat,mat_mask);
            matrnk = mat1T_ones * DR[k*oversampling+f][1];
            matrnk.copyTo(dmat_rank, mat_mask);
        }
    }

    //for(int i=0; i<DepthToRank.size(); i++){
    //    CPING2("DR0",DepthToRank[i][0]);
    //    CPING2("DR1",DepthToRank[i][1]);
    //}
    if(!vmat_sharp_i_set)
    {
        this->vmat_sharp_i = vectmat;
        vmat_sharp_i_set=1;
    }
    return true;
}


bool DepthClass::d_poly(const tdfp_depth & dparam, Mat1T & dmat, Mat1T & dmat_rank, Mat1T & dmat_score, cv::Mat1i & dmat_label){
    // Créer un vecteur des points en lesquels évaluer chaque profondeur
    // D abscisses des profondeurs
    // X formule w/ oversampling ( k*v[i+1]+(o-k)*v[i] )/o;
    // et stocker également X², X³, etc....
    // taille degree*oversampling*set_dim[2]

    // parcourir (i,j)
    //     parcourir X(k)
    //     stocker le maximum et son argument.

    // si mat_sharp_interpol_set == false
    //     enregistrer dans sharp_interpol_vmat;

    //int dim[] = {set_dim[0],set_dim[1]};
    //cv::Mat1d tmpmat(set_dim[0],set_dim[1],CV_TF);//(2,dim,CV_32F);
    vector<Mat1T> vectmat;
    if(!vmat_sharp_i_set)
    {
        vectmat.resize(oversampling*(set_dim[2]-1));
        for(int i=0; i<oversampling*(set_dim[2]-1); i++)
        {
            vectmat[i] = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
        }
    }
    
    dmat = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    dmat = dmat + 20.0;
    dmat_score = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    dmat_rank = cv::Mat::zeros(set_dim[0],set_dim[1],CV_TF);
    dmat_label = cv::Mat::zeros(set_dim[0],set_dim[1],CV_32S);

    //tmpmat = cv::Mat1d::Mat(set_dim[0],set_dim[1],CV_TF);
    CPING2(dmat.rows,dmat.cols);
    vector<fType> X( (set_dim[2]-1)*oversampling*degree);
    vector<vector<fType> > DR( (set_dim[2]-1)*oversampling+1, vector<fType>(2));
    
    
    for(int k = 0; k<set_dim[2]-1; k++){
    for(int j = 0; j<oversampling; j++){
        DR[oversampling*k+j][0] = (fType)(j*focus[k+1]+(oversampling-j)*focus[k]) / oversampling;
        if(j<=oversampling/2) DR[oversampling*k+j][1] = k;
        else DR[oversampling*k+j][1] = k+1;
        for(int i=0; i<degree; i++){
            X[k*oversampling*degree+j*degree+i] = (fType)pow(DR[oversampling*k+j][0],i+1);
        }
    }}
    DR[oversampling*(set_dim[2]-1)+0][0] = (fType)focus[ (set_dim[2]-1) ];
    DR[oversampling*(set_dim[2]-1)+0][1] = (fType) (set_dim[2]-1) ;
    

    fType tmp;
    fType max=0;
    fType arg=0;
    fType rnk=0;
    fType linarg=0;
    for(int i=0; i<set_dim[0]; i++) for(int j=0; j<set_dim[1]; j++)
    {
        max = -1;
        arg = 20;

        for(int f=0; f<oversampling*(set_dim[2]-1);f++){

        // profondeur D[f], coeffs X[f*oversampling*(set_dim[2]-1) +k]
            tmp = dparam.vmat[0].at<fType>(i,j);
            for(int k=0;k<degree;k++){
                tmp += dparam.vmat[k+1].at<fType>(i,j)*X[f*degree + k];
            }
            // now store the data there
            if(!vmat_sharp_i_set)
            {
                vectmat[f].at<fType>(i,j) = tmp;
            }
            
            if(tmp>max){
                max = tmp;
                arg = DR[f][0];
                rnk = DR[f][1];
                dmat_label.at<int>(i,j)=f;
            }
        }
        
        dmat.at<fType>(i,j)= (fType)arg; 
        dmat_score.at<fType>(i,j)= (fType)max; // usefull for scaling
        if(max>1000)//why alert ?
        {
            //CPING2("alleeeeert",max);
        }
        dmat_rank.at<fType>(i,j) = (fType)rnk;
    }
    this->DepthToRank = DR;
    build_HalfDepthToRank();
    CPING2("DR SET : ", this->DepthToRank.size());
    //for(int i=0; i<DepthToRank.size(); i++){
    //    CPING2("DR0",DepthToRank[i][0]);
    //    CPING2("DR1",DepthToRank[i][1]);
    //}
    if(!vmat_sharp_i_set)
    {
        this->vmat_sharp_i = vectmat;
        vmat_sharp_i_set=1;
    }
    return true;
}

bool DepthClass::d_polymod(Mat1T & dmat, Mat1T & dmat_rank, Mat1T & dmat_score, cv::Mat1i & dmat_label){
///Explores vectmat around dmat_rank to determine if there exists a local maxima near. Breaks the oversampling, though. 
    for(int i=0; i<set_dim[0]; i++) for(int j=0; j<set_dim[1]; j++){
        int tmprank = dmat_rank.at<fType>(i,j);
        int tmprank2 = tmprank+2;
        int tmprank3;
        fType scoremax=sharpSetStored[tmprank].ivmat[0].at<fType>(i,j);
        while(tmprank!=tmprank2){
            tmprank3=tmprank;
            tmprank2=tmprank;
            for(int k=max(0,tmprank-1);k<min(set_dim[2],tmprank+2);k++){
                if(sharpSetStored[k].ivmat[0].at<fType>(i,j) > scoremax){
                    tmprank3=k;
                    scoremax = sharpSetStored[k].ivmat[0].at<fType>(i,j);
                }
            }
            tmprank=tmprank3;
        }
        dmat_rank.at<fType>(i,j)=(fType)tmprank;
        dmat_label.at<int>(i,j)=min(tmprank*oversampling,
				      oversampling*(set_dim[2]-1)+1 );
        dmat_score.at<fType>(i,j) = scoremax;
        dmat.at<fType>(i,j)=focus[tmprank];
    }
    return true;
}

bool DepthClass::interpolate(const vector<fType> & x, const vector<fType> & y, int n, int N, vector<fType> X, vector<fType> & z){
    vector<fType> x2(x.size()+2),y2(y.size()+2);
    if(MiscClass::optional_features[0]){
        x2[0]=x[0]-0.01f; y2[0]=y[0];
        for(int l=1; l<x2.size()-1; ++l){
            x2[l]=x[l-1];
            y2[l]=y[l-1];
        }
        x2[x2.size()-1]=x[x.size()-1]+0.01f; y2[y2.size()-1]=y[y.size()-1];
        N=N+2;
        //myLog->as("using updated interpolation \n");
    } else {
        x2=x; y2=y;
    }

    int i,j,k;
    
    //myLog->as("n == "+to_string2(n)+" \n");
    //COUT2("out : n ==",n);
    //COUT2("out : N ==",N);


    vector<fType> tmp;
// n is the degree of Polynomial, B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients
    fType B[n+1][n+2],a[n+1];
    for (i=0;i<=n;i++)
        for (j=0;j<=n;j++)
            B[i][j]=X[i+j];
    fType Y[n+1];    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)

    for (i=0;i<n+1;i++)
    {    
        Y[i]=0;
        for (j=0;j<N;j++)
        Y[i]=Y[i]+pow(x2[j],i)*y2[j];
    }
    for (i=0;i<=n;i++)
        B[i][n+1]=Y[i];
    n=n+1;                //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations

    for (i=0;i<n;i++)                    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
        for (k=i+1;k<n;k++)
            if (B[i][i]<B[k][i])
                for (j=0;j<=n;j++)
                {
                    fType temp=B[i][j];
                    B[i][j]=B[k][j];
                    B[k][j]=temp;
                }
    
    for (i=0;i<n-1;i++)            //loop to perform the gauss elimination
        for (k=i+1;k<n;k++)
            {
                fType t=B[k][i]/B[i][i];
                for (j=0;j<=n;j++)
                    B[k][j]=B[k][j]-t*B[i][j];
            }
    for (i=n-1;i>=0;i--)                //back-substitution
    {//x is an array whose values correspond to the values of x,y,z..
        a[i]=B[i][n];                //make the variable to be calculated equal to the rhs of the last equation
        for (j=0;j<n;j++)
            if (j!=i)            //then subtract all the lhs values except the coefficient of the variable whose value is being calculated
                a[i]=a[i]-B[i][j]*a[j];
        a[i]=a[i]/B[i][i];            //now finally divide the rhs by the coefficient of the variable to be calculated
    }
    for (i=0;i<n;i++){
        tmp.push_back(a[i]);
    } 

    z = tmp;
    return 0;
    return true;
}

bool DepthClass::s_poly_ij(const vector<cv::Point> & vP, const tdf_imgset & sharpSet, const tdfp_depth & dparam,const Mat1T & dmat, const string & folder){
    // calls gnuplot in order to sketch the evolution of sharpness and sharpness interpolated on one pixel.
    
    // build a vector holding sharpness vs depth (actually not that will just be the call)
    // build a vector holding polynomial interpolation X D(oversampled)
    //     peut on se contenter du vecteur DepthToRank OUI devrait être renommer depthlist
    FILE *gnuplot = popen("gnuplot", "w");
    for(int i=0; i<vP.size(); i++)
    {
        fprintf(gnuplot, "set term 'pngcairo' \n");
        fprintf(gnuplot, "set output '%s/sharpness%ix%i.png'\n",folder.c_str(), vP[i].y,vP[i].x);
        fprintf(gnuplot, "plot '-' with lines, '-' with lines, '-' with lines\n");
        
        
        for(int k=0;k<sharpSet.size();k++){
            fprintf(gnuplot, "%g %g\n", (fType)sharpSet[k].focus,
			(fType)sharpSet[k].ivmat[0].at<fType>(vP[i].y,vP[i].x));
            } 
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");
        
        for(int k=0;k<DepthToRank.size();k++){
            fType tmp3 =0;
            tmp3 =0;
            for(int l=0;l<dparam.degree+1;l++){
                tmp3+=dparam.vmat[l].at<fType>(vP[i].y,vP[i].x)*pow(DepthToRank[k][0],l);
                }
            fprintf(gnuplot, "%g %g\n", DepthToRank[k][0],tmp3);
            }
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");

        // chosen rank
        for(int k=0;k<10;k+=3){
            fType tmp3 = dmat.at<fType>(vP[i].y,vP[i].x);
            fprintf(gnuplot, "%g %g\n", tmp3, (fType) k);
            }
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");

    

        fprintf(gnuplot,"unset output \n");
        // exit gnuplot

    }  
    fprintf(gnuplot,"exit \n"); 
    pclose(gnuplot);
    
    return true;
}
bool DepthClass::s_argmax_ij(const vector<cv::Point> & vP, const tdf_imgset & sharpSet, const tdfp_depth & dparam,const Mat1T & dmat, const string & folder){
    // calls gnuplot in order to sketch the evolution of sharpness and sharpness interpolated on one pixel.
    
    // build a vector holding sharpness vs depth (actually not that will just be the call)
    // build a vector holding polynomial interpolation X D(oversampled)
    //     peut on se contenter du vecteur DepthToRank OUI devrait être renommer depthlist
    FILE *gnuplot = popen("gnuplot", "w");
    for(int i=0; i<vP.size(); i++)
    {
        fprintf(gnuplot, "set term 'pngcairo' \n");
        fprintf(gnuplot, "set output '%s/sharpness%ix%i.png'\n",folder.c_str(), vP[i].y,vP[i].x);
        fprintf(gnuplot, "plot '-' with lines, '-' with lines, '-' with lines\n");
        
        
        for(int k=0;k<sharpSet.size();k++){
            fprintf(gnuplot, "%g %g\n", (fType)sharpSet[k].focus,
			(fType)sharpSet[k].ivmat[0].at<fType>(vP[i].y,vP[i].x));
            } 
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");
        
        for(int k=0;k<DepthToRank.size();k++){
            fType tmp3 =0;
            tmp3 = dparam.vmat[k].at<fType>(vP[i].y,vP[i].x);
            fprintf(gnuplot, "%g %g\n", DepthToRank[k][0],tmp3);
            }
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");

        // chosen rank
        for(int k=0;k<10;k+=3){
            fType tmp3 = dmat.at<fType>(vP[i].y,vP[i].x);
            fprintf(gnuplot, "%g %g\n", tmp3, (fType) k);
            }
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");

    

        fprintf(gnuplot,"unset output \n");
        // exit gnuplot

    }  
    fprintf(gnuplot,"exit \n"); 
    pclose(gnuplot);
    
    return true;
}




/////////////////////////////////////////////////////////////////
// Visualisation features
/////////////////////////////////////////////////////////////////
bool DepthClass::show_interpol_generic(const vector<cv::Point> & vP, const tdf_imgset & sharpSet, const tdfp_depth & dparam,const Mat1T & dmat, const string & folder){
    // calls gnuplot in order to sketch the evolution of sharpness and sharpness interpolated on one pixel.
    
    // build a vector holding sharpness vs depth (actually not that will just be the call)
    // build a vector holding polynomial interpolation X D(oversampled)
    //     peut on se contenter du vecteur DepthToRank OUI devrait être renommer depthlist
    FILE *gnuplot = popen("gnuplot", "w");
    for(int i=0; i<vP.size(); i++)
    {
        fprintf(gnuplot, "set term 'pngcairo' \n");
        fprintf(gnuplot, "set output '%s/sharpness%ix%i.png'\n",folder.c_str(), vP[i].y,vP[i].x);
        fprintf(gnuplot, "plot '-' with lines, '-' with lines, '-' with lines\n");
        
        
        for(int k=0;k<sharpSet.size();k++){
            fprintf(gnuplot, "%g %g\n", (fType)sharpSet[k].focus,
			(fType)sharpSet[k].ivmat[0].at<fType>(vP[i].y,vP[i].x));
            } 
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");
        
        for(int k=0;k<DepthToRank.size();k++){
            fType tmp3 =0;
            tmp3 = vmat_sharp_i[k].at<fType>(vP[i].y,vP[i].x);
            fprintf(gnuplot, "%g %g\n", DepthToRank[k][0],tmp3);
            }
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");

        // chosen rank
        for(int k=0;k<10;k+=3){
            fType tmp3 = dmat.at<fType>(vP[i].y,vP[i].x);
            fprintf(gnuplot, "%g %g\n", tmp3, (fType) k);
            }
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");

    

        fprintf(gnuplot,"unset output \n");
        // exit gnuplot

    }  
    fprintf(gnuplot,"exit \n"); 
    pclose(gnuplot);
    
    return true;
}




