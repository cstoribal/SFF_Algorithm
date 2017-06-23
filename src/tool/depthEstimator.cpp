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




DepthClass::DepthClass(){ this->set = false; }
DepthClass::~DepthClass(){}

bool DepthClass::setlogs(MyLog* mylog){
    this->myLog = mylog;
    return true;
}

bool DepthClass::set_param(string settype){
    this->type = settype;
    this->degree = 8;
    this->oversampling = 2;
    cout<< this->type <<endl;
    this->set = true;
    return true;
}

bool DepthClass::buildEstimation(void){
    return true;
}

bool DepthClass::buildEstimation(const tdf_imgset & sharpSet, tdfp_depth & pmat){
    this->set_dim[0] = sharpSet[0].ivmat[0].rows;
    this->set_dim[1] = sharpSet[0].ivmat[0].cols;
    this->set_dim[2] = sharpSet.size();

    for(int k=0;k<set_dim[2];k++)
    {
        this->focus.push_back(sharpSet[k].focus);
    }
    if(type=="polynome")
        f_poly(sharpSet, pmat);
    return true;
}

bool DepthClass::buildDepthmat(const tdfp_depth & dparam, cv::Mat1d & dmat,cv::Mat1d & dmat_rank, cv::Mat1d & dmat_score ){

    if(type=="polynome")
        d_poly(dparam, dmat, dmat_rank, dmat_score);
    cout << "depthmap built" << endl;
    return true;
}

int DepthClass::getRankFromDepth(double input){
    int rank = -1;
    int i=0;
    while(rank == -1 and i<DepthToRank.size() ){
        if(DepthToRank[i][0]>=input)
            rank = (int)round(DepthToRank[i][1]);
        i++;
    }
    if(rank == -1)
    {
        cout << "could not find rank, value set to last image" << endl;
        rank = (int)round(DepthToRank[DepthToRank.size()-1][1]);
    }
    return rank;
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
        

bool DepthClass::showInterpolationAt(const vector<cv::Point> & vP, const tdf_imgset & SharpSet, const tdfp_depth & dparam, const cv::Mat1d & dmat, const string & folder){
    if(type=="polynome") s_poly_ij(vP, SharpSet, dparam, dmat, folder);
    
    return true ;
}





//////////////////////////
// PRIVATE //
//////////////////////////





bool DepthClass::f_poly(const tdf_imgset & sharpSet, tdfp_depth & dparam){
    // gives the parameters of each i*j polynom (coefficients) from 
    if(sharpSet.size() == 0) return false;
    int N = set_dim[2];
    int dim[] = {set_dim[0],set_dim[1],degree+1};
    //vector<Mat1d> tmpmat(degree+1);
    vector<double> x(N); 
    vector<double> y(N);
    vector<double> z(degree+1);
    dparam.vmat = vector<cv::Mat1d>(degree+1);
    for(int k=0;k<degree+1;k++){
        dparam.vmat[k] = cv::Mat::zeros(set_dim[0],set_dim[1], CV_64F);}
    for(int k=0;k<N; k++)
        {x[k]=sharpSet[k].focus;}
    

    vector<double> X(2*degree+1);
    for(int i=0;i<2*degree+1;i++)
    {
        X[i]=0;
        for(int j=0;j<N;j++)
            X[i]=X[i]+pow(x[j],i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    }
    
    
    for(int i=0; i<set_dim[0]; i++){
        for(int j=0; j<set_dim[1]; j++){

            for(int k=0; k<set_dim[2]; k++){
                y[k]=sharpSet[k].ivmat[0].at<double>(i,j);
            }
            this->interpolate(x,y,degree,set_dim[2],X,z);
            
            for(int k=0; k<degree+1; k++){
                dparam.vmat[k].at<double>(i,j) = z[k];
            }
        }
    }
    dparam.type   = this->type;
    dparam.degree = this->degree;
    
    return true;
}

bool DepthClass::d_poly(const tdfp_depth & dparam, cv::Mat1d & dmat, cv::Mat1d & dmat_rank, cv::Mat1d & dmat_score){
    // Créer un vecteur des points en lesquels évaluer chaque profondeur
    // D abscisses des profondeurs
    // X formule w/ oversampling ( k*v[i+1]+(o-k)*v[i] )/o;
    // et stocker également X², X³, etc....
    // taille degree*oversampling*set_dim[2]

    // parcourir (i,j)
    //     parcourir X(k)
    //     stocker le maximum et son argument.

    //ctrlf pdmat

    //int dim[] = {set_dim[0],set_dim[1]};
    //cv::Mat1d tmpmat(set_dim[0],set_dim[1],CV_64F);//(2,dim,CV_32F);
    dmat = cv::Mat::zeros(set_dim[0],set_dim[1],CV_64F);
    dmat_score = cv::Mat::zeros(set_dim[0],set_dim[1],CV_64F);
    dmat_rank = cv::Mat::zeros(set_dim[0],set_dim[1],CV_64F);

    //tmpmat = cv::Mat1d::Mat(set_dim[0],set_dim[1],CV_64F);
    CPING2(dmat.rows,dmat.cols);
    vector<double> X( (set_dim[2]-1)*oversampling*degree);
    vector<vector<double> > DR( (set_dim[2]-1)*oversampling, vector<double>(2));
    
    
    for(int k = 0; k<set_dim[2]-1; k++){
    for(int j = 0; j<oversampling; j++){
        DR[oversampling*k+j][0] = (double)(j*focus[k+1]+(oversampling-j)*focus[k]) / oversampling;
        if(j<oversampling/2) DR[oversampling*k+j][1] = k;
        else DR[oversampling*k+j][1] = k+1;
        for(int i=0; i<degree; i++){
            X[k*oversampling*degree+j*degree+i] = (double)pow(DR[oversampling*k+j][0],i+1);
        }
    }}

    double tmp;
    double max=0;
    double arg=0;
    double rnk=0;
    double linarg=0;
    for(int i=0; i<set_dim[0]; i++){
    for(int j=0; j<set_dim[1]; j++){
        max = -1;
	//if(i==712 and j==1493) cout << "maximum is " << max << endl;
        arg = 0;

        for(int f=0; f<oversampling*(set_dim[2]-1);f++){

        // profondeur D[f], coeffs X[f*oversampling*(set_dim[2]-1) +k]
            tmp = dparam.vmat[0].at<double>(i,j);
            for(int k=0;k<degree;k++){
                tmp += dparam.vmat[k+1].at<double>(i,j)*X[f*degree + k];
            }
            //if(i==712 and j==1493) cout << "score" << tmp << " arg " << DR[f][0] << " bytheway " << DR[f][1] << endl;
            if(tmp>max){
                max = tmp;
                arg = DR[f][0];
                rnk = DR[f][1];
		//if(i==712 and j==1493) cout << "hit !" << endl;
            }
        }
        
        dmat.at<double>(i,j)= (double)arg; //TODO uncomment
        dmat_score.at<double>(i,j)= (double)max;
        dmat_rank.at<double>(i,j) = (double)rnk;
    }
    }
    this->DepthToRank = DR;

    cout<<"ping2"<<endl;
    return true;
}



bool DepthClass::interpolate(const vector<double> & x, const vector<double> & y, int n, int N, vector<double> X, vector<double> & z){
    int i,j,k;
    vector<double> tmp;
// n is the degree of Polynomial, B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients
    double B[n+1][n+2],a[n+1];
    for (i=0;i<=n;i++)
        for (j=0;j<=n;j++)
            B[i][j]=X[i+j];
    double Y[n+1];    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)

    for (i=0;i<n+1;i++)
    {    
        Y[i]=0;
        for (j=0;j<N;j++)
        Y[i]=Y[i]+pow(x[j],i)*y[j];
    }
    for (i=0;i<=n;i++)
        B[i][n+1]=Y[i];
    n=n+1;                //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations

    for (i=0;i<n;i++)                    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
        for (k=i+1;k<n;k++)
            if (B[i][i]<B[k][i])
                for (j=0;j<=n;j++)
                {
                    double temp=B[i][j];
                    B[i][j]=B[k][j];
                    B[k][j]=temp;
                }
    
    for (i=0;i<n-1;i++)            //loop to perform the gauss elimination
        for (k=i+1;k<n;k++)
            {
                double t=B[k][i]/B[i][i];
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

bool DepthClass::s_poly_ij(const vector<cv::Point> & vP, const tdf_imgset & sharpSet, const tdfp_depth & dparam,const cv::Mat1d & dmat, const string & folder){
    // calls gnuplot in order to sketch the evolution of sharpness and sharpness interpolated on one pixel.
    
    // build a vector holding sharpness vs depth (actually not that will just be the call)
    // build a vector holding polynomial interpolation X D(oversampled)
    //     peut on se contenter du vecteur DepthToRank OUI devrait être renommer depthlist
/*    for(int k=0; k<DepthToRank.size();k++)
    {
        cout << DepthToRank[k][0] << endl;
    }
*/
    FILE *gnuplot = popen("gnuplot", "w");
    for(int i=0; i<vP.size(); i++)
    {
        fprintf(gnuplot, "set term 'pngcairo' \n");
        fprintf(gnuplot, "set output '%s/sharpness%ix%i.png'\n",folder.c_str(), vP[i].y,vP[i].x);
        fprintf(gnuplot, "plot '-' with lines, '-' with lines, '-' with lines\n");
        
        
        for(int k=0;k<sharpSet.size();k++){
            fprintf(gnuplot, "%g %g\n", (double)sharpSet[k].focus,
			(double)sharpSet[k].ivmat[0].at<double>(vP[i].y,vP[i].x));
            } 
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");
        
        for(int k=0;k<DepthToRank.size();k++){
            double tmp3 =0;
            tmp3 =0;
            for(int l=0;l<dparam.degree+1;l++){
                tmp3+=dparam.vmat[l].at<double>(vP[i].y,vP[i].x)*pow(DepthToRank[k][0],l);
                }
            fprintf(gnuplot, "%g %g\n", DepthToRank[k][0],tmp3);
            }
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");

        // chosen rank
        for(int k=0;k<10;k+=3){
            double tmp3 = dmat.at<double>(vP[i].y,vP[i].x);
            fprintf(gnuplot, "%g %g\n", tmp3, (double) k);
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





