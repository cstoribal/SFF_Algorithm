/****************************
*	testClass.cpp
*	cstoribal
*	09-04-17
****************************/

/*
Aimed at holding test functions
*/



#include "testClass.h"

using namespace std;
using namespace cv;


TestClass::TestClass() {

}

TestClass::~TestClass(){

}

bool TestClass::coutMatType(Mat matrix){

    cout << matrix << endl;
    return true;
}


bool TestClass::matrixOp1(Mat& matrix){
    matrix = 2*matrix;

    return true;
}

bool TestClass::matrixOp2(Mat& matrix){
    matrix = (matrix>125)*256;

    return true;
}

bool TestClass::matrixOp3(Mat& M){
    //M -= M;
    // NB: gestion de l'overflow OK
    M += Mat::eye(M.rows, M.cols, CV_8UC3)*255 - M;
    M -= Vec3b (65,255,255);
    // M += Mat::eye(M.rows, M.cols, CV_8UC3(255,255,255) );

    return true;
}




bool TestClass::rgbToGray(Mat3f M, Mat1f& N){

    cvtColor(M, N, CV_BGR2GRAY);
    
    
    return true;
}

bool TestClass::cmatToArraymat(Mat3f M, vector<Mat1f>& Vm){
    split(M,Vm);
    
    
    return true;
}


bool TestClass::fillSharpPoly(typedef_imgset & sharpSet)
{
    int A[] = {3,100,-25,1};
    double tmp = 0;
    for(int k=0; k<sharpSet.size(); k++){
        tmp = A[0]+A[1]*k+A[2]*k*k+A[3]*k*k*k;
        sharpSet[k].ivmat[0] = Mat(sharpSet[0].ivmat[0].rows,sharpSet[0].ivmat[0].cols, CV_64F, tmp);
    cout << k << endl;
    }
    cout << "ss1" << sharpSet[0].ivmat[0].cols << " & " << sharpSet[0].ivmat[0].rows << endl;
    
    return true;
}


bool TestClass::fillAMatrix(Mat1d & imat)
{
    //image_MF = Mat::zeros(Size(dim1,dim2), CV_64FC3);
    //vector<Mat1d> imatBGR(imageSet[0].dim);
    int dim1 = 960;
    int dim2 = 540;
    cout<< "pong1" << endl;
    for(int d=0; d<3; d++){
        imat = Mat::zeros(dim1, dim2, CV_64F);
        for(int i=0; i<dim1; i++){
        for(int j=0; j<dim2; j++){
            cout << imat.at<double>(i,j) << " , at " << i << " " << j << " for d=" << d << endl;
            imat.at<double>(i,j) = 1; 
            cout << imat.cols << "&row&" << imat.rows << endl ;
        }}
    }
    return true;
}

// c'est donc Rows - Cols. Sorry.





int TestClass::polyfit()
{
    return 0;
}

