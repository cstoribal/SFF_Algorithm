/****************************
*	pretreat.cpp
*	cstoribal
*	15-06-17
****************************/

/*
Outputs at the same time a precise log of each simulation encountered
and the results under the format of a csv
*/


#include "pretreat.h"



PreTreatment::PreTreatment(){
    this->blurwin = 0;
    this->scale = 1.0;
    this->noise_a = 0.0;
    this->noise_b = 0.0;
    this->noise_ca = 0.0;
    this->noise_cs = 0.0;
    this->myLog = new MyLog();
    srand(time(NULL));
}
PreTreatment::~PreTreatment(){}

bool PreTreatment::setlogs(MyLog* mylog){
    this->myLog = mylog;
    return true;
}

bool PreTreatment::set_param(const tdf_input & prts){
    this->blurwin  = prts.gauss;
    this->scale    = prts.scale;
    this->noise_a  = prts.noise_a;
    this->noise_b  = prts.noise_b;
    this->noise_ca = prts.noise_ca;
    this->noise_cs = prts.noise_cs; 
    return true;
}
bool PreTreatment::compute_scale(Mat_<fType> & image){
    //if(scale!=1||scale!=0){
        resize(image,image,Size(),1,1);
        // déjà fait au chargement des images.
        return true;
    //}

    return false;
}

bool PreTreatment::compute_blur(Mat_<fType> & image){
    if(blurwin!=1||blurwin!=0)
    {
        GaussianBlur(image,image, Size(blurwin,blurwin),0,0);
        return true;
    }
    
    
    
    return true;
}

bool PreTreatment::compute_noises(Mat_<fType> & image){
    this->dimy = image.rows;
    this->dimx = image.cols;
    this->dimc = image.channels();

    if(this->dimx==0||this->dimy==0){
        //myLog->a("Pretreatment noise impossible, image dimention is null");
        return false;
    }

    if(noise_a){
        //myLog->a( "Computing noise_mult "+to_string(noise_a)+"\n");
        this->compute_nmult(image);
    }

    if(noise_b){
        //myLog->a( "Computing noise_add "+to_string(noise_b)+"\n");
        this->compute_nunif(image);
    }

    if(noise_ca){
        //myLog->a( "Computing noise_gauss "+to_string(noise_ca)+ "sigma " +to_string2(noise_cs)+"\n");
        this->compute_ngauss(image);
    }
    
    //image.convertTo(image, CV_32F);
    threshold( image, image, 1.0, 1.0,THRESH_TRUNC);
    image = -image;
    threshold( image, image, 0.0, 0.0,THRESH_TRUNC);
    image = -image;
    //image.convertTo(image, CV_TF);
    
    //image = image.mul(Mat::zeros(dimy,dimx,CV_TF),(image<0));
    
    return true;
}

bool PreTreatment::compute_nmult(Mat_<fType> & image){
    fType coef1,coef2;
    coef1=1+abs(noise_a);
    coef2=1-abs(noise_a);
    Mat_<fType> matnoise(this->dimy,this->dimx, image.type());
    if(this->dimc == 1)
        randu(matnoise,Scalar(coef1),Scalar(coef2) );
    if(this->dimc == 2)
        randu(matnoise,Scalar(coef1,coef1),Scalar(coef2,coef2) );
    if(this->dimc == 3)
        randu(matnoise,Scalar(coef1,coef1,coef1),
			Scalar(coef2,coef2,coef2) );
    if(this->dimc == 4)
        randu(matnoise,Scalar(coef1,coef1,coef1,coef1),
			Scalar(coef2,coef2,coef2,coef2) );
    image = image.mul(matnoise);

    return true;
}

bool PreTreatment::compute_nunif(Mat_<fType> & image){
    fType coef1,coef2;
    coef1=+abs(noise_b);
    coef2=-abs(noise_b);

    Mat_<fType> matnoise(this->dimy,this->dimx, image.type());
    if(this->dimc == 1)
        randu(matnoise,Scalar(coef1),Scalar(coef2) );
    if(this->dimc == 2)
        randu(matnoise,Scalar(coef1,coef1),Scalar(coef2,coef2) );
    if(this->dimc == 3)
        randu(matnoise,Scalar(coef1,coef1,coef1),
			Scalar(coef2,coef2,coef2) );
    if(this->dimc == 4)
        randu(matnoise,Scalar(coef1,coef1,coef1,coef1),
			Scalar(coef2,coef2,coef2,coef2) );

    image = image + matnoise;

    return true;
}

bool PreTreatment::compute_ngauss(Mat_<fType> & image){
    fType coef1,coef2;
    coef1=0;
    coef2=this->noise_cs;

    Mat_<fType> matnoise(this->dimy,this->dimx, image.type());
    if(this->dimc == 1)
        randn(matnoise,Scalar(coef1),Scalar(coef2) );
    if(this->dimc == 2)
        randn(matnoise,Scalar(coef1,coef1),Scalar(coef2,coef2) );
    if(this->dimc == 3)
        randn(matnoise,Scalar(coef1,coef1,coef1),
			Scalar(coef2,coef2,coef2) );
    if(this->dimc == 4)
        randn(matnoise,Scalar(coef1,coef1,coef1,coef1),
			Scalar(coef2,coef2,coef2,coef2) );

    image = image + noise_ca*matnoise;
    return true;
}














