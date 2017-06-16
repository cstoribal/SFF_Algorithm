/****************************
*	ProjSFF
*	sffClass.cpp
*	cstoribal
*	10-04-17
****************************/

/*
Holds datas
*/



#include "sffClass.h"

using namespace std;
using namespace cv;

MySFF::MySFF(){
    A_tst = Point(5,5);
}
MySFF::~MySFF(){ }



bool MySFF::loadProblem(int argc, char** argv){
    ioWizard.parseArgs(argc, argv);
    if(!ioWizard.checkArgs())
    {
        cout<<"Arguments checking failed"<<endl;
        ioWizard.displayHelp();
        return false;
    }


    ioWizard.setArgs(input_prts); //andget
    // Parsing done
    ioWizard.buildImageSet(imageSet);
    ioWizard.loadGroundTruth(gt_dmat,"");
    ioWizard.autosetImsetParameters(imageSet);
    ioWizard.mkdir(input_prts.outputfolder);
    ioWizard.set_auto_directory(input_prts.outputfolder);
    ioWizard.storeParameters();
    dim1 = imageSet[0].ivmat[0].rows;
    dim2 = imageSet[0].ivmat[0].cols;
    


    
    return true;
}

bool MySFF::preTreat(void){
    //for(int i=0; i<imageSet.size(); i++){
    //for(int j=0; j<imageSet[0].dim; j++){
        //resize(imageSet[i].ivmat[j], imageSet[i].ivmat[j], Size(), 0.2, 0.2);
        //GaussianBlur(imageSet[i].ivmat[j],imageSet[i].ivmat[j], Size(3,3),0,0);
    //}
    //}
    dim1 = imageSet[0].ivmat[0].rows;
    dim2 = imageSet[0].ivmat[0].cols;
    cout << "- d1 " << dim1 << " d2 " << dim2 << endl;

    // Gaussian blur


    return true;
}


bool MySFF::doSharpness(void){
    sharpOP.optypeSelector(input_prts.sharp);
    sharpOP.compute(this->imageSet,this->sharpSet);
    cout << "dsharp1 " << sharpSet[0].ivmat[0].cols << endl;
    //done
    
    
}

bool MySFF::doDepth(void){
    depthEst.set_param(input_prts.depth);
    // test 
    //testClass.fillSharpPoly(sharpSet);
    depthEst.buildEstimation(sharpSet, depth_parameters);
    depthEst.buildDepthmat(depth_parameters, dmat, dmat_rank, dmat_score);
    this->nb_labels = depthEst.getNbLabels();
    assert(this->nb_labels == (this->sharpSet.size()-1)*2); // Oversampling
    dmat.copyTo(rmat);
    return true;
}

bool MySFF::showInterpolationAt(void){
    //this->A_tst = ioWizard.clicked; // updates last clicked point and shows interpolation at that point 
    //TODO version vectorpoint
    //cout << "point: " << A_tst.y << " " << A_tst.x << endl;
    //depthEst.showInterpolationAt(A_tst.y,A_tst.x,sharpSet,depth_parameters, input_prts.outputfolder);
    depthEst.showInterpolationAt(ioWizard.clicked,sharpSet,depth_parameters,dmat, input_prts.outputfolder);
    ioWizard.clicked = vector<Point>();
    return true;
}

bool MySFF::clickInterpolation(Mat image, int timer){ //TODO interpolation non fonctionnelle
    ioWizard.clickImage("scale", image, timer, &forwarder, this );
    return true;
}

bool MySFF::showInterpolation(Point A){
    vector<Point> v = vector<Point>();
    v.push_back(A);
    depthEst.showInterpolationAt(v,sharpSet,depth_parameters,dmat, input_prts.outputfolder);
}

bool MySFF::setMultifocus(void){
    image_MF = Mat::zeros(Size(dim1,dim2), CV_TFC3);
    vector<Mat1T> imatBGR(imageSet[0].dim);
    //Mat1T imatBGR;
    for(int d=0; d<imageSet[0].dim; d++){
        imatBGR[d] = Mat::zeros(dim1, dim2, CV_TF);
        cout << "Multifocus " << d << " in progress..." << endl;
        for(int i=0; i<dim1; i++){
        for(int j=0; j<dim2; j++){
            // imatBGR[d].at<fType>(i,j) = imageSet[depthEst.getRankFromDepth(dmat.at<fType>(i,j))].ivmat[d].at<fType>(i,j);
            imatBGR[d].at<fType>(i,j) = imageSet[(int)dmat_rank.at<fType>(i,j)].ivmat[d].at<fType>(i,j);
        }}
    }

    // set a cross on the pointed pixel
    
    fType imin, imax;
    minMaxLoc(imatBGR[2], &imin, &imax);
    for(int i=-2; i<3; i++){
    for(int j=-2; j<3; j++){
    if(i*j==0){
        imatBGR[2].at<fType>(A_tst.y+i,A_tst.x+j)=imax;
        imatBGR[1].at<fType>(A_tst.y+i,A_tst.x+j)=(fType)0;
        imatBGR[0].at<fType>(A_tst.y+i,A_tst.x+j)=(fType)0;
    }}}

    merge(imatBGR,this->image_MF);
    
    return true;
}

bool MySFF::setMultifocusRmat(void){
    image_MF = Mat::zeros(Size(dim1,dim2), CV_TFC3);
    vector<Mat1T> imatBGR(imageSet[0].dim);
    //Mat1T imatBGR;
    for(int d=0; d<imageSet[0].dim; d++){
        imatBGR[d] = Mat::zeros(dim1, dim2, CV_TF);
        cout << "Multifocus " << d << " in progress..." << endl;
        for(int i=0; i<dim1; i++){
        for(int j=0; j<dim2; j++){
            imatBGR[d].at<fType>(i,j) = imageSet[depthEst.getRankFromDepth(rmat.at<fType>(i,j))].ivmat[d].at<fType>(i,j);
            //imatBGR[d].at<fType>(i,j) = imageSet[(int)dmat_rank.at<fType>(i,j)].ivmat[d].at<fType>(i,j);
        }}
    }

    // set a cross on the pointed pixel
    
    fType imin, imax;
    minMaxLoc(imatBGR[2], &imin, &imax);
    for(int i=-2; i<3; i++){
    for(int j=-2; j<3; j++){
    if(i*j==0){
        imatBGR[2].at<fType>(A_tst.y+i,A_tst.x+j)=imax;
        imatBGR[1].at<fType>(A_tst.y+i,A_tst.x+j)=(fType)0;
        imatBGR[0].at<fType>(A_tst.y+i,A_tst.x+j)=(fType)0;
    }}}

    merge(imatBGR,this->image_MF);
    
    return true;
}





bool MySFF::testEnergy(void){

    energyClass.set_parameters(input_prts.nrj_d, input_prts.nrj_r,sharpSet, dmat, this->nb_labels, 1);
    Mat1T rmattmp;
    dmat.copyTo(rmattmp);

    energyClass.computeMatEnergy(E_BOTH,rmattmp,vector<Point>() );
    energyClass.updateEnergy(E_BOTH);
    
    ioWizard.showImage("scale",energyClass.ed_mat,1000);
    ioWizard.showImage("scale",energyClass.er_mat,1000);
    
}
    


bool MySFF::optimize(void){
    //for(fType lambda=1;lambda<8;lambda+=0.5)
    for(int i=0;i<input_prts.optirange.size();i++)
    {
        fType lambda = input_prts.optirange[i];
        COUT2("Starting optimization at lambda= ",lambda);
        energyClass.set_parameters(input_prts.nrj_d, input_prts.nrj_r, sharpSet, dmat, this->nb_labels, lambda);
        //set
        
        opti_prts.type = input_prts.opti;
        opti_prts.energyclass = & this->energyClass;
        opti_prts.nb_pixels = dim1*dim2;
        opti_prts.nb_labels = nb_labels;
        opti_prts.height     = dim1;
        opti_prts.width      = dim2;
        opti_prts.labels    = depthEst.getLabels();
    
        optiClass.set_param(opti_prts);
        optiClass.set_optimization();
        try{
            optiClass.compute_optimization();
        }
        catch(const GCException & error)
        {
            COUT2("\nGCException encounterd at lambda = ",lambda);
            printf("\n%s\n",error.message);
            COUT("\n");
            continue;
        } 
       catch(const char* message)
        {
            COUT2("\n Error encounterd at lambda = ",lambda);
            printf("\n%s\n",message);
            COUT("\n");
            continue;
        }
        optiClass.writebackmatrix(rmat);
        string tmp = "I_depth_lambda"; tmp += to_string2(lambda); tmp += ".png" ;
        ioWizard.writeImage(tmp,this->rmat);
        ioWizard.showImage("scale",this->rmat,100);
        ioWizard.show3DImage("3D_rmat_l" + to_string2(lambda) +".png",this->rmat);

        evaluate();
        
        fType tmpmean = cv::sum(this->rmat)[0]/(dim1*dim2);
        CPING2("rmat mean is ",tmpmean);

        tmpmean = cv::sum(this->gt_dmat)[0]/(dim1*dim2);
        CPING2("gtmat mean is ",tmpmean);
        
    }
    
    
}


bool MySFF::evaluate(void){
    //dmat.copyTo(gt_dmat);
    //dmat.copyTo(rmat);
    evalClass.set_parameters(gt_dmat,depthEst.getLabels());
    evalClass.compute_RMSE(rmat,rmse,q_rmse);
}


bool forwarder(void* context, Point A) {
    static_cast<MySFF*>(context)->showInterpolation(A);
}




