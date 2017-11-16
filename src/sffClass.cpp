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
MySFF::~MySFF(){
    COUT("starting to write logs");
    myLog.write();
    //CPING("this is MySFF destructor's test");
}

bool MySFF::setlogs(void){
    // redirects all child classes pointers to the log output class
    ioWizard.setlogs(&myLog);
    pretreatClass.setlogs(&myLog);
    sharpOP.setlogs(&myLog);
    depthEst.setlogs(&myLog,&ioWizard);
    energyClass.setlogs(&myLog);
    optiPlan.set_logs(&ioWizard,&myLog);
    optiClass.setlogs(&myLog);
    evalClass.setup(&myLog,&depthEst);
    energyClass.set_class(&depthEst, &ioWizard);
}

bool MySFF::loadProblem(int argc, char** argv){
    ioWizard.parseArgs(argc, argv);
    if(!ioWizard.checkArgs())
    {
        cout<<"Arguments checking failed"<<endl;
        ioWizard.displayHelp();
        return false;
    }


    ioWizard.setArgs(this->input_prts); // get args
    myLog.log_data_out->setup(this->input_prts); // set to logs.

    // Parsing done
    ioWizard.buildImageSet(imageSet);
    ioWizard.loadGroundTruth(gt_dmat,"");
    ioWizard.autosetImsetParameters(imageSet);
    ioWizard.mkdir(input_prts.outputfolder);
    ioWizard.set_auto_directory(input_prts.outputfolder);
    ioWizard.storeParameters();
    MiscClass::optional_features = ioWizard.get_optional_features();
    return true;
}

bool MySFF::preTreat(void){
    for(int i=0; i<imageSet.size(); i++){
    for(int j=0; j<imageSet[0].dim; j++){
        pretreatClass.set_param(input_prts);
        pretreatClass.compute_scale(imageSet[i].ivmat[j]);
        pretreatClass.compute_noises(imageSet[i].ivmat[j]);
        pretreatClass.compute_blur(imageSet[i].ivmat[j]);
    }
    }
    
    pretreatClass.compute_scale(gt_dmat);
    
    dim1 = imageSet[0].ivmat[0].rows;
    dim2 = imageSet[0].ivmat[0].cols;
    myLog.a("Image Dimensions : "+to_string2(dim1)+"px rows per "+to_string2(dim2)+"px cols\n");
    //TODO Gaussian blur


    return true;
}


bool MySFF::doSharpness(void){
    sharpOP.optypeSelector(input_prts.sharp);
    sharpOP.compute(this->imageSet,this->sharpSet);
    //done
    return true;
}

bool MySFF::doDepth(void){
    depthEst.set_param(input_prts.depth);
    // test 
    //testClass.fillSharpPoly(sharpSet);
    depthEst.buildEstimation(sharpSet, depth_parameters);
    depthEst.buildDepthmat(depth_parameters, dmat, dmat_rank, dmat_score);
    vector<fType> labels = depthEst.getLabels();
    this->nb_labels = labels.size();
    // set dmat next outputs
    CPING("test depth");
    assert(this->nb_labels == (this->sharpSet.size())*depthEst.getOversampling()); // Check Oversampling
    fType d_min = labels[0];
    fType d_max = labels[this->nb_labels-1];
    ioWizard.img_setscale(d_min,d_max,1);
    ioWizard.img_setscale(10.0*(d_min-d_max),10.0*(d_max-d_min),2);
    ioWizard.img_setscale(0,nb_labels,3);
    
    depthEst.getMLabelFromMDepth(dmat,dmat_label);
    depthEst.getMLabelFromMDepth(gt_dmat,gt_label_mat);
    Utils::build_int_histogram(dmat_label,nb_labels,histogram_dmat);
    Utils::build_int_histogram(gt_label_mat,nb_labels,histogram_gtmat);
    CPING2("histogram size",histogram_dmat.size());
    ioWizard.draw_histogram(histogram_gtmat,histogram_dmat);

    
    
    
    ioWizard.mksubdir("data");
    ioWizard.img_setscale(1);
    ioWizard.writeImage("data/groundtruth.png",this->gt_dmat);
    ioWizard.img_setscale(3);
    ioWizard.writeImage("data/dmat_label.png",this->dmat_label);
    ioWizard.img_unsetscale();
    Mat1T tmpmat;
    log(dmat_score+1.0f,tmpmat);
    ioWizard.writeImage("data/Iscore.png",tmpmat);
    
    debug_check_all("end doDepth");
    dmat.copyTo(rmat);
    CPING2("nblabels", nb_labels);
    CPINGIF4("End Depth, scale1 : ", d_min, d_max, " !",true);
    return true;
}



bool MySFF::prepare_optimization_plan(void){
    optiPlan.set_groundtruth(gt_label_mat);
    std::vector<std::string> typelist(5);
    typelist[0]="binary";
    typelist[1]="binary_v2";
    typelist[2]="otsu";
    //typelist[3]="otsu_v0";
    typelist[3]="median";
    typelist[4]="2means";
    COUT("starting optimization plan");
    for(int i=0; i<typelist.size(); i++){
        optiPlan.set_param(typelist[i],nb_labels,dim1*dim2,histogram_dmat,1,1);
    }
    //COUT("optiplan set");
    CPING("showallRMSE");
    optiPlan.show_all_RMSE("RMSEall");
    Mat1i gt_dmat_label;
    
    COUT("optiplan RMSE computed");
    optiPlan.show_all_thresh_plans("Threshplan_");
    optiPlan.write_all_ThreshedMatrix(dmat_label); // rmse_gt set !
    optiPlan.computeCrossRMSEperf_andLog();

    optiPlan.show_all_RMSE2("RMSEall2");
    optiPlan.addToLog();
    
    return true;
}

bool MySFF::showInterpolationAt(void){
    depthEst.showInterpolationAt(ioWizard.clicked,sharpSet,depth_parameters,dmat, input_prts.outputfolder);
    ioWizard.clicked = vector<Point>();
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
        }}
    }

    
    merge(imatBGR,this->image_MF);
    
    CPING("pong");
    return true;
}

    


bool MySFF::optimize(void){

    ioWizard.mksubdir("optimized");
    
    int timernk = 5;
    //for(fType lambda=1;lambda<8;lambda+=0.5)
    tdf_input& ip = input_prts; //alias. too complex
    if(ip.vect_lambda_d.size()<ip.vect_lambda_r.size()){
        int start =ip.vect_lambda_d.size();
        fType lambdareplicate = ip.vect_lambda_d[start-1];
        ip.vect_lambda_d.resize(ip.vect_lambda_r.size());
        for(int i=start;i<ip.vect_lambda_r.size();i++){
            ip.vect_lambda_d[i]=lambdareplicate;
        }
    }
    if(ip.vect_lambda_d.size()>ip.vect_lambda_r.size()){
        int start = ip.vect_lambda_r.size();
        fType lambdareplicate = ip.vect_lambda_r[start-1];
        ip.vect_lambda_r.resize(ip.vect_lambda_d.size());
        for(int i=start;i<ip.vect_lambda_d.size();i++){
            ip.vect_lambda_r[i]=lambdareplicate;
        }
    } 
    // fill in differences.
    opti_prts.type = input_prts.opti;
    opti_prts.energyclass = & this->energyClass;
    opti_prts.nb_pixels = dim1*dim2;
    opti_prts.nb_labels = nb_labels;
    opti_prts.height    = dim1;
    opti_prts.width     = dim2;
    opti_prts.connexity = input_prts.connexity;
    opti_prts.labels    = depthEst.getLabels();
    opti_prts.histogram = histogram_dmat;
    COUT2("nb_pixels", opti_prts.nb_pixels);
    
    optiClass.set_param(opti_prts);
    optiClass.set_allneighbors();
    myLog.a("Optimization class set");
    
    vector<int>   vectindex(3);
    vector<fType> vectrmse(3);
    vector<fType> vectlamb_d(3);
    vectrmse[0]=1000000;
    vectrmse[1]=1000000;
    vectrmse[2]=1000000;

    int recorded_states=3;
    int higherRMSE_index=0;
    int secondRMSE_index=1;
    fType denom = 1;
    fType precrmse = 0;
    
    int i=0;
    int loop=1;
    int delay_end=0;
    fType lambda_r;
    fType lambda_d;

    int startiterationmax = input_prts.maxiteration<1?(int)floor(log2(this->nb_labels-1)-input_prts.maxiteration+10):input_prts.maxiteration;

    myLog.time_r(5);
    while((i<input_prts.vect_lambda_r.size()+20) & (loop>0)) 
    {// 20 max iterations search
    for(int maxiter=startiterationmax; maxiter>0; maxiter--)
    {
        if(i<input_prts.vect_lambda_r.size())
        {
            lambda_r = input_prts.vect_lambda_r[i];
            lambda_d = input_prts.vect_lambda_d[i];
        }
        else if(i==input_prts.vect_lambda_r.size())
        {
            if(input_prts.vect_lambda_r[vectindex[higherRMSE_index]]==0){
                input_prts.vect_lambda_r[vectindex[higherRMSE_index]]=1;
                input_prts.vect_lambda_d[vectindex[higherRMSE_index]]=input_prts.vect_lambda_d[vectindex[higherRMSE_index]]*10;
            }
            if(input_prts.vect_lambda_r[vectindex[secondRMSE_index]]==0){
                input_prts.vect_lambda_r[vectindex[secondRMSE_index]]==1;
                input_prts.vect_lambda_d[vectindex[secondRMSE_index]]=input_prts.vect_lambda_d[vectindex[higherRMSE_index]]*10;
            }
            denom = input_prts.vect_lambda_r[vectindex[higherRMSE_index]]*
                    input_prts.vect_lambda_r[vectindex[secondRMSE_index]];
            vectlamb_d[higherRMSE_index]=
                    input_prts.vect_lambda_d[vectindex[higherRMSE_index]]*
                    input_prts.vect_lambda_r[vectindex[secondRMSE_index]];
            vectlamb_d[secondRMSE_index]=
                    input_prts.vect_lambda_d[vectindex[secondRMSE_index]]*
                    input_prts.vect_lambda_r[vectindex[higherRMSE_index]];
            vectlamb_d[2*(higherRMSE_index + secondRMSE_index)%3]=
                    vectlamb_d[higherRMSE_index]*5.0/8.0+
                    vectlamb_d[secondRMSE_index]*3.0/8.0;
            lambda_r = denom;
            lambda_d = vectlamb_d[2*(higherRMSE_index + secondRMSE_index)%3];
            loop=2;
            loop=0; // No alpha search ! 
        }
        else
        {
            vectlamb_d[2*(higherRMSE_index + secondRMSE_index)%3]=
                    vectlamb_d[higherRMSE_index]*5.0/8.0+
                    vectlamb_d[secondRMSE_index]*3.0/8.0;
            lambda_r = denom;
            if(lambda_d ==vectlamb_d[2*(higherRMSE_index + secondRMSE_index)%3]) loop=0;
            lambda_d = vectlamb_d[2*(higherRMSE_index + secondRMSE_index)%3];
        }


        string tmp = "Dep-D" + to_string2(lambda_d) + "R" + to_string2(lambda_r)+"It"+to_string2(maxiter);
        
        if(!this->launch_optimization(lambda_d,lambda_r,maxiter, rmat) ) continue;
        evaluate();
        
        if(rmse<=vectrmse[higherRMSE_index])
        {
            CPING2("fattest RMSE : ",rmse);
            higherRMSE_index = 2*(higherRMSE_index + secondRMSE_index)%3;
            vectrmse[higherRMSE_index]=rmse;
            vectindex[higherRMSE_index]=i;
            secondRMSE_index = 2*(higherRMSE_index + secondRMSE_index)%3;
        }
        else if(rmse<=vectrmse[secondRMSE_index])
        {
            CPING2("second RMSE : ",rmse);
            secondRMSE_index = 2*(higherRMSE_index + secondRMSE_index)%3;
            vectrmse[secondRMSE_index]=rmse;
            vectindex[secondRMSE_index]=i;
        }
        else if(rmse<=vectrmse[2*(higherRMSE_index + secondRMSE_index)%3])
        {
            vectrmse[2*(higherRMSE_index + secondRMSE_index)%3]=rmse;
            vectindex[2*(higherRMSE_index + secondRMSE_index)%3]=i;
            if(loop==2) loop=0; //break;
        }
        else 
        {
            if(loop==2) loop=0;
            CPING2("skipped, RMSE : ", rmse);
        }
        if(loop==2 & vectrmse[higherRMSE_index]<vectrmse[secondRMSE_index]*1.0001)
        {
            if(delay_end++>0) loop=0;
        }
        // significates that rmse don't go much higher.
        precrmse=rmse;

        ioWizard.img_setscale(1);
        ioWizard.writeImage("optimized/2D-"+tmp+ ".png",this->rmat);
        //ioWizard.showImage("scale",tmp,this->rmat,100);
        ioWizard.write3DImage("optimized/3D-"+tmp+ ".png",this->rmat);

        ioWizard.img_setscale(2);
        ioWizard.writeImage("optimized/2Ddiff-"+tmp+ ".png",10*(this->rmat-this->gt_dmat));
        ioWizard.write3DImage("optimized/3Ddiff-"+tmp+ ".png",10*(this->rmat-this->gt_dmat) );
        ioWizard.img_unsetscale();
        myLog.a(tmp+"\n");
        myLog.time_r(7);
        myLog.set_state(lambda_r,lambda_d,maxiter);
        myLog.set_eval(rmse,psnr);
        myLog.write();
        maxiter =0;
               
    }
    i++;
    }
    
    


    
    debug_check_all("end");
}




bool MySFF::launch_optimization(fType l_d, fType l_r, int maxiter, Mat1T & output){

    COUT2("Starting optimization at lambda_r = ",l_r);
    COUT2("Starting optimization at lambda_d = ",l_d);

    energyClass.set_parameters(input_prts.nrj_d, input_prts.nrj_r, sharpSet, dmat, l_d, l_r);
    //set
    optiClass.reset(maxiter);


    try{
        if(!optiClass.do_optimization() ) 
        {
            CPING("continue");
            return false;
        }
    }
    catch(const GCException & error)
    {
        COUT2("\nGCException encounterd at lambda_r = ",l_r);
        printf("\n%s\n",error.message);
        COUT("\n");
        return false;
    } 
   catch(const char* message)
    {
        COUT2("\n Error encounterd at lambda_r = ",l_r);
        printf("\n%s\n",message);
        COUT("\n");
        return false;
    }
    optiClass.writebackmatrix(output);
    
    return true;
}


bool MySFF::evaluate(void){
    evalClass.set_parameters(gt_dmat,depthEst.getLabels());
    evalClass.compute_RMSE(rmat,rmse,q_rmse);
    evalClass.compute_RMSE_label(rmat,rmse,q_rmse);
    evalClass.compute_PSNR(rmat,psnr);
}

bool MySFF::setNewProblem(void){
    // if the 1st loading ioWizard returns -M as an option (macro) then
    // store the string vector with problems adresses indexed by iPRB,
    // destroy all classes, reset variables
    // iPRB++
    // logout the iprb (for debug)
    // reset classes, relink to logs,
    // call loadProblem with a new int argc, char** argv
    // basically argc = 2, (?), argv = [-D,../../Samples/.....]
    
    
    return false;
}

bool MySFF::debug_check_all(std::string context){
    debug_MMCheck(this->gt_dmat,"gt_dmatrix "+context);
    debug_MMCheck(this->gt_label_mat,"gt_label_mat "+context);
    debug_MMCheck(this->dmat,"dmatrix "+context);
    debug_MMCheck(this->dmat_rank,"drankmatrix "+context);
    debug_MMCheck(this->dmat_score,"dscorematrix "+context);
    debug_MMCheck(this->dmat_label,"dlabelmatrix "+context);
    //debug_MMCheck(this->dmat_labelfloat,"dlabelfloatmatrix "+context);
    debug_MMCheck(this->rmat,"rmatrix "+context);
    //debug_MMCheck(this->image_MF,"MF_image "+context);
    
}


bool MySFF::debug_MMCheck(const cv::Mat & matrix, std::string name){
    double _min,_max;
    minMaxLoc(matrix,&_min,&_max);
    CPINGIF4("matrix min max : ",name, _min, _max, true);
    return true;
}
