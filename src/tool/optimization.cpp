/****************************
*	ProjSFF
*	optimization.cpp
*	cstoribal
*	03-05-17
****************************/

#include "optimization.h"


OptiClass::OptiClass(){ 
    this->set = false;
    this->nbs_set = false;
    this->connexity = 4; //TODO
    //set all arrows to NULL pointers
    this->nbs_nb   = NULL;
    this->nbs_nbk  = NULL;
    this->nbs_n1D  = NULL;
    this->nbs_nk1D = NULL;
    this->nbs_N    = NULL;
    this->nbs_w1D  = NULL;
    this->nbs_wk1D = NULL;
    this->nbs_W    = NULL;
    
    //this->sorted_label_img = NULL;
    //this->adapt_index1D = NULL;
    //this->adapt_index = NULL;
}
OptiClass::~OptiClass(){
    if(this->nbs_set)
    {
        delete [] nbs_nb;
        delete [] nbs_n1D;
        delete [] nbs_nk1D;
        delete [] nbs_w1D;
        delete [] nbs_wk1D;
        
        delete [] nbs_nbk;
        delete [] nbs_W;
        delete [] nbs_N;
        this->nbs_set=false;
        //delete [] sorted_label_img;
        //delete [] adapt_index1D;
        //delete [] adapt_index;
        //delete adapt_Iterator;
    }
}

bool OptiClass::setlogs(IOWizard* _ioW, MyLog* mylog){
    this->ioW   = _ioW;
    this->myLog = mylog;
    return true;
}

bool OptiClass::set_optiplan(OptiPlan* _p_OptiPlanner){
    p_OptiPlanner= _p_OptiPlanner;
    return true;
}

bool OptiClass::set_blindestimation(const Mat1T & _blindmat){
    blindmat=_blindmat; //warning !! not to be changed !
    return true;
}

bool OptiClass::set_param(tdfp_opti popti, const Mat1T & _gt_dmat){

// TODO be sure to check that we are working with the good number of pixels, rows, cols...
    this->energyClass = popti.energyclass;
    this->depthClass  = popti.depthClass;
    this->evalClass  = popti.evalClass;
    this->name_opti = popti.type;
    this->nb_pixels = popti.nb_pixels;
    this->nb_labels = popti.nb_labels;
    this->width     = popti.width;
    this->height    = popti.height;
    this->labels    = popti.labels;  //ce sont des focus
    this->connexity = popti.connexity;
//TODO should use a different word for labels & focus
    _gt_dmat.copyTo(gt_dmat);

    build_rank4xy();
    ioW->mksubdir("optimized");
    newfolders=true;
    new_EdataMatrix=true;
    
    
    this->set = true;
}

bool OptiClass::do_all_optimizations(void){
    size_t nb_storedplans;
    p_OptiPlanner->get_nb_storedmethods(nb_storedplans);
    v_selected_method.resize(nb_storedplans);
    vv_rmse.resize(nb_storedplans);
    vv_psnr.resize(nb_storedplans);
    v_types.resize(nb_storedplans);
    for(size_t tmp=0; tmp<nb_storedplans; tmp++){
        v_selected_method[tmp]=tmp;
    }
    
    for(size_t tmp=0; tmp<nb_storedplans; tmp++){
        myLog->clear_iteration_times();
        actual_idx_method=tmp;
        reset(1,lambda,false);
        compute_opt_custom();
        ioW->img_unsetscale();
        ioW->writeImage("optimized/"+selected_typename+"/"+to_string2(lambda)+"/"+"BlindDiff.png",abs(blindmat-regularized_depthmat));
        myLog->write();
        myLog->clear_log();
    }
    newfolders=false;
    show_all_rmse();
    compute_write_cross_RMSE();
    
    return true;

}

bool OptiClass::do_optimization(void){
    CPING(this->name_opti);
    if(this->name_opti == "gco_grid"){
        set_optimization_gco_grid();
        return compute_gco_grid();
    }
    if(this->name_opti == "gco_custom"){ 
        return compute_opt_custom();
    }
    return error("Warning, gco" + this->name_opti + "not computed\n");
}

bool OptiClass::reset(fType l_d, fType l_r, bool _new_EdataMatrix){
    //copy weight matrix to reset it at initial status
    for(int k=0; k<nb_pixels*connexity; k++)// for(int c=0; c<connexity; c++)
    {
        nbs_wk1D[k] = nbs_w1D[k];
        nbs_nk1D[k] = nbs_n1D[k];
    }
    for(int k=0; k<nb_pixels; k++){
        nbs_nbk[k]   = nbs_nb[k];
    }
    if(l_d>0)lambda = l_r/l_d;

    regularized_labelmat = Mat::zeros(height,width,CV_32S);
    regularized_depthmat = Mat::zeros(height,width,CV_TF);
    if(_new_EdataMatrix){new_EdataMatrix=true;}
    //this->maxiteration=maxiter;
    return true;
}


bool OptiClass::writebackmatrix(Mat1T & do_mat){
// according to optimization type, stores the regularised matrix from the vectors to the sff class
    do_mat = Mat::zeros(height,width,CV_TF);
    convert_vec2mat(data_out,do_mat);
    return true;
}


///////////////////////////
//		privates //
///////////////////////////


/////////////////////////// Common


bool OptiClass::build_rank4xy(void){
    // knowing width & height, we can build the mat1i and vector<point> for conversions
    this->getrank = Mat::zeros(height,width,CV_32S);
    this->getxy   = vector<Point>(width*height);
    int k=0;
    for(int i=0;i<height;i++){
    for(int j=0;j<width;j++){
        getxy[k]=Point(j,i); 
        getrank.at<int>(i,j) = k;
        k++;        
    }}
    assert(k==this->nb_pixels);
    return true;
}

bool OptiClass::convert_mat2labvec(const vector<Mat1E> & vmat, vector<eType> & vect){
    // knowing height*width*labels we can construct the vector from the matrix vector.
    int nb_lab = vmat.size();
    vect.resize(height*width*nb_lab);
    for(int i=0; i<height;   i++){
    for(int j=0; j<width;    j++){
    for(int l=0; l<nb_lab; l++){
        vect[i*width*nb_lab+j*nb_lab+l]=vmat[l].at<eType>(i,j);        
    }}}

    return true;
}


bool OptiClass::convert_vec2mat(const vector<int> & vect, Mat1T & vmat){
    // if we know height*width (check it's nb_pixels) we can get the matrix from the vector
    // in addition, we convert it immediately to the old label format. so that's done.
    for(int i=0; i<height; i++){
    for(int j=0; j<width; j++){
        vmat.at<fType>(i,j)=labels[vect[i*width+j]];
    }}
    return true;
}

bool OptiClass::convert_vec2mat(const vector<int> & vect, Mat1i & vmat){
    // no typecast
    for(int i=0; i<height; i++){
    for(int j=0; j<width; j++){
        vmat.at<int>(i,j)=vect[i*width+j];
    }}
    return true;
}



bool OptiClass::set_allneighbors(void){
    // needed if optimisation type general graph
    if(nbs_set == true){
        myLog->a("Warning, re-setting neighbors without freeing memory");
        delete [] nbs_nb;
        delete [] nbs_n1D;
        delete [] nbs_nk1D;
        delete [] nbs_w1D;
        delete [] nbs_wk1D;

        delete [] nbs_nbk;
        delete [] nbs_W;
        delete [] nbs_N;
        this->nbs_set=false;
    }
    // sets up neighbornumarray, neighborarray, weightarray !
    // aka nbs_nb, nbs_n, nbs_w;
    Mat1E borders = Mat::zeros(height,width,CV_TE)+1;
    this->neighbor_mat = Mat::zeros(5,5,CV_TE);
    if(connexity==4 ||
       connexity==8)
    {
        neighbor_mat.at<eType>(2,1) = 0.25f;
        neighbor_mat.at<eType>(1,2) = 0.25f;
        neighbor_mat.at<eType>(3,2) = 0.25f;
        neighbor_mat.at<eType>(2,3) = 0.25f;
    }
    if(connexity==8)
    {
        neighbor_mat.at<eType>(1,1) = (eType)0.25f/(eType)sqrt(2.0);
        neighbor_mat.at<eType>(3,1) = (eType)0.25f/(eType)sqrt(2.0);
        neighbor_mat.at<eType>(1,3) = (eType)0.25f/(eType)sqrt(2.0);
        neighbor_mat.at<eType>(3,3) = (eType)0.25f/(eType)sqrt(2.0);
    }

    filter2D(borders, borders, -1, this->neighbor_mat, Point(2,2), 0, BORDER_CONSTANT );
    
    nbs_nb  = new int[nb_pixels];
    nbs_n1D = new int[nb_pixels*connexity];
    nbs_nk1D = new int[nb_pixels*connexity];
    nbs_w1D = new eType[nb_pixels*connexity];
    nbs_wk1D = new eType[nb_pixels*connexity];
    
    nbs_nbk = new int[nb_pixels];             // Nombre de voisins
    nbs_N   = new int*[nb_pixels];            // Noms des voisins
    nbs_W   = new eType*[nb_pixels];          // Poids associé au voisin

    for(int k=0; k<nb_pixels; k++){
        nbs_N[k] = &nbs_nk1D[k*connexity];
        nbs_W[k] = &nbs_wk1D[k*connexity];
    }


    int nbs_tmp_nb;

    for(int k=0; k<nb_pixels; k++)
    {
        nbs_tmp_nb = 0; // Number of neighbors
        for(int i=-2; i<3; i++) for(int j=-2; j<3; j++)
        {
            if( neighbor_mat.at<eType>(2+i,2+j) != 0 )
            if( getxy[k].y+i >= 0 && getxy[k].y+i < height )
            if( getxy[k].x+j >= 0 && getxy[k].x+j < width  )
            {
                nbs_N[k][nbs_tmp_nb]=getrank.at<int>
				(getxy[k].y+i,getxy[k].x+j);
                nbs_W[k][nbs_tmp_nb]=neighbor_mat.at<eType>(2+i,2+j);
                nbs_tmp_nb++;
            }
        }
        nbs_nbk[k]=nbs_tmp_nb;
        nbs_nb[k] = nbs_nbk[k];    //storing
    }
    
    for(int k=0; k<nb_pixels*connexity; k++) //storing
    {
        nbs_n1D[k]=nbs_nk1D[k];
        nbs_w1D[k]=nbs_wk1D[k];
    }
    
    nbs_set = true;
    return true;
}

/////////////// graph cuts
bool OptiClass::set_optimization_gco_grid(void){
// also stores the vector of matrix (depth*mat) to the data_in 
// which will call some functions from param to gco to initialise the gco
    // get vector of matrix;
    Mat1E tmp_mat = Mat::zeros(height,width,CV_TE);
    vector<Mat1E> tmp_data_energy(nb_labels);
    energyClass->getDataEnergy_3DMatrix(this->labels,tmp_data_energy);

    convert_mat2labvec(tmp_data_energy,data_in);
    tmp_mat = Mat::zeros(nb_labels,nb_labels,CV_TE);
    if(!energyClass->getCrossLabelMatrix(labels,tmp_mat)) return false;
    smoothvect.resize(nb_labels*nb_labels);
    for(int i=0; i<nb_labels; i++){
    for(int j=0; j<nb_labels; j++){
        smoothvect[i*nb_labels+j]=tmp_mat.at<eType>(i,j);
    }}
    
        myLog->time_r(6);
    try{
        gco = new GCoptimizationGridGraph(width,height,nb_labels);
        gco->setDataCost(&data_in[0]);
        gco->setSmoothCost(&smoothvect[0]);
        if(is_same<eType,int>::value) 
            printf("\nBefore optgrdnrj is %lli",gco->compute_energy());
        else
            printf("\nBefore optgrdnrj energy is %f",gco->compute_energy());
	}
    catch (GCException e){
		e.Report();
	}
    return true;
}



bool OptiClass::compute_gco_grid(void){
    gco->expansion(10);// For swap use gc->swap(num_iterations);
    data_out.resize(nb_pixels);
    for(int i=0; i<nb_pixels; i++){data_out[i] = gco->whatLabel(i);}
    if(is_same<eType,int>::value)
        printf("\nAfter optimization energy is %lli \n",gco->compute_energy());
    else
        printf("\nAfter optimization energy is %f \n",gco->compute_energy());
    return true;
    
}



bool OptiClass::compute_opt_custom(void){
    /////
    myLog->time_r(6);
    Mat1i M = Mat::zeros(height,width,CV_32S); //Global path output
    Mat1i N = Mat::zeros(height,width,CV_32S); //step   path output
    Mat1E Dp = Mat::zeros(height,width,CV_TE); //Datacost Offset
    Mat1E Dn = Mat::zeros(height,width,CV_TE); //Datacost Offset

    vector<Mat1E> Edata_slice(2);
    Edata_slice[0] = Mat::zeros(height,width,CV_TE);    //Datacost custom slice -
    Edata_slice[1] = Mat::zeros(height,width,CV_TE);    //Datacost custom slice +


    Mat1E lmat = Mat::zeros(nb_labels,nb_labels,CV_TE); // label-to-label matrix.
    if(!energyClass->getCrossLabelMatrix(labels,lmat))  // Smoothcosts
        return error("getCrossLabelMatrix");
    smoothvect.resize(2*2);
    for(int i=0; i<2; i++) for(int j=0; j<2; j++)
        smoothvect[i*2+j]=lmat.at<eType>(0+i,0+j);      // we only take 4 values 
 
    //vector<Mat1E> Edata(nb_labels);               // Datacosts
    if(new_EdataMatrix){
        Edata.resize(nb_labels);
        energyClass->getDataEnergy_3DMatrix(this->labels,Edata);
        new_EdataMatrix=false;
        ioW->img_unsetscale();
        for(int i=0; i<nb_labels; i++){
            //ioW->writeImage("edata"+to_string(i)+".png",Edata[i]);
        }
    }

    vector<int> v_nbs_reup_N(connexity);
    vector<eType> v_nbs_reup_W(connexity);
    size_t max_iterations=1;
    vector<vector<size_t> > vv_thresh(0);
    vector<vector<size_t> > vv_centro(0);
    p_OptiPlanner->get_thresholds_n_centroids(v_selected_method[actual_idx_method],
		 selected_typename, max_iterations, vv_thresh, vv_centro);
    v_types[actual_idx_method]=selected_typename;
    
    myLog->time_r(7);
    if(newfolders)ioW->mksubdir("optimized/"+selected_typename);
    ioW->mksubdir("optimized/"+selected_typename+"/"+to_string2(lambda) );

    std::string imagename;
    std::string foldername="optimized/"+selected_typename+"/"+to_string2(lambda)+"/";
    // Starting optimization  
    int n=1;
    COUT2(selected_typename," passe N° ");
    for(n=1; (n<max_iterations) ; n++)
    {    
        myLog->time_i();
        M = 2*M;
        
        std::cout << n << "..." << endl;
        for(int i=0; i<height; i++) for(int j=0; j<width; j++)
        {
            if(vv_thresh[n][M.at<int>(i,j)+1]<1)return error("indice vv_thresh 0");
            Edata_slice[0].at<eType>(i,j) = Edata[vv_thresh[n]
		[M.at<int>(i,j)+1]-1 ].at<eType>(i,j)+Dn.at<eType>(i,j);
            Edata_slice[1].at<eType>(i,j) = Edata[vv_thresh[n]
		[M.at<int>(i,j)+1]   ].at<eType>(i,j)+Dp.at<eType>(i,j);
        }
        this->convert_mat2labvec(Edata_slice,data_in);
        
        try{
            
            gco = new GCoptimizationGeneralGraph(this->nb_pixels,2);
            ((GCoptimizationGeneralGraph*)gco)
			->setAllNeighbors(nbs_nbk,nbs_N,nbs_W);
            gco->setDataCost(&data_in[0]);
            gco->setSmoothCost(&smoothvect[0]);
            gco->expansion(2);
            data_out.resize(nb_pixels);
            for(int i=0; i<nb_pixels; i++){
                data_out[i] = gco->whatLabel(i);
            }
            delete gco;
            
	}
        catch (GCException e){
	    e.Report();
            delete gco;
            return error("GCException failure in "+selected_typename+"!\n");
	}
        //MAJ VOISINAGES
        for(int i=0; i<nb_pixels; i++) for(int j=0; j<nbs_nbk[i]; j++)
        {
            if(data_out[i]!=data_out[ nbs_N[i][j]])
            {
                int l=0;
                for(int k=0; k<nbs_nbk[i]; k++)
                {
                    if(data_out[i]==data_out[ nbs_N[i][k]])
                    {
                        v_nbs_reup_N[l]=nbs_N[i][k];
                        v_nbs_reup_W[l]=nbs_W[i][k];
                        l++;
                    }
                    else{
                        if(data_out[i]>data_out[ nbs_N[i][k] ])
                            Dp.at<eType>(getxy[i]) += smoothvect[1];
                        else
                            Dn.at<eType>(getxy[i]) += smoothvect[1];
                    }
                }
                nbs_nbk[i]=l;
                for(int k=0; k<l; k++){
                    nbs_N[i][k]=v_nbs_reup_N[k];
                    nbs_W[i][k]=v_nbs_reup_W[k];
                }
            }

        }
        // OUTPUT - mise en forme
        convert_vec2mat(data_out,N); //version labels
        M = M+N;
        for(int i=0; i<height; i++) for(int j=0; j<width; j++)
        {
            size_t tmp = vv_centro[n][M.at<int>(i,j)];
            data_out[i*width+j] = tmp;
            regularized_labelmat.at<int>(i,j) = tmp;
            regularized_depthmat.at<fType>(i,j)=labels[tmp];
        }

        myLog->time_r(7+n);
        // OUTPUT - Plotting 
        imagename = to_string2(n)+"it-"+
		selected_typename+".png";
        ioW->img_setscale(1);
        ioW->writeImage(foldername+"2D-"+imagename,regularized_depthmat);
        ioW->write3DImage(foldername+"3D-"+imagename,regularized_depthmat);
        ioW->img_setscale(4);
        ioW->writeImage(foldername+"2Ddiff-"+imagename,abs(10*(regularized_depthmat-gt_dmat)) );
        ioW->img_setscale(2);
        ioW->write3DImage(foldername+"3Ddiff-"+imagename,10*(regularized_depthmat-gt_dmat) );
        //
        fType rmse,psnr;
        myLog->as("**** At iteration "+to_string2(n)+" ****\n");
        evalClass->compute_RMSE_label(regularized_depthmat,rmse);
        evalClass->compute_PSNR(regularized_depthmat,psnr);
        
        if(n==1){myLog->set_state(selected_typename,lambda);}
        myLog->set_eval_at(rmse,psnr,n);
        vv_rmse[actual_idx_method].resize(n);
        vv_psnr[actual_idx_method].resize(n);
        vv_rmse[actual_idx_method][n-1] = rmse;
        vv_psnr[actual_idx_method][n-1] = psnr;
    }
    CPING(" over");
    
    return true;
}


/////////////////////
// VISUALISATION   //
/////////////////////

bool OptiClass::show_all_rmse(void){
        FILE* gnuplot = popen("gnuplot","w");
    fprintf(gnuplot, "set term 'pngcairo' size 950,600 enhanced font 'Verdana,10'\n");
    ioW->set_gnuplot_output(gnuplot,"optimized/rmse-lambda-"+to_string2(lambda)+".png");

    std::string initplot_msg = "plot ";
    size_t nb_iter;
    initplot_msg += "'-' title '"+v_types[0]+"' with lines";
    for(int k=1;k<v_selected_method.size(); k++){
        //_tmptype="-";
        initplot_msg += ", '-' title '"+v_types[k]+"' with lines";
    }
    initplot_msg += "\n";
    fprintf(gnuplot, "%s",initplot_msg.c_str());


    for(int k=0;k<v_selected_method.size();k++){
        gnuplot_vect(gnuplot,vv_rmse[k]);
    }
    
    fprintf(gnuplot,"unset output \n");
    fprintf(gnuplot,"exit \n"); 
    pclose(gnuplot);
    return true;
}

bool OptiClass::compute_write_cross_RMSE(void){
    size_t nb_computed_plans=vv_rmse.size();
    if(nb_computed_plans<2) return error("not enough computed optimizations!");
    vector<vector<vector<std::string> > > vvv_deltaRMSEmatrix;
    vector<vector<std::string> >    vv_types;
    std::string typeA;
    std::string typeB;
    vector<double> v_rmseA, v_rmseB;
    vvv_deltaRMSEmatrix.resize(nb_computed_plans);
    vv_types.resize(nb_computed_plans);
    size_t maxiter;
    p_OptiPlanner->get_upperbound_iterations(maxiter);
    size_t miniter;
    for(int m=0; m<nb_computed_plans; m++){
        vvv_deltaRMSEmatrix[m].resize(m);
        vv_types[m].resize(m);
        typeA   = v_types[m];
        v_rmseA = vv_rmse[m];
        for(int n=0; n<m; n++){
            vvv_deltaRMSEmatrix[m][n].resize(maxiter);
            typeB   = v_types[n];
            v_rmseB = vv_rmse[n];
            vv_types[m][n]=typeA+"-"+typeB;
            miniter=min(v_rmseA.size(),v_rmseB.size());
            for(int i=0; i<min(miniter,maxiter); i++){
                vvv_deltaRMSEmatrix[m][n][i]=to_string2(v_rmseA[i]-v_rmseB[i]);
                }
            for(int i=min(miniter,maxiter); 
			i<maxiter; ++i){
                vvv_deltaRMSEmatrix[m][n][i]="";
                }
            }
        }
    if(!myLog->write_deltaRMSE("./deltaRMSE-post",vvv_deltaRMSEmatrix,vv_types)){
        return error("failure in XRMSE, size m<=1 ???");
    }
    
}

bool OptiClass::gnuplot_vect(FILE* gnuplot, vector<fType> vect){
    for(int i=0; i<vect.size(); i++){
        fprintf(gnuplot, "%i %g\n", i+1, vect[i]);
    }
    fflush(gnuplot);
    fprintf(gnuplot, "e\n");
    return true;
}

///
// ERRORS
///

bool OptiClass::error(std::string text){
    myLog->av("failure in opticlass!\n"+text+"\n");
    return false;
}

