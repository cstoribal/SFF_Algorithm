/****************************
*	ProjSFF
*	optiplan.cpp
*	cstoribal
*	13-10-17
****************************/

#include "../tool/optiplan.h"


/////////////////////////////
/////////////////////////////
/////////////////////////////

OptiStored::OptiStored():type(""),stored_label(0),nb_iterations(0),h_thresh(false),h_centroid(false),h_rmse(false),h_rmsegt(false){}
OptiStored::~OptiStored(){
    if(h_thresh){
        //delete [] a_thresh;
        //delete [] aa_thresh;
    }
    if(h_centroid){
        //delete [] a_centroid;
        //delete [] aa_centroid;
    }
    //if(h_rmse){delete [] a_rmse;}
}
bool OptiStored::copy_from_plan(std::string& _type, 
	size_t _label, size_t _nb_iter, 
	const vector<vector<size_t> > 	& _vv_thresh, 
	const vector<vector<size_t> > 	& _vv_centroid, 
	const vector<float> 		& _v_rmse){
    if(type!=""){return false;}
    if(h_thresh||h_rmse||h_centroid){return false;}
    type=_type;
    stored_label=_label;
    nb_iterations=_nb_iter;
    /// init
    //size_t arraysize = (size_t)pow(2,nb_iterations)+nb_iterations;
    //a_thresh = new size_t [arraysize];
    //aa_thresh= new size_t* [nb_iterations]; h_thresh=true;
    //a_centroid = new size_t [arraysize];
    //aa_centroid= new size_t* [nb_iterations]; h_centroid=true;
    //a_rmse   = new float [nb_iterations]; h_rmse  =true;
    /// fill
    vv_thresh = _vv_thresh;
    vv_centroid = _vv_centroid;
    v_rmse = _v_rmse;

    h_rmse= true;
    h_centroid=true;
    h_thresh=true;
    /*
    for(size_t i=0; i<nb_iterations; i++){
        a_rmse[i]    = _a_rmse[i];
        aa_thresh[i] = &a_thresh[(int)(pow(2,i)+i)-1];
        aa_centroid[i] = &a_centroid[(int)(pow(2,i)+i)-1];
        for(size_t j=0; j<=pow(2,nb_iterations);j++){
            a_thresh[(int)(pow(2,i)+i)-1+j] = 
		_a_thresh[(int)(pow(2,i)+i)-1+j];
            a_centroid[(int)(pow(2,i)+i)-1+j] = 
		_a_centroid[(int)(pow(2,i)+i)-1+j];
        }
        
    }
    */
    return true;
}


bool OptiStored::get_data_pointers(std::string& _type, size_t& _nb_iter,
	vector<vector<size_t> >	& _vv_thresh, 
	vector<vector<size_t> >	& _vv_centroid, 
	vector<float>		& _v_rmse){
    if(!h_thresh||!h_rmse||!h_centroid){return false;}
    if(type==""){return false;}
    
    _type = type; _nb_iter=nb_iterations;
    _vv_thresh = vv_thresh; 
    _vv_centroid = vv_centroid;
    _v_rmse = v_rmse;
    return true;
}

bool OptiStored::get_data_pointers(std::string& _type, size_t& _nb_iter,
	vector<vector<size_t> >	& _vv_thresh, 
	vector<vector<size_t> >	& _vv_centroid, 
	vector<float>		& _v_rmse,
        vector<float>		& _v_rmsegt){
    if(!h_thresh||!h_rmse||!h_centroid||!h_rmsegt){return false;}
    if(type==""){return false;}
    
    _type = type; _nb_iter=nb_iterations;
    _vv_thresh = vv_thresh; 
    _vv_centroid = vv_centroid;
    _v_rmse = v_rmse;
    _v_rmsegt = v_rmsegt;
    return true;
}


bool OptiStored::get_rmse_at_iter(const size_t & iter, float & rmse){
    if(iter>=nb_iterations||!h_rmse||!h_thresh||!h_centroid){return false;}
    rmse = v_rmse[iter];
    return true;
}

bool OptiStored::get_rmse(vector<float> & _v_rmse){
    if(!h_rmse){return false;}
    _v_rmse = v_rmse;
    return true;
}

bool OptiStored::get_type(std::string & _type){
    _type = type;
    return true;
}

size_t OptiStored::get_iter(void){
    return nb_iterations;
}

bool OptiStored::set_rmsegt(vector<float> _v_rmsegt){
    if(h_rmsegt) return false;
    h_rmsegt=true;
    v_rmsegt=_v_rmsegt;
    return true;
}




/////////////////////////////
/////////////////////////////
/////////////////////////////
OptiPlan::OptiPlan():type(""),nb_labels(0),nb_pixels(0),h_thresh(false),h_centroid(false),h_histogram(false),h_histogram0(false),h_rmse(false),h_rmsegt(false),nb_storedplans(0), stored_set(){}
// initialisation

OptiPlan::~OptiPlan(){
    if(h_thresh){
        //delete [] a_thresh;
        //delete [] aa_thresh;
    }
    if(h_centroid){
        //delete [] a_centroid;
        //delete [] aa_centroid;
    }
    //if(h_histogram){delete [] a_histogram;}
    //if(h_rmse){delete [] a_rmse;}
}
bool OptiPlan::set_logs(IOWizard* _ioWizard, MyLog* _myLog){
    ioWizard = _ioWizard;
    myLog = _myLog;
}

bool OptiPlan::set_groundtruth(cv::Mat1i gt_label_mat){
    gt_label_mat.copyTo(gtl_mat);
    return true;
}

bool OptiPlan::set_param(string _type, size_t _labels, size_t _pixels, vector<size_t> histogram, bool store, bool reset){
    if(type == ""){type=_type;}
    else {return error("type already set");}
    nb_labels = _labels;
    upperbound_iterations = floor(log2(nb_labels-1)+4);
    nb_pixels = _pixels;
    if(histogram.size()!=nb_labels){return error("histogram size  "+to_string2(histogram.size())+" , label size "+to_string2(nb_labels) );}
    if(!h_histogram0){
        v_histogram0 = histogram;
        h_histogram0  =  true;
    }
    //else {return error("histogram already set");}
    
    if(type == "binary")  
        pIT_search_thresh=&OptiPlan::IT_search_binary;
    if(type == "binary_v2")  
        pIT_search_thresh=&OptiPlan::IT_search_binary_v2;
    if(type == "otsu") 
        pIT_search_thresh=&OptiPlan::IT_search_otsu;
    if(type == "otsu_v0") 
        pIT_search_thresh=&OptiPlan::IT_search_otsu_v0;
    if(type == "median")
        pIT_search_thresh=&OptiPlan::IT_search_median;
    if(type == "2means")
        pIT_search_thresh=&OptiPlan::IT_search_2means;
    if(type == "1surX")
        pIT_search_thresh=&OptiPlan::IT_search_1surX;

    IT_sets_centroids = type=="otsu"||type=="binary_v2"||type=="2means";

    if(pIT_search_thresh==NULL){return error("pIT_search init failed");}
    
    set_histomod();
    CPING("set_thresh "+type);
    set_thresh();
    CPING("ComputeRMSE");
    compute_RMSE();
    if(store)store_setting(reset);
    
    
}

bool OptiPlan::set_histomod(void){
    v_histogram=v_histogram0;
    if(type == type){ //TODO implement various methods
        h_histogram=true;
        for(int i=0; i<nb_labels; i++){
            v_histogram[i] += 1;
        }
    }
    return true; 
}


bool OptiPlan::set_thresh(void){
    // must set nb_iterations somewhere
    if(type == type) // algo recherche itérative //TODO choose
    {
        // initialisation
        vv_thresh = vector<vector<size_t> >();
        vv_centroid = vector<vector<size_t> >();
        int iter=0;
        vv_thresh.resize(iter+1);          vv_centroid.resize(iter+1);
        vv_thresh[0].resize(pow(2,iter)+1);
        vv_centroid[0].resize(pow(2,iter)+1);//last component unused
        vv_thresh[0][0]=0;
        vv_thresh[0][1]=nb_labels;//n'est jamais atteint
        size_t _min,_max;
        bool fin = false;
        bool full_active = true;
        if(IT_sets_centroids){
            (this->*pIT_search_thresh)
	      (0, nb_labels-1, vv_centroid[iter][0],_min,_max);
        }


        //COUT("starting set thresh");
        while(!fin)
        {
            fin = true;
            ++iter;

            vv_thresh.resize(iter+1);   vv_centroid.resize(iter+1);
            vv_thresh[iter].resize(pow(2,iter)+1);
            vv_centroid[iter].resize(pow(2,iter)+1);
            vv_thresh[iter][pow(2,iter)]=
			vv_thresh[iter-1][pow(2,iter-1)];
            for(int i = 0; i<vv_thresh[iter-1].size()-1; i++)
            {
                if(i==1&&iter==4)verbose=true;else verbose=false;
                //CPINGIF2("accès à ", i, iter, 4);
                // fill in known values for v_threshold
                vv_thresh[iter][2*i]=vv_thresh[iter-1][i];
                // fill v_threshold[iter][2*i+1]
                _min = vv_thresh[iter-1][i];
                _max = vv_thresh[iter-1][i+1]-1;
                //CPINGIF2("acbis à ", i, iter, 4);
                //CPINGIF4("i, min,max",i,_min,_max,(iter<=4&&i<=1&&type=="otsu"));
                if(_min>=_max||
                  vv_thresh[iter-1][i]==vv_thresh[iter-1][i+1]){
                    //CPINGIF2("full_active set false",iter,iter,4);
                    full_active=false;
                    if(!IT_sets_centroids){
                        vv_centroid[iter-1][i]=_min;
                    }else {
	                vv_centroid[iter][2*i]=_min;//vv_centroid[iter-1][i];
                        vv_centroid[iter][2*i+1]=_min;//vv_centroid[iter-1][i];
                    }
                    //vv_thresh[iter][2*i+1]=_min; //WARNING
                    vv_thresh[iter][2*i+1]=min(vv_thresh[iter-1][i+1],nb_labels-1);
                }
                else{
                if(_min!=_max){
                    fin=false;
                    //CPINGIF4(min,max,iter,4);
                    (this->*pIT_search_thresh)
		(_min, _max, vv_thresh[iter][2*i+1],
		vv_centroid[iter][2*i],vv_centroid[iter][2*i+1]);
                    if(!IT_sets_centroids)
                        {vv_centroid[iter-1][i]=vv_thresh[iter][2*i+1];}
                }}
            }
            vv_centroid[iter-1][vv_thresh[iter-1].size()-1]=66666;
            if(iter>=upperbound_iterations){fin=true;}
        }
        //if(!IT_sets_centroids){
        //    for(int i=0; i<vv_thresh[iter].size()-1;i++){
        //        vv_centroid[iter][i]=vv_thresh[iter][i];
        //    }
        //}
        this->nb_iterations=iter;
    }

    h_thresh=true;
    h_centroid=true;
    return true;
}

/////// Outils de stockage - évaluation
bool OptiPlan::compute_RMSE(void){
    CPING("in rmse");
    if(h_rmse){return error("rmse_handler_already_set");}
    if(nb_iterations==0){return error("nb_iterations==0");}
    //a_rmse = new float[nb_iterations]; 
    v_rmse.resize(nb_iterations);
    h_rmse = true;
    CPING("rmse ok start");
    for(size_t iter=0; iter<nb_iterations; iter++){ //mod
        v_rmse[iter]=0;
        for(size_t path=0; path<vv_thresh[iter].size()-1; path++)
        {
            for(size_t label=vv_thresh[iter][path]; label<vv_thresh[iter][path+1]; label++)
            {
                v_rmse[iter]+=v_histogram[label]*
			pow((long int)label-(long int)vv_centroid[iter][path],2);
            }
        }
        v_rmse[iter]=sqrt(v_rmse[iter])/(float)nb_pixels;
    }
    
    return true;
}


bool OptiPlan::store_setting(bool reset){
    if(!h_thresh||!h_histogram||!h_rmse)
		{return error("store_setting_error");}
    stored_set.resize(nb_storedplans+1);
    if(!stored_set[nb_storedplans].copy_from_plan(type,  nb_storedplans,
		nb_iterations, vv_thresh, vv_centroid, v_rmse))
         {return error("copy_from_plan returned false");}

    nb_storedplans++;
    if(!reset){return true;}
    return reset_memory();
}

bool OptiPlan::load_method(size_t idx_storage){
    if(h_thresh||h_histogram||h_rmse){
        myLog->a("Warning, loading a method over a pre_existing one");
        reset_memory();
    }
    if(idx_storage>=nb_storedplans){
        return error("loading idx>nb_stored");}
    
    stored_set[idx_storage].get_data_pointers
	(type,nb_iterations,vv_thresh,vv_centroid,v_rmse);
    h_thresh = true; h_centroid=true; h_rmse=true;
    
    return true;
}

bool OptiPlan::reset_memory(void){
    bool output;
    if(!h_thresh && !h_histogram && !h_rmse)
	{output=error("Non critical : nothing to erase");}
    if(h_thresh) { vv_thresh.resize(0); h_thresh = false; }
    if(h_centroid){ vv_centroid.resize(0); h_centroid = false; }
    if(h_histogram){v_histogram.resize(0); h_histogram = false;}
    if(h_rmse){v_rmse.resize(0); h_rmse = false;}
    pIT_search_thresh = NULL;
    type = ""; nb_iterations = 0; 
    return true;
}

/////////////////////////////////////////////////
/////// Algorithmes de recherche itératifs //////
/////////////////////////////////////////////////
bool OptiPlan::IT_search_binary(size_t _min, size_t _max, size_t& thresh, size_t& c1, size_t& c2){

    thresh = (size_t)((((float)(_min+_max+1))/2.0f));
    thresh = max(_min+1,thresh);
    return true;
}

bool OptiPlan::IT_search_binary_v2(size_t _min, size_t _max, size_t& thresh, size_t& c1, size_t& c2){
    
    thresh = (size_t)((((float)(_min+_max+1))/2.0f));
    thresh = max(_min+1,thresh);
    size_t nbA=0,nbB=0;
    float meanA=0.0f;
    float meanB=0.0f;
    for(size_t i=_min;i<thresh;i++){
        meanA+=v_histogram[i]*i;
        nbA+=v_histogram[i];
    } 
    for(size_t i=thresh;i<=_max;i++){
        meanB+=v_histogram[i]*i;
        nbB+=v_histogram[i];
    }
    c1=(meanA/(float) nbA)+0.5f;
    c2=(meanB/(float) nbB)+0.5f;
    return true;
}



bool OptiPlan::IT_search_median(size_t _min, size_t _max, size_t& thresh, size_t& c1, size_t& c2){
    //trouver la valeur médiane :
    // get nb_px dans l'intervalle [_min,max];
    // trouver label médian dans ]min,max]
    size_t nb_AB=0;
    for(size_t arg=_min;arg<=_max;arg++){
        nb_AB+=v_histogram[arg];
    }
    size_t target = ((float)(nb_AB)+0.5f)/2.0f;
    size_t rank = 0;
    size_t arg;
    for(arg=_min; arg<=_max && rank<target; arg++){
        rank+=v_histogram[arg];
    }
    thresh = max(_min+1,arg-1);
    return true;
}



bool OptiPlan::IT_search_otsu(size_t _min, size_t _max, size_t& thresh, size_t& c1, size_t& c2){
    //trouver le seuil d'otsu pour le segment 
    if(_min==_max){return error("in otsu, min=max");}
    float mumax;
    float mu;
    float sumA,sumB,meanA,meanB,nb_all_inv_2;
    size_t nb_A,nb_B;
    size_t  tmp_thresh=_min+1;
    nb_A  = v_histogram[_min];
    sumA  = v_histogram[_min]*_min;
    meanA = _min;
    nb_B  = 0;
    sumB  = 0;
    for(size_t i=_min+1; i<=_max; i++){
        nb_B += v_histogram[i];
        sumB += v_histogram[i]*i;
    }
    if(nb_A+nb_B ==0){
        thresh = max(_min+1,(size_t)((( (float)(_min+_max+1))/2.0f)));
        c1 = max(_min,(size_t)((( (float)(_min+thresh+1))/2.0f)));
        c2 = max(_min+1,(size_t)((( (float)(_max+thresh+1))/2.0f)));
        return error("*not critical* - Otsu void interval, thresh set binary");
    }
    if(nb_B ==0){
        thresh = tmp_thresh;
        c1 = _min;
        c2 = max(_min+1,(size_t)((( (float)(_max+thresh+1))/2.0f)));
        return true;
    }

    nb_all_inv_2  = (float)1/(float)(nb_A+nb_B);
    nb_all_inv_2 *= nb_all_inv_2;
    meanB = sumB/nb_B;  
    
    mumax = (float)nb_A*(float)nb_B*nb_all_inv_2*pow(meanA-meanB,2);
    thresh=tmp_thresh;
    c1=(meanA+0.5f);c2=(meanB+0.5f);
    //CPINGIF4("\n\n\n_min,_max,nb_all_inv_2",_min,_max,nb_all_inv_2,verbose);
    //CPINGIF4("meanA,meanB,mumax",meanA,meanB,mumax,verbose);
    for(tmp_thresh=_min+2; tmp_thresh<=_max; tmp_thresh++){
        nb_A += v_histogram[tmp_thresh-1];
        nb_B -= v_histogram[tmp_thresh-1];
        sumA += v_histogram[tmp_thresh-1]*(tmp_thresh-1);
        sumB -= v_histogram[tmp_thresh-1]*(tmp_thresh-1);
        
        meanA = (nb_A!=0) ? sumA/(float)nb_A : tmp_thresh-1;
        meanB = (nb_B!=0) ? sumB/(float)nb_B : tmp_thresh;
        mu = (float)nb_A*(float)nb_B*nb_all_inv_2*pow(meanA-meanB,2);

        //CPINGIF4("tmp_thresh ,nb_A,nb_B",tmp_thresh ,nb_A,nb_B,verbose);
        //CPINGIF4("meanA,meanB,mu",meanA,meanB,mumax,verbose);
        if(mu>mumax){
            thresh = tmp_thresh;
            mumax  = mu        ;
            c1     = (meanA+0.5f);
            c2     = (meanB+0.5f);
        }
    }
    return true;
}

bool OptiPlan::IT_search_otsu_v0(size_t _min, size_t _max, size_t& thresh, size_t& c1, size_t& c2){
    size_t fc1,fc2;
    return IT_search_otsu(_min, _max, thresh, fc1, fc2);
}


bool OptiPlan::IT_search_2means(size_t _min, size_t _max, size_t& thresh, size_t& c1, size_t& c2){
    float tmpthresh =((((float)(_min+1+_max))/2.0f));
    long int nb_A(0),nb_B(0);
    float sumA(0.0f), sumB(0.0f),meanA,meanB;
    float newthresh=tmpthresh;
    bool new_is_sup;
    //CPING("in");

    for(size_t i=_min; i<tmpthresh; i++){
        nb_A+=v_histogram[i];
        sumA+=v_histogram[i]*i;
        meanA=sumA/(float)nb_A;
    }
    for(size_t i=tmpthresh; i<=_max; i++){
        nb_B+=v_histogram[i];
        sumB+=v_histogram[i]*i;
        meanB=sumB/(float)nb_B;
    }
    newthresh =  ((meanA+meanB+1.0f)/2.0f);
    //ComputeCentroid
    while(newthresh!=tmpthresh){
        new_is_sup = (newthresh>tmpthresh);
        for(long int i=(long int)min(tmpthresh,newthresh);i<(long int)max(tmpthresh,newthresh);i++){
            long int tmp=(long int)v_histogram[i]; //TODO check consistency.
            nb_A+=  new_is_sup?tmp   :-tmp;
            nb_B+= !new_is_sup?tmp   :-tmp;
            sumA+=  new_is_sup?tmp*i :-tmp*i;
            sumB+= !new_is_sup?tmp*i :-tmp*i;
        }
        meanA=sumA/(float)nb_A;
        meanB=sumB/(float)nb_B;
        tmpthresh = newthresh;
        newthresh = ((meanA+meanB+1.0f)/2.0f);
    }
    c1=meanA;
    c2=meanB;
    thresh = newthresh;
    return true;
}


bool OptiPlan::IT_search_1surX(size_t _min, size_t _max, size_t& thresh, size_t& c1, size_t& c2){
    
    
    thresh = _min+1;
    return error("1 sur X not implemented yet");
    return true;
}


////////////////////////////
//      CONSULTATION      //
////////////////////////////
bool OptiPlan::get_best_method_at_it(size_t iter, size_t & idx_storage, std::string& type, bool& active){
    type="";
    if(nb_storedplans==0){return error("No storedplans to compare");}
    bool rmseinit=false;
    active=true;
    float rmse,rmse_best;
    size_t idx_plan=0;
    idx_storage=0;
    while(!rmseinit && idx_plan<nb_storedplans){
        if(stored_set[idx_plan].get_rmse_at_iter(iter,rmse)){
            rmseinit=true;
            rmse_best=rmse;
            idx_storage=idx_plan;
        }
        else { active=false; }
        idx_plan++;
    }
    while(idx_plan<nb_storedplans){
        if(stored_set[idx_plan].get_rmse_at_iter(iter,rmse)){
            if(rmse<rmse_best){
                rmse_best=rmse;
                idx_storage=idx_plan;
            }
        }
        else { active=false; }
        idx_plan++;
    }
    
    stored_set[idx_storage].get_type(type);
    if(!rmseinit){return error("No method at it, to get RMSE from");}
    
    return true;
}


bool OptiPlan::get_best_method_at_it_pointers(size_t iter,
	std::string& _type, size_t& _nb_iter,
	vector<vector<size_t> >	& _vv_thresh, 
	vector<vector<size_t> >	& _vv_centroid, 
	vector<float>		& _v_rmse){
    if(nb_storedplans==0){return error("no plan stored for bestmatit");}
    size_t idx_bestplan;
    bool active;
    if(!get_best_method_at_it(iter,idx_bestplan,_type,active))
        {return error("Get Best Method failed");}
    if(!active){return false;}
    return stored_set[idx_bestplan].get_data_pointers
	(_type,_nb_iter,_vv_thresh,_vv_centroid,_v_rmse);
}

bool OptiPlan::get_thresholds_n_centroids(int method,
			std::string & _type,
			size_t & _nb_iterations, 
			vector<vector<size_t> > & _vv_thresh, 
			vector<vector<size_t> > & _vv_centro){
    vector<float> _v_rmse;
    if(method>(int)nb_storedplans){return error("no such storedplans (get_t&c)");}
    if(method<-1){return error("asked -2 in get_tnc");}
    if(method!=-1){
        return(stored_set[method].get_data_pointers(
		_type,_nb_iterations,_vv_thresh,_vv_centro,_v_rmse));
    }
    if(h_thresh && h_centroid){
        _vv_thresh = vv_thresh;
        _vv_centro = vv_centroid;
        _type      = type;
        _nb_iterations = nb_iterations;
        return true;
    }
    return error("nothing to do");
}

bool OptiPlan::get_nb_storedmethods(size_t & _nb_storedplans){
    _nb_storedplans = nb_storedplans;
    return true;    
}

////////////////////////////
//      TRAITEMENT        //
////////////////////////////

bool OptiPlan::get_ThreshedMatrix(const cv::Mat1i & mat_in, cv::Mat1i & mat_out, size_t iteration, int method){
    // if method=-1, on utilise la méthode stockée en mémoire, si elle existe.
    std::string _type;
    size_t _nb_iterations;
    vector<vector<size_t> > _vv_thresh;
    vector<vector<size_t> > _vv_centroid;
    vector<float>   _v_rmse;
    //gestion des erreurs et chargement des données
    if(method==-1){
        if(!h_thresh||!h_centroid||!h_rmse){
            return error("ThreshMatrix: failure, no plan stored");
        } else {
            _type = type;_nb_iterations=nb_iterations; _vv_thresh = vv_thresh;
            _vv_centroid = vv_centroid; _v_rmse = v_rmse;
    }} else {
        if(method>=nb_storedplans||method<-1){
            return error("ThreshedMatrix: invprt method |"+to_string2(method));
        } else {
            stored_set[method].get_data_pointers(
              _type, _nb_iterations, _vv_thresh, _vv_centroid, _v_rmse);
    }   }
    if(iteration>=_nb_iterations){
        return error("ThreshedMatrix : invprt iteration |"+to_string2(iteration));}
    // renvoie la matrice des labels (parmi les labels interpolés)
    mat_out = cv::Mat::zeros(mat_in.rows,mat_in.cols,CV_32S);
    cv::Mat1i m_zero = cv::Mat::zeros(mat_in.rows,mat_in.cols,CV_32S);
    cv::Mat1i m_hold = cv::Mat::zeros(mat_in.rows,mat_in.cols,CV_32S);
    cv::Mat1i m_mask = cv::Mat::zeros(mat_in.rows,mat_in.cols,CV_32S);
    size_t tempvalA,tempvalB;
    //int   temprank;
    for(int i=1; i<_vv_thresh[iteration].size();i++){
        m_mask = m_mask*0;
        //temprank = _vv_centroid[iteration][i-1];
        
        tempvalA=_vv_thresh[iteration][i-1];
        tempvalB=_vv_thresh[iteration][i];
        m_hold  = m_zero + _vv_centroid[iteration][i-1];
        m_hold.copyTo(m_mask,  (mat_in<tempvalB) );
        m_zero.copyTo(m_mask,  (mat_in<tempvalA) );
        m_mask.copyTo(mat_out, (m_mask>0) );
    }
    return true;
}

////////////////////////////
//      VISUALISATION     //
////////////////////////////

bool OptiPlan::show_RMSE(const std::string& filename){
    // s'appuie sur a_rmse pour les données stockées
    FILE* gnuplot = popen("gnuplot","w");
    fprintf(gnuplot, "set term 'pngcairo' size 950,600 enhanced font 'Verdana,10'\n");
    ioWizard->set_gnuplot_output(gnuplot,filename+".png");
    //
    fprintf(gnuplot, "plot '-' with lines\n");
    //COUT(nb_iterations);
    //CPING2("rmsesize",v_rmse.size());
    for(int k=1; k<nb_iterations; k++)
    {

        fprintf(gnuplot, "%i %g\n", k, v_rmse[k]);
    }
    fflush(gnuplot);fprintf(gnuplot, "e\n");
    fprintf(gnuplot,"unset output \n");
    //
    fprintf(gnuplot,"exit \n"); 
    pclose(gnuplot);

    return true;
}

///// BEGIN TODO recoder show_rmse en plus simple
bool OptiPlan::show_all_RMSE(const std::string& filename){
    //même chose que show mais pour toutes les méthodes. 
    FILE* gnuplot = popen("gnuplot","w");
    fprintf(gnuplot, "set term 'pngcairo' size 950,600 enhanced font 'Verdana,10'\n");
    ioWizard->set_gnuplot_output(gnuplot,filename+".png");

    std::string initplot_msg = "plot ";
    std::string _tmptype;
    size_t nb_iter;
    vector<vector<size_t> > vvs1;
    vector<vector<size_t> > vvs2;
    vector<float> vf;
    if(!stored_set[0].get_data_pointers(_tmptype, nb_iter, 
      vvs1, vvs2, vf))
        {return error("get_data_pointers at step showall0");}
    initplot_msg += "'-' title '"+_tmptype+"' with lines";
    for(int k=1;k<nb_storedplans; k++){
        if(!stored_set[k].get_data_pointers(_tmptype, nb_iter,
	  vvs1, vvs2, vf) )
            {return error("get_data_pointers at step showall");}
        //_tmptype="-";
        initplot_msg += ", '-' title '"+_tmptype+"' with lines";
    }
    initplot_msg += "\n";
    fprintf(gnuplot, "%s",initplot_msg.c_str());


    for(int k=0;k<nb_storedplans;k++){
        show_RMSE_elt_n(gnuplot,k);
    }
    
    fprintf(gnuplot,"unset output \n");
    fprintf(gnuplot,"exit \n"); 
    pclose(gnuplot);
    
    return true;
}

bool OptiPlan::show_RMSE_elt_n(FILE* gnuplot, size_t idx){
    //
    std::string type;
    size_t nb_iter;
    vector<float> vector_rmse;
    vector<vector<size_t> > vvs1, vvs2;
    
    if(!stored_set[idx].get_data_pointers(type, nb_iter, vvs1, vvs2, vector_rmse)){return error("get_data_pointers fail at show_rmse_elt_" + to_string2(idx) );}
    for(int i=1; i<nb_iter; i++){
        fprintf(gnuplot, "%i %g\n", i, vector_rmse[i]);
    }
    fflush(gnuplot);
    fprintf(gnuplot, "e\n");
    return true;
}
//////END TODO
//////Start TODO.ing

bool OptiPlan::show_all_RMSE2(const std::string& filename){
    FILE* gnuplot = popen("gnuplot","w");
    fprintf(gnuplot, "set term 'pngcairo' size 950,600 enhanced font 'Verdana,10'\n");
    ioWizard->set_gnuplot_output(gnuplot,filename+".png");

    std::string initplot_msg = "plot ";
    std::string _tmptype;
    size_t nb_iter;
    vector<vector<size_t> > vvs1;
    vector<vector<size_t> > vvs2;
    vector<float> _v_rmse;
    vector<float> _v_rmsegt;
    if(!stored_set[0].get_data_pointers(_tmptype, nb_iter, 
      vvs1, vvs2, _v_rmse, _v_rmsegt))
        {return error("get_data_pointers at step showall0");}
    initplot_msg += "'-' title '"+_tmptype+"' with lines";
    initplot_msg += ",'-' title '"+_tmptype+" vs gt' with lines";
    for(int k=1;k<nb_storedplans; k++){
        if(!stored_set[k].get_data_pointers(_tmptype, nb_iter,
	  vvs1, vvs2, _v_rmse, _v_rmsegt) )
            {return error("get_data_pointers at step showall");}
        //_tmptype="-";
        initplot_msg += ", '-' title '"+_tmptype+"' with lines";
        initplot_msg += ", '-' title '"+_tmptype+" vs gt' with lines";
    }
    initplot_msg += "\n";
    fprintf(gnuplot, "%s",initplot_msg.c_str());


    for(int k=0;k<nb_storedplans;k++){
        stored_set[k].get_data_pointers(_tmptype, nb_iter,
	  vvs1, vvs2, _v_rmse, _v_rmsegt);        
        gnuplot_vect(gnuplot,_v_rmse);   
        gnuplot_vect(gnuplot,_v_rmsegt);
    }
    
    fprintf(gnuplot,"unset output \n");
    fprintf(gnuplot,"exit \n"); 
    pclose(gnuplot);
    return true;
}

bool OptiPlan::gnuplot_vect(FILE* gnuplot, vector<float> vect){
    for(int i=1; i<vect.size(); i++){
        fprintf(gnuplot, "%i %g\n", i, vect[i]);
    }
    fflush(gnuplot);
    fprintf(gnuplot, "e\n");
    return true;
}

bool OptiPlan::show_all_thresh_plans(std::string filename){
    ioWizard->mksubdir("threshplans");
    filename=string("threshplans/")+filename;
    if(nb_storedplans<1){return error("no plan set! no thresh to show");}
    for(int k = 0; k<nb_storedplans; k++){
        show_thresh_plan(filename,k);
    }
}

bool OptiPlan::show_thresh_plan(std::string filename, int kplan){
    // if kplan =-1, then loads internal data (if existing)
    bool see_arrows = true;
    std::string t_type;
    size_t t_nb_iter;
    vector<vector<size_t> > t_vv_threshold,t_vv_centroid;
    vector<float> t_v_rmse;
    if(!(kplan>=-1&&kplan<nb_storedplans)){
        return error("calling a inexistant plan n°"+to_string2(kplan));
    }
    else if(kplan!=-1){
        stored_set[kplan].get_data_pointers(t_type,t_nb_iter,
	t_vv_threshold,t_vv_centroid,t_v_rmse);
    }
    else if(!h_thresh||!h_centroid||!h_histogram0){
        return error("handler not set at show thresh plan");
    }
    else{
        t_type        = type;
        t_nb_iter     = nb_iterations;
        t_vv_threshold= vv_thresh;
        t_vv_centroid = vv_centroid;
        t_v_rmse      = v_rmse;
    }
    
    FILE* gnuplot = popen("gnuplot","w");
    fprintf(gnuplot,"set terminal pngcairo size 950,600 enhanced font 'Verdana,10'\n");
    ioWizard->set_gnuplot_output(gnuplot,filename+t_type+".png");
    std::string tmp_string="set border linewidth 1.5\n"+
  to_string2("")+
  "set style line 1 linecolor rgb '#dd181f' linetype 1 linewidth 2\n"+
  "set style line 2 linecolor rgb '#0060ad' linetype 1 linewidth 2 pt 4\n"+
  "set style line 3 linecolor rgb '#dd181f' linetype 1 linewidth 2 pt 6\n"+
  "set style line 4 linecolor rgb '#aaffaa' linetype 1 linewidth 1\n"+
  "set key at "+to_string2(nb_labels-10)+",3\n"+
  "set title 'Optimization with "+to_string2(t_type)+"' font 'Helvetica,15' enhanced\n"+
  "set xlabel 'label'\n"+
  "set ylabel 'iteration'\n"+
  "set xrange[0:"+to_string2(nb_labels)+"]\n"+
  "set yrange["+to_string2(-(int)(t_nb_iter))+":"+to_string2(t_nb_iter/2+1)+"]\n"+
  "set xtics 4\n"+
  "set ytics 1\n"+
  "set tics scale 0.75\n"+
  "plot '-' title 'histogram' with lines ls 1, "+
       "'-' title 'threshold' with p ls 2, "    +
       "'-' title 'centroids' with p ls 3";
    if(see_arrows){
        tmp_string+=
       ", '-' notitle with lines ls 4";
    }
    tmp_string+="\n";

    fprintf(gnuplot,tmp_string.c_str());
    float scale_histogram;
    for(size_t i=0; i<nb_labels; i++)
    {
        scale_histogram = v_histogram0[i]>scale_histogram?
			  v_histogram0[i]:scale_histogram;
    }
    scale_histogram = ((float)t_nb_iter)*0.50f/scale_histogram;
    for(size_t i=0; i<nb_labels; i++)
    {
        fprintf(gnuplot, "%lu %g \n", i, 
			((float)v_histogram0[i])*scale_histogram );
    }
    fflush(gnuplot);
    fprintf(gnuplot, "e\n");

    for(int i=0; i<t_nb_iter; i++)
    {
        for(size_t path=0; path<pow(2,i)+1; path++)
        {
            fprintf(gnuplot,"%lu %i \n", t_vv_threshold[i][path], -i);
        }
    }
    fflush(gnuplot);
    fprintf(gnuplot, "e\n");

    for(int i=0; i<t_nb_iter; i++)
    {
        for(size_t path=0; path<pow(2,i); path++)
        {
            fprintf(gnuplot,"%lu %i \n", t_vv_centroid[i][path], -i);
        }
    }
    fflush(gnuplot);
    fprintf(gnuplot, "e\n");

    if(see_arrows){
    for(int i=0; i<t_nb_iter-1; i++)
        {
            for(size_t path=0; path<pow(2,i); path++)
            {
                fprintf(gnuplot,"%lu %i \n", t_vv_centroid[i][path], -i);
                fprintf(gnuplot,"%lu %i \n", t_vv_centroid[i+1][2*path], -i-1);
                fprintf(gnuplot,"\n");
                fprintf(gnuplot,"%lu %i \n", t_vv_centroid[i][path], -i);
                fprintf(gnuplot,"%lu %i \n", t_vv_centroid[i+1][2*path+1], -i-1);
                fprintf(gnuplot,"\n");
            }
        }
        fflush(gnuplot);
        fprintf(gnuplot, "e\n");
    }
    fprintf(gnuplot,"unset output \n");
    fprintf(gnuplot,"exit \n"); 
    pclose(gnuplot);
    return true;
}



bool OptiPlan::write_all_ThreshedMatrix(cv::Mat1i & mat_in, std::string folder, bool set_rmsegt){
    // pour tout plan enregistré
    // pour toutes les itérations
    // appeler la fonction getThreshedMatrix
    cv::Mat1i mat_tmp = cv::Mat::zeros(mat_in.rows,mat_in.cols,CV_32S);
    if(set_rmsegt&& (gtl_mat.cols!=mat_in.cols||gtl_mat.rows!=mat_in.rows) ){
        return error("Matrix dimension must agree for rmse (gt & threshedmatrix)");
    }
    bool plan_is_stored = h_thresh && h_centroid; //disabled.
    vector<float> _v_rmse_gt; //if set_rmsegt
    ioWizard->mksubdir(folder);
    std::string _type;
    std::string str_tmp = "";
    
    for(size_t method = plan_is_stored?0:0 ; method<nb_storedplans; method++){
        if(set_rmsegt)_v_rmse_gt.resize(stored_set[method].get_iter());
        for(size_t i=1; i<stored_set[method].get_iter();i++){
            if(!get_ThreshedMatrix(mat_in,mat_tmp,i,method)){
                error("Unexpected failure in get_ThreshedMatrix");}
            if(set_rmsegt){
                _v_rmse_gt[i]=Utils::compute_rmse_label(mat_tmp, gtl_mat,nb_labels);
            }
            stored_set[method].get_type(_type);
            str_tmp = _type + "-" + to_string2(i);
            ioWizard->img_setscale(3);
            ioWizard->writeImage(folder+"/2D-"+str_tmp+ ".png",mat_tmp);
            ioWizard->write3DImage(folder+"/3D-"+str_tmp+ ".png",mat_tmp);
        }
        if(set_rmsegt)stored_set[method].set_rmsegt(_v_rmse_gt);
    }
    
    ioWizard->img_unsetscale();
    return true;
}

bool OptiPlan::computeCrossRMSEperf_andLog(void){
    if(nb_storedplans==0||upperbound_iterations==0){
        return error("No storedplans or iterations for XRMSE");}
    vector<vector<vector<std::string> > > vvv_deltaRMSEmatrix;
    vector<vector<std::string> >    vv_types;
    std::string typeA;
    std::string typeB;
    vector<float> v_rmseA, v_rmseB;
    vvv_deltaRMSEmatrix.resize(nb_storedplans);
    vv_types.resize(nb_storedplans);
    size_t miniter;
    for(int m=0; m<nb_storedplans; m++){
        vvv_deltaRMSEmatrix[m].resize(m);
        vv_types[m].resize(m);
        stored_set[m].get_type(typeA);
        stored_set[m].get_rmse(v_rmseA);
        for(int n=0; n<m; n++){
            vvv_deltaRMSEmatrix[m][n].resize(upperbound_iterations);
            stored_set[n].get_type(typeB);
            stored_set[n].get_rmse(v_rmseB);
            vv_types[m][n]=typeA+"-"+typeB;
            miniter=min(v_rmseA.size(),v_rmseB.size());
            for(int i=0; i<min(miniter,upperbound_iterations); i++){
                vvv_deltaRMSEmatrix[m][n][i]=to_string2(v_rmseA[i]-v_rmseB[i]);
                }
            for(int i=min(miniter,upperbound_iterations); 
			i<upperbound_iterations; ++i){
                vvv_deltaRMSEmatrix[m][n][i]="";
                }
            }
        }
    if(!myLog->write_deltaRMSEtoHistogram(vvv_deltaRMSEmatrix,vv_types)){
        return error("failure in XRMSE, size m<=1 ???");
    }
    return true;
}


bool OptiPlan::addToLog(void){
    std::vector<std::string> v_types;
    v_types.resize(0);
    size_t iter=0;
    size_t idx;
    std::string _type;
    bool active=true;
    while(active)
    {
        get_best_method_at_it(iter,idx,_type,active);
        iter++;
        v_types.push_back(_type);
    }
    
    myLog->set_bestplans(v_types);
    return true;
}



bool OptiPlan::error(void){
    myLog->av("failure in optiplan!\n");
    return false;
}
bool OptiPlan::error(const std::string& text){
    myLog->av("failure in optiplan!\n"+text+"\n");
    return false;
}

