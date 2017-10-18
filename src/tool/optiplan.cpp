/****************************
*	ProjSFF
*	optiplan.cpp
*	cstoribal
*	13-10-17
****************************/

#include "../tool/optiplan.h"



OptiPlan::OptiPlan():type(NULL),nb_labels(0),nb_pixels(0),h_thresh(false),h_histogram(false),a_thresh(NULL),a_histogram(NULL){}
// initialisation

OptiPlan::~OptiPlan(){
    if(h_thresh){
        delete [] a_thresh;
        delete [] aa_thresh;
    }
    if(h_histogram){delete [] a_histogram;}
}

bool OptiPlan::set_param(string _type, IOWizard* _ioWizard, MyLog* _myLog, size_t _labels, size_t _pixels, vector<unsigned int> histogram){
    if(type == NULL){type=_type;}
    else {return error();}
    ioWizard = _ioWizard;
    myLog = _myLog;
    nb_labels = _labels;
    nb_pixels = _pixels;
    if(histogram.size()!=nb_labels){return error();}
    if(!h_histogram){
        a_histogram = new size_t[nb_labels];
        for(int i=0; i<nb_labels; i++){
            a_histogram[i]=histogram[i];
        }
        h_histogram = true;
    }
    else {return error();}
    
    if(type == type) //TODO choice
        pf_search=&OptiPlan::search_kmean_in;
    
    set_histomod();
    set_thresh();
    
    
    
}

bool OptiPlan::set_histomod(void){
    if(type == type){ //TODO implement various methods
        for(int i=0; i<nb_labels; i++){
            a_histogram[i] += 1;
        }
    }
    return true; 
}


bool OptiPlan::set_thresh(void){
    if(type == "k2means") // algo recherche itérative
    {
        // initialisation
        vector<vector<size_t> > v_threshold;
        v_threshold.resize(1);
        v_threshold[0].resize(1);
        size_t min, max;
        get_cluster(vector<size_t>(),0,min,max);
        search_kmean_in(1,0,min,max,v_threshold[0]);
        bool fin = false;
        bool full_active = true;
        int iter=1;
        
        while(!fin)
        {
            fin = true;
            v_threshold.resize(iter+1);
            v_threshold[iter].resize(pow(2,iter)-1);
            for(int i = 0; i<v_threshold[iter-1].size(); i++)
                // fill in known values for v_threshold
                v_threshold[iter][2*i+1]=v_threshold[iter-1][i];
            for(int path=0; path<pow(2,iter); path++)
            {
                get_cluster(v_threshold[iter-1],path, min, max);
                (this->*pf_search)(1,path, min, max, v_threshold[iter]);
                if(min==max)full_active=false;
                if(min!=max)fin=false;
            }
        }
    }
    
}

bool OptiPlan::get_cluster(vector<size_t> thresh, int path, size_t& min, size_t& max){
    if(thresh.size()==0){
        min=0; max=nb_labels-1; 
        return true;
    }
    if(path==0) min=0;
    else min=thresh[path-1];
    if(thresh.size()==path) max = nb_labels-1;
    else max=thresh[path]-1;
    return true;
}



bool OptiPlan::search_kmean_in(int k, int path, size_t min, size_t max, vector<size_t> & means){
    
    
    
    

    //write output in means[path*2]
    return true;
}  


bool OptiPlan::search_otsu(int k, int path, size_t min, size_t max, vector<size_t> & means){
    int argmax=min+1;
    int sigmax=0;
    
    
    
    
    
    

    //write output in means[path*2]
    return true;
}  



bool OptiPlan::error(void){
    myLog->av("failure in optiplan!\n");
    return false;
}

