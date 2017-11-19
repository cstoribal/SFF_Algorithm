/****************************
*	logs.cpp
*	cstoribal
*	15-06-17
****************************/

/*
Outputs at the same time a precise log of each simulation encountered
and the results under the format of a csv
*/


#include "logs.h"

MyLog::MyLog(){
    logversion = "v0.5";
    logs = "Starting log file \n "+logversion+" \n";
}
MyLog::~MyLog(){
    delete this->log_data_out;
}

bool MyLog::set_param(const string & file, bool verb){
    this->verbose = false;
    this->filepath = file;
    cout<<filepath<<endl;
    this->print_eftypes();
    this->log_data_out = new MyLogOut;
    this->log_data_out->logversion=logversion;
    return true;
}

bool MyLog::print_eftypes(void){
    logs += "*****\nfType is : ";
    if(is_same<fType,float>::value) logs += "float";
    if(is_same<fType,double>::value) logs += "double";
    logs += "\neType is : ";
    if(is_same<eType,int>::value) logs += "int";
    if(is_same<eType,float>::value) logs += "float";
    if(is_same<eType,double>::value) logs += "double";
    logs += "\n*****\n";
    return true;
}

bool MyLog::a(const string & addstring){
    logs += addstring;
    if(this->verbose)
    {
        cout<<addstring;
    }
    return true;
}

bool MyLog::av(const string & addstring){
    logs += addstring;
    cout<<addstring;
    return true;
}

bool MyLog::as(const string & addstring){
    logs += addstring;
    return true;
}

bool MyLog::write(void){
    ofstream out(this->filepath.c_str());
    out << this->logs;
    out.close();

    log_data_out->write();
    return true;
}


double MyLog::time_r(void){
    double tmp = timer.Time();
    timer.Init();
    return tmp;
}

double MyLog::time_r(int rank){
    double tmp = timer.Time();
    timer.Init();
    if(rank>=this->log_data_out->output.time.size()){
        CPING("error time_r");
        return tmp;
    }
    this->log_data_out->output.time[rank] = tmp;
    return tmp;
}

double MyLog::time(void){
    return timer.Time();
}

double MyLog::time(int rank){
    double tmp = timer.Time();
    if(rank>=this->log_data_out->output.time.size()){
        CPING("error time_r");
        return tmp;
    }
    this->log_data_out->output.time[rank] = tmp;
    return tmp;
}

void MyLog::time_i(void){
    return timer.Init();
}

bool MyLog::clear_log(void){
//ensures data are not overwritten
    this->log_data_out->output.lambda         = -999;
    this->log_data_out->output.iterationlvl   = -999;
    this->log_data_out->output.opti           = "none";
    this->log_data_out->output.rmse           = -999;
    this->log_data_out->output.psnr           = -999;
    this->log_data_out->output.types.resize(0);
    //for(int i=6; i<this->log_data_out->output.time.size(); i++){
    //    this->log_data_out->output.time[i] = -1;
    //}
    
}

bool MyLog::set_state(fType lr, fType ld, int iter){
    this->log_data_out->output.lambda = lr/ld;
    this->log_data_out->output.iterationlvl = iter;
    return true;
}

bool MyLog::set_state(std::string opti, fType lambda, int iter){
    this->log_data_out->output.lambda         = lambda;
    this->log_data_out->output.opti           = opti;
    this->log_data_out->output.iterationlvl   = iter;
    return true;
}
    

bool MyLog::set_eval(fType rmse, fType psnr){
    this->log_data_out->output.rmse = rmse;
    this->log_data_out->output.psnr = psnr;
    return true;
}

bool MyLog::set_bestplans(std::vector<std::string> types){
    this->log_data_out->output.types = types;
    this->log_data_out->output.types.resize(8);
    for(int i=types.size(); i<8; i++){
        this->log_data_out->output.types[i] = "end";
    }
    return true;
}

bool MyLog::write_deltaRMSEtoHistogram(vector<vector<vector<std::string> > > & vvv_deltarmse, vector<vector<std::string> > & vv_type12){
    return this->log_data_out->write_deltaRMSEtoHistogram(vvv_deltarmse, vv_type12);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

MyLogOut::MyLogOut(){};MyLogOut::~MyLogOut(){};

bool MyLogOut::setup(tdf_input & input){
    this->output.settings = &input;
    this->output.time.resize(15);
    this->output.type_best_theorical="";
    this->output.type_best_regularized="";
    return true;
}

bool MyLogOut::Format_txt(void){
    tdf_log & o = output;
    tdf_input* & in = o.settings;
    int tmp_nimg = (int)floor((double)(in->file1_lasti-in->file1_firsti+1) / (double)in->file1_deltai);

    strdata = logversion + "; ";
    strdata += 
	  (in->outputf_set?to_string2(in->outputfolder)	:"") + "; "
        + (in->file1_set?  to_string2(in->file1_path)	:"") + "; " 
	+ (in->file1_set?  to_string2(in->file1_ext)	:"") + "; " 
	+ (in->file1_set?  to_string2(tmp_nimg) 	:"") + "; "
	+ (in->file2_set?  to_string2(in->file2[0])	:"") + "; "
	+ (in->file2_set?  to_string2(in->file2.size())	:"") + "; "
	+ (in->groundt_set?to_string2(in->gtpath)	:"") + "; "
	+ (in->groundt_set?to_string2(in->gtb)  	:"") + "; "
	+ (in->groundt_set?to_string2(in->gta+in->gtb) 	:"") + "; "
	+ (in->focus_set?  to_string2(in->focus[0])	:"") + "; "
	+ (in->focus_set?to_string2(in->focus[in->focus.size()-1]):"") + "; "
	+ (in->preproc_set?to_string2(in->scale)	:"") + "; "
 	+ (in->preproc_set?to_string2(in->gauss)	:"") + "; "
 	+ (in->preproc_set?to_string2(in->noise_a)	:"") + "; "
 	+ (in->preproc_set?to_string2(in->noise_b)	:"") + "; "
 	+ (in->preproc_set?to_string2(in->noise_ca)	:"") + "; "
 	+ (in->preproc_set?to_string2(in->noise_cs)	:"") + "; "
 	+ (in->sharp_set?  to_string2(in->sharp)	:"") + "; "
	+ (in->depth_set?  to_string2(in->depth)	:"") + "; "
	+ (in->nrj_set?    to_string2(in->nrj_d)	:"") + "; "
	+ (in->nrj_set?    to_string2(in->nrj_r)	:"") + "; "
	//+ (in->opti_set?   to_string2(in->opti) 	:"") + "; "
	+ (in->opti_set?   to_string2(in->connexity)	:"") + "; ";
    
    strdata +=
	  to_string2(o.opti) + "; "
	+ to_string2(o.lambda) + "; "
	+ to_string2(o.iterationlvl) + "; "
	+ to_string2(o.rmse) + "; "
	+ to_string2(o.psnr);

    for(int k=0;k<o.time.size();k++)
        strdata += "; " + to_string2(o.time[k]);


    for(int k=0;k<o.types.size();k++)
        strdata += ";" + to_string2(o.types[k]);

    
    strdata += "\n";
    return true;
}


bool MyLogOut::write_deltaRMSEtoHistogram(vector<vector<vector<std::string> > > & vvv_deltarmse, vector<vector<std::string> > & vv_type12){
    //si le fichier n'existe pas, le créer
    size_t M  = vvv_deltarmse.size();
    if(M<=1){return false;}
    size_t IT = vvv_deltarmse[1][0].size();
    std::string image_path = this->output.settings->file1_set==1?this->output.settings->file1_path:this->output.settings->file2[0];
    bool file_exists=false;
    {
        ifstream fichier("./deltaRMSE.csv");
        file_exists = !fichier.fail();
        fichier.close();
    }
    std::ofstream outfile;
    std::string text = "v0.3 ; ";
    outfile.open("./deltaRMSE.csv", std::ios_base::app);
    if(!file_exists){
        COUT("création d'un nouveau fichier deltaRMSE");
        text += "filepath ; ";
        for(int m=0; m<M;  m++) for(int n=0; n<m; n++) for(int it=1;it<IT; it++){
            text+=vv_type12[m][n]+"-it-"+to_string2(it)+" ; ";
        }
        text+="\nv0.3 ;";
    }
    text+= image_path+" ;";
    for(int m=0; m<M;  m++) for(int n=0; n<m; n++) for(int it=1;it<IT; it++){
            text+=vvv_deltarmse[m][n][it]+";";
    }
    text+="\n";
    outfile << text;
    outfile.close();
    return true;
}


bool MyLogOut::write(void){
    this->Format_txt();

    bool file_exists=false;
    {
        ifstream fichier("../logout.csv");
        file_exists = !fichier.fail();
        fichier.close();
    }
    std::ofstream outfile;
    outfile.open("../logout.csv", std::ios_base::app);
    if(!file_exists){
        COUT("création d'un nouveau fichier de logs");
        create_new_logfile_header(outfile);
    }
    // TODO si pas de fichier alors créer l'en-tête !!
    outfile << strdata;
    outfile.close();
    return true;
}

bool MyLogOut::create_new_logfile_header(std::ofstream & outfile){
    std::string initdata;
    initdata = logversion + " ; outputfolder ; file1_p ; extension ; nbimg ; file2_p ; nb_img ; gt_path ; dmin ; dmax ; fmin ; fmax ; scale ; blur ; noiseA ; noiseB ; noiseCA ; noiseCS ; sharpOp ; depthOp ; nrj_d ; nrj_r ; opti ; connexity ; opti ; lambda ; iterationlvl ; rmse ; psnr ; t ; t ; t ; t ; t ; t ;  topti ; \n";
    outfile << initdata ;
    return true; //TODO
}




