/****************************
*	logs.cpp
*	cstoribal
*	15-06-17
****************************/

/*
Outputs at the same time a precise log of each simulation encountered
and the results under the format of a csv
*/


// 0. remove psnr, sets logs rmse in line in a vector.

#include "logs.h"

MyLog::MyLog(){
    logversion = "v0.7";
    logs = "Starting log file \n "+logversion+" \n";
}
MyLog::~MyLog(){
    ofstream out(this->filepath.c_str());
    out << this->logs;
    out.close();
    
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

bool MyLog::clear_iteration_times(void){
    for(int i=6; i<this->log_data_out->output.time.size(); i++){
        this->log_data_out->output.time[i] = -1;
    }
    for(int i=1; i<this->log_data_out->output.v_rmse.size(); i++){
        this->log_data_out->output.v_rmse[i] = -1;
    }
    return true;
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


bool MyLog::set_state(std::string opti, fType lambda){
    this->log_data_out->output.lambda         = lambda;
    this->log_data_out->output.opti           = opti;
    return true;
}

bool MyLog::set_eval_at(fType rmse, int iter){
    if(iter<10){
        this->log_data_out->output.v_rmse[iter] = rmse;
        return true;
    }
    return false;
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
    this->output.time.resize(20);
    this->output.type_best_theorical="";
    this->output.type_best_regularized="";
    this->output.v_rmse.resize(10);
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
	+ to_string2(o.lambda);
	//+ to_string2(o.iterationlvl) + "; "
	//+ to_string2(o.rmse) + "; "
	//+ to_string2(o.psnr);
    
    for(int k=0; k<o.v_rmse.size(); k++){
        strdata += "; " + to_string2(o.v_rmse[k]);
    }

    for(int k=0;k<o.time.size();k++)
        strdata += "; " + to_string2(o.time[k]);


    for(int k=0;k<o.types.size();k++)
        strdata += "; " + to_string2(o.types[k]);

    
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
            text+=vvv_deltarmse[m][n][it]+" ;";
    }
    text+="\n";
    outfile << text;
    outfile.close();
    return true;
}

bool MyLogOut::write(void){
    write("../logout.csv",1);
    CPING2("outputfolder : ",output.settings->outputfolder);
    if(output.settings->outputf_set) write("./"+output.settings->outputfolder+"/logout.csv",0);
    return true;
}

bool MyLogOut::write(std::string filename, bool verbose){
    this->Format_txt();

    bool file_exists=false;
    {
        ifstream fichier(filename);
        file_exists = !fichier.fail();
        fichier.close();
    }
    std::ofstream outfile;
    outfile.open(filename, std::ios_base::app);
    if(!file_exists){
        if(verbose) {COUT("création d'un nouveau fichier de logs"); }
        create_new_logfile_header(outfile);
    }
    // TODO si pas de fichier alors créer l'en-tête !!
    outfile << strdata;
    outfile.close();
    return true;
}

bool MyLogOut::create_new_logfile_header(std::ofstream & outfile){
    std::string initdata;
    initdata = logversion + " ; outputfolder ; file1_p ; extension ; nbimg ; file2_p ; nb_img ; gt_path ; dmin ; dmax ; fmin ; fmax ; scale ; blur ; noiseA ; noiseB ; noiseCA ; noiseCS ; sharpOp ; depthOp ; nrj_d ; nrj_r ; connexity ; opti ; lambda ; rmse0 ; rmse1; rmse2; rmse3; r4 ; r5 ; r6 ; r7 ; r8 ; r9 ; t0 ; t1 ; t2 ; t3 ; t4 ; t5 ; t6_charg ; t7_opti ; t8_o2 ; t9_o3 ; t10_o4 ; t11_o5 ; t12_o6 ; t13 ; t14 ; t15 ; t16 ; t17 ; t18 ; t19 ; type 1; type 2; type 3 ; type 4 ; type 5 ; type 6 ; type 7 ; type 8 \n";
    outfile << initdata ;
    return true; //TODO
}




