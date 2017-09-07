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
    logversion = "v0.1";
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
    this->log_data_out->output.time[rank] = tmp;
    return tmp;
}

double MyLog::time(void){
    return timer.Time();
}

double MyLog::time(int rank){
    double tmp = timer.Time();
    this->log_data_out->output.time[rank] = tmp;
    return tmp;
}

void MyLog::time_i(void){
    return timer.Init();
}

bool MyLog::set_state(fType lr, fType ld, int iter){
    this->log_data_out->output.lambda_r = lr;
    this->log_data_out->output.lambda_d = ld;
    this->log_data_out->output.iterationlvl = iter;
    return true;
}

bool MyLog::set_eval(fType rmse, fType psnr){
    this->log_data_out->output.rmse = rmse;
    this->log_data_out->output.psnr = psnr;
    return true;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

MyLogOut::MyLogOut(){};MyLogOut::~MyLogOut(){};

bool MyLogOut::setup(tdf_input & input){
    this->output.settings = &input;
    this->output.time.resize(10);
    return true;
}

bool MyLogOut::Format_txt(void){
    tdf_log & o = output;
    tdf_input* & in = o.settings;
    int tmp_nimg = (int)floor((double)(in->file1_lasti-in->file1_firsti+1) / (double)in->file1_deltai);

    strdata = logversion + ";";
    strdata += 
	  (in->outputf_set?to_string2(in->outputfolder)	:"") + ";"
        + (in->file1_set?  to_string2(in->file1_path)	:"") + ";" 
	+ (in->file1_set?  to_string2(in->file1_ext)	:"") + ";" 
	+ (in->file1_set?  to_string2(tmp_nimg) 	:"") + ";"
	+ (in->file2_set?  to_string2(in->file2[0])	:"") + ";"
	+ (in->file2_set?  to_string2(in->file2.size())	:"") + ";"
	+ (in->groundt_set?to_string2(in->gtpath)	:"") + ";"
	+ (in->groundt_set?to_string2(in->gtb)  	:"") + ";"
	+ (in->groundt_set?to_string2(in->gta+in->gtb) 	:"") + ";"
	+ (in->focus_set?  to_string2(in->focus[0])	:"") + ";"
	+ (in->focus_set?to_string2(in->focus[in->focus.size()-1]):"") + ";"
	+ (in->preproc_set?to_string2(in->scale)	:"") + ";"
 	+ (in->preproc_set?to_string2(in->gauss)	:"") + ";"
 	+ (in->preproc_set?to_string2(in->noise_a)	:"") + ";"
 	+ (in->preproc_set?to_string2(in->noise_b)	:"") + ";"
 	+ (in->preproc_set?to_string2(in->noise_ca)	:"") + ";"
 	+ (in->preproc_set?to_string2(in->noise_cs)	:"") + ";"
 	+ (in->sharp_set?  to_string2(in->sharp)	:"") + ";"
	+ (in->depth_set?  to_string2(in->depth)	:"") + ";"
	+ (in->nrj_set?    to_string2(in->nrj_d)	:"") + ";"
	+ (in->nrj_set?    to_string2(in->nrj_r)	:"") + ";"
	+ (in->opti_set?   to_string2(in->opti) 	:"") + ";"
	+ (in->opti_set?   to_string2(in->connexity)	:"") + ";";
    
    strdata +=
	  to_string2(o.lambda_d) + ";"
	+ to_string2(o.lambda_r) + ";"
	+ to_string2(o.iterationlvl) + ";"
	+ to_string2(o.rmse) + ";"
	+ to_string2(o.psnr);

    for(int k=0;k<o.time.size();k++)
        strdata += ";" + to_string2(o.time[k]);

    
    strdata += "\n";
    return true;
}

bool MyLogOut::write(void){
    this->Format_txt();
    std::ofstream outfile;
    outfile.open("1data/logout.csv", std::ios_base::app);
    outfile << strdata;
    outfile.close();
}




