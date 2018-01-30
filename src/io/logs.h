/****************************
*	logs.h
*	cstoribal
*	15-06-17
****************************/

/*
Outputs at the same time a precise log of each simulation encountered
and the results under the format of a csv
*/


#ifndef LOGS_H_ 
#define LOGS_H_

#include <string>
#include "../misc/miscdef.h"

using namespace std;

class MyLogOut{
public:
    MyLogOut(); ~MyLogOut();
    
    bool setup(tdf_input & input);
    bool write(void);
    bool write(std::string filename, bool verbose);
    bool write_deltaRMSE(std::string filename, vector<vector<vector<std::string> > > & vvv_deltarmse, vector<vector<std::string> > & vv_type12);
    string logversion;
    tdf_log output;
    
private:
    string strdata;
    bool create_new_logfile_header(std::ofstream & outfile);
    bool Format_txt(void);
    bool priv_write_deltaRMSE(std::string filename, vector<vector<vector<std::string> > > & vvv_deltarmse, vector<vector<std::string> > & vv_type12);
    
};


class MyLog{
public:
    MyLog(); ~MyLog();
    bool set_param(const string & filepath, bool verb); //setting path and filename for output;
    MyLogOut* log_data_out;

    bool a(const string & addstring); //add
    bool av(const string & addstring); //addverbose
    bool as(const string & addstring); //addsilent
    bool write(void);

    double time_r(void); 	// returns timer and reinit timer;
    double time_r(int rank); // returns timer and reinit timer;
    double time(void);   	// calls time()
    double time(int rank);   // calls time()
    void   time_i(void); // calls init timer.
    bool   clear_log(void);
    bool   clear_iteration_times(void); //
    bool   set_state(fType lr, fType ld, int iter);
    bool   set_state(std::string opti, fType lambda, int iter);
    bool   set_state(std::string opti, fType lambda);
    bool   set_eval(fType rmse, fType psnr);
    bool   set_eval_at(fType rmse, int iter);
    bool   set_eval_at(fType rmse, fType psnr, int iter);
    bool   set_bestplans(std::vector<std::string> types);
    bool   write_deltaRMSE(std::string filename, vector<vector<vector<std::string> > > &  vvv_deltarmse, vector<vector<std::string> > & vv_type12);
    
private:
    CTimer timer;
    string logs;
    string filepath; //file name+path, ready
    string logversion;
    bool verbose;

    bool print_eftypes(void);
    
};

#endif // LOGS_H_
