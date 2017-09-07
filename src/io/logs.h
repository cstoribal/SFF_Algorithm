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
    string logversion;
    tdf_log output;
    
private:
    string strdata;
    
    bool Format_txt(void);
    
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
    bool   set_state(fType lr, fType ld, int iter);
    bool   set_eval(fType rmse, fType psnr);
    
private:
    CTimer timer;
    string logs;
    string filepath; //file name+path, ready
    string logversion;
    bool verbose;

    bool print_eftypes(void);
    
};

#endif // LOGS_H_
