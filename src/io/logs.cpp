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
    logs = "starting log file \n  \n";
}
MyLog::~MyLog(){}

bool MyLog::set_param(const string & file, bool verb){
    this->verbose = verb;
    this->filepath = file;
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

bool MyLog::write(void){
    ofstream out(this->filepath.c_str());
    out << this->logs;
    out.close();
    return true;
}


