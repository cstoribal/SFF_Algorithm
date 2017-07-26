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
    logs = "Starting log file \n v0.0 \n";
}
MyLog::~MyLog(){}

bool MyLog::set_param(const string & file, bool verb){
    this->verbose = false;
    this->filepath = file;
    cout<<filepath<<endl;
    this->print_eftypes();
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
    return true;
}


