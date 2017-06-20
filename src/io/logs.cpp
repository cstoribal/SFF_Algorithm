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
    logs = "starting logs file \n  \n";
}
MyLog::~MyLog(){}

bool MyLog::set_param(string file){
    this->filepath = file;
    return true;
}

bool MyLog::a(string addstring){
    logs += addstring;
    cout << addstring << endl;
    return true;
}
