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


class MyLog{
public:
    MyLog(); ~MyLog();
    bool set_param(const string & filepath, bool verb); //setting path and filename for output;
    bool a(const string & addstring); //add
    bool av(const string & addstring); //addverbose
    bool as(const string & addstring); //addsilent
    bool write(void);
    
private:
    string logs;
    string filepath; //file name+path, ready
    bool verbose;




    bool print_eftypes(void);
    
    
};

#endif // LOGS_H_
