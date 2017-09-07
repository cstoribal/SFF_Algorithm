/**************************
*	ProjSFF
*	mainSFF.cpp
*	cstoribal
*	07-04-17
**************************/

#include <iostream>
#include <string>

#include <ctime> 
#include <stdio.h>
//#include <mgl2/qt.h>
#include <mgl2/mgl.h>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>

#include "../misc/miscdef.h"


using namespace cv;
using namespace std;


int main( int argc, char** argv )
{
   string pathmodel;
   string pathreplace;
   string output;
   int startreplace = 5;
   for (int i = 1; i < argc; i++) { 

        if (std::string(argv[i]) == "-M") {
        // The next argument should be present, model
            assert( i+1 < argc);
            pathmodel = string(argv[i+1]);
            
        }

        if (std::string(argv[i]) == "-f") {
        // The next argument should be present, filename
            assert( i+1 < argc);
            pathreplace = string(argv[i+1]);
        } 
    }
    
    vector<string> modeldata;
    vector<string> replacedata;
    vector<string> finaldata;
    fstream filemodel;
    fstream filereplace;

    string windows_eol_remove = "\r";
    string line;
    int lineidx=0;
    int start, end, weolrmv;

    filemodel.open(pathmodel.c_str());

    if (!filemodel.good())
    {
        COUT2("error opening file",pathmodel);
        return false;
    }

    // read each line of the file
    while (!filemodel.eof())
    {
        int tmplinestop =0 ;
        getline(filemodel,line);
        start = 0U;
        weolrmv = line.find(windows_eol_remove);
        if(weolrmv != string::npos){
            line = line.substr(start,weolrmv);
        }
        modeldata.push_back(line.c_str());
        lineidx+=1;
    }
    filemodel.close();


////////////////
    //lineidx = 0;
    filereplace.open(pathreplace.c_str());

    if (!filereplace.good())
    {
        COUT2("error opening file",pathreplace);
        return false;
    }

    // read each line of the file
    while (!filereplace.eof())
    {
        int tmplinestop =0 ;
        getline(filereplace,line);
        start = 0U;
        weolrmv = line.find(windows_eol_remove);
        if(weolrmv != string::npos){
            line = line.substr(start,weolrmv);
        }
        replacedata.push_back(line.c_str());
        //lineidx+=1;
    }
    filereplace.close();
    
    for(int j=0; j<lineidx; j++)
    {
        if(j<startreplace){
            finaldata.push_back(replacedata[j]);
            }
        else{
            finaldata.push_back(modeldata[j]);
            }
    }
 
    output = "";
    for(int j=0;j<lineidx; j++)
    {
        output += finaldata[j] + "\n";
    }
    ofstream out(pathreplace.c_str());
    out << output;
    out.close();
    
    
    
    
    



    return 0;

}




