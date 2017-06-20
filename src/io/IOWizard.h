/****************************
*	ProjSFF
*	IO_Wizard.h
*	cstoribal
*	06-04-17
****************************/

/*
Aimed at managing all the IO data tranfers
loading saving images, depths, and visualisations
*/

//updates updvl0.1

#ifndef IOWIZARD_H_
#define IOWIZARD_H_


#include <iostream>
#include <string>
#include <sstream>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <mgl2/mgl.h>
#include <mgl2/qt.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fstream>
#include <vector>

#include "../misc/miscdef.h"
#include "../io/mgldrw.h"



using namespace cv;
using namespace std;


class IOWizard{


public:
    IOWizard();
    ~IOWizard();

    // Parsing
    bool parseArgs(int argc, char** argv);
    bool setArgs(tdf_input & clonedinput);
    bool checkArgs(void);
    bool displayHelp(void);

    bool parsefile2vect(string filename, vector<vector<string> > & filedata); //ok
    bool parsevect2struct(const vector<vector<string> > & fd, tdf_input & inp); //ok
    bool storeParameters(void);


    int help;			// check if help mode is on   //TODO checking and default values if needed
    string arg_datapath;	// dataname 
    string arg_savedatapath;


    tdf_input input;
    string    input_folder;     // to add to input paths

    vector<Point> clicked; // last clicked point on window


    //bool addbuffer(const string & text); //updvl0.1
    //bool writebuffer(void); //updvl0.1

    // Image Processing
    bool readImage(const string filepath, Mat & image);
    bool nameImage(const string filepath, const string extension, int indice, string& samplepath);
    bool buildImageSet(tdf_imgset& imageSet);
    	//when parameters are written on a txt file along with samples
    bool autosetImsetParameters(tdf_imgset & imageSet);

    bool loadGroundTruth(Mat & gtmat, string filepath);



    bool showImage(const string param, const Mat & image, int timer);
    bool clickImage(const string param, const Mat & image, int timer,
			 bool (*fptr)(void*,Point), void* context) ;
    bool writeImage(const string filename, const Mat & image);

    bool write3DImage(const string filename, const Mat & image);
    bool show3DImage(const string filename, const Mat & image);
    bool mat3Dplot(void);


    bool mkdir(const string directory);
    bool set_auto_directory(string foldername);

    


private:
    string autofolder;          // automatically added to filenames
    // string program_logs;  //Upd VL0.1



    MyDrawer2D myDrawer2D;

    /////////////// Loading variables...
    string arg_filename;	// filename - used for help choice
    string arg_separator;       // default "-"
    string arg_extension;	// extension
    string outputfolder;        //
    string arg_groundtruthfile; // groundtruth filename
    int arg_groundtruth;        // is a groundtruth provided ?
    int arg_gta,arg_gtb;
    int arg_firstindice;	//
    int arg_deltaindice;        // sous-échantillonnage
    int arg_lastindice;		//
    int arg_sizeindice;		//
    double arg_scale;		// scales the image
    int arg_gauss;		// windows size
    double arg_noise_a,arg_noise_b,arg_noise_ca,arg_noise_cs;//
    int arg_sharpness;		// set focus measurement type
    int arg_depth;		// set local depth estimator type
    int arg_nrjdata;		// set energy for data related term
    int arg_nrjreg;		// set regularisation functional term
    int arg_opti;		// set optimisation algorithm
    ////////////////

};

void CallBackFunc(int event, int x, int y, int flags, void* userdata);
void CallBackFunc2(int event, int x, int y, int flags, void* userdata);

void CallBack3DQT(mglGraph *gr);
bool forwarder2(void* context);

#endif // IOWIZARD_H_
