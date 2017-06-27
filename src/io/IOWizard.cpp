/****************************
*	IO_Wizard.cpp
*	cstoribal
*	06-04-17
****************************/

/*
Aimed at managing all the IO data tranfers
loading saving images, depths, and visualisations
*/



#include "IOWizard.h"



using namespace std;
using namespace cv;



IOWizard::IOWizard() {

    ///// Parsing arguments - default - unset
    this->arg_datapath = "";
    this->arg_savedatapath = "data.txt";
    this->arg_filename = "";
    this->help=0;
    this->arg_firstindice=-1;
    this->arg_deltaindice=1;
    this->autofolder="";
    this->outputfolder="";
    this->arg_sizeindice=4;
    this->arg_groundtruth = 0;
    this->arg_groundtruthfile = "";
    this->arg_gta = 1;
    this->arg_gtb = 0;
    this->arg_scale = -1;
    this->arg_gauss = -1;
    this->arg_noise_a = -1;
    this->arg_noise_b = -1;
    this->arg_noise_ca = -1;
    this->arg_noise_cs = -1;
    this->arg_sharpness = -1;
    this->arg_depth = -1;
    this->arg_nrjdata = -1;
    this->arg_nrjreg = -1;
    this->arg_opti = -1;

    this->clicked = vector<Point>() ;
    ///// tdf_input default values
    this->input =(tdf_input) {
		.file1_set=0,.file1_path="",.file1_sizeind=0,
			.file1_sep="",.file1_ext="",.file1_firsti=-1,
			.file1_deltai=1,.file1_lasti=0,
		.file2_set=0,.file2=vector<string>(),
		.groundt_set=0,.gtpath="",.gta=0,.gtb=0,
		.outputf_set=0,.outputfolder="./",
		.focus_set=0,.focus=vector<fType>(),
                .preproc_set=0,.scale=1,.gauss=1,
                        .noise_a=0.05,.noise_b=0.05,.noise_ca=0.05,.noise_cs=0.05,
		.sharp_set=0,.sharp="",
                .depth_set=0,.depth="",
		.nrj_set=0,.nrj_d="",.nrj_r="",
		.opti_set=0,.opti="",
                .lambda_r_set=0,.vect_lambda_r={1,2,3},
                .lambda_d_set=0,.vect_lambda_d={3,2,1.5}};

    ///// IOWizard behaviorial set/unset parameters
    this->img_dout_scaleselect=0;

}

IOWizard::~IOWizard(){

}

bool IOWizard::setlogs(MyLog* mylog){
    this->myLog = mylog;
    return true;
}
bool IOWizard::setlogsoutput(void){
    // to be called when parsing arguments is done. currently called at parseArgs, seems okay
    this->myLog->set_param(this->input.outputfolder + "/logs.txt",true);
    return true;
}
/////////////////////////////////////////////////////////////////
// Parsing
/////////////////////////////////////////////////////////////////


bool IOWizard::parseArgs(int argc, char** argv) {

    if (argc == 2 or argc == 3){if(std::string(argv[1]) == "help"){
            cout<<"......................................"<<endl;
            this->help=1;
            if(argc==3) this->arg_filename = string(argv[2]);
            return true;
        }
    }

    for (int i = 1; i < argc; i++) { 

        if (std::string(argv[i]) == "-D") {
        // The next argument should be present, datapath .txt
            assert( i+1 < argc);
            this->arg_datapath = string(argv[i+1]);
            
        }

        if (std::string(argv[i]) == "-f") {
        // The next argument should be present, and it should be filename
            assert( i+1 < argc);
            this->arg_filename = string(argv[i+1]);
        }

        if (std::string(argv[i]) == "-fs") {
        // The next argument should be present, and it should be file separator (if needed)
            assert( i+1 < argc);
            this->arg_filename = string(argv[i+1]);
        }

        if (std::string(argv[i]) == "-fsi") {
        // The next argument should be present, and it should be the number of figures for filename
            assert( i+1 < argc);
            this->arg_sizeindice = atoi(argv[i+1]);
        }

        if (std::string(argv[i]) == "-fdi") {
        // The next argument should be present, and it should be step between two indices
            assert( i+1 < argc);
            this->arg_deltaindice = atoi(argv[i+1]);
        }

        if (std::string(argv[i]) == "-ffi") {
        // The next argument should be present, and it should be the first indice of the sequence
            assert( i+1 < argc);
            this->arg_firstindice = atoi(argv[i+1]);
        }

        if (std::string(argv[i]) == "-fli") {
        // The next argument should be present, and it should be the lst indice of the sequence
            assert( i+1 < argc);
            this->arg_lastindice = atoi(argv[i+1]);
        }

        if (std::string(argv[i]) == "-gt") {
        // The next argument should be present, and it should be the lst indice of the sequence
            assert( i+1 < argc);
            this->arg_groundtruth = 1;
            this->arg_groundtruthfile = string(argv[i+1]);
        }

        if (std::string(argv[i]) == "-ps") {
        // The next argument should be present, and it should be sharpness
            assert( i+1 < argc);
            this->arg_sharpness = atoi(argv[i+1]);
        }
        
        if (std::string(argv[i]) == "-pd") {
        // The next argument should be present, and it should be depth measure type
            assert( i+1 < argc);
            this->arg_depth = atoi(argv[i+1]);
        }
        
        if (std::string(argv[i]) == "-pnd") {
        // The next argument should be present, and it should be data energy term
            assert( i+1 < argc);
            this->arg_nrjdata = atoi(argv[i+1]);
        }
        
        if (std::string(argv[i]) == "-pnr") {
        // The next argument should be present, and it should be regularisation energy term
            assert( i+1 < argc);
            this->arg_nrjreg = atoi(argv[i+1]);
        }
        
        if (std::string(argv[i]) == "-po") {
        // The next argument should be present, and it should be optimisation type
            assert( i+1 < argc);
            this->arg_opti = atoi(argv[i+1]);
        }

        if (std::string(argv[i]) == "-o") {
        // The next argument should be present, and it should be a valid output folder
            assert( i+1 < argc);
            this->outputfolder = string(argv[i+1]);
            cout << this->outputfolder << endl;
        }        
    }

    return true;
}


bool IOWizard::setArgs(tdf_input & clonedinput){
    // Evaluate the arguments litterally (for the user)
    // and now compares the C++ input modified data to the file data
    // passes the argument thanks to line case -1: break;
    if(arg_datapath !="")
        input_folder = arg_datapath.substr(0,arg_datapath.find_last_of("/"));
    COUT2("input folder is",input_folder);
    if(arg_filename !=""){
        input_folder = arg_filename.substr(0,arg_filename.find_last_of("/"));
        COUT2("input folder set to",input_folder);
        if(input.file1_set==1)COUT("replacing file1 input");
        input.file1_set=1;
        input.file2_set=0;
        input.file1_path    = arg_filename.substr(1+arg_filename.find_last_of("/"));
        input.file1_sizeind = arg_sizeindice;
        input.file1_sep     = arg_separator;
        input.file1_ext     = arg_extension;
        input.file1_firsti  = arg_firstindice;
        input.file1_deltai  = arg_deltaindice;
        input.file1_lasti   = arg_lastindice;
    }


    if(arg_groundtruth)
    {
        if(input.groundt_set==1)COUT("replacing groundtruth");
        input.groundt_set=1;
        input.gtpath= arg_groundtruthfile;
        input.gta   = arg_gta;
        input.gtb   = arg_gtb;
    }
    
    if(outputfolder!="")
    {
        if(input.outputf_set==1)COUT("replacing outputfolder");
        input.outputf_set=1;
        input.outputfolder= outputfolder;
    }

    if((arg_scale!=-1 || arg_gauss!=-1 || arg_noise_a!=-1 || arg_noise_b!=-1 || arg_noise_ca !=-1 || arg_noise_cs !=-1) && input.preproc_set==1 )
    {
        COUT("replacing preprocessing");
    }
    if(arg_scale!=-1) input.scale = arg_scale;
    if(arg_gauss!=-1) input.gauss = arg_gauss;
    if(arg_noise_a!=-1) input.noise_a = arg_noise_a;
    if(arg_noise_b!=-1) input.noise_b = arg_noise_b;
    if(arg_noise_ca!=-1) input.noise_ca = arg_noise_ca;
    if(arg_noise_cs!=-1) input.noise_cs = arg_noise_cs;
    

    
    if(this->arg_sharpness!=-1 && input.sharp_set==1)
    {
        COUT("replacing sharp");
    }
    switch(this->arg_sharpness){
        case 1:
            input.sharp="SMLAP";
            break;
        case 0:
            input.sharp="";
            break;
        case -1:break;
        default:
            input.sharp="NONE";
            break;
    }
    
    if(this->arg_depth!=-1 && input.depth_set ==1)
    {
        COUT("replacing depth");
        input.depth_set=1;
    }
    switch(this->arg_depth){
        case 1:
            input.depth="polynome";
            break;
        case 0:
            input.depth="";
            break;
        case -1: break;
        default:
            input.depth="None";
            break;
    }
    
    if((this->arg_nrjdata!=-1 || this->arg_nrjreg!=-1) && input.nrj_set ==1)
    {
        COUT("replacing energies");
    }
    switch(this->arg_nrjdata){
        case 1:
            input.nrj_d="";
            break;
        case 0:
            input.nrj_d="";
            break;
        case -1: break;
        default:
            input.nrj_d="None";
            break;
    }

    switch(this->arg_nrjreg){
        case 1:
            input.nrj_r="absdiff";
            break;
        case 0:
            input.nrj_r="";
            break;
        case -1: break;
        default:
            input.nrj_r="None";
            break;
    }

    if(this->arg_opti!=-1 && input.opti_set ==1)
    {
        COUT("replacing optimisation");
        input.opti_set=1;
    }
    switch(this->arg_opti){
        case 1:
            input.opti="";
            break;
        case 0:
            input.opti="";
            break;
        case -1: break;
        default:
            input.opti="None";
            break;
    }
    
    clonedinput = input;
    this->setlogsoutput();
    return true;
}

bool IOWizard::checkArgs(void){  // TODO more checks like focus.size == images.size... but after check1 & check2 ?
    if(this->help==1)
    {
        return false;
    }
    
    // FILENAME
    // Check argument
    if(this->arg_filename == "" && this->arg_datapath == "")return false;
    if(this->arg_filename != ""){
        try{
            arg_extension = arg_filename.substr(
			arg_filename.find_last_of("."));
            arg_filename = arg_filename.substr(
			0,arg_filename.find_last_of("."));
        }
        catch(int i)
        {
            cout << "Filename is " << arg_filename << 
			", impossible to find the extension" << endl;
            return false;
        }
    }

    if(arg_datapath!=""){
        try{
            
            vector<vector<string> > tmp;
            parsefile2vect(arg_datapath,tmp);
            parsevect2struct(tmp,input);
        }
        catch(int i)
        {

            COUT("echec de lecture du dataset .txt");
            return false;
        }
    }
    // 
    return true;
}


bool IOWizard::displayHelp(void){
    int all = 0;
    if(help == 1)
    {
        cout<<this->arg_filename<<endl;
        cout<<"*********** ProjSFF - Help ***********"<<endl;
        if(this->arg_filename=="all") all=1;

        if(this->arg_filename=="" or all==1)
        { //Message générique
		cout<<"Computes the SFF algorithm on a set of images,"
		" depending on various parameters. To get a more precise"
		" help, type help [parameter] with one of the parameters"
		" amongst the follonwing ones, or help all." << endl;
		cout<<"Syntax is ProjSFF {[param] [arg] }"<<endl;
		COUT("Parameters: -D -S -f -fsi -fs -ffi -fdi -fli");
                COUT("            -gt -sc -gw -ps -pnd -pnr -po -o");
		COUT("**************************************");
        }
        if(this->arg_filename=="-D" or all==1)
        {
        COUT("[-D] parameter sets the name of the .txt data input.");
        COUT("     the format of this .txt file is the following one:");
        COUT("     1/0; -f; -fsi; -fs; ext; -ffi; -fdi; -fli \\n");
        COUT("     1/0; -filenames(a set of image paths);;;; \\n");
        COUT("     1/0; -gt; gta; gtb; \\n");
        COUT("     1/0; -o; \\n");
        COUT("     1/0; -focusdepths(a vector of imagefocus); \\n");
        COUT("     1/0; -sc; -gw; -n1; -n2; -nga; -ngs\\n");
        COUT("     1/0; -psharp; \\n");
        COUT("     1/0; -pdepth; \\n");
        COUT("     1/0; -pnd energy data; -pnr energyregul; \\n");
        COUT("     1/0; -poptimization; \\n");
        COUT("     1/0; -lambda_r vector; \\n");
        COUT("     1/0; -lambda_d vector; \\n");
        COUT("     note : called parameters overwrite those ones");
        COUT("     nb : data.txt should lie near images.png");
        }
        if(this->arg_filename=="-D" or all==1)
        {
        COUT("[-S] argument : path of the .txt data output");
        COUT("     will store all commands & parameters on a .txt file");
        }
        if(all ==1)
        {
		cout<<endl;
		cout<<"----------------------fileprocessing---------"<<endl;
        } 
        if(this->arg_filename=="-f" or all==1)
        {
	cout<<"[-f] parameter, required, is the filename for input data with "<<endl;
	cout<<"       extension but no number. For instance, sample-0003.jpg "<<endl;
	cout<<"       will be called sample.jpg. When working on a set, arguments "<<endl;
	cout<<"       ffi and fli will determine the sample autonaming." <<endl;
	cout<<"       Otherwise, use -F flag (not implemented yet) for specific sets"<<endl;
        }
        if(this->arg_filename=="-fs" or all==1)
        {
		cout<<"[-fs] parameter sets the separator type between" << endl;
		cout<<"      name & number. Warning, this cannot be another 'parameter' " <<endl;
		cout<<"      like -f for instance. set it to None for no separator"<<endl;
        }
        if(this->arg_filename=="-fsi" or all==1)
        {
		cout<<"[-fsi] parameter sets size of the index of the file (int)" << endl;
		cout<<"      for instance '3' sets filenameseparator001.ext " <<endl;
        }
        if(this->arg_filename=="-fdi" or all==1)
        {
		cout<<"[-fdi] specifies the delta between two indices"<<endl;
        }
        if(this->arg_filename=="-ffi" or all==1)
        {
		cout<<"[-ffi] parameter is first sample indice"<<endl;
        }
        if(this->arg_filename=="-fli" or all==1)
        {
		cout<<"[-fli] parameter is last sample indice"<<endl;
        }
        if(this->arg_filename=="-gt" or all==1)
        {
		cout<<"[-gt] parameter is groundtruth image file"<<endl;
        }
        if(this->arg_filename=="-o" or all==1)
        {
		cout<<"[-o] parameter is output folder. Better try './output' for instance "<<endl;
        }
        if(all==1)
        {
		cout<<endl;
		cout<<"---------------------preprocessing---------"<<endl;
        }
        if(this->arg_filename=="-sc" or all==1)
        {
		COUT("[-sc] argument: scale of img");
        }
        if(this->arg_filename=="-gw" or all==1)
        {
		COUT("[-gw] argument: gausswindow");
        }
        if(all==1)
        {
		cout<<endl;
		cout<<"----------------------sffprocessing---------"<<endl;
        }
        if(this->arg_filename=="-ps" or all==1)
        {
		cout<<"[-ps] parameter, integer, is the sharpness operator type (int)"<<endl;
		cout<<"      x 1 is SMLAP"<<endl; 
        }
        if(this->arg_filename=="-pnd" or all==1)
        {
		cout<<"[-pnd] parameter is depth energy type (int)"<<endl;
        }
        if(this->arg_filename=="pnr" or all==1)
        {
		cout<<"[-pnr] parameter is regularisation energy type (int)"<<endl;
        }
        if(this->arg_filename=="-po" or all==1)
        {
		cout<<"[-po] parameter is optimisation type (int)"<<endl;
        }

        if(this->arg_filename=="NONE" or all==1)
        {
		cout<<"[] parameter is "<<endl;
        }
        return true;
    }
    else
    {
        cout<<"Failure, somewhere, likely on the argument parsing"<<endl;
        cout<<"Try 'ProjSFF help all'"<<endl;
    }
    return true;
}


bool IOWizard::parsefile2vect(string filename, vector<vector<string> > & filedata){
    // Create a file reading object
    // filedata.clear();
    vector<vector<string> >().swap(filedata);
    string line;
    string word;
    fstream file;
    string delimiter = ";";
    string windows_eol_remove = "\r";
    string rmv_spaces = " ";
    int lineidx=0;
    int wrdidx=0;
    int start, end, weolrmv;
    int rmv;

    file.open(filename.c_str());
    if (!file.good())
    {
        COUT2("error opening file",filename);
        return false;
    }

    // read each line of the file
    while (!file.eof())
    {
        int tmplinestop =0 ;
        filedata.push_back(vector<string>() );
        getline(file,line);
        wrdidx=0;
        start = 0U;
        end = line.find(delimiter, start);
        weolrmv = line.find(windows_eol_remove);
        if(weolrmv != string::npos){
            line = line.substr(start,weolrmv);
        }
        while(end != string::npos ){// read each ";" word until eof line
            rmv = line.substr(start,end-start).find_first_not_of(rmv_spaces);
            if(rmv == string::npos){ // nothing in word except spaces
                filedata[lineidx].push_back("");
                start = end+delimiter.length();
                end = line.find(delimiter, start);
                
            } else {   // crop word to content 
                start = start +rmv;
                rmv = line.substr(start,end-start).find_last_not_of(rmv_spaces);
                filedata[lineidx].push_back(line.substr(start,rmv + 1));
                start = end+delimiter.length();
                end = line.find(delimiter, start);
                
            }
            wrdidx+=1;
        }
        rmv = line.substr(start).find_first_not_of(rmv_spaces);
        if(rmv == string::npos)
            filedata[lineidx].push_back("");
        else {
            start = start +rmv;
            rmv = line.substr(start).find_last_not_of(rmv_spaces);
            filedata[lineidx].push_back(line.substr(start,rmv+1));
        }
        lineidx+=1;
    }
    return true;
    //filedata op en théorie
}

bool IOWizard::parsevect2struct(const vector<vector<string> > & fd, tdf_input & inp){
    // Follows the exact structure of tdf_input
    // fd est un vecteur<lignes de vecteurs<mots> > 
    int l_idx =0;
    int w_idx =0;
    while( l_idx<fd.size() )
    {
    w_idx=0;
    while( w_idx<fd[l_idx].size() )
    {
        switch(l_idx){
        case 0:
            switch(w_idx){
            case 0:
                inp.file1_set=atoi(fd[l_idx][w_idx].c_str());
                if(inp.file1_set!=1)w_idx = fd[l_idx].size();
                break;
            case 1: inp.file1_path = fd[l_idx][w_idx];break;
            case 2: inp.file1_sizeind = atoi(fd[l_idx][w_idx].c_str());break;
            case 3: inp.file1_sep = fd[l_idx][w_idx];break;
            case 4: inp.file1_ext  = fd[l_idx][w_idx].c_str();break;
            case 5: inp.file1_firsti  = atoi(fd[l_idx][w_idx].c_str());break;
            case 6: inp.file1_deltai  = atoi(fd[l_idx][w_idx].c_str());break;
            case 7: inp.file1_lasti   = atoi(fd[l_idx][w_idx].c_str());break;
            default:
                break;
            }
            break;
        case 1:
            switch(w_idx){
            case 0:
                inp.file2_set=atoi(fd[l_idx][w_idx].c_str());
                if(inp.file2_set!=1)w_idx = fd[l_idx].size();
                else inp.file2.resize( fd[l_idx].size()-1 );
                break;
            default:
                inp.file2[w_idx-1]= fd[l_idx][w_idx] ;
                break;
            }
            break;
        case 2:
            switch(w_idx){
            case 0:
                inp.groundt_set=atoi(fd[l_idx][w_idx].c_str());
                if(inp.groundt_set!=1)w_idx = fd[l_idx].size();
                break;
            case 1:
                inp.gtpath =fd[l_idx][w_idx];
                break;
            case 2:
                inp.gta=atof(fd[l_idx][w_idx].c_str());
                break;
            case 3:
                inp.gtb=atof(fd[l_idx][w_idx].c_str());
                break;
            default:
                break;
            }
            break;
        case 3:
            switch(w_idx){
            case 0:
                inp.outputf_set=atoi(fd[l_idx][w_idx].c_str());
                if(inp.outputf_set!=1)w_idx = fd[l_idx].size();
                break;
            case 1:
                inp.outputfolder = fd[l_idx][w_idx];
                break;
            default:
                break;
            }
            break;
        case 4:
            switch(w_idx){
            case 0:
                inp.focus_set=atoi(fd[l_idx][w_idx].c_str());
                CPING("done");
                if(inp.focus_set!=1)w_idx = fd[l_idx].size();
                else inp.focus.resize( fd[l_idx].size()-1 );
                break;
            default:
                inp.focus[w_idx-1]= atof( fd[l_idx][w_idx].c_str() ) ;
                break;
            }
            break;
        case 5:
            switch(w_idx){
            case 0:
                inp.preproc_set=atoi(fd[l_idx][w_idx].c_str());
                if(inp.preproc_set!=1)w_idx = fd[l_idx].size();
                break;
            case 1:
                inp.scale = atof( fd[l_idx][w_idx].c_str() );
                break;
            case 2:
                inp.gauss = atoi( fd[l_idx][w_idx].c_str() );
                break;
            case 3:
                inp.noise_a = atof( fd[l_idx][w_idx].c_str() );
                break;
            case 4:
                inp.noise_b = atof( fd[l_idx][w_idx].c_str() );
                break;
            case 5:
                inp.noise_ca = atof( fd[l_idx][w_idx].c_str() );
                break;
            case 6:
                inp.noise_cs = atof( fd[l_idx][w_idx].c_str() );
                break;
            default:
                break;
            }
            break;
        case 6:
            switch(w_idx){
            case 0:
                inp.sharp_set=atoi(fd[l_idx][w_idx].c_str());
                if(inp.sharp_set!=1)w_idx = fd[l_idx].size();
                break;
            case 1:
                inp.sharp = fd[l_idx][w_idx];
                break;
            default:
                break;
            }
            break;
        case 7:
            switch(w_idx){
            case 0:
                inp.depth_set=atoi(fd[l_idx][w_idx].c_str());
                if(inp.depth_set!=1)w_idx = fd[l_idx].size();
                break;
            case 1:
                inp.depth = fd[l_idx][w_idx];
                break;
            default:
                break;
            }
            break;
        case 8:
            switch(w_idx){
            case 0:
                inp.nrj_set=atoi(fd[l_idx][w_idx].c_str());
                if(inp.nrj_set!=1)w_idx = fd[l_idx].size();
                break;
            case 1:
                inp.nrj_d = fd[l_idx][w_idx];
                break;
            case 2:
                inp.nrj_r = fd[l_idx][w_idx];
                break;
            default:
                break;
            }
            break;
        case 9:
            switch(w_idx){
            case 0:
                inp.opti_set=atoi(fd[l_idx][w_idx].c_str());
                if(inp.opti_set!=1)w_idx = fd[l_idx].size();
                break;
            case 1:
                inp.opti = fd[l_idx][w_idx];
                break;
            default:break;
            }
            break;
        case 10:
            switch(w_idx){
            case 0:
                inp.lambda_r_set=atoi(fd[l_idx][w_idx].c_str());
                if(inp.lambda_r_set!=1)w_idx = fd[l_idx].size();
                else inp.vect_lambda_r.resize( fd[l_idx].size()-1 );
                break;
            default:
                inp.vect_lambda_r[w_idx-1] = atof(fd[l_idx][w_idx].c_str());
                break;
            }
            break;
        case 11:
            switch(w_idx){
            case 0:
                inp.lambda_d_set=atoi(fd[l_idx][w_idx].c_str());
                if(inp.lambda_d_set!=1)w_idx = fd[l_idx].size();
                else inp.vect_lambda_d.resize( fd[l_idx].size()-1 );
                break;
            default:
                inp.vect_lambda_d[w_idx-1] = atof(fd[l_idx][w_idx].c_str());
                break;
            }
            break;
        default:break;
        }
        w_idx += 1;
    }
    l_idx += 1;
    }
    
    return true;
}


bool IOWizard::storeParameters(void){

    tdf_input i = this->input;
    string textdata;
    textdata = to_string2(i.file1_set) + " ; " + i.file1_path + " ; " + to_string2(i.file1_sizeind) + " ; " + i.file1_sep + " ; " + i.file1_ext + " ; " + to_string2(i.file1_firsti) + " ; " + to_string2(i.file1_deltai) + " ; " + to_string2(i.file1_lasti) + "\n" ;

    textdata += to_string2(i.file2_set);
    for(int j=0; j<i.file2.size();j++)
        textdata +=  " ; " + i.file2[j];
    textdata += "\n";

    textdata += to_string2(i.groundt_set) + " ; " + i.gtpath + " ; " + to_string2(i.gta) + " ; " + to_string2(i.gtb) + " \n";

    textdata += to_string2(i.outputf_set) + " ; " + i.outputfolder + "\n";

    textdata += to_string2(i.focus_set);
    for(int j=0; j<i.focus.size();j++)
        textdata +=  " ; " + to_string2(i.focus[j]);
    textdata += "\n";

    textdata += to_string2(i.preproc_set) + " ; " + to_string2(i.scale) + " ; " + to_string2(i.gauss) + " ; " + to_string2(i.noise_a) + " ; " + to_string2(i.noise_b) + " ; " + to_string2(i.noise_ca) + " ; " + to_string2(i.noise_cs) + "\n";

    textdata += to_string2(i.sharp_set) + " ; " + i.sharp + "\n";

    textdata += to_string2(i.depth_set) + " ; " + i.depth + "\n";

    textdata += to_string2(i.nrj_set) + " ; " + i.nrj_d+" ; "+i.nrj_r+"\n";

    textdata += to_string2(i.opti_set) +" ; "+ i.opti + "\n";

    textdata += to_string2(i.lambda_r_set);
    for(int j=0; j<i.vect_lambda_r.size();j++)
        textdata +=  " ; " + to_string2(i.vect_lambda_r[j]);
    textdata += "\n";

    textdata += to_string2(i.lambda_d_set);
    for(int j=0; j<i.vect_lambda_d.size();j++)
        textdata +=  " ; " + to_string2(i.vect_lambda_d[j]);
    textdata += "\n";

    // COUT(textdata);
    string filepath = i.outputfolder + "/" + this->arg_savedatapath;
    ofstream out(filepath.c_str());
    out << textdata;
    out.close();
    return true;
}

/*
bool IOWizard::addbuffer(const string & text){
    this->program_logs += text;    
    return true;
}

bool IOWizard::writebuffer(void){
    string filepath = input.outputfolder + "/" + "logs.txt";
    ofstream out(filepath.c_str() );
    out << this->program_logs;
    out.close();
    return true;    
}
*/ //Update updVL0.1

////////////////////////////////////////////////////////////////////////
////////////////////////////Load - Store////////////////////////////////
////////////////////////////Load - Store////////////////////////////////
////////////////////////////Load - Store////////////////////////////////
////////////////////////////////////////////////////////////////////////




bool IOWizard::readImage(const string filepath, Mat & image) {
    image = imread( filepath, IMREAD_COLOR );
    if( image.empty() )
    {
        cout << "Could not open or find the image :" << filepath << endl ;
        return false;
    }
    image.convertTo(image,CV_TF,1.0/255);
    return true;
}

bool IOWizard::nameImage(const string filepath, const string extension, int indice, string& samplepath){
    samplepath = "NONE";
    stringstream ss;
    ss << indice;
    string str= ss.str();
    
    if(str=="0") str="";
    if(indice<0) return false;
    if(indice>=pow(10,input.file1_sizeind) ) return false;
    samplepath = filepath;
    
    int k = input.file1_sizeind-1;
    while( indice<pow(10,k) )
    {
        samplepath += string("0");
        k--;
    }
    samplepath += str + extension;
    return true;
    
    return false;
}

bool IOWizard::buildImageSet(tdf_imgset& imageSet) {
    //Depending on the data builds the imageSet
    imageSet = vector<struct tagged_img>(); //vecteur
    int tmprank = 0;
    CPING("building dataset");
    if(input.file1_set==1 && input.file2_set==0)
    {
        if(input.file1_firsti ==-1){
            CPING("option1");
            struct tagged_img timg;
            timg.rank=-1;
            timg.name=input.file1_path+input.file1_ext;
            Mat imat;
            if(!readImage(input_folder+"/"+timg.name, imat)) return false;
            if(input.scale!=1)
                resize(imat, imat, Size(), input.scale, input.scale);
            split(imat,timg.ivmat);
            timg.dim = timg.ivmat.size();
            imageSet.push_back(timg);
            return true;
        }
        CPING("option2");
        for(int i=input.file1_firsti; i<input.file1_lasti+1; i+=input.file1_deltai){
            struct tagged_img timg;
            timg.rank = tmprank;
            tmprank++;
            timg.focus = (float)i;
            nameImage(input.file1_path, input.file1_ext, i, timg.name);
            // get name
            Mat imat = Mat();
            if(!readImage(input_folder+"/"+timg.name, imat)) return false;
            if(input.scale!=1)
                resize(imat, imat, Size(), input.scale, input.scale);
            // read mat
            split(imat,timg.ivmat);
            // cvt into vector dimension dim
            timg.dim = timg.ivmat.size();
            // copyTo timg.ivmat
            imageSet.push_back(timg);
        }
        return true;
    }
    if(input.file2_set==1)
    {
        CPING("option3");
        for(tmprank=0;tmprank<input.file2.size();tmprank++)
        {
            struct tagged_img timg;
            timg.rank = tmprank;
            if(input.focus_set==1)
                timg.focus = input.focus[tmprank];
            else
                timg.focus = (fType)tmprank;
            //nameImage(input.file1_path, input.file1_ext, i, timg.name);
            // get name
            Mat imat = Mat();

            
            if(!readImage(input_folder+"/"+input.file2[tmprank],imat))return false;
            if(input.scale!=1)
                resize(imat, imat, Size(), input.scale, input.scale);
            // read mat
            split(imat,timg.ivmat);
            // cvt into vector dimension dim
            timg.dim = timg.ivmat.size();
            // copyTo timg.ivmat
            imageSet.push_back(timg);
        }
        return true;   
        
    }
    return true;
}


bool IOWizard::autosetImsetParameters(tdf_imgset & imageSet) {
    for(int i=0; i<imageSet.size(); i++)
    {
        imageSet[i].dpth = 0; //TODO TEST //(fType)imageSet[i].rank;
        //imageSet[i].focus = (fType)imageSet[i].rank;
    }
    return true;
}


bool IOWizard::loadGroundTruth(Mat & gtmat,string filepath = "") {
    if(input.groundt_set == 0){
        cout << "Error : calling groundtruth while no path set" << endl;
        return false;
    }
    if(filepath == "") filepath = this->input_folder+"/"+this->input.gtpath;
    readImage(filepath, gtmat);
    vector<Mat1T> tmpmatv;
    split(gtmat,tmpmatv);
    gtmat = tmpmatv[0];
    gtmat = gtmat*(input.gta)+input.gtb;
    if(input.scale!=1)
        resize(gtmat, gtmat, Size(), input.scale, input.scale);
    // cvtColor(gtmat, gtmat, CV_BGR2GRAY); ok pour version float mfo
    
    //TODO check if 1C or 3C matrix 

    return true;
}

////SCALING////
////SCALING////
////SCALING////
bool IOWizard::img_setscale(fType min, fType max, int select){
    //select=0 -> erase all
    //select=1 -> Depthmap

    if(max == min) throw("error omax == omin in setscale");
    if(select<=0)
    {
        select = 0;
        this->img_dout_scaleset.resize(select,false);
        this->img_dout_sca.resize(select,1);
        this->img_dout_scb.resize(select,0);
        myLog->a("IOWizard: unset all scales");
        return true;
    }
        
    if(this->img_dout_sca.size()<select){
        this->img_dout_scaleset.resize(select,false);
        this->img_dout_sca.resize(select,1);
        this->img_dout_scb.resize(select,0);
    }
    this->img_dout_scaleset[select-1] = true;
    this->img_dout_sca[select-1] = 1/(max - min);
    this->img_dout_scb[select-1] = min;
    this->img_dout_scaleselect = select;
    return true;
}
bool IOWizard::img_setscale(int select){// 0 to unset.
    if(this->img_dout_scaleset.size() < select)
    {
        myLog->a("Failure in iowizard output scale selector too big");
        this->img_dout_scaleselect = 0;
        return false;
    }
    if(!this->img_dout_scaleset[select-1])
    {
        myLog->a("Falure in ioWizard output scale selected not set");
        this->img_dout_scaleselect = 0;
        return false;
    }
    this->img_dout_scaleselect = select;
    return true;
}

bool IOWizard::img_unsetscale(void){
    this->img_dout_scaleselect = 0;
    return true;
}
//// IMG DISPLAY & SAVE ////
//// IMG DISPLAY & SAVE ////
//// IMG DISPLAY & SAVE ////

bool IOWizard::showImage(const string param, const Mat & image, int timer) {
    if(0){return false;}
    fType imin, imax, scale;
    Mat imat;
    image.copyTo(imat);
    if(img_dout_scaleselect)
    {
        imat = (imat - this->img_dout_scb[img_dout_scaleselect-1])*this->img_dout_sca[img_dout_scaleselect-1];
    }
    else
    {
        cv::minMaxLoc(imat, &imin, &imax);
        scale = 1;
        if(imin != imax) scale = 1/(imax-imin);
        imat = (imat-imin)*scale;
    }
    imat.convertTo(imat,CV_8UC1,255);
    //normalize(imat, imat, 0, 255, NORM_MINMAX);
    if(param == "scale") namedWindow( "View", WINDOW_NORMAL );
    else namedWindow( "View", WINDOW_AUTOSIZE);
    
    setMouseCallback("View", CallBackFunc, &(this->clicked) );

    imshow( "View", imat );
    resizeWindow( "View" , 900,600);
    waitKey(timer);
    
    return true;
}

bool IOWizard::writeImage(const string filename, const Mat & image){
    if(0){return false;}
    fType imin, imax, scale;
    fType omin, omax;
    Mat imat;
    image.copyTo(imat);
    if(img_dout_scaleselect)
    {
        imat = (imat - this->img_dout_scb[img_dout_scaleselect-1])*this->img_dout_sca[img_dout_scaleselect-1];
    }
    else
    {
        cv::minMaxLoc(imat, &imin, &imax);
        scale = 1;
        if(imin != imax) scale = 1/(imax-imin);
        imat = (imat-imin)*scale;
    }
    imat.convertTo(imat, CV_8UC1, 255);
    //normalize(image, imat, 0, 255, NORM_MINMAX);
    imwrite( autofolder+filename, imat );

    return true;
}


bool IOWizard::write3DImage(const string filename, const Mat & image){
    
    //TODO perpare vector evolution (vector mat &)

    fType imin, imax;
    int resolution;
    int rows = image.rows;
    int cols = image.cols;
    mglGraph gr;//(NULL,filename);
    //Mat imat;
    //image copyTo(imat);
    mglData a;
    if(img_dout_scaleselect)
    {
        imin = this->img_dout_scb[img_dout_scaleselect-1];
        imax = 1/this->img_dout_sca[img_dout_scaleselect-1]+imin;
    }
    else
    {
        cv::minMaxLoc(image,&imin,&imax);
        
    }
    // Data computation for sketching
    a.Create(image.rows,image.cols);
    for(int i=0; i<rows; i++)    for(int j=0; j<cols; j++)
    {
        a.a[i+rows*j] = image.at<fType>(i,j);
    }
    // data import
    gr.Title(filename.c_str());
    gr.Light(true);  gr.Rotate(50,60);
    gr.SetRange('z',imax,imin); //gr.SetOrigin(-1,-1,-1); 
    gr.Axis();
    gr.Label('z',"depth",0); 
    gr.Box();  gr.Surf(a);//"kw|:");

    gr.WritePNG((autofolder + filename).c_str()); 
    

    return true;
    
    
}

bool IOWizard::show3DImage(const string filename, const Mat & image){
    //TODO prepare vector evolution
    fType imin, imax;
    int resolution;
    int cols = image.cols;
    vector<Mat> vmat_data;
    vmat_data.push_back(image);

    if(img_dout_scaleselect)
    {
        imin = this->img_dout_scb[img_dout_scaleselect-1];
        imax = 1/this->img_dout_sca[img_dout_scaleselect-1]+imin;
    }
    else
    {
        cv::minMaxLoc(image,&imin,&imax);
        
    }

    this->myDrawer2D.vmat_data = vmat_data;
    this->myDrawer2D.imin = imin;
    this->myDrawer2D.imax = imax;
    this->myDrawer2D.rows = image.rows;
    this->myDrawer2D.cols = image.cols;
    this->myDrawer2D.filename = filename;



    mglQT gr(&(this->myDrawer2D),filename.c_str());

    //gr.WritePNG((autofolder + filename).c_str());
    gr.Run();    
    
    

    return true;
    
    
}



bool IOWizard::mkdir(const string directory){
    string var = string("mkdir \"") + directory + string("\"");
    const char * c = var.c_str();
    system(c);
    
    return true;
}

bool IOWizard::set_auto_directory(string foldername){
    //this function impacts every loaded image by adding a reference folder.
    this->autofolder = foldername + "/";
    return true;
}






///////////////////////////////////////////////////////////
////	Live Input
///////////////////////////////////////////////////////////

void CallBackFunc(int event, int x, int y, int flags, void* userdata){

    if  ( event == EVENT_RBUTTONDOWN )
    {
        cout << "Right button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
        //vector<Point> tmpvect = *(vector<Point>*)userdata; 
        //tmpvect.push_back(Point(x,y));//push_back(Point(x,y));
        //*(vector<Point>*)userdata = tmpvect;
        (*(vector<Point>*)userdata).push_back(Point(x,y));
    }
    else if  ( event == EVENT_MBUTTONDOWN )
    {
        cout << "Middle button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
    }
    else if ( event == EVENT_MOUSEMOVE )
    {
        //cout << "Mouse move over the window - position (" << x << ", " << y << ")" << endl;
    }

}

/*
void CallBackFunc2(int event, int x, int y, int flags, void* userdata){
    // TODO debug here
    if  ( event == EVENT_RBUTTONDOWN )
    {
        cout << "Right button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
        //(*)(Point) userdata;
        (*(interpol_func*) userdata).fptr( (*(interpol_func*) userdata).context, Point(x,y));
    }
}

*/


/*
bool IOWizard::clickImage(const string param, const Mat & image, int timer, bool (*fptr)(void*,Point), void* context) {
    if(0){return false;}
    fType imin, imax, scale;
    Mat imat;
    image.copyTo(imat);
    cv::minMaxLoc(imat, &imin, &imax);
    scale = 1;
    if(imin != imax) scale = 1/(imax-imin);
    imat = (imat-imin)*scale;
    imat.convertTo(imat,CV_8UC1,255);
    normalize(imat, imat, 0, 255, NORM_MINMAX);
    if(param == "scale") namedWindow( "View", WINDOW_NORMAL );
    else namedWindow( "View", WINDOW_AUTOSIZE);
    
    interpol_func tmp;
    tmp.fptr = fptr;
    tmp.context = context;
    setMouseCallback("View", CallBackFunc2, &(tmp) );

    imshow( "View", imat );
    //resizeWindow( "View" , 900,600);
    waitKey(timer);
    
    return true;
}
*/



