/****************************
*	ProjSFF
*	optimization.cpp
*	cstoribal
*	03-05-17
****************************/

#include "optimization.h"


OptiClass::OptiClass(){ this->set = false; };
OptiClass::~OptiClass(){};

bool OptiClass::setlogs(MyLog* mylog){
    this->myLog = mylog;
    return true;
}

bool OptiClass::set_param(tdfp_opti popti){ //const void* param){
//sets name_opti, 
// calls a set_param_gcut
// will set new *gco, nb_pixels, nb_labels, etc., will malloc what's needed
// warning not to funk the system with const void* TODO
// We'll need void* to hold dmat AND the initialised instance of energy so that we can start it up. We'll also need a private function to convert matrix in vectors and reciproquely
// TODO be sure to check that we are working with the good number of pixels, rows, cols...
    this->energyClass = popti.energyclass;
    this->name_opti = popti.type;
    this->nb_pixels = popti.nb_pixels;
    this->nb_labels = popti.nb_labels;
    this->width     = popti.width;
    this->height    = popti.height;
    this->labels    = popti.labels;


    build_rank4xy();
    
    this->set = true;
}

bool OptiClass::do_optimization(void){
    if(0){
        set_optimization_gco_grid();
        return compute_gco_grid();
    }
    if(0){
        set_optimization_gco_gen();
        return compute_gco_gen();
    }
    if(1){
        return compute_opt_binary();
    }
    if(0){
        return compute_opt_multiscale();
    }
    return false;

}


/*
bool OptiClass::set_optimization(void){
// calls set_optimization_gcotype()
// also stores the vector of matrix (depth*mat) to the data_in 
// which will call some functions from param to gco to initialise the gco
    if(1) return set_optimization_gco_grid();
    if(0) return set_optimization_gco_gen();
    return false;
}

bool OptiClass::compute_optimization(void){
// basically computes optimization according to the type.
// calls compute_gco 
    if(1) return compute_gco_grid();
    if(0) return compute_gco_gen();
    return false;
}
*/

bool OptiClass::writebackmatrix(Mat1T & do_mat){
// according to optimization type, stores the regularised matrix from the vectors to the sff class
// 
    do_mat = Mat::zeros(height,width,CV_TF);
    convert_vec2mat(data_out,do_mat);
    return true;

}


///////////////////////////
//		privates //
///////////////////////////


/////////////////////////// Common


bool OptiClass::build_rank4xy(void){
    // knowing width & height, we can build the mat1i and vector<point> for conversions
    this->getrank = Mat::zeros(height,width,CV_32S); //TODO check h, w swap
    this->getxy   = vector<Point>(width*height);
    int k=0;
    for(int i=0;i<height;i++){
    for(int j=0;j<width;j++){
        getxy[k]=Point(j,i); //TODO warning 
        getrank.at<int>(i,j) = k;
        k++;        
    }}
    assert(k==this->nb_pixels);
    return true;
}

bool OptiClass::convert_mat2labvec(const vector<Mat1E> & vmat, vector<eType> & vect){
    // knowing height*width*labels we can construct the vector from the matrix vector.
    int nb_lab = vmat.size();
    vect.resize(height*width*nb_lab);
    for(int i=0; i<height;   i++){
    for(int j=0; j<width;    j++){
    for(int l=0; l<nb_lab; l++){
        vect[i*width*nb_lab+j*nb_lab+l]=vmat[l].at<eType>(i,j);        
    }}}

    return true;
}


bool OptiClass::convert_vec2mat(const vector<int> & vect, Mat1T & vmat){
    // if we know height*width (check it's nb_pixels) we can get the matrix from the vector
    // in addition, we convert it immediately to the old label format. so that's done.
    for(int i=0; i<height; i++){
    for(int j=0; j<width; j++){
        vmat.at<fType>(i,j)=labels[vect[i*width+j]];
    }}
    return true;
}

bool OptiClass::convert_vec2mat(const vector<int> & vect, Mat1i & vmat){
    // no typecast
    for(int i=0; i<height; i++){
    for(int j=0; j<width; j++){
        vmat.at<int>(i,j)=vect[i*width+j];
    }}
    return true;
}


/////////////// graph cuts
bool OptiClass::set_optimization_gco_grid(void){
// calls set_optimization_gcotype()
// also stores the vector of matrix (depth*mat) to the data_in 
// which will call some functions from param to gco to initialise the gco
    // get vector of matrix;
    Mat1E tmp_mat = Mat::zeros(height,width,CV_TE);
    vector<Mat1E> tmp_data_energy(nb_labels);
    energyClass->getDataEnergy_3DMatrix(this->labels,tmp_data_energy);

    convert_mat2labvec(tmp_data_energy,data_in);
    tmp_mat = Mat::zeros(nb_labels,nb_labels,CV_TE);
    energyClass->getCrossLabelMatrix(labels,tmp_mat);
    smoothvect.resize(nb_labels*nb_labels);
    for(int i=0; i<nb_labels; i++){
    for(int j=0; j<nb_labels; j++){
        smoothvect[i*nb_labels+j]=tmp_mat.at<eType>(i,j);
    }}
    
    try{
        gco = new GCoptimizationGridGraph(width,height,nb_labels);
        gco->setDataCost(&data_in[0]);
        gco->setSmoothCost(&smoothvect[0]);
        if(is_same<eType,int>::value)
            printf("\nBefore optimization energy is %lli",
			gco->compute_energy());
        else
            printf("\nBefore optimization energy is %f",
			gco->compute_energy());
	}
    catch (GCException e){
		e.Report();
	}



    return true;
}



bool OptiClass::set_optimization_gco_gen(void){
// first arranges the data
// calls set_optimization_gco_general()
// also stores the vector of matrix (depth*mat) to the data_in 
// which will call some functions from param to gco to initialise the gco
    // get vector of matrix;
}





bool OptiClass::compute_gco_grid(void){
    gco->expansion(10);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
    data_out.resize(nb_pixels);
    for(int i=0; i<nb_pixels; i++){
        data_out[i] = gco->whatLabel(i);
        }

    if(is_same<eType,int>::value)
        printf("\nAfter optimization energy is %lli \n",
		gco->compute_energy());
    else
        printf("\nAfter optimization energy is %f \n",
		gco->compute_energy());

    

    return true;
    
}

bool OptiClass::compute_gco_gen(void){

    gco->expansion(10);// run expansion for 2 iterations. For swap use gc->swap(num_iterations);
    data_out.resize(nb_pixels);
    for(int i=0; i<nb_pixels; i++){
        data_out[i] = gco->whatLabel(i);
        }

    if(is_same<eType,int>::value)
        printf("\nAfter optimization energy is %lli \n",
		gco->compute_energy());
    else
        printf("\nAfter optimization energy is %f \n",
		gco->compute_energy()); 

    return true;
    
}





bool OptiClass::compute_opt_binary(void){
// Parcourt tous les labels les uns après les autres.
// en déduit une matrice comparant les labels deux à deux. 
//     - M matrice des labels en construction (+1)
//     - E vectmat energie attache aux données
//     - L mat energie labels
//     - n indice d'itération ( n=0, n++, n<nblabels-1 )
//     - start, loop.
    Mat1i M = Mat::zeros(height,width,CV_32S);
    Mat1i N = Mat::zeros(height,width,CV_32S);
        
    Mat1E lmat = Mat::zeros(nb_labels,nb_labels,CV_TE);
    energyClass->getCrossLabelMatrix(labels,lmat);
    smoothvect.resize(2*2);
    
    vector<Mat1E> Edata(nb_labels);
    vector<Mat1E> Edata_slice(2);
    energyClass->getDataEnergy_3DMatrix(this->labels,Edata);
    for(int n=0; n<nb_labels-1; n++)
    {
        Edata_slice[0] = Edata[n];
        Edata_slice[1] = Edata[n+1];
        this->convert_mat2labvec(Edata_slice,data_in);

        for(int i=0; i<2; i++) for(int j=0; j<2; j++)
        {
            smoothvect[i*2+j]=lmat.at<eType>(n+i,n+j);
        }
        try{
            gco = new GCoptimizationGridGraph(width,height,2);
            gco->setDataCost(&data_in[0]);
            gco->setSmoothCost(&smoothvect[0]);
            //if(is_same<eType,int>::value)
            //    printf("\nBefore optimization energy is %lli",
		//	gco->compute_energy());
            //else
            //    printf("\nBefore optimization energy is %f",
		//	gco->compute_energy());
            //CPING("here");
            gco->expansion(10);
            data_out.resize(nb_pixels);
            for(int i=0; i<nb_pixels; i++){
                data_out[i] = gco->whatLabel(i);
            }
            delete gco;
            
	}
        catch (GCException e){
	    e.Report();
            return false;
	}
        convert_vec2mat(data_out,N);
        M = M+N;
    }
    //output = Mat::zeros(height,width,CV_TF);
    // go back to original output.
    for(int i=0; i<height; i++) for(int j=0; j<width; j++)
    {
        data_out[i*width+j] = M.at<int>(i,j);
    }
    return true;
}

bool OptiClass::compute_opt_multiscale(void){
// Parcourt les labels en opérant une dichotomie. 
//     - know label <-> focus ?
//     - M matrice masque des labels en construction. (+I)
//     - n indice d'itération ( n~log2(Nlabels)-1 -> 0 )
//     - I labelshift rang n  ( '00001' << n )
//     - E vectmat energie attache aux données 
//     - L mat energie labels.
//     - start !
// Loop

    
    
    
    return true;
}
















