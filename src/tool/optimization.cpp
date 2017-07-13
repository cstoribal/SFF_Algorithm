/****************************
*	ProjSFF
*	optimization.cpp
*	cstoribal
*	03-05-17
****************************/

#include "optimization.h"


OptiClass::OptiClass(){ 
    this->set = false;
    this->nbs_set = 0;
    this->connexity = 4; //TODO
};
OptiClass::~OptiClass(){
    if(this->nbs_set)
    {
        delete [] nbs_nb;
        delete [] nbs_n1D;
        delete [] nbs_n;
        delete [] nbs_w1D;
        delete [] nbs_w;
        this->nbs_set=0;
    }
};

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
    if(0){
        return compute_opt_binary();
    }
    if(1){
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
bool OptiClass::reset(void){
    //copy weight matrix to reset it at initial status
    for(int k=0; k<nb_pixels; k++) for(int c=0; c<connexity; c++)
    {
        nbs_wk[k][c] = nbs_w[k][c];
    }
    return true;
}


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
        getxy[k]=Point(j,i); 
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



bool OptiClass::set_allneighbors(void){
    // needed if optimisation type general graph
    if(nbs_set == 1){
        myLog->a("Warning, re-setting neighbors without freeing memory");
        delete [] nbs_nb;
        delete [] nbs_n1D;
        delete [] nbs_n;
        delete [] nbs_w1D;
        delete [] nbs_w;
        this->nbs_set=0;
    }
    // sets up neighbornumarray, neighborarray, weightarray !
    // aka nbs_nb, nbs_n, nbs_w;
    Mat1E borders = Mat::zeros(height,width,CV_TE)+1;
    this->neighbor_mat = Mat::zeros(5,5,CV_TE);
    if(connexity=4)
    {
        neighbor_mat.at<eType>(2,1) = 1;
        neighbor_mat.at<eType>(1,2) = 1;
        neighbor_mat.at<eType>(3,2) = 1;
        neighbor_mat.at<eType>(2,3) = 1;
    }

    filter2D(borders, borders, -1, this->neighbor_mat, Point(2,2), 0, BORDER_CONSTANT );
    
    nbs_nb  = new int[nb_pixels];
    nbs_n1D = new int[nb_pixels*connexity];
    nbs_w1D = new eType[nb_pixels*connexity];
    nbs_wk1D = new eType[nb_pixels*connexity];
    
    nbs_n   = new int*[nb_pixels];
    nbs_w   = new eType*[nb_pixels];
    nbs_wk   = new eType*[nb_pixels];


    for(int k=0; k<nb_pixels; k++){
        nbs_n[k] = &nbs_n1D[k*connexity];
        nbs_w[k] = &nbs_w1D[k*connexity];
        nbs_wk[k] = &nbs_wk1D[k*connexity];
    }


    int nbs_tmp_nb;

    for(int k=0; k<nb_pixels; k++)
    {
        nbs_tmp_nb = 0; // Number of neighbors
        for(int i=-2; i<3; i++) for(int j=-2; j<3; j++)
        {
            if( neighbor_mat.at<eType>(2+i,2+j) != 0 )
            if( getxy[k].y+i >= 0 && getxy[k].y+i < height )
            if( getxy[k].x+j >= 0 && getxy[k].x+j < width  )
            {
                nbs_n[k][nbs_tmp_nb]=getrank.at<int>
				(getxy[k].y+i,getxy[k].x+j);
                nbs_w[k][nbs_tmp_nb]=neighbor_mat.at<eType>(2+i,2+j);
                nbs_tmp_nb++;
            }
        }
        nbs_nb[k]=nbs_tmp_nb;
    }
    
/*
    for (int y = 0; y < nb_pixels; ++y)
    {
        for (int x = 0; x < connexity; ++x)
        {
            std::cout << std::hex << &(nbs_n[y][x]) << ' ';
        }
    }
    std::cout << "\n\n";

    // Print the array
    for (int y = 0; y < nb_pixels; y++)
    {
        std::cout << std::hex << &(nbs_n[y][0]) << std::dec;
        std::cout <<  "; " << (nbs_nb[y]);
        std::cout << ": ";
        for (int x = 0; x < connexity; x++)
        {
            std::cout << nbs_n[y][x] << ' ';
        }
        std::cout << std::endl;
    }
*/
    
    
    
    nbs_set = 1;
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
    if(!energyClass->getCrossLabelMatrix(labels,tmp_mat)) return false;
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
    Mat1i M = Mat::zeros(height,width,CV_32S);
    Mat1i N = Mat::zeros(height,width,CV_32S);
    Mat1E Dp = Mat::zeros(height,width,CV_TE);
    Mat1E Dn = Mat::zeros(height,width,CV_TE);
    int I = 1; int J;
    double range = floor(log2(nb_labels-1)+1); //nb of digits necessary
    Mat1E lmat = Mat::zeros(nb_labels,nb_labels,CV_TE);
    energyClass->getCrossLabelMatrix(labels,lmat);
    smoothvect.resize(2*2);
    
    vector<Mat1E> Edata(nb_labels);
    energyClass->getDataEnergy_3DMatrix(this->labels,Edata);
    vector<Mat1E> Edata_slice(2);
    Edata_slice[0] = Mat::zeros(height,width,CV_TE);
    Edata_slice[1] = Mat::zeros(height,width,CV_TE);

    for(int i=0; i<2; i++) for(int j=0; j<2; j++)
    {
        smoothvect[i*2+j]=lmat.at<eType>(0+i,0+j);
    }
    for(int n=(int)range-1; n>=0; n--)
    {
        I = 1<<n; J = I-1;
        //CPING2("passe n",n);
        //CPING2("I(n)",I);
        for(int i=0; i<height; i++) for(int j=0; j<width; j++)
        {
            // TODO check dp or dn (orientation du delta)
            if(M.at<int>(i,j)+I >= nb_labels)
            {
                Edata_slice[0].at<eType>(i,j)=Edata[M.at<int>(i,j)-1]
			.at<eType>(i,j)+Dn.at<eType>(i,j);
                Edata_slice[1].at<eType>(i,j)=Edata[M.at<int>(i,j)]
			.at<eType>(i,j)+Dp.at<eType>(i,j);
            }
            else
            {
                Edata_slice[0].at<eType>(i,j) = Edata[M.at<int>(i,j)+J]
			.at<eType>(i,j)+Dn.at<eType>(i,j);
                Edata_slice[1].at<eType>(i,j) = Edata[M.at<int>(i,j)+I]
			.at<eType>(i,j)+Dp.at<eType>(i,j);
            }
        }
        this->convert_mat2labvec(Edata_slice,data_in);

        try{
            gco = new GCoptimizationGeneralGraph(this->nb_pixels,2);
            ((GCoptimizationGeneralGraph*)gco)
			->setAllNeighbors(nbs_nb,nbs_n,nbs_wk);
            gco->setDataCost(&data_in[0]);
            gco->setSmoothCost(&smoothvect[0]);
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

        for(int k=0; k<nb_pixels; k++) for(int i=0; i<nbs_nb[k]; i++)
        {
            if(nbs_wk[k][i]!=0)
            if(data_out[k]!=data_out[ nbs_n[k][i]])
            {//i°th neighbor j of k°th pixel is different => break w_jk
                nbs_wk[k][i] = 0;
                if(data_out[k]>data_out[ nbs_n[k][i] ])
                    Dp.at<eType>(getxy[k]) += smoothvect[1];
                else
                    Dn.at<eType>(getxy[k]) += smoothvect[1];
            }

        }
        
        convert_vec2mat(data_out,N);
        N = N*I;
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
















