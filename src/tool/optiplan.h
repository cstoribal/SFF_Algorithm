/****************************
*	ProjSFF
*	optiplan.h
*	cstoribal
*	13-10-17
****************************/


#ifndef OPTI_P_H_
#define OPTI_P_H_





#include "../misc/miscdef.h"
#include "../io/IOWizard.h"
#include "../tool/utils.h"

using namespace cv;

class OptiStored {
private:
    std::string type;
    size_t   stored_label;
    size_t   nb_iterations;
    bool     h_thresh;
    //size_t*  a_thresh;
    vector<vector<size_t> > vv_thresh;
    bool     h_centroid;
    //size_t*  a_centroid;
    vector<vector<size_t> > vv_centroid;
    bool     h_rmse;
    vector<float>   v_rmse;
    bool     h_rmsegt;
    vector<float>   v_rmsegt;
public:
    OptiStored(); ~OptiStored();
    bool copy_from_plan(std::string& _type, size_t _label, 
	size_t _nb_iter, 
	const vector<vector<size_t> > 	& _vv_thresh, 
	const vector<vector<size_t> > 	& _vv_centroid, 
	const vector<float> 		& _v_rmse);
    bool get_data_pointers(std::string& _type, size_t& _nb_iter,
	vector<vector<size_t> >	& _vv_thresh, 
	vector<vector<size_t> >	& _vv_centroid, 
	vector<float>		& _v_rmse);
    bool get_data_pointers(std::string& _type, size_t& _nb_iter,
	vector<vector<size_t> >	& _vv_thresh, 
	vector<vector<size_t> >	& _vv_centroid, 
	vector<float>		& _v_rmse, 
	vector<float>		& _v_rmsegt);
    bool get_rmse_at_iter(const size_t & iter, float & rmse);
    bool get_rmse(vector<float> & _v_rmse);
    bool get_type(std::string & _type);
    size_t get_iter(void);
    bool set_rmsegt(vector<float> _v_rmsegt);
};


class OptiPlan{
private:
    IOWizard* ioWizard;
    MyLog* myLog;
    std::string type;
    size_t nb_labels;
    size_t upperbound_iterations;
    size_t nb_pixels;
    size_t nb_iterations; // fixé après le set_thresh. Dans les getters, créer une règle pour jeter une alerte en cas de tentative d'accès incorrect. ça carrive
    size_t nb_storedplans;
    size_t idx_bestplan;
    cv::Mat1i gtl_mat;
    bool verbose; //debug;

///////////////////////////////////
    bool            h_thresh; //handler
    //vector<size_t>  v_thresh; // full set of thresholds
    vector<vector<size_t> > vv_thresh;
    bool            h_centroid;
    //vector<size_t>  v_centroid;
    vector<vector<size_t> > vv_centroid;
    //size_t** aa_thresh; //pointeurs vers tableaux. A initialiser aussi.
    bool            h_histogram;
    bool            h_histogram0;
    vector<size_t>  v_histogram0;
    vector<size_t>  v_histogram; //modified histogram for optimization
    bool            h_rmse;
    vector<float>   v_rmse;   // RMSE théorique % iteration pour une méthode donnée
    bool            h_rmsegt;
    vector<float>   v_rmsegt; // RMSE réel. Doit être set avec une mat ext
///////////////////////////////////
    
    std::vector<OptiStored> stored_set;
    
    bool error(void); //shouts handler state : where dit it fail.
    bool error(const std::string& text); 
    bool set_histomod(void); // pour l'instant, +1 tout le temps, if type
    bool set_thresh(void); // fill in the threshold values
    bool compute_RMSE(void); // okay !
    bool show_RMSE_elt_n(FILE* gnuplot, size_t idx);
    bool gnuplot_vect(FILE* gnuplot, vector<float> vect);


    bool generic_2Dvector_to_1Darray(const vector<vector<size_t> >& vect, bool& handler, size_t*& array);
    bool generic_2Dvector_to_2Darray(const vector<vector<size_t> >& vect, bool& handler, size_t*& array, size_t**& aarray);


    // Fonctions de recherche de threshold : 
    // appartenant au cluster ]min;max]
    // renvoyer 
    bool (OptiPlan::* pIT_search_thresh)(
	size_t min, size_t max, size_t& thresh, size_t& c1, size_t& c2);
    bool IT_search_binary(
	size_t min, size_t max, size_t& thresh, size_t& c1, size_t& c2);
    bool IT_search_binary_v2(
	size_t min, size_t max, size_t& thresh, size_t& c1, size_t& c2);
    bool IT_search_adapt( // définition ?
	size_t min, size_t max, size_t& thresh, size_t& c1, size_t& c2);
    bool IT_search_median(
	size_t min, size_t max, size_t& thresh, size_t& c1, size_t& c2);
    bool IT_search_otsu(
	size_t min, size_t max, size_t& thresh, size_t& c1, size_t& c2);
    bool IT_search_otsu_v0(
	size_t min, size_t max, size_t& thresh, size_t& c1, size_t& c2);
    bool IT_search_2means(
	size_t min, size_t max, size_t& thresh, size_t& c1, size_t& c2);
    

    bool IT_search_1surX( // non itératif ?
	size_t min, size_t max, size_t& thresh, size_t& c1, size_t& c2);
    
    
    bool IT_sets_centroids; //Check if centroids are set or not

    // print that in a file.
    // log that

public:
    OptiPlan(); ~OptiPlan();
    bool set_logs(IOWizard* _ioWizard, MyLog* _myLog);
    bool set_groundtruth(cv::Mat1i gt_label_mat);
    bool set_param(string _type, size_t _labels, size_t _pixels, std::vector<size_t> histogram, bool store=false, bool reset=true);
    bool store_setting(bool reset=true);
    bool load_method(size_t idx_storage);
    bool reset_memory(void);
    bool get_best_method_at_it(size_t iter, size_t & idx_storage, std::string& type, bool& active);
    bool get_best_method_at_it_pointers(size_t iter,
	std::string& _type, size_t& _nb_iter,
	vector<vector<size_t> >	& _vv_thresh, 
	vector<vector<size_t> >	& _vv_centroid, 
	vector<float>		& _v_rmse);

    //process external data 
    bool get_ThreshedMatrix(const cv::Mat1i & mat_in, cv::Mat1i & mat_out, size_t iteration, int method=-1);

    bool get_thresholds_n_centroids(int method, // -1 current, or stored
			std::string & _type,
			size_t & _nb_iterations, 
			vector<vector<size_t> > & _vv_thresh, 
			vector<vector<size_t> > & _vv_centro);

    bool get_nb_storedmethods(size_t & _nb_storedplans);

    // Visualisation
    bool show_RMSE(const string& filename);
    bool show_all_RMSE(const string& filename);
    bool show_all_RMSE2(const string& filename);

    bool show_thresh_plan(std::string filename, int kplan);
    bool show_all_thresh_plans(std::string filename);

    bool write_all_ThreshedMatrix(
	cv::Mat1i & mat_in, std::string folder="threshed", bool set_rgt=true);

    bool computeCrossRMSEperf_andLog(void);
    bool addToLog(void);
    
    
    
    
    
};

#endif
