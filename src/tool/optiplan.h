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

using namespace std;
using namespace cv;




class OptiPlan{
private:
    IOWizard* ioWizard;
    MyLog* myLog;
    string type;
    size_t nb_labels;
    size_t nb_pixels;
    size_t nb_iterations; // à fixer. Dans les getters, créer une règle pour jeter une alerte en cas de tentative d'accès incorrecte.
    bool     h_thresh; //handler
    size_t*  a_thresh; // full set of thresholds
    size_t** aa_thresh; //pointeurs vers tableaux. A initialiser aussi.
    bool     h_histogram;
    size_t*  a_histogram; //modified histogram for optimization


    bool error(void); //shouts handler state : where dit it fail.
    bool set_histomod(void); // pour l'instant, +1 tout le temps, if type
    bool set_thresh(void); // fill in the threshold values

    bool get_cluster(vector<size_t> thresh, int path, size_t& min, size_t& max);

    bool search_kmean_in(int k, int path, size_t min, size_t max, vector<size_t>& means);
    bool search_otsu(int k, int path, size_t min, size_t max, vector<size_t>& means);



    bool (OptiPlan::* pf_search)(int k, int path, size_t min, size_t max, vector<size_t>& means) ;
    
    // kmeans, 
    // k2means,
    // adaptative...

    // calculer l'erreur rmse théorique.
    // print that in a file.
    // log that

public:
    OptiPlan(); ~OptiPlan();
    bool set_param(string _type, IOWizard* _ioWizard, MyLog* _myLog, size_t _labels, size_t _pixels, vector<unsigned int> histogram);
    
};

#endif
