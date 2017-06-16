/****************************
*	mgldrw.h
*	cstoribal
*	15-06-17
****************************/

/*
Just holds the class inherrited from mglDraw
*/


#include <mgl2/wnd.h>
#include <vector>
#include <string>
#include <opencv2/core.hpp>

using namespace cv;
using namespace std;


class MyDrawer2D: public mglDraw {
    public:
        int rows,cols,resolution;
        fType imin, imax;
        string filename;
        vector<Mat> vmat_data;
        

        int Draw (mglGraph *gr){
            mglData a;
            a.Create(rows,cols);
            for(int i=0; i<rows; i++)    for(int j=0; j<cols; j++)
            {
                a.a[i+rows*j] = vmat_data[0].at<fType>(i,j);
            }
            
            gr->Title(filename.c_str());
            gr->Light(true);  gr->Rotate(50,60);
            gr->SetRange('z',imin,imax); //gr.SetOrigin(-1,-1,-1); 
            gr->Axis();
            gr->Label('z',"depth",0); 
            gr->Box();  gr->Surf(a);
            return 0;
        }

};
        
