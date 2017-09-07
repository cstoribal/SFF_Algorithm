cmake_minimum_required(VERSION 2.8)

project(ProjSFF)
SET(CMAKE_MODULE_PATH ${ProjSFF_SOURCE_DIR})
##################################################
##OPENCV

find_package(OpenCV REQUIRED)
message(" -- OPENCV include: " ${OpenCV_INCLUDE_DIRS})
message(" -- OPENCV libs dir: " ${OpenCV_LIB_DIR})
message(" -- OPENCV libs: " ${OpenCV_LIBS})
# semble inclure toute la librairie opencv

include_directories (${OPENCV_INCLUDE_DIRS})
message("Cmake running...")
##################################################
##MathGL
##find_package(MathGL)## 

## include_directories (${MathGL_INCLUDE_DIRS})
## message(" -- MathGL include " ${MathGL_INCLUDE_DIRS})
## message(" -- MathGL libs dir: " ${MathGL_LIB_DIR})
## message(" -- MathGL libs: " ${MathGL_LIBS})
find_package(MathGL)
include_directories(${MathGL_INCLUDE_DIRS})

if(MathGL_FOUND)
    #set true for mathgl 

    SET(GCC_COVERAGE_COMPILE_FLAGS "-std=gnu++11")
    message("MathGL found... switching to -std=gnu++11")
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

endif(MathGL_FOUND)




##################################################
file(GLOB LIB_CPP  ./misc/miscdef.cpp ./io/IOWizard.cpp ./io/logs.cpp ./tool/evaluation.cpp ./tool/optimization.cpp ./tool/sharpnessOperator.cpp ./tool/depthEstimator.cpp ./tool/energy2.cpp ./tool/pretreat.cpp ./tool/timer.cpp ./sffClass.cpp  ./misc/testClass.cpp)
ADD_LIBRARY(myLIB SHARED ${LIB_CPP} )
file(GLOB GCO_CPP  ./gco/GCoptimization.cpp ./gco/LinkedBlockList.cpp)
ADD_LIBRARY(gcoLIB SHARED ${GCO_CPP})
##ADD_LIBRARY(mglLIB SHARED ${MathGL_LIBRAIRIES})



add_executable(ProjSFF mainSFF.cpp)
target_link_libraries(ProjSFF myLIB  gcoLIB  ${OpenCV_LIBS} ${MathGL_LIBRARIES} )

