#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <signal.h>

double globalP[10];
std::string file;

#include "LevMarOptimization.h"
#include "ModelOutput.h"
#include "boost/filesystem.hpp"


void writeComparison(int param)
{
    std::cout << "Writing comparison data " << std::endl << std::flush;
    bouc_wen_model_output modelMeasurementComparison(std::string("ModelFit.txt"), file);
    modelMeasurementComparison.setParameters(globalP);
    modelMeasurementComparison.output();
    exit(1);
}

int main ( int argc, char *argv[] )
{
    if (argc < 2)
    {
	std::cout << "Wrong number of arguments " << std::endl;
	std::cout << "Correct usage: CalibrateCoupling [INPUT_FILE1] [INPUT_FILE2] ..." << std::endl;
	return EXIT_FAILURE;
    }

    std::vector<std::string> strDataFiles;
    for(int i=1; i < argc; i++)
    {
	if(boost::filesystem::exists(argv[i]))
	{
	    strDataFiles.push_back(std::string(argv[i]));
	}
	else
	{
	    std::cout << "File not found : " << argv[i] << std::endl;
	    return EXIT_FAILURE;
	}
    }

//                  A      	     beta 		gamma  	n    		a     		ki   		gp      deflectionOffset  velocitySmooth
//    double initialValue[] = {2.360888     ,1.647663       , 1.555485      , 1.000000      , 0.413352      , 1.578969      , 0.640117      , 0.000000      , -0.017962, 0.01};
    double initialValue[] = {0.00    ,0.1       ,0.0       ,1.0       ,0.0       ,0.001         ,0.0            ,0.0      ,0.0     ,0.01};
    double lowerLimit[] =   {0.001   ,0.1	,-10.0     ,1.0	      ,0.0  	 ,0.001		,0.0	 	,0.0	  ,-5.0};
    double upperLimit[] =   {10      ,10.0	, 10.0     ,20.0      ,1.0	 ,10.0		,10.0	 	,3.0  	  ,+5.0};
    lev_mar_optim oLMO(std::string("CouplingCalibrationLog.txt"), strDataFiles, (int)strDataFiles.size());
    for(int i=0; i<10; i++)
	globalP[i] = initialValue[i];
    file = std::string(argv[1]);

    void (*prev_fn)(int);
    prev_fn = signal (SIGINT,writeComparison);
    if (prev_fn==SIG_IGN) signal (SIGINT,SIG_IGN);

    oLMO.optimize(initialValue, lowerLimit, upperLimit, globalP);

    for(int i=0; i<10; i++)
	globalP[i] = initialValue[i];
    writeComparison(0);

    return EXIT_SUCCESS;
}


