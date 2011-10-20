/*
 * =====================================================================================
 *
 *       Filename:  lev_mar_optim.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  23/10/10 11:10:25
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ajish Babu (), ajish.babu@dfki.de
 *        Company:  DFKI
 *
 * =====================================================================================
 */
#include "LevMarOptimization.h"

using namespace std;

#define GAMMA_BETA_WEIGHT 10

lev_mar_optim::lev_mar_optim(const string &strLogFile, vector<std::string> strDataFile,  int _nDataFiles)
    : mesTime		(_nDataFiles, vector<double>(0,0)),
      mesTorque		(_nDataFiles, vector<double>(0,0)),
      mesDeflection	(_nDataFiles, vector<double>(0,0)),
      mesDeflection_dot	(_nDataFiles, vector<double>(0,0))
{
    firstRun = true;
    if ( ! ( fileLog = fopen ( strLogFile.c_str(), "w" ) ) )
    {
	printf ( "Can't open file.\n" );
	exit ( 1 );
    }
    nDataFiles = _nDataFiles;

    for(int i = 0; i < nDataFiles ; i++)
    { 
	if ( ! ( fileMeaData = fopen ( strDataFile[i].c_str(),"rb" ) ) )
	{
	    printf ( "Can't open file.\n" );
	    exit ( 1 );
	}	

	readMeaInput(i);
	fclose ( fileMeaData );
    }
}

lev_mar_optim::~lev_mar_optim()
{
    fclose ( fileLog );
}

void lev_mar_optim::readMeaInput(int n)
{
    char str[400], sub[100];

    while(1)
    {
	fgets(str, 400, fileMeaData);
	if(feof(fileMeaData))
	    break;
    
	istringstream iss(str);
	iss >> sub;
	mesTime[n].push_back(atof(sub));	
	iss >> sub;
	mesTorque[n].push_back(atof(sub));	
	iss >> sub;
	mesDeflection[n].push_back(atof(sub) * 180.0/M_PI);	
	iss >> sub;
	mesDeflection_dot[n].push_back(atof(sub) * 180.0/M_PI);	
    }
}

void lev_mar_optim::optimize(double *x, double *lp, double *up)
{
    int m = 0, nP = 9;
    for(int i=0; i < nDataFiles; i++)
	m += (mesTime[i].size());

    std::vector<double> fvec(m, 0.0);
    optimumParameter.resize(nP);
    double info[10];

    dlevmar_bc_dif(&torque_error, x, &fvec[0], nP, m, lp, up ,1000, 
	    NULL, info, NULL, NULL, this);

    printf ( "\n\n RESULT \n A: %f beta: %f gamma: %f n: %f a: %f ki: %f gp: %f dO: %f da: %f\n", 
	    		       x[0],    x[1],     x[2], x[3], x[4],  x[5],  x[6], x[7], x[8] );
    printf ( "\n(%f,%f,%f,%f,%f,%f,%f,%f,%f)\n", 
	       x[0],    x[1],     x[2], x[3], x[4],  x[5],  x[6], x[7], x[8]);
    printf ( "\n INFORMATION : \n ");
    printf ( " 		Error at initial Parameter : %f\n " 	, info[1]);
    printf ( " 		||J^T e||_inf : %f\n " 			, info[2]);
    printf ( " 		||Dp||_2 : %f\n " 			, info[3]);
    printf ( " 		\\mu/max[J^T J]_ii ] : %f\n " 		, info[4]);
    printf ( " 		Iterations : %f\n " 			, info[5]);
    printf ( " 		Reason for termination : %f\n "		, info[6]);
    printf ( " 		Function evaluations : %f\n " 		, info[7]);
    printf ( " 		Jacobian evaluations : %f\n" 		, info[8]);
    printf ( " 		Number of linear systems solved : %f\n ", info[9]);
    fprintf(fileLog,     "\n\n RESULT \n A: %f beta: %f gamma: %f n: %f a: %f ki: %f, gp:%f dO: %f da: %f\n",
	    		                   x[0],    x[1],     x[2], x[3], x[4],  x[5],  x[6], x[7], x[8]);
    fflush(fileLog);
}

void torque_error(double *x, double *fvec, int np, int m, void *_adata)
{
    lev_mar_optim *adata = (lev_mar_optim*) _adata;
    adata->torqueError 		= 0.0;

    for(int i=0; i< adata->nDataFiles; i++)
    {
	hysteresis_model::BoucWenModel hysModel;
	hysModel.setParameters(x);
	adata->mIndex = 0;

	while (  hysModel.current_time < adata->mesTime[i].back()  )
	{
	    hysModel.getStress(
		    adata->mesTime[i][adata->mIndex],
		    adata->mesDeflection[i][adata->mIndex],
		    adata->mesDeflection_dot[i][adata->mIndex],
		    adata->modelTorque  );

	    adata->torqueError += pow((adata->modelTorque - adata->mesTorque[i][adata->mIndex] ),2);

	    fvec[adata->mIndex] = adata->modelTorque - adata->mesTorque[i][adata->mIndex] ;

	    if(fabs(x[2]) > x[1])
	    {
		fvec[adata->mIndex] += (fabs(x[2]) - x[1]) * GAMMA_BETA_WEIGHT;
		adata->torqueError += pow(((fabs(x[2]) - x[1]) * GAMMA_BETA_WEIGHT ),2);
	    }
	    
	    if( hysModel.current_time > adata->mesTime[i][adata->mIndex] )
		adata->mIndex++;
	}
    }

    if(adata->firstRun)
    {
	adata->firstRun = false;
	adata->minError = adata->torqueError;
	for(int i=0; i<np; i++)
	    adata->optimumParameter[i] = x[i];
    }
    else
    {
	if(adata->torqueError < adata->minError)
	{
	    adata->minError = adata->torqueError;
	    for(int i=0; i<np; i++)
		adata->optimumParameter[i] = x[i];

	    printf ( "  - A: %f\n    beta: %f\n    gamma: %f\n    n: %f\n    alpha: %f\n    ki: %f\n    gearPlay: %f\n    deflectionOffset: %f\n    dampingConst: %f\n    velocitySmoothFactor: 0.01\n     Err:%.2f\n", 
		    adata->optimumParameter[0], adata->optimumParameter[1], adata->optimumParameter[2], 
		    adata->optimumParameter[3], adata->optimumParameter[4], adata->optimumParameter[5],  
		    adata->optimumParameter[6], adata->optimumParameter[7], adata->optimumParameter[8], adata->minError);

	    fprintf ( adata->fileLog,  "  - A: %f\n    beta: %f\n    gamma: %f\n    n: %f\n    alpha: %f\n    ki: %f\n    gearPlay: %f\n    deflectionOffset: %f\n    dampingConst: %f\n    velocitySmoothFactor: 0.01\n     Err:%.2f\n", 
		    adata->optimumParameter[0], adata->optimumParameter[1], adata->optimumParameter[2], 
		    adata->optimumParameter[3], adata->optimumParameter[4], adata->optimumParameter[5],  
		    adata->optimumParameter[6], adata->optimumParameter[7], adata->optimumParameter[8], adata->minError);
	}
    }

    printf ( "A: %f beta: %f gamma: %f n: %f alpha: %f ki: %f gearPlay: %f deflectionOffset: %f dampingConst: %f velocitySmoothFactor: 0.01  Err:%.2f\n", 
	        x[0],    x[1],     x[2], x[3], x[4],  x[5],  x[6], x[7], x[8], adata->torqueError);

    fflush(adata->fileLog);
}

