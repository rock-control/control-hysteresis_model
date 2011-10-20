/*
 * =====================================================================================
 *
 *       Filename:  lev_mar_optim.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  23/10/10 11:11:09
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ajish Babu (), ajish.babu@dfki.de
 *        Company:  DFKI
 *
 * =====================================================================================
 */

#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <levmar.h>

#include "BoucWenHysteresisModel.hpp"

#ifndef _LEV_MAR_OPTIM_H
#define _LEV_MAR_OPTIM_H

void torque_error(double *x, double *fvec, int m, int np, void *adata);

class lev_mar_optim 
{
    private:
	void readMeaInput(int);

    public:
	double minError;
	bool firstRun;
	std::vector<double> optimumParameter;
	double mesIntpDeflection, mesIntpDeflection_dot, mesIntpTorque, timeFactor, torqueError;
	double prevDeflection, prevTime, modelTorque;

	int mIndex;
	int nDataFiles;

	FILE *fileLog, *fileMeaData;

	std::vector< std::vector<double> > mesTime, mesTorque, mesDeflection,mesDeflection_dot ; 

	lev_mar_optim(const std::string &strLogFile,std::vector<std::string> strDataFile,  int _nDataFiles);
	~lev_mar_optim();

	void optimize(double *p, double *lp, double *up);
};
#endif
