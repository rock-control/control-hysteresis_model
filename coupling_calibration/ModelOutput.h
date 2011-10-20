/*
 * =====================================================================================
 *
 *       Filename:  bouc_wen_model_output.h
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

#include "BoucWenHysteresisModel.hpp"

#ifndef _BOUC_WEN_MODEL_OUTPUT_H
#define _BOUC_WEN_MODEL_OUTPUT_H

class bouc_wen_model_output 
{
    private:
	hysteresis_model::BoucWenModel hysModel;
	double tmax;              // Total Simulation period

	double mesIntpDeflection, mesIntpDeflection_dot, mesIntpTorque, timeFactor;
	double prevDeflection, prevTime, modelTorque;

	int index;

	FILE *fileLog, *fileMeaData;
	std::vector<double> mesTime, mesTorque, mesDeflection, mesDeflection_dot; 

	void readMeaInput();

    public:
	bouc_wen_model_output(const std::string &strLogFile, const std::string &strDataFile);
	~bouc_wen_model_output();
	void setParameters(double *p);
	void output();
};
#endif
