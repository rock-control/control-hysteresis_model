/*
 * =====================================================================================
 *
 *       Filename:  brute_force_optim.cpp
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
#include "ModelOutput.h"

using namespace std;

bouc_wen_model_output::bouc_wen_model_output(const string &strLogFile, const string &strDataFile)
{
    if ( ! ( fileLog = fopen ( strLogFile.c_str(), "w" ) ) )
    {
	printf ( "Can't open LOG file.\n" );
	exit ( 1 );
    }

    if ( ! ( fileMeaData = fopen ( strDataFile.c_str(),"rb" ) ) )
    {
	printf ( "Can't open DATA file.\n" );
	exit ( 1 );
    }	
   
    readMeaInput();
    tmax = mesTime.back();              // Total Simulation period
}

bouc_wen_model_output::~bouc_wen_model_output()
{
    fclose ( fileLog );

}

void bouc_wen_model_output::readMeaInput()
{
    char str[400], sub[100];
    while(1)
    {
	fgets(str, 100, fileMeaData);
	if(feof(fileMeaData))
	    break;

	istringstream iss(str);
	iss >> sub;
	mesTime.push_back(atof(sub));	
	iss >> sub;
	mesTorque.push_back(atof(sub));	// kNm
	iss >> sub;
	mesDeflection.push_back(atof(sub) * 180.0/M_PI );	//Degrees
	iss >> sub;
	mesDeflection_dot.push_back(atof(sub) * 180.0/M_PI );	//Degrees/sec
    }
    fclose ( fileMeaData );
}

void bouc_wen_model_output::setParameters(double *p)
{
    hysModel.setParameters(p);
}

void bouc_wen_model_output::output()
{
	index = 0;
	while (  hysModel.current_time < tmax  )
	{
	    hysModel.getStress(
		    mesTime[index],
                    mesDeflection[index],
		    mesDeflection_dot[index],
		    modelTorque  );
	    
	    if( hysModel.current_time > mesTime[index] )
		index++;

	    fprintf ( fileLog,"%f %f %f %f %f \n"
		    , hysModel.current_time
		    , mesDeflection[index]
		    , modelTorque
		    , mesTorque[index]
		    , mesDeflection_dot[index]
		     ); 
	}
}

