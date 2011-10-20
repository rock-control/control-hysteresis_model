/*
 * =====================================================================================
 *
 *       Filename:  bouc_wen_model.cpp
 *
 *    Description:  Class implementing the hysteresis model
 *
 *        Version:  1.0
 *        Created:  10/16/09 15:36:43
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ajish Babu (), ajish.babu@dfki.de
 *        Company:  DFKI
 *
 * =====================================================================================
 */

#include "BoucWenHysteresisModel.hpp"

using namespace hysteresis_model;

BoucWenModel::BoucWenModel(double _sampling_period, int _simulation_per_cycle,
	double _initial_time, double *_initial_state) 
: RK4_SIM(1, 1, (_sampling_period/((double)_simulation_per_cycle)),
	_initial_time, _initial_state)
{
    // Initialize parameters
    sampling_period = _sampling_period;
    simulation_per_cycle = _simulation_per_cycle;

    firstRun = true;
}

inline void BoucWenModel::DERIV(const double t, const double *x, 
	const double *u, double *xdot)
{
    // The actual Bouc-Wen-Baber-Noori model which calculates the 
    // hysteretic defection using this intergral
    xdot[0] =
	( A * u[0]
	  - (
	        beta * fabs(u[0]) * pow(fabs(x[0]), n-1) * x[0]
	      + gamma*      u[0]  * pow(fabs(x[0]), n  ) 
	    )
	);
}

bool BoucWenModel::getStress(double currTime, double currStrain, double& strainVel, double& stress)
{
//    if(currTime == 0.0)
//    {
//      reset(currTime);
//      prevTime 		= currTime;
//      prevStrain 	= currStrain;
//      prevStrainVel	= 0; 
//
//      firstRun = false;
//      strainVel = 0;
//      stress = a*ki*(currStrain + deflectionOffset);
//      return true;
//    }
//
//    // deflection velocity without smoothing
//    strainVel = (currStrain - prevStrain) / (currTime - prevTime);
//
//    // First order smoothing
//    strainVel = 
//      velocitySmoothFactor*strainVel + (1.0 - velocitySmoothFactor) * prevStrainVel;
//
//    prevTime 		= currTime;
//    prevStrain 	        = currStrain;
//    prevStrainVel	= strainVel; 

    ctrl_input[0] = strainVel;
    for (int ii=0; ii < simulation_per_cycle; ii++)
    {
	solve();
    }

    // Returns the stress which is weighted sum of the actual strain part and
    // the hysteretic strain part
    torque = a*ki*(currStrain + deflectionOffset) + (1-a)*ki*plant_state[0];

 
    // Damping effect of the coupling
    torque -= dampingConstant * strainVel;
    // Gear play is applied when the torque changes
    if(torque <= torqueGearPlay && torque >= -torqueGearPlay)
    {
      torque = 0.0;
    }
    else
    {
      torque -= (torque/fabs(torque)) * torqueGearPlay; 
    }

    // resets when the torque is computed to be NaN
    if(torque != torque)
    {
// std::cout << " NaN @ hysteresis model" << std::endl;
      firstRun = true;
      return false;
    }
    stress = torque;
    return true; 
}

void BoucWenModel::getParameters(double *p) const
{
    p[0] = A;  		
    p[1] = beta;
    p[2] = gamma;
    p[3] = n;
    p[4] = a;
    p[5] = ki;

    p[6] = gearPlay;	
    p[7] = deflectionOffset;	
    p[8] = dampingConstant;	

    p[9] = velocitySmoothFactor;	
}

void BoucWenModel::printParameters() const
{
    std::cout <<"A                :"<<A        << std::endl
              <<"beta             :"<<beta     << std::endl
              <<"gamma            :"<<gamma    << std::endl
              <<"n                :"<<n        << std::endl
              <<"a                :"<<a        << std::endl
              <<"ki               :"<<ki       << std::endl
              <<"gearPlay         :"<<gearPlay          << std::endl
              <<"deflectionOffset :"<<deflectionOffset  << std::endl
              <<"dampingConstant  :"<<dampingConstant	<< std::endl
              <<"velocitySmoothFactor  :"<<velocitySmoothFactor	<< std::endl;
}

void BoucWenModel::setParameters(double* const p)
{
    setParameters(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9]);	
}

void BoucWenModel::setParameters(
		    double _A, 
		    double _beta,
		    double _gamma,
		    double _n,
		    double _a,
		    double _ki,
		    double _gearPlay,
		    double _deflectionOffset,
                    double _dampingConstant,
                    double _velocitySmoothFactor)
{
    A	 = _A;  		
    beta = _beta;
    gamma= _gamma;
    n    = _n;
    a    = _a;
    ki   = _ki;

    gearPlay = _gearPlay;
    deflectionOffset = _deflectionOffset;
    dampingConstant = _dampingConstant;

    velocitySmoothFactor = _velocitySmoothFactor;

    torqueGearPlay =  a*ki*gearPlay/2.0; 
}
