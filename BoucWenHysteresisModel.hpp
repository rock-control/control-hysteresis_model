/*
 * =====================================================================================
 *
 *       Filename:  bouc_wen_model.h
 *
 *    Description:  Class definition of bouc wen hysterisis model
 *
 *        Version:  1.0
 *        Created:  10/16/09 15:37:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ajish Babu (), ajish.babu@dfki.de
 *        Company:  DFKI
 *
 * =====================================================================================
 */

#ifndef BOUC_WEN_HYSTERESIS_MODEL_H
#define BOUC_WEN_HYSTERESIS_MODEL_H

#include <math.h>
#include "motor_controller/RK4Integrator.hpp"
#include "iostream"


namespace hysteresis_model
{
    class BoucWenModel : public motor_controller::RK4_SIM
    {
	public:
	    BoucWenModel(double samp_time = 0.001, int _simulation_per_cycle = 1,
		    double _initial_time = 0.0, double *_initial_state = NULL);
	    ~BoucWenModel() {};

	    double sampling_period;
	    int simulation_per_cycle;

	    inline void DERIV(const double t, const double *x, 
		    const double *u, double *xdot);

	    bool getStress (double currTime, double strain, double& strainVel, double& stress);
	    void reset(double initTime = 0.0)    {
		init_param( (sampling_period/((double)simulation_per_cycle)),
			initTime, NULL);
	    };
	    void setParameters(double* const p);
	    void setParameters(
		    double _A, 
		    double _beta,
		    double _gamma,
		    double _n,
		    double _a,
		    double _ki,
		    double _gearPlay,
		    double _deflectionOffset,
                    double _dampingConstant,
                    double _velocitySmoothFactor);

	    void getParameters(double *p) const;
            void printParameters() const;

	private:
	    double A;  		
	    double beta;
	    double gamma;
	    double n;
	    double a;
	    double ki;

            double gearPlay;
            double deflectionOffset;
            double dampingConstant;
            double velocitySmoothFactor;

            double torqueGearPlay;
            double torque;

            double prevTime;
            double prevStrain;
            double prevStrainVel;

            bool firstRun;
    };
}
#endif 
