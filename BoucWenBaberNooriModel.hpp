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

#ifndef BOUC_WEN_BABER_NOORI_MODEL_H
#define BOUC_WEN_BABER_NOORI_MODEL_H

#include <math.h>
#include "RK4Integrator.hpp"
#include "iostream"


namespace hysteresis_model
{
    class BoucWenBaberNooriModel : public RK4_SIM
    {
	public:
	    BoucWenBaberNooriModel(double samp_time = 0.001, int _simulation_per_cycle = 1,
		    double _initial_time = 0.0, double *_initial_state = NULL);
	    ~BoucWenBaberNooriModel() {};

	    double sampling_period;
	    int simulation_per_cycle;

	    inline void DERIV(const double t, const double *x, 
		    const double *u, double *xdot);

	    double getStress (double strain, double strainVel);
	    void reset()    {
		init_param( (sampling_period/((double)simulation_per_cycle)),
			0.0, NULL);
	    };
	    void setParameters(double* const p);
	    void setParameters(
		    double _A, 
		    double _beta,
		    double _gamma,
		    double _n,
		    double _a,
		    double _ki,
		    double _nu,
		    double _eta,
		    double _h,
                    double _gearPlay);

	    void getParameters(double *p) const;
            void printParameters() const;

	private:
	    double A;  		
	    double beta;
	    double gamma;
	    double n;
	    double a;
	    double ki;
	    double nu;
	    double eta;

	    double h;
            double gearPlay;
            double torqueGearPlay;
            double torque;
    };
}
#endif 


