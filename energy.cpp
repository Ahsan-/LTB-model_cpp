#include <iostream>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <cmath>
#include <stdexcept>
#include "sole_params.h"
#include "tools.h"

struct cosmic_time_param{
					/** Specialized struct to be used with cosmic_time. Do not use this 
					     struct outside energy.cpp
					**/ 
                   double E;
                   double M;
                   double Lambda;
                   double age;
                   };
                   
double cosmic_time(double R, void *p_)
{
	/** For convenience and optimization the following symbols have been redefined 
		R/R0     --> R
		2E/r^2   --> E
		2M/r^3   --> M
		Lambda/3 -->Lambda 
		
		Eq. (2.2) in ``Structures in the Universe by Exact Method'' by Krzystof 
		Bolejko etal. rearranged to give time as the independent variable, as the 
		normalized Scale factor runs from 0 to 1 in the integrand.
	**/
	struct cosmic_time_param *p = (struct cosmic_time_param*) p_;
	return sqrt(R/(p->E*R + p->M+p->Lambda*R*R*R));
}

double cosmic_time_minus_age(double E, void *p_)
{
	struct cosmic_time_param *p = (struct cosmic_time_param*) p_;
	
	//generously allocate workspace 1000 doubles
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	double result, error;
	
	gsl_function Integrand;
	Integrand.function = &cosmic_time;
	p->E = E;
	Integrand.params   = p;
//Note that since the age of the universe is of the order 13.7*306 in Mpc
//the relative and absolute errors can be set to acheive accuracy up to only
//a few decimal places.
//int gsl_integration_qag (const gsl_function * f, double a, double b, double epsabs, 
//                         double epsrel, size_t limit, int key, 
//                        gsl_integration_workspace * workspace, double * result, 
//                         double * abserr)
	gsl_integration_qag (&Integrand, 0.,1.,1.49e-6,1.49e-4,1000, //0.,1.,1.49e-4,1.49e-4,1000//1.49e-16,1.49e-12,1000,//
	                     6,//1 to 6
	                     w,&result,&error);
	gsl_integration_workspace_free(w);
//size_t neval;
//gsl_integration_qng (&Integrand, 0.,1.,1.49e-4,1.49e-4,//1.49e-16,1.49e-12
//	                     &result,&error,&neval);

return p->age-result;
}

double get_E(struct cosmic_time_param *p)
{
 int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double E = 0.;
  double E_low = 0., E_high = 0.1;//1e-6;
  //double E_low = -1e3, E_high = 1e3;//1e-6;
  gsl_function F;

  F.function = &cosmic_time_minus_age; 
  F.params = p;

  T = gsl_root_fsolver_brent;//bisection;//;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, E_low, E_high);

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      E = gsl_root_fsolver_root (s);
      E_low = gsl_root_fsolver_x_lower (s);
      E_high = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (E_low, E_high,  
                                       4.44e-16, 4.44e-15);
    
    }
  while (status == GSL_CONTINUE && iter < max_iter);
   if (status != GSL_SUCCESS)
      {
        std::cout <<  "Energy.cpp, get_E did not converge " << std::endl;
      }
  
  gsl_root_fsolver_free (s);
  
  return E;
}


void set_E_Er_vec(double *E, double *E_r, double *r_vec,int size,
				  spline_1d *spM, struct Params *model_params)
{
	//cosmic_time_param p;
	//p.Lambda = model_params->Lambda/3.;
	//p.age = model_params->age;
	
	#pragma omp parallel for
	for (int i=0; i<size; i++)
	{
	cosmic_time_param p;
	p.Lambda = model_params->Lambda/3.;

	p.age = model_params->age;

	double r = r_vec[i];
	p.M   = 2.*((*spM)(r,0))/(r*r*r); //see cosmic_time_param definition
	
	/**
		For E=0, the time integral has the solution:
		2/3*ln(sqrt(Lambda)R^1.5+sqrt(Lambda*R^3+M))/sqrt(Lambda) - 1/3*ln(M)/sqrt(Lambda)
		Note that the integral is performed in the notation
		R/R0     --> R
		2E/r^2   --> E
		2M/r^3   --> M
		Lambda/3 -->Lambda
		Furthermore since M<<1 avoid log(M) which can cause round off errors
	**/
	//double max_age = 2./3.*log(sqrt(p.Lambda)+sqrt(p.Lambda*+p.M))/sqrt(p.Lambda) 
	//             - 1./3.*log(p.M)/sqrt(p.Lambda);
	double max_age =log((2.*p.Lambda+p.M+2.*sqrt(p.Lambda*(p.Lambda+p.M)))/
	                     p.M)/(3.*sqrt(p.Lambda));
	if (  max_age < p.age)
	{
		std::cout << "At r = " << r << " " <<
		"Current parameters allow a maximum age of " << max_age << std::endl;
		std::cout << "You specified age of the universe to be " << p.age << std::endl;
		//std::cout << "Decrease Omega_Lambda and or age " << std::endl;
		throw std::invalid_argument("Decrease H_out and or age ");
	}
	else
	{
		E[i] = get_E(&p)*(r*r/2.); //see cosmic_time_param definition
	}
		
	}

	spline_1d spE_r(r_vec, E, size);
	for (int i=0; i<size; i++)
	{
		E_r[i] = spE_r(r_vec[i],1);
	}
	
}


