#include <iostream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <stdexcept>
#include <cmath>
#include "mass_function.h"
#include "sole_params.h"
using namespace std;

double dMdr(double r, Params *p)
{
	return 0.75*p->OmegaM_out*p->Hoverc_out*p->Hoverc_out*r*r*(2.+p->delta_w 
	       -p->delta_w*tanh((r-p->r0)/p->delta_r));
}

int dMdr_gsl (double r, const double y[], double f[], void *p_)
{
	struct Params *p = (struct Params *) p_;
	f[0] = 0.75*p->OmegaM_out*p->Hoverc_out*p->Hoverc_out*r*r*(2.+p->delta_w 
	       -p->delta_w*tanh((r-p->r0)/p->delta_r));
  return GSL_SUCCESS;
}

void set_M_Mr_vec(double *M_vec, double *Mr_vec, const double *r_vec,int size,
				  struct Params *p,
				  int (*f)(double, const double [], double [], void *))
{
  	int (*ptr) (double t, const double y[], double *dfdy, 
                      double dfdt[], void *args);
    ptr= NULL;                   
	gsl_odeiv2_system sys = { f, ptr, 1, p};

	gsl_odeiv2_driver *d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msadams,
                                   1e-4, 1e-12, 1e-15);//gsl_odeiv2_step_msadams,gsl_odeiv2_step_rkf45,rkck,rk8pd

  double r_init = 0.;
  double y_init[] = {0.0};
	
  for (int i=0; i < size; i++)
    {
      int status = gsl_odeiv2_driver_apply (d, &r_init, r_vec[i], y_init);

      if (status != GSL_SUCCESS)
	{
	  std::cout << "error Massfunction, return value= "<< status << std::endl;
	  break;
	}
      std::cout.precision(5);
      std::cout.width(10);
     // std::cout<< std::setw(4) << r_init<< std::setw(12) << y_init[0]  << std::endl;
      M_vec[i] = y_init[0];
      (*f)(r_init,y_init,&Mr_vec[i],p); 
    }
	gsl_odeiv2_driver_free (d);
}


