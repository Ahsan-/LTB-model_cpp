#ifndef MASS_FUNCTION
#define MASS_FUNCTION
#include "sole_params.h"

typedef int (*gsl_deriv)(double, const double [], double [], void *);
//void set_M_Mr_vec(double *M, double *M_r, double *r_vec,int size,
//				  struct Params *model_params, gsl_deriv f);

int dMdr_gsl (double r, const double y[], double f[], void *p_);
double dMdr(double r, Params *p);

void set_M_Mr_vec(double *M, double *M_r, const double *r_vec,int size,
				  Params *model_params, gsl_deriv f); 
				  //int (*gsl_deriv)(double, const double [], double [], void *));




#endif
