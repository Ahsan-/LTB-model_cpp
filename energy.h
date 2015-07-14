#ifndef ENERGY
#define ENERGY
#include "tools.h"
#include "sole_params.h"
void set_E_Er_vec(double *E, double *E_r, double *r_vec,int size,
				  spline_1d *spM, struct Params *model_params);
#endif
