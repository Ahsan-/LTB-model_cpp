#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>
#include <iomanip>
#include <cmath>
#include <iostream>
//#include "model_params.h"
#include "sole_params.h"
#include "LTB_background.h"
using namespace std;



/**
The idea behind the series solutions is to first write 
diff(t,R) = sqrt(3 R) / sqrt(Lambda R^3 + 6ER + 6M) 
          = sqrt(R/2M)/ sqrt((Lambda R^3 / 6M) + (ER/M) + 1)  
In all solutions of interest M>0 i.e., strictly positive, therefore the 
series expansion for (1+x)^(-1/2) can be used where x = (Lambda R^3 / 6M) + (ER/M).
R can always be chosen to be sufficiently small so the expansion 
(1+x)^(-1/2) = \sum_{n=0}^{k} \binom{k, n} x^n 
will be accurate even for small k. The series fails at late times, when R is large 
and x can not be gauranteed to be small.
**/  

double t_series(double R,double E,double M,double Lambda)
{
	/**
	Return R.H.S of the t series
	**/
	
double	t = (-0.12155e5 / 0.1376256e7 * sqrt(0.2e1) * pow(E, 0.9e1) * pow(M, -0.19e2 / 0.2e1) - 0.143e3 / 0.12288e5 * sqrt(0.2e1) * pow(E, 0.6e1) * Lambda * pow(M, -0.15e2 / 0.2e1) - 0.5e1 / 0.1536e4 * sqrt(0.2e1) * pow(E, 0.3e1) * Lambda * Lambda * pow(M, -0.11e2 / 0.2e1) - 0.5e1 / 0.72576e5 * sqrt(0.2e1) * pow(Lambda, 0.3e1) * pow(M, -0.7e1 / 0.2e1)) * pow(R, 0.21e2 / 0.2e1) + (0.6435e4 / 0.622592e6 * sqrt(0.2e1) * pow(E, 0.8e1) * pow(M, -0.17e2 / 0.2e1) + 0.231e3 / 0.19456e5 * sqrt(0.2e1) * pow(E, 0.5e1) * Lambda * pow(M, -0.13e2 / 0.2e1) + 0.35e2 / 0.14592e5 * sqrt(0.2e1) * E * E * Lambda * Lambda * pow(M, -0.9e1 / 0.2e1)) * pow(R, 0.19e2 / 0.2e1) + (-0.429e3 / 0.34816e5 * sqrt(0.2e1) * pow(E, 0.7e1) * pow(M, -0.15e2 / 0.2e1) - 0.105e3 / 0.8704e4 * sqrt(0.2e1) * pow(E, 0.4e1) * Lambda * pow(M, -0.11e2 / 0.2e1) - 0.5e1 / 0.3264e4 * sqrt(0.2e1) * E * Lambda * Lambda * pow(M, -0.7e1 / 0.2e1)) * pow(R, 0.17e2 / 0.2e1) + (0.77e2 / 0.5120e4 * sqrt(0.2e1) * pow(E, 0.6e1) * pow(M, -0.13e2 / 0.2e1) + 0.7e1 / 0.576e3 * sqrt(0.2e1) * pow(E, 0.3e1) * Lambda * pow(M, -0.9e1 / 0.2e1) + sqrt(0.2e1) * Lambda * Lambda * pow(M, -0.5e1 / 0.2e1) / 0.1440e4) * pow(R, 0.15e2 / 0.2e1) + (-0.63e2 / 0.3328e4 * sqrt(0.2e1) * pow(E, 0.5e1) * pow(M, -0.11e2 / 0.2e1) - 0.5e1 / 0.416e3 * sqrt(0.2e1) * E * E * Lambda * pow(M, -0.7e1 / 0.2e1)) * pow(R, 0.13e2 / 0.2e1) + (0.35e2 / 0.1408e4 * sqrt(0.2e1) * pow(E, 0.4e1) * pow(M, -0.9e1 / 0.2e1) + sqrt(0.2e1) * E * Lambda * pow(M, -0.5e1 / 0.2e1) / 0.88e2) * pow(R, 0.11e2 / 0.2e1) + (-0.5e1 / 0.144e3 * sqrt(0.2e1) * pow(E, 0.3e1) * pow(M, -0.7e1 / 0.2e1) - sqrt(0.2e1) * Lambda * pow(M, -0.3e1 / 0.2e1) / 0.108e3) * pow(R, 0.9e1 / 0.2e1) + 0.3e1 / 0.56e2 * sqrt(0.2e1) * E * E * pow(R, 0.7e1 / 0.2e1) * pow(M, -0.5e1 / 0.2e1) - sqrt(0.2e1) * E * pow(R, 0.5e1 / 0.2e1) * pow(M, -0.3e1 / 0.2e1) / 0.10e2 + sqrt(0.2e1) * pow(R, 0.3e1 / 0.2e1) * pow(M, -0.1e1 / 0.2e1) / 0.3e1;

	
	return t;
}

double dt_dR_series(double R,double E,double M,double Lambda)
{
	/**
	Returns R.H.S of the t series when it is differentiated w.r.t R. Partial 
	derivative.
	**/	
double	dt_dR = 0.21e2 / 0.2e1 * (-0.12155e5 / 0.1376256e7 * sqrt(0.2e1) * pow(E, 0.9e1) * pow(M, -0.19e2 / 0.2e1) - 0.143e3 / 0.12288e5 * sqrt(0.2e1) * pow(E, 0.6e1) * Lambda * pow(M, -0.15e2 / 0.2e1) - 0.5e1 / 0.1536e4 * sqrt(0.2e1) * pow(E, 0.3e1) * Lambda * Lambda * pow(M, -0.11e2 / 0.2e1) - 0.5e1 / 0.72576e5 * sqrt(0.2e1) * pow(Lambda, 0.3e1) * pow(M, -0.7e1 / 0.2e1)) * pow(R, 0.19e2 / 0.2e1) + 0.19e2 / 0.2e1 * (0.6435e4 / 0.622592e6 * sqrt(0.2e1) * pow(E, 0.8e1) * pow(M, -0.17e2 / 0.2e1) + 0.231e3 / 0.19456e5 * sqrt(0.2e1) * pow(E, 0.5e1) * Lambda * pow(M, -0.13e2 / 0.2e1) + 0.35e2 / 0.14592e5 * sqrt(0.2e1) * E * E * Lambda * Lambda * pow(M, -0.9e1 / 0.2e1)) * pow(R, 0.17e2 / 0.2e1) + 0.17e2 / 0.2e1 * (-0.429e3 / 0.34816e5 * sqrt(0.2e1) * pow(E, 0.7e1) * pow(M, -0.15e2 / 0.2e1) - 0.105e3 / 0.8704e4 * sqrt(0.2e1) * pow(E, 0.4e1) * Lambda * pow(M, -0.11e2 / 0.2e1) - 0.5e1 / 0.3264e4 * sqrt(0.2e1) * E * Lambda * Lambda * pow(M, -0.7e1 / 0.2e1)) * pow(R, 0.15e2 / 0.2e1) + 0.15e2 / 0.2e1 * (0.77e2 / 0.5120e4 * sqrt(0.2e1) * pow(E, 0.6e1) * pow(M, -0.13e2 / 0.2e1) + 0.7e1 / 0.576e3 * sqrt(0.2e1) * pow(E, 0.3e1) * Lambda * pow(M, -0.9e1 / 0.2e1) + sqrt(0.2e1) * Lambda * Lambda * pow(M, -0.5e1 / 0.2e1) / 0.1440e4) * pow(R, 0.13e2 / 0.2e1) + 0.13e2 / 0.2e1 * (-0.63e2 / 0.3328e4 * sqrt(0.2e1) * pow(E, 0.5e1) * pow(M, -0.11e2 / 0.2e1) - 0.5e1 / 0.416e3 * sqrt(0.2e1) * E * E * Lambda * pow(M, -0.7e1 / 0.2e1)) * pow(R, 0.11e2 / 0.2e1) + 0.11e2 / 0.2e1 * (0.35e2 / 0.1408e4 * sqrt(0.2e1) * pow(E, 0.4e1) * pow(M, -0.9e1 / 0.2e1) + sqrt(0.2e1) * E * Lambda * pow(M, -0.5e1 / 0.2e1) / 0.88e2) * pow(R, 0.9e1 / 0.2e1) + 0.9e1 / 0.2e1 * (-0.5e1 / 0.144e3 * sqrt(0.2e1) * pow(E, 0.3e1) * pow(M, -0.7e1 / 0.2e1) - sqrt(0.2e1) * Lambda * pow(M, -0.3e1 / 0.2e1) / 0.108e3) * pow(R, 0.7e1 / 0.2e1) + 0.3e1 / 0.16e2 * sqrt(0.2e1) * E * E * pow(R, 0.5e1 / 0.2e1) * pow(M, -0.5e1 / 0.2e1) - sqrt(0.2e1) * E * pow(R, 0.3e1 / 0.2e1) * pow(M, -0.3e1 / 0.2e1) / 0.4e1 + sqrt(0.2e1) * sqrt(R) * pow(M, -0.1e1 / 0.2e1) / 0.2e1;

	return dt_dR;
}

double dt_dE_series(double R,double E,double M,double Lambda)
{
	/**
	Return R.H.S of the t series when it is differentiated w.r.t E
	**/
double	dt_dE = (-0.36465e5 / 0.458752e6 * sqrt(0.2e1) * pow(E, 0.8e1) * pow(M, -0.19e2 / 0.2e1) - 0.143e3 / 0.2048e4 * sqrt(0.2e1) * pow(E, 0.5e1) * Lambda * pow(M, -0.15e2 / 0.2e1) - 0.5e1 / 0.512e3 * sqrt(0.2e1) * E * E * Lambda * Lambda * pow(M, -0.11e2 / 0.2e1)) * pow(R, 0.21e2 / 0.2e1) + (0.6435e4 / 0.77824e5 * sqrt(0.2e1) * pow(E, 0.7e1) * pow(M, -0.17e2 / 0.2e1) + 0.1155e4 / 0.19456e5 * sqrt(0.2e1) * pow(E, 0.4e1) * Lambda * pow(M, -0.13e2 / 0.2e1) + 0.35e2 / 0.7296e4 * sqrt(0.2e1) * E * Lambda * Lambda * pow(M, -0.9e1 / 0.2e1)) * pow(R, 0.19e2 / 0.2e1) + (-0.3003e4 / 0.34816e5 * sqrt(0.2e1) * pow(E, 0.6e1) * pow(M, -0.15e2 / 0.2e1) - 0.105e3 / 0.2176e4 * sqrt(0.2e1) * pow(E, 0.3e1) * Lambda * pow(M, -0.11e2 / 0.2e1) - 0.5e1 / 0.3264e4 * sqrt(0.2e1) * Lambda * Lambda * pow(M, -0.7e1 / 0.2e1)) * pow(R, 0.17e2 / 0.2e1) + (0.231e3 / 0.2560e4 * sqrt(0.2e1) * pow(E, 0.5e1) * pow(M, -0.13e2 / 0.2e1) + 0.7e1 / 0.192e3 * sqrt(0.2e1) * E * E * Lambda * pow(M, -0.9e1 / 0.2e1)) * pow(R, 0.15e2 / 0.2e1) + (-0.315e3 / 0.3328e4 * sqrt(0.2e1) * pow(E, 0.4e1) * pow(M, -0.11e2 / 0.2e1) - 0.5e1 / 0.208e3 * sqrt(0.2e1) * E * Lambda * pow(M, -0.7e1 / 0.2e1)) * pow(R, 0.13e2 / 0.2e1) + (0.35e2 / 0.352e3 * sqrt(0.2e1) * pow(E, 0.3e1) * pow(M, -0.9e1 / 0.2e1) + sqrt(0.2e1) * Lambda * pow(M, -0.5e1 / 0.2e1) / 0.88e2) * pow(R, 0.11e2 / 0.2e1) - 0.5e1 / 0.48e2 * sqrt(0.2e1) * E * E * pow(M, -0.7e1 / 0.2e1) * pow(R, 0.9e1 / 0.2e1) + 0.3e1 / 0.28e2 * sqrt(0.2e1) * E * pow(R, 0.7e1 / 0.2e1) * pow(M, -0.5e1 / 0.2e1) - sqrt(0.2e1) * pow(R, 0.5e1 / 0.2e1) * pow(M, -0.3e1 / 0.2e1) / 0.10e2;

	
	return dt_dE;
}

double dt_dM_series(double R,double E,double M,double Lambda)
{
	/**
	Return R.H.S of the t series when it is differentiated w.r.t M
	**/
double	dt_dM = (0.230945e6 / 0.2752512e7 * sqrt(0.2e1) * pow(E, 0.9e1) * pow(M, -0.21e2 / 0.2e1) + 0.715e3 / 0.8192e4 * sqrt(0.2e1) * pow(E, 0.6e1) * Lambda * pow(M, -0.17e2 / 0.2e1) + 0.55e2 / 0.3072e4 * sqrt(0.2e1) * pow(E, 0.3e1) * Lambda * Lambda * pow(M, -0.13e2 / 0.2e1) + 0.5e1 / 0.20736e5 * sqrt(0.2e1) * pow(Lambda, 0.3e1) * pow(M, -0.9e1 / 0.2e1)) * pow(R, 0.21e2 / 0.2e1) + (-0.109395e6 / 0.1245184e7 * sqrt(0.2e1) * pow(E, 0.8e1) * pow(M, -0.19e2 / 0.2e1) - 0.3003e4 / 0.38912e5 * sqrt(0.2e1) * pow(E, 0.5e1) * Lambda * pow(M, -0.15e2 / 0.2e1) - 0.105e3 / 0.9728e4 * sqrt(0.2e1) * E * E * Lambda * Lambda * pow(M, -0.11e2 / 0.2e1)) * pow(R, 0.19e2 / 0.2e1) + (0.6435e4 / 0.69632e5 * sqrt(0.2e1) * pow(E, 0.7e1) * pow(M, -0.17e2 / 0.2e1) + 0.1155e4 / 0.17408e5 * sqrt(0.2e1) * pow(E, 0.4e1) * Lambda * pow(M, -0.13e2 / 0.2e1) + 0.35e2 / 0.6528e4 * sqrt(0.2e1) * E * Lambda * Lambda * pow(M, -0.9e1 / 0.2e1)) * pow(R, 0.17e2 / 0.2e1) + (-0.1001e4 / 0.10240e5 * sqrt(0.2e1) * pow(E, 0.6e1) * pow(M, -0.15e2 / 0.2e1) - 0.7e1 / 0.128e3 * sqrt(0.2e1) * pow(E, 0.3e1) * Lambda * pow(M, -0.11e2 / 0.2e1) - sqrt(0.2e1) * Lambda * Lambda * pow(M, -0.7e1 / 0.2e1) / 0.576e3) * pow(R, 0.15e2 / 0.2e1) + (0.693e3 / 0.6656e4 * sqrt(0.2e1) * pow(E, 0.5e1) * pow(M, -0.13e2 / 0.2e1) + 0.35e2 / 0.832e3 * sqrt(0.2e1) * E * E * Lambda * pow(M, -0.9e1 / 0.2e1)) * pow(R, 0.13e2 / 0.2e1) + (-0.315e3 / 0.2816e4 * sqrt(0.2e1) * pow(E, 0.4e1) * pow(M, -0.11e2 / 0.2e1) - 0.5e1 / 0.176e3 * sqrt(0.2e1) * E * Lambda * pow(M, -0.7e1 / 0.2e1)) * pow(R, 0.11e2 / 0.2e1) + (0.35e2 / 0.288e3 * sqrt(0.2e1) * pow(E, 0.3e1) * pow(M, -0.9e1 / 0.2e1) + sqrt(0.2e1) * Lambda * pow(M, -0.5e1 / 0.2e1) / 0.72e2) * pow(R, 0.9e1 / 0.2e1) - 0.15e2 / 0.112e3 * sqrt(0.2e1) * E * E * pow(R, 0.7e1 / 0.2e1) * pow(M, -0.7e1 / 0.2e1) + 0.3e1 / 0.20e2 * sqrt(0.2e1) * E * pow(R, 0.5e1 / 0.2e1) * pow(M, -0.5e1 / 0.2e1) - sqrt(0.2e1) * pow(R, 0.3e1 / 0.2e1) * pow(M, -0.3e1 / 0.2e1) / 0.6e1;

	return dt_dM;
}








struct parameters_ { double E;
	                   double dEdr;
	                   double Lambda;
	                   double M;
	                   double dMdr;
	                   double r_loc;};
		
int
func_1 (double t, const double y[], double f[],
      void *p)//params_)
{
	double E, E_r, Lambda, M, M_r;
	//double * parameters = (double *) params_;
	struct parameters_ *params_ = (struct parameters_*) p;
	E = params_->E;
	E_r = params_->dEdr;
	Lambda = params_->Lambda;
	M = params_->M;
	M_r = params_->dMdr;
  
	f[0] = sqrt(2.*E + 2.*M/y[0] + Lambda/3.*y[0]*y[0]);
		
	f[1] = -M/(y[0]*y[0]) + Lambda/3.*y[0];
		
	f[2] = 1./(f[0])* (E_r + M_r/y[0] + f[1]*y[2]);
	
  return GSL_SUCCESS;
}

int
jac_1 (double t, const double y[], double *dfdy, 
     double dfdt[], void *p)//parameters_)
{
  	double E, E_r, Lambda, M, M_r;
	//double * parameters = (double *) parameters_;
	//E = parameters[0];
	//E_r = parameters[1];
	//Lambda = parameters[2];
	//M = parameters[3];
	//M_r = parameters[4];
		struct parameters_ *params_ = (struct parameters_*) p;
	E = params_->E;
	E_r = params_->dEdr;
	Lambda = params_->Lambda;
	M = params_->M;
	M_r = params_->dMdr;
	
	double f0 = sqrt(2.*E + 2.*M/y[0] + Lambda/3.*y[0]*y[0]);
	
	double f1 = -M/(y[0]*y[0]) + Lambda/3.*y[0];
		
	double f2 = 1./(f0)* (E_r + M_r/y[0] + f1*y[2]);
	
	double df0_dR= (-M/(y[0]*y[0])+Lambda*y[0]/3.)/f0; 
	double df1_dR= 2.*M/(y[0]*y[0]*y[0])+Lambda/3.; 
	double df2_dR= (-M_r/(y[0]*y[0]) + df1_dR*y[2]-df0_dR*f2)/f0;

  gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy, 3, 3);
  gsl_matrix * m = &dfdy_mat.matrix; 
  gsl_matrix_set (m, 0, 0, df0_dR);
  gsl_matrix_set (m, 0, 1, 0.);
  gsl_matrix_set (m, 0, 2, 0.);
  
  gsl_matrix_set (m, 1, 0, df1_dR);
  gsl_matrix_set (m, 1, 1, 0.);
  gsl_matrix_set (m, 0, 2, 0.);
  
  gsl_matrix_set (m, 2, 0, df2_dR);
  gsl_matrix_set (m, 2, 1, 0.);
  gsl_matrix_set (m, 2, 2,  f1/f0);
  
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  return GSL_SUCCESS;
}

double bigbang(double R_init, void *p)
{
  double E, E_r, Lambda, M, M_r;
  struct parameters_ *params_ = (struct parameters_*) p;
  E = params_->E;
  E_r = params_->dEdr;
  Lambda = params_->Lambda;
  M = params_->M;
  M_r = params_->dMdr;
  //cout << E <<" "<<" "<< M<<" "<<  Lambda<<" "<<  " param " << endl;
 // cout << "R_init, t_val, t_ls "<< R_init << " "<< t_series(R_init,E,M,Lambda)
 //      << " " << 1e-6 << endl;//1e-4*ageMpc  << endl;

  // time of last scattering is taken to be 100k years then in Mpc it is 1e-4*ageMpc
  double ans = t_series(R_init,E,M,Lambda)-1e-6;//1e-4*ageMpc;
  return ans;
}
void get_t_init(double *t_init, double *R_init , struct parameters_ *p)
{
  
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  double R_low = 0., R_high = 1e-5*p->r_loc;//1e-6;
  gsl_function F;

  F.function = &bigbang;
  F.params = p;
 // cout << " get_t_init "<< p->E << " " << p->M << " " << p->Lambda << endl;

  T = gsl_root_fsolver_brent;//_bisection;//
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, R_low, R_high);
/**
  printf ("using %s method\n", 
          gsl_root_fsolver_name (s));

 printf ("%5s [%9s, %9s] %9s %10s %9s\n",
         "iter", "lower", "upper", "root", 
        "err", "err(est)");
**/
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      R_low = gsl_root_fsolver_x_lower (s);
      R_high = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (R_low, R_high,
                                       0, 0.001);
/**
      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%d [%e, %e] %e  %e\n",
              iter, R_low, R_high,
              r, 
         
             R_high - R_low);
**/
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  //t_init = &r;
  //R_init = &x_lo;
  *R_init = r;

  gsl_root_fsolver_free (s);
}


ScaleFactor::ScaleFactor(struct Params *p, spline_1d *pM, spline_1d *pM_r,
                         spline_1d *pE, spline_1d *pE_r,double t_init=1e-6,
                         double t_max = 30.*ageMpc,int num_pt = 20000)
{
	ScaleFactor::p = p;
	M = pM; M_r = pM_r; E = pE; E_r = pE_r;
	
	ScaleFactor::t_init = t_init;
	ScaleFactor::t_max = t_max;
	size_t_vec = num_pt;
	
	t_vec = new double [num_pt];
    //create t_vec with log distribution
    linspace(t_vec, log10(t_init), log10(t_max), num_pt, true);
    for (int i=0; i<size_t_vec; i++)
    {
    	t_vec[i] = pow(10,t_vec[i]);
    }
}

void ScaleFactor::operator()(struct Bg_sol *sol,int index,double r_loc,
                             double atol=1e-12,double rtol=1e-8)//atol=1e-12,double rtol=1e-10
{
	double E, dEdr, M, dMdr, Lambda;
	E = (*(this->E))(r_loc,0);//(*(ScaleFactor::E))(r_loc);//(*ScaleFactor::E)(r_loc);
	dEdr = (*(this->E_r))(r_loc,0);
	M = (*(this->M))(r_loc,0);
	dMdr = (*(this->M_r))(r_loc,0);
	Lambda = p->Lambda;
	
	struct parameters_ parameters = {E,dEdr, Lambda, M, dMdr,r_loc};
	double t_init;// = t_series(R_init,E,M,Lambda);
	//t_init = 0.;
	double R_init = 0.;
	get_t_init(&t_init,&R_init,&parameters);
	t_init = 1e-6;//1e-4*ageMpc;
	double y_init[3];
	double y_deriv[3];
	y_init[0] = R_init;
	y_init[1] = 1./dt_dR_series(R_init,E,M,Lambda);
	y_init[2] = -(dt_dE_series(R_init,E,M,Lambda)*dEdr + dt_dM_series(R_init,E,M,Lambda)*dMdr)/ 
                dt_dR_series(R_init,E,M,Lambda);
	func_1(t_init,y_init,y_deriv,&parameters); 
	int num_pt = size_t_vec;
	sol->R[(num_pt)*index+0]    = y_init[0];
	sol->R_t[(num_pt)*index+0]  = y_init[1];
	sol->R_r[(num_pt)*index+0]  = y_init[2];
	sol->R_tt[(num_pt)*index+0]  = y_deriv[1];
	sol->R_rt[(num_pt)*index+0]  = y_deriv[2];
     
	//double parameters[5] = {E,dEdr, Lambda, M, dMdr};
	//parameters_ parameters = {E,dEdr, Lambda, M, dMdr};
   //gsl_odeiv2_step_msadams;// bsimp;//;//;//;//rkf45;//rk8pd/rkck;
   //double *ptr = NULL;
  	int (*ptr) (double t, const double y[], double *dfdy, 
                      double dfdt[], void *args);
   ptr= NULL;                   
  gsl_odeiv2_system sys = { func_1, ptr, 3, &parameters};

  gsl_odeiv2_driver *d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
                                   1e-4, rtol, atol);
//without jacobian msadams,rkf45,rkck,rk8pd
//with jacobian msbdf, bsimp, rk4imp, rk2imp, rk1imp
                       

  //double t[num_pt];
  //cout << "t_init " << t_init << " " << t_init/ageMpc << " R_init "<< R_init<< endl;
	
  for (int j = 1; j < num_pt; j++)
    {
      int status = gsl_odeiv2_driver_apply (d, &t_init, t_vec[j], y_init);

      if (status != GSL_SUCCESS)
	{
	  //printf ("error, return value=%d\n", status);
	  cout << "error ScaleFactor, return value= "<< status << endl;
	  break;
	}
	func_1(t_init,y_init,y_deriv,&parameters); 
	sol->R[(num_pt)*index+j]    = y_init[0];
	sol->R_t[(num_pt)*index+j]  = y_init[1];
	sol->R_r[(num_pt)*index+j]  = y_init[2];
	sol->R_tt[(num_pt)*index+j]  = y_deriv[1];
	sol->R_rt[(num_pt)*index+j]  = y_deriv[2];
      //printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
      cout.precision(5);
      //cout.width(10);
      //cout<< setw(4) << t_init<< setw(12) << y_init[0] << setw(12) << y_init[1] << endl;
    }
    gsl_odeiv2_driver_free (d);

}


