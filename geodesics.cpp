#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <typeinfo>
#include "geodesics.h"
using namespace std;


struct func_wrapper
{
	spline_1d *E;
	spline_1d *E_r;
	scipy_spline_2d *R;
	scipy_spline_2d *R_r;
	scipy_spline_2d *R_rt;
	scipy_spline_2d *R_t;
};

void set_LTB_init_conds(double *P_obs,double *Dir, double *y_init, struct func_wrapper *p)
{
	double t, r, theta, phi;
		t = P_obs[0]; r = P_obs[1]; theta = P_obs[2]; phi = P_obs[3];
	double a, b;
		a = Dir[0]; b = Dir[1];
		y_init[0] = t; 
		y_init[1] = r;
		y_init[2] = theta;
		y_init[3] = phi; 
		
		b = b-phi;
		// y_init[4] = dr_ds; y_init[5] = dtheta_ds;  y_init[6] = dphi_ds
		//std::cout << typeid(a).name() << '\n';
		double E;
		E = (*(p->E))(r,0);
		y_init[4] = -sin(a)*cos(b)*sqrt(1.+2.*E);
		y_init[5] = sin(a)*sin(b)/r;
		y_init[6] = cos(a)/(sin(theta)*r);
		//y_init[7] = affine_parameter = 0. when z=0.
		y_init[7] = 0.;
		
		//cout << "norm  " << y_init[4]*y_init[4]/sqrt(1.+2.*E)+
		//        r*r*(y_init[5]*y_init[5]+ sin(theta)*sin(theta)*y_init[6]*y_init[6])
		//        << " r " << r << endl;
}
		
int LTB_geodesic_derivs (double z, const double y[], double f[],
      void *p_)
{
	double E, E_r, R, R_r, R_rr, R_rt, R_t;
	struct func_wrapper *p = (struct func_wrapper*) p_;
  	
  	double t, r, theta, phi;
		t=y[0]; r=y[1]; theta=y[2]; phi=y[3];
	double sin_theta, cos_theta;
		sin_theta = sin(theta); cos_theta = cos(theta);
	
	double r_pos = abs(r);
	if ( r_pos > 20*1e3)
	{
		cout << "Gentlmen, we got a problem. Change r_pos limits" << endl;
		cout << "r_pos " << r_pos << " rmax "<< 20*1e3 << endl;
	}
	E    = (*(p->E))(r_pos,0);
	E_r  = (*(p->E))(r_pos,0);
	R    = (*(p->R))(r_pos,t,0,0);
	R_r  = (*(p->R_r))(r_pos,t,0,0);
	R_rr = (*(p->R_r))(r_pos,t,1,0);
	R_rt = (*(p->R_rt))(r_pos,t,0,0);
	R_t  = (*(p->R_t))(r_pos,t,0,0);

	double dt_ds = (1.+z);
	double dr_ds = y[4];
	double dtheta_ds = y[5];
	double dphi_ds = y[6];
	//double DA = y[7];
	
	double ddt_dss, ddr_dss, ds_dz, ddr_dsz, ddtheta_dss, ddtheta_dsz;
	double ddphi_dss, ddphi_dsz, dt_dz, dr_dz, dtheta_dz, dphi_dz;
	
	ddt_dss = -R_r*R_rt/(1.+2.*E)*dr_ds*dr_ds - R*R_t*( dtheta_ds*dtheta_ds + 
	                                      sin_theta*sin_theta*dphi_ds*dphi_ds );
	
	ds_dz = 1./ddt_dss;
	
	ddr_dss = -(2.*R_rr*E-R_r*E_r+R_rr)/(1.+2.*E)/R_r*dr_ds*dr_ds 
	          - 2.*R_rt/R_r*dr_ds*dt_ds + 
	          (1.+2.*E)*R/R_r*(dtheta_ds*dtheta_ds + sin_theta*sin_theta*dphi_ds*dphi_ds);
	ddr_dsz = ddr_dss*ds_dz;
	
	ddtheta_dss = -2./R*dtheta_ds*(R_r*dr_ds + R_t*dt_ds) + sin_theta*cos_theta*dphi_ds*dphi_ds;
	ddtheta_dsz = ddtheta_dss*ds_dz;
	
	ddphi_dss = -2./R*dphi_ds*(R_r*dr_ds + R_t*dt_ds) - 2.*cos_theta/sin_theta*dtheta_ds*dphi_ds;
	ddphi_dsz = ddphi_dss*ds_dz;
	
	//dDA_ds = ;
	//dDA_dz = -dDA_ds*ds_dz
	
	dt_dz = dt_ds*ds_dz;
	dr_dz = dr_ds*ds_dz;
	dtheta_dz = dtheta_ds*ds_dz;
	dphi_dz = dphi_ds*ds_dz;
	
	f[0] = dt_dz; f[1] = dr_dz; f[2] = dtheta_dz; f[3] = dphi_dz;
	f[4] = ddr_dsz; f[5] = ddtheta_dsz; f[6] = ddphi_dsz; 
	f[7] = ds_dz; //f[8] = dDA_dz;
	
  return GSL_SUCCESS;
}


LTB_Geodesics::LTB_Geodesics(struct Params *Bg_param_, 
                         spline_1d *pE, spline_1d *pE_r,
                         scipy_spline_2d* pR,scipy_spline_2d* pR_r,
                         scipy_spline_2d* pR_rt,scipy_spline_2d* pR_t,
                         double z_max = 3000., double z_init = 1e-6, int num_pt = 1600)
{
	Bg_param = Bg_param_;
	E    = pE; 
	E_r  = pE_r;
	R    = pR; 
	R_r  = pR_r;
	R_rt = pR_rt;
	R_t  = pR_t;
	
	LTB_Geodesics::z_max = z_max;
	LTB_Geodesics::z_init = z_init;
	LTB_Geodesics::size_z_vec = num_pt;
	
	z_vec = new double [num_pt];
	z_vec[0] = 0.;
    //create t_vec with log distribution
   /** linspace(&z_vec[1], log10(z_init), log10(z_max), size_z_vec-1, true);
    for (int i=1; i<size_z_vec; i++)
    {
    	z_vec[i] = pow(10,z_vec[i]);
    }
    **/
    
        linspace(z_vec,0., 1e-6,100,false);
        linspace(&z_vec[100],1e-6, 1.,100,false);
        linspace(&z_vec[200],1.,10,200,false);
        linspace(&z_vec[400],10,z_max, size_z_vec-400, true);
     
        // if z_vec is increasing too slowly then t will not be strictly decreasing
        //but if it is not fine enough then r increasing too rapidly as the ode solver fails
        //linspace(z_vec,0.,z_max, size_z_vec, true);
 
}

void LTB_Geodesics::operator()(struct Geo_sol *sol,double *P_obs,double * Dir,
                             double atol=1e-12,double rtol=1e-8)
{
	/**
	P_obs: Array identifying the position of the observer
		   [t_obs, r_obs, theta_obs, phi_obs)
	Dir: Array of angular direction in which the geodesic propagates.
		 [theta_star , phi_star]
	**/
	sol->z_vec = z_vec;
	
	struct func_wrapper parameters;
	parameters.E    = E; 
	parameters.E_r  = E_r;
	parameters.R    = R; 
	parameters.R_r  = R_r;
	parameters.R_rt = R_rt;
	parameters.R_t  = R_t;
	
	double z_init = 0.;
	double y_init[8];
	double y_deriv[8];	
	set_LTB_init_conds(P_obs, Dir, y_init, &parameters);


	sol->t[0]  = y_init[0];
	sol->r[0]  = y_init[1];
	sol->theta[0]  = y_init[2];
	sol->phi[0]  = y_init[3];
	sol->dr_ds[0]  = y_init[4];
	sol->dtheta_ds[0] = y_init[5];
	sol->dphi_ds[0] = y_init[6];
	sol->s[0] = y_init[7];
    
    LTB_geodesic_derivs(0, y_init, y_deriv,&parameters);
    sol->dt_dz[0] = y_deriv[0];
    sol->ddr_dsz[0]= y_deriv[4];
    sol->ddtheta_dsz[0]= y_deriv[5];
    sol->ddphi_dsz[0]= y_deriv[6];
    sol->ds_dz[0]= y_deriv[7];
    
    sol->dt_ds[0] = 1.; 
                
 	gsl_odeiv2_system sys = { LTB_geodesic_derivs, NULL, 8, &parameters};

  	gsl_odeiv2_driver *d =
                     gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                                    z_vec[1], rtol, atol);
	//without jacobian msadams,rkf45,rkck,rk8pd
	//with jacobian msbdf, bsimp, rk4imp, rk2imp, rk1imp
	
  for (int j = 1; j < size_z_vec; j++)
    {
      int status = gsl_odeiv2_driver_apply (d, &z_init, z_vec[j], y_init);

      if (status != GSL_SUCCESS)
	{
		cout <<"Geodesics " << Bg_param->r0 << " " 
		<< Bg_param->delta_w << " "
	    << Bg_param->H_out << endl;
	  cout << "error Geodesics, return value= "<< status << endl;
	  break;
	}
	LTB_geodesic_derivs(z_init, y_init, y_deriv,&parameters);
	sol->t[j]  = y_init[0];
	sol->r[j]  = y_init[1];
	sol->theta[j]  = y_init[2];
	sol->phi[j]  = y_init[3];
	sol->dr_ds[j]  = y_init[4];
	sol->dtheta_ds[j] = y_init[5];
	sol->dphi_ds[j] = y_init[6];
	sol->s[j] = y_init[7];
    
    sol->dt_dz[j]  = y_deriv[0];
    sol->ddr_dsz[j]= y_deriv[4];
    sol->ddtheta_dsz[j]= y_deriv[5];
    sol->ddphi_dsz[j]= y_deriv[6];
    sol->ds_dz[j]= y_deriv[7]; 
 	sol->dt_ds[j] = 1.+z_init;
    }
    
    gsl_odeiv2_driver_free (d);

}


/**
int Szekeres_geodesic_derivs (double t, const double y[], double f[],
      void *p)
{
	
  return GSL_SUCCESS;
}
**/
