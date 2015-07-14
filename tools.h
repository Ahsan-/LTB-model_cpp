#ifndef TOOLS
#define TOOLS
#include <iostream>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <stdexcept>
#include <iomanip>
#include <Python.h>
#include "sole_params.h"
#include "splev_bispeu.h"
#include "scipy_tools.h"


class spline_1d
{	
	public:
		gsl_interp *interp;
		const gsl_interp_type *T;
		size_t size;
		const double *xa;
		const double *ya;
		//double *xa;
		//double *ya;
		gsl_interp_accel * acc;
		

	spline_1d(const double *xa_, const double *ya_, int size_)
	//spline_1d(double *xa_, double *ya_, int size_)
	{
		size = (size_t) size_;
		xa = xa_;
		ya = ya_;
		T = gsl_interp_cspline;
		acc = gsl_interp_accel_alloc(); 
		interp = gsl_interp_alloc(T, size);
		gsl_interp_init(interp,xa, ya,size);
		
		
	}
			
	double operator()(double x,int nu)
	{	
		switch(nu){
			case 0:
				return gsl_interp_eval(interp,xa,ya,x,acc);
			case 1:
				return gsl_interp_eval_deriv(interp,xa,ya,x,acc);
			case 2:
				return gsl_interp_eval_deriv2(interp,xa,ya,x,acc);
			default:
				throw std::invalid_argument( "nu can only take on the values {0,1,2}" );
		}
	}
	double operator()(double a, double b)
	{	
		
		return gsl_interp_eval_integ(interp,xa,ya,a,b,acc);
		
	}
	~spline_1d()
	{
		gsl_interp_accel_free(acc);
		gsl_interp_free(interp);
		//std::cout << "yo free as can be" << std::endl;
	}
};


class scipy_spline_1d
{	// set_1d_spline(char *fun_name, double *x_, double *y_, int size,int smoothing)
	public:
    struct spline_1d_param sp_p;
    struct spline_1d_param dsp_p;
	scipy_spline_1d(double *x, double *y,int size_x, int smoothing=0)
	{
		set_1d_spline(&sp_p, &dsp_p, &x[0],&y[0],size_x,smoothing);
		
	}
	double operator()(double x, int nu)
	{
		int ier =0;
		int size_x_val = 1;
		double x_val[1] = {x};
		double y_val[1];
		double ans=0.;
		
		switch(nu){
			case 0:
				splev_cpp_(sp_p.t,&sp_p.n,sp_p.c,&sp_p.k,&x_val[0],&y_val[0],
				           &size_x_val,&ier);
				ans = y_val[0];
				break;
			case 1:
				splev_cpp_(dsp_p.t,&dsp_p.n,dsp_p.c,&dsp_p.k,&x_val[0],
				           &y_val[0],&size_x_val,&ier);
				ans = y_val[0];
				break;
			default:
				throw std::invalid_argument( "nu can only take on the values {0,1}" );
			
			if ( ier !=0 )
			{
				throw std::invalid_argument( "Spline evaluation failed, possibly "
				"the tried x value is out of range" );
			}
		 }
		 
		return ans;
	  }
	~scipy_spline_1d()
	{
		//nothing for now
	}
};

class scipy_spline_2d
{	// set_1d_spline(char *fun_name, double *x_, double *y_, int size,int smoothing)
	public:
    struct spline_2d_param sp;
    struct spline_2d_param dspx;
    struct spline_2d_param dspy;
	scipy_spline_2d(double *x, double *y,double *z,int size_x, int size_y, int smoothing=0)
	{
		set_2d_spline(&sp, &dspx,&dspy, &x[0],&y[0],&z[0],
		              size_x,size_y,smoothing);
		
	}
	double operator()(double x,double y, int dx, int dy)
	{
		int ier =0;
		int size_x_val = 1;
		double x_val[1] = {x};
		double y_val[1] = {y};
		double z_val[1]; 
		double ans=0.;
		
		if ( dx==0 && dy ==0)
		{
			double wrk[sp.lwrk];
			bispeu_cpp_(sp.tx,&sp.nx,sp.ty,&sp.ny,sp.c,&sp.nc,&sp.kx,&sp.ky,
                         &x_val[0],&y_val[0],&z_val[0],&size_x_val,&wrk[0],&sp.lwrk,&ier);
				ans = z_val[0];
		}
		else if (dx==1 && dy ==0)
		{
			double wrk[dspx.lwrk];
			bispeu_cpp_(dspx.tx,&dspx.nx,dspx.ty,&dspx.ny,dspx.c,&dspx.nc,
			             &dspx.kx,&dspx.ky,
                         &x_val[0],&y_val[0],&z_val[0],&size_x_val,&wrk[0],&dspx.lwrk,&ier);
				ans = z_val[0];
		}
		else if (dx==0 && dy ==1)
		{
			double wrk[dspy.lwrk];
			bispeu_cpp_(dspy.tx,&dspy.nx,dspy.ty,&dspy.ny,dspy.c,&dspy.nc,
			             &dspy.kx,&dspy.ky,
                         &x_val[0],&y_val[0],&z_val[0],&size_x_val,&wrk[0],&dspy.lwrk,&ier);
				ans = z_val[0];
		}
		else
		{
				throw std::invalid_argument( "dx and dy can only take on the values {0,1}" );
		}
			
		if ( ier !=0 )
			{
				throw std::invalid_argument( "Spline evaluation failed, possibly "
				"the tried x or y value is out of range" );
			}
		 
		return ans;
	  }
	~scipy_spline_2d()
	{
		//nothing for now
	}
};

void linspace(double *t, double t_init, double t_end, int num_pt,
              bool end_point);

void sample_radial_coord(double *r_vec, struct Params *p, double r_min,
                         double r_max, int num_pt1, int num_pt2);
#endif


















