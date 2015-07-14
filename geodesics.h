#ifndef GEODESICS
#define GEODESICS
#include <iostream>
#include "sole_params.h"
#include "tools.h"

struct Geo_sol
{
	// contains arrays of all the geodesics and their derivatives
	double *t;
	double *r;
	double *theta;
	double *phi;
	double *s; // affine parameter
	double *dt_dz;
	double *dr_ds;
	double *dtheta_ds;
	double *dphi_ds;
	double *ds_dz;
	double *ddr_dsz;
	double *ddtheta_dsz;
	double *ddphi_dsz;
	double *dt_ds;
	//double *DA; angular diameter distance
	
	int size_z_vec;
	double *z_vec; //set to z_vec in LTB_Geodesics function operator
	
	Geo_sol(int size_z_vec)
	{
		Geo_sol::size_z_vec = size_z_vec;
		t           = new double [size_z_vec];
		r           = new double [size_z_vec];
		theta       = new double [size_z_vec];
		phi         = new double [size_z_vec];
		s           = new double [size_z_vec];
		dt_dz       = new double [size_z_vec];
		dr_ds       = new double [size_z_vec];
		dtheta_ds   = new double [size_z_vec];
		dphi_ds     = new double [size_z_vec];
		ds_dz       = new double [size_z_vec];
		ddr_dsz     = new double [size_z_vec];
		ddtheta_dsz = new double [size_z_vec];
		ddphi_dsz   = new double [size_z_vec];
		dt_ds       = new double [size_z_vec]; 
		
	};
	~Geo_sol()
	{
		delete [] t;
		delete [] r;
		delete [] theta;
		delete [] phi;
		delete [] s;
		delete [] dt_dz;
		delete [] dr_ds;
		delete [] dtheta_ds;
		delete [] dphi_ds;
		delete [] ds_dz;
		delete [] ddr_dsz;
		delete [] ddtheta_dsz;
		delete [] ddphi_dsz;
		delete [] dt_ds;
		//delete [] DA;
	
	}
	
};

class LTB_Geodesics
{
	/**
	Solves the Null geodesics in LTB model. The independent variable is
	redshift and not the affine parameter. The most common approach is to set 
	phi to zero and remove its differential equation from the system of PDEs 
	describing the LTB model. This is possible because of the inherent symmetry in 
	the LTB metric. Here I have not made any simplifications and both theta and 
	phi coordiantes are evolved.
	**/
	public:
		struct Params *Bg_param;
		spline_1d *E, *E_r;
		scipy_spline_2d *R, *R_r, *R_rt, *R_t;
		
		double *z_vec;
		int size_z_vec;
		double z_max;
		double z_init;
		
	    
	    void operator()(struct Geo_sol *,double *,double *, double , double );
	    
	    LTB_Geodesics(struct Params *, 
                         spline_1d *, spline_1d *,
                         scipy_spline_2d*, scipy_spline_2d*,
                         scipy_spline_2d*, scipy_spline_2d*,
                         double, double, int );
	   
       ~LTB_Geodesics()
       {
        	delete [] z_vec;
       }
	    
};


#endif
