#ifndef LTB_BACKGROUND
#define LTB_BACKGROUND
#include <iostream>
#include "sole_params.h"
#include "tools.h"

struct Bg_sol
{
	// contains arrays of all the background solutions
	double *R;
	double *R_r;
	double *R_rr;
	double *R_rt;
	double *R_t;
	double *R_tt;
	
	int size_t_vec;
	int size_r_vec;
	
	//Bg_sol(){};
	
	Bg_sol(int size_t_vec,int size_r_vec)
	{
		Bg_sol::size_t_vec = size_t_vec;
		Bg_sol::size_r_vec = size_r_vec;
		
		R = new double [size_t_vec*size_r_vec];
		R_r = new double [size_t_vec*size_r_vec];
		R_rr = new double [size_t_vec*size_r_vec];
		R_rt = new double [size_t_vec*size_r_vec];
		R_t = new double [size_t_vec*size_r_vec];
		R_tt = new double [size_t_vec*size_r_vec];
	};
	~Bg_sol()
	{
		delete [] R;
		delete [] R_r;
		delete [] R_rr;
		delete [] R_rt;
		delete [] R_t;
		delete [] R_tt;
		//std::cout << "Bg_sol freed" << std::endl;
	}
	
};

class ScaleFactor
{
	/**
	A class for solving Eq. (2.2) in ``Structures in the Universe by Exact Method``
	by Krzystof Bolejko etal. The user provides Lambda, E(r), M(r), diff(E(r),r)
	and diff(M(r),r). Once solved they can be used in the geodesic equations for 
	both the LTB and Szekeres models.
	**/
	public:
		struct Params *p;
		spline_1d *M, *M_r, *E, *E_r;
		
		double *t_vec;
		int size_t_vec;
		double t_init;
		double t_max;
		
	    
	    void operator()(struct Bg_sol *,int,double,double,double);
	    
	    ScaleFactor(struct Params *, spline_1d *, spline_1d *, spline_1d *,spline_1d *,
	                double , double , int);
	   
	   int (*func) (double t, const double y[], double f[], void *args);
	   int (*jac) (double t, const double y[], double *dfdy, 
                      double dfdt[], void *args);
       ~ScaleFactor()
       {
        	delete [] t_vec;
       }
	    
};


#endif
