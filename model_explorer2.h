#ifndef MODEL_EXPLORER
#define MODEL_EXPLORER
#include <Python.h>
#include "splev_bispeu.h"
#include "sole_params.h"
#include "scipy_tools.h"
#include "tools.h"
#include "mass_function.h"
#include "LTB_background.h"
#include "energy.h"
#include "geodesics.h"
//ModelExplorer
class bg_1dspline_vecs 
{
	private:
		bool allocated;
	public:
		int size_bgr;
		int size_rdense;
		int size_bgt;
		int size_geoz;
		
		double *r_vec ;
		double *M_vec;
		double *Mr_vec;
		double *E_vec;
		double *Er_vec;
		
		double *t_vec; //not allocated set to that from Scalefactor
		
	
	void set_bg_1dspline_vecs(int size_bgr, int size_bgt)
	{
		bg_1dspline_vecs::size_bgr    = size_bgr;
		bg_1dspline_vecs::size_bgt    = size_bgt;
		
		r_vec = new double[size_bgr];			
		M_vec = new double[size_bgr];
		Mr_vec = new double[size_bgr];
		E_vec = new double[size_bgr];
		Er_vec = new double[size_bgr];	
		
		allocated =true;
	}
	

  ~bg_1dspline_vecs()
  {
		if (allocated)
		{
			delete [] r_vec;
			delete [] M_vec;
			delete [] Mr_vec;
			delete [] E_vec;
			delete [] Er_vec;
		}
  }
};

class ModelExplorer
{
	private:
	bool bg_allocated;
	bool geo_allocated;
	void release_mem(void)
	{
		if (bg_allocated == true)
		{
			delete spvecs1d;
			delete spM;
			delete spMr;
			delete spE;
			delete spEr;
			delete spR;
			delete spR_r;
			delete spR_rt;
			delete spR_t;
			delete SF;
			delete solbg;
			
			bg_allocated = false;
		}
		if (geo_allocated == true)
		{
			delete LTBGeo;
			delete solgeo;
			
			geo_allocated = false;
		}
	}
	public:
	Params *p;
	
	int size_bgr;
	int size_bgt;
	double r_min;
	double r_max;
	double t_init;
	double t_max;
	//int (*dMdr_gsl)(double, const double, double, void);
	gsl_deriv deriv_MassFunc;
	spline_1d * spM;
	spline_1d * spMr;
	spline_1d * spE;
	spline_1d * spEr;
	scipy_spline_2d *spR;
	scipy_spline_2d *spR_r;
	scipy_spline_2d *spR_rt;
	scipy_spline_2d *spR_t;	
	bg_1dspline_vecs *spvecs1d;
	Bg_sol *solbg;
	ScaleFactor *SF;
	
	int size_geoz;
	double z_min;
	double z_max;
	Geo_sol *solgeo;
	LTB_Geodesics *LTBGeo;
	
	void initialize_bg(Params *,
                       int (*deriv_MassFunc_)(double, const double *, double*, void*) ,
                       int , int, double, double,double, double);
    void initialize_bg(Params *p_,
                  int (*deriv_MassFunc_)(double, const double*, double*, void*));
    void initialize_geo(double , double, int);
	void evolve_geo(double *, double *, double, double);
	void reinitialize_bg(Params *, int, int, double,double, double, double);
	void reset_bg_and_geo(void)
	{
		release_mem();
	}
	~ModelExplorer()
	{
		release_mem();
	}
};

struct type_MCMC_param
{
	double lower;
	double upper;
	double value;
	double sigma;
};

struct MCMC_params
{
	type_MCMC_param r_loc; //location of observer
	type_MCMC_param r0; //void size
	type_MCMC_param delta_w; //void steepness
	// declare derived parameters as mcmc parameters, they can be provided 
	// which can be used a priors
	type_MCMC_param dA; //angular diameter distance to last scattering surface
	type_MCMC_param DM; //dipole magnitude
	type_MCMC_param zdec; //last scattering redshift
	type_MCMC_param tdec;
	type_MCMC_param rdec;

	bool reset_bg;

	MCMC_params(bool reset_bg = false)
	{
		MCMC_params::reset_bg = reset_bg;
	}
};

class MCMC
{
	public:
	ModelExplorer *explorer;
	Params *p;
	MCMC_params *mu;
	MCMC(Params *p, ModelExplorer *explorer, MCMC_params *mu)
	{
		MCMC::p = p;
		MCMC::explorer = explorer;
		MCMC::mu = mu;
	}
	
	double get_likelihood(MCMC_params *mu);
	void set_derived_parameters(MCMC_params *, MCMC_params *);

};

#endif












