#include <assert.h>
#include "model_explorer2.h"

void ModelExplorer::initialize_bg(Params *p_,
                                  int (*deriv_MassFunc_)(double, const double*, double*, void*),
                                  int size_bgr_ = 200, 
                                  int size_bgt_=2000, double r_min_ = 1e-4,
                                  double r_max_ = 20*1e3,
                                  double t_init_ = 1e-6, double t_max_ = 30.*ageMpc)
{
	size_bgr = size_bgr_;
	size_bgt = size_bgt_;
	p = p_;
	deriv_MassFunc = deriv_MassFunc_;
	r_min = r_min_;
	r_max = r_max_;
	t_init = t_init_;
	t_max = t_max;
	
	assert( (size_bgr_ % 2 == 0) && "size_bgr must be divisible by 2");
	spvecs1d = new bg_1dspline_vecs;
	spvecs1d->set_bg_1dspline_vecs(size_bgr, size_bgt);
  	sample_radial_coord(spvecs1d->r_vec, p, r_min, r_max, size_bgr/2, size_bgr/2);
	double *r_vec = spvecs1d->r_vec;
	
	set_M_Mr_vec(spvecs1d->M_vec, spvecs1d->Mr_vec, r_vec, size_bgr, 
	             p, deriv_MassFunc);
	spM = new spline_1d(r_vec, spvecs1d->M_vec, size_bgr);
	spMr= new spline_1d(r_vec, spvecs1d->Mr_vec,size_bgr);

	set_E_Er_vec(spvecs1d->E_vec,spvecs1d->Er_vec, r_vec,size_bgr, spM, p);
	spE= new spline_1d(r_vec, spvecs1d->E_vec, size_bgr);
	spEr= new spline_1d(r_vec, spvecs1d->Er_vec,size_bgr);
	
	SF = new ScaleFactor(p,spM,spMr,spE,spEr,
                                 t_init,p->age,size_bgt);
	solbg = new Bg_sol(size_bgt,size_bgr);
	#pragma omp parallel for
	for (int i = 0; i<size_bgr; i++)
	{
  		(*SF)(solbg, i, r_vec[i],1e-15,1e-12);//1e-12,1e-10);
	}
	
	double *t_vec = SF->t_vec;
	spR = new scipy_spline_2d(r_vec, t_vec,solbg->R,size_bgr,size_bgt,0);
	spR_r = new scipy_spline_2d(r_vec, t_vec,solbg->R_r,size_bgr,size_bgt,0);
	spR_rt= new scipy_spline_2d(r_vec, t_vec,solbg->R_rt,size_bgr,size_bgt,0);
	spR_t= new scipy_spline_2d(r_vec, t_vec,solbg->R_t,size_bgr,size_bgt,0);
	//scipy_spline_2d spR_tt(r_vec, t_vec,sol->R_tt,size_bgr,size_bgt,0);
	bg_allocated = true;
		
}

void ModelExplorer::initialize_bg(Params *p_,
                                  int (*deriv_MassFunc_)(double, const double*, double*, void*))
{
initialize_bg(p_, deriv_MassFunc_, 200, 2000, 1e-4,20*1e3,1e-6, 30.*ageMpc);
}

void ModelExplorer::initialize_geo(double z_max_= 1200., double z_min_=0.25,
                                   int size_geoz_ = 1600)
{
	size_geoz = size_geoz_;
	z_min = z_min_;
	z_max = z_max_;
	LTBGeo = new LTB_Geodesics(p,spE,spEr,spR,spR_r,spR_rt,spR_t,
                            z_max,z_min,size_geoz);//1e-4,1600);
    solgeo = new Geo_sol(size_geoz);
    geo_allocated = true;
    //std::cout << "ModelExplorer::initialize_geo " << z_min << " " << z_max << std::endl;
}

void ModelExplorer::evolve_geo(double *P_obs, double *Dir,
                               double atol = 1e-12, double rtol = 1e-10)
{
	if (geo_allocated == true)
	{
		(*LTBGeo)(solgeo,P_obs,Dir,atol,rtol);
	}
	else
	{
		std::cout << "You got the blues!. geo_allocated = false" << std::endl;
	}
}

void ModelExplorer::reinitialize_bg(Params *p_,
                                  int size_bgr_ = 200, 
                                  int size_bgt_=2000, double r_min_ = 1e-4,
                                  double r_max_ = 20*1e3,
                                  double t_init_ = 1e-6, double t_max_ = 30.*ageMpc)
{
	assert(bg_allocated == false && 
	       "You got the blues!!. Destroy current background instance before making a new one"); 
	initialize_bg(p_,deriv_MassFunc,size_bgr_,size_bgt_, r_min_, r_max_, 
	              t_init_,t_max_);
}









