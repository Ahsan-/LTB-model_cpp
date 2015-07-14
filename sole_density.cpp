#include <iomanip>
#include <omp.h>
#include <typeinfo>
#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <Python.h>
#include "splev_bispeu.h"
#include "sole_params.h"
#include "scipy_tools.h"
#include "tools.h"
#include "mass_function.h"
#include "LTB_background.h"
#include "energy.h"
#include "geodesics.h"
#include "model_explorer2.h"
#include "cmb.h"


using namespace std;


void perform_MCMC(void);

int main()
{
	Py_Initialize();
    initscipy_tools();
    initemcee_wrapper();

	Params WILT_params(0.73, 0.673,0.315, 53./0.673, 15./0.673, -0.99);
	//r_loc = 75.22 , r0=78.87 -0.98072 0.673
	//Params WILT_params(0.73, 0.673,0.315, 78.87, 15./0.673, -0.98072);
	ModelExplorer universe;
	universe.initialize_bg(&WILT_params,dMdr_gsl);
	universe.initialize_geo(1500.,0.,1600);
	//check that all went well
	/**************************************************************************/
	double P_obs[] = {WILT_params.age, 27.9/WILT_params.H_out, 
	                  PI/2.+0.*(90.+29.3)*PI/180., 276.4*PI/180.};
	//P_obs[1] = 75.22;//1e-3;
	double Dir[] = {PI/2., 276.4*PI/180.};

	universe.evolve_geo(P_obs,Dir,1e-15,1e-12);//1e-12,1e-10);
	
	
	ofstream f("geodesics_sol.dat");
	f << "z t r theta phi" << endl;
	for (int i=0; i<universe.size_geoz; i++)
	{
		f << boost::format("%.6e %.6e %.6e %.6e %.6e\n")
		%universe.solgeo->z_vec[i] %universe.solgeo->t[i] %universe.solgeo->r[i] 
		%universe.solgeo->theta[i] %universe.solgeo->phi[i];
	}
	f.close();

	/**************************************************************************/
	//Params WILT_params1(0.5, 0.673,0.315, 23./0.673, 15./0.673, -0.96);
	Params WILT_params1(0.73, 0.673,0.315, 78.87, 15./0.673, -0.98072);

	universe.reset_bg_and_geo();
	universe.reinitialize_bg(&WILT_params1,
                             200, 2000, 1e-4,20*1e3,1e-6, 30.*ageMpc);
                            
    universe.initialize_geo(1200.,1e-6,1600);                        
	P_obs[1] = 1e-3;//44./WILT_params.H_out;

	universe.evolve_geo(P_obs,Dir,1e-15,1e-12);
	
	ofstream f2("geodesics_sol.dat",std::ios_base::app);
	f2 << "z t r theta phi" << endl;
	for (int i=0; i<universe.size_geoz; i++)
	{
		f2 << boost::format("%.6e %.6e %.6e %.6e %.6e\n")
		%universe.solgeo->z_vec[i] %universe.solgeo->t[i] %universe.solgeo->r[i] 
		%universe.solgeo->theta[i] %universe.solgeo->phi[i];
	}
	f2.close();

	/**
	ofstream f3("Helmholtz.dat");
	f3 << "comoving_r       K(r)       \\dot{R}     K>(\\dot{R}^2 - 1)/2" << endl;
	double age = WILT_params.age;
	for (int i=0; i<universe.size_bgr; i++)
	{
		double rval = universe.spvecs1d->r_vec[i];
		double K = (*(universe.spE))(rval,0);
		double Rdot = (*(universe.spR_t))(rval,age,0,0); 
		
		f3 << boost::format("%.6e  %.6e   %.6e   %d\n")
		%rval   %K   %Rdot  %(K>(Rdot*Rdot - 1.)/2.);
	}
	f3.close();
	**/
	perform_MCMC();

  Py_Finalize(); 
	return 0;
}





double MCMC::get_likelihood(MCMC_params *mu_)
{
	double lnlkh;
	//rmin = 1e-4
	double P_obs[] = {p->age, 1e-3, PI/2.+0.*(90.+29.3)*PI/180., 276.4*PI/180.};
	double Dir[] = {PI/2., 276.4*PI/180.};
	if ( mu_->reset_bg == true)
	{
		//cout << "reinitalized bg" << endl;
		p->r0 = mu_->r0.value;
		p->delta_w = mu_->delta_w.value;
		explorer->reset_bg_and_geo();
		explorer->reinitialize_bg(p,200, 2000, 1e-4,20*1e3,1e-6, 30.*ageMpc);
    	explorer->initialize_geo(1200.,0.,1600);//1e-6,1600);   
		mu_->reset_bg = false;
	}
	
	MCMC_params mu_central_1;
    //fill ltbSol for the central observer
	explorer->evolve_geo(P_obs,Dir,1e-15,1e-12);
	set_derived_parameters(&mu_central_1, NULL);
	//fill ltbSol for the current position of observer
	P_obs[1] = mu_->r_loc.value;
	explorer->evolve_geo(P_obs,Dir,1e-15,1e-12);
	set_derived_parameters(mu_, &mu_central_1);
	// - lnL = chi^2/2
	double sigma_DM = 0.9*1e3/c_;
	lnlkh = -0.5*(mu_->DM.value - 1.23e-3)*(mu_->DM.value - 1.23e-3)/ (sigma_DM*sigma_DM); 
	return lnlkh;
	//http://arxiv.org/pdf/1303.5087v2.pdf
	//http://users-phys.au.dk/haugboel/pdf/0802.1523.pdf
}

void MCMC::set_derived_parameters(MCMC_params *mu_, MCMC_params *mu_central)
{
	if (mu_central == NULL)
	{// assume that calculations are for the central observer
		spline_1d spt(explorer->solgeo->z_vec, explorer->solgeo->t,
		explorer->solgeo->size_z_vec);
		mu_->zdec.value = 1100.;
		mu_->tdec.value = spt(mu_->zdec.value,0);
		//cout << " mu_ " << mu_->tdec.value << endl;
		spline_1d spr(explorer->solgeo->z_vec, explorer->solgeo->r,
		explorer->solgeo->size_z_vec);
		mu_->rdec.value = abs(spr(mu_->zdec.value,0));
		mu_->dA.value = (*(explorer->spR))(mu_->rdec.value,mu_->tdec.value,0,0);
		mu_->DM.value = 0.; //DM can not be set.
	}
	else
	{// assume mu_central is filled and can be used for off center calculations
		mu_->tdec.value = mu_central->tdec.value;
		double t_neg_temp[explorer->size_geoz];
		for (int i = 0; i < explorer->size_geoz; i++)
		{
			t_neg_temp[i] = - explorer->solgeo->t[i];
		}
 
		//for (int i=0; i< explorer->size_geoz-1; i++)
		//{
			//assert(explorer->solgeo->t[i+1]- explorer->solgeo->t[i] < 0 && 
			//      "time not strictly decreasing");
			//cout << "time not strictly decreasing" << endl;
		//}
		spline_1d spz(&t_neg_temp[200], &explorer->solgeo->z_vec[200],
		explorer->solgeo->size_z_vec-200);
		mu_->zdec.value = spz(-mu_->tdec.value,0);
		spline_1d spr(explorer->solgeo->z_vec, explorer->solgeo->r,
		explorer->solgeo->size_z_vec);
		mu_->rdec.value = abs(spr(mu_->zdec.value,0));
		mu_->dA.value = (*(explorer->spR))(mu_->rdec.value,mu_->tdec.value,0,0);
		mu_->DM.value = (mu_central->zdec.value - mu_->zdec.value)/
		(1. + mu_central->zdec.value);
	}


}

void perform_MCMC(void)
{
	Params WILT_params(0.73, 0.673,0.315, 53./0.673, 15./0.673, -0.99);
	//MCMC_params communicated to emcee
	MCMC_params mu;
	
	//position of observer
	mu.r_loc.lower = 1e-2*Mpc;
	mu.r_loc.upper = 500*Mpc;//2e3*Mpc;
	mu.r_loc.value = 50./WILT_params.H_out;//27.9/WILT_params.H_out;
	mu.r_loc.sigma = 1.*Mpc;//0.1*Mpc;//0.1*Mpc;
	//void size
	mu.r0.lower = 50.*Mpc/WILT_params.H_out;
	mu.r0.upper = 150.*Mpc/WILT_params.H_out;
	mu.r0.value = 53.*Mpc/WILT_params.H_out;
	mu.r0.sigma = 1.*Mpc;//0.001*Mpc;//1.*Mpc/WILT_params.H_out;
	//void steepness
	mu.delta_w.lower = -0.995;
	mu.delta_w.upper = -0.95;
	mu.delta_w.value = -0.98;
	mu.delta_w.sigma = -0.001;//-1e-10;//-0.01;
	//void steepness
	//derived parameters which may or may not be used as a prior
	mu.dA.lower = 5.*Mpc;
	mu.dA.upper = 50.*Mpc;
	
	mu.DM.lower = 1e-10;
	mu.DM.upper = 10.;
	
	cout << "performing mcmc ... " << endl;
	//pass Background and Geodesics instances to MCMC class
	ModelExplorer universe;
	universe.initialize_bg(&WILT_params,dMdr_gsl);
	universe.initialize_geo(1200.,1e-6,1600);

	class MCMC myMCMC(&WILT_params, &universe, &mu);
	int ndim = 3;
	int nwalkers = 30;//100;//50;//18;
	int num_samples = 3;//200;//10;//200;
	int burn_in = 1;
	
	cout << "get_likelihood " << myMCMC.get_likelihood(&mu) << endl;
	set_forth_EmceeInstance( ndim, nwalkers, burn_in, num_samples,
                            &myMCMC, &mu);
    
    
}









