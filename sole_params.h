#ifndef SOLE_PARAMS
#define SOLE_PARAMS

const double PI = 3.14159265358979323846264338327950288419716939937510;
const double c_ = 299792458.; // ms^-1
const double Mpc= 1.;
const double Gpc= 1e3*Mpc;
const double km_over_c= 1e5/c_; 
//e.g H_in;                    // 0.5-0.85 units km s^-1 Mpc^-1
//Hoverc_in = H_in*km_over_c;  // units of Mpc^-1
const double ageMpc =306.60139383811764; //age in Mega Parsecs
// age = 15. billion years
//     = 15. *10**9*(365.25*24.*60.*60.) s
//     = 15. *10**9*(365.25*24.*60.*60.)* 299792458. m
//     = 15. *10**9*(365.25*24.*60.*60.)* 299792458. *3.24077929*10**(-23) Mpc
//     = 15. * 306.60139383811764 Mpc

struct Params
{
	/**
	Most LTB solutions modeling a void have a density profile defined interms of 
	the void radius and steepness. The Hubble expansion rate in these models is 
	characterized by its value at the location of the observer and its asymptotic 
	value. Here we collect these parameters.
	**/
	
	double H_in;  // 0.5-0.85 units km s^-1 Mpc^-1
	double Hoverc_in;// =H_in*1e5/c;  // units of Mpc^-1
	double H_out;  // 0.3-0.7 units km s^-1 Mpc^-1
	double Hoverc_out;// = H_out*1e5/c;  // units of Mpc^-1
	//double H_not;  // 0.5-0.95 units km s^-1 Mpc^-1
	//double Hoverc_not;// = H_not*1e5/c;  // units of Mpc^-1
	//double Omega_in;  // 0.05-0.35
	//if Lambda is nonzero check equations for correct units. [Lambda]=[R]^-2
	//if Lambda is nonzero check equations for correct units. [Lambda]=[R]^-2
	double OmegaM_out;// = 0.99999 - Lambda;  //
	double OmegaR_out;
	double Omega_Lambda;
	double Lambda;  //  0.7
	
	double r0;  //  2.5*Gpc 3.5 0.33 0.3-4.5 units Gpc
	double delta_r;  // 0.2*r0  0.1r0-0.9r0
	double delta_w;
	
	double age;//in Mpc
	
	
	Params(double H_in,double H_out, double OmegaM_out,
	       double r0, double delta_r, double delta_w)
	{
	  Params::H_in = H_in;
	  Params::H_out = H_out;
	  Hoverc_in = H_in*km_over_c;
	  Hoverc_out = H_out*km_over_c;
	  //Hoverc_not = H_not*km_over_c;
	  
	  OmegaR_out = 4.1834e-5/(H_out*H_out);
	  Params::OmegaM_out=OmegaM_out;
	  Omega_Lambda=(1. - OmegaM_out - OmegaR_out)*0.999;
	  Lambda= Omega_Lambda*3.*Hoverc_out*Hoverc_out;
	  //Params::Omega_in = Omega_in;
	  Params::r0 = r0;
	  Params::delta_r=delta_r;
	  Params::delta_w=delta_w;
	  
	  age = ageMpc*13.819;
	};
	
};

#endif
