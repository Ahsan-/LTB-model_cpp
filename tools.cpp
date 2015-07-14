#include <cmath>
#include <assert.h>
#include "tools.h"
void linspace(double *t, double t_init, double t_end, int num_pt,
                 bool end_point = true)
{
//t_init = t_series(R=R_init,E=E,M=M,Lambda=Lambda)
	double step;
	if (end_point)
		step = (t_end-t_init)/ double(num_pt-1);
	else
		step = (t_end-t_init)/ double(num_pt);
	
	for (int i=0; i<num_pt; i++)
		//t[i] = t_init + double(i)*step;
		*(t+i) = t_init + double(i)*step;
	
}

void sample_radial_coord(double *r_vec, struct Params *p, double r_min=1e-4,
                         double r_max=20*1e3, int num_pt1=100, int num_pt2=100)
{
	/**The comoving radial coordinate r is assumed to be in mega parsecs
	num_pt1: 
			must be an integer divisible by 5. Then 3/5th of num_pt1 is linearly
	distributed from r0-1.5delta_r to r0+1.5delta_r1 and 1/5th is linearly
	distributed from r0-3delta_r to r0-1.5delta_r and 1/5th is linearly
	distributed from r0+1.5delta_r to r0+3delta_r
	num_pt2:
			must be divisible by 2. Half of it is logrithmicly distributed
	between r_init and r0-3delta_r and the other half between r0+3deltar_r
	and r_max
	**/
	
	assert( (num_pt1 % 5) ==0 &&"num_pt1 must be an integer divisible by 5");
	assert( (num_pt2 % 2) ==0 && "num_pt2 must be even");

	if ( (p->r0 - 3.*p->delta_r ) > 0. )
	{
		linspace(r_vec,
		         log(r_min),log(p->r0-3.*p->delta_r),num_pt2/2,false);
		for (int i = 0; i<num_pt2/2; i++)
			r_vec[i] = exp(r_vec[i]);

		linspace(&r_vec[num_pt2/2],
	         	p->r0-3.*p->delta_r,p->r0-1.5*p->delta_r,num_pt1/5,false);
	
		linspace(&r_vec[num_pt2/2+num_pt1/5],
	         	p->r0-1.5*p->delta_r,p->r0+1.5*p->delta_r,(num_pt1/5)*3,false);
	
		linspace(&r_vec[num_pt2/2+(num_pt1/5)*4],
	    	    p->r0+1.5*p->delta_r,p->r0+3.*p->delta_r,num_pt1/5,false);
	
		linspace(&r_vec[num_pt2/2+num_pt1],
	    	     log(p->r0+3.*p->delta_r),log(r_max),num_pt2/2,true);
		for(int i=num_pt2/2+num_pt1; i<num_pt2+num_pt1; i++)
			r_vec[i] = exp(r_vec[i]);
		}
	else{
		linspace(r_vec,
		         log(r_min),log(1.),num_pt2/2,false);
		for (int i=0; i<num_pt2/2; i++)
			r_vec[i] = exp(r_vec[i]);
		
		linspace(&r_vec[num_pt2/2],
		         1.,p->r0+3.*p->delta_r,num_pt1,false);
		
		linspace(&r_vec[num_pt2/2+num_pt1],
		         log(p->r0+3.*p->delta_r),log(r_max),num_pt2/2,true);
		for(int i=num_pt2/2+num_pt1; i<num_pt2+num_pt1; i++)
			r_vec[i] = exp(r_vec[i]);
		}

}


