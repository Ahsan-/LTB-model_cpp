#!python2.7
#cython: boundscheck=False, wraparound=False
import numpy as np
cimport numpy as np
from cython.operator cimport dereference as deref
from libcpp cimport bool
from libc.string cimport memcpy
from cython.view cimport array as cvarray
import emcee
from scipy import multiply


#cpdef lnprob(x, ivar):
#	return -0.5*np.sum(ivar*x**2)
#
#ndim, nwalkers = 4, 100
#ivar = 1./np.random.rand(ndim)
#p0 = [np.random.rand(ndim) for i in range(nwalkers)]
#
#sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[ivar])
#sampler.run_mcmc(p0,100)

class EmceeWrapper:
	def __init__(self, ndim, llh):
		self.ndim = ndim
		self.llh = llh
	
	def plot(self):
		from matplotlib import pylab as plt
		for i in range(self.ndim):
			plt.figure()
			plt.hist(self.sampler.flatchain[:,i],100,color="k",histtype="step")
			plt.title("Dimension {0:d}".format(i))
		plt.show()
		
		samples = self.sampler.chain[:,:, :].reshape((-1, self.ndim))
		import emcee.triangle as triangle
		fig = triangle.corner(samples)
		fig.savefig("triangle.png")
		
	def lnprob(self,x,mu,icov):
		diff = x-mu
		return -np.dot(diff,np.dot(icov,diff))/2.
	
	def __call__(self,means, cov,nwalkers=100, burn_in=100, num_samples=1000):
		print " ndim nwalkers, num_samples ", self.ndim, nwalkers, num_samples
		self.nwalkers = nwalkers
		#p0 = np.random.rand(self.ndim * nwalkers ).reshape((nwalkers,self.ndim))
		p0 = np.random.multivariate_normal(means,cov,nwalkers)
		print " means and cov ", means, cov
		sampler = emcee.EnsembleSampler(nwalkers, self.ndim, self.llh )
		pos, prob, state, blobs = sampler.run_mcmc(p0,burn_in)
		#pos, prob, state, blobs = sampler.run_mcmc(p0,1)
		#sampler.reset()
		#sampler.run_mcmc(pos,num_samples)
		self.sampler = sampler
		fc = open("fchain.dat","w")
		fp = open("fprob.dat","w")
		fb = open("fblob.dat","w")
		fc.close()
		fp.close()
		fb.close()
		for pos, prob, state, blobs in sampler.sample(pos,prob,state,blobs,\
        iterations=num_samples,storechain=True):#500
			fc = open("fchain.dat","a")
			fc.write("\n".join(["\t".join([str(q) for q in p]) for p in pos]))
			fc.write("\n")
			fc.close()
			fb = open("fblob.dat","a")
			fb.write("\n".join(["\t".join([str(b) for b in blob]) for blob in blobs]))
			fb.write("\n")
			fb.close()
			fp = open("fprob.dat","a")
			fp.write("\n".join( str(pr) for pr in prob))
			fp.write("\n")
			fp.close()
		
		print ("Mean acceptance fraction: {0:.3f}".format(
		       np.mean(sampler.acceptance_fraction))) 
		
		self.plot()
		

cdef extern from "model_explorer2.h":
	cdef cppclass type_MCMC_param:
		double lower
		double upper
		double value
		double sigma
	cdef cppclass MCMC_params:
		type_MCMC_param r_loc
		type_MCMC_param r0
		type_MCMC_param delta_w
		type_MCMC_param dA
		type_MCMC_param DM
		type_MCMC_param zdec
		type_MCMC_param tdec
		type_MCMC_param rdec
		bool reset_bg
	
	cdef cppclass MCMC:
		double get_likelihood(MCMC_params *)
		void set_derived_parameters(MCMC_params *, MCMC_params * = NULL)


cdef class Likelihood:
	cdef MCMC *mcmc
	cdef MCMC_params *mu
	cdef void set_init(self, MCMC *mcmc_, MCMC_params *mu_):
		self.mcmc = mcmc_
		self.mu = mu_
	
	cpdef lnprob(self,x):
		#check priors are satisfied if not reject proposed parameters
		if ( x[0] < deref(self.mu).r_loc.lower or x[0] > deref(self.mu).r_loc.upper or
		     x[1] < deref(self.mu).r0.lower or x[1] > deref(self.mu).r0.upper or 
		     x[2] < deref(self.mu).delta_w.lower or x[2] > deref(self.mu).delta_w.upper):
			return -np.inf, np.zeros(5)
		
		cdef bool reset_bg = True
		deref(self.mu).reset_bg = reset_bg
		deref(self.mu).r_loc.value = x[0]
		deref(self.mu).r0.value = x[1]
		deref(self.mu).delta_w.value = x[2]
		
		cdef double lnk = deref(self.mcmc).get_likelihood( self.mu )
		blobs = np.array([deref(self.mu).dA.value, deref(self.mu).DM.value,
		                  deref(self.mu).zdec.value, deref(self.mu).tdec.value,
		                  deref(self.mu).rdec.value])
		#return derived parameters as a blob
		return lnk, blobs

	
cdef public void set_forth_EmceeInstance(int ndim, int nwalkers, 
                                         int burn_in, int num_samples,
                                         MCMC *mcmc, MCMC_params *mu):
#{
	llh = Likelihood()
	llh.set_init(mcmc, mu)
	means = np.array([deref(mu).r_loc.value, deref(mu).r0.value, 
	                 deref(mu).delta_w.value])
	sigmas = np.array([deref(mu).r_loc.sigma, deref(mu).r0.sigma, 
	                 deref(mu).delta_w.sigma])
	cov = np.diag(sigmas)
	cov = multiply(cov,cov)
	run = EmceeWrapper(ndim, llh.lnprob)
	run(means, cov, nwalkers,burn_in,num_samples)
	
#}

























