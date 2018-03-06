from pylab import *

import numpy
import numpy.random
from scipy.special import gamma as gammaF
from scipy.misc import factorial
import numpy.random


class NegativeBinomial(object):

	# negative binomial p,r parameterization
	@classmethod
	def pdfPR(self,x,p,r):
		"""Negative binomial distribution, p is the probability of a FAILURE. Result is
			the number of failures until the rth success."""
		return gammaF(x + r) / (factorial(x)*gammaF(r)) * p**r * (1-p)**x


	@classmethod
	def convertParameters(self,mu,sigma=None,theta=None,alpha=None):
		"""Convert between mu/[sigma|theta|alpha] parameterization and p/r"""
		if sigma is not None:
			r = mu**2 / (sigma**2-mu)
		elif theta is not None:
			r = theta
		elif alpha is not None:
			r = 1./alpha
		p = r / (r + mu)
		return r,p


	@classmethod
	def pdf(self,x,mu,sigma=None,theta=None,alpha=None):
		"""Regression style negative binomial.  Set one of sigma, theta or alpha.
			Theta is an anti-dispersion or clumping parameter (NB -> poisson as theta -> infinity)
			Alpha is 1/theta, a dispersion parameter (NB -> poisson as alpha -> 0)"""
	
		r,p = self.convertParameters(mu,sigma,theta,alpha)
		return self.pdfPR(x,p,r)
	
	
	@classmethod
	def sample(self,mu,sigma=None,theta=None,alpha=None,size=None):
		"""Draw random samples from a negative binomial distribution"""
		r,p = self.convertParameters(mu,sigma,theta,alpha)
		return numpy.random.negative_binomial(n=r,p=p,size=size)


if __name__ == '__main__':
	mx = 30
	x = arange(0,mx,0.1); p = 0.2; r = theta = 1.5; mu = 6.0
	clf(); 
	h,bins = histogram(NegativeBinomial.sample(mu=mu,theta=theta,size=100000),normed=True,bins=mx,range=(0,mx))
	width = 0.7 * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	bar(bins[:-1], h, align='center', width=width)
	plot(x,NegativeBinomial.pdf(x,mu,theta=theta))
	xlim(xmin=-0.5)


	# mu = 2; clf();
	# for theta in arange(0.3,0.5+0.05,0.05):
	# 	sigma = sqrt(mu**2 / theta + mu)
	# 	plot(x,negativeBinomial(x,mu,sigma),label=str(theta))
