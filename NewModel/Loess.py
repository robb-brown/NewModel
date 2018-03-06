import numpy
import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects.packages import importr
import copy
import traceback
from time import time
from math import *
import scipy.stats


class Loess(object):
	
	def __init__(self,x=None,y=None,span=0.75):
		self.x = x
		self.y = y
		self.R = robjects.r
		self.loess = None
		self.weights = None
		self.span = span
		self.std = None
		self.xs = None
		self.prediction = None
		self.estimation = None
		self.ci = None
		
	
	def setData(self,x=None,y=None,weights=None):
		if x:
			self.x = x
		if y:
			self.y = y
		if weights:
			self.weights = weights
		self.loess = None
			
	def predict(self,xs):
		if not self.loess:
			data = self.R['data.frame'](x=robjects.FloatVector(self.x),y=robjects.FloatVector(self.y))
			if self.weights:
				w = robjects.FloatVector(self.weights)
				self.loess = self.R.loess('y~x',data,weights=w,span=self.span)
			else:
				self.loess = self.R.loess('y~x',data,span=self.span)
		pred = self.R.predict(self.loess,self.R['data.frame'](x=robjects.FloatVector(xs)),se=True)
		prediction = numpy.ma.array(pred[0])
		se = numpy.ma.array(pred[1])
		df = float(pred[-1][0])
		ci = se * scipy.stats.t(df).ppf(0.975)
		prediction.mask = numpy.where(numpy.isnan(prediction),1,0)
		ci.mask = prediction.mask
		self.std = float(pred[2][0])
		self.xs = xs
		self.prediction = prediction
		self.ci = ci
		return prediction,ci
		
		
	def getNormalizedResiduals(self):
		if self.estimation == None:
			self.estimation,_ = self.predict(self.x)
		residuals = numpy.abs(self.estimation-numpy.array(self.y))
		return residuals / self.std