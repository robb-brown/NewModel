from __future__ import print_function
import numpy
import rpy2.robjects as robjects
import rpy2.rinterface as rinterface
from rpy2.robjects.packages import importr
import copy
import traceback
from time import time
from math import *
import scipy.stats
import os.path as pth
import pandas

try:
	from .ModelData import ModelData, formatFloat
except:
	from .ModelData import ModelData, formatFloat
	

from rpy2.robjects.vectors import DataFrame, FloatVector, IntVector, StrVector, ListVector, FactorVector, Vector, Matrix
from collections import OrderedDict


class bcolors:
	BLUE = '\033[94m'
	GREEN = '\033[92m'
	YELLOW = '\033[93m'
	RED = '\033[91m'
	END = '\033[0m'





# New LME
p = pth.join(pth.dirname(__file__),'Rlibs','new')
lme = importr("lme4",lib_loc='"%s"' % p,suppress_messages=True)
lmerTest = importr('lmerTest')
mass = importr("MASS")


def convertVal(val):
	try:
		ret = int(val)
	except:
		try:
			ret = float(val)
		except:
			ret = val
	return ret
	
	
def rToDict(data):
	try:
		named = not (data.names == rinterface.NULL or numpy.all(numpy.array(data.names) == rinterface.NULL))
	except (rinterface.RRuntimeError):
		return None
	except:
		named = None
	rDictTypes = [ DataFrame]
	rArrayTypes = [FloatVector,IntVector,Vector]
	rListTypes=[ListVector,StrVector]
	rMatrixTypes = [Matrix]
	
	kind = None
	if type(data) in rMatrixTypes:
		kind = 'matrix'
	elif hasattr(data,'slotnames') and 'Dim' in list(data.slotnames()) and len(list(data.do_slot('Dim')))>1:
		kind = 'matrixObject'
	elif type(data) in rDictTypes or named:
		kind = 'dict'
	elif type(data) in rListTypes:
		kind = 'list'
	elif type(data) in rArrayTypes:
		kind = 'array'
	elif type(data) == FactorVector:
		kind = 'factorVector'
	try:
		if kind == 'dict':
			sub = [rToDict(elt) for elt in data]
			if not data.names == rinterface.NULL:
				return OrderedDict(zip(data.names, sub))
			else:
				return OrderedDict()
		elif kind == 'list':
			return data[0] if len(data) == 1 else [rToDict(elt) for elt in data]
		elif kind == 'array':
			return data[0] if len(data) == 1 else numpy.array(data)
		elif kind == 'matrix':
			rN,cN = numpy.shape(data)
			if not data.rownames == rinterface.NULL:
				return OrderedDict([(data.rownames[j],OrderedDict([(data.colnames[i],data[i*rN+j]) for i in range(0,cN)])) for j in range(0,rN)])
			elif not data.colnames == rinterface.NULL:
				return OrderedDict([(data.colnames[c],[data[c*rN+r] for r in range(0,rN)]) for c in range(0,cN)])
			else:
				return numpy.array(data)
		elif kind == 'matrixObject':
			mat = numpy.array(data.do_slot('x')).reshape(list(data.do_slot('Dim')))
			dimNames = [i if isinstance(i,list) else [i] for i in rToDict(data.do_slot('Dimnames'))]
			return pandas.DataFrame(mat,columns=dimNames[1],index=dimNames[0])
		elif kind == 'factorVector':
			return [data.levels[i-1] for i in data]
		else:
			if hasattr(data, "rclass"): # An unsupported r class
				return data
			else:
				return data # We reached the end of recursion	
	except:
		print('Exception converting R datatype. Thought it was a {}'.format(kind))
		traceback.print_exc()
		return data
			
			
def sortkeypicker(keynames):
	def getit(adict):
	   composite = [adict[k] for k in keynames]
	   return composite
	return getit



def formatTable(table,pCol=None,**args):	
	if isinstance(table.index,pandas.MultiIndex):
		headerRows = 2
	else:
		headerRows = 1 if table.index.name is None else 2
	index = args.get('index',False if table.index.name is None and isinstance(table.index,pandas.RangeIndex) else True)
	s = table.to_string(index=index).split('\n')
	if not pCol is None:
		for r in range(len(s[headerRows:])):
			if table[pCol][r] < 0.05:
				s[r+headerRows] = bcolors.BLUE+s[r+headerRows]+bcolors.END
	return '\n'.join(s)
	
# def formatTable(coefficients,colnames,rownames=None,pCol='p',colHeaders=None,colSpace0=None,colSpace=None,quiet=True):
# 	rownames = rownames if not rownames is None else ['']*len(coefficients)
# 	s = ''
# 	if not colSpace0:
# 		colSpace0 = int(numpy.max(numpy.array([len(i) for i in rownames])))+1
# 	if not colSpace:
# 		colSpace = 15
# 	if not colHeaders:
# 		colHeaders = colnames
# 	colSpace0S = "%%%ds" % colSpace0
# 	colSpaceS = "%%%ds" % colSpace
# 	s += colSpace0S % ''
# 	s += ''.join([colSpaceS % col for col in colHeaders])
# 	s += '\n'
# 	for r,row in enumerate(rownames):
# 		if not pCol == None:
# 			try:
# 				p = coefficients.iloc[r][pCol]
# 				if p < 0.05:
# 					s += bcolors.BLUE
# 			except:
# 				if not quiet:
# 					traceback.print_exc()
# 					print()
# 		s += colSpace0S % row
# 		for c,col in enumerate(colnames):
# 			try:
# 				try:
# 					s += colSpaceS % formatFloat(float(coefficients.iloc[r][col]),colSpace-2)
# 				except:
# 					s += colSpaceS % coefficients.iloc[r][col]
# 			except:
# 				s += colSpaceS % ' '
# 				if not quiet:
# 					traceback.print_exc()
# 					print()
# 		s += bcolors.END + '\n'
# 	return s




class LinearModel(object):

	def __init__(self,data=None,modelSpec=None,weights=None,title=None,quiet=False):
		self.data = None
		self.modelSpec = None
		self.specifyModel(modelSpec)
		self.model = None
		self.summary = None
		self.R = robjects.r
		self.weights = weights
		self.timings = {}
		self.title = title
		self.setData(data)
		self.quiet = quiet


	def makeNullModelSpec(self):
		null = self.modelSpec
		responding = null.split('~')[0]
		null = responding + '~ 1'
		self.nullModelSpec = null


	def setData(self,data):
		self.model = None
		self.data = copy.deepcopy(data)
		if self.data and self.modelSpec:
			self.execute()


	def specifyModel(self,model):
		self.model = None
		self.modelSpec = model
		self.makeNullModelSpec()
		if self.data and self.modelSpec:
			self.execute()


	def fitModel(self,modelSpec,weights=None):
		if (weights):
			model = self.R.lm(modelSpec,weights=weights,model=True,x=True,y=True,qr=True)
		else:
			model = self.R.lm(modelSpec,model=True,x=True,y=True,qr=True)
		return model


	def execute(self):
		robjects.globalenv["data"] = self.data.getRData()
		if self.weights:
			rdata = self.data.getRData()
			weights = rdata[rdata.names.index(self.weights)]
		else:
			weights = None
		try:
			self.R('attach(data)')
			t1 = time()

			self.model = self.fitModel(self.modelSpec,weights)

			t2 = time()
			self.timings['MainModelFit'] = t2-t1
			self.R('detach(data)')
		except:
			self.R('detach(data)')
			traceback.print_exc()
		self.ModelData = rToDict(self.model)
		self.summary = rToDict(self.R.summary(self.model))
		self.coefficients = self.summary['coefficients']
		self.R['rm']("data")


	def getEstimates(self,estimateData,groups,y,groupValues={},fixedValues=None):
		""" Construct estimates given a ModelData containing a sample timecourse and a list of grouping factors """
		# Make a copy of the estimateData and calculate derived values based on what was done to the real data
		estimateDataT = copy.deepcopy(estimateData)
		estimateDataT.calculateDerivedValues(self.data.derivedValues)
		samples = copy.deepcopy(estimateDataT.samples)
		# set the fixed values
		for key in fixedValues.keys():
			for sample in samples:
				sample[key] = fixedValues[key]

		# Make copies of the basic timecourse for all the groups levels, including the supplied groupValues (if any)
		for group in set(groups + groupValues.keys()):
			newSamples = []
			for sample in samples:
				if group in groupValues.keys():
					for value in groupValues[group]:
						sample[group] = value
						newSamples.append(copy.deepcopy(sample))
				else:
					for value in self.data.getRData()[list(self.data.getRData().names).index(group)].levels:
						sample[group] = value
						newSamples.append(copy.deepcopy(sample))

			samples = newSamples

		# Put the samples back into the ModelData
		estimateDataT.setData(samples,style='stats')

		r = self.R

		# Use the R predict funciton
		(fit,se,df,resdualscale) = r.predict(**{'object':self.model,'newdata':estimateDataT.getRData(),'se.fit':True,'interval':'confidence','level':0.95})
		fit = numpy.array(fit)
		Y = fit[:,0]
		CIs = fit[:,2]-Y

		# Insert the Y values and CIs into the samples
		for s,sample in enumerate(samples):
			sample[y] = Y[s]
			sample['_CI'] = CIs[s]
		estimateDataT.setData(samples,style='stats')
		return estimateDataT


	def getResiduals(self):
		summary = self.R.summary(self.model)
		return numpy.array(self.R['residuals'](summary))


	def getNormalizedResiduals(self):
		residuals = self.getResiduals()
		return numpy.abs(residuals) / numpy.std(numpy.abs(residuals))


	def printSummary(self):
		if self.title:
			print(bcolors.GREEN + '---------------------------------------------------' + bcolors.END)
			print(bcolors.GREEN + '%s' % self.title + bcolors.END)
			print(bcolors.GREEN + '---------------------------------------------------\n' + bcolors.END)

		summary = self.R.summary(self.model)
		# fit summary
		try:
			pvalue = 1-self.R.pf(self.summary['fstatistic'][0],self.summary['fstatistic'][1],self.summary['fstatistic'][2])[0]
			print("Residual standard error: %f on %d degrees of freedom" % (self.summary['sigma'][0],self.summary['df'][1]))
			print("Multiple R-squared: %f, Adjusted R-squared: %f" % (self.summary['r.squared'][0],self.summary['adj.r.squared'][0]))
			print('F-statistic: %f on %d and %d DF, p = ' % (self.summary['fstatistic'][0],self.summary['fstatistic'][1],self.summary['fstatistic'][2]),end='')
			if pvalue < 0.05:
				print(bcolors.BLUE + '%f' % pvalue + bcolors.END)
			else:
				print(bcolors.RED + '%f' % pvalue + bcolors.END)
		except:
			pass
		print('AIC: %f' % (self.R.AIC(self.model)[0]))
		print()

		# coefficients
		coef = self.summary['coefficients']
		space0 = int(numpy.max(numpy.array([len(i) for i in self.summary['coefficients'].keys()])))+1
		space = 15
		print("Coefficients:")
		s = ''
		s += '%*s' % (space0,'') +''.join(['%*s' % (space,i) for i in self.summary['coefficients'].values()[0].keys()])
		s += '\n'
		for r,row in enumerate(self.summary['coefficients'].keys()):
			if coef[row].values()[-1] < 0.05:
				s += bcolors.BLUE
			s += '%*s' % (space0,self.summary['coefficients'].keys()[r])
			s += ''.join(["%*s" % (space,formatFloat(i,space-2)) for i in coef[row].values()])
			s += bcolors.END + '\n'
		print(s)
		print()
		
		
	def testContrast(self,contrast={},printResult=True):
		mcomp = importr("multcomp")
		c1 = [contrast.get(k,0.) for k in self.summary['coefficients'].keys()]
		mc = mcomp.glht(self.model,self.R.matrix(robjects.FloatVector(numpy.array([c1]).ravel()),nrow=1,byrow=True))
		s = self.R.summary(mc)
		estimate = numpy.array(s[-1][2])
		sigma = numpy.array(s[-1][3])
		t = numpy.array(s[-1][4])
		p = numpy.array(s[-1][5])

		if printResult:
			s = 'Test of hypothesis:  '
			s2 = ''.join([u' %+d\u22C5%s' % (c1[k],name) for k,name in enumerate(self.summary['coefficients'].keys()) if not c1[k] == 0])[1:].replace('-','- ').replace('+','+ ').replace(u'1\u22C5',u'')
			if s2.startswith('+'):
				s2 = s2[2:]
			s += s2 + ' == 0'
			print(s)
			print(formatTable({'':dict(Estimate=estimate,Sigma=sigma,t=t,p=p)},rownames=[''],colnames=['Estimate','Sigma','t','p'],colSpace=10))

		if len(p) > 0:
			return dict(estimate=estimate,sigma=sigma,t=t,p=p)
		else:
			return dict(estimate=estimate[0],sigma=sigma[0],t=t[0],p=p[0])
		
		
		
class GeneralizedLinearModel(LinearModel):
	""" Generalized Linear Model.  Works just like the linear model except it takes a family parameter:
		family=	'gaussian': gaussian distribution, identity link, just like the regular LM
				'binomial': 		link = logit
				'Gamma': 			link = inverse
				'inverse.gaussian':	link = 1/mu^2
				'poisson'	: 		link = log
				'quasi'		: 		link = identity, constant variance
				'quasibinomial' : 	link = logit
				'quasipoisson' : 	link = log
	"""

	def __init__(self,data=None,modelSpec=None,weights=None,title=None,quiet=False, family='gaussian'):
		self.family = family
		super(GeneralizedLinearModel,self).__init__(data=data,modelSpec=modelSpec,weights=weights,title=title,quiet=quiet)


	def fitModel(self,modelSpec,weights=None):
		if (weights):
			model = self.R.glm(modelSpec,weights=weights,family=self.family,model=True,x=True,y=True)
		else:
			model = self.R.glm(modelSpec,model=True,family=self.family,x=True,y=True)
		return model


	def printSummary(self):
		summary = rToDict(self.R.summary(self.model))
		if self.title:
			print(bcolors.GREEN + '---------------------------------------------------' + bcolors.END)
			print(bcolors.GREEN + '%s' % self.title + bcolors.END)
			print(bcolors.GREEN + '---------------------------------------------------\n' + bcolors.END)
		print("Generalized %s Linear Model with %s Link" % (summary['family']['family'].title(),summary['family']['link'].title()))

		# fit summary
		try:
			if 'dispersion' in summary:
				print("Dispersion factor for %s: %0.2f" % (summary['family']['family'],summary['dispersion']))
			print("Null Deviance: %0.2f on %d degrees of freedom" % (summary['null.deviance'],summary['df.null']))
			print("Residual Deviance: %0.2f on %d degrees of freedom" % (summary['deviance'],summary['df.residual']))
		except:
			pass
		print("AIC: %0.2f" % (summary['aic']))
		print()

		# coefficients
		coef = summary['coefficients']

		# If we've got a link function, insert the transformed estimate
		if not summary['family']['link'] == 'gaussian':
			for row in coef.keys():
				untransformedVal = coef[row]['Estimate']
				if not row == '(Intercept)':
					untransformedVal += coef['(Intercept)']['Estimate']

				val = summary['family']['linkinv'](untransformedVal)[0]
				coef[row] = OrderedDict([('%s^-1(Est)' % (summary['family']['link']),val)]+list(coef[row].items()))


		space0 = int(numpy.max(numpy.array([len(i) for i in coef.keys()])))+1
		space = 15
		print("Coefficients:")
		s = ''
		s += '%*s' % (space0,'') +''.join(['%*s' % (space,i) for i in coef[coef.keys()[0]].keys()])
		s += '\n'
		for row in coef.keys():
			if coef[row].get('Pr(>|z|)',coef[row].get('Pr(>|z|)',None)) < 0.05:
				s += bcolors.BLUE
			s += '%*s' % (space0,row)
			s += ''.join(["%*s" % (space,formatFloat(i,space-2)) for i in coef[row].values()])
			s += bcolors.END + '\n'
		print(s)
		print()



class NegativeBinomialModel(LinearModel):

	def execute(self):
#		import code; code.interact(local=dict(globals(), **locals()))
#		robjects.globalenv["data"] = self.data.getRData()
		rdata = self.data.getRData()
		if self.weights:
			weights = rdata[rdata.names.index(self.weights)]
		else:
			weights = None
		try:
#			self.R('attach(data)')
			t1 = time()
			if (weights):
				self.model = mass.glm_nb(data=rdata,formula=self.modelSpec,weights=weights,model=True,x=True,y=True)
			else:
				self.model = mass.glm_nb(data=rdata,formula=self.modelSpec,model=True,x=True,y=True)
			t2 = time()
			self.timings['MainModelFit'] = t2-t1

			# NULL model
			if weights:
				self.nullModel = mass.glm_nb(data=rdata,formula=self.nullModelSpec,weights=weights,model=True,x=True,y=True)
			else:
				self.nullModel = mass.glm_nb(data=rdata,formula=self.nullModelSpec,model=True,x=True,y=True)
			robjects.globalenv['model'] = self.model
			robjects.globalenv['nullModel'] = self.nullModel
			self.modelComparison = self.R('anova(model,nullModel)')

#			self.R('detach(data)')
		except:
#			self.R('detach(data)')
			traceback.print_exc()


		# self.ModelData = {}
		# for k,v in self.model.iteritems():
		# 	self.ModelData[k] = v
		# self.summary = {}
		# for k,v in self.R.summary(self.model).iteritems():
		# 	self.summary[k] = v

		self.ModelData = rToDict(self.model)
		self.summary = rToDict(self.R.summary(self.model))


		self.coefficients = self.summary['coefficients']
		self.R['rm']("data")



	def printSummary(self):
		if self.modelComparison:
			print("Comparing to the Null Model (evaluating the fit of the fixed effects):")
			print("	Chi Sq: %f; Chi df: %f;" % (self.modelComparison[4][1],self.modelComparison[5][1]),end='')
			p = self.modelComparison[-1][-1]
			if p < 0.05:
				modelSignificant = True	
				print(bcolors.BLUE + "p = %f" % p + bcolors.END)
			else:
				modelSignificant = False
				print(bcolors.RED + "p = %f" % p + bcolors.END)
			print()

		super(NegativeBinomialModel,self).printSummary()
		print("Family: Negative Binomial, log link function")
		print("Theta = %f" % (self.summary['theta'][0]))






class MixedModel(LinearModel):
	# Wraps R lmer

	def __init__(self,data=None,modelSpec=None,nullModelSpec=None,weights=None,REML=True,mcmc=None,title=None,quiet=False,**args):		
		self.args = args
		self.modelSpec = modelSpec
		self.nullModelSpec = nullModelSpec
		if self.nullModelSpec == None:
			self.makeNullModelSpec()
		self.specifiedFixedEffects = [i for i in self.modelSpec.split('~')[1].split('(')[0].replace('+',' ').replace('*',' ').replace(':',' ').split(' ') if not i == '']
		self.modelComparison = None
		self.REML = REML
		self.mcmc = mcmc
		self.LRResults = None
		super(MixedModel,self).__init__(data,modelSpec,weights,title,quiet)
#		if self.data and self.modelSpec:
#			self.execute()


	def setData(self,data):
		self.data = copy.deepcopy(data)
		if self.data and self.modelSpec:
			try:
				self.execute()
			except:
				pass


	def makeNullModelSpec(self):
		null = self.modelSpec
		responding = null.split('~')[0]
		random = [i for i in null.split('~')[-1].split('+') if i.find('(') >= 0]
		null = responding + '~ 1 +' + '+'.join(random)
		self.nullModelSpec = null


	def specifyModel(self,model,nullModel=None):
		super(MixedModel,self).specifyModel(model)
		if nullModel:
			self.nullModelSpec = nullModel
		else:
			self.makeNullModelSpec()
		self.specifiedFixedEffects = [i for i in self.modelSpec.split('~')[1].split('(')[0].replace('+',' ').replace('*',' ').replace(':',' ').split(' ') if not i == '']


	def makeLevelTerms(self,levelTerms=None,levels=[]):
		if levelTerms == None:
			levelTerms = ['']
		nextTerm = levels[0][0]
		nextLevels = levels[0][1]
		lTerms = []
		for term in levelTerms:
			for level in nextLevels:
				if term == '':
					lTerms.append(nextTerm+level)
				else:
					lTerms.append(term+':'+nextTerm+level)
		if len(levels) > 1:
			lTerms = self.makeLevelTerms(lTerms,levels[1:])
		return lTerms



	@property
	def X(self):
		return numpy.matrix(self.R['getME'](self.model,'X'))

	@property
	def Z(self):
		return numpy.matrix(self.R['as.matrix'](self.R['getME'](self.model,'Z')))

	@property
	def B(self):
		return numpy.matrix(self.R['fixef'](self.model)).T

	@property
	def fixedEffectNames(self):
		return list(self.R['fixef'](self.model).names)

	@property
	def b(self):
		return numpy.array(self.R['ranef'](self.model))

	@property
	def randomEffectNames(self):
		return list(self.R['ranef'](self.model)[0].names)

	@property
	def sigmas(self):
		return [float(self.randomCoefficients.loc[g,l]['variance']) for g in self.sigmaGroups for l in self.sigmaLevels(g)]

	@property
	def sigmaGroups(self):
		re = self.randomCoefficients.reset_index()
		return list(re['group'][re['group'] != 'Residual'].unique())
		
	def sigmaLevels(self,sigmaGroup):
		re = self.randomCoefficients.reset_index()
		return list(re[re['group'] == sigmaGroup]['level'])

	@property
	def sigma(self):
		return float(self.randomCoefficients.loc['Residual','']['variance'])

	@property
	def fixedPredictions(self):
		return self.X*self.B

	@property
	def sigmaF(self):
		return numpy.var(self.fixedPredictions)

	@property	
	def mR2(self):
		return self.sigmaF / (self.sigmaF + sum(self.sigmas) + self.sigma)

	@property
	def cR2(self):
		return (self.sigmaF + sum(self.sigmas)) / (self.sigmaF + sum(self.sigmas) + self.sigma)


	def fitModel(self,data,modelSpec,nullModelSpec=None,weights=None,REML=True,**args):
		if weights:
			model = lmerTest.lmer(self.modelSpec,weights=weights,REML=REML,data=data)
		else:
			model = lmerTest.lmer(self.modelSpec,REML=REML,data=data)

		if nullModelSpec is not None:
			if weights:
				nullModel = lmerTest.lmer(nullModelSpec,REML=REML,weights=weights,data=data)
			else:
				nullModel = lmerTest.lmer(nullModelSpec,REML=REML,data=data)
		else:
			nullModel = None

		return model,nullModel


	def execute(self,REML=None,**args):
		self.summary = None
		self.nullSummary = None
		self.modelComparison = None
		self.LRResults = None
		rdata = self.data.getRData()
		if self.weights:
			weights = rdata[rdata.names.index(self.weights)]
			weights = numpy.array(weights)
			weights /= numpy.average(weights)
			weights = robjects.FloatVector(weights)
		else:
			weights = None
		if REML:
			self.REML = REML
		try:
			t1 = time()
			self.model,self.nullModel = self.fitModel(data=self.data.getRData(mustHave=self.specifiedFixedEffects),modelSpec=self.modelSpec,nullModelSpec=self.nullModelSpec,weights=weights,REML=self.REML,**args)
			t2 = time()

			self.timings['MainModelFit'] = t2-t1
			self.summary = rToDict(self.R.summary(self.model))
			if self.nullModel is not None:
				self.nullSummary = rToDict(self.R.summary(self.nullModel))
				self.modelComparison = rToDict(self.R.anova(self.model,self.nullModel))
		except:
			if not self.quiet:
				traceback.print_exc()

		self.groups = OrderedDict()
		ngrps = self.summary['ngrps']
		for g,group in enumerate(ngrps.keys()):
			self.groups[group] = ngrps[group]
		
		self.totalSamples = self.summary['devcomp']['dims']['N']
		self.coefficients = pandas.DataFrame(self.summary['coefficients']).transpose()
		self.coefficients = self.coefficients.reset_index().rename(columns=dict(index='Effect')).set_index('Effect')

		# Random effects
		self.randomCoefficients = OrderedDict()
		for group in self.summary['varcor'].keys():
			self.randomCoefficients[group] = OrderedDict()
			for level in self.summary['varcor'][group].keys():
				self.randomCoefficients[group][level] = OrderedDict()
				self.randomCoefficients[group][level]['variance'] = self.summary['varcor'][group][level][level]

		self.randomCoefficients['Residual'] = OrderedDict({'':OrderedDict({'variance':self.summary['sigma']**2})})
		re = []
		for group in self.randomCoefficients.keys():
			re.append(pandas.DataFrame(self.randomCoefficients[group]).transpose().reset_index()); 
			re[-1]['group'] = group
			re[-1] = re[-1].rename(columns=dict(index='level'))
			re[-1]['std. dev.'] = numpy.sqrt(re[-1]['variance'])
		self.randomCoefficients = pandas.concat(re).set_index(['group','level'])
		
		# Calculate model terms for use by other methods
		terms = list(self.R.terms(self.R.formula(self.modelSpec)).do_slot('term.labels'))
		fixed = [t for t in terms if t.find('|') < 0]
		atomicTerms = [t.split(':') for t in fixed]
		atomicTerms = set([item for sublist in atomicTerms for item in sublist])
		randomTerms = [t for t in terms if t.find('|') > 0][0]
		randomTerms = set([t.strip() for t in randomTerms.replace('+','|').replace('-','|').replace('*','|').replace(':','|').replace('\\','|').replace('/','|').split('|')])
		formula = self.R.formula('~ ' + '+'.join(fixed))
		y = self.modelSpec.split('~')[0].strip()

		self.terms = terms
		self.atomicTerms = atomicTerms
		self.randomTerms = randomTerms
		self.formula = formula
		self.y = y
		
		return True


	def removeBadSamples(self,data):
		toRemove = []
		terms = list(self.atomicTerms.union(self.randomTerms)) + [self.y]
		terms = [r.replace('.','-') for r in terms if not r == '1']
		for s,sample in enumerate(data.samples):
			for t,term in enumerate(terms):
				try:
					if numpy.isnan(sample[term]):
						toRemove.append(sample)
						break
				except:
					pass

				if sample[term] in [None,'na','nan','None'] or sample.__class__ == numpy.ma.masked.__class__:
					toRemove.append(sample)
					break
		for s in toRemove:
			data.samples.remove(s)
		data.convertToRData()



	def getIndividualEstimates(self,estimateData=None,fixedTerms={},excludeRandomEffects = {},addResiduals=False,primaryKeys=None,allData=False,y=None):
		""" Construct individual estimates.  
			- If you do not supply estimateData the model's data will be used.  
			- Values in the fixedTerms dictionary will be substituted into the data so that you can hold some values constant.
			- excludeRandomEffects is a dictionary where the key is the random effect group (like Subject) and the value
			is a list of specific effects relating to that group to exclude (like '(Intercept)', 'elapsed').  Random effects
			associated with these will not be taken into account.  If excludeRandomEffects == 'all', no random effects will be
			used.  
			- set addResiduals to True if you want to add in the calculated residuals.  This is handy if you want to graph
			the data "corrected" for various effects.
			- primary keys are the set of properties that should be checked to uniquely identify a sample.  These can be set if
			you're using addResiduals.  Unless you know better, use all the items the fixed and random effects depend upon.
			If you don't set this, the effects from the model will be used.  This should work in general.
			- alternately, if you're using all the samples from the data in this model, set allData = True and the indices of your
			samples will be assumed to correspond to the indices of the source data.
			- You can specify y, but I'm not sure what this would do.  It's automatically determined from the model if
			unspecified.

			The return is a ModelData object that can be graphed with ModelPlotter.
		"""

		r = self.R;

		terms = self.terms
		atomicTerms = self.atomicTerms
		randomTerms = self.randomTerms
		y = self.y.replace('.','-')
		formula = self.formula

		modelData = copy.deepcopy(self.data); self.removeBadSamples(modelData)
		estimateData = copy.deepcopy(estimateData) if (not estimateData is None and not allData)  else modelData
		self.removeBadSamples(estimateData)		

		estimateDataT = copy.deepcopy(estimateData)	

		# Override any set fixed effects
		for s,sample in enumerate(estimateDataT.samples):
			for term in fixedTerms.keys():
				if not term in atomicTerms:
					print("The model doesn't have key %s.  Spelling mistake?" % term)
				sample[term] = fixedTerms[term]
		estimateDataT.convertToRData()
		rdata = estimateDataT.getRData()

		# calculation of fixed effects (from model.getEstimates)
#		if self.debug:
#			import code; code.interact(local=dict(globals(), **locals()))

		# try:						# This code uses R's model.matrix.  Unfortunately it doesn't work without all contrast levels represented
		# 	X = r['model.matrix'](formula,data=rdata)
		# except:
		# 	# R's model.matrix is being stupid.  Use ours instead
		# 	print "Ignore the preceeding error.  R is being stupid.  We'll do it ourselves."

		X = self.makeModelMatrix(estimateDataT)
		B = r.fixef(self.model); Bnames = list(B.names); B2 = list(B);
		pp = self.model.do_slot('pp')		# Need to do this separately because R is stupid
		Z = numpy.array(R['as.matrix'](pp.slots['.xData']['Zt']))#self.model.do_slot('Zt')		
		Y = numpy.dot(numpy.array(X),numpy.array(B))
		
		u = pandas.DataFrame(rToDict(r.ranef(self.model)))

		
		# ROBB FIX LATER 
		if 0:
			# Calculation of random effects
			# Add in the random effects
			if not excludeRandomEffects == 'all':
				for s,sample in enumerate(estimateDataT.samples):
					# Find the random effects
					for groupName in u.columns:
						exclude = excludeRandomEffects.get(groupName,[])
						for effectName in u[groupName].index:
							if effectName in exclude:
								continue
						
							individual = individuals.index(sample[groupName])
							effectValue = 1 if effectName == '(Intercept)' else sample[effectName]
							adjustment = effectValue * effect[individual]
							Y[s] += adjustment

		primaryKeys = primaryKeys if primaryKeys is not None else list(atomicTerms.union(randomTerms))
		residuals = numpy.array(r.resid(self.model))
		for s,sample in enumerate(estimateDataT.samples):
			if addResiduals:
				if allData:
					Y[s] += residuals[s]
				else:
					Y[s] += residuals[modelData.getSampleIndex(estimateData.samples[s],primaryKeys)]		# Important that this looks up the unmodified samples

			sample[y] = Y[s]
			estimateDataT.samples[s] = sample

		estimateDataT.convertToRData()
		return estimateDataT




	def getEstimates(self,estimateData,groups=[],y=None,groupValues={},fixedValues={}):
		""" Construct estimates given a ModelData containing a sample timecourse and a list of grouping factors """

		r = self.R


		# Get the terms and remake the formula without the random effects
		terms = list(r.terms(r.formula(self.modelSpec)).do_slot('term.labels'))
		terms = [t for t in terms if t.find('|') < 0]
		atomicTerms = [t.split(':') for t in terms]
		atomicTerms = set([item for sublist in atomicTerms for item in sublist])

		if not y:
			y = self.modelSpec.split('~')[0].strip()

		# Omit terms that do not have multiple values
		#		newTerms = []
		#		for term in terms:
		#			if term.find(':') >= 0:
		#			if len(numpy.unique(estimateDataT.getRData()[list(estimateDataT.getRData().names).index(term)])) > 1:
		#				newTerms.append(term)
		formula = r.formula('~ ' + '+'.join(terms))

		rdata = self.data.getRData()

	# Might want to finish writing this.  Not sure if it's needed.
#		if not estimateData:		# Make estimates to cover the appropriate groups
#			for term in terms:
#				if 'levels' in dir(rdata[list(rdata.names).index(term)]):		# It's a discrete factor

#			mn = numpy.min([numpy.min(self.groupedData[group]['x']) for group in self.groupedData.keys()])
#			mx = numpy.max([numpy.max(self.groupedData[group]['x']) for group in self.groupedData.keys()])
#			self.estimateData = ModelData()
#			self.estimateData.makeEstimate(self.y,self.x,[mn,mx],extraKeys=self.data.samples[0].keys())

#		else:

		# Make a copy of the estimateData and calculate derived values based on what was done to the real data
		estimateDataT = copy.deepcopy(estimateData)
		newDerived = []; keysToSet = groupValues.keys() + fixedValues.keys()
		for d in self.data.derivedValues:
			if not d['key'] in keysToSet:
				newDerived.append(d)

		estimateDataT.calculateDerivedValues(newDerived)
		samples = copy.deepcopy(estimateDataT.samples)
		for sample in samples:
			for k in fixedValues.keys():
				sample[k] = fixedValues[k]	


		# Make copies of the basic timecourse for all the groups levels
		if len(groups) > 0:
			for group in groups:
				newSamples = []
				for sample in samples:
					values = groupValues.get(group,None)
					if values == None:
						try:
							values = list(rdata[list(rdata.names).index(group)].levels)
						except:
							traceback.print_exc()
							print("\nYou must provide values to plot %s as a group\n" % group)
							#return None

					for value in values:
						sample[group] = value	
						newSamples.append(copy.deepcopy(sample))

				samples = newSamples

		# Put the samples back into the ModelData
		estimateDataT.setData(samples,style='stats')

		# Run the derived values again
		estimateDataT.calculateDerivedValues(newDerived)

		# Check to make sure we have all the levels for all the terms
		rdata = estimateDataT.getRData()
		modelData = self.data.getRData()
		samplesToAdd = []
		for term in atomicTerms:
			if term in rdata.colnames:
				i = rdata.colnames.index(term)
				if 'levels' in dir(rdata[i]) and rdata[i].nlevels > 0:
					needLevels = list(modelData[modelData.colnames.index(term)].levels)
					haveLevels = list(rdata[i].levels)
					missingLevels = [l for l in needLevels if not l in haveLevels]
					for l in missingLevels:
						s = dict(estimateDataT.samples[0])
						s[term] = l
						s['_temp'] = 'delete'
						samplesToAdd.append(s)

		if len(samplesToAdd) > 0:
			estimateDataT.addSamples(samplesToAdd,extras={'_temp':'delete'})

		# Check to make sure our multilevel terms have the same reference order as in the model
		rdata = estimateDataT.getRData()
		modelData = self.data.getRData()
		for term in atomicTerms:
			if term in rdata.colnames:
				i = rdata.colnames.index(term)
				if 'levels' in dir(rdata[i]) and rdata[i].nlevels > 0:
					rdata[i] = r.factor(rdata[i],levels=list(modelData[modelData.colnames.index(term)].levels))

		# Grab the X and B
		X = r['model.matrix'](formula,data=rdata)
		B = r.fixef(self.model); Bnames = list(B.names); B2 = list(B);
		V = r.vcov(self.model); vshape = V.do_slot('Dim')
		V = numpy.array(V.do_slot('x')); V = numpy.reshape(V,vshape)
		vkeep = numpy.ones(numpy.shape(V))

		for i,n in enumerate(B.names):
			if not n in X.colnames:
				Bnames.remove(n)
				B2.remove(B[i])
				vkeep[i] = 0; vkeep[:,i] = 0

		B = B2
		V = V[numpy.where(vkeep>0)]; V.shape = (len(B),len(B))

		for i,n in enumerate(X.colnames):
			if not Bnames[i] == n:
				if not self.quiet:
					print("Mismatch in level references (%s, %s)!" % (B.names[i],n))
				return None

		# Calculate the estimates based on the X and B
		Y = numpy.dot(numpy.array(X),numpy.array(B))

		# Calculate the variances based on the X and V

		# Get the between-subject variances		 	DEBUG
#		Vr = float(self.randomCoefficients['subject']['Variance'])

		# Calculate the CIs
		variances = numpy.diag(numpy.dot(numpy.dot(X,V),numpy.array(X).T))# + Vr
		stds = numpy.sqrt(variances)
		try:
			ddf = min([i['dDF'] for i in self.coefficients.values() if i.has_key('dDF')])
			CIs = stds * scipy.stats.t(ddf).ppf(0.975)
		except:
			CIs = None

		# Insert the Y values and CIs into the samples
		samples = estimateDataT.samples
		for s,sample in enumerate(samples):
			sample[y] = Y[s]
			if not CIs == None:
				sample['_CI'] = CIs[s]
			else:
				sample['_CI'] = None
		finalSamples = [s for s in samples if not s.get('_temp')]
		estimateDataT.setData(finalSamples,style='stats')
		return estimateDataT


	def getResiduals(self):
		return numpy.array(self.R.resid(self.model))


	def getNormalizedResiduals(self):
		residuals = self.getResiduals()
		return numpy.abs(residuals) / numpy.std(numpy.abs(residuals))


	def getTermValue(self,sample,term):
		value = None
		if term == '(Intercept)':
			value = 1
		elif term in sample.keys():		# continuous variables
			value = sample[term]
		elif term.find(':')>0:			# Interaction
			value = numpy.prod([self.getTermValue(sample,t) for t in term.split(':')])
		else:			# Dummy coded categorical variables
			for t in sample.keys():
				if term.startswith(t):
					value = int(sample[t] == term.replace(t,''))
		return value


	def makeModelMatrix(self,data=None):
		data = data if data is not None else self.data
		Xt = numpy.zeros((len(data.samples),len(list(self.coefficients.index))),numpy.float64)
		for s,sample in enumerate(data.samples):
			for i,term in enumerate(list(self.coefficients.index)):
				value = self.getTermValue(sample,term)
				if value is None:
					print("Could not find a value for %s" % term)

				Xt[s,i] = value
		return Xt



	def testContrast(self,contrast={},printResult=True):
		mcomp = importr("multcomp")
		c1 = [contrast.get(k,0.) for k in list(self.coefficients.index)]
		mc = mcomp.glht(self.model,self.R.matrix(robjects.FloatVector(numpy.array([c1]).ravel()),nrow=1,byrow=True))
		s = self.R.summary(mc)
		estimate = numpy.array(s[-1][2])
		sigma = numpy.array(s[-1][3])
		t = numpy.array(s[-1][4])
		p = numpy.array(s[-1][5])

		if printResult:
			s = 'Test of hypothesis:  '
			[contrast.get(k,0.) for k in list(self.coefficients.index)]
			s2 = ''.join([u' %+d\u22C5%s' % (c1[k],name) for k,name in enumerate(list(self.coefficients.index)) if not c1[k] == 0])[1:].replace('-','- ').replace('+','+ ').replace(u'1\u22C5',u'')
			if s2.startswith('+'):
				s2 = s2[2:]
			s += s2 + ' == 0'
			print(s)
			d = pandas.DataFrame([OrderedDict(Estimate=estimate[0],Sigma=sigma[0],t=t[0],p=p[0])])[['Estimate','Sigma','t','p']]
			print(formatTable(d,pCol='p'))

		if len(p) > 0:
			return dict(estimate=estimate,sigma=sigma,t=t,p=p)
		else:
			return dict(estimate=estimate[0],sigma=sigma[0],t=t[0],p=p[0])


	def printSummary(self,null=False):
		if self.title:
			print(bcolors.GREEN + '---------------------------------------------------' + bcolors.END)
			print(bcolors.GREEN + '%s' % self.title + bcolors.END)
			print(bcolors.GREEN + '---------------------------------------------------\n' + bcolors.END)
		summary = self.summary
		print(summary['methTitle']); print()
		print("Number of samples: %d in %d total groups" % (self.totalSamples,summary['devcomp']['dims']['q']))
		print("Groups:	")
		ngrps = summary['ngrps']
		for g,group in enumerate(ngrps.keys()):
			print('%s: %d; ' % (group,ngrps[group]),end='')
		print(); print()
		if self.modelComparison:
			print("Comparing to the Null Model (evaluating the fit of the fixed effects):")
			print("	Chi Sq: {chisq}; Chi df: {df}; ".format(chisq=self.modelComparison['Chisq'][-1],df=self.modelComparison['Chi Df'][-1]),end='')
			p = self.modelComparison['Pr(>Chisq)'][-1]
			if p < 0.05:
				modelSignificant = True	
				print(bcolors.BLUE + "p = %f" % p + bcolors.END)
			else:
				modelSignificant = False
				print(bcolors.RED + "p = %f" % p + bcolors.END)
			print()

		print('Log Likelihood: {}'.format(self.summary['logLik']))
		if self.summary['devcomp']['cmp'].get('REML',False) and self.modelComparison:
			print("Warning, fit with REML:")
			print("	AIC: {AIC};  BIC: {BIC}   (lower is better)".format(AIC=self.modelComparison['AIC'][-1],BIC=self.modelComparison['BIC'][-1]))
			print()
			
		print("Marginal R^2 (fixed effects):               %0.4f" % (self.mR2))
		print("Conditional R^2 (fixed+random effects):     %0.4f" % (self.cR2))
		print()

		print("Random effects:")
		colnames = ['variance','std. dev.'];
		print(formatTable(self.randomCoefficients[colnames]))

		print("Fixed effects:")
		colnames = ['Estimate','Std. Error','df','t value','Pr(>|t|)']
		colHeaders = ['Estimate','Std. Error', 'df', 't', 'p']
		print(formatTable(self.coefficients.rename(columns=dict(zip(colnames,colHeaders)))[colHeaders],pCol='p'))
		print()

		# Correlation of fixed effects
		print('Correlation of Fixed Effects:')
		print(formatTable(summary['vcov']))



	def printSummaryNew(self,null=False):
		if null:
			summary = merSummary(self.R.summary(self.nullModel))
		else:
			summary = merSummary(self.R.summary(self.model))
		print(summary.methTitle[0])
		print("Number of samples: %d in %d total groups" % (summary.dims[1],summary.dims[3]))
		print("Groups:	")
		ngrps = summary.ngrps
		for g,group in enumerate(ngrps.names):
			print('%s: %d;' % (group,ngrps[g]),end='')
		print(); print()
		if self.modelComparison:
			print("Comparing to the Null Model:")
			print("Chi Sq: %f; Chi df: %f;" % (self.modelComparison[4][1],self.modelComparison[5][1]),end='')
			p = self.modelComparison[6][1]
			if p < 0.05:
				modelSignificant = True	
				print(bcolors.BLUE + "p = %f" % p + bcolors.END)
			else:
				modelSignificant = False
				print(bcolors.RED + "p = %f" % p + bcolors.END)
			print()
		print(summary.AICtab)
		print("Random effects:")
		print(str(summary.REmat).replace('"',' '))
		print("Fixed effects:")
		d = summary.coefs
		dd = numpy.array(d)
		n = self.summary.dims[1]
		p = self.summary.dims[2]
		q = self.summary.dims[3]
		colnames = ['estimate','stderr','nDF','dDF','f value','p value']
		rownames = self.coefficients.keys()
		s = ''
		colSpace0 = int(numpy.max(numpy.array([len(i) for i in rownames])))+1
		colSpace = 15
		colSpace0S = "%%%ds" % colSpace0
		colSpaceS = "%%%ds" % colSpace
		s += colSpace0S % ''
		s += ''.join([colSpaceS % col for col in colnames])
		if self.LRResults:
			s += colSpaceS % '	LR p value'
#		s += colSpaceS % 'Est. p Range'
		s += '\n'
		lrp = None


		for r,row in enumerate(rownames):
			factor = row
			try:
				p = self.coefficients[factor]['p']
				if p < 0.05:
					s += bcolors.BLUE
			except:
				pass
			s += colSpace0S % row
			for c,col in enumerate(colnames):
				try:
					s += colSpaceS % formatFloat(float(self.coefficients[factor][col]),colSpace-2)
				except:
					s += colSpaceS % '%*s' % (colSpace,'')  
			s += bcolors.END + '\n'
		print(s)
		print()
		if (self.mcmc):
			print("MCMC Results:")
			print("Fixed effects:")

			# Print fixed effects
			d = self.mcmc[0]
			colnames = d.colnames
			rownames = d.rownames
			rowNameLengths = []
			for r in rownames:
				rowNameLengths.append(len(r))
			s = ''
			colSpace0 = max(rowNameLengths)
			colSpace = 11
			colSpace0S = "%%%ds" % colSpace0
			colSpaceS = "%%%ds" % colSpace
			s += colSpace0S % ''
			for col in colnames:
				s += colSpaceS % col
			s += '\n'
			for r,row in enumerate(rownames):
				significant = float(d[4][r]) < 0.05
				if significant: s += bcolors.BLUE
				s += colSpace0S % row
				for c,col in enumerate(colnames):
					s += colSpaceS % formatFloat(float(d[c][r]),colSpace-2)
				if significant: s += bcolors.END
				s += '\n'
			print(s)

			print()
			print("Random effects:")
			print(self.mcmc[1])
		# Correlation of fixed effects
		s = str(summary)
		print(s[s.find('Correlation of Fixed Effects:'):])




class GeneralizedMixedModel(MixedModel):

	""" Generalized Linear Mixed Model.  Works just like the mixed model except it takes a family parameter:
		family=	'gaussian': gaussian distribution, identity link, just like the regular GLMM
				'binomial': 		link = logit
				'Gamma': 			link = inverse
				'inverse.gaussian':	link = 1/mu^2
				'poisson'	: 		link = log
				'quasi'		: 		link = identity, constant variance
				'quasibinomial' : 	link = logit
				'quasipoisson' : 	link = log
	"""

	def __init__(self,data=None,modelSpec=None,nullModelSpec=None,weights=None,REML=True,mcmc=None,title=None,quiet=False, family='gaussian',**args):
		self.family = family
		super(GeneralizedMixedModel,self).__init__(data=data,modelSpec=modelSpec,nullModelSpec=nullModelSpec,weights=weights,REML=REML,mcmc=mcmc,title=title,quiet=quiet,**args)


	def fitModel(self,modelSpec,nullModelSpec=None,weights=None,REML=True,**args):
		if weights:
			model = lme.glmer(self.modelSpec,family=self.family,weights=weights,REML=REML,model=True,x=True)
		else:
			model = lme.glmer(self.modelSpec,family=self.family,REML=REML,model=True,x=True)

		if nullModelSpec is not None:
			if weights:
				nullModel = lme.lmer(nullModelSpec,REML=REML,weights=weights,model=True,x=True)
			else:
				nullModel = lme.lmer(nullModelSpec,REML=REML,model=True,x=True)
		else:
			nullModel = None

		return model,nullModel


	def printSummary(self):
		super(GeneralizedMixedModel,self).printSummary()




class NegativeBinomialMixedModel(GeneralizedMixedModel):

	""" Generalized Linear Mixed Model with a negative binomial distribution and log link. 
		Appears not to be implemented in my version of LME4.  Have to get on the latest version...."""

	def __init__(self,data=None,modelSpec=None,nullModelSpec=None,weights=None,REML=True,mcmc=None,title=None,quiet=False):
		self.family = family
		super(GeneralizedMixedModel,self).__init__(data=data,modelSpec=modelSpec,nullModelSpec=nullModelSpec,weights=weights,REML=REML,mcmc=mcmc,title=title,quiet=quiet,family='negative binomial')


	def fitModel(self,modelSpec,nullModelSpec=None,weights=None,REML=True,**args):
		if weights:
			model = lme.glmer(self.modelSpec,family=self.family,weights=weights,REML=REML,model=True,x=True)
		else:
			model = lme.glmer(self.modelSpec,family=self.family,REML=REML,model=True,x=True)

		if nullModelSpec is not None:
			if weights:
				nullModel = lme.lmer(nullModelSpec,REML=REML,weights=weights,model=True,x=True)
			else:
				nullModel = lme.lmer(nullModelSpec,REML=REML,model=True,x=True)
		else:
			nullModel = None

		return model,nullModel


	def printSummary(self):
		super(GeneralizedMixedModel,self).printSummary()



R = robjects.r


# MAIN PART
if __name__ == '__main__':
	from ModelData import ModelData
	import random
	data = ModelData([dict(x=i,y=i*2+random.random(),sub=int(i/2)) for i in range(0,10)])

	# Mixed Model
	# 
	# robjects.globalenv["data"] = data.getRData(mustHave=['x'])
	# R('attach(data)')
	# model = lmerTest.lmer('y ~ x + (1|sub)',model=True,x=True)
	# summary = R.summary(model)
	# randSummary = lmerTest.rand(model)


	# Negative Binomial Mixed
	from lingo.useDB2 import *
	db = labDB
	counts = ModelData(db.getValues('subject,timepoint,lesionCount',table='MaskMetadata',contrast='NewT2Lesions'))

	intervals = dict([(s['subject'],dict(zip(*db.getOrderedTimepoints(subject=s['subject'],includeIntervals=True)))) for s in counts.samples])

	counts.calculateDerivedValues([ {	'key' : 'age', 'criteria' : 'db.getAge(subject,timepoint)', 'variables' : dict(db=db), 'doNotStore' : True },
									{	'key' : 'diagnosis', 'criteria' : 'db.getDiagnosis(subject,timepoint)', 'variables' : dict(db=db), 'doNotStore' : True },
									{	'key' : 'treatment', 'criteria' : 'db.getTreatment(subject,timepoint)', 'variables' : dict(db=db), 'doNotStore' : True },
									{	'key' : 'ageC', 'criteria' : "age-16"},
									{	'key' : 'treated', 'criteria' : "0 if treatment == 'none' else 1 if treatment in ['copaxone','betaseron','avonex','rebif'] else None", 'variables' : dict(db=db), 'doNotStore' : True },
									{	'key' : 'interval', 'criteria' : 'intervals[subject][timepoint]', 'variables':dict(intervals=intervals),'doNotStore':True},
									{	'key' : 'lesionRate', 'criteria' : 'lesionCount / (interval / 365.)'}
									])

	#counts.filter("[s for s in samples if s['diagnosis'] in ['MS','monoADS']]")
	#formula='lesionCount ~ ageC*diagnosis + (1|subject)'

	counts.filter("[s for s in samples if s['diagnosis'] in ['MS']]")
	formula='lesionCount ~ 1 + (1|subject)'
	counts2 = counts.filtered("[s for s in samples if 300 < s['interval'] < 400]")



	model = lme.glmer_nb(formula=formula,data=counts2.getRData(),nAGQ=2)

	print(R.summary(model))

	# Plot
	mx = max([s['lesionCount'] for s in counts.samples])
	h,bins = histogram([s['lesionCount'] for s in counts.samples],normed=True,bins=mx,range=(0,mx))
	width = 0.7 * (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	bar(bins[:-1], h, align='center', width=width)
	x = arange(0,mx)
	from NegativeBinomialDistribution import NegativeBinomial
	theta = float(R.summary(model)[6][0].split('(')[1].split(')')[0])
	mu = R.fixef(model)[0]
	plot(x,NegativeBinomial.pdf(x,e**mu,theta=theta))
	# End plot



	# # Test contrasts
	# mcomp = importr("multcomp")
	# 
	# # Can't just do rToDict on summary(model) because the coefficients section is a weird matrix
	# coefficients = rToDict(lme.fixef(model)).keys()	
	# 
	# # This is an improvement on the existing testcontrast that lets you have two vectors
	# # so you don't have to do the subtraction yourself
	# contrast1 = {'age': 10,'diagnosisMS': 1, 'age:diagnosisMS':10}
	# contrast2 = {'age': 11,'diagnosisMS':1, 'age:diagnosisMS' : 11}
	# c1 = array([contrast1.get(k,0.) for k in coefficients])
	# c2 = array([contrast2.get(k,0.) for k in coefficients])
	# cc = c2-c1
	# mc = mcomp.glht(model,R.matrix(robjects.FloatVector(numpy.array([cc]).ravel()),nrow=1,byrow=True))
	# print R.summary(mc)





	#model = MixedModel(data,'y ~ x + (1|sub)')
	#model.printSummary()





