from __future__ import print_function

try:
	from .ModelData import *
	from .Loess import Loess
	from .NewModel import MixedModel
except:
	from ModelData import *
	from Loess import Loess
	from NewModel import MixedModel
		
	
import pylab
import numpy
import numpy.ma as ma
try:
	import seaborn as sns
except:
	sns = None	

class ModelPlotter(object):
	
	def __init__(self,data=None,model=None,x=None,y=None,weights=None,groups=[],specialGroups=[],silentGroups=[],
					individuals=[],mean=True,plotIndividuals=True,span=0.6,ialpha=None,malpha=1.0,cialpha=None,meanStyle=None,errorBars='95',
					estimateData=None,binFunction=None,iLinewidth=1,mLinewidth=1,linestyle='-',
					colors=None,styles=None,derivedValues=None,groupValues={},fixedValues={},xUnitConversion=None,
					yUnitConversion=None, drawLegend=True,groupOrdering=None,includeN=True,**args):
		self.model = model
		if self.model:
			self.data = self.model.data
			self.estimateData = estimateData
			self.meanStyle = 'model'
		else:
			self.data = data
			self.meanStyle = 'loess'
			self.estimateData = None
		if data:
			self.data = data
		if meanStyle:
			self.meanStyle = meanStyle
		self.errorBars = errorBars
		self.x = x
		self.y = y
		self.weights = weights
		self.groups = groups
		self.specialGroups = specialGroups
		self.silentGroups = silentGroups		
		self.individuals = individuals
		self.plotIndividuals = plotIndividuals
		self.mean = mean
		self.span = span
		self.ialpha = ialpha
		self.malpha = malpha
		self.cialpha = cialpha if cialpha is not None else self.malpha/10.
		
		self.connect = args.get('connect',True if individuals and len(individuals) > 0 else False)
		self.marker = args.get('marker',None if self.connect else 'o')
		self.markerSize = args.get('markerSize',5)
		self.markerAlpha = args.get('markerAlpha',0.5)
		
		self.label = args.get('label',None)
		self.includeN = includeN
		
		self.groupColorKey = {}
		self.groupStyleKey = {}
		self.groupFormatKey = {}
		self.groupedData = {}
		self.colorCounter = -1
		self.styleCounter = -1
		if colors:
			self.colors = colors
		else:
			if sns is not None:
				self.colors = sns.color_palette()
			else:
				self.colors = ['g','b','c','m','y','k','r']
		if styles:
			self.styles = styles
		else:
			self.styles = ['-','--',':','-.']
		self.estimatedTimecourses = None
		self.binFunction = binFunction
		self.iLinewidth = iLinewidth
		self.mLinewidth = mLinewidth
		self.linestyle = linestyle
		self.derivedValues = None
		self.groupValues = groupValues
		self.fixedValues = fixedValues
		self.xUnitConversion = xUnitConversion
		self.yUnitConversion = yUnitConversion
		self.drawLegend = drawLegend
		self.groupOrdering = groupOrdering
			
		
		
	def setData(self,data):
		self.data = data
		self.estimatedTimecourses = None
		
	def setModel(self,model):
		self.model = model
		self.estimatedTimecourses = None
		
	def setGroups(self,groups):
		self.groups = groups
		

	def getNestedVal(self,data,keys,create=True):
		if create and not keys[0] in data: data[keys[0]] = {}
		if len(keys) == 1:
			return data[keys[0]]
		else:
			return self.getNestedVal(data[keys[0]],keys[1:])

	def getGroupChain(self,groups,sample):
		groupChain = []
		for group in groups:
			groupChain.append(sample[group])
		return groupChain


	def nestedPlot(self,data,groups=None,chain='',color=0):
		for k,key in enumerate(data.keys()):
			if len(groups) > 1:
				color += 1
				try:
					self.nestedPlot(data[key],groups=groups[1:],chain=chain+'%s-'%key,color=color)
				except:
					print()
					traceback.print_exc()
					print()
			else:
				xD,yD = zip(*sorted(zip(data[key]['x'],data[key]['y'])))
				if k == 0: 
					label = chain[0:-1]
				else: label = ''
				if (self.connect):
#					print "Plotting %s %s %s" % (xD,key,label)
					pylab.plot(xD,yD,alpha=self.ialpha,color=self.colors[color%len(self.colors)],label=label)
				else:
					pylab.scatter(xD,yD,alpha=self.ialpha,color=self.colors[color%len(self.colors)],label=label)
		

	def getTimecourses(self,data,**kwargs):
		x = kwargs.get('x',self.x); y = kwargs.get('y',self.y)
		x = x if x in data.samples[0] else x.replace('.','-')
		y = y if y in data.samples[0] else y.replace('.','-')
		groups = kwargs.get('groups',self.groups); silentGroups = kwargs.get('silentGroups',self.silentGroups)
		individuals = kwargs.get('individuals',self.individuals)
		overKey=[x]; accKey=[x,y];
		if '_CI' in data.samples[0]:
			accKey.append('_CI')
		constant=groups+silentGroups+individuals  
#		timecourses, sampleIndex = data.accumulate(overKey=overKey,accKey=accKey+['dataTimepoint'],constant=constant,returnSampleIndex=True)
		timecourses, sampleIndex = data.accumulate(overKey=overKey,accKey=accKey,constant=constant,returnSampleIndex=True)
		return timecourses, sampleIndex
		
		
	def getGroupedSample(self,sample):
		newSample = {}
		for g, group in enumerate(self.groups):
			samplevalue = sample[group]
			if group in self.groupValues.keys():
				try:
					if not samplevalue in self.groupValues[group]:
						possibilities = numpy.array(self.groupValues[group])
						diffs = abs(possibilities - samplevalue)
						samplevalue = possibilities[numpy.argmin(diffs)]
				except:
					traceback.print_exc()
			
			newSample[group] = samplevalue
		return newSample
		
		
	def makeGroupDescriptor(self,sample):
		newSample = self.getGroupedSample(sample)
		try:
			group = '-'.join([str(newSample[k]).replace("'",'') for k in self.groups])
		except:
			print("Could not make group descriptor out of: %s" % ([str(newSample[k]).replace("'",'') for k in self.groups]))
			return None
		return group
		
		
	def getColorAndStyle(self,sample):
		style = self.linestyle
		color = self.colors[0]

		newSample = self.getGroupedSample(sample)
		groupDescription = self.makeGroupDescriptor(newSample)
		if groupDescription in self.groupFormatKey:
			color = self.groupFormatKey[groupDescription]['color']
			style = self.groupFormatKey[groupDescription]['style']
		else:
			for g,group in enumerate(self.groups[0:2]):			
				samplevalue = newSample[group]
				if g == 0:			# first group chooses colour
					color = self.groupColorKey.get(samplevalue,None)
					style = self.groupStyleKey.get(samplevalue,self.linestyle)
					if not color:
						self.colorCounter = (self.colorCounter + 1) % len(self.colors)
						color = self.colors[self.colorCounter]
						if color.__class__ == dict:
							self.groupStyleKey[samplevalue] = color.get('style',self.linestyle); style = self.groupStyleKey[samplevalue]
							self.groupColorKey[samplevalue] = color['color']; color = color['color']
						else:
							self.groupColorKey[samplevalue] = self.colors[self.colorCounter]
							color = self.groupColorKey[samplevalue]
				if g == 1:
					style = self.groupStyleKey.get(samplevalue,None)
					if style is None:		# Try composite group
						try:
							style = self.groupStyleKey.get('-'.join(self.groups[1:])).get('-'.join(groupDescription.split('-')[1:]))						
						except:
							pass
					if style is None:		# give up and make a new one
						self.styleCounter = (self.styleCounter + 1) % len(self.styles)
						style = self.styles[self.styleCounter]
						if style.__class__ == dict:
							self.groupColorKey[samplevalue] = style.get('color',self.colors[0]); self.groupColorKey[samplevalue]
							self.groupStyleKey[samplevalue] = style['style']; style = style['style']
						else:
							self.groupStyleKey[samplevalue] = self.styles[self.styleCounter]
							style = self.groupStyleKey.get(samplevalue,None)
			
		if not groupDescription in self.groupFormatKey:
			self.groupFormatKey[groupDescription] = {'color':color,'style':style}
		return color,style
		
		
	def substituteValues(self,alternateModel,estimateData=None):
		""" Replaces values for matching x's in this plotter with ones from the alternate model"""
		estimateData = estimateData if estimateData is not None else self.estimateData
		aplotter = ModelPlotter(model=alternateModel,x=self.x,y=self.y,groups=self.groups,individuals=self.individuals,silentGroups=self.silentGroups,fixedValues=self.fixedValues,groupValues=self.groupValues,mean=True,estimateData=estimateData)
		aplotter.plot(prepareOnly=True)
		self.plot(prepareOnly=True)
		et1 = aplotter.estimatedTimecourses
		et2 = self.estimatedTimecourses
		for t,tc in enumerate(et2):
			for i,x in enumerate(tc[self.x]):
				if x in et1[t][aplotter.x]:
					idx = et1[t]['elapsed'].index(x)
					for k in ['_CI',self.y]:
						tc[k][i] = et1[t][k][idx]
		
		
	def sorter(self,item):
		value = 0
		for k,key in enumerate(self.groupOrdering[::-1]):
			group = key[0]; order = key[1]
			try:
				value += k*100 + order.index(item[group])
			except:
				pass
			return value
	
	
	def plotEstimatedTimecourses(self,estimatedTimecourses,**kwargs):
		x = kwargs.get('x',self.x); y = kwargs.get('y',self.y)
		colors = kwargs.get('colors',None); styles = kwargs.get('styles',None)
		
		if self.groupOrdering is not None:
			estimatedTimecourses.sort(key=self.sorter)				
		
		for n,i in enumerate(estimatedTimecourses):
			group = self.makeGroupDescriptor(i)
			if colors and styles:
				color = colors[n]
				style = styles[n]
			else:
				color,style = self.getColorAndStyle(i)
			try:
				groupLabel = group
				if self.includeN:
					groupLabel += ' N=%d' % self.groupedData[group]['N']
				label = self.label if self.label is not None else groupLabel
			except:
				label = group
			xs = i[x]; ys = i[y]
			if self.xUnitConversion:
				xs = eval('x ' + self.xUnitConversion,{'x':numpy.array(xs)},globals())
			if self.yUnitConversion:
				ys = eval('y ' + self.yUnitConversion,{'y':numpy.array(ys)},globals())						
			if not i['_CI'][0] == None:
				pred = numpy.array(ys); ci = i['_CI']
				if self.yUnitConversion:
					ci = eval('x ' + self.yUnitConversion,{'x':numpy.array(ci)},globals())
				pylab.fill_between(xs,pred-ci,pred+ci,color=color,alpha=self.cialpha)
			pylab.plot(xs,ys,color=color,alpha=self.malpha,label=label,linewidth=self.mLinewidth,linestyle=style)
			
			
	def plotEstimates(self,estimates,**kwargs):
		estimatedTimecourses,_ = self.getTimecourses(estimates,**kwargs)
		self.plotEstimatedTimecourses(estimatedTimecourses,**kwargs)
		
		
	def plotData(self,prepareOnly=False,correctFixed={},correctRandom={},addResiduals=True,**kwargs):
		if self.model is not None and not 'data' in kwargs:			# If we were supplied data, use that
			self.data = self.model.getIndividualEstimates(fixedTerms=correctFixed, excludeRandomEffects=correctRandom,addResiduals=addResiduals, allData=True)
		
		s,self.sampleIndex = self.getTimecourses(self.data,**kwargs)
		self.timecourses = s
		self.groupedData = {}
		if self.ialpha == None:
			ialpha = min(1.,1. / len(s) * 4.5)
		else:
			ialpha = self.ialpha
		for i in s:
			group = self.makeGroupDescriptor(i)
			if kwargs.get('individualColor',None):
				color = kwargs['individualColor']
				style = kwargs.get('individualStyle','-')
			else:
				color,style = self.getColorAndStyle(i)

			if not group in self.groupedData:
				self.groupedData[group] = {'x':[],'y':[],'weights':[],'N':0}
				if not self.mean:
					label = group
				else:
					label = None
				label = None
			else:
				label = None
			if self.plotIndividuals and not prepareOnly:
				try:
					xs = i[self.x]; ys = i[self.y]
					if self.xUnitConversion:
						xs = eval('x ' + self.xUnitConversion,{'x':numpy.array(xs)},globals())
					if self.yUnitConversion:
						ys = eval('x ' + self.yUnitConversion,{'x':numpy.array(ys)},globals())
					if not self.connect:
						style = ''
					pylab.plot(xs,ys,color=color,label=label,alpha=ialpha,marker=self.marker,markersize=self.markerSize,linewidth=self.iLinewidth,linestyle=style)
				except:
					traceback.print_exc()
			ix = ma.array(i.get(self.x,[None])); ix.mask = numpy.equal(ix,None)
			iy = ma.array(i.get(self.x,[None])); iy.mask = numpy.equal(ix,None)
			combinedMask = numpy.logical_or(ix.mask,iy.mask)
			ix.mask = combinedMask; iy.mask = combinedMask
			if (~combinedMask).any():
				self.groupedData[group]['x'] += i[self.x]
				self.groupedData[group]['y'] += i[self.y]
				if isinstance(self.model,MixedModel):
					self.groupedData[group]['N'] += 1
				else:
					self.groupedData[group]['N'] = len(i[self.y])
				if self.weights:
					self.groupedData[group]['weights'] += i[self.weights]
		self.individualLines = pylab.gca().lines
		
		
	def plotResiduals(self,x,**kwargs):
		residuals = self.model.getResiduals()		
		samples = []; counter = 0
		for s,sample in enumerate(self.model.data.samples):
			if all([sample.get(term,None) is not None for term in list(self.model.atomicTerms) + [self.model.y]]):
				samples.append({x:sample[x],'residual':residuals[counter]})
				counter += 1
		
		data = ModelData(samples)
		plotter = ModelPlotter(data=data,x=x,y='residual',marker='o',connect=False,**kwargs)
		plotter.plot()
		pylab.ylabel('Residual');
		pylab.gca().legend_ = None
		pylab.draw()

		
		
			
	def plot(self,**kwargs):
		prepareOnly = kwargs.get('prepareOnly',False)
		for k in list(kwargs.keys()):
			if hasattr(self,k):
				setattr(self,k,kwargs.pop(k))
		if 'figure' in kwargs:
			figure(figure)
		
		# Plot the individuals
		self.plotData(**kwargs)

		# Plot the mean
		if self.mean:
			if self.meanStyle == 'model' or (self.meanStyle == None and self.model and self.estimateData):
				# Model mean
				if not self.estimatedTimecourses:
					if self.estimateData == None:
						knots = kwargs.get('knots',None)
						if knots is None:
							mn = numpy.ma.min([numpy.ma.min(self.groupedData[group]['x']) for group in self.groupedData.keys()])
							mx = numpy.ma.max([numpy.ma.max(self.groupedData[group]['x']) for group in self.groupedData.keys()])
							knots = numpy.arange(mn,mx,(mx-mn)/10.)
						self.estimateData = ModelData()
						individual = (self.individuals[0],'subject1') if not self.individuals is None and len(self.individuals) > 0 else None
						self.estimateData.makeEstimate(self.y,self.x,knots,individual=individual,extraKeys=self.data.samples[0].keys())

					if self.derivedValues:
						self.estimateData.calculateDerivedValues(self.derivedValues)

					self.estimates = self.model.getEstimates(self.estimateData,self.groups,self.y,self.groupValues,self.fixedValues)
					y = [sample[self.y] for sample in self.estimates.samples]
					if numpy.any(numpy.isnan(y)):
						print("There are nans in plotter's estimate data. Probably you forgot to specify a group or group value")
					if not self.estimates is None:
						self.estimatedTimecourses,_ = self.getTimecourses(self.estimates,**kwargs)
					else:
						print('Plotter could not create estimates. Did you forget to specify something?')
						return
				if prepareOnly:
					return
				
				self.plotEstimatedTimecourses(self.estimatedTimecourses,**kwargs)
				
			if (self.meanStyle.lower() == 'loess' or (self.meanStyle == None and self.model == None)):
				# LOESS mean
				self.loess = {}
				for group in self.groupedData.keys():
					try:
						loess = Loess(None,None,span=self.span)
						loess.setData(x=self.groupedData[group]['x'],y=self.groupedData[group]['y'],weights=self.groupedData[group]['weights'])
						mn = min(self.groupedData[group]['x']); mx = max(self.groupedData[group]['x'])
						xs = numpy.arange(mn,mx,(mx-mn)/1000.)
						pred,ci = loess.predict(xs)

						if self.xUnitConversion:
							xs = eval('x ' + self.xUnitConversion,{'x':numpy.array(xs)},globals())
						if self.yUnitConversion:
							pred = eval('x ' + self.yUnitConversion,{'x':numpy.array(pred)},globals())						
							ci = eval('x ' + self.yUnitConversion,{'x':numpy.array(ci)},globals())						
						
						label = group;
						if self.includeN:
							label += ' N=%d' % self.groupedData[group]['N']
						pylab.plot(xs,pred,color=self.groupFormatKey[group]['color'],alpha=self.malpha,label=label,linewidth=self.mLinewidth,linestyle=self.groupFormatKey[group]['style'])
						pylab.fill_between(xs,pred-ci,pred+ci,color=self.groupFormatKey[group]['color'],alpha=self.cialpha)
						self.loess[group] = loess
					except:
						print("Error building loess curve for %s" % group)
						traceback.print_exc()
						print()

			if self.meanStyle == 'binned':
				# Binned mean
				self.binnedMean = {}
				if not self.binFunction:
					delta = 1000000000000000
					for group in self.groupedData.keys():
						x = numpy.array(self.groupedData[group]['x'])
						delta = min([1. / ((max(x)-min(x)) / 100.),delta])
					self.binFunction = "numpy.around(x*%f).astype(numpy.int16) / %f" % (delta,delta)
				for group in self.groupedData.keys():
					x = numpy.array(self.groupedData[group]['x'])
					y = numpy.array(self.groupedData[group]['y'])
					binnedX = eval(self.binFunction,{'x':x},globals())
					newY = []
					newX = []
					ci = []
					N = []
					for i in numpy.unique(binnedX):
						vals = numpy.ma.masked_invalid(y[numpy.where(binnedX==i)[0]])
						newX.append(i)
						newY.append(numpy.ma.average(vals))
						if self.errorBars == '95':
							ci.append(numpy.ma.std(vals)/numpy.ma.sqrt(numpy.ma.count(vals))*1.96)
						elif self.errorBars == 'sd':
							ci.append(numpy.ma.std(vals))
						N.append(numpy.ma.count(vals))
					newY = numpy.array(newY); ci = numpy.array(ci)
					label = group 
					if self.includeN:
						label += ' N=%d' % self.groupedData[group]['N']

					if self.xUnitConversion:
						newX = eval('x ' + self.xUnitConversion,{'x':numpy.array(newX)},globals())
					if self.yUnitConversion:
						newY = eval('x ' + self.yUnitConversion,{'x':numpy.array(newY)},globals())						
						ci = eval('x ' + self.yUnitConversion,{'x':numpy.array(ci)},globals())						

					pylab.plot(newX,newY,color=self.groupFormatKey[group]['color'],alpha=self.malpha,label=label,linewidth=self.mLinewidth,linestyle=self.groupFormatKey[group]['style'])
					pylab.fill_between(newX,newY-ci,newY+ci,color=self.groupFormatKey[group]['color'],alpha=self.cialpha)
					self.binnedMean[group] = dict(x=newX,y=newY,ci=ci)
			
		pylab.xlabel(self.x)
		pylab.ylabel(self.y)
		if self.drawLegend:
			pylab.legend()
		

StatPlotter = ModelPlotter

