import sys
import os
import subprocess
import numpy as np
import pandas as pd
import pylab
import matplotlib.pyplot as plt
import rpy2
import rpy2.robjects as robjects
import scipy.stats as stats
from scipy.optimize import curve_fit


class Mk:
	
	def __init__(self,sfs,div,table=False):
	
		self.sfs = sfs
		self.div = div
		self.table = table
		self.p0 = sfs.p0.sum()
		self.pi = sfs.pi.sum()
		self.d0 = div.d0[0]
		self.di = div.di[0]
		self.m0 = div.m0[0]
		self.mi = div.mi[0]

	data = np.zeros([4,1], dtype={'names':('name', 'age', 'weight'),'formats':('U10', 'i4', 'f8')})

	def createTableFull(self):
		## Create MKT table 
		output = pd.DataFrame([[self.pi,self.di,self.mi],[self.p0,self.d0,self.m0]],columns = ['polymorphism','divergence','sites'],index=['selected','neutral'])
		return(output)

	def createTableFull(self):
		## Create MKT table 
		output = pd.DataFrame([[self.pi,self.di,self.mi],[self.p0,self.d0,self.m0]],columns = ['polymorphism','divergence','sites'],index=['selected','neutral'])
		return(output)

	def standard(self):
		alpha = 1 - ((self.d0*self.pi)/(self.p0*self.di))
		oddsRatio, pvalue = stats.fisher_exact([[self.p0, self.pi ], [self.d0,self.di]])
		output = np.array((alpha,pvalue), dtype={'names':('alpha', 'pvalue',),'formats':('f8','f8')})
		# output = pd.DataFrame(alpha,pvalue, columns=['alpha','pvalue'],index=1)
		return(output)

	def emkt(self,listCutoffs,plot=False):

		## Declare output lists and data frames
		output     = list()
				
		## Estimation of alpha
		alphaCorrected = list()
		fractions = list()

		for c in listCutoffs:
			
			## Estimating alpha with pi/p0 ratio 
			piMinus     = self.sfs[self.sfs['freq'] <= c].pi.sum() 
			piGreater   = self.sfs[(self.sfs['freq'] > c) & (self.sfs['freq'] < 0.9)].pi.sum() 
			p0Minus     = self.sfs[self.sfs['freq'] <= c].p0.sum() 
			p0Greater   = self.sfs[(self.sfs['freq'] > c) & (self.sfs['freq'] < 0.9)].p0.sum() 
			
			ratiop0     = p0Minus/p0Greater
			deleterious = piMinus - (piGreater * ratiop0)
			piNeutral   = self.pi - deleterious
			
			alphaC      = 1 - (((self.pi - deleterious)/self.p0)*(self.d0/self.d0))

			## Estimation of b: weakly deleterious
			b    = (deleterious/self.d0)*(self.m0/self.mi)
			
			## Estimation of f: neutral sites
			f = (self.m0*piNeutral)/(self.mi*self.p0)
			
			## Estimation of d, strongly deleterious sites
			d = 1 - (f+b)
			
			## Fisher exact test p-value from the MKT
			oddsRatio, pvalue = stats.fisher_exact([[self.p0, (self.pi -deleterious)], [self.d0,self.di]])
						
			## Store output  
			alphaCorrected.append([alphaC, pvalue])
			fractions.append([d, f, b])

		## Fractions
		# names(fraction) = 'Correction'
		
		## Store output 
		## Output format
		alphaCorrected        = pd.DataFrame(alphaCorrected,columns = ['alphaCorrected','pvalue'],index=[listCutoffs])
		# colnames(output[['alphaCorrected']]) = c('cutoff', 'alphaCorrected', 'pvalue')
		fractions = pd.DataFrame(fractions,columns = ['d','f','b'],index=[listCutoffs])
		# names(output[['fractions']])                     = c('cutoff', 'd','f','b')
		
		## Render plot
		if plot:

			## Cut-offs graph
		
			## Melt fractions data
			def plot(cutoffs,alphaValues):
				names = [str(x) for x in cutoffs]    
				fig = plt.figure(figsize=(8, 6))
				plt.style.use('seaborn')
				plt.plot(names, alphaValues, 'o',color='black');
				plt.plot(names, alphaValues, color='black');
				plt.ylim((alphaValues[0] - alphaValues[0]*50/100), (alphaValues[-1] + alphaValues[-1]*50/100));
				plt.xlabel('Frequency cut-offs')
				plt.ylabel('Adaptation rate')
				return(fig)

			ax = plot(listCutoffs,alphaCorrected.alphaCorrected)
			
			return(alphaCorrected,fractions,ax)

		## If no plot to render
		else:

			return(alphaCorrected,fractions)

	def amkt(self,xlow=0.1,xhigh=0.9,plot=False):

		## Compute alpha values and trim
		alpha         = 1 - ((self.sfs.pi*self.d0)/(self.sfs.p0*self.di))
		cutoff1       = xlow
		cutoff2       = xhigh
		trim          = ((self.sfs['freq'].values >= cutoff1) & (self.sfs['freq'].values <= cutoff2))
		fTrimmed     = self.sfs['freq'][trim]
		alphaTrimmed = alpha[trim]
		
		## Compute the original MK alpha using trimmed values
		alphaNonasymp = (1 - ((self.sfs.pi.sum()*self.d0)/(self.sfs.p0.sum()*self.di)))

		fitMKmodel="""
			function(alpha_trimmed, f_trimmed, res) {
				suppressMessages(library(nls2))

				## First fitting using starting values (st)
				mod = tryCatch({
					
					## Starting values to fit the model  
					st = expand.grid(const_a=seq(-1,1,length.out=res + 1), const_b=seq(-1,1,length.out=res), const_c=seq(1,10,length.out=res + 1))
					
					## Fitting
					nls2::nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start=st, algorithm='brute-force', control=nls.control(maxiter=NROW(st)))
					
				}, error=function(cond) {}) ## Return condition of error when unable to fit
				
				## If mod fails...
				if (length(mod) == 0) { return(NULL) }
				
				## Second fitting, starting from previous fit (mod)
				mod2 = tryCatch({
					nls2::nls2(alpha_trimmed ~ const_a + const_b * exp(-const_c* f_trimmed), start = mod, control=nls.control(maxiter=200))
					
				}, error=function(cond) {}) ## Same error handling than the previous step
				
				## If mod2 fails...
				if (length(mod2) == 0) { return(NULL) }
				
				## Return mod2 if fitted
				return(mod2)
			}
		"""
		predictNLS = """
			## Compute confidence intervals of alpha using predictNLS 
			## Get a CI using Monte Carlo simulation based upon a fitted model.  
			## Thanks to Andrej-Nikolai Spiess (http://www.dr-spiess.de) for this code.
			function(object, newdata, level = 0.95, nsim = 5000) {
				suppressMessages(library(MASS))
				## get right-hand side of formula
				RHS  = as.list(object$call$formula)[[3]]
				EXPR = as.expression(RHS)
				
				## all variables in model
				VARS = all.vars(EXPR)
				
				## coefficients
				COEF = coef(object)
				
				## extract predictor variable    
				predNAME = setdiff(VARS, names(COEF))  
				
				## take fitted values, if 'newdata' is missing
				if (missing(newdata)) {
					newdata           = eval(object$data)[predNAME]
					colnames(newdata) = predNAME
				}
				  
				## check that 'newdata' has same name as predVAR
				if (names(newdata)[1] != predNAME) stop("newdata should have name '", predNAME, "'!")

				## get parameter coefficients
				COEF = coef(object)

				## get variance-covariance matrix
				VCOV = vcov(object)

				## augment variance-covariance matrix for 'mvrnorm' 
				## by adding a column/row for 'error in x'
				NCOL = ncol(VCOV)
				ADD1 = c(rep(0, NCOL))
				ADD1 = matrix(ADD1, ncol = 1)
				colnames(ADD1) = predNAME
				VCOV = cbind(VCOV, ADD1)
				ADD2 = c(rep(0, NCOL + 1))
				ADD2 = matrix(ADD2, nrow = 1)
				rownames(ADD2) = predNAME
				VCOV = rbind(VCOV, ADD2) 

				## iterate over all entries in 'newdata' as in usual 'predict.' functions
				NR = nrow(newdata)
				respVEC = numeric(NR)
				seVEC = numeric(NR)
				varPLACE = ncol(VCOV)   

				## define counter function
				counter = function(i) {
					if (i%%10 == 0) { cat(i) 
				} else { cat(".") }
					if (i%%50 == 0) { cat("\n") }
					flush.console()
				}
				  
				## create output matrix (df)
				outMAT = NULL 
				
				for (i in 1:NR) {
					
					## get predictor values and optional errors
					predVAL = newdata[i, 1]
					if (ncol(newdata) == 2) predERROR = newdata[i, 2] else predERROR = 0
					names(predVAL) = predNAME  
					names(predERROR) = predNAME  
					
					## create mean vector for 'mvrnorm'
					MU = c(COEF, predVAL)
					
					## create variance-covariance matrix for 'mvrnorm'
					## by putting error^2 in lower-right position of VCOV
					newVCOV = VCOV
					newVCOV[varPLACE, varPLACE] = predERROR^2
					
					## create MC simulation matrix
					simMAT = mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, empirical = TRUE)
					
					## evaluate expression on rows of simMAT
					EVAL = try(eval(EXPR, envir = as.data.frame(simMAT)), silent = TRUE)
					if (inherits(EVAL, "try-error")) stop("There was an error evaluating the simulations!")
					
					## collect statistics
					PRED = data.frame(predVAL)
					colnames(PRED) = predNAME   
					FITTED = predict(object, newdata = data.frame(PRED))
					MEAN.sim = mean(EVAL, na.rm = TRUE)
					SD.sim = sd(EVAL, na.rm = TRUE)
					MEDIAN.sim = median(EVAL, na.rm = TRUE)
					MAD.sim = mad(EVAL, na.rm = TRUE)
					QUANT = quantile(EVAL, c((1 - level)/2, level + (1 - level)/2))
					RES = c(FITTED, MEAN.sim, SD.sim, MEDIAN.sim, MAD.sim, QUANT[1], QUANT[2])
					outMAT = rbind(outMAT, RES)
				}
				  
				colnames(outMAT) = c("fit", "mean", "sd", "median", "mad", names(QUANT[1]), names(QUANT[2]))
				rownames(outMAT) = NULL   
				return(outMAT)      
			}
		"""

		fitMKmodel=robjects.r(fitMKmodel)
		predictNLS=robjects.r(predictNLS)

		predict=robjects.r('stats::predict')
		coef=robjects.r('stats::coef')
		rdf=robjects.r('base::data.frame')
		rdata1 = robjects.FloatVector(alphaTrimmed.values)
		rdata2 = robjects.FloatVector(fTrimmed.values)

		## Two-step nls2() model fit at a given level of precision (res)
		mod1 = fitMKmodel(rdata1,rdata2, 10)
		## If mod1 did not work, try a deeper scan for a decent fit (res=20)
		try:
			alpha_1_est = predict(mod1, newdata=rdf(f_trimmed=1))
			const = coef(mod1)
			const = dict(zip(const.names, list(const)))
			const_a = const['const_a']
			const_b = const['const_b']
			const_c = const['const_c']
			
			ci_pred = list(predictNLS(mod1, newdata=rdf(f_trimmed=1)))
			alpha_1_low = ci_pred[5]
			alpha_1_high = ci_pred[6]

			asympDf = pd.DataFrame({'model':'exponential', 'a':const_a, 'b':const_b, 'c':const_c, 'alphaAsymptotic':alpha_1_est, 'ciLow':alpha_1_low, 'ciHigh':alpha_1_high, 'alphaOriginal':alphaNonasymp},index=[0])

			if plot:
				# Plotting fit and alpha data
				def asympPlot(x):
					return const_a + const_b * np.exp(-const_c*x)
				vectorizeFunc = np.vectorize(asympPlot)
				x = np.arange(0.0, 1, 0.00001)
				y = vectorizeFunc(x)

				def plot(x,y,freq,alpha,ciLow,ciHigh):
					
					fig = plt.figure(figsize=(8, 6))
					plt.style.use('seaborn')
					plt.plot(x,y, color='black')
					plt.plot(x, y, color='black');
					plt.xlim(0 , 1);
					plt.ylim(y[0], 1);
					plt.plot(x, np.array([ciLow]*x.shape[0]), color='black');
					plt.plot(x, np.array([ciHigh]*x.shape[0]), color='black');
					plt.plot(freq, alpha, 'o', color='black');
					plt.xlabel('Allele Frequency')
					plt.ylabel('Rate of adaptation')
					plt.fill_between(x,np.array([ciLow]*x.shape[0]), np.array([ciHigh]*x.shape[0]),facecolor='gray', alpha=0.2)
					return(fig)

				ax = plot(x,y,fTrimmed,alphaTrimmed,alpha_1_low,alpha_1_high)

		# Force another fit
		except rpy2.rinterface.RRuntimeError:
			mod1 = fitMKmodel(rdata1,rdata2, 20)
			
			try:
				alpha_1_est = predict(mod1, newdata=rdf(f_trimmed=1))
				const = coef(mod1)
				const = dict(zip(const.names, list(const)))
				const_a = const['const_a']
				const_b = const['const_b']
				const_c = const['const_c']
				
				ci_pred = list(predictNLS(mod1, newdata=rdf(f_trimmed=1)))
				alpha_1_low = ci_pred[5]
				alpha_1_high = ci_pred[6]

				asympDf = pd.DataFrame({'model':'exponential', 'a':const_a, 'b':const_b, 'c':const_c, 'alphaAsymptotic':alpha_1_est, 'ciLow':alpha_1_low, 'ciHigh':alpha_1_high, 'alphaOriginal':alphaNonasymp},index=[0])
			except:
				asympDf = pd.DataFrame({'model':'exponential', 'a':np.nan, 'b':np.nan, 'c':np.nan, 'alphaAsymptotic':np.nan, 'ciLow':np.nan, 'ciHigh':np.nan, 'alphaOriginal':alphaNonasymp},index=[0])
		except:

			asympDf = pd.DataFrame({'model':'exponential', 'a':np.nan, 'b':np.nan, 'c':np.nan, 'alphaAsymptotic':np.nan, 'ciLow':np.nan, 'ciHigh':np.nan, 'alphaOriginal':alphaNonasymp},index=[0])

		if plot:
			return(asympDf,ax)
		else:
			return(asympDf,np.nan)

