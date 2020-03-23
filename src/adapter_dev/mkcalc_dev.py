# created by LHU
# 11/28/2016
# modified by jmurga
# 03/2020


from scipy.stats import binom
from scipy.optimize import fsolve
# from scipy.optimize import minimize
from scipy.integrate import quad
from scipy.special import gamma as gammaFunc
from mpmath import gammainc
# from scipy.stats import gamma as gamInt
from scipy.special import polygamma
import subprocess
import sys
import os
import mpmath
import math
import numpy as np
from numba import *

import julia
jl = julia.Julia(compiled_modules=False)
fn = jl.include('/home/jmurga/mktest/scripts/fullNeg.jl')
fn2 = jl.include('/home/jmurga/mktest/scripts/SfsPos.jl')

@njit
def binomOp_numba(NN,B,nn):

	NN = int(round(NN*B))
	samps2 = np.empty((51,501))
	sampFreqs2 = np.empty((51,501))

	for i in range(0,nn+1):
		for j in range(0,NN+1):
			samps2[i,j] = i
			sampFreqs2[i,j] = j/(NN+0.)

	return samps2,nn,sampFreqs2


class AsympMK:

	def __init__(self,gam_neg=-40,N=250,theta_f=1e-5,theta_mid_neutral=1e-3,
				 alLow=0.2,alTot=0.2,gL=10,gH=200,al=0.184,be=0.000402,B=1.,
				 pposL=0.001,pposH=0.0,n=10,Lf=10**6,rho=0.001,neut_mid=False,
				 L_mid=150,pref="",nsim=100,TE=5.,demog=False,ABC=False,al2=0.0415,
				 be2=0.00515625,gF=False,skip=0,expan=1.98,scratch=False):

		 self.skip = skip
		 self.gam_neg = gam_neg
		 self.N = N
		 self.NN = 2*N
		 self.theta_f = theta_f
		 self.theta_mid_neutral = theta_mid_neutral
		 self.alLow = alLow
		 self.alTot = alTot
		 self.gL = gL
		 self.gH = gH
		 self.al = al
		 self.be = be
		 self.al2 = al2
		 self.be2 = be2
		 self.B = B
		 self.task_id = 1
		 self.pposL = pposL
		 self.pposH = pposH
		 self.n = n
		 self.nn = 2*n
		 self.alpha_x = np.zeros(self.nn-1)
		 self.Lf = Lf 
		 self.L_mid = L_mid
		 self.pref = pref
		 self.nsim=nsim
		 self.gF = gF
		 self.cumuSfsNeut= []
		 self.cumuSfsSel = []
		 self.TE = TE
		 self.demog = demog
		 self.ABC=ABC
		 self.expan=expan
		 self.scratch=scratch
		 while self.L_mid % 3 != 0:
			 self.L_mid +=1 
		 self.rho = rho
		 self.neut_mid = neut_mid
		 if 'SGE_TASK_ID' in os.environ:
			 self.task_id = int(os.environ['SGE_TASK_ID'])
		 if 'SLURM_ARRAY_TASK_ID' in os.environ:
			 self.task_id = int(os.environ['SLURM_ARRAY_TASK_ID'])
		 self.task_id = self.task_id + skip

	def GamSfsNeg(self, x):
		beta = self.be/(1.*self.B)
		return (2.**-self.al)*(beta**self.al)*(-mpmath.zeta(self.al,x+beta/2.) + mpmath.zeta(self.al,(2+beta)/2.))/((-1.+x)*x)

	def SfsPos(self, gamma, x):
		gam = gamma*self.B
		return 0.5*(mpmath.exp(2*gam)*(1-mpmath.exp(-2.*gam*(1.-x)))/((mpmath.exp(2*gam)-1.)*x*(1.-x)))

	def FullSfs(self, gamma, ppos, x):
		if x > 0 and x < 1.:
			#S = abs(self.gam_neg/(1.*self.NN))
			#r = self.rho/(2.*self.NN)
			#u = self.theta_f/(2.*self.NN)
			#s = gamma/(self.NN*1.)
			#p0 = polygamma(1,(s+S)/r)
			#p1 = polygamma(1,1.+(r*self.Lf+s+S)/r)
			#pposplus = ppos*mpmath.exp(-2.*S*u*(p0-p1)/r**2)
			return ppos*self.SfsPos(gamma, x) + (1.-ppos)*self.GamSfsNeg(x)
		return 0.
 
	def FullPos(self, gamma, ppos, x):
		if x > 0 and x < 1.:
			return ppos*self.SfsPos(gamma, x)
		return 0.
 
	def FullNeg(self, ppos, x):
		beta = self.be/(1.*self.B)
		if x > 0 and x < 1.:
			return (1.-ppos)*(2.**-self.al)*(beta**self.al)*(-mpmath.zeta(self.al,x+beta/2.) + mpmath.zeta(self.al,(2+beta)/2.))/((-1.+x)*x)

		return 0.

	def binomOp(self):
		NN = int(round(self.NN*self.B))
		samps = [[i for j in range(0,NN+1)] for i in range(0,self.nn+1)]
		sampFreqs = [[j/(NN+0.) for j in range(0,NN+1)] for i in range(0,self.nn+1)]
		return binom.pmf(samps,self.nn,sampFreqs)

	def binomOp_numba(self):

		NN = int(round(self.NN*self.B))
		samps2 = np.empty((self.nn+1,self.NN+1))
		sampFreqs2 = np.empty((self.nn+1,self.NN+1))

		for i in range(0,self.nn+1):
			for j in range(0,self.NN+1):
				samps2[i,j] = i
				sampFreqs2[i,j] = j/(self.NN+0.)

		return samps2,self.nn,sampFreqs2


	def DiscSFSSel(self,gamma,ppos):
		NN = int(round(self.NN*self.B))
		# dFunc = np.vectorize(self.FullSfs)
		# return np.multiply(1./(NN+0.),dFunc(gamma,ppos,[i/(NN+0.) for i in range(0,NN+1)]))
		n = np.array([i/(NN+0.) for i in range(0,NN+1)])
		return np.multiply(1./(NN+0.),dFunc(gamma,ppos,[i/(NN+0.) for i in range(0,NN+1)]))
	
	def DiscSFSSelPos(self,gamma,ppos):
		NN = int(round(self.NN*self.B))
		# dFunc = np.vectorize(self.FullPos)
		# return np.multiply(1./(NN+0.),dFunc(gamma,ppos,[i/(NN+0.) for i in range(0,NN+1)]))
		n = np.array([i/(NN+0.) for i in range(0,NN+1)])
		return np.multiply(1./(NN+0.),fn2(gamma,self.B,n))
	
	def DiscSFSSelNeg(self,ppos):
		NN = int(round(self.NN*self.B))
		# dFunc = np.vectorize(self.FullNeg)
		beta = self.be/(1.*self.B)
		n = np.array([i/(NN+0.) for i in range(0,NN+1)])
		return np.multiply(1./(NN+0.),fn(ppos,beta,self.al,n))

	def DiscSFSNeut(self):
		NN = int(round(self.NN*self.B))
		def sfs(i):
			if i > 0 and i < NN:
				 return 1./(i+0.)
			return 0.
		sfs = np.vectorize(sfs)
		return sfs([(i+0.) for i in range(0,NN+1)])

	def DiscSFSSelDown(self, gamma, ppos):
		return self.DiscSFSSelPosDown(gamma,ppos)+self.DiscSFSSelNegDown(ppos)
	
	def DiscSFSSelPosDown(self, gamma, ppos):
		S = abs(self.gam_neg/(1.*self.NN))
		r = self.rho/(2.*self.NN)
		u = self.theta_f/(2.*self.NN)
		s = gamma/(self.NN*1.)
		p0 = polygamma(1,(s+S)/r)
		p1 = polygamma(1,1.+(r*self.Lf+s+S)/r)
		red_plus = mpmath.exp(-2.*S*u*(p0-p1)/r**2)

		return (self.theta_mid_neutral)*red_plus*0.745*(np.dot(self.binomOp(),self.DiscSFSSelPos(gamma,ppos)))[1:-1]
	
	def DiscSFSSelNegDown(self, ppos):
		return self.B*(self.theta_mid_neutral)*0.745*(np.dot(self.binomOp(),self.DiscSFSSelNeg(ppos)))[1:-1]

	def DiscSFSNeutDown(self):
		return self.B*(self.theta_mid_neutral)*0.255*(np.dot(self.binomOp(),self.DiscSFSNeut()))[1:-1]

	def fixNeut(self):
		return 0.255*(1./(self.B*self.NN))

	#def fixNeg(self,pfpos):
	#    return 0.745*(1-ppos)*(2**(-self.al))*(self.be**self.al)*(-mpmath.zeta(self.al,(2.+self.be)/2.)+mpmath.zeta(self.al,0.5*(2-1./self.NN+self.be)))
	#  0.745*(1-ppos)*(2**(-alpha))*(B**-alpha)*(beta**alpha)*(-mpmath.zeta(alpha,1.+beta/(2.*B))+mpmath.zeta(alpha,0.5*(2-1./NN+beta/B)))
	def fixNegB(self,ppos):
		return 0.745*(1-ppos)*(2**(-self.al))*(self.B**(-self.al))*(self.be**self.al)*(-mpmath.zeta(self.al,1.+self.be/(2.*self.B))+mpmath.zeta(self.al,0.5*(2-1./(self.N*self.B)+self.be/self.B)))

	def pFix(self,gamma):
		s = gamma/(self.NN+0.)
		pfix = (1.-mpmath.exp(-2.*s))/(1.-mpmath.exp(-2.*gamma))
		if s >= 0.1:
			pfix = mpmath.exp(-(1.+s))
			lim = 0
			while(lim < 200):
				pfix = mpmath.exp((1.+s)*(pfix-1.))
				lim +=1
			pfix = 1-pfix
		return pfix
		
	def fixPosSim(self,gamma,ppos):

		S = abs(self.gam_neg/(1.*self.NN))
		r = self.rho/(2.*self.NN)
		u = self.theta_f/(2.*self.NN)
		s = gamma/(self.NN*1.)
		p0 = polygamma(1,(s+S)/r)
		p1 = polygamma(1,(r+self.Lf*r+s+S)/r)
		#CC = 2*s/self.pFix(gamma)
		CC = 1.
		#print mpmath.exp(-2.*S*u*(p0-p1)/r**2)
		return 0.745*ppos*mpmath.exp(-2.*S*u*(p0-p1)*CC**2/r**2)*self.pFix(gamma)

	def alphaExpSimLow(self,pposL,pposH):
		return self.fixPosSim(self.gL,0.5*pposL)/(self.fixPosSim(self.gL,0.5*pposL)+self.fixPosSim(self.gH,0.5*pposH)+self.fixNegB(0.5*pposL+0.5*pposH))

	def alphaExpSimTot(self,pposL,pposH):
		return (self.fixPosSim(self.gL,0.5*pposL)+self.fixPosSim(self.gH,0.5*pposH))/(self.fixPosSim(self.gL,0.5*pposL)+self.fixPosSim(self.gH,0.5*pposH)+self.fixNegB(0.5*pposL+0.5*pposH))

	def solvEqns(self,params):
		pposL,pposH = params
		return (self.alphaExpSimTot(pposL,pposH)-self.alTot,self.alphaExpSimLow(pposL,pposH)-self.alLow)

	def setPpos(self):
		pposL,pposH =  fsolve(self.solvEqns,(0.001,0.001))
		#print self.alphaExpSimLow(pposL,pposH)
		#print self.alphaExpSimTot(pposL,pposH)
		if pposL < 0.:
			pposL = 0.
		if pposH < 0.:
			pposH = 0.
		self.pposH = pposH
		self.pposL = pposL

	def GammaDist(self,gamma):
		return ((self.be**self.al)/gammaFunc(self.al))*(gamma**(self.al-1))*mpmath.exp(-self.be*gamma)

	def PiP0(self,gamma):
		U = 4*self.theta_f*self.Lf/(2.*self.NN)
		R = 2*self.Lf*self.rho/(2.*self.NN)
		return self.GammaDist(gamma)*mpmath.exp(-(self.GammaDist(gamma)*U/(2.*self.NN))/(gamma/(self.NN+0.)+R/(2.*self.NN)))

	def intPiP0(self):
		ret = lambda gam: self.PiP0(gam)
		return quad(ret,0.,1000)[0]

	def calcBGam(self,L,alpha,beta,theta):
		# not sure whats going on here, seems to overestimate or under sometimes
		u = 2*theta/(2.*self.NN)
		r = self.rho/(2.*self.NN)

		a = -1.*(2**(alpha+1))*np.exp(L*self.NN*r*beta)*L*self.N*((1./(L*self.N*r))**(1-alpha))*u*(beta**alpha)
		b =  float(gammainc(1-alpha,L*self.NN*r*beta))
		#fudge = 1-gamInt.cdf(0.2,a=self.al2,scale=1./self.be2)
		#fudge = 1.
		fudge = 0.25

		c = np.exp(a*b*fudge)
 
		return c
		#return np.frompyfunc(mpmath.exp(-1.*(2**(1+alpha))*mpmath.exp(L*self.NN*r*beta)*L*self.N*((1./(L*self.N*r))**(1-alpha))*u*(beta**alpha)*float(gammainc(1-alpha,L*self.NN*r*beta))),1,1)

	def calcB(self,L,theta):
		t = -1.*self.gam_neg/(self.NN+0.)
		u = 2.*theta/(2.*self.NN)
		r = self.rho/(2.*self.NN)

		#return u*t/((t+r*L)**2)
		#return u*t/(2*(t+(1.-np.exp(-2*r*L))/2.)**2)
		
		# Nordborg models
		#####return u/(2*r*t*(1+((1.-np.exp(-2*r*L))/2.)*(1-t)/t)**2)
		#return u/(t*(1+((1-np.exp(-2.*r*L))/2.)*(1-t)/t)**2)
		#return u/(t*(1+r*L*(1-t)/t)**2)
		#####return u/(t*(1+(np.exp(-2*r*L)/2)/t)**2)

	def Br(self,Lmax,theta):
		t = -1.*self.gam_neg/(self.NN+0.)
		u = theta/(2.*self.NN)
		r = self.rho/(2.*self.NN)
		#return np.exp(-2.*quad(lambda L: self.calcB(L,theta), 0., Lmax)[0])
		#print np.exp(-2.*t*u*(polygamma(1,(r+t)/r)-polygamma(1,1+Lmax+t/r))/r**2), np.exp(-u*2.*Lmax/(Lmax*r+2*t))
		return np.exp(-4*u*Lmax/(2*Lmax*r+t))

	def get_B_vals(self):
		ret = []
		for i in range(20,50):
			L = int(round(1.3**i))
			ret.append(self.Br(L),L)
		return ret
	
	def set_Lf(self):
		Lf  = fsolve(lambda L: self.Br(L,self.theta_f)-self.B,100)
		self.Lf = int(round(Lf[0]))

	def set_theta_f(self):
		theta_f  = fsolve(lambda theta: self.Br(self.Lf,theta)-self.B,0.00001)
		self.theta_f = theta_f[0]
	
	def set_theta_f_gam(self):
		theta_f  = fsolve(lambda theta: self.calcBGam(self.Lf,self.al2,self.be2,theta)-self.B,0.0000001)
		self.theta_f = theta_f[0]

	def cumuSfs(self,sfsTemp):
		out = [np.sum(sfsTemp)]
		for i in range(0,len(sfsTemp)):
			app = out[i]-sfsTemp[i]
			if app > 0.:
				out.append(out[i]-sfsTemp[i])
			else:
				out.append(0.) 
		return out

	def alx(self,gammaL,gammaH,pposL,pposH):
		ret = []

		#Fixation
		fN = self.B*self.fixNeut()
		fNeg = self.B*self.fixNegB(0.5*pposH+0.5*pposL)
		fPosL = self.fixPosSim(gammaL,0.5*pposL)
		fPosH = self.fixPosSim(gammaH,0.5*pposH)

		#Pol
		neut = self.cumuSfs(self.DiscSFSNeutDown())
		selH = self.cumuSfs(self.DiscSFSSelPosDown(gammaH,pposH))
		selL = self.cumuSfs(self.DiscSFSSelPosDown(gammaL,pposL))
		selN = self.cumuSfs(self.DiscSFSSelNegDown(pposH+pposL))
		
		sel = []
		for i in range(0,len(selH)):
			sel.append((selH[i]+selL[i])+selN[i])
		for i in range(0,self.nn-1):
			ret.append(float(1. - (fN/(fPosL + fPosH+  fNeg+0.))* sel[i]/neut[i]))
		return (ret,sel,neut,selH,selL,selN,fN,fNeg,fPosL,fPosH)
	
	def alx_noCumu(self,gammaL,gammaH,pposL,pposH):
		ret = []
		fN = self.B*self.fixNeut()
		fNeg = self.B*self.fixNegB(0.5*pposH+0.5*pposL)
		fPosL = self.fixPosSim(gammaL,0.5*pposL)
		fPosH = self.fixPosSim(gammaH,0.5*pposH)
		neut = self.DiscSFSNeutDown()
		selH = self.DiscSFSSelPosDown(gammaH,pposH)
		selL = self.DiscSFSSelPosDown(gammaL,pposL)
		selN = self.DiscSFSSelNegDown(pposH+pposL)        

		sel = []
		for i in range(0,len(selH)):
			sel.append((selH[i]+selL[i])+selN[i])
		for i in range(0,self.nn-1):
			ret.append(float(1. - (fN/(fPosL + fPosH+  fNeg+0.))* sel[i]/neut[i]))
		return ret

	def alx_nopos(self,gammaL,gammaH,pposL,pposH):
		# mean = (alpha/beta)
		# sdv  = alpha/beta^2
		# nal/nbe  = (alpha/beta)/2
		# nal/nbe^2 = sqrt(alpha)/beta^2/2
		ret = []
		fN = self.B*self.fixNeut()*(self.theta_mid_neutral/2.)*self.TE*self.NN
		fNeg = self.B*(self.theta_mid_neutral/2.)*self.TE*self.NN*self.fixNegB(0.5*pposH+0.5*pposL)
		fPosL = self.fixPosSim(gammaL,0.5*pposL)*(self.theta_mid_neutral/2.)*self.TE*self.NN
		fPosH = self.fixPosSim(gammaH,0.5*pposH)*(self.theta_mid_neutral/2.)*self.TE*self.NN
		
		neut = self.cumuSfs(self.DiscSFSNeutDown())
		selN = self.cumuSfs(self.DiscSFSSelNegDown(pposL+pposH))

		sel = []
		for i in range(0,len(selN)):
			sel.append(selN[i])
		for i in range(0,self.nn-1):
			ret.append(float(1. - (fN/(fPosL + fPosH+  fNeg+0.))* sel[i]/neut[i]))
		return ret
	
	def alx_nopos_noCumu(self,gammaL,gammaH,pposL,pposH):
		# mean = (alpha/beta)
		# sdv  = alpha/beta^2
		# nal/nbe  = (alpha/beta)/2
		# nal/nbe^2 = sqrt(alpha)/beta^2/2
		ret = []
		fN = self.B*self.fixNeut()*(self.theta_mid_neutral/2.)*self.TE*self.NN
		fNeg = self.B*(self.theta_mid_neutral/2.)*self.TE*self.NN*self.fixNegB(0.5*pposH+0.5*pposL)
		fPosL = self.fixPosSim(gammaL,0.5*pposL)*(self.theta_mid_neutral/2.)*self.TE*self.NN
		fPosH = self.fixPosSim(gammaH,0.5*pposH)*(self.theta_mid_neutral/2.)*self.TE*self.NN
		
		neut = self.DiscSFSNeutDown()
		selN = self.DiscSFSSelNegDown(pposL+pposH)

		sel = []
		for i in range(0,len(selN)):
			sel.append(selN[i])
		for i in range(0,self.nn-1):
			ret.append(float(1. - (fN/(fPosL + fPosH+  fNeg+0.))* sel[i]/neut[i]))
		return ret

	def makeDf(self,alphaPos,alphaNonPos):
		data = {'alphaPos':alphaPos,'alphaNonPos':alphaNonPos}
		df = pd.DataFrame(data)
		df = df.reset_index().rename(columns={'index':'counts'})
		df.counts = list(range(1,self.nn))
		return(df)



