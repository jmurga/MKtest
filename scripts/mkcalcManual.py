import sys
import os
import mpmath
import math
import numpy as np
from scipy.stats import binom
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.integrate import quad
from scipy.special import gamma as gammaFunc
from mpmath import gammainc
from scipy.stats import gamma as gamInt
from scipy.special import polygamma
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
from numba import *
from numba.extending import overload

B=0.999
gam_neg=-83
theta_f=5.417709306354643e-07
alLow=0.2
alTot=0.2
neut_mid=False
L_mid=501
Lf=10**6
N=500
NN=1000
n=250
nn=500
gL=10
gH=500
theta_mid_neutral=0.001
al=0.184
be=0.000402
rho=0.001
alpha_x = np.zeros(nn-1)
cumuSfsNeut= []
cumuSfsSel = []
pposL=0.005841064021461007
pposH=0.0

############################

def Br_manual(Lmax,theta):
	t = -1.*gam_neg/(NN+0.)
	u = theta/(2.*NN)
	r = rho/(2.*NN)
	return np.exp(-4*u*Lmax/(2*Lmax*r+t))

def set_theta_f_manual():
	theta_f  = fsolve(lambda theta: Br_manual(Lf,theta)-B,0.00001)
	theta_f = theta_f[0]

def fixNeut_manual():
	return 0.255*(1./(B*NN))

def fixNegB_manual(ppos):
	return 0.745*(1-ppos)*(2**(-al))*(B**(-al))*(be**al)*(-mpmath.zeta(al,1.+be/(2.*B))+mpmath.zeta(al,0.5*(2-1./(N*B)+be/B)))

def fixPosSim_manual(gamma,ppos):

	S = abs(gam_neg/(1.*NN))
	r = rho/(2.*NN)
	u = theta_f/(2.*NN)
	s = gamma/(NN*1.)
	p0 = polygamma(1,(s+S)/r)
	p1 = polygamma(1,(r+Lf*r+s+S)/r)
	#CC = 2*s/pFix(gamma)
	CC = 1.
	#print mpmath.exp(-2.*S*u*(p0-p1)/r**2)
	return 0.745*ppos*mpmath.exp(-2.*S*u*(p0-p1)*CC**2/r**2)*pFix_manual(gamma)
def pFix_manual(gamma):
	s = gamma/(NN+0.)
	pfix = (1.-mpmath.exp(-2.*s))/(1.-mpmath.exp(-2.*gamma))
	if s >= 0.1:
		pfix = mpmath.exp(-(1.+s))
		lim = 0
		while(lim < 200):
			pfix = mpmath.exp((1.+s)*(pfix-1.))
			lim +=1
		pfix = 1-pfix
	return pfix

def cumuSfs_manual(sfsTemp):
	out = [np.sum(sfsTemp)]
	for i in range(0,len(sfsTemp)):
		app = out[i]-sfsTemp[i]
		if app > 0.:
			out.append(out[i]-sfsTemp[i])
		else:
			out.append(0.) 
	return out

def DiscSFSNeutDown_manual():
	return B*(theta_mid_neutral)*0.255*(np.dot(binomOp_manual(),DiscSFSNeut_manual()))[1:-1]

def binomOp_manual():
	NN2 = int(round(NN*B))
	samps = [[i for j in range(0,NN2+1)] for i in range(0,nn+1)]
	sampFreqs = [[j/(NN2+0.) for j in range(0,NN2+1)] for i in range(0,nn+1)]
	return binom.pmf(samps,nn,sampFreqs)

def DiscSFSNeut_manual():
	NN2 = int(round(NN*B))
	def sfs(i):
		if i > 0 and i < NN2:
			 return 1./(i+0.)
		return 0.
	sfs = np.vectorize(sfs)
	return sfs([(i+0.) for i in range(0,NN2+1)])

def DiscSFSSelPosDown_manual( gamma, ppos):
	S = abs(gam_neg/(1.*NN))
	r = rho/(2.*NN)
	u = theta_f/(2.*NN)
	s = gamma/(NN*1.)
	p0 = polygamma(1,(s+S)/r)
	p1 = polygamma(1,1.+(r*Lf+s+S)/r)
	red_plus = mpmath.exp(-2.*S*u*(p0-p1)/r**2)

	return (theta_mid_neutral)*red_plus*0.745*(np.dot(binomOp_manual(),DiscSFSSelPos_manual(gamma,ppos)))[1:-1]

def DiscSFSSelPos_manual(gamma,ppos):
	NN2 = int(round(NN*B))
	dFunc = np.vectorize(FullPos_manual)
	return np.multiply(1./(NN2+0.),dFunc(gamma,ppos,[i/(NN2+0.) for i in range(0,NN2+1)]))

def SfsPos_manual( gamma, x):
	gam = gamma*B
	return 0.5*(mpmath.exp(2*gam)*(1-mpmath.exp(-2.*gam*(1.-x)))/((mpmath.exp(2*gam)-1.)*x*(1.-x)))
def FullPos_manual( gamma, ppos, x):
	if x > 0 and x < 1.:
		return ppos*SfsPos_manual(gamma, x)
	return 0.

def DiscSFSSelNegDown_manual(ppos):
		return B*(theta_mid_neutral)*0.745*(np.dot(binomOp_manual(),DiscSFSSelNeg_manual(ppos)))[1:-1]

def DiscSFSSelNegDown_manual2(ppos):
		return B*(theta_mid_neutral)*0.745*(np.dot(binomOp_manual(),DiscSFSSelNeg_manual2(ppos)))[1:-1]

def FullNeg_manual( ppos, x):
	beta = be/(1.*B)
	if x > 0 and x < 1.:
		return (1.-ppos)*(2.**-al)*(beta**al)*(-mpmath.zeta(al,x+beta/2.) + mpmath.zeta(al,(2+beta)/2.))/((-1.+x)*x)

	return 0.

def DiscSFSSelNeg_manual(ppos):
	NN2 = int(round(NN*B))
	dFunc = np.vectorize(FullNeg_manual)
	return np.multiply(1./(NN2+0.),dFunc(ppos,[i/(NN2+0.) for i in range(0,NN2+1)]))

def DiscSFSSelNeg_manual2(ppos):
	NN2 = int(round(NN*B))
	# dFunc = np.vectorize(FullNeg_manual)
	beta = be/(1.*B)
	n = np.array([i/(NN2+0.) for i in range(0,NN2+1)])
	return np.multiply(1./(NN2+0.),fn(ppos,beta,al,n))

def alx_manual(gammaL,gammaH,pposL,pposH):
	ret = []

	#Fixation
	fN = B*fixNeut_manual()
	fNeg = B*fixNegB_manual(0.5*pposH+0.5*pposL)
	fPosL = fixPosSim_manual(gammaL,0.5*pposL)
	fPosH = fixPosSim_manual(gammaH,0.5*pposH)

	#Pol
	neut = cumuSfs_manual(DiscSFSNeutDown_manual())
	selH = cumuSfs_manual(DiscSFSSelPosDown_manual(gammaH,pposH))
	selL = cumuSfs_manual(DiscSFSSelPosDown_manual(gammaL,pposL))
	selN = cumuSfs_manual(DiscSFSSelNegDown_manual(pposH+pposL))
	
	sel = []
	for i in range(0,len(selH)):
		sel.append((selH[i]+selL[i])+selN[i])
	for i in range(0,nn-1):
		ret.append(float(1. - (fN/(fPosL + fPosH+  fNeg+0.))* sel[i]/neut[i]))
	return (ret,sel,neut,selH,selL,selN,fN,fNeg,fPosL,fPosH)

############################
# fN = B*fixNeut()
# fNeg = B*fixNegB(0.5*pposH+0.5*pposL)
# fPosL = fixPosSim(gL,0.5*pposL)
# fPosH = fixPosSim(gH,0.5*pposH)

# #Pol
# neut = cumuSfs(DiscSFSNeutDown())
# selH = cumuSfs(DiscSFSSelPosDown(gH,pposH))
# selL = cumuSfs(DiscSFSSelPosDown(gL,pposL))
# selN = cumuSfs(DiscSFSSelNegDown(pposH+pposL))


#%timeit fix = pFix(gamma=gH,NN=500)
#%timeit fix = pFix_numba(gamma=gH,NN=500)


# fPosH= fixPosSim_numba(gamma=gH,ppos=pposH,gam_neg=gam_neg,NN=NN,rho=rho,Lf=Lf)
#%timeit fPosH = fixPosSim(gamma=gH,ppos=pposH,gam_neg=gam_neg,NN=NN,rho=rho,Lf=Lf)
#%timeit fPosH = fixPosSim_numba(gamma=gH,ppos=pposH,gam_neg=gam_neg,NN=NN,rho=rho,Lf=Lf)

##############Neutral fixation##############


# fN    = B*fixNeut_numba(B,NN=1000)
#%time fN    = B*fixNeut(B,NN=1000)
#%time fN    = B*fixNeut_numba(B,NN=1000)


##############Negative fixation##############



# fNeg          = B*fixNegB_numba(0.5*pposH+0.5*pposL,al,B,N,be)
#%timeit fNeg  = B*fixNegB(0.5*pposH+0.5*pposL,al,B,N,be)
#%timeit fNeg  = B*fixNegB_numba(0.5*pposH+0.5*pposL,al,B,N,be)


############################SFS############################


# for i in range(0,10):
	#%time alpha = alx(gL,gH,pposL,pposH,NN,B,gam_neg,rho,Lf,theta_mid_neutral)
# print('++++++++++++++++')
# for i in range(0,10):
	#%time alpha = alx_numba(gL,gH,pposL,pposH,NN,B,gam_neg,rho,Lf,theta_mid_neutral)

#%lprun -f alx alpha = alx(gL,gH,pposL,pposH,NN,B,gam_neg,rho,Lf,theta_mid_neutral)
#%lprun -f alx_numba alpha = alx_numba(gL,gH,pposL,pposH,NN,B,gam_neg,rho,Lf,theta_mid_neutral)

#########################NUMBA#########################

@njit
def pFix_numba(gamma):
	s = gamma/(NN+0.)
	pfix = (1.-np.exp(-2.*s))/(1.-np.exp(-2.*gamma))
	if s >= 0.1:
		pfix = np.exp(-(1.+s))
		lim = 0
		# for i in prange(1,200):
		# 	pfix = np.exp((1.+s)*(pfix-1.))
		# 	lim +=1
		while(lim < 200):
			pfix = np.exp((1.+s)*(pfix-1.))
			lim +=1
		pfix = 1-pfix
	return pfix	

@overload(polygamma)
def fixPosSim_numba(gamma,ppos):
	
	S = abs(gam_neg/(1.*NN))
	r = rho/(2.*NN)
	u = theta_f/(2.*NN)
	s = gamma/(NN*1.)
	p0 = polygamma(1,(s+S)/r)
	p1 = polygamma(1,(r+Lf*r+s+S)/r)
	CC = 1.
	#print mpmath.exp(-2.*S*u*(p0-p1)/r**2)
	return 0.745*ppos*np.exp(-2.*S*u*(p0-p1)*CC**2/r**2)*pFix_numba(gamma)

@njit(fastmath=True)
def fixNeut_numba():
		return 0.255*(1./(float(B)*float(NN)))
@overload(mpmath.zeta)
def fixNegB_numba(ppos):
	return 0.745*(1-ppos)*(2**(-al))*(B**(-al))*(be**al)*(-mpmath.zeta(al,1.+be/(2.*B))+mpmath.zeta(al,0.5*(2-1./(N*B)+be/B)))

@overload(mpmath.zeta)
def FullNeg_numba( ppos, x):
	beta = be/(1.*B)	
	if x > 0 and x < 1.:
		return (1.-ppos)*(2.**-al)*(beta**al)*(-mpmath.zeta(al,np.array(x)+beta/2.) + mpmath.zeta(al,(2+beta)/2.))/((-1.+x)*x)
	else:
		return 0.

@njit
def returnAlpha():
	ret = np.zeros((nn-1,))
	for i in range(0,nn-1):
		ret[i] = (float(fN/(fNeg+0.)))
	return ret

def alx_numba(gammaL,gammaH,pposL,pposH):
	
	fN    = B*fixNeut_numba(B,NN=NN)
	fNeg  = B*fixNegB_numba(0.5*pposH+0.5*pposL,al,B,N,be)
	fPosL = fixPosSim_numba(gamma=gH,ppos=pposH,gam_neg=gam_neg,NN=NN,rho=rho,Lf=Lf)
	fPosH = fixPosSim_numba(gamma=gH,ppos=pposH,gam_neg=gam_neg,NN=NN,rho=rho,Lf=Lf)	
	ret = returnAlpha(fN,float(fNeg),50)
	ret2 = returnAlpha(fN,float(fNeg),50)

	return ret
	# neut  = cumuSfs(DiscSFSNeutDown())
	# selH  = cumuSfs(DiscSFSSelPosDown(gammaH,pposH))
	# selL  = cumuSfs(DiscSFSSelPosDown(gammaL,pposL))

	# #%time selN  = cumuSfs_numba(DiscSFSSelNegDown_numba(pposH+pposL))
	# sel   = np.zeros_like(selH)
	# ret   = []
	# for i in range(0,len(selH)):
	# 	sel.append((selH[i]+selL[i])+selN[i])
	# for i in range(0,nn-1):
	# 	ret.append(float(1. - (fN/(fPosL + fPosH+  fNeg+0.))* sel[i]/neut[i]))


@njit
def binomOp_numba():

	NN = int(round(NN*B))
	samps2 = np.zeros((nn+1,NN+1))
	sampFreqs2 = np.zeros((nn+1,NN+1))

	for i in range(0,nn+1):
		for j in range(0,NN+1):
			samps2[i,j] = i
			sampFreqs2[i,j] = j/(NN+0.)

	return samps2,nn,sampFreqs2


@njit
def cumuSfs_numba(sfsTemp):
	
	out = np.empty((len(sfsTemp)+1))
	out[0] = np.sum(sfsTemp)

	for i in range(0,len(sfsTemp)):
		app = out[i]-sfsTemp[i]
		if app > 0.:
			out[i+1] = app
		else:
			out[i+1] = (0.) 
	return out

def DiscSFSSelNeg_numba(ppos):
	NN2 = int(round(NN*B))
	dFunc = np.vectorize(FullNeg_manual)
	return np.multiply(1./(NN2+0.),dFunc(ppos,[i/(NN2+0.) for i in range(0,NN2+1)]))


def DiscSFSSelNegDown_numba(ppos):
	a,b,c = binomOp_numba()
	bn = binom.pmf(a,b,c)
	return B*(theta_mid_neutral)*0.745*(np.dot(bn,DiscSFSSelNeg(ppos)))[1:-1]



# import time 
# s = time.time()
# a=cumuSfs()
# print(time.time() -s )

# import time 
# s = time.time()
# a= cumuSfs_numba()
# print(time.time() -s )

#%time a=binomOp(B,500,50)
#%time a=binomOp_numba(B,500,50)
#%time a=binomOp(B,5000,500)
#%time a=binomOp_numba(B,5000,500)

#%time a = binomOp(B,50000,5000)
#%time a = 




# import numpy as N
# def z(r,i,E=1e-24):
# 	R=0;I=0;n=0;
# 	while(True):
# 		a=0;b=0;m=2(-n-1)
# 		for k in range(0,n+1):
# 			M=(-1)kN.product([x/(x-(n-k))for x in range(n-k+1,n+1)]);A=(k+1)**-r;t=-iN.log(k+1);a+=MAN.cos(t);b+=MAN.sin(t)
# 			a=m;b=m;R+=a;I+=b;n+=1
# 			if aa+bb<E:
# 				break
# 			A=2*(1-r);t=-iN.log(2);a=1-AN.cos(t);b=-AN.sin(t);d=aa+bb;a=a/d;b=-b/d
# 	print(Ra-Ib,Rb+Ia)

# 	return()