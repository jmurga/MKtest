#!/usr/bin/python
#$ -e sim.div.log
#$ -o sim.div.log
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=48:00:00
#$ -t 1-100

# initialize parameters
# B is desired fraction of diversity remaining (on scale from 0 to 1, although if you give it exactly 0 or 1 it will choke)
# al is the desired total alpha
# alLow is the proportion of alpha due to weakly beneficial fixations
# N is population size, but the results only very weakly depend on N
# n is sample size
# gL is 2Ns for weakly beneficial mutations
# gam_neg is 2Ns for deleterious alleles driving Background selection

# several other parameters make no difference in the calculations (only matter when calling the software for simulations), including:
	# L_mid is the length of the central coding locus (is effectively doubled in the simulations though since two such loci are simulated per simulation)
	# Lf is the length of flanking sequence harboring deleterious mutations on each side of the coding loci 

from adapter import mkcalc
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import sys



B  = float(0.999)
al = float(0.2)
alLow = [0,0.1,0.2]
weak=0.2
###########################
f = []
# for weak in alLow:
import time
s = time.time()

adap = mkcalc.AsympMK(B=B,gam_neg=-83,theta_f=0.001,alLow=weak,alTot=al,neut_mid=False,L_mid=501,Lf=10**6,N=500,nsim=2500,pref="unc",gL=10,n=250,gH=500)

# here the software calculates the mutation rates corresponding to the desired terms
adap.set_theta_f()
theta_f = adap.theta_f
adap.B = 0.999
adap.set_theta_f()
adap.setPpos()
adap.theta_f = theta_f
adap.B = B


	# run the calucation
pos = adap.alx(adap.gL,adap.gH,adap.pposL,adap.pposH)
nopos = adap.alx_nopos(adap.gL,adap.gH,adap.pposL,adap.pposH)

def asympPlot(dfAlpha,originalAsymp):
	fig    = plt.figure(figsize=(8, 6))
	plt.style.use('seaborn')
	
	x      = dfAlpha.counts.values
	y1     = dfAlpha.alphaPos.values
	y2     = dfAlpha.alphaNonPos.values

	plt.plot(x, y1,'-',color='#91d5db')
	plt.plot(x, y1,'o',color='#91d5db')
	plt.plot(x, y2,'-',color='#fab4be')
	plt.plot(x, y2,'o',color='#fab4be')
	plt.plot(x, [originalAsymp]*x.shape[0],'-.',color='gray')
	plt.grid()
	plt.xscale('log')
	plt.xlim(0,60)
	# plt.ylim((y1.min() - (y1.min()*50/100)),0.25)
	plt.xticks([1,2,5,10,20,50],[1,2,5,10,20,50])
	plt.grid()

	return(fig)

def bgsPlot(freq,bgs1,bgs2,bgs3,bgs4,trueAlpha):
	fig    = plt.figure(figsize=(8, 6))
	plt.style.use('seaborn')
	
	x      = freq
	y1     = bgs1
	y2     = bgs2
	y3     = bgs3
	y4     = bgs4

	plt.plot(x, y1,'-',color='blue')
	plt.plot(x, y2,'-',color='orange')
	plt.plot(x, y3,'-',color='yellow')
	plt.plot(x, y4,'-',color='red')
	plt.plot(x, [trueAlpha]*x.shape[0],'-.',color='gray')
	plt.grid()
	plt.xscale('log')
	plt.xticks([1,2,5,10,20,50],[1,2,5,10,20,50])
	plt.grid()

	return(fig)











	# df.to_csv('/home/jmurga/alphaWeakly'+str(weak)+'.tsv',sep='\t',index=False,header=True)
###########################
B  = [0.4,0.6,0.8,0.999]
al = float(0.2)
alLow = float(0)
###########################
dfList = list()
for bgs in B:
	print(bgs)
	adap = mkcalc.AsympMK(B=bgs,gam_neg=-83,theta_f=0.001,alLow=alLow,alTot=al,neut_mid=False,L_mid=501,Lf=10**6,N=500,nsim=2500,pref="unc",gL=10,n=25,gH=500)

	# here the software calculates the mutation rates corresponding to the desired terms
	adap.set_theta_f()
	theta_f = adap.theta_f
	adap.B = 0.999
	adap.set_theta_f()
	adap.setPpos()
	adap.theta_f = theta_f
	adap.B = bgs
	# run the calucation

	pos = adap.alx(adap.gL,adap.gH,adap.pposL,adap.pposH)
	nopos = adap.alx_nopos(adap.gL,adap.gH,adap.pposL,adap.pposH)
	dfList.append(adap.makeDf(pos,nopos))


bgsPlot(dfList[0].counts.values,dfList[0].alphaPos.values,dfList[1].alphaPos.values,dfList[2].alphaPos.values,dfList[3].alphaPos.values,0.2).show()
bgsPlot(dfList[0].counts.values,dfList[0].alphaNonPos.values,dfList[1].alphaNonPos.values,dfList[2].alphaNonPos.values,dfList[3].alphaNonPos.values,0.2).show()
# library(ggplot2)
# library(data.table)
# df <- fread('/home/jmurga/bgs.tsv')
# df <- as.data.frame(df)
# df$bgs <- as.factor(df$bgs)
# ggplot(df) + geom_line(aes(x=counts,y=nopos,colour=bgs)) + geom_point(aes(x=counts,y=nopos,colour=bgs)) +  scale_x_continuous(trans="log10",breaks=c(1,2,5,10,20,50,seq(1, 50, by = 39)))