
import subprocess
import sys
import os
import numpy as np
import pandas as pd


class Mk:
	
	def __init__(self,daf,div,table=False):
	
		self.daf = daf
		self.div = div
		self.table = table

	def createTable(self):
		## Create MKT table 
		output = np.array([np.sum(self.daf['P0']),np.sum(self.daf['Pi']),self.div['D0'],self.div['Di']])
		return(output)

	# def standard(self):
	# 	## Declare output data frame
	# 	# output = pd.DataFrame({'alpha': float, 'pvalue' : float},index=[1])
	# 	output = np.zeros((1),dtype=[('alpha', '<i4'), ('pvalue', '<i4')])
		

	# 	## Estimation of alpha and fisher exact test p-value
		# alpha <- 1-(mktTable[2,1]/mktTable[1,1])*(mktTable[1,2]/mktTable[2,2])
		# pvalue <- fisher.test(mktTable)$p.value
		
	# 	## Ka, Ks, omega, omegaA, omegaD
	# 	Ka     = divergence$Di/divergence$mi
	# 	Ks     = divergence$D0/divergence$m0
	# 	omega  = Ka/Ks
	# 	omegaA = omega*alpha
	# 	omegaD = omega-omegaA
	# 	divergenceMetrics = data.frame(Ka, Ks, omega, omegaA, omegaD)
		
	# 	## Output  

	# 	output <- list('alpha'= data.frame('alpha' = alpha,'Fisher Test'= pvalue),
	# 		'mktTable'= mktTable,'divMetrics'=divergenceMetrics)
		
	# 	return(output)
	# def emkt:

	# def amkt:

