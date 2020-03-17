#!/usr/bin/python
#$ -e sim.div.log
#$ -o sim.div.log
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=48:00:00
#$ -t 1-100

from adapter import mkcalc
from scipy.optimize import fsolve
import sys

B  = float(sys.argv[1])
al = float(sys.argv[2])
alLow = float(sys.argv[3])

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

adap1 = mkcalc.AsympMK(B=0.4,gam_neg=-83,theta_f=0.001,alLow=0,alTot=0.2,neut_mid=False,L_mid=501,Lf=10**6,N=500,nsim=2500,pref="unc",gL=10,n=25,gH=500)
df1 = makeDf(adap1.alx(adap.gL,adap.gH,adap.pposL,adap.pposH),adap1.alx_nopos(adap.gL,adap.gH,adap.pposL,adap.pposH))

adap2 = mkcalc.AsympMK(B=0.6,gam_neg=-83,theta_f=0.001,alLow=0.2,alTot=0.2,neut_mid=False,L_mid=501,Lf=10**6,N=500,nsim=2500,pref="unc",gL=10,n=25,gH=500)

df2 = makeDf(adap2.alx(adap.gL,adap.gH,adap.pposL,adap.pposH),adap2.alx_nopos(adap.gL,adap.gH,adap.pposL,adap.pposH))

adap3 = mkcalc.AsympMK(B=0.8,gam_neg=-83,theta_f=0.001,alLow=0.2,alTot=0.2,neut_mid=False,L_mid=501,Lf=10**6,N=500,nsim=2500,pref="unc",gL=10,n=25,gH=500)
df3 = makeDf(adap3.alx(adap.gL,adap.gH,adap.pposL,adap.pposH),adap3.alx_nopos(adap.gL,adap.gH,adap.pposL,adap.pposH))

adap4 = mkcalc.AsympMK(B=0.999,gam_neg=-83,theta_f=0.001,alLow=0,alTot=0.2,neut_mid=False,L_mid=501,Lf=10**6,N=500,nsim=2500,pref="unc",gL=10,n=25,gH=500)
df4 = makeDf(adap4.alx(adap.gL,adap.gH,adap.pposL,adap.pposH),adap4.alx_nopos(adap.gL,adap.gH,adap.pposL,adap.pposH))

plt.plot(df1['counts'].values[1:],df1['pos'].values[1:], linestyle='-',color='blue')
plt.plot(df2['counts'].values[1:], df2['pos'].values[1:], linestyle='-', color='orange');
plt.plot(df3['counts'].values[1:], df3['pos'].values[1:], linestyle='-', color='yellow');
plt.plot(df4['counts'].values[1:], df4['pos'].values[1:], linestyle='-', color='red');
plt.plot(df4['counts'].values[1:], df4['pos'].values[1:], linestyle='-', color='red');
plt.plot(df4['counts'].values[1:], [0.2]*48, linestyle='--', color='black');
plt.ylim([-0.4,0.2])
plt.show(block=False)


fig = plt.figure(figsize=(8, 6))
plt.style.use('seaborn')
plt.plot(df['counts'].values,df['pos'].values, linestyle='-', marker='o',color='#efb2ba')
plt.plot(df['counts'].values, df['nopos'].values, linestyle='-', marker='o', color='#87d1e6');
plt.show(block=False)

# run the calucation and print out
# the out format looks like

#1 -0.221219196467 c a
#1 -0.20387225102 c np
#2 -0.0986506815748 c a
#2 -0.079336980885 c np

# the first column is allele count x
# the second column is alpha(x)
# the third column you can ignore 
# the fourth column indicates whether the alpha(x) includes segregating beneficial alleles

adap.print_alx(adap.gL,adap.gH,adap.pposL,adap.pposH)
pos = adap.alx(adap.gL,adap.gH,adap.pposL,adap.pposH)
nopos=adap.alx_nopos(adap.gL,adap.gH,adap.pposL,adap.pposH)
data = {'pos':pos,'nopos':nopos}
df = pd.DataFrame(data)
df = df.reset_index().rename(columns={'index':'counts'})

fig = plt.figure(figsize=(8, 6))
plt.style.use('seaborn')
plt.plot(df['counts'].values,df['pos'].values, linestyle='-', marker='o',color='#efb2ba')
plt.plot(df['counts'].values, df['nopos'].values, linestyle='-', marker='o', color='#87d1e6');
plt.show(block=False)

###########################
adap = mkcalc.AsympMK(B=0.999,gam_neg=-83,theta_f=0.001,alLow=0.1,alTot=0.074,neut_mid=False,L_mid=501,Lf=10**6,N=500,nsim=2500,pref="unc",gL=10,n=25,gH=500)

# here the software calculates the mutation rates corresponding to the desired terms
adap.set_theta_f()
adap.setPpos()
adap.theta_f = theta_f

pos = adap.alx(adap.gL,adap.gH,adap.pposL,adap.pposH)
nopos = adap.alx_nopos(adap.gL,adap.gH,adap.pposL,adap.pposH)

def makeDf(alx,alxnopos):
	data = {'pos':alx,'nopos':alxnopos}
	df = pd.DataFrame(data)
	df = df.reset_index().rename(columns={'index':'counts'})
	return(df)

fig = plt.figure(figsize=(8, 6))
plt.style.use('seaborn')
plt.plot(df['counts'].values,df['pos'].values, linestyle='-', marker='o',color='#87d1e6')
plt.plot(df['counts'].values, df['nopos'].values, linestyle='-', marker='o', color='#efb2ba');
plt.plot(df['counts'].values, np.array([0.074]*df.shape[0]), color='gray');

plt.show(block=False)