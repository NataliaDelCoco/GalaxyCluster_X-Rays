import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import pymc3 as pm
from scipy import optimize
import math
#====================================================================
#    FUNCOES
#====================================================================

def kT_func(R_med, Rc, Rs, eta,beta,alpha,eps,T0, T1):
  num1=(R_med/Rc)**(-eta)
  div=(1+((R_med/Rc)**2))**(beta - eta/2)
  num2=(1+(R_med/Rs)**alpha)**(eps/alpha)
  
  aux=(num1/div)*num2

  T=T0 + T1*aux
  return T
  #------------------------------------------------------------------

def chi2red(obs,exp,par_model, obs_err):
  chi_vet=[]
  soma=0
  for i in range(len(obs)):
    aux=((((obs[i]-exp[i])**2))/(obs_err[i]**2))
    soma = soma+aux
    chi_vet.append(aux)

  chi=np.sum(soma)
  num=(len(obs)) - par_model
  chired=chi/num
  return [chired, chi_vet]


#====================================================================
#    ENTRADAS
#====================================================================
# in_data=sys.argv[1]

in_data_T=sys.argv[1]
sem_ultimo=int(sys.argv[2]) #sem_ultimo = 1 > exclui ultimo ponto
# in_data_Ab='metal-profile-err.txt'

# nom,ext=in_data_.split('.')
# nom=nom.split('_')
# nome=nom[0]

data_T=pd.read_csv(in_data_T, delimiter=" ", header=None)

if (sem_ultimo == 1):
  data_T = data_T[:-1]
# data_Ab=pd.read_csv(in_data_Ab, delimiter=" ", header=None)

R_arc=data_T[0]
kT=data_T[1]
kTinf=data_T[2]
kTerInf=kT-kTinf
kTsup=data_T[3]
kTerSup=kTsup-kT

# R_arc=data_T[0]
# Ab=data_Ab[1]
# Abinf=data_Ab[2]
# AberInf=Ab-Abinf
# Absup=data_Ab[3]
# AberSup=Absup-Ab



#erro kT medio
kTEr_med=[]
for i in range(len(kTerInf)) :
  aux=(kTerInf[i]+kTerSup[i])/2
  if math.isnan(aux) == True:
    aux2=kT[i]*0.25
    kTEr_med.append(aux2)
  else:
    kTEr_med.append(aux)

#=====================================================================
#    AJUSTE FUNCAO
#=====================================================================

var=0.000001
eta= 0 
beta = 1
alpha = 1
eps= 0 
T1=3
T0 = 0
Rc= 1.
Rs=3.

#solta Rc e beta
params, params_cov= optimize.curve_fit(kT_func, R_arc, kT, sigma=kTEr_med ,p0=[Rc,Rs,eta, beta, alpha, eps, T0,T1],\
  bounds=((0.,Rs-var,eta-var,0,alpha-var,eps-var,T0-var,T1-var),\
  (2.9,Rs+var,eta+var,np.inf,alpha+var,eps+var,T0+var,T1+var)))


Rc_0=params[0]
beta_0=params[3]

#solta eta
params, params_cov= optimize.curve_fit(kT_func, R_arc, kT, sigma=kTEr_med ,p0=[Rc_0,Rs, eta,beta_0,alpha, eps,T0,T1],\
  bounds=((0.,Rs-var,0,0.,alpha-var,eps-var,T0-var,T1-var),\
  (3,Rs+var,np.inf,np.inf,alpha+var,eps+var, T0+var,T1+var)))


Rc_0=params[0]
eta_0=params[2]
beta_0=params[3]

#solta Rs, alpha, eps
params, params_cov= optimize.curve_fit(kT_func, R_arc, kT, sigma=kTEr_med ,p0=[Rc_0,Rs,eta_0,beta_0,alpha, eps,T0,T1],\
  bounds=((0.,0,0,0, 0., 0.,T0-var,T1-var),\
  (3,10,np.inf,np.inf,np.inf, np.inf,T0+var,T1+var)))

Rc_0=params[0]
Rs_0=params[1]
eta_0=params[2]
beta_0=params[3]
alpha_0=params[4]
eps_0=params[5]

#solta T0 e T1
params, params_cov= optimize.curve_fit(kT_func, R_arc, kT, sigma=kTEr_med ,p0=[Rc_0,Rs_0,eta_0,beta_0,alpha_0,eps_0,T0,T1],\
  bounds=((0., 0,0,0,0,0,0,0),\
  (3,10,np.inf, np.inf,np.inf,np.inf, np.inf,np.inf)))

Rc_0=params[0]
Rs_0=params[1]
eta_0=params[2]
beta_0=params[3]
alpha_0=params[4]
eps_0=params[5]
T0_0=params[6]
T1_0=params[7]

Rc_0_sd=np.sqrt(params_cov[0][0])
Rs_0_sd=np.sqrt(params_cov[1][1])
eta_0_sd=np.sqrt(params_cov[2][2])
beta_0_sd=np.sqrt(params_cov[3][3])
alpha_0_sd=np.sqrt(params_cov[4][4])
eps_0_sd=np.sqrt(params_cov[5][5])
T0_0_sd=np.sqrt(params_cov[6][6])
T1_0_sd=np.sqrt(params_cov[7][7])



kT_ajuste=kT_func(R_arc, Rc_0, Rs_0, eta_0, beta_0,alpha_0,eps_0, T0_0, T1_0)


chired,chivet=chi2red(kT,kT_ajuste,8,kTEr_med)
resid=kT_ajuste - kT

# SALVA AS SAIDAS EM FILES
nome_sai='kT_fit_params.txt'
if (sem_ultimo ==1):
  nome_sai='kT_fit_params_sem_ultimo.txt'

head1=('Parametros beta model modificado - Temperature')
head2=('rc e rs em arcmin')
head3=('T1, T1Err, rc, rcErr, rs, rsErr, eta, etaErr,beta, betaErr, epsilon, epsilonErr, alpha,\
 alphaErr, T0, T0Err')
val=[T1_0,T1_0_sd,Rc_0,Rc_0_sd,Rs_0,Rs_0_sd,eta_0,eta_0_sd,beta_0,beta_0_sd,\
alpha_0,alpha_0_sd,eps_0,eps_0_sd,T0_0,T0_0_sd]

with open(nome_sai,'w') as f:
  f.write(head1)
  f.write("\n")
  f.write(head2)
  f.write("\n")
  f.write(head3)
  f.write("\n")
  for i in val:
    u=str(i)
    f.write(u)
    f.write(", ")
f.close()

# columns='R (arcmin), R (kpc), kT (keV), kTErr, kT_ajustado (keV), chi2, Residuo'
# cols=[R_med, R_kpc, kT, kTEr_med, kT_ajuste, chivet, resid]

# nome_sai2=nome+'_BrilhoSup_fit.txt'
# np.savetxt(nome_sai2, np.transpose(cols), header=columns, delimiter=',', fmt='%s')


nome_sai='temp_prof_fit.png'
if (sem_ultimo ==1):
  nome_sai='temp_prof_fit_sem_ultimo.png'


plt.errorbar(R_arc,kT, yerr=kTEr_med, capsize=5, linestyle='',mec='dimgrey', mfc='dimgrey', \
ms=6, elinewidth=1, ecolor='dimgrey')
plt.plot(R_arc, kT_ajuste, color='indianred')
plt.xlabel('RA (arcmin)')
plt.ylabel('kT (keV)')


plt.savefig(nome_sai)
























