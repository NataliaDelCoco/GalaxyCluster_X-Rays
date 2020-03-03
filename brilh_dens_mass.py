import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import pymc3 as pm
import scipy
from scipy import optimize
from scipy import stats
import matplotlib.ticker
from scipy.integrate import quad
import math
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import astropy.units as u
import csv
import astropy.constants as ctes
#===================================================================
#    FUNCOES
#===================================================================
def brilho(R,sig0,rc,rs,beta,gamma, eps, alpha, cte):
  num_b1 = (R/rc)**(-gamma)
  pot_b1 = 3*beta - (1/2) - (gamma/2)
  div_b1 = (1+(R/rc)**2)**pot_b1
  b1=num_b1/div_b1

  pot_b2 = eps/alpha
  b2 = (1+(R/rs)**alpha)**pot_b2

  return sig0*b1*b2 + cte

#------------------------------------------------------------------

def integrand(R,rc,rs,beta,gamma, eps, alpha, r):
  p1_1 = 2*(1/2 - 3*beta + gamma/2)*R*(1 + (R/rc)**2)**(-(1/2) - 3*beta + gamma/2 )
  p1_2 = ((R/rc)**(-gamma))*(1 + (R/rs)**alpha)**(eps/alpha)
  p1 = p1_1*p1_2/(rc**2)

  p2_1 = gamma*(1 + (R/rc)**2)**(1/2 - 3*beta + gamma/2)*((R/rc)**(-1 - gamma))
  p2_2 = (1 + (R/rs)**alpha)**(eps/alpha)
  p2 = p2_1*p2_2/rc

  p3_1 = eps*(1 + (R/rc)**2)**(1/2 - 3*beta + gamma/2)*((R/rc)**(-gamma))
  p3_2 = ((1 + (R/rs)**alpha)**(-1 + eps/alpha))*((R/rs)**(-1 + alpha))
  p3 = p3_1*p3_2/rs

  der = (p1 - p2 + p3)

  aux=np.sqrt(R**2 - r**2)

  fim=der/aux

  return fim

#------------------------------------------------------------------

def densi(r,rc,rs,beta,gamma,eps,alpha):
  num_d1 = (r/rc)**(-gamma/2)
  pot_d1 = (3*beta - (gamma/2))/2
  div_d1 = (1+(r/rc)**2)**pot_d1
  d1 = num_d1/div_d1

  pot_d2 = (-eps/alpha)/2
  d2 = (1+(r/rs)**alpha)**pot_d2
  return d1*d2

#------------------------------------------------------------------

def densi_n0(r,rc,rs,beta,gamma,eps,alpha,n0):
  num_d1 = (r/rc)**(-gamma/2)
  pot_d1 = (3*beta - (gamma/2))/2
  div_d1 = (1+(r/rc)**2)**pot_d1
  d1 = num_d1/div_d1

  pot_d2 = (-eps/alpha)/2
  d2 = (1+(r/rs)**alpha)**pot_d2
  return n0*d1*d2


#------------------------------------------------------------------

def integ_densi_n0(r,rc,rs,beta,gamma,eps,alpha,n0):
  num_d1 = (r/rc)**(-gamma/2)
  pot_d1 = (3*beta - (gamma/2))/2
  div_d1 = (1+(r/rc)**2)**pot_d1
  d1 = num_d1/div_d1

  pot_d2 = (-eps/alpha)/2
  d2 = (1+(r/rs)**alpha)**pot_d2

  densi = n0*d1*d2*(r**2)

  return densi
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

#------------------------------------------------------------------

def integ_norm(r,rc,rs,beta,gamma,eps,alpha):
  num_d1 = (r/rc)**(-gamma/2)
  pot_d1 = (3*beta - (gamma/2))/2
  div_d1 = (1+(r/rc)**2)**pot_d1
  d1 = num_d1/div_d1

  pot_d2 = (-eps/alpha)/2
  d2 = (1+(r/rs)**alpha)**pot_d2

  ngas=d1*d2
  ngas2=ngas**2

  integrand=ngas2*4*np.pi*(r**2)

  return integrand

#------------------------------------------------------------------
def der_ln_pho(r, rc, rs, beta, gamma, alpha, eps):
  p1=-gamma/2
  p2=-eps/(1+(rs/r)**alpha)
  p3=-(3*beta - gamma/2)/(1+(rc/r)**2)

  return p1+p2+p3

#------------------------------------------------------------------
#MASSA TOTAL MAUGHAN 2006
def M500_maug(Yx,z):
  OmM=0.3
  OmL=0.7
  h=0.7
  By=0.57
  Ez=np.sqrt(OmM*(1+z)**3 + (1-OmM-OmL)*(1+z)**2 + OmL)
  Ay=5.77*1e14*np.sqrt(h) #np.sqrt(h) #Solar masses
  

  aux_M = Ez**(-2/5) * Ay * (Yx**By) / (3*1e14)
  return aux_M

#===================================================================
#    ENTRADAS
#===================================================================
in_norm='norm_vals.csv'
in_brilho='comb-obj-im-400-7200-ps_brilhoSup.csv'
in_T='avrgT.txt'
out_name='A773'
z=0.217

#===================================================================
#LE ENTRADAS
#===================================================================
#Dados da normalizacao
dataN=pd.read_csv(in_norm,delimiter=',')

ct1=dataN['ct1'].values
ct1_er=dataN['ct1_er'].values
ct2=dataN['ct2'].values
norm=dataN['norm'].values
norm_er=dataN['norm_er'].values
r_extr_arc=dataN['r_extracao(arcmin)'].values

#-------------------------------------------------------
# DADOS DO BRILHO SUPERFICIAL

dataS=pd.read_csv(in_brilho, delimiter=",")

dataS = dataS.drop(dataS.index[0])

S_obs=dataS['INTENS'].values
S_obsEr=dataS['INT_ERR'].values
R_kpc=dataS['RAIO (kpc)'].values
R_arc=dataS['RAIO (arcmin)'].values


for i, item in enumerate(S_obs):
  if item == 'INDEF':
    data=data.drop(data.index[i])


for i, item in enumerate(S_obsEr):
  if item == 'INDEF':
    aux1=np.where(S_obsEr != 'INDEF')[0]
    a=min(aux1)
    perc=float(S_obsEr[a])/float(S_obs[a])
    aux2=S_obs[i]*perc
    S_obsEr[i] = aux2


S_obs=S_obs.astype(float)
S_obsEr=S_obsEr.astype(float)


R_med=[] #ARCMIN
for i in range(len(R_arc)):
  if i == 0:
    aux=R_arc[i]**(3/2)
    aux2= (aux/2.)**(2/3)
    R_med.append(aux2)
  else:
    aux=R_arc[i]**(3/2) + R_arc[i-1]**(3/2)
    aux2= (aux/2)**(2/3)
    R_med.append(aux2)

S_obs_red=[]
S_obsEr_red=[]
R_med_red=[]
for i in range(len(S_obs)):
  if i>3:
    S_obs_red.append(S_obs[i])
    S_obsEr_red.append(S_obsEr[i])
    R_med_red.append(R_med[i])

#---------------------------------------------------------
#Dados T media e R500
dataT=pd.read_csv(in_T, delimiter=' ')

R500_arc=dataT['r500(arcmin)'].values
kT=dataT['kT(keV)'].values

R500_arc = R500_arc[-1]
kT=kT[-1]


#===================================================================
#CONVERSAO PARA CM
#===================================================================
conv=cosmo.angular_diameter_distance(z) #equivale a ang de 1 rad na dist z, em Mpc

Da_cm=(conv.to(u.cm)).value
Rext_cm=r_extr_arc[0]*Da_cm*np.pi/(180*60)
R500_cm=R500_arc*Da_cm*np.pi/(180*60)

R_med_cm=[]
R_med_kpc=[]

for item in R_med:
  aux=item*Da_cm*np.pi/(180*60)
  R_med_cm.append(aux)
  aux= aux * u.cm
  aux_mpc = aux.to(u.Mpc).value
  R_med_kpc.append(aux_mpc*1e3)

R500_aux = R500_cm * u.cm
Rext_aux = Rext_cm * u.cm
# R500_ErrInf = R500_ErrInf * u.Mpc
# R500_ErrSup = R500_ErrSup * u.Mpc

R500_kpc=(R500_aux.to(u.Mpc).value)*1e3
Rext_mpc=Rext_aux.to(u.Mpc).value
# R500_ErrInf_cm=R500_ErrInf.to(u.cm).value
# R500_ErrSup_cm=R500_ErrSup.to(u.cm).value

#===================================================================
#    AJUSTA BRILHO SUPERFICIAL
#===================================================================

#----OPTMIZE------------------------------------

# começa com beta simples:
#    b2=1 => eps=0  alpha=1
#    num_b1=1 => gamma=0
sig0=1.
rc=0.5
rs=5
beta=0.6
gamma=0.00000001
eps=0.00000001
alpha=3.

cte=0.00000001
var=0.00000001

params, params_cov=optimize.curve_fit(brilho, R_med, S_obs, sigma=S_obsEr, p0=[sig0,rc,rs,beta,gamma,eps,alpha,cte],\
  bounds=((0.,0.,rs-var,0.3,gamma-var, eps-var, alpha-var, cte-var),\
  (np.inf,np.inf,rs+var,2, gamma+var, eps+var, alpha+var,cte+var)))

sig0=params[0]
rc=params[1]
beta=params[3]

# agora solto gamma

params, params_cov=optimize.curve_fit(brilho, R_med, S_obs, sigma=S_obsEr, p0=[sig0,rc,rs,beta,1.,eps,alpha,cte],\
  bounds=((0.,0.,rs-var,0.3, 0., eps-var, alpha-var, cte-var),\
  (np.inf,np.inf,rs+var,1, np.inf, eps+var, alpha+var,cte+var)))

sig0=params[0]
rc=params[1]
beta=params[3]
gamma=params[4]

#agora solta rs, eps e alpha

params, params_cov=optimize.curve_fit(brilho, R_med, S_obs, sigma=S_obsEr, p0=[sig0,rc,rs,beta,gamma,1.,alpha,cte],\
  bounds=((0.,0.,rc,0.3, 0., 0, alpha-var, cte-var),\
  (np.inf,np.inf,10.,1,np.inf,5.,alpha+var,cte+var)))


sig0=params[0]
rc=params[1]
rs=params[2]
beta=params[3]
gamma=params[4]
eps=params[5]
alpha=params[6]

# solta cte

params, params_cov=optimize.curve_fit(brilho, R_med, S_obs, sigma=S_obsEr, p0=[sig0,rc,rs,beta,gamma,eps,alpha,0.05],\
  bounds=((0.,0.,rc,0.3, 0., 0, alpha-var, 0.),\
  (np.inf,np.inf,10.,1,np.inf,5.,alpha+var,np.inf)))

sig0=params[0]
rc=params[1]
rs=params[2]
beta=params[3]
gamma=params[4]
eps=params[5]
alpha=params[6]
cte=params[7]

#refaz

params, params_cov=optimize.curve_fit(brilho, R_med, S_obs, sigma=S_obsEr, p0=[sig0,rc,rs,beta,gamma,eps,alpha,cte],\
  bounds=((0.,0.,rc,0.3, 0., 0, alpha-var, 0.),\
  (np.inf,np.inf,10.,1,np.inf,5.,alpha+var,np.inf)))

sig0_bri=params[0]
rc_bri=params[1]
rs_bri=params[2]
beta_bri=params[3]
gamma_bri=params[4]
eps_bri=params[5]
alpha_bri=params[6]
cte_bri=params[7]

sig0_bri_sd=np.sqrt(params_cov[0][0])
rc_bri_sd=np.sqrt(params_cov[1][1])
rs_bri_sd=np.sqrt(params_cov[2][2])
beta_bri_sd=np.sqrt(params_cov[3][3])
gamma_bri_sd=np.sqrt(params_cov[4][4])
eps_bri_sd=np.sqrt(params_cov[5][5])
alpha_bri_sd=np.sqrt(params_cov[6][6])
cte_bri_sd=np.sqrt(params_cov[7][7])


S_opt=brilho(R_med,sig0_bri,rc_bri,rs_bri,beta_bri,gamma_bri,eps_bri,alpha_bri,cte_bri)

chired_S,chivet_S=chi2red(S_obs,S_opt,8,S_obsEr)
resid_S=S_opt-S_obs

# SALVA AS SAIDAS EM FILES
nome_sai=out_name+'_BrilhoSup_fit_params.txt'
head1=('Parametros beta model modificado - BRILHO SUPERFICIAL')
head2=('rc e rs em arcmin')
head3=('sig0, sig0Err, rc, rcErr, rs, rsErr, beta, betaErr,gamma, gammaErr, epsilon, epsilonErr, alpha,\
 alphaErr, bkg, bkgErr')
val=[sig0_bri,sig0_bri_sd,rc_bri,rc_bri_sd,rs_bri,rs_bri_sd,beta_bri,beta_bri_sd,gamma_bri,gamma_bri_sd,\
eps_bri,eps_bri_sd,alpha_bri,alpha_bri_sd,cte_bri,cte_bri_sd]

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

columns='R (arcmin), R (kpc), S_observado (counts/arcmin²), S_obs_err, S_ajustado, chi2, Residuo'
cols=[R_med, R_kpc, S_obs, S_obsEr, S_opt, chivet_S, resid_S]

nome_sai2=out_name+'_BrilhoSup_fit.txt'
np.savetxt(nome_sai2, np.transpose(cols), header=columns, delimiter=',', fmt='%1.6e')


#===================================================================
#    PLOTA BRILHO SUPERFICIAL
#===================================================================


heights = [6, 1,1]
gs_kw = dict(height_ratios=heights)

f, (ax1,ax2,ax3)=plt.subplots(ncols=1, nrows=3, sharex=True,gridspec_kw=gs_kw )
# f.suptitle('Perfil radial de brilho superficial - A2142')
f.suptitle('Surface brightness radial profile - Abell 1795')
#data and fit


ax1.errorbar(R_med, S_obs, yerr=S_obsEr, fmt='.', capsize=5, mec='dimgrey', mfc='dimgrey', \
 ms=6, elinewidth=1, ecolor='dimgrey' )
ax1.plot(R_med, S_opt, label='adapted-β-model', color='indianred')
ax1.set_yscale('log')
ax1.set_ylabel('Surface brightness (counts/arcmin²)')
ax1.legend(loc='best')

ax2.plot(R_med,chivet_S, linestyle='',marker='.', color='indianred')
ax2.axhline(y=0, linestyle='--',marker='', color='dimgrey')
ax2.set_title('$\chi$²', pad=-10., fontsize=8)

ax3.plot(R_med,resid_S, linestyle='',marker='.', color='indianred')
ax3.axhline(y=0, linestyle='--',marker='', color='dimgrey')
ax3.set_xscale('log')
ax3.set_title('Residue', pad=-10., fontsize=8)
ax3.set_xlabel('Radius (arcmin)')

ax3.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax3.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

nome_fig=out_name+'_BrilhoSup_fit.png'
plt.savefig(nome_fig)

#===================================================================
#    CALCULA DENSIDADE
#===================================================================
#apos achar a derivada da funcao do brilho sup, temos que avalia-la pros
#parametros ajustados

#cria um array log de raio com r_max = R_med_max
r_max=np.log10(R_med[-1])
r = np.logspace(-1.5, r_max, 100)


#conv para cm
r_max_cm=np.log10(R500_cm)
r_min_cm=np.log10(R_med_cm[0])
r_cm = np.logspace(r_min_cm, r_max_cm, 100)
rc_bri_cm=rc_bri*Da_cm*np.pi/(180*60)
rs_bri_cm=rs_bri*Da_cm*np.pi/(180*60)
rc_bri_sd_cm=rc_bri_sd*Da_cm*np.pi/(180*60)
rs_bri_sd_cm=rs_bri_sd*Da_cm*np.pi/(180*60)


#nabs = [-(1/Pi)*(rasultado da integral)]^(1/2)
#TEM QUE SER CALCULADO EM ARCMIN
# "nabs" sai em unidades de (arcmin)^(-1/2)

n = []
for raio in r:
  aux = quad(integrand,raio, np.inf, args=(rc_bri,rs_bri,beta_bri,gamma_bri,eps_bri,alpha_bri, raio))
  n.append(aux[0])

nabs=[]

for i in range(len(n)):
  aux=((-1/np.pi)*n[i])
  aux2=np.sqrt(aux)
  nabs.append(aux)

#-----------------------------------------------------------
# Acha erro da densidade
# NAO DA PRA FAZER EM CM, TEM QUE SER EM ARCMIN
rerun=25
nerr=[]


for i in range(len(r)):
  raio = r[i]
  print('RAIO')
  print(raio)
  aux_nerr=[]
  count=0
  porc=nabs[i]*0.15
  n=2
  erros=0
  while count <= rerun:
    rc_re = stats.truncnorm.rvs(0, np.inf,rc_bri, rc_bri_sd/n)
    rs_re = stats.truncnorm.rvs(0, np.inf,rs_bri, rs_bri_sd/n)
    beta_re = stats.truncnorm.rvs(0, 10, beta_bri, beta_bri_sd/n)
    gamma_re = np.random.normal(gamma_bri, gamma_bri_sd/n)
    eps_re = np.random.normal(eps_bri, eps_bri_sd/n)
    alpha_re = np.random.normal(alpha_bri, alpha_bri_sd/n)

    integral,ig=quad(integrand,raio, np.inf, args=(rc_re,rs_re,beta_re,gamma_re,eps_re,alpha_re, raio))
    val=np.sqrt((-1/np.pi)*integral)
    # print('val')
    # print(val)

    if math.isnan(val) == True:
      erros=erros+1
    elif val > (nabs[i] + 1*nabs[i]):
      erros=erros+1
    elif val <(nabs[i] - 1*nabs[i]):
      erros=erros+1
    else:
      aux_nerr.append(val)
      count=count+1
      erros=0

    if erros >20:
      count=rerun+1
      if i > 0:
        std=nerr[-1]


  std=np.std(aux_nerr)
  if std < porc:
    nerr.append(std)
  else:
    nerr.append(porc)

#Antes de achar a normalizacao, vamos fazer um primeiro ajuste dos valores da densidade,
# SEM CONSIDERAR "n0"

rc=1.
rs=5.
beta=0.65
gamma=0.00000001
eps=0.00000001
alpha=3.

var=0.00000001

#solta rc e beta
params, params_cov=optimize.curve_fit(densi, r, nabs,sigma=nerr , p0=[rc,rs,beta,gamma,eps,alpha],\
  bounds=((0.,rs-var,0.,gamma-var, eps-var, alpha-var),\
  (np.inf,rs+var,1., gamma+var, eps+var, alpha+var)))

rc=params[0]
beta=params[2]

# agora solto gamma
params, params_cov=optimize.curve_fit(densi, r, nabs,sigma=nerr , p0=[rc,rs,beta,gamma_bri,eps,alpha],\
  bounds=((0.2,rs-var,0.,0., eps-var, alpha-var),\
  (np.inf,rs+var,2., np.inf, eps+var, alpha+var)))


rc=params[0]
beta=params[2]
gamma=params[3]

#agora solta rs e eps

params, params_cov=optimize.curve_fit(densi, r, nabs, sigma=nerr , p0=[rc,rs_bri,beta,gamma,eps_bri,alpha],\
  bounds=((0.,rc,0., 0., 0., 0.),\
  (rs,np.inf,2.,np.inf,np.inf,np.inf)))

rc=params[0]
rs=params[1]
beta=params[2]
gamma=params[3]
eps=params[4]
alpha=params[5]

#refaz
params, params_cov=optimize.curve_fit(densi, r, nabs, sigma=nerr , p0=[rc,rs,beta,gamma,eps,alpha],\
  bounds=((0.,rc,0.3, 0., 0., alpha-var),\
  (rs,10.,2.,np.inf,np.inf,alpha+var)))


rc_den=params[0] 
rs_den=params[1]
beta_den=params[2]
gamma_den=params[3]
eps_den=params[4]
alpha_den=params[5]


rc_den_sd=np.sqrt(params_cov[0][0])
rs_den_sd=np.sqrt(params_cov[1][1])
beta_den_sd=np.sqrt(params_cov[2][2])
gamma_den_sd=np.sqrt(params_cov[3][3])
eps_den_sd=np.sqrt(params_cov[4][4])
alpha_den_sd=np.sqrt(params_cov[5][5])

#-----------------------------------------------------------
#Com o valor dos parametros, vamos tentar encontrar n0, usando 'norm'

rc_den_cm=rc_den*Da_cm*np.pi/(180*60)
rc_den_sd_cm = rc_den_sd*Da_cm*np.pi/(180*60)
rs_den_cm=rs_den*Da_cm*np.pi/(180*60)
rs_den_sd_cm = rs_den_sd*Da_cm*np.pi/(180*60)
val_int_norm, ig =quad(integ_norm, 0, Rext_cm,args=(rc_den_cm,rs_den_cm,beta_den,gamma_den,eps_den,alpha_den))

#nenh = k*ngas^2 = k * no^2 *pho(r)^2
#NORM EH DADA EM CM


#vamos fazer 3 opcoes: so norm, norm*ct1, norm*ct1*ct2
k=1.15
num_n0=norm*(Da_cm**2)*((1+z)**2)
div_n0 = (10**(-14))*k*val_int_norm
n0=np.sqrt(num_n0/div_n0) # UNIDADE DE CM^(-3)

num_n01=norm*ct1*(Da_cm**2)*((1+z)**2)
div_n01 = (10**(-14))*k*val_int_norm
n01=np.sqrt(num_n01/div_n0) # UNIDADE DE CM^(-3)

num_n012=norm*ct1*ct2*(Da_cm**2)*((1+z)**2)
div_n0 = (10**(-14))*k*val_int_norm
n012=np.sqrt(num_n012/div_n0) # UNIDADE DE CM^(-3)

#-----------------------------------------------------------
#N0 FINAL, RESIDUOS E CHI2

n_obs=[]
n_obs_err=[]

for i in range(len(nabs)):
  aux1=nabs[i]*n012
  n_obs.append(aux1[0])
  aux2=nerr[i]*n012
  n_obs_err.append(aux2)

n_opt=densi_n0(r,rc_den,rs_den,beta_den,gamma_den,eps_den,alpha_den,n012)

chired_n,chivet_n=chi2red(n_obs,n_opt,7,n_obs_err)
resid_n=n_opt-n_obs




# # SALVA AS SAIDAS EM FILES
nome_sai=out_name+'_densi_fit_params.txt'
head1=('Parametros beta model modificado - DENSIDADE')
head2=('n0 (cm-3), rc (arcmin), rcErr, rs (arcmin), rsErr, rc (cm), rcErr, rs (cm), rsErr,beta, betaErr,gamma, gammaErr, epsilon, epsilonErr, alpha,\
 alphaErr')
val=[n0,rc_den,rc_den_sd,rs_den,rs_den_sd, rc_den_cm,rc_den_sd_cm,rs_den_cm,rs_den_sd_cm,beta_den,beta_den_sd,gamma_den,gamma_den_sd,\
eps_den,eps_den_sd,alpha_den,alpha_den_sd]

with open(nome_sai,'w') as f:
  f.write(head1)
  f.write("\n")
  f.write(head2)
  f.write("\n")

  for i in val:
    u=str(i)
    f.write(u)
    f.write(", ")

columns='R (arcmin), R (cm), n_obs (cm^-3), n_obs_err, n_ajustado, chi2, Residuo'
cols=[r, r_cm, n_obs, n_obs_err, n_opt, chivet_n, resid_n]

nome_sai2=out_name+'_BrilhoSup_fit.txt'
np.savetxt(nome_sai2, np.transpose(cols), header=columns, delimiter=',')



#===================================================================
#    PLOTA DENSIDADE
#===================================================================


heights = [6, 1,1]
gs_kw = dict(height_ratios=heights)

f, (ax1,ax2,ax3)=plt.subplots(ncols=1, nrows=3, sharex=True,gridspec_kw=gs_kw )
# f.suptitle('Perfil radial de brilho superficial - A2142')
f.suptitle('Eletronic Density radial profile - Abell 1795')
#data and fit


ax1.errorbar(r, n_obs, yerr=n_obs_err, fmt='.', capsize=5, mec='dimgrey', mfc='dimgrey', \
 ms=6, elinewidth=1, ecolor='dimgrey' )
ax1.plot(r, n_opt, label='adapted-β-model', color='indianred')
ax1.set_yscale('log')
ax1.set_ylabel('Eletronic density ($cm^{-3}$')
ax1.legend(loc='best')

ax2.plot(r,chivet_n, linestyle=' ',marker='.', color='indianred')
ax2.axhline(y=0, linestyle='--',marker='', color='dimgrey')
ax2.set_title('$\chi$²', pad=-10., fontsize=8)

ax3.plot(r,resid_n, linestyle='',marker='.', color='indianred')
ax3.axhline(y=0, linestyle='--',marker='', color='dimgrey')
ax3.set_xscale('log')
ax3.set_title('Residue', pad=-10., fontsize=8)
ax3.set_xlabel('Radius (arcmin)')

ax3.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax3.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

nome_fig=out_name+'_densi_fit.png'
plt.savefig(nome_fig)

#===================================================================
#    MASSA DO GAS
#===================================================================
int_densi, ig=quad(integ_densi_n0,0, R500_cm, args=(rc_den_cm,rs_den_cm,beta_den,gamma_den,eps_den,alpha_den, n012))

rc_cm=rc*arcmin2cm
rs_cm=rs*arcmin2cm
int_densi, ig=quad(integ_densi_n0,0, R500_cm, args=(rc_cm,rs_cm,b,a,e,g,n012))


import astropy.units as u
mA=1.89210622e-27 #kg
Mgas=float(4*np.pi*mA*int_densi) #KG
Mgas = (Mgas * u.kg).to(u.M_sun).value #SOLAR MASS





#perfil massa do gas
int_densi_prof=[]
for i in range(len(r_cm)-1):
  r1=r_cm[i]
  r2=r_cm[i+1]
  aux,ig=quad(integ_densi_n0,r1, r2, args=(rc_den_cm,rs_den_cm,beta_den,gamma_den,eps_den,alpha_den, n0))
  int_densi_prof.append(aux)

Mgas_prof=[]
for item in int_densi_prof:
  aux=4*np.pi*mA*item
  Mgas_prof.append(aux)

Mgas_prof = (Mgas_prof *u.kg).to(u.M_sun).value

#===================================================================
#    MASSA TOTAL
#===================================================================
#calculo da massa------------------------------------------------------------------------

def Mtot_isot(params,kT,r,arcmin2cm):
  rcc,rss,alpha,beta,epsilon,gamma,_,bkg = pars
  rcc=rcc*arcmin2cm
  rss=rss*arcmin2cm
  r=r*arcmin2cm

  #Para o calculo da massa total, preciso deixar tudo em SI
  # passando kT de keV -> eV -> J
  conv_eV=1e3*((1 *u.eV).si.value)

  kT_si=kT*conv_eV #kg m^2 s^-2
# kT_ErrSup_si=kT_ErrSup*conv_eV
# kT_ErrInf_si=kT_ErrInf*conv_eV

  r_m = r/100
  mu=0.6
#--------------------------------
  #tanto faz a unidade, corta tudo
  del_lnpho=der_ln_pho(r, rcc, rss, alpha, beta, epsilon, alpha)
  aux=(kT_si*r_m)/(ctes.G.value*mu*ctes.m_p.value)

  #ISOTERMICA
  Mtot=-aux*del_lnpho #kg
  Mtot = ((Mtot *u.kg).to(u.M_sun)).value/1e13
  return Mtot


#perfil massa total

del_lnpho_prof=der_ln_pho(r_cm, rc_den_cm, rs_den_cm, beta_den, gamma_den, alpha_den,eps_den)
# del_lnT_prof=der_ln_T(r_cm, rc_tem_cm, rs_tem_cm, beta_tem, eta_tem, alpha_tem,eps_tem)

r_m = (r_cm *u.cm).to(u.m).value
aux_prof=(kT_si*r_m)/(ctes.G.value*mu*ctes.m_p.value)

Mtot_prof=-aux_prof*del_lnpho_prof #kg
Mtot_prof = ((Mtot_prof *u.kg).to(u.M_sun)).value

# Mtot_prof2=-aux_prof*(del_lnpho_prof + del_lnT_prof) #kg
# Mtot_prof2 = ((Mtot_prof2 *u.kg).to(u.M_sun)).value



Yx_tot=kT*Mgas
Mmaug_tot = M500_maug(Yx_tot,z)





#===================================================================
#    PLOTA MASSA
#===================================================================
r_kpc= (r_cm *u.cm).to(u.kpc).value
r_r500=r_cm/R500_cm

plt.plot(r_cm, Mgas_prof, marker='.', linestyle='--',label='Massa do gás', color='indianred')
plt.plot(r_cm, Mtot_prof, marker='.', linestyle='--',label='Massa total (T constante)', color='deepskyblue')
# plt.plot(r_r500, Mtot_prof2, marker='.', linestyle='--',label='Massa total (T varia)', color='olivedrab')

plt.suptitle('Perfil de massa - Abell 1795')
plt.legend(loc='best')
plt.yscale('log')
plt.ylabel('Massa [$M_{\odot}$]')
plt.xlabel('R (arcmin)')
# plt.xlabel('R/$R_{500}$ (kpc)')
plt.xscale('log')

plt.show()

nome_fig=nome+'_massa_prof.png'
plt.savefig(nome_fig)


