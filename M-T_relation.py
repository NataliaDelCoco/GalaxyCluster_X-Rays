import sys
import numpy as np
import astropy.units as un
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70 * un.km / un.s / un.Mpc, Om0=0.3)

kT=float(sys.argv[1])
kT_err1=float(sys.argv[2])
kT_err2=float(sys.argv[3])
z=float(sys.argv[4])

#fósmula tirada da apostila do gastao, que veio de Evrard et al 1996

cte=2.2158*1e15/1.4

mass_aux=kT/((1+z)*10)
mass=mass_aux**(3/2)

M500=cte*mass/1e13

#Erro
aux_err = cte*1.5*(mass_aux**(1/2))*(mass_aux/kT)
M500_err1=np.sqrt((aux_err**2)*(kT_err1**2))/1e13
M500_err2=np.sqrt((aux_err**2)*(kT_err2**2))/1e13

print ('M500 = %.5f -+ (%.5f, %.5f) 1e13Msun' %(M500, M500_err1, M500_err2))


#====================
densi_c=cosmo.critical_density(z)
Msun_Mpc3=un.M_sun/((1e6*un.pc)**3)
densi_c = (densi_c.to(Msun_Mpc3)).value

cte2=3.993*1e3
cte2_err=0.934*1e3

mass_aux=kT/((1+z)*10)
mass=mass_aux**(3/2)
M500_2 = (cte2*densi_c*(mass_aux**1.5)*((0.7*(1+z))**(-3)))

#erro
M500_err1=M500_2*np.sqrt((cte2_err/cte2)**2 + (1.5*kT_err1/kT)**2)/1e13
M500_err2=M500_2*np.sqrt((cte2_err/cte2)**2 + (1.5*kT_err2/kT)**2)/1e13

M500_2=M500_2/1e13

print ('M5002 = %.5f -+ (%.5f, %.5f) 1e13Msun' %(M500_2, M500_err1, M500_err2))