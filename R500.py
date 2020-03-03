#usando do valore de kT, acha R500

import sys
import numpy as np


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass


#ENTRADAS
if len(sys.argv) < 4:
  kT=float(sys.argv[1])
  z=float(sys.argv[2])
  kTneg=0.05*kT
  kTpos=0.05*kT
else:
  kT=float(sys.argv[1])
  z=float(sys.argv[2])
  kTneg=float(sys.argv[3])
  kTpos=float(sys.argv[4])



cteEr=0.09

#FUNCAO R500

def r500(kT,z):

  z_aux=(1+z)**(-3/2)
  kT_aux=np.sqrt(kT/10)

  return 1.24*kT_aux*z_aux

def r500_Er(kT, z, kTneg, kTpos, cteEr):
  z_aux=(1+z)**(-3./2.)
  kT_aux=np.sqrt(kT/10)

  deriv_cte=kT_aux*z_aux
  deriv_kT=1.24*kT_aux*z_aux/(2*(kT**2))

  aux_err_cte=(deriv_cte**2)*(cteEr**2)
  aux_err_kTneg=(deriv_kT**2)*(kTneg**2)
  aux_err_kTpos=(deriv_kT**2)*(kTpos**2)

  r500_neg=np.sqrt(aux_err_cte + aux_err_kTneg)
  r500_pos=np.sqrt(aux_err_cte + aux_err_kTpos)

  return [r500_neg,r500_pos]

#chama as funcoes

R500 = r500(kT,z)
R500Neg,R500Pos = r500_Er(kT, z, kTneg, kTpos, cteEr)

R500_015=0.15*R500
R500_015Neg=0.15*R500Neg
R500_015Pos=0.15*R500Pos

sai=open('R500.txt', 'w')
sai.write('R500_1 = %s'%R500)
sai.write(' (-%s'%R500Neg)
sai.write(',%s'%R500Pos)
sai.write(")")
sai.write('h-1 Mpc')
sai.write('\n')
sai.write('R500_0.15 = %s'%R500_015)
sai.write(' (-%s'%R500_015Neg)
sai.write(',%s'%R500_015Pos)
sai.write(")")
sai.write('h-1 Mpc')
sai.write('\n')


sai.close()

