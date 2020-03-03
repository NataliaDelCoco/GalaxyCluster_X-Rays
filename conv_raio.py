import numpy as np
import sys

conv=float(sys.argv[1]) #kpc/arcsec
r_init=float(sys.argv[2]) #arcmin ou Mpc
convers=sys.argv[3] #(arcmin2Mpc ou Mpc2arcmin)

#acha quantos arcmin equivalem a 50 Mpc

#Mpc=50000
#Mpc=3

def Mpc2arc(conv,Mpc):
  kpc=Mpc*1000
  arcmin=(kpc/60)/conv
  return arcmin

def arc2Mpc(conv,arc):
  Mpc=arc*conv*60/1000
  return Mpc


#+++++++++++++++++++++++++++++++

if convers == 'arcmin2Mpc' :
  r_final = arc2Mpc(conv,r_init)
else:
  r_final = Mpc2arc(conv,r_init)

print(r_final)
