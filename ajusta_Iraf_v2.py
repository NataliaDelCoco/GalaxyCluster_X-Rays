
#Natalia Crepaldi


#cria uma imagem corrigida do exposure map
#Tira o brilho superficial em elipses usando o iraf
#ajusta modelo bets 2D

#entradas:
#      1)imagem do aglomerado
#      2)exposure map
#      4)imagem cheese
#      5) mascara cheese

import os
import sys
from astropy.io import fits
import numpy as np
import csv
import pandas as pd
import stsci.tools

local=os.getcwd()

#entra no home/natalia pra inicializar
os.chdir('/home/natalia')

from pyraf import iraf
from iraf import stsdas
from iraf import analysis
from iraf import isophote
from iraf import images
from iraf import listpixels
from iraf import noao
from iraf import imred
from iraf import ccdred
from iraf import badpiximage
from iraf import ccdmask
from iraf import imutil
from iraf import tv
from iraf import imedit 
from iraf import proto
from iraf import fixpix
from iraf import mskexpr
from iraf import mskregions

iraf.ellipse.unlearn()
iraf.magpar.unlearn()
iraf.geompar.unlearn()
iraf.samplepar.unlearn()
iraf.imarith.unlearn()
# iraf.display.unlearn()
iraf.isoplot.unlearn()


#vai pro diretorio 'analysis'
os.chdir(local)

#======================================================
#  ENTRADAS

fich_nps_m1=sys.argv[1] #imagem mos1 sem fontes pontuais
fich_ps_m1=sys.argv[2] #imagem mos1 com fontes pontuais
exp_m1 = sys.argv[3] #exp_map mos1
chee_m1= sys.argv[4]# cheese mos1
fich_nps_m2=sys.argv[5] #imagem mos2 sem fontes pontuais
fich_ps_m2=sys.argv[6] #imagem mos2 com fontes pontuais
exp_m2 = sys.argv[7] #exp_map mos2
chee_m2= sys.argv[8]# cheese mos2

x0=int(sys.argv[9]) #centro x da iamgem
y0=int(sys.argv[10]) #centro y da iamgem

nome='m1m2_expcor_norm_cor'

fich_comb = nome+".fit"
fich_ch="cheese_99.fits"
fich_mask = "cheese&.pl"
fich4 = "blu.txt"
tab   = nome+".tab"
tabcsv = nome+".csv"
model = nome+"_Model.fits"
resid = nome+"_Resid.fits"


arqs=[tab,tabcsv,model,resid, fich4, fich_mask, fich_ch]

for item in arqs:
  if os.path.exists(item):
    os.remove(item)

#======================================================
#  IMAGEM CORRIGIDA DO EXPMAP

# MOS1------------------------------------------------------------
#compara valores maximos da imagem sem point sources e do exposure map
hdul_nps=fits.open(fich_nps_m1)
hdul_exp=fits.open(exp_m1)
rat_ps_exp=np.amax(hdul_nps[0].data)/np.amax(hdul_exp[0].data)

#divide imagem com point sources pelo exposure map

iraf.imarith(operand1=fich_ps_m1, op="/",  operand2=exp_m1, result='mos1_expcor.fits')
#tem que normalizar as contagens da imagem corrigida pelo exp_map.
#vamos normalizar usando a razao rat_ps_exp

hdul_expcor=fits.open('mos1_expcor.fits')
expcor_data=hdul_expcor[0].data


rat_ps_expcor=np.amax(hdul_nps[0].data)/np.amax(expcor_data)
norm_cor=expcor_data*rat_ps_expcor

norm=expcor_data/rat_ps_exp

hdul_norm=fits.PrimaryHDU(norm_cor) 
hdul_norm.writeto('mos1_expcor_norm.fits', overwrite=True)

# MOS2------------------------------------------------------------
#compara valores maximos da imagem sem point sources e do exposure map
hdul_nps=fits.open(fich_nps_m2)
hdul_exp=fits.open(exp_m2)
rat_ps_exp=np.amax(hdul_nps[0].data)/np.amax(hdul_exp[0].data)

#divide imagem com point sources pelo exposure map

iraf.imarith(operand1=fich_ps_m2, op="/",  operand2=exp_m2, result='mos2_expcor.fits')
#tem que normalizar as contagens da imagem corrigida pelo exp_map.
#vamos normalizar usando a razao rat_ps_exp

hdul_expcor=fits.open('mos2_expcor.fits')
expcor_data=hdul_expcor[0].data


rat_ps_expcor=np.amax(hdul_nps[0].data)/np.amax(expcor_data)
norm_cor=expcor_data*rat_ps_expcor


norm=expcor_data/rat_ps_exp

hdul_norm=fits.PrimaryHDU(norm_cor) 
hdul_norm.writeto('mos2_expcor_norm.fits', overwrite=True)

# SOMA IMAGENS EXPCOR_NORM MOS1 E MOS2-------------------------------
hdul_norm_m1=fits.open('mos1_expcor_norm.fits')
hdul_norm_m2=fits.open('mos2_expcor_norm.fits')
m1_data = hdul_norm_m1[0].data
m2_data = hdul_norm_m2[0].data

#soma
final=m1_data+m2_data
hdul_final = fits.PrimaryHDU(final)
hdul_final.writeto(fich_comb, overwrite=True)
#======================================================
# ARRUMA MASCARA

iraf.imarith(operand1=chee_m1, op="*", operand2=chee_m2, result="chm1m2.fits")
# iraf.imarith(operand1="chm1m2.fits", op="*", operand2=chpn, result="chtot.fits")

hd_orig=fits.open(fich_ps_m2)
hd_chm=fits.open("chm1m2.fits")

orig_data=hd_orig[0].data
chm_data=hd_chm[0].data


#cria uma nova imagem, colocando 99999 onde eh observacao e 0 onde eh mask
lin=0
for row in chm_data:
  col=0 
  for elem in row:
     if elem==1:      
         orig_data[lin,col]=int(999999)
     col=col+1                       
  lin=lin+1 



hd_orig[0].data=orig_data
hd_orig.writeto(fich_ch)
hd_orig.close()

#cria mascara.pl
# bad: 0 -> 1
# god: 999999 -> 0
iraf.imexpr(expr='(a == 0) ? 1 : 0',output=fich_mask, a=fich_ch)

#======================================================

##  number of sigma-clip iterations
iraf.samplepar.nclip=12

## Parametros da geometria EM IMG
iraf.geompar.x0 = x0
iraf.geompar.y0 = y0
iraf.geompar.ellip0 = 0.1
iraf.geompar.sma0 = 15

## -----  nao ajustamos em menos de 1 pixels  ##
iraf.geompar.minsma = 1

## -----  este valor pode ser maior ou menor; depende da imagem  ##
iraf.geompar.maxsma=1000

## Muda parametros default do display para podermos ver o centro da galaxia
#so descomentar quando quiser visualizar

# !ds9&
# iraf.display.zscale='no'
# iraf.display.ztrans="log"

## =====>>>>  Roda o stsdas/ellipse
print(x0)
print(y0)
iraf.ellipse(input=fich_comb, output=tab, interactive='no', masksz=1, region='yes')

## ====>>>>  Constroi um modelo baseado no ellipse
iraf.bmodel(table=tab,output=model)
iraf.imarith(operand1=fich_comb,  op="-", operand2=model,  result=resid)

## ====>>>> Faz o grafico mag X R^(1/4) (perfil de de Vaucouleurs)
# iraf.isoplot(input=tab, xaxis="RSMA", yaxis="INTENS")

### ------  Saida de uma tabela ascii (.txt) para ajuste em outro programa


colunas="SMA,RSMA,INTENS,INT_ERR,RMS,MAG,MAG_LERR,MAG_UERR"
colunas_2=['SMA','RSMA','INTENS','INT_ERR','RMS','MAG','MAG_LERR','MAG_UERR']

iraf.tdump(table=tab, datafile=fich4, columns=colunas)

#ACABA AQUI 


data=pd.read_csv(fich4, sep="\s+", header=None)
data.columns=colunas_2
data.to_csv(tabcsv)



### Apaga arquivos temporarios

# os.remove(fich2)
os.remove(fich4)
os.remove(tab)

print('Terminou. Veja os arquivos ',nome,' Model e Resid.fits')


# transforma de .txt pra .csv

