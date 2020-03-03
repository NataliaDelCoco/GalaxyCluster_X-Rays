#!/bin/bash

#ACHA O CENTRO 
#usando o centro de referencia, que a gente lê no arquivo de eventos, 
#fazemos 6 iteracoes pra achar o centro real do aglomerado

# centro_ref=$cluster'_centro_ref.txt'
# x_cen=$(cat $centro_ref | grep X_ref | cut -d ' ' -f 2)
# y_cen=$(cat $centro_ref | grep Y_ref | cut -d ' ' -f 2)
# r_cen=$(cat $centro_ref | grep R_ref | cut -d ' ' -f 2)


cam=$1
echo $cam
x_cen=$(fkeyprint $cam'_filt.fits'"[0]" REFXCRPX | grep WCS | awk '{print $2}')
y_cen=$(fkeyprint $cam'_filt.fits'"[0]" REFYCRPX | grep WCS | awk '{print $2}')
echo $x_cen $y_cen
pix_ref=$(fkeyprint $cam'_filt.fits'"[0]" REFYCDLT | grep WCS | awk '{print $2}')
pix_ref=$(echo $pix_ref | awk '{print $1 * 3600.}') #segundo de arco

rin=10.0 #minutos de arco
rinPhy_ref=$(echo $rin $pix_ref | awk '{print $1*60/$2}')
r_cen=$rinPhy_ref
i=1
#-------------LOOP---------------
while [ $i -le 6 ]; do

#echo $i
regiao_cen="((X,Y) in CIRCLE($x_cen,$y_cen,$r_cen))"
echo $regiao_cen
eregionanalyse imageset=$pn'_img_sembin'.fits srcexp="$regiao_cen" \
  exposuremap=$pn'_expmap_sembin'.fits > saida.txt 

x_cen=$(cat saida.txt | grep -v Executing | grep xcentroid | cut -d ' ' -f 3)
y_cen=$(cat saida.txt | grep -v Executing | grep ycentroid | cut -d ' ' -f 3)
r_cen=$(echo "$r_cen/1.75" | bc)


let i=i+1

done #while

arq_sai=$cam'_centro_PHY.reg'

echo $x_cen $y_cen $r_cen
echo "(X,Y) in circle("$x_cen","$y_cen","$r_cen")" >>  $arq_sai
rm saida.txt
