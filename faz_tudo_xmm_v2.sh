#!/bin/bash

#----------------------
#ENTRADAS
echo 'cluster name:'
read cluster
#limites de energia
echo "Emin (eV):"
read PImin #500
echo "Emax (eV):"
read PImax #10000
echo "Bin size:"
read bin #128



#VALORES DOS RAIOS [ARCMIN]
r_ext_agl=8. #raio externo da reg usada pra spec do agl
r_cen_agl=50 #kpc
r_int_bkg=12. #raio interno do anel pra analise do bkg
r_ext_bkg=14. #raio externo do anel para analise do bkg
#de arcmin pra physical r*1200

r_ext_agl=$(bc <<< "scale=1; $r_ext_agl*1200.")
#r_cen_agl=$(bc <<< "scale=1; $r_cen_agl*1200.")
r_int_bkg=$(bc <<< "scale=1; $r_int_bkg*1200.")
r_ext_bkg=$(bc <<< "scale=1; $r_ext_bkg*1200.")

#------------------------
#DOWNLOAD E PREPARACAO dos arquivos e dos diretorios

path="/home/natalia/Dados"
mkdir -p $path/$cluster

catalogo=/home/natalia/Dados/XMM_catalogues/Catalogue_xmm_ned.csv
obsid=$(cat $catalogo |grep $cluster | cut -d ',' -f 1)
echo "obsid="$obsid
len=${#obsid}
dif=$(echo $len | awk '{print(10-$1)}')
if [ "$dif" -ne "0" ]; then
  while [ $dif -gt 0 ]; do
    obsid=0$obsid
    dif=$(echo $dif | awk '{print($1-1)}')
  done
fi



outdir=/home/natalia/Dados/XMM/clusters
cd $outdir
nome_file=$obsid'.tar'
curl -o $nome_file 'http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno='$obsid

tar -xvf $nome_file

mv $obsid $path/$cluster
cd ../../$cluster

arq=$(ls)


export SAS_CCFPATH=/home/natalia/Dados/ccf

local=$(pwd)


dir=${arq%.tar.gz}
# mv $arq $local/XMM/clusters/ $path/$cluster/
# cd $path/$cluster/
mkdir -p $dir
cd $dir
mkdir odf analysis
cd ..
mv $arq $dir/odf/
cd $dir/odf

tar -xzf $arq.tar.gz
tar_file=$(ls *.TAR)
tar -xvf $tar_file
mv $arq.tar.gz $outdir
rm $tar_file

local=$(pwd)
export SAS_ODF=$local


cd ../analysis/
cifbuild
local=$(pwd)
export SAS_CCF=$local/ccf.cif

cd ../odf/
odfingest
SUM=$(ls *SUM.SAS)
local=$(pwd)
export SAS_ODF=$local/$SUM

cd ../analysis/
path_anl=$(pwd)
#-----------------------------------------------------
#REDUCAO

m1=$cluster'_mos1'
m2=$cluster'_mos2'
pn=$cluster'_pn'



emchain
epchain
epchain withoutoftime=yes

#curva de luz
pn-filter
mos-filter


#filtragens-----------------------
m1_clean=$(ls mos1*clean.*)
m2_clean=$(ls mos2*clean.*)
pn_clean=$(ls pn*clean.*)

evselect table=$m1_clean filteredset=$m1'_filt'.fits withfilteredset=yes filtertype=expression \
expression="(PATTERN<=12)&&(FLAG==0)" \
keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=$m2_clean filteredset=$m2'_filt'.fits withfilteredset=yes filtertype=expression \
expression="(PATTERN<=12)&&(FLAG==0)" \
keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=$pn_clean filteredset=$pn'_filt'.fits withfilteredset=yes filtertype=expression \
expression="(PATTERN==0)&&(FLAG==0)" \
keepfilteroutput=yes updateexposure=yes filterexposure=yes

#imagens filtradas
#VAMOS FAZER TUDO EM COORDENADAS DO DETECTOR

evselect table=$m1'_filt'.fits  imageset=$m1'_img'.fits imagebinning='binSize' withimageset=yes \
xcolumn='DETX' ycolumn='DETY' ximagebinsize=$bin yimagebinsize=$bin \
expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN in [0:4])"

evselect table=$m2'_filt'.fits imageset=$m2'_img'.fits imagebinning='binSize' withimageset=yes \
xcolumn='DETX' ycolumn='DETY' ximagebinsize=$bin yimagebinsize=$bin \
expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN in [0:4])"


evselect table=$pn'_filt'.fits imageset=$pn'_img'.fits imagebinning='binSize' withimageset=yes \
xcolumn='DETX' ycolumn='DETY' ximagebinsize=$bin yimagebinsize=$bin \
expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN==0)"

#atitude do satelite
atthkgen atthkset=atti.dat timestep=1
#mapa de exposição

eexpmap imageset=$m1'_img'.fits eventset=$m1'_filt'.fits expimageset=$m1'_expmap'.fits  attitudeset='atti.dat' \
pimin=$PImin pimax=$PImax withdetcoords=true withvignetting=true

eexpmap imageset=$m2'_img'.fits eventset=$m2'_filt'.fits expimageset=$m2'_expmap'.fits  attitudeset='atti.dat' \
pimin=$PImin pimax=$PImax withdetcoords=true withvignetting=true

eexpmap imageset=$pn'_img'.fits eventset=$pn'_filt'.fits expimageset=$pn'_expmap'.fits  attitudeset='atti.dat' \
pimin=$PImin pimax=$PImax withdetcoords=true withvignetting=true


#fontes pontuais-------------------
#TEM QUE GERAR IMAGENS =>NAO<= BINADAS
#ESPECIFICAMENTE AQUI USAREMOS =>COORDS PHY<=
#POR CAUSA DO ECOORDCONV


#cria imagem
#limites de energia das imagens:
#m1 e m2 = [500:10000] 
evselect table=$m1'_filt'.fits  imageset=$m1'_img_sembin'.fits withimageset=yes \
xcolumn='X' ycolumn='Y' expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN in [0:4])"

evselect table=$m2'_filt'.fits imageset=$m2'_img_sembin'.fits withimageset=yes \
xcolumn='X' ycolumn='Y' expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN in [0:4])"


evselect table=$pn'_filt'.fits imageset=$pn'_img_sembin'.fits withimageset=yes \
xcolumn='X' ycolumn='Y' expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN==0)"

#mapa de exposicao
eexpmap imageset=$m1'_img_sembin'.fits eventset=$m1'_filt'.fits expimageset=$m1'_expmap_sembin'.fits  attitudeset='atti.dat' \
pimin=$PImin pimax=$PImax

eexpmap imageset=$m2'_img_sembin'.fits eventset=$m2'_filt'.fits expimageset=$m2'_expmap_sembin'.fits  attitudeset='atti.dat' \
pimin=$PImin pimax=$PImax

eexpmap imageset=$pn'_img_sembin'.fits eventset=$pn'_filt'.fits expimageset=$pn'_expmap_sembin'.fits  attitudeset='atti.dat' \
pimin=$PImin pimax=$PImax

#ewavelet
ewavelet imageset=$m1'_img_sembin'.fits expmapset=$m1'_expmap_sembin'.fits srclistset=$m1'_wave'.fits savebkgrnd=yes \
bkgrndset=$m1'_wave_back'.fits makerecon=yes recimageset=$m1'_wave_rec'.fits

ewavelet imageset=$m2'_img_sembin'.fits expmapset=$m2'_expmap_sembin'.fits srclistset=$m2'_wave'.fits savebkgrnd=yes \
bkgrndset=$m2'_wave_back'.fits makerecon=yes recimageset=$m2'_wave_rec'.fits

ewavelet imageset=$pn'_img_sembin'.fits expmapset=$pn'_expmap_sembin'.fits srclistset=$pn'_wave'.fits savebkgrnd=yes \
bkgrndset=$pn'_wave_back'.fits makerecon=yes recimageset=$pn'_wave_rec'.fits \
threshold=6 minscale=3

#esse programa verifica quais fontes pontuais realmente existem
#e verifica o centro de referencia de cada observacao
#saidas:
#cluster_mos1/mos2/pn_pointsources.fits/txt
#cluster_centro_ref.txt

#UNIDADES DOS ARQES DE SAIDA EM "IMAGEM"
python /home/natalia/Codigos/compara_fontes.py 



#ACHA CENTRO => usando pn-------------------------------------------

source /home/natalia/Codigos/acha_centro.sh $pn

bash /home/natalia/Codigos/transforma_coord.sh $pn'_filt.fits' centro_PHY.reg DET 128 coord centro_DET.reg

x_cen=$(cat centro_DET.reg | grep circle| cut -d '(' -f 2 | cut -d ',' -f 1)
y_cen=$(cat centro_DET.reg | grep circle| cut -d '(' -f 2 | cut -d ',' -f 2)
r_cen=$(cat centro_DET.reg | grep circle| cut -d '(' -f 2 | cut -d ',' -f 3 | cut -d ')' -f 1)
#-----------------------------


#chama programa que tranforma as coordenadas das fontes pontuais a serem mascaradas
#saidas: cluster_camera_coord_ps.reg

bash /home/natalia/Codigos/transforma_coord.sh $m1'_filt.fits' $m1'_pointsources.txt' DET 128 tab $m1'_DET_ps.reg'
bash /home/natalia/Codigos/transforma_coord.sh $m2'_filt.fits' $m2'_pointsources.txt' DET 128 tab $m2'_DET_ps.reg'
bash /home/natalia/Codigos/transforma_coord.sh $pn'_filt.fits' $pn'_pointsources.txt' DET 128 tab $pn'_DET_ps.reg'



#cria anel de [14,12] min 

echo "((DETX,DETY) in circle("$x_cen","$y_cen","$r_ext_bkg"))&&!((DETX,DETY) \
in circle("$x_cen","$y_cen","$r_int_bkg"))" > anel_bkg.reg

#----------------------------------
#cria imagem da total e da regiao a ser analisada, MASCARANDO AS PS


ps_m1=$(< $m1'_DET_ps_excl.reg')
ps_m2=$(< $m2'_DET_ps_excl.reg')
ps_pn=$(< $pn'_DET_ps_excl.reg')


#TOTAL SEM PS
evselect table=$m1'_filt'.fits  imageset=$m1'_img_semps'.fits imagebinning='binSize' withimageset=yes \
xcolumn='DETX' ycolumn='DETY' ximagebinsize=$bin yimagebinsize=$bin \
expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN in [0:4])$ps_m1"

evselect table=$m2'_filt'.fits  imageset=$m2'_img_semps'.fits imagebinning='binSize' withimageset=yes \
xcolumn='DETX' ycolumn='DETY' ximagebinsize=$bin yimagebinsize=$bin \
expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN in [0:4])$ps_m2"

evselect table=$pn'_filt'.fits  imageset=$pn'_img_semps'.fits imagebinning='binSize' withimageset=yes \
xcolumn='DETX' ycolumn='DETY' ximagebinsize=$bin yimagebinsize=$bin \
expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN == 0)$ps_pn"


#background------------------------
#verificar no arquivo de eventos a data e hora & o tipo de filtro (thin, medium, thick)


path_back=/home/natalia/Dados/XMM/blank_sky

date=$(fkeyprint $m1'_filt.fits'"[0]" DATE-OBS | grep Start | cut -d ' ' -f 2)
filter=$(fkeyprint $m1'_filt.fits'"[0]" FILTER | grep ID | awk '{print $3}' | cut -d "'" -f 2)
rev=$(fkeyprint $m1'_filt.fits'"[0]" REVOL | grep Revolution | cut -d '=' -f 2 | cut -d '/' -f 1) # revolucao


if [ "$filter" == "Thick" ]
then
  f=k
elif [ "$filter" == "Thin" ]
then
  f=t
else
  f=m
fi

if [ $rev -gt 0961 ] #rev=0961 => quebrou ccd6
then #se nao tem a ccd6
  b_m1='m1'$f'ffg_nccd6_events'
else #se tem a ccd6
  b_m1='m1'$f'ffg_wccd6_events'
fi

b_m2='m2'$f'ffg_events'
b_pn='pn'$f'ffg_events'

back_m1=$b_m1'.fits'
back_m2=$b_m2'.fits'
back_pn=$b_pn'.fits'

#apos achar as blanck_files certas, copia elas pra dentro da pasta Analysis

cp $path_back'/'$back_m1 $path_anl
cp $path_back'/'$back_m2 $path_anl
cp $path_back'/'$back_pn $path_anl

mv $back_m1 $m1'_BLANK.fits'
mv $back_m2 $m2'_BLANK.fits'
mv $back_pn $pn'_BLANK.fits'


fparkey $date $m1'_BLANK.fits'+0 DATE-OBS add=yes
fparkey $date $m1'_BLANK.fits'+1 DATE-OBS add=yes
fparkey $date $m2'_BLANK.fits'+0 DATE-OBS add=yes
fparkey $date $m2'_BLANK.fits'+1 DATE-OBS add=yes
fparkey $date $pn'_BLANK.fits'+0 DATE-OBS add=yes
fparkey $date $pn'_BLANK.fits'+1 DATE-OBS add=yes

#filtra as blanck_files do mesmo jeito que filtrou os arquivos de eventos 

evselect table=$m1'_BLANK.fits' filteredset=$m1'_BLANK_filt.fits' withfilteredset=yes filtertype=expression \
expression="(PATTERN<=12)&&(FLAG==0)" \
keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=$m2'_BLANK.fits' filteredset=$m2'_BLANK_filt.fits' withfilteredset=yes filtertype=expression \
expression="(PATTERN<=12)&&(FLAG==0)" \
keepfilteroutput=yes updateexposure=yes filterexposure=yes

evselect table=$pn'_BLANK.fits' filteredset=$pn'_BLANK_filt.fits' withfilteredset=yes filtertype=expression \
expression="(PATTERN==0)&&(FLAG==0)" \
keepfilteroutput=yes updateexposure=yes filterexposure=yes

#coloca as blank_files na mesma posição e coords que as observações
#MOS1
file=$m1'_filt'.fits
ra_pnt=$(fkeyprint $file"[0]" RA_PNT | grep = | cut -d ' ' -f 4)
dec_pnt=$(fkeyprint $file"[0]" DEC_PNT | grep = | cut -d ' ' -f 3)
pa_pnt=$(fkeyprint $file"[0]" PA_PNT | grep = | cut -d ' ' -f 4)
ranom=$(fkeyprint $file"[0]" REFXCRVL | grep = | cut -d ' ' -f 2)
decnom=$(fkeyprint $file"[0]" REFYCRVL | grep = | cut -d ' ' -f 2)

outfile=$m1'_BLANK_filt_correctpos.fits'
cp $m1'_BLANK_filt.fits' $outfile
attcalc -w 0 -V 0 eventset=$outfile attitudelabel=fixed fixedra=$ra_pnt fixeddec=$dec_pnt \
fixedposangle=$pa_pnt withatthkset=N refpointlabel=user nominalra=$ranom nominaldec=$decnom
fparkey $ra_pnt $outfile+0 RA_PNT add=yes
fparkey $ra_pnt $outfile+1 RA_PNT add=yes
fparkey $dec_pnt $outfile+0 DEC_PNT add=yes
fparkey $dec_pnt $outfile+1 DEC_PNT add=yes
fparkey $pa_pnt $outfile+0 PA_PNT add=yes
fparkey $pa_pnt $outfile+1 PA_PNT add=yes

unset file outfile

#MOS2
file=$m2'_filt'.fits
ra_pnt=$(fkeyprint $file"[0]" RA_PNT | grep = | cut -d ' ' -f 4)
dec_pnt=$(fkeyprint $file"[0]" DEC_PNT | grep = | cut -d ' ' -f 3)
pa_pnt=$(fkeyprint $file"[0]" PA_PNT | grep = | cut -d ' ' -f 4)
ranom=$(fkeyprint $file"[0]" REFXCRVL | grep = | cut -d ' ' -f 2)
decnom=$(fkeyprint $file"[0]" REFYCRVL | grep = | cut -d ' ' -f 2)

outfile=$m2'_BLANK_filt_correctpos.fits'
cp $m2'_BLANK_filt.fits' $outfile
attcalc -w 0 -V 0 eventset=$outfile attitudelabel=fixed fixedra=$ra_pnt fixeddec=$dec_pnt \
fixedposangle=$pa_pnt withatthkset=N refpointlabel=user nominalra=$ranom nominaldec=$decnom
fparkey $ra_pnt $outfile+0 RA_PNT add=yes
fparkey $ra_pnt $outfile+1 RA_PNT add=yes
fparkey $dec_pnt $outfile+0 DEC_PNT add=yes
fparkey $dec_pnt $outfile+1 DEC_PNT add=yes
fparkey $pa_pnt $outfile+0 PA_PNT add=yes
fparkey $pa_pnt $outfile+1 PA_PNT add=yes

unset file outfile

#pn
file=$pn'_filt'.fits
ra_pnt=$(fkeyprint $file"[0]" RA_PNT | grep = | cut -d ' ' -f 4)
dec_pnt=$(fkeyprint $file"[0]" DEC_PNT | grep = | cut -d ' ' -f 3)
pa_pnt=$(fkeyprint $file"[0]" PA_PNT | grep = | cut -d ' ' -f 4)
ranom=$(fkeyprint $file"[0]" REFXCRVL | grep = | cut -d ' ' -f 2)
decnom=$(fkeyprint $file"[0]" REFYCRVL | grep = | cut -d ' ' -f 2)

outfile=$pn'_BLANK_filt_correctpos.fits'
cp $pn'_BLANK_filt.fits' $outfile
attcalc -w 0 -V 0 eventset=$outfile attitudelabel=fixed fixedra=$ra_pnt fixeddec=$dec_pnt \
fixedposangle=$pa_pnt withatthkset=N refpointlabel=user nominalra=$ranom nominaldec=$decnom
fparkey $ra_pnt $outfile+0 RA_PNT add=yes
fparkey $ra_pnt $outfile+1 RA_PNT add=yes
fparkey $dec_pnt $outfile+0 DEC_PNT add=yes
fparkey $dec_pnt $outfile+1 DEC_PNT add=yes
fparkey $pa_pnt $outfile+0 PA_PNT add=yes
fparkey $pa_pnt $outfile+1 PA_PNT add=yes

unset file outfile



aux_bkg=$(<anel_bkg.reg)
reg_m1_bkg=$(echo $aux_bkg $ps_m1 )
reg_m2_bkg=$(echo $aux_bkg $ps_m2 )
reg_pn_bkg=$(echo $aux_bkg $ps_pn )


# echo 'Qual o corte inferior de energia para ajuste espectral?'
# read E_min

#    EXTRAI ESPECTROS E ACHA R500-----------------------------------
#
E_min=0.5
nh 2000.0 $ranom $decnom >> $cluster'_nh.txt'
nhh=$( cat $cluster'_nh.txt' | grep Weighted | awk '{print $7}')
nh_n=$(echo $nhh |cut -d E -f 1)
rm $cluster'_nh.txt'
#colocando em unidades de E+22

nh_22=$(bc <<< "scale=4; $nh_n/100")

xs_m1=$m1'_spec_v2_gr.fits'
xs_m2=$m2'_spec_v2_gr.fits'
xs_pn=$pn'_spec_v2_gr.fits'

#REDSHIFT
file_z=/home/natalia/Dados/XMM_catalogues/Catalogue_xmm_ned.csv
z=$(cat $file_z | grep $cluster  |cut -d ',' -f 12)

#arrumando o raio (chute inicial)
R_kpc=$(echo $r_cen_agl)
r0_Mpc=$(echo $R_kpc | awk '{print($1/1000)}')


continua=sim
conta=0

if [ "$continua" == "sim" ]; then

  #-----AGLOMERADOS---------------------------------
  #   refaz as regioes
  rm r_arcmin.txt -f
  ipython << ULA
import pandas as pd
import numpy as np
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

conv=cosmo.angular_diameter_distance($z)
conv=conv.value

div=conv*0.00029088820866572
R_arcmin=$r0_Mpc/div

sai=open('r_arcmin.txt', 'w')

sai.write("r_arcmin=%s" % R_arcmin)
sai.close()
ULA
  r_arcmin=$(cat r_arcmin.txt | cut -d '=' -f 2)
  #transforma de arcmin pra DET
  r0_DET=$(bc  <<<  "scale=8; $r_arcmin*1200.")


  echo "((DETX,DETY) in circle("$x_cen","$y_cen","$r_ext_agl"))&&!((DETX,DETY) \
  in circle("$x_cen","$y_cen","$r0_DET"))" > aglomerado.reg

  aux_ag=$(<aglomerado.reg)
  reg_m1=$(echo $aux_ag $ps_m1)
  reg_m2=$(echo $aux_ag $ps_m2)
  reg_pn=$(echo $aux_ag $ps_pn)

  #   espectro
  evselect table=$m1'_filt'.fits':EVENTS' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$reg_m1" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$m1'_spec_v2.fits' spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999

  evselect table=$m2'_filt'.fits destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$reg_m2" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$m2'_spec_v2.fits' spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999

  evselect table=$pn'_filt'.fits destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$reg_pn" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$pn'_spec_v2.fits' spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479
 
  #-----BLANKS---------------------------------
  #   espectro
  evselect table=$m1'_BLANK_filt_correctpos.fits' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$reg_m1" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$m1'_BLANK_spec_v2.fits' spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999

  evselect table=$m2'_BLANK_filt_correctpos.fits' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$reg_m2" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$m2'_BLANK_spec_v2.fits' spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999

  evselect table=$pn'_BLANK_filt_correctpos.fits' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$reg_pn" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$pn'_BLANK_spec_v2.fits' spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479

  #backscale
  backscale spectrumset=$m1'_spec_v2.fits' badpixlocation=$m1'_filt.fits'
  backscale spectrumset=$m2'_spec_v2.fits' badpixlocation=$m2'_filt.fits'
  backscale spectrumset=$pn'_spec_v2.fits' badpixlocation=$pn'_filt.fits'

  backscale spectrumset=$m1'_BLANK_spec_v2.fits' badpixlocation=$m1'_filt.fits'
  backscale spectrumset=$m2'_BLANK_spec_v2.fits' badpixlocation=$m2'_filt.fits'
  backscale spectrumset=$pn'_BLANK_spec_v2.fits' badpixlocation=$pn'_filt.fits'

  #   rmf
  rmfgen rmfset=$m1'v2.rmf' withenergybins=no  spectrumset=$m1'_spec_v2.fits' format='var' detmaptype='flat'
    rmf_m1=$m1'v2.rmf'

  rmfgen rmfset=$m2'v2.rmf' withenergybins=no  spectrumset=$m2'_spec_v2.fits' format='var' detmaptype='flat'
    rmf_m2=$m2'v2.rmf'

  rmfgen rmfset=$pn'v2.rmf' withenergybins=no  spectrumset=$pn'_spec_v2.fits' format='var' detmaptype='flat'
    rmf_pn=$pn'v2.rmf'

  #   arf
  #         mos1
  evselect table=$m1'_filt'.fits  imageset=$m1'_img_arf'.fits imagebinning='binSize' \
    withimageset=yes xcolumn='DETX' ycolumn='DETY' ximagebinsize=128 yimagebinsize=128 \
    filtertype='expression' expression="(PATTERN<=12)&&(FLAG==0)"

  arfgen spectrumset=$m1'_spec_v2.fits' arfset=$m1'v2.arf' withrmfset=true \
    rmfset=$m1'v2.rmf' detmaptype='dataset' detmaparray=$m1'_img_arf.fits' \
    extendedsource=yes badpixlocation=$m1'_filt.fits' useodfatt=yes

  arf_m1=$m1'v2.arf'
  rm $m1'_img_arf.fits'

  #        mos2
  evselect table=$m2'_filt'.fits  imageset=$m2'_img_arf'.fits imagebinning='binSize' \
    withimageset=yes xcolumn='DETX' ycolumn='DETY' ximagebinsize=128 yimagebinsize=128 \
    filtertype='expression' expression="(PATTERN<=12)&&(FLAG==0)"

  arfgen spectrumset=$m2'_spec_v2.fits' arfset=$m2'v2.arf' withrmfset=true \
    rmfset=$m2'v2.rmf' detmaptype='dataset' detmaparray=$m2'_img_arf.fits' \
    extendedsource=yes badpixlocation=$m2'_filt.fits' useodfatt=no

  arf_m2=$m2'v2.arf'
  rm $m2'_img_arf.fits'

  #          pn
  evselect table=$pn'_filt'.fits  imageset=$pn'_img_arf'.fits imagebinning='binSize' \
    withimageset=yes xcolumn='DETX' ycolumn='DETY' ximagebinsize=128 yimagebinsize=128 \
    filtertype='expression' expression="(PATTERN==0)&&(FLAG==0)"

  arfgen spectrumset=$pn'_spec_v2.fits' arfset=$pn'v2.arf' withrmfset=true \
    rmfset=$pn'v2.rmf' detmaptype='dataset' detmaparray=$pn'_img_arf.fits' \
    extendedsource=yes badpixlocation=$pn'_filt.fits' useodfatt=yes

  arf_pn=$pn'v2.arf'
  rm $pn'_img_arf.fits'



  spec_blank_m1=$m1'_BLANK_spec_v2.fits'
  spec_blank_m2=$m2'_BLANK_spec_v2.fits'
  spec_blank_pn=$pn'_BLANK_spec_v2.fits'


  #-----GRPPHA---------------------------------
  #COLOCANDO O BLANK, RMF E ARF NO ESPEC DO AGL
  # E binando com minimo de counts = 8
  nbin=25

grppha $m1'_spec_v2.fits' \!$m1'_spec_v2_gr.fits' << LALA
chkey RESPFILE $rmf_m1
chkey ANCRFILE $arf_m1
chkey BACKFILE $spec_blank_m1
group min $nbin
exit
LALA

grppha $m2'_spec_v2.fits' \!$m2'_spec_v2_gr.fits' << LALA
chkey RESPFILE $rmf_m2
chkey ANCRFILE $arf_m2
chkey BACKFILE $spec_blank_m2
group min $nbin
exit
LALA

grppha $pn'_spec_v2.fits' \!$pn'_spec_v2_gr.fits' << LALA
chkey RESPFILE $rmf_pn
chkey ANCRFILE $arf_pn
chkey BACKFILE $spec_blank_pn
group min $nbin
exit
LALA

    if [ $conta -le 2 ]; then

      rm kT.txt -f
      bash /home/natalia/Codigos/spec_med.sh $xs_m1 $xs_m2 $xs_pn $nh_22 $z $E_min
      kT_val=kT2.txt

      kT=$( cat $kT_val | tr -d ' ' | cut -d '=' -f 2 | cut -d '(' -f 1 )
      kTErNeg=$( cat $kT_val | cut -d '-' -f 2 | cut -d ',' -f 1 )
      kTErPos=$( cat $kT_val | cut -d ',' -f 2 | cut -d ')' -f 1 )

      rm R500.txt -f
      python /home/natalia/Codigos/R500.py $kT $kTErNeg $kTErPos $z
      #esses valores saem em Mpc
      aux_r500=R500.txt

      r500=$(cat $aux_r500 | awk -F 'R500_1' '{print $2}' | tr -d ' ' | cut -d '=' -f 2 | cut -d '(' -f 1)
      r500Neg=$(cat $aux_r500 | awk -F 'R500_1' '{print $2}' | tr -d ' ' | cut -d '(' -f 2 | cut -d '-' -f 2 | cut -d ',' -f 1 )
      r500Pos=$(cat $aux_r500 | awk -F 'R500_1' '{print $2}' | tr -d ' ' | cut -d ',' -f 2  | cut -d ')' -f 1)

      r500_015=$(cat $aux_r500 | awk -F 'R500_0.15' '{print $2}' |  tr -d ' ' | cut -d '=' -f 2 | cut -d '(' -f 1)
      r500_015Neg=$(cat $aux_r500 | awk -F 'R500_0.15' '{print $2}' | tr -d ' ' | cut -d '(' -f 2 | cut -d '-' -f 2 | cut -d ',' -f 1 )
      r500_015Pos=$(cat $aux_r500 | awk -F 'R500_0.15' '{print $2}' | tr -d ' ' | cut -d ',' -f 2  | cut -d ')' -f 1)

      conta=$(bc <<< "$conta +1")
    fi

    if [ $conta -le 2 ]; then
      continua=sim

    else
      r_mais=$(echo $r500_015  $r500_015Pos | awk '{print($1+ ($2/2))}')
      r_menos=$(echo $r500_015  $r500_015Neg | awk '{print($1- ($2/2))}')

      #pra comparar, precisa transformar r0 pra Mpc
      if [ "$r0_Mpc" \> "$r_menos" ] && [ "$r0_Mpc" \< $r_mais ]; then
        continua=nao

      else

        rm kT.txt -f
        kT_val=$(bash /home/natalia/Codigos/spec_med.sh $xs_m1 $xs_m2 $xs_pn $nh_22 $z $E_min)
        kT_val=kT.txt

        kT=$( echo $kT_val | tr -d ' ' | cut -d '=' -f 2 | cut -d '(' -f 1 )
        kTErNeg=$( echo $kT_val | cut -d '-' -f 2 | cut -d ',' -f 1 )
        kTErPos=$( echo $kT_val | cut -d ',' -f 2 | cut -d ')' -f 1 )


        rm R500.txt -f
        python /home/natalia/Codigos/R500.py $kT $kTErNeg $kTErPos $z
        aux_r500=R500.txt

        r500=$(cat $aux_r500 | awk -F 'R500_1' '{print $2}' | tr -d ' ' | cut -d '=' -f 2 | cut -d '(' -f 1)
        r500Neg=$(cat $aux_r500 | awk -F 'R500_1' '{print $2}' | tr -d ' ' | cut -d '(' -f 2 | cut -d '-' -f 2 | cut -d ',' -f 1 )
        r500Pos=$(cat $aux_r500 | awk -F 'R500_1' '{print $2}' | tr -d ' ' | cut -d ',' -f 2  | cut -d ')' -f 1)

        r500_015=$(cat $aux_r500 | awk -F 'R500_0.15' '{print $2}' |  tr -d ' ' | cut -d '=' -f 2 | cut -d '(' -f 1)
        r500_015Neg=$(cat $aux_r500 | awk -F 'R500_0.15' '{print $2}' | tr -d ' ' | cut -d '(' -f 2 | cut -d '-' -f 2 | cut -d ',' -f 1 )
        r500_015Pos=$(cat $aux_r500 | awk -F 'R500_0.15' '{print $2}' | tr -d ' ' | cut -d ',' -f 2  | cut -d ')' -f 1)

        continua=sim

      fi
    fi
  r0_Mpc=$(echo $r500_015)

fi

r_cen_agl_Mpc=$(echo $r0_Mpc)

rm R_arcmin.txt -f
ipython << ULA
import pandas as pd
import numpy as np
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

conv=cosmo.angular_diameter_distance($z)
conv=conv.value

div=conv*0.00029088820866572
R_arcmin=$r0_Mpc/div

sai=open('r_arcmin.txt', 'w')

sai.write("r_arcmin=%s" % R_arcmin)
sai.close()
ULA
r_arcmin=$(cat r_arcmin.txt | cut -d '=' -f 2)
#transforma de arcmin pra DET
r0_DET=$(bc  <<<  "scale=8; $r_arcmin*1200.")


r_cen_agl=$(echo $r0_DET)

echo "((DETX,DETY) in circle("$x_cen","$y_cen","$r_ext_agl"))&&!((DETX,DETY) \
in circle("$x_cen","$y_cen","$r_cen_agl"))" > aglomerado.reg

aux_ag=$(<aglomerado.reg)
reg_m1=$(echo $aux_ag $ps_m1)
reg_m2=$(echo $aux_ag $ps_m2)
reg_pn=$(echo $aux_ag $ps_pn)


#cria imagem da total e da regiao a ser analisada, MASCARANDO O CENTRO E AS PS

#TOTAL SEM PS
evselect table=$m1'_filt'.fits  imageset=$m1'_img_regsemps'.fits imagebinning='binSize' withimageset=yes \
xcolumn='DETX' ycolumn='DETY' ximagebinsize=$bin yimagebinsize=$bin \
expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN in [0:4])&& ($reg_m1)"

evselect table=$m2'_filt'.fits  imageset=$m2'_img_regsemps'.fits imagebinning='binSize' withimageset=yes \
xcolumn='DETX' ycolumn='DETY' ximagebinsize=$bin yimagebinsize=$bin \
expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN in [0:4])&&($reg_m2)"

evselect table=$pn'_filt'.fits  imageset=$pn'_img_regsemps'.fits imagebinning='binSize' withimageset=yes \
xcolumn='DETX' ycolumn='DETY' ximagebinsize=$bin yimagebinsize=$bin \
expression="(FLAG==0)&&(PI in [$PImin:$PImax])&&(PATTERN == 0)&&($reg_pn)"

# exptime=$(fkeyprint A2142_mos1_expmap.fits+0 EXPOSURE | grep '=' | cut -d ' ' -f 2)



#========================================================
#ANEIS PARA PERFIL RADIAL
#vamos usar o MOS2 pra tirar os aneis

#calcula os aneis pra analise de perfil
#coords em =>PHY<=

bash /home/natalia/Codigos/transforma_coord.sh $m2'_filt.fits' $m2'_pointsources.txt' PHY 128 tab $m2'_PHY_ps.reg'

echo ' '
echo '======================================================'
echo ' '
echo 'Qual numero minimo de contagens para cada anel?'
read counts_min
echo ' '
echo '======================================================'

bash /home/natalia/Codigos/faz_aneis.sh $m2'_filt.fits' $counts_min 8 $m2'_PHY_ps_excl.reg' centro_PHY.reg 128 1
#saida: aneis_circ.reg

#agora que já tem as regioes, precisa tirar o espectro e o bkg de cada anel. Pro bkg, vamos usar o valor da 
#escala ja calculado

#tem que fazer pra linhas-1 aneis, + 1 regiao, que eh a central

#transforma pra DET coords
bash /home/natalia/Codigos/transforma_coord.sh $m2'_filt.fits' aneis_circ.reg DET 128 coord aneis_circ_DET.reg

arq_aneis='aneis_circ_DET.reg'
linhas=$(wc -l $arq_aneis | awk '{print $1}')
naneis=$(bc <<<"$linhas-1")



#MOS1----------------------


lin=$(echo $linhas)

ExcReg=$(<$m1'_DET_ps_excl.reg')

escala=$(cat $cluster'_BLANK_escalas.txt'| grep EMOS1 | grep escala | awk '{print $4}')

while [ $lin -gt 0 ]; do
  lin_ext=$(echo $lin)
  lin_int=$(bc <<< "$lin_ext -1")
  anel_ext=$(sed "${lin_ext}q;d" $arq_aneis)
  anel_ext="((DETX,DETY) in $anel_ext)"
  anel_int=$(sed "${lin_int}q;d" $arq_aneis)
  anel_int="((DETX,DETY) in $anel_int)"
  anel="$anel_ext&&! $anel_int"
  anel_ps="$anel_ext &&! $anel_int $ExcReg"
  echo "tirando espectro da regiao:   $anel_ext &&! $anel_int"

  #determina qual regiao user
  #pros anies intenos (<5) desconsidera as fontes pontuais
  #pra $lin=1, eh um circulo e nao um anel
  if [ $lin -gt 5 ]; then
    usa_reg=$anel_ps
  elif [[ $lin -le 5 && $lin -gt 1 ]]; then
    usa_reg=$anel
  else
    circ=$(sed "${lin}q;d" $arq_aneis)
    usa_reg="((DETX,DETY) in $circ)"
  fi

  #espec do anel
  spec_anel=$m1'_spec_v2_anel_'$lin'.fits'

  evselect table=$m1'_filt'.fits':EVENTS' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$usa_reg" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$spec_anel spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999

    
  rmf_anel=$m1'_anel_'$lin'.rmf'

  rmfgen rmfset=$rmf_anel withenergybins=no  spectrumset=$spec_anel format='var' detmaptype='flat' 
  
  img=$m1'_img_arf'.fits 
  evselect table=$m1'_filt'.fits  imageset=$img imagebinning='binSize' \
    withimageset=yes xcolumn='DETX' ycolumn='DETY' ximagebinsize=128 yimagebinsize=128 \
    filtertype='expression' expression="($usa_reg)&&(PATTERN<=12)&&(FLAG==0)"

  arf_anel=$m1'_anel_'$lin'.arf'

  arfgen spectrumset=$spec_anel arfset=$arf_anel withrmfset=true \
    rmfset=$rmf_anel detmaptype='dataset' detmaparray=$img \
    extendedsource=yes badpixlocation=$m1'_filt.fits' useodfatt=yes

  
  rm $m1'_img_arf.fits'

  #espec do bkg
  evselect table=$m1'_BLANK_filt_correctpos.fits' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$usa_reg" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$m1'_BLANK_spec_v2_anel_'$lin'.fits' spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999


grppha $m1'_BLANK_spec_v2_anel_'$lin'.fits' \!$m1'_BLANK_spec_v2_anel_'$lin'_gr.fits' << LALA
chkey AREASCAL $escala
exit
LALA

  spec_bkg=$m1'_BLANK_spec_v2_anel_'$lin'_gr.fits'
  spec_bkg=$m1'_BLANK_spec_v2_anel_'$lin'.fits'


grppha $m1'_spec_v2_anel_'$lin'.fits' \!$m1'_spec_v2_anel_'$lin'_gr.fits' << LALA
chkey RESPFILE $rmf_anel
chkey ANCRFILE $arf_anel
chkey BACKFILE $spec_bkg
group min $nbin
exit
LALA

  unset rmf_anel arf_anel spec_bkg anel_ext anel_int lin_ext lin_int usa_reg

  lin=$(bc <<< "$lin -1")

done

unset escala


#MOS2----------------------


lin=$(echo $linhas)

ExcReg=$(<$m2'_DET_ps_excl.reg')

escala=$(cat $cluster'_BLANK_escalas.txt'| grep EMOS2 | grep escala | awk '{print $4}')

while [ $lin -gt 0 ]; do
  lin_ext=$(echo $lin)
  lin_int=$(bc <<< "$lin_ext -1")
  anel_ext=$(sed "${lin_ext}q;d" $arq_aneis)
  anel_ext="((DETX,DETY) in $anel_ext)"
  anel_int=$(sed "${lin_int}q;d" $arq_aneis)
  anel_int="((DETX,DETY) in $anel_int)"
  anel="$anel_ext&&! $anel_int"
  anel_ps="$anel_ext &&! $anel_int $ExcReg"
  echo "tirando espectro da regiao:   $anel_ext &&! $anel_int"

  #determina qual regiao user
  #pros anies intenos (<5) desconsidera as fontes pontuais
  #pra $lin=1, eh um circulo e nao um anel
  if [ $lin -gt 5 ]; then
    usa_reg=$anel_ps
  elif [[ $lin -le 5 && $lin -gt 1 ]]; then
    usa_reg=$anel
  else
    circ=$(sed "${lin}q;d" $arq_aneis)
    usa_reg="((DETX,DETY) in $circ)"
  fi

  #espec do anel
  spec_anel=$m2'_spec_v2_anel_'$lin'.fits'

  evselect table=$m2'_filt'.fits':EVENTS' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$usa_reg" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$spec_anel spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999

    
  rmf_anel=$m2'_anel_'$lin'.rmf'

  rmfgen rmfset=$rmf_anel withenergybins=no  spectrumset=$spec_anel format='var' detmaptype='flat' 


  img=$m2'_img_arf'.fits 
  evselect table=$m2'_filt'.fits  imageset=$img imagebinning='binSize' \
    withimageset=yes xcolumn='DETX' ycolumn='DETY' ximagebinsize=128 yimagebinsize=128 \
    filtertype='expression' expression="($usa_reg)&&(PATTERN<=12)&&(FLAG==0)"

  arf_anel=$m2'_anel_'$lin'.arf'

  arfgen spectrumset=$spec_anel arfset=$arf_anel withrmfset=true \
    rmfset=$rmf_anel detmaptype='dataset' detmaparray=$img \
    extendedsource=yes badpixlocation=$m2'_filt.fits' useodfatt=yes

  
  rm $m2'_img_arf.fits'

  #espec do bkg
  evselect table=$m2'_BLANK_filt_correctpos.fits' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$usa_reg" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$m2'_BLANK_spec_v2_anel_'$lin'.fits' spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999


grppha $m2'_BLANK_spec_v2_anel_'$lin'.fits' \!$m2'_BLANK_spec_v2_anel_'$lin'_gr.fits' << LALA
chkey AREASCAL $escala
exit
LALA

  spec_bkg=$m2'_BLANK_spec_v2_anel_'$lin'_gr.fits'



grppha $m2'_spec_v2_anel_'$lin'.fits' \!$m2'_spec_v2_anel_'$lin'_gr.fits' << LALA
chkey RESPFILE $rmf_anel
chkey ANCRFILE $arf_anel
chkey BACKFILE $spec_bkg
group min $nbin
exit
LALA

  unset rmf_anel arf_anel spec_bkg anel_ext anel_int lin_ext lin_int usa_reg

  lin=$(bc <<< "$lin -1")

done

unset escala

#PN-----------------------

lin=$(echo $linhas)

ExcReg=$(< $pn'_DET_ps_excl.reg')

escala=$(cat $cluster'_BLANK_escalas.txt'| grep EPN | grep escala | awk '{print $4}')

#para os aneis internos, nao mascara as fontes pontuais
while [ $lin -gt 0 ]; do
  lin_ext=$(echo $lin)
  lin_int=$(bc <<< "$lin_ext -1")
  anel_ext=$(sed "${lin_ext}q;d" $arq_aneis)
  anel_ext="((DETX,DETY) in $anel_ext)"
  anel_int=$(sed "${lin_int}q;d" $arq_aneis)
  anel_int="((DETX,DETY) in $anel_int)"
  anel="$anel_ext&&! $anel_int"
  anel_ps="$anel_ext &&! $anel_int $ExcReg"
  echo "tirando espectro da regiao:   $anel_ext &&! $anel_int"

  #determina qual regiao user
  #pros anies intenos (<5) desconsidera as fontes pontuais
  #pra $lin=1, eh um circulo e nao um anel
  if [ $lin -gt 5 ]; then
    usa_reg=$anel_ps
  elif [[ $lin -le 5 && $lin -gt 1 ]]; then
    usa_reg=$anel
  else
    circ=$(sed "${lin}q;d" $arq_aneis)
    usa_reg="((DETX,DETY) in $circ)"
  fi


  #espec do anel
  spec_anel=$pn'_spec_v2_anel_'$lin'.fits'

  evselect table=$pn'_filt'.fits':EVENTS' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$usa_reg" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$spec_anel spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479
    
  rmf_anel=$pn'_anel_'$lin'.rmf'

  rmfgen rmfset=$rmf_anel withenergybins=no  spectrumset=$spec_anel format='var' detmaptype='flat' 
  
  img=$pn'_img_arf'.fits 
  evselect table=$pn'_filt'.fits  imageset=$img imagebinning='binSize' \
    withimageset=yes xcolumn='DETX' ycolumn='DETY' ximagebinsize=128 yimagebinsize=128 \
    filtertype='expression' expression="($usa_reg)&&(PATTERN==0)&&(FLAG==0)"

  arf_anel=$pn'_anel_'$lin'.arf'

  arfgen spectrumset=$spec_anel arfset=$arf_anel withrmfset=true \
    rmfset=$rmf_anel detmaptype='dataset' detmaparray=$img \
    extendedsource=yes badpixlocation=$pn'_filt.fits' useodfatt=yes
  
  rm $pn'_img_arf.fits'

  #espec do bkg
  evselect table=$pn'_BLANK_filt_correctpos.fits' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' \
    expression="$usa_reg" dssblock='' writedss=yes cleandss=no updateexposure=yes \
    filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' \
    withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no withspectrumset=yes \
    spectrumset=$pn'_BLANK_spec_v2_anel_'$lin'.fits' spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479

grppha $pn'_BLANK_spec_v2_anel_'$lin'.fits' \!$pn'_BLANK_spec_v2_anel_'$lin'_gr.fits' << LALA
chkey AREASCAL $escala
exit
LALA

  spec_bkg=$pn'_BLANK_spec_v2_anel_'$lin'_gr.fits'

grppha $pn'_spec_v2_anel_'$lin'.fits' \!$pn'_spec_v2_anel_'$lin'_gr.fits' << LALA
chkey RESPFILE $rmf_anel
chkey ANCRFILE $arf_anel
chkey BACKFILE $spec_bkg
group min $nbin
exit
LALA

  unset rmf_anel arf_anel spec_bkg anel_ext anel_int lin_ext lin_int usa_reg

  lin=$(bc <<< "$lin -1")

done

unset escala

#=== AJUSTE ANEIS ==========

lin=$(echo $linhas)
E_min=0.5

while [ $lin -gt 0 ]; do

  anel_mos1=$m1"_spec_v2_anel_"$lin"_gr.fits"
  anel_mos2=$m2"_spec_v2_anel_"$lin"_gr.fits"
  anel_pn=$pn"_spec_v2_anel_"$lin"_gr.fits"
  out=$cluster'_saidaXspec_anel_'$lin'.txt'
  out_plot=$cluster'_ajusteXspec_anel_'$lin'.ps'


xspec<< AIAI

data 1:1 $anel_mos2 2:2 $anel_mos1 3:3 $anel_pn
ignore bad
ignore 1:1 **-$E_min 8.-** 2:2 **-$E_min 8.-** 3:3 **-$E_min 8.-**
setplot energy
cpd /xs
setplot rebin 4 100
cosmo ,,0.7

model wabs(apec) & &
newpar 1 $nh_22 0.05 
newpar 2 8. -0.1 
newpar 3 0.3 .1 
newpar 4 $z 
newpar 5 0.1 
newpar 10 0.1 
newpar 15 0.1 

fit 2000 0.01
log $out
fit
plot ldata chi
hardcopy $out_plot color

error 1 2 3 4

flux 0.5 8.0

newpar 1 0. -1
lum 0.5 8.0 $z

dummy 0.01 100. 2000 log
lum 0.01 100. $z

resp

label top $cluster Anel $lin
plot ldata chi




exit
y
AIAI

  lin=$(bc <<< "$lin -1")

done

#pegando as informações relevantes pra fazer graficos

out=$cluster'_valores_saidaXspec_aneis.csv'


#HEADER
echo "RAIO (arcmin),nH,nHErPos,nHErNeg,ABUNDANCE,ABErPos,ABErNeg,kT,kTErPo,kTErNeg,CHI2_RED,\
FLUX (ergs/cm^2/s) EMOS1,FLUX (ergs/cm^2/s) EMOS2,FLUX (ergs/cm^2/s) EPN,\
LUM EMOS1 (erg/s),LUM EMOS2 (erg/s),LUM EPN (erg/s),LUM_BOL EMOS1,LUM_BOL EMOS2,LUM_BOL EPN" >>$out

#RAIO => eh bom deixar em RA DEC pra ter mais nocao
# precisa converter
bash /home/natalia/Codigos/transforma_coord.sh $m2'_filt.fits' aneis_circ.reg RA 128 coord aneis_circ_RA.reg

lin=$(echo $linhas)

while [ $lin -gt 0 ]; do


  reg_ext=$(sed "${lin}q;d" aneis_circ_RA.reg)
  ra_ext=$( echo $reg_ext | grep circle | cut -d '(' -f 2 | cut -d ',' -f 3 | cut -d ')' -f 1)

  lin_int=$(bc <<< "$lin -1")

  if [ $lin -gt 1 ]; then
    reg_int=$(sed "${lin_int}q;d" aneis_circ_RA.reg)
    ra_int=$( echo $reg_int | cut -d ',' -f 3 | cut -d ')' -f 1 )
  else
    ra_int=0
  fi

  ra_anel=$(bc <<< "scale=2; ($ra_ext+$ra_int)/2")

  file_in=$cluster'_saidaXspec_anel_'$lin'.txt'

  #nH
  nH_xs=$(cat $file_in | grep 'wabs' | grep 'nH' | head -1 | awk -F '22 ' '{print $2}' | awk -F '+' '{print $1}' |tr -d ' ')
  nHErNeg=$(cat $file_in | grep '(-' | grep '#     1' | cut -d '-' -f 2 | cut -d ',' -f 1)
  nHErPos=$(cat $file_in | grep '(-' | grep '#     1' | cut -d ',' -f 2 | cut -d ')' -f 1)


  #Abund
  Abun=$(cat $file_in | grep 'apec' | grep 'Abundanc' | head -1 | awk -F 'Abundanc' '{print $2}' | awk -F '+' '{print $1}' |tr -d ' ')
  AbunErNeg=$(cat $file_in | grep '(-' | grep '#     3' | cut -d '-' -f 2 | cut -d ',' -f 1)
  AbunErPos=$(cat $file_in | grep '(-' | grep '#     3' | cut -d ',' -f 2 | cut -d ')' -f 1)

  #kT
  kT_xs=$(cat $file_in | grep 'apec' | grep 'kT' | head -1 | awk -F 'keV' '{print $2}' | awk -F '+' '{print $1}' |tr -d ' ')
  kTErNeg=$(cat $file_in | grep '(-' | grep '#     2' | cut -d '-' -f 2 | cut -d ',' -f 1)
  kTErPos=$(cat $file_in | grep '(-' | grep '#     2' | cut -d ',' -f 2 | cut -d ')' -f 1)


  #CHI_RED
  chi_red=$(cat $file_in | grep 'Reduced chi-squared' | head -1 | cut -d '=' -f 2 | awk -F 'for' '{print $1}' | tr -d ' ')

  #FLUX
  flux_m2=$(cat $file_in | grep 'Model Flux' | head -1 | cut -d '(' -f 2 | awk -F 'erg' '{print $1}' | tr -d ' ')
  flux_m1=$(cat $file_in | grep 'Model Flux' | head -2 | tail -1 | cut -d '(' -f 2 | awk -F 'erg' '{print $1}' | tr -d ' ')
  flux_pn=$(cat $file_in | grep 'Model Flux' | tail -1 | cut -d '(' -f 2 | awk -F 'erg' '{print $1}' | tr -d ' ')

  #LUM
  lum_m2=$(cat $file_in | grep 'Model Luminosity' | grep '0.5'| head -1 | cut -d 'y' -f 2 | awk -F 'erg' '{print $1}' | tr -d ' ')
  lum_m1=$(cat $file_in | grep 'Model Luminosity' | grep '0.5'| head -2 | tail -1 | cut -d 'y' -f 2 | awk -F 'erg' '{print $1}' | tr -d ' ')
  lum_pn=$(cat $file_in | grep 'Model Luminosity' | grep '0.5'| tail -1 | cut -d 'y' -f 2 | awk -F 'erg' '{print $1}' | tr -d ' ')

  #LUM_BOL

  lumbol_m2=$(cat $file_in | grep 'Model Luminosity' | grep '0.01'| head -1 | cut -d 'y' -f 2 | awk -F 'erg' '{print $1}' | tr -d ' ')
  lumbol_m1=$(cat $file_in | grep 'Model Luminosity' | grep '0.01'| head -2 | tail -1 | cut -d 'y' -f 2 | awk -F 'erg' '{print $1}' | tr -d ' ')
  lumbol_pn=$(cat $file_in | grep 'Model Luminosity' | grep '0.01'| tail -1 | cut -d 'y' -f 2 | awk -F 'erg' '{print $1}' | tr -d ' ')

  echo $ra_anel','$nH_xs','$nHErPos','$nHErNeg','$Abun','$AbunErPos','$AbunErNeg','\
  $kT_xs','$kTErPos','$kTErNeg','$chi_red','$flux_m1','$flux_m2','$flux_pn','$lum_m1','$lum_m2','\
  $lum_pn','$lumbol_m1','$lumbol_m2','$lumbol_pn >> $out

  lin=$(bc <<< "$lin -1")
done
#adiciona uma coluna com os valores do raio em kpc
#nome arq_saida = $out'_kpc.csv'
python /home/natalia/Codigos/arcminTOkpc.py $out $z
rm $out
out=$cluster'_valores_saidaXspec_aneis_kpc.csv'



#PLOTANDO------------------------

ipython <<AIAI

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

data=pd.read_csv("$out", delimiter=",")

R_min=data['RAIO (arcmin)']
R_pc=data['RAIO (kpc)']  
kT=data['kT']
kTErPo=data['kTErPo']
kTErNeg=data['kTErNeg']

Ab=data['ABUNDANCE']
AbErPo=data['ABErPos']
AbErNeg=data['ABErNeg']


#define eixos
x1=R_min
x2=R_pc
y=kT
yerr=[kTErNeg,kTErPo]


#plota
ax1 = plt.subplot(1,1,1)
ax1.errorbar(x1,y, yerr=yerr, fmt='.', capsize=5, mec='dimgrey', mfc='dimgrey', \
 ms=6, elinewidth=1, ecolor='dimgrey' )

ax2 = ax1.twiny()
ax2.plot(x2,y,marker='', linestyle='')

#labels
ax1.set_ylabel(r'kT [keV]')
ax1.set_xlabel(u'Raio [arcmin]')
ax2.set_xlabel(u'Raio [kpc]')


plt.title("Perfil radial de Temperatura - Abell 2142", \
  fontsize='large', y=1.10)

#limite superior original
sup=plt.xlim()[1]


#ktmedio
#medio
plt.axhline(y=$kT_med, color='indianred', linestyle='--')


errsup=$kT_med+$kTErPos_med
errinf=$kT_med-$kTErNeg_med

plt.fill_between(2*x2-1, errinf, errsup, facecolor='lightcoral', alpha=0.5)


plt.xlim(0.,sup)

#plt.show()
plt.savefig('Radial_kT_profile.jpeg')


y=None
yerr=None
x1=None
x2=None
sup=None
errsup=None
errinf=None
plt.cla()
plt.clf()
#----------

#define eixos
x1=R_min
x2=R_pc
y=Ab
yerr=[AbErNeg,AbErPo]



#plota
ax1 = plt.subplot(1,1,1)
ax1.errorbar(x1,y, yerr=yerr, fmt='.', capsize=5, mec='dimgrey', mfc='dimgrey', \
 ms=6, elinewidth=1, ecolor='dimgrey' )
ax2 = ax1.twiny()
ax2.plot(x2,y,marker='', linestyle='')



#labels
ax1.set_ylabel(r'Metalicidade')
ax1.set_xlabel(u'Raio [arcmin]')
ax2.set_xlabel(u'Raio [kpc]')




#limite superior original
sup=plt.xlim()[1]

plt.title("Perfil radial de Metalicidade - Abell 2142", \
  fontsize='large', y=1.10)

plt.axhline(y=$Abun_med, color='indianred', linestyle='--')

errsup=$Abun_med+$AbunErPos_med
errinf=$Abun_med-$AbunErNeg_med

plt.fill_between(2*x2-1, errinf, errsup, facecolor='lightcoral', alpha=0.5)


plt.xlim(0.,sup)

plt.yscale('log')

plt.savefig('Radial_abund_profile.jpeg')
AIAI


#==============================================
#    PERFIL DE BRILHO SUPERFICIAL
#==============================================

#vamos usar o iraf
#primeiro, precisa de um arquivo que mascara os GAPS e as FONTES


emask expimageset=$m1'_expmap.fits' detmaskset=$m1'_emask.fits' 
emask expimageset=$m2'_expmap.fits' detmaskset=$m2'_emask.fits'
emask expimageset=$pn'_expmap.fits' detmaskset=$pn'_emask.fits'

python /home/natalia/Codigos/troca_emask.py $m1'_emask.fits'
python /home/natalia/Codigos/troca_emask.py $m2'_emask.fits'
python /home/natalia/Codigos/troca_emask.py $pn'_emask.fits'


#acha o brilho superficial
#saidas: $m1/$m2/$pn_counts_aneis.csv
bash /home/natalia/Codigos/brilho.sh $cluster

#adiciona a coluna com os raios em arcmin e kpc
in_ra=$cluster'_valores_saidaXspec_aneis_kpc.csv'
in_cts_m1=$m1'_counts_aneis.csv'
in_cts_m2=$m2'_counts_aneis.csv'
in_cts_pn=$pn'_counts_aneis.csv'

ipython <<AIAI
import pandas as pd
import csv

in_ra='$in_ra'
in_cts_m1='$in_cts_m1'
in_cts_m2='$in_cts_m2'
in_cts_pn='$in_cts_pn'


ra = pd.read_csv(in_ra,delimiter=',')
cts_m1= pd.read_csv(in_cts_m1, delimiter=',')
cts_m2= pd.read_csv(in_cts_m2, delimiter=',')
cts_pn= pd.read_csv(in_cts_pn, delimiter=',')

tam=len(cts_m1)

ra_r_arc=ra['RAIO (arcmin)']
ra_r_kpc=ra['RAIO (kpc)']

cts_r_arc_m1=[]
cts_r_kpc_m1=[]
cts_r_arc_m2=[]
cts_r_kpc_m2=[]
cts_r_arc_pn=[]
cts_r_kpc_pn=[]

for i in range(tam):
  cts_r_arc_m1.append(ra_r_arc[i])
  cts_r_kpc_m1.append(ra_r_kpc[i])
  cts_r_arc_m2.append(ra_r_arc[i])
  cts_r_kpc_m2.append(ra_r_kpc[i])
  cts_r_arc_pn.append(ra_r_arc[i])
  cts_r_kpc_pn.append(ra_r_kpc[i])

cts_m1['raio (arcmin)'] = cts_r_arc_m1
cts_m1['raio (kpc)'] = cts_r_kpc_m1
cts_m2['raio (arcmin)'] = cts_r_arc_m2
cts_m2['raio (kpc)'] = cts_r_kpc_m2
cts_pn['raio (arcmin)'] = cts_r_arc_pn
cts_pn['raio (kpc)'] = cts_r_kpc_pn


nome_arq_m1, ext = in_cts_m1.split('.')
cts_sai_m1=nome_arq_m1+'_kpc.csv'
nome_arq_m2, ext = in_cts_m2.split('.')
cts_sai_m2=nome_arq_m2+'_kpc.csv'
nome_arq_pn, ext = in_cts_pn.split('.')
cts_sai_pn=nome_arq_pn+'_kpc.csv'

cts_m1.to_csv(cts_sai_m1)
cts_m2.to_csv(cts_sai_m2)
cts_pn.to_csv(cts_sai_pn)

AIAI

