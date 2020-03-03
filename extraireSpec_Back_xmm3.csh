qui#!/bin/csh -f
#
# Script pour l'extraction d'un spectre et background dans une region
# Associa uma RMF e uma ARF a cada espectro
# XMM - PN e MOS
#
# 19/04/2007 acrescenta opcao "useodfatt=no" no arfgen se necessario
# 20/04/2007 saida opcional no xspec (xspecops = y ou t ou s)
# 24/04/2007 Versao para blank field de A. Read
# 13/01/2009 Versao 3, apropriada para ser rodada pelo pipeline
# 16/01/2009 O rmfgen nao roda bem para o PN no MAC: a grid de energia tem que ser colocada na mao
# 24/02/2011 Passa para xspec12
# 11/05/2011 Leva em consideracao o CCD6 do MOS 1 morto depois de 9/marco/2005
# 01/06/2011 O pipeline ja' produz o arquivo corrigido do 
#            MOS1 depois do micrometeoroide
# 08/09/2011 deixa energia minima variavel para ajuste (normalizacao)
# 20/02/2013 chama xspec12, um wrapper de xspec para rodar junto com SAS

if ( ${#argv} != 6 ) then
    echo
    echo "Uso: extraireSpec_Back_xmm   EventosClean   EventosCleanBack .regExtracao"
    echo "    .regSrcPts    binspec   xspec?"
    echo "regioes formato ds9/acis, uma por linha com + ou -"
    echo "se nao ha exclusao, regSrcPts = null"
    echo "binspec binagem do grppha."
    echo "xspec? = y ou s ou t ==> liga xspec automaticamente"
    echo
   exit
endif

## PARA O HPC64 apenas
## setenv TMPDIR '/home/tmp'

if ($?SAS_ODF) then
  echo "Este e'  script extraireSpec_Back_xmm3"
  echo 'Estamos usando SAS_ODF:' $SAS_ODF
  echo 
else
  echo 'A variavel SAS_ODF nao esta definida. To fora...'
  exit
endif

set comp = `uname`

set eminPN = 5.0
set eminMOS = 5.0

# Parametres d'entree:
set fich_in = $1      # Fichier d'events "cleaned" des flares, pattern, etc...
set fich_back_in = $2 # Arquivo de eventos Clean e alinhado  do Background
set fich_reg = $3     # Fichier .reg en coordonnees DET (phys no ds9)
set excReg = $4       # Fichier .reg des sources exclues
set nbin = $5         # numero minimo de counts/canal de energia (grppha)
set xspecops = $6     # y ou s ou t ==> usa 'useodfatt=no' no arfgen"

## fich_reg e' a regiao onde extraimos o espectro, p.ex., 
## ellipseMaior + rotbox - ellipseMenor
## O "+" significa um "&&" ("and" logico)
## excReg contem (em geral) as fontes puntuais
## O arquivo .reg deve ter UMA REGIAO POR LINHA com sinal de + ou -.

## Determina se e' MOS ou PN

set inst = `dmkeypar $fich_in INSTRUME echo+`

set obs = "P"`dmkeypar $fich_in OBS_ID echo+`$inst
## set obs = `echo $fich_in:r | awk '{print substr($1,7,13)}'`
set reg = `echo $fich_reg:r:t`


## Determina o numero da orbita e set a orbita do impacto do 
#  micrometeoroid no MOS1
## set orbita = `dmkeypar $fich_in REVOLUT echo+`
## set impact = "961"


# Fichiers de sortie:
set <<="spec_"$obs$reg".fits"            # Spectre non-binne
set fich_spec_bin="spec_bin_"$obs$reg".fits"    # Spectre binne par grppha
set fich_rmf=$obs$reg".rmf"                     # RMF
set fich_ext_img="img_"$obs$reg".fits"          # Image de la region
set fich_arf=$obs$reg".arf"                     # ARF
set fich_spec_back="spec_back_"$obs$reg".fits"  # Spectre Background (meme region DET)
set fich_spec_back_bin="spec_back_bin_"$obs$reg".fits" # Spectre Back binado
set fich_spec_back_bin_esca = "spec_back_bin_esca"$obs$reg".fits" # Spectre Back escalonado
set fich_spec_back_soft = "spec_soft"$obs$reg".fits"   # espectro de fundo soft


### O background (blankfield ja' deve estar alinhado !!

set fich_back =  $fich_back_in
set fich_back_alin = $fich_back


## Extraire la region avec la bonne syntaxe
## set regiao=`cat $fich_reg |grep circle | sed '1,$s/image\;//'`
##set regiao=`cat $fich_reg |grep -v "#" `

#set re = "((X,Y) IN "`cat $fich_reg | grep -v "#"`
#set regiao = `echo $re | sed 's/-/\) \&\& \(\!\(X,Y\) IN /g'`")"
#unset re

###  A regiao "principal". Se houver "rotbox", transforma-lo em sintaxe XMM
###   rotbox  --> box(xcen,ycen,xwidth/2,width/2,rotacao)

\rm -f totoTMP1.reg totoTMP.reg

set xwid = `cat $fich_reg | grep rotbox | sed 's/,/ /g'| awk '{print $3}'`

if ( $xwid == "" ) then

  cat $fich_reg > totoTMP1.reg

else

  set yhei = `cat $fich_reg | grep rotbox | sed 's/,/ /g'| awk '{print $4}'`

  set xwid2 = `echo $xwid | awk '{print $1/2}'`
  set ywid2 = `echo $yhei | awk '{print $1/2}'`

  cat $fich_reg | sed 's/'$xwid'/'$xwid2'/' | sed 's/'$yhei'/'$ywid2'/' > totoTMP1.reg

endif

## Concatena com fontes puntuais `a excluir

if ($excReg == null || $excReg == NULL || $excReg == Null || $excReg == nul || $excReg == NUL || $excReg == Nul) then

  cat totoTMP1.reg  | grep -v "#" | sed 's/rotbox/box/g' > totoTMP.reg

  set regiao = `echo "(X,Y) in ";cat totoTMP.reg | sed 's/+/\&\& (X,Y) in /g' | sed 's/-/\&\& \!(X,Y) in /g'`

else

  cat totoTMP1.reg $excReg | grep -v "#" | sed 's/rotbox/box/g' > totoTMP.reg

  set regiao = `echo "(X,Y) in ";cat totoTMP.reg | sed 's/+/\&\& (X,Y) in /g' | sed 's/-/\&\& \!(X,Y) in /g'`

endif


## Extraire le spectre dans la region "regiao"


if ($inst == EMOS1 || $inst == EMOS2) then
evselect table=$fich_in':EVENTS' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' expression="$regiao" dssblock='' writedss=yes cleandss=no updateexposure=yes filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no  withspectrumset=yes spectrumset=$fich_spec spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999

else

evselect table=$fich_in':EVENTS' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' expression="$regiao" dssblock='' writedss=yes cleandss=no updateexposure=yes filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no  withspectrumset=yes spectrumset=$fich_spec spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479

endif

# A linha de baixo funciona se o xmgrace estiver instalado
# dsplot table=$fich_spec &


if (-e $fich_rmf) then
    echo "RMF ja existe"
else
    echo "Cree le RMF"

    if ($inst == EPN && $comp == Darwin) then

    ####### o rmfgen nao funciona alem de 13.5 keV para o PN no Mac...Gastao 16/01/2009

        rmfgen rmfset=$fich_rmf withenergybins=yes energymin=0.05 energymax=13.55 nenergybins=1000 spectrumset=$fich_spec format='var' detmaptype='flat'

    else
        rmfgen rmfset=$fich_rmf withenergybins=no  spectrumset=$fich_spec format='var' detmaptype='flat'
    endif

endif

if (-e $fich_ext_img) then
    echo "Imagem ja existe"
else
    echo "Extraire image de la region"

    evselect table=$fich_in':EVENTS' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' expression="$regiao" dssblock='' writedss=yes cleandss=no updateexposure=yes filterexposure=yes blockstocopy='' attributestocopy=''  ignorelegallimits=no withimageset=yes imageset=$fich_ext_img xcolumn='DETX' ycolumn='DETY' imagebinning='binSize' squarepixels=yes ximagebinsize=64 yimagebinsize=64 

endif

if (-e $fich_arf) then
    echo "ARF ja existe"
else
    echo "Cree ARF"

    arfgen spectrumset=$fich_spec withrmfset=yes rmfset=$fich_rmf arfset=$fich_arf detmaptype='dataset' detmaparray=$fich_ext_img withdetbounds=no detxoffset=1200 detxbins=5 detyoffset=1200 detybins=5 psfenergy=5 filterdss=yes extendedsource=yes modeleffarea=yes modelquantumeff=yes modelfiltertrans=yes modelee=yes modelootcorr=yes eegridfactor=100 withbadpixcorr=yes badpixlocation=$fich_in setbackscale=no keeparfset=yes useodfatt=yes

    if (! -e $fich_arf) then
 
       echo "Opa, nao deu, vou tentar com 'useodfatt=no'"

       arfgen spectrumset=$fich_spec withrmfset=yes rmfset=$fich_rmf arfset=$fich_arf detmaptype='dataset' detmaparray=$fich_ext_img withdetbounds=no detxoffset=1200 detxbins=5 detyoffset=1200 detybins=5 psfenergy=5 filterdss=yes extendedsource=yes modeleffarea=yes modelquantumeff=yes modelfiltertrans=yes modelee=yes modelootcorr=yes eegridfactor=100 withbadpixcorr=yes badpixlocation=$fich_in setbackscale=no keeparfset=yes useodfatt=no

   endif

endif

#==================   BACKGROUND  =================================

echo "Extraire le spectre BACK dans la region "$regiao

if ($inst == EMOS1 || $inst == EMOS2) then

## if ($inst == EMOS1) then
##    if ($orbita > $impact) then
##    echo "CCD 6 do MOS1 kaput"
##
##    set regiao = `echo $regiao "&& CCDNR != 6"`
##    
##    endif
## endif

evselect table=$fich_back_alin':EVENTS' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' expression="$regiao" dssblock='' writedss=yes cleandss=no updateexposure=yes filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no  withspectrumset=yes spectrumset=$fich_spec_back spectralbinsize=15 withspecranges=yes specchannelmin=0 specchannelmax=11999

else

## setenv SAS_MEMORY_MODEL low

evselect table=$fich_back_alin':EVENTS' destruct=yes flagcolumn='EVFLAG' flagbit=-1 filtertype='expression' expression="$regiao" dssblock='' writedss=yes cleandss=no updateexposure=yes filterexposure=yes blockstocopy='' attributestocopy='' energycolumn='PI' withzcolumn=no zcolumn='WEIGHT' withzerrorcolumn=yes zerrorcolumn='EWEIGHT' ignorelegallimits=no  withspectrumset=yes spectrumset=$fich_spec_back spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479

## setenv SAS_MEMORY_MODEL high

endif

## unsetenv TMPDIR

echo "Binning avec  GRPPHA"

#===================    GRPPHA   para data  =================================

grppha $fich_spec \!$fich_spec_bin << LALA
chkey RESPFILE $fich_rmf
chkey ANCRFILE $fich_arf
exit
LALA

#===================    GRPPHA   para Back  =================================

grppha $fich_spec_back \!$fich_spec_back_bin << LALA
chkey RESPFILE $fich_rmf
chkey ANCRFILE $fich_arf
exit
LALA

## Determina escalonamento do fundo

if ($inst == EPN) then
##  set ig1 = "**-1.3 8.0-**"  <== original
  set ig1 = "**-5. 11.-**"
  set ig2 = "7.0-9.5"

else
##  set ig1 = "**-2.0 8.0-**"  <== original
  set ig1 = "**-5.0 8.0-**"
  set ig2 = "1.0-2.0"

endif

if (! -e $fich_rmf) then
   echo
   echo "Erro: nao tem arquivo rmf"
   echo
   exit
else if (! -e $fich_arf) then
   echo
   echo "Erro: nao tem arquivo arf"
   echo
   exit
endif

xspec12<<EOF
data 1:1 $fich_spec_bin 2:2 $fich_spec_back_bin
query no
ignore 1:1 $ig1 2:2 $ig1
ignore 1:1 $ig2 2:2 $ig2
ignore bad
log xspecTEMPO.txt
show rates
exit
y
EOF

set ctsfonte = `cat xspecTEMPO.txt | grep Net | head -1 | awk '{print $7}'`
set ctsback = `cat xspecTEMPO.txt | grep Net | tail -1 | awk '{print $7}'`

set escala = `echo $ctsback $ctsfonte | awk '{print $1/$2}'`

set expFonte = `fkeyprint $fich_spec_bin+1 EXPOSURE | grep -v \# | grep EXPOSURE | awk '{print $2}'`
set expBlank = `fkeyprint $fich_spec_back_bin+1 EXPOSURE | grep -v \# | grep EXPOSURE | awk '{print $2}'`

set ratioT = `echo $expBlank $expFonte | awk '{print $1/$2}'`

echo
echo 'Razao de contagens (blank/fonte):' $escala
echo 'Razao de exposure  (blank/fonte):' $ratioT
echo

echo
echo "Espectro da obs., anel externo:   " $fich_spec_bin
echo "Espectro Background Escalonado:   " $fich_spec_back_bin
echo "Espectro Back Escalonado e binado:" $fich_spec_back_bin_esca
echo "Espectro de residuo (soft excess):" $fich_spec_back_soft
echo

#===================    GRPPHA   para Back ESCALONADO  =================================

grppha $fich_spec_back_bin \!$fich_spec_back_bin_esca << LALA
chkey AREASCAL $escala
exit
LALA


#======== Calcula espectro soft-excess


cp $fich_spec_bin specSoft.fits
fcalc specSoft.fits specSoft2.fits FOCOUNTS "COUNTS"
fdelcol specSoft2.fits+1 COUNTS N Y
faddcol specSoft2.fits $fich_spec_back_bin COUNTS
fcalc specSoft2.fits specSoft3.fits FOCOUNTS "FOCOUNTS-(COUNTS/$ratioT)/$escala" clobber=yes
fdelcol specSoft3.fits+1 COUNTS N Y
fcalc specSoft3.fits specSoft4.fits COUNTS "FOCOUNTS"

mv specSoft4.fits $fich_spec_back_soft

\rm specSoft.fits specSoft2.fits specSoft3.fits xspecTEMPO.txt totoTMP1.reg totoTMP.reg


# ==================   XSPEC      ==================================
if ($xspecops == t || $xspecops == s || $xspecops == y) then

set eps = $inst$reg".eps/cps"

xspec12 << TOTO
data 1:1 $fich_spec_bin 2:2 $fich_spec_back_bin_esca 3:3 $fich_spec_back_soft
cpd $eps
setplot energy
ignore 1:1 **-0.15 13.-** 2:2 **-0.15 13.-** 3:3 **-0.15 13.-**
ignore bad
setplot rebin 5 30
plot data
exit
y
TOTO

endif

unset xspecops eps fich_spec_bin fich_spec_back_bin fich_spec_back_bin_esca fich_spec_back_soft
## unset orbita impact

\rm OFFSET*.ds

exit

#
#
#===================    GRPPHA   para data  =================================

grppha $fich_spec \!$fich_spec_bin << LALA
chkey RESPFILE $fich_rmf
chkey ANCRFILE $fich_arf
group min $nbin
exit
LALA

#===================    GRPPHA   para Back  =================================

grppha $fich_spec_back \!$fich_spec_back_bin << LALA
chkey RESPFILE $fich_rmf
chkey ANCRFILE $fich_arf
group min $nbin
exit
LALA
########
