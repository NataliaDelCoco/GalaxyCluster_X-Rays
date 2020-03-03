#!/bin/bash

## Cria uma serie de aneis concentricos com uma contagem liquida (fonte - fundo) dada
## Usa o script  def_centro_regs.csh
## O fundo e' estimado na propria observacao, em um anel na borda

## Natalia 18/08/18




if [ $# -ne 8 ]
then
  echo
  echo "uso: bash faz_aneis.sh   event_file  counts_minimos  r_max(arcmin) fontes_pontuais  pos_centro\
  bin_imagem  quiet nome_output"
  echo "quieto = 0 => verbose minima"
  echo "quieto = 1 => verbose tagarela"
  echo

  return
fi


evt=$1  ## nome do arquivo de eventos
netcountLim=$(echo $2 | awk '{print int($1)}')  ## numero de contagens net desejado para cada anel
rLimMaxArcMin=$3        ## raio externo maximo em arcmin
aux_fontes=$4
centro=$5
binimg=$6               ## binagem da imagem  128 ou 64 parecem OK
quieto=$7
saida=$8


inst=$( fkeyprint $evt+0 INSTRUME | grep '=' | cut -d ' ' -f 2 | cut -c 2- )

if [ "$inst" = "EMOS1" ]
then
  nome='mos1'
elif [ "$inst" = "EMOS2" ]
then
  nome='mos2'
else
  nome='pn'
fi

regExc=$(<$aux_fontes)

emin=500   #  em eV
emax=8000  #  em eV
rbkg2=14.0   ## arcmin raio ext fundo. Provavelmente OK
rbkg1=12.   ## arcmin raio int fundo. Idem

rStep=160    ## pixel physical. Se for menor OK, mas demora mais

## se o valor final de contagens for até 80% do limite, OK
limTolera=$(echo $netcountLim | awk '{print int($1 * 0.8)}')


filtroEnergy="(PI in ("$emin":"$emax"))"

pix=$(fkeyprint $evt"[0]" REFYCDLT | grep WCS |awk '{print $2}')
pix=$(echo $pix | awk '{print $1 * 3600.}')

## raio maximo em pixel physical inteiro (nao real)
rLimMax=$(echo $rLimMaxArcMin $pix |  awk '{print int($1 * 60. / $2)}')

# define regiao para estimar o bkgrnd
xcenExp=$(fkeyprint $evt"[0]" REFXCRPX | grep WCS | awk '{print $2}')
ycenExp=$(fkeyprint $evt"[0]" REFYCRPX | grep WCS | awk '{print $2}')

rbkg1phy=$(echo $rbkg1 $pix | awk '{print $1 * 60. / $2}')
rbkg2phy=$(echo $rbkg2 $pix | awk '{print $1 * 60. / $2}')

backReg="((X,Y) IN circle($xcenExp,$ycenExp,$rbkg2phy)) &&! ((X,Y) IN circle($xcenExp,$ycenExp,$rbkg1phy))"

toto=bkgIMG_tmp.fits

\rm -f $toto
# extrai a imagem
evselect table=$evt':EVENTS' filtertype='expression' expression="$backReg && $filtroEnergy" imageset=$toto imagedatatype=Real32 xcolumn='X' ycolumn='Y' imagebinning='binSize' squarepixels=yes ximagebinsize=$binimg yimagebinsize=$binimg -V 0

## Usa sigma clip com 3.5 sigma para tentar garantir convergencia nao nula
# punlearn dmstat
# dmstat $toto  centroid=no clip=yes nsigma=3.5 > /dev/null

# set med = `pget dmstat out_mean`
# set sig = `pget dmstat out_sigma`

eregionanalyse imageset=$toto srcexp="$backReg" > saida_img_bkg.txt

med=$( cat saida_img_bkg.txt | grep "cnts per pixel" | cut -d ' ' -f 7 )  #num de counts por pixel

# if ($med == 0) then  ## tenta de novo com sigma maior
#   dmstat $toto  centroid=no clip=yes nsigma=5.5 > /dev/null

#   set med = `pget dmstat out_mean`
#   set sig = `pget dmstat out_sigma`

# endif

# echo
# echo "Fundo = " $med "+/-" $sig "cnt/pixel (imagem binada, bin= "$binimg"). Rmax(phys)= "$rLimMax



## Acha o centro do obj de interesse -----------------------------------------------------------

cenx=$(cat $centro | grep circle | cut -d '(' -f 3 | cut -d ',' -f 1)
ceny=$(cat $centro | grep circle | cut -d '(' -f 3 | cut -d ',' -f 2)
echo "Centro (phys): " $cenx" , "$ceny
## ----------------------------------------------------------------------------------------------


## Acha as fontes que serao mascaradas  ---------------------------------------------------------
## point sources e alguma subestrutura

# source $raiz/def_ptsrcs.csh $evt $emin $emax $binimg

# \rm -f totoTMP.reg
# cat  pts_wav_list.reg | grep -v "#" > totoTMP.reg
# set  regExc = `cat totoTMP.reg | sed 's/-/\&\& !(X,Y) IN /g'`

## ----------------------------------------------------------------------------------------------



## Raios iniciais dos aneis : o primeiro e' um circulo (r interno = 0)
rINTphy=0.0
rOUTphy=300.0
rOUTphyINT=300

\rm -f $saida
touch $saida  ##  cria arquivo com os circulos (aneis)

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - LOOP EXTERNO
while [ $rOUTphyINT -le $rLimMax ];do

  netcoutINT=0

#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . LOOP INTERNO
  while [ $netcoutINT -le $netcountLim -a $rOUTphyINT -le $rLimMax ];do

    # aumenta r_externo e tenta
    rOUTphy=$(echo $rOUTphy $rStep | awk '{print $1 + $2}')
    rOUTphyINT=$(echo $rOUTphy | awk '{print int($1)}')

    region="((X,Y) in circle($cenx,$ceny,$rOUTphy)) &&! ((X,Y) in circle($cenx,$ceny,$rINTphy))"
    region=$(echo $region $regExc)

    #echo "regiao=" $region


    totoA=anel_tmp.fits   ## arquivo temporario com a imagem do anel

    \rm -rf $totoA
    # extrai imagem
    evselect table=$evt':EVENTS' filtertype='expression' expression="$region && $filtroEnergy" imageset=$totoA imagedatatype=Real32 xcolumn='X' ycolumn='Y' imagebinning='binSize' squarepixels=yes ximagebinsize=$binimg yimagebinsize=$binimg -V 0

    # punlearn dmstat
    # dmstat $totoA centroid=no  > /dev/null
    # set medAnel = `pget dmstat out_mean`
    # set npix = `pget dmstat out_good`

    eregionanalyse imageset=$totoA srcexp="$region" > saida_img_anel.txt 

    medAnel=$( cat saida_img_anel.txt | grep "cnts per pixel" | cut -d ' ' -f 7 ) #counts per pixel
    ncount=$( cat saida_img_anel.txt | grep "source region" | cut -d ' ' -f 6 ) #num total counts
    npix=$(echo $ncount $medAnel | awk '{print int($1/$2)}')

    netcount=$(echo $medAnel $med $npix |awk '{print ($1 - $2) * $3}')
    netcoutINT=$(echo $netcount | awk '{print int($1)}')   ## arredonda p/ valor inteiro (p/ o loop)

    if [ $quieto == 1 ]
    then
      echo "Total counts: " $netcount "entre "$rINTphy" e "$rOUTphy" (phys)"
    fi

    done   ### while : ate' chegar no numero de counts necessario
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if [ $netcoutINT -ge $limTolera ]
     then     ####    OK, temos um novo circulo

     echo " Convergiu total counts: " $netcount "entre "$rINTphy" e "$rOUTphy" (phys)"

     echo "circle("$cenx,$ceny,$rOUTphy")" >> $saida

     rINTphy=$rOUTphy

     rOUTphyINT=$(echo $rOUTphy | awk '{print int($1)}')

    fi

done    ## while
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




# ds9 $evt -bin factor $binimg -scale sqrt -regions -format ciao $regiaoAneis \
#     -regions -format ciao pts_wav_list.reg &

\rm -f $toto totoTMP.reg

 unset  rINTphy rOUTphyINT netcount netcoutINT npix medAnel totoA toto

rm saida_img_bkg.txt anel_tmp.fits saida_img_anel.txt

# exit