#!/bin/csh -f

## Cria uma serie de aneis concentricos com uma contagem liquida (fonte - fundo) dada
## Usa o script  def_centro_regs.csh
## O fundo e' estimado na propria observacao, em um anel na borda

## Gastao 16/04/2018

if (${#argv} != 5) then
    echo
    echo "Uso: def_centro_regs.csh  arquivo_eventos.fits  counts_minimos   r_max(arcmin)  bin_imagem  quieto"
    echo "quieto = 0 => verbose minima"
    echo "quieto = 1 => verbose tagarela"
    echo
    exit
endif

set bla = $0
if ($bla == "tcsh" || $bla == "csh") then      ## usou source
    set raiz = `history 1 |awk '{print $4}'`
    set raiz = `eval "echo $raiz:h"`
else                                           ## e' um executavel  
    set raiz = $bla:h
endif
unset bla


set T0 = `date +%s` ## tempo em segundos no incio do script


set evt = $1  ## nome do arquivo de eventos
set netcountLim = `echo $2 | awk '{print int($1)}'`  ## numero de contagens net desejado para cada anel
set rLimMaxArcMin = $3        ## raio externo maximo em arcmin
set binimg = $4               ## binagem da imagem  128 ou 64 parecem OK
set quieto = $5

set emin = 500   #  em eV
set emax = 7000  #  em eV
set rbkg2 = 12.0   ## arcmin raio ext fundo. Provavelmente OK
set rbkg1 = 10.5   ## arcmin raio int fundo. Idem

set rStep = 160    ## pixel physical. Se for menor OK, mas demora mais

## se o valor final de contagens for até 80% do limite, OK
set limTolera = `echo $netcountLim | awk '{print int($1 * 0.8)}'`

set regiaoAneis = "aneis_circ.reg"

set filtroEnergy="(PI in ("$emin":"$emax"))"

set pix =  `fkeyprint $evt"[0]" REFYCDLT | grep WCS |awk '{print $2}'`
set pix = `echo $pix | awk '{print $1 * 3600.}'`

## raio maximo em pixel physical inteiro (nao real)
set rLimMax = `echo $rLimMaxArcMin $pix |  awk '{print int($1 * 60. / $2)}'`

# define regiao para estimar o bkgrnd
set xcenExp = `fkeyprint $evt"[0]" REFXCRPX | grep WCS | awk '{print $2}'`
set ycenExp = `fkeyprint $evt"[0]" REFYCRPX | grep WCS | awk '{print $2}'`

set rbkg1phy = `echo $rbkg1 $pix | awk '{print $1 * 60. / $2}'`
set rbkg2phy = `echo $rbkg2 $pix | awk '{print $1 * 60. / $2}'`

set backReg = "((X,Y) IN circle($xcenExp,$ycenExp,$rbkg2phy) && (!(X,Y) IN circle($xcenExp,$ycenExp,$rbkg1phy)))"

set toto = bkgIMG_tmp.fits

\rm -f $toto
# extrai a imagem
evselect table=$evt':EVENTS' filtertype='expression' expression="$backReg && $filtroEnergy" imageset=$toto imagedatatype=Real32 xcolumn='X' ycolumn='Y' imagebinning='binSize' squarepixels=yes ximagebinsize=$binimg yimagebinsize=$binimg -V 0

## Usa sigma clip com 3.5 sigma para tentar garantir convergencia nao nula
punlearn dmstat
dmstat $toto  centroid=no clip=yes nsigma=3.5 > /dev/null

set med = `pget dmstat out_mean`
set sig = `pget dmstat out_sigma`

if ($med == 0) then  ## tenta de novo com sigma maior
  dmstat $toto  centroid=no clip=yes nsigma=5.5 > /dev/null

  set med = `pget dmstat out_mean`
  set sig = `pget dmstat out_sigma`

endif

echo
echo "Fundo = " $med "+/-" $sig "cnt/pixel (imagem binada, bin= "$binimg"). Rmax(phys)= "$rLimMax

set T1 =  `date +%s`
echo $T1 $T0 | awk '{print $1-$2,"segs"}'


## Acha o centro do obj de interesse -----------------------------------------------------------

# O parametro 'quito' tem que ser igual a zero aqui
set bla = `source $raiz/def_centro_regs.csh $evt $emin $emax $binimg 0`

set cenx = `echo $bla | awk '{print $3}'`
set ceny = `echo $bla | awk '{print $7}'`

echo "Centro (phys): " $cenx" , "$ceny

## ----------------------------------------------------------------------------------------------

set T1 =  `date +%s`
echo $T1 $T0 | awk '{print $1-$2,"segs"}'


## Acha as fontes que serao mascaradas  ---------------------------------------------------------
## point sources e alguma subestrutura

source $raiz/def_ptsrcs.csh $evt $emin $emax $binimg

\rm -f totoTMP.reg
cat  pts_wav_list.reg | grep -v "#" > totoTMP.reg
set  regExc = `cat totoTMP.reg | sed 's/-/\&\& !(X,Y) IN /g'`

## ----------------------------------------------------------------------------------------------

set T1 =  `date +%s`
echo $T1 $T0 | awk '{print $1-$2,"segs"}'


## Raios iniciais dos aneis : o primeiro e' um circulo (r interno = 0)
set rINTphy = 0.0
set rOUTphy = 300.0
set rOUTphyINT = 300

\rm -f $regiaoAneis
touch $regiaoAneis  ##  cria arquivo com os circulos (aneis)

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - LOOP EXTERNO
while ($rOUTphyINT <= $rLimMax)

set netcoutINT = 0

#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . LOOP INTERNO
while ($netcoutINT <= $netcountLim && $rOUTphyINT <= $rLimMax)

# aumenta r_externo e tenta
set rOUTphy = `echo $rOUTphy $rStep | awk '{print $1 + $2}'`
set rOUTphyINT = `echo $rOUTphy | awk '{print int($1)}'`

set region =  "((X,Y) IN circle($cenx,$ceny,$rOUTphy) && (!(X,Y) IN circle($cenx,$ceny,$rINTphy)))"
set region = `echo $region $regExc`


set totoA = anel_tmp.fits   ## arquivo temporario com a imagem do anel

\rm -rf $totoA
# extrai imagem
evselect table=$evt':EVENTS' filtertype='expression' expression="$region && $filtroEnergy" imageset=$totoA imagedatatype=Real32 xcolumn='X' ycolumn='Y' imagebinning='binSize' squarepixels=yes ximagebinsize=$binimg yimagebinsize=$binimg -V 0

punlearn dmstat
dmstat $totoA centroid=no  > /dev/null
set medAnel = `pget dmstat out_mean`
set npix = `pget dmstat out_good`

set netcount = `echo $medAnel $med $npix |awk '{print ($1 - $2) * $3}'`
set netcoutINT = `echo $netcount | awk '{print int($1)}'`   ## arredonda p/ valor inteiro (p/ o loop)

if ($quieto == 1)  echo "Total counts: " $netcount "entre "$rINTphy" e "$rOUTphy" (phys)"

end   ### while : ate' chegar no numero de counts necessario
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if ($netcoutINT >= $limTolera) then     ####    OK, temos um novo circulo

 echo " Convergiu total counts: " $netcount "entre "$rINTphy" e "$rOUTphy" (phys)"

 echo "circle("$cenx,$ceny,$rOUTphy")" >> $regiaoAneis

 set rINTphy = $rOUTphy

 set rOUTphyINT = `echo $rOUTphy | awk '{print int($1)}'`

endif

end    ## while
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

echo
set T1 =  `date +%s`
echo $T1 $T0 | awk '{print $1-$2,"segs"}'


ds9 $evt -bin factor $binimg -scale sqrt -regions -format ciao $regiaoAneis \
    -regions -format ciao pts_wav_list.reg &

\rm -f $toto totoTMP.reg

unset T0, T1, rINTphy, rOUTphyINT, netcount, netcoutINT, npix, medAnel, totoA, toto

exit
