#!/bin/csh -f
 
## Localisa as fontes pontuais e cria um .reg ciao/physical com elipses 
## sobre estas regioes.
## Usa wavelets (ciao) mas podemos comparar com  Voronoi

## Gastao 16/04/2018
 
if (${#argv} != 4) then
    echo
    echo "Uso: def_ptsrcs  arquivo_de_eventos  emin(eV)  emax(eV)  bin_imagem"
    echo "Numero de parametros errado... "
    exit
endif

set evt = $1     #  arquivo de eventos
set emin = $2    #  em eV
set emax = $3    #  em eV
set binimg = $4  #  como a imagem e' binada : o resultado e' em funcao desta binagem

set imag = "imagem_orig.fits"
\rm -f $imag pts_vor_list.fits

# extrai a imagem
evselect table=$evt':EVENTS' filtertype='expression' expression="$filtroEnergy" imageset=$imag xcolumn='X' ycolumn='Y' imagebinning='binSize' squarepixels=yes ximagebinsize=$binimg yimagebinsize=$binimg -V 0

### ----------------------------------
###  Usa Voronoi para detectar fontes
### ----------------------------------

punlearn vtpdetect
vtpdetect infile=$imag expfile=none outfile=pts_vor_list.fits scale=1 limit=1e-06 coarse=5 maxiter=25

##  ds9 $imag -regions load pts_vor_list.fits &

### ----------------------------------
###  Usa wavelet para detectar fontes
### ----------------------------------

punlearn wavdetect
wavdetect infile=$imag expfile=none outfile=pts_wav_list.fits scellfile=pts_wav_img.fits imagefile=pts_wav_rec.fits defnbkgfile=pts_wav_bkg.fits scale='2 4 8' psffile="" clobber=yes verbose=0


\rm -f  toto_wav_list2.txt toto_wav_list3.txt pts_wav_list.reg
dmsort pts_wav_list.fits pts_wav_listSort.fits keys=-net_counts copyall=yes clobber=yes
dmlist pts_wav_listSort.fits"[cols POS,NET_COUNTS,NET_RATE,R,ROTANG]" data,clean > toto_wav_list2.txt

set nlin = `wc -l toto_wav_list2.txt | awk '{print $1 - 2}'`
tail -$nlin toto_wav_list2.txt > toto_wav_list3.txt

cat toto_wav_list3.txt | grep -v \# | awk '{print "-ellipse(",$1,",",$2,","$5,",",$6,",",$7")"}' | sed 's/ //g' > pts_wav_list.reg

\rm -f  toto_wav_list2.txt toto_wav_list3.txt pts_wav_list.fits

# visualiza imagem com regioes de fontes, exceto pela fonte principal (o aglomerado, esperamos)

ds9 $imag -regions load pts_wav_list.reg -scale log \
    pts_wav_img.fits -regions load pts_vor_list.fits -cmap b&

exit
