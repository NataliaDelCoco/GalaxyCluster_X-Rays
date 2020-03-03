lista_cls=/home/natalia/Dados/XMM_catalogues/cls.txt
cls_file=/home/natalia/Dados/XMM_catalogues/clusters_usados.csv
salva_file=/home/natalia/Dados/Aglomerados/Xdata_total.csv

rm -f $salva_file

echo "Cluster, obs, RA(X-ray), DEC(X-ray), redshift, kTc, kTc+, kTc-,\
 kTm, kTm+, kTm-, Zc, Zc+, Zc-, Zm, Zm+, Zm-,\
 n0fit, n0fitErr, rc, rcErr, alpha, alphaErr, beta, betaErr, bkg, bkgErr,\
 chi_red, Cuspiness, Cuspiness_err, CSB, CSB_err, CSB4, CSB4_err, kT_ratio, kT_ratioErr,\
 R500(arcmin), R500(kpc), n0 (count/arcmin³), n0Err, Mgas_tot (1e13Msun),\
 MgasErr, Mtot_isoT(1e13Msun), MtotErr, Mmaug, Yx" >> $salva_file

while IFS= read -r line
do
  ag=$(echo "$line" )
  echo $ag

  path=/home/natalia/Dados/Aglomerados/$ag/*
  cd $path
  obs=$(pwd |rev | cut -d '/' -f 1 |rev)
  cd analysis

  # INFOS GERAIS
  RA=$(cat $cls_file | grep $ag | cut -d ',' -f 2 | xargs)
  DEC=$(cat $cls_file | grep $ag | cut -d ',' -f 3 | xargs)
  z=$(cat $cls_file | grep $ag | cut -d ',' -f 4 | xargs)
  #conv=$(cat $cls_file | grep $ag | rev | cut -d ',' -f 1 | rev | xargs)


  # Tlines=$(wc -l < avrgT.txt)
  # T=$(head -$Tlines avrgT.txt | tail -1 | awk '{print $2}')
  # r500_arc=$(head -$(($Tlines-2)) avrgT.txt | tail -1 | awk '{print $5}')
  # r500_Mpc=$(python /home/natalia/Codigos/conv_raio.py $conv $r500_arc 'arcmin2Mpc')

  # AJUSTE ESPECTRAL
  tc=$(head -56 fitavrgfinal.log | tail -1 | awk '{print $7}')     #Fetches the fitted temperature and metalicity from the bulk (R500-0.15 R500) and core (0.15 R500)
  metalc=$(head -57 fitavrgfinal.log | tail -1 | awk '{print $6}')

  errorTcl=$(head -1 errorout-avrg.dat | tail -1 | awk '{print $1}')
  errorTch=$(head -1 errorout-avrg.dat | tail -1 | awk '{print $2}')
  errorMcl=$(head -2 errorout-avrg.dat | tail -1 | awk '{print $1}')
  errorMch=$(head -2 errorout-avrg.dat | tail -1 | awk '{print $2}')

  t=$(head -206 fitavrgfinal.log | tail -1 | awk '{print $7}')     
  metal=$(head -207 fitavrgfinal.log | tail -1 | awk '{print $6}')

  errorTl=$(head -3 errorout-avrg.dat | tail -1 | awk '{print $1}')
  errorTh=$(head -3 errorout-avrg.dat | tail -1 | awk '{print $2}')
  errorMl=$(head -4 errorout-avrg.dat | tail -1 | awk '{print $1}')
  errorMh=$(head -4 errorout-avrg.dat | tail -1 | awk '{print $2}')

  # #BRILHO SUPERFICIAL
  BS=$(cat BrilhoSup_fit_params.txt | tail -1)

  # #MASSAS
  MASS=$(cat Massas_tot_fit.txt |tail -1)


  echo "$ag, $obs, $RA, $DEC, $z, $tc, $errorTch, $errorTcl,\
  $t, $errorTh, $errorTl, $metalc, $errorMch, $errorMcl, $metal, $errorMh, $errorMl, $BS $MASS" >> $salva_file 


  # echo "$ag, $obs, $RA, $DEC, $z, $tc, $errorTch, $errorTcl,\
  # $t, $errorTh, $errorTl, $metalc, $errorMch, $errorMcl, $metal, $errorMh, $errorMl\
  # $BS, $MASS" >> $salva_file 


done < "$lista_cls"
