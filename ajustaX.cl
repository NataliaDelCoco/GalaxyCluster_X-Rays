#
#
# pra rodar, no iraf: cl < ajustaX.cl
# talvez tenha que abrir um termnal e um ds9 antes
#--------------------------------------------------
#Quando estiver no modo interativo do ellipse/iraf, você pode usar as seguintes 
#teclas:

#x -> para mudar a posição do centro
#f -> para refazer o 1o ajuste

#h -> para fazer todos os ajustes

#Quando acabar, abra os arquivos *model.fits e *resid.fits
#-----------------------------------------------------------------


stsdas
analys
isoph
unlearn ellipse
unlearn magpar
unlearn geompar
unlearn samplepar
unlearn imarith
unlearn display
unlearn isoplot

string nome, fich, fich2, fich3, fich4, fich5, tab, tabasc, model, resid

## ------------------------------
## Coloque apenas o prefixo do nome como no exemplo abaixo
## Para rodar, dentro do iraf, fazer cl < ajustaX.cl (nome deste arquivo)

cd /home/natalia/Dados/A1795/0097820101/analysis

nome = "A1795_mos2_img_semps"

## ------------------------------


fich  = nome//".fits"
fich2 = nome//"_100x.fits"
fich3 = "blo.fits"
fich4 = "blu.txt"
tab   = nome//"_acis.tab"
tabasc= nome//"_acis.txt"
model = nome//"_Model.fits"
resid = nome//"_Resid.fits"

del(tab); del(tabasc); del(model); del(resid)
del(fich2); del(fich3); del(fich4)


imarith(operand1=fich, op="/",  operand2="10.", result=fich2)

##  number of sigma-clip iterations
samplepar.nclip=12


## Parametros da geometria EM IMG

geompar.x0= 155
geompar.y0= 156

geompar.ellip0 = 0.1
geompar.sma0 = 15

## -----  nao ajustamos em menos de 1 pixels  ##
geompar.minsma = 1

## -----  este valor pode ser maior ou menor; depende da imagem  ##
geompar.maxsma=200

## Muda parametros default do display para podermos ver o centro da galaxia
display.zscale=no
display.ztrans="log"

## =====>>>>  Roda o stsdas/ellipse
ellipse(input=fich2, output=tab, interac=yes, masksz=1, region=yes)

## ====>>>>  Constroi um modelo baseado no ellipse
bmodel(table=tab,output=fich3)

imarith(operand1=fich3, op="*", operand2="10.", result=model)
imarith(operand1=fich,  op="-", operand2=model,  result=resid)

## ====>>>> Faz o grafico mag X R^(1/4) (perfil de de Vaucouleurs)
isoplot(input=tab, xaxis="RSMA", yaxis="MAG")

### ------  Saida de uma tabela ascii (.txt) para ajuste em outro programa
string colunas
colunas="SMA,RSMA,INTENS,RMS,MAG,MAG_LERR,MAG_UERR"

ttools.tdump(table=tab, datafil=fich4, columns=colunas)

string linha
linha = "!echo " // nome //  " > toto.txt"
print(linha)|cl()

linha = "!echo '   ' " // colunas // " >> toto.txt"
print(linha)|cl()

linha = "!cat " // fich4 // " >> toto.txt"
print(linha)|cl()

linha = "!cat toto.txt | sed 's/,/ /g' > " // tabasc
print(linha)|cl()

### Apaga arquivos temporarios

!\rm -f toto.txt

del(fich3) ; del(fich4)

print("Terminou. Veja os arquivos "//nome//" Model e Resid.fits")

