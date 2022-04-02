#Calculo da amplitude de espalhamento de uma onda em um cristal 1D
#Segundo caso: posicoes com vibracoes termicas aproximadas por gaussianas
# F(K) = soma de x(exp(-i*K*x))

from numpy import*
import matplotlib.pyplot as plt
from random import*
from time import*

print "começo :", clock()
print localtime()

n = 1000          #Número de átomos no cristal
m = 1000          #Número de cristais 
cristal = []      #Lista com as posições de cada átomo no cristal
eixoX = []        #Eixo x do gráfico
eixoY = []        #Eixo y do gráfico

#Construindo os eixos x e y do gráfico (vetor de onda)

for K in range (-10000, 10001):
   eixoX.append(K*pi/1000)
   eixoY.append(0)

print clock()

#Contruindo o cristal

for k in range (m):
   for i in range (1, n+1):
      cristal.append(gauss(i, 0.25))
         
###"Jogando" a onda no cristal
##
##for K in eixoX :  #Um certo vetor de onda é escolhido
##   F = 0   #Amplitude de espalhamento
##   for x in cristal:     #O vetor percorre o cristal
##      F += e**(-K*x*1j)
##   eixoY[eixoX.index(K)] = (abs(F)/m)

###Construindo o cristal (2) e jogando a onda!!
##
##for k in range (m):
##   cristal = []
##   for i in range (1, n+1):
##      cristal.append(gauss(i, 0.25))
##   for K in eixoX :
##      F = 0
##      for x in cristal :
##         F += e**(K*x*1j)
##      eixoY[eixoX.index(K)] += abs(F)/m
   
#Contantdo o tempo...

##print "fim :", clock()/60, 'minutos ou ', clock()/(3600),'horas'
##print localtime()
##
###Montando o gráfico
##
##plt.plot(eixoX, eixoY)
##plt.xlabel("k'-k")
##plt.ylabel("F(k'-k)")
##plt.show()

plt.hist(cristal, bins = 5000)
plt.show()
