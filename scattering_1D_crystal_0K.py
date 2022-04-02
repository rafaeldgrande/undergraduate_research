#Calculo da amplitude de espalhamento de uma onda em um cristal 1D
#Primeiro caso : posi��es bem definidas
# F(K) = soma de x(exp(-i*K*x))

from math import*
import matplotlib.pyplot as plt
from visual.graph import*
from time import*

print "come�o :", localtime()

n = 1000          #N�mero de �tomos no cristal
a = 1.            #Parametro de rede
h = pi/(1000*a)   #Passo da onda
cristal = []      #Lista com as posi��es de cada �tomo no cristal
eixoX = []        #Eixo x do gr�fico
eixoY = []        #Eixo y do gr�fico

#Contruindo o cristal

for i in range (-n, n+1):
   cristal.append(i)

#"Jogando" a onda no cristal

for K in arange (-10*pi/a, 10*pi/a + h, h) :  #Um certo vetor de onda � escolhido
   F = 0   #Amplitude de espalhamento
   for x in arange (len(cristal)):     #O vetor percorre o cristal
      F += e**(-K*cristal[x]*1j)
   if eixoX.count(K) == 0 :
      eixoX.append(K)
      eixoY.append(abs(F))
   else :
      eixoY[eixoX.index(K)] += abs(F)

#Contantdo o tempo...

print "fim :", localtime()

#Montando o gr�fico

plt.plot(eixoX, eixoY)
plt.xlabel("k'-k")
plt.ylabel("F(k'-k)")
plt.show()

         
      
