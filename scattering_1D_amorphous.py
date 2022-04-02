#Calculo da amplitude de espalhamento de uma onda em um cristal 1D
#Terceiro caso : cristal amorfo (Xn = X(n-1) + deltaX
# F(K) = soma de x(exp(-i*K*x))

from numpy import*
import matplotlib.pyplot as plt
from random import*
from time import*

print "come�o :", clock()

n = 1000          #N�mero de �tomos no cristal
m = 1000          #N�mero de cristais 
cristal = []      #Lista com as posi��es de cada �tomo no cristal
eixoX = []        #Eixo x do gr�fico
eixoY = []        #Eixo y do gr�fico

#Contruindo o cristal

for k in range (m):
   x = 0
   for i in range (1, n+1):
      x += (gauss(1, 0.01))
      cristal.append(x)

#Construindo o eixo x do gr�fico (vetor de onda)

for K in range (-10000, 10001):
   eixoX.append(K*pi/1000)

print clock()
         
#"Jogando" a onda no cristal

for K in eixoX :  #Um certo vetor de onda � escolhido
   F = 0   #Amplitude de espalhamento
   for x in cristal:     #O vetor percorre o cristal
      F += e**(-K*x*1j)
   eixoY.append(abs(F)/m)
   
#Contantdo o tempo...

print "fim :", clock()/60, 'minutos ou ', clock()/(3600),'horas'

#Montando o gr�fico

plt.plot(eixoX, eixoY)
plt.xlabel("k'-k")
plt.ylabel("F(k'-k)")
plt.show()

         
      
