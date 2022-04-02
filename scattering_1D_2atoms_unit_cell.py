#Calculo da amplitude de espalhamento de uma onda em um cristal 1D
#Quarto caso : liga bin�ria de dois tipos de �tomos A e B distribuidos
#aleatoriamente. Cada �tomo contribui para a densidade com uma fun��o
#delta multiplicada por um fator de forma. Fa = 1 e Fb = 2
# F(K) = soma de x(exp(-j*K*x))

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
Fa = 1            #Fator de forma do �tomo A
Fb = 2            #Fator de forma do �tomo B
Pa = 0.25          #Probabilidade do �tomo A
Pb = 1 - Pa       #Probabilidade do �tomo B

#Contruindo o cristal

for k in range (m):
   for i in range (1, n+1):
      p = random()
      if p <= Pa :
         cristal.append("A")
      else :
         cristal.append("B")

#Construindo o eixo x do gr�fico (vetor de onda)

for K in range (-10000, 10001):
   eixoX.append(K*pi/1000)

print clock()
         
#"Jogando" a onda no cristal

for K in eixoX :  #Um certo vetor de onda � escolhido
   F = 0   #Amplitude de espalhamento
   for x in range(len(cristal)) :     #O vetor percorre o cristal
      if cristal[x] == "A" :
         F += Fa*e**(-K*(x%n)*1j)
      else :
         F += Fb*e**(-K*(x%n)*1j)
   eixoY.append(abs(F)/m)
   
#Contantdo o tempo...

print "fim :", clock()/60, 'minutos ou ', clock()/(3600),'horas'

#Montando o gr�fico

plt.plot(eixoX, eixoY)
plt.xlabel("k'-k")
plt.ylabel("F(k'-k)")
plt.show()

         
      
