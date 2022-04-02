#Calculo da amplitude de espalhamento de uma onda em um cristal 1D
#Quarto caso : quasi-cristal (cadeia de Fibonacci)
# F(K) = soma de x(exp(-i*K*x))

from numpy import*
import matplotlib.pyplot as plt
from random import*
from time import*

print "começo :", clock()
print localtime()

n = 1000          #Número de átomos no cristal
m = 1000          #Número de cristais
cristal2 = []
cristal = []      #Lista com as posições de cada átomo no cristal
eixoX = []        #Eixo x do gráfico
eixoY = []        #Eixo y do gráfico
p1 = "A"
p2 = "AB"          #Passos na cadeia de Fibbonacci
p = 0

#Contruindo o cristal

while len(p2) <= 2000 :
   p1 = p2 + p1
   p2 = p1 + p2

for ch in p2 :
   if ch == "A" :
     p += 1
   elif ch == "B":
      p += 2
   cristal2.append(p)

for i in range (m):
   for j in range (n) :
      cristal.append(cristal2[j+i] - cristal2[i])

#Construindo o eixo x do gráfico (vetor de onda)

for K in range (-10000, 10001):
   eixoX.append(K*pi/1000)

print clock()
         
#"Jogando" a onda no cristal

for K in eixoX :  #Um certo vetor de onda é escolhido
   F = 0   #Amplitude de espalhamento
   for x in cristal:     #O vetor percorre o cristal
      F += e**(-K*x*1j)
   eixoY.append(abs(F)/m)
   
#Contantdo o tempo...

print "fim :", clock()/60, 'minutos ou ', clock()/(3600),'horas'
print localtime()

#Montando o gráfico

plt.plot(eixoX, eixoY)
plt.xlabel("k'-k")
plt.ylabel("F(k'-k)")
plt.show()


         
      
