from math import*
import matplotlib.pyplot as plt
distancias = []
contador = []
camada = 0
alfa = 0
eixoX = []
eixoY = []
eixoY2 = []
T = 0    #Total de atomos 
n = 40  #"Dimensões do cristal"
for i in range (-n, n+1):
    for j in range (-n, n+1):
        for k in range (-n, n+1):
            d = sqrt(i**2 + j**2 + k**2)
            p = distancias.count(d)
            if p == 0 :
               distancias.append(d)
               contador.append(d)
            else :
               q = distancias.index(d)
               contador[q] += 1
for D in range (1, n**2 + 1) :
   if distancias.count(D**(0.5)) != 0 :
      camada += 1
      p = distancias.index(D**(0.5))
      if camada % 2 == 1 :
         alfa += contador[p]/(distancias[p])
      else :
         alfa -= contador[p]/(distancias[p])
      T += contador[p]   #Total de átomos
      eixoX.append(camada)
      eixoY.append(alfa)
plt.plot(eixoX, eixoY)
plt.xlabel('Camada')
plt.ylabel('Alfa')
plt.show()

a = input('Aperte enter para finalizar')
