#Programa de dinâmica molecular (2D)

import numpy as np
import matplotlib.pyplot as plt
from time import*

print "START :", localtime()
clock()     #Começa a contar o tempo da simulação

#------------------------------------------------------------------------------
#------------------Constantes(em unidades reduzidas) e listas------------------
#------------------------------------------------------------------------------

sigma = 1                # = 3.4e-10 metros
epsilon = 1              # = 1.65e-21 Joules
m = 1                    # = 6.69e-26 Kg
Lx = 6                  #Comprimento e largura da caixa em unidades de sigma
Ly = 6
k = 1                    #Constante de Boltzman
R = 0                    #Módulo da distância entre duas particulas
r = [0,0]                #Vetor unitário que aponta de uma particula até outra
Particulas = []          #Lista contendo as paticulas do sistema
t = 0                    #Tempo inicial 
dt = 0.01                #Intervalo dt (unidade = 2.17e-12 = sigma*(sqrt(m/epsilon))
limite = 1              #Limite de tempo da simulação
U = []                   #Energia potencial total do sistema
K = []                   #Energia cinética total do sistema
E = []                   #Energia total do sistema
P = []                   #Pressão do sistema
Temp = []                #Temperatura do sistema
DeltaEdet = []           #Variação da energia total na simulação
TrajetoriasX = []        #Trajetórias das partículas no eixo x
TrajetoriasY = []        #Trajetórias das particulas no eixo y
Nt = []                  #Número de particulas de um lado da célula
Particulas_tempo_t = []  #Posições e velocidades no final da simulação t = 0.5
T = []                   #Lista com os valores do tempo
Dt = []                  #Lista contendo valores diferentes para o passo
Px = []                  #Momento total na direção x
Py = []                  #Momento total na direção y
ErroP = 0.001            #Erro máximo aceitável do momento
  
#------------------------------------------------------------------------------
#--------------Definindo a classe particula e seus atributos-------------------
#------------------------------------------------------------------------------

class Particula :
   rx = float
   ry = float            #Posições x e y
   vx = float
   vy = float            #Velocidades nas direções x e y
   ax = float
   ay = float            #Acelerações mas direções x e y
   m = float             #Massa da particula

   def ModuloDistancia (self):
      return np.sqrt(self.rx**2 + self.ry**2)

   def ModuloVelocidade (self):
      return np.sqrt(self.vx**2 + self.vy**2)

   def Momento (self):
      return [self.vx*self.m, self.vy*self.m]

   def EnergiaCinetica (self):
      return (self.m/2.)*self.ModuloVelocidade()**2


#--------------------------------------------------------------------------
#------Esta parte do prgrama lê o arquivo 'Estado inicial.txt'-------------
#------que contem os valores a serem atribuidos a cada particula-----------
#--------------------------------------------------------------------------
   
def ConverteStringLista (String) :
   x = list()
   P = str()
   for ch in String :
      if not(ch == '[' or ch == ',' or ch == ' ' or ch == ']'):
         P += ch
      elif ch == ',' or ch == ']' :
         x.append(float(P))
         P = str()
   return x

def Ler_arquivo (Arquivo_in) :
   arq = open(str(Arquivo_in)+'.txt')

   for i in arq :
      if i[0] != '#' and i != '\n' :
         x = ConverteStringLista(i)
         A = Particula()
         A.rx, A.ry, A.vx, A.vy, A.m = x[0], x[1], x[2], x[3], x[4]
         Particulas.append(A)

   for A in Particulas :
      TrajetoriasX.append([A.rx])
      TrajetoriasY.append([A.ry])

   arq.close()

#-----------------------------------------------------------------------------
#----------------Escrevendo um arquivo txt com as posições e velocidades------
#---------------------------------no tempo t----------------------------------
#-----------------------------------------------------------------------------

def Estado_no_tempo_t (Particulas_tempo_t, Arquivo_out) :
   arq = open(Arquivo_out+'.txt', 'w')

   arq.write('#Dados com das posições e velocidades no tempo t\n'
          '#para o programa de dinâmica molecular\n'
          '#Os dados estão organizados em strings da seguinte maneira :\n'
          '#[rx, ry, vx, vy, m], onde : \n'
          "# 'r' é o vetor posição, 'v' e velocidade e 'm' é a massa.\n"
          '\n')

   for i in range(len(Particulas_tempo_t)) :
      arq.write(str(Particulas_tempo_t[i])+'\n')

   arq.close()

#-----------------------------------------------------------------------------
#-------------Funções que determinam características do sistema---------------
#-----------------------------------------------------------------------------

def Densidade (Particulas, Lx, Ly) :
   return float(len(Particulas))/(Lx*Ly)
         
def EnergiaCineticaTotal (Particulas):
   K = 0
   for A in Particulas:
      K += A.EnergiaCinetica()
   return K

def CentroDeMassa (Particulas):
   CM = [0,0]
   M = 0.
   for A in Particulas :
      CM[0], CM[1] = CM[0]+A.rx*A.m, CM[1]+A.ry*A.m
      M += A.m
   CM[0], CM[1] = CM[0]/M, CM[1]/M
   return CM
   
def MomentoTotal (Particulas):
   P = [0,0]
   for A in Particulas :
      P[0] += A.Momento()[0]
      P[1] += A.Momento()[1]
   return P

def PotencialLJ (A, B) :
   R = DeltaR(A, B)[0]
   return 4*epsilon*((sigma/R)**12 - (sigma/R)**6)

def EnergiaPotSist (Particulas) :  #Calcula a energia potencial entre as particulas 
   U = 0                      #do sistema
   for i in range (len(Particulas)) :
      for j in range (i+1, len(Particulas)) :
         U += PotencialLJ(Particulas[i], Particulas[j])
   return U

def Temperatura (Particulas) :
   return EnergiaCineticaTotal(Particulas)/(len(Particulas))  

def Pressao (Particulas) :
   P = len(Particulas)*Temperatura(Particulas)/(Lx*Ly)
   for i in range (len(Particulas)) :
      for j in range (i+1, len(Particulas)) :
         P += (np.dot(Forca(Particulas[i], Particulas[j]),
                     DeltaR(Particulas[i], Particulas[j])[1])*
               DeltaR(Particulas[i], Particulas[j])[0]/(2*Lx*Ly))
   return P

#-----------------------------------------------------------------------------
#--------------------------Medidas estatísticas-------------------------------
#-----------------------------------------------------------------------------

def Media(Lista) :
   return float(sum(Lista))/len(Lista)

def Varianca (Lista) :
   v = 0
   m = Media(Lista)
   for i in Lista :
      v += pow(i - m, 2)
   return v/len(Lista)

def DesvioPadrao (Lista) :
   return np.sqrt(Varianca(Lista))

#-----------------------------------------------------------------------------
#-------------------------------Ajuste linear---------------------------------
#-----------------------------------------------------------------------------

def AjusteLinear (x, y) :
   A = np.vstack([x, np.ones(len(x))]).T

   return np.linalg.lstsq(A, y)[0]    #Primeiro item da lista é o coeficiente
                                      #angular e o segundo o linear

#-----------------------------------------------------------------------------
#-------------Funções necessárias ao programa---------------------------------
#-----------------------------------------------------------------------------

def DeltaR (A, B) : #Mostra a menor distância no sistema com condições de contorno
   x1, x2, y1, y2 = A.rx, B.rx, A.ry, B.ry
   dx = x2 - x1
   dy = y2 - y1
   if abs(dx) > Lx/2. :
      if dx < 0 :
         dx += Lx
      else :
         dx -= Lx
   if abs(dy) > Lx/2. :
      if dy < 0 :
         dy += Ly
      else :
         dy -= Ly
   R = np.sqrt(dx**2 + dy**2)
   if R == 0 :
      r = [0,0]
   else :
      r = [dx/R, dy/R]
   return [R, r]

def Forca (A, B) :   #Devolve o módulo da força de A em B
   R = DeltaR(A, B)[0]
   if R <= 2.5 :
      f = (24*epsilon/R)*(2*pow(sigma/R,12) - pow(sigma/R,6))
      fx, fy = DeltaR(A,B)[1][0]*f, DeltaR(A,B)[1][1]*f 
      return [fx, fy]
   else :
      return [0, 0]

def CondicoesContorno (A, Lx, Ly) :         #Não deixa que as particulas saiam
   if A.rx > Lx :                        # do sistema (caixa)
      A.rx -= Lx
   elif A.rx < 0 :
      A.rx += Lx
   if A.ry > Ly :
      A.ry -= Ly
   elif A.ry < 0 :
      A.ry += Ly

def Aceleracao (Particulas, sigma, epsilon) :   #Calcula a aceleração sobre todas as
   for A in Particulas :
      A.ax, A.ay = 0, 0
   for i in range (len(Particulas)) :
      for j in range (i+1, len(Particulas)) :
         Particulas[i].ax += Forca(Particulas[j], Particulas[i])[0]
         Particulas[i].ay += Forca(Particulas[j], Particulas[i])[1]
         Particulas[j].ax -= Particulas[i].ax
         Particulas[j].ay -= Particulas[i].ay


#-----------------------------------------------------------------------------
#---------------------------Algoritmo de Verlet-------------------------------
#-----X(n+1) = Xn + Vn*dt + (An)*dt**2/2   e V(n+1) = (An + An+1)*dt/2--------
#-----------------------------------------------------------------------------


def Verlet (Particulas,t, dt, sigma, epsilon, Lx, Ly, limite, Particulas_tempo_t, ErroP) :
   p = True
   while (t <= limite and p == True) :
      for A in Particulas :
         A.rx += A.vx*dt + A.ax*dt**2/2.
         A.ry += A.vy*dt + A.ay*dt**2/2.
         A.vx += A.ax*dt/2
         A.vy += A.ay*dt/2
      Aceleracao (Particulas, sigma, epsilon)    #Atualiza aceleração
      for A in Particulas :
         A.vx += A.ax*dt/2
         A.vy += A.ay*dt/2
         CondicoesContorno(A, Lx, Ly)     
         if (0 <= A.rx <= Lx and 0 <= A.ry <= Ly) == False :
            print "Erro! Alguma particula está com velocidade muito alta!"
            p  = False     #Verificando se há algum erro grave
            
#----------------------------Algoritmo de Verlet acaba aqui.----------------------

      U.append(EnergiaPotSist(Particulas))   #Energias potencial e cinética
      K.append(EnergiaCineticaTotal(Particulas))
      E.append(EnergiaPotSist(Particulas)
               +EnergiaCineticaTotal(Particulas))
      P.append(Pressao(Particulas))
      Temp.append(Temperatura(Particulas))

      n = 0                                     #Número de particulas de um lado da caixa      
      for A in Particulas :                     #Plotando as trajetórias das particulas
         TrajetoriasX[Particulas.index(A)].append(A.rx)
         TrajetoriasY[Particulas.index(A)].append(A.ry)
         if A.rx >= Lx/2. :      #Plotando o número de particulas de um lado 
            n += 1               # da caixa
      Nt.append(n)

      Px.append(MomentoTotal(Particulas)[0]) #Verificando variação do momento total
      Py.append(MomentoTotal(Particulas)[1])
      
      if (t/dt) % 1000 == 0 :                #Resetar o momento total
         for A in Particulas :
            A.vx -= MomentoTotal(Particulas)[0]/(m*len(Particulas))
            A.vy -= MomentoTotal(Particulas)[1]/(m*len(Particulas))
                            
      T.append(t)
      t += dt
      
   for A in Particulas :
      Particulas_tempo_t.append([A.rx, A.ry, -A.vx, -A.vy, m])  #O sinal de '-' é para
                                                          #ver se o sistema volta para
                                                          #o estado inicial

   
   
#-----------------------------------------------------------------------------
#-----------------------Executando o programa em si --------------------------
#-----------------------------------------------------------------------------        

Arquivo_in = 'Estado inicial'
Arquivo_out = 'Estado no tempo t'

Ler_arquivo (Arquivo_in)

Aceleracao (Particulas, sigma, epsilon)    #Atualiza aceleração
Verlet (Particulas, t, dt, sigma, epsilon, Lx, Ly, limite, Particulas_tempo_t, ErroP)   

Estado_no_tempo_t (Particulas_tempo_t, Arquivo_out)


print "FIM :", localtime()
print clock(), "segundos =", clock()/60.,"minutos =", clock()/3600., "horas"

#-----------------------------------------------------------------------------
#---------------------------Plotando gráficos---------------------------------
#-----------------------------------------------------------------------------

plt.figure(1)

plt.plot(T, Temp)
plt.xlabel("t")
plt.ylabel("T(t)")
plt.grid(True)
plt.text(limite*0.75, 1.43, '$\mu$ ='+str(Media(Temp)))
plt.text(limite*0.75, 1.40, '$\sigma$ ='+str(DesvioPadrao(Temp)))




#----------------------Gráfico com ajuste linear-----------------------------

def Energia_Ajuste_Linear (T, E) :

   plt.figure(2)

   a, b = AjusteLinear(T, E)    #Fazendo o ajuste linear

   E2 = []
   Ruido = 0
   for i in range (len(T)) :
      E2.append(T[i]*a + b)
      Ruido += (E2[i] - E[i])**2

   Ruido = np.sqrt(Ruido)/len(E)
      

   plt.plot(T, E, label = "E(t)")
   plt.plot(T, E2, label = "Ajuste linear")
   plt.text(limite*0.75,Media(E),'Ruido ='+str(Ruido))
   plt.xlabel("t")
   plt.ylabel("E(t)")

#-----------------------------------------------------------------------------

Energia_Ajuste_Linear (T, E)
plt.show()


