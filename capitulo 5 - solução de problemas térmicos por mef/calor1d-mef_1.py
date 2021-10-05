

import numpy as np
import matplotlib.pyplot as plt
import mesh as ms

# parametros da simulacao
L = 1.0
ne = 4 # numero de elementos da malha quadratica (igual linear)
mode_N=3 # Ordem de N
mode_X=1 # distribuicao dos pontos da malha (linear(1), quadratico(2), gauss("gauss"))
Te = 100.0
Td = 0.0
q = 180.0
rho = 1.0
cv = 1.0
kappa = 1.0

nx=(ne+1)+(ne*mode_N-ne) #numero de pontos da malha no eixo x
alpha = kappa/(rho*cv)

# geracao dos pontos e IEN

msh=ms.mesh_1d(nx,L,mode_X) #usando gerador de malha
X=msh[1]
IEN=np.array(msh[0]) #IEN para quando N é linear ja e dada pelo gerador
if mode_N!=1:
    IEN=[] # se N nao for linear zera a malha
    for i in range(ne):
        if mode_N==2:
            IEN.append([i*mode_N,i*mode_N+1,i*mode_N+2]) # [0,1,2],[2,3,4]...
        if mode_N==3:
            IEN.append([i*mode_N,i*mode_N+1,i*mode_N+2,i*mode_N+3]) # [0,1,2,3]...   
IEN=np.array(IEN)  #transformando lista em array 

# plot da condicao de contorno de T
T = np.zeros( (nx),dtype='float' )
T[0]=Te
T[-1]=Td #condições de temperatura na ponta
plt.plot(X,T,'r-') #plota condição de temperatura inicial

# inicializacao das matrizes e vetores
K = np.zeros( (nx,nx),dtype='float' )
M = np.zeros( (nx,nx),dtype='float' )
f = np.zeros( (nx),dtype='float' )
# loop nos elementos
#-------------------------------------------------- 
for e in range(0,ne):
 if mode_N==1: #N linear
  [v1,v2] = IEN[e] 
  # calculo do comprimento do elemento
  h = X[v2]-X[v1]
  # matrizes e vetores do elemento mestre
  kelem = (1.0/h)*np.array([ [1,-1],
                             [-1,1] ])
  felem = (h/2.0)*np.array( [1, 1] )
 if mode_N==2: # N quadrado
  [v1,v2,v3] = IEN[e] 
  # calculo do comprimento do elemento quadratico
  h = X[v3]-X[v1]
  # matrizes e vetores do elemento mestre
  kelem = (1.0/(3*h))*np.array([ [7,-8,1],
                                 [-8,16,-8], 
                                 [1,-8,7] ]) 
  felem = (h/6.0)*np.array( [1,4,1] )
 if mode_N==3: # N cubico
  [v1,v2,v3,v4] = IEN[e] 
  # calculo do comprimento do elemento quadratico
  h = X[v4]-X[v1]
  # matrizes e vetores do elemento mestre
  kelem = (1.0/(40*h))*np.array([ [148,-189,54,-13],
                                  [-189,432,-297,54],
                                  [54,-297,432,-189],
                                  [-13,54,-189,148]]) 
  felem = (h/8.0)*np.array( [1,3,3,1] )

 # loops nos pontos do elemento e
 for i in range(0,mode_N+1):
  iglobal = IEN[e,i]
  f[iglobal] += q/(rho*cv)*felem[i]
  for j in range(0,mode_N+1):
   jglobal = IEN[e,j]

   K[iglobal,jglobal] += alpha*kelem[i,j]


A = K.copy()
# formulacao classica (com vetor f)
b = f.copy()


# condicao de contorno de Dirichlet
A[0,:] = 0.0
A[0,0] = 1.0
b[0] = Te
A[-1,:] = 0.0
A[-1,-1] = 1.0
b[-1] = Td

# solucao do sistema linear

T = np.linalg.solve(A,b)

#Plot da temperatura
plt.plot(X,T,'ko-')
plt.xlabel('comprimento da barra [m]')
plt.ylabel('temperatura [oC]')
plt.title("Calor 1d mef com geração de calor") # adicionando titulo ao plot
plt.show()

