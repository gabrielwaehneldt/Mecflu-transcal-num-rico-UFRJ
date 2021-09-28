

import numpy as np
import matplotlib.pyplot as plt

# parametros da simulacao
L = 1.0
npoints = 30
ne = npoints-1
dx = L/ne
Q = 100.0 # fonte de calor
k = 1.0  # condutividade termica do material
cv = 1.0 # capacidade termica
rho = 1.0 # densidade
dt = 0.05 # passo de tempo

time = 0.0 # tempo de simulacao

# condicao de contorno de Dirichlet (espacial)
Te = 100.0
Td = 0.0

# condicao de inicial (temporal)
T = np.zeros( (npoints),dtype='float' )
T[0] = Te
T[-1] = Td

# geracao dos pontos
X = np.linspace(0,L,npoints)

# geracao da matriz de conectividade
IEN = np.zeros( (ne,2),dtype='int' )
for e in range(0,ne):
 IEN[e] = [e,e+1]

# vetor de indices de contorno
cc = [0,npoints-1]
# vetor dos valores do contorno
bcc = np.zeros( (npoints),dtype='float' )
bcc[0] = Te
bcc[-1] = Td

# plot da condicao inicial (time=0.00)
plt.plot(X,T,'r-')

# Eq. do ponto i (pontos do miolo) para time=0.01

# inicializar a matriz A e o vetor b
A = np.zeros( (npoints,npoints),dtype='float' )
b = np.zeros( (npoints),dtype='float' )
# populando os valores dos pontos de contorno
for i in cc:
  A[i,i] = 1.0
  b[i]   = bcc[i]
for i in range(1,npoints-1):
  A[i,i]   =  1+ 2*( dt/(dx*dx) ) * (k/(rho*cv)) # diagonal principal
  A[i,i-1] =  -( dt/(dx*dx) ) * (k/(rho*cv)) # diagonal inferior
  A[i,i+1] =  -( dt/(dx*dx) ) * (k/(rho*cv)) # diagonal superior

for n in range(20):

 # populando os valores dos pontos internos de b
 for i in range(1,npoints-1):
  b[i]     =(( dt/(rho*cv) ) * Q)+T[i] #so precisa do T[i] dentro do for
 T = np.linalg.solve(A,b)

plt.plot(X,T,'ko-')
plt.xlabel('comprimento da barra [m]')
plt.ylabel('temperatura [oC]')
plt.title("Calor 1d transiente implicito") # adicionando titulo ao plot
plt.show()
