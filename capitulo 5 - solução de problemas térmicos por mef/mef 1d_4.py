import numpy as np
import sympy as sp
import mesh as ms
import matplotlib.pyplot as plt

# parametros da simulacao
L = 1.0
npoints = 5
ne = npoints-1
#x = sp.symbols("x")
mode_X=2 #1 linear; 2 quadratico


# condicao de contorno de Dirichlet (espacial)
Te = 100.0
Td = 0.0

# geracao dos pontos e IEN por mesh

msh=ms.mesh_1d(npoints,L,mode_X)
X=msh[1] ; IEN=np.array(msh[0])
print (IEN)



#Montando a matriz k e A

A=np.zeros((npoints,npoints),dtype='float')
k=np.zeros((2,2),dtype='float')
for i in range(ne):
    v1=IEN[i,0]
    v2=IEN[i,1]
    h=X[v2]-X[v1]
    k=1/h*np.array([[1,-1],
                   [-1,1]])
    A[i:i+(2),i:i+(2)]+=-k #vai somando a matriz k na matriz A ne vezes percorrendo na diagonal


A[0]=A[0]*0 #zerando a linha de baixo e cima e colocando 1 nas 2 pontas da diagonal principal de A
A[-1]=A[-1]*0
A[0][0]=1
A[-1][-1]=1

#montando o vetor b
b = np.zeros( (npoints),dtype='float' )
b[0]=Te
b[-1]=Td #condições de temperatura na ponta

# plot da condicao de contorno de T
plt.plot(X,b,'r-')
#print (A)
#print (b)
# resolvendo a matriz e achando T
T = np.linalg.solve(A,b)

#plotando as temperaturas
plt.plot(X,T,'ko-')
plt.xlabel('comprimento da barra [m]')
plt.ylabel('temperatura [oC]')
plt.title("Calor 1d mef") # adicionando titulo ao plot
plt.show()
