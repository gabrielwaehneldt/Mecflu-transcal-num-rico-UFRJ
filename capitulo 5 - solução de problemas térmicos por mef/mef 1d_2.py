import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# parametros da simulacao
L = 1.0
npoints = 15
ne = npoints-1
x = sp.symbols("x")
mode=1 #1-N linear; 2-N quadratico


# condicao de contorno de Dirichlet (espacial)
Te = 100.0
Td = 0.0

# geracao dos pontos
X = np.linspace(0,L,npoints)

# geracao da matriz de conectividade
IEN = np.zeros( (ne,2),dtype='int' )
for e in range(0,ne):
 IEN[e] = [e,e+1]

#geração do N (codigo de aluno)
def Ni(X, mode=1):

    x = sp.symbols("x")

    N=[]

    for k in range(len(X)-1):

        e_pos = np.linspace(X[k],X[k+1],mode+1)

        A=np.ones((mode+1,mode+1))
        for i in range(len(A)):
            for j in range(len(A)):
                if j!=(len(A)-1):
                    A[i,j]=e_pos[i]**(mode-j)
        for p in range(len(e_pos)):
            b=np.zeros(mode+1)
            b[p]=1

            n=list(np.linalg.solve(A,b))
            f=0
            for i in range(len(n)):
                f+=n[i]*x**(mode-i)
            N.append(f)
    return N


#Montando a matriz k
N=(Ni(X, mode))
k=np.zeros((mode+1,mode+1),dtype='float')
for i in range(mode+1):
    for j in range(mode+1):
        der_0=sp.diff(N[i],x) #derivando N
        der_1=sp.diff(N[j],x) #derivando N
        inte=sp.integrate(der_0*der_1, (x, 0, L/ne)) #integrando as derivadas
        k[i,j]=inte #inserindo as integrais na matriz k


#montando matriz A

A=np.zeros((npoints,npoints),dtype='float')
for i in range(npoints-mode): #somando a matriz k na A
    A[i:i+(mode+1),i:i+(mode+1)]+=k

A[0]=A[0]*0 #zerando a linha de baixo e cima e colocando 1 nas quinas de A
A[-1]=A[-1]*0
A[0][0]=1
A[-1][-1]=1 

#montando o vetor b
b = np.zeros( (npoints),dtype='float' )
b[0]=Te
b[-1]=Td #condições de temperatura na ponta

# plot da condicao de contorno de T
plt.plot(X,b,'r-')

# resolvendo a matriz e achando T
T = np.linalg.solve(A,b)

#plotando as temperaturas
plt.plot(X,T,'ko-')
plt.xlabel('comprimento da barra [m]')
plt.ylabel('temperatura [oC]')
plt.title("Calor 1d mef") # adicionando titulo ao plot
plt.show()
