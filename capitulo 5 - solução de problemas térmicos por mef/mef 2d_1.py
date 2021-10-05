import numpy as np
import sympy as sp
import mesh as ms
import matplotlib.pyplot as plt

# parametros da simulacao

nx = 5
ny = 5
npoints = nx*ny
ne = (nx-1)*(ny-1)
Q=100.0 #fonte de calor [W/m^3]
rho=1.0 #densidade do material [kg/m^3]
cv=1.0 #capacidade termica [J/kgK]
k=1.0 #condutividade termica
alpha=k/(rho*cv) #difusidade termica
theta=1.0
dt = 0.3 # passo de tempo
it= 10 #numero de iteracoes



# condicao de contorno de Dirichlet (espacial)


# geracao dos pontos e IEN por mesh

msh=ms.mesh_2d(nx,ny)
IEN=np.array(msh[0]) ; X=msh[2] ; Y=msh[3]
#print (X)



#Montando a matriz K e M

K=np.zeros((npoints,npoints),dtype='float')
M=np.zeros((npoints,npoints),dtype='float')
for i in range(ne):
    #vertices do triangulo
    v1=IEN[i,0]
    v2=IEN[i,1]
    v3=IEN[i,2]

    #coeficientes b e c
    bi=Y[v2]-Y[v3]
    bj=Y[v3]-Y[v1]
    bk=Y[v1]-Y[v2]
    ci=X[v3]-X[v2]
    cj=X[v1]-X[v3]
    ck=X[v2]-X[v1]
    

    area=(1.0/2.0)*np.linalg.det([[1,X[v1],Y[v1]],
                                  [1,X[v2],Y[v2]],
                                   [1,X[v3],Y[v3]]])
    #print (area)
    kxel=(1.0/(4*area))*np.array([[bi*bi,bi*bj,bi*bk],
                                 [bj*bi,bj*bj,bj*bk],
                                 [bk*bi,bk*bj,bk*bk]])
    kyel=(1.0/(4*area))*np.array([[ci*ci,ci*cj,ci*ck],
                                 [cj*ci,cj*cj,cj*ck],
                                 [ck*ci,ck*cj,ck*ck]])
    kel=alpha*(kxel+kyel)
    mel=(area/12.0)*np.array([[2,1,1],
                              [1,2,1],
                              [1,1,2]])
    for ilocal in range(3):
        iglobal=IEN[i,ilocal]
        for jlocal in range(3):
            jglobal=IEN[i,jlocal]

            K[iglobal,jglobal]+=kel[ilocal,jlocal]
            M[iglobal,jglobal]+=mel[ilocal,jlocal]
    #K[i:i+(3),i:i+(3)]+=kel #vai somando a matriz kel na matriz K ne vezes percorrendo na diagonal
    #M[i:i+(3),i:i+(3)]+=mel #vai somando a matriz mel na matriz M ne vezes percorrendo na diagonal
print (K)
#print (M)


#montando o vetor b
f=np.dot(M,(Q/(rho*cv))*np.ones((npoints),dtype="float"))
print (f)

# plot da condicao de contorno de T
T = np.zeros( (npoints),dtype='float' )
plt.plot(X,T,'r-') #plota condição de temperatura inicial


# resolvendo a matriz e achando T
T = np.linalg.solve(K,f)

#plotando as temperaturas
plt.plot(X,T,'ko-')
plt.xlabel('comprimento da barra [m]')
plt.ylabel('temperatura [oC]')
plt.title("Calor 1d mef com geração de calor") # adicionando titulo ao plot
plt.show()

