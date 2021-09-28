
import numpy as np
import matplotlib.pyplot as plt
import mesh
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import time
start = time.time() #inicia contagem de tempo do programa

#inputs
nx = 100
ny = 100
npoints = nx*ny
ne = (nx-1)*(ny-1)
dx = 1.0/nx
dy = 1.0/ny
dt = 0.3 # passo de tempo
it= 10 #numero de iteracoes

mesh_info=mesh.mesh_2d(nx,ny)

#posição dos pontos

X = np.array(mesh_info[2]) 
Y = np.array(mesh_info[3])

#pontos de contorno

cc1=np.array([],dtype="int");cc2=np.array([],dtype="int");cc3=np.array([],dtype="int");cc4=np.array([],dtype="int")
for i in range(nx):
   cc1=np.append(cc1,i)
for i in range(nx):
   cc3=np.append(cc3,(nx*(ny-1))+i)
for i in range(1,(ny-1)):
   cc2=np.append(cc2,(nx-1)+nx*i)
for i in range(1,(ny-1)):
   cc4=np.append(cc4,nx*i)
cc = np.concatenate((cc1,cc2,cc3,cc4))
#pontos internos
inner = []
for j in range(1,ny-1):
    for i in range(1,nx-1):
       inner.append(i+nx*j)


# condicoes de contorno espaciais
bc = 1*np.ones( (npoints),dtype='float' )
for i in cc1:
 bc[i] = X[i] 
for i in cc2:
 bc[i] = Y[i] + 1
for i in cc3:
 bc[i] = X[i] + 1
for i in cc4:
 bc[i] = Y[i]


# condicao inicial de temperatura incluindo as bordas
T = 0*np.ones( (npoints),dtype='float' )
for i in cc:
 T[i] = bc[i]# temperaturas de borda no array


# distribuicao de fonte de calor
Q = np.zeros( (npoints),dtype='float' )
for i in range(0,npoints):
 # cpu 1
 if X[i]<=0.4 and X[i] >= 0.2 and Y[i] <= 0.4 and Y[i]>=0.2: 
  Q[i] = 1000.0
 # cpu 2
 if X[i]<=0.8 and X[i] >= 0.6 and Y[i] <= 0.4 and Y[i]>=0.2: 
  Q[i] = 1000.0
 # cpu 3
 if X[i]<=0.8 and X[i] >= 0.6 and Y[i] <= 0.8 and Y[i]>=0.6: 
  Q[i] = 1000.0
 # cpu 4
 if X[i]<=0.4 and X[i] >= 0.2 and Y[i] <= 0.8 and Y[i]>=0.6: 
  Q[i] = 1000.0

# inicializar a matriz A e o vetor b
A = np.zeros( (5,npoints),dtype='float' )
b = np.zeros( (npoints),dtype='float' )


# populando os valores dos pontos de contorno em A em b
for i in cc:
  A[0,i] = 1.0
  b[i]   = bc[i]
# montando a matriz para o miolo
for i in inner:
  A[0,i]   = 1+2*( dt/(dx*dx) )+2*( dt/(dy*dy)) # diagonal principal
  A[1,i-1] =  -( dt/(dx*dx) )  # diagonal inferior -1 (porque i-1?)
  A[2,i] =  -( dt/(dx*dx) )  # diagonal superior +1
  A[3,i]=  -( dt/(dy*dy) ) #diagonal superior_2 +nx
  A[4,i-nx]=  -( dt/(dy*dy) ) #diagonal inferior_2 -nx (porque i-nx?)
A=diags(A, [0,-1,1,nx,-nx],shape=(npoints, npoints),format="csr",dtype="float") #criando matriz esparsa A

for n in range(it):
 print ("Ciclo:"+str(n+1))

 # populando os valores dos pontos internos de b no tempo 
 for i in inner:
  b[i]     =T[i] +(dt*Q[i])

 T = spsolve(A,b) #resolvendo matriz esparsa


# # plot 2D cor (quadrilatero)
T = T.reshape(ny,nx) #armando a matriz da Temperatura no espaço
surf = plt.imshow(T, interpolation='quadric', origin='lower',
                  cmap=plt.cm.jet, extent=(X.min(),
                  X.max(), Y.min(), Y.max()))
plt.colorbar(surf,shrink=1.0, aspect=20)
plt.grid(color='black', linestyle='solid', linewidth=0.5)
labx = np.linspace(X.min(),X.max(),nx)
laby = np.linspace(Y.min(),Y.max(),ny)
plt.xticks(labx)
plt.yticks(laby)

end = time.time() #termina contagem de tempo do programa
print ("Exucution time of the program is- ",end-start) #printa tempo gasto no intervalo

plt.title("Calor 2d transiente implicito (sparse)") # adicionando titulo ao plot
plt.show()
