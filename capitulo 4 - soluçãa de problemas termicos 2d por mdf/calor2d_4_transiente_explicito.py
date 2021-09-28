
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mesh


#inputs
nx = 25
ny = 25
npoints = nx*ny
ne = (nx-1)*(ny-1)
dx = 1.0/nx
dy = 1.0/ny
time=0.0 # tempo de simulacao
dt = 0.0006 # passo de tempo

mesh_info=mesh.mesh_2d(nx,ny)

#posição dos pontos

X = np.array(mesh_info[2]) 
Y = np.array(mesh_info[3])

#pontos de contorno

cc1=[];cc2=[];cc3=[];cc4=[]
for i in range(nx):
   cc1.append(i)
for i in range(nx):
   cc3.append((nx*(ny-1))+i)
for i in range(1,(ny-1)):
    cc2.append((nx-1)+nx*i)
for i in range(1,(ny-1)):
    cc4.append(nx*i)
cc = cc1 + cc2 + cc3 + cc4

#pontos internos

inner = []
for i in range(npoints):
    inner.append(i)
for i in cc:
    inner.remove(i)


# condicoes de contorno espaciais
bc = 0*np.ones( (npoints),dtype='float' )
for i in cc1:
 bc[i] = X[i]
for i in cc2:
 bc[i] = Y[i]*Y[i] + 1
for i in cc3:
 bc[i] = X[i]*X[i] + 1
for i in cc4:
 bc[i] = Y[i]

# condicao inicial (temporal)
T = 0*np.ones( (npoints),dtype='float' )
for i in cc:
 T[i] = bc[i]

# distribuicao de fonte de calor
Q = 0*np.ones( (npoints),dtype='float' )
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

T = T.reshape(ny,nx) #armando a matriz da Temperatura no espaço
Q = Q.reshape(ny,nx) #armando a matriz da geração de calor no espaço


# Calculo de T no tempo
for h in range(100000):
 Tp = T[1,1]
 for j in range(1,ny-1):
    for i in range(1,nx-1):
       T[i,j] = T[i,j] + \
         dt*(T[i+1,j]-2*T[i,j]+T[i-1,j])/dx**2 + \
         dt*(T[i,j+1]-2*T[i,j]+T[i,j-1])/dy**2 + dt*Q[i,j]     
 Tf = T[1,1]
 
 if (Tf-Tp) < 1E-5:
  print ("break")
  break


# # plot 2D cor (quadrilatero)
surf = plt.imshow(T, interpolation='quadric', origin='lower',
                  cmap=matplotlib.cm.jet, extent=(X.min(),
                  X.max(), Y.min(), Y.max()))
plt.colorbar(surf,shrink=1.0, aspect=20)
plt.grid(color='black', linestyle='solid', linewidth=0.5)
labx = np.linspace(X.min(),X.max(),nx)
laby = np.linspace(Y.min(),Y.max(),ny)
plt.xticks(labx)
plt.yticks(laby)


plt.title("Calor 2d transiente explicito") # adicionando titulo ao plot
plt.show()
