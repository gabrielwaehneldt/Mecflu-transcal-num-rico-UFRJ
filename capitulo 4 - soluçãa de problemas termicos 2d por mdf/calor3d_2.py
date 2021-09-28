import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mesh



Lx = 1.0
Ly = 1.0
Lz = 1.0
nx = 10
ny = 10
nz = 20
npoints = nx*ny*nz
ne = (nx-1)*(ny-1)*(nz-1)
dx = Lx/(nx-1)
dy = Ly/(ny-1)
dz = Lz/(nz-1)

# geracao de malha

mesh_info=mesh.mesh_3d(nx,ny,nz)
X=mesh_info[0] ; Y=mesh_info[1] ; Z=mesh_info[2] #criando as coordanadas dos pontos


# pontos de contorno
cc1=[];cc2=[];cc3=[];cc4=[];cc5=[];cc6=[];inner=[]
for r in range (1,nz-1):
    for i in range(nx):
        cc1.append(nx*ny*r+i)               #baixo

for r in range(1,nz-1):
    for i in range(1,ny-1):
        cc2.append(nx*ny*r+ny*i+(ny-1))     #direita
                                   
for r in range(1,nz-1):
    for i in range(nx):
        cc3.append(nx*ny*r+(nx*(ny-1))+i)   #topo

for r in range(1,nz-1):
    for i in range(1,ny-1):
        cc4.append(nx*ny*r+ny*i)            #esquerda

for i in range(nx*ny):
    cc5.append(i)                           #frente
for i in range(nx*ny*(nz-1),nx*ny*nz):
    cc6.append(i)                           #costas
cc = cc1 + cc2 + cc3 + cc4 + cc5 + cc6

# pontos internos
for i in range(npoints):
    inner.append(i)
for i in cc:
    inner.remove(i)



#-------------------------------------------------- 
# plt.plot(X,Y,'bo')
# plt.plot(X[cc],Y[cc],'ko')
# plt.show()
#-------------------------------------------------- 

bc = 0*np.ones( (npoints),dtype='float' )
for i in cc1:
 bc[i] = X[i]+Z[i]
for i in cc2:
 bc[i] = 3
for i in cc3:
 bc[i] = 4
for i in cc4:
 bc[i] = Z[i]+Y[i]
for i in cc5:
 bc[i] = 3
for i in cc6:
 bc[i] = 2+Y[i]
#--------------------------------------------------
# # # plot 2D cor (quadrilatero)
# Z = bc.reshape(ny,nx)
# surf = plt.imshow(Z, interpolation='quadric', origin='lower',
#                   cmap=matplotlib.cm.jet, extent=(X.min(),
#                   X.max(), Y.min(), Y.max()))
# plt.colorbar(surf,shrink=1.0, aspect=20)
# plt.grid(color='black', linestyle='solid', linewidth=0.5)
# labx = np.linspace(X.min(),X.max(),nx)
# laby = np.linspace(Y.min(),Y.max(),ny)
# plt.xticks(labx)
# plt.yticks(laby)
# plt.show()
#-------------------------------------------------- 

A = np.zeros( (npoints,npoints),dtype='float' )
b = np.zeros( (npoints),dtype='float' )

# Eq. da cond. de contorno de Dirichlet
for i in cc:
 A[i,i] = 1.0
 b[i]   = bc[i]

# Eqs. pontos internos (inner):
for i in inner:
 A[i,i-1] = 1/(dx*dx)                       # direita
 A[i,i] = -2/(dx*dx) -2/(dy*dy) -2/(dz*dz)  # ele mesmo
 A[i,i+1] = 1/(dx*dx)                       # esquerda
 A[i,i-nx] = 1/(dy*dy)                      # baixo
 A[i,i+nx] = 1/(dy*dy)                      # cima
 A[i,i-(nx*ny)] = 1/(dz*dz)                 # frente
 A[i,i+(nx*ny)] = 1/(dz*dz)                 # tras
print (A)
# solucao do sistema linear Ax=b
start = time.time() #inicia contagem de tempo do programa
T = np.linalg.solve(A,b)
end = time.time() #termina contagem de tempo do programa
# plot 3D (plano) cor (cubo)
Z = T.reshape(ny*nz,nx)
surf = plt.imshow(Z, interpolation='quadric', origin='lower',
                  cmap=matplotlib.cm.jet, extent=(X.min(),
                  X.max(), Y.min(), Y.max()))
plt.colorbar(surf,shrink=1.0, aspect=20)
plt.grid(color='black', linestyle='solid', linewidth=0.5)
labx = np.linspace(X.min(),X.max(),nx)
laby = np.linspace(Y.min(),Y.max(),ny)
plt.xticks(labx)
plt.yticks(laby)


print ("Exucution time of the program is- ",end-start) #printa tempo gasto no intervalo
plt.title("Calor 3d") # adicionando titulo ao plot
plt.show()
