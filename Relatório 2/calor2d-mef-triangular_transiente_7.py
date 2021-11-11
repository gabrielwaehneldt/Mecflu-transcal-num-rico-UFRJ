#organizar o programa

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mesh as ms
import meshio

# parametros da simulacao

k_FR4 = 0.9
K_copper = 400 #W/mK
k = k_FR4*0.9 + K_copper*0.1
rho = 8800*0.1 + 1100*0.9 #g/cm^3
cv_copper  = 380 #J/gK
cv_FR4 = 1100
cv = cv_copper*0.1 + cv_FR4*0.9
kappa_x = 1.0
kappa_y = 1.0
alpha_x = kappa_x/(rho*cv)
alpha_y = kappa_y/(rho*cv)
q = 1000.0 # fonte de calor
teta=1 #parametro para controle do mdf temporal
dt=0.01 #passo do tempo
it=6 #numero de iterações
q2 = 9306.0 # calor cpu W/m^2
q3 = 5114.0 # calor do vcore
q4 = 4573.0 # calor do vigpu
q5 = 262.0  # calor da ram
q6 = 2624.0 # calor vddr
q7 = 14588.0# calor placa de video
q8 = 7031.0 # calor do chipset


#importando e lendo a malha do gmsh

msh = meshio.read('Placa_mae_teste - Copia.msh')
X = msh.points[:,0]
Y = msh.points[:,1]
IEN = msh.cells['triangle']
npoints = len(X) 
ne = IEN.shape[0]
regions = msh.cell_data['triangle']['gmsh:geometrical']


# captura os pontos de contorno e gera cc
cc=[]
for i in range(npoints):
    if X[i]==0.185 or X[i]==0.0 or Y[i]==0.226 or Y[i]==0.0:
        cc.append(i)

# inicializacao das matrizes e vetores
K = np.zeros( (npoints,npoints),dtype='float' )
M = np.zeros( (npoints,npoints),dtype='float' )

for e in range(0,ne):
 # definicao das matrizes elementares
 [v1,v2,v3] = IEN[e]
 x1 = X[v1]; x2 = X[v2]; x3 = X[v3] #coordenada x dos 3 vertices
 y1 = Y[v1]; y2 = Y[v2]; y3 = Y[v3] #coordenada y dos 3 vertices
 a=np.array([x2-x1,y2-y1]) ; b=np.array([x3-x2,y3-y2]) #vetores
 Area=(-b[0]*a[1]+(a[0]*b[1]))/2.0
 c1=x3-x2 ; c2=x1-x3 ; c3=x2-x1
 b1=y2-y3 ; b2=y3-y1 ; b3=y1-y2
 melem = (Area/12.0)*np.array([ [2,1,1],
                           [1,2,1],
                           [1,1,2] ])


 kx = (1/(4.0*Area))*np.array([ [b1**2,b1*b2,b1*b3],
                              [b2*b1,b2**2,b2*b3],
                              [b3*b1,b3*b2,b3**2]])

 ky = (1/(4.0*Area))*np.array([ [c1**2,c1*c2,c1*c3],
                              [c2*c1,c2**2,c2*c3],
                              [c3*c1,c3*c2,c3**2]])

 kelem=kx+ky
 
 for ilocal in range(0,3):
  iglobal = IEN[e,ilocal]
  for jlocal in range(0,3):
   jglobal = IEN[e,jlocal]

   K[iglobal,jglobal] += kelem[ilocal,jlocal]
   M[iglobal,jglobal] += melem[ilocal,jlocal]

print (kelem)
print (melem)

b = np.zeros( (npoints),dtype='float' )
T = np.zeros( (npoints),dtype='float' )
A = M/dt + teta*K.copy()
#Montando o vetor Q
q=100*np.ones(npoints)
for i in range(ne):
    if regions[i]!=1:
        if regions[i]==2:
            q[IEN[i][0]]=q2
            q[IEN[i][1]]=q2
            q[IEN[i][2]]=q2
        if regions[i]==8:
            q[IEN[i][0]]=q8
            q[IEN[i][1]]=q8
            q[IEN[i][2]]=q8
        if regions[i]==3:
            q[IEN[i][0]]=q3
            q[IEN[i][1]]=q3
            q[IEN[i][2]]=q3
        if regions[i]==7:
            q[IEN[i][0]]=q7
            q[IEN[i][1]]=q7
            q[IEN[i][2]]=q7
        if regions[i]==6:
            q[IEN[i][0]]=q6
            q[IEN[i][1]]=q6
            q[IEN[i][2]]=q6
        if regions[i]==5:
            q[IEN[i][0]]=q5
            q[IEN[i][1]]=q5
            q[IEN[i][2]]=q5
        if regions[i]==4:
            q[IEN[i][0]]=q4
            q[IEN[i][1]]=q4
            q[IEN[i][2]]=q4
Q = (q/(rho*cv))
MQ= M@Q

# imposicao das c.c.s de Dirichlet
for i in cc:
 A[i,:] = 0.0
 A[i,i] = 1.0
 b[i] = 40

 


for n in range(it):
 print ("Ciclo:"+str(n+1))

 # populando os valores dos pontos internos de b no tempo 

 
 
 b=M@(T/dt) + MQ - K@((1-teta)*T)

# imposicao das c.c.s de Dirichlet
 for i in cc:
  b[i] = 40
 
  

 T = np.linalg.solve(A,b) #resolvendo matriz


# # plot 2D cor (triangular)
triang = matplotlib.tri.Triangulation(X,Y,IEN)
ax = plt.axes()
ax.set_aspect("equal")
ff = ax.tricontourf(triang,T,cmap="jet")
plt.colorbar(ff)
plt.title("Calor 2d MEF") # adicionando titulo ao plot
plt.show()
