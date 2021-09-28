import matplotlib.pyplot as plt
import numpy as np

#dados
x=15 #posicao x
y=-25 #posicao y
z=-20 #posicao z
t=0.0 #tempo inicial
dt=0.001 #passo
p=28
teta=10
beta=8.0/3.0
graph=1 #0:x-y, 1:x-z, 2:y-z

while t<100:
    x1=x+dt*teta*(y-x)
    y1=y+dt*x*(p-z)-y*dt
    z1=z+dt*x*y-beta*z*dt
    t+=dt
    if graph==0:
        #print (x1)
        plt.plot(x,y, color='black', marker='o', linestyle='dashed',linewidth=2, markersize=0.2)
        #plt.plot(t,ang1,"bo-")
    if graph==1:
        plt.plot(x,z,color='black', marker='o', linestyle='dashed',linewidth=2, markersize=1)
    if graph==2:
        plt.plot(y,z,color='black', marker='o', linestyle='dashed',linewidth=2, markersize=1)
    x=x1
    y=y1
    z=z1
if graph==0: 
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
if graph==1:
    plt.xlabel("x [m]")
    plt.ylabel("z [m]")
if graph==2:
    plt.xlabel("y [m]")
    plt.ylabel("z [m]")
plt.title("Borboleta de Lorenz")
#plt.rc("text",usetex=True)
plt.show()
