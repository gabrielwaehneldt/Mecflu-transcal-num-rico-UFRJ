import matplotlib.pyplot as plt
import numpy as np

#dados
m=1.0 #massa do carrinho
g=9.81 #aceleracao gravidade
dt=0.01 #passo
l=1.0 #comprimento da haste
t=0 #tempo inicial
vx0=0 #velocidade inicial
ang0=0.175 #posição inicial [rad]
mode=0 #-0 ang/t -1 v/ang

while t<10:
    ang1=ang0+dt*vx0
    #print (ang1)
    vx1=vx0-dt*(g/l)*np.sin(ang1)
    vx0=vx1
    ang0=ang1
    t+=dt
    if mode==0:
        plt.plot(t,ang1,color='red', marker='o', linestyle='dashed',linewidth=2, markersize=1)
    else:
        plt.plot(ang1,vx1,color='blue', marker='o', linestyle='dashed',linewidth=2, markersize=1)
if mode==0:
    
    plt.xlabel("tempo[s]")
    plt.ylabel("angulo [rad]")
else:
    plt.xlabel("ang [rad]")
    plt.ylabel("velocidade [rad/s]")
plt.title("Pendulo simples")
plt.show()
