import matplotlib.pyplot as plt
import numpy as np

#dados
y=0.2
dt=0.01 #passo
w=6.28 #frequencia motriz
w0=w*1.5    #frequencia natural
beta2=w0/4.0
t=0.0 #tempo inicial
vx0=0.0 #velocidade inicial
ang0=0.0 #posição inicial [rad]
mode=0 #-0 ang/t -1 v/ang

while t<6:
    ang1=ang0+dt*vx0
    print (vx0)
    vx1=vx0-dt*((w0**2)*np.sin(ang1)+beta2*vx0-y*(w0**2)*np.cos(w*t))
    vx0=vx1
    ang0=ang1
    t+=dt
    if mode==0:
        #plt.plot(t,ang1,color='black', marker='o', linestyle='dashed',linewidth=2, markersize=1)
        plt.plot(t,ang1,"bo-")
    else:
        plt.plot(ang1,vx1,color='black', marker='o', linestyle='dashed',linewidth=2, markersize=1)
if mode==0:
    
    plt.xlabel("tempo[s]")
    plt.ylabel("angulo [rad]")
else:
    plt.xlabel("ang [rad]")
    plt.ylabel("velocidade [rad/s]")
plt.title("Pendulo forcado amortecido")
plt.rc("text",usetex=True)
plt.show()
