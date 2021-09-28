
import matplotlib.pyplot as plt
import numpy as np

#dados iniciais

t=0 # tempo
v=30 # velocidade em metros por segundo
ang=50 # angulo da direção da velocidade em graus
x=0 # posição do projetil
y=0.4 # altura do projetil
d=0.07 # diametro projetil
V= (3.1415*d**3)/6 # volume do projetil
m=0.15 #massa do projetil
g=9.81 # aceleracao da gravidade
dt= 0.01 #passo
mode=0#0-Vácuo ; 1-atrito linear ; 2-atrito quadrático
colr=["green","red","black"]

#atrito linear dados
beta=0.16
b=beta*d

#atrito quadratico dados
eps=0.25
c=eps*(d**2)
vx0=v*np.cos(ang*3.1416/180)
vy0=v*np.sin(ang*3.1416/180)


def drag(vx0,vy0,eixo,mode):
    if mode==0:
        return 0
    if mode==1:
        if eixo=="x":
            return -b*vx0
        if eixo=="y":
            return -b*vy0
    if mode==2:
        F=-c*((vx0**2)+(vy0*2))
        ang=np.arctan(vy0/vx0)
        if eixo=="x":
            return F*np.cos(ang)
        if eixo=="y":
            return F*np.sin(ang)
        
while x<192:
    #print (x)
    vx1=vx0 +(drag(vx0,vy0,"x",mode)*dt/m)
    vy1=-g*dt + vy0 +(drag(vx0,vy0,"y",mode)*dt/m)
    x+=vx1*dt
    y+=vy1*dt
    t+=t
    if y<0:
        y=-y
        vy1=-vy1
    ang=np.arctan(vy1/vx1)
    vx0=vx1
    vy0=vy1
    if x>=192 and mode<2:
        x=0
        y=0
        vx0=v*np.cos(50*3.1416/180)
        vy0=v*np.sin(50*3.1416/180)
        mode+=1
    plt.plot(x,y,color=colr[mode], marker='o', markersize=4, label='cubic')
plt.xlabel("x[m]")
plt.ylabel("y[m]")
plt.title("Lançamento de projétil")
plt.show()
