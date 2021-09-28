import matplotlib.pyplot as plt
#dados

m=1.0 #massa do carrinho
mi=0.5 # amortecimento
dt=0.02 #passo
t=0 #tempo inicial
vx0=0.1 #velocidade inicial
x=0.2 #posição inicial
mode=0 #-0 v/t -1 v/x

while t<100:
    k1=(dt*mi*(1-x**2)*vx0)-dt*x
    k2=dt*mi*(1-x**2)*(vx0+k1/2)-dt*x
    vx1=vx0+(0.5*(k1+k2))
    x+=vx1*dt
    vx0=vx1
    t+=dt
    if mode==0:
        plt.plot(t,vx1,color='red', marker='o', linestyle='dashed',linewidth=2, markersize=1)
    else:
        plt.plot(x,vx1,color='blue', marker='o', linestyle='dashed',linewidth=2, markersize=1)
if mode==0:
    plt.xlabel("tempo[s]")
    plt.ylabel("velocidade [m/s]")
else:
    plt.xlabel("x[m]")
    plt.ylabel("velocidade [m/s]")
plt.title("Oscilador de Van Der Pol")
plt.show()
