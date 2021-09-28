import matplotlib.pyplot as plt
#dados

m=0.15 #massa do carrinho
k=0.1 # coeficiente da mola
dt=0.01 #passo
g=9.81 #aceleracao gravidade
t=0 #tempo inicial
vy0=10 #velocidade inicial
y=0 #posição inicial
graph=0 #0 v-t, 1 v-y
while t<30:
    vy1=(-k*y*dt/m)+vy0+g*dt
    y+=vy1*dt
    vy0=vy1
    t+=dt
    if graph==0:
        plt.plot(t,vy1,"bo")
        plt.xlabel("tempo[s]")
        plt.ylabel("velocidade [m/s]")
        plt.title("Sistema massa mola")
    if graph==1:
        plt.plot(y,vy1,"bo")
        plt.xlabel("y[m]")
        plt.ylabel("velocidade [m/s]")
        plt.title("Sistema massa mola vertical")
plt.show()
