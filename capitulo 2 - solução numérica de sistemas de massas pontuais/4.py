import matplotlib.pyplot as plt
#dados

m=1.0 #massa do carrinho
k=0.1 # coeficiente da mola
w=(k/m)**0.5 #frequencia
dt=0.1 #passo
t=0 #tempo inicial
vx0=10 #velocidade inicial
vx1=9.0
v_old=vx0
x=0 #posição inicial
graph=0 #0 v-t, 1 v-x
while t<100:
    v_old=vx0
    vx1=(-k*x*dt/m)+vx0
    x+=vx1*dt
    vx0=vx1
    t+=dt
    if graph==0:
        plt.plot(t,vx1,"bo")
        plt.xlabel("tempo[s]")
        plt.ylabel("velocidade [m/s]")
        plt.title("Sistema massa mola")
    if graph==1:
        plt.plot(x,vx1,"bo")
        plt.xlabel("x[m]")
        plt.ylabel("velocidade [m/s]")
        plt.title("Sistema massa mola")
plt.show()
