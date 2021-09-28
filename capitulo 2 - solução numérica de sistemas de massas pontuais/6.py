import matplotlib.pyplot as plt
#dados

m=1.0 #massa do carrinho
beta=0.5 # amortecimento
w0=0.5 #frequencia
dt=0.1 #passo
t=0 #tempo inicial
vx0=10 #velocidade inicial
x=0 #posição inicial
while t<34:
    vx1=-2*beta*vx0*dt-(w0**2)*x*dt+vx0
    x+=vx1*dt
    vx0=vx1
    t+=dt
    plt.plot(t,x,"bo")
    plt.xlabel("tempo[s]")
    plt.ylabel("posicao [m]")
    plt.title("Sistema massa mola amortecedor")
plt.show()
