import matplotlib.pyplot as plt

#Dados
m=1.0 # massa do carrinho
b=0.1 #arrasto
passo=0.01 # passo do tempo
vx=10 # Velocidade inicial
vx1=9
v_old=10
s=0 # distância percorrida
t=0 # tempo
while v_old-vx1>10**-6:
    v_old=vx
    vx1=vx-(passo*b*vx/m)
    s+=vx1*passo
    vx=vx1
    t+=passo
    S=str(s)
    plt.plot(t,vx1,"bo")
print (S+" metros de distancia percorrida.")
plt.xlabel('Distância percorrida [m]')
plt.ylabel('Velocidade [m/s]')
plt.title("Movimento Horizontal de um Carrinho") # adicionando titulo ao plot
plt.show()
