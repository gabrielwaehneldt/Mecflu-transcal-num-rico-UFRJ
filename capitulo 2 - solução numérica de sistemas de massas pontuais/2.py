import matplotlib.pyplot as plt

#dados iniciais
p=840 # massa especifica
d=(1.5*(10**-6))# diametro da gota
t=0 # tempo inicial
x=0 # posição inicial
beta=1.6*(10**-4)
dt=0.0000001
vx0=0 # velocidade inicial
vx1=9.0
v_old=vx0
m=(p*3.1415*((d)**3)/6.0)
b=beta*d

while abs(v_old-vx1)>10**-9.5:
    v_old=vx0
    vx1=vx0-(dt*b*vx0/m) + dt*9.81
    x+=vx1*dt
    vx0=vx1
    t+=dt
    plt.plot(t,vx1,"bo")
print (vx1)
plt.xlabel('Distância percorrida [m]')
plt.ylabel('tempo [s]')
plt.title("Velocidade terminal de uma gota") # adicionando titulo ao plot
plt.show()
        
