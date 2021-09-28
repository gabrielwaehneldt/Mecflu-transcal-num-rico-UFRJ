#Gera malha 2d aleatória com cor baseada em 2 critérios


from math import acos, degrees
from iteration_utilities import unique_everseen
from iteration_utilities import duplicates
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import timing

#Dados
nx = 3 ; ny = 2 #numero de pontos em x e y
xf = 3 ; yf = 2 #tamanho do dominio em x e y
npoints = nx*ny #Numero de pontos
ne = 2*((nx-1)*(ny-1)) #Numero de elementos para triangulo
mode=0 #1 - normal; 0 - equilateralidade
complemento=""

#funcoes

def normal_func(X,Y,IEN):
    ne = len(IEN)
    t_area=0
    areas=[]
    big_area=0
    flag2=0
    for i in range(0,ne):
        [v1,v2,v3]=IEN[i]
        flag=0
        x1 = X[v1]; x2 = X[v2]; x3 = X[v3] #coordenada x dos 3 vertices
        y1 = Y[v1]; y2 = Y[v2]; y3 = Y[v3] #coordenada y dos 3 vertices
        a=np.array([x2-x1,y2-y1]) ; b=np.array([x3-x2,y3-y2]) ; c=np.array([x1-x3,y1-y3]) #vetores
        n1=(-b[0]*a[1]+(a[0]*b[1]))/2.0 ; n2=(-c[0]*b[1]+(b[0]*c[1]))/2.0 ; n3=(-a[0]*c[1]+(c[0]*a[1]))/2.0 # produto vetorial
        if (n1>0 and n2>0 and n3>0) or (n1<0 and n2<0 and n3<0):
            flag=1
        if flag==0:
            IEN[i]=[v3,v2,v1]
            flag2=1
        if n3>big_area:
            big_area=n3
        t_area+=(n3)
        areas.append(abs(n3))
    t_area=abs(t_area)
    return (t_area,areas,big_area,flag2)

def mesh_quality_tri(X,Y,IEN,ne):
	crit_list=[]
	ne = IEN.shape[0]
	for i in range(0,ne):
		[v1,v2,v3]=IEN[i]
		x1 = X[v1]; x2 = X[v2]; x3 = X[v3] #coordenada x dos 3 vertices
		y1 = Y[v1]; y2 = Y[v2]; y3 = Y[v3] #coordenada y dos 3 vertices
		l1=((x2-x1)**2+(y2-y1)**2)**0.5
		l2=((x3-x2)**2+(y3-y2)**2)**0.5
		l3=((x1-x3)**2+(y1-y3)**2)**0.5
		a1=degrees(acos((l2 * l2 + l3 * l3 - l1 * l1)/(2.0 * l2 * l3)))
		a2=degrees(acos((l3 * l3 + l1 * l1 - l2 * l2)/(2.0 * l3 * l1)))
		a3=degrees(acos((l1 * l1 + l2 * l2 - l3 * l3)/(2.0 * l1 * l2)))
		crit=(((abs(a2-60))+ (abs(a3-60))+ (abs(a1-60))))
		crit_list.append(crit)
	return crit_list

def IEN_bound(IEN): #retorna a IEN dos vertices de fronteira
    IEN_bound=[]
    org_list=[]
    for y in IEN:
        org_list.append([y[0],y[1]])
        org_list.append([y[1],y[2]])
        org_list.append([y[2],y[0]])
    count=0
    for i in org_list:
        n=org_list[count:].count(i)
        n+=org_list[count:].count([i[1],i[0]])
        if n==1:
            IEN_bound.append(i)
        if n==2:
            org_list.remove([i[1],i[0]])
        count+=1
    IEN_bound.append([IEN_bound[-1][1],IEN_bound[0][0]])
    return IEN_bound

# geracao aleatoria de pontos
X = np.random.uniform(0.0,xf,npoints)
Y = np.random.uniform(0.0,yf,npoints)

# gerando malha de triangulos
Tri = matplotlib.tri.Triangulation(X,Y)
IEN = Tri.triangles

# gerando a IEN_bound                
IEN_bound=IEN_bound(IEN)

 
# plot da malha de triangulos
fig, ax = plt.subplots()
ax.set_aspect('equal')

for i in range(0,npoints):
 plt.text(X[i]+0.02,Y[i]+0.03,str(i),color='b')
if mode==0: #plota equilateralidade dos triangulos
    crit_list=mesh_quality_tri(X,Y,IEN,ne)
    plt.tripcolor(Tri, crit_list, edgecolors='k',cmap='Blues',vmin=0,vmax=200)
    complemento+=" com qualidade da malha"

if mode==1: #plota area dos triangulos e informacoes da normal
    normal_info=normal_func(X,Y,IEN)
    plt.tripcolor(X,Y,IEN, normal_info[1], edgecolors='k',cmap='Greens',vmin=0,vmax=normal_info[2])
    complemento+=" com indicador de área"
    if normal_info[3]==1:
        print ("-Orientacao da malha teve de ser corrigida")
    print ("-Area coberta pela malha de " + str(normal_info[0])+ " unidades")

count=0
while len(IEN_bound)-1>count: #plota a borda linha a linha
    xx=([X[IEN_bound[count][0]],X[IEN_bound[count][1]]])
    yy=([Y[IEN_bound[count][0]],Y[IEN_bound[count][1]]])
    count+=1    
    plt.plot(xx,yy,"-ko",color="red")
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Malha aleatoria " + complemento+ " e borda")
plt.show()
