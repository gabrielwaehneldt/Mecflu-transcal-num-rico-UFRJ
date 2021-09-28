
import time ; from math import acos, degrees ; import numpy as np
import matplotlib.pyplot as plt ; import matplotlib ; import meshio

#Funcoes secundárias ####################################################################
def normal_func(X,Y,IEN,mesh):
        ne = len(IEN)
        t_area=0 #somatorio das areas dos triangulos
        areas=[] #lista das areas dos triangulos
        flag2=0 #flag para saber se houve necessidade de correcao do sentido da malha
        big_area=0 #area do maior triangulo
        for i in range(0,ne):
            if mesh=="quad":
                [v1,v2,v3,v4]=IEN[i]
                flag=0
                x1 = X[v1]; x2 = X[v2]; x3 = X[v3]; x4 = X[v4] #coordenada x dos 3 vertices
                y1 = Y[v1]; y2 = Y[v2]; y3 = Y[v3]; y4 = Y[v4] #coordenada y dos 3 vertices
                a=np.array([x2-x1,y2-y1]) ; b=np.array([x3-x2,y3-y2]) ; c=np.array([x4-x3,y4-y3]) ; d=np.array([x1-x4,y1-y4]) #vetores
                n1=(-b[0]*a[1]+(a[0]*b[1])) ; n2=(-c[0]*b[1]+(b[0]*c[1])) ; n3=(-d[0]*c[1]+(c[0]*d[1])) ; n4=(-a[0]*d[1]+(d[0]*a[1])) # produto vetorial
                if (n1>0 and n2>0 and n3>0 and n4>0) or (n1<0 and n2<0 and n3<0 and n4<0): #se a malha tiver sentido condizente ira entrar
                    flag=1
                if flag==0:
                    IEN[i]=[v3,v4,v1,v2]
                    flag2=1
            if mesh=="tri":
                [v1,v2,v3]=IEN[i]
                flag=0
                x1 = X[v1]; x2 = X[v2]; x3 = X[v3] #coordenada x dos 3 vertices
                y1 = Y[v1]; y2 = Y[v2]; y3 = Y[v3] #coordenada y dos 3 vertices
                a=np.array([x2-x1,y2-y1]) ; b=np.array([x3-x2,y3-y2]) ; c=np.array([x1-x3,y1-y3]) #vetores
                n1=(-b[0]*a[1]+(a[0]*b[1]))/2.0 ; n2=(-c[0]*b[1]+(b[0]*c[1]))/2.0 ; n3=(-a[0]*c[1]+(c[0]*a[1]))/2.0 # produto vetorial
                if (n1>0 and n2>0 and n3>0) or (n1<0 and n2<0 and n3<0): #se a malha tiver sentido condizente ira entrar
                    flag=1
                if flag==0:
                    IEN[i]=[v3,v2,v1]
                    flag2=1
            if n3>big_area:
                big_area=n3
            t_area+=(n3)
            areas.append([0.1,0.6,0.3,abs(n3)])
        t_area=abs(t_area)
        for e in range(len(areas)):
            areas[e][3]=areas[e][3]/big_area
        return (t_area,areas,flag2) #retorna lista com area total, lista com area de cada triangulo, maior area de um triangulo, flag para necessidade de correcao da orientacao da malha    

def mesh_quality_tri(X,Y,IEN): #analisa qualidade de malhas triangulares
        crit_list=[]
        ne = IEN.shape[0]
        for i in range(0,ne):
            [v1,v2,v3]=IEN[i]
            x1 = X[v1]; x2 = X[v2]; x3 = X[v3] #coordenada x dos 3 vertices
            y1 = Y[v1]; y2 = Y[v2]; y3 = Y[v3] #coordenada y dos 3 vertices
            l1=((x2-x1)**2+(y2-y1)**2)**0.5 #comprimentos dos lados
            l2=((x3-x2)**2+(y3-y2)**2)**0.5
            l3=((x1-x3)**2+(y1-y3)**2)**0.5
            a1=degrees(acos((l2 * l2 + l3 * l3 - l1 * l1)/(2.0 * l2 * l3))) #calculo dos angulos
            a2=degrees(acos((l3 * l3 + l1 * l1 - l2 * l2)/(2.0 * l3 * l1)))
            a3=degrees(acos((l1 * l1 + l2 * l2 - l3 * l3)/(2.0 * l1 * l2)))
            crit=(((abs(a2-60))+ (abs(a3-60))+ (abs(a1-60))))/200#criterio de equilateralidade
            crit_list.append([0.0,0.1,0.6,crit])
        return crit_list

def mesh_quality_quad(X,Y,IEN): #analisa qualidade de malhas triangulares
    crit_list=[]
    crit_large=0
    ne = IEN.shape[0]
    for i in range(0,ne):
        [v1,v2,v3,v4]=IEN[i]
        x1 = X[v1]; x2 = X[v2]; x3 = X[v3]; x4 = X[v4] #coordenada x dos 3 vertices
        y1 = Y[v1]; y2 = Y[v2]; y3 = Y[v3]; y4 = Y[v4] #coordenada y dos 3 vertices
        l1=((x2-x1)**2+(y2-y1)**2)**0.5 #comprimentos dos lados
        l2=((x3-x2)**2+(y3-y2)**2)**0.5
        l3=((x4-x3)**2+(y4-y3)**2)**0.5
        l4=((x1-x4)**2+(y1-y4)**2)**0.5
        l_medio=(l1+l2+l3+l4)/4
        crit=((abs(l1-l_medio))+ (abs(l2-l_medio))+ (abs(l3-l_medio)) + (abs(l4-l_medio)))#criterio de equilateralidade
        if crit_large<crit:
            crit_large=crit
        crit_list.append([0.0,0.1,0.6,crit])
    for i in range(len(crit_list)):
        crit_list[i][3]=(crit_list[i][3]/crit_large)
    return crit_list


def gaussian(x,x0,sigma):
        a=1/(sigma*(2*np.pi)**0.5)
        return a*np.exp(-np.power((x - x0)/sigma, 2.)/2.)

def tri_IEN(nx,ne): #cria IEN para triangulos
        I=[]
        c0=0 ; c1=1 ; c2=nx ; d0=nx+1
        count=0 ; count2=0
        for i in range(ne):
            I.append([c0+count,c1+count,c2+count])
            I.append([c0+count+1,d0+count,c2+count])
            count+=1 ; count2+=1
            if count2==(nx-1):
                count+=1 ; count2=0
        IEN=np.array(I)
        return (IEN)   

def quad_IEN(nx,ne): #cria IEN para quadrados
        I=[]
        c0=0 ; c1=1 ; c2=nx+1 ; c3=nx
        count=0 ; count2=0
        for i in range(ne):
            I.append([c0+count,c1+count,c2+count,c3+count])
            count+=1 ; count2+=1
            if count2==(nx-1):
                count+=1 ; count2=0 
        IEN=np.array(I)
        return (IEN)



def perturbacao(X,Y,A,teta,lamb,xf,yf,nx,ny): # perturbacao dos pontos da malha
        for e in range(ny-1):
            for i in range(nx):
                Y[nx*(ny-1-e)+i]=Y[nx*(ny-1-e)+i]+(A*np.sin((2*3.1415*X[i]/lamb)-teta))*((ny-1)-e)/(ny-1)
        return "Perturbacao realizada"

def IEN_bound_gen(IEN,mesh): #retorna a IEN dos vertices de fronteira
        IEN_bound=[]
        org_list=[]
        for y in IEN:
            org_list.append([y[0],y[1]])
            org_list.append([y[1],y[2]])
            if mesh=="tri":
                org_list.append([y[2],y[0]])
            if mesh=="quad":
                org_list.append([y[2],y[3]])
                org_list.append([y[3],y[0]])
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
        print (IEN_bound)
        return IEN_bound

#############################################################################################
#funcoes primarias ##########################################################################
def mesh_1d(nx,xf=1,mode=1):
    IEN=[]
    # geracao da malha 1D uniforme
    if type(mode)==int or type(mode)==float  or mode=="exp":    
            if mode=="exp":
                mode=2.718281828
            X = np.linspace(0,xf**(1/mode),nx)
            X = (X**mode)
    if mode=="gauss":
        count=0
        total=0
        passo=xf/(nx*10)
        x=0
        X = np.linspace(0,xf,nx)
        count2=0
        while count <nx*10:
            total+=xf*gaussian(x,xf/2.0,xf/6.0)*passo
            x+=passo
            count+=1
            if total>=X[count2]:
                X[count2]=x
                count2+=1
    plt.plot(X,0*X,'bo-')
    plt.show()
    # Criando a IEN
    for v in range(nx-1):
        IEN.append([v,v+1])
    #print (IEN)
    return [IEN,X]
#mesh_2d(numero de pontos em x,igual para y,tamanho do dominio em x,igual para y,Amplitude perturbacao,fase perturbacao,comprimento de onda, distribuição de pontos em x("exp", "gauss", qualquer expoente), igual para y,
#tipo der malha, 1-ligado ; 0-desligado visualiza a qualidade da malha triangular por cor /tem prioridade quanto a normal/,1-ligado ; 0-desligado ; liga uma borda visualizando os vertices de fronteira,
#1-ligado ; 0- desligado ; organiza malhas desorganizadas e visualiza areas dos elementos,1-ligado ; 0- desligado ; gera e propaga perturbacao senoide na malha,1-ligado ; 0- desligado ; retira pontos e IEN de arquivo externo
#nome do arquivo com a malha

def mesh_2d(nx,ny,plot=0,xf=1,yf=1,A=0.5,teta=1.57075,lamb=0.2,modex=1,modey=1,mesh="quad",mesh_quality=0,borda=0,normal=0,pertub=0,ext_mesh=0,arquivo='malha-furo.msh'):

    npoints = nx*ny #Numero de pontos
    ne = (nx-1)*(ny-1) #Numero de elementos para quadrilatero
    lisx=[] ; lisy=[]
    complemento="" #informacoes no titulo do plot

    if mesh=="tri":
        complemento+=" triangular"
    if mesh=="quad":
        complemento+=" retangular"
    
    
    

    # distribuicao de pontos #########################################################
    if ext_mesh==1: #pega os pontos e IEN de arquivo externo
        msh = meshio.read(arquivo)
        X = msh.points[:,0]
        Y = msh.points[:,1]
        IEN =(msh.cells_dict["triangle"])
        IEN=np.array(IEN)
        ne = IEN.shape[0]
        npoints = len(X)

    if ext_mesh==0: #so entra caso nao tenha uma malha externa/ gera os pontos
        if type(modex)==int or type(modex)==float  or modex=="exp":    
            if modex=="exp":
                modex=2.718281828
            X = np.linspace(0,xf**(1/modex),nx)
            X = (X**modex)
        if modex=="gauss" or modex=="inv_gauss":
            exp=1
            if modex=="inv_gauss":
                exp=-1
            count=0 ; count2=0
            total=0
            passo=xf/(nx*10)
            x=0
            X = np.linspace(0,xf,nx)
            while count <nx*10:
                total+=xf*gaussian(x,xf/2.0,xf/6.0)*passo
                x+=passo
                count+=1
                if total>=X[count2]: 
                    X[count2]=x
                    count2+=1
        for i in range(0,ny): #transformando em lista
            for i in range(0,nx):
                lisx.append(X[i])
        if type(modey)==int or type(modey)==float or modey=="exp":    
            if modey=="exp":
                modey=2.718281828
            Y = np.linspace(0,yf**(1/modey),ny)
            Y = (Y**modey)
        if modey=="gauss" or modey=="inv_gauss":
            exp=1
            if modey=="inv_gauss":
                exp=-1
            count=0 ; count2=0
            total=0
            passo=yf/(ny*10)
            y=0
            Y = np.linspace(0,yf,ny)
            while count <ny*10:
                total+=yf*gaussian(y,yf/2.0,yf/6.0,exp)*passo
                y+=passo
                count+=1
                if total>=(Y[count2]): 
                    Y[count2]=y
                    count2+=1
        for r in range(0,ny): #transformando em lista
            for i in range(0,nx):
                lisy.append(Y[r])

        X=lisx
        Y=lisy

    #perturbacao dos pontos
    if pertub==1 and ext_mesh==0:
        p=perturbacao(X,Y,A,teta,lamb,xf,yf,nx,ny)

    # IEN para quadrilateros #########################################################

    if mesh== "quad" and ext_mesh==0:
        IEN=quad_IEN(nx,ne)
    # IEN para triangulos #########################################################

    if mesh=="tri" and ext_mesh==0:
        IEN=tri_IEN(nx,ne)
    # plota malha de quadrilateros e triangulos ################################################
    if plot==1:
            cor="ko"
            if mesh=="tri":
                cor="o"
            xy = np.stack((X, Y), axis=-1)
            verts = xy[IEN]
            ax=plt.gca()
            pc = matplotlib.collections.PolyCollection(verts,edgecolors=('black',),facecolors='pink',linewidths=(0.7,))                                                                                   
            ax.add_collection(pc)
            plt.plot(X,Y,cor)
            for i in range(0,npoints):
             plt.text(X[i]+0.02,Y[i]+0.03,str(i),color='b')
                 
            plt.gca().set_aspect('equal')

            if mesh_quality==1: #visualiza qualidade de malha 
                if mesh=="tri":
                    crit_list=mesh_quality_tri(X,Y,IEN)
                if mesh=="quad":
                    crit_list=mesh_quality_quad(X,Y,IEN)
                pc.set_color(crit_list)
                complemento+=" com qualidade da malha"
            if normal==1 and mesh_quality==0: #visualiza e da informacoes da normal
                normal_info=normal_func(X,Y,IEN,mesh)
                pc.set_color(normal_info[1])
                complemento+=" com indicador de área"
                if normal_info[2]==1:
                    print ("-Orientacao da malha teve de ser corrigida")
                print ("-Area coberta pela malha de " + str(normal_info[0])+ " unidades")

            org_list=[] #plota por cima as linhas dos vertices
            for y in IEN:
                org_list.append([y[0],y[1]])
                org_list.append([y[1],y[2]])
                if mesh=="tri":
                    org_list.append([y[2],y[0]])
                if mesh=="quad":
                    org_list.append([y[2],y[3]])
                    org_list.append([y[3],y[0]])
            for i in org_list:
                xx=([X[i[0]],X[i[1]]])
                yy=([Y[i[0]],Y[i[1]]])   
                plt.plot(xx,yy,"-ko",linewidth=0.5,color="black")

            if borda==1: #plota por cima a borda
                IEN_bound=IEN_bound_gen(IEN,mesh)
                count=0
                complemento+=" e borda"
                while len(IEN_bound)-1>count: #plota a borda linha a linha
                    xx=([X[IEN_bound[count][0]],X[IEN_bound[count][1]]])
                    yy=([Y[IEN_bound[count][0]],Y[IEN_bound[count][1]]])
                    count+=1    
                    plt.plot(xx,yy,"-ko",color="red")
            plt.xlabel("X") ; plt.ylabel("Y") # adicionando legenda para eixo x e y
            plt.title("Malha" + complemento) # adicionando titulo ao plot
            plt.show()
    if borda==1: #so retorna ien_bound se pedir borda
            return ([IEN,IEN_bound,X,Y])
    return ([IEN,0,X,Y])



def mesh_3d(nx,ny,nz,xf=1,yf=1,zf=1,modex=1,modey=1,modez=1):
    lisx=[];lisy=[];lisz=[]
    X = np.linspace(0,xf**(1/modex),nx)
    X = (X**modex)
    Y = np.linspace(0,yf**(1/modey),ny)
    Y = (Y**modey)
    Z = np.linspace(0,zf**(1/modez),nz)
    Z = (Z**modez)
    for i in range(0,nz):#transformando em lista
        for i in range(0,ny): 
                for i in range(0,nx):
                    lisx.append(X[i])
    for i in range(0,nz):#transformando em lista
        for r in range(0,ny): 
                for i in range(0,nx):
                    lisy.append(Y[r])
    for r in range(0,nz):#transformando em lista
        for i in range(0,ny): 
                for i in range(0,nx):
                    lisz.append(Z[r])
    X=np.array(lisx);Y=np.array(lisy);Z=np.array(lisz)
    return [X,Y,Z]



#print (mesh_2d(3,3,1,1,1,0.5,1.57075,0.2,1,1,"quad",0,1,0,0,0,'malha-furo.msh'))





