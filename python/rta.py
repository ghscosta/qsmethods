#PACKAGES

import numpy as np
import pandas as pd

#INPUT PARAMETERS

file = 'rrnhub.dat'
lb = 0.18
trlx = 10000
tave = 10000

#RECONSTRUCT THE NETWORK
edges = pd.read_csv(file,header=None)

#IF VERTICES ARE LABELLED 0,1,...N, UNCOMMENT THE LINE BELOW
edges = edges - 1
############################################################

n_edges = len(edges)
n_vert = edges[0].unique()
n_vert = len(n_vert)

ini = np.zeros(n_vert).astype(int)
deg = np.zeros(n_vert).astype(int)
deg2= np.zeros(n_vert).astype(int)
adj = np.zeros(n_edges).astype(int)
n_boxes = 100
boxes = np.zeros((n_boxes,int(n_vert/n_boxes))).astype(int)
lab_boxes = np.zeros(n_vert).astype(int)

####BUILDING THE DEGREE MATRIX
for i in range(0,n_edges):
    deg[edges[0][i]] += 1
    deg2[edges[0][i]]+= 1

###############################

#########BUILDING THE INI MATRIX
ini[0] = 0
for i in range(1,n_vert):
    ini[i] = ini[i-1]+deg[i-1]
################################    
#########BUILDING THE ADJACENCY LIST
for i in range(0,n_edges):
    vt_in = edges[0][i]
    vt_out= edges[1][i]

    if(vt_out > vt_in):
        adj[ini[vt_in]+deg2[vt_in]-1] = vt_out
        adj[ini[vt_out]+deg2[vt_out]-1] = vt_in
        deg2[vt_in] -= 1
        deg2[vt_out]-= 1
    
deg_min = np.amin(deg)
deg_max = np.amax(deg)

for i in range (0,n_boxes):
		for j in range (0,int(n_vert/n_boxes)):
			boxes[j][i] = (n_vert/n_boxes)*(j) + i
			lab_boxes[int(n_vert/n_boxes)*(j)+i] = j
print(lab_boxes[0:200])


#EVOLUTION OF SIS MODEL

#############################INITIAL CONDITION

his = np.zeros(n_vert).astype(float)
sigma = np.zeros(n_vert).astype(int)
lis = np.zeros(n_vert).astype(int)

t = 0.0 ; tser = 0
rho = 0.0 ; rho2 = 0.0 ; 
cont = 0 ; edg_inf = 0.0 ;
sinft = 0 ;
    
for i in range(0,n_vert):
    sigma[i] = 1
    cont += 1
    lis[cont-1] = i
    edg_inf += deg[i]

#################################################

#IMPORTANT VARIABLES
t_i = np.zeros(n_vert).astype(int)
tn_i= np.zeros(n_vert).astype(int)
tbox= np.zeros(n_boxes).astype(int)
sinf= np.zeros(n_boxes).astype(int)
tmax= np.zeros(n_boxes).astype(int)

def qs_method_rta():
    global edg_inf
    global cont
    
    maxtbox = np.amax(sinf)
    propo = 1.0*sinft/total_t
    pdeci = propo - int(propo)
    n_reativ = int(propo)
    z1 = np.random.rand()

    if( z1 < pdeci):
        n_reativ = n_reativ + 1

    cont = 0

    for nativ in range (1,n_reativ+1):

        while True:
            z = np.random.randint(0,n_boxes)
            z1 = np.random.rand()
            if(z1 < 1.0*sinf[z]/maxtbox):
                break
                
        while True:
            z1 = np.random.randint(0,int(n_vert/n_boxes))
            z2 = np.random.rand()
            tgt = boxes[z][z1]
            if(z2 < 1.0*t_i[tgt]/tmax[z] and sigma[tgt] == 0):
                break
        sigma[tgt] = 1
        cont = cont+1
        lis[cont-1] = tgt
        edg_inf = edg_inf + deg[tgt]
        


#RELAXATION

while (t < trlx):
   
    p = 1.0*cont/(1.0*cont+1.0*lb*edg_inf)
    z = np.random.randint(0,cont)
    z1 = np.random.rand()

    if(z1 < p):
        sigma[lis[z]] = 0
        edg_inf -= deg[lis[z]]
        
        t_i[lis[z]] += 1
        box = lab_boxes[lis[z]]
        tmax[box] = max(tmax[box],t_i[lis[z]])
        sinft += 1
        sinf[box] += + 1
        
        lis[z] = lis[cont-1]
        cont -= 1
    
        if(cont == 0):
            qs_method_rta()
    else:
       
        while True:
            z = np.random.randint(0,cont)
            z2 = np.random.rand()
            if(z2 < 1.0*deg[lis[z]]/deg_max):
                break
        rnd = np.random.randint(0,deg[lis[z]])
        
        if(sigma[adj[ini[lis[z]]+rnd]] == 0):
            sigma[adj[ini[lis[z]]+rnd]] = 1
            cont += 1
            lis[cont-1] = adj[ini[lis[z]]+rnd]
            edg_inf += deg[lis[cont-1]]

    dt = 1.0*p/cont
    t = t+dt
    total_t = t
    if(t > tser):
        print("Time elapsed (Relaxation):", t,cont,t_i[0],t_i[1])
        tser = tser+trlx/10
        
t = 0 ; tser = 1

while (t < tave):
   
    p = 1.0*cont/(1.0*cont+1.0*lb*edg_inf)
    z = np.random.randint(0,cont)
    z1 = np.random.rand()
    
    if(z1 < p):
        sigma[lis[z]] = 0
        edg_inf -= deg[lis[z]]
        
        t_i[lis[z]] = t_i[lis[z]] + 1
        box = lab_boxes[lis[z]]
        tmax[box] = max(tmax[box],t_i[lis[z]])
        sinft = sinft+1
        sinf[box] = sinf[box] + 1
        
        lis[z] = lis[cont-1]
        cont -= 1
    
        if(cont == 0):
            qs_method_rta()
    else:
       
        while True:
            z = np.random.randint(0,cont)
            z2 = np.random.rand()
            if(z2 < 1.0*deg[lis[z]]/deg_max):
                break
        rnd = np.random.randint(0,deg[lis[z]])
        
        if(sigma[adj[ini[lis[z]]+rnd]] == 0):
            sigma[adj[ini[lis[z]]+rnd]] = 1
            cont += 1
            lis[cont-1] = adj[ini[lis[z]]+rnd]
            edg_inf += deg[lis[cont-1]]

    dt = 1.0*p/cont
    t = t+dt
    total_t = total_t + dt
    his[cont] = his[cont] + dt
    if(t > tser):
        print("Time elapsed (Average):", t,cont)
        tser = tser+tave/10
         
norm = sum(his)
his = 1.0*his/norm
print(his)  
for i in range(1,n_vert):
    rho =  rho + 1.0*i*his[i]
    rho2 = rho2 + 1.0*i*i*his[i]


tau = 1.0/his[1]
rho = 1.0*rho/n_vert
rho2 = 1.0*rho2/(1.0*n_vert)/(1.0*n_vert)
chi = 1.0*n_vert*(rho2-rho*rho)/rho

if(tau > 1e10):
    tau = t_ave

np.savetxt("pn.dat",his)
print("Density of infected vertices:",rho)
print("Dynamical susceptibility:",chi)
print("Lifespam:",tau)
