import numpy as np

def findI(v,T,G):


    k=1.380649e-23
    G_ref=1000
    T_ref=298.13
    q=1.60217663e-19
    E_g=1.1
    k_i=0.005254

    i_sc=9.06
    v_oc=46.22
    N_s=72
    n=1.3
    
    i_rs=[]
    i_o=[]
    i_ph=[]
    i=[]
    for j in range(len(T)):
        t=T[j]
        g=G[j]
        i_rs.append(i_sc/(np.exp((q*v_oc)/(n*N_s*k*t))-1))
        i_o.append(i_rs[j]*((t/T_ref)**3)*np.exp(((q*E_g)/(n*k))*((1/T_ref)-(1/t))))
        i_ph.append((i_sc+k_i*(t-T_ref))*(g/G_ref))
   
        i.append(i_ph[j]-i_o[j]*(np.exp((q*(v)/(n*k*N_s*t)))-1))
                
    return i

def findP(i,v):
    p=0
    for j in range(len(i)):
        if p < i[j]*v[j]:
            p=i[j]*v[j]
    return p

def power(T,G):
    V=np.arange(-5,46.22,0.1)
    I=findI(V,T,G)
    P=[]
    for item in I:
        p=findP(item,V)
        P.append(p)

    return P