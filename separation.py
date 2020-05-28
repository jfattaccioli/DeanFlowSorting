import matplotlib.pyplot as plt 
import numpy as np
import math


'''
Forces en présence :

net lift force : Fl=rho*(Umax/D)²*Cl*r⁴ 
    rho : densité du fluide (s u)
    Umax : vitesse maximum du fluide m/s = 2*Uf
    D : diametre hydraulique m
    r : rayon de la particule m
    Cl : lift coefiecient (0.5 s u)
'''


'''
Dean force : Fd = 5.4*10^(-4)*pi*mu*De^(1.63)*r
        De : nombre de Dean = Re*sqrt(D/(2*R))
    Re : nombre de Reynolds
    D : diametre hydraulique (m)
    R : rayon de courbure (m)
    r : rayon de la particule (m)
    mu : viscosité (kg/ms)
    
Constantes :
rho = 1000 kg/m³
mu =10^-3 kg/ms

Essais
De1=0.23  U=5  microL/min => Uf = 33 mm/s
De2=0.47  U=10 microL/min => Uf = 66 mm/s
De3=0.94  U=20 microL/min => Uf = 130 mm/s

r1 = 1.9 micrometre
r2 = 7.32 micrometre 

Tube : 
Large                       = 100 micrometre
Hauteur                     = 50  micrometre
Longeur                     = 13  cm
Rayon de courbure initial   = 3   mm
espace inter-boucle         = 250 micrometre
5 boucles


'''
def Fl(rho,Umax,H,r,pr):
    if pr>H/5+0.00000001:
        Cl=np.abs(np.sin(10*np.pi/3/H*pr-5*np.pi/3))
    else :
        Cl =-np.abs(np.sin(10*np.pi/3/H*pr-5*np.pi/3))
    return rho*Umax*Umax*Cl*r*r*r*r/H/H


'''H=50e-6
x=np.linspace(H/5,H,10000)
y=np.abs(np.sin(x))
y2=np.abs(np.sin(10*np.pi/3/H*x-5*np.pi/3))
plt.plot(x,y2)'''

def Fd(mu,r,Rc,rho,Uf,H):
    De=rho*Uf*H/mu*np.sqrt(H/(2*Rc))
    #De=0.94
    return 5.4e-4*np.pi*mu*De**(1.63)*r
    
'''t_min=0.00001
t_max=10*math.pi
n_t=10000
tabt=np.linspace(t_min,t_max,n_t)
tableposition=np.zeros(len(tabt))'''

'''tabrint=(a*tabt)
tabrext=tabrint+100e-6
tabxI=tabrint*np.cos(-tabt)
tabyI=tabrint*np.sin(-tabt)
tabxE=tabrext*np.cos(-tabt)
tabyE=tabrext*np.sin(-tabt)
tabrm=tabrext+tabrint/2
rc= (a*a+tabrm*tabrm)**(3/2)/(2*a*a+tabrm*tabrm)'''



'''for i in range(100):
    tab_r=np.append(tab_r,pr)
    tab_r2=np.append(tab_r2,pr2)
    tab_phi=np.append(tab_phi,pphi)
    #np.append(tab_phi2,pphi2)
    ar=(Fl(rho,Umax,H,r,pr)+Fd(mu,r,D,rc,rho,Uf,H))/m
    #ar2=Fl(rho,Umax,H,r2,m2,pr2)+Fd(mu,r2,D,rc,Re,m2)
    #print (Fl(rho,Umax,H,r,pr))      
    #print(Fd(mu,r,D,rc,rho,Uf,H))
    vr=ar*dt
    #vr2=ar2*dt
    pr=vr*dt+pr
    #print(pr)
    #pr2=vr2*dt+pr2
    pphi=vphi*dt+pphi
    #pphi2=r2*vphi*dt+pphi2
    rc= (a*a+pr*pr)**(3/2)/(2*a*a+pr*pr)
print(len(tab_r))'''

H=50e-6
mu = 0.001
r=20e-6
r2=15e-6
a=250e-6/2/np.pi



m=1050*4/3*np.pi*(r/2)**3
rho=1000
Uf=130e-3
Umax=2*Uf
m2=1050*4/3*np.pi*(r2/2)**3

pphi=10*np.pi
rc0=(a)*((pphi*pphi+1)**(3/2))/(pphi*pphi+2)
vphi=0.27
t=0
dt=0.00001
vr=0
vr2=0
pr0=4*H/5
pr02=H/3
rc=rc0
pr=rc+pr0
pr2=rc+pr02
tabr=[pr]
tabphi=[pphi]
tabrc=[rc0]
tabT=[0]
tabr2=[pr2]
tab_vr=[0]
tab_vr2=[0]
tab_Fl=[0]
tab_Fd=[0]
tabr0=[pr0]
while pphi>0.001 : 
    vr= vr-Fl(rho,Umax,H,r,pr0)/m*dt+Fd(mu,r,rc,rho,Uf,H)/m*dt
    vr2=vr2-Fl(rho,Umax,H,r2,pr02)/m2*dt+Fd(mu,r2,rc,rho,Uf,H)/m2*dt
    pr0=pr0+vr*dt
    pr02=pr02+vr2*dt
    rc= (a)*((pphi*pphi+1)**(3/2))/(pphi*pphi+2)
    pr=pr0+rc
    pr2=pr02+rc
    pphi=pphi-0.27/rc*dt
    #rc= (a*a+pr*pr)**(3/2)/(2*a*a+pr*pr)
    tabrc.append(rc)
    tabphi.append(pphi)
    tabr.append(pr)
    tabr2.append(pr2)
    t=t+dt
    tabT.append(t)
    tabr0.append(pr0)
    tab_vr.append(vr)
    tab_vr2.append(vr2)
    tab_Fl.append(Fl(rho,Umax,H,r,pr0))
    tab_Fd.append(Fd(mu,r,rc,rho,Uf,H))

#plt.plot(tabT,tab_Fd)
#plt.plot(tabphi,tabrc)
#plt.plot(tabT,tabrc)
axes=plt.gca()
axes.set_xlim(-0.0013,0.0013)
axes.set_ylim(-0.0013,0.0013)
x=tabr*np.cos(tabphi)
y=tabr*np.sin(tabphi)
x2=tabr2*np.cos(tabphi)
#plt.plot(tabT,tabr0)
plt.plot(x,y)
plt.plot(x2,y)
#plt.plot(tabphi,tabr)
#plt.plot(tabT,tab_vr,'b')
#plt.plot(tabT,tab_vr2,'r')

#plt.plot(tabphi,tabr)
#print (tab_r)
#print (tab_phi)

#x=tab_r*np.cos(tab_phi)
#y=tab_r*np.sin(tab_phi)
#plt.close()
#plt.plot(tab_phi,tab_r,'r')
#plt.plot(tab_phi,tab_r2)
#plt.xlabel('phi')
#plt.plot(x,y)
#plt.show()

'''tabrint=(a*tabt)
tabrext=tabrint+100e-6
tabxI=tabrint*np.cos(tabt)
tabyI=tabrint*np.sin(tabt)
tabxE=tabrext*np.cos(tabt)
tabyE=tabrext*np.sin(tabt)
tabrm=tabrext+tabrint/2
plt.plot(tabxI,tabyI)
plt.plot(tabxE,tabyE)
plt.show()'''
    
    
    


tabt=np.linspace(0,10*np.pi,10001)
tabrint=(a*tabt)
tabrext=tabrint+50e-6
tabmil=tabrint+25e-6
tabxI=tabrint*np.cos(tabt)
tabyI=tabrint*np.sin(tabt)
tabxE=tabrext*np.cos(tabt)
tabyE=tabrext*np.sin(tabt)
tabxM=tabmil*np.cos(tabt)
tabyM=tabmil*np.sin(tabt)
plt.plot(tabxI,tabyI,'k')
plt.plot(tabxE,tabyE,'k')
plt.plot(tabxM,tabyM,'k')

#print(rc)
#plt.plot(tabr,rc)


'''M=22e-6
M2=20e-6
vitr=0
vitr2=0
vity=0.27
y=0
dt=0.00001
rc= (a*a+M*M)**(3/2)/(2*a*a+M*M)
tabM=[M]
tabM2=[M2]
taby=[y]
while y<0.0013 :
    vitr=vitr+(Fl(rho,Umax,H,r,M)+Fd(mu,r,H,rc,rho,Uf,H))/m*dt
    vitr2=vitr2+(Fl(rho,Umax,H,r2,M2)+Fd(mu,r2,H,rc,rho,Uf,H))/m2*dt
    print(Fl(rho,Umax,H,r2,M2))
    print('Fd=',Fd(mu,r2,H,rc,rho,Uf,H))
    M=M+vitr*dt
    M2=M2+vitr2*dt
    y=y+vity*dt
    rc= (a*a+M*M)**(3/2)/(2*a*a+M*M)
    tabM.append(M)
    taby.append(y)
    tabM2.append(M2)
axes=plt.gca()
axes.set_ylim(0,0.000050)  
plt.plot(taby,tabM)
plt.plot(taby,tabM2)'''
#ar=Fd(0.001,1.9e-6,150e-6,rc,0.47,m)#+Fl(1000,260e-3,150e-6,0.5,7e-6,m)
#print (ar*m)
#vr=ar*tabt+0
#pr=vr*tabt
#plt.plot(tabt,pr)
plt.show()

