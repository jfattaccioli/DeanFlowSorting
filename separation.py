import matplotlib.pyplot as plt 
import numpy as np
import math
from scipy.integrate import odeint





H=50e-6                                          #Hauteur du tube
mu = 0.001 #
r=20e-6                                           #Taille particule 1
r2=15e-6                                          #Taille particule 2
a=250e-6/2/np.pi                                  #Espace entre chaque tour de la spirale
m=1050*4/3*np.pi*(r/2)**3                         #Masse particule
rho=1000
Uf=130e-3
Umax=2*Uf
m2=1050*4/3*np.pi*(r2/2)**3

pphi=10*np.pi                                     #Position initale de la particule selon phi
rc=(a)*((pphi*pphi+1)**(3/2))/(pphi*pphi+2)       #Rayon de courbure inital de la spirale
vphi=0.27                                         #Vitesse constante longitudinale
dt=0.00001                                        #Pas de temps
vr=0                                              #Vitesse latérale initiale de la particule 1
vr2=0                                             #Vitesse latérale initiale de la particule 2
pr0=4*H/5                                         #Position initale de la particule 1 au sein du tube
pr02=H/3                                          #Position initale de la particule 2 au sein du tube
pr=rc+pr0                                         #Position initale de la particule 1 
pr2=rc+pr02                                       #Position initale de la particule 2







#Net Lift Force
def Fl(rho,Umax,H,r,pr):
    if pr>H/5+0.00000001:
        Cl=np.abs(np.sin(10*np.pi/3/H*pr-5*np.pi/3))
    else :
        Cl =-np.abs(np.sin(10*np.pi/3/H*pr-5*np.pi/3))
    return rho*Umax*Umax*Cl*r*r*r*r/H/H



#Dean Force
def Fd(mu,r,Rc,rho,Uf,H):
    De=rho*Uf*H/mu*np.sqrt(H/(2*Rc))
    #De=0.94
    return 5.4e-4*np.pi*mu*De**(1.63)*r
    



#Première formulation du code - ! les tableaux ici sont des arrays
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


#Creation de tableaux des positions,vitesses,forces et rayons de courbures


#Méthode en prenant un pas de temps arbitraire et en passant de l'accélération à la vitesse en utilisant ce pas

def methodepas(pr,pphi,pr2,pr0,pr02,vr,vr2,rc):
    tabr=np.array([pr])
    tabphi=np.array([pphi])
    tabrc=np.array([rc])
    tabr2=np.array([pr2])
    tab_vr=np.array([0])
    tab_vr2=np.array([0])
    tab_Fl=np.array([0])
    tab_Fd=np.array([0])
    tabr0=np.array([pr0])
    while pphi>0.001 : 
        vr= vr-Fl(rho,Umax,H,r,pr0)/m*dt+Fd(mu,r,rc,rho,Uf,H)/m*dt              #Utilisation du PFD sur la particule 1
        vr2=vr2-Fl(rho,Umax,H,r2,pr02)/m2*dt+Fd(mu,r2,rc,rho,Uf,H)/m2*dt        #Utilisation du PFD sur la particule 2
        pr0=pr0+vr*dt
        pr02=pr02+vr2*dt
        rc= (a)*((pphi*pphi+1)**(3/2))/(pphi*pphi+2)
        pr=pr0+rc
        pr2=pr02+rc
        pphi=pphi-0.27/rc*dt
        
        #Mise à jour des tableaux à chaque passage
        tabrc=np.append(tabrc,rc)
        tabphi=np.append(tabphi,pphi)
        tabr=np.append(tabr,pr)
        tabr2=np.append(tabr2,pr2)
        tab_vr=np.append(tab_vr,vr)
        tab_vr2=np.append(tab_vr2,vr2)
        tab_Fl=np.append(tab_Fl,Fl(rho,Umax,H,r,pr0))
        tab_Fd=np.append(tab_Fd,Fd(mu,r,rc,rho,Uf,H))
    return(tabr,tabr2,tabphi)



#Tracé la spirale d'archimède de 5 tours à bonne échelle

def spirale():
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



#Méthode en utilisant odeint qui permet d'intégrer notre équation différentielle : pr''(t)+Fl-Fd=0 où pr est la position de la particule
'''pr''(t)+Fl-Fd=0
pr'(t)=p(t)
p'(t)=Fd-Fl'''
def integrale(y,t,rho,Umax,H,r,mu,Uf):
    rc=(a)*((t*t+1)**(3/2))/(t*t+2)
    pr,p=y
    dydt=[p,Fd(mu,r,rc,rho,Uf,H)-Fl(rho,Umax,H,r,pr0)]
    print(rc)
    print(dydt)
    return dydt
   




def main(n):
    if n==1 :
        spirale()
        tab1,tab2,tab3=methodepas(pr,pphi,pr2,pr0,pr02,vr,vr2,rc)
        axes=plt.gca()
        axes.set_xlim(-0.0013,0.0013)
        axes.set_ylim(-0.0013,0.0013)
        x=tab1*np.cos(tab3)
        y=tab1*np.sin(tab3)
        x2=tab2*np.cos(tab3)
        #plt.plot(tabT,tabr0)
        plt.plot(x,y)
    if n==2 : 
        t=np.linspace(10*np.pi,0,num=1000)
        y0=[rc+pr0,0.0]
        sol=odeint(integrale,y0,t,args=(rho,Umax,H,r,mu,Uf))
        x=sol[:,0]*np.cos(t)
        y=sol[:,0]*np.sin(t)
        plt.plot(x,y)
    plt.show()

main(2)

