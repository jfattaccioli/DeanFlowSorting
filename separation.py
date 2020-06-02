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
tabr=[pr]
tabphi=[pphi]
tabrc=[rc]
tabr2=[pr2]
tab_vr=[0]
tab_vr2=[0]
tab_Fl=[0]
tab_Fd=[0]
tabr0=[pr0]

#Méthode en prenant un pas de temps arbitraire et en passant de l'accélération à la vitesse en utilisant ce pas
'''while pphi>0.001 : 
    vr= vr-Fl(rho,Umax,H,r,pr0)/m*dt+Fd(mu,r,rc,rho,Uf,H)/m*dt              #Utilisation du PFD sur la particule 1
    vr2=vr2-Fl(rho,Umax,H,r2,pr02)/m2*dt+Fd(mu,r2,rc,rho,Uf,H)/m2*dt        #Utilisation du PFD sur la particule 2
    pr0=pr0+vr*dt
    pr02=pr02+vr2*dt
    rc= (a)*((pphi*pphi+1)**(3/2))/(pphi*pphi+2)
    pr=pr0+rc
    pr2=pr02+rc
    pphi=pphi-0.27/rc*dt
    
    #Mise à jour des tableaux à chaque passage
    tabrc.append(rc)
    tabphi.append(pphi)
    tabr.append(pr)
    tabr2.append(pr2)
    tab_vr.append(vr)
    tab_vr2.append(vr2)
    tab_Fl.append(Fl(rho,Umax,H,r,pr0))
    tab_Fd.append(Fd(mu,r,rc,rho,Uf,H))'''



#Tracés des axes et de la position des particules en coordonnées cartésiennes 
'''axes=plt.gca()
axes.set_xlim(-0.0013,0.0013)
axes.set_ylim(-0.0013,0.0013)
x=tabr*np.cos(tabphi)
y=tabr*np.sin(tabphi)
x2=tabr2*np.cos(tabphi)
#plt.plot(tabT,tabr0)
plt.plot(x,y)
plt.plot(x2,y)'''


#Tracé la spirale d'archimède de 5 tours à bonne échelle
'''tabt=np.linspace(0,10*np.pi,10001)
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
plt.plot(tabxM,tabyM,'k')'''



#Méthode en utilisant odeint qui permet d'intégrer notre équation différentielle : pr''(t)+Fl-Fd=0 où pr est la position de la particule
'''pr''(t)+Fl-Fd=0
pr'(t)=p(t)
p'(t)=Fd-Fl'''
def integrale(y,t,rho,Umax,H,r,mu,rc,Uf):
    pr,p=y
    dydt=[p,Fd(mu,r,rc,rho,Uf,H)/m-Fl(rho,Umax,H,r,pr0)/m]
    return dydt
y0=[rc+pr0,0.0]
t=np.linspace(0,10*np.pi,1000)
sol=odeint(integrale,y0,t,args=(rho,Umax,H,r,mu,rc,Uf))
x=sol[:,0]*np.cos(t)
y=sol[:,0]*np.sin(t)

plt.plot(x,y)    


plt.show()

