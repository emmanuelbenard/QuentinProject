# -*- coding: utf-8 -*-
import numpy as np;
import math as m;
import matplotlib.pyplot as plt
class Propu():
    """ Classe qui représente le système propulsif de l'avion.
    Caractérisé au travers des arguments:
        - xp :np.array  x position des centres des moteurs, par rapport wing root C/4
            positive si en arrière
        - yp : np.array y position des centres des moteurs, par rapport wing root C/4
            positive si sur l'aile droite
        - zp : np.array z position des centres des moteurs, par rapport wing root C/4
            positive si au dessus
        - D : np.array qui donne les diamètres des hélices, en partant de l'extrème
            gauche vers l'extrème droite
        - rHub : np.array qui donne les rayons des centres d'hélice
        - Tc : np.array qui donne les coefficients de poussée de chaque moteurs;
            de l'extrême gauche vers droite. Tc_i = 2*T_i/(rho*V_0^2 * S_wing)
        - T np.array qui donne les poussées de chaque moteur (Newton);
            de l'extrême gauche vers la droite
        - Omega : np.array qui donne la vitesse de rotation (rad/sec) de chaque moteur;
            de l'extrême gauche vers la droite.
        - OWU : np.array qui donne le sens de rotation : si le moteur tourne de manière
            à ce que l'extérieur de l'aile voit un vent upward : OWU_i = 1
        - bool : booléen qui dit si l'élément est présent sur l'avion
    """
    def __init__(self):
        self.xp = np.empty(0,dtype = float);
        self.yp = np.empty(0,dtype = float);
        self.zp = np.empty(0,dtype = float);
        self.D = np.empty(0,dtype = float);
        self.rHub=  np.empty(0,dtype = float);
        self.Tc = np.empty(0,dtype = float);
        self.T = np.empty(0,dtype = float);
        self.Omega = np.empty(0,dtype = float);
        self.OWU = np.empty(0,dtype = bool);
        self.bool = False;
        self.vix = np.empty(0,dtype = float);
        self.viy = np.empty(0,dtype = float);
        self.viz = np.empty(0,dtype = float);
        
    def setXp(self,x):
        """ Falps deflection """
        self.xp = x;
    def getXp(self,ii='all'):
        """ Return the flaps deflection on the ii-th segment"""
        if ii == 'all':
            return self.xp[:];
        else:
            return self.xp[ii];

    def setYp(self,y):
        self.yp = y;
    def getYp(self,ii='all'):
        if ii == 'all':
            return self.yp[:];
        else:
            return self.yp[ii];

    def setZp(self,z):
        self.zp = z;
    def getZp(self,ii='all'):
        if ii == 'all':
            return self.zp[:];
        else:
            return self.zp[ii];

    def setD(self,d):
        """ Définit le diamètre des hélices"""
        self.D = d;
    def getD(self,ii='all'):
        if ii == 'all':
            return self.D;
        else:
            return self.D[ii];

    def setRh(self,r):
        """ Définit le rayon du capot central de l'hélice """
        self.rHub = r;
    def getRh(self,ii='all'):
        if ii == 'all':
            return self.rHub;
        else:
            return self.rHub[ii];

    def setTc(self,Tc):
        """ Définit le coefficient de traction des hélices """
        self.bool = True;
        self.Tc = Tc;
    def getTc(self,ii = 'all'):
        if ii == 'all':
            return self.Tc;
        else:
            return self.Tc[ii]
    def setT(self,T):
        """ Définit la traction de rotation des hélices """
        self.T = T;
    def getT(self,ii = 'all'):
        if ii == 'all':
            return self.T;
        else:
            return self.T[ii]
        
    def setOmega(self,Omega):
        """ Définit la vitesse de rotation des hélices """
        self.Omega = Omega;
    def getOmega(self,ii = 'all'):
        if ii == 'all':
            return self.Omega;
        else:
            return self.Omega[ii]
        
    def setOWU(self,OWU):
        """ Définit le sens de rotation des hélices OWU : 1 = Outboard Wind Up"""
        for i in range(len(OWU)):
            o = OWU[i];
            if o == 0:
                OWU[i] = 0;
            else:
                OWU[i] == 1;
        self.OWU = OWU;
    def getOWU(self,ii = 'all'):
        if ii == 'all':
            return self.OWU;
        else:
            return self.OWU[ii];
    
    def setV(self,vix,viy,viz):
        """ Définit les vitesses induites par l'hélice"""
        self.vix = np.concatenate([self.vix,vix]);
        self.viy = np.concatenate([self.viy,viy]);
        self.viz = np.concatenate([self.viz,viz]);
    def getV(self):
        """ Retourne un tableau reprennant les vitesses induites par l'hélice"""
        return self.vix,self.viy;self.viz;
        
def PropManager_SEI(pI,W,ht,vt,flow):
    """ Met en forme les données liées à la propulsion:
        - Format position hélices;
        - Vitesses induites sur l'aile et l'empennage horizontal"""
    p = Propu();
    p.setXp(pI.getXp());
    p.setYp(pI.getYp());
    p.setZp(pI.getZp());
    p.setD(pI.getD());
    p.setRh(pI.getRh());
    p.setTc(pI.getTc());
    for i in range(len(p.Tc)):
        if p.Tc[i] != 0 :
            p.Tc[i] -= np.sign(p.Tc[i]) * 0.15 * m.pi * p.rHub[i]**2/(W.getS());
        else:
            p.Tc[i] = -0.015;
    p.setT(pI.getT());
    p.setOmega(pI.getOmega());
    p.setOWU(pI.getOWU());
    vix,viy,viz = EIVelocities_SEI(p,W,flow);
    p.setV(vix,viy,viz);
    if ht.bool:
        vixt,viyt,vizt = EIVelocities_SEI(p,ht,flow)
        p.setV(vixt,viyt,vizt);
    return p;

def PropManager_DEI(pI,W,ht,vt,Vx,Vy,Vz):
    """ Met en forme les données liées à la propulsion:
        - Format position hélices;
        - Vitesses induites sur l'aile et l'empennage horizontal"""
    p = Propu();
    p.setXp(pI.getXp()-W.getChordDist(W.getR()/2)*0.25);
    p.setYp(pI.getYp());
    p.setZp(pI.getZp());
    p.setD(pI.getD());
    p.setRh(pI.getRh());
    p.setTc(pI.getTc());
    for i in range(len(p.Tc)):
        if p.Tc[i] == 0. :
            p.Tc[i] = -0.015;
    p.setT(pI.getT());
    p.setOmega(pI.getOmega());
    p.setOWU(pI.getOWU());
    # Function not implemented for computing the EIV: don't know how to compute incidence on prop
    # EIV should normally be computed at each correction interation, taking into account the 
    # current velocity triangle based on wingtip vortex induced velocity
    return p;
    
def EIVelocities_SEI(p,W,flow):
    """
        ATTENTION : LA FONCTION QUI PERMET DE CALCULER L'INFLUENCE DES MOTEURS
        EN PRESENCE D'UN ECOULEMENT NON UNIFORME N'EST PAS IMPLEMENTEE.
    """
    """ EIVelocities(p,W,flow) retourne la distribution de vitesses induites par le moteur sur
    les panneaux de W.
     p : les données du système propulsif
     W : les données de la surface portante
     flow : les données de l'écoulement
     Neidp is the number of elements in which the wing will be divided.
     V0: aircraft velocity                   [m/s]
     h: altitude study                       [Km]
     b: Semi-Span                         [m]
     l: Corde en chaque section           [m]
     Rhub: radius du hub                     [m]
     D: Diametre helice                      [m]
     Omega: vitesse de rotation de l'helice  [rad/s]
     T: Traction de l'helice             [N]
     yp: Position helices. (y=0 always)      [m]
         i.e.: 
               none: yp=0[]
               1 in the middle: yp=[0]
               2: yp=[-3 3]"""

    #Calcul masse volumique
    rho0=1.225; #masse volumique à niveau de la mer           [kg/m^3]
    dT=-6.5;    #gradiente de temperature dans la troposphere [K/km]
    T0=288.15;  #Temperature à niveau de la mer               [K]
    g=9.80665;  #gravité                                      [m/s^2]
    Rair=287.1;   #Constante de l'air                           [m^2/(s^2*K)]
    h = flow.getH();
    V0 = flow.getV0();
    rho = rho0 * (1. + dT*h/T0)**(- g/(Rair*dT*10**(-3)) - 1.);
    nbinter = 100;
    Neidp = nbinter*W.getR()+1;
    ##
    # Variables helice
    Sh = m.pi * p.getD()**2 *0.25;#Surface disque actuator          [m^2]
    N = len(p.D);
    ##
    # Discretitation 
    x = np.interp(np.linspace(0,W.getR()+1,Neidp),np.linspace(0,W.getR()+1,W.getR()+1),W.getX());
    y = np.interp(np.linspace(0,W.getR()+1,Neidp),np.linspace(0,W.getR()+1,W.getR()+1),W.getY());
    z = np.interp(np.linspace(0,W.getR()+1,Neidp),np.linspace(0,W.getR()+1,W.getR()+1),W.getZ());
    vix = np.zeros(W.getR(),dtype = float);
    viy = np.zeros(W.getR(),dtype = float);
    viz = np.zeros(W.getR(),dtype = float);
    dvix = np.zeros(Neidp,dtype = float);
    dviy = np.zeros(Neidp,dtype = float);
    dviz = np.zeros(Neidp,dtype = float);
    dvitheta = np.zeros(Neidp,dtype = float);
    signe = [-1., 1.];
    for i in range(Neidp):
        #Calcul de l'increment de vitesse, si l'on est dans la zone d'influence 
        #de l'helice.
        for j in range(N):
            d = m.sqrt((y[i] - p.yp[j])**2 + (z[i] - p.zp[j])**2);
            yd = abs((y[i] - p.yp[j]));
            if abs((p.zp[j]-W.vDist)) < p.rHub[j]:
                rP = abs(p.rHub[j]*m.cos(m.asin((p.zp[j]-W.vDist)/p.rHub[j])));
            else:
                rP = 0.;
            if abs((p.zp[j]-W.vDist)) < p.D[j]*0.5:
                D = abs(p.D[j]*m.cos(m.asin(2.*(p.zp[j]-W.vDist)/p.D[j])));
            else:
                D = 0;
            if ((yd >= rP) and (yd < D * 0.5) and p.Omega[j] != 0.):
                dvix[i] += 0.5*V0*(m.sqrt(1.+2.*p.T[j]/(rho*Sh[j]*V0**2))-1.);
                vix2 = 0.5*V0*(m.sqrt(1.+2.*p.T[j]/(rho*Sh[j]*V0**2))-1.);
                a = vix2/V0;
                aprim = 0.5 * (1. - m.sqrt(abs(1.-4.*a*(1.+a)*(V0/(p.Omega[j] * d))**2)));
                dvitheta[i] = abs((aprim * 2. * p.Omega[j] * d));
                angle = m.atan2(z[i] - p.zp[j], y[i] - p.yp[j]);
                dviy[i] -= signe[p.OWU[j]] * dvitheta[i] * m.sin(angle);
                dviz[i] += signe[p.OWU[j]] * dvitheta[i] * m.cos(angle);
    if N != 1:
        dviy[:int((Neidp-1.)/2.)] *= -1;
        dviz[:int((Neidp-1.)/2.)] *= -1;
    chord = np.interp(y,W.y,W.chordDistrib);
    ds = np.sqrt((x-W.getX(W.r/2))**2+(y-W.getY(W.r/2))**2+(z-W.getZ(W.r/2))**2);
    for i in range(W.getR()):
        debut = nbinter*i+1;
        fin = nbinter*(i+1);
        vix[i] = np.trapz(dvix[debut:fin]*chord[debut:fin],ds[debut:fin])/((ds[fin-1]-ds[debut])*(chord[debut]+chord[fin-1])*0.5);
        viy[i] = np.trapz(dviy[debut:fin]*chord[debut:fin],ds[debut:fin])/((ds[fin-1]-ds[debut])*(chord[debut]+chord[fin-1])*0.5);
        viz[i] = np.trapz(dviz[debut:fin]*chord[debut:fin],ds[debut:fin])/((ds[fin-1]-ds[debut])*(chord[debut]+chord[fin-1])*0.5);
    return vix,viy,viz;
        