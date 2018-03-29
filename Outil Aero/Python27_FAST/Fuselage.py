# -*- coding: utf-8 -*-
import numpy as np;
import math as m;
import utilitaire as u;
import matplotlib.pyplot as plt;
class Fuselage():
    """Classe définissant le fuselage dans son ensemble reprenant les diférentes
    charactéristiques de celles ci:
        nL : fuselage nose length
        cL : fuselage cylinder length
        bL : fuselage bottom length
        
        cD : fuselage cylinder diameter
        bD : fuselage bottom diameter
        
        hDist : horizontale distance between fuselage nose and wing root section C/4
        vDist : vertical distance between fuselage nose and wing root section C/4
            positive : high-wing.
        bool : un booléen qui dit si l'élément est présent sur l'avion
        yFus : y position of the fuselage-wing intersection
    """
    def __init__(self):
        self.nL = 0;
        self.cL = 0;
        self.bL = 0;
        self.cD = 0;
        self.bD = 0;
        self.hDist = 0;
        self.vDist = 0;
        self.yFus = np.zeros(2,dtype = float);
        self.lEq = 0.;
        self.bool = False;
        self.CorrFact = np.empty(0,dtype = float);
        self.generatrice = np.zeros((20,2),dtype = float);
        self.L = 0.;
        self.D = 0.;
        self.M = 0.;
        self.Y = 0.;
        self.N = 0.;
        
    
    def setNL(self,l):
        """ Définit la longueur du nez du fuselage"""
        self.nL = l;
    def getNL(self):
        """ Retourne la longueur du nez du fuselage"""
        return self.nL;
    
    def setCL(self,l):
        """ Définit la longueur du cylindre du fuselage"""
        self.bool = True;
        self.cL = l;
    def getCL(self):
        """ Retourne la longueur du cylindre du fuselage"""
        return self.cL;
    
    def setBL(self,l):
        """ Définit la longueur du culot du fuselage"""
        self.bL = l;
    def getBL(self):
        """ Retourne la longueur du culot du fuselage"""
        return self.bL;
    
    def setCD(self,d):
        """ Définit le diamètre du cylindre du fuselage"""
        self.cD = d;
    def getCD(self):
        """ Retourne le diamètre du cylindre du fuselage"""
        return self.cD;
    
    def setBD(self,d):
        """ Définit le diamètre du culot du fuselage"""
        self.bD = d;
    def getBD(self):
        """ Retourne le diamètre du culot du fuselage"""
        return self.bD;
    
    def setHD(self,d):
        """ Définit la distance horizontale entre le nez du fuselage
        et la wing root C/4"""
        self.hDist = d;
    def getHD(self):
        """ Définit la distance horizontale entre le nez du fuselage
        et la wing root C/4"""
        return self.hDist;
    
    def setVD(self,d):
        """ Définit la distance verticale entre le nez du fuselage
        et la wing root C/4"""
        self.vDist = d;
    def getVD(self):
        """ Définit la distance verticale entre le nez du fuselage
        et la wing root C/4"""
        return self.vDist;
    
    def setYF(self,yFus):
        """ Définit la position selon y de l'intersection fuselage-wing """
        self.yFus = yFus;
    def getYF(self):
        """ Retourne la position selon y de l'intersection fuselage-wing """
        return self.yFus;
    
    def setLEq(self,l):
        """ Définit la longueur equivatente telle que lEq*D_cl = int_0^L r(x) dx
        Où L est la ongueur totale du fuselage, r(x) la distance de la peau du fuselage
        à l'axe de révolution"""
        self.lEq = l;
    def getLEq(self):
        """ Retourne la longueur equivalente du fuselage"""
        return self.lEq;
    
    def setCF(self,l):
        """ Définit le facteur de correction pour l'angle d'attaque vu sur chaque panneau."""
        self.lEq = l;
    def getCF(self):
        """ Retourne le facteur de correction"""
        return self.lEq;
    
    def setGen(self,gen):
        """ Définit la génératrice du fuselage r(x) """
        self.generatrice = gen;
    
    def getGen(self):
        """ Retourne la génératrice du fuselage r(x) """
        return self.generatrice;
    
    def getL(self):
        """ Retourne le lift du fuselage"""
        return self.L
    def getD(self):
        """ Retourne le drag du fuselage"""
        return self.D
    def getM(self):
        """ Retourne le moment tangage du fuselage"""
        return self.M
    def getY(self):
        """ Retourne la force latérale du fuselage"""
        return self.Y
    def getN(self):
        """ Retourne le moment de lacet du fuselage"""
        return self.N
    

def FusManager(FI,W,vtail,flow,rp):
    F = Fuselage();
    F.setBD(FI.getBD());
    F.setBL(FI.getBL());
    F.setCD(FI.getCD());
    F.setCL(FI.getCL());
    F.setNL(FI.getNL());
    F.setVD(FI.getVD());
    F.setHD(FI.getHD());
    F.setYF(W.getYF());
    generatrice = np.zeros((20,2),dtype = float);
    Angle = np.linspace(0.,m.pi*0.5,18);
    for i in range(18):
        generatrice[i] = np.array([F.getNL()*(1.-m.cos(Angle[i])),0.5*F.getCD()*m.sin(Angle[i])]);
    generatrice[18] = generatrice[17] + np.array([F.getCL(),0.]);
    generatrice[19] = generatrice[18] + np.array([F.getBL(),0.5*(F.getBD()-F.getCD())]);
    F.setGen(generatrice);
    S = 2 * np.trapz(generatrice[:,1],generatrice[:,0]);
    F.setLEq(S/F.getCD());
    CF = np.ones(W.getR(),dtype = float);
    for i in range(W.getR()):
        if F.vDist > 0:
            if not (W.getYP(i) > F.getYF()[0] and W.getYP(i) < F.getYF()[1]):
                if not (W.getYP(i) > -F.getCD()*0.5 and W.getYP(i) < F.getCD()*0.5):
                    CF[i] += cFact(F.getCD()*0.5,W.getY(i),W.getY(i+1),W.getChordDist(i),W.getChordDist(i+1));
                else:
                    CF[i] += cFact2(F.getCD()*0.5,W.getY(i),W.getY(i+1),W.getChordDist(i),W.getChordDist(i+1));
    [F.L,F.D,F.M,F.Y,F.N] = aeroFus(F,flow.getAlphas(),flow.beta,rp,flow.V0);
    return F;

def cFact(rf,y1,y2,c1,c2):
    """ Source : propeller-wing interaction prediction for Early Design,
    Giovanna Ferraro, Timoleon Kipouros and A. Mark Savill"""
    y = np.linspace(y1,y2,20);
    c = np.linspace(c1,c2,20);
    return rf**2/(np.trapz(y,c*y**2))*(0.5*(y2-y1)*(c1+c2));

def cFact2(rf,y1,y2,c1,c2):
    """ Source : propeller-wing interaction prediction for Early Design, 
    Giovanna Ferraro, Timoleon Kipouros and A. Mark Savill"""
    y = np.linspace(y1,y2,20);
    c = np.linspace(c1,c2,20);
    return np.trapz(y,-c*(1.-2*(y/rf)**2))/(0.5*(y2-y1)*(c1+c2));

def aeroFus(F,alpha,beta,rp,V0=50.):
    """ Results based on slender bodies and Allen's theories"""
    a = alpha * m.pi/180;
    S_Cyl = 0.25 * m.pi * F.cD ** 2;
    f_F = (F.cL + F.bL + F.nL) / F.cD;
    FF = 1. + 2.2/(f_F ** (1.5)) - 0.9/(f_F ** (3.));
    gen = F.getGen();
    x = np.concatenate([np.linspace(0,gen[17,0],50), np.linspace(gen[17,0],gen[18,0],100),np.linspace(gen[18,0],gen[19,0],100)]);
    ReX = V0*x/(1.57e-5);
    delta = np.concatenate([[0.],4.92*x[1:50]/(ReX[1:50]**(0.5)),4.92*x[50]/(ReX[50]**(0.5))+0.35*(x[50:]-x[50])/(ReX[50:]**(0.2))]);
    rayon = np.interp(x,gen[:,0],gen[:,1])+delta;
    S = m.pi * rayon ** 2;
    cp_X = (x[1:]+x[:-1])*0.5;
    dS = S[1:]-S[:-1];
    
    CN_lin = np.sin(2.*a) * S[-1];
    CX_lin = np.sin(a)**2 * S[-1];
    L_lin = CN_lin * np.cos(a) - CX_lin * np.sin(a);
    D_lin = CN_lin * np.sin(a) + CX_lin * np.cos(a);
    M_lin = -np.sin(2.*a) * np.sum(dS*cp_X);
    ReF = V0 * cp_X /(1.57e-5);
    CF = np.concatenate([1.328/(ReF[:49]**0.5), 0.0442/(ReF[49:]**(1./6.))]);
    dX = x[1:]-x[:-1];
    r2 = np.interp(x,gen[:,0],gen[:,1]);
    ds = ((r2[1:]-r2[:-1])**2 + dX ** 2) ** (0.5);
    rMoy = (r2[1:]+r2[:-1])*0.5;
    dSW = 2*m.pi*ds*rMoy;
    Frot = np.sum(CF*dSW)*FF;
    S_Culot = 0.25 * m.pi*F.bD**2;
    D_Culot = 0.14*S_Culot; # Lecture Aerodynamics 2A : slender body
    if F.bL/F.cD < 0.8:
        D_con = 1.4*S_Cyl;
    else:
        D_con = 1.4 * m.exp(-((F.bL/F.cD)-0.8)*3./0.8) * S_Cyl;
    D_windscreen = S_Cyl * 2.e-3;
    D_par = (Frot+D_Culot + D_con + D_windscreen);
    
    Cx = 1.2; # Allen ensures that for transversal flow ok if no compressible effect
    L_visc = np.sin(a) * np.abs(np.sin(a)) * np.cos(a) * Cx * np.sum(2.*dX*rMoy);
    D_visc = np.abs(np.sin(a) ** 3) * Cx * np.sum(2.*dX*rMoy);
    M_visc = - Cx * np.abs(np.sin(a)) * np.sin(a) * np.sum(2.*dX*rMoy*cp_X);
    
    L = L_lin + L_visc - D_par * np.sin(a);
    D = D_lin + D_visc + D_par * np.cos(a);
    Moment = M_lin + M_visc;
    d = ((rp[0] + F.hDist) ** 2 + (rp[2]+F.vDist)**2)**(0.5);
    TP = np.arctan2(-(rp[2]+F.vDist),(rp[0]+F.hDist));
    M = Moment + d * L * np.cos(TP+a) + d * D * np.sin(TP + a);
    
    beta *= -m.pi/180;
    CN_lin = np.sin(2.*beta) * S[-1];
    CX_lin = np.sin(beta)**2 * S[-1];
    L_lin = CN_lin * np.cos(beta) - CX_lin * np.sin(beta);
    D_lin = CN_lin * np.sin(beta) + CX_lin * np.cos(beta);
    M_lin = -np.sin(2.*beta) * np.sum(dS*cp_X);
    L_visc = np.sin(beta) * np.abs(np.sin(beta)) * np.cos(beta) * Cx * np.sum(2.*dX*rMoy);
    D_visc = np.abs(np.sin(beta) ** 3) * Cx * np.sum(2.*dX*rMoy);
    M_visc = - Cx * np.abs(np.sin(beta)) * np.sin(beta) * np.sum(2.*dX*rMoy*cp_X);
    N = M_lin + M_visc;
    Y = L_lin + L_visc - D_par * np.sin(beta);
    Dt = D_lin + D_visc + D_par * np.cos(beta);
    dv = ((rp[0] + F.hDist) ** 2 + (rp[1])**2)**(0.5);
    TPv = np.arctan2(-(rp[1]),(rp[0]+F.hDist));
    N += dv * Y * m.cos(TPv - beta*m.pi/180) + dv * Dt * m.sin(TPv-beta*m.pi/180);
    return L,D,M,Y,N