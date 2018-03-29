# -*- coding: utf-8 -*-
import numpy as np;
import math as m;
class Wing():
    """ Classe définissant une aile dans son ensemble reprenant les diférentes
    charactéristiques de celles ci:

        b : span
        chord : vertical array with the chord on any pannel;
        flapsDiscY : vertical array with the spanwise coordinate of the flaps
           discontinuity sections
        afDiscY : vertical array with the spanwise coordinate of the
           airfoil discontinuity sections
        polarDiscY : vertical array with the spanwise coordinate of the
           polar discontinuity sections (2D curves)
        taperDiscY : vertical array with the spanwise coordinate of the
           taper ratio discontinuity sections
        sweepDiscY : vertical array with the spanwise coordinate of the
           sweep angle discontinuity sections
        dihDiscY : vertical array with the spanwise coordinate of the
           dihedral angle discontinuity sections
        twistDiscY : vertical array with the spanwise coordinate of the
           twist rate discontinuity sections
        FTC : the ratio between the flaps length and the local chord
           length
        polar : a cell-array with for each cell a "polar" object with field
           data (the alpha-cl-cd data), clmax, alphal0 (zero-lift angle).
        airfoil : a cell-array with each cell gives the airfoil naca number
           representation, cell 1 correspond to first panel after root.
        sweepC4 : vertical array with wing.sweep(i) is the sweep angle of 
           the panel from wing.y(i) to wing.y(i+1) at c/4 (rad)
        dih : vertical array with wing.dih(i) is the dihedral angle of 
           the panel from wing.y(i) to wing.y(i+1) (rad)
        twist : vertical array with wing.twist(i) is the twist angle of 
           the section at wing.y(i) (rad)
        DF : vertical array with the different flaps defection
           along the span (deg)
        DFR : Flaps deflection along the right hand side wing
        DFL : flaps deflection along the LHS wing.
        r : number of spanwise panel along the wing;
        p : number of chordwise panel along the airfoil;
        vix : vertical array with wing.vix(i) is the mean (surface integrated
           normalized by the surface of the panel) of the tangencial engine 
           induced velocity on the panel from wing.y(i) to wing.y(i+1)
        viy : vertical array with wing.viy(i) is the mean (surface integrated
           normalized by the surface of the panel) of the lateral engine 
           induced velocity on the panel from wing.y(i) to wing.y(i+1)
        viz : vertical array with wing.viz(i) is the mean (surface integrated
           normalized by the surface of the panel) of the normal engine 
           induced velocity on the panel from wing.y(i) to wing.y(i+1)
        theta : the angle (0 -> pi) spanwise location of the limits of the panels
        eta : the normalized (-1 -> 1) spanwise location of the limits of the panels
        y : the  spanwise location of (-b/2 -> b/2) the limits of the panels
        discY : vertical array of the complete set of the spanwise location
           of the discontinuity
        discTheta : same as above but in angle coordinates;
        airfoilIndex : vertical array with wing.airfoilIndex(i) is the index of 
           the airfoil (wing.airfoil) to use for the section at wing.y(i)
        polarIndex : vertical array with wing.polarIndex(i) is the index of 
           the airfoil (wing.polar) to use for the section at wing.y(i)
        chordDistrib : vertical array with wing.chordDistrib(i) is the chord length of 
           the section at wing.y(i)
        bool : un booléen qui dit si l'élément est présent sur l'avion
     """
    def __init__(self):
        """ Définition des attributs pour le premier constructeur """
        # general properties
        self.bool = True;
        self.b = 0.;
        self.S = 0.;
        self.hDist = 0.; #Only for htail
        self.vDist = 0.;
        self.r = 40;
        self.mC = 8;
        self.x = np.empty(0,dtype = float);
        self.y = np.empty(0,dtype = float);
        self.z = np.empty(0,dtype = float);
        self.yf = 0;
        self.xP = np.empty(0,dtype = float);
        self.yP = np.empty(0,dtype = float);
        self.zP = np.empty(0,dtype = float);
        self.theta = np.empty(0,dtype = float);
        self.eta = np.empty(0,dtype = float);
        self.FTC = 0.;
        self.MAC = 0.;

        # Geometric distribution along the span
        self.chordDistrib = np.empty(0,dtype = float);
        self.sweepc4 = np.array(0,dtype = float);
        self.dih = np.array(0,dtype = float);
        self.twist = np.array(0,dtype = float);
        self.airfoilIndex = np.empty(0,dtype=int);
        self.polarIndex = np.empty(0,dtype=int);
        self.DF = np.array(0,dtype=float);
        self.chord = np.empty(0,dtype=float);

        # Discontinuity properties
        self.flapsDiscY = np.empty(0,dtype = float);
        self.afDiscY = np.empty(0,dtype = float);
        self.polarDiscY = np.empty(0,dtype = float);
        self.taperDiscY = np.empty(0,dtype = float);
        self.sweepDiscY = np.empty(0,dtype = float);
        self.dihDiscY = np.empty(0,dtype = float);
        self.twistDiscY = np.empty(0,dtype = float);
        self.discTheta = np.empty(0,dtype = float);
        self.discThetaPol = np.empty(0,dtype = float);

        self.polarR = np.empty(0,dtype = object);
        self.polarL = np.empty(0,dtype = object);
        self.airfoil = np.empty(0,dtype = str);
        self.DFL = np.empty(0,dtype = float);
        self.DFR = np.empty(0,dtype = float);
        
        # Engines induced velocities
        self.vix = np.zeros(0,dtype = float);
        self.viy = np.zeros(0,dtype = float);
        self.viz = np.zeros(0,dtype = float);
        

        

    def incSpan(self,b):
        """ Incrémente le span de l'aile de b"""
        self.b += b;
    def decSpan(self,b):
        """ Diminue le span de l'aile de b"""
        self.b -= b;
    def getSpan(self):
        """ Methode qui retourne le span"""
        return self.b;

    def setR(self,r):
        """ Définit le nombre de panneaux le long du span """
        self.bool = True;
        self.r = r;
    def getR(self):
        """ Retourne le nombre de panneaux le long du span """
        return self.r;

    def getmC(self):
        """ Retourne le nombre de point de controle le long de la corde """
        return self.mC;
    
    def setX(self,x):
        """ Définit la position (x-comp) des quarts de corde des 
        sections de controle le long du span, negative : left wing """
        self.x = x;
    def getX(self,ii='all'):
        """ Retourne la position (x-comp) des quarts de corde des 
        sections de controle le long du span, negative : left wing """
        if ii == 'all':
            return self.x;
        else:
            return self.x[ii];
    
    def setY(self,y):
        """ Définit la position (y-comp) des quarts de corde des 
        sections de controle le long du span, negative : left wing """
        self.y = y;
    def getY(self,ii='all'):
        """ Retourne la position (y-comp) des quarts de corde des 
        sections de controle le long du span, negative : left wing """
        if ii == 'all':
            return self.y;
        else:
            return self.y[ii];

    def setZ(self,z):
        """ Définit la position (z-comp) des quarts de corde des 
        sections de controle le long du span, negative : left wing """
        self.z = z;
    def getZ(self,ii='all'):
        """ Retourne la position (z-comp) des quarts de corde des 
        sections de controle le long du span, negative : left wing """
        if ii == 'all':
            return self.z;
        else:
            return self.z[ii];
    
    def getYF(self,ii='all'):
        """ Retourne la position (y-comp) de l'intersection Wing-Fuselage """
        if ii == 'all':
            return self.yf;
        else:
            return self.yf[ii];
        
    def setXP(self,x):
        """ Définit la position (x-comp) des quarts de corde des 
        panneaux le long du span, negative : left wing """
        self.xP = x;
    def getXP(self,ii='all'):
        """ Retourne la position (x-comp) des quarts de corde des 
        panneaux le long du span, negative : left wing """
        if ii == 'all':
            return self.xP;
        else:
            return self.xP[ii];
    
    def setYP(self,y):
        """ Définit la position (y-comp) des quarts de corde des 
        panneaux le long du span, negative : left wing """
        self.yP = y;
    def getYP(self,ii='all'):
        """ Retourne la position (y-comp) des quarts de corde des 
        panneaux le long du span, negative : left wing """
        if ii == 'all':
            return self.yP;
        else:
            return self.yP[ii];

    def setZP(self,z):
        """ Définit la position (z-comp) des quarts de corde des 
        panneaux le long du span, negative : left wing """
        self.zP = z;
    def getZP(self,ii='all'):
        """ Retourne la position (z-comp) des quarts de corde des 
        panneaux le long du span, negative : left wing """
        if ii == 'all':
            return self.zP;
        else:
            return self.zP[ii];
    
    def setTheta(self,theta):
        """ Définit la position (y-comp) en radians telque theta = 0 : right tip;
        theta = pi : left tip;
        des quarts de corde des sections de controle le long du span"""
        self.theta = theta;
    def getTheta(self,ii='all'):
        """ Retourne la position (y-comp) en radians telque theta = 0 : right tip;
        theta = pi : left tip;
        des quarts de corde des sections de controle le long du span"""
        if ii == 'all':
            return self.theta;
        else:
            return self.theta[ii];
        
    def setEta(self,eta):
        """ Définit la position (y-comp) normalisée des quarts de corde des 
        sections de controle le long du span, -1 : left left tip, 1 right tip """
        self.eta = eta;
    def getEta(self,ii='all'):
        """ Retourne la position (y-comp) normalisée des quarts de corde des 
        sections de controle le long du span, -1 : left left tip, 1 right tip """
        if ii == 'all':
            return self.eta;
        else:
            return self.eta[ii];
    

    def setCF(self,xv):
        """ Définit le ratio longueur du volet sur longueur de la corde locale"""
        self.FTC = xv;
    def getCF(self):
        """ Retourne le ratio longueur du volet sur longueur de la corde locale"""
        return self.FTC;


    def addChordDist(self,c):
        """ Définit la longueur de corde des sections de controle"""
        self.chordDistrib = np.append(self.chordDistrib,c);
    def getChordDist(self,ii='all'):
        """ Retourne la longueur de corde des sections de controle"""
        if ii == 'all':
            return self.chordDistrib[:];
        else:
            return self.chordDistrib[ii];
    def resetSweepC4(self,phi):
        """ Redéfinit l'angle de flèche à l'emplanture"""
        self.sweepc4 = np.empty(0,dtype = float);
        self.sweepc4 = np.append(self.sweepc4,phi);
    def addSweepC4(self,phi):
        """ Ajoute l'angle de flèche au quart de corde pour chaque panneau"""
        self.sweepc4 = np.append(self.sweepc4,phi);
    def setSweepC4(self,phi):
        """ Déinit l'angle de flèche au quart de corde pour chaque panneau"""
        self.sweepc4 = phi;
    def getSweepC4(self,ii='all'):
        """ Retourne l'angle de flèche au quart de corde pour chaque panneau"""
        if ii == 'all':
            return self.sweepc4;
        else:
            return self.sweepc4[ii];

    def resetDih(self,psi):
        """ Redéfinit l'angle de dihèdre à l'emplanture"""
        self.dih = np.append(np.empty(0,dtype = float),psi);
    def addDih(self,psi):
        """ Ajoute l'angle de dihèdre au quart de corde pour chaque panneau"""
        self.dih = np.append(self.dih,psi);
    def setDih(self,phi):
        """ Déinit l'angle de dihèdre au quart de corde pour chaque panneau"""
        self.dih = phi;
    def getDih(self,ii='all'):
        """ Retourne l'angle de dihèdre au quart de corde pour chaque panneau"""
        if ii == 'all':
            return self.dih[:];
        else:
            return self.dih[ii];
    
    def resetTwist(self,a):
        """ Redéfinit l'angle de twist à l'emplanture = incidence aile"""
        self.twist = np.append(np.empty(0,dtype = float),a);
    def addTwist(self,theta):
        """ Ajoute l'angle de twist au quart de corde pour chaque panneau"""
        self.twist = np.append(self.twist,theta);
    def setTwist(self,phi):
        """ Déinit l'angle de twist au quart de corde pour chaque panneau"""
        self.twist = phi;
    def getTwist(self,ii='all'):
        """ Retourne l'angle de twist au quart de corde pour chaque panneau"""
        if ii == 'all':
            return self.twist[:];
        else:
            return self.twist[ii];

    def addAFI(self,index):
        """ Ajoute l'index référant le profil utilisé sur chaque panneau"""
        self.airfoilIndex = np.append(self.airfoilIndex,index);
    def getAFI(self,ii='all'):
        """ Retourne l'index référant le profil utilisé sur chaque panneau"""
        if ii == 'all':
            return self.airfoilIndex;
        else:
            return self.airfoilIndex[ii];

    def addPolI(self,index):
        """ Ajoute l'index référant le fichier polaire utilisé sur chaque panneau"""
        self.polarIndex = np.append(self.polarIndex,index);
    def getPolI(self,ii='all'):
        """ Retourne l'index référant le fichier polaire utilisé sur chaque panneau"""
        if ii == 'all':
            return self.polarIndex;
        else:
            return self.polarIndex[ii];

    def addDsF(self,delta):
        """ Ajoute la déflection du volet sur chaque panneau """
        self.deltasFlaps = np.append(self.deltasFlaps,delta);
    def getDsF(self,ii='all'):
        """ Retourne la déflection du volet sur chaque panneau """
        if ii == 'all':
            return self.deltasFlaps;
        else:
            return self.deltasFlaps[ii];

    def addChord(self,c):
        """ Ajoute la corde de chaque panneau """
        self.chord = np.append(self.chord,c);
    def getChord(self,ii='all'):
        """ Retourne la corde de chaque panneau"""
        if ii == 'all':
            return self.chord;
        else:
            return self.chord[ii];




    def addFlapsD(self,y):
        """ Ajoute une position de discontinuité de déflection de volet
        le long du span de l'aile"""
        self.flapsDiscY = np.append(self.flapsDiscY,y);
    def setFlapsD(self,y):
        """ Définit les positions de discontinuité de déflection de volet
        le long du span de l'aile"""
        self.flapsDiscY = y;
    def rmFlapsD(self):
        """ Enlève la dernière position de discontinuité de déflection de volet
        le long du span de l'aile"""
        self.flapsDiscY = self.flapsDiscY[:-1];
    def getFlapsD(self,ii='all'):
        """ Retourne les positions de discontinuité de déflection de volet
        le long du span de l'aile"""
        if ii =='all':
            return self.flapsDiscY;
        else:
            return self.flapsDiscY[ii];

    def addAfD(self,y):
        """ Ajoute une position de discontinuité de profil
        le long du span de l'aile"""
        self.afDiscY = np.append(self.afDiscY,y);
    def setAfD(self,y):
        """ Définit les positions de discontinuité de profil
        le long du span de l'aile"""
        self.afDiscY = y;
    def rmAfD(self):
        """ Supprime la dernière position de discontinuité de profil
        le long du span de l'aile"""
        self.afDiscY = self.afDiscY[:-1];
    def getAfD(self,ii='all'):
        """ Retourne les positions de discontinuité de profil
        le long du span de l'aile"""
        if ii =='all':
            return self.afDiscY;
        else:
            return self.afDiscY[ii];

    def addPolD(self,y):
        """ Ajoute une position de discontinuité de fichier polar
        le long du span de l'aile"""
        self.polarDiscY = np.append(self.polarDiscY,y);
    def setPolD(self,y):
        """ Définit les positions de discontinuité de fichier polar
        le long du span de l'aile"""
        self.polarDiscY = y;
    def rmPolD(self):
        """ Supprime la dernière position de discontinuité de fichier polar
        le long du span de l'aile"""
        self.polarDiscY = self.polarDiscY[:-1];
    def getPolD(self,ii='all'):
        """ Retourne les positions de discontinuité de fichier polar
        le long du span de l'aile"""
        if ii =='all':
            return self.polarDiscY;
        else:
            return self.polarDiscY[ii];

    def addTapD(self,y):
        """ Ajoute une position de discontinuité de taper
        le long du span de l'aile"""
        self.taperDiscY = np.append(self.taperDiscY,y);
    def setTapD(self,y):
        """ Définit les positions de discontinuité de taper
        le long du span de l'aile"""
        self.taperDiscY = y;
    def rmTapD(self):
        """ Supprime la dernère position de discontinuité de taper
        le long du span de l'aile"""
        self.taperDiscY = self.taperDiscY[:-1];
    def getTapD(self,ii='all'):
        """ Retourne les positions de discontinuité de taper
        le long du span de l'aile"""
        if ii =='all':
            return self.taperDiscY;
        else:
            return self.taperDiscY[ii];

    def addSweepD(self,y):
        """ Ajoute une position de discontinuité d'angle de flèche
        le long du span de l'aile"""
        self.sweepDiscY = np.append(self.sweepDiscY,y);
    def setSweepD(self,y):
        """ Définit les positions de discontinuité d'angle de flèche
        le long du span de l'aile"""
        self.sweepDiscY = y;
    def rmSweepD(self):
        """ Supprime la dernière position de discontinuité d'angle de flèche
        le long du span de l'aile"""
        self.sweepDiscY = self.sweepDiscY[:-1];
    def getSweepD(self,ii='all'):
        """ Retourne les positions de discontinuité d'angle de flèche
        le long du span de l'aile"""
        if ii =='all':
            return self.sweepDiscY;
        else:
            return self.sweepDiscY[ii];

    def addDihD(self,y):
        """ Ajoute une position de discontinuité d'angle de dihèdre
        le long du span de l'aile"""
        self.dihDiscY = np.append(self.dihDiscY,y);
    def setDihD(self,y):
        """ Définit les positions de discontinuité d'angle de dihèdre
        le long du span de l'aile"""
        self.dihDiscY = y;
    def rmDihD(self):
        """ Supprime la dernière position de discontinuité d'angle de dihèdre
        le long du span de l'aile"""
        self.dihDiscY = self.dihDiscY[:-1];
    def getDihD(self,ii='all'):
        """ Retourne les positions de discontinuité d'angle de dihèdre
        le long du span de l'aile"""
        if ii =='all':
            return self.dihDiscY;
        else:
            return self.dihDiscY[ii];

    def addTwistD(self,y):
        """ Ajoute une position de discontinuité d'angle de twist
        le long du span de l'aile"""
        self.twistDiscY = np.append(self.twistDiscY,y);
    def setTwistD(self,y):
        """ Définit les positions de discontinuité d'angle de twist
        le long du span de l'aile"""
        self.twistDiscY = y;
    def rmTwistD(self):
        """ Supprime la dernière position de discontinuité d'angle de twist
        le long du span de l'aile"""
        self.twistDiscY = self.twistDiscY[:-1];
    def getTwistD(self,ii='all'):
        """ Retourne les positions de discontinuité d'angle de twist
        le long du span de l'aile"""
        if ii =='all':
            return self.twistDiscY;
        else:
            return self.twistDiscY[ii];

    def setDiscTheta(self,theta):
        """ Définit les positions de discontinuité en coordonnées theta
        le long du span de l'aile"""
        self.discTheta = theta;
    def getThetaD(self):
        """ Retourne les positions de discontinuité en coordonnées theta
        le long du span de l'aile"""
        return self.discTheta;

    def setDiscThetaPol(self,theta):
        """ Définit les positions de discontinuité de polair en coordonnées theta
        le long du span de l'aile"""
        self.discThetaPol = theta;
    def getThetaDPol(self):
        """ Retourne les positions de discontinuité de polair en coordonnées theta
        le long du span de l'aile"""
        return self.discThetaPol;

    def addPolarR(self,polar):
        """ Ajoute un object polar """
        self.polarR = np.append(self.polarR,polar);
    def addPolarL(self,polar):
        """ Ajoute un object polar """
        self.polarL = np.append(self.polarL,polar);
    def getPolar(self,side,ii='all'):
        """ Retourne les objects polar """
        if side == 0:
            if ii == 'all':
                return self.polarL;
            else:
                return self.polarL[ii];
        else:
            if ii == 'all':
                return self.polarR;
            else:
                return self.polarR[ii];

    def addAF(self,af):
        """ Ajoute un profil"""
        self.airfoil = np.append(self.airfoil,af);
    def getAF(self,ii='all'):
        """ Retourne les profils"""
        if ii == 'all':
            return self.airfoil[:];
        else:
            return self.airfoil[ii];

    def addDFL(self,delta):
        """ Ajoute une déflection de volet pour l'aile gauche """
        self.DFL = np.append(self.DFL,delta);
    def getDFL(self,ii='all'):
        """ Retourne les déflections de volet pour l'aile gauche"""
        if ii == 'all':
            return self.DFL[:];
        else:
            return self.DFL[ii];
    
    def addDFR(self,delta):
        """ Ajoute une déflection de volet pour l'aile droite """
        self.DFR = np.append(self.DFR,delta);
    def getDFR(self,ii='all'):
        """ Retourne les déflections de volet pour l'aile droite """
        if ii == 'all':
            return self.DFR[:];
        else:
            return self.DFR[ii];
    def getDF(self,i):
        if i < self.r/2:
            return self.DFL[self.r/2 - i -1];
        else:
            return self.DFR[i - self.r/2];

    def setViX(self,vix):
        """ Définit les vitesses induites par le moteur
        selon l'axe longitudinal sur chaque panneau """
        self.vix = vix;
    def getViX(self,ii='all'):
        """ Retourne les vitesses induites par le moteur
        selon l'axe longitudinal sur chaque panneau """
        if ii == 'all':
            return self.vix;
        else:
            return self.vix[ii];

    def setViY(self,viy):
        """ Définit les vitesses induites par le moteur
        selon l'axe longitudinal sur chaque panneau """
        self.viy = viy;
    def getViY(self,ii='all'):
        """ Retourne les vitesses induites par le moteur
        selon l'axe longitudinal sur chaque panneau """
        if ii == 'all':
            return self.viy;
        else:
            return self.viy[ii];
        
    def setViZ(self,viz):
        """ Définit les vitesses induites par le moteur
        selon l'axe longitudinal sur chaque panneau """
        self.viz = viz;
    def getViZ(self,ii='all'):
        """ Retourne les vitesses induites par le moteur
        selon l'axe longitudinal sur chaque panneau """
        if ii == 'all':
            return self.viz;
        else:
            return self.viz[ii];
    def getMac(self):
        y = self.y;
        c = self.chordDistrib;
        dA = 0.5*(y[1:]-y[:-1])*(c[1:]+c[:-1]);
        S = sum(dA);
        cMean = 0.5*(c[1:]+c[:-1]);
        return np.sum(cMean*dA)/S;
    
    def getSI(self):
        Y = np.concatenate([[0.],self.taperDiscY,[self.b*0.5]]);
        dY = Y[1:]-Y[0:-1];
        c = (self.chord[1:]+self.chord[0:-1]);
        return np.sum(dY*c);
    
    def getS(self):
        dY = self.y[1:]-self.y[0:-1];
        c = (self.chordDistrib[1:]+self.chordDistrib[0:-1])*0.5;
        return np.sum(dY*c);
    
def WingManager(WI,prop,fus):
    W = Wing();
    for i in range(len(WI.getPolar(0))):
        W.addPolarL(WI.getPolar(0,i));
        W.addPolarR(WI.getPolar(1,i))
    W.addAF(WI.getAF());
    W.S = WI.getSI();
    r = WI.r;
    #theta = np.linspace(0.,m.pi,W.getR());
    #W.setTheta(theta[::-1]);
    #W.setEta(np.cos(W.getTheta()));
    #theta = np.linspace(-1.,1.,W.getR());
    W.incSpan(WI.getSpan());
    #W.setY(theta*W.getSpan()*0.5);
    W.addTapD(np.sort(np.concatenate([WI.getTapD(),-WI.getTapD()])));
    W.addDihD(np.sort(np.concatenate([WI.getDihD(),-WI.getDihD()])));
    W.addSweepD(np.sort(np.concatenate([WI.getSweepD(),-WI.getSweepD()])));
    W.addTwistD(np.sort(np.concatenate([WI.getTwistD(),-WI.getTwistD()])));
    W.addAfD(np.sort(np.concatenate([WI.getAfD(),-WI.getAfD()])));
    W.addFlapsD(np.sort(np.concatenate([WI.getFlapsD(),-WI.getFlapsD()])));
    W.addPolD(np.sort(np.concatenate([WI.getPolD(),-WI.getPolD()])));
    YEng = getYEng(WI,prop);
    YFus = getYFus(WI,fus);
    W.yf = YFus;
    YTip = W.getSpan()* np.array([-0.5,0.,0.5]);
    discY = np.sort(np.unique(np.concatenate([W.getTapD(),W.getDihD(),W.getSweepD(),\
        W.getTwistD(),W.getAfD(),W.getFlapsD(),W.getPolD(),YEng,YFus,YTip])));
    dY = discY[1:] - discY[:-1];
    nbPan = np.round(r * dY / W.b);
    nbPan[(nbPan == 0.)] = 1;
    nbY =  int(sum(nbPan)+1);
    i = 0;
    y = np.zeros(nbY,dtype = float);
    if len(discY) != 3:
        theta_i = np.linspace(0.,m.pi,nbPan[i]+1);
        eta_i = 0.5 * (np.cos(theta_i) + 1.);
        y_i = discY[i] * eta_i + (1. - eta_i) * discY[i+1];
        y[:int(nbPan[i])] = y_i[:-1];
        i = 1;
        while discY[i+1] != 0. and discY[i+1] != YTip[1]:
            theta_i = np.linspace(0.,m.pi,nbPan[i]+1);
            eta_i = 0.5 * (np.cos(theta_i) + 1.);
            y_i = discY[i] * eta_i + (1. - eta_i) * discY[i+1];
            y[int(sum(nbPan[:i])):int(sum(nbPan[:i+1]))] = y_i[:-1];
            i += 1;
        if discY[i+1] == 0.:
            theta_i = np.linspace(0.,m.pi,nbPan[i]+nbPan[i+1]+1);
            eta_i = 0.5 * (np.cos(theta_i) + 1.);
            y_i = discY[i] * eta_i + (1. - eta_i) * discY[i+2];
            y[int(sum(nbPan[:i])):int(sum(nbPan[:i+2]))+1] = y_i;
            i += 2;
        while discY[i] != YTip[1] and i+1 < len(discY):
            theta_i = np.linspace(0.,m.pi,nbPan[i]+1);
            eta_i = 0.5 * (np.cos(theta_i) + 1.);
            y_i = discY[i] * eta_i + (1. - eta_i) * discY[i+1];
            y[int(sum(nbPan[:i]))+1:int(sum(nbPan[:i+1]))+1] = y_i[1:];
            i += 1;
    else:
        theta_i = np.linspace(0.,m.pi,nbPan[i]+nbPan[i+1]+1);
        eta_i = 0.5 * (np.cos(theta_i) + 1.);
        y_i = discY[i] * eta_i + (1. - eta_i) * discY[i+2];
        y[int(sum(nbPan[:i])):int(sum(nbPan[:i+2]))+1] = y_i;
#    tol = W.getSpan()*0.25/W.getR();
#    for i in range(W.getR()):
#        for j in range(len(discY)):
#            if np.abs(W.getY(i)-discY[j]) < tol:
#                W.y[i] = 0.;
#    tol = W.getSpan()/W.getR();
#    for i in range(W.getR()):
#        for j in range(len(YEng)):
#            if np.abs(W.getY(i)-YEng[j]) < tol:
#                W.y[i] = 0.;
    W.setY(y);
    W.setR(len(W.getY())-1);
    W.setYP((W.getY(range(1,1+W.getR()))+W.getY(range(0,W.getR())))*0.5);
    W.setEta(W.getY()*2./W.getSpan());
    W.setTheta(np.round(np.arccos(W.getEta()),12));
    W.setDiscThetaPol(np.round(np.arccos(W.getPolD()/(WI.getSpan()*0.5)),12));
    W.setCF(WI.getCF());

    flapsZone = np.zeros(W.getR()/2,dtype = int);
    dfr = WI.getDFR();
    dfl = WI.getDFL();
    if np.any(dfr) or np.any(dfl):
        for i in range(len(WI.getFlapsD())):
            indice = np.array([W.getYP(ii) >= WI.getFlapsD(i) for ii in range(W.getR()/2, W.getR())],dtype = bool);
            flapsZone += indice;
    W.addDFL(dfl[flapsZone]);
    W.addDFR(dfr[flapsZone]);
    
    afZone = np.zeros(W.getR(),dtype = int);
    if np.any(WI.getAfD()):
        for i in range(len(WI.getAfD())):
            indice = np.array([W.getYP(ii) >= WI.getAfD(i) for ii in range(W.getR()/2, W.getR())],dtype = bool);
            incr = np.concatenate([np.flipud(indice),indice]);
            afZone += incr;
    W.addAFI(afZone);
    
    polZone = np.zeros(W.getR(),dtype = int);
    if np.any(WI.getPolD()):
        for i in range(len(WI.getPolD())):
            indice = np.array([W.getYP(ii) >= WI.getPolD(i) for ii in range(W.getR()/2, W.getR())],dtype = bool);
            incr = np.concatenate([np.flipud(indice),indice]);
            polZone += incr;
    W.addPolI(polZone);
    
    taper =  WI.getTapD();
    if np.any(taper):
        x = np.round(np.concatenate([[0.],taper,[W.getSpan()*0.5]]),12);
    else:
        x = np.round(np.array([0.,W.getSpan()*0.5]),12);
    W.addChordDist(np.interp(np.absolute(W.getY()), x, WI.getChord()));
    W.addChord(np.interp(np.absolute(W.getYP()),W.getY(),W.getChordDist()));
    
    
    sweepZone = np.zeros(W.getR(),dtype = int);
    sweep = WI.getSweepC4()*m.pi/180.;
    if np.any(WI.getSweepD()):
        for i in range(len(WI.getSweepD())):
            indice = np.array([W.getYP(ii) > WI.getSweepD(i) for ii in range(W.getR()/2, W.getR())],dtype = bool);
            incr = np.concatenate([np.flipud(indice),indice]);
            sweepZone += incr;
    W.setSweepC4(sweep[sweepZone]);
    span = W.getY(range(W.getR()/2+1,W.getR()+1))-W.getY(range(W.getR()/2,W.getR()));
    x = np.zeros(W.getR()/2+1);
    for i in range(W.getR()/2):
        x[i+1] = x[i] + span[i] * m.tan(W.getSweepC4(W.getR()/2+i));
    W.setX(np.concatenate([np.flipud(x[1:]),x]));
    W.setXP((W.getX(range(1,W.getR()+1))+W.getX(range(0,W.getR())))*0.5);
    
    dihZone = np.zeros(W.getR(),dtype = int);
    dih = WI.getDih()*m.pi/180.;
    if np.any(WI.getDihD()):
        for i in range(len(WI.getDihD())):
            indice = np.array([W.getYP(ii) > WI.getDihD(i) for ii in range(W.getR()/2, W.getR())],dtype = bool);
            incr = np.concatenate([np.flipud(indice),indice]);
            dihZone += incr;
    W.setDih(dih[dihZone]);
    
    z = np.zeros(W.getR()/2+1);
    for i in range(W.getR()/2):
        z[i+1] = z[i] + span[i] * m.tan(W.getDih(W.getR()/2+i));
    W.setZ(np.concatenate([np.flipud(z[1:]),z]));
    W.setZP((W.getZ(range(1,W.getR()+1))+W.getZ(range(0,W.getR())))*0.5);
    
    if np.any(WI.getTwistD()):
        x = np.round(np.concatenate([[0.],WI.getTwistD(),[WI.getSpan()*0.5]]),12);
    else:
        x = np.round(np.array([0.,WI.getSpan()*0.5]),12);
    twist = WI.getTwist()*m.pi/180;
    W.setTwist(np.interp(np.absolute(W.getYP()), x, twist));
    W.twSec = np.interp(np.absolute(y), x, twist)
    return W;

def getYEng(W,prop):
    if prop.bool:
        bDemi = W.getSpan()*0.5;
        nbE = len(prop.getTc());
        rP = np.zeros(nbE,dtype = float);
        dP = np.zeros(nbE,dtype = float);
        yP = prop.getYp();
        zP = prop.getZp();
        Rh = prop.getRh();
        D = prop.getD();
        nInter = 5 + 2 * 0;
        for ii in range(nbE):
            if abs((zP[ii]-W.vDist)) < Rh[ii]:
                rP[ii] = abs(Rh[ii]*m.cos(m.asin((zP[ii]-W.vDist)/Rh[ii])));
            if abs((zP[ii]-W.vDist)) < D[ii]*0.5:
                dP[ii] = abs(D[ii]*m.cos(m.asin(2.*(zP[ii]-W.vDist)/D[ii])));
        Y = np.zeros(nInter*nbE,dtype = float);
        for ii in range(nbE):
            if (yP[ii] > -bDemi) and (yP[ii] < bDemi):
                Y[nInter*ii] = yP[ii];
            if (yP[ii]-dP[ii]*0.5 > -bDemi) and (yP[ii] - dP[ii]*0.5 < bDemi):
                Y[nInter*ii+1] = yP[ii] - dP[ii]*0.5;
#            if (yP[ii] - rP[ii] - (dP[ii]*0.5 - rP[ii]) * 0.5 > -bDemi) and (yP[ii] - rP[ii] - (dP[ii]*0.5 - rP[ii]) * 0.5 < bDemi):
#                Y[nInter*ii+5] = yP[ii] - rP[ii] - (dP[ii]*0.5 - rP[ii]) * 0.5; 
#            if (yP[ii] - rP[ii] - (dP[ii]*0.5 - rP[ii])/3. > -bDemi) and (yP[ii] - rP[ii] - (dP[ii]*0.5 - rP[ii]) /3. < bDemi):
#                Y[nInter*ii+3] = yP[ii] - rP[ii] - (dP[ii]*0.5 - rP[ii])/3.; 
            if (yP[ii]-rP[ii] > -bDemi) and (yP[ii] - rP[ii] < bDemi):
                Y[nInter*ii+2] = yP[ii] - rP[ii];
            if (yP[ii]+rP[ii] > -bDemi) and (yP[ii] + rP[ii] < bDemi):
                Y[nInter*ii+3] = yP[ii] + rP[ii];
#            if (yP[ii] + rP[ii] + (dP[ii]*0.5 - rP[ii]) * 0.5 > -bDemi) and (yP[ii] + rP[ii] + (dP[ii]*0.5 - rP[ii]) * 0.5 < bDemi):
#                Y[nInter*ii+6] = yP[ii] + rP[ii] + (dP[ii]*0.5 - rP[ii]) * 0.5;
#            if (yP[ii] + rP[ii] + (dP[ii]*0.5 - rP[ii]) * 2./3. > -bDemi) and (yP[ii] + rP[ii] + (dP[ii]*0.5 - rP[ii]) * 2./3. < bDemi):
#                Y[nInter*ii+7] = yP[ii] + rP[ii] + (dP[ii]*0.5 - rP[ii]) * 2./3.;
            if (yP[ii]+dP[ii]*0.5 > -bDemi) and (yP[ii] + dP[ii]*0.5 < bDemi):
                Y[nInter*ii+4] = yP[ii] + dP[ii]*0.5;
        Y = np.unique(Y);
    else:
        Y = np.zeros(1,dtype = float);
    return Y;

def getYFus(W,fus):
    if fus.bool and abs(m.acos(fus.getVD()*2./fus.getCD())) < 1.:
        Angle = abs(m.acos(fus.getVD()*2./fus.getCD()));
        y = m.sin(Angle) * fus.getCD()*0.5 * np.array([-1.,1.],dtype = float)
    else:
        y = np.empty(0,dtype = float);
    return y;