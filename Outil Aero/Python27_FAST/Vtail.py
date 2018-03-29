# -*- coding: utf-8 -*-
import numpy as np;
import math as m;
class Vtail():
    """ Classe définissant une dérive verticale dans son ensemble reprenant les diférentes
    charactéristiques de celles ci:

        b : span
        chord : vertical array with the chord at the sections;
        flapsDiscZ : vertical array with the spanwise coordinate of the flaps
           discontinuity sections
        afDiscZ : vertical array with the spanwise coordinate of the
           airfoil discontinuity sections
        polarDiscZ : vertical array with the spanwise coordinate of the
           polar discontinuity sections (2D curves)
        taperdiscZ : vertical array with the spanwise coordinate of the
           taper ratio discontinuity sections
        sweepDiscZ : vertical array with the spanwise coordinate of the
           sweep angle discontinuity sections
        FTC : the ratio between the flaps length and the local chord
           length
        polar : a cell-array with for each cell a "polar" object with field
           data (the alpha-cl-cd data), clmax, alphal0 (zero-lift angle).
        airfoil : a cell-array with each cell gives the airfoil naca number
           representation, cell 1 correspond to first panel after root.
        sweepC4 : vertical array with vtail.sweep(i) is the sweep angle of 
           the panel from vtail.y(i) to vtail.y(i+1) at c/4 (rad)
        DF : vertical array with the different flaps defection
           along the span (deg)
        r : number of spanwise panel along the vtail;
        p : number of chordwise panel along the airfoil;
        cFlaps_cLoc : vertical array with vtail.cFlaps_cLocs(i) is the 
           local flaps to chord ratio
        theta : the angle (0 -> pi/2) spanwise location of the limits of the panels
        eta : the normalized (0 -> 1) spanwise location of the limits of the panels
        y : the  spanwise location of (-b/2 -> b/2) the limits of the panels
        discZ : vertical array of the complete set of the spanwise location
           of the discontinuity
        discTheta : same as above but in angle coordinates;
        airfoilIndex : vertical array with vtail.airfoilIndex(i) is the index of 
           the airfoil (vtail.airfoil) to use for the section at vtail.y(i)
        polarIndex : vertical array with vtail.polarIndex(i) is the index of 
           the airfoil (vtail.polar) to use for the section at vtail.y(i)
        chordDistrib : vertical array with vtail.chordDistrib(i) is the chord length of 
           the section at vtail.z(i)
         bool : un booléen qui dit si élément présent sur l'avion

     """
    def __init__(self):
        """ Définition des attributs pour le premier constructeur """
        # general properties
        self.b = 0;
        self.mC = 8;
        self.bool = False;
        self.p = 10;
        self.x = np.empty(0,dtype = float);
        self.y = np.empty(0,dtype = float);
        self.z = np.empty(0,dtype = float);
        self.xP = np.empty(0,dtype = float);
        self.yP = np.empty(0,dtype = float);
        self.zP = np.empty(0,dtype = float);
        self.theta = np.empty(0,dtype = float);
        self.eta = np.empty(0,dtype = float);
        self.FTC = 0;

        # Geometric distribution along the span
        self.chordDistrib = np.empty(0,dtype = float);
        self.sweepc4 = np.zeros(1,dtype = float);
        self.airfoilIndex = np.empty(0,dtype=int);
        self.polarIndex = np.empty(0,dtype=int);
        self.chord = np.empty(0,dtype=float);

        # Discontinuity properties
        self.flapsDiscZ = np.empty(0,dtype = float);
        self.afDiscZ = np.empty(0,dtype = float);
        self.polarDiscZ = np.empty(0,dtype = float);
        self.taperdiscZ = np.empty(0,dtype = float);
        self.sweepDiscZ = np.empty(0,dtype = float);
        self.dihDiscZ = np.empty(0,dtype = float);
        self.twistdiscZ = np.empty(0,dtype = float);
        self.discTheta = np.empty(0,dtype = float);
        self.discZPol = np.empty(0,dtype = float);

        self.polar = np.empty(0,dtype = object);
        self.airfoil = np.empty(0,dtype = str);
        self.DF = np.empty(0,dtype = float);
        
        self.hDist = 0.;
        self.vDist = 0.;
        

        

    def incSpan(self,b):
        """ Incrémente le span de la dérive verticale de b"""
        self.b += b;
    def decSpan(self,b):
        """ Diminue le span de la dérive verticale de b"""
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

    def getP(self):
        """ Retourne le nombre de point de controle le long de la corde """
        return self.p;

    
    def setX(self,x):
        """ Définit la position (x-comp) des quarts de corde des 
        sections de controle le long du span, negative : left vtail """
        self.x = x;
    def getX(self,ii='all'):
        """ Retourne la position (x-comp) des quarts de corde des 
        sections de controle le long du span, negative : left vtail """
        if ii == 'all':
            return self.x;
        else:
            return self.x[ii];
    
    def setY(self,y):
        """ Définit la position (y-comp) des quarts de corde des 
        sections de controle le long du span, negative : left vtail """
        self.y = y;
    def getY(self,ii='all'):
        """ Retourne la position (y-comp) des quarts de corde des 
        sections de controle le long du span, negative : left vtail """
        if ii == 'all':
            return self.y;
        else:
            return self.y[ii];

    def setZ(self,z):
        """ Définit la position (z-comp) des quarts de corde des 
        sections de controle le long du span, negative : left vtail """
        self.z = z;
    def getZ(self,ii='all'):
        """ Retourne la position (z-comp) des quarts de corde des 
        sections de controle le long du span, negative : left vtail """
        if ii == 'all':
            return self.z;
        else:
            return self.z[ii];
        
        
    def setXP(self,x):
        """ Définit la position (x-comp) des quarts de corde des 
        panneaux le long du span, negative : left vtail """
        self.xP = x;
    def getXP(self,ii='all'):
        """ Retourne la position (x-comp) des quarts de corde des 
        panneaux le long du span, negative : left vtail """
        if ii == 'all':
            return self.xP;
        else:
            return self.xP[ii];
    
    def setYP(self,y):
        """ Définit la position (y-comp) des quarts de corde des 
        panneaux le long du span, negative : left vtail """
        self.yP = y;
    def getYP(self,ii='all'):
        """ Retourne la position (y-comp) des quarts de corde des 
        panneaux le long du span, negative : left vtail """
        if ii == 'all':
            return self.yP;
        else:
            return self.yP[ii];

    def setZP(self,z):
        """ Définit la position (z-comp) des quarts de corde des 
        panneaux le long du span, negative : left vtail """
        self.zP = z;
    def getZP(self,ii='all'):
        """ Retourne la position (z-comp) des quarts de corde des 
        panneaux le long du span, negative : left vtail """
        if ii == 'all':
            return self.zP;
        else:
            return self.zP[ii];
    def setHDist(self,l):
        """ Définit la distance horizontale entre la position c/4 de 
        l'aile et la position c/4 de la dérive vertical"""
        self.hDist = l;
    def getHDist(self):
        """ Retourne la distance horizontale entre la position c/4 de 
        l'aile et la position c/4 de la dérive vertical"""
        return self.hDist;
    def setVDist(self,l):
        """ Définit la distance verticale entre la position c/4 de 
        l'aile et la position c/4 de la dérive vertical"""
        self.vDist = l;
    def getVDist(self):
        """ Retourne la distance verticale entre la position c/4 de 
        l'aile et la position c/4 de la dérive vertical"""
        return self.vDist;
    
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
        self.sweepc4 = np.append(np.empty(0,dtype = float),phi);
    def addSweepC4(self,phi):
        """ Ajoute l'angle de flèche au quart de corde pour chaque panneau"""
        self.sweepc4 = np.append(self.sweepc4,phi);
    def setSweepC4(self,phi):
        """ Definit l'angle de flèche au quart de corde pour chaque panneau"""
        self.sweepc4 = phi;
    def getSweepC4(self,ii='all'):
        """ Retourne l'angle de flèche au quart de corde pour chaque panneau"""
        if ii == 'all':
            return self.sweepc4;
        else:
            return self.sweepc4[ii];

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

    def addDF(self,delta):
        """ Ajoute la déflection du volet sur chaque panneau """
        self.DF = np.append(self.DF,delta);
    def getDF(self,ii='all'):
        """ Retourne la déflection du volet sur chaque panneau """
        if ii == 'all':
            return self.DF;
        else:
            return self.DF[ii];

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
        le long du span de la dérive verticale"""
        self.flapsDiscZ = np.append(self.flapsDiscZ,y);
    def setFlapsD(self,y):
        """ Définit les positions de discontinuité de déflection de volet
        le long du span de la dérive verticale"""
        self.flapsDiscZ = y;
    def rmFlapsD(self):
        """ Enlève la dernière position de discontinuité de déflection de volet
        le long du span de la dérive verticale"""
        self.flapsDiscZ = self.flapsDiscZ[:-1];
    def getFlapsD(self,ii='all'):
        """ Retourne les positions de discontinuité de déflection de volet
        le long du span de la dérive verticale"""
        if ii =='all':
            return self.flapsDiscZ;
        else:
            return self.flapsDiscZ[ii];

    def addAfD(self,y):
        """ Ajoute une position de discontinuité de profil
        le long du span de la dérive verticale"""
        self.afDiscZ = np.append(self.afDiscZ,y);
    def setAfD(self,y):
        """ Définit les positions de discontinuité de profil
        le long du span de la dérive verticale"""
        self.afDiscZ = y;
    def rmAfD(self):
        """ Supprime la dernière position de discontinuité de profil
        le long du span de la dérive verticale"""
        self.afDiscZ = self.afDiscZ[:-1];
    def getAfD(self,ii='all'):
        """ Retourne les positions de discontinuité de profil
        le long du span de la dérive verticale"""
        if ii =='all':
            return self.afDiscZ;
        else:
            return self.afDiscZ[ii];

    def addPolD(self,y):
        """ Ajoute une position de discontinuité de fichier polar
        le long du span de la dérive verticale"""
        self.polarDiscZ = np.append(self.polarDiscZ,y);
    def setPolD(self,y):
        """ Définit les positions de discontinuité de fichier polar
        le long du span de la dérive verticale"""
        self.polarDiscZ = y;
    def rmPolD(self):
        """ Supprime la dernière position de discontinuité de fichier polar
        le long du span de la dérive verticale"""
        self.polarDiscZ = self.polarDiscZ[:-1];
    def getPolD(self,ii='all'):
        """ Retourne les positions de discontinuité de fichier polar
        le long du span de la dérive verticale"""
        if ii =='all':
            return self.polarDiscZ;
        else:
            return self.polarDiscZ[ii];

    def addTapD(self,y):
        """ Ajoute une position de discontinuité de taper
        le long du span de la dérive verticale"""
        self.taperdiscZ = np.append(self.taperdiscZ,y);
    def setTapD(self,y):
        """ Définit les positions de discontinuité de taper
        le long du span de la dérive verticale"""
        self.taperdiscZ = y;
    def rmTapD(self):
        """ Supprime la dernère position de discontinuité de taper
        le long du span de la dérive verticale"""
        self.taperdiscZ = self.taperdiscZ[:-1];
    def getTapD(self,ii='all'):
        """ Retourne les positions de discontinuité de taper
        le long du span de la dérive verticale"""
        if ii =='all':
            return self.taperdiscZ;
        else:
            return self.taperdiscZ[ii];

    def addSweepD(self,y):
        """ Ajoute une position de discontinuité d'angle de flèche
        le long du span de la dérive verticale"""
        self.sweepDiscZ = np.append(self.sweepDiscZ,y);
    def setSweepD(self,y):
        """ Définit les positions de discontinuité d'angle de flèche
        le long du span de la dérive verticale"""
        self.sweepDiscZ = y;
    def rmSweepD(self):
        """ Supprime la dernière position de discontinuité d'angle de flèche
        le long du span de la dérive verticale"""
        self.sweepDiscZ = self.sweepDiscZ[:-1];
    def getSweepD(self,ii='all'):
        """ Retourne les positions de discontinuité d'angle de flèche
        le long du span de la dérive verticale"""
        if ii =='all':
            return self.sweepDiscZ;
        else:
            return self.sweepDiscZ[ii];

    def addDihD(self,y):
        """ Ajoute une position de discontinuité d'angle de dihèdre
        le long du span de la dérive verticale"""
        self.dihDiscZ = np.append(self.dihDiscZ,y);
    def setDihD(self,y):
        """ Définit les positions de discontinuité d'angle de dihèdre
        le long du span de la dérive verticale"""
        self.dihDiscZ = y;
    def rmDihD(self):
        """ Supprime la dernière position de discontinuité d'angle de dihèdre
        le long du span de la dérive verticale"""
        self.dihDiscZ = self.dihDiscZ[:-1];
    def getDihD(self,ii='all'):
        """ Retourne les positions de discontinuité d'angle de dihèdre
        le long du span de la dérive verticale"""
        if ii =='all':
            return self.dihDiscZ;
        else:
            return self.dihDiscZ[ii];

    def addTwistD(self,y):
        """ Ajoute une position de discontinuité d'angle de twist
        le long du span de la dérive verticale"""
        self.twistdiscZ = np.append(self.twistdiscZ,y);
    def setTwistD(self,y):
        """ Définit les positions de discontinuité d'angle de twist
        le long du span de la dérive verticale"""
        self.twistdiscZ = y;
    def rmTwistD(self):
        """ Supprime la dernière position de discontinuité d'angle de twist
        le long du span de la dérive verticale"""
        self.twistdiscZ = self.twistdiscZ[:-1];
    def getTwistD(self,ii='all'):
        """ Retourne les positions de discontinuité d'angle de twist
        le long du span de la dérive verticale"""
        if ii =='all':
            return self.twistdiscZ;
        else:
            return self.twistdiscZ[ii];

    def setDiscTheta(self,theta):
        """ Définit les positions de discontinuité en coordonnées theta
        le long du span de la dérive verticale"""
        self.discTheta = theta;
    def getThetaD(self):
        """ Retourne les positions de discontinuité en coordonnées theta
        le long du span de la dérive verticale"""
        return self.discTheta;

    def setDiscThetaPol(self,theta):
        """ Définit les positions de discontinuité de polair en coordonnées theta
        le long du span de la dérive verticale"""
        self.discThetaPol = theta;
    def getThetaDPol(self):
        """ Retourne les positions de discontinuité de polair en coordonnées theta
        le long du span de la dérive verticale"""
        return self.discThetaPol;

    def addPolar(self,polar):
        """ Ajoute un object polar """
        self.polar = np.append(self.polar,polar);
    def getPolar(self,ii='all'):
        """ Retourne les objects polar """
        if ii == 'all':
            return self.polar;
        else:
            return self.polar[ii];

    def addAF(self,af):
        """ Ajoute un profil"""
        self.airfoil = np.append(self.airfoil,af);
    def getAF(self,ii='all'):
        """ Retourne les profils"""
        if ii == 'all':
            return self.airfoil[:];
        else:
            return self.airfoil[ii];

    def addDF(self,delta):
        """ Ajoute une déflection de volet pour la dérive verticale """
        self.DF = np.append(self.DF,delta);
    def getDF(self,ii='all'):
        """ Retourne les déflections de volet pour la dérive verticale"""
        if ii == 'all':
            return self.DF[:];
        else:
            return self.DF[ii];

    

    def addSeg(self,panel):
        """ Ajoute un segment d'aile """
        self.section = np.append(self.section,panel);
    def getSeg(self,ii):
        """ Retourne un segment d'aile """
        temp = np.empty([6,65],dtype = float);
        temp[:3,:] = self.section[ii].getBeg();
        temp[3:,] = self.section[ii].getEnd();
        return np.transpose(temp);

def VtailManager(VTI):
    VT = Vtail();
    for i in range(len(VTI.getPolar())):
        VT.addPolar(VTI.polar[i]);
    VT.addAF(VTI.getAF());
    VT.hDist = VTI.hDist;
    VT.vDist = VTI.vDist;
    r = 16;
    VT.setR(r);
    
    theta = np.linspace(0.,1.,r);
    VT.setEta(theta);
    
    
    VT.incSpan(VTI.getSpan());
    VT.setZ(VT.getEta()*VT.getSpan());
    VT.addTapD(np.sort(VTI.getTapD()));
    VT.addSweepD(VTI.getSweepD());
    VT.addAfD(VTI.getAfD());
    VT.addFlapsD(VTI.getFlapsD());
    VT.addPolD(VTI.getPolD());
    discZ = np.unique(np.concatenate([VT.getTapD(),VT.getSweepD(),\
        VT.getAfD(),VT.getFlapsD(),VT.getPolD(),[0.,VT.getSpan()]]));
    
    dZ = discZ[1:] - discZ[:-1];
    nbPan = np.round(r * dZ / VT.b);
    nbPan[(nbPan == 0.)] = 1;
    nbZ =  int(sum(nbPan)+1);
    i = 0;
    z = np.zeros(nbZ,dtype = float);
    
    if len(discZ) != 2:
        theta_i = np.linspace(0.,m.pi,nbPan[i]+1);
        eta_i = 0.5 * (np.cos(theta_i) + 1.);
        z_i = discZ[i] * eta_i + (1. - eta_i) * discZ[i+1];
        z[:int(nbPan[i])] = z_i[:-1];
        i = 1;
        while discZ[i+1] != VT.getSpan():
            theta_i = np.linspace(0.,m.pi,nbPan[i]+1);
            eta_i = 0.5 * (np.cos(theta_i) + 1.);
            z_i = discZ[i] * eta_i + (1. - eta_i) * discZ[i+1];
            z[int(sum(nbPan[:i])):int(sum(nbPan[:i+1]))] = z_i[:-1];
            i += 1;
    else:
        theta_i = np.linspace(0.,m.pi,nbPan[i]+1);
        eta_i = 0.5 * (np.cos(theta_i) + 1.);
        z_i = discZ[i] * eta_i + (1. - eta_i) * discZ[i+1];
        z[int(sum(nbPan[:i])):int(sum(nbPan[:i+1]))+1] = z_i;
    
    VT.setZ(z);
    VT.setR(len(VT.getZ())-1);
    VT.setZP((VT.getZ(range(1,1+VT.getR()))+VT.getZ(range(0,VT.getR())))*0.5);
    VT.setEta(VT.getZ()/VT.getSpan());
    VT.setTheta(np.round(np.arccos(VT.getEta()),12));
    VT.setDiscThetaPol(np.round(np.arccos(VT.getPolD()/(VTI.getSpan())),12));
    VT.setCF(VTI.getCF());
    
    VT.setY(np.zeros(VT.getR()+1,dtype = float));
    VT.setYP(np.zeros(VT.getR(),dtype = float));
    
    flapsZone = np.zeros(VT.getR(),dtype = int);
    df = VTI.getDF();
    if np.any(df):
        for i in range(len(VTI.getFlapsD())):
            indice = np.array([VT.getZP(ii) >= VTI.getFlapsD(i) for ii in range(VT.getR())],dtype = bool);
            flapsZone += indice;
    VT.addDF(df[flapsZone]);
    
    afZone = np.zeros(VT.getR(),dtype = int);
    if np.any(VTI.getAfD()):
        for i in range(len(VTI.getAfD())):
            indice = np.array([VT.getZP(ii) >= VTI.getAfD(i) for ii in range(VT.getR())],dtype = bool);
            afZone += indice;
    VT.addAFI(afZone);
    
    polZone = np.zeros(VT.getR(),dtype = int);
    if np.any(VTI.getPolD()):
        for i in range(len(VTI.getPolD())):
            indice = np.array([VT.getZP(ii) >= VTI.getPolD(i) for ii in range(VT.getR())],dtype = bool);
            polZone += indice;
    VT.addPolI(polZone);
    
    taper =  VTI.getTapD();
    if np.any(taper):
        x = np.round(np.concatenate([[0.],taper,[VT.getSpan()*0.5]]),12);
    else:
        x = np.round(np.array([0.,VT.getSpan()]),12);
    VT.addChordDist(np.interp(np.absolute(VT.getZ()), x, VTI.getChord()));
    VT.addChord(np.interp(np.absolute(VT.getZP()),VT.getZ(),VT.getChordDist()));
    
    
    sweepZone = np.zeros(VT.getR(),dtype = int);
    sweep = VTI.getSweepC4()*m.pi/180.;
    if np.any(VTI.getSweepD()):
        for i in range(len(VTI.getSweepD())):
            indice = np.array([VT.getYP(ii) > VTI.getSweepD(i) for ii in range(VT.getR()/2, VT.getR())],dtype = bool);
            incr = np.concatenate([np.flipud(indice),indice]);
            sweepZone += incr;
    VT.setSweepC4(sweep[sweepZone]);
    span = VT.getZ(range(1,VT.getR()+1))-VT.getZ(range(VT.getR()));
    x = np.zeros(VT.getR()+1);
    for i in range(VT.getR()):
        x[i+1] = x[i] + span[i] * m.tan(VT.getSweepC4(i));
    VT.setX(x);
    VT.setXP((VT.getX(range(1,VT.getR()+1))+VT.getX(range(0,VT.getR())))*0.5);
    
    VT.x += VT.hDist;
    VT.xP += VT.hDist;
    
    VT.z += VT.vDist;
    VT.zP += VT.vDist;
    VT.taperdiscZ += VT.vDist;
    VT.sweepDiscZ += VT.vDist;
    VT.afDiscZ += VT.vDist;
    VT.flapsDiscZ += VT.vDist;
    VT.polarDiscZ += VT.vDist;
    return VT;