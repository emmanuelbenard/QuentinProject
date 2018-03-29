# -*- coding: utf-8 -*-
import Wing;
import Htail;
import Vtail;
import Fuselage;
import Propu;
import numpy as np;
class Aircraft_aero():
    """ Classe qui reprend l'avion en son entier:
            wing : object qui représente l'aile
            htail : objet qui représente l'empennage horizontale
            vtail : objet qui représente la dérive verticale
            fus : objet qui représente le fuselage
            prop : objet qui représente le système propulsif
            rp : np.array le point de référence pour le calcul des moment. Par rapport
                au wing root quater chord section point (0,0,0).
    """
    def __init__(self):
        self.wing = Wing.Wing();
        self.htail = Htail.Htail();
        self.vtail = Vtail.Vtail();
        self.fus = Fuselage.Fuselage();
        self.prop = Propu.Propu();
        self.rp = np.zeros(3,dtype=float);
    
    def setWing(self,wing):
        """ Attribuer une aile à l'avion"""
        self.wing = wing;
    def getWing(self):
        """ Retourne l'aile de l'avion"""
        return self.wing;
    
    def setHT(self,htail):
        """ Attribuer un empennage à l'avion"""
        self.htail = htail;
    def getHT(self):
        """ Retourne l'empennage de l'avion"""
        return self.htail;
    
    def setVT(self,vtail):
        """ Attribuer une dérive verticale à l'avion"""
        self.vtail = vtail;
    def getVT(self):
        """ Retourne la dérive de l'avion de l'avion"""
        return self.vtail;
    
    
    def setFus(self,fus):
        """ Attribuer un fuselage à l'avion"""
        self.fus = fus;
    def getFus(self):
        """ Retourne le fuselage de l'avion"""
        return self.fus;
    
    def setProp(self,prop):
        """ Attribuer un système propulsif à l'avion"""
        self.prop = prop;
    def getProp(self):
        """ Retourne le système propulsif de l'avion"""
        return self.prop;
    
    def setRP(self,rp):
        """ Attribuer un point de référence à l'avion pour le calcul des moments"""
        self.rp = rp;
    def getRP(self):
        """ Retourne le point de référence de l'avion"""
        return self.rp;

def AircraftManager_SEI(ACI,flow):
    AC = Aircraft_aero();
    AC.rp = ACI.rp;
    if ACI.wing.bool:
        AC.wing = Wing.WingManager(ACI.wing,ACI.prop,ACI.fus);
    if ACI.htail.bool:
        fus = Fuselage.Fuselage();
        fus.bool = False;
        AC.htail = Htail.HtailManager(ACI.htail,ACI.prop,fus);
    else:
        AC.htail.r = 0;
    if ACI.vtail.bool:
        AC.vtail = Vtail.VtailManager(ACI.vtail);
    if ACI.fus.bool:
        AC.fus = Fuselage.FusManager(ACI.fus,AC.wing,AC.vtail,flow,AC.rp);
    if ACI.prop.bool:
        AC.prop = Propu.PropManager_SEI(ACI.prop,AC.wing,AC.htail,AC.vtail,flow);
    else:
        AC.prop.vix = np.zeros(AC.wing.getR()+AC.htail.getR(),dtype = float);
        AC.prop.viy = np.zeros(AC.wing.getR()+AC.htail.getR(),dtype = float);
        AC.prop.viz = np.zeros(AC.wing.getR()+AC.htail.getR(),dtype = float);
    return AC;

def AircraftManager_DEI(ACI,flow):
    AC = Aircraft_aero();
    AC.rp = ACI.rp;
    if ACI.wing.bool:
        AC.wing = Wing.WingManager(ACI.wing,ACI.prop,ACI.fus);
    if ACI.htail.bool:
        fus = Fuselage();
        fus.bool = False;
        AC.htail = Htail.HtailManager(ACI.htail,ACI.prop,fus);
    else:
        print AC.htail.r
    if ACI.vtail.bool:
        AC.vtail = Vtail.VtailManager(ACI.vtail);
    if ACI.fus.bool:
        AC.fus = Fuselage.FusManager(ACI.fus,AC.wing,flow,AC.rp);
    if ACI.prop.bool:
        AC.prop = Propu.PropManager_DEI(ACI.prop,AC.wing,AC.htail,AC.vtail,flow);
    return AC;