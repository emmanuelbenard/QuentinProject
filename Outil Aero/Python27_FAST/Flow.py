# -*- coding: utf-8 -*-
import utilitaire as u;
class Flow():
    """ Classe définissant l'écoulement reprenant les diférentes
    charactéristiques de celles ci:
        V0 : flow velocity (m/s)
        aMin : minimum a.o.a. of the study (deg)
        self.DA : interval between two a.o.a. for the study (deg)
        aMax : maximum a.o.a. of the study (deg)
        beta : sideslip angle of the study (deg)
        h : altitude of the flight conditions (m)
        M : Mach number of the flight conditions
    """
    def __init__(self):
        self.V0 = 1.;
        self.aMin = 0.;
        self.DA = 0.25;
        self.aMax = 10.;
        self.beta = 0.;
        self.h = 0.;
        self.M = 0.;
        self.mW = 8;
        self.at = 0.;
        
        
    def setV0(self,V):
        """ Définit la vitesse de vol pour l'étude"""
        self.V0 = V;
    def getV0(self):
        """ Retourne la vitesse de vol pour l'étude"""
        return self.V0;
    
    def setAMin(self,a):
        """ Définit l'a.o.a. min de l'étude en degrés"""
        self.aMin = a;
    def getAMin(self):
        """ Retourne l'a.o.a. min de l'étude en degrés"""
        return self.aMin;
    def setAMax(self,a):
        """ Définit l'a.o.a. max de l'étude en degrés"""
        self.aMax = a;
    def getAMax(self):
        """ Retourne l'a.o.a. max de l'étude en degrés"""
        return self.aMax;
    def setDA(self,a):
        """ Définit l'intervalle d'a.o.a. pour l'étude en degrés"""
        self.DA = a;
    def getDA(self):
        """ Retourne l'intervalle d'a.o.a. pour l'étude en degrés"""
        return self.DA;
    def setBeta(self,b):
        """ Définit l'angle de sideslipe pour l'étude en degrés"""
        self.beta = b;
    def getBeta(self):
        """ Retourne l'angle de sideslipe pour l'étude en degrés"""
        return self.beta;
    def setH(self,h):
        """ Définit l'altitude de vol"""
        self.h = h;
    def getH(self):
        """Retourne l'altitude de vol"""
        return self.h;
    def setMach(self,M):
        self.M = M;
    def getMach(self):
        return self.M;
    
    def getAlphas(self):
        """Retourne less a.o.a. de l'étude en degrés"""
        return u.rangeF(self.aMin, self.aMax, self.DA);