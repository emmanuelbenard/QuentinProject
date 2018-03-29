# -*- coding: utf-8 -*-
import Wing;
class Htail(Wing.Wing):

    def __init__(self):
        Wing.Wing.__init__(self);
        self.hDist = 0.; #Only for htail
        self.vDist = 0.;
        self.bool = False;

    def setHDist(self,l):
        """ Définit la distance horizontale entre la position c/4 de 
        l'aile et la position c/4 de l'empennage"""
        self.hDist = l;
    def getHDist(self):
        """ Retourne la distance horizontale entre la position c/4 de 
        l'aile et la position c/4 de l'empennage"""
        return self.hDist;
    def setVDist(self,l):
        """ Définit la distance verticale entre la position c/4 de 
        l'aile et la position c/4 de l'empennage"""
        self.vDist = l;
    def getVDist(self):
        """ Retourne la distance verticale entre la position c/4 de 
        l'aile et la position c/4 de l'empennage"""
        return self.vDist;
def HtailManager(HTI,prop,fus):
    HT = Htail();
    HTI.r = 20;
    W = Wing.WingManager(HTI,prop,fus);
    HT.__dict__ = W.__dict__;
    HT.setHDist(HTI.getHDist());
    HT.setVDist(HTI.getVDist());
    HT.x += HT.hDist;
    HT.xP += HT.hDist;
    HT.z += HT.vDist;
    HT.zP += HT.vDist;
    return HT;