# -*- coding: utf-8 -*-
import math as m
import numpy as np
import numpy.matlib
import copy
import matplotlib.pyplot as plt
import utilitaire as u
class Panel:
	""" Classe définissant un segment de span de l'aile.
		Classe munie des champs suivants:
			- begin section
			- end section

	 """
	def __init__(self):
		""" Définition des attributs pour le premier constructeur """
		self.begin = np.zeros([65,3],dtype = float);
		self.end = np.zeros([65,3],dtype = float);

	def setBeg(self,begX,begY,begZ):
		self.begin[:,0] = begX;
		self.begin[:,1] = begY;
		self.begin[:,2] = begZ;
	def getBeg(self):
		return self.begin[:,0],self.begin[:,1],self.begin[:,2];

	def setEnd(self,endX,endY,endZ):
		self.end[:,0] = endX;
		self.end[:,1] = endY;
		self.end[:,2] = endZ;
	def getEnd(self):
		return self.end[:,0],self.end[:,1],self.end[:,2];