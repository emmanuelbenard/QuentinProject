# -*- coding: utf-8 -*-
import math as m
import numpy as np
import numpy.matlib
import numpy.linalg
import copy
import matplotlib.pyplot as plt
import utilitaire as u
import Polar
import Panel
import sys
import os
def An(V0,h,b,D,Rhub,Omega,T,zh,xh):
	"""Coefficiets of Prandtl lifting line and Blade element theory
	An(Neidp) gives the 3 coeficients of the EIDP.
	 Neidp is the number of elements in which the wing will be divided.
	 V0: aircraft velocity                   [m/s]
	 h: altitude study                       [Km]
	 b: Semi-Span                         [m]
	 k: Clalpha/2 en chaque section       [adim]
	 l: Corde en chaque section           [m]
	 alpha0: angle d'ataque de portance nule [degrees]
	 Rhub: radius du hub                     [m]
	 D: Diametre helice                      [m]
	 Omega: vitesse de rotation de l'helice  [rad/s]
	 T: Traction de l'helice             [N]
	 zh: Position helices. (y=0 always)      [m]
		 i.e.: 
			   none: zh=[]
			   1 in the middle: zh=[0]
			   2: zh=[-3 3]
	
	RESULTS
	 an is the influence of the wing with angle of attack null
	 bn is the influence of the velocity tangential of the propeller
	 cn is the influence of des incidences de portances nulles"""

	#Calcul masse volumique
	rho0=1.225; #masse volumique à niveau de la mer           [kg/m^3]
	dT=-6.5;    #gradiente de temperature dans la troposphere [K/km]
	T0=288.15;  #Temperature à niveau de la mer               [K]
	g=9.80665;  #gravité                                      [m/s^2]
	Rair=287.;   #Constante de l'air                           [m^2/(s^2*K)]
	rho = rho0 * (1. + dT*h/T0)**(- g/(Rair*dT*10**(-3)) - 1.);
	Neidp = 100;
	##
	# Variables helice
	Sh = m.pi * D**2 *0.25;#Surface disque actuator          [m^2]
	##
	# Discretitation 
	theta = np.zeros(Neidp+1,dtype = float);
	z = np.zeros(Neidp+1,dtype = float);
	vix = np.zeros(Neidp+1,dtype = float);
	vitheta = np.zeros(Neidp+1,dtype = float);
	for i in range(Neidp+1):
		#initialization de calcul en section -b
		theta[i]=float(i)*m.pi/Neidp;
		z[i]=b*m.cos(theta[i]);
		#Calcul de l'increment de vitesse, si l'on est dans la zone d'influence 
		#de l'helice.
		if ((z[i]>=(zh-D/2) and z[i]<=(zh-Rhub))or (z[i]<=(zh+D/2) and z[i]>=(zh+Rhub))):
			vix[i] += V0*0.5 * (m.sqrt(1.+2*T/(rho*Sh*V0**2))-1.);
		if ((z[i]>=(zh-D/2) and z[i]<=(zh-Rhub))or (z[i]<=(zh+D/2) and z[i]>=(zh+Rhub))):
			vix2 = V0*0.5 * (m.sqrt(1.+2*T/(rho*Sh*V0**2))-1.);
			a = vix2/V0;
			aprim = np.nan_to_num(0.5 * ( 1. - np.sqrt( np.absolute(1. - 4*a*(1.+a)*(V0/(Omega*np.abs(z[i]-zh)))**2))));
			vitheta[i] += aprim*2.*Omega*(z[i]-zh);
	vitheta[:Neidp/2+1] *= -1;
	return vix,vitheta,z;