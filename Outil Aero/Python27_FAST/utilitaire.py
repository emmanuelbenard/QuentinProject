# -*- coding: utf-8 -*-
import math as m
import numpy as np
import sys
def naca456(n,xF=0.,delta=0.):
    """ NACA Airfoil Generator
     This function generates a set of points containing the coordinates of a
     NACA airfoil from the 4 Digit Series, 5 Digit Series and 6 Series given
     its number and, as additional features, the chordt, the number of points
     to be calculated, spacing type (between linear and cosine spacing),
     opened or closed trailing edge and the angle of attack of the airfoil.
     It also plots the airfoil for further comprovation if it is the required
     one by the user.
     
     -------------------------------------------------------------------------

     MIT License
     
     Copyright (c) 2016 Alejandro de Haro

     -------------------------------------------------------------------------

     INPUT DATA
       n --> NACA number (4, 5 or 6 digits)

     OPTIONAL INPUT DATA
       xF --> flaps to chord ratio
       delta --> flaps deflection"""
    #----------------------- Numerical parameters -------------------
    n = float(n);
    xF = float(xF);
    delta = float(delta);
    s = 33;
    #----------------------- COMPROVATION OF AIRFOIL SERIES -------------------
    if m.floor(n/1e7)==0:
        if m.floor(n/1e6)==0:
            if m.floor(n/1e5)==0:
                if m.floor(n/1e4)==0:
                    nc=4;   # NACA 4 digit series
                else:
                    nc=5;   # NACA 5 digit series
            else:
                nc=6;   # NACA 6 digit series
        else:
            nc=7;   # NACA 7 digit series
    else:
        nc=8;   # NACA 8 digit series

    beta=np.linspace(0,m.pi,s);  # Angle for cosine spacing
    x=(1-np.cos(beta))/2;  # X coordinate of airfoil (cosine spacing)

    t=float(n%100)/100;   # Maximum thickness as fraction of chord (two last digits)
    sym=0;  # Symetric airfoil variable
    #----------------------- VARIABLE PRELOCATION -----------------------------
    y_c=np.zeros(s); # Mean camber vector prelocation
    dyc_dx=np.zeros(s);  # Mean camber fisrt derivative vector prelocation
    #----------------------- THICKNESS CALCULATION ----------------------------

    y_t=t/0.2*(0.2969*np.sqrt(x)-0.126*x-0.3516*x**2+0.2843*x**3-0.1036*x**4); # Thickness y coordinate with closed trailing edge
        
    if nc==4:
    #----------------------- MEAN CAMBER 4 DIGIT SERIES CALCULATION -----------
        #----------------------- CONSTANTS ------------------------------------
        maxt=m.floor(n/1e3)/100;    # Maximum camber (1st digit)
        p=(m.floor(n/100)%10)/10;  # Location of maximum camber (2nd digit)
        if maxt==0:
            if p==0:
                sym=1;  # Comprovation of symetric airfoil with two 0
            else:
                sym=2;  # Comprovation of symetric airfoil with one 0
        #----------------------- CAMBER ---------------------------------------
        for i in range(s):
            if x[i]<p:
                y_c[i]=maxt*x[i]/p**2*(2*p-x[i]);    # Mean camber y coordinate
                dyc_dx[i]=2*maxt/p**2*(p-x[i]);    # Mean camber first derivative
            else:
                y_c[i]=maxt*(1-x[i])/(1-p)**2*(1+x[i]-2*p);    # Mean camber y coordinate
                dyc_dx[i]=2*maxt/(1-p)**2*(p-x[i]);    # Mean camber first derivative
    elif nc==5:
    #----------------------- MEAN CAMBER 5 DIGIT SERIES CALCULATION -----------
        #----------------------- CONSTANTS ------------------------------------
        p=(m.floor(n/1000)%10)/20;  # Location of maximum camber (2nd digit)
        rn=(m.floor(n/100)%10);    # Type of camber (3rd digit)
        if rn==0:
        #----------------------- STANDARD CAMBER ------------------------------
            #----------------------- CONSTANTS --------------------------------
            r=3.33333333333212*p**3+0.700000000000909*p**2+1.19666666666638*p-0.00399999999996247;    # R constant calculation by interpolation
            k1=1514933.33335235*p**4-1087744.00001147*p**3+286455.266669048*p**2-32968.4700001967*p+1420.18500000524;  # K1 constant calculation by interpolation
            #----------------------- CAMBER -----------------------------------
            for i in range(s):
                if x[i]<r:
                    y_c[i]=k1/6*(x[i]**3-3*r*x[i]**2+r**2*(3-r)*x[i]); # Mean camber y coordinate
                    dyc_dx[i]=k1/6*(3*x[i]**2-6*r*x[i]+r**2*(3-r)); # Mean camber first derivative
                else:
                    y_c[i]=k1*r**3/6*(1-x[i]);   # Mean camber y coordinate
                    dyc_dx[i]=-k1*r**3/(6);    # Mean camber first derivative
        elif rn==1:
        #----------------------- REFLEXED CAMBER ------------------------------
            #----------------------- CONSTANTS --------------------------------
            r=10.6666666666861*p**3-2.00000000001601*p**2+1.73333333333684*p-0.0340000000002413;  # R constant calculation by interponation
            k1=-27973.3333333385*p**3+17972.8000000027*p**2-3888.40666666711*p+289.076000000022;  # K1 constant calculation by interpolation
            k2_k1=85.5279999999984*p**3-34.9828000000004*p**2+4.80324000000028*p-0.21526000000003;    # K1/K2 constant calculation by interpolation
            #----------------------- CAMBER -----------------------------------
            for i in range(s):
                if x[i]<r:
                    y_c[i]=k1/6*((x[i]-r)**3-k2_k1*(1-r)**3*x[i]-r**3*x[i]+r**3);   # Mean camber y coordinate
                    dyc_dx[i]=k1/6*(3*(x[i]-r)**2-k2_k1*(1-r)**3-r**3);    # Mean camber first derivative
                else:
                    y_c[i]=k1/6*(k2_k1*(x[i]-r)**3-k2_k1*(1-r)**3*x[i]-r**3*x[i]+r**3); # Mean camber y coordinate
                    dyc_dx[i]=k1/6*(3*k2_k1*(x[i]-r)**2-k2_k1*(1-r)**3-r**3);  # Mean camber first derivative
    elif nc==6:
        #----------------------- MEAN CAMBER 6 DIGIT SERIES CALCULATION -----------
        #----------------------- CONSTANTS ------------------------------------
        ser=m.floor(n/100000);    # Number of series (1st digit)
        a=(m.floor(n/10000)%10)/10;  # Chordwise position of minimum pressure (2nd digit)
        c_li=(m.floor(n/100)%10)/10;  # Design lift coefficient (4th digit)
        g=-1./(1-a)*(a**2*(0.5*m.log(a)-0.25)+0.25);  # G constant calculation
        h=1./(1-a)*(0.5*(1-a)**2*m.log(1-a)-0.25*(1-a)**2)+g; # H constant calculation
        #----------------------- CAMBER ---------------------------------------
        indice = np.array([not ((x[i] == 0. or x[i] == 1. or x[i] == a)) for i in range(len(x))],dtype = bool);
        for i in range(len(x)):
            if indice[i]:
                y_c[i]=c_li/(2*m.pi*(a+1))*(1./(1-a)*(0.5*(a-x[i])**2*np.log(np.abs(a-x[i]))-0.5*(1-x[i])**2*np.log(1-x[i])+0.25*(1-x[i])**2-0.25*(a-x[i])**2)-x[i]*np.log(x[i])+g-h*x[i]); # Mean camber y coordinate
                dyc_dx[i]=-(c_li*(h+np.log(x[i])-(x[i]*0.5-a*0.5+(np.log(1-x[i])*(2*x[i]-2))*0.5+(np.log(np.abs(a-x[i]))*(2*a-2*x[i]))*0.5+(np.sign(a-x[i])*(a-x[i])**2)/(2*np.abs(a-x[i])))/(a-1)+1))/(2*m.pi*(a+1));    # Mean camber first derivative
    #----------------------- FINAL CALCULATIONS -------------------------------
    theta=np.arctan(dyc_dx); # Angle for modifying x coordinate
    #----------------------- COORDINATE ASSIGNATION ---------------------------
    x_e=x;  # X extrados coordinate
    x_i=np.flipud(x);  # X intrados coordinate
    y_e=(y_c+y_t*np.cos(theta));    # Y extrados coordinate
    y_i=np.flipud((y_c-y_t*np.cos(theta)));    # Y intrados coordinate
    X = np.concatenate([x_i,x_e[1:]]);    # from LE counter clocjwise
    Y = np.concatenate([y_i,y_e[1:]]);    # from LE counter clockwise
    #----------------------- Flaps deflection ----------------------------------------

    if xF != 0:
        X2 = np.concatenate([np.flipud(np.unique(np.concatenate([X[0:33],[1-xF]]))),np.unique((np.concatenate([X[34:],[1-xF]])))]);
        Y2 = np.concatenate([np.flipud(np.interp(np.flipud(X2[0:34]),np.flipud(X[0:33]),np.flipud(Y[0:33]))),np.interp(X2[34:],X[34:],Y[34:])]);
        index = np.where(X2 == 1-xF);
        index = index[0];
        indexIntra = index[0];
        indexExtra = index[1]
        middle = np.array([X2[indexIntra]+X2[indexExtra],0,Y2[indexIntra]+Y2[indexExtra]])/2;
        Rot = roty(delta);
        point = [X2[0:indexIntra],0,Y2[0:indexIntra]]-middle;
        point = np.dot(Rot,point)+middle;
        X2[0:indexIntra] = point[0];
        Y2[0:indexIntra] = point[2];
        point2 = [X2[indexExtra+1:],0,Y2[indexExtra+1:]]-middle;
        point2 = np.dot(Rot,point2)+middle;
        X2[indexExtra+1:] = point2[0];
        Y2[indexExtra+1:] = point2[2];
    else:
        X2 = X;
        Y2 = Y;

    section = np.array([X2,np.zeros(len(X2)),Y2]);
    return np.transpose(section);

def rotx(delta):
    """ Function that return the rotation matrix along the Y-axis for an angle of delta in degrees.
    INPUT:
        - delta : angle in degrees
    OUTPUT:
        - Rot : [3 x 3] rotation matrix (np.array([3,3]))"""
    deltaRad = m.pi*delta/180;
    return np.array([[1.,0.,0.],[0.,m.cos(deltaRad),-m.sin(deltaRad)],[0.,m.sin(deltaRad),m.cos(deltaRad)]]);

def roty(delta):
    """ Function that return the rotation matrix along the Y-axis for an angle of delta in degrees.
    INPUT:
        - delta : angle in degrees
    OUTPUT:
        - Rot : [3 x 3] rotation matrix (np.array([3,3]))"""
    deltaRad = m.pi*delta/180;
    return np.array([[m.cos(deltaRad),0,m.sin(deltaRad)],[0,1,0],[-m.sin(deltaRad),0,m.cos(deltaRad)]]);

def rotz(delta):
    """ Function that returns the rotation matrix along the Z-axis for an angle of delta in degrees.
    INPUT:
        - delta : angle in degrees
    OUTPUT:
        - Rot : [3 x 3] rotation matrix (np.array([3,3]))"""
    deltaRad = m.pi*delta/180;
    return np.array([[m.cos(deltaRad),-m.sin(deltaRad),0.],[m.sin(deltaRad),m.cos(deltaRad),0.],[0.,0.,1.]]);

def localVelTri(Vx,Vy,Vz,tp,sweep,dih):
    """ Fucntion that returns the local velocity triangle for any panel taking
    into account the twist, sweep and dihedre. 
    INPUT:
        - Vx,Vy,Vz : velocity triangle in the absolute referential
    OUTPUT:
        - Vxl, Vyl, Vzl : local velocity triangle for the panel"""
    
    Vxz = Vx * np.cos(sweep) - Vy * np.sin(sweep);
    Vyz = Vx * np.sin(sweep) + Vy * np.cos(sweep);
    Vzz = Vz;
    
    Vxx = Vxz * np.cos(tp) - Vzz * np.sin(tp);
    Vyx = Vyz;
    Vzx = Vxz * np.sin(tp) + Vzz * np.cos(tp);
    
    Vxl = Vxx;
    Vyl = Vyx * np.cos(dih) + Vzx * np.sin(dih);
    Vzl = - Vyx * np.sin(dih) + Vzx * np.cos(dih);
    return Vxl,Vyl,Vzl;

def localVelTriVT(Vx,Vy,Vz,sweep):
    """ Fucntion that returns the local velocity triangle for any panel taking
    into account the sweep for the vertical tail. 
    INPUT:
        - Vx,Vy,Vz : velocity triangle in the absolute referential
    OUTPUT:
        - Vxl, Vyl, Vzl : local velocity triangle for the panel"""
        
    Vxl = Vx * np.cos(sweep) - Vz * np.sin(sweep);
    Vyl = Vy;
    Vzl = Vx * np.sin(sweep) + Vz * np.cos(sweep);
    
    return Vxl,Vyl,Vzl;

def rangeF(a,b,step=1.):
    """ Cette fonction retourne la série d valeur comprises entre "a" et "b"
    avec comme intervalle entre valeurs "step"."""
    a = float(a);
    b = float(b);
    Delta = np.abs(b-a);
    l = int(Delta/step+1);
    R = np.zeros(l,dtype = float);
    for i in range(l):
        R[i] = a + step*i;
    return R;

def justLoad(fp,sk=0):
    try :
        af = np.loadtxt(fp,skiprows=sk);
    except ValueError:
        af = justLoad(fp,sk+1);
    except IOError:
        print('Document not found, verify the path and the existence of the file.');
        sys.exit();
    return af;