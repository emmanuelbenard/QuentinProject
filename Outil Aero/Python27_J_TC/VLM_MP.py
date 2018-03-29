  # -*- coding: utf-8 -*-
import math as m
import numpy as np
import utilitaire as u
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import test

def ICMatrix(ac,cla,flow):
    """ Prediction of aerodynamic characteristics of the wing
       Autor : Quentin borlon 
       Date : 5 mai 2017
    
       Function that predicts the aerodynamic coefficients for a given wing.
       Based on the wing geometry and the sectional 2D aerodynamic datas.
    
       Function initially based, and notation preserved from :
           2004  Mihai Pruna, Alberto Davila;
           free software; under the terms of the GNU General Public License
          
     INPUT:
       clAlpha : vertical array with clAlphas(i) is the lift curve slope of the
           panel from wing.y(i) to wing.y(i+1);
       wing : a structral object with as fields:
    
           b : span
           chord : vertical array with the chord at the root (1) any
               discontinuity of taper ratio (2:end-1) and at the tip (end);
           flapsDiscY : vertical array with the spanwise coordinate of the flaps
               discontinuity sections
           afDiscY : vertical array with the spanwise coordinate of the
               airfoil discontinuity sections
           airfoil : a cell-array with each cell gives the airfoil naca number
               representation, cell 1 correspond to first panel after root.
           sweep : vertical array with wing.sweep(i) is the sweep angle of 
               the panel from wing.y(i) to wing.y(i+1) (rad)
           dih : vertical array with wing.dih(i) is the dihedral angle of 
               the panel from wing.y(i) to wing.y(i+1) (rad)
           twist : vertical array with wing.twist(i) is the twist angle of 
               the section at wing.y(i) (rad)
           deltasFlaps : vertical array with wing.deltasFlaps(i) is the 
               flaps defection  of the panel from wing.y(i) to wing.y(i+1)
               (deg)
           r : number of spanwise panel along the wing;
           m : number of chordwise panel along the airfoil;
           Mach : flight mach number
           cFlaps_cLoc : vertical array with wing.cFlaps_cLocs(i) is the 
               local flaps to chord ratio
           y : the  spanwise location of (-b/2 -> b/2) the limits of the panels
           discY : vertical array of the complete set of the spanwise location
           airfoilIndex : vertical array with wing.airfoilIndex(i) is the index of 
               the airfoil (wing.airfoil) to use for the section at wing.y(i)
           chordDistrib : vertical array with wing.chordDistrib(i) is the chord length of 
               the section at wing.y(i)
     
     OUTPUT:
       A : the influence coefficient matrix [n x n] such that A*{GAMMA/2} + {Q}*{normal} = 0
       normal : a [3 x (wing.getR()/2+1)] matrix that provides the normal downward
           of the panel."""        

    prop = ac.prop;
    wing = ac.wing;
    cf = wing.getCF();
    # Generate grid coordinates
    # Generate collocation points and normal : where tangency condition is
    # satisfied. Distance from bound vortex depends on the sectional lift
    # curve slope : (dist/localChord) = clAlphas/(4*pi), clAlphas assumed to be 2 *pi
    if prop.bool and cf != 0.:
        return getGridF_Engines(flow,ac,cla);
    elif prop.bool:
        return getGrid_Engines(flow,ac,cla);
    elif cf !=0.:
        return getGridF_NOEngines(flow,ac,cla);
    else:
        return getGrid_NOEngines(flow,ac,cla);

def getGrid_NOEngines(flow,ac,cla):
    # flow sideslip and aoa angles
    beta = - flow.getBeta()*m.pi/180.;
    aoa = m.pi * (flow.getAMax()+flow.getAMin())/360.;
    
    # Main lifting surfaces
    wing = ac.wing;
    htail = ac.htail;
    vtail = ac.vtail;
    # Numerical parameters for  discretization
    
    mC = wing.mC;                                                               # chordwise discretisation number of control point for the chord
    mW = 8;                                # discretisation of the wake, get correct direction of it behind parts
    n = wing.getR()+htail.getR()+vtail.getR();                                               # spanwise discretisation number of panel
    
    # Recover the wing parameters
    
    
    
    
    X,Y,Z,COLOCX,COLOCY,COLOCZ,ds,dS,normal = test.HorSurfaces_NOEngines(wing,cla[:wing.r],aoa,beta,mC,mW);
    if htail.bool:
        XHT,YHT,ZHT,COLOCXHT,COLOCYHT,COLOCZHT,dsHT,dSHT,normalHT = test.HorSurfaces_NOEngines(htail,cla[wing.r:wing.r+htail.r],aoa,beta,mC,mW);
        X = np.concatenate([X,XHT]);
        Y = np.concatenate([Y,YHT]);
        Z = np.concatenate([Z,ZHT]);
        COLOCX = np.concatenate([COLOCX,COLOCXHT]);
        COLOCY = np.concatenate([COLOCY,COLOCYHT]);
        COLOCZ = np.concatenate([COLOCZ,COLOCZHT]);
        ds = np.concatenate([ds,dsHT]);
        dS = np.concatenate([dS,dSHT]);
        normal = np.concatenate([normal,normalHT],axis = 1);
    if vtail.bool:
        XVT,YVT,ZVT,COLOCXVT,COLOCYVT,COLOCZVT,dsVT,dSVT,normalVT = test.VerSurfaces_NOEngines(vtail,cla[wing.r+htail.r:],aoa,beta,mC,mW);
        X = np.concatenate([X,XVT]);
        Y = np.concatenate([Y,YVT]);
        Z = np.concatenate([Z,ZVT]);
        COLOCX = np.concatenate([COLOCX,COLOCXVT]);
        COLOCY = np.concatenate([COLOCY,COLOCYVT]);
        COLOCZ = np.concatenate([COLOCZ,COLOCZVT]);
        ds = np.concatenate([ds,dsVT]);
        dS = np.concatenate([dS,dSVT]);
        normal = np.concatenate([normal,normalVT],axis = 1);

    select = np.zeros([n,n * (mC-1)]); # rechercher intensité du dernier vortex uniquement
    select2 = np.zeros([n * (mC-1),n]); # pour chaque paneau sur même section y, même velocity triangle
    select3 = np.zeros([n + len(ac.prop.D),n * (mC-1) + len(ac.prop.D)]); # 
    for i in range(n):
        select[i,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i] = 1.;
        select3[i,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i];
    if ac.prop.bool:
        select3[-len(ac.prop.D):,-len(ac.prop.D):] = np.eye(len(ac.prop.D));
    Ao,Vxo,Vyo,Vzo = ICM(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac,n,mC,mW);
    invA = np.linalg.inv(Ao);
    A = invA;
    Vx = np.dot(select3,Vxo);
    Vy = np.dot(select3,Vyo);
    Vz = np.dot(select3,Vzo);
    return A,normal,Vx,Vy,Vz,select,select2;


def getGrid_Engines(flow,ac,cla):
    # flow sideslip and aoa angles
    beta = - flow.getBeta()*m.pi/180.;
    aoa = m.pi * (flow.getAMax()+flow.getAMin())/360.;
    # Main lifting surfaces
    wing = ac.wing;
    htail = ac.htail;
    vtail = ac.vtail;
    prop = ac.prop;
    
    rho0=1.225; #masse volumique à niveau de la mer           [kg/m^3]
    dT=-6.5;    #gradiente de temperature dans la troposphere [K/km]
    T0=288.15;  #Temperature à niveau de la mer               [K]
    g=9.80665;  #gravité                                      [m/s^2]
    Rair=287.1;   #Constante de l'air                           [m^2/(s^2*K)]
    h = flow.getH();    # flight altitude [km]
    V0 = 1.; # freestream velocity [m/s]
    rho = rho0 * (1. + dT*h/T0)**(- g/(Rair*dT*10**(-3)) - 1.); # air density
    Sh = m.pi * prop.getD()**2 *0.25;#Surface disque actuator          [m^2]
    nbE = len(prop.getD());
    OWU = prop.getYp()/np.abs(prop.getYp());
    for ii in range(nbE):
        if not(prop.OWU[ii]):
            OWU[ii] *= -1.;
    # Numerical parameters for  discretization
    lF = 2.*wing.getSpan();                                                     # temps caractéristique de convection des vortex
    nbEch = 3.*m.ceil(max(1./(prop.J*prop.D)));                                    # Nombre minimal de points de controle par rotation des lignes de courant/tourbillons
    mW = int(lF*nbEch);                                                         # discretisation of the wake, get correct direction of it behind parts
    times = np.linspace(0.,lF,mW);
    mC = wing.mC;                                                               # chordwise discretisation number of control point for the chord
    n = wing.getR()+htail.getR()+vtail.getR();                                  # spanwise discretisation number of panel
    
    X,Y,Z,COLOCX,COLOCY,COLOCZ,ds,dS,normal = test.HorSurfaces_Engines(wing,cla[:wing.r],prop,Sh,OWU,nbE,V0,rho,aoa,beta,mC,mW,times);
    if htail.bool:
        XHT,YHT,ZHT,COLOCXHT,COLOCYHT,COLOCZHT,dsHT,dSHT,normalHT = test.HorSurfaces_Engines(htail,cla[wing.r:wing.r+htail.r],prop,Sh,OWU,nbE,V0,rho,aoa,beta,mC,mW,times);
        X = np.concatenate([X,XHT]);
        Y = np.concatenate([Y,YHT]);
        Z = np.concatenate([Z,ZHT]);
        COLOCX = np.concatenate([COLOCX,COLOCXHT]);
        COLOCY = np.concatenate([COLOCY,COLOCYHT]);
        COLOCZ = np.concatenate([COLOCZ,COLOCZHT]);
        ds = np.concatenate([ds,dsHT]);
        dS = np.concatenate([dS,dSHT]);
        normal = np.concatenate([normal,normalHT],axis = 1);
    if vtail.bool:
        XVT,YVT,ZVT,COLOCXVT,COLOCYVT,COLOCZVT,dsVT,dSVT,normalVT = test.VerSurfaces_Engines(vtail,cla[wing.r+htail.r:],prop,Sh,OWU,nbE,V0,rho,aoa,beta,mC,mW,times);
        X = np.concatenate([X,XVT]);
        Y = np.concatenate([Y,YVT]);
        Z = np.concatenate([Z,ZVT]);
        COLOCX = np.concatenate([COLOCX,COLOCXVT]);
        COLOCY = np.concatenate([COLOCY,COLOCYVT]);
        COLOCZ = np.concatenate([COLOCZ,COLOCZVT]);
        ds = np.concatenate([ds,dsVT]);
        dS = np.concatenate([dS,dSVT]);
        normal = np.concatenate([normal,normalVT],axis = 1);
    select = np.zeros([n,n * (mC-1)]); # rechercher intensité du dernier vortex uniquement
    select2 = np.zeros([n * (mC-1),n]); # pour chaque paneau sur même section y, même velocity triangle
    select3 = np.zeros([n + len(ac.prop.D),n * (mC-1) + len(ac.prop.D)]); # 
    for i in range(n):
        select[i,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i] = 1.;
        select3[i,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i];
    if ac.prop.bool:
        select3[n:,(mC-1)*n:] = np.eye(len(ac.prop.D));
    Ao,Vxo,Vyo,Vzo = ICM(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac,n,mC,mW);
    invA = np.linalg.inv(Ao);
    A = invA;
    Vx = np.dot(select3,Vxo);
    Vy = np.dot(select3,Vyo);
    Vz = np.dot(select3,Vzo);
    return A,normal,Vx,Vy,Vz,select,select2;
        
def getGridF_NOEngines(flow,ac,cla):
    # flow sideslip and aoa angles
    beta = - flow.getBeta()*m.pi/180.;
    aoa = m.pi * (flow.getAMax()+flow.getAMin())/360.;
    # Main lifting surfaces
    wing = ac.wing;
    htail = ac.htail;
    vtail = ac.vtail;
    
    # Numerical parameters for  discretization
    
    mC = wing.mC;                                                               # chordwise discretisation number of control point for the chord
    mW = flow.mW;                                                               # discretisation of the wake, get correct direction of it behind parts
    n = 2*wing.getR()+htail.getR() + vtail.getR();                              # spanwise discretisation number of panel# Recover the wing parameters
    
    X,Y,Z,COLOCX,COLOCY,COLOCZ,ds,dS,normal = test.HorSurfacesF_NOEngines(wing,cla[:wing.r],aoa,beta,mC,mW);
    if htail.bool:
        XHT,YHT,ZHT,COLOCXHT,COLOCYHT,COLOCZHT,dsHT,dSHT,normalHT = test.HorSurfaces_NOEngines(htail,cla[wing.r:wing.r+htail.r],aoa,beta,mC,mW);
        X = np.concatenate([X,XHT]);
        Y = np.concatenate([Y,YHT]);
        Z = np.concatenate([Z,ZHT]);
        COLOCX = np.concatenate([COLOCX,COLOCXHT]);
        COLOCY = np.concatenate([COLOCY,COLOCYHT]);
        COLOCZ = np.concatenate([COLOCZ,COLOCZHT]);
        ds = np.concatenate([ds,dsHT]);
        dS = np.concatenate([dS,dSHT]);
        normal = np.concatenate([normal,normalHT],axis = 1);
    if vtail.bool:
        XVT,YVT,ZVT,COLOCXVT,COLOCYVT,COLOCZVT,dsVT,dSVT,normalVT = test.VerSurfaces_NOEngines(vtail,cla[wing.r+htail.r:],aoa,beta,mC,mW);
        X = np.concatenate([X,XVT]);
        Y = np.concatenate([Y,YVT]);
        Z = np.concatenate([Z,ZVT]);
        COLOCX = np.concatenate([COLOCX,COLOCXVT]);
        COLOCY = np.concatenate([COLOCY,COLOCYVT]);
        COLOCZ = np.concatenate([COLOCZ,COLOCZVT]);
        ds = np.concatenate([ds,dsVT]);
        dS = np.concatenate([dS,dSVT]);
        normal = np.concatenate([normal,normalVT],axis = 1);
        
#    plt.plot(Y[:wing.r*(2*(mC+mW)+1)],X[:wing.r*(2*(mC+mW)+1)]);
#    plt.plot(Y[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)],X[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)]);
#    plt.plot(Y[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)],X[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)]);
#    plt.plot(Y[(2*wing.r+htail.r)*(2*(mC+mW)+1):],X[(2*wing.r+htail.r)*(2*(mC+mW)+1):]);
#    plt.axis([-7,7,-1,13])
#    plt.show()
#    
#    plt.plot(Y[:wing.r*(2*(mC+mW)+1)],Z[:wing.r*(2*(mC+mW)+1)]);
#    plt.plot(Y[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)],Z[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)]);
#    plt.plot(Y[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)],Z[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)]);
#    plt.plot(Y[(2*wing.r+htail.r)*(2*(mC+mW)+1):],Z[(2*wing.r+htail.r)*(2*(mC+mW)+1):]);
#    plt.axis([-7,7,-3,7])
#    plt.show()
#    
#    plt.plot(X[:wing.r*(2*(mC+mW)+1)],Z[:wing.r*(2*(mC+mW)+1)]);
#    plt.plot(X[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)],Z[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)]);
#    plt.plot(X[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)],Z[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)]);
#    plt.plot(X[(2*wing.r+htail.r)*(2*(mC+mW)+1):],Z[(2*wing.r+htail.r)*(2*(mC+mW)+1):]);
#    plt.axis([-1,13,-3,7])
#    plt.show()
#    return
    
    
    select = np.zeros([wing.r + htail.r + vtail.r,n * (mC-1)]); # rechercher intensité du dernier vortex uniquement
    select2 = np.zeros([n * (mC-1),wing.r + htail.r + vtail.r]); # pour chaque paneau sur même section y, même velocity triangle
    select3 = np.zeros([wing.r + htail.r + vtail.r + len(ac.prop.D),n * (mC-1) + len(ac.prop.D)]); # 
    for i in range(wing.r):
        select[i,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i] = 1.;
        select3[i,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i];
    for i in range(wing.r,n):
        select[i-wing.r,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i - wing.r] = 1.;
        select3[i - wing.r,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i-wing.r];
    if ac.prop.bool:
        select3[wing.r + htail.r + vtail.r:,n * (mC-1):] = np.eye(len(ac.prop.D));
    Ao,Vxo,Vyo,Vzo = ICM_F(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac,n,mC,mW);
    invA = np.linalg.inv(Ao);
    A = invA;
    Vx = np.dot(select3,Vxo);
    Vy = np.dot(select3,Vyo);
    Vz = np.dot(select3,Vzo);
    return A,normal,Vx,Vy,Vz,select,select2;


def getGridF_Engines(flow,ac,cla):
    # flow sideslip and aoa angles
    beta = - flow.getBeta()*m.pi/180.;
    aoa = m.pi * (flow.getAMax()+flow.getAMin())/360.;
    # Main lifting surfaces
    wing = ac.wing;
    htail = ac.htail;
    vtail = ac.vtail;
    prop = ac.prop;
    
    rho0=1.225; #masse volumique à niveau de la mer           [kg/m^3]
    dT=-6.5;    #gradiente de temperature dans la troposphere [K/km]
    T0=288.15;  #Temperature à niveau de la mer               [K]
    g=9.80665;  #gravité                                      [m/s^2]
    Rair=287.1;   #Constante de l'air                           [m^2/(s^2*K)]
    h = flow.getH();    # flight altitude [km]
    V0 = 1.; # freestream velocity [m/s]
    rho = rho0 * (1. + dT*h/T0)**(- g/(Rair*dT*10**(-3)) - 1.); # air density
    Sh = m.pi * prop.getD()**2 *0.25;#Surface disque actuator          [m^2]
    nbE = len(prop.getD());
    OWU = prop.getYp()/np.abs(prop.getYp());
    for ii in range(nbE):
        if not(prop.OWU[ii]):
            OWU[ii] *= -1.;
    # Numerical parameters for  discretization
    lF = 2.*wing.getSpan();                                                     # temps caractéristique de convection des vortex
    nbEch = 3.*m.ceil(max(1./(prop.J*prop.D)));                                    # Nombre minimal de points de controle par rotation des lignes de courant/tourbillons
    mW = int(lF*nbEch);                                                         # discretisation of the wake, get correct direction of it behind parts
    times = np.linspace(0.,lF,mW);
    mC = wing.mC;                                                               # chordwise discretisation number of control point for the chord
    n = 2*wing.getR()+htail.getR()+vtail.getR();                                  # spanwise discretisation number of panel
    
    X,Y,Z,COLOCX,COLOCY,COLOCZ,ds,dS,normal = test.HorSurfacesF_Engines(wing,cla[:wing.r],prop,Sh,OWU,nbE,V0,rho,aoa,beta,mC,mW,times);
    if htail.bool:
        XHT,YHT,ZHT,COLOCXHT,COLOCYHT,COLOCZHT,dsHT,dSHT,normalHT = test.HorSurfaces_Engines(htail,cla[wing.r:wing.r+htail.r],prop,Sh,OWU,nbE,V0,rho,aoa,beta,mC,mW,times);
        X = np.concatenate([X,XHT]);
        Y = np.concatenate([Y,YHT]);
        Z = np.concatenate([Z,ZHT]);
        COLOCX = np.concatenate([COLOCX,COLOCXHT]);
        COLOCY = np.concatenate([COLOCY,COLOCYHT]);
        COLOCZ = np.concatenate([COLOCZ,COLOCZHT]);
        ds = np.concatenate([ds,dsHT]);
        dS = np.concatenate([dS,dSHT]);
        normal = np.concatenate([normal,normalHT],axis = 1);
    if vtail.bool:
        XVT,YVT,ZVT,COLOCXVT,COLOCYVT,COLOCZVT,dsVT,dSVT,normalVT = test.VerSurfaces_Engines(vtail,cla[wing.r+htail.r:],prop,Sh,OWU,nbE,V0,rho,aoa,beta,mC,mW,times);
        X = np.concatenate([X,XVT]);
        Y = np.concatenate([Y,YVT]);
        Z = np.concatenate([Z,ZVT]);
        COLOCX = np.concatenate([COLOCX,COLOCXVT]);
        COLOCY = np.concatenate([COLOCY,COLOCYVT]);
        COLOCZ = np.concatenate([COLOCZ,COLOCZVT]);
        ds = np.concatenate([ds,dsVT]);
        dS = np.concatenate([dS,dSVT]);
        normal = np.concatenate([normal,normalVT],axis = 1);
#    plt.plot(Y[:wing.r*(2*(mC+mW)+1)],X[:wing.r*(2*(mC+mW)+1)]);
#    plt.plot(Y[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)],X[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)]);
#    plt.plot(Y[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)],X[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)]);
#    plt.plot(Y[(2*wing.r+htail.r)*(2*(mC+mW)+1):],X[(2*wing.r+htail.r)*(2*(mC+mW)+1):]);
#    plt.axis([-7,7,-1,13])
#    plt.show()
#    
#    plt.plot(Y[:wing.r*(2*(mC+mW)+1)],Z[:wing.r*(2*(mC+mW)+1)]);
#    plt.plot(Y[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)],Z[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)]);
#    plt.plot(Y[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)],Z[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)]);
#    plt.plot(Y[(2*wing.r+htail.r)*(2*(mC+mW)+1):],Z[(2*wing.r+htail.r)*(2*(mC+mW)+1):]);
#    plt.axis([-7,7,-3,7])
#    plt.show()
#    
#    plt.plot(X[:wing.r*(2*(mC+mW)+1)],Z[:wing.r*(2*(mC+mW)+1)]);
#    plt.plot(X[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)],Z[wing.r*(2*(mC+mW)+1):2*wing.r*(2*(mC+mW)+1)]);
#    plt.plot(X[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)],Z[2*wing.r*(2*(mC+mW)+1):(2*wing.r+htail.r)*(2*(mC+mW)+1)]);
#    plt.plot(X[(2*wing.r+htail.r)*(2*(mC+mW)+1):],Z[(2*wing.r+htail.r)*(2*(mC+mW)+1):]);
#    plt.axis([-1,13,-3,7])
#    plt.show()
#    return
    select = np.zeros([wing.r + htail.r + vtail.r,n * (mC-1)]); # rechercher intensité du dernier vortex uniquement
    select2 = np.zeros([n * (mC-1),wing.r + htail.r + vtail.r]); # pour chaque paneau sur même section y, même velocity triangle
    select3 = np.zeros([wing.r + htail.r + +vtail.r + len(ac.prop.D),n * (mC-1) + len(ac.prop.D)]); # 
    for i in range(wing.r):
        select[i,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i] = 1.;
        select3[i,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i];
    for i in range(wing.r,n):
        select[i-wing.r,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i - wing.r] = 1.;
        select3[i - wing.r,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i-wing.r];
    if ac.prop.bool:
        select3[wing.r + htail.r + +vtail.r:,n * (mC-1):] = np.eye(len(ac.prop.D));
    Ao,Vxo,Vyo,Vzo = ICM_F(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac,n,mC,mW);
    invA = np.linalg.inv(Ao);
    A = invA;
    Vx = np.dot(select3,Vxo);
    Vy = np.dot(select3,Vyo);
    Vz = np.dot(select3,Vzo);
    return A,normal,Vx,Vy,Vz,select,select2;

def ICM(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac,n,mC,mW):
    if ac.fus.bool:
        HWing = ac.fus.vDist > 0;
    if ac.htail.bool and ac.vtail.bool:
        HTail = ac.htail.z[ac.htail.getR()/2] > ((ac.vtail.z[-1]-ac.vtail.z[0]) * 0.66) + ac.vtail.z[0];
    if not(ac.fus.bool):
        if (not(ac.vtail.bool) or not(ac.htail.bool) or HTail) :
            A,Vx,Vy,Vz = OnlyWing(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
#        elif not HWing:
#            A,Vx,Vy,Vz = halfWing_Wall(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
        else:
            A,Vx,Vy,Vz = BothWingOneTailVtail(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
    else:
        if not(ac.vtail.bool):
            if HWing:
                A,Vx,Vy,Vz = OnlyWing(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            else:
                A,Vx,Vy,Vz = OneWingBothTail(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
        else:
            if HWing and HTail:
                A,Vx,Vy,Vz = OnlyWing(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            elif HTail:
                A,Vx,Vy,Vz = OneWingBothTail(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            elif HWing:
                A,Vx,Vy,Vz = BothWingOneTailVtail(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            else:
                A,Vx,Vy,Vz = OneWingOneTailVtail(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
    return A,Vx,Vy,Vz;

def ICM_F(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac,n,mC,mW):
    if ac.fus.bool:
        HWing = ac.fus.vDist > 0;
    if ac.htail.bool and ac.vtail.bool:
        HTail = ac.htail.z[ac.htail.getR()/2] > ((ac.vtail.z[-1]-ac.vtail.z[0]) * 0.66) + ac.vtail.z[0];
    if not(ac.fus.bool):
        if (not(ac.vtail.bool) or not(ac.htail.bool) or HTail):
            A,Vx,Vy,Vz = OnlyWing_F(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
        else:
            A,Vx,Vy,Vz = BothWingOneTailVtailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
    else:
        if not(ac.vtail.bool):
            if HWing:
                A,Vx,Vy,Vz = OnlyWing_F(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            else:
                A,Vx,Vy,Vz = OneWingBothTailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
        else:
            if HWing and HTail:
                A,Vx,Vy,Vz = OnlyWing_F(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            elif HTail:
                A,Vx,Vy,Vz = OneWingBothTailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            elif HWing:
                A,Vx,Vy,Vz = BothWingOneTailVtailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            else:
                A,Vx,Vy,Vz = OneWingOneTailVtailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
    return A,Vx,Vy,Vz;

def halfWing_Wall(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    m = n * (mC-1);
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vx = np.zeros([m,n*(mC-1)],dtype = float);
    Vy = np.zeros([m,n*(mC-1)],dtype = float);
    Vz = np.zeros([m,n*(mC-1)],dtype = float);
    # Influence de la partie gauche de l'aile sur elle-même et les autres surfaces
    index = np.array(range(ac.wing.r/2*(mC-1)));
    for j in range(ac.wing.r/2-1):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    j = ac.wing.r/2-1;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NR_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie droite de l'aile sur elle-même et les autres surfaces
    index = np.array(range(ac.wing.r/2*(mC-1),ac.wing.r*(mC-1)));
    j = ac.wing.r/2;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    for j in range(ac.wing.r/2+1,ac.wing.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    if ac.prop.bool:
        test.On_PropInfluence(ac,X,Y,Z,n,m,mC,mW,A,Vx,Vy,Vz);
    return A,Vx,Vy,Vz;

def OnlyWing(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    m = n * (mC-1);
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vx = np.zeros([m,n*(mC-1)],dtype = float);
    Vy = np.zeros([m,n*(mC-1)],dtype = float);
    Vz = np.zeros([m,n*(mC-1)],dtype = float);
    # Influence de l'aile sur tous les points
    index = np.array(range(n*(mC-1)));
    for j in range(ac.wing.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la HTail uniquement sur aile et elle-même
    index = np.array(range((ac.wing.r+ac.htail.r)*(mC-1)));
    for j in range(ac.wing.r,ac.wing.r+ac.htail.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la VTail sur tous les points
    index = np.array(range(n*(mC-1)));
    if ac.vtail.bool:
        j = ac.wing.r+ac.htail.r;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
        for j in range(ac.wing.r+ac.htail.r+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    if ac.prop.bool:
        test.On_PropInfluence(ac,X,Y,Z,n,m,mC,mW,A,Vx,Vy,Vz);
    return A,Vx,Vy,Vz;

def OnlyWing_F(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    m = n * (mC-1);
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vx = np.zeros([m,n*(mC-1)],dtype = float);
    Vy = np.zeros([m,n*(mC-1)],dtype = float);
    Vz = np.zeros([m,n*(mC-1)],dtype = float);
    # Influence de l'aile sur tous les points
    index = np.array(range(n*(mC-1)));
    for j in range(ac.wing.r*2):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la HTail uniquement sur aile, flaps et elle-même
    index = np.array(range((2*ac.wing.r+ac.htail.r)*(mC-1)));
    for j in range(2*ac.wing.r,2*ac.wing.r+ac.htail.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la VTail sur tous les points
    index = np.array(range(n*(mC-1)));
    if ac.vtail.bool:
        j = 2*ac.wing.r+ac.htail.r;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
        for j in range(2*ac.wing.r+ac.htail.r+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    if ac.prop.bool:
        test.On_PropInfluenceF(ac,X,Y,Z,n,m,mC,mW,A,Vx,Vy,Vz);
    return A,Vx,Vy,Vz;

def BothWingOneTailVtail(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    m = n * (mC-1);
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vx = np.zeros([m,n*(mC-1)],dtype = float);
    Vy = np.zeros([m,n*(mC-1)],dtype = float);
    Vz = np.zeros([m,n*(mC-1)],dtype = float);
    # Influence de l'aile sur tous les points
    index = np.array(range(n*(mC-1)));
    for j in range(ac.wing.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie gauche de la HTail uniquement sur aile et elle-même
    index = np.array(range((ac.wing.r+ac.htail.r/2)*(mC-1)));
    for j in range(ac.wing.r,ac.wing.r+ac.htail.r/2-1):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    j = ac.wing.r+ac.htail.r/2-1;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NR_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie droite de la HTail uniquement sur aile et elle-même
    index = np.concatenate((range(ac.wing.r*(mC-1)),range((ac.wing.r+ac.htail.r/2)*(mC-1),(ac.wing.r+ac.htail.r)*(mC-1))));
    j = ac.wing.r+ac.htail.r/2;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    for j in range(ac.wing.r+ac.htail.r/2+1,ac.wing.r+ac.htail.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la VTail sur tous les points
    index = np.array(range(n*(mC-1)));
    if ac.vtail.bool:
        j = ac.wing.r+ac.htail.r;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
        for j in range(ac.wing.r+ac.htail.r+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    if ac.prop.bool:
        test.On_PropInfluence(ac,X,Y,Z,n,m,mC,mW,A,Vx,Vy,Vz);
    return A,Vx,Vy,Vz;

def BothWingOneTailVtailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    m = n * (mC-1);
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vx = np.zeros([m,n*(mC-1)],dtype = float);
    Vy = np.zeros([m,n*(mC-1)],dtype = float);
    Vz = np.zeros([m,n*(mC-1)],dtype = float);
    # Influence de l'aile sur tous les points
    index = np.array(range(n*(mC-1)));
    for j in range(2*ac.wing.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie gauche de la HTail uniquement sur aile et elle-même
    index = np.array(range((2*ac.wing.r+ac.htail.r/2)*(mC-1)));
    for j in range(2*ac.wing.r,2*ac.wing.r+ac.htail.r/2-1):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    j = 2*ac.wing.r+ac.htail.r/2-1;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NR_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie droite de la HTail uniquement sur aile et elle-même
    index = np.concatenate((range(2*ac.wing.r*(mC-1)),range((2*ac.wing.r+ac.htail.r/2)*(mC-1),(2*ac.wing.r+ac.htail.r)*(mC-1))));
    j = 2*ac.wing.r+ac.htail.r/2;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    for j in range(2*ac.wing.r+ac.htail.r/2+1,2*ac.wing.r+ac.htail.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la VTail sur tous les points
    index = np.array(range(n*(mC-1)));
    if ac.vtail.bool:
        j = 2*ac.wing.r+ac.htail.r;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
        for j in range(2*ac.wing.r+ac.htail.r+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    if ac.prop.bool:
        test.On_PropInfluenceF(ac,X,Y,Z,n,m,mC,mW,A,Vx,Vy,Vz);
    return A,Vx,Vy,Vz;

def OneWingBothTail(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    m = n * (mC-1);
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vx = np.zeros([m,n*(mC-1)],dtype = float);
    Vy = np.zeros([m,n*(mC-1)],dtype = float);
    Vz = np.zeros([m,n*(mC-1)],dtype = float);
    # Influence de la partie gauche de l'aile sur elle-même et les autres surfaces
    index = np.concatenate((range(ac.wing.r/2*(mC-1)),range(ac.wing.r*(mC-1),n*(mC-1))));
    for j in range(ac.wing.r/2-1):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    j = ac.wing.r/2-1;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NR_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie droite de l'aile sur elle-même et les autres surfaces
    index = np.array(range(ac.wing.r/2*(mC-1),n*(mC-1)));
    j = ac.wing.r/2;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    for j in range(ac.wing.r/2+1,ac.wing.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la HTail uniquement sur aile et elle-même
    index = np.array(range((ac.wing.r+ac.htail.r)*(mC-1)));
    for j in range(ac.wing.r,ac.wing.r+ac.htail.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la VTail sur tous les points
    index = np.array(range(n*(mC-1)));
    if ac.vtail.bool:
        j = ac.wing.r+ac.htail.r;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
        for j in range(ac.wing.r+ac.htail.r+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    if ac.prop.bool:
        test.On_PropInfluence(ac,X,Y,Z,n,m,mC,mW,A,Vx,Vy,Vz);
    return A,Vx,Vy,Vz;

def OneWingBothTailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    m = n * (mC-1);
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vx = np.zeros([m,n*(mC-1)],dtype = float);
    Vy = np.zeros([m,n*(mC-1)],dtype = float);
    Vz = np.zeros([m,n*(mC-1)],dtype = float);
    # Influence de la partie gauche de l'aile sur elle-même, le volet gauche et les autres surfaces
    index = np.concatenate((range(ac.wing.r/2*(mC-1)),range(ac.wing.r*(mC-1),3*ac.wing.r/2*(mC-1)),range(ac.wing.r*2*(mC-1),n*(mC-1))));
    for j in range(ac.wing.r/2-1):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    j = ac.wing.r/2-1;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NR_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie gauche du volet sur lui-même, l'aile gauche et les autres surfaces
    for j in range(ac.wing.r,ac.wing.r*3/2-1):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    j = ac.wing.r*3/2-1;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NR_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie droite de l'aile sur elle-même et les autres surfaces
    index = np.concatenate((range(ac.wing.r/2*(mC-1),ac.wing.r*(mC-1)),range(3*ac.wing.r/2*(mC-1),n*(mC-1))));
    j = ac.wing.r/2;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    for j in range(ac.wing.r/2+1,ac.wing.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie droite du volet sur lui-même, l'aile droite et les autres surfaces
    j = 3*ac.wing.r/2;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    for j in range(3*ac.wing.r/2+1,2*ac.wing.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la HTail uniquement sur aile et elle-même
    index = np.array(range((2*ac.wing.r+ac.htail.r)*(mC-1)));
    for j in range(2*ac.wing.r,ac.wing.r+ac.htail.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la VTail sur tous les points
    index = np.array(range(n*(mC-1)));
    if ac.vtail.bool:
        j = 2*ac.wing.r+ac.htail.r;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
        for j in range(2*ac.wing.r+ac.htail.r+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    if ac.prop.bool:
        test.On_PropInfluenceF(ac,X,Y,Z,n,m,mC,mW,A,Vx,Vy,Vz);
    return A,Vx,Vy,Vz;

def OneWingOneTailVtail(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    m = n * (mC-1);
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vx = np.zeros([m,n*(mC-1)],dtype = float);
    Vy = np.zeros([m,n*(mC-1)],dtype = float);
    Vz = np.zeros([m,n*(mC-1)],dtype = float);
    # Influence de la partie gauche de l'aile sur elle-même et les autres surfaces
    index = np.concatenate((range(ac.wing.r/2*(mC-1)),range(ac.wing.r*(mC-1),n*(mC-1))));
    for j in range(ac.wing.r/2-1):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    j = ac.wing.r/2-1;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NR_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie droite de l'aile sur elle-même et les autres surfaces
    index = np.array(range(ac.wing.r/2*(mC-1),n*(mC-1)));
    j = ac.wing.r/2;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    for j in range(ac.wing.r/2+1,ac.wing.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie gauche de la HTail uniquement sur aile et elle-même
    index = np.array(range((ac.wing.r+ac.htail.r/2)*(mC-1)));
    for j in range(ac.wing.r,ac.wing.r+ac.htail.r/2-1):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    j = ac.wing.r+ac.htail.r/2-1;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NR_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie droite de la HTail uniquement sur aile et elle-même
    index = np.concatenate((range(ac.wing.r*(mC-1)),range((ac.wing.r+ac.htail.r/2)*(mC-1),(ac.wing.r+ac.htail.r)*(mC-1))));
    j = ac.wing.r+ac.htail.r/2;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    for j in range(ac.wing.r+ac.htail.r/2+1,ac.wing.r+ac.htail.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la VTail sur tous les points
    index = np.array(range(n*(mC-1)));
    if ac.vtail.bool:
        j = ac.wing.r+ac.htail.r;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
        for j in range(ac.wing.r+ac.htail.r+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    if ac.prop.bool:
        test.On_PropInfluence(ac,X,Y,Z,n,m,mC,mW,A,Vx,Vy,Vz);
    return A,Vx,Vy,Vz;

def OneWingOneTailVtailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    m = n * (mC-1);
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vx = np.zeros([m,n*(mC-1)],dtype = float);
    Vy = np.zeros([m,n*(mC-1)],dtype = float);
    Vz = np.zeros([m,n*(mC-1)],dtype = float);
    # Influence de la partie gauche de l'aile sur elle-même, le volet gauche et les autres surfaces
    index = np.concatenate((range(ac.wing.r/2*(mC-1)),range(ac.wing.r*(mC-1),3*ac.wing.r/2*(mC-1)),range(ac.wing.r*2*(mC-1),n*(mC-1))));
    for j in range(ac.wing.r/2-1):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    j = ac.wing.r/2-1;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NR_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie gauche du volet sur lui-même, l'aile gauche et les autres surfaces
    for j in range(ac.wing.r,ac.wing.r*3/2-1):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    j = ac.wing.r*3/2-1;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NR_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie droite de l'aile sur elle-même et les autres surfaces
    index = np.concatenate((range(ac.wing.r/2*(mC-1),ac.wing.r*(mC-1)),range(3*ac.wing.r/2*(mC-1),n*(mC-1))));
    j = ac.wing.r/2;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    for j in range(ac.wing.r/2+1,ac.wing.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie droite du volet sur lui-même, l'aile droite et les autres surfaces
    j = 3*ac.wing.r/2;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    for j in range(3*ac.wing.r/2+1,2*ac.wing.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie gauche de la HTail uniquement sur aile et elle-même
    index = np.array(range((2*ac.wing.r+ac.htail.r/2)*(mC-1)));
    for j in range(2*ac.wing.r,2*ac.wing.r+ac.htail.r/2-1):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    j = 2*ac.wing.r+ac.htail.r/2-1;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NR_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la partie droite de la HTail uniquement sur aile et elle-même
    index = np.concatenate((range(2*ac.wing.r*(mC-1)),range((2*ac.wing.r+ac.htail.r/2)*(mC-1),(2*ac.wing.r+ac.htail.r)*(mC-1))));
    j = 2*ac.wing.r+ac.htail.r/2;
    pathX = getVal(X,2*(mW+mC)+1,j);
    pathY = getVal(Y,2*(mW+mC)+1,j);
    pathZ = getVal(Z,2*(mW+mC)+1,j);
    test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    for j in range(2*ac.wing.r+ac.htail.r/2+1,2*ac.wing.r+ac.htail.r):
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    # Influence de la VTail sur tous les points
    index = np.array(range(n*(mC-1)));
    if ac.vtail.bool:
        j = 2*ac.wing.r+ac.htail.r;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        test.NL_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
        for j in range(2*ac.wing.r+ac.htail.r+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            test.Full_Influence(COLOCX,COLOCY,COLOCZ,normal,index,pathX,pathY,pathZ,A,Vx,Vy,Vz,mC,mW,j);
    if ac.prop.bool:
        test.On_PropInfluenceF(ac,X,Y,Z,n,m,mC,mW,A,Vx,Vy,Vz);
    return A,Vx,Vy,Vz;


def vortxl(x,y,z,normal,pathX,pathY,pathZ,mC,mW):
    """ Computing of the unit influence of the vortex on the colocation point
    
    
       Initially Copyright (C) 2004  Mihai Pruna, Alberto Davila
    
       Modified by Quentin Borlon (5 mai 2017)
    
       Same as proposed by Mondher Yahyaoui
       ( International Journal of Mechanical, Aerospace, Industrial, 
       Mechatronic and Manufacturing Engineering Vol:8, No:10, 2014 ).
    
       Exception : the influence of the vortex that goes to infinity."""
    
    nbRing = len(pathX) -1;                                                     # number of wake elements on the outer ring
    nbLine = nbRing + mC - 1;
    r1r2x = np.zeros(nbLine,dtype = float);
    r1r2y = np.zeros(nbLine,dtype = float);
    r1r2z = np.zeros(nbLine,dtype = float);
    square = np.zeros(nbLine,dtype = float);
    r1 = np.zeros(nbLine,dtype = float);
    r2 = np.zeros(nbLine,dtype = float);
    ror1 = np.zeros(nbLine,dtype = float);
    ror2 = np.zeros(nbLine,dtype = float);
    coeff = np.zeros(nbLine,dtype = float);
    
    # Contribution of the outer ring
    x1 = pathX[:-1];
    y1 = pathY[:-1];
    z1 = pathZ[:-1];
    x2 = pathX[1:];
    y2 = pathY[1:];
    z2 = pathZ[1:];
    rcut = 1e-15;
    r1r2x[:nbRing] = (y-y1)*(z-z2)-(z-z1)*(y-y2);
    r1r2y[:nbRing] = -((x-x1)*(z-z2)-(z-z1)*(x-x2));
    r1r2z[:nbRing] = (x-x1)*(y-y2)-(y-y1)*(x-x2);
    square[:nbRing] = r1r2x[:nbRing]*r1r2x[:nbRing]+r1r2y[:nbRing]*r1r2y[:nbRing]+r1r2z[:nbRing]*r1r2z[:nbRing];
    r1[:nbRing] = np.sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1));
    r2[:nbRing] = np.sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2));
    ror1[:nbRing] = (x2-x1)*(x-x1)+(y2-y1)*(y-y1)+(z2-z1)*(z-z1);
    ror2[:nbRing] = (x2-x1)*(x-x2)+(y2-y1)*(y-y2)+(z2-z1)*(z-z2);
    x1T = pathX[2:mC+1];
    y1T = pathY[2:mC+1];
    z1T = pathZ[2:mC+1];
    x2T = pathX[-2:-mC-1:-1];
    y2T = pathY[-2:-mC-1:-1];
    z2T = pathZ[-2:-mC-1:-1];
    r1r2x[nbRing:] = (y-y1T)*(z-z2T)-(z-z1T)*(y-y2T);
    r1r2y[nbRing:] = -((x-x1T)*(z-z2T)-(z-z1T)*(x-x2T));
    r1r2z[nbRing:] = (x-x1T)*(y-y2T)-(y-y1T)*(x-x2T);
    square[nbRing:] = r1r2x[nbRing:]*r1r2x[nbRing:]+r1r2y[nbRing:]*r1r2y[nbRing:]+r1r2z[nbRing:]*r1r2z[nbRing:];
    r1[nbRing:] = np.sqrt((x-x1T)*(x-x1T) + (y-y1T)*(y-y1T) + (z-z1T)*(z-z1T));
    r2[nbRing:] = np.sqrt((x-x2T)*(x-x2T) + (y-y2T)*(y-y2T) + (z-z2T)*(z-z2T));
    ror1[nbRing:] = (x2T-x1T)*(x-x1T)+(y2T-y1T)*(y-y1T)+(z2T-z1T)*(z-z1T);
    ror2[nbRing:] = (x2T-x1T)*(x-x2T)+(y2T-y1T)*(y-y2T)+(z2T-z1T)*(z-z2T);
    indice = np.array([not ((r1[i]<rcut) or (r2[i]<rcut) or (square[i]<rcut) ) for i in range(nbLine)],dtype = bool);
    coeff[indice] = 0.25/(m.pi*square[indice])*(ror1[indice]/r1[indice]-ror2[indice]/r2[indice]);
    ax = r1r2x * coeff;
    ay = r1r2y * coeff;
    az = r1r2z * coeff;
    a = np.zeros(mC-1,dtype = float);
    a[0] = (ax[0] + ax[1] + ax[nbRing] + ax[nbRing-1]) * normal[0] + \
        (ay[0] + ay[1] + ay[nbRing-1] + ay[nbRing]) * normal[1] + \
        (az[0] + az[1] + az[nbRing-1] + az[nbRing]) * normal[2];
    a[1:-1] = (ax[2:mC-1] + ax[nbRing+1:nbRing+mC-2] + ax[nbRing-2:nbRing-mC+1:-1] - ax[nbRing : nbRing + mC - 3]) * normal[0] + \
            (ay[2:mC-1] + ay[nbRing+1:nbRing+mC-2] + ay[nbRing-2:nbRing-mC+1:-1] - ay[nbRing : nbRing + mC - 3]) * normal[1] + \
            (az[2:mC-1] + az[nbRing+1:nbRing+mC-2] + az[nbRing-2:nbRing-mC+1:-1] - az[nbRing : nbRing + mC - 3]) * normal[2];
    a[-1] = (np.dot(r1r2x[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]) - ax[nbRing + mC - 3]) * normal[0] + \
            (np.dot(r1r2y[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]) - ay[nbRing + mC - 3]) * normal[1] + \
            (np.dot(r1r2z[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]) - az[nbRing + mC - 3]) * normal[2];
#    vi = np.array([np.dot(r1r2x[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1]),np.dot(r1r2y[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1]),np.dot(r1r2z[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1])]);
    vix = np.zeros(mC-1,dtype = float);
    viy = np.zeros(mC-1,dtype = float);
    viz = np.zeros(mC-1,dtype = float);
    vix[:-1] = (ax[1:mC-1] + ax[nbRing-1:nbRing-mC+1:-1]);
    viy[:-1] = (ay[1:mC-1] + ay[nbRing-1:nbRing-mC+1:-1]);
    viz[:-1] = (az[1:mC-1] + az[nbRing-1:nbRing-mC+1:-1]);
    vix[-1] = np.dot(r1r2x[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]);
    viy[-1] = np.dot(r1r2y[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]);
    viz[-1] = np.dot(r1r2z[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]);
    return a,vix,viy,viz;

def vortxl_NL(x,y,z,normal,pathX,pathY,pathZ,mC,mW):
    """ Computing of the unit influence of the vortex on the colocation point
    
    
       Initially Copyright (C) 2004  Mihai Pruna, Alberto Davila
    
       Modified by Quentin Borlon (5 mai 2017)
    
       Same as proposed by Mondher Yahyaoui
       ( International Journal of Mechanical, Aerospace, Industrial, 
       Mechatronic and Manufacturing Engineering Vol:8, No:10, 2014 ).
    
       Exception : the influence of the vortex that goes to infinity."""
    
#    nbRing = len(pathX) - (mC + mW);                                                     # number of wake elements on the outer ring
#    nbLine = nbRing + mC - 1 ;
#    r1r2x = np.zeros(nbLine,dtype = float);
#    r1r2y = np.zeros(nbLine,dtype = float);
#    r1r2z = np.zeros(nbLine,dtype = float);
#    square = np.zeros(nbLine,dtype = float);
#    r1 = np.zeros(nbLine,dtype = float);
#    r2 = np.zeros(nbLine,dtype = float);
#    ror1 = np.zeros(nbLine,dtype = float);
#    ror2 = np.zeros(nbLine,dtype = float);
#    coeff = np.zeros(nbLine,dtype = float);
#    
#    # Contribution of the outer ring
#    x1 = pathX[:-(mC + mW)];
#    y1 = pathY[:-(mC + mW)];
#    z1 = pathZ[:-(mC + mW)];
#    x2 = pathX[1:-(mC + mW)+1];
#    y2 = pathY[1:-(mC + mW)+1];
#    z2 = pathZ[1:-(mC + mW)+1];
#    rcut = 1e-15;
#    r1r2x[:nbRing] = (y-y1)*(z-z2)-(z-z1)*(y-y2);
#    r1r2y[:nbRing] = -((x-x1)*(z-z2)-(z-z1)*(x-x2));
#    r1r2z[:nbRing] = (x-x1)*(y-y2)-(y-y1)*(x-x2);
#    square[:nbRing] = r1r2x[:nbRing]*r1r2x[:nbRing]+r1r2y[:nbRing]*r1r2y[:nbRing]+r1r2z[:nbRing]*r1r2z[:nbRing];
#    r1[:nbRing] = np.sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1));
#    r2[:nbRing] = np.sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2));
#    ror1[:nbRing] = (x2-x1)*(x-x1)+(y2-y1)*(y-y1)+(z2-z1)*(z-z1);
#    ror2[:nbRing] = (x2-x1)*(x-x2)+(y2-y1)*(y-y2)+(z2-z1)*(z-z2);
#    x1T = pathX[2:mC+1];
#    y1T = pathY[2:mC+1];
#    z1T = pathZ[2:mC+1];
#    x2T = pathX[-2:-mC-1:-1];
#    y2T = pathY[-2:-mC-1:-1];
#    z2T = pathZ[-2:-mC-1:-1];
#    r1r2x[nbRing:] = (y-y1T)*(z-z2T)-(z-z1T)*(y-y2T);
#    r1r2y[nbRing:] = -((x-x1T)*(z-z2T)-(z-z1T)*(x-x2T));
#    r1r2z[nbRing:] = (x-x1T)*(y-y2T)-(y-y1T)*(x-x2T);
#    square[nbRing:] = r1r2x[nbRing:]*r1r2x[nbRing:]+r1r2y[nbRing:]*r1r2y[nbRing:]+r1r2z[nbRing:]*r1r2z[nbRing:];
#    r1[nbRing:] = np.sqrt((x-x1T)*(x-x1T) + (y-y1T)*(y-y1T) + (z-z1T)*(z-z1T));
#    r2[nbRing:] = np.sqrt((x-x2T)*(x-x2T) + (y-y2T)*(y-y2T) + (z-z2T)*(z-z2T));
#    ror1[nbRing:] = (x2T-x1T)*(x-x1T)+(y2T-y1T)*(y-y1T)+(z2T-z1T)*(z-z1T);
#    ror2[nbRing:] = (x2T-x1T)*(x-x2T)+(y2T-y1T)*(y-y2T)+(z2T-z1T)*(z-z2T);
#    indice = np.array([not ((r1[i]<rcut) or (r2[i]<rcut) or (square[i]<rcut) ) for i in range(nbLine)],dtype = bool);
#    coeff[indice] = 0.25/(m.pi*square[indice])*(ror1[indice]/r1[indice]-ror2[indice]/r2[indice]);
#    coeff[mC+mW+1:nbRing] *= 0;
#    ax = r1r2x * coeff;
#    ay = r1r2y * coeff;
#    az = r1r2z * coeff;
#    a = np.zeros(mC-1,dtype = float);
#    a[0] = (ax[0] + ax[1] + ax[nbRing]) * normal[0] + \
#        (ay[0] + ay[1] + ay[nbRing]) * normal[1] + \
#        (az[0] + az[1] + az[nbRing]) * normal[2];
#    a[1:-1] = (ax[2:mC-1] + ax[nbRing+1:nbRing+mC-2] - ax[nbRing : nbRing + mC - 3]) * normal[0] + \
#            (ay[2:mC-1] + ay[nbRing+1:nbRing+mC-2] - ay[nbRing : nbRing + mC - 3]) * normal[1] + \
#            (az[2:mC-1] + az[nbRing+1:nbRing+mC-2] - az[nbRing : nbRing + mC - 3]) * normal[2];
#    a[-1] = (np.dot(r1r2x[mC-1:mC+mW+1],coeff[mC-1:mC+mW+1]) - ax[nbRing + mC - 3]) * normal[0] + \
#            (np.dot(r1r2y[mC-1:mC+mW+1],coeff[mC-1:mC+mW+1]) - ay[nbRing + mC - 3]) * normal[1] + \
#            (np.dot(r1r2z[mC-1:mC+mW+1],coeff[mC-1:mC+mW+1]) - az[nbRing + mC - 3]) * normal[2];
##    vi = np.array([np.dot(r1r2x[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1]),np.dot(r1r2y[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1]),np.dot(r1r2z[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1])]);
#    vix = np.zeros(mC-1,dtype = float);
#    viy = np.zeros(mC-1,dtype = float);
#    viz = np.zeros(mC-1,dtype = float);
#    vix[:-1] = ax[1:mC-1];
#    viy[:-1] = ay[1:mC-1];
#    viz[:-1] = az[1:mC-1];
#    vix[-1] = np.dot(r1r2x[mC-1:mC+mW+1],coeff[mC-1:mC+mW+1]);
#    viy[-1] = np.dot(r1r2y[mC-1:mC+mW+1],coeff[mC-1:mC+mW+1]);
#    viz[-1] = np.dot(r1r2z[mC-1:mC+mW+1],coeff[mC-1:mC+mW+1]);
    nbRing = len(pathX) -1;                                                     # number of wake elements on the outer ring
    nbLine = nbRing + mC - 1;
    r1r2x = np.zeros(nbLine,dtype = float);
    r1r2y = np.zeros(nbLine,dtype = float);
    r1r2z = np.zeros(nbLine,dtype = float);
    square = np.zeros(nbLine,dtype = float);
    r1 = np.zeros(nbLine,dtype = float);
    r2 = np.zeros(nbLine,dtype = float);
    ror1 = np.zeros(nbLine,dtype = float);
    ror2 = np.zeros(nbLine,dtype = float);
    coeff = np.zeros(nbLine,dtype = float);
    
    # Contribution of the outer ring
    x1 = pathX[:-1];
    y1 = pathY[:-1];
    z1 = pathZ[:-1];
    x2 = pathX[1:];
    y2 = pathY[1:];
    z2 = pathZ[1:];
    rcut = 1e-15;
    r1r2x[:nbRing] = (y-y1)*(z-z2)-(z-z1)*(y-y2);
    r1r2y[:nbRing] = -((x-x1)*(z-z2)-(z-z1)*(x-x2));
    r1r2z[:nbRing] = (x-x1)*(y-y2)-(y-y1)*(x-x2);
    square[:nbRing] = r1r2x[:nbRing]*r1r2x[:nbRing]+r1r2y[:nbRing]*r1r2y[:nbRing]+r1r2z[:nbRing]*r1r2z[:nbRing];
    r1[:nbRing] = np.sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1));
    r2[:nbRing] = np.sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2));
    ror1[:nbRing] = (x2-x1)*(x-x1)+(y2-y1)*(y-y1)+(z2-z1)*(z-z1);
    ror2[:nbRing] = (x2-x1)*(x-x2)+(y2-y1)*(y-y2)+(z2-z1)*(z-z2);
    x1T = pathX[2:mC+1];
    y1T = pathY[2:mC+1];
    z1T = pathZ[2:mC+1];
    x2T = pathX[-2:-mC-1:-1];
    y2T = pathY[-2:-mC-1:-1];
    z2T = pathZ[-2:-mC-1:-1];
    r1r2x[nbRing:] = (y-y1T)*(z-z2T)-(z-z1T)*(y-y2T);
    r1r2y[nbRing:] = -((x-x1T)*(z-z2T)-(z-z1T)*(x-x2T));
    r1r2z[nbRing:] = (x-x1T)*(y-y2T)-(y-y1T)*(x-x2T);
    square[nbRing:] = r1r2x[nbRing:]*r1r2x[nbRing:]+r1r2y[nbRing:]*r1r2y[nbRing:]+r1r2z[nbRing:]*r1r2z[nbRing:];
    r1[nbRing:] = np.sqrt((x-x1T)*(x-x1T) + (y-y1T)*(y-y1T) + (z-z1T)*(z-z1T));
    r2[nbRing:] = np.sqrt((x-x2T)*(x-x2T) + (y-y2T)*(y-y2T) + (z-z2T)*(z-z2T));
    ror1[nbRing:] = (x2T-x1T)*(x-x1T)+(y2T-y1T)*(y-y1T)+(z2T-z1T)*(z-z1T);
    ror2[nbRing:] = (x2T-x1T)*(x-x2T)+(y2T-y1T)*(y-y2T)+(z2T-z1T)*(z-z2T);
    indice = np.array([not ((r1[i]<rcut) or (r2[i]<rcut) or (square[i]<rcut) ) for i in range(nbLine)],dtype = bool);
    coeff[indice] = 0.25/(m.pi*square[indice])*(ror1[indice]/r1[indice]-ror2[indice]/r2[indice]);
    coeff[mC+mW+1:nbRing] *=0.3;
    ax = r1r2x * coeff;
    ay = r1r2y * coeff;
    az = r1r2z * coeff;
    a = np.zeros(mC-1,dtype = float);
    a[0] = (ax[0] + ax[1] + ax[nbRing] + ax[nbRing-1]) * normal[0] + \
        (ay[0] + ay[1] + ay[nbRing-1] + ay[nbRing]) * normal[1] + \
        (az[0] + az[1] + az[nbRing-1] + az[nbRing]) * normal[2];
    a[1:-1] = (ax[2:mC-1] + ax[nbRing+1:nbRing+mC-2] + ax[nbRing-2:nbRing-mC+1:-1] - ax[nbRing : nbRing + mC - 3]) * normal[0] + \
            (ay[2:mC-1] + ay[nbRing+1:nbRing+mC-2] + ay[nbRing-2:nbRing-mC+1:-1] - ay[nbRing : nbRing + mC - 3]) * normal[1] + \
            (az[2:mC-1] + az[nbRing+1:nbRing+mC-2] + az[nbRing-2:nbRing-mC+1:-1] - az[nbRing : nbRing + mC - 3]) * normal[2];
    a[-1] = (np.dot(r1r2x[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]) - ax[nbRing + mC - 3]) * normal[0] + \
            (np.dot(r1r2y[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]) - ay[nbRing + mC - 3]) * normal[1] + \
            (np.dot(r1r2z[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]) - az[nbRing + mC - 3]) * normal[2];
#    vi = np.array([np.dot(r1r2x[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1]),np.dot(r1r2y[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1]),np.dot(r1r2z[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1])]);
    vix = np.zeros(mC-1,dtype = float);
    viy = np.zeros(mC-1,dtype = float);
    viz = np.zeros(mC-1,dtype = float);
    vix[:-1] = (ax[1:mC-1] + ax[nbRing-1:nbRing-mC+1:-1]);
    viy[:-1] = (ay[1:mC-1] + ay[nbRing-1:nbRing-mC+1:-1]);
    viz[:-1] = (az[1:mC-1] + az[nbRing-1:nbRing-mC+1:-1]);
    vix[-1] = np.dot(r1r2x[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]);
    viy[-1] = np.dot(r1r2y[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]);
    viz[-1] = np.dot(r1r2z[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]);
    return a,vix,viy,viz;
def vortxl_NR(x,y,z,normal,pathX,pathY,pathZ,mC,mW):
    """ Computing of the unit influence of the vortex on the colocation point
    
    
       Initially Copyright (C) 2004  Mihai Pruna, Alberto Davila
    
       Modified by Quentin Borlon (5 mai 2017)
    
       Same as proposed by Mondher Yahyaoui
       ( International Journal of Mechanical, Aerospace, Industrial, 
       Mechatronic and Manufacturing Engineering Vol:8, No:10, 2014 ).
    
       Exception : the influence of the vortex that goes to infinity."""
    
#    nbRing = len(pathX) - (mC + mW);                                                     # number of wake elements on the outer ring
#    nbLine = nbRing + mC - 1 ;
#    r1r2x = np.zeros(nbLine,dtype = float);
#    r1r2y = np.zeros(nbLine,dtype = float);
#    r1r2z = np.zeros(nbLine,dtype = float);
#    square = np.zeros(nbLine,dtype = float);
#    r1 = np.zeros(nbLine,dtype = float);
#    r2 = np.zeros(nbLine,dtype = float);
#    ror1 = np.zeros(nbLine,dtype = float);
#    ror2 = np.zeros(nbLine,dtype = float);
#    coeff = np.zeros(nbLine,dtype = float);
#    
#    # Contribution of the outer ring
#    x1 = pathX[range((mC + mW),len(pathX))];
#    y1 = pathY[range((mC + mW),len(pathX))];
#    z1 = pathZ[range((mC + mW),len(pathX))];
#    x2 = pathX[np.concatenate([range((mC + mW)+1,len(pathX)),[1]])];
#    y2 = pathY[np.concatenate([range((mC + mW)+1,len(pathX)),[1]])];
#    z2 = pathZ[np.concatenate([range((mC + mW)+1,len(pathX)),[1]])];
#    rcut = 1e-15;
#    r1r2x[:nbRing] = (y-y1)*(z-z2)-(z-z1)*(y-y2);
#    r1r2y[:nbRing] = -((x-x1)*(z-z2)-(z-z1)*(x-x2));
#    r1r2z[:nbRing] = (x-x1)*(y-y2)-(y-y1)*(x-x2);
#    square[:nbRing] = r1r2x[:nbRing]*r1r2x[:nbRing]+r1r2y[:nbRing]*r1r2y[:nbRing]+r1r2z[:nbRing]*r1r2z[:nbRing];
#    r1[:nbRing] = np.sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1));
#    r2[:nbRing] = np.sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2));
#    ror1[:nbRing] = (x2-x1)*(x-x1)+(y2-y1)*(y-y1)+(z2-z1)*(z-z1);
#    ror2[:nbRing] = (x2-x1)*(x-x2)+(y2-y1)*(y-y2)+(z2-z1)*(z-z2);
#    x1T = pathX[2:mC+1];
#    y1T = pathY[2:mC+1];
#    z1T = pathZ[2:mC+1];
#    x2T = pathX[-2:-mC-1:-1];
#    y2T = pathY[-2:-mC-1:-1];
#    z2T = pathZ[-2:-mC-1:-1];
#    r1r2x[nbRing:] = (y-y1T)*(z-z2T)-(z-z1T)*(y-y2T);
#    r1r2y[nbRing:] = -((x-x1T)*(z-z2T)-(z-z1T)*(x-x2T));
#    r1r2z[nbRing:] = (x-x1T)*(y-y2T)-(y-y1T)*(x-x2T);
#    square[nbRing:] = r1r2x[nbRing:]*r1r2x[nbRing:]+r1r2y[nbRing:]*r1r2y[nbRing:]+r1r2z[nbRing:]*r1r2z[nbRing:];
#    r1[nbRing:] = np.sqrt((x-x1T)*(x-x1T) + (y-y1T)*(y-y1T) + (z-z1T)*(z-z1T));
#    r2[nbRing:] = np.sqrt((x-x2T)*(x-x2T) + (y-y2T)*(y-y2T) + (z-z2T)*(z-z2T));
#    ror1[nbRing:] = (x2T-x1T)*(x-x1T)+(y2T-y1T)*(y-y1T)+(z2T-z1T)*(z-z1T);
#    ror2[nbRing:] = (x2T-x1T)*(x-x2T)+(y2T-y1T)*(y-y2T)+(z2T-z1T)*(z-z2T);
#    indice = np.array([not ((r1[i]<rcut) or (r2[i]<rcut) or (square[i]<rcut) ) for i in range(nbLine)],dtype = bool);
#    coeff[indice] = 0.25/(m.pi*square[indice])*(ror1[indice]/r1[indice]-ror2[indice]/r2[indice]);
#    ax = r1r2x * coeff;
#    ay = r1r2y * coeff;
#    az = r1r2z * coeff;
#    a = np.zeros(mC-1,dtype = float);
#    a[0] = (ax[nbRing-2] + ax[nbRing] + ax[nbRing-1]) * normal[0] + \
#        (ay[nbRing-2] + ay[nbRing] + ay[nbRing-1]) * normal[1] + \
#        (az[nbRing-2] + az[nbRing] + az[nbRing-1]) * normal[2];
#    a[1:-1] = (ax[nbRing-3:nbRing-mC:-1] - ax[nbRing : nbRing + mC - 3]  + ax[nbRing +1 : nbRing + mC - 2]) * normal[0] + \
#            (ay[nbRing-3:nbRing-mC:-1] - ay[nbRing : nbRing + mC - 3]  + ay[nbRing +1 : nbRing + mC - 2]) * normal[1] + \
#            (az[nbRing-3:nbRing-mC:-1] - az[nbRing : nbRing + mC - 3]  + az[nbRing +1 : nbRing + mC - 2]) * normal[2];
#    a[-1] = (np.dot(r1r2x[0:nbRing-mC+1],coeff[0:nbRing-mC+1]) - ax[nbRing + mC - 3]) * normal[0] + \
#            (np.dot(r1r2y[0:nbRing-mC+1],coeff[0:nbRing-mC+1]) - ay[nbRing + mC - 3]) * normal[1] + \
#            (np.dot(r1r2z[0:nbRing-mC+1],coeff[0:nbRing-mC+1]) - az[nbRing + mC - 3]) * normal[2];
##    vi = np.array([np.dot(r1r2x[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1]),np.dot(r1r2y[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1]),np.dot(r1r2z[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1])]);
#    vix = np.zeros(mC-1,dtype = float);
#    viy = np.zeros(mC-1,dtype = float);
#    viz = np.zeros(mC-1,dtype = float);
#    vix[:-1] = ax[nbRing-2:nbRing-mC:-1];
#    viy[:-1] = ay[nbRing-2:nbRing-mC:-1];
#    viz[:-1] = az[nbRing-2:nbRing-mC:-1];
#    vix[-1] = np.dot(r1r2x[0:nbRing-mC+1],coeff[0:nbRing-mC+1]);
#    viy[-1] = np.dot(r1r2y[0:nbRing-mC+1],coeff[0:nbRing-mC+1]);
#    viz[-1] = np.dot(r1r2z[0:nbRing-mC+1],coeff[0:nbRing-mC+1]);
    nbRing = len(pathX) -1;                                                     # number of wake elements on the outer ring
    nbLine = nbRing + mC - 1;
    r1r2x = np.zeros(nbLine,dtype = float);
    r1r2y = np.zeros(nbLine,dtype = float);
    r1r2z = np.zeros(nbLine,dtype = float);
    square = np.zeros(nbLine,dtype = float);
    r1 = np.zeros(nbLine,dtype = float);
    r2 = np.zeros(nbLine,dtype = float);
    ror1 = np.zeros(nbLine,dtype = float);
    ror2 = np.zeros(nbLine,dtype = float);
    coeff = np.zeros(nbLine,dtype = float);
    
    # Contribution of the outer ring
    x1 = pathX[:-1];
    y1 = pathY[:-1];
    z1 = pathZ[:-1];
    x2 = pathX[1:];
    y2 = pathY[1:];
    z2 = pathZ[1:];
    rcut = 1e-15;
    r1r2x[:nbRing] = (y-y1)*(z-z2)-(z-z1)*(y-y2);
    r1r2y[:nbRing] = -((x-x1)*(z-z2)-(z-z1)*(x-x2));
    r1r2z[:nbRing] = (x-x1)*(y-y2)-(y-y1)*(x-x2);
    square[:nbRing] = r1r2x[:nbRing]*r1r2x[:nbRing]+r1r2y[:nbRing]*r1r2y[:nbRing]+r1r2z[:nbRing]*r1r2z[:nbRing];
    r1[:nbRing] = np.sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1));
    r2[:nbRing] = np.sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2));
    ror1[:nbRing] = (x2-x1)*(x-x1)+(y2-y1)*(y-y1)+(z2-z1)*(z-z1);
    ror2[:nbRing] = (x2-x1)*(x-x2)+(y2-y1)*(y-y2)+(z2-z1)*(z-z2);
    x1T = pathX[2:mC+1];
    y1T = pathY[2:mC+1];
    z1T = pathZ[2:mC+1];
    x2T = pathX[-2:-mC-1:-1];
    y2T = pathY[-2:-mC-1:-1];
    z2T = pathZ[-2:-mC-1:-1];
    r1r2x[nbRing:] = (y-y1T)*(z-z2T)-(z-z1T)*(y-y2T);
    r1r2y[nbRing:] = -((x-x1T)*(z-z2T)-(z-z1T)*(x-x2T));
    r1r2z[nbRing:] = (x-x1T)*(y-y2T)-(y-y1T)*(x-x2T);
    square[nbRing:] = r1r2x[nbRing:]*r1r2x[nbRing:]+r1r2y[nbRing:]*r1r2y[nbRing:]+r1r2z[nbRing:]*r1r2z[nbRing:];
    r1[nbRing:] = np.sqrt((x-x1T)*(x-x1T) + (y-y1T)*(y-y1T) + (z-z1T)*(z-z1T));
    r2[nbRing:] = np.sqrt((x-x2T)*(x-x2T) + (y-y2T)*(y-y2T) + (z-z2T)*(z-z2T));
    ror1[nbRing:] = (x2T-x1T)*(x-x1T)+(y2T-y1T)*(y-y1T)+(z2T-z1T)*(z-z1T);
    ror2[nbRing:] = (x2T-x1T)*(x-x2T)+(y2T-y1T)*(y-y2T)+(z2T-z1T)*(z-z2T);
    indice = np.array([not ((r1[i]<rcut) or (r2[i]<rcut) or (square[i]<rcut) ) for i in range(nbLine)],dtype = bool);
    coeff[indice] = 0.25/(m.pi*square[indice])*(ror1[indice]/r1[indice]-ror2[indice]/r2[indice]);
    coeff[1:mW+mC] *= 0.3;
    ax = r1r2x * coeff;
    ay = r1r2y * coeff;
    az = r1r2z * coeff;
    a = np.zeros(mC-1,dtype = float);
    a[0] = (ax[0] + ax[1] + ax[nbRing] + ax[nbRing-1]) * normal[0] + \
        (ay[0] + ay[1] + ay[nbRing-1] + ay[nbRing]) * normal[1] + \
        (az[0] + az[1] + az[nbRing-1] + az[nbRing]) * normal[2];
    a[1:-1] = (ax[2:mC-1] + ax[nbRing+1:nbRing+mC-2] + ax[nbRing-2:nbRing-mC+1:-1] - ax[nbRing : nbRing + mC - 3]) * normal[0] + \
            (ay[2:mC-1] + ay[nbRing+1:nbRing+mC-2] + ay[nbRing-2:nbRing-mC+1:-1] - ay[nbRing : nbRing + mC - 3]) * normal[1] + \
            (az[2:mC-1] + az[nbRing+1:nbRing+mC-2] + az[nbRing-2:nbRing-mC+1:-1] - az[nbRing : nbRing + mC - 3]) * normal[2];
    a[-1] = (np.dot(r1r2x[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]) - ax[nbRing + mC - 3]) * normal[0] + \
            (np.dot(r1r2y[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]) - ay[nbRing + mC - 3]) * normal[1] + \
            (np.dot(r1r2z[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]) - az[nbRing + mC - 3]) * normal[2];
#    vi = np.array([np.dot(r1r2x[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1]),np.dot(r1r2y[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1]),np.dot(r1r2z[mC:mC+2*mW+1],coeff[mC:mC+2*mW+1])]);
    vix = np.zeros(mC-1,dtype = float);
    viy = np.zeros(mC-1,dtype = float);
    viz = np.zeros(mC-1,dtype = float);
    vix[:-1] = (ax[1:mC-1] + ax[nbRing-1:nbRing-mC+1:-1]);
    viy[:-1] = (ay[1:mC-1] + ay[nbRing-1:nbRing-mC+1:-1]);
    viz[:-1] = (az[1:mC-1] + az[nbRing-1:nbRing-mC+1:-1]);
    vix[-1] = np.dot(r1r2x[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]);
    viy[-1] = np.dot(r1r2y[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]);
    viz[-1] = np.dot(r1r2z[mC-1:mC+2*mW+2],coeff[mC-1:mC+2*mW+2]);
    return a,vix,viy,viz;
def camber(naca,x):
    """ Compute the camber of the naca 4-,5- and 6-digits.
                #   Taken over and lightly adapted by Quentin Borlon
                # NACA Airfoil Generator
                # This function generates a set of points containing the coordinates of a
                # NACA airfoil from the 4 Digit Series, 5 Digit Series and 6 Series given
                # its number and, as additional features, the chordt, the number of points
                # to be calculated, spacing type (between linear and cosine spacing),
                # opened or closed trailing edge and the angle of attack of the airfoil.
                # It also plots the airfoil for further comprovation if it is the required
                # one by the user.
                # 
                # -------------------------------------------------------------------------
                #
                # MIT License
                # 
                # Copyright (c) 2016 Alejandro de Haro"""
    try:
        naca = float(naca);
        Cam = np.zeros(len(x));
        # 6-digits
        if m.floor(naca/(1e5)):
            a=(m.floor(naca/10000)%10)/10;  # Chordwise position of minimum pressure (2nd digit)
            c_li=(m.floor(naca/100)%10)/10;  # Design lift coefficient (4th digit)
            g=-1./(1-a)*(a**2*(0.5*m.log(a)-0.25)+0.25);  # G constant calculation
            h=1./(1-a)*(0.5*(1-a)**2*m.log(1-a)-0.25*(1-a)**2)+g; # H constant calculation
            #----------------------- CAMBER ---------------------------------------
            indice = np.array([not ((x[i] == 0. or x[i] == 1. or x[i] == a)) for i in range(len(x))],dtype = bool);
            for i in range(len(x)):
                if indice[i]:
                    Cam[i]=c_li/(2*m.pi*(a+1))*(1./(1-a)*(0.5*(a-x[i])**2*np.log(np.abs(a-x[i]))-0.5*(1-x[i])**2*np.log(1-x[i])+0.25*(1-x[i])**2-0.25*(a-x[i])**2)-x[i]*np.log(x[i])+g-h*x[i]); # Mean camber y coordinate
        # 5-digits
        elif m.floor(naca/(1e4)):
            p=(m.floor(naca/1000)%10)/20;  # Location of maximum camber (2nd digit)
            rn=(m.floor(naca/100)%10);    # Type of camber (3rd digit)
            if rn==0:
            #----------------------- STANDARD CAMBER ------------------------------
                #----------------------- CONSTANTS --------------------------------
                r=3.33333333333212*p**3+0.700000000000909*p**2+1.19666666666638*p-0.00399999999996247;    # R constant calculation by interpolation
                k1=1514933.33335235*p**4-1087744.00001147*p**3+286455.266669048*p**2-32968.4700001967*p+1420.18500000524;  # K1 constant calculation by interpolation
                #----------------------- CAMBER -----------------------------------
                for i in range(len(x)):
                    if x[i]<r:
                        Cam[i]=k1/6*(x[i]**3-3*r*x[i]**2+r**2*(3-r)*x[i]); # Mean camber y coordinate
                    else:
                        Cam[i]=k1*r**3/6*(1-x[i]);   # Mean camber y coordinate
            elif rn==1:
            #----------------------- REFLEXED CAMBER ------------------------------
                #----------------------- CONSTANTS --------------------------------
                r=10.6666666666861*p**3-2.00000000001601*p**2+1.73333333333684*p-0.0340000000002413;  # R constant calculation by interponation
                k1=-27973.3333333385*p**3+17972.8000000027*p**2-3888.40666666711*p+289.076000000022;  # K1 constant calculation by interpolation
                k2_k1=85.5279999999984*p**3-34.9828000000004*p**2+4.80324000000028*p-0.21526000000003;    # K1/K2 constant calculation by interpolation
                #----------------------- CAMBER -----------------------------------
                for i in range(len(x)):
                    if x[i]<r:
                        Cam[i]=k1/6*((x[i]-r)**3-k2_k1*(1-r)**3*x[i]-r**3*x[i]+r**3);   # Mean camber y coordinate
                    else:
                        Cam[i]=k1/6*(k2_k1*(x[i]-r)**3-k2_k1*(1-r)**3*x[i]-r**3*x[i]+r**3); # Mean camber y coordinate
        # 4-digits
        else:
            maxt=m.floor(naca/1e3)/100;    # Maximum camber (1st digit)
            p=(m.floor(naca/100)%10)/10; 
            #----------------------- CAMBER ---------------------------------------
            for i in range(len(x)):
                if x[i]<p:
                    Cam[i]=maxt*x[i]/p**2*(2*p-x[i]);    # Mean camber y coordinate
                else:
                    Cam[i]=maxt*(1-x[i])/(1-p)**2*(1+x[i]-2*p);    # Mean camber y coordinate
    except ValueError:
        Cam  = getCamFromDataFile(naca,x);
    return Cam;

def getCamFromDataFile(filePath,x):
    """ Function that loads the 2D polar data from the path"""
    section = u.justLoad(filePath,0);
    

    if (section[0,0] !=1. or section[0,1] !=0.):
            section[0,0]  = 1.;
            section[0,1] = 0.;
    if (section[-1,0] !=1 or section[-1,1] !=0):
        section[-1,0] = 1.;
        section[-1,1] = 0.;
    n = np.where(np.logical_and(section[:,0] ==0,section[:,1] == 0.));
    if not(np.any(n)):
        n = m.floor(np.size(section,axis=0)/2);
        section = np.concatenate([section[:n,:],[np.array([0.,0.])],section[n:,:]],0);
    else:
        n = n[0][0];
    inf = np.interp(x,np.flipud(section[:n+1,0]),np.flipud(section[0:n+1,1]));
    sup = np.interp(x,section[n:,0],section[n:,1]);
    Cam = (inf+sup)*0.5
    Cam[0] = 0;
    return Cam;

def getCamF(xF):
    section = u.justLoad('./PolarFiles/flaps.dat',0);
    Cam = np.interp(xF,section[:,0],section[:,1]);
    return Cam

def setTable(table,dim2,pan,val):
    i0 = pan*dim2;
    for i in range(len(val)):
        table[i0+i] = val[i];

def getVal(table,dim2,pan):
    i0 = pan*dim2;
    return table[i0:i0+dim2];

def vortxlV(x,y,z,x1,y1,z1,x2,y2,z2):
    """ Computing of the unit influence of the vortex on the colocation point
    
    
       Initially Copyright (C) 2004  Mihai Pruna, Alberto Davila
    
       Modified by Quentin Borlon (5 mai 2017)
    
       Same as proposed by Mondher Yahyaoui
       ( International Journal of Mechanical, Aerospace, Industrial, 
       Mechatronic and Manufacturing Engineering Vol:8, No:10, 2014 ).
    
       Exception : the influence of the vortex that goes to infinity."""

    nbPoint = len(x1);
    rcutSq=1e-8;
    rcut = 1e-8;
    r1r2x = (y-y1)*(z-z2)-(z-z1)*(y-y2);
    r1r2y = -((x-x1)*(z-z2)-(z-z1)*(x-x2));
    r1r2z = (x-x1)*(y-y2)-(y-y1)*(x-x2);
    square = r1r2x*r1r2x+r1r2y*r1r2y+r1r2z*r1r2z;
    r1 = np.sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1));
    r2 = np.sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2));
    indice = np.array([not ((r1[i]<rcut) or (r2[i]<rcut) or (square[i]<rcutSq) ) for i in range(nbPoint)],dtype = bool);
    ror1 = np.zeros(len(r1));
    ror2 = np.zeros(len(r1));
    ror1[indice] = (x2[indice]-x1[indice])*(x-x1[indice])+(y2[indice]-y1[indice])*(y-y1[indice])+(z2[indice]-z1[indice])*(z-z1[indice]);
    ror2[indice] = (x2[indice]-x1[indice])*(x-x2[indice])+(y2[indice]-y1[indice])*(y-y2[indice])+(z2[indice]-z1[indice])*(z-z2[indice]);
    coeff = np.zeros(len(r1));
    coeff[indice] = 0.25/(m.pi*square[indice])*(ror1[indice]/r1[indice]-ror2[indice]/r2[indice]);
    a = np.array([np.dot(r1r2x,coeff),np.dot(r1r2y,coeff),np.dot(r1r2z,coeff)]);
    coeff[0] = 0.;
    vi = np.array([np.dot(r1r2x,coeff),np.dot(r1r2y,coeff),np.dot(r1r2z,coeff)]);
    return a,vi;

def ICMatrixV(vtail,cla,flow):
    """ Prediction of aerodynamic characteristics of the vertical tail. 
        Assumed to be independant of the flow on the lifting surfaces to avoid
        too strong coupling with vortex of the horizontal tail. If interactions
        with HTP must be neglected to avoid infinite values, not necessary to 
        compute interaction because of the little influence of the wing on it.
        Allows to have smaller matrix and reduces a lot the cpu costs.
        
       Autor : Quentin borlon 
       Date : 28 october 2017
    
       Function that predicts the aerodynamic coefficients for a given vtail.
       Based on the vtail geometry and the sectional 2D aerodynamic datas.
    
          
     INPUT:
       clAlpha : vertical array with clAlphas(i) is the lift curve slope of the
           panel from wing.y(i) to wing.y(i+1);
       vtail : a structral object with as fields:
    
           b : span
           chord : vertical array with the chord at the root (1) any
               discontinuity of taper ratio (2:end-1) and at the tip (end)
           airfoil : a cell-array with each cell gives the airfoil naca number
               representation, cell 1 correspond to first panel after root.
           sweep : vertical array with wing.sweep(i) is the sweep angle of 
               the panel from wing.y(i) to wing.y(i+1) (rad)
           deltasFlaps : vertical array with wing.deltasFlaps(i) is the 
               flaps defection  of the panel from wing.y(i) to wing.y(i+1)
               (deg)
           r : number of spanwise panel along the vtail;
           cFlaps_cLoc : vertical array with wing.cFlaps_cLocs(i) is the 
               local flaps to chord ratio
           z : the  spanwise location of the limits of the panels
           discY : vertical array of the complete set of the spanwise location
           airfoilIndex : vertical array with wing.airfoilIndex(i) is the index of 
               the airfoil (wing.airfoil) to use for the section at wing.y(i)
           chordDistrib : vertical array with wing.chordDistrib(i) is the chord length of 
               the section at wing.y(i)
     
     OUTPUT:
       A : the influence coefficient matrix [n x n] such that A*{GAMMA/2} + {Q}*{normal} = 0
       normal : a [3 x (wing.getR()/2+1)] matrix that provides the normal downward
           of the panel."""
    # Recover the numerical parameters
    n = vtail.getR();                                  # spanwise discretisation number of panel
    mC = vtail.mC;                                                                   # chordwise discretisation number of checkpoint for the wake
    mW = flow.mW;
    beta = -flow.beta * m.pi/180;
    aoa = flow.at;
    # Recover the vtail parameters
    
    c = vtail.getChordDist();
    cf = vtail.getCF();
    
    x = vtail.getX();
    z = vtail.getZ();
    # Rudder, Assumed to be as plain flaps
    cf = vtail.getCF();
        
    if cf != 0: 
        xT = np.unique(np.concatenate([(1.-cf)*0.5*(np.cos(np.linspace(m.pi,0.,mC))+1.),[0.25]]));
        mC = len(xT);
    else:
        xT = np.unique(np.concatenate([0.5*(np.cos(np.linspace(m.pi,0.,mC))+1.),[0.25]]));
        mC = len(xT);
    yT = np.zeros([mC,len(vtail.getAF())],dtype = float);
    for ii in range(len(vtail.getAF())):
        yT[:,ii-1]= camber(vtail.getAF(ii),xT);
    
    

    X = np.zeros(n * (2 * (mC + mW)+1),dtype = float);
    Y = np.zeros(n * (2 * (mC + mW)+1),dtype = float);                                   # initialization
    Z = np.zeros(n * (2 * (mC + mW)+1),dtype = float);
    COLOCX=np.zeros((mC-1)*n);
    COLOCY=np.zeros((mC-1)*n);
    COLOCZ=np.zeros((mC-1)*n);
    normal = np.zeros([3,(mC-1)*n]);
    coef = 0.25+cla*0.25/m.pi;
    ds = np.zeros((mC-1)*n);                                                    # vector of area of any panel
    dS = np.zeros(n);                                                           # vector of area of a spanwise section
    
    xvl = np.zeros(mC + mW,dtype = float);
    yvl = np.zeros(mC + mW,dtype = float);
    zvl = np.zeros(mC + mW,dtype = float);
    xvt = np.zeros(mC + mW,dtype = float);
    yvt = np.zeros(mC + mW,dtype = float);
    zvt = np.zeros(mC + mW,dtype = float);
    dydx = np.zeros(mW-1,dtype = float);
    dzdx = np.zeros(mW-1,dtype = float);
    for i in range(n):
        camb = yT[:,vtail.getAFI(i)]
        
        il = i;
        cl = c[il];
        
        xl = (xT - 0.25) * cl + x[il];
        yl = camb * cl;
        zl = z[il] * np.ones(mC);
        if vtail.getDF(i) != 0.:
            delta = vtail.getDF(i);
            RotF = u.rotz(delta);
            center = np.array([xl[-2],yl[-2],zl[-2]]);
            point = np.array([xl[-1],yl[-1],zl[-1]])-center;
            point = np.dot(RotF,point) + center;
            xl[-1] = point[0];
            yl[-1] = point[1];
            zl[-1] = point[2];
            
        xvl[:mC-1] = 0.75 * xl[:-1] + 0.25 * xl[1:];
        yvl[:mC-1] = 0.75 * yl[:-1] + 0.25 * yl[1:];
        zvl[:mC-1] = 0.75 * zl[:-1] + 0.25 * zl[1:];
        xvl[mC-1] = xvl[mC-2] + (xl[-1]-xl[-2]);
        yvl[mC-1] = yvl[mC-2] + (yl[-1]-yl[-2]);
        zvl[mC-1] = zvl[mC-2] + (zl[-1]-zl[-2]);
        # End of chord vortex = begining of wake vortex
        xvl[mC:-1] = xvl[mC-1] + 2.5 * cl * (1.+np.array(range(mW-1),dtype = float))/mW;
        xvl[-1] = 50. * vtail.b;
        dydxl = (yl[mC-1]-yl[mC-2])/(xl[mC-1]-xl[mC-2]);
        dydx = dydxl * np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC]))\
            + m.tan(beta) * (1.-np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])));
        dzdx = m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])));
        for ii in range(mW-1):
            zvl[mC+ii] = zvl[mC+(ii-1)] + dzdx[ii] * (xvl[mC+ii] - xvl[mC+(ii-1)]);
            yvl[mC+ii] = yvl[mC+(ii-1)] + dydx[ii] * (xvl[mC+ii] - xvl[mC+(ii-1)]);
        zvl[-1] = zvl[-2] + m.tan(aoa) * (xvl[-1] - xvl[-2]);
        yvl[-1] = yvl[-2] + m.tan(beta) * (xvl[-1] - xvl[-2]);
        
        it = i+1;
        ct = c[it];
        
        xt = (xT - 0.25) * ct + x[it];
        yt = camb * ct;
        zt = z[it] * np.ones(mC);
        if vtail.getDF(i) != 0.:
            delta = vtail.getDF(i);
            RotF = u.rotz(-delta);
            center = np.array([xt[-2],yt[-2],zt[-2]]);
            point = np.array([xt[-1],yt[-1],zt[-1]])-center;
            point = np.dot(RotF,point) + center;
            xt[-1] = point[0];
            yt[-1] = point[1];
            zt[-1] = point[2];
            
        xvt[:mC-1] = 0.75 * xt[:-1] + 0.25 * xt[1:];
        yvt[:mC-1] = 0.75 * yt[:-1] + 0.25 * yt[1:];
        zvt[:mC-1] = 0.75 * zt[:-1] + 0.25 * zt[1:];
        xvt[mC-1] = xvt[mC-2] + (xt[-1]-xt[-2]);
        yvt[mC-1] = yvt[mC-2] + (yt[-1]-yt[-2]);
        zvt[mC-1] = zvt[mC-2] + (zt[-1]-zt[-2]);
        # End of chord vortex = begining of wake vortex
        xvt[mC:-1] = xvt[mC-1] + 2.5 * ct * (1.+np.array(range(mW-1),dtype = float))/mW;
        xvt[-1] = 50. * vtail.b;
        dydxt = (yt[mC-1]-yt[mC-2])/(xt[mC-1]-xt[mC-2]);
        dydx = dydxt * np.exp(-3.*(np.array(xvt[mC:-1] - xvt[mC]))/(xvt[-2] - xvt[mC]))\
            + m.tan(beta) * (1.-np.exp(-3.*(np.array(xvt[mC:-1] - xvl[mC]))/(xvt[-2] - xvt[mC])));
        dzdx = m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvt[mC:-1] - xvl[mC]))/(xvt[-2] - xvt[mC])));
        for ii in range(mW-1):
            zvt[mC+ii] = zvt[mC+(ii-1)] + dzdx[ii] * (xvt[mC+ii] - xvt[mC+(ii-1)]);
            yvt[mC+ii] = yvt[mC+(ii-1)] + dydx[ii] * (xvt[mC+ii] - xvt[mC+(ii-1)]);
        zvt[-1] = zvt[-2] + m.tan(aoa) * (xvt[-1] - xvt[-2]);
        yvt[-1] = yvt[-2] + m.tan(beta) * (xvt[-1] - xvt[-2]);
        
        setTable(X,2*(mC+mW)+1,i,np.concatenate([[xvl[0]],xvt,xvl[::-1]]));
        setTable(Y,2*(mC+mW)+1,i,np.concatenate([[yvl[0]],yvt,yvl[::-1]]));
        setTable(Z,2*(mC+mW)+1,i,np.concatenate([[zvl[0]],zvt,zvl[::-1]]));
        for j in range(mC-1):
            
            val = [xvl[j],xvt[j], 0.5* (xl[j+1] + xt[j+1]),0.5* (xl[j] + xt[j])];
            COLOCX[i * (mC-1) + j] = val[3] * (1.-coef[i]) + val[2] * coef[i];
            cpx1 = val[1] - val[0];
            cpx2 = val[3] - val[2];
            
            val = [yvl[j],yvt[j], 0.5* (yl[j+1] + yt[j+1]),0.5* (yl[j] + yt[j])];
            COLOCY[i * (mC-1) + j] = val[3] * (1.-coef[i]) + val[2] * coef[i];
            cpy1 = val[1] - val[0];
            cpy2 = val[3] - val[2];
            
            val = [zvl[j],zvt[j], 0.5* (zl[j+1] + zt[j+1]),0.5* (zl[j] + zt[j])];
            COLOCZ[i * (mC-1) + j] = val[3] * (1.-coef[i]) + val[2] * coef[i];
            cpz1 = val[1] - val[0];
            cpz2 = val[3] - val[2];
            
            cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
            cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
            ds[i * (mC-1) + j] = cpmag;
            normal[:, i * (mC-1) + j] = cp/cpmag;
        dS[i] = sum(ds[i * (mC-1):(i+1) * (mC-1)]);
    select = np.zeros([vtail.r,n * (mC-1)]); # rechercher intensité du dernier vortex uniquement
    select2 = np.zeros([n * (mC-1),vtail.r]); # pour chaque paneau sur même section y, même velocity triangle
    select3 = np.zeros([vtail.r,n * (mC-1)]); # 
    for i in range(vtail.r):
        select[i,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i] = 1.;
        select3[i,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i];
    ##
    Ao,Vxo,Vyo,Vzo = ICM_V(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,n,mC,mW);
    A = np.linalg.inv(Ao);
    Vx = np.dot(select3,Vxo);
    Vy = np.dot(select3,Vyo);
    Vz = np.dot(select3,Vzo);
    return A,normal,Vx,Vy,Vz,select,select2;

def ICM_V(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,n,mC,mW):
    A = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vx = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vy = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    Vz = np.zeros([n*(mC-1),n*(mC-1)],dtype = float);
    for b in range(n * (mC - 1)):
        j = 0;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    return A,Vx,Vy,Vz;