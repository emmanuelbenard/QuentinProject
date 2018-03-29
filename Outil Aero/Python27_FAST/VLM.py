# -*- coding: utf-8 -*-
import math as m
import numpy as np
import utilitaire as u
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
def ICMatrix(ac,cla):
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

    wing = ac.wing;
    htail = ac.htail;
    # Recover the numerical parameters
    n = wing.getR()+htail.getR();                                  # spanwise discretisation number of panel
    mC =  15;                                                                   # chordwise discretisation number of checkpoint for the wake
    mW = 6;
    
    # Recover the wing parameters
    
    c = wing.getChordDist();
    cf = wing.getCF();
    tw = wing.getTwist();
    dih = wing.getDih();
    sw = wing.getSweepC4();
    
    x = wing.getX();
    y = wing.getY();
    z = wing.getZ();
    xp = wing.getXP();
    yp = wing.getYP();
    zp = wing.getZP();
    cP = wing.getChord();
    
    twSec = wing.twSec;
    if cf != 0: 
        xW = np.unique(np.concatenate([np.linspace(1,1.-cf,4),np.linspace(1.-cf,0,mC-3)]));
    else:
        xW = np.linspace(0,1,mC);
        
    zW = np.zeros([mC,len(wing.getAF())],dtype = float);
    for ii in range(len(wing.getAF())):
        zW[:,ii-1]= camber(wing.getAF(ii),xW);
    
    if htail.bool:
        c = np.concatenate([c,htail.getChordDist()]);
        cfT = htail.getCF();
        tw = np.concatenate([tw,htail.getTwist()]);
        twSec = np.concatenate([wing.twSec,htail.twSec]);
        dih = np.concatenate([dih,htail.getDih()]);
        sw = np.concatenate([sw,htail.getSweepC4()]);
    
        x = np.concatenate([x,htail.getX()]);
        y = np.concatenate([y,htail.getY()]);
        z = np.concatenate([z,htail.getZ()]);
        xp = np.concatenate([xp,htail.getXP()]);
        yp = np.concatenate([yp,htail.getYP()]);
        zp = np.concatenate([zp,htail.getZP()]);
        cP = np.concatenate([cP,htail.getChord()]);
        if cfT != 0: 
            xT = np.unique(np.concatenate([np.linspace(1,1.-cfT,2),np.linspace(1.-cfT,0,mC-1)]));
        else:
            xT = np.linspace(0,1,mC);
            
        zT = np.zeros([mC,len(htail.getAF())],dtype = float);
        for ii in range(len(htail.getAF())):
            zT[:,ii-1]= camber(htail.getAF(ii),xT);
            

    #generate grid corner coordinates
    # generate collocation points and normal : where tangency condition is
    # satisfied. Distance from bound vortex depends on the sectional lift
    # curve slope : (dist/localChord) = clAlphas/(4*pi)
    
    Ypanelgrid = np.zeros(4*n,dtype = float);
    Xpanelgrid = np.zeros(4*n,dtype = float);                                   # initialization
    Zpanelgrid = np.zeros(4*n,dtype = float);
    COLOCX=np.zeros(n);
    COLOCY=np.zeros(n);
    COLOCZ=np.zeros(n);
    normal = np.zeros([3,n]);
    coef = 0.25+cla*0.25/m.pi;
    for i in range(wing.getR()):
        val = [x[i],x[i+1],xp[i] - 0.25 * cP[i] * m.cos(tw[i]), xp[i] + 0.75 * cP[i] * m.cos(tw[i])];
        COLOCX[i] = val[2] * (1.-coef[i]) + val[3] * coef[i];
        cpx1 = val[1] - val[0];
        cpx2 = val[3] - val[2];
        setTable(Xpanelgrid,4,i,val);
        val = [y[i],y[i+1],yp[i],yp[i]];
        cpy1 = val[1] - val[0];
        cpy2 = val[3] - val[2];
        COLOCY[i] = val[2] * (1.-coef[i]) + val[3] * coef[i];
        setTable(Ypanelgrid,4,i,val);
        val = [z[i],z[i+1],zp[i] + 0.25 * cP[i] * m.sin(tw[i]), zp[i] - 0.75 * cP[i] * m.sin(tw[i])];
        cpz1 = val[1] - val[0];
        cpz2 = val[3] - val[2];
        COLOCZ[i] = val[2] * (1.-coef[i]) + val[3] * coef[i];
        setTable(Zpanelgrid,4,i,val);
        cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
        cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
        normal[:,i] = cp/cpmag;
    for i in range(wing.getR(),wing.getR()+htail.getR()):
        i0 = i+1;
        val = [x[i0],x[i0+1],xp[i] - 0.25 * cP[i] * m.cos(tw[i]), xp[i] + 0.75 * cP[i] * m.cos(tw[i])];
        COLOCX[i] = val[2] * (1.-coef[i]) + val[3] * coef[i];
        cpx1 = val[1] - val[0];
        cpx2 = val[3] - val[2];
        setTable(Xpanelgrid,4,i,val);
        val = [y[i0],y[i0+1],yp[i],yp[i]];
        cpy1 = val[1] - val[0];
        cpy2 = val[3] - val[2];
        COLOCY[i] = val[2] * (1.-coef[i]) + val[3] * coef[i];
        setTable(Ypanelgrid,4,i,val);
        val = [z[i0],z[i0+1],zp[i] + 0.25 * cP[i] * m.sin(tw[i]), zp[i] - 0.75 * cP[i] * m.sin(tw[i])];
        cpz1 = val[1] - val[0];
        cpz2 = val[3] - val[2];
        COLOCZ[i] = val[2] * (1.-coef[i]) + val[3] * coef[i];
        setTable(Zpanelgrid,4,i,val);
        cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
        cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
        normal[:,i] = cp/cpmag;

    # generate vortex ring corner coordinates
    # The bound vortex is not at the c/4 of the whole chord but at the 1/4
    # of the first panel ! If changes also modifiates the location (coef)
    # of the colocation points such that the distance between both is
    # conserved.
    
    X = np.zeros(((mC+mW)*2+1)*n);
    Y = np.zeros(((mC+mW)*2+1)*n);
    Z = np.zeros(((mC+mW)*2+1)*n);
    
    wxL = np.zeros(mW,dtype = float);
    wzL = np.zeros(mW,dtype = float);
    wxR = np.zeros(mW,dtype = float);
    wzR = np.zeros(mW,dtype = float);
    RotF = u.roty(180.);
    for i in range(wing.getR()):
        il = i;
        ir = i+1;
        pxL = xW * c[il] + x[il];
        pyL = np.ones(mC) * y[il];
        pzL = zW[:,wing.getAFI(i)]*c[il]+ z[il];
        center = np.array([pxL[0],pyL[0],pzL[0]]);
        alpha = 180./m.pi*twSec[il];
        Rot = u.roty(alpha);
        for ii in range(1,mC):
            point = np.array([pxL[ii],0,pzL[ii]])-center;
            point = np.dot(Rot,point) + center;
            pxL[ii] = point[0];
            pzL[ii] = point[2];
        delta = wing.getDF(i);
        RotD = u.roty(delta)
        centerD = np.array([pxL[-2],pyL[-2],pzL[-2]]);
        point = np.array([pxL[-1],0,pzL[-1]])-centerD;
        point = np.dot(RotD,point) + centerD;
        pxL[-1] = point[0];
        pzL[-1] = point[2];
        centreF = np.array([pxL[-1],pyL[-1],pzL[-1]]);
        for ii in range(mW):
            point = np.array([pxL[mC-2-ii],0,pzL[mC-2-ii]])-centreF;
            point = np.dot(RotF,point) + centreF;
            wxL[ii] = point[0];
            wzL[ii] = point[2];
        wxL[-1] = 10*wing.getSpan();
        wzL[-1] = wzL[-2];
        wyL = np.ones(mW) * y[il];
        
        pxR = xW * c[ir] + x[ir];
        pyR = np.ones(mC) * y[ir];
        pzR = zW[:,wing.getAFI(i)]*c[ir]+ z[ir];
        alpha = 180./m.pi*twSec[ir];
        Rot = u.roty(alpha);
        center = np.array([pxR[0],pyR[0],pzR[0]]);
        for ii in range(1,mC):
            point = np.array([pxR[ii],0,pzR[ii]])-center;
            point = np.dot(Rot,point) + center;
            pxR[ii] = point[0];
            pzR[ii] = point[2];
        centerD = np.array([pxR[-2],pyR[-2],pzR[-2]]);
        point = np.array([pxR[-1],0,pzR[-1]])-centerD;
        point = np.dot(RotD,point) + centerD;
        pxR[-1] = point[0];
        pzR[-1] = point[2];
        centreF = np.array([pxR[-1],pyR[-1],pzR[-1]]);
        for ii in range(mW):
            point = np.array([pxR[mC-2-ii],0,pzR[mC-2-ii]])-centreF;
            point = np.dot(RotF,point) + centreF;
            wxR[ii] = point[0];
            wzR[ii] = point[2];
        wxR[-1] = 10.*wing.getSpan();
        wzR[-1] = wzR[-2];
        wyR = np.ones(mW) * y[ir];
        val = np.concatenate([[pxL[0]],pxR,wxR,wxL[::-1],pxL[::-1]]);
        setTable(X,2*(mW+mC)+1,i,val);
        val = np.concatenate([[pyL[0]],pyR,wyR,wyL[::-1],pyL[::-1]]);
        setTable(Y,2*(mW+mC)+1,i,val);
        val = np.concatenate([[pzL[0]],pzR,wzR,wzL[::-1],pzL[::-1]]);
        setTable(Z,2*(mW+mC)+1,i,val);
    
    for i in range(wing.getR(),wing.getR()+htail.getR()):
        il = i+1;
        ir = i+2;
        pxL = xT * c[il] + x[il];
        pyL = np.ones(mC) * y[il];
        pzL = zT[:,htail.getAFI(i-wing.getR())]*c[il]+ z[il];
        center = np.array([pxL[0],pyL[0],pzL[0]]);
        alpha = 180./m.pi*twSec[il];
        Rot = u.roty(alpha);
        for ii in range(1,mC):
            point = np.array([pxL[ii],0,pzL[ii]])-center;
            point = np.dot(Rot,point) + center;
            pxL[ii] = point[0];
            pzL[ii] = point[2];
        delta = htail.getDF(i-wing.getR());
        RotD = u.roty(delta);
        centerD = np.array([pxL[-2],pyL[-2],pzL[-2]]);
        point = np.array([pxL[-1],0,pzL[-1]])-centerD;
        point = np.dot(RotD,point) + centerD;
        pxL[-1] = point[0];
        pzL[-1] = point[2];
        centreF = np.array([pxL[-1],pyL[-1],pzL[-1]]);
        for ii in range(mW):
            point = np.array([pxL[mC-2-ii],0,pzL[mC-2-ii]])-centreF;
            point = np.dot(RotF,point) + centreF;
            wxL[ii] = point[0];
            wzL[ii] = point[2];
        wxL[-1] = 10*wing.getSpan();
        wzL[-1] = wzL[-2];
        wyL = np.ones(mW) * y[il];
        pxR = xT * c[ir] + x[ir];
        pyR = np.ones(mC) * y[ir];
        pzR = zT[:,htail.getAFI(i-wing.getR())]*c[ir]+ z[ir];
        center = np.array([pxR[0],pyR[0],pzR[0]]);
        alpha = 180./m.pi*twSec[ir];
        Rot = u.roty(alpha);
        for ii in range(1,mC):
            point = np.array([pxR[ii],0,pzR[ii]])-center;
            point = np.dot(Rot,point) + center;
            pxR[ii] = point[0];
            pzR[ii] = point[2];
        centerD = np.array([pxR[-2],pyR[-2],pzR[-2]]);
        point = np.array([pxR[-1],0,pzR[-1]])-centerD;
        point = np.dot(RotD,point) + centerD;
        pxR[-1] = point[0];
        pzR[-1] = point[2];
        wyR = np.ones(mW) * y[ir];
        centreF = np.array([pxR[-1],pyR[-1],pzR[-1]]);
        for ii in range(mW):
            point = np.array([pxR[mC-2-ii],0,pzR[mC-2-ii]])-centreF;
            point = np.dot(RotF,point) + centreF;
            wxR[ii] = point[0];
            wzR[ii] = point[2];
        wxR[-1] = 10*wing.getSpan();
        wzR[-1] = wzR[-2];
        val = np.concatenate([[pxL[0]],pxR,wxR,wxL[::-1],pxL[::-1]]);
        setTable(X,2*(mW+mC)+1,i,val);
        val = np.concatenate([[pyL[0]],pyR,wyR,wyL[::-1],pyL[::-1]]);
        setTable(Y,2*(mW+mC)+1,i,val);
        val = np.concatenate([[pzL[0]],pzR,wzR,wzL[::-1],pzL[::-1]]);
        setTable(Z,2*(mW+mC)+1,i,val);
    A,Vx,Vy,Vz = ICM(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac,n,mC,mW);
    return A,normal,Vx,Vy,Vz;
        
def ICM(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac,n,mC,mW):
    if ac.fus.bool:
        HWing = ac.fus.vDist > 0;
    if ac.htail.bool and ac.vtail.bool:
        HTail = ac.htail.z[ac.htail.getR()/2] > ((ac.vtail.z[-1]-ac.vtail.z[0]) * 0.66) + ac.vtail.z[0];
    if not(ac.fus.bool):
        if not(ac.vtail.bool) or not(ac.htail.bool) or HTail:
            A,Vx,Vy,Vz = OnlyWing(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
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

def OnlyWing(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    velinduced=np.zeros(3,dtype = float);
    m = n;
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n,n],dtype = float);
    Vx = np.zeros([m,n],dtype = float);
    Vy = np.zeros([m,n],dtype = float);
    Vz = np.zeros([m,n],dtype = float);
    for b in range(n):
        for j in range(n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j] = np.dot(a,normal[:,b]);
            Vx[b,j] = velinduced[0];
            Vy[b,j] = velinduced[1];
            Vz[b,j] = velinduced[2];
    if ac.prop.bool:
        for b in range(n,m):
            x = ac.prop.xp[b-n];
            y = ac.prop.yp[b-n];
            z = ac.prop.zp[b-n];
            for j in range(n):
                pathX = getVal(X,2*(mW+mC)+1,j);
                pathY = getVal(Y,2*(mW+mC)+1,j);
                pathZ = getVal(Z,2*(mW+mC)+1,j);
                a,velinduced = vortxl(x,y,z,pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
                Vx[b,j] = velinduced[0];
                Vy[b,j] = velinduced[1];
                Vz[b,j] = velinduced[2];
    return A,Vx,Vy,Vz;

def BothWingOneTailVtail(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    velinduced=np.zeros(3,dtype = float);
    A = np.zeros([n,n],dtype = float);
    m = n;
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n,n],dtype = float);
    Vx = np.zeros([m,n],dtype = float);
    Vy = np.zeros([m,n],dtype = float);
    Vz = np.zeros([m,n],dtype = float);
    for b in range(ac.wing.getR()):
        for j in range(n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
    for b in range(ac.wing.getR(),ac.htail.getR()/2+ac.wing.getR()):
        for j in range(ac.htail.getR()/2+ac.wing.getR()-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathX = pathX[::-1];
        pathX = np.concatenate([[pathX[-2]],pathX]);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathY = pathY[::-1];
        pathY = np.concatenate([[pathY[-2]],pathY]);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        pathZ = pathZ[::-1];
        pathZ = np.concatenate([[pathZ[-2]],pathZ]);
        a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-(mC+mW+1)],pathY[:-(mC+mW+1)],pathZ[:-(mC+mW+1)],pathX[1:-(mC+mW)],pathY[1:-(mC+mW)],pathZ[1:-(mC+mW)]);
        A[b,j]=np.dot(-a,normal[:,b]);
        Vx[b,j]=-velinduced[0];
        Vy[b,j]=-velinduced[1];
        Vz[b,j]=-velinduced[2];
    for b in range(ac.wing.getR()+ac.htail.getR()/2,n):
        for j in range(ac.wing.getR()):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
        j = ac.wing.getR()+ac.htail.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:mC+mW+1],pathY[:mC+mW+1],pathZ[:mC+mW+1],pathX[1:mC+mW+2],pathY[1:mC+mW+2],pathZ[1:mC+mW+2]);
        A[b,j]=np.dot(a,normal[:,b]);
        Vx[b,j]=velinduced[0];
        Vy[b,j]=velinduced[1];
        Vz[b,j]=velinduced[2];
        for j in range(ac.wing.getR()+ac.htail.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
    if ac.prop.bool:
        for b in range(n,m):
            x = ac.prop.xp[b-n];
            y = ac.prop.yp[b-n];
            z = ac.prop.zp[b-n];
            for j in range(n):
                pathX = getVal(X,2*(mW+mC)+1,j);
                pathY = getVal(Y,2*(mW+mC)+1,j);
                pathZ = getVal(Z,2*(mW+mC)+1,j);
                a,velinduced = vortxl(x,y,z,pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
                Vx[b,j] = velinduced[0];
                Vy[b,j] = velinduced[1];
                Vz[b,j] = velinduced[2];
    return A,Vx,Vy,Vz;

def OneWingBothTail(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    velinduced=np.zeros(3,dtype = float);
    A = np.zeros([n,n],dtype = float);
    m = n;
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n,n],dtype = float);
    Vx = np.zeros([m,n],dtype = float);
    Vy = np.zeros([m,n],dtype = float);
    Vz = np.zeros([m,n],dtype = float);
    for b in range(ac.wing.getR()/2):
        for j in range(ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathX = pathX[::-1];
        pathX = np.concatenate([[pathX[-2]],pathX]);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathY = pathY[::-1];
        pathY = np.concatenate([[pathY[-2]],pathY]);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        pathZ = pathZ[::-1];
        pathZ = np.concatenate([[pathZ[-2]],pathZ]);
        a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-(mC+mW+1)],pathY[:-(mC+mW+1)],pathZ[:-(mC+mW+1)],pathX[1:-(mC+mW)],pathY[1:-(mC+mW)],pathZ[1:-(mC+mW)]);
        A[b,j]=np.dot(-a,normal[:,b]);
        Vx[b,j]=-velinduced[0];
        Vy[b,j]=-velinduced[1];
        Vz[b,j]=-velinduced[2];
        for j in range(ac.wing.getR(),n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
    for b in range(ac.wing.getR()/2,ac.wing.getR()):
        j = ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:mC+mW+1],pathY[:mC+mW+1],pathZ[:mC+mW+1],pathX[1:mC+mW+2],pathY[1:mC+mW+2],pathZ[1:mC+mW+2]);
        A[b,j]=np.dot(a,normal[:,b]);
        Vx[b,j]=velinduced[0];
        Vy[b,j]=velinduced[1];
        Vz[b,j]=velinduced[2];
        for j in range(ac.wing.getR()+ac.htail.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
    for b in range(ac.wing.getR(),n):
        for j in range(n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
    if ac.prop.bool:
        for b in range(n,m):
            for j in range(n):
                x = ac.prop.xp[b-n];
                y = ac.prop.yp[b-n];
                z = ac.prop.zp[b-n];
                pathX = getVal(X,2*(mW+mC)+1,j);
                pathY = getVal(Y,2*(mW+mC)+1,j);
                pathZ = getVal(Z,2*(mW+mC)+1,j);
                a,velinduced = vortxl(x,y,z,pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
                Vx[b,j] = velinduced[0];
                Vy[b,j] = velinduced[1];
                Vz[b,j] = velinduced[2];
    return A,Vx,Vy,Vz;

def OneWingOneTailVtail(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac):
    velinduced=np.zeros(3,dtype = float);
    A = np.zeros([n,n],dtype = float);
    m = n;
    if ac.prop.bool:
        nbE = len(ac.prop.D);
        m += nbE;
    A = np.zeros([n,n],dtype = float);
    Vx = np.zeros([m,n],dtype = float);
    Vy = np.zeros([m,n],dtype = float);
    Vz = np.zeros([m,n],dtype = float);
    for b in range(ac.wing.getR()/2):
        for j in range(ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathX = pathX[::-1];
        pathX = np.concatenate([[pathX[-2]],pathX]);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathY = pathY[::-1];
        pathY = np.concatenate([[pathY[-2]],pathY]);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        pathZ = pathZ[::-1];
        pathZ = np.concatenate([[pathZ[-2]],pathZ]);
        a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-(mC+mW+1)],pathY[:-(mC+mW+1)],pathZ[:-(mC+mW+1)],pathX[1:-(mC+mW)],pathY[1:-(mC+mW)],pathZ[1:-(mC+mW)]);
        A[b,j]=np.dot(-a,normal[:,b]);
        Vx[b,j]=-velinduced[0];
        Vy[b,j]=-velinduced[1];
        Vz[b,j]=-velinduced[2];
        for j in range(ac.wing.getR(),n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
    for b in range(ac.wing.getR()/2,ac.wing.getR()):
        j = ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:mC+mW+1],pathY[:mC+mW+1],pathZ[:mC+mW+1],pathX[1:mC+mW+2],pathY[1:mC+mW+2],pathZ[1:mC+mW+2]);
        A[b,j]=np.dot(a,normal[:,b]);
        Vx[b,j]=velinduced[0];
        Vy[b,j]=velinduced[1];
        Vz[b,j]=velinduced[2];
        for j in range(ac.wing.getR()+ac.htail.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
    for b in range(ac.wing.getR(),ac.htail.getR()/2+ac.wing.getR()):
        for j in range(ac.htail.getR()/2+ac.wing.getR()-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathX = pathX[::-1];
        pathX = np.concatenate([[pathX[-2]],pathX]);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathY = pathY[::-1];
        pathY = np.concatenate([[pathY[-2]],pathY]);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        pathZ = pathZ[::-1];
        pathZ = np.concatenate([[pathZ[-2]],pathZ]);
        a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-(mC+mW+1)],pathY[:-(mC+mW+1)],pathZ[:-(mC+mW+1)],pathX[1:-(mC+mW)],pathY[1:-(mC+mW)],pathZ[1:-(mC+mW)]);
        A[b,j]=np.dot(-a,normal[:,b]);
        Vx[b,j]=-velinduced[0];
        Vy[b,j]=-velinduced[1];
        Vz[b,j]=-velinduced[2];
    for b in range(ac.wing.getR()+ac.htail.getR()/2,n):
        for j in range(ac.wing.getR()):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
        j = ac.wing.getR()+ac.htail.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:mC+mW+1],pathY[:mC+mW+1],pathZ[:mC+mW+1],pathX[1:mC+mW+2],pathY[1:mC+mW+2],pathZ[1:mC+mW+2]);
        A[b,j]=np.dot(a,normal[:,b]);
        Vx[b,j]=velinduced[0];
        Vy[b,j]=velinduced[1];
        Vz[b,j]=velinduced[2];
        for j in range(ac.wing.getR()+ac.htail.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
    if ac.prop.bool:
        for b in range(n,m):
            for j in range(n):
                x = ac.prop.xp[b-n];
                y = ac.prop.yp[b-n];
                z = ac.prop.zp[b-n];
                pathX = getVal(X,2*(mW+mC)+1,j);
                pathY = getVal(Y,2*(mW+mC)+1,j);
                pathZ = getVal(Z,2*(mW+mC)+1,j);
                a,velinduced = vortxl(x,y,z,pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
                Vx[b,j] = velinduced[0];
                Vy[b,j] = velinduced[1];
                Vz[b,j] = velinduced[2];
    return A,Vx,Vy,Vz;


def vortxl(x,y,z,x1,y1,z1,x2,y2,z2):
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
    except ValueError as e:
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

def setTable(table,dim2,pan,val):
    i0 = pan*dim2;
    for i in range(len(val)):
        table[i0+i] = val[i];

def getVal(table,dim2,pan):
    i0 = pan*dim2;
    return table[i0:i0+dim2];

def ICMatrixV(vtail,cla):
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
    mC =  11;                                                                   # chordwise discretisation number of checkpoint for the wake
    mW = 6;
    # Recover the vtail parameters
    
    c = vtail.getChordDist();
    cf = vtail.getCF();
    
    x = vtail.getX();
    z = vtail.getZ();
    xp = vtail.getXP();
    yp = vtail.getYP();
    zp = vtail.getZP();
    cP = vtail.getChord();
    if cf != 0: 
        xVT = np.unique(np.concatenate([np.linspace(1,1.-cf,2),np.linspace(1.-cf,0,mC-1)]));
    else:
        xVT = np.linspace(0,1,mC);
    
    

    #generate grid corner coordinates
    # generate collocation points and normal : where tangency condition is
    # satisfied. Distance from bound vortex depends on the sectional lift
    # curve slope : (dist/localChord) = clAlphas/(4*pi)

    Xpanelgrid = np.zeros(4*n,dtype = float);                                   # initialization
    Zpanelgrid = np.zeros(4*n,dtype = float);
    COLOCX=np.zeros(n);
    COLOCY=np.zeros(n);
    COLOCZ=np.zeros(n);
    normal = np.zeros([3,n]);
    coef = 0.25+cla*0.25/m.pi;
    
    for i in range(n):
        i0 = i;
        val = [xp[i] - 0.25 * cP[i], xp[i] + 0.75 * cP[i], x[i0],x[i0+1]];
        cpx1 = val[1] - val[0];
        cpx2 = val[3] - val[2];
        COLOCX[i] = val[0] * (1.-coef[i]) + val[1] * coef[i];
        setTable(Xpanelgrid,4,i,val);
        val = [zp[i],zp[i],z[i0], z[i0+1]];
        cpz1 = val[1] - val[0];
        cpz2 = val[3] - val[2];
        COLOCZ[i] = zp[i];
        setTable(Zpanelgrid,4,i,val);
        cp= np.cross(np.array([cpx1,0.,cpz1]),np.array([cpx2,0.,cpz2]));
        cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
        normal[:,i] = cp/cpmag;
        
    # generate vortex ring corner coordinates
    # The bound vortex is not at the c/4 of the whole chord but at the 1/4
    # of the first panel ! If changes also modifiates the location (coef)
    # of the colocation points such that the distance between both is
    # conserved.
    
    X = np.zeros(((mC+mW)*2+1)*n);
    Y = np.zeros(((mC+mW)*2+1)*n);
    Z = np.zeros(((mC+mW)*2+1)*n);
    
    wxL = np.zeros(mW,dtype = float);
    wyL = np.zeros(mW,dtype = float);
    wzL = np.zeros(mW,dtype = float);
    wxR = np.zeros(mW,dtype = float);
    wyR = np.zeros(mW,dtype = float);
    wzR = np.zeros(mW,dtype = float);
    RotFV = u.rotz(180.);
    for i in range(n):
        il = i;
        ir = i+1;
        pxL = xVT * c[il] + x[il];
        pyL = np.zeros(mC);
        pzL = np.ones(mC) * z[il];
        delta = vtail.getDF(i);
        RotD = u.rotz(delta);
        centerD = np.array([pxL[-2],pyL[-2],pzL[-2]]);
        point = np.array([pxL[-1],0.,pzL[-1]])-centerD;
        point = np.dot(RotD,point) + centerD;
        pxL[-1] = point[0];
        pyL[-1] = point[1];
        centreF = np.array([pxL[-1],pyL[-1],pzL[-1]]);
        wzL = np.ones(mW) * z[il];
        for ii in range(mW):
            point = np.array([pxL[mC-2-ii],0.,pzL[mC-2-ii]])-centreF;
            point = np.dot(RotFV,point) + centreF;
            wxL[ii] = point[0];
            wyL[ii] = point[1];
        wxL[-1] = 100*vtail.getSpan();
        wyL[-1] = wyL[-2];
        pxR = xVT * c[ir] + x[ir];
        pyR = np.zeros(mC,dtype = float);
        pzR = np.ones(mC) * z[ir];
        centerD = np.array([pxR[-2],pyR[-2],pzR[-2]]);
        point = np.array([pxR[-1],0.,pzR[-1]])-centerD;
        point = np.dot(RotD,point) + centerD;
        pxR[-1] = point[0];
        pyR[-1] = point[1];
        centreF = np.array([pxR[-1],pyR[-1],pzR[-1]]);
        for ii in range(mW):
            point = np.array([pxR[mC-2-ii],0.,pzR[mC-2-ii]])-centreF;
            point = np.dot(RotFV,point) + centreF;
            wxR[ii] = point[0];
            wyR[ii] = point[1];
        wxR[-1] = 100*vtail.getSpan();
        wyR[-1] = wyR[-2];
        wzR = np.ones(mW) * z[ir];
        val = np.concatenate([[pxL[0]],pxR,wxR,wxL[::-1],pxL[::-1]]);
        setTable(X,2*(mW+mC)+1,i,val);
        val = np.concatenate([[pyL[0]],pyR,wyR,wyL[::-1],pyL[::-1]]);
        setTable(Y,2*(mW+mC)+1,i,val);
        val = np.concatenate([[pzL[0]],pzR,wzR,wzL[::-1],pzL[::-1]]);
        setTable(Z,2*(mW+mC)+1,i,val);
    [A,Vx,Vy,Vz] = ICM_V(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal);
    return A,normal,Vx,Vy,Vz;

def ICM_V(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal):
    velinduced=np.zeros(3,dtype = float);
    A = np.zeros([n,n],dtype = float);
    Vx = np.zeros([n,n],dtype = float);
    Vy = np.zeros([n,n],dtype = float);
    Vz = np.zeros([n,n],dtype = float);
    for b in range(n):
        j = 0;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:mW+mC+1],pathY[:mW+mC+1],pathZ[:mW+mC+1],pathX[1:mW+mC+2],pathY[1:mW+mC+2],pathZ[1:mW+mC+2]);
        A[b,j]=np.dot(a,normal[:,b]);
        Vx[b,j]=velinduced[0];
        Vy[b,j]=velinduced[1];
        Vz[b,j]=velinduced[2];
        for j in range(1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,velinduced = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],pathX[:-1],pathY[:-1],pathZ[:-1],pathX[1:],pathY[1:],pathZ[1:]);
            A[b,j]=np.dot(a,normal[:,b]);
            Vx[b,j]=velinduced[0];
            Vy[b,j]=velinduced[1];
            Vz[b,j]=velinduced[2];
    return A,Vx,Vy,Vz;