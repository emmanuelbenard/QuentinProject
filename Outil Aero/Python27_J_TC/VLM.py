  # -*- coding: utf-8 -*-
import math as m
import numpy as np
import utilitaire as u
import Polar as p
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Flow
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
    aoa = m.pi * (flow.getAMax()+flow.getAMin())/180.;
    
    # Main lifting surfaces
    wing = ac.wing;
    htail = ac.htail;
    
    # Numerical parameters for  discretization
    
    mC = wing.mC;                                                               # chordwise discretisation number of control point for the chord
    mW = max([8,int(3.*flow.V0/wing.getMac())]);                                # discretisation of the wake, get correct direction of it behind parts
    n = wing.getR()+htail.getR();                                               # spanwise discretisation number of panel
    
    # Recover the wing parameters
    
    # Panels' coordinates and main parameters (at c/4)
    xp = wing.getXP();
    yp = wing.getYP();
    zp = wing.getZP();
    
    cP = wing.getChord();
    tw = wing.getTwist();
    dih = wing.getDih();
    sw = wing.getSweepC4();
    
    # Panel bordes' coordinate and main parameters (at c/4)    
    x = wing.getX();
    y = wing.getY();
    z = wing.getZ();
    c = wing.getChordDist();
    twSec = wing.twSec;
    
    xW = np.unique(np.concatenate([0.5*(np.cos(np.linspace(m.pi,0.,mC))+1.),[0.25]]));
    mC = len(xW);
    iC4W = np.where(xW == 0.25)[0][0];
    
    zW = np.zeros([mC,len(wing.getAF())],dtype = float);
    for ii in range(len(wing.getAF())):
        zW[:,ii]= camber(wing.getAF(ii),xW);
    
    if htail.bool:
        # Panel bordes' coordinate and main parameters (at c/4)
        x = np.concatenate([x,htail.getX()]);
        y = np.concatenate([y,htail.getY()]);
        z = np.concatenate([z,htail.getZ()]);
        c = np.concatenate([c,htail.getChordDist()]);
        twSec = np.concatenate([wing.twSec,htail.twSec]);
        
        # Panels' coordinates and main parameters (at c/4)
        xp = np.concatenate([xp,htail.getXP()]);
        yp = np.concatenate([yp,htail.getYP()]);
        zp = np.concatenate([zp,htail.getZP()]);
        
        cP = np.concatenate([cP,htail.getChord()]);
        tw = np.concatenate([tw,htail.getTwist()]);
        dih = np.concatenate([dih,htail.getDih()]);
        sw = np.concatenate([sw,htail.getSweepC4()]);
        
        # Elevator, Assumed to be as plain flaps
        cfT = htail.getCF();
        
        if cfT != 0: 
            xT = np.unique(np.concatenate([np.linspace(1.,1.-cfT,2),(1.-cfT)*0.5*(np.cos(np.linspace(m.pi,0.,mC-1))+1.)]));
            xT[abs((xT-0.25)) == np.min(abs(xT-0.25))] = 0.25;
        else:
            xT = 0.5*(np.cos(np.linspace(m.pi,0.,mC))+1.);
            xT[abs((xT-0.25)) == np.min(abs(xT-0.25))] = 0.25;
        iC4T = np.where(xT == 0.25)[0][0];
        zT = np.zeros([mC,len(htail.getAF())],dtype = float);
        for ii in range(len(htail.getAF())):
            zT[:,ii-1]= camber(htail.getAF(ii),xT);
    
    
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
    xvr = np.zeros(mC + mW,dtype = float);
    yvr = np.zeros(mC + mW,dtype = float);
    zvr = np.zeros(mC + mW,dtype = float);
    dzdx = np.zeros(mW-1,dtype = float);                       
    dydx = np.zeros(mW-1,dtype = float);                       
    for i in range(wing.getR()):
        camb = zW[:,wing.getAFI(i)]
        
        il = i;
        cl = c[il];
        twl = twSec[il];
        
        xl = (xW - 0.25) * cl + x[il];
        yl = y[il] * np.ones(mC);
        zl = camb * cl + z[il];
        center = np.array([xl[iC4W],yl[iC4W],zl[iC4W]]);
        alpha = 180./m.pi*twl;
        Rot = u.roty(alpha);
        for ii in range(mC):
            point = np.array([xl[ii],yl[ii],zl[ii]])-center;
            point = np.dot(Rot,point) + center;
            xl[ii] = point[0];
            yl[ii] = point[1];
            zl[ii] = point[2];
        
            
        xvl[:mC-1] = 0.75 * xl[:-1] + 0.25 * xl[1:];
        yvl[:mC-1] = 0.75 * yl[:-1] + 0.25 * yl[1:];
        zvl[:mC-1] = 0.75 * zl[:-1] + 0.25 * zl[1:];
        xvl[mC-1] = xvl[mC-2] + (xl[-1]-xl[-2]);
        yvl[mC-1] = yvl[mC-2] + (yl[-1]-yl[-2]);
        zvl[mC-1] = zvl[mC-2] + (zl[-1]-zl[-2]);
        # End of chord vortex = begining of wake vortex
        xvl[mC:-1] = xvl[mC-1] + 2.5 * cl * (1.+np.array(range(mW-1),dtype = float))/mW;
        xvl[-1] = 10. * wing.b;
        dzdxl = (zl[mC-1]-zl[mC-2])/(xl[mC-1]-xl[mC-2]);
        dydx = m.tan(beta) * (1.-np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])));
        dzdx = dzdxl * np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])) \
            + m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])));
        for ii in range(mW-1):
            zvl[mC+ii] = zvl[mC+(ii-1)] + dzdx[ii] * (xvl[mC+ii] - xvl[mC+(ii-1)]);
            yvl[mC+ii] = yvl[mC+(ii-1)] + dydx[ii] * (xvl[mC+ii] - xvl[mC+(ii-1)]);
        zvl[-1] = zvl[-2] + m.tan(aoa) * (xvl[-1] - xvl[-2]);
        yvl[-1] = yvl[-2] + m.tan(beta) * (xvl[-1] - xvl[-2]);
        
        ir = i+1;
        cr = c[ir];
        twr = twSec[ir];
        
        xr = (xW - 0.25) * cr + x[ir];
        yr = y[ir] * np.ones(mC);
        zr = camb * cr + z[ir];
        center = np.array([xr[iC4W],yr[iC4W],zr[iC4W]]);
        alpha = 180./m.pi*twr;
        Rot = u.roty(alpha);
        for ii in range(0,mC):
            point = np.array([xr[ii],yr[ii],zr[ii]])-center;
            point = np.dot(Rot,point) + center;
            xr[ii] = point[0];
            yr[ii] = point[1];
            zr[ii] = point[2];
        
        xvr[:mC-1] = 0.75 * xr[:-1] + 0.25 * xr[1:];
        yvr[:mC-1] = 0.75 * yr[:-1] + 0.25 * yr[1:];
        zvr[:mC-1] = 0.75 * zr[:-1] + 0.25 * zr[1:];
        xvr[mC-1] = xvr[mC-2] + (xr[-1]-xr[-2]);
        yvr[mC-1] = yvr[mC-2] + (yr[-1]-yr[-2]);
        zvr[mC-1] = zvr[mC-2] + (zr[-1]-zr[-2]);
        # End of chord vortex = begining of wake vortex
        xvr[mC:-1] = xvr[mC-1] + 2.5 * cr * (1.+np.array(range(mW-1),dtype = float))/mW;
        xvr[-1] = 10. * wing.b;
        dzdxr = (zr[mC-1]-zr[mC-2])/(xr[mC-1]-xr[mC-2]);
        dydx = m.tan(beta) * (1.-np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])));
        dzdx = dzdxr * np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])) \
            + m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])));
        for ii in range(mW-1):
            zvr[mC+ii] = zvr[mC+(ii-1)] + dzdx[ii] * (xvr[mC+ii] - xvr[mC+(ii-1)]);
            yvr[mC+ii] = yvr[mC+(ii-1)] + dydx[ii] * (xvr[mC+ii] - xvr[mC+(ii-1)]);
        zvr[-1] = zvr[-2] + m.tan(aoa) * (xvr[-1] - xvr[-2]);
        yvr[-1] = yvr[-2] + m.tan(beta) * (xvr[-1] - xvr[-2]);
        
        setTable(X,2*(mC+mW)+1,i,np.concatenate([[xvl[0]],xvr,xvl[::-1]]));
        setTable(Y,2*(mC+mW)+1,i,np.concatenate([[yvl[0]],yvr,yvl[::-1]]));
        setTable(Z,2*(mC+mW)+1,i,np.concatenate([[zvl[0]],zvr,zvl[::-1]]));
        
        for j in range(mC-1):
            
            val = [xvl[j],xvr[j],0.5* (xl[j] + xr[j]), 0.5* (xl[j+1] + xr[j+1])];
            COLOCX[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpx1 = val[1] - val[0];
            cpx2 = val[3] - val[2];
            
            val = [yvl[j],yvr[j],0.5* (yl[j] + yr[j]), 0.5* (yl[j+1] + yr[j+1])];
            COLOCY[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpy1 = val[1] - val[0];
            cpy2 = val[3] - val[2];
            
            val = [zvl[j],zvr[j],0.5* (zl[j] + zr[j]), 0.5* (zl[j+1] + zr[j+1])];
            COLOCZ[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpz1 = val[1] - val[0];
            cpz2 = val[3] - val[2];
            
            cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
            cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
            ds[i * (mC-1) + j] = cpmag;
            normal[:, i * (mC-1) + j] = cp/cpmag;
        dS[i] = sum(ds[i * (mC-1):(i+1) * (mC-1)]);
    
    for i in range(wing.getR(),wing.getR()+htail.getR()):
        iPT = i-wing.getR();
        camb = zT[:,htail.getAFI(iPT)]
        
        il = i+1;
        cl = c[il];
        twl = twSec[il];
        
        xl = (xT - 0.25) * cl + x[il];
        yl = y[il] * np.ones(mC);
        zl = camb * cl + z[il];
        center = np.array([xl[iC4T],yl[iC4T],zl[iC4T]]);
        alpha = 180./m.pi*twl;
        Rot = u.roty(alpha);
        for ii in range(mC):
            point = np.array([xl[ii],yl[ii],zl[ii]])-center;
            point = np.dot(Rot,point) + center;
            xl[ii] = point[0];
            yl[ii] = point[1];
            zl[ii] = point[2];
        if htail.getDF(iPT) != 0.:
            delta = htail.getDF(iPT);
            RotF = u.roty(delta);
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
        xvl[-1] = 10. * wing.b;
        dzdxl = (zl[mC-1]-zl[mC-2])/(xl[mC-1]-xl[mC-2]);
        dydx = m.tan(beta) * (1.-np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])));
        dzdx = dzdxl * np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])) \
            + m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])));
        for ii in range(mW-1):
            zvl[mC+ii] = zvl[mC+(ii-1)] + dzdx[ii] * (xvl[mC+ii] - xvl[mC+(ii-1)]);
            yvl[mC+ii] = yvl[mC+(ii-1)] + dydx[ii] * (xvl[mC+ii] - xvl[mC+(ii-1)]);
        zvl[-1] = zvl[-2] + m.tan(aoa) * (xvl[-1] - xvl[-2]);
        yvl[-1] = yvl[-2] + m.tan(beta) * (xvl[-1] - xvl[-2]);
        
        ir = i+2;
        cr = c[ir];
        twr = twSec[ir];
        
        xr = (xT - 0.25) * cr + x[ir];
        yr = y[ir] * np.ones(mC);
        zr = camb * cr + z[ir];
        center = np.array([xr[iC4T],yr[iC4T],zr[iC4T]]);
        alpha = 180./m.pi*twr;
        Rot = u.roty(alpha);
        for ii in range(0,mC):
            point = np.array([xr[ii],yr[ii],zr[ii]])-center;
            point = np.dot(Rot,point) + center;
            xr[ii] = point[0];
            yr[ii] = point[1];
            zr[ii] = point[2];
        if htail.getDF(iPT) != 0.:
            delta = htail.getDF(iPT);
            RotF = u.roty(delta);
            center = np.array([xr[-2],yr[-2],zr[-2]]);
            point = np.array([xr[-1],yr[-1],zr[-1]])-center;
            point = np.dot(RotF,point) + center;
            xr[-1] = point[0];
            yr[-1] = point[1];
            zr[-1] = point[2];
            
        xvr[:mC-1] = 0.75 * xr[:-1] + 0.25 * xr[1:];
        yvr[:mC-1] = 0.75 * yr[:-1] + 0.25 * yr[1:];
        zvr[:mC-1] = 0.75 * zr[:-1] + 0.25 * zr[1:];
        xvr[mC-1] = xvr[mC-2] + (xr[-1]-xr[-2]);
        yvr[mC-1] = yvr[mC-2] + (yr[-1]-yr[-2]);
        zvr[mC-1] = zvr[mC-2] + (zr[-1]-zr[-2]);
        # End of chord vortex = begining of wake vortex
        xvr[mC:-1] = xvr[mC-1] + 2.5 * cr * (1.+np.array(range(mW-1),dtype = float))/mW;
        xvr[-1] = 10. * wing.b;
        dzdxr = (zr[mC-1]-zr[mC-2])/(xr[mC-1]-xr[mC-2]);
        dydx = m.tan(beta) * (1.-np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])));
        dzdx = dzdxr * np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])) \
            + m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])));
        for ii in range(mW-1):
            zvr[mC+ii] = zvr[mC+(ii-1)] + dzdx[ii] * (xvr[mC+ii] - xvr[mC+(ii-1)]);
            yvr[mC+ii] = yvr[mC+(ii-1)] + dydx[ii] * (xvr[mC+ii] - xvr[mC+(ii-1)]);
        zvr[-1] = zvr[-2] + m.tan(aoa) * (xvr[-1] - xvr[-2]);
        yvr[-1] = yvr[-2] + m.tan(beta) * (xvr[-1] - xvr[-2]);
        setTable(X,2*(mC+mW)+1,i,np.concatenate([[xvl[0]],xvr,xvl[::-1]]));
        setTable(Y,2*(mC+mW)+1,i,np.concatenate([[yvl[0]],yvr,yvl[::-1]]));
        setTable(Z,2*(mC+mW)+1,i,np.concatenate([[zvl[0]],zvr,zvl[::-1]]));
        
        for j in range(mC-1):
            
            val = [xvl[j],xvr[j],0.5* (xl[j] + xr[j]), 0.5* (xl[j+1] + xr[j+1])];
            COLOCX[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpx1 = val[1] - val[0];
            cpx2 = val[3] - val[2];
            
            val = [yvl[j],yvr[j],0.5* (yl[j] + yr[j]), 0.5* (yl[j+1] + yr[j+1])];
            COLOCY[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpy1 = val[1] - val[0];
            cpy2 = val[3] - val[2];
            
            val = [zvl[j],zvr[j],0.5* (zl[j] + zr[j]), 0.5* (zl[j+1] + zr[j+1])];
            COLOCZ[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpz1 = val[1] - val[0];
            cpz2 = val[3] - val[2];
            
            cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
            cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
            ds[i * (mC-1) + j] = cpmag;
            normal[:, i * (mC-1) + j] = cp/cpmag;
        dS[i] = sum(ds[i * (mC-1):(i+1) * (mC-1)]);
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
    aoa = m.pi * (flow.getAMax()+flow.getAMin())/180.;
    # Main lifting surfaces
    wing = ac.wing;
    htail = ac.htail;
    prop = ac.prop;
    rho0=1.225; #masse volumique à niveau de la mer           [kg/m^3]
    dT=-6.5;    #gradiente de temperature dans la troposphere [K/km]
    T0=288.15;  #Temperature à niveau de la mer               [K]
    g=9.80665;  #gravité                                      [m/s^2]
    Rair=287.1;   #Constante de l'air                           [m^2/(s^2*K)]
    h = flow.getH();    # flight altitude [km]
    V0 = flow.getV0(); # freestream velocity [m/s]
    rho = rho0 * (1. + dT*h/T0)**(- g/(Rair*dT*10**(-3)) - 1.); # air density
    Sh = m.pi * prop.getD()**2 *0.25;#Surface disque actuator          [m^2]
    nbE = len(prop.getD());
    OWU = prop.getYp()/np.abs(prop.getYp());
    for ii in range(nbE):
        if not(prop.OWU[ii]):
            OWU *= -1.;
    # Numerical parameters for  discretization
    tF = 2.;                                                                    # temps caractéristique de convection des vortex
    nbEch = 1.;                                                                # Nombre minimal de points de controle par rotation des lignes de courant/tourbillons
    mW = int(tF*nbEch/(2*m.pi)*max(prop.getOmega()));                           # discretisation of the wake, get correct direction of it behind parts
    times = np.linspace(0.,tF,mW);
    mC = wing.mC;                                                               # chordwise discretisation number of control point for the chord
    n = wing.getR()+htail.getR();                                               # spanwise discretisation number of panel
    
    # Recover the wing parameters
    
    # Panels' coordinates and main parameters (at c/4)
    xp = wing.getXP();
    yp = wing.getYP();
    zp = wing.getZP();
    
    cP = wing.getChord();
    tw = wing.getTwist();
    dih = wing.getDih();
    sw = wing.getSweepC4();
    
    # Panel bordes' coordinate and main parameters (at c/4)    
    x = wing.getX();
    y = wing.getY();
    z = wing.getZ();
    c = wing.getChordDist();
    twSec = wing.twSec;
    
    # Flaps
    
    xW = np.unique(np.concatenate([0.5*(np.cos(np.linspace(m.pi,0.,mC))+1.),[0.25]]));
    mC = len(xW);
    iC4W = np.where(xW == 0.25)[0][0];
    
    zW = np.zeros([mC,len(wing.getAF())],dtype = float);
    for ii in range(len(wing.getAF())):
        zW[:,ii]= camber(wing.getAF(ii),xW);
    
    if htail.bool:
        # Panel bordes' coordinate and main parameters (at c/4)
        x = np.concatenate([x,htail.getX()]);
        y = np.concatenate([y,htail.getY()]);
        z = np.concatenate([z,htail.getZ()]);
        c = np.concatenate([c,htail.getChordDist()]);
        twSec = np.concatenate([wing.twSec,htail.twSec]);
        
        # Panels' coordinates and main parameters (at c/4)
        xp = np.concatenate([xp,htail.getXP()]);
        yp = np.concatenate([yp,htail.getYP()]);
        zp = np.concatenate([zp,htail.getZP()]);
        
        cP = np.concatenate([cP,htail.getChord()]);
        tw = np.concatenate([tw,htail.getTwist()]);
        dih = np.concatenate([dih,htail.getDih()]);
        sw = np.concatenate([sw,htail.getSweepC4()]);
        
        # Elevator, Assumed to be as plain flaps
        cfT = htail.getCF();
        
        if cfT != 0: 
            xT = np.unique(np.concatenate([np.linspace(1.,1.-cfT,2),(1.-cfT)*0.5*(np.cos(np.linspace(m.pi,0.,mC-1))+1.)]));
            xT[abs((xT-0.25)) == np.min(abs(xT-0.25))] = 0.25;
        else:
            xT = 0.5*(np.cos(np.linspace(m.pi,0.,mC))+1.);
            xT[abs((xT-0.25)) == np.min(abs(xT-0.25))] = 0.25;
        iC4T = np.where(xT == 0.25)[0][0];
        zT = np.zeros([mC,len(htail.getAF())],dtype = float);
        for ii in range(len(htail.getAF())):
            zT[:,ii-1]= camber(htail.getAF(ii),xT);
    
    
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
    xvr = np.zeros(mC + mW,dtype = float);
    yvr = np.zeros(mC + mW,dtype = float);
    zvr = np.zeros(mC + mW,dtype = float);
    dzdx = np.zeros(mW-1,dtype = float);                       
    dydx = np.zeros(mW-1,dtype = float);                       
    for i in range(wing.getR()):
        camb = zW[:,wing.getAFI(i)]
        
        il = i;
        cl = c[il];
        twl = twSec[il];
        
        xl = (xW - 0.25) * cl + x[il];
        yl = y[il] * np.ones(mC);
        zl = camb * cl + z[il];
        center = np.array([xl[iC4W],yl[iC4W],zl[iC4W]]);
        alpha = 180./m.pi*twl;
        Rot = u.roty(alpha);
        for ii in range(mC):
            point = np.array([xl[ii],yl[ii],zl[ii]])-center;
            point = np.dot(Rot,point) + center;
            xl[ii] = point[0];
            yl[ii] = point[1];
            zl[ii] = point[2];
        
            
        xvl[:mC-1] = 0.75 * xl[:-1] + 0.25 * xl[1:];
        yvl[:mC-1] = 0.75 * yl[:-1] + 0.25 * yl[1:];
        zvl[:mC-1] = 0.75 * zl[:-1] + 0.25 * zl[1:];
        xvl[mC-1] = xvl[mC-2] + (xl[-1]-xl[-2]);
        yvl[mC-1:] = yvl[mC-2] + (yl[-1]-yl[-2]);                               # initial guess : stay straight at the end of the wing
        zvl[mC-1:] = zvl[mC-2] + (zl[-1]-zl[-2]);
        # End of chord vortex = begining of wake vortex
        
        
        # introduced effect of prop
        # attention : prendre en compte le fait que certaines ldc sont sous l'influence de 2 moteurs!
        
        centerPropY = prop.getYp() + (xvl[mC-1] - prop.getXp()) * m.tan(beta);
        centerPropZ = prop.getZp() + (xvl[mC-1] - prop.getXp()) * m.tan(aoa);
        vix = 0.;
        for j in range(nbE):
            d = m.sqrt((yvl[mC-1] - centerPropY[j])**2 + (zvl[mC-1] - centerPropZ[j])**2);
            rP = prop.rHub[j];
            D = prop.D[j];
            vitheta = 0.;
            theta0 = np.arctan2(zvl[mC-1] - centerPropZ[j],yvl[mC-1] - centerPropY[j]);
            if ((d >= rP) and (d <= D * 0.5) and prop.Omega[j] != 0.):
                vix += 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                vix2 = 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                a = vix2/V0;
                aprim = 0.5 * (1. - m.sqrt(abs(1.-4.*a*(1.+a)*(V0/(prop.Omega[j] * d))**2)));
                vitheta = OWU[j]*abs((aprim * 2. * prop.Omega[j] * d));
                Theta = times*vitheta/d + theta0;
                dY = np.cos(Theta[1:]) * d + centerPropY[j] - yvl[mC-1];
                dZ = np.sin(Theta[1:]) * d + centerPropZ[j] - zvl[mC-1] ;
                yvl[mC:-1] += dY;
                zvl[mC:-1] += dZ;
        xvl[mC-1:-1] = xvl[mC-1] + times * (V0+vix);
        xvl[-1] = 10. * wing.b;
        
        
        indiceFinLocalEffectCamber = np.where(xvl >= xvl[mC-1] + 2.5 * cl)[0][1]; 
        
        
        dzdxl = (zl[mC-1]-zl[mC-2])/(xl[mC-1]-xl[mC-2]);
        
        # Vérifie ça!
        dydx = V0/(V0+vix) *m.tan(beta) * (1.-np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])));
        dzdx = V0/(V0+vix) * dzdxl * np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])) \
            + V0/(V0+vix) * m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])));
        dY = np.zeros(mW+1);
        dZ = np.zeros(mW+1);
        for ii in range(1,indiceFinLocalEffectCamber-mC+1):
            dZ[ii] = dZ[(ii-1)] + dzdx[ii-1] * (xvl[mC-1+ii] - xvl[(mC-1+ii-1)]);
            dY[ii] = dY[(ii-1)] + dydx[ii-1] * (xvl[mC-1+ii] - xvl[(mC-1+ii-1)]);
        dZ[indiceFinLocalEffectCamber-mC+1:] = dZ[indiceFinLocalEffectCamber-mC] + m.tan(aoa) * (xvl[indiceFinLocalEffectCamber:] - xvl[indiceFinLocalEffectCamber-1]);
        dY[indiceFinLocalEffectCamber-mC+1:] = dY[indiceFinLocalEffectCamber-mC] + m.tan(beta) * (xvl[indiceFinLocalEffectCamber:] - xvl[indiceFinLocalEffectCamber-1]);
        yvl[mC-1:] += dY;
        zvl[mC-1:] += dZ;
        
        ir = i+1;
        cr = c[ir];
        twr = twSec[ir];
        
        xr = (xW - 0.25) * cr + x[ir];
        yr = y[ir] * np.ones(mC);
        zr = camb * cr + z[ir];
        center = np.array([xr[iC4W],yr[iC4W],zr[iC4W]]);
        alpha = 180./m.pi*twr;
        Rot = u.roty(alpha);
        for ii in range(0,mC):
            point = np.array([xr[ii],yr[ii],zr[ii]])-center;
            point = np.dot(Rot,point) + center;
            xr[ii] = point[0];
            yr[ii] = point[1];
            zr[ii] = point[2];
        
        xvr[:mC-1] = 0.75 * xr[:-1] + 0.25 * xr[1:];
        yvr[:mC-1] = 0.75 * yr[:-1] + 0.25 * yr[1:];
        zvr[:mC-1] = 0.75 * zr[:-1] + 0.25 * zr[1:];
        xvr[mC-1] = xvr[mC-2] + (xr[-1]-xr[-2]);
        yvr[mC-1:] = yvr[mC-2] + (yr[-1]-yr[-2]);
        zvr[mC-1:] = zvr[mC-2] + (zr[-1]-zr[-2]);
        # End of chord vortex = begining of wake vortex
        centerPropY = prop.getYp() + (xvr[mC-1] - prop.getXp()) * m.tan(beta);
        centerPropZ = prop.getZp() + (xvr[mC-1] - prop.getXp()) * m.tan(aoa);
        vix = 0.;
        for j in range(nbE):
            d = m.sqrt((yvr[mC-1] - centerPropY[j])**2 + (zvr[mC-1] - centerPropZ[j])**2);
            rP = prop.rHub[j];
            D = prop.D[j];
            vitheta = 0.;
            theta0 = np.arctan2(zvr[mC-1] - centerPropZ[j],yvr[mC-1] - centerPropY[j]);
            if ((d >= rP) and (d <= D * 0.5) and prop.Omega[j] != 0.):
                vix += 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                vix2 = 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                a = vix2/V0;
                aprim = 0.5 * (1. - m.sqrt(abs(1.-4.*a*(1.+a)*(V0/(prop.Omega[j] * d))**2)));
                vitheta = OWU[j]*abs((aprim * 2. * prop.Omega[j] * d));
                Theta = times*vitheta/d + theta0;
                dY = np.cos(Theta[1:]) * d + centerPropY[j] - yvr[mC-1];
                dZ = np.sin(Theta[1:]) * d + centerPropZ[j] - zvr[mC-1] ;
                yvr[mC:-1] += dY;
                zvr[mC:-1] += dZ;
        xvr[mC-1:-1] = xvr[mC-1] + times * (V0+vix);
        xvr[-1] = 10. * wing.b;
        
        
        indiceFinLocalEffectCamber = np.where(xvr >= xvr[mC-1] + 2.5 * cr)[0][1]; 
        
        
        dzdxr = (zr[mC-1]-zr[mC-2])/(xr[mC-1]-xr[mC-2]);
        
        # Vérifie ça!
        dydx = V0/(V0+vix) * m.tan(beta) * (1.-np.exp(-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1])));
        dzdx = V0/(V0+vix) * dzdxr * np.exp(-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1])) \
            + V0/(V0+vix) * m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1])));
        dY = np.zeros(mW+1);
        dZ = np.zeros(mW+1);
        for ii in range(1,indiceFinLocalEffectCamber-mC+1):
            dZ[ii] = dZ[(ii-1)] + dzdx[ii-1] * (xvr[mC-1+ii] - xvr[(mC-1+ii-1)]);
            dY[ii] = dY[(ii-1)] + dydx[ii-1] * (xvr[mC-1+ii] - xvr[(mC-1+ii-1)]);
        dZ[indiceFinLocalEffectCamber-mC+1:] = dZ[indiceFinLocalEffectCamber-mC] + m.tan(aoa) * (xvr[indiceFinLocalEffectCamber:] - xvr[indiceFinLocalEffectCamber-1]);
        dY[indiceFinLocalEffectCamber-mC+1:] = dY[indiceFinLocalEffectCamber-mC] + m.tan(beta) * (xvr[indiceFinLocalEffectCamber:] - xvr[indiceFinLocalEffectCamber-1]);
        yvr[mC-1:] += dY;
        zvr[mC-1:] += dZ;
        
        setTable(X,2*(mC+mW)+1,i,np.concatenate([[xvl[0]],xvr,xvl[::-1]]));
        setTable(Y,2*(mC+mW)+1,i,np.concatenate([[yvl[0]],yvr,yvl[::-1]]));
        setTable(Z,2*(mC+mW)+1,i,np.concatenate([[zvl[0]],zvr,zvl[::-1]]));
        
        for j in range(mC-1):
            
            val = [xvl[j],xvr[j],0.5* (xl[j] + xr[j]), 0.5* (xl[j+1] + xr[j+1])];
            COLOCX[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpx1 = val[1] - val[0];
            cpx2 = val[3] - val[2];
            
            val = [yvl[j],yvr[j],0.5* (yl[j] + yr[j]), 0.5* (yl[j+1] + yr[j+1])];
            COLOCY[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpy1 = val[1] - val[0];
            cpy2 = val[3] - val[2];
            
            val = [zvl[j],zvr[j],0.5* (zl[j] + zr[j]), 0.5* (zl[j+1] + zr[j+1])];
            COLOCZ[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpz1 = val[1] - val[0];
            cpz2 = val[3] - val[2];
            
            cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
            cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
            ds[i * (mC-1) + j] = cpmag;
            normal[:, i * (mC-1) + j] = cp/cpmag;
        dS[i] = sum(ds[i * (mC-1):(i+1) * (mC-1)]);
    
    for i in range(wing.getR(),wing.getR()+htail.getR()):
        iPT = i-wing.getR();
        camb = zT[:,htail.getAFI(iPT)]
        
        il = i+1;
        cl = c[il];
        twl = twSec[il];
        
        xl = (xT - 0.25) * cl + x[il];
        yl = y[il] * np.ones(mC);
        zl = camb * cl + z[il];
        center = np.array([xl[iC4T],yl[iC4T],zl[iC4T]]);
        alpha = 180./m.pi*twl;
        Rot = u.roty(alpha);
        for ii in range(mC):
            point = np.array([xl[ii],yl[ii],zl[ii]])-center;
            point = np.dot(Rot,point) + center;
            xl[ii] = point[0];
            yl[ii] = point[1];
            zl[ii] = point[2];
        if htail.getDF(iPT) != 0.:
            delta = htail.getDF(iPT);
            RotF = u.roty(delta);
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
        yvl[mC-1:] = yvl[mC-2] + (yl[-1]-yl[-2]);
        zvl[mC-1:] = zvl[mC-2] + (zl[-1]-zl[-2]);
        # End of chord vortex = begining of wake vortex
        centerPropY = prop.getYp() + (xvl[mC-1] - prop.getXp()) * m.tan(beta);
        centerPropZ = prop.getZp() + (xvl[mC-1] - prop.getXp()) * m.tan(aoa);
        vix = 0.;
        for j in range(nbE):
            d = m.sqrt((yvl[mC-1] - centerPropY[j])**2 + (zvl[mC-1] - centerPropZ[j])**2);
            rP = prop.rHub[j];
            D = prop.D[j];
            vitheta = 0.;
            theta0 = np.arctan2(zvl[mC-1] - centerPropZ[j],yvl[mC-1] - centerPropY[j]);
            if ((d >= rP) and (d <= D * 0.5) and prop.Omega[j] != 0.):
                vix += 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                vix2 = 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                a = vix2/V0;
                aprim = 0.5 * (1. - m.sqrt(abs(1.-4.*a*(1.+a)*(V0/(prop.Omega[j] * d))**2)));
                vitheta = OWU[j]*abs((aprim * 2. * prop.Omega[j] * d));
                Theta = times*vitheta/d + theta0;
                dY = np.cos(Theta[1:]) * d + centerPropY[j] - yvl[mC-1];
                dZ = np.sin(Theta[1:]) * d + centerPropZ[j] - zvl[mC-1] ;
                yvl[mC:-1] += dY;
                zvl[mC:-1] += dZ;
        xvl[mC-1:-1] = xvl[mC-1] + times * (V0+vix);
        xvl[-1] = 10. * wing.b;
        
        
        indiceFinLocalEffectCamber = np.where(xvl >= xvl[mC-1] + 2.5 * cl)[0][1]; 
        
        
        dzdxl = (zl[mC-1]-zl[mC-2])/(xl[mC-1]-xl[mC-2]);
        
        # Vérifie ça!
        dydx = V0/(V0+vix) * m.tan(beta) * (1.-np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])));
        dzdx = V0/(V0+vix)*dzdxl * np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])) \
            + V0/(V0+vix)*m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])));
        dY = np.zeros(mW+1);
        dZ = np.zeros(mW+1);
        for ii in range(1,indiceFinLocalEffectCamber-mC+1):
            dZ[ii] = dZ[(ii-1)] + dzdx[ii-1] * (xvl[mC-1+ii] - xvl[(mC-1+ii-1)]);
            dY[ii] = dY[(ii-1)] + dydx[ii-1] * (xvl[mC-1+ii] - xvl[(mC-1+ii-1)]);
        dZ[indiceFinLocalEffectCamber-mC+1:] = dZ[indiceFinLocalEffectCamber-mC] + m.tan(aoa) * (xvl[indiceFinLocalEffectCamber:] - xvl[indiceFinLocalEffectCamber-1]);
        dY[indiceFinLocalEffectCamber-mC+1:] = dY[indiceFinLocalEffectCamber-mC] + m.tan(beta) * (xvl[indiceFinLocalEffectCamber:] - xvl[indiceFinLocalEffectCamber-1]);
        yvl[mC-1:] += dY;
        zvl[mC-1:] += dZ;
        
        ir = i+2;
        cr = c[ir];
        twr = twSec[ir];
        
        xr = (xT - 0.25) * cr + x[ir];
        yr = y[ir] * np.ones(mC);
        zr = camb * cr + z[ir];
        center = np.array([xr[iC4T],yr[iC4T],zr[iC4T]]);
        alpha = 180./m.pi*twr;
        Rot = u.roty(alpha);
        for ii in range(0,mC):
            point = np.array([xr[ii],yr[ii],zr[ii]])-center;
            point = np.dot(Rot,point) + center;
            xr[ii] = point[0];
            yr[ii] = point[1];
            zr[ii] = point[2];
        if htail.getDF(iPT) != 0.:
            delta = htail.getDF(iPT);
            RotF = u.roty(delta);
            center = np.array([xr[-2],yr[-2],zr[-2]]);
            point = np.array([xr[-1],yr[-1],zr[-1]])-center;
            point = np.dot(RotF,point) + center;
            xr[-1] = point[0];
            yr[-1] = point[1];
            zr[-1] = point[2];
            
        xvr[:mC-1] = 0.75 * xr[:-1] + 0.25 * xr[1:];
        yvr[:mC-1] = 0.75 * yr[:-1] + 0.25 * yr[1:];
        zvr[:mC-1] = 0.75 * zr[:-1] + 0.25 * zr[1:];
        xvr[mC-1] = xvr[mC-2] + (xr[-1]-xr[-2]);
        yvr[mC-1:] = yvr[mC-2] + (yr[-1]-yr[-2]);
        zvr[mC-1:] = zvr[mC-2] + (zr[-1]-zr[-2]);
        # End of chord vortex = begining of wake vortex
        centerPropY = prop.getYp() + (xvr[mC-1] - prop.getXp()) * m.tan(beta);
        centerPropZ = prop.getZp() + (xvr[mC-1] - prop.getXp()) * m.tan(aoa);
        vix = 0.;
        for j in range(nbE):
            d = m.sqrt((yvr[mC-1] - centerPropY[j])**2 + (zvr[mC-1] - centerPropZ[j])**2);
            rP = prop.rHub[j];
            D = prop.D[j];
            vitheta = 0.;
            theta0 = np.arctan2(zvr[mC-1] - centerPropZ[j],yvr[mC-1] - centerPropY[j]);
            if ((d >= rP) and (d <= D * 0.5) and prop.Omega[j] != 0.):
                vix += 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                vix2 = 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                a = vix2/V0;
                aprim = 0.5 * (1. - m.sqrt(abs(1.-4.*a*(1.+a)*(V0/(prop.Omega[j] * d))**2)));
                vitheta = OWU[j]*abs((aprim * 2. * prop.Omega[j] * d));
                Theta = times*vitheta/d + theta0;
                dY = np.cos(Theta[1:]) * d + centerPropY[j] - yvr[mC-1];
                dZ = np.sin(Theta[1:]) * d + centerPropZ[j] - zvr[mC-1] ;
                yvr[mC:-1] += dY;
                zvr[mC:-1] += dZ;
        xvr[mC-1:-1] = xvr[mC-1] + times * (V0+vix);
        xvr[-1] = 10. * wing.b;
        
        
        indiceFinLocalEffectCamber = np.where(xvr >= xvr[mC-1] + 2.5 * cr)[0][1]; 
        
        
        dzdxr = (zr[mC-1]-zr[mC-2])/(xr[mC-1]-xr[mC-2]);
        
#        # Vérifie ça!
        dydx = V0/(V0+vix)*m.tan(beta) * (1.-np.exp(-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1])));
        dzdx = V0/(V0+vix)*dzdxr * np.exp(-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1])) \
            + V0/(V0+vix)*m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1])));
        dY = np.zeros(mW+1);
        dZ = np.zeros(mW+1);
        for ii in range(1,indiceFinLocalEffectCamber-mC+1):
            dZ[ii] = dZ[(ii-1)] + dzdx[ii-1] * (xvr[mC-1+ii] - xvr[(mC-1+ii-1)]);
            dY[ii] = dY[(ii-1)] + dydx[ii-1] * (xvr[mC-1+ii] - xvr[(mC-1+ii-1)]);
        dZ[indiceFinLocalEffectCamber-mC+1:] = dZ[indiceFinLocalEffectCamber-mC] + m.tan(aoa) * (xvr[indiceFinLocalEffectCamber:] - xvr[indiceFinLocalEffectCamber-1]);
        dY[indiceFinLocalEffectCamber-mC+1:] = dY[indiceFinLocalEffectCamber-mC] + m.tan(beta) * (xvr[indiceFinLocalEffectCamber:] - xvr[indiceFinLocalEffectCamber-1]);
        yvr[mC-1:] += dY;
        zvr[mC-1:] += dZ;
        setTable(X,2*(mC+mW)+1,i,np.concatenate([[xvl[0]],xvr,xvl[::-1]]));
        setTable(Y,2*(mC+mW)+1,i,np.concatenate([[yvl[0]],yvr,yvl[::-1]]));
        setTable(Z,2*(mC+mW)+1,i,np.concatenate([[zvl[0]],zvr,zvl[::-1]]));
        
        for j in range(mC-1):
            
            val = [xvl[j],xvr[j],0.5* (xl[j] + xr[j]), 0.5* (xl[j+1] + xr[j+1])];
            COLOCX[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpx1 = val[1] - val[0];
            cpx2 = val[3] - val[2];
            
            val = [yvl[j],yvr[j],0.5* (yl[j] + yr[j]), 0.5* (yl[j+1] + yr[j+1])];
            COLOCY[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpy1 = val[1] - val[0];
            cpy2 = val[3] - val[2];
            
            val = [zvl[j],zvr[j],0.5* (zl[j] + zr[j]), 0.5* (zl[j+1] + zr[j+1])];
            COLOCZ[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpz1 = val[1] - val[0];
            cpz2 = val[3] - val[2];
            
            cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
            cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
            ds[i * (mC-1) + j] = cpmag;
            normal[:, i * (mC-1) + j] = cp/cpmag;
        dS[i] = sum(ds[i * (mC-1):(i+1) * (mC-1)]);
    select = np.zeros([n,n * (mC-1)]); # rechercher intensité du dernier vortex uniquement
    select2 = np.zeros([n * (mC-1),n]); # pour chaque paneau sur même section y, même velocity triangle
    select3 = np.zeros([n + len(ac.prop.D),n * (mC-1) + len(ac.prop.D)]); # 
    for i in range(n):
        select[i,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i] = 1.;
        select3[i,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i];
    if ac.prop.bool:
        select3[-len(ac.prop.D):,-len(ac.prop.D):] = np.eye(len(ac.prop.D));
#    plt.plot(Y,X),plt.axis([-8,8,-1,15]),plt.show()
#    plt.plot(Y,Z),plt.axis([-8,8,-1,15]),plt.show()
#    plt.plot(X,Z),plt.axis([-1,15,-1,15]),plt.show()
#    return
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
    aoa = m.pi * (flow.getAMax()+flow.getAMin())/180.;
    # Main lifting surfaces
    wing = ac.wing;
    htail = ac.htail;
    cf = wing.getCF();
    
    # Numerical parameters for  discretization
    
    mC = wing.mC;                                                               # chordwise discretisation number of control point for the chord
    mW = flow.mW;                                                               # discretisation of the wake, get correct direction of it behind parts
    n = 2*wing.getR()+htail.getR();                                               # spanwise discretisation number of panel# Recover the wing parameters
    
    # Panels' coordinates and main parameters (at c/4)
    xp = wing.getXP();
    yp = wing.getYP();
    zp = wing.getZP();
    
    cP = wing.getChord();
    tw = wing.getTwist();
    dih = wing.getDih();
    sw = wing.getSweepC4();
    
    # Panel bordes' coordinate and main parameters (at c/4)    
    x = wing.getX();
    y = wing.getY();
    z = wing.getZ();
    c = wing.getChordDist();
    twSec = wing.twSec;
    
    xW = np.unique(np.concatenate([(1.-cf)*0.5*(np.cos(np.linspace(m.pi,0.,mC))+1.),[0.25]]));
    mC = len(xW);
    iC4W = np.where(xW == 0.25)[0][0];                                       # Indice du début du flaps / fin de main element
    zW = np.zeros([mC,len(wing.getAF())],dtype = float);
    for ii in range(len(wing.getAF())):
        zW[:,ii]= camber(wing.getAF(ii),xW);

    xF = 1. - cf * np.linspace(1.,0,mC);
    zF = np.zeros([mC,len(wing.getAF())],dtype = float);
    for ii in range(len(wing.getAF())):
        zF[:,ii] = camber(wing.getAF(ii),xF);
    
    if htail.bool:
        # Panel bordes' coordinate and main parameters (at c/4)
        x = np.concatenate([x,htail.getX()]);
        y = np.concatenate([y,htail.getY()]);
        z = np.concatenate([z,htail.getZ()]);
        c = np.concatenate([c,htail.getChordDist()]);
        twSec = np.concatenate([wing.twSec,htail.twSec]);
        
        # Panels' coordinates and main parameters (at c/4)
        xp = np.concatenate([xp,htail.getXP()]);
        yp = np.concatenate([yp,htail.getYP()]);
        zp = np.concatenate([zp,htail.getZP()]);
        
        cP = np.concatenate([cP,htail.getChord()]);
        tw = np.concatenate([tw,htail.getTwist()]);
        dih = np.concatenate([dih,htail.getDih()]);
        sw = np.concatenate([sw,htail.getSweepC4()]);
        
        # Elevator, Assumed to be as plain flaps
        cfT = htail.getCF();
        
        if cfT != 0: 
            xT = np.unique(np.concatenate([np.linspace(1.,1.-cfT,2),(1.-cfT)*0.5*(np.cos(np.linspace(m.pi,0.,mC-1))+1.)]));
            xT[abs((xT-0.25)) == np.min(abs(xT-0.25))] = 0.25;
        else:
            xT = 0.5*(np.cos(np.linspace(m.pi,0.,mC))+1.);
            xT[abs((xT-0.25)) == np.min(abs(xT-0.25))] = 0.25;
        iC4T = np.where(xT == 0.25)[0][0];
        zT = np.zeros([mC,len(htail.getAF())],dtype = float);
        for ii in range(len(htail.getAF())):
            zT[:,ii-1]= camber(htail.getAF(ii),xT);
            

    #generate grid corner coordinates
    # generate collocation points and normal : where tangency condition is
    # satisfied. Distance from bound vortex depends on the sectional lift
    # curve slope : (dist/localChord) = clAlphas/(4*pi)
    
    X = np.zeros(n * (2 * (mC + mW)+1),dtype = float);
    Y = np.zeros(n * (2 * (mC + mW)+1),dtype = float);                                   # initialization
    Z = np.zeros(n * (2 * (mC + mW)+1),dtype = float);
    COLOCX=np.zeros((mC-1)*n);
    COLOCY=np.zeros((mC-1)*n);
    COLOCZ=np.zeros((mC-1)*n);
    normal = np.zeros([3,(mC-1)*n]);
    coef = 0.25+cla*0.25/m.pi;
    ds = np.zeros((mC-1)*n);                                                    # vector of area of any panel
    dS = np.zeros(wing.r+htail.r);                                                           # vector of area of a spanwise section
    
    xvl = np.zeros(mC + mW,dtype = float);
    yvl = np.zeros(mC + mW,dtype = float);
    zvl = np.zeros(mC + mW,dtype = float);
    xvr = np.zeros(mC + mW,dtype = float);
    yvr = np.zeros(mC + mW,dtype = float);
    zvr = np.zeros(mC + mW,dtype = float);
    xvlf = np.zeros(mC + mW,dtype = float);
    yvlf = np.zeros(mC + mW,dtype = float);
    zvlf = np.zeros(mC + mW,dtype = float);
    xvrf = np.zeros(mC + mW,dtype = float);
    yvrf = np.zeros(mC + mW,dtype = float);
    zvrf = np.zeros(mC + mW,dtype = float);
    dzdx = np.zeros(mW-1,dtype = float);       
    dzdxf = np.zeros(mW-1,dtype = float);                       
    for i in range(wing.getR()):
        camb = zW[:,wing.getAFI(i)]
        cambF = zF[:,wing.getAFI(i)];
        il = i;
        cl = c[il];
        twl = twSec[il];
        
        xl = (xW - 0.25) * cl + x[il];
        yl = y[il] * np.ones(mC);
        zl = camb * cl + z[il];
        
        xlf = (xF - 0.25) * cl + x[il];
        ylf = y[il] * np.ones(mC);
        zlf = cambF * cl + z[il];
        center = np.array([xl[iC4W],yl[iC4W],zl[iC4W]]);
        alpha = 180./m.pi*twl;
        Rot = u.roty(alpha);
        for ii in range(mC):
            point = np.array([xl[ii],yl[ii],zl[ii]])-center;
            point = np.dot(Rot,point) + center;
            xl[ii] = point[0];
            yl[ii] = point[1];
            zl[ii] = point[2];
            pointf = np.array([xlf[ii],ylf[ii],zlf[ii]])-center;
            pointf = np.dot(Rot,pointf) + center;
            xlf[ii] = pointf[0] - 0.02 * cl;
            ylf[ii] = pointf[1];
            zlf[ii] = pointf[2] - 0.02 * cl;

        centerf = np.array([xlf[0],ylf[0],zlf[0]]);
        delta = wing.getDF(i);
        Rotf = u.roty(delta);
        for ii in range(mC):
            pointf = np.array([xlf[ii],ylf[ii],zlf[ii]])-centerf;
            pointf = np.dot(Rotf,pointf) + centerf;
            xlf[ii] = pointf[0];
            ylf[ii] = pointf[1];
            zlf[ii] = pointf[2];  
        xvl[:mC-1] = 0.75 * xl[:-1] + 0.25 * xl[1:];
        yvl[:mC-1] = 0.75 * yl[:-1] + 0.25 * yl[1:];
        zvl[:mC-1] = 0.75 * zl[:-1] + 0.25 * zl[1:];
        xvl[mC-1] = xvl[mC-2] + (xl[-1]-xl[-2]);
        yvl[mC-1] = yvl[mC-2] + (yl[-1]-yl[-2]);
        zvl[mC-1] = zvl[mC-2] + (zl[-1]-zl[-2]);
        
        xvlf[:mC-1] = 0.75 * xlf[:-1] + 0.25 * xlf[1:];
        yvlf[:mC-1] = 0.75 * ylf[:-1] + 0.25 * ylf[1:];
        zvlf[:mC-1] = 0.75 * zlf[:-1] + 0.25 * zlf[1:];
        xvlf[mC-1] = xvlf[mC-2] + (xlf[-1]-xlf[-2]);
        yvlf[mC-1] = yvlf[mC-2] + (ylf[-1]-ylf[-2]);
        zvlf[mC-1] = zvlf[mC-2] + (zlf[-1]-zlf[-2]);
        
        
        # End of chord vortex = begining of wake vortex
        Wake = 1.;
        xvl[mC:-1] = xvl[mC-1] + Wake * 2.5 * cl * (1.+np.array(range(mW-1),dtype = float))/mW;
        xvl[-1] = 10. * wing.b * Wake + (1.- Wake) * xvl[mC-1];
        dzdxl = (zl[mC-1]-zl[mC-2])/(xl[mC-1]-xl[mC-2]);
        dydx = m.tan(beta) * (1.-np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])));
        dzdx = dzdxl * np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])) \
            + m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])));
        for ii in range(mW-1):
            zvl[mC+ii] = zvl[mC+(ii-1)] + dzdx[ii] * (xvl[mC+ii] - xvl[mC+(ii-1)]);
            yvl[mC+ii] = yvl[mC+(ii-1)] + dydx[ii] * (xvl[mC+ii] - xvl[mC+(ii-1)]);
        zvl[-1] = zvl[-2] + m.tan(aoa) * (xvl[-1] - xvl[-2]);
        yvl[-1] = yvl[-2] + m.tan(beta) * (xvl[-1] - xvl[-2]);
        
        xvlf[mC:-1] = xvlf[mC-1] + 2.5 * cl * (1.+np.array(range(mW-1),dtype = float))/mW;
        xvlf[-1] = 10. * wing.b;
        dzdxlf = (zlf[mC-1]-zlf[mC-2])/(xlf[mC-1]-xlf[mC-2]);
        dzdxf = dzdxlf * np.exp(-3.*(np.array(xvlf[mC:-1] - xvlf[mC]))/(xvlf[-2] - xvlf[mC])) \
            + m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvlf[mC:-1] - xvlf[mC]))/(xvlf[-2] - xvlf[mC])));
        dydxf = m.tan(beta) * (1.-np.exp(-3.*(np.array(xvlf[mC:-1] - xvlf[mC]))/(xvlf[-2] - xvlf[mC])));
        for ii in range(mW-1):
            zvlf[mC+ii] = zvlf[mC+(ii-1)] + dzdxf[ii] * (xvlf[mC+ii] - xvlf[mC+(ii-1)]);
            yvlf[mC+ii] = yvlf[mC+(ii-1)] + dydxf[ii] * (xvlf[mC+ii] - xvlf[mC+(ii-1)]);
        zvlf[-1] = zvlf[-2] + m.tan(aoa) * (xvlf[-1] - xvlf[-2]);
        yvlf[-1] = yvlf[-2] + m.tan(beta) * (xvlf[-1] - xvlf[-2]);
        
        
        
        ## Right Part
        
        ir = i+1;
        cr = c[ir];
        twr = twSec[ir];
        
        xr = (xW - 0.25) * cr + x[ir];
        yr = y[ir] * np.ones(mC);
        zr = camb * cr + z[ir];
        
        xrf = (xF - 0.25) * cr + x[ir];
        yrf = y[ir] * np.ones(mC);
        zrf = cambF * cr + z[ir];
        
        center = np.array([xr[iC4W],yr[iC4W],zr[iC4W]]);
        alpha = 180./m.pi*twr;
        Rot = u.roty(alpha);
        for ii in range(0,mC):
            point = np.array([xr[ii],yr[ii],zr[ii]])-center;
            point = np.dot(Rot,point) + center;
            xr[ii] = point[0];
            yr[ii] = point[1];
            zr[ii] = point[2];
            pointf = np.array([xrf[ii],yrf[ii],zrf[ii]])-center;
            pointf = np.dot(Rot,pointf) + center;
            xrf[ii] = pointf[0] - 0.02 * cr;
            yrf[ii] = pointf[1];
            zrf[ii] = pointf[2] - 0.02 * cr;
        
        centerf = np.array([xrf[0],yrf[0],zrf[0]]);
        for ii in range(mC):
            pointf = np.array([xrf[ii],yrf[ii],zrf[ii]])-centerf;
            pointf = np.dot(Rotf,pointf) + centerf;
            xrf[ii] = pointf[0];
            yrf[ii] = pointf[1];
            zrf[ii] = pointf[2];  
        
        xvr[:mC-1] = 0.75 * xr[:-1] + 0.25 * xr[1:];
        yvr[:mC-1] = 0.75 * yr[:-1] + 0.25 * yr[1:];
        zvr[:mC-1] = 0.75 * zr[:-1] + 0.25 * zr[1:];
        xvr[mC-1] = xvr[mC-2] + (xr[-1]-xr[-2]);
        yvr[mC-1] = yvr[mC-2] + (yr[-1]-yr[-2]);
        zvr[mC-1] = zvr[mC-2] + (zr[-1]-zr[-2]);
        
        xvrf[:mC-1] = 0.75 * xrf[:-1] + 0.25 * xrf[1:];
        yvrf[:mC-1] = 0.75 * yrf[:-1] + 0.25 * yrf[1:];
        zvrf[:mC-1] = 0.75 * zrf[:-1] + 0.25 * zrf[1:];
        xvrf[mC-1] = xvrf[mC-2] + (xrf[-1]-xrf[-2]);
        yvrf[mC-1] = yvrf[mC-2] + (yrf[-1]-yrf[-2]);
        zvrf[mC-1] = zvrf[mC-2] + (zrf[-1]-zrf[-2]);
        
        # End of chord vortex = begining of wake vortex
        xvr[mC:-1] = xvr[mC-1] + Wake * 2.5 * cr * (1.+np.array(range(mW-1),dtype = float))/mW;
        xvr[-1] = 10. * wing.b * Wake + (1.- Wake) * xvr[mC-1];
        dzdxr = (zr[mC-1]-zr[mC-2])/(xr[mC-1]-xr[mC-2]);
        dydx = m.tan(beta) * (1.-np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])));
        dzdx = dzdxr * np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])) \
            + m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])));
        for ii in range(mW-1):
            zvr[mC+ii] = zvr[mC+(ii-1)] + dzdx[ii] * (xvr[mC+ii] - xvr[mC+(ii-1)]);
            yvr[mC+ii] = yvr[mC+(ii-1)] + dydx[ii] * (xvr[mC+ii] - xvr[mC+(ii-1)]);
        zvr[-1] = zvr[-2] + m.tan(aoa) * (xvr[-1] - xvr[-2]);
        yvr[-1] = yvr[-2] + m.tan(beta) * (xvr[-1] - xvr[-2]);
        
        
        xvrf[mC:-1] = xvrf[mC-1] + 2.5 * cr * (1.+np.array(range(mW-1),dtype = float))/mW;
        xvrf[-1] = 10. * wing.b;
        yvrf[mC:] = yvrf[mC-1] + m.tan(beta) * (xvrf[mC:] - xvrf[mC-1]);
        dzdxrf = (zrf[mC-1]-zrf[mC-2])/(xrf[mC-1]-xrf[mC-2]);
        dydx = m.tan(beta) * (1.-np.exp(-3.*(np.array(xvrf[mC:-1] - xvrf[mC]))/(xvrf[-2] - xvrf[mC])));
        dzdx = dzdxrf * np.exp(-3.*(np.array(xvrf[mC:-1] - xvrf[mC]))/(xvrf[-2] - xvrf[mC])) \
            + m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvrf[mC:-1] - xvrf[mC]))/(xvrf[-2] - xvrf[mC])));
        for ii in range(mW-1):
            zvrf[mC+ii] = zvrf[mC+(ii-1)] + dzdxf[ii] * (xvrf[mC+ii] - xvrf[mC+(ii-1)]);
            yvrf[mC+ii] = yvrf[mC+(ii-1)] + dydxf[ii] * (xvrf[mC+ii] - xvrf[mC+(ii-1)]);
        zvrf[-1] = zvrf[-2] + m.tan(aoa) * (xvrf[-1] - xvrf[-2]);
        yvrf[-1] = yvrf[-2] + m.tan(beta) * (xvrf[-1] - xvrf[-2]);
        
        
        setTable(X,2*(mC+mW)+1,i,np.concatenate([[xvl[0]],xvr,xvl[::-1]]));
        setTable(Y,2*(mC+mW)+1,i,np.concatenate([[yvl[0]],yvr,yvl[::-1]]));
        setTable(Z,2*(mC+mW)+1,i,np.concatenate([[zvl[0]],zvr,zvl[::-1]]));
        
        setTable(X,2*(mC+mW)+1,wing.r+i,np.concatenate([[xvlf[0]],xvrf,xvlf[::-1]]));
        setTable(Y,2*(mC+mW)+1,wing.r+i,np.concatenate([[yvlf[0]],yvrf,yvlf[::-1]]));
        setTable(Z,2*(mC+mW)+1,wing.r+i,np.concatenate([[zvlf[0]],zvrf,zvlf[::-1]]));
        
        
        for j in range(mC-1):
            
            val = [xvl[j],xvr[j],0.5* (xl[j] + xr[j]), 0.5* (xl[j+1] + xr[j+1])];
            COLOCX[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpx1 = val[1] - val[0];
            cpx2 = val[3] - val[2];
            
            val = [yvl[j],yvr[j],0.5* (yl[j] + yr[j]), 0.5* (yl[j+1] + yr[j+1])];
            COLOCY[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpy1 = val[1] - val[0];
            cpy2 = val[3] - val[2];
            
            val = [zvl[j],zvr[j],0.5* (zl[j] + zr[j]), 0.5* (zl[j+1] + zr[j+1])];
            COLOCZ[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpz1 = val[1] - val[0];
            cpz2 = val[3] - val[2];
            
            cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
            cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
            ds[i * (mC-1) + j] = cpmag;
            normal[:, i * (mC-1) + j] = cp/cpmag;
            
            val = [xvlf[j],xvrf[j],0.5* (xlf[j] + xrf[j]), 0.5* (xlf[j+1] + xrf[j+1])];
            COLOCX[(i+wing.r) * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpx1 = val[1] - val[0];
            cpx2 = val[3] - val[2];
            
            val = [yvlf[j],yvrf[j],0.5* (ylf[j] + yrf[j]), 0.5* (ylf[j+1] + yrf[j+1])];
            COLOCY[(i+wing.r) * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpy1 = val[1] - val[0];
            cpy2 = val[3] - val[2];
            
            val = [zvlf[j],zvrf[j],0.5* (zlf[j] + zrf[j]), 0.5* (zlf[j+1] + zrf[j+1])];
            COLOCZ[(i+wing.r) * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpz1 = val[1] - val[0];
            cpz2 = val[3] - val[2];
            
            cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
            cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
            ds[(i + wing.r) * (mC-1) + j] = cpmag;
            normal[:, (i + wing.r) * (mC-1) + j] = cp/cpmag;
            
        dS[i] = sum(ds[i * (mC-1):(i+1) * (mC-1)]) + sum(ds[(i+wing.r) * (mC-1):(i+wing.r+1) * (mC-1)]);
    for i in range(2*wing.getR(),2*wing.getR()+htail.getR()):
        iPT = i- 2 * wing.getR();
        camb = zT[:,htail.getAFI(iPT)]
        
        il = i+1 - wing.r;
        cl = c[il];
        twl = twSec[il];
        
        xl = (xT - 0.25) * cl + x[il];
        yl = y[il] * np.ones(mC);
        zl = camb * cl + z[il];
        center = np.array([xl[iC4T],yl[iC4T],zl[iC4T]]);
        alpha = 180./m.pi*twl;
        Rot = u.roty(alpha);
        for ii in range(mC):
            point = np.array([xl[ii],yl[ii],zl[ii]])-center;
            point = np.dot(Rot,point) + center;
            xl[ii] = point[0];
            yl[ii] = point[1];
            zl[ii] = point[2];
        if htail.getDF(iPT) != 0.:
            delta = htail.getDF(iPT);
            RotF = u.roty(delta);
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
        xvl[-1] = 10. * wing.b;
        dzdxl = (zl[mC-1]-zl[mC-2])/(xl[mC-1]-xl[mC-2]);
        dydx = m.tan(beta) * (1.-np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])));
        dzdx = dzdxl * np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])) \
            + m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvl[mC:-1] - xvl[mC]))/(xvl[-2] - xvl[mC])));
        for ii in range(mW-1):
            zvl[mC+ii] = zvl[mC+(ii-1)] + dzdx[ii] * (xvl[mC+ii] - xvl[mC+(ii-1)]);
            yvl[mC+ii] = yvl[mC+(ii-1)] + dydx[ii] * (xvl[mC+ii] - xvl[mC+(ii-1)]);
        zvl[-1] = zvl[-2] + m.tan(aoa) * (xvl[-1] - xvl[-2]);
        yvl[-1] = yvl[-2] + m.tan(beta) * (xvl[-1] - xvl[-2]);
        
        ir = i+2 - wing.r;
        cr = c[ir];
        twr = twSec[ir];
        
        xr = (xT - 0.25) * cr + x[ir];
        yr = y[ir] * np.ones(mC);
        zr = camb * cr + z[ir];
        center = np.array([xr[iC4T],yr[iC4T],zr[iC4T]]);
        alpha = 180./m.pi*twr;
        Rot = u.roty(alpha);
        for ii in range(0,mC):
            point = np.array([xr[ii],yr[ii],zr[ii]])-center;
            point = np.dot(Rot,point) + center;
            xr[ii] = point[0];
            yr[ii] = point[1];
            zr[ii] = point[2];
        if htail.getDF(iPT) != 0.:
            delta = htail.getDF(iPT);
            RotF = u.roty(delta);
            center = np.array([xr[-2],yr[-2],zr[-2]]);
            point = np.array([xr[-1],yr[-1],zr[-1]])-center;
            point = np.dot(RotF,point) + center;
            xr[-1] = point[0];
            yr[-1] = point[1];
            zr[-1] = point[2];
            
        xvr[:mC-1] = 0.75 * xr[:-1] + 0.25 * xr[1:];
        yvr[:mC-1] = 0.75 * yr[:-1] + 0.25 * yr[1:];
        zvr[:mC-1] = 0.75 * zr[:-1] + 0.25 * zr[1:];
        xvr[mC-1] = xvr[mC-2] + (xr[-1]-xr[-2]);
        yvr[mC-1] = yvr[mC-2] + (yr[-1]-yr[-2]);
        zvr[mC-1] = zvr[mC-2] + (zr[-1]-zr[-2]);
        # End of chord vortex = begining of wake vortex
        xvr[mC:-1] = xvr[mC-1] + 2.5 * cr * (1.+np.array(range(mW-1),dtype = float))/mW;
        xvr[-1] = 10. * wing.b;
        dzdxr = (zr[mC-1]-zr[mC-2])/(xr[mC-1]-xr[mC-2]);
        dydx = m.tan(beta) * (1.-np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])));
        dzdx = dzdxr * np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])) \
            + m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvr[mC:-1] - xvr[mC]))/(xvr[-2] - xvr[mC])));
        for ii in range(mW-1):
            zvr[mC+ii] = zvr[mC+(ii-1)] + dzdx[ii] * (xvr[mC+ii] - xvr[mC+(ii-1)]);
            yvr[mC+ii] = yvr[mC+(ii-1)] + dydx[ii] * (xvr[mC+ii] - xvr[mC+(ii-1)]);
        zvr[-1] = zvr[-2] + m.tan(aoa) * (xvr[-1] - xvr[-2]);
        yvr[-1] = yvr[-2] + m.tan(beta) * (xvr[-1] - xvr[-2]);
        setTable(X,2*(mC+mW)+1,i,np.concatenate([[xvl[0]],xvr,xvl[::-1]]));
        setTable(Y,2*(mC+mW)+1,i,np.concatenate([[yvl[0]],yvr,yvl[::-1]]));
        setTable(Z,2*(mC+mW)+1,i,np.concatenate([[zvl[0]],zvr,zvl[::-1]]));
        
        for j in range(mC-1):
            
            val = [xvl[j],xvr[j],0.5* (xl[j] + xr[j]), 0.5* (xl[j+1] + xr[j+1])];
            COLOCX[i * (mC-1) + j] = val[2] * (1.-coef[i - wing.r]) + val[3] * coef[i - wing.r];
            cpx1 = val[1] - val[0];
            cpx2 = val[3] - val[2];
            
            val = [yvl[j],yvr[j],0.5* (yl[j] + yr[j]), 0.5* (yl[j+1] + yr[j+1])];
            COLOCY[i * (mC-1) + j] = val[2] * (1.-coef[i - wing.r]) + val[3] * coef[i - wing.r];
            cpy1 = val[1] - val[0];
            cpy2 = val[3] - val[2];
            
            val = [zvl[j],zvr[j],0.5* (zl[j] + zr[j]), 0.5* (zl[j+1] + zr[j+1])];
            COLOCZ[i * (mC-1) + j] = val[2] * (1.-coef[i - wing.r]) + val[3] * coef[i - wing.r];
            cpz1 = val[1] - val[0];
            cpz2 = val[3] - val[2];
            
            cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
            cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
            ds[i * (mC-1) + j] = cpmag;
            normal[:, i * (mC-1) + j] = cp/cpmag;
        dS[i-wing.r] = sum(ds[i * (mC-1):(i+1) * (mC-1)]);
    select = np.zeros([wing.r + htail.r,n * (mC-1)]); # rechercher intensité du dernier vortex uniquement
    select2 = np.zeros([n * (mC-1),wing.r + htail.r]); # pour chaque paneau sur même section y, même velocity triangle
    select3 = np.zeros([wing.r + htail.r + len(ac.prop.D),n * (mC-1) + len(ac.prop.D)]); # 
    for i in range(wing.r):
        select[i,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i] = 1.;
        select3[i,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i];
    for i in range(wing.r,n):
        select[i-wing.r,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i - wing.r] = 1.;
        select3[i - wing.r,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i-wing.r];
    if ac.prop.bool:
        select3[-len(ac.prop.D):,-len(ac.prop.D):] = np.eye(len(ac.prop.D));
    Ao,Vxo,Vyo,Vzo = ICM_F(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac,n,mC,mW);
    invA = np.linalg.inv(Ao);
    A = invA;
    Vx = np.dot(select3,Vxo);
    Vy = np.dot(select3,Vyo);
    Vz = np.dot(select3,Vzo);
    return A,normal,Vx,Vy,Vz,select,select2;


def getGridF_Engines(flow,ac,cla):
    
    beta = - flow.getBeta()*m.pi/180.;
    aoa = m.pi * (flow.getAMax()+flow.getAMin())/180.;
    # Main lifting surfaces
    wing = ac.wing;
    htail = ac.htail;
    prop = ac.prop;
    cf = wing.getCF();
    rho0=1.225; #masse volumique à niveau de la mer           [kg/m^3]
    dT=-6.5;    #gradiente de temperature dans la troposphere [K/km]
    T0=288.15;  #Temperature à niveau de la mer               [K]
    g=9.80665;  #gravité                                      [m/s^2]
    Rair=287.1;   #Constante de l'air                           [m^2/(s^2*K)]
    h = flow.getH();    # flight altitude [km]
    V0 = flow.getV0(); # freestream velocity [m/s]
    rho = rho0 * (1. + dT*h/T0)**(- g/(Rair*dT*10**(-3)) - 1.); # air density
    Sh = m.pi * prop.getD()**2 *0.25;#Surface disque actuator          [m^2]
    nbE = len(prop.getD());
    OWU = prop.getYp()/np.abs(prop.getYp());
    for ii in range(nbE):
        if not(prop.OWU[ii]):
            OWU *= -1.;
    # Numerical parameters for  discretization
    tF = 2.;                                                                    # temps caractéristique de convection des vortex
    nbEch = 1.;                                                                # Nombre minimal de points de controle par rotation des lignes de courant/tourbillons
    mW = int(tF*nbEch/(2*m.pi)*max(prop.getOmega()));                           # discretisation of the wake, get correct direction of it behind parts
    times = np.linspace(0.,tF,mW);
    mC = wing.mC;                                                               # chordwise discretisation number of control point for the chord
    n = 2*wing.getR()+htail.getR();                                               # spanwise discretisation number of panel

     
    # Panels' coordinates and main parameters (at c/4)
    xp = wing.getXP();
    yp = wing.getYP();
    zp = wing.getZP();
    
    cP = wing.getChord();
    tw = wing.getTwist();
    dih = wing.getDih();
    sw = wing.getSweepC4();
    
    # Panel bordes' coordinate and main parameters (at c/4)    
    x = wing.getX();
    y = wing.getY();
    z = wing.getZ();
    c = wing.getChordDist();
    twSec = wing.twSec;
    
    xW = np.unique(np.concatenate([(1.-cf)*0.5*(np.cos(np.linspace(m.pi,0.,mC))+1.),[0.25]]));
    mC = len(xW);
    iC4W = np.where(xW == 0.25)[0][0];                                       # Indice du début du flaps / fin de main element
    zW = np.zeros([mC,len(wing.getAF())],dtype = float);
    for ii in range(len(wing.getAF())):
        zW[:,ii]= camber(wing.getAF(ii),xW);

    xF = 1. - cf * np.linspace(1.,0,mC);
    zF = np.zeros([mC,len(wing.getAF())],dtype = float);
    for ii in range(len(wing.getAF())):
        zF[:,ii] = camber(wing.getAF(ii),xF);
    
    if htail.bool:
        # Panel bordes' coordinate and main parameters (at c/4)
        x = np.concatenate([x,htail.getX()]);
        y = np.concatenate([y,htail.getY()]);
        z = np.concatenate([z,htail.getZ()]);
        c = np.concatenate([c,htail.getChordDist()]);
        twSec = np.concatenate([wing.twSec,htail.twSec]);
        
        # Panels' coordinates and main parameters (at c/4)
        xp = np.concatenate([xp,htail.getXP()]);
        yp = np.concatenate([yp,htail.getYP()]);
        zp = np.concatenate([zp,htail.getZP()]);
        
        cP = np.concatenate([cP,htail.getChord()]);
        tw = np.concatenate([tw,htail.getTwist()]);
        dih = np.concatenate([dih,htail.getDih()]);
        sw = np.concatenate([sw,htail.getSweepC4()]);
        
        # Elevator, Assumed to be as plain flaps
        cfT = htail.getCF();
        
        if cfT != 0: 
            xT = np.unique(np.concatenate([np.linspace(1.,1.-cfT,2),(1.-cfT)*0.5*(np.cos(np.linspace(m.pi,0.,mC-1))+1.)]));
            xT[abs((xT-0.25)) == np.min(abs(xT-0.25))] = 0.25;
        else:
            xT = 0.5*(np.cos(np.linspace(m.pi,0.,mC))+1.);
            xT[abs((xT-0.25)) == np.min(abs(xT-0.25))] = 0.25;
        iC4T = np.where(xT == 0.25)[0][0];
        zT = np.zeros([mC,len(htail.getAF())],dtype = float);
        for ii in range(len(htail.getAF())):
            zT[:,ii-1]= camber(htail.getAF(ii),xT);
            

    #generate grid corner coordinates
    # generate collocation points and normal : where tangency condition is
    # satisfied. Distance from bound vortex depends on the sectional lift
    # curve slope : (dist/localChord) = clAlphas/(4*pi)
    
    X = np.zeros(n * (2 * (mC + mW)+1),dtype = float);
    Y = np.zeros(n * (2 * (mC + mW)+1),dtype = float);                                   # initialization
    Z = np.zeros(n * (2 * (mC + mW)+1),dtype = float);
    COLOCX=np.zeros((mC-1)*n);
    COLOCY=np.zeros((mC-1)*n);
    COLOCZ=np.zeros((mC-1)*n);
    normal = np.zeros([3,(mC-1)*n]);
    coef = 0.25+cla*0.25/m.pi;
    ds = np.zeros((mC-1)*n);                                                    # vector of area of any panel
    dS = np.zeros(wing.r+htail.r);                                                           # vector of area of a spanwise section
    
    xvl = np.zeros(mC + mW,dtype = float);
    yvl = np.zeros(mC + mW,dtype = float);
    zvl = np.zeros(mC + mW,dtype = float);
    xvr = np.zeros(mC + mW,dtype = float);
    yvr = np.zeros(mC + mW,dtype = float);
    zvr = np.zeros(mC + mW,dtype = float);
    xvlf = np.zeros(mC + mW,dtype = float);
    yvlf = np.zeros(mC + mW,dtype = float);
    zvlf = np.zeros(mC + mW,dtype = float);
    xvrf = np.zeros(mC + mW,dtype = float);
    yvrf = np.zeros(mC + mW,dtype = float);
    zvrf = np.zeros(mC + mW,dtype = float);
    dzdx = np.zeros(mW-1,dtype = float);       
    for i in range(wing.getR()):
        camb = zW[:,wing.getAFI(i)]
        cambF = zF[:,wing.getAFI(i)];
        il = i;
        cl = c[il];
        twl = twSec[il];
        
        xl = (xW - 0.25) * cl + x[il];
        yl = y[il] * np.ones(mC);
        zl = camb * cl + z[il];
        
        xlf = (xF - 0.25) * cl + x[il];
        ylf = y[il] * np.ones(mC);
        zlf = cambF * cl + z[il];
        center = np.array([xl[iC4W],yl[iC4W],zl[iC4W]]);
        alpha = 180./m.pi*twl;
        Rot = u.roty(alpha);
        for ii in range(mC):
            point = np.array([xl[ii],yl[ii],zl[ii]])-center;
            point = np.dot(Rot,point) + center;
            xl[ii] = point[0];
            yl[ii] = point[1];
            zl[ii] = point[2];
            pointf = np.array([xlf[ii],ylf[ii],zlf[ii]])-center;
            pointf = np.dot(Rot,pointf) + center;
            xlf[ii] = pointf[0] - 0.02 * cl;
            ylf[ii] = pointf[1];
            zlf[ii] = pointf[2] - 0.02 * cl;

        centerf = np.array([xlf[0],ylf[0],zlf[0]]);
        delta = wing.getDF(i);
        Rotf = u.roty(delta);
        for ii in range(mC):
            pointf = np.array([xlf[ii],ylf[ii],zlf[ii]])-centerf;
            pointf = np.dot(Rotf,pointf) + centerf;
            xlf[ii] = pointf[0];
            ylf[ii] = pointf[1];
            zlf[ii] = pointf[2];  
        xvl[:mC-1] = 0.75 * xl[:-1] + 0.25 * xl[1:];
        yvl[:mC-1] = 0.75 * yl[:-1] + 0.25 * yl[1:];
        zvl[:mC-1] = 0.75 * zl[:-1] + 0.25 * zl[1:];
        xvl[mC-1] = xvl[mC-2] + (xl[-1]-xl[-2]);
        yvl[mC-1:] = yvl[mC-2] + (yl[-1]-yl[-2]);
        zvl[mC-1:] = zvl[mC-2] + (zl[-1]-zl[-2]);
        
        xvlf[:mC-1] = 0.75 * xlf[:-1] + 0.25 * xlf[1:];
        yvlf[:mC-1] = 0.75 * ylf[:-1] + 0.25 * ylf[1:];
        zvlf[:mC-1] = 0.75 * zlf[:-1] + 0.25 * zlf[1:];
        xvlf[mC-1] = xvlf[mC-2] + (xlf[-1]-xlf[-2]);
        yvlf[mC-1:] = yvlf[mC-2] + (ylf[-1]-ylf[-2]);
        zvlf[mC-1:] = zvlf[mC-2] + (zlf[-1]-zlf[-2]);
        
        
        
        centerPropY = prop.getYp() + (xvl[mC-1] - prop.getXp()) * m.tan(beta);
        centerPropZ = prop.getZp() + (xvl[mC-1] - prop.getXp()) * m.tan(aoa);
        vix = 0.;
        for j in range(nbE):
            d = m.sqrt((yvl[mC-1] - centerPropY[j])**2 + (zvl[mC-1] - centerPropZ[j])**2);
            rP = prop.rHub[j];
            D = prop.D[j];
            vitheta = 0.;
            theta0 = np.arctan2(zvl[mC-1] - centerPropZ[j],yvl[mC-1] - centerPropY[j]);
            if ((d >= rP) and (d <= D * 0.5) and prop.Omega[j] != 0.):
                vix += 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                vix2 = 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                a = vix2/V0;
                aprim = 0.5 * (1. - m.sqrt(abs(1.-4.*a*(1.+a)*(V0/(prop.Omega[j] * d))**2)));
                vitheta = OWU[j]*abs((aprim * 2. * prop.Omega[j] * d));
                Theta = times*vitheta/d + theta0;
                dY = np.cos(Theta[1:]) * d + centerPropY[j] - yvl[mC-1];
                dZ = np.sin(Theta[1:]) * d + centerPropZ[j] - zvl[mC-1] ;
                yvl[mC:-1] += dY;
                zvl[mC:-1] += dZ;
        xvl[mC-1:-1] = xvl[mC-1] + times * (V0+vix);
        xvl[-1] = 10. * wing.b;
        
        
        indiceFinLocalEffectCamber = np.where(xvl >= xvl[mC-1] + 2.5 * cl)[0][1]; 
        
        
        dzdxl = (zl[mC-1]-zl[mC-2])/(xl[mC-1]-xl[mC-2]);
        
        # Vérifie ça!
        dydx = V0/(V0+vix) *m.tan(beta) * (1.-np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])));
        dzdx = V0/(V0+vix) * dzdxl * np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])) \
            + V0/(V0+vix) * m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])));
        dY = np.zeros(mW+1);
        dZ = np.zeros(mW+1);
        for ii in range(1,indiceFinLocalEffectCamber-mC+1):
            dZ[ii] = dZ[(ii-1)] + dzdx[ii-1] * (xvl[mC-1+ii] - xvl[(mC-1+ii-1)]);
            dY[ii] = dY[(ii-1)] + dydx[ii-1] * (xvl[mC-1+ii] - xvl[(mC-1+ii-1)]);
        dZ[indiceFinLocalEffectCamber-mC+1:] = dZ[indiceFinLocalEffectCamber-mC] + m.tan(aoa) * (xvl[indiceFinLocalEffectCamber:] - xvl[indiceFinLocalEffectCamber-1]);
        dY[indiceFinLocalEffectCamber-mC+1:] = dY[indiceFinLocalEffectCamber-mC] + m.tan(beta) * (xvl[indiceFinLocalEffectCamber:] - xvl[indiceFinLocalEffectCamber-1]);
        yvl[mC-1:] += dY;
        zvl[mC-1:] += dZ;
        
        
        centerPropY = prop.getYp() + (xvlf[mC-1] - prop.getXp()) * m.tan(beta);
        centerPropZ = prop.getZp() + (xvlf[mC-1] - prop.getXp()) * m.tan(aoa);
        vix = 0.;
        for j in range(nbE):
            d = m.sqrt((yvlf[mC-1] - centerPropY[j])**2 + (zvlf[mC-1] - centerPropZ[j])**2);
            rP = prop.rHub[j];
            D = prop.D[j];
            vitheta = 0.;
            theta0 = np.arctan2(zvlf[mC-1] - centerPropZ[j],yvlf[mC-1] - centerPropY[j]);
            if ((d >= rP) and (d <= D * 0.5) and prop.Omega[j] != 0.):
                vix += 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                vix2 = 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                a = vix2/V0;
                aprim = 0.5 * (1. - m.sqrt(abs(1.-4.*a*(1.+a)*(V0/(prop.Omega[j] * d))**2)));
                vitheta = OWU[j]*abs((aprim * 2. * prop.Omega[j] * d));
                Theta = times*vitheta/d + theta0;
                dY = np.cos(Theta[1:]) * d + centerPropY[j] - yvlf[mC-1];
                dZ = np.sin(Theta[1:]) * d + centerPropZ[j] - zvlf[mC-1] ;
                yvlf[mC:-1] += dY;
                zvlf[mC:-1] += dZ;
        xvlf[mC-1:-1] = xvlf[mC-1] + times * (V0+vix);
        xvlf[-1] = 10. * wing.b;
        
        
        indiceFinLocalEffectCamber = np.where(xvlf >= xvlf[mC-1] + 2.5 * cl)[0][1]; 
        
        
        dzdxl = (zlf[mC-1]-zlf[mC-2])/(xlf[mC-1]-xlf[mC-2]);
        
        # Vérifie ça!
        dydx = V0/(V0+vix) *m.tan(beta) * (1.-np.exp(-3.*(np.array(xvlf[mC:indiceFinLocalEffectCamber] - xvlf[mC-1]))/(xvlf[indiceFinLocalEffectCamber-1] - xvlf[mC-1])));
        dzdx = V0/(V0+vix) * dzdxl * np.exp(-3.*(np.array(xvlf[mC:indiceFinLocalEffectCamber] - xvlf[mC-1]))/(xvlf[indiceFinLocalEffectCamber-1] - xvlf[mC-1])) \
            + V0/(V0+vix) * m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvlf[mC:indiceFinLocalEffectCamber] - xvlf[mC-1]))/(xvlf[indiceFinLocalEffectCamber-1] - xvlf[mC-1])));
        dY = np.zeros(mW+1);
        dZ = np.zeros(mW+1);
        for ii in range(1,indiceFinLocalEffectCamber-mC+1):
            dZ[ii] = dZ[(ii-1)] + dzdx[ii-1] * (xvlf[mC-1+ii] - xvlf[(mC-1+ii-1)]);
            dY[ii] = dY[(ii-1)] + dydx[ii-1] * (xvlf[mC-1+ii] - xvlf[(mC-1+ii-1)]);
        dZ[indiceFinLocalEffectCamber-mC+1:] = dZ[indiceFinLocalEffectCamber-mC] + m.tan(aoa) * (xvlf[indiceFinLocalEffectCamber:] - xvlf[indiceFinLocalEffectCamber-1]);
        dY[indiceFinLocalEffectCamber-mC+1:] = dY[indiceFinLocalEffectCamber-mC] + m.tan(beta) * (xvlf[indiceFinLocalEffectCamber:] - xvlf[indiceFinLocalEffectCamber-1]);
        yvlf[mC-1:] += dY;
        zvlf[mC-1:] += dZ;
        
        
        
        
        ## Right Part
        
        ir = i+1;
        cr = c[ir];
        twr = twSec[ir];
        
        xr = (xW - 0.25) * cr + x[ir];
        yr = y[ir] * np.ones(mC);
        zr = camb * cr + z[ir];
        
        xrf = (xF - 0.25) * cr + x[ir];
        yrf = y[ir] * np.ones(mC);
        zrf = cambF * cr + z[ir];
        
        center = np.array([xr[iC4W],yr[iC4W],zr[iC4W]]);
        alpha = 180./m.pi*twr;
        Rot = u.roty(alpha);
        for ii in range(0,mC):
            point = np.array([xr[ii],yr[ii],zr[ii]])-center;
            point = np.dot(Rot,point) + center;
            xr[ii] = point[0];
            yr[ii] = point[1];
            zr[ii] = point[2];
            pointf = np.array([xrf[ii],yrf[ii],zrf[ii]])-center;
            pointf = np.dot(Rot,pointf) + center;
            xrf[ii] = pointf[0] - 0.02 * cr;
            yrf[ii] = pointf[1];
            zrf[ii] = pointf[2] - 0.02 * cr;
        
        centerf = np.array([xrf[0],yrf[0],zrf[0]]);
        for ii in range(mC):
            pointf = np.array([xrf[ii],yrf[ii],zrf[ii]])-centerf;
            pointf = np.dot(Rotf,pointf) + centerf;
            xrf[ii] = pointf[0];
            yrf[ii] = pointf[1];
            zrf[ii] = pointf[2];  
        
        xvr[:mC-1] = 0.75 * xr[:-1] + 0.25 * xr[1:];
        yvr[:mC-1] = 0.75 * yr[:-1] + 0.25 * yr[1:];
        zvr[:mC-1] = 0.75 * zr[:-1] + 0.25 * zr[1:];
        xvr[mC-1] = xvr[mC-2] + (xr[-1]-xr[-2]);
        yvr[mC-1:] = yvr[mC-2] + (yr[-1]-yr[-2]);
        zvr[mC-1:] = zvr[mC-2] + (zr[-1]-zr[-2]);
        
        xvrf[:mC-1] = 0.75 * xrf[:-1] + 0.25 * xrf[1:];
        yvrf[:mC-1] = 0.75 * yrf[:-1] + 0.25 * yrf[1:];
        zvrf[:mC-1] = 0.75 * zrf[:-1] + 0.25 * zrf[1:];
        xvrf[mC-1] = xvrf[mC-2] + (xrf[-1]-xrf[-2]);
        yvrf[mC-1:] = yvrf[mC-2] + (yrf[-1]-yrf[-2]);
        zvrf[mC-1:] = zvrf[mC-2] + (zrf[-1]-zrf[-2]);
        
        centerPropY = prop.getYp() + (xvr[mC-1] - prop.getXp()) * m.tan(beta);
        centerPropZ = prop.getZp() + (xvr[mC-1] - prop.getXp()) * m.tan(aoa);
        vix = 0.;
        for j in range(nbE):
            d = m.sqrt((yvr[mC-1] - centerPropY[j])**2 + (zvr[mC-1] - centerPropZ[j])**2);
            rP = prop.rHub[j];
            D = prop.D[j];
            vitheta = 0.;
            theta0 = np.arctan2(zvr[mC-1] - centerPropZ[j],yvr[mC-1] - centerPropY[j]);
            if ((d >= rP) and (d <= D * 0.5) and prop.Omega[j] != 0.):
                vix += 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                vix2 = 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                a = vix2/V0;
                aprim = 0.5 * (1. - m.sqrt(abs(1.-4.*a*(1.+a)*(V0/(prop.Omega[j] * d))**2)));
                vitheta = OWU[j]*abs((aprim * 2. * prop.Omega[j] * d));
                Theta = times*vitheta/d + theta0;
                dY = np.cos(Theta[1:]) * d + centerPropY[j] - yvr[mC-1];
                dZ = np.sin(Theta[1:]) * d + centerPropZ[j] - zvr[mC-1] ;
                yvr[mC:-1] += dY;
                zvr[mC:-1] += dZ;
        xvr[mC-1:-1] = xvr[mC-1] + times * (V0+vix);
        xvr[-1] = 10. * wing.b;
        
        
        indiceFinLocalEffectCamber = np.where(xvr >= xvr[mC-1] + 2.5 * cr)[0][1]; 
        
        
        dzdxr = (zr[mC-1]-zr[mC-2])/(xr[mC-1]-xr[mC-2]);
        
        # Vérifie ça!
        dydx = V0/(V0+vix) *m.tan(beta) * (1.-np.exp(-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1])));
        dzdx = V0/(V0+vix) * dzdxr * np.exp(-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1])) \
            + V0/(V0+vix) * m.tan(aoa) * (1.-np.exp((-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1]))));
        dY = np.zeros(mW+1);
        dZ = np.zeros(mW+1);
        for ii in range(1,indiceFinLocalEffectCamber-mC+1):
            dZ[ii] = dZ[(ii-1)] + dzdx[ii-1] * (xvr[mC-1+ii] - xvr[(mC-1+ii-1)]);
            dY[ii] = dY[(ii-1)] + dydx[ii-1] * (xvr[mC-1+ii] - xvr[(mC-1+ii-1)]);
        dZ[indiceFinLocalEffectCamber-mC+1:] = dZ[indiceFinLocalEffectCamber-mC] + m.tan(aoa) * (xvr[indiceFinLocalEffectCamber:] - xvr[indiceFinLocalEffectCamber-1]);
        dY[indiceFinLocalEffectCamber-mC+1:] = dY[indiceFinLocalEffectCamber-mC] + m.tan(beta) * (xvr[indiceFinLocalEffectCamber:] - xvr[indiceFinLocalEffectCamber-1]);
        yvr[mC-1:] += dY;
        zvr[mC-1:] += dZ;
        
        centerPropY = prop.getYp() + (xvrf[mC-1] - prop.getXp()) * m.tan(beta);
        centerPropZ = prop.getZp() + (xvrf[mC-1] - prop.getXp()) * m.tan(aoa);
        vix = 0.;
        for j in range(nbE):
            d = m.sqrt((yvrf[mC-1] - centerPropY[j])**2 + (zvrf[mC-1] - centerPropZ[j])**2);
            rP = prop.rHub[j];
            D = prop.D[j];
            vitheta = 0.;
            theta0 = np.arctan2(zvrf[mC-1] - centerPropZ[j],yvrf[mC-1] - centerPropY[j]);
            if ((d >= rP) and (d <= D * 0.5) and prop.Omega[j] != 0.):
                vix += 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                vix2 = 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                a = vix2/V0;
                aprim = 0.5 * (1. - m.sqrt(abs(1.-4.*a*(1.+a)*(V0/(prop.Omega[j] * d))**2)));
                vitheta = OWU[j]*abs((aprim * 2. * prop.Omega[j] * d));
                Theta = times*vitheta/d + theta0;
                dY = np.cos(Theta[1:]) * d + centerPropY[j] - yvrf[mC-1];
                dZ = np.sin(Theta[1:]) * d + centerPropZ[j] - zvrf[mC-1] ;
                yvrf[mC:-1] += dY;
                zvrf[mC:-1] += dZ;
        xvrf[mC-1:-1] = xvrf[mC-1] + times * (V0+vix);
        xvrf[-1] = 10. * wing.b;
        
        
        indiceFinLocalEffectCamber = np.where(xvrf >= xvrf[mC-1] + 2.5 * cr)[0][1]; 
        
        
        dzdxr = (zrf[mC-1]-zrf[mC-2])/(xrf[mC-1]-xrf[mC-2]);
        
        # Vérifie ça!
        dydx = V0/(V0+vix) *m.tan(beta) * (1.-np.exp(-3.*(np.array(xvrf[mC:indiceFinLocalEffectCamber] - xvrf[mC-1]))/(xvrf[indiceFinLocalEffectCamber-1] - xvrf[mC-1])));
        dzdx = V0/(V0+vix) * dzdxl * np.exp(-3.*(np.array(xvrf[mC:indiceFinLocalEffectCamber] - xvrf[mC-1]))/(xvrf[indiceFinLocalEffectCamber-1] - xvrf[mC-1])) \
            + V0/(V0+vix) * m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvrf[mC:indiceFinLocalEffectCamber] - xvrf[mC-1]))/(xvrf[indiceFinLocalEffectCamber-1] - xvrf[mC-1])));
        dY = np.zeros(mW+1);
        dZ = np.zeros(mW+1);
        for ii in range(1,indiceFinLocalEffectCamber-mC+1):
            dZ[ii] = dZ[(ii-1)] + dzdx[ii-1] * (xvrf[mC-1+ii] - xvrf[(mC-1+ii-1)]);
            dY[ii] = dY[(ii-1)] + dydx[ii-1] * (xvrf[mC-1+ii] - xvrf[(mC-1+ii-1)]);
        dZ[indiceFinLocalEffectCamber-mC+1:] = dZ[indiceFinLocalEffectCamber-mC] + m.tan(aoa) * (xvrf[indiceFinLocalEffectCamber:] - xvrf[indiceFinLocalEffectCamber-1]);
        dY[indiceFinLocalEffectCamber-mC+1:] = dY[indiceFinLocalEffectCamber-mC] + m.tan(beta) * (xvrf[indiceFinLocalEffectCamber:] - xvrf[indiceFinLocalEffectCamber-1]);
        yvrf[mC-1:] += dY;
        zvrf[mC-1:] += dZ;
        
        
        setTable(X,2*(mC+mW)+1,i,np.concatenate([[xvl[0]],xvr,xvl[::-1]]));
        setTable(Y,2*(mC+mW)+1,i,np.concatenate([[yvl[0]],yvr,yvl[::-1]]));
        setTable(Z,2*(mC+mW)+1,i,np.concatenate([[zvl[0]],zvr,zvl[::-1]]));
        
        setTable(X,2*(mC+mW)+1,wing.r+i,np.concatenate([[xvlf[0]],xvrf,xvlf[::-1]]));
        setTable(Y,2*(mC+mW)+1,wing.r+i,np.concatenate([[yvlf[0]],yvrf,yvlf[::-1]]));
        setTable(Z,2*(mC+mW)+1,wing.r+i,np.concatenate([[zvlf[0]],zvrf,zvlf[::-1]]));
        
        
        for j in range(mC-1):
            
            val = [xvl[j],xvr[j],0.5* (xl[j] + xr[j]), 0.5* (xl[j+1] + xr[j+1])];
            COLOCX[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpx1 = val[1] - val[0];
            cpx2 = val[3] - val[2];
            
            val = [yvl[j],yvr[j],0.5* (yl[j] + yr[j]), 0.5* (yl[j+1] + yr[j+1])];
            COLOCY[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpy1 = val[1] - val[0];
            cpy2 = val[3] - val[2];
            
            val = [zvl[j],zvr[j],0.5* (zl[j] + zr[j]), 0.5* (zl[j+1] + zr[j+1])];
            COLOCZ[i * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpz1 = val[1] - val[0];
            cpz2 = val[3] - val[2];
            
            cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
            cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
            ds[i * (mC-1) + j] = cpmag;
            normal[:, i * (mC-1) + j] = cp/cpmag;
            
            val = [xvlf[j],xvrf[j],0.5* (xlf[j] + xrf[j]), 0.5* (xlf[j+1] + xrf[j+1])];
            COLOCX[(i+wing.r) * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpx1 = val[1] - val[0];
            cpx2 = val[3] - val[2];
            
            val = [yvlf[j],yvrf[j],0.5* (ylf[j] + yrf[j]), 0.5* (ylf[j+1] + yrf[j+1])];
            COLOCY[(i+wing.r) * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpy1 = val[1] - val[0];
            cpy2 = val[3] - val[2];
            
            val = [zvlf[j],zvrf[j],0.5* (zlf[j] + zrf[j]), 0.5* (zlf[j+1] + zrf[j+1])];
            COLOCZ[(i+wing.r) * (mC-1) + j] = val[2] * (1.-coef[i]) + val[3] * coef[i];
            cpz1 = val[1] - val[0];
            cpz2 = val[3] - val[2];
            
            cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
            cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
            ds[(i + wing.r) * (mC-1) + j] = cpmag;
            normal[:, (i + wing.r) * (mC-1) + j] = cp/cpmag;
            
        dS[i] = sum(ds[i * (mC-1):(i+1) * (mC-1)]) + sum(ds[(i+wing.r) * (mC-1):(i+wing.r+1) * (mC-1)]);
    for i in range(2*wing.getR(),2*wing.getR()+htail.getR()):
        iPT = i- 2 * wing.getR();
        camb = zT[:,htail.getAFI(iPT)]
        
        il = i+1 - wing.r;
        cl = c[il];
        twl = twSec[il];
        
        xl = (xT - 0.25) * cl + x[il];
        yl = y[il] * np.ones(mC);
        zl = camb * cl + z[il];
        center = np.array([xl[iC4T],yl[iC4T],zl[iC4T]]);
        alpha = 180./m.pi*twl;
        Rot = u.roty(alpha);
        for ii in range(mC):
            point = np.array([xl[ii],yl[ii],zl[ii]])-center;
            point = np.dot(Rot,point) + center;
            xl[ii] = point[0];
            yl[ii] = point[1];
            zl[ii] = point[2];
        if htail.getDF(iPT) != 0.:
            delta = htail.getDF(iPT);
            RotF = u.roty(delta);
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
        yvl[mC-1:] = yvl[mC-2] + (yl[-1]-yl[-2]);
        zvl[mC-1:] = zvl[mC-2] + (zl[-1]-zl[-2]);
        # End of chord vortex = begining of wake vortex
        centerPropY = prop.getYp() + (xvl[mC-1] - prop.getXp()) * m.tan(beta);
        centerPropZ = prop.getZp() + (xvl[mC-1] - prop.getXp()) * m.tan(aoa);
        vix = 0.;
        for j in range(nbE):
            d = m.sqrt((yvl[mC-1] - centerPropY[j])**2 + (zvl[mC-1] - centerPropZ[j])**2);
            rP = prop.rHub[j];
            D = prop.D[j];
            vitheta = 0.;
            theta0 = np.arctan2(zvl[mC-1] - centerPropZ[j],yvl[mC-1] - centerPropY[j]);
            if ((d >= rP) and (d <= D * 0.5) and prop.Omega[j] != 0.):
                vix += 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                vix2 = 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                a = vix2/V0;
                aprim = 0.5 * (1. - m.sqrt(abs(1.-4.*a*(1.+a)*(V0/(prop.Omega[j] * d))**2)));
                vitheta = OWU[j]*abs((aprim * 2. * prop.Omega[j] * d));
                Theta = times*vitheta/d + theta0;
                dY = np.cos(Theta[1:]) * d + centerPropY[j] - yvl[mC-1];
                dZ = np.sin(Theta[1:]) * d + centerPropZ[j] - zvl[mC-1] ;
                yvl[mC:-1] += dY;
                zvl[mC:-1] += dZ;
        xvl[mC-1:-1] = xvl[mC-1] + times * (V0+vix);
        xvl[-1] = 10. * wing.b;
        
        
        indiceFinLocalEffectCamber = np.where(xvl >= xvl[mC-1] + 2.5 * cl)[0][1]; 
        
        
        dzdxl = (zl[mC-1]-zl[mC-2])/(xl[mC-1]-xl[mC-2]);
        
        # Vérifie ça!
        dydx = V0/(V0+vix) *m.tan(beta) * (1.-np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])));
        dzdx = V0/(V0+vix) * dzdxl * np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])) \
            + V0/(V0+vix) * m.tan(aoa) * (1.-np.exp(-3.*(np.array(xvl[mC:indiceFinLocalEffectCamber] - xvl[mC-1]))/(xvl[indiceFinLocalEffectCamber-1] - xvl[mC-1])));
        dY = np.zeros(mW+1);
        dZ = np.zeros(mW+1);
        for ii in range(1,indiceFinLocalEffectCamber-mC+1):
            dZ[ii] = dZ[(ii-1)] + dzdx[ii-1] * (xvl[mC-1+ii] - xvl[(mC-1+ii-1)]);
            dY[ii] = dY[(ii-1)] + dydx[ii-1] * (xvl[mC-1+ii] - xvl[(mC-1+ii-1)]);
        dZ[indiceFinLocalEffectCamber-mC+1:] = dZ[indiceFinLocalEffectCamber-mC] + m.tan(aoa) * (xvl[indiceFinLocalEffectCamber:] - xvl[indiceFinLocalEffectCamber-1]);
        dY[indiceFinLocalEffectCamber-mC+1:] = dY[indiceFinLocalEffectCamber-mC] + m.tan(beta) * (xvl[indiceFinLocalEffectCamber:] - xvl[indiceFinLocalEffectCamber-1]);
        yvl[mC-1:] += dY;
        zvl[mC-1:] += dZ;
        
        ir = i+2 - wing.r;
        cr = c[ir];
        twr = twSec[ir];
        
        xr = (xT - 0.25) * cr + x[ir];
        yr = y[ir] * np.ones(mC);
        zr = camb * cr + z[ir];
        center = np.array([xr[iC4T],yr[iC4T],zr[iC4T]]);
        alpha = 180./m.pi*twr;
        Rot = u.roty(alpha);
        for ii in range(0,mC):
            point = np.array([xr[ii],yr[ii],zr[ii]])-center;
            point = np.dot(Rot,point) + center;
            xr[ii] = point[0];
            yr[ii] = point[1];
            zr[ii] = point[2];
        if htail.getDF(iPT) != 0.:
            delta = htail.getDF(iPT);
            RotF = u.roty(delta);
            center = np.array([xr[-2],yr[-2],zr[-2]]);
            point = np.array([xr[-1],yr[-1],zr[-1]])-center;
            point = np.dot(RotF,point) + center;
            xr[-1] = point[0];
            yr[-1] = point[1];
            zr[-1] = point[2];
            
        xvr[:mC-1] = 0.75 * xr[:-1] + 0.25 * xr[1:];
        yvr[:mC-1] = 0.75 * yr[:-1] + 0.25 * yr[1:];
        zvr[:mC-1] = 0.75 * zr[:-1] + 0.25 * zr[1:];
        xvr[mC-1] = xvr[mC-2] + (xr[-1]-xr[-2]);
        yvr[mC-1:] = yvr[mC-2] + (yr[-1]-yr[-2]);
        zvr[mC-1:] = zvr[mC-2] + (zr[-1]-zr[-2]);
        # End of chord vortex = begining of wake vortex
        centerPropY = prop.getYp() + (xvr[mC-1] - prop.getXp()) * m.tan(beta);
        centerPropZ = prop.getZp() + (xvr[mC-1] - prop.getXp()) * m.tan(aoa);
        vix = 0.;
        for j in range(nbE):
            d = m.sqrt((yvr[mC-1] - centerPropY[j])**2 + (zvr[mC-1] - centerPropZ[j])**2);
            rP = prop.rHub[j];
            D = prop.D[j];
            vitheta = 0.;
            theta0 = np.arctan2(zvr[mC-1] - centerPropZ[j],yvr[mC-1] - centerPropY[j]);
            if ((d >= rP) and (d <= D * 0.5) and prop.Omega[j] != 0.):
                vix += 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                vix2 = 0.5*V0*(m.sqrt(1.+2.*prop.T[j]/(rho*Sh[j]*V0**2))-1.);
                a = vix2/V0;
                aprim = 0.5 * (1. - m.sqrt(abs(1.-4.*a*(1.+a)*(V0/(prop.Omega[j] * d))**2)));
                vitheta = OWU[j]*abs((aprim * 2. * prop.Omega[j] * d));
                Theta = times*vitheta/d + theta0;
                dY = np.cos(Theta[1:]) * d + centerPropY[j] - yvr[mC-1];
                dZ = np.sin(Theta[1:]) * d + centerPropZ[j] - zvr[mC-1] ;
                yvr[mC:-1] += dY;
                zvr[mC:-1] += dZ;
        xvr[mC-1:-1] = xvr[mC-1] + times * (V0+vix);
        xvr[-1] = 10. * wing.b;
        
        
        indiceFinLocalEffectCamber = np.where(xvr >= xvr[mC-1] + 2.5 * cr)[0][1]; 
        
        
        dzdxr = (zr[mC-1]-zr[mC-2])/(xr[mC-1]-xr[mC-2]);
        
        # Vérifie ça!
        dydx = V0/(V0+vix) *m.tan(beta) * (1.-np.exp(-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1])));
        dzdx = V0/(V0+vix) * dzdxr * np.exp(-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1])) \
            + V0/(V0+vix) * m.tan(aoa) * (1.-np.exp((-3.*(np.array(xvr[mC:indiceFinLocalEffectCamber] - xvr[mC-1]))/(xvr[indiceFinLocalEffectCamber-1] - xvr[mC-1]))));
        dY = np.zeros(mW+1);
        dZ = np.zeros(mW+1);
        for ii in range(1,indiceFinLocalEffectCamber-mC+1):
            dZ[ii] = dZ[(ii-1)] + dzdx[ii-1] * (xvr[mC-1+ii] - xvr[(mC-1+ii-1)]);
            dY[ii] = dY[(ii-1)] + dydx[ii-1] * (xvr[mC-1+ii] - xvr[(mC-1+ii-1)]);
        dZ[indiceFinLocalEffectCamber-mC+1:] = dZ[indiceFinLocalEffectCamber-mC] + m.tan(aoa) * (xvr[indiceFinLocalEffectCamber:] - xvr[indiceFinLocalEffectCamber-1]);
        dY[indiceFinLocalEffectCamber-mC+1:] = dY[indiceFinLocalEffectCamber-mC] + m.tan(beta) * (xvr[indiceFinLocalEffectCamber:] - xvr[indiceFinLocalEffectCamber-1]);
        yvr[mC-1:] += dY;
        zvr[mC-1:] += dZ;
        setTable(X,2*(mC+mW)+1,i,np.concatenate([[xvl[0]],xvr,xvl[::-1]]));
        setTable(Y,2*(mC+mW)+1,i,np.concatenate([[yvl[0]],yvr,yvl[::-1]]));
        setTable(Z,2*(mC+mW)+1,i,np.concatenate([[zvl[0]],zvr,zvl[::-1]]));
        
        for j in range(mC-1):
            
            val = [xvl[j],xvr[j],0.5* (xl[j] + xr[j]), 0.5* (xl[j+1] + xr[j+1])];
            COLOCX[i * (mC-1) + j] = val[2] * (1.-coef[i - wing.r]) + val[3] * coef[i - wing.r];
            cpx1 = val[1] - val[0];
            cpx2 = val[3] - val[2];
            
            val = [yvl[j],yvr[j],0.5* (yl[j] + yr[j]), 0.5* (yl[j+1] + yr[j+1])];
            COLOCY[i * (mC-1) + j] = val[2] * (1.-coef[i - wing.r]) + val[3] * coef[i - wing.r];
            cpy1 = val[1] - val[0];
            cpy2 = val[3] - val[2];
            
            val = [zvl[j],zvr[j],0.5* (zl[j] + zr[j]), 0.5* (zl[j+1] + zr[j+1])];
            COLOCZ[i * (mC-1) + j] = val[2] * (1.-coef[i - wing.r]) + val[3] * coef[i - wing.r];
            cpz1 = val[1] - val[0];
            cpz2 = val[3] - val[2];
            
            cp= np.cross(np.array([cpx1,cpy1,cpz1]),np.array([cpx2,cpy2,cpz2]));
            cpmag= m.sqrt(cp[1]*cp[1]+cp[2]*cp[2]+cp[0]*cp[0]);
            ds[i * (mC-1) + j] = cpmag;
            normal[:, i * (mC-1) + j] = cp/cpmag;
        dS[i-wing.r] = sum(ds[i * (mC-1):(i+1) * (mC-1)]);
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
    select = np.zeros([wing.r + htail.r,n * (mC-1)]); # rechercher intensité du dernier vortex uniquement
    select2 = np.zeros([n * (mC-1),wing.r + htail.r]); # pour chaque paneau sur même section y, même velocity triangle
    select3 = np.zeros([wing.r + htail.r + len(ac.prop.D),n * (mC-1) + len(ac.prop.D)]); # 
    for i in range(wing.r):
        select[i,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i] = 1.;
        select3[i,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i];
    for i in range(wing.r,n):
        select[i-wing.r,(mC-2) + (mC-1)*i] = 1.;
        select2[(mC-1)*i:(mC-1)*(i+1),i - wing.r] = 1.;
        select3[i - wing.r,(mC-1)*i:(mC-1)*(i+1)] = ds[(mC-1)*i:(mC-1)*(i+1)]/dS[i-wing.r];
    if ac.prop.bool:
        select3[-len(ac.prop.D):,-len(ac.prop.D):] = np.eye(len(ac.prop.D));
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

def ICM_F(X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac,n,mC,mW):
    if ac.fus.bool:
        HWing = ac.fus.vDist > 0;
    if ac.htail.bool and ac.vtail.bool:
        HTail = ac.htail.z[ac.htail.getR()/2] > ((ac.vtail.z[-1]-ac.vtail.z[0]) * 0.66) + ac.vtail.z[0];
    if not(ac.fus.bool):
        if not(ac.vtail.bool) or not(ac.htail.bool) or HTail:
            A,Vx,Vy,Vz = OnlyWing(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
        else:
            A,Vx,Vy,Vz = BothWingOneTailVtailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
    else:
        if not(ac.vtail.bool):
            if HWing:
                A,Vx,Vy,Vz = OnlyWing(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            else:
                A,Vx,Vy,Vz = OneWingBothTailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
        else:
            if HWing and HTail:
                A,Vx,Vy,Vz = OnlyWing(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            elif HTail:
                A,Vx,Vy,Vz = OneWingBothTailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            elif HWing:
                A,Vx,Vy,Vz = BothWingOneTailVtailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
            else:
                A,Vx,Vy,Vz = OneWingOneTailVtailF(n,mC,mW,X,Y,Z,COLOCX,COLOCY,COLOCZ,normal,ac);
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
    for b in range(n * (mC - 1)):
        for j in range(n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    if ac.prop.bool:
        for b in range(n * (mC - 1),m):
            x = ac.prop.xp[b-n* (mC - 1)];
            y = ac.prop.yp[b-n* (mC - 1)];
            z = ac.prop.zp[b-n* (mC - 1)];
            for j in range(n):
                pathX = getVal(X,2*(mW+mC)+1,j);
                pathY = getVal(Y,2*(mW+mC)+1,j);
                pathZ = getVal(Z,2*(mW+mC)+1,j);
                a,vix,viy,viz = vortxl(x,y,z,np.array([-1.,0.,0.]),pathX,pathY,pathZ,mC,mW);
                Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
                Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
                Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
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
    for b in range(ac.wing.getR()*(mC-1)):
        for j in range(n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(ac.wing.getR()*(mC-1),(ac.htail.getR()/2+ac.wing.getR())*(mC-1)):
        for j in range(ac.htail.getR()/2+ac.wing.getR()-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range((ac.wing.getR()+ac.htail.getR()/2)*(mC-1),n*(mC-1)):
        for j in range(ac.wing.getR()):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j = ac.wing.getR()+ac.htail.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()+ac.htail.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    if ac.prop.bool:
        for b in range(n * (mC - 1),m):
            x = ac.prop.xp[b-n* (mC - 1)];
            y = ac.prop.yp[b-n* (mC - 1)];
            z = ac.prop.zp[b-n* (mC - 1)];
            for j in range(n):
                pathX = getVal(X,2*(mW+mC)+1,j);
                pathY = getVal(Y,2*(mW+mC)+1,j);
                pathZ = getVal(Z,2*(mW+mC)+1,j);
                a,vix,viy,viz = vortxl(x,y,z,np.array([-1.,0.,0.]),pathX,pathY,pathZ,mC,mW);
                Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
                Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
                Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
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
    for b in range(2*ac.wing.getR()*(mC-1)):
        for j in range(n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(2*ac.wing.getR()*(mC-1),(ac.htail.getR()/2+2*ac.wing.getR())*(mC-1)):
        for j in range(ac.htail.getR()/2+2*ac.wing.getR()-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range((2*ac.wing.getR()+ac.htail.getR()/2)*(mC-1),n*(mC-1)):
        for j in range(2*ac.wing.getR()):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j = 2*ac.wing.getR()+ac.htail.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(2*ac.wing.getR()+ac.htail.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    if ac.prop.bool:
        for b in range(n * (mC - 1),m):
            x = ac.prop.xp[b-n* (mC - 1)];
            y = ac.prop.yp[b-n* (mC - 1)];
            z = ac.prop.zp[b-n* (mC - 1)];
            for j in range(n):
                pathX = getVal(X,2*(mW+mC)+1,j);
                pathY = getVal(Y,2*(mW+mC)+1,j);
                pathZ = getVal(Z,2*(mW+mC)+1,j);
                a,vix,viy,viz = vortxl(x,y,z,np.array([-1.,0.,0.]),pathX,pathY,pathZ,mC,mW);
                Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
                Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
                Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
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
    for b in range(ac.wing.getR()/2*(mC-1)):
        for j in range(ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR(),n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(ac.wing.getR()/2*(mC-1),ac.wing.getR()*(mC-1)):
        j = ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(ac.wing.getR()*(mC-1),n*(mC-1)):
        for j in range(n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    if ac.prop.bool:
        for b in range(n * (mC - 1),m):
            x = ac.prop.xp[b-n* (mC - 1)];
            y = ac.prop.yp[b-n* (mC - 1)];
            z = ac.prop.zp[b-n* (mC - 1)];
            for j in range(n):
                pathX = getVal(X,2*(mW+mC)+1,j);
                pathY = getVal(Y,2*(mW+mC)+1,j);
                pathZ = getVal(Z,2*(mW+mC)+1,j);
                a,vix,viy,viz = vortxl(x,y,z,np.array([-1.,0.,0.]),pathX,pathY,pathZ,mC,mW);
                Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
                Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
                Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
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
    for b in range(ac.wing.getR()/2*(mC-1)):
        for j in range(ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR(),ac.wing.getR()+ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(2*ac.wing.getR(),n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(ac.wing.getR()/2*(mC-1),ac.wing.getR()*(mC-1)):
        j = ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()/2+1,ac.wing.getR()):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j = ac.wing.r+ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()+ac.wing.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(ac.wing.getR()*(mC-1),ac.wing.getR()*(mC-1)+ac.wing.getR()/2*(mC-1)):
        for j in range(ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR(),ac.wing.getR()+ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(2*ac.wing.getR(),n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(ac.wing.getR()*(mC-1)+ac.wing.getR()/2*(mC-1),2*ac.wing.getR()*(mC-1)):
        j = ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()/2+1,ac.wing.getR()):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j = ac.wing.r+ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()+ac.wing.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(2*ac.wing.getR()*(mC-1),n*(mC-1)):
        for j in range(n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    if ac.prop.bool:
        for b in range(n * (mC - 1),m):
            x = ac.prop.xp[b-n* (mC - 1)];
            y = ac.prop.yp[b-n* (mC - 1)];
            z = ac.prop.zp[b-n* (mC - 1)];
            for j in range(n):
                pathX = getVal(X,2*(mW+mC)+1,j);
                pathY = getVal(Y,2*(mW+mC)+1,j);
                pathZ = getVal(Z,2*(mW+mC)+1,j);
                a,vix,viy,viz = vortxl(x,y,z,np.array([-1.,0.,0.]),pathX,pathY,pathZ,mC,mW);
                Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
                Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
                Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
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
    for b in range(ac.wing.getR()/2*(mC-1)):
        for j in range(ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR(),n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(ac.wing.getR()/2*(mC-1),ac.wing.getR()*(mC-1)):
        j = ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(ac.wing.getR()*(mC-1),(ac.htail.getR()/2+ac.wing.getR())*(mC-1)):
        for j in range(ac.htail.getR()/2+ac.wing.getR()-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range((ac.wing.getR()+ac.htail.getR()/2)*(mC-1),n*(mC-1)):
        for j in range(ac.wing.getR()):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j = ac.wing.getR()+ac.htail.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()+ac.htail.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    if ac.prop.bool:
        for b in range(n * (mC - 1),m):
            x = ac.prop.xp[b-n* (mC - 1)];
            y = ac.prop.yp[b-n* (mC - 1)];
            z = ac.prop.zp[b-n* (mC - 1)];
            for j in range(n):
                pathX = getVal(X,2*(mW+mC)+1,j);
                pathY = getVal(Y,2*(mW+mC)+1,j);
                pathZ = getVal(Z,2*(mW+mC)+1,j);
                a,vix,viy,viz = vortxl(x,y,z,np.array([-1.,0.,0.]),pathX,pathY,pathZ,mC,mW);
                Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
                Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
                Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
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
    for b in range(ac.wing.getR()/2*(mC-1)):
        for j in range(ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR(),ac.wing.getR()+ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(2*ac.wing.getR(),n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(ac.wing.getR()/2*(mC-1),ac.wing.getR()*(mC-1)):
        j = ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()/2+1,ac.wing.getR()):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j = ac.wing.r+ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()+ac.wing.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(ac.wing.getR()*(mC-1),ac.wing.getR()*(mC-1)+ac.wing.getR()/2*(mC-1)):
        for j in range(ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR(),ac.wing.getR()+ac.wing.getR()/2-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(2*ac.wing.getR(),n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(ac.wing.getR()*(mC-1)+ac.wing.getR()/2*(mC-1),2*ac.wing.getR()*(mC-1)):
        j = ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()/2+1,ac.wing.getR()):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j = ac.wing.r+ac.wing.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()+ac.wing.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range(2*ac.wing.getR()*(mC-1),(ac.htail.getR()/2+2*ac.wing.getR())*(mC-1)):
        for j in range(ac.htail.getR()/2+2*ac.wing.getR()-1):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j += 1;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NR(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range((2*ac.wing.getR()+ac.htail.getR()/2)*(mC-1),n*(mC-1)):
        for j in range(2*ac.wing.getR()):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j = 2*ac.wing.getR()+ac.htail.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(2*ac.wing.getR()+ac.htail.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    for b in range((ac.wing.getR()+ac.htail.getR()/2)*(mC-1),n*(mC-1)):
        for j in range(ac.wing.getR()):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        j = ac.wing.getR()+ac.htail.getR()/2;
        pathX = getVal(X,2*(mW+mC)+1,j);
        pathY = getVal(Y,2*(mW+mC)+1,j);
        pathZ = getVal(Z,2*(mW+mC)+1,j);
        a,vix,viy,viz = vortxl_NL(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
        A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
        Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
        Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
        Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
        for j in range(ac.wing.getR()+ac.htail.getR()/2+1,n):
            pathX = getVal(X,2*(mW+mC)+1,j);
            pathY = getVal(Y,2*(mW+mC)+1,j);
            pathZ = getVal(Z,2*(mW+mC)+1,j);
            a,vix,viy,viz = vortxl(COLOCX[b],COLOCY[b],COLOCZ[b],normal[:,b],pathX,pathY,pathZ,mC,mW);
            A[b,j*(mC-1) : (j+1) *(mC-1)] = a;
            Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
            Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
            Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
    if ac.prop.bool:
        for b in range(n * (mC - 1),m):
            x = ac.prop.xp[b-n* (mC - 1)];
            y = ac.prop.yp[b-n* (mC - 1)];
            z = ac.prop.zp[b-n* (mC - 1)];
            for j in range(n):
                pathX = getVal(X,2*(mW+mC)+1,j);
                pathY = getVal(Y,2*(mW+mC)+1,j);
                pathZ = getVal(Z,2*(mW+mC)+1,j);
                a,vix,viy,viz = vortxl(x,y,z,np.array([-1.,0.,0.]),pathX,pathY,pathZ,mC,mW);
                Vx[b,j*(mC-1) : (j+1) *(mC-1)] = vix;
                Vy[b,j*(mC-1) : (j+1) *(mC-1)] = viy;
                Vz[b,j*(mC-1) : (j+1) *(mC-1)] = viz;
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
    coeff[mC+mW+1:nbRing] *= 0.3;
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