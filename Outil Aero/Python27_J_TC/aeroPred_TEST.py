# -*- coding: utf-8 -*-
import math as m;
import numpy as np;
import Polar2 as p;
import utilitaire as u
import VLM_MP
import copy
import Propu as propu
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
def main_SEI(ac,flow,lin):
    """function CL,CD,CY,CM,Cl,Cn,CDi,CD0,Alphas = main(ac,flow,lin);
    Prediction of aerodynamic characteristics of the wing
    Autor : Quentin borlon
    Date : 17-10-2017

    Function that predicts the aerodynamic coefficients for a given aircraft.
    Based on the aircraft geometry and the sectional 2D aerodynamic datas.

    INPUT:
    ac : the full aircraft configuration : wing, htail, vtail, fuselage and prop
    flow : flow configuration
    lin : boolean for linear study

    OUTPUT:
        Alphas : a.o.a. where convergence is reached
        CL : lift coefficient for entire ac
        CD : total drag coefficient for the entire ac
        CY : lateral force coefficient for entire ac
        CM : pitching moment coefficient for entire ac
        Cl : rolling moment coefficient for entire ac
        Cn : yawing moment coefficient for entire ac
        CDi : induced drag coefficient for entire ac
        CD0 : "airfoil" drag coefficient for entire ac"""

    Strongly = False
    if Strongly:
        return StronglyCoupled(ac,flow,lin);
    else:
        return LooselyCoupled(ac,flow,lin);

def VLM_Solve(invA,select,select2,Vx,Vy,Vz,velTriX,velTriY,velTriZ,normal):
    norme = np.sqrt(velTriX*velTriX+velTriY*velTriY+velTriZ*velTriZ);
    RHS = -(normal[0,:] * np.dot(select2,velTriX) + normal[1,:] * np.dot(select2,velTriY) + normal[2,:] * np.dot(select2,velTriZ));
    gamma = np.dot(invA,RHS);
    ViX = np.dot(Vx,gamma);
    ViY = np.dot(Vy,gamma);
    ViZ = np.dot(Vz,gamma);
    ccn = 2. * np.dot(select,gamma)/norme;
    return ViX,ViY,ViZ,ccn;

def VLM_SolveV(invA,Vx,Vy,Vz,velTriX,velTriY,velTriZ,normal):
    norme = np.sqrt(velTriX*velTriX+velTriY*velTriY+velTriZ*velTriZ);
    RHS = -(normal[0,:]* velTriX + normal[1,:] * velTriY + normal[2,:] * velTriZ);
    gamma = np.dot(invA,RHS);
    ViX = np.dot(Vx,gamma);
    ViY = np.dot(Vy,gamma);
    ViZ = np.dot(Vz,gamma);
    ccn = 2. * gamma/norme;
    return ViX,ViY,ViZ,ccn;

    
def main_DEI(ac,flow,lin):
    """
    /!\ NOT YET IMPLEMENTED, should make iterations btw propellers and induced
    velociy by wake on them. 
    
    Personnaly, think it is too time consuming for this stage of the design
    process as we don"t know the exact geometry of the prop/airfoils..
    
    function CL,CD,CY,CM,Cl,Cn,CDi,CD0,Alphas = main(ac,flow,lin);
    Prediction of aerodynamic characteristics of the wing
    Autor : Quentin borlon
    Date : 17-10-2017

    Function that predicts the aerodynamic coefficients for a given aircraft.
    Based on the aircraft geometry and the sectional 2D aerodynamic datas.

    INPUT:
    ac : the full aircraft configuration : wing, htail, vtail, fuselage and prop
    flow : flow configuration
    lin : boolean for linear study

    OUTPUT:
        Alphas : a.o.a. where convergence is reached
        CL : lift coefficient for entire ac
        CD : total drag coefficient for the entire ac
        CY : lateral force coefficient for entire ac
        CM : pitching moment coefficient for entire ac
        Cl : rolling moment coefficient for entire ac
        Cn : yawing moment coefficient for entire ac
        CDi : induced drag coefficient for entire ac
        CD0 : "airfoil" drag coefficient for entire ac"""

    return main_SEI(ac,flow,lin);

def LooselyCoupled(ac,flow,lin):
    # AC configuration
    wing = ac.wing;
    htail = ac.htail;
    vtail = ac.vtail;
    fus = ac.fus;
    prop = ac.prop;
    
    # Numerical parameters initialisation
    tol = 1e-3;                                                                 # Tolerance for the a.o.a. correction loop
    nIterMAX = 100;                                                              # Maximum number of iteration in the correction loop before introducing the numerical damping
    tolWarning = 2e-3;                                                          # Tolerance above which a warning is send to the promp because of lack of accuracy
    tolStop = 5e-2;                                                             # Tolerance above which the software is stopped because of lack of accuracy
    
    # Flight conditions
    Alphas = flow.getAlphas() * m.pi/180;
    Beta = m.pi/180.*flow.getBeta();
    V0 = 1.;
    nbA = len(Alphas);
    # Geometric characteristics
    # Wing
    sweepC4 = wing.getSweepC4();
    for i in range(len(sweepC4)/2):
        sweepC4[i] *= -1.;
    dih = wing.getDih();
    for i in range(len(dih)/2):
        dih[i] *= -1.;
    tp = wing.getTwist();
    x = wing.getXP();
    y = wing.getYP();
    z = wing.getZP();
    c = wing.getChord();
    dA = 0.5*(wing.y[1:]-wing.y[:-1])*(wing.chordDistrib[1:]+wing.chordDistrib[:-1]);
    S = wing.getS();
    MAC = wing.getMac();
    [waPanel,wclPanel, wcdPanel, wcmPanel, alphal0PanelVLM, aMax, cl_alphaPanel] = p.panelPolar(wing);
    nbPan = len(y)
    #Fuselage
    yFus = fus.getYF();
    index = np.array(range(nbPan));
    # Aerodynamic contribution of fuselage
    if fus.bool:
        for i in range(nbPan):
            if (y[i] > yFus[0] and y[i] < yFus[1]):
                index[i] = -1;
        index = index[index >= 0];
        L_F= fus.getL();
        D_F= fus.getD();
        M_F= fus.getM();
        N_F = fus.getN();
        Y_F = fus.getY();
        c_r = np.interp(yFus[1],y,c);
        Rec_r = c_r*V0/(0.00001568);
        dih_r = np.interp(yFus[1],y,dih);
        tc = 0.075;
        C_Dint = (c_r**2)/S*(0.1112 - 0.2572*m.sin(dih_r) + 3.44 * tc - 0.02097 * m.log10(Rec_r) + 0.09009 * m.sin(dih_r) * m.sin(dih_r) - 2.549 * tc * m.sin(dih_r) + 0.03010 * m.log10(Rec_r) * m.sin(dih_r) -  0.1462 * tc * m.log10(Rec_r));
    else:
        C_Dint = 0.;
        L_F = np.zeros(nbA,dtype = float);
        Y_F = 0.;
        N_F = 0.;
        D_F = np.zeros(nbA,dtype = float);
        M_F = np.zeros(nbA,dtype = float);
    
    # Engines
    # Induced velocity by the engine
    # tangencial increases the local speed and then the local lift
    # contribution 
    if prop.bool:
        dyp = prop.yp - ac.rp[1];
        dzp = prop.zp - ac.rp[2];
    else:
        dyp = 0;
        dzp = 0;
        prop.Tc = 0;
    
    # Htail
    if htail.bool:
        yTail = htail.getY();
        cTail = htail.getChordDist();
        dihT = htail.getDih();
        index = np.concatenate([index,np.array(range(index[-1]+1,index[-1]+1+htail.getR()))]);
        dA = np.concatenate([dA , 0.5 * (cTail[1:] + cTail[:-1]) * (yTail[1:] - yTail[:-1])]);
        tp = np.concatenate([tp,htail.getTwist()]);
        
        for i in range(htail.getR()/2):
            dihT[i] *= -1.;
        dih = np.concatenate([dih,dihT]);
        sweepC4T = htail.getSweepC4();
        for i in range(htail.getR()/2):
            sweepC4T[i] *= -1.;
        sweepC4 = np.concatenate([sweepC4,sweepC4T]);
        c = np.concatenate([c,htail.getChord()]);
        x = np.concatenate([x,htail.getXP()]);
        y = np.concatenate([y,htail.getYP()]);
        z = np.concatenate([z,htail.getZP()]);
        nbPan = len(y)
        [waPanelT,wclPanelT, wcdPanelT, wcmPanelT, alphal0PanelVLMT, alphaMaxT, cl_alphaPanelT] = p.panelPolar(htail);
        waPanel = np.concatenate([waPanel,waPanelT],0);
        wclPanel = np.concatenate([wclPanel,wclPanelT],0);
        wcdPanel = np.concatenate([wcdPanel,wcdPanelT],0);
        wcmPanel = np.concatenate([wcmPanel,wcmPanelT],0);
        alphal0PanelVLM = np.concatenate([alphal0PanelVLM,alphal0PanelVLMT]);
        aMax = np.concatenate([aMax,alphaMaxT]);
        cl_alphaPanel= np.concatenate([cl_alphaPanel,cl_alphaPanelT]);
        
    # Vtail
    if vtail.bool:
        zTail = vtail.getZ();
        cTail = vtail.getChordDist();
        indexV = np.array(range(index[-1]+1,wing.getR()+htail.getR()+vtail.getR()),dtype = int);
        
        dA = np.concatenate([dA , 0.5 * (cTail[1:] + cTail[:-1]) * (zTail[1:] - zTail[:-1])]);
        
        sweepC4VT = vtail.getSweepC4();
        sweepC4 = np.concatenate([sweepC4,sweepC4VT]);
        c = np.concatenate([c,vtail.getChord()]);
        x = np.concatenate([x,vtail.getXP()]);
        y = np.concatenate([y,vtail.getYP()]);
        z = np.concatenate([z,vtail.getZP()]);
        nbPan = len(y)
        [waPanelVT,wclPanelVT, wcdPanelVT, wcmPanelVT, alphal0PanelVLMVT, alphaMaxVT, cl_alphaPanelVT] = p.vtailPolar(vtail);
        waPanel = np.concatenate([waPanel,waPanelVT],0);
        wclPanel = np.concatenate([wclPanel,wclPanelVT],0);
        wcdPanel = np.concatenate([wcdPanel,wcdPanelVT],0);
        wcmPanel = np.concatenate([wcmPanel,wcmPanelVT],0);
        alphal0PanelVLM = np.concatenate([alphal0PanelVLM,alphal0PanelVLMVT]);
        aMax = np.concatenate([aMax,alphaMaxVT]);
        cl_alphaPanel= np.concatenate([cl_alphaPanel,cl_alphaPanelVT]);
        if fus.bool:
            c_r = vtail.getChord(0);
            Rec_r = c_r*V0/(0.00001568);
            dih_r = 0.;
            C_Dint += 0.5*(c_r**2)/S*(0.1112 - 0.2572*m.sin(dih_r) + 3.44 * tc - 0.02097 * m.log10(Rec_r) + 0.09009 * m.sin(dih_r) * m.sin(dih_r) - 2.549 * tc * m.sin(dih_r) + 0.03010 * m.log10(Rec_r) * m.sin(dih_r) -  0.1462 * tc * m.log10(Rec_r));
    else:
        CYV= np.zeros(nbA,float);
        CMV= np.zeros(nbA,float);
        ClV= np.zeros(nbA,float);
        CnV= np.zeros(nbA,float);
        CDiV= np.zeros(nbA,float);
        CD0V= np.zeros(nbA,float);
        indexV = np.empty(0,dtype = int);

    # Initialisation and beginning of the computation
    deltaAlpha = np.zeros(nbPan);
    al_i = np.zeros([nbPan,nbA],dtype = float);
    CL = np.zeros(nbA,dtype = float);
    CY = np.zeros(nbA,dtype = float);
    CM = np.zeros(nbA,dtype = float);
    Cl = np.zeros(nbA,dtype = float);
    Cn = np.zeros(nbA,dtype = float);
    CDi = np.zeros(nbA,dtype = float);
    CD0 = np.zeros(nbA,dtype = float);

    lonDist = np.sqrt((x-ac.rp[0])**2 + (z-ac.rp[2])**2);
    angle = np.arctan2(z-ac.rp[2],x-ac.rp[0])

    yDist = y - ac.rp[1];
    stop = False;
    nCorrection = 0;
    K = np.ones(nbPan);
    # Construction of the A matrix
    [A,normal,Vx,Vy,Vz,select,select2] = VLM_MP.ICMatrix(ac,cl_alphaPanel,flow);
    artVisc = np.eye(nbPan);
    indexRes = np.ones(nbA,dtype = bool);
    if lin:
        for ii in range(nbA):
            # Construction of velocity triangles
            AOA = Alphas[ii];
            velTriX = V0 * np.cos(AOA) * np.cos(-Beta) + vix;
            velTriY = V0 * np.sin(-Beta) + viy;
            velTriZ = V0 * np.sin(AOA) * np.cos(-Beta) + viz;
            AOA_VLM = np.arctan2(velTriZ,velTriX);
            
            
            VxA = V0 * np.cos(Alphas[ii]) * np.cos(-Beta)+ vix;
            VyA = V0 * np.sin(-Beta) + viy;
            VzA = V0 * np.sin(Alphas[ii]) * np.cos(-Beta)+ viz;
            velTriXLoc,velTriYLoc,velTriZLoc = u.localVelTri(VxA,VyA,VzA,tp,sweepC4,dih);
            
            Vx0 = V0 * np.cos(Alphas[ii]) * np.cos(-Beta);
            Vy0 = V0 * np.sin(-Beta);
            Vz0 = V0 * np.sin(Alphas[ii]) * np.cos(-Beta);
            velTriX0,velTriY0,velTriZ0 = u.localVelTri(Vx0,Vy0,Vz0,tp,sweepC4,dih);
            af0 = np.arctan2(velTriZ0,velTriX0);
            
            [Vix,Viy,Viz,ccl] = VLM_Solve(A,Vx,Vy,select,select2,Vz,velTriX,velTriY,velTriZ,normal);
            clw = ccl / c;
            velTriXVisc,velTriYVisc,velTriZVisc = u.localVelTri(VxA + Vix[:nbPan] * m.cos(Alphas[ii]) - Viz[:nbPan] * m.sin(Alphas[ii]),VyA + Viy[:nbPan],VzA + Viz[:nbPan]* m.cos(Alphas[ii]) + Vix[:nbPan] * m.sin(Alphas[ii]),tp,sweepC4,dih);
            alphaLoc = np.arctan2(velTriZVisc,velTriXVisc);
            [clVisc,cd,cm] = p.AeroCoef(waPanel,wclPanel,wcdPanel,wcmPanel,alphaLoc);
            
            speedFactor = (velTriXVisc*velTriXVisc+velTriZVisc*velTriZVisc)/(V0*V0);
            al_i[:,ii] = af0 - alphaLoc;
            
            dF = clVisc *  dA * speedFactor;
            
            dCDi = np.sin(al_i[index,ii]) * dF[index] * np.cos(sweepC4);
            dCL = dF[index]* np.cos(al_i[index,ii])*np.cos(dih[index])**2 - cd[index] *  dA[index]*speedFactor[index] * np.sin(al_i[index,ii]);
            dCD0 = cd[index] *  dA[index]*speedFactor[index] * np.cos(al_i[index,ii]);
            dCM = - ( dCL * lonDist[index] * np.cos(angle[index]-Alphas[ii])) + \
                (dCDi + dCD0) * lonDist[index] * np.sin(angle[index] - Alphas[ii]) + \
                cm[index] * dA[index] * c[index] * speedFactor[index];
            Cl[ii] = -np.sum(dCL*yDist[index]) / (S*MAC) + ClV[ii];
            Cn[ii] = (np.sum((dCD0+dCDi)*yDist[index]) + N_F) / (S*MAC)\
              - np.sum(dyp*prop.Tc/MAC) + CnV[ii];
            CY[ii] = Y_F/S + CYV[ii];
            CL[ii] = (np.sum(dCL)+L_F[ii])/S + np.sum(prop.Tc) * m.sin(Alphas[ii]) ;
            CDi[ii] = np.sum(dCDi)/S + CDiV[ii];
            CD0[ii] = CD0V[ii] + (np.sum(dCD0)+D_F[ii])/S - np.sum(prop.Tc) * m.cos(Alphas[ii]);
            CM[ii] = CMV[ii] + (np.sum(dCM) + M_F[ii])/(S*MAC) - np.sum(dzp*prop.Tc/MAC);
            al_i[:,ii] = 180./m.pi*(al_i[:,ii]);
        CD = CD0+CDi+C_Dint;
        Alphas = Alphas*180./m.pi;
        return CL,CD,CY,CM,Cl,Cn,CDi,CD0,Alphas
    else:
        for ii in range(nbA):                                              # if possible no numerical damping to fasten the computation
            AOA = Alphas[ii];
            vix, viy, viz = propu.computeVelocity(ac, AOA, Beta);
            velTriX = V0 * np.cos(AOA) * np.cos(-Beta) + vix;
            velTriY = V0 * np.sin(-Beta) + viy;
            velTriZ = V0 * np.sin(AOA) * np.cos(-Beta) + viz;
            AOA_VLM = np.concatenate([np.arctan2(velTriZ[index],velTriX[index]),np.arctan2(velTriY[indexV],velTriX[indexV])]);
            AOA_VLM0 = AOA_VLM;
            AOA_VLM += deltaAlpha;
            Norme = np.concatenate([np.sqrt(velTriZ[index]*velTriZ[index]+velTriX[index]*velTriX[index]),np.sqrt(velTriY[indexV]*velTriY[indexV]+velTriX[indexV]*velTriX[indexV])]);
            velTriX = Norme * np.cos(AOA_VLM);
            velTriZ[index] = Norme[index] * np.sin(AOA_VLM[index]);
            velTriY[indexV] = Norme[indexV] * np.sin(AOA_VLM[indexV]);
            
            
            VxA = V0 * np.cos(Alphas[ii]) * np.cos(-Beta) + vix;
            VyA = V0 * np.sin(-Beta) + viy;
            VzA = V0 * np.sin(Alphas[ii]) * np.cos(-Beta) + viz;
            velTriXLoc,velTriYLoc,velTriZLoc = u.localVelTri(VxA[index],VyA[index],VzA[index],tp,sweepC4[index],dih);
            velTriXLocV,velTriYLocV,velTriZLocV = u.localVelTriVT(VxA[indexV],VyA[indexV],VzA[indexV],sweepC4[indexV]);
            
            
            Vx0 = V0 * np.cos(Alphas[ii]) * np.cos(-Beta);
            Vy0 = V0 * np.sin(-Beta);
            Vz0 = V0 * np.sin(Alphas[ii]) * np.cos(-Beta);
            velTriX0,velTriY0,velTriZ0 = u.localVelTri(Vx0,Vy0,Vz0,tp,sweepC4[index],dih);
            velTriX0V,velTriY0V,velTriZ0V = u.localVelTriVT(Vx0,Vy0,Vz0,sweepC4[indexV]);
            AOA_0 = np.arctan2(velTriZ0,velTriX0);
            AOA_0V = np.arctan2(velTriY0V,velTriX0V);
            AOA_0 = np.concatenate((AOA_0,AOA_0V));
            
            # Linear prediction of cl and induced velocities
            [Vix,Viy,Viz,ccl] = VLM_Solve(A,select,select2,Vx,Vy,Vz,velTriX,velTriY,velTriZ,normal);
            clw = ccl / c;
            clw[indexV] *= -1;
            # Viscous prediction based on prediction of induced velocity
            velTriXVisc,velTriYVisc,velTriZVisc = u.localVelTri(VxA[index] + Vix[index] * m.cos(Alphas[ii]) - Viz[index] * m.sin(Alphas[ii]),VyA[index] + Viy[index],VzA[index] + Viz[index] * m.cos(Alphas[ii]) + Vix[index] * m.sin(Alphas[ii]),tp,sweepC4[index],dih);
            velTriXViscV,velTriYViscV,velTriZViscV = u.localVelTriVT(VxA[indexV]+ Vix[indexV] * m.cos(Alphas[ii]) - Viz[indexV] * m.sin(Alphas[ii]),VyA[indexV] + Viy[indexV],VzA[indexV] + Viz[indexV] * m.cos(Alphas[ii]) + Vix[indexV] * m.sin(Alphas[ii]),sweepC4[indexV]);
            alphaLoc = np.arctan2(velTriZVisc,velTriXVisc);
            alphaLocV = np.arctan2(velTriYViscV,velTriXViscV);
            alphaLoc = np.concatenate((alphaLoc,alphaLocV));
            [clVisc,cd,cm] = p.AeroCoef(waPanel,wclPanel,wcdPanel,wcmPanel,alphaLoc);
            
            # Convergence criterion computation
            deltaAlpha = (clVisc-clw) / cl_alphaPanel; 
            nIter = 0;
            # apparent inviscid difference of a.o.a. such that the inviscid lift at AOA+deltaAlpha
            # Begining of the correction loop
            while np.amax(np.absolute(deltaAlpha))*2*m.pi > tol or np.any(np.isnan(deltaAlpha)):
                nIter = nIter +1;
                # If convergence too tough  : introduction of numerical damping
                if nIter == nIterMAX or np.any(np.isnan(deltaAlpha)):
                    nIter = 0;
                    K = 0.5*K;
                    if max(K) < 0.5:
                        K = 2.*K;
                        tol = 1e-3;
                        nCorrection = nCorrection +1;
#                        nIterMAX = 100;
                        if nCorrection == 1:
                            if ii == nbA-1:
                                Alphas[ii] = (Alphas[ii]+Alphas[ii-1])*0.5;
                            else:
                                Alphas[ii] = (Alphas[ii]+Alphas[ii+1])*0.5;
                            AOA = Alphas[ii];
                            
                            velTriX = V0 * np.cos(AOA) * np.cos(-Beta) + vix;
                            velTriY = V0 * np.sin(-Beta) + viy;
                            velTriZ = V0 * np.sin(AOA) * np.cos(-Beta) + viz;
                            AOA_VLM = np.concatenate([np.arctan2(velTriZ[index],velTriX[index]),np.arctan2(velTriY[indexV],velTriX[indexV])]);
                            Norme = np.concatenate([np.sqrt(velTriZ[index]*velTriZ[index]+velTriX[index]*velTriX[index]),np.sqrt(velTriY[indexV]*velTriY[indexV]+velTriX[indexV]*velTriX[indexV])]);
                            velTriX = Norme * np.cos(AOA_VLM);
                            velTriZ[index] = Norme[index] * np.sin(AOA_VLM[index]);
                            velTriY[indexV] = Norme[indexV] * np.sin(AOA_VLM[indexV]);
                            
                        elif nCorrection == 2:
                            stop = True;
                            break;
                    deltaAlpha = np.zeros(nbPan,dtype = float);
                artVisc = artViscLC(ac,alphaLoc,Alphas[1],Alphas[0],aMax,waPanel,wclPanel,cl_alphaPanel,Norme,c,y);
                AOA_VLM += np.dot(artVisc,deltaAlpha);
                velTriX = Norme * np.cos(AOA_VLM);
                velTriZ[index] = Norme[index] * np.sin(AOA_VLM[index]);
                velTriY[indexV] = Norme[indexV] * np.sin(AOA_VLM[indexV]);
                
                [Vix,Viy,Viz,ccl] = VLM_Solve(A,select,select2,Vx,Vy,Vz,velTriX,velTriY,velTriZ,normal);
                clw = ccl/(c);
                clw[indexV] *= -1;
                
                velTriXVisc,velTriYVisc,velTriZVisc = u.localVelTri(VxA[index] + Vix[index] * m.cos(Alphas[ii]) - Viz[index] * m.sin(Alphas[ii]),VyA[index] + Viy[index],VzA[index] + Viz[index] * m.cos(Alphas[ii]) + Vix[index] * m.sin(Alphas[ii]),tp,sweepC4[index],dih);
                velTriXViscV,velTriYViscV,velTriZViscV = u.localVelTriVT(VxA[indexV]+ Vix[indexV] * m.cos(Alphas[ii]) - Viz[indexV] * m.sin(Alphas[ii]),VyA[indexV] + Viy[indexV],VzA[indexV] + Viz[indexV] * m.cos(Alphas[ii]) + Vix[indexV] * m.sin(Alphas[ii]),sweepC4[indexV]);
                alphaLoc = np.arctan2(velTriZVisc,velTriXVisc);
                alphaLocV = np.arctan2(velTriYViscV,velTriXViscV);
                alphaLoc = np.concatenate((alphaLoc,alphaLocV));
                [clVisc,cd,cm] = p.AeroCoef(waPanel,wclPanel,wcdPanel,wcmPanel,alphaLoc);
                
                deltaAlpha = (clVisc-clw)/(cl_alphaPanel)*K;
            speedFactor = np.concatenate((velTriXVisc[index]*velTriXVisc[index] + velTriZVisc[index] *velTriZVisc[index] + velTriYVisc[index] *velTriYVisc[index], velTriXViscV*velTriXViscV +velTriZViscV*velTriZViscV + velTriYViscV *velTriYViscV ))/(V0*V0);
            al_i[:,ii] = (AOA_0 - alphaLoc);
            nCorrection = 0;
            K = np.ones(nbPan,dtype = float)
            nIterMAX = 25;
            tol = 1e-3;                  
            if stop:
#                print('WARNING : the tolerance on cl for alpha = '+str(180./m.pi*Alphas[ii])+' has raised up to '+str(np.amax(np.absolute(deltaAlpha)*cl_alphaPanel)));
                indexRes[ii] = False;
                deltaAlpha = np.zeros(nbPan);
                stop = False;
            else:
                deltaAlpha = AOA_VLM - AOA_VLM0;
                dF = clw * dA * speedFactor;
                dCDi = np.sin(al_i[:,ii]) * dF * np.cos(sweepC4);
                dCL = dF[index] * np.cos(al_i[index,ii]) * np.cos(dih[index])**2 - cd[index] *  dA[index] * speedFactor[index] * np.sin(al_i[index,ii]);
                dY = dF[indexV] * np.cos(al_i[indexV,ii]) - cd[indexV] *  dA[indexV] * speedFactor[indexV] * np.sin(al_i[indexV,ii]);
                dCD0 = cd *  dA * speedFactor * np.cos(al_i[:,ii]) ;
                dCM = np.sum(- dCL * lonDist[index] * np.cos(angle[index]-Alphas[ii])) + \
                    np.sum((dCDi + dCD0) * lonDist * np.sin(angle - Alphas[ii])) + \
                    np.sum(cm[index] * dA[index] * c[index] * speedFactor[index]);
                Cl[ii] = -np.sum(dCL*yDist[index]) / (S*MAC) + np.sum(dY * lonDist[indexV] *  np.sin(angle[indexV] - Alphas[ii])) / (S*MAC);
                Cn[ii] = (np.sum((dCD0+dCDi) * yDist) + N_F) / (S*MAC)\
                  - np.sum(dyp*prop.Tc/MAC) - np.sum(dY * lonDist[indexV] * np.cos(angle[indexV]-Alphas[ii])) / (S*MAC) + \
                  np.sum(cm[indexV] * dA[indexV] * c[indexV] * speedFactor[indexV])/ (S*MAC);
                CY[ii] = (Y_F + np.sum(dY))/S;
                CL[ii] = (np.sum(dCL)+L_F[ii])/S + np.sum(prop.Tc) * m.sin(Alphas[ii]) ;
                CDi[ii] = np.sum(dCDi)/S;
                CD0[ii] = (np.sum(dCD0)+D_F[ii])/S;
                CM[ii] = (dCM + M_F[ii])/(S*MAC) - np.sum(dzp*prop.Tc/MAC);
                al_i[:,ii] = 180./m.pi*(al_i[:,ii]);
                if ii < nbA-1:
                    artVisc = artViscLC(ac,alphaLoc,Alphas[ii],Alphas[ii+1],aMax,waPanel,wclPanel,cl_alphaPanel,Norme,c,y);
            
        CD = CDi+CD0+C_Dint;
        Alphas = 180./m.pi * Alphas;
        return CL[indexRes],CD[indexRes],CY[indexRes],CM[indexRes],Cl[indexRes],Cn[indexRes],CDi[indexRes],CD0[indexRes],Alphas[indexRes];
    
def artViscLC(ac,alphaLoc,AOA,AOA1,aMax,waPanel,wclPanel,cl_alphaPanel,Norme,c,y):
    nbPan = len(alphaLoc);
    artVisc = np.eye(nbPan);
    wing = ac.wing;
    htail = ac.htail;
    marge = 0.95;
    beta = m.log(0.5)*wing.r/wing.b;
    for kk in range(1,wing.r-1):
        if (alphaLoc[kk] + AOA1 - AOA) > aMax[kk]*marge:
            clM = np.interp(alphaLoc[kk] - AOA + AOA1,waPanel[kk],wclPanel[kk]);
            clm = np.interp(alphaLoc[kk],waPanel[kk],wclPanel[kk]);
            claLoc = (clM-clm)/(AOA1-AOA);
            artVisc[kk,:wing.r] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[:wing.r],beta);
    if (alphaLoc[0] + AOA1 - AOA) > aMax[0]*marge:
        kk = 0;
        clM = np.interp(alphaLoc[kk] + AOA1 - AOA,waPanel[kk],wclPanel[kk]);
        clm = np.interp(alphaLoc[kk],waPanel[kk],wclPanel[kk]);
        claLoc = (clM-clm)/( AOA1 - AOA);
        artVisc[kk,:wing.r] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[:wing.r],beta);
    if (alphaLoc[wing.r-1] + AOA1 - AOA) > aMax[wing.r-1]*marge:
        kk = wing.r-1
        clM = np.interp(alphaLoc[kk] + AOA1 - AOA,waPanel[kk],wclPanel[kk]);
        clm = np.interp(alphaLoc[kk],waPanel[kk],wclPanel[kk]);
        claLoc = (clM-clm)/( AOA1 - AOA);
        artVisc[kk,:wing.r] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[:wing.r],beta);
    if htail.bool:
        beta = m.log(0.5)*htail.r/htail.b;
        if (alphaLoc[wing.r] + AOA1 - AOA) > aMax[wing.r]*marge:
            kk = wing.r;
            clM = np.interp(alphaLoc[kk] + AOA1 - AOA,waPanel[kk],wclPanel[kk]);
            clm = np.interp(alphaLoc[kk],waPanel[kk],wclPanel[kk]);
            claLoc = (clM-clm)/( AOA1 - AOA);
            artVisc[kk,wing.r:] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[wing.r:],beta);
        if (alphaLoc[nbPan-1] + AOA1 - AOA) > aMax[nbPan-1]*marge:
            kk = nbPan-1;
            clM = np.interp(alphaLoc[kk] + AOA1 - AOA,waPanel[kk],wclPanel[kk]);
            clm = np.interp(alphaLoc[kk],waPanel[kk],wclPanel[kk]);
            claLoc = (clM-clm)/( AOA1 - AOA);
            artVisc[kk,wing.r:] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[wing.r:],beta);
        for kk in range(wing.r+1,wing.r+htail.r-1):
            if (alphaLoc[kk] + AOA1 - AOA) > aMax[kk]*marge:
                clM = np.interp(alphaLoc[kk]+ AOA1 - AOA,waPanel[kk],wclPanel[kk]);
                clm = np.interp(alphaLoc[kk],waPanel[kk],wclPanel[kk]);
                claLoc = (clM-clm)/(AOA1-AOA);
                artVisc[kk,wing.r:] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[wing.r:],beta);
    return artVisc;
def muCoef(claLoc,cla,U0,c,yStar,y,beta):
    """ Function that returns the multiplier coefficient from a decreasing exponential law"""
    mu = max([0,-0.5 * (claLoc/cla * (2./(cla*U0*c)+1.) - (2./(cla*U0*c)-1.))])*np.exp(beta*abs(yStar-y));
    mu[y==yStar] = - np.sum(mu[y!=yStar]);
    if -0.5 * (claLoc/cla * (2./(cla*U0*c)+1.) - (2./(cla*U0*c)-1.)) > 0:
        mu *= 2*(-0.5 * (claLoc/cla * (2./(cla*U0*c)+1.) - (2./(cla*U0*c)-1.)))/ mu[y==yStar];
    return mu
def StronglyCoupled(ac,flow,lin):
    # AC configuration
    wing = ac.wing;
    htail = ac.htail;
    vtail = ac.vtail;
    fus = ac.fus;
    prop = ac.prop;
    
    # Numerical parameters initialisation
    tol = 1e-3;                                                                 # Tolerance for the a.o.a. correction loop
    nIterMAX = 200;                                                              # Maximum number of iteration in the correction loop before introducing the numerical damping
    tolWarning = 2e-3;                                                          # Tolerance above which a warning is send to the promp because of lack of accuracy
    tolStop = 5e-2;                                                             # Tolerance above which the software is stopped because of lack of accuracy
    
    # Flight conditions
    Alphas = flow.getAlphas() * m.pi/180;
    Beta = m.pi/180.*flow.getBeta();
    V0 = flow.getV0();
    nbA = len(Alphas);
    # Geometric characteristics
    # Wing
    sweepC4 = wing.getSweepC4();
    for i in range(len(sweepC4)/2):
        sweepC4[i] *= -1.;
    dih = wing.getDih();
    for i in range(len(dih)/2):
        dih[i] *= -1.;
    tp = wing.getTwist();
    x = wing.getXP();
    y = wing.getYP();
    z = wing.getZP();
    c = wing.getChord();
    dA = 0.5*(wing.y[1:]-wing.y[:-1])*(wing.chordDistrib[1:]+wing.chordDistrib[:-1]);
    S = wing.getS();
    MAC = wing.getMac();
    [waPanel,wclPanel, wcdPanel, wcmPanel, alphal0PanelVLM, aMax, cl_alphaPanel] = p.panelPolar(wing);
    nbPan = len(y)
    #Fuselage
    yFus = fus.getYF();
    index = np.array(range(nbPan));
    # Aerodynamic contribution of fuselage
    if fus.bool:
        for i in range(nbPan):
            if (y[i] > yFus[0] and y[i] < yFus[1]):
                index[i] = -1;
        index = index[index >= 0];
        L_F= fus.getL();
        D_F= fus.getD();
        M_F= fus.getM();
        N_F = fus.getN();
        Y_F = fus.getY();
        c_r = np.interp(yFus[1],y,c);
        Rec_r = c_r*V0/(0.00001568);
        dih_r = np.interp(yFus[1],y,dih);
        tc = 0.075;
        C_Dint = (c_r**2)/S*(0.1112 - 0.2572*m.sin(dih_r) + 3.44 * tc - 0.02097 * m.log10(Rec_r) + 0.09009 * m.sin(dih_r) * m.sin(dih_r) - 2.549 * tc * m.sin(dih_r) + 0.03010 * m.log10(Rec_r) * m.sin(dih_r) -  0.1462 * tc * m.log10(Rec_r));
    else:
        C_Dint = 0.;
        L_F = np.zeros(nbA,dtype = float);
        Y_F = 0.;
        N_F = 0.;
        D_F = np.zeros(nbA,dtype = float);
        M_F = np.zeros(nbA,dtype = float);
    
    # Engines
    # Induced velocity by the engine
    # tangencial increases the local speed and then the local lift
    # contribution 
    vix = prop.vix;
    viy = prop.viy;
    viz = prop.viz;
    if prop.bool:
        dyp = prop.yp - ac.rp[1];
        dzp = prop.zp - ac.rp[2];
    else:
        dyp = 0;
        dzp = 0;
        prop.Tc = 0;
    
    # Htail
    if htail.bool:
        yTail = htail.getY();
        cTail = htail.getChordDist();
        dihT = htail.getDih();
        index = np.concatenate([index,np.array(range(index[-1]+1,index[-1]+1+htail.getR()))]);
        dA = np.concatenate([dA , 0.5 * (cTail[1:] + cTail[:-1]) * (yTail[1:] - yTail[:-1])]);
        tp = np.concatenate([tp,htail.getTwist()]);
        
        for i in range(htail.getR()/2):
            dihT[i] *= -1.;
        dih = np.concatenate([dih,dihT]);
        sweepC4T = htail.getSweepC4();
        for i in range(htail.getR()/2):
            sweepC4T[i] *= -1.;
        sweepC4 = np.concatenate([sweepC4,sweepC4T]);
        c = np.concatenate([c,htail.getChord()]);
        x = np.concatenate([x,htail.getXP()]);
        y = np.concatenate([y,htail.getYP()]);
        z = np.concatenate([z,htail.getZP()]);
        nbPan = len(y)
        [waPanelT,wclPanelT, wcdPanelT, wcmPanelT, alphal0PanelVLMT, alphaMaxT, cl_alphaPanelT] = p.panelPolar(htail);
        waPanel = np.concatenate([waPanel,waPanelT],0);
        wclPanel = np.concatenate([wclPanel,wclPanelT],0);
        wcdPanel = np.concatenate([wcdPanel,wcdPanelT],0);
        wcmPanel = np.concatenate([wcmPanel,wcmPanelT],0);
        alphal0PanelVLM = np.concatenate([alphal0PanelVLM,alphal0PanelVLMT]);
        aMax = np.concatenate([aMax,alphaMaxT]);
        cl_alphaPanel= np.concatenate([cl_alphaPanel,cl_alphaPanelT]);
        
    # Vtail
    if vtail.bool:
        CYV,CMV,ClV,CnV,CDiV,CD0V,Alphas = AeroPredVTail(vtail,flow,lin,S,MAC,ac.rp);
        L_F= np.interp(Alphas*180./m.pi, flow.getAlphas(), L_F);
        D_F= np.interp(Alphas*180./m.pi, flow.getAlphas(), D_F);
        M_F= np.interp(Alphas*180./m.pi, flow.getAlphas(), M_F);
        nbA = len(Alphas);
        if fus.bool:
            c_r = vtail.getChord(0);
            Rec_r = c_r*V0/(0.00001568);
            dih_r = 0.;
            C_Dint += (c_r**2)/S*(0.1112 - 0.2572*m.sin(dih_r) + 3.44 * tc - 0.02097 * m.log10(Rec_r) + 0.09009 * m.sin(dih_r) * m.sin(dih_r) - 2.549 * tc * m.sin(dih_r) + 0.03010 * m.log10(Rec_r) * m.sin(dih_r) -  0.1462 * tc * m.log10(Rec_r));
    else:
        CYV= np.zeros(nbA,float);
        CMV= np.zeros(nbA,float);
        ClV= np.zeros(nbA,float);
        CnV= np.zeros(nbA,float);
        CDiV= np.zeros(nbA,float);
        CD0V= np.zeros(nbA,float);
    
    # Initialisation and beginning of the computation
    deltaAlpha = np.zeros(nbPan);
    al_i = np.zeros([nbPan,nbA],dtype = float);
    CL = np.zeros(nbA,dtype = float);
    CY = np.zeros(nbA,dtype = float);
    CM = np.zeros(nbA,dtype = float);
    Cl = np.zeros(nbA,dtype = float);
    Cn = np.zeros(nbA,dtype = float);
    CDi = np.zeros(nbA,dtype = float);
    CD0 = np.zeros(nbA,dtype = float);

    lonDist = np.sqrt((x-ac.rp[0])**2 + (z-ac.rp[2])**2);
    angle = np.arctan2(z-ac.rp[2],x-ac.rp[0]);
    yDist = y - ac.rp[1];
    stop = False;
    # Construction of the A matrix
    [A,normal,Vx,Vy,Vz,select,select2] = VLM_MP.ICMatrix(ac,cl_alphaPanel,flow);
    nA = np.shape(A)[0];
    SA = np.zeros([nA + nbPan,nA+nbPan],dtype = float);
    SA[:nA,:nA] = np.linalg.inv(A);
    SA[nA:,nA:] = np.eye(nbPan);
    J = np.zeros([nA+nbPan,nA+nbPan]);
    J[:nA,:nA] = np.linalg.inv(A);
    claLoc = np.zeros(nbPan);
    mC = nA/nbPan;
    I = np.eye(nbPan);
    # Construction of the A matrix
    alphal0PanelVLM *= np.sum(select,1);
    if lin:
        for ii in range(nbA):
            # Construction of velocity triangles
            AOA = Alphas[ii];
            velTriX = V0 * np.cos(AOA) * np.cos(-Beta) + vix;
            velTriY = V0 * np.sin(-Beta) + viy;
            velTriZ = V0 * np.sin(AOA) * np.cos(-Beta) + viz;
            AOA_VLM = np.arctan2(velTriZ,velTriX);
            
            
            VxA = V0 * np.cos(Alphas[ii]) * np.cos(-Beta)+ vix;
            VyA = V0 * np.sin(-Beta) + viy;
            VzA = V0 * np.sin(Alphas[ii]) * np.cos(-Beta)+ viz;
            velTriXLoc,velTriYLoc,velTriZLoc = u.localVelTri(VxA,VyA,VzA,tp,sweepC4,dih);
            
            Vx0 = V0 * np.cos(Alphas[ii]) * np.cos(-Beta);
            Vy0 = V0 * np.sin(-Beta);
            Vz0 = V0 * np.sin(Alphas[ii]) * np.cos(-Beta);
            velTriX0,velTriY0,velTriZ0 = u.localVelTri(Vx0,Vy0,Vz0,tp,sweepC4,dih);
            af0 = np.arctan2(velTriZ0,velTriX0);
            
            [Vix,Viy,Viz,ccl] = VLM_Solve(A,Vx,Vy,select,select2,Vz,velTriX,velTriY,velTriZ,normal);
            clw = ccl / c;
            velTriXVisc,velTriYVisc,velTriZVisc = u.localVelTri(VxA + Vix[:nbPan] * m.cos(Alphas[ii]) - Viz[:nbPan] * m.sin(Alphas[ii]),VyA + Viy[:nbPan],VzA + Viz[:nbPan]* m.cos(Alphas[ii]) + Vix[:nbPan] * m.sin(Alphas[ii]),tp,sweepC4,dih);
            alphaLoc = np.arctan2(velTriZVisc,velTriXVisc);
            [clVisc,cd,cm] = p.AeroCoef(waPanel,wclPanel,wcdPanel,wcmPanel,alphaLoc);
            
            speedFactor = (velTriXVisc*velTriXVisc+velTriZVisc*velTriZVisc)/(V0*V0);
            al_i[:,ii] = af0 - alphaLoc;
            
            dF = clVisc *  dA * speedFactor;
            
            dCDi = np.sin(al_i[index,ii]) * dF[index] * np.cos(sweepC4);
            dCL = dF[index]* np.cos(al_i[index,ii])*np.cos(dih[index])**2 - cd[index] *  dA[index]*speedFactor[index] * np.sin(al_i[index,ii]);
            dCD0 = cd[index] *  dA[index]*speedFactor[index] * np.cos(al_i[index,ii]);
            dCM = - ( dCL * lonDist[index] * np.cos(angle[index]-Alphas[ii])) + \
                (dCDi + dCD0) * lonDist[index] * np.sin(angle[index] - Alphas[ii]) + \
                cm[index] * dA[index] * c[index] * speedFactor[index];
            Cl[ii] = -np.sum(dCL*yDist[index]) / (S*MAC) +ClV[ii];
            Cn[ii] = (np.sum((dCD0+dCDi)*yDist[index]) + N_F) / (S*MAC)\
              - np.sum(dyp*prop.Tc/MAC) + CnV[ii];
            CY[ii] = Y_F/S + CYV[ii];
            CL[ii] = (np.sum(dCL)+L_F[ii])/S + np.sum(prop.Tc) * m.sin(Alphas[ii]) ;
            CDi[ii] = np.sum(dCDi)/S + CDiV[ii];
            CD0[ii] = CD0V[ii] + (np.sum(dCD0)+D_F[ii])/S - np.sum(prop.Tc) * m.cos(Alphas[ii]);
            CM[ii] = CMV[ii] + (np.sum(dCM) + M_F[ii])/(S*MAC) - np.sum(dzp*prop.Tc/MAC);
            al_i[:,ii] = 180./m.pi*(al_i[:,ii]);
        CD = CD0+CDi+C_Dint;
        Alphas = Alphas*180./m.pi;
        return CL,CD,CY,CM,Cl,Cn,CDi,CD0,Alphas
    else:
        for ii in range(nbA):                                              # if possible no numerical damping to fasten the computation
            deltaAlpha0 = deltaAlpha;
            AOA = Alphas[ii];
            velTriX = V0 * np.cos(AOA) * np.cos(-Beta) + vix;
            velTriY = V0 * np.sin(-Beta) + viy;
            velTriZ = V0 * np.sin(AOA) * np.cos(-Beta) + viz;
            AOA_VLM = np.arctan2(velTriZ,velTriX);
            AOA_VLM0 = AOA_VLM;
            AOA_VLM += deltaAlpha0;
            Norme = np.sqrt(velTriZ*velTriZ+velTriX*velTriX);
            velTriX = Norme * np.cos(AOA_VLM);
            velTriZ = Norme * np.sin(AOA_VLM);
            norme = np.sqrt(velTriX*velTriX+velTriY*velTriY+velTriZ*velTriZ);
            
            VxA = V0 * np.cos(Alphas[ii]) * np.cos(-Beta) + vix;
            VyA = V0 * np.sin(-Beta) + viy;
            VzA = V0 * np.sin(Alphas[ii]) * np.cos(-Beta) + viz;
            velTriXLoc,velTriYLoc,velTriZLoc = u.localVelTri(VxA,VyA,VzA,tp,sweepC4,dih);
            alpha0 = np.arctan2(velTriZLoc,velTriXLoc);
            
            Vx0 = V0 * np.cos(Alphas[ii]) * np.cos(-Beta);
            Vy0 = V0 * np.sin(-Beta);
            Vz0 = V0 * np.sin(Alphas[ii]) * np.cos(-Beta);
            velTriX0,velTriY0,velTriZ0 = u.localVelTri(Vx0,Vy0,Vz0,tp,sweepC4,dih);
            AOA_0 = np.arctan2(velTriZ0,velTriX0);
            
            # Linear prediction of cl and induced velocities
            [Vix,Viy,Viz,ccl] = VLM_Solve(A,select,select2,Vx,Vy,Vz,velTriX,velTriY,velTriZ,normal);
            clw = ccl / c;
            
            # Viscous prediction based on prediction of induced velocity
            velTriXVisc,velTriYVisc,velTriZVisc = u.localVelTri(VxA + Vix[:nbPan],VyA + Viy[:nbPan],VzA + Viz[:nbPan],tp,sweepC4,dih);
            alphaLoc = np.arctan2(velTriZVisc,velTriXVisc);
            [clVisc,cd,cm] = p.AeroCoef(waPanel,wclPanel,wcdPanel,wcmPanel,alphaLoc);
            artVisc = artViscSC(ac,alphaLoc,0.,m.pi/360.,aMax,waPanel,wclPanel,cl_alphaPanel,Norme,c,y,nA);
            # Convergence criterion computation
            deltaAlpha = (clVisc-clw) / cl_alphaPanel; 
            # apparent inviscid difference of a.o.a. such that the inviscid lift at AOA+deltaAlpha
            
            ## Try to strong coupling
            
            RHS = -(normal[0,:] * np.dot(select2,velTriX) + normal[1,:] * np.dot(select2,velTriY) + normal[2,:] * np.dot(select2,velTriZ));
            gamma = np.dot(A,RHS);
            Vix = np.dot(Vx,gamma);
            Viy = np.dot(Vy,gamma);
            Viz = np.dot(Vz,gamma);
            X = np.concatenate([gamma,deltaAlpha0]);
            nbIter = 0;
            while(max(abs(deltaAlpha))>tol and nbIter < nIterMAX):
                b = np.concatenate([RHS,deltaAlpha]);
                F = np.dot(SA,X) - b - np.dot(artVisc,X);
                for kk in range(nbPan):
                    J[kk*mC:(kk+1)*mC,nA+kk] = - normal[0,kk*mC:(kk+1)*mC] * velTriZ[kk] + normal[2,kk*mC:(kk+1)*mC] * velTriX[kk];
                    clM = np.interp(alphaLoc[kk] + m.pi/360.,waPanel[kk],wclPanel[kk]);
                    clm = np.interp(alphaLoc[kk] - m.pi/360.,waPanel[kk],wclPanel[kk]);
                    claLoc[kk] = (clM-clm)*180./m.pi;
                    J[nA+kk,:nA] = -select[kk,:]/cl_alphaPanel[kk]*(claLoc[kk] * 2./(cl_alphaPanel[kk]*c[kk]*norme[kk]) - 2./(c[kk]*norme[kk]));
                    J[nA+kk,nA:] =  (1.+claLoc[kk]/cl_alphaPanel[kk])*I[kk,:] - artVisc[kk,nA:];
                X1 = X - np.dot(np.linalg.inv(J),F);
                RHS = -(normal[0,:]* np.dot(select2,velTriX) + normal[1,:] * np.dot(select2,velTriY) + normal[2,:] * np.dot(select2,velTriZ));
                clw = 2. * np.dot(select,gamma)/(norme*c);
                Vix = np.dot(Vx,gamma);
                Viy = np.dot(Vy,gamma);
                Viz = np.dot(Vz,gamma);
                velTriXVisc,velTriYVisc,velTriZVisc = u.localVelTri(VxA + Vix[:nbPan],VyA + Viy[:nbPan],VzA + Viz[:nbPan],tp,sweepC4,dih);
                alphaLoc2 = np.arctan2(velTriZVisc,velTriXVisc);
                [clVisc2,cd,cm] = p.AeroCoef(waPanel,wclPanel,wcdPanel,wcmPanel,alphaLoc2);
                alphaLoc = (clw/cl_alphaPanel + alphal0PanelVLM + alpha0 - (AOA_VLM + tp));
                [clVisc,cd,cm] = p.AeroCoef(waPanel,wclPanel,wcdPanel,wcmPanel,alphaLoc);
                deltaAlpha = (clVisc2-clw) / cl_alphaPanel;
#                print 'Alpha'
#                plt.plot(y,alphaLoc*180./m.pi)
#                plt.plot(y,alphaLoc2*180./m.pi)
#                plt.show()
#                print 'CL'
#                plt.plot(y,clw)
#                plt.plot(y,clVisc)
#                plt.plot(y,clVisc2)
#                plt.show()
#                print 'DeltaAlpha'
#                plt.plot(y,X[nA:]*180./m.pi)
#                plt.plot(y,X1[nA:]*180./m.pi)
#                plt.plot(y,X1[nA:]*180./m.pi-X[nA:]*180./m.pi)
#                plt.plot(y,deltaAlpha*180./m.pi)
#                plt.show()
                X = X1
                artVisc = artViscSC(ac,alphaLoc,0.,m.pi/360.,aMax,waPanel,wclPanel,cl_alphaPanel,Norme,c,y,nA);
                gamma = X[:nA];
                
                AOA_VLM = AOA_VLM0 + X[nA:];
                velTriX = Norme * np.cos(AOA_VLM);
                velTriZ = Norme * np.sin(AOA_VLM);
                nbIter += 1;
            speedFactor = (velTriXVisc*velTriXVisc+velTriZVisc*velTriZVisc)/(V0*V0);
            al_i[:,ii] = (AOA_0 - alphaLoc);
            if np.amax(np.absolute(deltaAlpha))*2.*m.pi >= tolStop or stop:
                print('WARNING : the tolerance on cl for alpha = '+str(180./m.pi*Alphas[ii])+' has raised up to '+str(np.amax(np.absolute(deltaAlpha)*cl_alphaPanel)));
                print('Because of the too low precison, the calcul is stopped here');
                Alphas = 180./m.pi*(Alphas[:ii]);
                CL = CL[:ii];
                CDi = CDi[:ii];
                CD0 = CD0[:ii];
                CY = CY[:ii];
                CM = CM[:ii];
                Cl = Cl[:ii];
                Cn = Cn[:ii];
                CD = CDi + CD0 + C_Dint;
                return CL,CD,CY,CM,Cl,Cn,CDi,CD0,Alphas;
            elif np.amax(abs(deltaAlpha)) >= tolWarning:
                print('WARNING : the tolerance on cl for alpha = '+str(180./m.pi*Alphas[ii])+' has raised up to '+str(np.amax(np.absolute(deltaAlpha)*cl_alphaPanel)));
#            nIterMAX = 50;
            tol = 1e-3;                                                         # reset the tolerance for next computation
            deltaAlpha = AOA_VLM - AOA_VLM0;
            dF = clw *  dA * speedFactor;
            dCDi = np.sin(al_i[index,ii]) * dF[index] * np.cos(sweepC4)[index];
            dCL = dF[index] * np.cos(al_i[index,ii]) * np.cos(dih[index])**2 - cd[index] *  dA[index] * speedFactor[index] * np.sin(al_i[index,ii]);
            dCD0 = (cd[index]) *  dA[index] * speedFactor[index] * np.cos(al_i[index,ii]) ;
            dCM = - dCL * lonDist[index] * np.cos(angle[index]-Alphas[ii]) + \
                (dCDi + dCD0) * lonDist[index] * np.sin(angle[index] - Alphas[ii]) + \
                cm[index] * dA[index] * c[index] * speedFactor[index];
            Cl[ii] = -np.sum(dCL*yDist[index]) / (S*MAC) + ClV[ii];
            Cn[ii] = (np.sum((dCD0+dCDi) * yDist[index]) + N_F) / (S*MAC)\
              - np.sum(dyp*prop.Tc/MAC) + CnV[ii];
            CY[ii] = Y_F/S + CYV[ii];
            CL[ii] = (np.sum(dCL)+L_F[ii])/S + np.sum(prop.Tc) * m.sin(Alphas[ii]) ;
            CDi[ii] = np.sum(dCDi)/S + CDiV[ii];
            CD0[ii] = CD0V[ii] + (np.sum(dCD0)+D_F[ii])/S;
            CM[ii] = CMV[ii] + (np.sum(dCM) + M_F[ii])/(S*MAC) - np.sum(dzp*prop.Tc/MAC);
            al_i[:,ii] = 180./m.pi*(al_i[:,ii]);
        CD = CDi+CD0+C_Dint;
        Alphas = 180./m.pi * Alphas;
        return CL,CD,CY,CM,Cl,Cn,CDi,CD0,Alphas;
def artViscSC(ac,alphaLoc,AOA,AOA1,aMax,waPanel,wclPanel,cl_alphaPanel,Norme,c,y,nA):
    nbPan = len(alphaLoc);
    artVisc = np.zeros([nbPan+nA,nbPan+nA]);
    wing = ac.wing;
    htail = ac.htail;
    marge = 0.9;
    beta = m.log(0.5)*wing.r/wing.b;
    for kk in range(1,wing.r-1):
        if (alphaLoc[kk] + AOA1 - AOA) > aMax[kk]*marge:
            clM = np.interp(alphaLoc[kk] + AOA1,waPanel[kk],wclPanel[kk]);
            clm = np.interp(alphaLoc[kk] + AOA,waPanel[kk],wclPanel[kk]);
            claLoc = (clM-clm)/(AOA1-AOA);
            artVisc[nA+kk,nA:nA+wing.r] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[:wing.r],beta);
    if (alphaLoc[0] + AOA1 - AOA) > aMax[0]*marge:
        kk = 0;
        clM = np.interp(alphaLoc[kk] + AOA1 - AOA,waPanel[kk],wclPanel[kk]);
        clm = np.interp(alphaLoc[kk],waPanel[kk],wclPanel[kk]);
        claLoc = (clM-clm)/( AOA1 - AOA);
        artVisc[nA+kk,nA:nA+wing.r] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[:wing.r],beta);
    if (alphaLoc[wing.r-1] + AOA1 - AOA) > aMax[wing.r-1]*marge:
        kk = wing.r-1
        clM = np.interp(alphaLoc[kk] + AOA1 - AOA,waPanel[kk],wclPanel[kk]);
        clm = np.interp(alphaLoc[kk],waPanel[kk],wclPanel[kk]);
        claLoc = (clM-clm)/( AOA1 - AOA);
        artVisc[nA+kk,nA:nA+wing.r] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[:wing.r],beta);
    if htail.bool:
        beta = m.log(0.5)*htail.r/htail.b;
        if (alphaLoc[wing.r] + AOA1 - AOA) > aMax[wing.r]*marge:
            kk = wing.r;
            clM = np.interp(alphaLoc[kk] + AOA1 - AOA,waPanel[kk],wclPanel[kk]);
            clm = np.interp(alphaLoc[kk],waPanel[kk],wclPanel[kk]);
            claLoc = (clM-clm)/( AOA1 - AOA);
            artVisc[nA+kk,nA+wing.r:] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[wing.r:],beta);
        if (alphaLoc[nbPan-1] + AOA1 - AOA) > aMax[nbPan-1]*marge:
            kk = nbPan-1;
            clM = np.interp(alphaLoc[kk] + AOA1 - AOA,waPanel[kk],wclPanel[kk]);
            clm = np.interp(alphaLoc[kk],waPanel[kk],wclPanel[kk]);
            claLoc = (clM-clm)/( AOA1 - AOA);
            artVisc[nA+kk,nA+wing.r:] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[wing.r:],beta);
        for kk in range(wing.r+1,wing.r+htail.r-1):
            if (alphaLoc[kk] + AOA1 - AOA) > aMax[kk]*marge:
                clM = np.interp(alphaLoc[kk]+ AOA1 - AOA,waPanel[kk],wclPanel[kk]);
                clm = np.interp(alphaLoc[kk],waPanel[kk],wclPanel[kk]);
                claLoc = (clM-clm)/(AOA1-AOA);
                artVisc[nA+kk,nA+wing.r:] += muCoef(claLoc,cl_alphaPanel[kk],Norme[kk],c[kk],y[kk],y[wing.r:],beta);
    return artVisc;