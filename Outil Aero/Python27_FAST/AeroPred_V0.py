# -*- coding: utf-8 -*-
import math as m;
import numpy as np;
import Polar as p;
import utilitaire as u
import VLM
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
def main(ac,flow,lin):
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

    wing = ac.wing;
    htail = ac.htail;
    vtail = ac.vtail;
    fus = ac.fus;
    prop = ac.prop;
    # Numerical parameters initialisation
    tol = 1e-3;                                                                 # Tolerance for the a.o.a. correction loop
    nIterMAX = 25;                                                              # Maximum number of iteration in the correction loop before introducing the numerical damping
    tolWarning = 2e-3;                                                          # Tolerance above which a warning is send to the promp because of lack of accuracy
    tolStop = 5e-2;                                                             # Tolerance above which the software is stopped because of lack of accuracy

    # Flight conditions
    Alphas = flow.getAlphas() * m.pi/180;
    Beta = flow.getBeta();
    V0 = flow.getV0();
    # Geometric characteristics
    sweepC4 = wing.getSweepC4();
    for i in range(len(sweepC4)/2):
        sweepC4[i] *= -1;
    dih = wing.getDih();
    tp = wing.getTwist();
    
    
    x = wing.getXP();
    y = wing.getYP();
    z = wing.getZP();
    c = wing.getChord();
    
    dA = 0.5*(wing.y[1:]-wing.y[:-1])*(wing.chordDistrib[1:]+wing.chordDistrib[:-1]);
    S = wing.getS();
    MAC = wing.getMac();
    
    yFus = fus.getYF();
    
    index = np.array(range(wing.getR()));
    indexV = [];
    
    [wclPanel, wcdPanel, wcmPanel, wcdaPanel,alphal0PanelVLM,cl_alphaPanel] = p.panelPolar(wing);
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
    
    # Aerodynamic contribution of fuselage
    L_F = np.zeros(len(Alphas),dtype = float);
    Y_F = 0.;
    N_F = 0.;
    D_F = np.zeros(len(Alphas),dtype = float);
    M_F = np.zeros(len(Alphas),dtype = float);
    
    if fus.bool:
        for i in range(wing.getR()):
            if (y[i] > yFus[0] and y[i] < yFus[1]):
                index[i] = -1;
        index = index[index >= 0];
        L_F= np.interp(flow.getAlphas(), u.rangeF(flow.aMin-20.,flow.aMax+20.,0.1), fus.getL());
        D_F= np.interp(flow.getAlphas(), u.rangeF(flow.aMin-20.,flow.aMax+20.,0.1), fus.getD());
        M_F= np.interp(flow.getAlphas(), u.rangeF(flow.aMin-20.,flow.aMax+20.,0.1), fus.getM());
        N_F = fus.getN();
        Y_F = fus.getY();
    
    if htail.bool:
        yTail = htail.getY();
        cTail = htail.getChordDist();
        index = np.concatenate([index,np.array(range(index[-1]+1,index[-1]+1+htail.getR()))]);
        dA = np.concatenate([dA , 0.5 * (cTail[1:] + cTail[:-1]) * (yTail[1:] - yTail[:-1])]);
        tp = np.concatenate([tp,htail.getTwist()]);
        dih = np.concatenate([dih,htail.getDih()]);
        sweepC4T = htail.getSweepC4();
        for i in range(htail.getR()/2):
            sweepC4T[i] *= -1;
        sweepC4 = np.concatenate([sweepC4,sweepC4T]);
        c = np.concatenate([c,htail.getChord()]);
        x = np.concatenate([x,htail.getXP()]);
        y = np.concatenate([y,htail.getYP()]);
        z = np.concatenate([z,htail.getZP()]);
        
        [wclPanelT, wcdPanelT, wcmPanelT, wcdaPanelT,alphal0PanelVLMT,cl_alphaPanelT] = p.panelPolar(htail);
        wclPanel = np.concatenate([wclPanel,wclPanelT],0);
        wcdPanel = np.concatenate([wcdPanel,wcdPanelT],0);
        wcmPanel = np.concatenate([wcmPanel,wcmPanelT],0);
        wcdaPanel = np.concatenate([wcdaPanel,wcdaPanelT],0);
        alphal0PanelVLM = np.concatenate([alphal0PanelVLM,alphal0PanelVLMT]);
        cl_alphaPanel= np.concatenate([cl_alphaPanel,cl_alphaPanelT]);
        
    if vtail.bool:
        CDV,CYV,CMV,ClV,CnV,CDiV,CD0V,Alphas = AeroPredVTail(vtail,flow,lin,S,MAC,ac.rp);
    # Construction of the A matrix
    
    [A,normal] = VLM.ICMatrix(ac,cl_alphaPanel);
    invA = np.linalg.inv(A);
    nbPan = len(y);
    nbA = len(Alphas);
    
    betaEff = sweepC4-Beta;
    Af = np.transpose(np.matlib.repmat(Alphas,nbPan,1));
    Af_eff = Af;
    
    twist_eff = tp;
    Af_twist = np.transpose(np.matlib.repmat(Alphas,nbPan,1)) + twist_eff;
    
    A0 = np.transpose(np.matlib.repmat(Alphas,nbPan,1)) - alphal0PanelVLM;
    
    # Physical flow velocities triangle
    vix0 = V0 * np.cos(betaEff) * np.cos(Af_eff);
    viy0 = V0 * np.sin(betaEff);
    viz0 = V0 * np.cos(betaEff) * np.sin(Af_eff);
    # Physical flow velocities triangle in the pannel frame
    vixT = V0 * np.cos(betaEff) * np.cos(Af_twist);
    viyT = viy0;
    vizT = V0 * np.cos(betaEff) * np.sin(Af_twist);
    # Flow velocities taken into account the zero-lift angle
    vixA0 = V0 * np.cos(betaEff) * np.cos(A0);
    vizA0 = V0 * np.cos(betaEff) * np.sin(A0);
    
    # Initialisation and beginning of the computation
    deltaAlpha = np.zeros(nbPan);
    # Intialization of data containers
    al_i = np.zeros([nbPan,nbA],dtype = float);
    CL = np.zeros(nbA,dtype = float);
    CY = np.zeros(nbA,dtype = float);
    CM = np.zeros(nbA,dtype = float);
    Cl = np.zeros(nbA,dtype = float);
    Cn = np.zeros(nbA,dtype = float);
    CDi = np.zeros(nbA,dtype = float);
    CD0 = np.zeros(nbA,dtype = float);
    
    #dV = np.zeros([3,nbPan],dtype = float);
    lonDist = np.sqrt((x-ac.rp[0])**2 + (z-ac.rp[2])**2);
    angle = np.arctan2(z-ac.rp[2],x-ac.rp[0])
    yDist = y - ac.rp[1];
    if lin:
        for ii in range(nbA):
            vin = viz+ viz0[ii,:];
            speedFactor = ((vix+vix0[ii])**2+vin**2)/(V0**2);
            vNref = np.concatenate([vizT[ii,index]+viz[index], viyT[indexV]+viy[indexV]]);
            vTref = np.concatenate([vixT[ii] + vix]);
            alphaRef = np.arctan2(vNref,vTref);
            a0z = viz + vizA0[ii];
            a0x = vix + vixA0[ii];
            a0 = np.arctan2(a0z,a0x);
            AOA = a0 + deltaAlpha - alphal0PanelVLM;
            ccl = VLM_Solve(invA, AOA ,normal , nbPan);
            clw = ccl/(c*np.cos(sweepC4));
            al_e = clw/(cl_alphaPanel);
            alphaLoc = alphaRef - (AOA - al_e);
            [clVisc,cd,cm,cda] = p.AeroCoef(wclPanel,wcdPanel,wcmPanel,wcdaPanel,alphaLoc);
            al_i[:,ii] = Alphas[ii] - alphaLoc + tp;
            dF = clw *  dA * speedFactor;
            dCDi = np.sin(al_i[index,ii]) * dF[index];
            dCL = dF[index]* np.cos(al_i[index,ii])*np.cos(dih[index])**2 - cd[index] *  dA[index]*speedFactor[index] * np.sin(al_i[index,ii]);
            dCD0 = cd[index] *  dA[index]*speedFactor[index] * np.cos(al_i[index,ii]);
            dCM = - ( dCL * lonDist[index] * np.cos(angle[index]-Alphas[ii])) + \
                (dCDi + dCD0) * lonDist[index] * np.sin(angle[index] - Alphas[ii]) + \
                cm[index] * dA[index] * c[index] * speedFactor[index];
            Cl[ii] = -np.sum(dCL*yDist[index]) / (S*wing.getMac()) +ClV[ii];
            Cn[ii] = (np.sum((dCD0+dCDi)*yDist[index]) + N_F) / (S*wing.getMac())\
              - np.sum(dyp*prop.Tc/wing.getMac()) + CnV[ii];
            CY[ii] = Y_F/S + CYV[ii];
            CL[ii] = (np.sum(dCL)+L_F[ii])/S + np.sum(prop.Tc) * m.sin(Alphas[ii]) ;
            CDi[ii] = np.sum(dCDi)/S + CDiV[ii];
            CD0[ii] = CD0V[ii] + (np.sum(dCD0)+D_F[ii])/S - np.sum(prop.Tc) * m.cos(Alphas[ii]);
            CM[ii] = CMV[ii] + (np.sum(dCM) + M_F[ii])/(S*wing.getMac()) - np.sum(dzp*prop.Tc/wing.getMac());
            al_i[:,ii] = 180./m.pi*(al_i[:,ii]);
        CD = CD0+CDi;
        Alphas = Alphas*180./m.pi;
        return CL,CD,CY,CM,Cl,Cn,CDi,CD0,Alphas
    else:
        for ii in range(nbA):
            K = np.ones(nbPan,dtype = float);                                   # if possible no numerical damping to fasten the computation
            vin = viz+ viz0[ii,:];
            speedFactor = ((vix+vix0[ii])**2+vin**2)/(V0**2);
            vNref = np.concatenate([vizT[ii,index]+viz[index], viyT[indexV]+viy[indexV]]);
            vTref = np.concatenate([vixT[ii] + vix]);
            alphaRef = np.arctan2(vNref,vTref);
            a0z = viz + vizA0[ii];
            a0x = vix + vixA0[ii];
            a0 = np.arctan2(a0z,a0x);
            AOA = a0 + deltaAlpha - alphal0PanelVLM;
            ccl = VLM_Solve(invA, AOA ,normal , nbPan);
            clw = ccl/(c*np.cos(sweepC4));
            al_e = clw/(cl_alphaPanel);
            alphaLoc = alphaRef - (AOA - al_e);
            [clVisc,cd,cm,cda] = p.AeroCoef(wclPanel,wcdPanel,wcmPanel,wcdaPanel,alphaLoc);
            deltaAlpha = (clVisc-clw)/(cl_alphaPanel*K);                       # apparent inviscid difference of a.o.a. such that the inviscid lift at AOA+deltaAlphas
            # is the viscous lift at alphaLoc
            nIter = 0;
            nCorrection = 0;
            # Begining of the correction loop
            while np.amax(np.absolute(deltaAlpha))*2*m.pi > tol or np.any(np.isnan(deltaAlpha)):
                nIter = nIter +1;
                # If convergence too tough  : introduction of numerical damping
                if nIter == nIterMAX or np.any(np.isnan(deltaAlpha)):
                    tol = tol*2;
                    nIter = 0;
                    deltaAlpha = np.zeros(nbPan,dtype = float);
                    K = 2*K;
                    if tol > tolStop:
                        tol = 1e-3;
                        nCorrection = nCorrection +1;
                        if nCorrection == 1:
                            nIterMAX = 2*nIterMAX;
                            K = 2*K;
                            if ii == nbA-1:
                                Alphas[ii] = (Alphas[ii]+Alphas[ii-1])/2;
                            else:
                                Alphas[ii] = (Alphas[ii]+Alphas[ii+1])/2;
                            #vNref = np.concatenate([vizT[ii,index]+viz[index], viyT[indexV]+viy[indexV]]);
                            #vTref = np.concatenate([vixT[ii] + vix]);
                            #alphaRef = np.arctan2(vNref,vTref);
                        elif nCorrection == 2:
                            K = np.ones(nbPan,dtype = float);
                            for kk in range(nbPan):
                                indice = 0;
                                for idi in range(ii):
                                    if np.isnan(al_i[kk,idi]):
                                        pass;
                                    else:
                                        indice = idi;
                                alpha = (Alphas[indice]+ tp[kk])*180./m.pi-al_i[kk,indice];
                                cl_alphaPanel[kk] = 180./m.pi*((np.polyval(np.polyder(wclPanel[kk,:]),alpha-0.2))+np.polyval(np.polyder(wclPanel[kk,:]),alpha-0.2))/(0.4);
                                cl = np.polyval(wclPanel[kk,:],alpha);
                                alphal0PanelVLM[kk] = -cl/cl_alphaPanel[kk] +alpha*m.pi/180;
                                while np.absolute(cl_alphaPanel[kk]) <4:
                                    alpha = alpha+0.5;
                                    cl_alphaPanel[kk] = 180./m.pi*((np.polyval(np.polyder(wclPanel[kk,:]),alpha-0.2))+np.polyval(np.polyder(wclPanel[kk,:]),alpha-0.2))/(0.4);
                                    cl = np.polyval(wclPanel[kk,:],alpha);
                                    alphal0PanelVLM[kk] = -cl/cl_alphaPanel[kk] +alpha*m.pi/180;
                            [A,normal] = VLM.ICMatrix(ac,cl_alphaPanel);
                            invA = np.linalg.inv(A);
                        elif nCorrection == 3:
                            break;
                AOA += deltaAlpha;
                ccl = VLM_Solve(invA, AOA,normal, nbPan);
                clw = ccl/(c*np.cos(sweepC4));
                al_e = clw/(cl_alphaPanel);

                alphaLoc = alphaRef - (AOA - al_e);
                [clVisc,cd,cm,cda] = p.AeroCoef(wclPanel,wcdPanel,wcmPanel,wcdaPanel,alphaLoc);
                deltaAlpha = (clVisc-clw)/(cl_alphaPanel*K);                # recovering lift and drag coeff
            al_i[:,ii] = Alphas[ii] - alphaLoc + tp;
            # Warning for user in case of bad coputation point
            if np.amax(np.absolute(deltaAlpha)) >= tolStop:
                print('Attention, the tolerance on cl for alpha = '+str(180./m.pi*Alphas[ii]-wing.getiW())+' has raised up to ',str(np.amax(np.absolute(deltaAlpha)*cl_alphaPanel)));
                print( 'Because of the too low precison, the calcul is stopped here');
                Alphas = 180./m.pi*(Alphas[:ii]);
                CL = CL[:ii];
                CDi = CDi[:ii];
                CD0 = CD0[:ii];
                CY = CY[:ii];
                CM = CM[:ii];
                Cl = Cl[:ii];
                Cn = Cn[:ii];
                CD = CDi + CD0;
                return CL,CD,CY,CM,Cl,Cn,CDi,CD0,Alphas;
            elif np.amax(abs(deltaAlpha)) >= tolWarning:
                print('Attention, the tolerance on cl for alpha = '+str(180./m.pi*Alphas[ii])+' has raised up to '+str(np.amax(np.absolute(deltaAlpha)*cl_alphaPanel)));
            nIterMAX = 100;
            tol = 1e-3;                    
            
            dF = clw *  dA * speedFactor;
            dCDi = np.sin(al_i[index,ii]) * dF[index];
            dCL = dF[index]* np.cos(al_i[index,ii])*np.cos(dih[index])**2 - cd[index] *  dA[index]*speedFactor[index] * np.sin(al_i[index,ii]);
            dCD0 = cd[index] *  dA[index]*speedFactor[index] * np.cos(al_i[index,ii]);
            dCM = - ( dCL * lonDist[index] * np.cos(angle[index]-Alphas[ii])) + \
                (dCDi + dCD0) * lonDist[index] * np.sin(angle[index] - Alphas[ii]) + \
                cm[index] * dA[index] * c[index] * speedFactor[index];
            Cl[ii] = -np.sum(dCL*yDist[index]) / (S*wing.getMac()) +ClV[ii];
            Cn[ii] = (np.sum((dCD0+dCDi)*yDist[index]) + N_F) / (S*wing.getMac())\
              - np.sum(dyp*prop.Tc/wing.getMac()) + CnV[ii];
            CY[ii] = Y_F/S + CYV[ii];
            CL[ii] = (np.sum(dCL)+L_F[ii])/S + np.sum(prop.Tc) * m.sin(Alphas[ii]) ;
            CDi[ii] = np.sum(dCDi)/S + CDiV[ii];
            CD0[ii] = CD0V[ii] + (np.sum(dCD0)+D_F[ii])/S - np.sum(prop.Tc) * m.cos(Alphas[ii]);
            CM[ii] = CMV[ii] + (np.sum(dCM) + M_F[ii])/(S*wing.getMac()) - np.sum(dzp*prop.Tc/wing.getMac());
            al_i[:,ii] = 180./m.pi*(al_i[:,ii]);
            
        CD = CDi+CD0;
        Alphas=180./m.pi*(Alphas);
        return CL,CD,CY,CM,Cl,Cn,CDi,CD0,Alphas;

def VLM_Solve(invA, AOA,normal, n):
    RHS = np.zeros(n,dtype = float);
    Qinf=np.array([np.cos(AOA),np.zeros(n),np.sin(AOA)]);
    for b in range(n):
        RHS[b]=-np.dot(Qinf[:,b],normal[:,b]);
    return 2* np.dot(invA,RHS);

def VLM_SolveV(invA, velTriX,velTriY,velTriZ,normal, n):
    RHS = np.zeros(n,dtype = float);
    #computing the right hand side
    for b in range(n):
        vel = np.array([velTriX[b], velTriY[b], velTriZ[b]]);
        norme = m.sqrt(velTriX[b]*velTriX[b]+velTriY[b]*velTriY[b]+velTriZ[b]*velTriZ[b]);
        RHS[b]=-np.dot(vel,normal[:,b])/norme;
    return 2* np.dot(invA,RHS);

def AeroPredVTail(vtail,flow,lin,S,MAC,rp):
    """ Function that computes independently the aerodynamics of the VTP.
    Considered independent for following reasons:
        - Interaction with HTP leads to infinite values because of the too
            short distance at the connnection. May not be considered to obtain
            finite values.
        - Interaction with the wing is really negligible.
    
    author : Quentin Borlon
    date : 28-10-2017
    
    INPUT:
            - vtail : object that contain all the data of the vertical tail;
            - flow : object that contains all the data of the flow;
            - lin : boolean that gives the information for linear study or not;
            - S : wing's area;
            - MAC : wing's mean aerodynamic chord;
            - rp : reference point for the computation of moment coefficients.
    OUTPUT:
            - CD,CY,CDi,CD0 : forces coefficients;
            - CM,Cl,Cn : moment coefficients;
            - Alphas : AOA at which the convergence was reached."""
    z = vtail.getZ();
    c = vtail.getChordDist();
    dA = 0.5 * (c[1:] + c[:-1]) * (z[1:] - z[:-1]);
    sweepC4 = vtail.getSweepC4();
    c = vtail.getChord();
    x = vtail.getXP();
    y = vtail.getYP();
    z = vtail.getZP();
    
    [wclPanel, wcdPanel, wcmPanel, wcdaPanel,alphal0PanelVLM,cl_alphaPanel] = p.vtailPolar(vtail);
    
    tol = 1e-3;                                                                 # Tolerance for the a.o.a. correction loop
    nIterMAX = 25;                                                              # Maximum number of iteration in the correction loop before introducing the numerical damping
    tolWarning = 2e-3;                                                          # Tolerance above which a warning is send to the promp because of lack of accuracy
    tolStop = 5e-2;                                                             # Tolerance above which the software is stopped because of lack of accuracy

    # Flight conditions
    Alphas = flow.getAlphas() * m.pi/180;
    Beta = flow.getBeta()*m.pi/180.;
    V0 = flow.getV0();
    
    [A,normal] = VLM.ICMatrixV(vtail,cl_alphaPanel);
    invA = np.linalg.inv(A);
    nbPan = len(y);
    nbA = len(Alphas);
    
    betaEff = -Beta*np.ones(vtail.getR(),dtype=float);
    betaVert = sweepC4;
    Af = np.transpose(np.matlib.repmat(Alphas,nbPan,1));
    Af_eff = Af+betaVert;
    
    A0 = np.transpose(np.matlib.repmat(Alphas,nbPan,1)) - alphal0PanelVLM;
    A0V = alphal0PanelVLM-Beta;
    
    # Physical flow velocities triangle
    vix0 = V0 * np.cos(betaEff) * np.cos(Af_eff);
    viy0 = V0 * np.sin(betaEff);
    # Flow velocities taken into account the zero-lift angle
    vixA0 = V0 * np.cos(betaEff) * np.cos(A0);
    viyA0 = V0 * np.sin(A0V);
    vizA0 = V0 * np.cos(betaEff) * np.sin(A0);
    
    # Initialisation and beginning of the computation
    deltaAlpha = np.zeros(nbPan);
    # Intialization of data containers
    al_i = np.zeros([nbPan,nbA],dtype = float);
    CY = np.zeros(nbA,dtype = float);
    CM = np.zeros(nbA,dtype = float);
    Cl = np.zeros(nbA,dtype = float);
    Cn = np.zeros(nbA,dtype = float);
    CDi = np.zeros(nbA,dtype = float);
    CD0 = np.zeros(nbA,dtype = float);
    
    dV = np.zeros([3,nbPan],dtype = float);
    lonDist = np.sqrt((x-rp[0])**2 + (z-rp[2])**2);
    angle = np.arctan2(z-rp[2],x-rp[0]);
    if lin:
        for ii in range(nbA):
            speedFactor = (vix0[ii]**2+viy0**2)/(V0**2);
#            vNref = np.concatenate([vizT[ii,index]+viz[index], viyT[indexV]+viy[indexV]]);
#            vTref = np.concatenate([vixT[ii] + vix]);
#            alphaRef = np.arctan2(vNref,vTref);
            velTriX = vixA0[ii] + dV[0];
            velTriY = viyA0 + dV[1];
            velTriZ = vizA0[ii] + dV[2];
#            velTriN = np.concatenate([velTriZ[index],velTriY[indexV]]);
            ccl = VLM_SolveV(invA, velTriX,velTriY,velTriZ,normal, nbPan);
            clw = ccl/(c*np.cos(sweepC4));
            al_e = clw/(cl_alphaPanel);
#            alI = al_e - np.arctan2(velTriN,velTriX) - tp;
#            alphaLoc = alphaRef + alI;
            alphaLoc = (al_e + alphal0PanelVLM);
            [clVisc,cd,cm,cda] = p.AeroCoef(wclPanel,wcdPanel,wcmPanel,wcdaPanel,alphaLoc);
            al_i[:,ii] = Alphas[ii] - alphaLoc;
            dF = clw *  dA * speedFactor;
            dCDi = np.sin(al_i[:,ii]) * dF;
            dCD0 = cd * dA * speedFactor * np.cos(al_i[:,ii]);
            dY = dF * np.cos(al_i[:,ii]) - cd *  dA * speedFactor * np.sin(al_i[:,ii]);
            dCM = (dCDi + dCD0) * lonDist * np.sin(angle - Alphas[ii]);
            Cl[ii] = np.sum(dY * lonDist *  np.sin(angle - Alphas[ii])) / (S * MAC);
            Cn[ii] = np.sum(dY * lonDist * np.cos(angle-Alphas[ii])) / (S * MAC);
            CY[ii] = np.sum(dY)/S;
            CDi[ii] = np.sum(dCDi)/S;
            CD0[ii] = np.sum(dCD0)/S;
            CM[ii] = np.sum(dCM)/(S * MAC);
        CD = CD0+CDi;
        return CD,CY,CM,Cl,Cn,CDi,CD0,Alphas
    else:
        for ii in range(nbA):
            K = np.ones(nbPan,dtype = float);  
            speedFactor = (vix0[ii]**2+viy0**2)/(V0**2);
#            vNref = np.concatenate([vizT[ii,index]+viz[index], viyT[indexV]+viy[indexV]]);
#            vTref = np.concatenate([vixT[ii] + vix]);
#            alphaRef = np.arctan2(vNref,vTref);
            velTriX = vixA0[ii] + dV[0];
            velTriY = viyA0 + dV[1];
            velTriZ = vizA0[ii] + dV[2];
            velTriXRef = vixA0[ii];
            velTriYRef = viyA0;
            velTriZRef = vizA0[ii];
#            velTriN = np.concatenate([velTriZ[index],velTriY[indexV]]);
            ccl = -VLM_SolveV(invA, velTriX,velTriY,velTriZ,normal, nbPan);
            clw = ccl/(c*np.cos(sweepC4));
            al_e = clw/(cl_alphaPanel);
#            alI = al_e - np.arctan2(velTriN,velTriX) - tp;
#            alphaLoc = alphaRef + alI;
            alphaLoc = (al_e + alphal0PanelVLM);
            [clVisc,cd,cm,cda] = p.AeroCoef(wclPanel,wcdPanel,wcmPanel,wcdaPanel,alphaLoc);
            deltaAlpha = (clVisc-clw)/(cl_alphaPanel*K);                       # apparent inviscid difference of a.o.a. such that the inviscid lift at AOA+deltaAlphas
            # is the viscous lift at alphaLoc
            nIter = 0;
            nCorrection = 0;
            # Begining of the correction loop
            while np.amax(np.absolute(deltaAlpha))*2*m.pi > tol or np.any(np.isnan(deltaAlpha)):
                nIter = nIter +1;
                # If convergence too tough  : introduction of numerical damping
                if nIter == nIterMAX or np.any(np.isnan(deltaAlpha)):
                    tol = tol*2;
                    nIter = 0;
                    deltaAlpha = np.zeros(nbPan,dtype = float);
                    K = 2*K;
                    if tol > tolStop:
                        tol = 1e-3;
                        nCorrection = nCorrection +1;
                        if nCorrection == 1:
                            nIterMAX = 2*nIterMAX;
                            K = 2*K;
                            if ii == nbA-1:
                                Alphas[ii] = (Alphas[ii]+Alphas[ii-1])/2;
                            else:
                                Alphas[ii] = (Alphas[ii]+Alphas[ii+1])/2;
                            #vNref = np.concatenate([vizT[ii,index]+viz[index], viyT[indexV]+viy[indexV]]);
                            #vTref = np.concatenate([vixT[ii] + vix]);
                            #alphaRef = np.arctan2(vNref,vTref);
                        elif nCorrection == 2:
                            K = np.ones(nbPan,dtype = float);
                            for kk in range(nbPan):
                                indice = 0;
                                for idi in range(ii):
                                    if np.isnan(al_i[kk,idi]):
                                        pass;
                                    else:
                                        indice = idi;
                                alpha = (Beta)*180./m.pi-al_i[kk,indice];
                                cl_alphaPanel[kk] = 180./m.pi*((np.polyval(np.polyder(wclPanel[kk,:]),alpha-0.2))+np.polyval(np.polyder(wclPanel[kk,:]),alpha-0.2))/(0.4);
                                cl = np.polyval(wclPanel[kk,:],alpha);
                                alphal0PanelVLM[kk] = -cl/cl_alphaPanel[kk] +alpha*m.pi/180;
                                while np.absolute(cl_alphaPanel[kk]) <4:
                                    alpha = alpha+0.5;
                                    cl_alphaPanel[kk] = 180./m.pi*((np.polyval(np.polyder(wclPanel[kk,:]),alpha-0.2))+np.polyval(np.polyder(wclPanel[kk,:]),alpha-0.2))/(0.4);
                                    cl = np.polyval(wclPanel[kk,:],alpha);
                                    alphal0PanelVLM[kk] = -cl/cl_alphaPanel[kk] +alpha*m.pi/180;
                            [A,normal] = VLM.ICMatrixV(vtail,cl_alphaPanel);
                            invA = np.linalg.inv(A);
                        elif nCorrection == 3:
                            break;
                    dV = np.zeros([3,nbPan],dtype = float);
                    velTriX = vixA0[ii] + dV[0];
                    velTriY = viyA0 + dV[1];
                    velTriZ = vizA0[ii] + dV[2];
                dV[0] = V0 * deltaAlpha * m.sin(Beta) * np.ones(nbPan);
                dV[1] = V0 * deltaAlpha * m.cos(Beta) * np.ones(nbPan);
                dV[2] = V0 * deltaAlpha * np.zeros(nbPan);
                velTriX = vixA0[ii] + dV[0];
                velTriY = viyA0 + dV[1];
                velTriZ = vizA0[ii] + dV[2];
                ccl = -VLM_SolveV(invA, velTriX,velTriY,velTriZ,normal, nbPan);
                clw = ccl/(c*np.cos(sweepC4));
                al_e = clw/(cl_alphaPanel);
                #            alI = al_e - np.arctan2(velTriN,velTriX) - tp;
                #            alphaLoc = alphaRef + alI;
                alphaLoc = (al_e + alphal0PanelVLM);
                [clVisc,cd,cm,cda] = p.AeroCoef(wclPanel,wcdPanel,wcmPanel,wcdaPanel,alphaLoc);
                deltaAlpha = (clVisc-clw)/(cl_alphaPanel*K);
            [clVisc,cd,cm,cda] = p.AeroCoef(wclPanel,wcdPanel,wcmPanel,wcdaPanel,alphaLoc);                # recovering lift and drag coeff
            al_i[:,ii] = -Beta - alphaLoc;
            # Warning for user in case of bad coputation point
            if np.amax(np.absolute(deltaAlpha)) >= tolStop:
                Alphas = 180./m.pi*(Alphas[:ii]);
                CDi = CDi[:ii];
                CD0 = CD0[:ii];
                CY = CY[:ii];
                CM = CM[:ii];
                Cl = Cl[:ii];
                Cn = Cn[:ii];
                CD = CDi + CD0;
                return CD,CY,CM,Cl,Cn,CDi,CD0,Alphas;
            elif np.amax(abs(deltaAlpha)) >= tolWarning:
                print('Attention, the tolerance on cl for alpha = '+str(180./m.pi*Alphas[ii])+' has raised up to '+str(np.amax(np.absolute(deltaAlpha)*cl_alphaPanel)));
            nIterMAX = 100;
            tol = 1e-3;                                                         # reset the tolerance for next computation
            dV[0] = velTriX - velTriXRef;
            dV[1] = velTriY - velTriYRef;
            dV[2] = velTriZ - velTriZRef;
            dF = clw *  dA * speedFactor;
            dCDi = np.sin(al_i[:,ii]) * dF;
            dCD0 = cd *  dA * speedFactor * np.cos(al_i[:,ii]);
            dY = dF * np.cos(al_i[:,ii]) - cd *  dA * speedFactor * np.sin(al_i[:,ii]);
            dCM = (dCDi + dCD0) * lonDist * np.sin(angle - Alphas[ii]);
            Cl[ii] = np.sum(dY * lonDist *  np.sin(angle - Alphas[ii])) / (S*MAC);
            Cn[ii] = np.sum(dY * lonDist * np.cos(angle-Alphas[ii])) / (S*MAC);
            CY[ii] = np.sum(dY)/S;
            CDi[ii] = np.sum(dCDi)/S;
            CD0[ii] = np.sum(dCD0)/S;
            CM[ii] = (np.sum(dCM))/(S*MAC);
            al_i[:,ii] = 180./m.pi*(al_i[:,ii]);
        CD = CDi+CD0;
        return CD,CY,CM,Cl,Cn,CDi,CD0,Alphas;
    
    