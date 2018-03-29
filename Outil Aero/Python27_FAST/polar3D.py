# -*- coding: utf-8 -*-
import math as m
import numpy as np
import numpy.matlib
import numpy.linalg
import copy
import Polar
import Panel
def panelPolar(wing,Mach):
    """
    
     Autor : Quentin Borlon
     This def provides the different polar data required for the rest pf
     the code based on the wing.NLPol files.
     
     INPUT:
       wing.getR() : number of control station along the span +1
       wing.polarIndex : array of index. wing.polarIndex(i) indicate the index
       of the file in wing.NLPol to use for the i-th station along the span
       
     OUTPUT:
       [wcl, wcd] : [(wing.getR()+1) x 31] matrix. array i gives the weight
           coefficients for 30 degree legendre interpolation the lift/drag 
            curves for the i-th station along the span, taken the y<0 part into
           account!
       al0Section : (wing.getR()+1)-array giving the zero-lift-angle 
           of the station along the span, taken the y<0 part into
           account!
       cl_alphaSection : (wing.getR()+1)-array giving the lift curve slope 
           of the station along the span, taken the y<0 part into
           account!"""

    n = 40;                                                                 # Degre of the Legendre polynomial approximation
    if Mach < 1.:
        Beta = np.sqrt(1. - (Mach*np.cos(wing.sweepc4))**2);
    else:
        Beta = np.sqrt((Mach*np.cos(wing.sweepc4))**2-1.);
    if np.any(wing.getPolD()):
        wclPanel = np.zeros([wing.getR(),n+1]);
        wcdPanel = np.zeros([wing.getR(),n+1]);
        al0Panel = np.zeros(wing.getR(),dtype = float);
        alMaxPan = np.zeros(wing.getR(),dtype = float);
        al0PanelVLM = np.zeros(wing.getR(),dtype = float);
        cl_alphaPanel = np.zeros(wing.getR(),dtype = float);
        
        thetaDisc = wing.getThetaDPol();
        yPanel = (wing.getY(range(1,wing.getR()+1))+wing.getY(range(0,wing.getR())))/2;
        mu = np.zeros([len(thetaDisc),wing.getR()],dtype = float);
        muPos = np.zeros([len(thetaDisc),wing.getR()],dtype = float);
        muNeg = np.zeros([len(thetaDisc),wing.getR()],dtype = float);
        deltaclmaxStar = np.zeros(len(thetaDisc),dtype = float);
        # Apply the polar modification to ensure continuity of the clmax
        for ii in range(len(thetaDisc)):
            thetaStar1 = thetaDisc[ii];
            mStar = np.where(np.abs(wing.getTheta()-(thetaStar1-0.0001))<=0.0002)[0][0];
            # characteristic distance of influence of the discontinuity =
            # len of the flaps
            dist = wing.getChordDist(mStar);
            distneg = dist * (1.-wing.getCF());
            distpos = dist * (1.+wing.getCF());
            muNeg[ii,:] = exponentialCoef(wing.b/2*m.cos(thetaStar1),yPanel,distneg);
            muPos[ii,:] = exponentialCoef(wing.b/2*m.cos(thetaStar1),yPanel,distpos);
            
            # Compute the initial delta in cl_max at the discontinuity
            deltaclmaxStar[ii] = wing.getPolar(wing.getPolI(mStar)).getClMax()-\
                wing.getPolar(wing.getPolI(int(mStar-np.sign(m.cos(thetaStar1))))).getClMax(); # le plus extérieur moins le plus intérieur
            if deltaclmaxStar[ii]<0 and np.sign(m.cos(thetaStar1)) < 0:
                mu[ii,:mStar] = 0.6*muPos[ii,:mStar];
                mu[ii,mStar:] = -0.4*muNeg[ii,mStar:];
            elif deltaclmaxStar[ii] > 0 and np.sign(m.cos(thetaStar1)) < 0:
                mu[ii,:mStar] = -0.4*muNeg[ii,:mStar];
                mu[ii,mStar:] = 0.6*muPos[ii,mStar:];
            elif deltaclmaxStar[ii] > 0 and np.sign(m.cos(thetaStar1)) > 0:
                mu[ii,:mStar] = 0.6*muPos[ii,:mStar];
                mu[ii,mStar:] = -0.4*muNeg[ii,mStar:];
            else:
                mu[ii,:mStar] = -0.4*muNeg[ii,:mStar];
                mu[ii,mStar:] = 0.6*muPos[ii,mStar:];
        deltaclmaxStar = np.absolute(deltaclmaxStar);
        # Polar altering like in Sivells Westrick NACA TN 2283
        for jj in range(wing.getR()/2):
            Polars = copy.deepcopy(wing.getPolar(wing.getPolI(jj+1)));
            clmaxAlt = Polars.getClMax() + np.dot(mu[:,jj],deltaclmaxStar);
            coef = clmaxAlt/Polars.getClMax();
            Polars.setCl(coef*Polars.getCl()*Beta[jj+1]/(Beta[jj+1]+np.absolute(Polars.getCl())*0.5*(Mach**2)/(1.+Beta[jj+1])));
            Polars.setAlpha(coef*(Polars.getAlpha()-Polars.getAl0())+Polars.getAl0());
            Polars.setClMax(clmaxAlt);
            al0Panel[jj] = m.pi/180.*(Polars.getAl0());
            if wing.getDsF(jj+1):
                iiClmin1 = np.where(Polars.getAlpha() <= -6)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= -2)[0][0];
            else:
                iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
            if not(np.any(iiClmin1)):
                iiClmin1 = m.floor(len(Polars.getAlpha())*0.25);
            if not(np.any(iiClmax1)):
                iiClmax1 = m.floor(len(Polars.getAlpha())*0.5);
            cl_alphaPanel[jj] = 180./m.pi*(Polars.getCl(iiClmax1)-Polars.getCl(iiClmin1)) / (Polars.getAlpha(iiClmax1)-Polars.getAlpha(iiClmin1));
            al0PanelVLM[jj]= Polars.getAl0VLM()*m.pi/180.;
            angle = np.linspace(Polars.getAlphaMax()*0.8,Polars.getAlphaMax()*1.2,101)*coef;
            alMaxPan[jj] = m.pi/180*angle[np.polyval(wclPanel[jj,:],angle) == max(np.polyval(wclPanel[jj,:],angle))][0];
            [wclPanel[jj,:],wcdPanel[jj,:]] = intrapolLegend(Polars,n);
        for jj in range(wing.getR()/2,wing.getR()):
            Polars = copy.deepcopy(wing.getPolar(wing.getPolI(jj)));
            clmaxAlt = Polars.getClMax() + np.dot(mu[:,jj],deltaclmaxStar);
            coef = clmaxAlt/Polars.getClMax();
            Polars.setCl(coef*Polars.getCl()*Beta[jj]/(Beta[jj]+np.absolute(Polars.getCl())*0.5*(Mach**2)/(1.+Beta[jj])));
            Polars.setAlpha(coef*(Polars.getAlpha()-Polars.getAl0())+Polars.getAl0());
            Polars.setClMax(clmaxAlt);
            al0Panel[jj] = m.pi/180.*(Polars.getAl0());
            if wing.getDsF(jj):
                iiClmin1 = np.where(Polars.getAlpha() <= -6)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= -2)[0][0];
            else:
                iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
            if not(np.any(iiClmin1)):
                iiClmin1 = m.floor(len(Polars.getAlpha())*0.25);
            if not(np.any(iiClmax1)):
                iiClmax1 = m.floor(len(Polars.getAlpha())*0.5);
            cl_alphaPanel[jj] = 180./m.pi*(Polars.getCl(iiClmax1)-Polars.getCl(iiClmin1)) / (Polars.getAlpha(iiClmax1)-Polars.getAlpha(iiClmin1));
            al0PanelVLM[jj]= Polars.getAl0VLM()*m.pi/180.;
            [wclPanel[jj,:],wcdPanel[jj,:]] = intrapolLegend(Polars,n);
            angle = np.linspace(Polars.getAlphaMax()*0.8,Polars.getAlphaMax()*1.2,101)*coef;
            alMaxPan[jj] = m.pi/180*angle[np.polyval(wclPanel[jj,:],angle) == max(np.polyval(wclPanel[jj,:],angle))][0];
    else:
        wclPanel = np.zeros([wing.getR(),n+1]);
        wcdPanel = np.zeros([wing.getR(),n+1]);
        cl_alphaPanel = np.zeros(wing.getR(),dtype = float);
        for ii in range(wing.getR()/2):
            Polars = copy.deepcopy(wing.getPolar(0));
            Polars.setCl(Polars.getCl()*Beta[ii+1]/(Beta[ii+1]+np.absolute(Polars.getCl())*0.5*(Mach**2)/(1.+Beta[ii+1])));
            if wing.getDsF(ii+1):
                iiClmin1 = np.where(Polars.getAlpha() <= -6)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= -4)[0][0];
            else:
                iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
            if not(np.any(iiClmin1)):
                iiClmin1 = m.floor(len(Polars.getAlpha())*0.25);
            if not(np.any(iiClmax1)):
                iiClmax1 = m.floor(len(Polars.getAlpha())*0.5);
            cl_alphaPanel[ii] = 180./m.pi*(Polars.getCl(iiClmax1)-Polars.getCl(iiClmin1))/ (Polars.getAlpha(iiClmax1)-Polars.getAlpha(iiClmin1));
            [wclPanel[ii,:],wcdPanel[ii,:]] = intrapolLegend(Polars,n);
        for ii in range(wing.getR()/2,wing.getR()):
            Polars = copy.deepcopy(wing.getPolar(0));
            Polars.setCl(Polars.getCl()*Beta[ii]/(Beta[ii]+np.absolute(Polars.getCl())*0.5*(Mach**2)/(1.+Beta[ii])));
            if wing.getDsF(ii):
                iiClmin1 = np.where(Polars.getAlpha() <= -6)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= -2)[0][0];
            else:
                iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
            if not(np.any(iiClmin1)):
                iiClmin1 = m.floor(len(Polars.getAlpha())*0.25);
            if not(np.any(iiClmax1)):
                iiClmax1 = m.floor(len(Polars.getAlpha())*0.5);
            cl_alphaPanel[ii] = 180./m.pi*(Polars.getCl(iiClmax1)-Polars.getCl(iiClmin1))/ (Polars.getAlpha(iiClmax1)-Polars.getAlpha(iiClmin1));
            [wclPanel[ii,:],wcdPanel[ii,:]] = intrapolLegend(Polars,n);
        angle = np.linspace(Polars.getAlphaMax()*0.8,Polars.getAlphaMax()*1.2,101);
        alMaxPan = m.pi/180*angle[np.polyval(wclPanel[0,:],angle) == max(np.polyval(wclPanel[0,:],angle))][0]*np.ones(wing.getR(),dtype = float);
        al0Panel = m.pi/180*(Polars.getAl0())*np.ones(wing.getR(),dtype = float);
        al0PanelVLM = Polars.getAl0VLM()*m.pi/180*np.ones(wing.getR(),dtype = float);
    return wclPanel, wcdPanel, alMaxPan,al0PanelVLM,cl_alphaPanel;



def intrapolLegend(polar,n):
    """Autor : Quentin Borlon, based on the work of
    
           Copyright (c) 2001, Jordi / Esther Soler / Bonet
           All rights reserved.
     
       This def compute the weight w_i of the legendre approximation for
     the lift and drag section curve. Using the gauus points as control points
     during the approximation step. It allows to go faster during the
     recovering of the section data at a given a.o.a. and reduces the 
     induced numerical oscillations at high angle of attack.
    
    
     INPUT:
       polar : the section polar to interpolate.
       n : order of the legendre approximation
    
     Output:
       [wcl,wcd] : the weight of the legendre polynomes for lift/drag curve
       approximation."""
    i = polar.getAlpha(0);
    j = polar.getAlpha(-1);
    x=lgwt(n+1,i,j);
    #n=len(x)-1;
    a = x;
    #cl = polar.getCl(x);
    cl = np.interp(x,polar.getAlpha(),polar.getCl());
    alphas = np.zeros([n+1,n+1],dtype= float);
    expo = np.array(range(n,-1,-1),dtype = float);
    for i in range(n+1):
        alphas[i,:] = a[i]**expo;
    pol = polegende(n);
    mat = np.dot(alphas,np.transpose(pol));
    w = numpy.linalg.solve(mat,cl);
    wcl = np.dot(w,pol);
    cd = np.interp(x,polar.getAlpha(),polar.getCd());
    w = numpy.linalg.solve(mat,cd);
    wcd = np.dot(w,pol);
    return wcl,wcd;

def polegende(n):
    """Autor : Quentin Borlon, based on the work of
    
          Copyright (c) 2001, Jordi / Esther Soler / Bonet
          All rights reserved.
    
      This def returns a square matrix (n+1) containing the coefficients 
    of the legendre polynomes. Cell p(n+1-i,n+1-j) corresponds to the 
    the coefficient before the term in (j) order for the polynome of order
    (i)
      
    INPUT:
      n: order of the greatest polynome
        
    OUTPUT:
      p : the matrix containing the coefficients"""

    p = np.zeros([n+1,n+1],dtype = float);
    p[n,n]=1;
    p[n-1,n-1:]=[1,0];
    for k in range(2,n+1):
        p[n-k,n-k:]=((2.*k-1.)*np.concatenate([[0],p[n+1-k,n+1-k:]])-(k-1)*np.concatenate([p[n+2-k,n+2-k:],[0,0]]))/k;
    return p;

def lgwt(n,a,b):
    
    """This script is for computing definite integrals using Legendre-Gauss
         Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
         [a,b] with truncation order N
        
         Suppose you have a continuous def f(x) which is defined on [a,b]
         which you can evaluate at any x in [a,b]. Simply evaluate it at all of
         the values contained in the x vector to obtain a vector f. Then compute
         the definite integral using sum(f.*w);
        
         Written by Greg von Winckel - 02/25/2004
        
         INPUT:
           n: order of the greatest polynome
           [a,b] : bound of the interval for the gauss points
        
         OUTPUT:
           x : row of gauss point in the desired interval"""
    n=n-1;
    N1=n+1; 
    N2=n+2;

    xu=np.linspace(-1,1,N1);

    # Initial guess
    y=np.cos((2*(np.array(range(N1)))+1)*m.pi/(2*n+2))+(0.27/N1)*np.sin(m.pi*xu*n/N2);
    # Legendre-Gauss Vandermonde Matrix
    L=np.zeros([N1,N2],dtype = float);

    # Derivative of LGVM
    Lp=np.zeros([N1,N2],dtype = float);

    # Compute the np.zeros of the N+1 Legendre Polynomial
    # using the recursion relation and the Newton-Raphson method

    y0=2;

    # Iterate until new points are uniformly within epsilon of old points
    while np.amax(np.absolute(y-y0))>1e-3:


        L[:,0] = np.ones(N1,dtype = float);
        Lp[:,0] = np.zeros(N1,dtype = float);

        L[:,1] = y;
        Lp[:,1] = np.ones(N1,dtype = float);

        for k in range(2,N1+1):
            L[:,k] = ( (2*k-1)*y*L[:,k-1]-(k-1)*L[:,k-2] )/(k);

        der = N2*( L[:,N1-1]- y*L[:,N2-1] )/(1.-y**2);

        y0 = y;
        y = y0 - L[:,N2-1]/der;

    # Linear map from[-1,1] to [a,b]
    x=(a*(1.-y)+b*(1.+y))*0.5;
    return x;
    
        
    
def exponentialCoef(yStar,y,dist):
    """ Function that returns the multiplier coefficient from a decreasing exponential law"""
    yStar2 = np.round(yStar,3);
    y2 = np.round(y,3);
    mu = np.exp(-np.absolute(y2-yStar2)*3/dist);
    return mu;