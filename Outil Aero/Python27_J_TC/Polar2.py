# -*- coding: utf-8 -*-
import math as m
import numpy as np
from copy import deepcopy
import sys
import matplotlib.pyplot as plt
class Polar2:
    """ Classe définissant une polaire 2D.
        Classe munie des champs suivants:
            - clmax : maximum lift coefficient
            - clAlpha : lift curve slope
            - al0 : zero-lift angle
            - al0VLM : angle such that the projection of the linear part of the polar cross the cl=0 line
            - data : table with aerodynamic data

     """
    def __init__(self):
        """ Définition des attributs pour le premier constructeur """
        self.clAlpha = 2*m.pi;
        self.clmax = 1.5;
        self.al0 = 0;
        self.al0VLM = 0;
        self.alphaMax = 0;
        self.data = np.zeros([551,5],dtype = float);
        self.data[:,0] = np.linspace(-55,55,551);
        self.filePath = '';

    def setClMax(self,cl):
        self.clmax = cl;
    def getClMax(self):
        return self.clmax;

    def setClAlpha(self,clAlpha):
        self.clAlpha = clAlpha;
    def getClAlpha(self):
        return self.clAlpha;

    def setAl0(self,alpha):
        self.al0 = alpha;
    def getAl0(self):
        return self.al0;

    def setAl0VLM(self,alpha):
        self.al0VLM = alpha;
    def getAl0VLM(self):
        return self.al0VLM;

    def setAlphaMax(self,alpha):
        self.alphaMax = alpha;
    def getAlphaMax(self):
        return self.alphaMax;

    def setAlpha(self,alpha):
        self.data[:,0] = alpha;
    def getAlpha(self,ii='all'):
        if ii == 'all':
            return self.data[:,0];
        else:
            return self.data[ii,0];

    def setData(self,alpha,cl,cd,cm = False ,cda = False):
        self.data = np.zeros([len(alpha),5],dtype = float);
        self.data[:,0] = alpha;
        self.data[:,1] = cl;
        self.data[:,2] = cd;
        if type(cm) == bool:
            self.data[:,3] = np.zeros(len(alpha),dtype = float);
        else:
            self.data[:,3] = cm;
        if type(cda) == bool:
            self.data[:,4] = np.zeros(len(alpha),dtype = float);
        else:
            self.data[:,4] = cda;
        
    def getData(self,ii='all'):
        if ii == 'all':
            return self.data;
        else:
            return self.data[ii,:];

    def setCl(self,cl):
        self.data[:,1] = cl;
    def getCl(self,ii='all'):
        if ii == 'all':
            return self.data[:,1];
        else:
            return self.data[ii,1];

    def setCd(self,cd):
        self.data[:,2] = cd;
    def getCd(self,ii='all'):
        if ii == 'all':
            return self.data[:,2];
        else:
            return self.data[ii,2];
    def setCm(self,cm):
        self.data[:,3] = cm;
    def getCm(self,ii='all'):
        if ii == 'all':
            return self.data[:,3];
        else:
            return self.data[ii,3];

    def setFP(self,fp):
        self.filePath = fp;
    def getFP(self):
        return self.filePath;

def panelPolar(wing):
    """
    
     Autor : Quentin Borlon
     This def provides the altered polar data required for the reste of
     the code based on the 2D polar files.
     
     INPUT:
       wing.getR() : number of control station along the span +1
       wing.polarIndex : array of index. wing.polarIndex(i) indicate the index
       of the file in wing.NLPol to use for the i-th station along the span
       
     OUTPUT:
       [wcl, wcd, wcm, wcda] : [(wing.getR()) x 41] matrix.
           Array i gives the weight coefficients for 40 degree legendre 
           interpolation the lift/drag:piutching moment and additional drag 
           curves for the i-th panel along the span. Aktered data like 
           in Sivells
       al0 : (wing.getR())-array giving the zero-lift-angle 
           of the panel along the span linear approximation from the linear part
           of the polar curves
       cl_alpha : wing.getR()-array giving the lift curve slope 
           of the panel along the span"""

    n = 200;                   
    wa = np.zeros([wing.getR(),n]);                                            
    wcl = np.zeros([wing.getR(),n]);
    wcd = np.zeros([wing.getR(),n]);
    wcm = np.zeros([wing.getR(),n]);
    al0 = np.zeros(wing.getR(),dtype = float);
    alphaMax =  np.zeros(wing.getR(),dtype = float);
    cl_alpha = np.zeros(wing.getR(),dtype = float);
    if np.any(wing.getPolD()):
        yDisc = wing.getPolD();
        yPanel = wing.getYP();
        mu = np.zeros([len(yDisc),wing.getR()],dtype = float);
        muPos = np.zeros([len(yDisc),wing.getR()],dtype = float);
        muNeg = np.zeros([len(yDisc),wing.getR()],dtype = float);
        deltaclmaxStar = np.zeros(len(yDisc),dtype = float);
        # Apply the polar modification to ensure continuity of the clmax
        for ii in range(len(yDisc)):
            yDisc1 = yDisc[ii];
            mStar = np.where(np.abs(wing.getY()-(yDisc1-0.0001))<=0.0002)[0][0];
            # characteristic distance of influence of the discontinuity =
            # len of the flaps
            dist = wing.getChordDist(mStar);
            distneg = dist * (1.-wing.getCF());
            distpos = dist * (1.+wing.getCF());
            muNeg[ii,:] = exponentialCoef(yDisc1,yPanel,distneg);
            muPos[ii,:] = exponentialCoef(yDisc1,yPanel,distpos);
            side = 0;
            if yDisc1 > 0.:
                side = 1;
            # Compute the initial delta in cl_max at the discontinuity
            deltaclmaxStar[ii] = wing.getPolar(side,wing.getPolI(mStar)).getClMax()-\
                wing.getPolar(side,wing.getPolI(int(mStar - 1))).getClMax(); #  gauche - droit
            if deltaclmaxStar[ii]<0 and yDisc1 < 0:
                mu[ii,:mStar] = -0.4*muNeg[ii,:mStar];
                mu[ii,mStar:] = 0.6*muPos[ii,mStar:];
            elif deltaclmaxStar[ii] > 0 and yDisc1 < 0:
                mu[ii,:mStar] = 0.6*muPos[ii,:mStar];
                mu[ii,mStar:] = -0.4*muNeg[ii,mStar:];
            elif deltaclmaxStar[ii] > 0 and yDisc1 > 0:
                mu[ii,:mStar] = 0.6*muPos[ii,:mStar];
                mu[ii,mStar:] = -0.4*muNeg[ii,mStar:];
            else:
                mu[ii,:mStar] = -0.4*muNeg[ii,:mStar];
                mu[ii,mStar:] = 0.6*muPos[ii,mStar:];
        deltaclmaxStar = np.absolute(deltaclmaxStar);
        # Polar altering like in Sivells Westrick NACA TN 2283
        for jj in range(wing.getR()/2):
            Polars = deepcopy(wing.getPolar(0,wing.getPolI(jj)));
            if wing.getDFL(wing.getR()/2 - 1 - jj):
                iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
            else:
                iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
            if not(np.any(iiClmin1)):
                iiClmin1 = m.floor(len(Polars.getAlpha())*0.25);
            if not(np.any(iiClmax1)):
                iiClmax1 = m.floor(len(Polars.getAlpha())*0.5);
            cl_alpha[jj] = 180./m.pi*(Polars.getCl(iiClmax1)-Polars.getCl(iiClmin1)) / (Polars.getAlpha(iiClmax1)-Polars.getAlpha(iiClmin1));
            clmaxAlt = Polars.getClMax() + np.dot(mu[:,jj],deltaclmaxStar);
            coef = clmaxAlt/Polars.getClMax();
            Polars.setCl(coef*Polars.getCl());
            Polars.setAlpha(coef*(Polars.getAlpha()-Polars.getAl0())+Polars.getAl0());
            Polars.setClMax(clmaxAlt);
            al0[jj] = m.pi/180.*(Polars.getAl0());
            alphaMax[jj] = m.pi/180.*(Polars.getAlphaMax())*coef;
            wa[jj,:] = np.linspace(Polars.getAlpha(0),Polars.getAlpha(-1),n);
            wcl[jj,:] = np.interp(wa[jj,:],Polars.getAlpha(),Polars.getCl());
            wcd[jj,:] = np.interp(wa[jj,:],Polars.getAlpha(),Polars.getCd());
            wcm[jj,:] = np.interp(wa[jj,:],Polars.getAlpha(),Polars.getCm());
            wa[jj,:] *= m.pi/180.;
        for jj in range(wing.getR()/2,wing.getR()):
            Polars = deepcopy(wing.getPolar(1,wing.getPolI(jj)));
            if wing.getDFR(jj-wing.getR()/2):
                iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
            else:
                iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
                iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
            if not(np.any(iiClmin1)):
                iiClmin1 = m.floor(len(Polars.getAlpha())*0.25);
            if not(np.any(iiClmax1)):
                iiClmax1 = m.floor(len(Polars.getAlpha())*0.5);
            cl_alpha[jj] = 180./m.pi*(Polars.getCl(iiClmax1)-Polars.getCl(iiClmin1)) / (Polars.getAlpha(iiClmax1)-Polars.getAlpha(iiClmin1));
            clmaxAlt = Polars.getClMax() + np.dot(mu[:,jj],deltaclmaxStar);
            coef = clmaxAlt/Polars.getClMax();
            Polars.setCl(coef*Polars.getCl());
            Polars.setAlpha(coef*(Polars.getAlpha()-Polars.getAl0())+Polars.getAl0());
            Polars.setClMax(clmaxAlt);
            al0[jj] = m.pi/180.*(Polars.getAl0());
            alphaMax[jj] = m.pi/180.*(Polars.getAlphaMax())*coef;
            wa[jj,:] = np.linspace(Polars.getAlpha(0),Polars.getAlpha(-1),n);
            wcl[jj,:] = np.interp(wa[jj,:],Polars.getAlpha(),Polars.getCl());
            wcd[jj,:] = np.interp(wa[jj,:],Polars.getAlpha(),Polars.getCd());
            wcm[jj,:] = np.interp(wa[jj,:],Polars.getAlpha(),Polars.getCm());
            wa[jj,:] *= m.pi/180.;
    else:
        Polars = deepcopy(wing.getPolar(0,0));
        if wing.getDFR(0):
            iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
            iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
        else:
            iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
            iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
        if not(np.any(iiClmin1)):
            iiClmin1 = m.floor(len(Polars.getAlpha())*0.25);
        if not(np.any(iiClmax1)):
            iiClmax1 = m.floor(len(Polars.getAlpha())*0.5);
        cl_alpha = np.ones(wing.getR(),dtype = float) * 180./m.pi*(Polars.getCl(iiClmax1)-Polars.getCl(iiClmin1))/ (Polars.getAlpha(iiClmax1)-Polars.getAlpha(iiClmin1));
        wa[0,:] = np.linspace(Polars.getAlpha(0),Polars.getAlpha(-1),n);
        wcl[0,:] = np.interp(wa[0,:],Polars.getAlpha(),Polars.getCl());
        wcd[0,:] = np.interp(wa[0,:],Polars.getAlpha(),Polars.getCd());
        wcm[0,:] = np.interp(wa[0,:],Polars.getAlpha(),Polars.getCm());
        wa[0,:] *= m.pi/180.;
        alphaMax = m.pi/180.*(Polars.getAlphaMax())*np.ones(wing.getR(),dtype = float);
        for ii in range(1,wing.getR()):
            wa[ii,:] = wa[0,:];
            wcl[ii,:] = wcl[0,:];
            wcd[ii,:] = wcd[0,:];
            wcm[ii,:] = wcm[0,:];
        al0 = m.pi/180*(Polars.getAl0())*np.ones(wing.getR(),dtype = float);
    return wa,wcl, wcd, wcm, al0, alphaMax, cl_alpha;



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
    w = np.linalg.solve(mat,cl);
    wcl = np.dot(w,pol);
    cd = np.interp(x,polar.getAlpha(),polar.getCd());
    w = np.linalg.solve(mat,cd);
    wcd = np.dot(w,pol);
    cm = np.interp(x,polar.getAlpha(),polar.getCm());
    w = np.linalg.solve(mat,cm);
    wcm = np.dot(w,pol);
    return wcl,wcd,wcm;

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
    mu = np.exp(-np.absolute(y-yStar)*3./dist);
    return mu;

def justLoad(fp,sk=0):
    try :
        dataPolar = np.loadtxt(fp,skiprows=sk);
        data = np.zeros([len(dataPolar[:,0]),5],dtype = float);
        if sk == 13:
            addDrag = 0.25 * (dataPolar[:,5]-dataPolar[:,6]);
            data[:,0:4] = dataPolar[:,[0,1,2,3]];
            data[:,4] = addDrag;
        else:
            data[:,0:4] = dataPolar[:,[0,1,2,4]];
    except ValueError:
        data = justLoad(fp,sk+1);
    except IOError:
        print('Document not found, verify the path and the existence of the file.');
        sys.exit();
    return data;

def loadPol(fp,sk=0):
    """ Function that loads the 2D polar data from the path"""
    dataPolar = justLoad(fp,0);

    polar = Polar2();

    # find clMax
    iiClmax = np.where(dataPolar[:,1] == np.amax(dataPolar[:,1]))[0][0];
    polar.setClMax(dataPolar[iiClmax,1]);
    polar.setAlphaMax(dataPolar[iiClmax,0]);
    dA = (dataPolar[-1,0]-dataPolar[0,0])/len(dataPolar[:,0]);
    nbA = int(max([0.,m.ceil((35.-dataPolar[-1,0])/dA)]));
    # find the zero lift angle of the polar file + enlarge "field" of alphas for polar curves, linear extrapolation
    lowAlpha = dataPolar[0,0] - (dataPolar[1,0]-dataPolar[0,0]) * np.array(range(10,0,-1));
    highAlpha = dataPolar[-1,0] + dA*np.array(range(1,nbA+1));
    lowCL = dataPolar[0,1] - (dataPolar[1,1]-dataPolar[0,1]) * np.array(range(10,0,-1));
    highCL = dataPolar[-1,1] * np.exp(-(highAlpha-dataPolar[-1,0])*1.5) + 1.2 * (1. - np.exp(-(highAlpha-dataPolar[-1,0])*1.5));
#    for i in range(10):
#        if highCL[i] < 0.7:
#            highCL[i] = 0.7;
    lowCD = dataPolar[0,2] - (dataPolar[1,2]-dataPolar[0,2]) * np.array(range(10,0,-1));
    highCD = dataPolar[-1,2] + 2. *np.sin(highAlpha*m.pi/180.) - 2. *m.sin(dataPolar[-1,0]*m.pi/180.);
    lowCM = dataPolar[0,3] - (dataPolar[1,3]-dataPolar[0,3]) * np.array(range(10,0,-1));
    highCM = dataPolar[-1,3] * np.ones(nbA);
   

    data1 = np.zeros([len(np.concatenate([lowAlpha,dataPolar[:,0],highAlpha])),4],dtype = float);
    data1[:,0] = np.concatenate([lowAlpha,dataPolar[:,0],highAlpha]);
    data1[:,1] = np.concatenate([lowCL,dataPolar[:,1],highCL]);
    data1[:,2] = np.concatenate([lowCD,dataPolar[:,2],highCD]);
    data1[:,3] = np.concatenate([lowCM,dataPolar[:,3],highCM]);

#    data1 = np.zeros([len(np.concatenate([lowAlpha,dataPolar[:,0]])),4],dtype = float);
#    data1[:,0] = np.concatenate([lowAlpha,dataPolar[:,0]]);
#    data1[:,1] = np.concatenate([lowCL,dataPolar[:,1]]);
#    data1[:,2] = np.concatenate([lowCD,dataPolar[:,2]]);
#    data1[:,3] = np.concatenate([lowCM,dataPolar[:,3]]);
    
    iiClmin = np.where(data1[:,1] == min(data1[:,1]))[0][-1];
    iiClmax = np.where(data1[:,1] == max(data1[:,1]))[0][0];
    al0_Polar= np.interp(0,data1[iiClmin:iiClmax+1,1],data1[iiClmin:iiClmax+1,0],left=999, right=999);

    while al0_Polar == 999:
        lowAlpha = data1[0,0] - (data1[1,0]-data1[0,0]) * np.array(range(5,0,-1));
        lowCL = data1[0,1] - (data1[1,1]-data1[0,1]) * np.array(range(5,0,-1));
        lowCD = data1[0,2] - (data1[1,2]-data1[0,2]) * np.array(range(5,0,-1));
        lowCM = data1[0,3] - (data1[1,3]-data1[0,3]) * np.array(range(5,0,-1));
        data1 = np.concatenate([np.vstack(np.concatenate([lowAlpha,data1[:,0]])),\
            np.vstack(np.concatenate([lowCL,data1[:,1]])),\
            np.vstack(np.concatenate([lowCD,data1[:,2]])),\
            np.vstack(np.concatenate([lowCM,data1[:,3]]))],1);
        iiClmin = np.where(data1[:,1] == min(data1[:,1]))[0][-1];
        iiClmax = np.where(data1[:,1] == max(data1[:,1]))[0][0];
        al0_Polar= np.interp(0,data1[iiClmin:iiClmax+1,1],data1[iiClmin:iiClmax+1,0],left=999, right=999);
    polar.setAl0(al0_Polar);

    # put convergence limits at -0.7 and 0.7 for resp. low angles and high angles
    lowAlpha = data1[0,0] -  np.array(range(10,0,-1));
    highAlpha = data1[-1,0] + 0.2*np.array([1,2,3,4,5,5.2,5.4,5.6,5.8,6],dtype = float);
    lowCL = np.flipud(np.concatenate([np.linspace(data1[0,1]+(-1.2-data1[0,1])/5,-1.2,4),-1.2*np.ones(6,dtype=float)]));
    highCL = (np.concatenate([np.linspace(data1[-1,1]+(1.2-data1[-1,1])/6,0.7,5),1.2*np.ones(5,dtype=float)]));
    lowCD = data1[0,2] - (data1[1,2]-data1[0,2]) *np.array(range(10,0,-1));
    highCD = data1[-1,2] + 0.2*(data1[-1,2]-data1[-2,2])*np.array([1,2,3,4,5,5.2,5.4,5.6,5.8,6],dtype = float);
    lowCM = data1[0,3] - (data1[1,3]-data1[0,3]) *np.array(range(10,0,-1));
    highCM = data1[-1,3] + 0.2*(data1[-1,3]-data1[-2,3])*np.array([1,2,3,4,5,5.2,5.4,5.6,5.8,6],dtype = float);
    a =np.concatenate([lowAlpha,data1[:,0],highAlpha]);
    b = np.concatenate([lowCL,data1[:,1],highCL]);
    c = np.concatenate([lowCD,data1[:,2],highCD]);
    d = np.concatenate([lowCM,data1[:,3],highCM]);
    polar.setData(a,b,c,d);
    polar.setFP(fp);
    return polar;

def vtailPolar(vtail):
    """
    
     Autor : Quentin Borlon
     This def provides the altered polar data required for the reste of
     the code based on the 2D polar files.
     
     INPUT:
       vtail.getR() : number of control station along the span +1
       vtail.polarIndex : array of index. vtail.polarIndex(i) indicate the index
       of the file in vtail.NLPol to use for the i-th station along the span
       
     OUTPUT:
       [wcl, wcd, wcm, wcda] : [(vtail.getR()) x 41] matrix.
           Array i gives the weight coefficients for 40 degree legendre 
           interpolation the lift/drag:piutching moment and additional drag 
           curves for the i-th panel along the span. Aktered data like 
           in Sivells
       al0 : (vtail.getR())-array giving the zero-lift-angle 
           of the panel along the span linear approximation from the linear part
           of the polar curves
       cl_alpha : vtail.getR()-array giving the lift curve slope 
           of the panel along the span"""

    n = 200;                   
    wa = np.zeros([vtail.getR(),n]);                                            
    wcl = np.zeros([vtail.getR(),n]);
    wcd = np.zeros([vtail.getR(),n]);
    wcm = np.zeros([vtail.getR(),n]);
    al0 = np.zeros(vtail.getR(),dtype = float);
    cl_alpha = np.zeros(vtail.getR(),dtype = float);
    alphaMax = np.zeros(vtail.getR(),dtype = float);
    if np.any(vtail.getPolD()):
        zDisc = vtail.getPolD();
        zPanel = vtail.getZP();
        mu = np.zeros([len(zDisc),vtail.getR()],dtype = float);
        muPos = np.zeros([len(zDisc),vtail.getR()],dtype = float);
        muNeg = np.zeros([len(zDisc),vtail.getR()],dtype = float);
        deltaclmaxStar = np.zeros(len(zDisc),dtype = float);
        # Apply the polar modification to ensure continuity of the clmax
        for ii in range(len(zDisc)):
            zStar1 = zDisc[ii];
            mStar = np.where(np.abs(vtail.getZ()-(zStar1-0.0001))<=0.0002)[0][0];
            # characteristic distance of influence of the discontinuity =
            # len of the flaps
            dist = vtail.getChordDist(mStar);
            distneg = dist * (1. - vtail.getCF());
            distpos = dist * (1. + vtail.getCF());
            muNeg[ii,:] = exponentialCoef(zStar1,zPanel,distneg);
            muPos[ii,:] = exponentialCoef(zStar1,zPanel,distpos);
            # Compute the initial delta in cl_max at the discontinuity
            deltaclmaxStar[ii] = vtail.getPolar(vtail.getPolI(mStar)).getClMax()-\
                vtail.getPolar(vtail.getPolI(mStar-1)).getClMax(); # le plus haut - le plus bas
            if deltaclmaxStar[ii] > 0.:
                mu[ii,:mStar] = 0.6*muPos[ii,:mStar];
                mu[ii,mStar:] = -0.4*muNeg[ii,mStar:];
            else:
                mu[ii,:mStar] = -0.4*muNeg[ii,:mStar];
                mu[ii,mStar:] = 0.6*muPos[ii,mStar:];
        deltaclmaxStar = np.absolute(deltaclmaxStar);
        # Polar altering like in Sivells Westrick NACA TN 2283
        for jj in range(vtail.getR()):
            Polars = deepcopy(vtail.getPolar(vtail.getPolI(jj)));
            clmaxAlt = Polars.getClMax() + np.dot(mu[:,jj],deltaclmaxStar);
            coef = clmaxAlt/Polars.getClMax();
            Polars.setCl(coef*Polars.getCl());
            Polars.setAlpha(coef*(Polars.getAlpha()-Polars.getAl0())+Polars.getAl0());
            Polars.setClMax(clmaxAlt);
            al0[jj] = m.pi/180.*(Polars.getAl0());
            iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
            iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
            if not(np.any(iiClmin1)):
                iiClmin1 = m.floor(len(Polars.getAlpha())*0.25);
            if not(np.any(iiClmax1)):
                iiClmax1 = m.floor(len(Polars.getAlpha())*0.5);
            cl_alpha[jj] = 180./m.pi*(Polars.getCl(iiClmax1)-Polars.getCl(iiClmin1)) / (Polars.getAlpha(iiClmax1)-Polars.getAlpha(iiClmin1));
            alphaMax[jj] = m.pi/180.*(Polars.getAlphaMax())*coef;
            wa[jj,:] = np.linspace(Polars.getAlpha(0),Polars.getAlpha(-1),n);
            wcl[jj,:] = np.interp(wa[jj,:],Polars.getAlpha(),Polars.getCl());
            wcd[jj,:] = np.interp(wa[jj,:],Polars.getAlpha(),Polars.getCd());
            wcm[jj,:] = np.interp(wa[jj,:],Polars.getAlpha(),Polars.getCm());
            wa[jj,:] *= m.pi/180.;
    else:
        Polars = deepcopy(vtail.getPolar(0));
        iiClmin1 = np.where(Polars.getAlpha() <= -2)[0][-1];
        iiClmax1 = np.where(Polars.getAlpha() >= 2)[0][0];
        if not(np.any(iiClmin1)):
            iiClmin1 = m.floor(len(Polars.getAlpha())*0.25);
        if not(np.any(iiClmax1)):
            iiClmax1 = m.floor(len(Polars.getAlpha())*0.5);
        cl_alpha = np.ones(vtail.getR(),dtype = float) * 180./m.pi*(Polars.getCl(iiClmax1)-Polars.getCl(iiClmin1))/ (Polars.getAlpha(iiClmax1)-Polars.getAlpha(iiClmin1));
        wa[0,:] = np.linspace(Polars.getAlpha(0),Polars.getAlpha(-1),n);
        wcl[0,:] = np.interp(wa[0,:],Polars.getAlpha(),Polars.getCl());
        wcd[0,:] = np.interp(wa[0,:],Polars.getAlpha(),Polars.getCd());
        wcm[0,:] = np.interp(wa[0,:],Polars.getAlpha(),Polars.getCm());
        wa[0,:] *= m.pi/180.;
        alphaMax = m.pi/180.*(Polars.getAlphaMax())*np.ones(vtail.getR(),dtype = float);
        for ii in range(1,vtail.getR()):
            wa[ii,:] = wa[0,:];
            wcl[ii,:] = wcl[0,:];
            wcd[ii,:] = wcd[0,:];
            wcm[ii,:] = wcm[0,:];
        al0 = m.pi/180*(Polars.getAl0())*np.ones(vtail.getR(),dtype = float);
    return wa,wcl, wcd, wcm, al0, alphaMax, cl_alpha;

def StalledPolar():
    Polars = loadPol('./PolarFiles/Stalled.txt');
    [wcl,wcd,wcm] = intrapolLegend(Polars,40);
    return wcl, wcd, wcm;
def AeroCoef(wa,wcl,wcd,wcm,alpha):
    """
    Autor : Quentin Borlon

    This function compute the interpollation of the lift and drag coefficient
    based on the correspondin weight using the legendre approximation series.

    INPUT:
        [wcl,wcd] : weight for the lift and drag curves.
        alpha : angle of attack at which the coefficients need to be computed

    OUTPUT:
        [cl,cd] : lift and drag coefficient thanks to the legendre 
        approximations of the curves."""
    n = len(alpha)
    cl = np.zeros(n);
    cd = np.zeros(n);
    cm = np.zeros(n);
    for ii in range(n):
        cl[ii] = np.interp(alpha[ii],wa[ii],wcl[ii]);
        cd[ii] = np.interp(alpha[ii],wa[ii],wcd[ii]);
        cm[ii] = np.interp(alpha[ii],wa[ii],wcm[ii]);
    return cl,cd,cm;