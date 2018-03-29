# -*-coding:Latin-1 -*
import math as m
import numpy as np
import numpy.matlib
import numpy.linalg
import copy
import matplotlib.pyplot as plt
import utilitaire as u
from Wing import Wing
def wingManager(wingI):
    """ Creation of the wing object under the right format
     Authors : Quentin Borlon
     Date : 5 mai 2017
    
     These function handels the different discontinuities and manages the new
     y-coordinates of the wing.
     The output wing now has the next additional fields:
    
    INPUT:
        wingI : a structral object with as fields:
    
           b : span
           chord : vertical array with the chord at the root (1) any
               discontinuity of taper ratio (2:end-1) and at the tip (end);
           flapsDiscY : vertical array with the spanwise coordinate of the flaps
               discontinuity sections
           afDiscY : vertical array with the spanwise coordinate of the
               airfoil discontinuity sections
           polarDiscY : vertical array with the spanwise coordinate of the
               polar discontinuity sections (2D curves)
           taperDiscY : vertical array with the spanwise coordinate of the
               taper ratio discontinuity sections
           sweepDiscY : vertical array with the spanwise coordinate of the
               sweep angle discontinuity sections
           dihDiscY : vertical array with the spanwise coordinate of the
               dihedral angle discontinuity sections
           twistDiscY : vertical array with the spanwise coordinate of the
               twist rate discontinuity sections
           cfToCLoc : the ratio between the flaps length and the local chord
               length
           polar : a cell-array with for each cell a "polar" object with field
               data (the alpha-cl-cd data), clmax, alphal0 (zero-lift angle),
               and clAlpha : lift curve slope in linear domain.
           airfoil : a cell-array with each cell gives the airfoil naca number
               representation, cell 1 correspond to first panel after root.
           deltaFlaps : vertical array with the different flaps defection
               along the span (deg)
           r : number of spanwise panel along the wing;
           m : number of chordwise panel along the airfoil;
           z : vertical array that gives the spanwise location for a
               discretisation of 100 station along the span
           vix : vertical array with wing.vix(i) is the tangencial engine 
               induced velocity on the station located at wing.z(i)
           vitheta : vertical array with wing.vitheta(i) is the normal engine 
               induced velocity on the station located at wing.z(i)
           V0 : flight velocity
    
    OUTPUT:
        wing : a structral object with additional fields:
    
           chord : vertical array with the chord at the root (1) any
               discontinuity of taper ratio (2:end-1) and at the tip (end);
           sweep : vertical array with wing.sweep(i) is the sweep angle of 
               the panel from wing.y(i) to wing.y(i+1) (rad)
           dih : vertical array with wing.dih(i) is the dihedral angle of 
               the panel from wing.y(i) to wing.y(i+1) (rad)
           twist : vertical array with wing.twist(i) is the twist angle of 
               the section at wing.y(i) (rad)
           deltasFlaps : vertical array with wing.deltasFlaps(i) is the 
               flaps defection  of the panel from wing.y(i) to wing.y(i+1)
               (deg)
           Mach : flight mach number
           vix : vertical array with wing.vix(i) is the mean (surface integrated
               normalized by the surface of the panel) of the tangencial engine 
               induced velocity on the panel from wing.y(i) to wing.y(i+1)
           vitheta : vertical array with wing.vix(i) is the mean (surface 
               integrated normalized by the surface of the panel) of the 
               normal engine induced velocity on the panel 
               from wing.y(i) to wing.y(i+1)
           V0 : flight velocity
           cFlaps_cLoc : vertical array with wing.cFlaps_cLocs(i) is the 
               local flaps to chord ratio
           theta : the angle (0 -> m.pi) spanwise location of the limits of the panels
           eta : the normalized (-1 -> 1) spanwise location of the limits of the panels
           y : the  spanwise location of (-b/2 -> b/2) the limits of the panels
           discY : vertical array of the complete set of the spanwise location
               of the discontinuity
           discTheta : same as above but in angle coordinates;
           airfoilIndex : vertical array with wing.airfoilIndex(i) is the index of 
               the airfoil (wing.airfoil) to use for the section at wing.y(i)
           polarIndex : vertical array with wing.polarIndex(i) is the index of 
               the airfoil (wing.polar) to use for the section at wing.y(i)
           chordDistrib : vertical array with wing.chordDistrib(i) is the chord length of 
               the section at wing.y(i)
           taper : vertical array with wing.taper(i) is the taper ratio at
               the section at wing.y(i) :  wing.chordDistrib(i+1)/ wing.chordDistrib(i)"""

    wing = Wing();
    wing.setiW(wingI.getiW());
    if wingI.getTailBool():
        wing.setTailBool(1);
        wing.setFusDiam(wingI.getFusDiam());
        wing.setFusLength(wingI.getFusLength());
        wing.tail = wingI.tail;
        wing.tail.setD(wingI.getD());
        wing.tail.setRh(wingI.getRh());
    wing.incSpan(wingI.getSpan());
    wing.setCF(wingI.getCF());
    wing.addPolar(wingI.getPolar());
    wing.addAF(wingI.getAF());
    wing.addDF(wingI.getDF());
    wing.setM(wingI.getM());
    wing.setV0(wingI.getV0());
    wing.setMach(wingI.getMach());
    # wing spanwise stations placement
    wing.setTheta(np.array(range(wingI.getR(),-1,-1),dtype = float)*m.pi/(wingI.getR()));                               # theta angle of the border of each spanwise station firs lw then rw                        
    wing.setEta(np.cos(wing.getTheta()));                                            # y-normalised coordinates of borders of the spanwise stations

    wing.setY(np.round(wing.eta*wing.b/2,4));                                  # y coordinates of the spanwise stations

    # Asuming the wing is sym, the discontinuities are also present on the lhw
    # (y < 0)
    wing.setTapD(np.sort(np.concatenate([-wingI.getTapD(),wingI.getTapD()])));
    wing.setDihD(np.sort(np.concatenate([-wingI.getDihD(),wingI.getDihD()])));
    wing.setSweepD(np.sort(np.concatenate([-wingI.getSweepD(),wingI.getSweepD()])));
    wing.setTwistD(np.sort(np.concatenate([-wingI.getTwistD(),wingI.getTwistD()])));
    wing.setAfD(np.sort(np.concatenate([-wingI.getAfD(),wingI.getAfD()])));
    wing.setFlapsD(np.sort(np.concatenate([-wingI.getFlapsD(),wingI.getFlapsD()])));
    wing.setPolD(np.sort(np.concatenate([-wingI.getPolD(),wingI.getPolD()])));
    discY = np.unique(np.concatenate([wing.getTapD(),wing.getDihD(),wing.getSweepD(),\
        wing.getTwistD(),wing.getAfD(),wing.getFlapsD(),wing.getPolD()]));
    if wing.getTailBool():
        discY = np.unique(np.concatenate([wing.getTapD(),wing.getDihD(),wing.getSweepD(),\
            wing.getTwistD(),wing.getAfD(),wing.getFlapsD(),wing.getPolD(),[-wing.getFusDiam()/2,wing.getFusDiam()/2]]));
    discY = np.round(discY,4);
    if wingI.getEIS():
        engineSection = np.zeros(4*len(wingI.getYh()),dtype = float);
        for i in range(len(wingI.getZh())):
            if wingI.yh[i] + wingI.getRh(i) <= wingI.getSpan()/2:
                engineSection[4*i] = wingI.yh[i] + wingI.getRh(i);
            if wingI.yh[i] + wingI.getD(i)/2 <= wingI.getSpan()/2:
                engineSection[4*i+1] = wingI.yh[i] + wingI.getD(i)/2;
            if wingI.yh[i] - wingI.getRh(i) >= -wingI.getSpan()/2:
                engineSection[4*i+3] = wingI.yh[i] - wingI.getRh(i);
            if wingI.yh[i] - wingI.getRh(i)/2 >= -wingI.getSpan()/2:
                engineSection[4*i+2] = wingI.yh[i] - wingI.getD(i)/2;
    else:
        engineSection = np.empty(0,dtype=float);
    Y = np.concatenate([wing.y,discY,engineSection]);
    tol = 1e-4;
    ii_unique = np.unique(Y,True)[1];
    Y = Y[ii_unique];
    AfOrFlaps = np.concatenate([wing.getAfD(),wing.getFlapsD()]);
    for ii in range(2,len(Y)-2):
        jj = ii-1;
        if Y[ii] in AfOrFlaps:
            Y[ii-1] = (Y[ii]+Y[ii-2])*0.5;
        if Y[jj] in  AfOrFlaps:
            Y[jj+1] = (Y[jj]+Y[jj+2])*0.5;

    wing.setY(Y);
    wing.setEta(wing.getY()*2/wingI.getSpan());
    wing.setTheta(np.round(np.arccos(np.round(wing.getEta(),4)),4));
    wing.setR(len(Y)-1);
    wing.setDiscTheta(np.round(np.arccos(discY/(wingI.getSpan()*0.5)),4));   
    wing.setDiscThetaPol(np.round(np.arccos(wing.getPolD()/(wingI.getSpan()*0.5)),4));     



    # the value is the value for the panel to come after, from the wing root to
    # the wing tip.

    # flaps discontinuity
    flapsZone = np.zeros(wing.getR()+1,dtype = int);
    df = wingI.getDF();
    if np.any(df):
        for i in range(len(wing.getFlapsD())/2,len(wing.getFlapsD())):
            indice = np.array([wing.getY(ii) >= wing.getFlapsD(i) for ii in range(wing.getR()/2, wing.getR()+1)],dtype = bool);
            incr = np.concatenate([np.flipud(indice[1:]),indice]);
            flapsZone += incr;
    wing.addDsF(df[flapsZone]);

    # airfoil discontinuity

    afZone = np.zeros(wing.getR()+1,dtype = int);
    if np.any(wing.getAfD()):
        for i in range(len(wing.getAfD())/2,len(wing.getAfD())):
            indice = np.array([wing.getY(ii) >= wing.getAfD(i) for ii in range(wing.getR()/2, wing.getR()+1)],dtype = bool);
            incr = np.concatenate([np.flipud(indice[1:]),indice]);
            afZone += incr;
    wing.addAFI(afZone);

    """ Normalement pas nécessaire
    for ii =1:length(wing.airfoil)
        wing.airfoil{ii} = str2double(wing.airfoil{ii});
    end"""
        


    # polar file allocation
    polZone = np.zeros(wing.getR()+1,dtype = int);
    if np.any(wing.getPolD()):
        for i in range(len(wing.getPolD())/2,len(wing.getPolD())):
            indice = np.array([wing.getY(ii) >= wing.getPolD(i) for ii in range(wing.getR()/2, wing.getR()+1)],dtype = bool);
            incr = np.concatenate([np.flipud(indice[1:]),indice]);
            polZone += incr;
    wing.addPolI(polZone);
    
    # assume symmetric taper distribution
    taper =  wing.getTapD();
    taper = taper[len(taper)/2:];
    if np.any(taper):
        x = np.round(np.concatenate([[0],taper,[wing.getSpan()/2]]),4);
    else:
        x = np.round(np.array([0.,wing.getSpan()/2]),4);
    wing.addChordDist(numpy.interp(np.absolute(wing.getY()), x, wingI.getChord()));
    taper = wing.getChordDist();
    taper = taper[wing.getR()/2+1:]/taper[wing.getR()/2:-1];
    taper = np.concatenate([np.flipud(taper[1:]),taper]);
    wing.addTaper(taper);

    # assume symmetric sweep distribution
    sweepZone = np.zeros(wing.getR()+1,dtype = int);
    sweep = wingI.getSweep()*m.pi/180;
    if np.any(wing.getSweepD()):
        for i in range(len(wing.getSweepD())/2,len(wing.getSweepD())):
            indice = np.array([wing.getY(ii) >= np.round(wing.getSweepD(i),4) for ii in range(wing.getR()/2, wing.getR()+1)],dtype = bool);
            incr = np.concatenate([np.flipud(indice[1:]),indice]);
            sweepZone += incr;
    wing.addSweep(sweep[sweepZone]);

    sweep = wing.getSweep();
    sweepC_4 = copy.deepcopy(sweep[wing.getR()/2:]);

    for i in range(wing.getR()/2):
        j = wing.getR()/2+i;
        d = wing.getY(j+1)-wing.getY(j);
        sweepC_4[i]= m.atan(m.tan(sweep[j])+0.25*(wing.chordDistrib[j+1]-wing.chordDistrib[j])/d);

    sweepC_4[-1] = sweepC_4[-2];
    sweepC_4 = np.round(sweepC_4,4);
    wing.addSweepC4(np.concatenate([np.flipud(sweepC_4[1:]),sweepC_4]));

    
    # assume symmetric dihedral distribution
    dihZone = np.zeros(wing.getR()+1,dtype = int);
    dih = wingI.getDih()*m.pi/180;
    if np.any(wing.getDihD()):
        for i in range(len(wing.getDihD())/2,len(wing.getDihD())):
            indice = np.array([wing.getY(ii) >= wing.getDihD(i) for ii in range(wing.getR()/2, wing.getR()+1)],dtype = bool);
            incr = np.concatenate([np.flipud(indice[1:]),indice]);
            dihZone += incr;
    wing.addDih(dih[dihZone]);             # local section dih

    # assume symmetric twist distribution and linear evolution btw sections
    if np.any(wing.getTwistD()):
        x = np.round(np.concatenate([[0],wingI.getTwistD(),[wingI.getSpan()/2]]),4);
    else:
        x = np.round(np.array([0.,wingI.getSpan()/2]),4);
    twist = wingI.getTwist()*m.pi/180;
    wing.addTwist(numpy.interp(np.absolute(wing.getY()), x, twist));
    
    # Engine velocties
    chord = numpy.interp(wingI.getY(), wing.getY(), wing.getChordDist());
    vix = np.zeros(wing.getR(),dtype = float);
    vitheta = np.zeros(wing.getR(),dtype = float);
    wingI.setY(np.flipud(np.round(wingI.getY(),4)));
    for jj in range(wing.getR()):
        debut = np.where(wingI.getY() <= wing.getY(jj))[0][-1];
        fin = np.where(wingI.getY() >= wing.getY(jj+1))[0][0]+1;
        vix[jj] = np.trapz(wingI.getViX(range(debut,fin))*chord[debut:fin],x=wingI.getY(range(debut,fin)))/((wingI.getY(fin-1)-wingI.getY(debut))*(chord[debut]+chord[fin-1])*0.5);
        vitheta[jj] = np.trapz(wingI.getViTheta(range(debut,fin))*chord[debut:fin],wingI.getY(range(debut,fin)))/((wingI.getY(fin-1)-wingI.getY(debut))*(chord[debut]+chord[fin-1])*0.5);
    if wing.getTailBool():
        wing.tail.vix = copy.deepcopy(wingI.vix);
        wing.tail.vitheta = copy.deepcopy(wingI.vitheta)*0.2;
        wing.tail.setY(np.flipud(wingI.getY()));
    wing.setViX(vix*m.cos(wing.getiW()*m.pi/180) - vitheta*m.sin(wing.iW*m.pi/180));                                     
    wing.setViTheta(vitheta*m.cos(wing.getiW()*m.pi/180) + vix*m.sin(wing.iW*m.pi/180));
    
    return wing;
