# -*- coding: utf-8 -*-
import numpy as np;
import Aircraft;
import Flow;
import Polar as p;
def readParam(pFile):
    """ This function read the input parameters and return them for the main function
    INPUT:
        pFile : the path and name of the input parameter file
    OUTPUT:
        ac : an object that contain the aircraft parameters
        flow : an object that containt the flow parameters
        resFile : the path of the results' file
    """
    ac = Aircraft.Aircraft();
    flow = Flow.Flow();
    iW = 0.;
    iT = 0.;
    lin = False;
    DEI = False;
    with  open(pFile,'r') as files:
        contenu = files.readlines();
    for line in contenu:
        mots = line.split();
        if len(mots)>2:
            if mots[0] == 'output':
                resFile = mots[2];
            elif mots[0] == 'EIS':
                if int(mots[2]) != 0:
                    ac.prop.bool = True;
            elif mots[0] == 'DEI':
                if int(mots[2]) != 0:
                    DEI = True;
            elif mots[0] == 'linear':
                if int(mots[2]) != 0:
                    lin = True;
            elif mots[0] == 'aMin':
                flow.setAMin(float(mots[2]));
            elif mots[0] == 'aMax':
                flow.setAMax(float(mots[2]));
            elif mots[0] == 'deltaA':
                flow.setDA(float(mots[2]));
            elif mots[0] == 'beta':
                flow.setBeta(float(mots[2]));
            elif mots[0] == 'Mach':
                flow.setMach(float(mots[2]));
            elif mots[0] == 'V0':
                flow.setV0(float(mots[2]));
            elif mots[0] == 'h':
                flow.setH(float(mots[2]));
            elif mots[0] == 'b':
                ac.wing.incSpan(float(mots[2]));
            elif mots[0] == 'cRoot':
                ac.wing.addChord(float(mots[2]));
            elif mots[0] == 'cTip':
                cTip = float(mots[2]);
            elif mots[0] == 'twistTip':
                tTip = float(mots[2]);
            elif mots[0] == 'cfToCLoc':
                ac.wing.setCF(float(mots[2]));
            elif mots[0] == 'iW':
                iW = float(mots[2]);
            elif mots[0] == 'htail':
                if float(mots[2]) != 0:
                    ac.htail.bool = True;
            elif mots[0] == 'vtail':
                if float(mots[2]) != 0:
                    ac.vtail.bool = True;
            elif mots[0] == 'fus':
                if float(mots[2]) != 0:
                    ac.fus.bool = True;
            elif mots[0] == 'refPoint':
                rp = np.zeros(3,dtype = float);
                rp[0] = float(mots[2]);
                rp[2] = float(mots[3]);
                ac.setRP(rp);
            elif mots[0] == 'yTaperDisc':
                for ii in range(2,len(mots)):
                    try:
                        y = float(mots[ii]);
                        ac.wing.addTapD(y);
                    except ValueError:
                        break;
            elif mots[0] == 'chord':
                for ii in range(2,len(mots)):
                    try:
                        c = float(mots[ii]);
                        ac.wing.addChord(c);
                    except ValueError:
                        break;
            elif mots[0] == 'ySweepDisc':
                for ii in range(2,len(mots)):
                    try:
                        y = float(mots[ii]);
                        ac.wing.addSweepD(y);
                    except ValueError:
                        break;
            elif mots[0] == 'sweep':
                sw = float(mots[2]);
                ac.wing.resetSweepC4(sw);
                for ii in range(3,len(mots)):
                    try:
                        sw = float(mots[ii]);
                        ac.wing.addSweepC4(sw);
                    except ValueError:
                        break;
            elif mots[0] == 'yTwistDisc':
                for ii in range(2,len(mots)):
                    try:
                        y = float(mots[ii]);
                        ac.wing.addTwistD(y);
                    except ValueError:
                        break;
            elif mots[0] == 'twist':
                ac.wing.resetTwist(iW);
                for ii in range(2,len(mots)):
                    try:
                        tw = float(mots[ii])+iW;
                        ac.wing.addTwist(tw);
                    except ValueError:
                        break;
            elif mots[0] == 'yDihDisc':
                for ii in range(2,len(mots)):
                    try:
                        y = float(mots[ii]);
                        ac.wing.addDihD(y);
                    except ValueError:
                        break;
            elif mots[0] == 'Dih':
                dih = float(mots[2]);
                ac.wing.resetDih(dih);
                for ii in range(3,len(mots)):
                    try:
                        dih = float(mots[ii]);
                        ac.wing.addDih(dih);
                    except ValueError:
                        break;
            elif mots[0] == 'yDiscAf':
                for ii in range(2,len(mots)):
                    try:
                        y = float(mots[ii]);
                        ac.wing.addAfD(y);
                    except ValueError:
                        break;
            elif mots[0] == 'airfoil':
                for ii in range(2,3+len(ac.wing.getAfD())):
                    try:
                        af = float(mots[ii]);
                        ac.wing.addAF(mots[ii]);
                    except ValueError:
                        ac.wing.addAF(mots[ii]);
            elif mots[0] == 'yDiscPol':
                for ii in range(2,len(mots)):
                    try:
                        y = float(mots[ii]);
                        ac.wing.addPolD(y);
                    except ValueError:
                        break;
            elif mots[0] == 'polarL':
                for ii in range(2,3+len(ac.wing.getPolD())):
                    polPath = mots[ii];
                    polar = p.loadPol(polPath);
                    ac.wing.addPolarL(polar);
            elif mots[0] == 'polarR':
                for ii in range(2,3+len(ac.wing.getPolD())):
                    polPath = mots[ii];
                    polar = p.loadPol(polPath);
                    ac.wing.addPolarR(polar);
            elif mots[0] == 'yDiscFlaps':
                for ii in range(2,len(mots)):
                    try:
                        y = float(mots[ii]);
                        ac.wing.addFlapsD(y);
                    except ValueError:
                        break;
            elif mots[0] == 'deltaL':
                for ii in range(2,len(mots)):
                    try:
                        delta = float(mots[ii]);
                        ac.wing.addDFL(delta);
                    except ValueError:
                        break;
            elif mots[0] == 'deltaR':
                for ii in range(2,len(mots)):
                    try:
                        delta = float(mots[ii]);
                        ac.wing.addDFR(delta);
                    except ValueError:
                        break;
            if ac.prop.bool:
                if mots[0] == 'nbE':
                    nbE = int(mots[2]);
                elif mots[0] == 'Tc':
                    Tc = np.zeros(nbE,dtype = float);
                    for ii in range(2,nbE+2):
                        try:
                            Tc[ii-2] = float(mots[ii]);
                        except ValueError:
                            break;
                    ac.prop.setTc(Tc);
                elif mots[0] == 'D':
                    D = np.zeros(nbE,dtype = float);
                    for ii in range(2,nbE+2):
                        try:
                            D[ii-2] = float(mots[ii]);
                        except ValueError:
                            break;
                    ac.prop.setD(D);
                elif mots[0] == 'Rhub':
                    Rhub = np.zeros(nbE,dtype = float);
                    for ii in range(2,nbE+2):
                        try:
                            Rhub[ii-2] = float(mots[ii]);
                        except ValueError:
                            break;
                    ac.prop.setRh(Rhub);
                elif mots[0] == 'Omega':
                    Omega = np.zeros(nbE,dtype = float);
                    for ii in range(2,nbE+2):
                        try:
                            Omega[ii-2] = float(mots[ii]);
                        except ValueError:
                            break;
                    ac.prop.setOmega(Omega);
                elif mots[0] == 'T':
                    T = np.zeros(nbE,dtype = float);
                    for ii in range(2,nbE+2):
                        try:
                            T[ii-2] = float(mots[ii]);
                        except ValueError:
                            break;
                    ac.prop.setT(T);
                elif mots[0] == 'yh':
                    yh = np.zeros(nbE,dtype = float);
                    for ii in range(2,nbE+2):
                        try:
                            y = float(mots[ii]);
                            yh[ii-2] = y;
                        except ValueError:
                            break;
                    ac.prop.setYp(yh);
                elif mots[0] == 'xh':
                    xh = np.zeros(nbE,dtype = float);
                    for ii in range(2,nbE+2):
                        try:
                            x = float(mots[ii]);
                            xh[ii-2] = x;
                        except ValueError:
                            break;
                    ac.prop.setXp(xh);
                elif mots[0] == 'zh':
                    zh = np.zeros(nbE,dtype = float);
                    for ii in range(2,nbE+2):
                        try:
                            z = float(mots[ii]);
                            zh[ii-2] = z;
                        except ValueError:
                            break;
                    ac.prop.setZp(zh);
                elif mots[0] == 'OWU':
                    OWU = np.empty(nbE,dtype = bool);
                    for ii in range(2,nbE+2):
                        try:
                            inv = float(mots[ii]);
                            if inv == 1:
                                OWU[ii-2] = True;
                            else:
                                OWU[ii-2] = False;
                        except ValueError:
                            break;
                    ac.prop.setOWU(OWU);
            if ac.htail.bool:
                if mots[0] == 'bT':
                    ac.htail.incSpan(float(mots[2]));
                elif mots[0] == 'cRootT':
                    ac.htail.addChord(float(mots[2]));
                elif mots[0] == 'cTipT':
                    cTipT = float(mots[2]);
                elif mots[0] == 'twistTipT':
                    tTipT = float(mots[2]);
                elif mots[0] == 'cfToCLocT':
                    ac.htail.setCF(float(mots[2]));
                elif mots[0] == 'iT':
                    iT = float(mots[2]);
                elif mots[0] == 'hDistT':
                    ac.htail.setHDist(float(mots[2]));
                elif mots[0] == 'vDistT':
                    ac.htail.setVDist(float(mots[2]));
                elif mots[0] == 'yTaperDiscT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.htail.addTapD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'chordT':
                    for ii in range(2,len(mots)):
                        try:
                            c = float(mots[ii]);
                            ac.htail.addChord(c);
                        except ValueError:
                            break;
                elif mots[0] == 'ySweepDiscT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.htail.addSweepD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'sweepT':
                    sw = float(mots[2]);
                    ac.htail.resetSweepC4(sw);
                    for ii in range(3,len(mots)):
                        try:
                            sw = float(mots[ii]);
                            ac.htail.addSweep(sw);
                        except ValueError:
                            break;
                elif mots[0] == 'yTwistDiscT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.htail.addTwistD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'twistT':
                    ac.htail.resetTwist(iT);
                    for ii in range(2,len(mots)):
                        try:
                            tw = float(mots[ii])+iT;
                            ac.htail.addTwist(tw);
                        except ValueError:
                            break;
                elif mots[0] == 'yDihDiscT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.htail.addDihD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'DihT':
                    dih = float(mots[2]);
                    ac.htail.resetDih(dih);
                    for ii in range(3,len(mots)):
                        try:
                            dih = float(mots[ii]);
                            ac.htail.addDih(dih);
                        except ValueError:
                            break;
                elif mots[0] == 'yDiscAfT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.htail.addAfD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'airfoilT':
                    for ii in range(2,3+len(ac.htail.getAfD())):
                        try:
                            af = float(mots[ii]);
                            ac.htail.addAF(mots[ii]);
                        except ValueError:
                            ac.htail.addAF(mots[ii]);
                elif mots[0] == 'yDiscPolT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.htail.addPolD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'polarTL':
                    for ii in range(2,3+len(ac.htail.getPolD())):
                        polPath = mots[ii];
                        polar = p.loadPol(polPath);
                        ac.htail.addPolarL(polar);
                elif mots[0] == 'polarTR':
                    for ii in range(2,3+len(ac.htail.getPolD())):
                        polPath = mots[ii];
                        polar = p.loadPol(polPath);
                        ac.htail.addPolarR(polar);
                elif mots[0] == 'yDiscFlapsT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.htail.addFlapsD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'deltaT':
                    for ii in range(2,len(mots)):
                        try:
                            delta = float(mots[ii]);
                            ac.htail.addDFL(delta);
                            ac.htail.addDFR(delta);
                        except ValueError:
                            break;
            if ac.vtail.bool:
                if mots[0] == 'bVT':
                    ac.vtail.incSpan(float(mots[2]));
                elif mots[0] == 'cRootVT':
                    ac.vtail.addChord(float(mots[2]));
                elif mots[0] == 'cTipVT':
                    cTipVT = float(mots[2]);
                elif mots[0] == 'cfToCLocVT':
                    ac.vtail.setCF(float(mots[2]));
                elif mots[0] == 'hDistVT':
                    ac.vtail.setHDist(float(mots[2]));
                elif mots[0] == 'vDistVT':
                    ac.vtail.setVDist(float(mots[2]));
                elif mots[0] == 'zTaperDiscVT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.vtail.addTapD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'chordVT':
                    for ii in range(2,len(mots)):
                        try:
                            c = float(mots[ii]);
                            ac.vtail.addChord(c);
                        except ValueError:
                            break;
                elif mots[0] == 'zSweepDiscVT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.vtail.addSweepD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'sweepVT':
                    sw = float(mots[2]);
                    ac.vtail.resetSweepC4(sw);
                    for ii in range(3,len(mots)):
                        try:
                            sw = float(mots[ii]);
                            ac.vtail.addSweep(sw);
                        except ValueError:
                            break;
                elif mots[0] == 'zDiscAfVT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.vtail.addAfD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'airfoilVT':
                    for ii in range(2,3+len(ac.vtail.getAfD())):
                        try:
                            af = float(mots[ii]);
                            ac.vtail.addAF(mots[ii]);
                        except ValueError:
                            ac.vtail.addAF(mots[ii]);
                elif mots[0] == 'zDiscPolVT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.vtail.addPolD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'polarVT':
                    for ii in range(2,3+len(ac.vtail.getPolD())):
                        polPath = mots[ii];
                        polar = p.loadPol(polPath);
                        ac.vtail.addPolar(polar);
                elif mots[0] == 'zDiscFlapsVT':
                    for ii in range(2,len(mots)):
                        try:
                            y = float(mots[ii]);
                            ac.vtail.addFlapsD(y);
                        except ValueError:
                            break;
                elif mots[0] == 'deltaVT':
                    for ii in range(2,len(mots)):
                        try:
                            delta = float(mots[ii]);
                            ac.vtail.addDF(delta);
                        except ValueError:
                            break;
            if ac.fus.bool:
                if mots[0] == 'nL':
                    ac.fus.setNL(float(mots[2]));
                elif mots[0] == 'cL':
                    ac.fus.setCL(float(mots[2]));
                elif mots[0] == 'bL':
                    ac.fus.setBL(float(mots[2]));
                elif mots[0] == 'cD':
                    ac.fus.setCD(float(mots[2]));
                elif mots[0] == 'bD':
                    ac.fus.setBD(float(mots[2]));
                elif mots[0] == 'hDistF':
                    ac.fus.setHD(float(mots[2]));
                elif mots[0] == 'vDistF':
                    ac.fus.setVD(float(mots[2]));
    ac.wing.addChord(cTip);
    ac.wing.addTwist(tTip+iW);
    if ac.htail.bool:
        ac.htail.addChord(cTipT);
        ac.htail.addTwist(tTipT+iT);
    if ac.vtail.bool:
        ac.vtail.addChord(cTipVT);
    return ac,flow,lin,resFile,DEI;