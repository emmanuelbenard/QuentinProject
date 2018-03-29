# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os
import AerodynamicPredictor as AP;
from datetime import date
import InputReader as IR;
import Aircraft_aero as AC;

if __name__ == '__main__':
    pFile = '/Users/Quentin/Documents/Toulouse/Cours/PIE/highlift/Python27_FAST/Inputs/AeroInputFile_Optim.txt';
    acI,flow,lin,resFile,DEI = IR.readParam(pFile);
    (filepath, filename) = os.path.split(resFile);
    if not(os.path.isdir(filepath)):
        os.makedirs(filepath);
    if DEI and acI.prop.bool:
        ac = AC.AircraftManager_DEI(acI,flow);
        CL,CD,CY,CM,Cl,Cn,CDi,CD0,Alphas = AP.main_DEI(ac,flow,lin);
    else:
        ac = AC.AircraftManager_SEI(acI,flow);
        CL,CD,CY,CM,Cl,Cn,CDi,CD0,Alphas = AP.main_SEI(ac,flow,lin);
    t= date.today();
    with open(resFile,'w') as files:
        files.write('    alpha        CL         CD         CDi        CD0        CY         CM         Cl         Cn        L/D   \n');
        files.write('   -------     ------     ------     ------     ------     ------     ------     ------     ------    ------- \n');
        for ii in range(len(CL)):
            if not np.isnan(CL[ii]):
                files.write('%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n' % (Alphas[ii], CL[ii], CD[ii],CDi[ii],CD0[ii],CY[ii],CM[ii], Cl[ii],Cn[ii],CL[ii]/CD[ii]));
    print 'CL(alpha)';
    plt.plot(Alphas,CL);
    plt.show();
    print 'CL(CD)';
    plt.plot(CD,CL);
    plt.plot(CDi,CL);
    plt.plot(CD0,CL);
    plt.show();
    print 'CL(CM)';
    plt.plot(CM,CL);
    plt.show();
    print 'Cl(alpha)';
    plt.plot(Alphas,Cl);
    plt.show();
    print 'Cn(alpha)';
    plt.plot(Alphas,Cn);
    plt.show();
    print 'CL/CD (alpha)';
    plt.plot(Alphas,CL/CD);
    plt.axis([0,20,0,15])
    plt.show();
    print max(CL/CD)