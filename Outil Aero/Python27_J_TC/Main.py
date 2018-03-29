# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import os
import aeroPred_TEST as AP;
from datetime import date
import InputReader as IR;
import Aircraft as AC;

if __name__ == '__main__':
    pFile = '/Users/Quentin/Documents/Toulouse/Cours/PIE/highlift/Python27_J_TC/Inputs/tn4448.txt';
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
        files.write('    alpha        CL         CD         CDi        CD0        CY         CM         Cl         Cn        L/D  \n');
        files.write('   -------     ------     ------     ------     ------     ------     ------     ------     ------    ------ \n');
        for ii in range(len(CL)):
            if not np.isnan(CL[ii]):
                files.write('%10.2f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n' % (Alphas[ii], CL[ii], CD[ii],CDi[ii],CD0[ii],CY[ii],CM[ii], Cl[ii],Cn[ii],CL[ii]/CD[ii]));
    a = np.array([-8.3613, -4.2042, 0.058605,4.4243,8.7895,12.946,17.415,19.598,21.575,23.551,25.840,27.923],dtype = float);
    cl = np.array([-0.65621, -0.19877, 0.18694, 0.63439, 1.1024, 1.5701, 2.0998, 2.3081, 2.3827, 2.5497, 2.6863, 2.6994],dtype = float);
    print 'CL(alpha)';
    plt.plot(Alphas,CL);
    plt.scatter(a,cl)
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
    print 'CL/CD (alpha), max = ', max(CL/CD);
    plt.plot(Alphas,CL/CD);
    plt.show();