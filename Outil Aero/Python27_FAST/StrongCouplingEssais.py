#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 11:52:19 2018

@author: Quentin
"""

while(max(abs(deltaAlpha))>tol and nbIter < nIterMAX):
                RHS = -(normal[0,:]* np.dot(select2,velTriX) + normal[1,:] * np.dot(select2,velTriY) + normal[2,:] * np.dot(select2,velTriZ));
                b = np.concatenate([RHS,deltaAlpha]);
                F = np.dot(SA,X) - b + np.dot(artVisc,b);
                URC = - normal[0,:]*np.dot(select2,velTriZ) + normal[2,:]*np.dot(select2,velTriX);
                for kk in range(nbPan):
                    J[kk*mC:(kk+1)*mC,nA+kk] = URC[kk*mC:(kk+1)*mC];
                    clM = np.interp(alphaLoc[kk] + 0.01,waPanel[kk],wclPanel[kk]);
                    clm = np.interp(alphaLoc[kk] - 0.01,waPanel[kk],wclPanel[kk]);
                    claLoc[kk] = (clM-clm)/(0.02);
                    VzPrimeij = (m.cos(dih[kk]) * m.sin(tp[ii]) - m.sin(dih[kk]) * m.sin(sweepC4[kk])) * Vx[kk,:] + m.cos(dih[kk]) * m.cos(tp[kk]) * Vz[kk,:];
                    VxPrimeij = m.cos(dih[kk]) * m.cos(tp[kk]) * Vx[kk,:] - m.sin(tp[kk]) * Vz[kk,:];
                    dClVdGamm = claLoc[kk]/(Norme[kk]*Norme[kk]) * (VzPrimeij * velTriXVisc[kk] - velTriZVisc[kk] * VxPrimeij);
                    dClIndGamm = select[kk,:] * 2./(c[kk]*norme[kk]);
                    J[nA+kk,:nA] = -1./(cl_alphaPanel[kk])*(dClVdGamm - dClIndGamm);
#                LRC = 1.+( claLoc)/cl_alphaPanel;
#                for i in range(nbPan):
#                    J[nA+i,nA+i] =  LRC[i];
                X1 = X - np.dot(np.linalg.inv(J),F);
                clw = 2. * np.dot(select,gamma)/(norme*c);
                Vix = np.dot(Vx,gamma);
                Viy = np.dot(Vy,gamma);
                Viz = np.dot(Vz,gamma);
#                velTriXVisc,velTriYVisc,velTriZVisc = u.localVelTri(VxA + Vix[:nbPan] * m.cos(Alphas[ii]) - Viz[:nbPan] * m.sin(Alphas[ii]),VyA + Viy[:nbPan],VzA + Viz[:nbPan] * m.cos(Alphas[ii]) + Vix[:nbPan] * m.sin(Alphas[ii]),tp,sweepC4,dih);
                velTriXVisc,velTriYVisc,velTriZVisc = u.localVelTri(VxA + Vix[:nbPan],VyA + Viy[:nbPan],VzA + Viz[:nbPan],tp,sweepC4,dih);
                alphaLoc = np.arctan2(velTriZVisc,velTriXVisc);
#                alphaLoc = alpha0 + clw/cl_alphaPanel + alphal0PanelVLM - (AOA_VLM0 + X[nA:]) - tp;
                [clVisc,cd,cm] = p.AeroCoef(waPanel,wclPanel,wcdPanel,wcmPanel,alphaLoc);
                X = X1
                gamma = X[:nA];
                deltaAlpha = (clVisc-clw) / cl_alphaPanel;
                AOA_VLM = AOA_VLM0 + X[nA:];
                velTriX = Norme * np.cos(AOA_VLM);
                velTriZ = Norme * np.sin(AOA_VLM);
                nbIter += 1
                plt.plot(y,clw)
                plt.plot(y,clVisc)
                plt.show();
            plt.plot(y,(clw/cl_alphaPanel + alphal0PanelVLM - (AOA_VLM+tp) + alpha0)*180./m.pi)
            plt.plot(y,180./m.pi*alphaLoc)
            plt.show();
            plt.plot(y,180./m.pi*deltaAlpha)
            plt.show();
            return