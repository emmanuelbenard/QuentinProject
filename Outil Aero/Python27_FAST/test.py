# -*- coding: utf-8 -*-
import InputReader as IR;
import Aircraft as AC;
ac,flow,lin,resFile = IR.readParam('./Inputs/Exemple_INPUTS_PIE.txt');
ac2 = AC.AircraftManager(ac,flow);