close all
clear all
clc
%% amplitude influence

% choice the coefficient
coefficient='Cl';
% Amplitude of test
Amp_list =[0.25 0.5 1 2 3 4 15];
run('amplitude_influence.m')

%% Confidence interval (only for non-null coefficients)

% choice the coefficient
coefficient='Cn';
run('Confidence_interval')


%% model validation
% two simulation for the coeffcient
coefficient='Cn';
% Simulation 1 for estimation of the coefficient
% elevator
Amp_elevator1 = 1;
elevator1 =4;
% rudder
Amp_rudder1 =3;
rudder1=1;
% aileron
Amp_aileron1 =6 ;
aileron1=1;
% Simulation 2 for the validation by camparison with datas
% elevator
Amp_elevator2= 1;
elevator2 =4;
% rudder
Amp_rudder2 = 8;
rudder2 = 2;
% aileron
Amp_aileron2 = -3;
aileron2 = 3;

run('model_validation')
