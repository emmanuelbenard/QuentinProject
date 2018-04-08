% Datas script

%% Inputs

% elevator
% start time of the deflection (elevator)
Te=4; 

% amplitude of the deflection
Ae = 1;

test_elevator = 4;
% test_elevator = 1 => input11
% test_elevator = 2 => input211
% test_elevator = 3 => input3211
% test_elevator = 4 => no action on the elevator

% rudder
% start time of the deflection (rudder)
Tr=2; 

% amplitude of the deflection
Ar = 1;

test_rudder = 1;
% test_rudder = 1 => input11
% test_rudder = 2 => input211
% test_rudder = 3 => input3211
% test_rudder = 4 => no action on the elevator

% aileron

% start time of the deflection (aileron)
Ta=6; 

% amplitude of the deflection
Aa = 1;

test_aileron = 1;
% test_aileron = 1 => input11
% test_aileron = 2 => input211
% test_aileron = 3 => input3211
% test_aileron = 4 => no action on the elevator


%% Noise 

% STATE NOISE

Noise_power_Accelerometer = 1*10^-9;
Offset_Accelerometer=57*10^-6;

Noise_power_gyro = 1*10^-7;
Offset_gyro=0.15/180*pi;

Noise_power_velocity_angles = 1*10^-10;
Offset_velocity_angles=0.1/180*pi;

Noise_power_airspeed = 1*10^-4;

% external disturbance

% wind
K_wind = [0 0 0 0 0 0]';



%% inertial coefficients
Ix = 24675560;
Iy = 44876980;  
Iz = 67383260;  
Jxy = 0;  
Jxz = 1315126;
Jyz = 0;
% wing reference area
S = 510.95;
% reference wing span
b = 59.74;
% aircraft mass
m = 288772;
V = 150;
cbar= 8.32;

%% True coefficient

% true lift coefficient
CL_0_true = 0.21;
CL_alpha_true = 4.4;
CL_q_true = 6.6;
CL_deltae_true = 0.32;
CLih_true = 0.7;
CL_true = [CL_0_true CL_alpha_true CL_q_true CL_deltae_true CLih_true];

% true drag coefficient

Cd_0_true = 0.0164;
Cd_alpha_true = 0.2;
Cd_q_true = 0;
Cd_deltae_true = 0;
Cdih_true = 0;
Cd_true = [Cd_0_true Cd_alpha_true Cd_q_true Cd_deltae_true Cdih_true];

% true Cy coefficient
Cy_0_true = 0;
Cy_b_true = -0.9;
Cy_p_true = 0;
Cy_r_true = 0;
Cy_deltaa_true = 0;
Cy_deltar_true = 0.12;
Cy_true = [Cy_b_true Cy_deltar_true];

% true Y moment
Cm_0_true = 0;
Cm_alpha_true = -1;
Cm_q_true = -20.5;
Cm_deltae_true = -1.3;
Cm_ih_true = -2.7;
Cm_true = [Cm_0_true Cm_alpha_true Cm_q_true Cm_deltae_true Cm_ih_true];

% true X moment
Cl_0_true = 0;
Cl_b_true = -0.16;
Cl_p_true = -0.34;
Cl_r_true = 0.13;
Cl_deltaa_true = -0.013;
Cl_deltar_true =  0.008;
Cl_true = [Cl_b_true Cl_p_true Cl_r_true Cl_deltaa_true Cl_deltar_true];

% true Z moment
Cn_0_true =0 ;
Cn_b_true =0.16 ;
Cn_p_true =-0.026  ;
Cn_r_true = -0.28;
Cn_deltaa_true =-0.0018 ;
Cn_deltar_true= -0.1;
Cn_true = [Cn_b_true Cn_p_true Cn_r_true Cn_deltaa_true Cn_deltar_true];
