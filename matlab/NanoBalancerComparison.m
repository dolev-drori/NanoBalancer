% The following code is intended to compare the findings 
% of our mathematical model of the robot with those of a 
% paper on a similar robot. The goal is to understand if 
% they are of the same order of magnitude and to identify 
% any physical parameters that may be incorrect due to our 
% measurements.

%%%% our phisical parameters: %%%%
g = 9.8; % Gravity [m/s^2]
b_f = 0; % Friction coefficient
m_b = 0.3064; % Body mass [kg]
l_b = 0.08; % The distance between the mass centers of the body and wheels [m]
I_b = 0.00261; % Body inertia torque [kgm^2]
m_w = 0.0043; % Wheel mass [kg]
l_w = 0.02; % Wheel radius [m]
I_w = 0.00000086; % Wheel inertia torque [kgm^2]
R_m = 20.83; % Electrical resistance of motor [ohm]
l_m = 0; % ?
b_m = 0; % ?
K_e = 0.108; % Motor Electrical constant [s/r]
K_t = 0.108; % Motor torque constant [Nm/A]

%state space equetions
gamma_11 = (I_w)/(l_w) + l_w * m_b + l_w * m_w;
gamma_12 = m_b * l_b * l_w;
alpha_11 = 0;
alpha_12 = -(((K_e * K_t) / R_m) + b_f) / l_w;
alpha_13 = 0;
alpha_14 = ((K_e * K_t) / R_m) + b_f;
beta_11 = K_t/R_m;

gamma_21 = m_b * l_b;
gamma_22 = I_b + (m_b * l_b^2);
alpha_21 = 0;
alpha_22 = (((K_e * K_t) / R_m) + b_f) / l_w;
alpha_23 = m_b * l_b * g;
alpha_24 = - (((K_e * K_t) / R_m) +b_f);
beta_21 = -K_t / R_m;

% with disturbance:
beta_12 = l_w;
beta_21 = l_b;

delta = gamma_11*gamma_22 - gamma_12*gamma_21;

a_22 = (gamma_22*alpha_12 - gamma_12*alpha_22) / delta;
% alpha_13 =0: a_23 = (gamma_22*alpha_13 - gamma_12*alpha_23) / delta;
a_23 = (- gamma_12*alpha_23) / delta;
a_24 = (gamma_22*alpha_14 - gamma_12*alpha_24) / delta;
a_42 = (-gamma_21*alpha_12 + gamma_11*alpha_22) / delta;
% alpha_13 =0: a_43 = (-gamma_21*alpha_13 + gamma_11*alpha_23) / delta;
a_43 = (gamma_11*alpha_23) / delta;
a_44 = (-gamma_21*alpha_14 + gamma_11*alpha_24) / delta;

b_21 = (gamma_22*beta_11 - gamma_12*beta_21)/delta;
b_41 = (-gamma_21*beta_11 + gamma_11*beta_21)/delta;

%%% thires phisical parameters: %%%
g_d = 9.8; % Gravity [m/s^2]
b_f_d = 0; % Friction coefficient
m_b_d = 0.381; % Body mass [kg]
l_b_d = 0.112; % The distance between the mass centers of the body and wheels [m]
I_b_d = 0.00616; % Body inertia torque [kgm^2]
m_w_d = 0.036; % Wheel mass [kg]
l_w_d = 0.021; % Wheel radius [m]
I_w_d = 0.00000746; % Wheel inertia torque [kgm^2]
R_m_d = 4.4; % Electrical resistance of motor [ohm]
l_m_d = 0; % ?
b_m_d = 0; % ?
K_e_d = 0.444; % Motor Electrical constant [s/r]
K_t_d = 0.470; % Motor torque constant [Nm/A]

gamma_11_d = (I_w_d)/(l_w_d) + l_w_d * m_b_d + l_w_d * m_w_d;
gamma_12_d = m_b_d * l_b_d * l_w_d;
alpha_11_d = 0;
alpha_12_d = -(((K_e_d * K_t_d) / R_m_d) + b_f_d) / l_w_d;
alpha_13_d = 0;
alpha_14_d = ((K_e_d * K_t_d) / R_m_d) + b_f_d;
beta_11_d = K_t/R_m;

gamma_21_d = m_b_d * l_b_d;
gamma_22_d = I_b_d + (m_b_d * l_b_d^2);
alpha_21_d = 0;
alpha_22_d = (((K_e_d * K_t_d) / R_m_d) + b_f_d) / l_w_d;
alpha_23_d = m_b_d * l_b_d * g_d;
alpha_24_d = - (((K_e_d * K_t_d) / R_m_d) +b_f_d);
beta_21_d = -K_t_d / R_m_d;

% with disturbance:
beta_12_d = l_w_d;
beta_21_d = l_b_d;

delta_d = gamma_11_d*gamma_22_d - gamma_12_d*gamma_21_d;

a_22_d = (gamma_22_d*alpha_12_d - gamma_12_d*alpha_22_d) / delta_d;
a_23_d = (gamma_22_d*alpha_13_d - gamma_12_d*alpha_23_d) / delta_d;
a_24_d = (gamma_22_d*alpha_14_d - gamma_12_d*alpha_24_d) / delta_d;
a_42_d = (-gamma_21_d*alpha_12_d + gamma_11_d*alpha_22_d) / delta_d;
a_43_d = (-gamma_21_d*alpha_13_d + gamma_11_d*alpha_23_d) / delta_d;
a_44_d = (-gamma_21_d*alpha_14_d + gamma_11_d*alpha_24_d) / delta_d;

b_21_d = (gamma_22_d*beta_11_d)/delta_d;
b_41_d = (-gamma_21_d*beta_11_d)/delta_d;

%b_21_d = (gamma_22_d*beta_11_d - gamma_12_d*beta_21_d)/delta_d;
%b_41_d = (-gamma_21_d*beta_11_d + gamma_11_d*beta_21_d)/delta_d;

%%% comparison: difference betwin our state variables to theirs %%%

gamma_11_res = gamma_11-gamma_11_d;
gamma_12_res = gamma_12-gamma_12_d;
alpha_11_res = alpha_11-alpha_11_d;
alpha_12_res = alpha_12- alpha_12_d;
alpha_13_res = alpha_13-alpha_13_d;
alpha_14_res = alpha_14-alpha_14_d;
beta_11_res = beta_11-beta_11_d;
gamma_21_res = gamma_21-gamma_21_d;
gamma_22_res = gamma_22-gamma_22_d;
alpha_21_res = alpha_21-alpha_21_d;
alpha_22_res = alpha_22-alpha_22_d;
alpha_23_res =alpha_23-alpha_23_d;
alpha_24_res = alpha_24-alpha_24_d;
beta_21_res = beta_21-beta_21_d;

% with disturbance:
beta_12_res = beta_12-beta_12_d;
beta_21_res = beta_21-beta_21_d;

delta_res = delta-delta_d;

a_22_res = a_22-a_22_d;
a_23_res = a_23-a_23_d;
a_24_res = a_24-a_24_d;
a_42_res = a_42-a_42_d;
a_43_res = a_43-a_43_d;
a_44_res = a_44-a_44_d;

b_21_res = b_21-b_21_d;
b_41_res = b_41-b_41_d;

% our state variables:
A = [0 1 0 0; 0 a_22 a_23 a_24; 0 0 0 1; 0 a_42 a_43 a_44];
B = [0; b_21; 0; b_41];

% their state variables 
A_d = [0 1 0 0; 0 a_22_d a_23_d a_24_d; 0 0 0 1; 0 a_42_d a_43_d a_44_d];
B_d = [0; b_21_d; 0; b_41_d];

% comperation betwin our and theirs A matrix and B matrix
A_res = A-A_d;
%A_res = [0 1 0 0; 0 a_22_res a_23_res a_24_res; 0 0 0 1; 0 a_42_res a_43_res a_44_res];
%B_res = [0; b_21_res; 0; b_41_res];
B_res = B-B_d;

% more defenitions
C = [0 0 1 0];
% our trunsfer function
sys = ss(A,B,C,0);
t_f = tf(sys);
% our trunsfer function poles and zeros
p = pole(t_f);
z = zero(t_f);

% theirs trunsfer function
sys_d = ss(A_d,B_d,C,0);
t_f_d = tf(sys_d);
% theirs trunsfer function poles and zeros
p_d = pole(t_f_d);
z_d = zero(t_f_d);

%sys_res = ss(A_res,B_res,C,0);
%t_f_res = tf(sys_res);
%p_res = pole(t_f_res);
%z_res = zero(t_f_res);

%Arkadi's PID
P = 1.75;
D = 0.15;
I=0.1;
N = 50;
s = tf('s');
Controller = P+(I/s)+D*(N/(1+(N/s)));
%open loop with Arkadi's PID
L = t_f*Controller;
%open loop parameters - Arkadi's PID
p_L = pole(L);
z_L = zero(L);
%close loop with Arkadi's PID
T =L*12/(1+L*12); %duplicate by 12 same as the controller in the simulation
pole_end = pole(T);


% controlability test
ctrl = [B ,A*B, A^2*B, A^3*B]; 

rank_control = rank(ctrl);

% observability test
obsv = [C; C*A; C*A^2; C*A^3];

rank_observ = rank(obsv);

% my controller by calculation
kd = 0.39;
kp = 6.93;
ki = 5.97;
controller_pid = (kd*s^2+kp*s+ki)/s;

%open loop
L_pid = t_f*controller_pid;
T_pid = L_pid/(1+L_pid);

% Bode for my PID
bode(T_pid);
margin(T_pid);
title('Bode for my PID');
grid on;

% Bode for my open loop
bode(L_pid);
margin(L_pid);
title('Bode for my open loop');
grid on;

% Step for my PID 
figure;
step(T_pid);
title('step Response for close loop');
grid on; 
stepinfo(T_pid);

% Nyquist for my PID
figure;
nyquist(L_pid);
grid on; 

%root locus for my PID
figure;
rlocus(L_pid);

%nichols for my PID
figure;
nichols(L_pid);
grid on; 

