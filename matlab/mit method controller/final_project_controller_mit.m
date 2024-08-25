s = tf("s");
%ksi = -log(20/100) / sqrt(pi^2 + (log(20/100))^2);
ts = 10; % settling time [sec] 
k=4; % the order of the system

% the poles according to the method:
p1 = -4.0156+5.0723j;
p2 = -4.0156-5.0723j;

p3 = -5.5281+1.6553j;
p4 = -5.5281-1.6553j;

% We will divide the pole locations by the settling time
p1 = p1/ts;
p2 = p2/ts;
p3 = p3/ts;
p4 = p4/ts;

phi = (s+p1)*(s+p2)*(s+p3);

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
beta_22 = l_b;

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

A = [0 1 0 0; 0 a_22 a_23 a_24; 0 0 0 1; 0 a_42 a_43 a_44];
B = [0; b_21; 0; b_41];

% more defenitions
C = [0 0 1 0];

% Calculate the feedback gain matrix K using acker
poles = [p1, p2, p3, p4];
K = acker(A, B, poles);

% Display the result
disp('The feedback gain matrix K is:');
disp(K);

% Define the state-space matrices for the closed-loop system
A_cl = A - B * K; % Closed-loop A matrix

% Define the D matrix (usually zero in state-space representations for feedback systems)
D_cl = 0;

% Define the state-space system for the closed-loop system
sys_cl = ss(A_cl, B, C, D_cl);

% Display the closed-loop system
disp('Closed-loop system:');
disp(sys_cl);

% Optionally, plot the step response of the closed-loop system
figure;
step(sys_cl);
title('Closed-Loop Step Response');

info = stepinfo(sys_cl);
disp('Step response characteristics:');
disp(info);

% Calculate poles of the closed-loop system
poles_cl = pole(sys_cl);
disp('Poles of the closed-loop system:');
disp(poles_cl);

% Check if all poles are in the left half-plane
if all(real(poles_cl) < 0)
    disp('The system is stable.');
else
    disp('The system is not stable.');
end