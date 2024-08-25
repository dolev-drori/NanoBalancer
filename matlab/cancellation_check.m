% check if there is pole-zero cancellation wich is not allowed
% and makes problems in real life control

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

% our state variables:
A = [0 1 0 0; 0 a_22 a_23 a_24; 0 0 0 1; 0 a_42 a_43 a_44];
B = [0; b_21; 0; b_41];

% more defenitions
C = [0 0 1 0];
% our trunsfer function

% Create the state-space system
sys = ss(A, B, C, 0);

% Compute the transfer function from the state-space system
G = tf(sys);

% Find the poles and zeros of the transfer function
poles = pole(G);
zeros = zero(G);

% Display the poles and zeros
disp('Poles:');
disp(poles);
disp('Zeros:');
disp(zeros);

% Check if there are common poles and zeros
common_poles_zeros = intersect(poles, zeros);

if ~isempty(common_poles_zeros)
    disp('There are cancelled poles and zeros:');
    disp(common_poles_zeros);
else
    disp('There are no cancelled poles and zeros.');
end

% rank A matrix
rank_A = rank(A);
disp('A rank:');
disp(rank_A);

% controlability test
ctrl = [B ,A*B, A^2*B, A^3*B]; 
rank_control = rank(ctrl);
disp('controlability rank:');
disp(rank_control);

% observability test
obsv = [C; C*A; C*A^2; C*A^3];
rank_observ = rank(obsv);
disp('observability rank:');
disp(rank_observ);
