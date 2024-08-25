% check if there is pole-zero cancellation wich is not allowed
% and makes problems in real life control
% This code was written by ChatGPT

% Define the state-space matrices based on your system parameters
A = [0 1 0 0;
     0 -((params.K_e * params.K_t) / params.R_m + params.b_f) / params.l_w params.m_b * params.l_b * params.g / params.I_b 0;
     0 0 0 1;
     0 ((params.K_e * params.K_t) / params.R_m + params.b_f) / params.l_w -params.m_b * params.g * params.l_b / params.I_b -params.b_f / params.I_b]; 

B = [0; 
    params.K_t / params.R_m; 
    0; 
    -params.K_t / params.R_m]; 

C = [1 0 0 0]; % Output matrix

D = 0; % Feedforward matrix (assumed to be zero)

% Create the state-space system
sys = ss(A, B, C, D);

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
