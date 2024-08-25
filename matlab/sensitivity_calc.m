%%%% our phisical parameters: %%%%
params = struct(...
'g', 9.8, ... % Gravity [m/s^2]
'b_f', 0, ... % Friction coefficient
'm_b', 0.3064, ... % Body mass [kg]
'l_b', 0.08, ... % The distance between the mass centers of the body and wheels [m]
'I_b', 0.00261, ... % Body inertia torque [kgm^2]
'm_w', 0.0043, ... % Wheel mass [kg]
'l_w', 0.02, ... % Wheel radius [m]
'I_w', 0.00000086, ... % Wheel inertia torque [kgm^2]
'R_m', 20.83, ... % Electrical resistance of motor [ohm]
'l_m', 0, ... % ?
'b_m', 0, ... % ?
'K_e', 0.108,  ...% Motor Electrical constant [s/r]
'K_t', 0.108); % Motor torque constant [Nm/A]

% Get the field names
paramNames = fieldnames(params);

% Function to calculate the state-space matrix A
function A = calculate_A(params)
    gamma_11 = (params.I_w)/(params.l_w) + params.l_w * params.m_b + params.l_w * params.m_w;
    gamma_12 = params.m_b * params.l_b * params.l_w;
    alpha_12 = -(((params.K_e * params.K_t) / params.R_m) + params.b_f) / params.l_w;
    alpha_14 = ((params.K_e * params.K_t) / params.R_m) + params.b_f;
    beta_11 = params.K_t / params.R_m;
    
    gamma_21 = params.m_b * params.l_b;
    gamma_22 = params.I_b + (params.m_b * params.l_b^2);
    alpha_22 = (((params.K_e * params.K_t) / params.R_m) + params.b_f) / params.l_w;
    alpha_23 = params.m_b * params.l_b * params.g;
    alpha_24 = - (((params.K_e * params.K_t) / params.R_m) + params.b_f);
    beta_21 = -params.K_t / params.R_m;
    
    delta = gamma_11 * gamma_22 - gamma_12 * gamma_21;
    
    a_22 = (gamma_22 * alpha_12 - gamma_12 * alpha_22) / delta;
    a_23 = (- gamma_12 * alpha_23) / delta;
    a_24 = (gamma_22 * alpha_14 - gamma_12 * alpha_24) / delta;
    a_42 = (-gamma_21 * alpha_12 + gamma_11 * alpha_22) / delta;
    a_43 = (gamma_11 * alpha_23) / delta;
    a_44 = (-gamma_21 * alpha_14 + gamma_11 * alpha_24) / delta;

    A = [0 1 0 0; 0 a_22 a_23 a_24; 0 0 0 1; 0 a_42 a_43 a_44];
end

% Calculate initial eigenvalues
A_initial = calculate_A(params);
eigenvalues_initial = eig(A_initial);

% Function to calculate percentage change in eigenvalues
function change = calculate_change(original, new)
    change = abs(new - original) ./ abs(original) * 100;
end

% Function to generate combinations of parameter indices
function combs = generate_combinations(n)
    combs = {};
    for k = 1:n
        combs = [combs; num2cell(nchoosek(1:n, k), 2)];
    end
end

% Generate all combinations of parameter indices
param_indices = 1:length(paramNames);
combinations = generate_combinations(length(param_indices));

% Check the sensitivity for each combination of parameters
for c = 1:length(combinations)
    comb = combinations{c};
    
    % Increase parameters by 10%
    params_changed = params;
    for k = 1:length(comb)
        param_name = paramNames{comb(k)};
        params_changed.(param_name) = params.(param_name) * 1.1;
    end
    A_changed = calculate_A(params_changed);
    eigenvalues_changed = eig(A_changed);
    
    % Calculate the change in eigenvalues
    change = calculate_change(eigenvalues_initial, eigenvalues_changed);
    
    if any(change > 15)
        fprintf('Change in parameters %s by +10%% results in eigenvalue change of %.2f%%\n', ...
            strjoin(paramNames(comb), ', '), max(change));
    end
    
    % Decrease parameters by 10%
    params_changed = params;
    for k = 1:length(comb)
        param_name = paramNames{comb(k)};
        params_changed.(param_name) = params.(param_name) * 0.9;
    end
    A_changed = calculate_A(params_changed);
    eigenvalues_changed = eig(A_changed);
    
    % Calculate the change in eigenvalues
    change = calculate_change(eigenvalues_initial, eigenvalues_changed);
    
    if any(change > 15)
        fprintf('Change in parameters %s by -10%% results in eigenvalue change of %.2f%%\n', ...
            strjoin(paramNames(comb), ', '), max(change));
    end
end