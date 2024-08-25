% Initial parameters
params = struct(...
    'g', 9.8, ...
    'b_f', 0, ...
    'm_b', 0.3064, ...
    'l_b', 0.08, ...
    'I_b', 0.00261, ...
    'm_w', 0.0043, ...
    'l_w', 0.02, ...
    'I_w', 0.00000086, ...
    'R_m', 20.83, ...
    'K_e', 0.108, ...
    'K_t', 0.108);

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

% Function to generate all possible combinations of +10% and -10% for given parameters
function variations = generate_variations(paramNames)
    n = length(paramNames);
    variations = dec2bin(0:(2^n-1)) - '0'; % Generate binary variations
    variations = (variations * 2 - 1) * 10; % Convert to +10% and -10%
end

% Generate all possible variations of +10% and -10% for all parameters
variations = generate_variations(paramNames);

% Check the sensitivity for each combination of parameters
for i = 1:size(variations, 1)
    params_changed = params;
    variation = variations(i, :);
    
    change_description = ''; % To hold the description of the changes
    for k = 1:length(paramNames)
        param_name = paramNames{k};
        params_changed.(param_name) = params.(param_name) * (1 + variation(k) / 100);
        change_description = strcat(change_description, sprintf('%s %+d%%, ', param_name, variation(k)));
    end
    
    A_changed = calculate_A(params_changed);
    eigenvalues_changed = eig(A_changed);
    
    % Calculate the change in eigenvalues
    change = calculate_change(eigenvalues_initial, eigenvalues_changed);
    
    if any(change > 15)
        fprintf('Change in parameters %s results in eigenvalue change of %.2f%%\n', ...
            change_description(1:end-2), ... % Remove the last comma and space
            max(change));
    end
end
