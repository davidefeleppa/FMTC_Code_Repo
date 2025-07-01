function omega_t = calculate_omega_t(t, qmax, phi, kappa, lambda_a, lambda_b, beta, a, b, gamma, T)
    % Initialize A_matrix and vector
    A_matrix = zeros(2*qmax+1, 2*qmax+1);
    vector = zeros(2*qmax+1, 1);
    
    for i = 0:(2*qmax)  % MATLAB uses 1-based indexing, but we adjust with qmax-i
        inventory = qmax - i;
        
        % Diagonal elements
        A_matrix(i+1, i+1) = -phi * kappa * inventory^2 + beta * kappa * (lambda_a - lambda_b) * inventory;
        
        % Vector elements
        vector(i+1) = exp(kappa * ((a-b)/2)*inventory - (gamma - beta/2)*inventory^2);
        
        % Upper diagonal (i+1 terms)
        if i < 2*qmax
            A_matrix(i+1, i+2) = lambda_b * exp(-1 + kappa*b - kappa*beta/2);
        end
        
        % Lower diagonal (i-1 terms)
        if i > 0
            A_matrix(i+1, i) = lambda_a * exp(-1 + kappa*a - kappa*beta/2);
        end
    end
    
    % Matrix exponential and multiplication
    omega_t = expm(A_matrix*(T-t)) * vector;
end

function gt = calculate_gt(t, kappa, qmax, phi, lambda_a, lambda_b, beta, a, b, gamma, T)
    % Calculate omega_t by calling the previously converted function
    omega_function = calculate_omega_t(t, qmax, phi, kappa, lambda_a, lambda_b, beta, a, b, gamma, T);
    
    % Compute the final result
    gt = (1 / kappa) * log(omega_function);
end

function [delta_a, delta_b] = calculate_deltas(t, T, q, q_tilde, z, qmax, kappa, phi, gamma, a, b, lambda_a, lambda_b, beta, take_maximum)
    % Calculate g_qs by calling the previously converted function
    % Note: Need to pass all required parameters to calculate_gt
    g_qs = calculate_gt(t, kappa, qmax, phi, lambda_a, lambda_b, beta, a, b, gamma, T);
    
    % Calculate indices with clipping
    indices = min(max(qmax - q, 0), 2 * qmax);
    indices = int32(indices);  % Convert to integer indices
    indices_minus_one = min(max(indices + 1, 0), 2 * qmax);  % Lower inventory
    indices_plus_one = min(max(indices - 1, 0), 2 * qmax);   % Higher inventory
    
    % MATLAB uses 1-based indexing, so we need to adjust
    g_qs = squeeze(g_qs);  % Remove singleton dimensions if any
    
    % Calculate c_hat values
    c_hat_a = 1/kappa - g_qs(indices_minus_one + 1) + g_qs(indices + 1);
    c_hat_b = 1/kappa - g_qs(indices_plus_one + 1) + g_qs(indices + 1);
    
    % Apply maximum condition if needed
    if take_maximum
        c_a = max(c_hat_a, a - beta/2);
        c_b = max(c_hat_b, b - beta/2);
    else
        c_a = c_hat_a;
        c_b = c_hat_b;
    end
    
    % Calculate final deltas
    delta_a = c_a - beta * q_tilde - z + beta/2;
    delta_b = c_b + beta * q_tilde + z + beta/2;
end


sigma = 1;
S0 = 100;

T = 1;
lambda_a = 10;
lambda_b = 10; 

qmax = 10;
qmin = -qmax;
sigma_Z = 1;
kappa = 2;

% Model for misses
eta = 0.01; % TICK SIZE
a = 0.1;
b = 0.1;
beta = 0.05;

% Penalties
phi = 0.1;
gamma = 0.03;

all_inventories = qmax:-1:qmin-1;

% Numerics

% Parameters
t = 0.5;
z = 0;

% Create vectors (note: MATLAB's range syntax is start:step:end)
q_vector = (qmax-1):-1:qmin;
q_tilde_vector = 9:-1:-10; % Equivalent to np.arange(10-1,-10,step=-1)

% Initialize matrices
deltas_a = zeros(length(q_vector), length(q_tilde_vector));
deltas_b = zeros(length(q_vector), length(q_tilde_vector));

% Calculate deltas
for j = 1:length(q_tilde_vector)
    for i = 1:length(q_vector)
        [delta_a, delta_b] = calculate_deltas(t, T, q_vector(i), q_tilde_vector(j), z, qmax, kappa, phi, gamma, a, b, lambda_a, lambda_b, beta, false);
        deltas_a(i, j) = delta_a;
        deltas_b(i, j) = delta_b;
    end
end

% Plotting
figure;
hold on;

% Create colormap
cmap = colormap("default");
norm_min = min(q_tilde_vector);
norm_max = max(q_tilde_vector);

% Plot each line
for i = 1:length(q_tilde_vector)
    % Normalize for color
    color_weight = (q_tilde_vector(i) - norm_min) / (norm_max - norm_min);
    color_idx = max(1, min(size(cmap,1), round(color_weight * size(cmap,1))));
    line_color = cmap(color_idx, :);
    
    % Plot ask and bid deltas
    plot(q_vector, deltas_a(:,i), '-', 'Color', line_color, 'LineWidth', 1.5);
    plot(q_vector, deltas_b(:,i), '--', 'Color', line_color, 'LineWidth', 1.5);
end

% Add grid and labels
grid on;
alpha(0.3);
xlabel('$Q_t$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\hat{\delta}^{a}_t, \hat{\delta}^{b}_t$', 'Interpreter', 'latex', 'FontSize', 12);

% Add colorbar
colormap("default");
cbar = colorbar;
clim([norm_min norm_max]);
cbar.Label.String = '$\tilde{Q}_t$';
cbar.Label.Interpreter = 'latex';
cbar.Label.FontSize = 12;

% Adjust layout
set(gcf, 'Position', [100 100 600 400]); % Width=600, Height=400
set(gcf, 'PaperPositionMode', 'auto');
hold off;
