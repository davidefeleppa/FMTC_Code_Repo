function omega_t = calculate_omega_t(t, qmax, phi, kappa, theta, lambda_a, lambda_b, beta, a, b, gamma, T)
    % Initialize Matrices 
    A_matrix = zeros(2*qmax+1, 2*qmax+1);
    B_matrix = zeros(2*qmax+1, 1);
    
    % For each inventory level [qmax:-1:qmin]
    for i = 0:(2*qmax) 
        inventory = qmax - i;
        
        % Diagonal elements
        A_matrix(i+1, i+1) = -phi * kappa * inventory^2 + beta * kappa * (lambda_a - lambda_b) * inventory;
        
        % Vector elements
        B_matrix(i+1) = exp(kappa * ((a-b)/2)*inventory - (gamma + 0.5*(theta-beta))*inventory^2);
        
        % Upper diagonal (i+1 terms)
        if i < 2*qmax
            A_matrix(i+1, i+2) = lambda_b * exp(-1 + kappa*b - kappa*0.5*(theta+beta));
        end
        
        % Lower diagonal (i-1 terms)
        if i > 0
            A_matrix(i+1, i) = lambda_a * exp(-1 + kappa*a - kappa*0.5*(theta+beta));
        end
    end
    
    % Matrix exponential and multiplication
    omega_t = expm(A_matrix*(T-t)) * B_matrix;

    if (sum(omega_t < 0) > 0)
        omega_t = 0;
        warning("Omega negative for given parameters");
    end

end

function gt = calculate_gt(t, kappa, theta, qmax, phi, lambda_a, lambda_b, beta, a, b, gamma, T)
    % Calculate omega(t)
    omega_function = calculate_omega_t(t, qmax, phi, kappa, theta, lambda_a, lambda_b, beta, a, b, gamma, T);
    
    % Compute g(t)
    gt = (1 / kappa) * log(omega_function);
end

function [delta_a, delta_b] = calculate_deltas(t, T, q, q_tilde, qmax, kappa, phi, gamma, a, b, lambda_a, lambda_b, beta, theta, take_maximum)
    % Calculate g(t) for each inventory level
    g_q = calculate_gt(t, kappa, theta, qmax, phi, lambda_a, lambda_b, beta, a, b, gamma, T);
    
    % Convert current inventory to 1-based index
    idx = qmax - q + 1;
    
    % Symmetric boundary handling
    if q == -qmax
        % Minimum inventory (max short position): can only increase inventory
        c_hat_a = 1/kappa - g_q(idx) + g_q(idx);  % No change for asks
        c_hat_b = 1/kappa - g_q(idx-1) + g_q(idx); % Can only buy
    elseif q == qmax
        % Maximum inventory (max long position): can only decrease inventory
        c_hat_a = 1/kappa - g_q(idx+1) + g_q(idx); % Can only sell
        c_hat_b = 1/kappa - g_q(idx) + g_q(idx);    % No change for bids
    else
        % Normal case
        c_hat_a = 1/kappa - g_q(idx+1) + g_q(idx);
        c_hat_b = 1/kappa - g_q(idx-1) + g_q(idx);
    end
    
    % Apply maximum constraints if needed
    if take_maximum
        c_a = max(c_hat_a, a - 0.5*(theta+beta));
        c_b = max(c_hat_b, b - 0.5*(theta+beta));
    else
        c_a = c_hat_a;
        c_b = c_hat_b;
    end
    
    % Final delta calculations
    delta_a = c_a - beta*q_tilde + 0.5*(theta+beta) - theta*q;
    delta_b = c_b + beta*q_tilde + 0.5*(theta+beta) + theta*q;
end

% Stock Price
sigma = 1;
S0 = 100;

% Inventory
qmax = 10;
qmin = -qmax;

% Fill probability sensitivity 
kappa = 2; 

% Arrival Probabilities ~ Poission(λ)
lambda_0 = 10;

% Penalties
phi = 0.1;
gamma = 0.03;

% Parameters
t = 0.1;
T = 1;
a = 0.3;
b = a;
lambda_a = lambda_0 * exp(-kappa*a);
lambda_b = lambda_0 * exp(-kappa*b);
beta = 0.05;

%% Analyze over time

% Parameters
T = 1;  % Terminal time
time_vector = linspace(0, T, 50);  % Time grid
theta_values = linspace(-0.1,0.1,3);     % Different theta values to analyze
q_values = -10:1:10;        % Inventory levels to plot (-10 to +10 in steps of 1)
q_colors = parula(length(q_values)); % Color map for inventory levels

% Create figure
figure('Position', [100 100 1000 800]);

% Time dependent parameters
a_func = @(t, a0, a1) a + 0.5*t;
b_func = @(t, b0, b1) b + 0.5*t;

% Loop through theta values (rows)
for theta_idx = 1:length(theta_values)
    theta = theta_values(theta_idx);
    
    % Pre-compute deltas for this theta
    deltas_a = zeros(length(time_vector), length(q_values));
    deltas_b = zeros(length(time_vector), length(q_values));
    
    for q_idx = 1:length(q_values)
        for t_idx = 1:length(time_vector)
                
            a_t = a_func(time_vector(t_idx));
            b_t = b_func(time_vector(t_idx));
            lambda_a_t = lambda_0 * exp(-kappa*a_t);
            lambda_b_t = lambda_0 * exp(-kappa*b_t); 
                
            [deltas_a(t_idx,q_idx), deltas_b(t_idx,q_idx)] = calculate_deltas(...
                time_vector(t_idx), T, q_values(q_idx), 0, qmax, kappa, phi, gamma, ...
                a, b, lambda_a, lambda_b, beta, theta, false);

        end
    end
    
    % ===== Ask Spread Subplot =====
    subplot(length(theta_values), 2, (theta_idx-1)*2 + 1);
    hold on; grid on;
    
    % Plot ask spreads (skip minimum inventory level)
    for q_idx = 1:length(q_values)
        if q_values(q_idx) ~= qmin  % Skip minimum inventory for ask side
            plot(time_vector, deltas_a(:,q_idx), ...
                'Color', q_colors(q_idx,:), ...
                'LineWidth', 1.5, ...
                'DisplayName', sprintf('q=%d', q_values(q_idx)));
        end
    end
    
    title(sprintf('$\\theta=%.2f$', theta), 'Interpreter', 'latex', 'FontSize', 12);
    if theta_idx == length(theta_values)
        xlabel('Time', 'FontSize', 10);
    end
    ylabel('$\delta^{a, F \star}$', 'Interpreter', 'latex', 'FontSize', 10);
    xlim([0 T]);
    
    % ===== Bid Spread Subplot =====
    subplot(length(theta_values), 2, (theta_idx-1)*2 + 2);
    hold on; grid on;
    
    % Plot bid spreads (skip maximum inventory level)
    for q_idx = 1:length(q_values)
        if q_values(q_idx) ~= qmax  % Skip maximum inventory for bid side
            plot(time_vector, deltas_b(:,q_idx), ...
                'Color', q_colors(q_idx,:), ...
                'LineWidth', 1.5, ...
                'DisplayName', sprintf('q=%d', q_values(q_idx)));
        end
    end
    
    title(sprintf('$\\theta=%.2f$', theta), 'Interpreter', 'latex', 'FontSize', 12);
    if theta_idx == length(theta_values)
        xlabel('Time', 'FontSize', 10);
    end
    ylabel('$\delta^{b, F \star}$', 'Interpreter', 'latex', 'FontSize', 10);
    xlim([0 T]);
end

% Add colorbar
colormap(q_colors);
cbar = colorbar;
clim([min(q_values) max(q_values)]);
cbar.Label.String = 'Inventory';
cbar.Label.FontSize = 10;
set(cbar, 'Position', [0.92 0.15 0.02 0.7]);

% Adjust all axes font sizes
set(findall(gcf, 'Type', 'axes'), 'FontSize', 9);