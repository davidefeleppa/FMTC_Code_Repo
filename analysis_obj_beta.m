%% 2D Plot: Leader Objective vs Theta for Different Beta Values
% Evaluates MM_Matrix for varying beta and theta values,
% then plots the leader objective values for each beta.

% For reproducibility
rng(0);

% Number of simulations
sims = 10000;

% Parameter ranges
theta_values = linspace(-0.25, 0.25, 50); % Finer resolution if needed
beta_values = linspace(0, 0.2, 10); % 10 beta values from 0 to 0.2

% Fixed variables
phi = 0.1;
gamma = 0.03;

% Time dependant variables
a_func = @(t, a0, a1) 0.1;
b_func = @(t, b0, b1) 0.1;

% Initialize results storage
num_theta = length(theta_values);
num_beta = length(beta_values);
leader_results = zeros(num_beta, num_theta);

fprintf('\n--- Running 2D Theta-Beta Sweep (Leader Only) ---\n');

% Main simulation loop
for i = 1:num_beta
    beta = beta_values(i);
    fprintf('Processing beta = %.3f (%d/%d)\n', beta, i, num_beta);
    
    for j = 1:num_theta
        theta = theta_values(j);
        
        % Run simulation
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, obj_follower, obj_leader] = MM_Matrix(a_func, b_func, beta, theta, phi, gamma, sims);
        
        % Store results (only leader)
        leader_results(i, j) = mean(obj_leader);
        
        fprintf('  Completed theta = %.3f (%d/%d)\n', theta, j, num_theta);
    end
end

% Plot results
figure;
hold on;

% Define colors for different beta lines
colors = parula(num_beta); % Generate distinct colors

% Smooth and plot the lines
for i = 1:num_beta
    % Smooth the data using moving average
    smoothed_results = smoothdata(leader_results(i, :), 'movmean', 5);
    
    plot(theta_values, smoothed_results, ...
        'Color', colors(i, :), ...
        'LineWidth', 2);
end

% Format plot
xlabel('\theta', 'FontWeight', 'bold');
grid on;

% Add colorbar
colormap(colors);
cbar = colorbar;
clim([min(beta_values) max(beta_values)]);
cbar.Label.String = '\beta';
cbar.Label.FontSize = 20;
set(cbar, 'Position', [0.92 0.15 0.02 0.7]);

hold off;