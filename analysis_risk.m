%% Evaluate MM_Matrix over a Grid of (phi, gamma) for a single theta value
% This script evaluates MM_Matrix over a grid of phi & gamma pairs for a single theta value
% and plots the leader objective contour plot.

% For reproducibility
rng(0);

% Simulations
sims = 1000;
% Fixed Variables
a_val = 0.1; % Fixed a value
b_val = 0.1; % Fixed b value (set equal to a)
beta_val = 0.05; % Fixed beta value
theta_val = -0.05; % Single theta value

% Parameter Ranges
phi_values = linspace(0, 0.2, 20); % phi values
gamma_values = linspace(0, 0.1, 20); % gamma values

% Initialize result matrix
num_phi = length(phi_values);
num_gamma = length(gamma_values);
leader_results = zeros(num_phi, num_gamma);

% Time dependent parameters - fixed values
a_func = @(t) a_val;
b_func = @(t) b_val;

% Loop through each phi and gamma combination
fprintf('\n--- Running 2D Parameter Sweep for θ = %.3f ---\n', theta_val);
total_iterations = num_phi * num_gamma;

for i = 1:num_phi

    for j = 1:num_gamma
        % Current parameter values
        phi_val = phi_values(i);
        gamma_val = gamma_values(j);
        % Run simulation
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, obj_follower, obj_leader] = ...
            MM_Matrix(a_func, b_func, beta_val, theta_val, phi_val, gamma_val, sims);
        % Store mean objective values
        leader_results(i, j) = mean(obj_leader);
        % Progress Bar
        current_progress = (i-1)*num_gamma + j;
        if mod(current_progress, 10) == 0 || current_progress == total_iterations
            fprintf('Progress: [%-20s] %d/%d (%.1f%%)\n', ...
            repmat('=',1,round(20*current_progress/total_iterations)), ...
            current_progress, total_iterations, 100*current_progress/total_iterations);
        end
    end
end

% Smooth the surface before plotting
leader_results_smooth = smoothdata(smoothdata(leader_results, 1, 'movmean', 3), 2, 'movmean', 3);

% Meshgrid for surface plotting
[PHI, GAMMA] = meshgrid(phi_values, gamma_values);

% Create figure for leader contour plot
figure('Position', [100, 100, 800, 600]);

% Create contour plot for leader objective function
contourf(PHI, GAMMA, leader_results_smooth', 20, 'LineColor', 'none');
colormap(jet);
colorbar;
xlabel('$\phi$', 'FontWeight', 'bold');
ylabel('$\gamma$', 'FontWeight', 'bold');
title(sprintf('Leader Objective Function (θ = %.3f, a = %.1f, b = %.1f, β = %.3f)', ...
 theta_val, a_val, b_val, beta_val), 'FontWeight', 'bold');
grid on;
hold on;
[C,h] = contour(PHI, GAMMA, leader_results_smooth', 10, 'k-');
clabel(C,h,'FontSize',8,'Color','k','LabelSpacing',500);

% Find and mark optimal point
%[max_val, max_idx] = max(leader_results, [], 'all');
%[opt_i, opt_j] = ind2sub(size(leader_results), max_idx);
% plot(phi_values(opt_i), gamma_values(opt_j), 'r', 'MarkerSize', 15, 'LineWidth', 3);
hold off;