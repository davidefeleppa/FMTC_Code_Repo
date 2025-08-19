%% Evaluate MM_Matrix over a Grid of (a, beta) for 3 different theta values
% This script evaluates MM_Matrix over a grid of a(t) & beta pairs for multiple theta values
% and plots the leader objective contour plots in 3 subplots on a single line.

% For reproducibility
rng(0);

% Simulations
sims = 10000;

% Fixed Variables
phi_float = 0.1;  % Penalty parameter
gamma_float = 0.03; % Penalty parameter

% Parameter Ranges
a_values = linspace(0, 1.0, 25);  % a values (b will equal a)
beta_values = linspace(0, 1.0, 25);  % beta values

% Different theta values to test (reduced to 3)
theta_values = [0.0, 0.10, 0.50];

% Initialize result matrices for each theta
num_a = length(a_values);
num_beta = length(beta_values);
num_theta = length(theta_values);
leader_results = zeros(num_a, num_beta, num_theta);

% Loop through each theta, a, and beta combination
fprintf('\n--- Running 3D Parameter Sweep for %d theta values ---\n', num_theta);
total_iterations = num_a * num_beta * num_theta;

for k = 1:num_theta
    theta = theta_values(k);
    fprintf('Processing theta = %.3f (%d/%d)\n', theta, k, num_theta);
    
    for i = 1:num_a
        for j = 1:num_beta
            
            % Current parameter values
            a_val = a_values(i);
            beta_val = beta_values(j);
            
            % Time dependent parameters - setting b = a
            a_func = @(t) a_val;  % Use the current a value
            b_func = @(t) a_val;  % Set b = a
            
            % Run simulation
            [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, obj_follower, obj_leader] = ...
                MM_Matrix(a_func, b_func, beta_val, theta, phi_float, gamma_float, sims);
            
            % Store mean objective values
            leader_results(i, j, k) = mean(obj_leader);
            
            % Progress Bar
            current_progress = (k-1)*num_a*num_beta + (i-1)*num_beta + j;
            if mod(current_progress, 20) == 0 || current_progress == total_iterations
                fprintf('Overall Progress: [%-20s] %d/%d (%.1f%%)\n', ...
                    repmat('=',1,round(20*current_progress/total_iterations)), ...
                    current_progress, total_iterations, 100*current_progress/total_iterations);
            end
        end
    end
end

% Meshgrid for surface plotting
[AA, BB] = meshgrid(a_values, beta_values);

% Create figure with 3 subplots in a single row for leader contour plots
figure('Position', [100, 100, 1800, 500]);

for k = 1:num_theta
    subplot(1, 3, k);
    
    % Create contour plot for leader objective function
    contourf(AA, BB, leader_results(:, :, k)', 20, 'LineColor', 'none');
    colormap(jet);
    colorbar;
    
    % Set labels with LaTeX interpreter
    xlabel('$a$ (= $b$)', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('$\beta$', 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold');
    % title(sprintf('$\\theta = %.2f$', theta_values(k)), 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
    
    grid on;
    hold on;
    [C,h] = contour(AA, BB, leader_results(:, :, k)', 10, 'k-');
    clabel(C,h,'FontSize',8,'Color','k','LabelSpacing',500);
    hold off;
    
    % Find and mark optimal point for this theta
    [max_val, max_idx] = max(leader_results(:, :, k), [], 'all');
    [opt_i, opt_j] = ind2sub(size(leader_results(:, :, k)), max_idx);
    
    % Uncomment to mark optimal points
    % hold on;
    % plot(a_values(opt_i), beta_values(opt_j), 'r*', 'MarkerSize', 15, 'LineWidth', 3);
    % hold off;
    
    % Print optimal values
    fprintf('θ=%.3f: Optimal a=%.3f, β=%.3f, objective=%.3f\n', ...
        theta_values(k), a_values(opt_i), beta_values(opt_j), max_val);
end

% Add overall title with LaTeX interpreter
% sgtitle('Leader Objective Function for Different $\theta$ Values', ...
%    'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

% Create summary table of optimal values
fprintf('\n=== SUMMARY OF OPTIMAL VALUES ===\n');
fprintf('θ\t\tOptimal a\tOptimal β\tMax Objective\n');
fprintf('-----------------------------------------------\n');
for k = 1:num_theta
    [max_val, max_idx] = max(leader_results(:, :, k), [], 'all');
    [opt_i, opt_j] = ind2sub(size(leader_results(:, :, k)), max_idx);
    fprintf('%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n', ...
        theta_values(k), a_values(opt_i), beta_values(opt_j), max_val);
end