%% 2D Plot: Objective vs Beta for Fixed a
% This script evaluates MM_Euler for a fixed a and a range of beta values,
% and plots the leader/follower objective values.

rng(0);  % For reproducibility

% Simulation parameters
steps = 1e6;
time_vals = linspace(0,1,steps);

% Define beta range
theta_values = linspace(-0.5, 0.5, 20);  % Finer resolution if desired
num_theta = length(theta_values);

beta = 0.15;

% Prepare a time dependand 'A'/'B' input
A = 0.1 * ones(steps, 1);
B = A;
beta = 0.05;

% Initialize result vectors
leader_results = zeros(1, num_theta);
follower_results = zeros(1, num_theta);

fprintf('\n--- Running 2D Theta Sweep ---\n');
for j = 1:num_theta
    theta = theta_values(j);
    
    % Run MM_Euler
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, obj_follower, obj_leader] = MM_Matrix(A, B, beta, theta, 1000);

    % Store results
    leader_results(j) = mean(obj_leader);
    follower_results(j) = mean(obj_follower);

    fprintf('Completed theta = %.3f\n', theta);
end

% Plot the results
figure;
plot(theta_values, leader_results, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r'); hold on;
plot(theta_values, follower_results, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('theta', 'FontWeight', 'bold');
ylabel('Objective Value', 'FontWeight', 'bold');
title(sprintf('Objective vs \\theta'), 'FontWeight', 'bold');
legend({'Leader', 'Follower'}, 'Location', 'best');
grid on;