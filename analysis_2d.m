%% 2D Plot: Objective vs Beta for Fixed a
% This script evaluates MM_Euler for a fixed a and a range of beta values,
% and plots the leader/follower objective values.

rng(0);  % For reproducibility

% Simulation parameters
steps = 1e6;
time_vals = linspace(0,1,steps);

% Fixed value of 'a'
a_fixed = 0.7;

% Define beta range
beta_values = linspace(0, 0.3, 20);  % Finer resolution if desired
num_beta = length(beta_values);

% Prepare constant 'A'/'B' input
%A_fixed = ones(steps, 1) * a_fixed;

% Prepare a time dependand 'A'/'B' input
ref_a = 0.1;
A_ref = ref_a * exp(-time_vals);
B_ref = A_ref;

% Initialize result vectors
leader_results = zeros(1, num_beta);
follower_results = zeros(1, num_beta);

fprintf('\n--- Running 2D Beta Sweep for a = %.3f ---\n', a_fixed);
for j = 1:num_beta
    beta = beta_values(j);
    
    % Run MM_Euler
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, obj_follower, obj_leader] = MM_Euler(A_ref, B_ref, beta, 1000);

    % Store results
    leader_results(j) = mean(obj_leader);
    follower_results(j) = mean(obj_follower);

    fprintf('Completed beta = %.3f\n', beta);
end

% Plot the results
figure;
plot(beta_values, leader_results, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r'); hold on;
plot(beta_values, follower_results, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('\beta', 'FontWeight', 'bold');
ylabel('Objective Value', 'FontWeight', 'bold');
title(sprintf('Objective vs \\beta (a = %.3f)', a_fixed), 'FontWeight', 'bold');
legend({'Leader', 'Follower'}, 'Location', 'best');
grid on;