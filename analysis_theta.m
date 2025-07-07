%% 2D Plot: Objective vs Theta for Fixed Beta
% Evaluates MM_Matrix for a fixed beta and varying theta values,
% then plots the leader/follower objective values.

% For reproducibility
rng(0);  

% Number of simulations
sims = 10000;

% Parameter ranges
theta_values = linspace(-0.5, 0.5, 10);  % Finer resolution if needed

% Fixed variables
beta = 0.05; 

% Time dependant variables
a_func = @(t, a0, a1) 0.1;
b_func = @(t, b0, b1) 0.1;

% Initialize results storage
num_theta = length(theta_values);
leader_results = zeros(1, num_theta);
follower_results = zeros(1, num_theta);

fprintf('\n--- Running 2D Theta Sweep ---\n');

% Main simulation loop
for j = 1:num_theta
    theta = theta_values(j);
    
    % Run simulation
    [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, obj_follower, obj_leader] = MM_Euler(a_func, b_func, beta, theta, sims);

    % Store results
    leader_results(j) = mean(obj_leader);
    follower_results(j) = mean(obj_follower);

    fprintf('Completed theta = %.3f (%d/%d)\n', theta, j, num_theta);
end

% Plot results
figure;
plot(theta_values, leader_results, 'r-o', ...
    'LineWidth', 2, ...
    'MarkerFaceColor', 'r', ...
    'DisplayName', 'Leader');
hold on;
plot(theta_values, follower_results, 'b-s', ...
    'LineWidth', 2, ...
    'MarkerFaceColor', 'b', ...
    'DisplayName', 'Follower');

% Format plot
xlabel('\theta', 'FontWeight', 'bold');
ylabel('Objective Value', 'FontWeight', 'bold');
title('Objective vs \theta', 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;