%% Evaluate MM_Euler over a Grid of (a, beta)
% This script evaluates MM_Euler over a grid of (a, beta) pairs
% and plots the mean objective values as surface plots.

rng(0);  % For reproducibility

% Simulation parameters
steps = 1e6;
time_vals = linspace(0,1,steps);

% Parameter ranges
a_values = linspace(0, 100, 5);
beta_values = linspace(0, 100, 5);
num_a = length(a_values);
num_beta = length(beta_values);

% Simulation
sims = 1000;
steps = 1e6; 

% Initialize result matrices
leader_results = zeros(num_a, num_beta);
follower_results = zeros(num_a, num_beta);

% Loop through each (a, beta) pair
fprintf('\n--- Running 3D Parameter Sweep ---\n');
for i = 1:num_a
    for j = 1:num_beta
        a = a_values(i);
        beta = beta_values(j);

        A = ones(steps, 1) * a;  % Constant a over time

        %A = a * exp(-time_vals); % Time dependant

        % Run simulation
        [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, obj_follower, obj_leader] = MM_Euler(A, A, beta, sims);

        % Store mean objective values
        leader_results(i, j) = mean(obj_leader);
        follower_results(i, j) = mean(obj_follower);

        fprintf('Completed a = %.3f, beta = %.3f\n', a, beta);
    end
end

% Meshgrid for surface plotting
[AA, BB] = meshgrid(a_values, beta_values);

% Plot: Leader Objective Surface
figure;
surf(AA, BB, leader_results');
xlabel('a', 'FontWeight', 'bold');
ylabel('\beta', 'FontWeight', 'bold');
zlabel('Objective Value', 'FontWeight', 'bold');
title('Leader Objective Surface', 'FontWeight', 'bold');
colorbar;
grid on;
view(135, 30);
shading interp;

% Plot: Follower Objective Surface
figure;
surf(AA, BB, follower_results');
xlabel('a', 'FontWeight', 'bold');
ylabel('\beta', 'FontWeight', 'bold');
zlabel('Objective Value', 'FontWeight', 'bold');
title('Follower Objective Surface', 'FontWeight', 'bold');
colorbar;
grid on;
view(135, 30);
shading interp;

% Combined Plot: Leader vs Follower
figure;
hold on;
surf1 = surf(AA, BB, leader_results');
set(surf1, 'FaceColor', 'r', 'FaceAlpha', 0.6, 'EdgeAlpha', 0.3);
surf2 = surf(AA, BB, follower_results');
set(surf2, 'FaceColor', 'b', 'FaceAlpha', 0.6, 'EdgeAlpha', 0.3);
xlabel('a', 'FontWeight', 'bold');
ylabel('\beta', 'FontWeight', 'bold');
zlabel('Objective Value', 'FontWeight', 'bold');
title('Leader (Red) vs Follower (Blue) Objectives', 'FontWeight', 'bold');
grid on;
view(135, 30);
legend({'Leader', 'Follower'}, 'Location', 'best');
hold off;

% Create contour plot for leader objective function
figure;
contourf(AA, BB, leader_results', 20, 'LineColor', 'none'); % 20 contour levels, no lines
colormap(jet); % Use jet colormap (or try 'parula', 'hot', etc.)
colorbar; % Show color scale
xlabel('a', 'FontWeight', 'bold');
ylabel('\beta', 'FontWeight', 'bold');
title('Leader Objective Function (Contour Plot)', 'FontWeight', 'bold');
grid on;
hold on;
[C,h] = contour(AA, BB, leader_results', 10, 'k-'); % 10 black contour lines
clabel(C,h,'FontSize',8,'Color','k','LabelSpacing',500);
hold off;