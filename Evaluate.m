%% Script 1
% Script to evaluate MM_Euler function over a range of a and beta values
% and plot the results as surfaces

% Define parameter ranges
a_values = 0:0.1:1;
beta_values = 0:0.1:1;
num_a = length(a_values);
num_beta = length(beta_values);

% Initialize result matrices
leader_results = zeros(num_a, num_beta);
follower_results = zeros(num_a, num_beta);

% Loop through all combinations of a and beta
for i = 1:num_a
    for j = 1:num_beta
        a = a_values(i);
        beta = beta_values(j);
        
        % Call the MM_Euler function (suppress display)
        [obj_follower, obj_leader] = MM_Euler(a, beta, false);
        
        % Store the mean results
        leader_results(i,j) = mean(obj_leader);
        follower_results(i,j) = mean(obj_follower);
        
        fprintf('Completed a=%.1f, beta=%.1f\n', a, beta);
    end
end

% Create meshgrid for surface plots
[A, BETA] = meshgrid(a_values, beta_values);

% Plot leader objective results
figure;
surf(A, BETA, leader_results');
xlabel('a');
ylabel('beta');
zlabel('Objective Value');
title('Leader Objective Function');
colorbar;
grid on;
view(3);

% Plot follower objective results
figure;
surf(A, BETA, follower_results');
xlabel('a');
ylabel('beta');
zlabel('Objective Value');
title('Follower Objective Function');
colorbar;
grid on;
view(3);

% Optional: Plot both surfaces together for comparison
figure;
hold on;
surf1 = surf(A, BETA, leader_results');
set(surf1, 'FaceColor', 'r', 'FaceAlpha', 0.7);
surf2 = surf(A, BETA, follower_results');
set(surf2, 'FaceColor', 'b', 'FaceAlpha', 0.7);
xlabel('a');
ylabel('beta');
zlabel('Objective Value');
title('Comparison of Leader (Red) and Follower (Blue) Objectives');
legend('Leader', 'Follower');
grid on;
view(3);
hold off;

%% Script 2
% Script to evaluate MM_Euler function for fixed a and varying beta
% and plot the results as 2D line plots

% Fixed parameter
a = 0.03;  % You can change this to any fixed value you want to analyze

% Define beta range
beta_values = 0:0.1:1;
num_beta = length(beta_values);

% Initialize result vectors
leader_results = zeros(1, num_beta);
follower_results = zeros(1, num_beta);

% Loop through all beta values
for j = 1:num_beta
    beta = beta_values(j);
    
    % Call the MM_Euler function (suppress display)
    [obj_follower, obj_leader] = MM_Euler(a, beta, false);
    
    % Store the mean results
    leader_results(j) = mean(obj_leader);
    follower_results(j) = mean(obj_follower);
    
    fprintf('Completed beta=%.1f\n', beta);
end

% Create 2D plots
figure;
hold on;
plot(beta_values, leader_results, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
plot(beta_values, follower_results, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('Beta');
ylabel('Objective Value');
title(sprintf('Objective Functions vs Beta (a = %.3f)', a));
legend('Leader Objective', 'Follower Objective', 'Location', 'best');
grid on;
hold off;

% Optional: Plot with error bars (if you want to show standard deviation)
% You would need to modify the loop to also store std values
% figure;
% errorbar(beta_values, leader_results, leader_stds, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
% hold on;
% errorbar(beta_values, follower_results, follower_stds, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b');
% xlabel('Beta');
% ylabel('Objective Value');
% title(sprintf('Objective Functions vs Beta (a = %.1f)', a));
% legend('Leader Objective', 'Follower Objective', 'Location', 'best');
% grid on;
% hold off;