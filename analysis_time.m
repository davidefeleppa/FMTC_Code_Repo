%% Evaluate Game Over Time
fprintf('\n--- Running Game---\n');

% Simulation parameters
steps = 1e6;
time_vals = linspace(0,1,steps);

% A and B functions
% ref_a = 0.1;
% ref_beta = 0.03;
% A_ref = ref_a * exp(-time_vals);
% B_ref = A_ref;

A = 0.1 * ones(steps, 1);
B = A;
beta = 0.05;
theta = 0.05;

[time, Q, Q_tilde, da, db, ta, tb, qmax, qmin, X, X_tilde, pnl, pnl_tilde, obj_follower, obj_leader] = MM_Matrix(A, B, beta, theta, 1);

% Fix terminal deltas/tildes for plotting
da(:,end) = da(:,end-1); db(:,end) = db(:,end-1);
ta(:,end) = ta(:,end-1); tb(:,end) = tb(:,end-1);

% Display summary
fprintf('--- Reference MM Results ---\n');
fprintf('Terminal Q:   mean = %.4f, std = %.4f\n', mean(Q(:,end)), std(Q(:,end)));
fprintf('Terminal X:   mean = %.4f, std = %.4f\n', mean(X(:,end)), std(X(:,end)));
fprintf('Terminal PnL: mean = %.4f, std = %.4f\n', mean(pnl), std(pnl));
fprintf('Objective:    mean = %.4f, std = %.4f\n', mean(obj_follower), std(obj_follower));

plot_time_series_results(time, Q, Q_tilde, da, db, ta, tb, qmax, qmin);