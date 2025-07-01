function g = euler_scheme(steps, T, qmax, qmin, a, b, gamma, beta, kappa, phi, lambda_a, lambda_b)
    % Initialize variables
    h = T / steps;
    all_inventories = qmax:-1:qmin; 
    g = zeros(steps, length(all_inventories));
    
    % Terminal condition
    for iq = 1:length(all_inventories)
        q = all_inventories(iq);
        g(end, iq) = ((a-b)/2)*q - (gamma-beta/2)*(q^2);
    end
    
    % Time stepping (backward)
    for it = 1:(steps-1)
        for iq = 1:length(all_inventories)
            q = all_inventories(iq);
            
            % Interior points
            if iq ~= 1 && iq ~= length(all_inventories)
                c_star_a = max(1/kappa - g(end-it+1, iq-1) + g(end-it+1, iq), a - beta/2);
                c_star_b = max(1/kappa - g(end-it+1, iq+1) + g(end-it+1, iq), b - beta/2);
                
                g(end-it, iq) = g(end-it+1, iq) + h * ( ...
                    -phi*(q^2) ...
                    + (lambda_a-lambda_b)*beta*q ...
                    + lambda_a*exp(-kappa*(c_star_a + beta/2 - a))*(c_star_a + g(end-it+1,iq-1) - g(end-it+1,iq)) ...
                    + lambda_b*exp(-kappa*(c_star_b + beta/2 - b))*(c_star_b + g(end-it+1,iq+1) - g(end-it+1,iq)));
            
            % Lower boundary (iq == 1)
            elseif iq == 1
                c_star_b = max(1/kappa - g(end-it+1, iq+1) + g(end-it+1, iq), b - beta/2);
                
                g(end-it, iq) = g(end-it+1, iq) + h * ( ...
                    -phi*(q^2) ...
                    + (lambda_a-lambda_b)*beta*q ...
                    + lambda_b*exp(-kappa*(c_star_b + beta/2 - b))*(c_star_b + g(end-it+1,iq+1) - g(end-it+1,iq)));
            
            % Upper boundary (iq == end)
            else
                c_star_a = max(1/kappa - g(end-it+1, iq-1) + g(end-it+1, iq), a - beta/2);
                
                g(end-it, iq) = g(end-it+1, iq) + h * ( ...
                    -phi*(q^2) ...
                    + (lambda_a-lambda_b)*beta*q ...
                    + lambda_a*exp(-kappa*(c_star_a + beta/2 - a))*(c_star_a + g(end-it+1,iq-1) - g(end-it+1,iq)));
            end
        end
    end
end

function [delta_a, delta_b] = calculate_deltas_euler(t, q, q_tilde, z, gt, gt_q_minus, gt_q_plus, kappa, a, b, beta, take_maximum)
    % Calculate c_hat values
    c_hat_a = 1/kappa - gt_q_minus + gt;
    c_hat_b = 1/kappa - gt_q_plus + gt;
    
    % Apply maximum condition if needed
    if take_maximum
        c_a = max(c_hat_a, a - beta/2);
        c_b = max(c_hat_b, b - beta/2);
    else
        c_a = c_hat_a;
        c_b = c_hat_b;
    end
    
    % Calculate final deltas
    delta_a = c_a - beta*q_tilde - z + beta/2;
    delta_b = c_b + beta*q_tilde + z + beta/2;
end

function [arrival_a, arrival_b] = get_arrival(sims, dt, lambda_a, lambda_b)
    % Generate uniform random numbers
    unif_a = rand(sims, 1);
    unif_b = rand(sims, 1);
    
    % Determine arrivals using Poisson process probabilities
    arrival_a = unif_a < (1 - exp(-lambda_a * dt));
    arrival_b = unif_b < (1 - exp(-lambda_b * dt));
    
end

% Parameters
sigma = 1;
S0 = 100;

T = 1;
lambda_a = 10;
lambda_b = 10; 

qmax = 10;
qmin = -qmax;
sigma_Z = 1;
kappa = 2;

% Model for misses
eta = 0.01; % TICK SIZE
a = 0.1;
b = 0.1;
beta = 0.05;

% Penalties
phi = 0.1;
gamma = 0.03;

all_inventories = qmax:-1:qmin-1;

t = 0.5;
z = 0;

Nt = 1000;
sims = 10000;
dt = T/Nt;
euler_steps = 1000000;

% Initialize arrays
X = zeros(sims, Nt+1);
Q = zeros(sims, Nt+1);
Q_tilde = zeros(sims, Nt+1);
S = zeros(sims, Nt+1);
Z = zeros(sims, Nt+1);
pnl_euler = zeros(sims, 1);
intQ = zeros(sims, 1);
objective_euler = zeros(sims, 1);

deltas_a = zeros(sims, Nt+1);
deltas_b = zeros(sims, Nt+1);
tildes_a = zeros(sims, Nt+1);
tildes_b = zeros(sims, Nt+1);

arrivaltimes_a = zeros(sims, Nt+1);
arrivaltimes_b = zeros(sims, Nt+1);

S(:,1) = S0;  % MATLAB uses 1-based indexing
time = linspace(0, T, Nt+1);

% Run Euler scheme
g = euler_scheme(euler_steps, T, qmax, qmin, a, b, gamma, beta, kappa, phi, lambda_a, lambda_b);

% Round all_inventories
all_inventories_list = round(all_inventories, 2);

% Main simulation loop
for it = 1:Nt
    t = time(it);
    
    % Get current g values
    gt_all_q = g(round(it*(euler_steps/Nt)), :);
    G = zeros(sims, 1);
    G_down = zeros(sims, 1);
    G_up = zeros(sims, 1);
    
    for sim = 1:sims
        q = Q(sim,it);
        [~, idx] = min(abs(all_inventories_list - q));
        G(sim) = gt_all_q(idx);
        
        if q ~= qmin
            [~, idx_down] = min(abs(all_inventories_list - (q-1)));
            G_down(sim) = gt_all_q(idx_down);
        end
        
        if q ~= qmax
            [~, idx_up] = min(abs(all_inventories_list - (q+1)));
            G_up(sim) = gt_all_q(idx_up);
        end
    end
    
    % Calculate deltas
    [delta_a, delta_b] = calculate_deltas_euler(t, Q(:,it), Q_tilde(:,it), Z(:,it), G, G_down, G_up, ...
                         kappa, a, b, beta, true);
    
    deltas_a(:, it) = delta_a;
    deltas_b(:, it) = delta_b;
    
    % Update Z and S (Brownian motions)
    Z(:,it+1) = Z(:,it) + sqrt(dt) * randn(sims, 1);
    S(:,it+1) = S(:,it) + sigma * sqrt(dt) * randn(sims, 1);
    
    % Simulate order arrivals
    [arrival_a, arrival_b] = get_arrival(sims, dt, lambda_a, lambda_b);
    arrivaltimes_a(:,it) = arrival_a * t;
    arrivaltimes_b(:,it) = arrival_b * t;
    
    % Calculate fill probabilities
    fill_a = rand(sims, 1);
    fill_b = rand(sims, 1);
    
    delta_tilde_a = a - beta * Q_tilde(:,it) - Z(:,it);
    delta_tilde_b = b + beta * Q_tilde(:,it) + Z(:,it);
    tildes_a(:,it) = delta_tilde_a;
    tildes_b(:,it) = delta_tilde_b;
    
    prob_a = min(exp(-kappa*(delta_a - delta_tilde_a)), 1);
    prob_b = min(exp(-kappa*(delta_b - delta_tilde_b)), 1);
    
    filled_trade_a = (fill_a <= prob_a) & (Q(:,it) ~= (qmin+1));
    filled_trade_b = (fill_b <= prob_b) & (Q(:,it) ~= (qmax-1));
    
    % Update inventories
    Q(:,it+1) = Q(:,it) - filled_trade_a.*arrival_a + filled_trade_b.*arrival_b;
    Q_tilde(:,it+1) = Q_tilde(:,it) - (1-filled_trade_a).*arrival_a + (1-filled_trade_b).*arrival_b;
    
    % Update cash
    X(:,it+1) = X(:,it) + filled_trade_a.*arrival_a.*(S(:,it) + delta_a) ...
                - filled_trade_b.*arrival_b.*(S(:,it) - delta_b);
end

% Calculate final PnL and objective
pnl_euler = X(:,end) + Q(:,end).*(S(:,end) + (a-b)/2 - beta*Q_tilde(:,end) - Z(:,end));
for s = 1:sims
    intQ(s) = sum(Q(s,:).^2)/Nt;
end
objective_euler = pnl_euler - a*Q(:,end).^2 - phi*intQ;

% Copy last values for plotting
deltas_a(:,end) = deltas_a(:,end-1);
deltas_b(:,end) = deltas_b(:,end-1);
tildes_a(:,end) = tildes_a(:,end-1);
tildes_b(:,end) = tildes_b(:,end-1);

% Display results
fprintf('Terminal Q mean = %.4f\tTerminal Q sd = %.4f\n', mean(Q(:,end)), std(Q(:,end)));
fprintf('Terminal $\\tilde{Q}$ mean = %.4f\tTerminal $\\tilde{Q}$ sd = %.4f\n', mean(Q_tilde(:,end)), std(Q_tilde(:,end)));
fprintf('Terminal X mean = %.4f\tTerminal X sd = %.4f\n', mean(X(:,end)), std(X(:,end)));
fprintf('Terminal PnL mean = %.4f\tTerminal PnL sd = %.4f\n', mean(pnl_euler), std(pnl_euler));
fprintf('Terminal objective mean = %.4f\tTerminal objective sd = %.4f\n', mean(objective_euler), std(objective_euler));

% Plot results
nmin = 45;
nmax = 48;
figure;
subplot(4,2,1);
plot(time, Q(nmin:nmax,:)');
hold on;
plot(time, repmat(qmax,1,length(time)), '--');
plot(time, repmat(qmin,1,length(time)), '--');
ylabel('$Q_t$', 'Interpreter', 'latex');
grid on;

subplot(4,2,2);
plot(time, Q_tilde(nmin:nmax,:)');
ylabel('$\tilde{Q}_t$', 'Interpreter', 'latex');
grid on;

% Plot delta a and tilde delta a
subplot(4,2,3);
plot(time, deltas_a(nmin:nmax,:)');
ylabel('$\delta^{*,a}_t$', 'Interpreter', 'latex');
grid on;
alpha(0.3);

subplot(4,2,4);
plot(time, tildes_a(nmin:nmax,:)');
ylabel('$\tilde{\delta}^{a}_t$', 'Interpreter', 'latex');
grid on;
alpha(0.3);

% Plot delta b and tilde delta b
subplot(4,2,5);
plot(time, deltas_b(nmin:nmax,:)');
ylabel('$\delta^{*,b}_t$', 'Interpreter', 'latex');
grid on;
alpha(0.3);

subplot(4,2,6);
plot(time, tildes_b(nmin:nmax,:)');
ylabel('$\tilde{\delta}^{b}_t$', 'Interpreter', 'latex');
grid on;
alpha(0.3);

% Plot spreads
subplot(4,2,7);
plot(time, deltas_a(nmin:nmax,:)' + deltas_b(nmin:nmax,:)');
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\delta^{*,a}_t+\delta^{*,b}_t$', 'Interpreter', 'latex');
grid on;
alpha(0.3);

subplot(4,2,8);
plot(time, tildes_a(nmin:nmax,:)' + tildes_b(nmin:nmax,:)');
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\tilde{\delta}^{a}_t+\tilde{\delta}^{b}_t$', 'Interpreter', 'latex');
grid on;
alpha(0.3);

% Adjust layout
set(gcf, 'Color', 'w'); % White background
set(gcf, 'Position', [100 100 900 900]); % Set figure size

% Add legends where needed
subplot(4,2,3);
legend('Sim 1', 'Sim 2', 'Sim 3', 'Location', 'best');

subplot(4,2,5);
legend('Sim 1', 'Sim 2', 'Sim 3', 'Location', 'best');

% Tight layout
h = findobj(gcf, 'type', 'axes');
set(h, 'FontSize', 12); % Consistent font size
