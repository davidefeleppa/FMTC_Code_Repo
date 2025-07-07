function [time, Q, Q_tilde, deltas_a, deltas_b, delta_tildes_a, delta_tildes_b, qmax, qmin, X, X_tilde, pnl, pnl_tilde, obj_follower, obj_leader] = MM_Euler(a_func, b_func, beta, theta, sims)    
    % Stock Price
    sigma = 1;
    S0 = 100;
    
    % Inventory
    qmax = 10;
    qmin = -qmax;

    % Fill probability sensitivity 
    kappa = 2; 

    % Arrival Probabilities Initial Intensity ~ Poission(λ)
    lambda_0 = 10;

    % Penalties
    phi = 0.1;
    gamma = 0.03;

    % All inventories
    all_inventories_list = qmax:-1:qmin;
    
    % Time scale
    T = 1;
    steps = 1000; % Number of time steps
    dt = T/steps; % Time step size
    time = linspace(0, T, steps); % Time axis
    
    % Initialize arrays
    X = zeros(sims, steps); % Cash
    X_tilde = zeros(sims, steps);
    Q = zeros(sims, steps); % Inventory
    Q_tilde = zeros(sims, steps);
    S = zeros(sims, steps); % Stock price path
    S(:,1) = S0; % Start stock at S0

    % Bid/ask spread
    deltas_a = zeros(sims, steps);
    deltas_b = zeros(sims, steps);
    delta_tildes_a = zeros(sims, steps);
    delta_tildes_b = zeros(sims, steps);
        
    % Order arrival probability intensities as functions of time 
    lambda_a_func = @(t) lambda_0 * exp(-kappa * a_func(t));
    lambda_b_func = @(t) lambda_0 * exp(-kappa * b_func(t));
    
    % Run Euler approximation to calculate G as a function of time and inventory: G(t,q)
    euler_steps = 1e6;
    G_approx = euler_approx(euler_steps, T, qmax, qmin, a_func, b_func, gamma, beta, kappa, phi, theta, lambda_a_func, lambda_b_func);
    
    % Optionally print the predited pnl value V
    %fprintf('V = %.3f \n', G_approx(1,11));
    
    % Main simulation loop
    % For each time step (i)
    for i = 1:steps-1
        % Current time step
        t = time(i);

        % Get current g values for this time step
        t_idx = floor(i*(euler_steps/steps));
        gt_all_q = G_approx(t_idx, :);
    
        % Pre-allocate space for G, G up and G down
        G = zeros(sims, 1);
        G_down = zeros(sims, 1);
        G_up = zeros(sims, 1);
        
        % For each simulation (sim)
        for sim = 1:sims
            % Get inventory level 
            q = Q(sim,i);

            % Get G for this inventory level
            [~, idx] = min(abs(all_inventories_list - q)); % Index of current inventory level
            G(sim) = gt_all_q(idx);
            
            % If inventory is not empty get G Down, for when we fill an ask 
            if q ~= qmin
                [~, idx_down] = min(abs(all_inventories_list - (q-1)));
                G_down(sim) = gt_all_q(idx_down);
            end
            
            % If inventory is not full get G Up, for when we fill a bid
            if q ~= qmax
                [~, idx_up] = min(abs(all_inventories_list - (q+1)));
                G_up(sim) = gt_all_q(idx_up);
            end
        end
        
        % Calculate deltas
        [delta_a, delta_b] = calculate_deltas_euler(Q(:,i), Q_tilde(:,i), G, G_down, G_up, ...
                             kappa, a_func(t), b_func(t), beta, theta, true);
        % Save deltas
        deltas_a(:, i) = delta_a;
        deltas_b(:, i) = delta_b;
        
        % Update S
        S(:,i+1) = S(:,i) + sigma * sqrt(dt) * randn(sims, 1);
        
        % Simulate order arrivals
        [arrival_a, arrival_b] = get_arrival(sims, dt, lambda_a_func(t), lambda_b_func(t));
        
        % Calculate fill probabilities
        fill_a = rand(sims, 1);
        fill_b = rand(sims, 1);
        
        % Get competitor bid ask quotes
        delta_tilde_a = a_func(t) - beta * Q_tilde(:,i) - theta*Q(:, i);
        delta_tilde_b = b_func(t) + beta * Q_tilde(:,i) + theta*Q(:, i);
        delta_tildes_a(:,i) = delta_tilde_a;
        delta_tildes_b(:,i) = delta_tilde_b;
        
        % Calculate prob of filling orders and respective binary operator 
        prob_a = min(exp(-kappa*(delta_a - delta_tilde_a)), 1);
        prob_b = min(exp(-kappa*(delta_b - delta_tilde_b)), 1);
        
        filled_trade_a = (fill_a <= prob_a) & (Q(:,i) ~= (qmin+1));
        filled_trade_b = (fill_b <= prob_b) & (Q(:,i) ~= (qmax-1));
        
        filled_trade_a_tilde = (fill_a <= 1-prob_a) & (Q_tilde(:,i) ~= (qmin+1));
        filled_trade_b_tilde = (fill_b <= 1-prob_b) & (Q_tilde(:,i) ~= (qmax-1));
        
        % Update inventories
        Q(:,i+1) = Q(:,i) - filled_trade_a.*arrival_a + filled_trade_b.*arrival_b;
        Q_tilde(:,i+1) = Q_tilde(:,i) - filled_trade_a_tilde.*arrival_a + (filled_trade_b_tilde).*arrival_b;
        
        % Update cashes
        X(:,i+1) = X(:,i) + filled_trade_a.*arrival_a.*(S(:,i) + delta_a) ...
                    - filled_trade_b.*arrival_b.*(S(:,i) - delta_b);
        
        X_tilde(:,i+1) = X_tilde(:,i) + filled_trade_a_tilde.*arrival_a.*(S(:,i) + delta_tilde_a) ...
                    - filled_trade_b_tilde.*arrival_b.*(S(:,i) - delta_tilde_b);
    end
    
    % Calculate final PnL and objectives
    
    % For follower
    final_price_follower = S(:,end) - beta*Q_tilde(:,end) - theta*Q(:,end) + 0.5*(a_func(T) - b_func(T));
    pnl = X(:,end) + Q(:,end).*final_price_follower;
    intQ = sum(Q.^2, 2)*dt/T; 
    obj_follower = pnl - gamma*Q(:,end).^2 - phi*intQ;
    
    % For leader
    final_price_leader = S(:,end);
    pnl_tilde = X_tilde(:,end) + Q_tilde(:,end).*final_price_leader;
    intQ_tilde = sum(Q_tilde.^2, 2)*dt/T;
    obj_leader = pnl_tilde - gamma*Q_tilde(:,end).^2 - phi*intQ_tilde;
    
end

function G = euler_approx(euler_steps, T, qmax, qmin, a_func, b_func, gamma, beta, kappa, phi, theta, lambda_a_func, lambda_b_func)
    % Initialize variables
    h = T / euler_steps;
    time = linspace(0, T, euler_steps); % Time axis
    all_inventories = qmax:-1:qmin; 
    num_inventories = length(all_inventories);
    
    % Initialize 3D matrix (time × Q × Q̃)
    G = zeros(euler_steps, num_inventories, num_inventories);
    
    % Terminal condition - G(T, Q, Q̃) doesn't depend on Q̃
    for iq = 1:num_inventories
        q = all_inventories(iq);
        terminal_value = ((a_func(T) - b_func(T))/2)*q - (gamma + 0.5*(theta-beta))*(q^2);
        G(end, iq, :) = terminal_value; % Same for all Q̃
    end
    
    % Time stepping (backward)
    for i = 1:(euler_steps-1)
        t_idx = euler_steps - i + 1;  % Current time index in backward scheme
        a = a_func(time(t_idx));
        b = b_func(time(t_idx));
        lambda_a = lambda_a_func(time(t_idx));
        lambda_b = lambda_b_func(time(t_idx));
        
        % Compute G(t-1, Q, Q̃) for all Q and Q̃
        % We have not yet implemented a time dependant beta/theta value so we only 
        % need G dependant on the Q inventory
        for iq_tilde = 1:1:1
            for iq = 1:num_inventories
                q = all_inventories(iq);
                
                % Interior points
                if iq ~= 1 && iq ~= num_inventories
                    c_star_a = max(1/kappa - G(t_idx, iq-1, iq_tilde) + G(t_idx, iq, iq_tilde), a - 0.5*(beta+theta));
                    c_star_b = max(1/kappa - G(t_idx, iq+1, iq_tilde) + G(t_idx, iq, iq_tilde), b - 0.5*(beta+theta));
                    
                    G(t_idx-1, iq, iq_tilde) = G(t_idx, iq, iq_tilde) + h * ( ...
                        -phi*(q^2) ...
                        + (lambda_a-lambda_b)*beta*q ...
                        + lambda_a*exp(-kappa*(c_star_a + 0.5*(theta+beta) - a))*(c_star_a + G(t_idx,iq-1,iq_tilde) - G(t_idx,iq,iq_tilde)) ...
                        + lambda_b*exp(-kappa*(c_star_b + 0.5*(theta+beta) - b))*(c_star_b + G(t_idx,iq+1,iq_tilde) - G(t_idx,iq,iq_tilde)));
                
                % Lower boundary (iq == 1)
                elseif iq == 1
                    c_star_b = max(1/kappa - G(t_idx, iq+1, iq_tilde) + G(t_idx, iq, iq_tilde), b - 0.5*(beta+theta));
                    
                    G(t_idx-1, iq, iq_tilde) = G(t_idx, iq, iq_tilde) + h * ( ...
                        -phi*(q^2) ...
                        + (lambda_a-lambda_b)*beta*q ...
                        + lambda_b*exp(-kappa*(c_star_b + 0.5*(beta+theta) - b))*(c_star_b + G(t_idx,iq+1,iq_tilde) - G(t_idx,iq,iq_tilde)));
                
                % Upper boundary (iq == end)
                else
                    c_star_a = max(1/kappa - G(t_idx, iq-1, iq_tilde) + G(t_idx, iq, iq_tilde), a - 0.5*(beta+theta));
                    
                    G(t_idx-1, iq, iq_tilde) = G(t_idx, iq, iq_tilde) + h * ( ...
                        -phi*(q^2) ...
                        + (lambda_a-lambda_b)*beta*q ...
                        + lambda_a*exp(-kappa*(c_star_a + 0.5*(beta+theta) - a))*(c_star_a + G(t_idx,iq-1,iq_tilde) - G(t_idx,iq,iq_tilde)));
                end
            end
        end
    end

    % We have not yet implemented a time dependant beta/theta value so we only
    % need G dependant on the Q inventory
    G = G(:,:,1);
end

function [delta_a, delta_b] = calculate_deltas_euler(q, q_tilde, gt, gt_q_minus, gt_q_plus, kappa, a, b, beta, theta, take_maximum)
    % Calculate c_hat values
    c_hat_a = 1/kappa - gt_q_minus + gt;
    c_hat_b = 1/kappa - gt_q_plus + gt;
    
    % Apply maximum condition if needed
    if take_maximum
        c_a = max(c_hat_a, a - 0.5*(beta+theta));
        c_b = max(c_hat_b, b - 0.5*(beta+theta));
    else
        c_a = c_hat_a;
        c_b = c_hat_b;
    end
    
    % Calculate final deltas
    delta_a = c_a - beta*q_tilde + 0.5*(beta+theta) - theta*q;
    delta_b = c_b + beta*q_tilde + 0.5*(beta+theta) + theta*q;
end

function [arrival_a, arrival_b] = get_arrival(sims, dt, lambda_a, lambda_b)
    % Generate uniform random numbers
    unif_a = rand(sims, 1);
    unif_b = rand(sims, 1);
    
    % Determine arrivals using Poisson process probabilities
    arrival_a = unif_a < (1 - exp(-lambda_a * dt));
    arrival_b = unif_b < (1 - exp(-lambda_b * dt));
end