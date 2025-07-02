function [time, Q, Q_tilde, deltas_a, deltas_b, tildes_a, tildes_b, qmax, qmin, X, X_tilde, pnl, pnl_tilde, obj_follower, obj_leader] = MM_Euler(A, B, beta, sims)    
    % Stock Price
    sigma = 1;
    S0 = 100;
    
    % Arrival Probabilities ~ Poission(λ)
    lambda_a = 10;
    lambda_b = 10; 
    
    % Inventory
    qmax = 10;
    qmin = -qmax;

    % Fill probability sensitivity 
    kappa = 2; 

    % Penalties
    phi = 0.01;
    gamma = 0.03;
    
    % All possible inventory levels 
    all_inventories = qmax:-1:qmin-1;
    
    % Time scale
    T = 1;
    Nt = 1000; % Number of steps
    dt = T/Nt; % Time step size
    euler_steps = 1000000; % Euler steps
    time = linspace(0, T, Nt+1); % Time axis
    
    % Initialize arrays
    X = zeros(sims, Nt+1); % Cash
    X_tilde = zeros(sims, Nt+1); % Cash
    Q = zeros(sims, Nt+1); % Inventory
    Q_tilde = zeros(sims, Nt+1); % Competitors Inventory
    S = zeros(sims, Nt+1); % Stock price path
    S(:,1) = S0; % Start stock at S0

    % Bid ask quotes of us and the competitor
    deltas_a = zeros(sims, Nt+1);
    deltas_b = zeros(sims, Nt+1);
    tildes_a = zeros(sims, Nt+1);
    tildes_b = zeros(sims, Nt+1);
    
    % Arrival of the asks and bids
    arrivaltimes_a = zeros(sims, Nt+1);
    arrivaltimes_b = zeros(sims, Nt+1);
    
    % Run Euler scheme to calculate g as a function of time step and inventory
    % g(t,q)
    g = euler_scheme(euler_steps, T, qmax, qmin, A, B, gamma, beta, kappa, phi, lambda_a, lambda_b);
    
    % Round all_inventories to the cent
    all_inventories_list = round(all_inventories, 2);
    
    % Main simulation loop
    % For each time step (i)
    for i = 1:Nt
        t = time(i);
        
        % Get current g values for this time step
        gt_all_q = g(round(i*(euler_steps/Nt)), :);
    
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
        [delta_a, delta_b] = calculate_deltas_euler(Q_tilde(:,i), G, G_down, G_up, ...
                             kappa, A(i), B(i), beta, true);
        % Save deltas
        deltas_a(:, i) = delta_a;
        deltas_b(:, i) = delta_b;
        
        % Update Z and S (Brownian motions)
        S(:,i+1) = S(:,i) + sigma * sqrt(dt) * randn(sims, 1);
        
        % Simulate order arrivals
        [arrival_a, arrival_b] = get_arrival(sims, dt, lambda_a, lambda_b);
        arrivaltimes_a(:,i) = arrival_a * t;
        arrivaltimes_b(:,i) = arrival_b * t;
        
        % Calculate fill probabilities
        fill_a = rand(sims, 1);
        fill_b = rand(sims, 1);
        
        % Get competitor bid ask quotes
        delta_tilde_a = A(i) - beta * Q_tilde(:,i);
        delta_tilde_b = B(i) + beta * Q_tilde(:,i);
        tildes_a(:,i) = delta_tilde_a;
        tildes_b(:,i) = delta_tilde_b;
        
        % Calculate prob of filling orders and and respective binary operator 
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
    
    % Calculate final PnL and objectives
    % For follower
    final_price_follower = S(:,end) - beta*Q_tilde(:,end);
    pnl = X(:,end) + Q(:,end).*final_price_follower;
    intQ = sum(Q.^2, 2)*dt/T; 
    obj_follower = pnl - A(end)*Q(:,end).^2 - phi*intQ;
    
    % For leader
    final_price_leader = S(:,end) - beta*Q_tilde(:,end);
    pnl_tilde = X_tilde(:,end) + Q_tilde(:,end).*final_price_leader;
    intQ_tilde = sum(Q_tilde.^2, 2)*dt/T;
    obj_leader = pnl_tilde - A(end)*Q_tilde(:,end).^2 - phi*intQ_tilde;
    
end

function g = euler_scheme(steps, T, qmax, qmin, A, B, gamma, beta, kappa, phi, lambda_a, lambda_b)
    % Initialize variables
    h = T / steps;
    all_inventories = qmax:-1:qmin; 
    g = zeros(steps, length(all_inventories));
    
    % Terminal condition
    for iq = 1:length(all_inventories)
        q = all_inventories(iq);
        g(end, iq) = ((A(end)-B(end))/2)*q - (gamma-beta/2)*(q^2);
    end
    
    % Time stepping (backward)
    for i = 1:(steps-1)
        t_idx = steps - i + 1;  % Current time index in backward scheme
        a = A(t_idx);
        b = B(t_idx);
        
        for iq = 1:length(all_inventories)
            q = all_inventories(iq);
            
            % Interior points
            if iq ~= 1 && iq ~= length(all_inventories)
                c_star_a = max(1/kappa - g(t_idx, iq-1) + g(t_idx, iq), a - beta/2);
                c_star_b = max(1/kappa - g(t_idx, iq+1) + g(t_idx, iq), b - beta/2);
                
                g(t_idx-1, iq) = g(t_idx, iq) + h * ( ...
                    -phi*(q^2) ...
                    + (lambda_a-lambda_b)*beta*q ...
                    + lambda_a*exp(-kappa*(c_star_a + beta/2 - a))*(c_star_a + g(t_idx,iq-1) - g(t_idx,iq)) ...
                    + lambda_b*exp(-kappa*(c_star_b + beta/2 - b))*(c_star_b + g(t_idx,iq+1) - g(t_idx,iq)));
            
            % Lower boundary (iq == 1)
            elseif iq == 1
                c_star_b = max(1/kappa - g(t_idx, iq+1) + g(t_idx, iq), b - beta/2);
                
                g(t_idx-1, iq) = g(t_idx, iq) + h * ( ...
                    -phi*(q^2) ...
                    + (lambda_a-lambda_b)*beta*q ...
                    + lambda_b*exp(-kappa*(c_star_b + beta/2 - b))*(c_star_b + g(t_idx,iq+1) - g(t_idx,iq)));
            
            % Upper boundary (iq == end)
            else
                c_star_a = max(1/kappa - g(t_idx, iq-1) + g(t_idx, iq), a - beta/2);
                
                g(t_idx-1, iq) = g(t_idx, iq) + h * ( ...
                    -phi*(q^2) ...
                    + (lambda_a-lambda_b)*beta*q ...
                    + lambda_a*exp(-kappa*(c_star_a + beta/2 - a))*(c_star_a + g(t_idx,iq-1) - g(t_idx,iq)));
            end
        end
    end
end

function [delta_a, delta_b] = calculate_deltas_euler(q_tilde, gt, gt_q_minus, gt_q_plus, kappa, a, b, beta, take_maximum)
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
    delta_a = c_a - beta*q_tilde + beta/2;
    delta_b = c_b + beta*q_tilde + beta/2;
end

function [arrival_a, arrival_b] = get_arrival(sims, dt, lambda_a, lambda_b)
    % Generate uniform random numbers
    unif_a = rand(sims, 1);
    unif_b = rand(sims, 1);
    
    % Determine arrivals using Poisson process probabilities
    arrival_a = unif_a < (1 - exp(-lambda_a * dt));
    arrival_b = unif_b < (1 - exp(-lambda_b * dt));
end