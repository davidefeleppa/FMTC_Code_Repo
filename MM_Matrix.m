function [time, Q, Q_tilde, deltas_a, deltas_b, delta_tildes_a, delta_tildes_b, qmax, qmin, X, X_tilde, pnl, pnl_tilde, obj_follower, obj_leader] = MM_Matrix(a_func, b_func, beta, theta, phi_float, gamma_float, sims)
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
    phi = 0.02; %0.1
    gamma = 0.01; % 0.03
    
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

    % Main simulation loop
    % For each time step (i)
    for i = 1:steps-1
        % Current time step
        t = time(i);

        % Adjust arrival probability to spread
        lambda_a = lambda_0 * exp(-kappa * a_func(t));
        lambda_b = lambda_0 * exp(-kappa * b_func(t));
        
        % Calculate optimal spreads of the follower
        [delta_a, delta_b] = calculate_deltas(t, T, Q(:,i), Q_tilde(:,i), qmax, kappa, phi_float, gamma_float, a_func(t), b_func(t), lambda_a, lambda_b, beta, theta, true);
        deltas_a(:, i) = delta_a;
        deltas_b(:, i) = delta_b;
        
        % Update S 
        S(:,i+1) = S(:,i) + sigma * sqrt(dt) * randn(sims, 1);
        
        % Simulate order arrivals
        [arrival_a, arrival_b] = get_arrival(sims, dt, lambda_a, lambda_b);
        
        % Calculate fill probabilities
        fill_a = rand(sims, 1);
        fill_b = rand(sims, 1);
        
        % Get leader bid ask quotes
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

function [arrival_a, arrival_b] = get_arrival(sims, dt, lambda_a, lambda_b)
    arrival_a = poissrnd(lambda_a*dt, sims, 1);
    arrival_b = poissrnd(lambda_b*dt, sims, 1);
    arrival_a(arrival_a > 1) = 1;
    arrival_b(arrival_b > 1) = 1;
end

function omega_t = calculate_omega_t(t, qmax, phi, kappa, theta, lambda_a, lambda_b, beta, a, b, gamma, T)
    % Initialize Matrices 
    A_matrix = zeros(2*qmax+1, 2*qmax+1);
    B_matrix = zeros(2*qmax+1, 1);
    
    % For each inventory level [qmax:-1:qmin]
    for i = 0:(2*qmax) 
        inventory = qmax - i;
        
        % Diagonal elements
        A_matrix(i+1, i+1) = -phi * kappa * inventory^2 + beta * kappa * (lambda_a - lambda_b) * inventory;
        
        % Vector elements
        B_matrix(i+1) = exp(kappa * ((a-b)/2)*inventory - (gamma + 0.5*(theta-beta))*inventory^2);
        
        % Upper diagonal (i+1 terms)
        if i < 2*qmax
            A_matrix(i+1, i+2) = lambda_b * exp(-1 + kappa*b - kappa*0.5*(theta+beta));
        end
        
        % Lower diagonal (i-1 terms)
        if i > 0
            A_matrix(i+1, i) = lambda_a * exp(-1 + kappa*a - kappa*0.5*(theta+beta));
        end
    end
    
    % Matrix exponential and multiplication
    omega_t = expm(A_matrix*(T-t)) * B_matrix;

    if (sum(omega_t < 0) > 0)
        omega_t = 0;
        warning("Omega negative for given parameters");
    end

end

function gt = calculate_gt(t, kappa, theta, qmax, phi, lambda_a, lambda_b, beta, a, b, gamma, T)
    % Calculate omega(t)
    omega_function = calculate_omega_t(t, qmax, phi, kappa, theta, lambda_a, lambda_b, beta, a, b, gamma, T);
    
    % Compute g(t)
    gt = (1 / kappa) * log(omega_function);
end

function [delta_a, delta_b] = calculate_deltas(t, T, q, q_tilde, qmax, kappa, phi, gamma, a, b, lambda_a, lambda_b, beta, theta, take_maximum)
    % Calculate g(t) for each inventory level
    g_q = calculate_gt(t, kappa, theta, qmax, phi, lambda_a, lambda_b, beta, a, b, gamma, T);
    
    % % Optionally print the predited pnl value V
    % if (t==0)
    %     fprintf('V = %.3f \n', g_q(11))
    % end

    % Convert current inventory to 1-based index
    idx = qmax - q + 1;
    
    % Symmetric boundary handling
    if q == -qmax
        % Minimum inventory (max short position): can only increase inventory
        c_hat_a = 1/kappa - g_q(idx) + g_q(idx);  % No change for asks
        c_hat_b = 1/kappa - g_q(idx-1) + g_q(idx); % Can only buy
    elseif q == qmax
        % Maximum inventory (max long position): can only decrease inventory
        c_hat_a = 1/kappa - g_q(idx+1) + g_q(idx); % Can only sell
        c_hat_b = 1/kappa - g_q(idx) + g_q(idx);    % No change for bids
    else
        % Normal case
        c_hat_a = 1/kappa - g_q(idx+1) + g_q(idx);
        c_hat_b = 1/kappa - g_q(idx-1) + g_q(idx);
    end
    
    % Apply maximum constraints if needed
    if take_maximum
        c_a = max(c_hat_a, a - 0.5*(theta+beta));
        c_b = max(c_hat_b, b - 0.5*(theta+beta));
    else
        c_a = c_hat_a;
        c_b = c_hat_b;
    end
    
    % Final delta calculations
    delta_a = c_a - beta*q_tilde + 0.5*(theta+beta) - theta*q;
    delta_b = c_b + beta*q_tilde + 0.5*(theta+beta) + theta*q;
end
