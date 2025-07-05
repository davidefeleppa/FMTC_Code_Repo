function plot_time_series_results(time, Q, Q_tilde, deltas_a, deltas_b, tildes_a, tildes_b, qmax, qmin)
    % PLOT_RESULTS Generate all visualization plots with trade markers
    % Inputs:
    %   time - Time vector
    %   Q, Q_tilde - Inventory trajectories
    %   deltas_a/b, tildes_a/b - Quote trajectories
    %   qmax, qmin - Inventory bounds
    
    % Detect trades from inventory changes 
    filled_trade_a = [false(size(Q,1),1), diff(Q,1,2) < 0];
    filled_trade_b = [false(size(Q,1),1), diff(Q,1,2) > 0];
    filled_trade_a_tilde = [false(size(Q_tilde,1),1), diff(Q_tilde,1,2) < 0];
    filled_trade_b_tilde = [false(size(Q_tilde,1),1), diff(Q_tilde,1,2) > 0];
    
    nmin = 1;
    disp_number = 1;
    nmax = nmin + disp_number - 1;
    
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [100 100 900 900]);
    
    % Get trade times for follower and leader
    trade_times_follower_a = time(any(filled_trade_a(nmin:nmax,:), 1));
    trade_times_follower_b = time(any(filled_trade_b(nmin:nmax,:), 1));
    trade_times_leader_a = time(any(filled_trade_a_tilde(nmin:nmax,:), 1));
    trade_times_leader_b = time(any(filled_trade_b_tilde(nmin:nmax,:), 1));
    
    % Convert trade times to column vectors
    trade_times_follower_a = trade_times_follower_a(:);
    trade_times_follower_b = trade_times_follower_b(:);
    trade_times_leader_a = trade_times_leader_a(:);
    trade_times_leader_b = trade_times_leader_b(:);
    
    % Inventory plots with trade markers
    subplot(4,2,1);
    plot(time, Q(nmin:nmax,:)');
    hold on;
    plot(time, repmat(qmax,1,length(time)), '--');
    plot(time, repmat(qmin,1,length(time)), '--');
    % Mark ask trades (inventory decreases)
    if ~isempty(trade_times_follower_a)
        scatter(trade_times_follower_a, interp1(time, Q(nmin,:), trade_times_follower_a), 'ro', 'filled');
    end
    % Mark bid trades (inventory increases)
    if ~isempty(trade_times_follower_b)
        scatter(trade_times_follower_b, interp1(time, Q(nmin,:), trade_times_follower_b), 'go', 'filled');
    end
    ylabel('$Q_t$', 'Interpreter', 'latex');
    grid on;
    ylim([qmin-1, qmax+1])  % Sets y-axis limits
    title('Follower Inventory');
    legend('Inventory', 'Max', 'Min', 'Ask Trade', 'Bid Trade', 'Location', 'best');
    
    subplot(4,2,2);
    plot(time, Q_tilde(nmin:nmax,:)');
    hold on;
    % Mark ask trades (inventory decreases)
    if ~isempty(trade_times_leader_a)
        scatter(trade_times_leader_a, interp1(time, Q_tilde(nmin,:), trade_times_leader_a), 'ro', 'filled');
    end
    % Mark bid trades (inventory increases)
    if ~isempty(trade_times_leader_b)
        scatter(trade_times_leader_b, interp1(time, Q_tilde(nmin,:), trade_times_leader_b), 'go', 'filled');
    end
    ylabel('$\tilde{Q}_t$', 'Interpreter', 'latex');
    grid on;
    title('Leader Inventory');
    legend('Inventory', 'Ask Trade', 'Bid Trade', 'Location', 'best');
    
    % Ask quotes with trade markers
    subplot(4,2,3);
    plot(time, deltas_a(nmin:nmax,:)');
    hold on;
    % Mark ask trades
    if ~isempty(trade_times_follower_a)
        scatter(trade_times_follower_a, interp1(time, deltas_a(nmin,:), trade_times_follower_a), 'ro', 'filled');
    end
    ylabel('$\delta^{*,a}_t$', 'Interpreter', 'latex');
    grid on;
    title('Follower Ask Quotes');
    legend('Quotes', 'Trades', 'Location', 'best');
    
    subplot(4,2,4);
    plot(time, tildes_a(nmin:nmax,:)');
    hold on;
    % Mark ask trades
    if ~isempty(trade_times_leader_a)
        scatter(trade_times_leader_a, interp1(time, tildes_a(nmin,:), trade_times_leader_a), 'ro', 'filled');
    end
    ylabel('$\tilde{\delta}^{a}_t$', 'Interpreter', 'latex');
    grid on;
    title('Leader Ask Quotes');
    legend('Quotes', 'Trades', 'Location', 'best');
    
    % Bid quotes with trade markers
    subplot(4,2,5);
    plot(time, deltas_b(nmin:nmax,:)');
    hold on;
    % Mark bid trades
    if ~isempty(trade_times_follower_b)
        scatter(trade_times_follower_b, interp1(time, deltas_b(nmin,:), trade_times_follower_b), 'go', 'filled');
    end
    ylabel('$\delta^{*,b}_t$', 'Interpreter', 'latex');
    grid on;
    title('Follower Bid Quotes');
    legend('Quotes', 'Trades', 'Location', 'best');
    
    subplot(4,2,6);
    plot(time, tildes_b(nmin:nmax,:)');
    hold on;
    % Mark bid trades
    if ~isempty(trade_times_leader_b)
        scatter(trade_times_leader_b, interp1(time, tildes_b(nmin,:), trade_times_leader_b), 'go', 'filled');
    end
    ylabel('$\tilde{\delta}^{b}_t$', 'Interpreter', 'latex');
    grid on;
    title('Leader Bid Quotes');
    legend('Quotes', 'Trades', 'Location', 'best');
    
    % Spreads with trade markers
    subplot(4,2,7);
    spread_follower = deltas_a(nmin:nmax,:)' + deltas_b(nmin:nmax,:)';
    plot(time, spread_follower);
    hold on;
    % Mark all trades
    if ~isempty(trade_times_follower_a) || ~isempty(trade_times_follower_b)
        all_times = [trade_times_follower_a; trade_times_follower_b];
        scatter(all_times, interp1(time, spread_follower, all_times), 'ko', 'filled');
    end
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$\delta^{*,a}_t+\delta^{*,b}_t$', 'Interpreter', 'latex');
    grid on;
    title('Follower Spread');
    legend('Spread', 'Trades', 'Location', 'best');
    
    subplot(4,2,8);
    spread_leader = tildes_a(nmin:nmax,:)' + tildes_b(nmin:nmax,:)';
    plot(time, spread_leader);
    hold on;
    % Mark all trades
    if ~isempty(trade_times_leader_a) || ~isempty(trade_times_leader_b)
        all_times = [trade_times_leader_a; trade_times_leader_b];
        scatter(all_times, interp1(time, spread_leader, all_times), 'ko', 'filled');
    end
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$\tilde{\delta}^{a}_t+\tilde{\delta}^{b}_t$', 'Interpreter', 'latex');
    grid on;
    title('Leader Spread');
    legend('Spread', 'Trades', 'Location', 'best');
    
    % Adjust all axes
    h = findobj(gcf, 'type', 'axes');
    set(h, 'FontSize', 12);
end