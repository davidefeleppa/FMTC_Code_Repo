function plot_time_series_results(time, Q_F, Q_L, delta_a_F, delta_b_F, delta_a_L, delta_b_L, qmax, qmin)
    % PLOT_RESULTS Generate combined visualization plots with trade markers
    % Inputs:
    %   time - Time vector
    %   Q_F, Q_L - Follower and Leader inventory trajectories
    %   delta_a/b_F, delta_a/b_L - Follower and Leader quote trajectories
    %   qmax, qmin - Inventory bounds
    
    % Detect trades from inventory changes 
    filled_trade_a_F = [false(size(Q_F,1),1), diff(Q_F,1,2) < 0];
    filled_trade_b_F = [false(size(Q_F,1),1), diff(Q_F,1,2) > 0];
    filled_trade_a_L = [false(size(Q_L,1),1), diff(Q_L,1,2) < 0];
    filled_trade_b_L = [false(size(Q_L,1),1), diff(Q_L,1,2) > 0];
    
    nmin = 1;
    disp_number = 1;
    nmax = nmin + disp_number - 1;
    
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [100 100 900 450]); % Adjusted for fewer subplots
    
    % Get trade times (properly closed brackets)
    trade_times_a_F = time(any(filled_trade_a_F(nmin:nmax,:), 1));
    trade_times_b_F = time(any(filled_trade_b_F(nmin:nmax,:), 1));
    trade_times_a_L = time(any(filled_trade_a_L(nmin:nmax,:), 1));
    trade_times_b_L = time(any(filled_trade_b_L(nmin:nmax,:), 1));
    
    % Convert to column vectors
    trade_times_a_F = trade_times_a_F(:);
    trade_times_b_F = trade_times_b_F(:);
    trade_times_a_L = trade_times_a_L(:);
    trade_times_b_L = trade_times_b_L(:);
    
    % Set colors and markers
    F_color = [0 0.4470 0.7410]; % Blue for Follower
    L_color = [0.8500 0.3250 0.0980]; % Orange for Leader
    
    % Combined Inventory Plot (properly closed brackets) %%%
    subplot(2,2,1);
    plot(time, Q_F(nmin:nmax,:)', 'Color', F_color, 'LineWidth', 1.5);
    hold on;
    plot(time, Q_L(nmin:nmax,:)', 'Color', L_color, 'LineWidth', 1.5);
    
    % Mark trades (properly closed brackets)
    if ~isempty(trade_times_a_F)
        scatter(trade_times_a_F, interp1(time, Q_F(nmin,:), trade_times_a_F), 'ro', 'filled');
    end
    if ~isempty(trade_times_b_F)
        scatter(trade_times_b_F, interp1(time, Q_F(nmin,:), trade_times_b_F), 'go', 'filled');
    end
    if ~isempty(trade_times_a_L)
        scatter(trade_times_a_L, interp1(time, Q_L(nmin,:), trade_times_a_L), 'ro', 'filled');
    end
    if ~isempty(trade_times_b_L)
        scatter(trade_times_b_L, interp1(time, Q_L(nmin,:), trade_times_b_L), 'go', 'filled');
    end

    plot(time, repmat(qmax,1,length(time)), 'k--');
    plot(time, repmat(qmin,1,length(time)), 'k--');
    
    ylabel('Inventory $Q_t$', 'Interpreter', 'latex');
    grid on;
    ylim([qmin-1, qmax+1]);
    title('Inventory Trajectories');
    legend('Follower', 'Leader', 'Ask', 'Bid', 'Location', 'best');

    % Combined Spreads (properly closed brackets) %%%
    subplot(2,2,2);
    spread_F = delta_a_F(nmin:nmax,:)' + delta_b_F(nmin:nmax,:)';
    spread_L = delta_a_L(nmin:nmax,:)' + delta_b_L(nmin:nmax,:)';
    
    plot(time, spread_F, 'Color', F_color, 'LineWidth', 1.5);
    hold on;
    plot(time, spread_L, 'Color', L_color, 'LineWidth', 1.5);
    
    % Mark all trades (properly closed brackets)
    all_times_F = [trade_times_a_F; trade_times_b_F];
    all_times_L = [trade_times_a_L; trade_times_b_L];
    
    if ~isempty(all_times_F)
        scatter(all_times_F, interp1(time, spread_F, all_times_F), 'ko', 'filled');
    end
    if ~isempty(all_times_L)
        scatter(all_times_L, interp1(time, spread_L, all_times_L), 'ko', 'filled');
    end
    
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('Spread $\delta^a_t+\delta^b_t$', 'Interpreter', 'latex');
    grid on;
    title('Spread Trajectories');
    legend('Follower', 'Leader', 'Trades', 'Location', 'best');
    
    % Combined Bid Quotes (properly closed brackets) %%%
    subplot(2,2,3);
    plot(time, delta_b_F(nmin:nmax,:)', 'Color', F_color, 'LineWidth', 1.5);
    hold on;
    plot(time, delta_b_L(nmin:nmax,:)', 'Color', L_color, 'LineWidth', 1.5);
    
    if ~isempty(trade_times_b_F)
        scatter(trade_times_b_F, interp1(time, delta_b_F(nmin,:), trade_times_b_F), 'go', 'filled');
    end
    if ~isempty(trade_times_b_L)
        scatter(trade_times_b_L, interp1(time, delta_b_L(nmin,:), trade_times_b_L), 'go', 'filled');
    end
    
    ylabel('Bid Quotes $\delta^b_t$', 'Interpreter', 'latex');
    grid on;
    title('Bid Quote Trajectories');
    legend('Follower', 'Leader', 'Trades', 'Location', 'best');

    % Combined Ask Quotes (properly closed brackets) %%%
    subplot(2,2,4);
    plot(time, delta_a_F(nmin:nmax,:)', 'Color', F_color, 'LineWidth', 1.5);
    hold on;
    plot(time, delta_a_L(nmin:nmax,:)', 'Color', L_color, 'LineWidth', 1.5);
    
    if ~isempty(trade_times_a_F)
        scatter(trade_times_a_F, interp1(time, delta_a_F(nmin,:), trade_times_a_F), 'ro', 'filled');
    end
    if ~isempty(trade_times_a_L)
        scatter(trade_times_a_L, interp1(time, delta_a_L(nmin,:), trade_times_a_L), 'ro', 'filled');
    end
    
    ylabel('Ask Quotes $\delta^a_t$', 'Interpreter', 'latex');
    grid on;
    title('Ask Quote Trajectories');
    legend('Follower', 'Leader', 'Trades', 'Location', 'best');

    
    % Adjust all axes
    h = findobj(gcf, 'type', 'axes');
    set(h, 'FontSize', 12);
    set(h, 'TickLabelInterpreter', 'latex');
end