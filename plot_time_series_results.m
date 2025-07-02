function plot_time_series_results(time, Q, Q_tilde, deltas_a, deltas_b, tildes_a, tildes_b, qmax, qmin)
    % PLOT_RESULTS Generate all visualization plots
    % Inputs:
    %   time - Time vector
    %   Q, Q_tilde - Inventory trajectories
    %   deltas_a/b, tildes_a/b - Quote trajectories
    %   qmax, qmin - Inventory bounds
    
    nmin = 50;
    nmax = 52;
    
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [100 100 900 900]);
    
    % Inventory plots
    subplot(4,2,1);
    plot(time, Q(nmin:nmax,:)');
    hold on;
    plot(time, repmat(qmax,1,length(time)), '--');
    plot(time, repmat(qmin,1,length(time)), '--');
    ylabel('$Q_t$', 'Interpreter', 'latex');
    grid on;
    title('Follower Inventory');
    
    subplot(4,2,2);
    plot(time, Q_tilde(nmin:nmax,:)');
    ylabel('$\tilde{Q}_t$', 'Interpreter', 'latex');
    grid on;
    title('Leader Inventory');
    
    % Ask quotes
    subplot(4,2,3);
    plot(time, deltas_a(nmin:nmax,:)');
    ylabel('$\delta^{*,a}_t$', 'Interpreter', 'latex');
    grid on;
    alpha(0.3);
    title('Follower Ask Quotes');
    legend('Sim 1', 'Sim 2', 'Sim 3', 'Location', 'best');
    
    subplot(4,2,4);
    plot(time, tildes_a(nmin:nmax,:)');
    ylabel('$\tilde{\delta}^{a}_t$', 'Interpreter', 'latex');
    grid on;
    alpha(0.3);
    title('Leader Ask Quotes');
    
    % Bid quotes
    subplot(4,2,5);
    plot(time, deltas_b(nmin:nmax,:)');
    ylabel('$\delta^{*,b}_t$', 'Interpreter', 'latex');
    grid on;
    alpha(0.3);
    title('Follower Bid Quotes');
    legend('Sim 1', 'Sim 2', 'Sim 3', 'Location', 'best');
    
    subplot(4,2,6);
    plot(time, tildes_b(nmin:nmax,:)');
    ylabel('$\tilde{\delta}^{b}_t$', 'Interpreter', 'latex');
    grid on;
    alpha(0.3);
    title('Leader Bid Quotes');
    
    % Spreads
    subplot(4,2,7);
    plot(time, deltas_a(nmin:nmax,:)' + deltas_b(nmin:nmax,:)');
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$\delta^{*,a}_t+\delta^{*,b}_t$', 'Interpreter', 'latex');
    grid on;
    alpha(0.3);
    title('Follower Spread');
    
    subplot(4,2,8);
    plot(time, tildes_a(nmin:nmax,:)' + tildes_b(nmin:nmax,:)');
    xlabel('$t$', 'Interpreter', 'latex');
    ylabel('$\tilde{\delta}^{a}_t+\tilde{\delta}^{b}_t$', 'Interpreter', 'latex');
    grid on;
    alpha(0.3);
    title('Leader Spread');
    
    % Adjust all axes
    h = findobj(gcf, 'type', 'axes');
    set(h, 'FontSize', 12);

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
end