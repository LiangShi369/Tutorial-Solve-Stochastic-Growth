clear

load sg_golden.mat

% Set simulation parameters
T = 100;          % Number of periods to simulate
burn_in = 50;     % Number of initial periods to discard
seed = 11;        % Random seed for reproducibility
rng(seed);

% Initialize variables
k_sim = zeros(T + burn_in, 1);   % Simulated capital path
a_sim = zeros(T + burn_in, 1);   % Simulated productivity shocks
c_sim = zeros(T + burn_in, 1);   % Simulated consumption path

% Initial conditions
k_sim(1) = mean(gridk);          % Start at the middle of the capital grid
a_sim(1) = ceil(na / 2);         % Start at the median productivity state

% Precompute cumulative distribution function (CDF) for shocks
cdf_a = cumsum(pdfa, 2);

% Simulation loop
for t = 1:(T + burn_in - 1)
   
    u = rand; % Draw random uniform number
    a_sim(t + 1) = find(cdf_a(a_sim(t), :) > u, 1) ;
    
    % Interpolate policy function for next-period capital
    k_next = interp1(gridk, kpol(:, a_sim(t)), k_sim(t), 'spline', 'extrap') ;
    k_sim(t + 1) = max(min(k_next, max(gridk)), min(gridk)); % Ensure within bounds
    
    % Compute consumption
    c_sim(t) = grida(a_sim(t)) * k_sim(t)^alpha + (1 - Del) * k_sim(t) - k_sim(t + 1) ;
end

% Discard burn-in periods
k_sim = k_sim((burn_in + 1):end-1);
a_sim = a_sim((burn_in + 1):end-1);
c_sim = c_sim((burn_in + 1):end-1);

% Map productivity states to actual values
a_sim_values = grida(a_sim);

%%

figure;
subplot(3, 1, 1);
plot(k_sim, 'LineWidth', 1.25);
xlabel('Time');
ylabel('Capital (k)');
title('Simulated Capital Path');

subplot(3, 1, 2);
plot(c_sim, 'LineWidth', 1.25);
xlabel('Time');
ylabel('Consumption (c)');
title('Simulated Consumption Path');

subplot(3, 1, 3);
plot(a_sim_values, 'LineWidth', 1.25);
xlabel('Time');
ylabel('Productivity Shock (A)');
title('Simulated Productivity Shocks');

%%
save('sg_simulation.mat', 'k_sim', 'c_sim', 'a_sim_values');

