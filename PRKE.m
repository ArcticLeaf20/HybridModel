% OSTR Hybrid Simulation Model, PRKE Solver Module
% Developed by Aidan Marsters
% Solves for power and precursor concentrations over time for a given
% reactivity insertion. Uses integrating factor approach to solve PRKE.

increment = 15 / 0.01; % Define the increment based on the time step
time_values = linspace(0, 15, increment); % Create an array of evenly spaced time values
% Parameters
generationTime = 1 * (10^-7); % Mean neutron generation time
decay_constants = [0.0128, 0.0301, 0.124, 0.325, 1.12, 2.69]; % Decay constants for each precursor group
beta_i = [0.000073, 0.000626, 0.000443, 0.000685, 0.000181, 0.000092]; % Beta values for precursor groups
beta = sum(beta_i); % Effective beta
concentration_values = zeros(1, 6); % Array for concentration values
power_values = zeros(1, length(time_values)); % Array for power values
initial_power = 3000 * 10^6; % Initial power in watts
reactivity = zeros(1, increment); % Array for reactivity values

% Initialize initial precursor concentrations
for j = 1:6
   concentration_values(j) = (beta_i(j) / (decay_constants(j) * generationTime)) * initial_power;
end
    
C_old = concentration_values;
time_step = 0.01;
    
% Main loop to calculate new concentration/power values
for k = 1:length(time_values)
    if time_values(k) < 5
       reactivity(k) = 0.08 * time_values(k) * beta;
    else
       reactivity(k) = (0.4 - 0.08 * (time_values(k) - 5)) * beta;
    end
        
    alpha = (reactivity(k) - beta) / generationTime;
        
    concentration_sum = 0;
    for x = 1:6
       concentration_sum = concentration_sum + decay_constants(x) * C_old(x);
    end
        
    power_values(k) = initial_power * exp(alpha * time_step) + ((1 / alpha) * (exp(alpha * time_step) - 1) * (concentration_sum));
        
    for z = 1:6
            concentration_values(z) = (C_old(z) * exp(-decay_constants(z) * time_step)) + ((beta_i(z) / generationTime) * 0.5 * (power_values(k) + initial_power) * (1 - exp(-decay_constants(z) * time_step)) * (1 / decay_constants(z)));
    end
        
    C_old = concentration_values;
    initial_power = power_values(k);
end

figure;
plot(time_values, power_values);
title('Solution to PRKE');
xlabel('Time (Seconds)');
ylabel('Power (watts)');
