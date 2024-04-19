% OSTR hybrid simulation, reactor kinetics module
% Developed by Aidan Marsters & Gabriel Eilhardt
% Solves for power and precursor concentrations over a given time interval
% Uses RK4 explicit numerical method to solve kinetics equations

% time configuration
h = 0.5; % set your h here to use as the time step for the time values
time_values = 0:h:15;  % use h to create the time values, so you only need to make one change

% Values
generationTime = 1 * (10^-7); % Mean neutron generation time
decay_constants = [0.0128, 0.0301, 0.124, 0.325, 1.12, 2.69]; % Decay constants for each precursor group
beta_i = [0.000073, 0.000626, 0.000443, 0.000685, 0.000181, 0.000092]; % Beta values for precursor groups
beta = sum(beta_i); % Effective beta
concentration_values = zeros(1, 6); % Array for concentration values
power_values = zeros(1, length(time_values)); % Array for power values  you may want to add room for the initial value -AMW
initial_power = 1 * 10^6; % Initial power in watts
reactivity = zeros(1, length(time_values)); % Array for reactivity values    this should be the same length as the time values, right? -AMW
% what are the units for everything above? reactivity? power? generation
% time? etc

for j = 1:6
   concentration_values(j) = (beta_i(j) / (decay_constants(j) * generationTime)) * initial_power;
end

C_old = concentration_values;
concentration_sum = 0;

for i = 1:length(time_values)
    for x = 1:6
        concentration_sum = concentration_sum + decay_constants(x) * C_old(x);  % this is blowing up, check that this is the correct math to be doing. -AMW
    end
    k1 = ((reactivity(i) - beta) / generationTime) * initial_power + concentration_sum;
    k2 = ((reactivity(i) - beta) / generationTime) * (initial_power + h*(k1/2)) + concentration_sum;  % I don't think the reactivity should be changing here because you do not know it for the next time step, right? -AMW
    k3 = ((reactivity(i) - beta) / generationTime) * (initial_power + h*(k2/2)) + concentration_sum;
    k4 = ((reactivity(i) - beta) / generationTime) * (initial_power + h*k3) + concentration_sum;
    new_power = initial_power + (h / 6) * (k1 + 2*k2 + 2*k3 + k4);
    initial_power = new_power;

    for z = 1:6
        concentration_values(z) = (beta / generationTime) * initial_power - decay_constants(z) * C_old(z);  % does this need to be a part of the RK4 scheme? or is using the RK4 solved power enough? -AMW
    end

    C_old = concentration_values;
    power_values(i) = initial_power; % similarily to above, this does not store the initial power value in the array. not a big deal, just be aware of it. -AMW
end

