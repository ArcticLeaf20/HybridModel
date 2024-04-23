%OSTR Hybrid Model PRKE Module

%Time
h = 0.5; %time step
time_values = 0:h:15; % matrix of time values 0 to 15 seconds at 0.5 second increment

%Parameters
generationTime = 1 * (10^-7); % Mean neutron generation time
decay_constants = [0.0128, 0.301, 0.124, 0.325, 1.12, 2.69]; % decay constants for each precursor group
beta_i = [0.000073, 0.000626, 0.000443, 0.000685, 0.000181, 0.000092]; % beta values for each precursor group
beta = sum(beta_i); % beta effective
concentration_values = zeros(1, 6); % empty array to store concentration values
changein_concentration = zeros(1, 6);
power_values = zeros(1, length(time_values)); % empty array to store power values
initial_power = 1 * 10^6; % initial power in watts

for j = 1:6
    concentration_values(j) = (beta_i(j) / (decay_constants(j) * generationTime)) * initial_power;
end

C_old = concentration_values;
reactivity = 0.04;

for i = 1:length(time_values)
    concentration_sum = 0;
    for x = 1:6
        concentration_sum = concentration_sum + decay_constants(x) * C_old(x);
    end
    
    k1 = ((reactivity - beta) / generationTime) * initial_power + concentration_sum;
    k2 = ((reactivity - beta) / generationTime) * (initial_power + h*(k1/2)) + concentration_sum;
    k3 = ((reactivity - beta) / generationTime) * (initial_power + h*(k2/2)) + concentration_sum;
    k4 = ((reactivity - beta) / generationTime) * (initial_power + h*k3) + concentration_sum;
    new_power = initial_power + (h / 6) * (k1 + 2*k2 + 2*k3 + k4);
    initial_power = new_power;

    for z = 1:6
        changein_concentration(z) = (beta / generationTime) * initial_power - decay_constants(z) * C_old(z);
        concentration_values(z) = changein_concentration(z) + concentration_values(z);
    end
    C_old = concentration_values;
    power_values(i) = new_power;
end