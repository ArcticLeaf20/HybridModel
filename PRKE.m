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
power_values = zeros(1, length(time_values)); % empty array to store power values
initial_power = 1 * 10^6; % initial power in watts

concentration_value2=zeros(1,6)
% where does the initial concentration value come from? it is very large to
% being with, which makes sense if you are counting atoms/molecules.
% However, this will be added to the power value so it seems like you
% should look into the units.
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
    % should there be a factor for converting from number of neutrons to
    % power? the PRKEs are for neutron population and precursor
    % concentration. if the orders of magnitude are for neutrons and
    % precursor concentration, it will be way too high for power. that
    % seems to be what is happening. the concentration sum is starting out
    % at 1E10 and the (reactivity-beta)/generationTime term is 1E5 so you
    % are multiplying the power by something on the order of 1E5, then
    % adding something on the order of 1E10, so it is going to be a huge
    % number right off the bat
    k1 = ((reactivity - beta) / generationTime) * initial_power + concentration_sum;
    k2 = ((reactivity - beta) / generationTime) * (initial_power + h*(k1/2)) + concentration_sum;
    k3 = ((reactivity - beta) / generationTime) * (initial_power + h*(k2/2)) + concentration_sum;
    k4 = ((reactivity - beta) / generationTime) * (initial_power + h*k3) + concentration_sum;
    new_power = initial_power + (h / 6) * (k1 + 2*k2 + 2*k3 + k4);
    initial_power = new_power;

    for z = 1:6
        concentration_value2(z) = (beta / generationTime) * initial_power - decay_constants(z) * C_old(z);
        concentration_values(z) = concentration_values(z) + concentration_value2(z)*h
        
    end
    C_old = concentration_values;
    power_values(i) = new_power;
end
