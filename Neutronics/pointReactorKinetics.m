% OSTR Hybrid Model, Point reactor kinetics module
% Solves power & precursor concentration over time using RK4 method.
% Known reactivity insertion.

% Parameters
% Rho will eventually become some defined function based on parameters such as xenon, control rod insertion, etc.
beta = 0.0065;
generationTime = 0.0001;
lambda = 0.08; % one-delayed group of neutrons (for now at least)

% initial conditions
time_step = 0.01; % in seconds
initial_time = 0; % in seconds
final_time = 15; % in seconds
initial_power = 1; % Megawatt
initial_concentration = 0;

% initialize vectors to store values from RK4 iterations
time_values = initial_time:time_step:final_time;
power_values = zeros(length(time_values));  % naked / magic numbers. use variables.
power_values(1) = initial_power;   % should power_values be initialized as a 2-d array?
concentration_values = zeros(length(time_values));  % should concentration_values be initialized as a 2-d array?
concentration_values(1) = initial_concentration;

% RK4 iteration
for i = 1:(length(power_values)-1)  % what you had before would not actually index through the array
    rho = (0.05 * i) - 0.2 ;  % naked / magic numbers here. what is being calculated?
    k1 = time_step * changeInPower(time_values(i), power_values(:,i), rho, concentration_values(i));  % is power_values 1-d or 2-d? you use a mix
    k2 = time_step * changeInPower(time_values(i) + time_step/2, power_values(:,i) + k1 / 2, rho, concentration_values(i));
    k3 = time_step * changeInPower(time_values(i) + time_step/2, power_values(:,i) + k2 / 2, rho, concentration_values(i));
    k4 = time_step * changeInPower(time_values(i) + time_step, power_values(:,i) + k3, rho, concentration_values(i));
    power_values(i+1) = power_values(i) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;

    z1 = time_step * changeInConcentration(time_values(i), power_values(:,i), concentration_values(i));
    z2 = time_step * changeInConcentration(time_values(i) + time_step/2, power_values(:,i) + z1 / 2, concentration_values(i));
    z3 = time_step * changeInConcentration(time_values(i) + time_step/2, power_values(:,i + z2 / 2), concentration_values(i));
    z4 = time_step * changeInConcentration(time_values(i) + time_step, power_values(:,i) + z3, concentration_values(i));
    concentration_values(:,i+1) = concentration_values(:,i) + (z1 + 2 * z2 + 2 * z3 + z4) / 6;   % same concern with concentration_values. 1-d or 2-d?

end

plot(time_values, power_values(1,:));
title('PRKE RK4 Method');
xlabel('Time');
ylabel('Power');

% I dont think these functions work how they are intended. I will
% investigate
% function dPdt = changeInPower(time, power, rho, concentration)
%     dPdt = ( (rho - beta) / generationTime) * power(time) + (lambda * concentration(time));
% 
% end
% 
% function dCdt = changeInConcentration(time, power, concentration)
%     dCdt = (beta - generationTime) * power(time) - lambda * concentration(time);
% 
% end
