% OSTR Hybrid Model, Point reactor kinetics module

function pointReactorKinetics
    % Define parameters
    beta = 0.007;   % Delayed neutron fraction
    lambda = [0.08, 0.12, 0.25, 0.3, 1.2, 3.0];  % Decay constants

    % Set up initial conditions
    P0 = 1.0;  % Initial power

    % Set up time span
    tspan = [0, 15];

    % Solve the system using ode15s
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    [t, results] = ode15s(@reactorODE, tspan, [P0, 0], options, beta, lambda);

    % Extract power and reactivity results
    P = results(:, 1);
    reactivity = results(:, 2);

    % Plot the results
    figure;

    % Plot Power
    subplot(2, 1, 1);
    plot(t, P, 'LineWidth', 2);
    title('Point Reactor Kinetics');
    xlabel('Time');
    ylabel('Power');
    grid on;

    % Plot Reactivity
    subplot(2, 1, 2);
    plot(t, reactivity, 'LineWidth', 2);
    xlabel('Time');
    ylabel('Reactivity');
    grid on;
end

function dXdt = reactorODE(t, X, beta, lambda)
    % Point reactor kinetics equations
    P = X(1);
    
    % Define piecewise reactivity function
    if t >= 0 && t < 5
        reactivity = 0.08 * t;
    elseif t >= 5 && t <= 15
        reactivity = 0.4 - 0.08 * (t-5);
    else
        reactivity = 0; % Reactivity is zero outside the specified intervals
    end

    lambda_eff = sum(lambda) - beta;
    dPdt = (reactivity * P + beta * sum(lambda .* P)) / lambda_eff;
    dRhodt = reactivity;
    dXdt = [dPdt; dRhodt];
end
