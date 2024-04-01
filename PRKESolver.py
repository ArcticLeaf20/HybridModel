import numpy as np
import matplotlib.pyplot as plt
import math

time_steps = np.array([3, 2, 1, 0.3, 0.2, 0.1, 0.01, 0.001]) # initialize an array of time steps we want to analyze

fig,ax1 = plt.subplots() # Create subplots for each time step using matplotlib. Most of the plotting stuff I grabbed from here: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html


# Create a for loop that encapsulates our entire code to make the code iterate through each time step we want it to analyze.
for i in range(len(time_steps)):
    increment = int(15/time_steps[i]) # Define our increment depending on the time step of focus
    time_values = np.linspace(0,15,increment) # create an array of evenly spaced time values that increases by our increment

    # Parameters
    generationTime = 10 * (10**-6) # Mean neutron generation time as defined in the project description
    decay_constants = np.array([0.0128, 0.0301, 0.124, 0.325, 1.12, 2.69]) # An array of decay constants that correspond to each precursor group
    beta_i = np.array([0.000073, 0.000626, 0.000443, 0.000685, 0.000181, 0.000092]) # an array of beta values that correspond to the 6 precursor groups
    beta = sum(beta_i) # Sum the beta of each precursor group to get our beta_eff
    concentration_values = np.zeros(6) # Create an empty arrray for concentration values, the length is 6 as we are dealing with 6 precursor groups
    power_values = np.zeros(len(time_values)) # an empty array of power values that depends on the time step we choose, smaller time step equals more calculations!
    initial_power = 3000*(10**6) # initial power converted to watts
    reactivity = np.zeros(increment) # create an empty array for reactivity values whose length is determined by the increment of focus

    # Initialize our initial precursor concentrations according to the equation described in the PRKE numerical solution PDF on canvas
    for j in range(6):
        concentration_values[j] += (beta_i[j]/(decay_constants[j] * generationTime)) * initial_power
    
    C_old = concentration_values # Set the old concentration to initial concentration

    time_step = time_steps[i] # set our time step to the one we are actually on

    # This for loop is the main part of our code, this calculates new concentration/power values and updates them.
    for k in range(len(time_values)):

        # We first need to figure out where in time we are, depending on the time we need to use the correct reactivity equation
        if time_values[k] < 5:
            reactivity[k] = 0.08 * time_values[k] * beta
        else:
            reactivity[k] = (0.4 - 0.08 * (time_values[k] - 5)) * beta

        alpha = (reactivity[k] - beta) / generationTime # update alpha depending on the reactivity

        concentration_sum = 0 # create a variable to hold the sum of each decay constant and precursor concentration respectively
        for x in range(6):
            concentration_sum += decay_constants[x] * C_old[x]
        
        # Update our power using equation (7) as described in the PRKE numerical solution PDF on canvas
        power_values[k] = initial_power * math.exp(alpha * time_step) + ((1 / alpha) * (math.exp(alpha*time_step) - 1) * (concentration_sum))

        # We need to create a for loop to update each precursor concentration separately.
        for z in range(6):
            concentration_values[z] = (C_old[z] * math.exp(-decay_constants[z] * time_step)) + ((beta_i[z]/generationTime)) * (0.5) * (power_values[k] + initial_power) * (1 - math.exp(-decay_constants[z] * time_step)) * (1 / decay_constants[z])

        # Set initial power/concentration to our newly calculated values as described in the iteration scheme.
        C_old = concentration_values
        initial_power = power_values[k]

    ax1.plot(time_values,power_values,label = f'Power,dt = {time_steps[i]}') # Plots each respective graph of each time step, also creates a legend with the label of power and the defined time step.

# Plot the results and create a legend
plt.title('Power over time')
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('Power (watts)')
plt.legend(loc='upper right')
plt.show()
