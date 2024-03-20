import numpy as np
import matplotlib.pyplot as plt

# Initialize an array of time steps we want to analyze
time_steps = np.array([3, 2, 1, 0.3, 0.2, 0.1, 0.01, 0.001])

# Create a figure to plot the results of different time steps
# We use ax1 because we are analyzing a single axis
# More information about matplotlib and subplots I found here: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html
fig,ax1 = plt.subplots()

# Create main for loop that will iterate through each time step we have selected and perform our calculations

for i in range(len(time_steps)):
    increment = int(15/time_steps[i]) # This initializes the number of elements in our array. A smaller time step will analyze more.
    time = np.linspace(0,15,increment) # This creates an array of numbers from 0 to 15 with increments determined by our time steps

    # Parameters
    generationTime = 0.000001 # Mean neutron generation time as given in the project description
    decay_constants = np.array([0.0128, 0.0301, 0.124, 0.325, 1.12, 2.69]) # An array that holds the decay constant values for each precursor group
    beta_i = np.array([0.000073, 0.000626, 0.000443, 0.000685, 0.000181, 0.000092]) # An array that holds the beta values for each precursor group
    beta = sum(beta_i) # Sum each precursor beta to get beta_eff
    initial_power = 3000 * (10**6) # Convert MW to W and define our initial power as given in the project description
    power_values = np.zeros(len(time)) # Create an empty array of power values that's length is defined by the length of the time array as they must match
    concentration_values = np.zeros(6) # Create an empty array of 6 values that will hold each precursor concentration
    reactivity = np.zeros(increment) # Create an empty array of reactivity values that will be as long as the increment defined by the time step we choose.

    # Initialize each precursor group with the equation described in the numerical solution PRKE PDF on canvas
    for k in range(6):
        concentration_values[k] += (beta_i[k] / (decay_constants[k] * generationTime)) * initial_power

    C_old = concentration_values # Set initial precursor concentration
    time_step = time_steps[i] # Set the time step to the actual time step we are on

    for j in range(len(time)):
        # Depending on the time frame we are in, we need to use the correct reactivity function. Each iteration, we will check what time we are currently at and update the reactivity accordingly
        if time[j] < 5:
            reactivity[j] = 0.08 * time[j] * beta # convert $ to reactivity
        else:
            reactivity[j] = (0.4 - 0.08 * (time[j] - 5)) * beta
        alpha = (reactivity[j] - beta) / generationTime

        concentration_sum = 0 # Create a variable called concentration_sum and use it to represent the sum of each concentration and decay constant respectively.
        for x in range(6):
            concentration_sum += decay_constants[x] * C_old[x]

        # Update our power array with our newly calculated value
        power_values[j] = initial_power * np.exp(alpha * time_step) + ((1 / alpha) * (np.exp(alpha * time_step) * (concentration_sum)))

        # Update our concentration array with newly calculated values
        for z in range(6):
            concentration_values[z] = (C_old[z] * np.exp(-decay_constants[z] * time_step)) + ((beta_i[z] / generationTime)) * (0.5) * (power_values[j] + initial_power) * (1 - np.exp(-decay_constants[z] * time_step)) * (1 / decay_constants[z])
        
        initial_power = power_values[j]
        C_old = concentration_values

    ax1.plot(time, power_values, label = f'Power,dt = {time_steps[k]}')

plt.title('Power over time')
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('Power (Watts)')
plt.legend(loc='lower right')
plt.show()




