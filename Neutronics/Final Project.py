import numpy as np
import matplotlib.pyplot as plt

time_steps = np.array([3, 2, 1, 0.3, 0.2, 0.1, 0.01, 0.001])

fig,ax1 = plt.subplots()

for i in range(len(time_steps)):
    increment = int(15/time_steps[i])
    time_values = np.linspace(0,15,increment)

    reactivity = np.zeros(increment)
    generationTime = 10 * (10**-6)

    concentration_values = np.zeros(6)
    beta = 0.0021
    beta_i = np.array([0.000073, 0.000626, 0.000443, 0.000685, 0.000181, 0.000092])
    initial_power = 3000*(10**6)
    power_values = np.zeros(len(time_values))
    decay_constants = np.array([0.0128, 0.0301, 0.124, 0.325, 1.12, 2.69])

    for j in range(6):
        concentration_values[j] += (beta_i[j]/(decay_constants[j] * generationTime)) * initial_power
    
    C_old = concentration_values

    t = time_steps[i]

    for k in range(len(time_values)):

        if time_values[k] < 5:
            reactivity[k] = 0.08 * time_values[k] * beta
        else:
            reactivity[k] = (0.4 - 0.08 * (time_values[k] - 5)) * beta

        alpha = (reactivity[k] - beta) / generationTime

        C_SUM = 0
        for x in range(6):
            C_SUM += decay_constants[x] * C_old[x]
        
        power_values[k] = initial_power * np.exp(alpha * t) + ((1 / alpha) * (np.exp(alpha*t) - 1) * (C_SUM))

        for z in range(6):
            concentration_values[z] = (C_old[z] * np.exp(-decay_constants[z] * t)) + ((beta_i[z]/generationTime)) * (0.5) * (power_values[k] + initial_power) * (1 - np.exp(-decay_constants[z] * t)) * (1 / decay_constants[z])

        initial_power = power_values[k]
        C_old = concentration_values

    ax1.plot(time_values,power_values,label = f'Power,dt = {time_steps[i]}')

plt.title('Power over time')
ax1.set_xlabel('Time (seconds)')
ax1.set_ylabel('Power (watts)')
plt.legend(loc='upper right')
plt.grid(True)
plt.show()




