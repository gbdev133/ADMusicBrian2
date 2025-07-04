from brian2 import *
import matplotlib.pyplot as plt
import sys
inp = sys.stdin.readline


dt = 0.1 * ms
# Parameters
tau = 60*ms
gamma = 0.641
I_0 = 0.334 * nA
a = 270 * Hz/nA
b = 108 * Hz
d = 0.154 * second
tau_ampa = 2 * ms
output_rates_A = []
print("Display noise with music (Y/N)?: ", end="")
x = inp().strip().lower()
sigma = 0.015 * nA if x == "y" else 0.03 * nA

# Synaptic strengths
J_s = 0.38  # Self-connection strength (<0.4182 in Murray)
J_t = 0.28387  # Cross-connection strength
# J_self = 0.316935 * nA
J_self = (J_s + J_t)/2 * nA
# J_cross = -0.033065 * nA  # derived by solving a sys of equations with Js = 0.35 nA, Jt = 0.28387 nA
J_cross = (J_t - J_s)/2 * nA

# Time-varying external input
dt_input = 1*ms
I_app = 0.0295 * nA
# I_ext_A_array = [I_app if 0 <= x <= 500 else 0 * nA for x in range(4001)]
I_ext_A_array = [I_app if x < 500 else 0 * nA for x in range(4001)]
I_ext_B_array = [I_app if 2000 <= x <= 2500 else 0 * nA for x in range(4001)]
I_ext_A = TimedArray(I_ext_A_array, dt=dt_input)
I_ext_B = TimedArray(I_ext_B_array, dt=dt_input)
# print(I_ext_A_array)
# print(I_ext_B_array)

# Define equations for two populations (A and B)
eqs = '''
dS_A/dt = -S_A/tau + (1 - S_A) * gamma * r_A : 1
dS_B/dt = -S_B/tau + (1 - S_B) * gamma * r_B : 1
r_A = (a*I_A - b) / (1 - exp(-d*(a*I_A - b))) : Hz
r_B = (a*I_B - b) / (1 - exp(-d*(a*I_B - b))) : Hz
I_A = J_self*S_A + J_cross*S_B + I_0 + I_ext_A(t) + I_noise : amp
# I_A = J_self*S_A + J_cross*S_B + I_0 + int(t < 500 * ms) * 0.075 * nA * sin(2 * pi * 50 * t/second) + I_noise : amp
I_B = J_self*S_B + J_cross*S_A + I_0 + I_ext_B(t) + I_noise : amp
dI_noise/dt = -I_noise / tau_ampa + sigma * tau_ampa**-0.5 * xi: amp
# I_ext_A: amp
# I_noise : amp
'''

# Create a dummy NeuronGroup containing both units
group = NeuronGroup(1, eqs, method='euler', dt=dt)

# Monitor
mon = StateMonitor(group, ['S_A', 'S_B', 'r_A', 'r_B'], record=0)

run(4000*ms)

# Plotting firing rates over time
plt.figure(figsize=(10, 4))
plt.plot(mon.t / ms, mon.r_A[0] / Hz, label='r_A (Hz)', color='tab:blue')
plt.plot(mon.t / ms, mon.r_B[0] / Hz, label='r_B (Hz)', color='tab:orange')
plt.xlabel('Time (ms)')
plt.ylabel('Firing rate (Hz)')
plt.title('Firing Rates Over Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
