from brian2 import *
import matplotlib.pyplot as plt
import sys
inp = sys.stdin.readline

dt = 0.1 * ms
# Parameters
tau = 60 * ms
gamma = 0.641
I_0 = 0.334 * nA
a = 270 * Hz / nA
b = 108 * Hz
d = 0.154 * second
tau_ampa = 2 * ms

print("Display noise with music (Y/N)?: ", end="")
x = inp().strip().lower()
sigma = 0.015 * nA if x == "y" else 0.03 * nA

# Synaptic strengths
J_t = 0.28387  # Cross-connection strength

# External input
dt_input = 1 * ms
I_app = 0.0295 * nA
I_ext_A_array = [I_app if t < 500 else 0 * nA for t in range(4001)]
I_ext_B_array = [I_app if 2000 <= t <= 2500 else 0 * nA for t in range(4001)]
I_ext_A = TimedArray(I_ext_A_array, dt=dt_input)
I_ext_B = TimedArray(I_ext_B_array, dt=dt_input)

eqs = '''
    dS_A/dt = -S_A / tau + (1 - S_A) * gamma * r_A : 1
    dS_B/dt = -S_B / tau + (1 - S_B) * gamma * r_B : 1
    r_A = (a * I_A - b) / (-expm1(-d * (a * I_A - b))) : Hz
    r_B = (a * I_B - b) / (-expm1(-d * (a * I_B - b))) : Hz
    I_A = J_self * S_A + J_cross * S_B + I_0 + I_ext_A(t) + I_noise : amp
    I_B = J_self * S_B + J_cross * S_A + I_0 + I_ext_B(t) + I_noise : amp
    dI_noise/dt = -I_noise / tau_ampa + sigma * tau_ampa**-0.5 * xi : amp
'''
group = NeuronGroup(1, eqs, method='euler', dt=dt)
mon = StateMonitor(group, ['r_A', 'r_B'], record=0)

# J_s values to sweep
J_s_values = linspace(0.35, 0.42, 250)  # 8 points between 0.35 and 0.42
output_rates_A = []
store()

for J_s in J_s_values:
    restore()
    J_self = (J_s + J_t) / 2 * nA
    J_cross = (J_t - J_s) / 2 * nA

    run(4000 * ms)

    # Average firing rate over last 1000 ms
    avg_r_A = np.mean(mon.r_A[0][30000:]) / Hz  # 30000 = 3000 ms / 0.1 ms
    output_rates_A.append(avg_r_A)

    print(f"J_s = {J_s}, Avg r_A = {avg_r_A:.2f} Hz")

print(*output_rates_A)
# Plot average firing rate as a function of J_s
plt.figure(figsize=(10, 4))
plt.plot(J_s_values, output_rates_A, marker='o', linestyle='None')
plt.xlabel('Self-connection strength J_s')
plt.ylabel('Average firing rate of A (Hz)')
plt.title('Output firing rate of A vs J_s')
plt.grid(True)
plt.tight_layout()
plt.show()
