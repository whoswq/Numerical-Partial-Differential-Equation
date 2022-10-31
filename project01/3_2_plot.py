import numpy as np
import matplotlib.pyplot as plt

T = 1


def U(x):
    return np.exp(-5 * (x - T - 1) * (x - T - 1))


U_nu_0_5 = np.fromfile("3_2_wave_nu_0_5.dat", dtype="float64")
U_nu_1_1 = np.fromfile("3_2_wave_nu_1_1.dat", dtype="float64")
U_nu_2_1 = np.fromfile("3_2_wave_nu_2_1.dat", dtype="float64")

x_array = np.array([1/32 * i - 1 for i in range(6 * 32 + 1)])
U_exact = U(x_array)

plt.plot(x_array, U_nu_0_5, label="$\\nu = 0.5$")
plt.plot(x_array, U_nu_1_1, label="$\\nu = 1$")
plt.plot(x_array, U_nu_2_1, label="$\\nu = 2$")
plt.plot(x_array, U_exact, label="exact")

plt.legend()
plt.xlabel("x")
plt.ylabel("u")
plt.savefig("3_2_upwind_compare.pdf")
plt.show()