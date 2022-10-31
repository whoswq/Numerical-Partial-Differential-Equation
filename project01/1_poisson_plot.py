import numpy as np
import matplotlib.pyplot as plt


def U(r):
    return -0.25 * r * r + 0.25


U_h = np.fromfile("1_poisson_Nr_16.dat", dtype="float64")
N = len(U_h) + 1
delta_r = 1 / (N - 0.5)
r_array = np.array([delta_r * (i + 0.5) for i in range(N - 1)])

U_exact = U(r_array)

plt.plot(r_array, U_h, label="Numerical results, N = 16")
plt.plot(r_array, U_exact, label="Exact results")
plt.legend()
plt.xlabel("r")
plt.ylabel("U(r)")
plt.savefig("1_poisson_Nr_16.pdf")
plt.show()
print(U_h - U_exact)
