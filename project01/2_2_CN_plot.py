import numpy as np
import matplotlib.pyplot as plt

# x_array_1 = np.array([1 / 32 * i for i in range(33)])
# u_1_32 = np.fromfile("2_1_diffusion_1_32.dat", dtype="float64")
x_array = np.array([1/50 * i for i in range(1, 51)])
u_T1 = np.fromfile("2_2_diffusion_2_T1.dat", dtype="float64")
u_T2 = np.fromfile("2_2_diffusion_2_T2.dat", dtype="float64")
u_T3 = np.fromfile("2_2_diffusion_2_T3.dat", dtype="float64")
u_T4 = np.fromfile("2_2_diffusion_2_T4.dat", dtype="float64")

plt.plot(x_array, u_T1, label="T = 1")
plt.plot(x_array, u_T2, label="T = 2")
plt.plot(x_array, u_T3, label="T = 3")
plt.plot(x_array, u_T4, label="T = 4")
plt.xlabel("x")
plt.ylabel("u")
plt.legend()
plt.savefig("2_2_CN_2.pdf")
plt.show()