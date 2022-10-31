import numpy as np
import matplotlib.pyplot as plt

u_1_4 = np.fromfile("3_1_wave_1_4.dat", dtype="float64")
u_1_8 = np.fromfile("3_1_wave_1_8.dat", dtype="float64")
u_1_16 = np.fromfile("3_1_wave_1_16.dat", dtype="float64")
u_1_32 = np.fromfile("3_1_wave_1_32.dat", dtype="float64")
u_1_64 = np.fromfile("3_1_wave_1_64.dat", dtype="float64")
u_1_128 = np.fromfile("3_1_wave_1_128.dat", dtype="float64")
u_1_256 = np.fromfile("3_1_wave_1_256.dat", dtype="float64")

error_2 = []
error_inf = []
x_array = [np.log(4), np.log(8), np.log(16), np.log(32), np.log(64), np.log(128)]
error_2.append(np.log(
    np.linalg.norm((u_1_8[::2] - u_1_4), ord=2) / np.sqrt(4)))
error_2.append(
    np.log(np.linalg.norm((u_1_16[::2] - u_1_8), ord=2) / np.sqrt(8)))
error_2.append(
    np.log(np.linalg.norm((u_1_32[::2] - u_1_16), ord=2) / np.sqrt(16)))
error_2.append(
    np.log(np.linalg.norm((u_1_64[::2] - u_1_32), ord=2) / np.sqrt(32)))
error_2.append(
    np.log(np.linalg.norm((u_1_128[::2] - u_1_64), ord=2) / np.sqrt(64)))
error_2.append(
    np.log(np.linalg.norm((u_1_256[::2] - u_1_128), ord=2) / np.sqrt(128)))
error_inf.append(np.log(np.linalg.norm(u_1_8[::2] - u_1_4, ord=np.inf)))
error_inf.append(np.log(np.linalg.norm(u_1_16[::2] - u_1_8, ord=np.inf)))
error_inf.append(np.log(np.linalg.norm(u_1_32[::2] - u_1_16, ord=np.inf)))
error_inf.append(np.log(np.linalg.norm(u_1_64[::2] - u_1_32, ord=np.inf)))
error_inf.append(np.log(np.linalg.norm(u_1_128[::2] - u_1_64, ord=np.inf)))
error_inf.append(np.log(np.linalg.norm(u_1_256[::2] - u_1_128, ord=np.inf)))

plt.plot(x_array, error_2, label="$L^2$ norm")
plt.plot(x_array, error_inf, label="$L^{\infty}$ norm")
plt.legend()
plt.xlabel("$\ln h^{-1}$")
plt.ylabel("$\ln||U_h - U_{h/2}||$")
plt.savefig("3_1_error_new.pdf")
plt.show()

print((error_2[-1] - error_2[-2]) / (x_array[-1] - x_array[-2]))
print((error_inf[-1] - error_inf[-2]) / (x_array[-1] - x_array[-2]))