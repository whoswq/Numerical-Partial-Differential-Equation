import numpy as np
import matplotlib.pyplot as plt


def U(x, T):
    return np.exp(-5 * (x - T - 1) * (x - T - 1))


U_0_5 = np.fromfile("3_3_LW_T_0_5.dat", dtype="float64")
U_1_0 = np.fromfile("3_3_LW_T_1_0.dat", dtype="float64")
U_1_5 = np.fromfile("3_3_LW_T_1_5.dat", dtype="float64")
U_2_0 = np.fromfile("3_3_LW_T_2_0.dat", dtype="float64")
U_4_5 = np.fromfile("3_3_LW_T_4_5.dat", dtype="float64")
U_4_0 = np.fromfile("3_3_LW_T_4_0.dat", dtype="float64")

x_array = np.array([1 / 100 * i - 1 for i in range(6 * 100 + 1)])

plt.plot(x_array, U_0_5, label="Lax-Wendroff")
plt.plot(x_array, U(x_array, 0.5), label="Exact")
plt.legend()
plt.xlabel("x")
plt.ylabel("u")
plt.title("T = 0.5")
plt.savefig("3_3_LW_0_5.pdf")
plt.show()
plt.close()

plt.plot(x_array, U_1_0, label="Lax-Wendroff")
plt.plot(x_array, U(x_array, 1), label="Exact")
plt.legend()
plt.xlabel("x")
plt.ylabel("u")
plt.title("T = 1")
plt.savefig("3_3_LW_1_0.pdf")
plt.show()
plt.close()

plt.plot(x_array, U_1_5, label="Lax-Wendroff")
plt.plot(x_array, U(x_array, 1.5), label="Exact")
plt.legend()
plt.xlabel("x")
plt.ylabel("u")
plt.title("T = 1.5")
plt.savefig("3_3_LW_1_5.pdf")
plt.show()
plt.close()

plt.plot(x_array, U_2_0, label="Lax-Wendroff")
plt.plot(x_array, U(x_array, 2), label="Exact")
plt.legend()
plt.xlabel("x")
plt.ylabel("u")
plt.title("T = 2")
plt.savefig("3_3_LW_2_0.pdf")
plt.show()
plt.close()

plt.plot(x_array, U_4_0, label="Lax-Wendroff")
plt.plot(x_array, U(x_array, 4), label="Exact")
plt.legend()
plt.xlabel("x")
plt.ylabel("u")
plt.title("T = 4.0")
plt.savefig("3_3_LW_4_0.pdf")
plt.show()
plt.close()