#include <math.h>
#include <stdlib.h>

#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "My_Matrix.h"

using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

int main() {
  double dx = 1.0 / 50.0;
  double mu = 2.5;  // grid ratio
  int N = 1.0 / dx;
  double dt = mu * dx * dx;
  double time = 0;
  double time_max = 4;
  double x_array[N + 1] = {0};
  for (int i = 0; i < N + 1; i++) {
    x_array[i] = (double)i / N;
  }

  // apply for memory
  My_Matrix U(N, 1);
  double save_u[N] = {0};
  My_Matrix F(N, 1);

  My_Matrix A_h(N, N);
  // construct matrix
  A_h.set_element(0, 0, mu + 1.0);
  A_h.set_element(0, 1, -0.5 * mu);
  for (int row = 1; row < N - 1; row++) {
    A_h.set_element(row, row - 1, -0.5 * mu);
    A_h.set_element(row, row, mu + 1.0);
    A_h.set_element(row, row + 1, -0.5 * mu);
  }
  A_h.set_element(N - 1, N - 2, -mu);
  A_h.set_element(N - 1, N - 1, mu + 1);

  // construct initial condition
  //   u(x, 0) = x^2 / 2
  //   for (int i = 0; i < N; i++) {
  //     U.set_element(i, 0, x_array[i + 1] * x_array[i + 1] / 2.0);
  //   }
  //   u(x, 0) = x^2 - x
  for (int i = 0; i < N; i++) {
    U.set_element(i, 0, x_array[i + 1] * (x_array[i + 1] - 1));
  }

  int cnt = 0;
  while (time < time_max) {
    // construct F
    F.set_element(
        0, 0,
        0.5 * mu * U.read_element(1, 0) + (1.0 - mu) * U.read_element(0, 0));
    for (int j = 1; j < N - 1; j++) {
      F.set_element(
          j, 0,
          0.5 * mu * (U.read_element(j - 1, 0) + U.read_element(j + 1, 0)) +
              (1.0 - mu) * U.read_element(j, 0));
    }
    F.set_element(N - 1, 0,
                  2.0 * mu * dx + mu * U.read_element(N - 2, 0) +
                      (1.0 - mu) * U.read_element(N - 1, 0));

    // update all inintial points
    U = A_h.LU_solvers(F);
    time += dt;
    cnt += 1;
    if (cnt == 999) {
      U.show_matrix();
      for (int i = 0; i < N; i++) {
        save_u[i] = (double)U.read_element(i, 0);
      }
      ofstream out;
      out.open("2_2_diffusion_2_T1.dat", std::ofstream::binary);
      out.write(reinterpret_cast<const char*>(save_u), sizeof(double) * (N));
      out.close();
    } else if (cnt - 1000 == 999) {
      U.show_matrix();
      for (int i = 0; i < N; i++) {
        save_u[i] = (double)U.read_element(i, 0);
      }
      ofstream out;
      out.open("2_2_diffusion_2_T2.dat", std::ofstream::binary);
      out.write(reinterpret_cast<const char*>(save_u), sizeof(double) * (N));
      out.close();
    } else if (cnt - 2000 == 999) {
      U.show_matrix();
      for (int i = 0; i < N; i++) {
        save_u[i] = (double)U.read_element(i, 0);
      }
      ofstream out;
      out.open("2_2_diffusion_2_T3.dat", std::ofstream::binary);
      out.write(reinterpret_cast<const char*>(save_u), sizeof(double) * (N));
      out.close();
    } else if (cnt - 3000 == 999) {
      U.show_matrix();
      for (int i = 0; i < N; i++) {
        save_u[i] = (double)U.read_element(i, 0);
      }
      ofstream out;
      out.open("2_2_diffusion_2_T4.dat", std::ofstream::binary);
      out.write(reinterpret_cast<const char*>(save_u), sizeof(double) * (N));
      out.close();
    }
  }
  cout << "total steps is " << cnt << endl;

  return 0;
}