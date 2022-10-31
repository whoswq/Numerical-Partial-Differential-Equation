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

using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

int main() {
  double dx = 1.0 / 100.0;
  double nu = 0.5;     // grid ratio
  int N = 6.0 / dx;  // grid number
  double dt = nu * dx;
  double time = 0;
  double time_max = 4;
  double x_array[N + 1] = {0};
  for (int i = 0; i < N + 1; i++) {
    x_array[i] = (double)i * dx - 1.0;
  }

  // apply for memory
  double* U_old = NULL;
  double* U_new = NULL;
  double* U_swp = NULL;
  U_new = new double[N + 1];
  U_old = new double[N + 1];

  // construct initial condition
  for (int i = 0; i < N + 1; i++) {
    U_old[i] = exp(-5.0 * (x_array[i] - 1.0) * (x_array[i] - 1.0));
  }
  int cnt = 0;
  while (time < time_max) {
    // update all inintial points
    U_new[0] = 0.0;
    for (int j = 1; j < N; j++) {
      U_new[j] = -0.5 * nu * (1.0 - nu) * U_old[j + 1] + (1 - nu * nu) * U_old[j] + 0.5 * nu * (1.0 + nu) * U_old[j - 1];
    }
    U_new[N] = (-0.5 * nu + 1.0 - 0.5 * nu * nu) * U_old[N] + 0.5 * nu * (1.0 + nu) * U_old[N - 1];
    // swap the pointer
    U_swp = U_new;
    U_new = U_old;
    U_old = U_swp;
    U_swp = NULL;
    // update time
    time += dt;
    cnt += 1;
  }
  cout << "total steps is " << cnt << endl;
  ofstream out;
  out.open("3_3_LW_T_4_0.dat", std::ofstream::binary);
  out.write(reinterpret_cast<const char*>(U_old), sizeof(double) * (N + 1));
  out.close();

  //  release the memory
  delete[] U_new;
  delete[] U_old;
  return 0;
}