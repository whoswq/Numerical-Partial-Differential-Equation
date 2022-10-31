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
// using namespace std;

using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

int main() {
  int N = 16;
  long double delta_r = 1.0 / (N - 0.5);
  My_Matrix A_h(N - 1, N - 1);
  My_Matrix F_h(N - 1, 1);
  My_Matrix U_h(N - 1, 1);
  for (int row = 0; row < N - 1; row++) {
    // construct the tridiagonal matrix
    if (row == 0) {
      A_h.set_element(0, 0, 2);
      A_h.set_element(0, 1, -2);
    } else if (row > 0 and row < N - 2) {
      A_h.set_element(row, row - 1, -(double)row / (0.5 + row));
      A_h.set_element(row, row, 2.0);
      A_h.set_element(row, row + 1, -(1.0 + row) / (0.5 + row));
    } else {
      A_h.set_element(row, row - 1, -(double)row / (0.5 + row));
      A_h.set_element(row, row, 2.0);
    }
    // construct F
    F_h.set_element(row, 0, delta_r * delta_r);
  }
  // this gives the results on 1, 2, ..., N-1
  U_h = A_h.LU_solvers(F_h);
  U_h.show_matrix();

  // save data
  // I want to use python to plot, but numpy only support
  // double precision floats
  double save_u[N - 1] = {0};
  for (int i = 0; i < N - 1; i++) {
    save_u[i] = (double)U_h.read_element(i, 0);
  }
  ofstream out;
  out.open("1_poisson_Nr_16.dat", std::ofstream::binary);
  out.write(reinterpret_cast<const char *>(save_u), sizeof(double) * (N - 1));
  out.close();

  return 0;
}