#include <cmath>
#include <iostream>
#include <cstdio>

#include "util.hpp"
#include "laplace.h"
#define OUTPUT_FORMAT "output/%07d.txt"

using namespace std;



int test_evolve_heat_equation_2d(){
  // Setup grid
  const double L = 1.0;
  const int nx = 100;
  const int ny = nx;

  const double dx = 1.0/nx;
  const double dy = dx;

  // Setup time stepping
  double dt =dx*dx/2 * 10;
  int nt = 100;
  int output_interval = 5;

  // Allocate arrays
  double *x0, *work;
  x0  = new double[(nx+2)*(ny+2)];
  int i, j;

  // Initialize x0
  int k = 2;
  int l = 4;

  for (i = 1; i < nx +1; i++) {
    for (j = 1; j < ny +1; j++) {
      x0[IJ(i,j,nx+2)] = sin(2 * PI / L *3 * (i-1) *dx) * sin(2*PI/L*(j-1)*dy);
    }
  }

 
  // Solve heat equation
  evolve_heat_equation_2d(x0, nx, dx, nt, dt, output_interval);

  free(x0);
  return 0;
}


int main(int argc, char *argv[])
{
  test_evolve_heat_equation_2d();
  return 0;
}
