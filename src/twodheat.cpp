#include <cmath>

#include "util.hpp"
#include "laplace.h"


void evolve_heat_equation_2d(double *x, double *y, int nx,
			     const double dt_outer, const double dt_inner){
  
}


int test_evolve_heat_equation_2d(){
  // Setup grid
  const double L = 1.0;
  const int nx = 100;
  const int ny = nx;

  const double dx = 1.0/nx;
  const double dy = dx;

  // Setup time stepping
  const double dt_outer = 1.0;
  const double dt_inner = .1;

  // Allocate arrays
  double *x0, *x;
  x0  = new double[(nx+2)*(ny+2)];
  x  = new double[(nx+2)*(ny+2)];
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
  setup_solvers(nx, ny);
  evolve_heat_equation_2d(x0, x, nx, dt_outer, dt_inner);
  print_state("t0.txt", nx,ny, x0);
  print_state("t1.txt", nx,ny, x);

  free_solvers();

  free(x);
  free(x0);
}
