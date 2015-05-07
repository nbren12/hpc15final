#include <cmath>

#include "util.hpp"
#include "laplace.h"



void evolve_heat_equation_2d(double *x, int n, const double tout,
			     const double dx, const double dt_min,
			     double* work){

  int i;
  double lambda = dt_min / dx / dx;


  setup_laplacian(n,n);
  set_lambda_cn(lambda);


  // forward step
  apply_laplacian(work, x);
  for (i = 0; i < (n+2)*(n+2); i++) {
    x[i] += work[i]/2.0;
    work[i] = x[i]; // copy back into work array
  }

  // backward step
  backward_solve(x, work);
  
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
  
  print_state("t0.txt", nx,ny, x0);
  evolve_heat_equation_2d(x0, nx, dt_outer, dx, dt_inner, x);
  print_state("t1.txt", nx,ny, x0);

  free_solvers();

  free(x);
  free(x0);
  return 0;
}
