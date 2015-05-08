#include <cmath>
#include <iostream>
#include <cstdio>

#include "util.hpp"
#include "laplace.h"
#define OUTPUT_FORMAT "output/%07d.txt"

using namespace std;


void evolve_heat_equation_2d(double *x, int n, double dx,
			     int nt, double dt, int output_interval){

  int i;

  // Setup time stepping
  const double lambda = dt / dx / dx;

  cout << "Running heat solver" << endl;
  cout << "Total time = " << nt*dt << endl;
  cout << "dt = " << dt << endl;
  cout << "num timesteps = " << nt << endl;

  // Allocate stuff
  LaplacianOp lapl(n, n);
  lapl.set_lambda(lambda);
  double * work = new double[(n+2) * (n+2)];

  int it;

  char output_filename[100];
  int count=0;

  sprintf(output_filename, OUTPUT_FORMAT, count++);
  cout << 0* dt << " " << output_filename << endl;
  print_state(output_filename, n, n, x);
  
  for (it = 1; it < nt+1; it++) {

    // forward step
    lapl.apply_laplacian(work, x);
    for (i = 0; i < (n+2)*(n+2); i++) {
      x[i] += work[i]*lambda/2.0;
      work[i] = x[i]; // copy back into work array
    }

    // backward step
    lapl.backward_solve(x, work);

    if (output_interval > 0) {
      if (it%output_interval == 0){
	sprintf(output_filename, OUTPUT_FORMAT, count++);
	cout << it * dt << " " << output_filename << endl;
	print_state(output_filename, n, n, x);
      
      }
    }
  }

  free(work);
}


int test_evolve_heat_equation_2d(){
  // Setup grid
  const double L = 1.0;
  const int nx = 100;
  const int ny = nx;

  const double dx = 1.0/nx;
  const double dy = dx;

  // Setup time stepping
  double dt =dx*dx/2;

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
  evolve_heat_equation_2d(x0, nx, 1.0/nx, 100, dt*10, 5);

  free(x0);
  return 0;
}
