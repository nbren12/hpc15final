#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <math.h>

#include "util.hpp"

#include "laplace.h"

#define PERIODIC_BC    1
#define DIRICHLET_BC   2



/*************************************************************
 *          Setup for sparse laplacian
 *************************************************************/

LaplacianOp::LaplacianOp(int nx, int ny) : nx(nx), ny(ny) {
  int i,j;

  // number of non zero entries
  // -4 for the corners
  n = (nx+2)*(ny+2);
  nz = 5 * nx * ny + 2 * ( nx + 2) + 2 * ( ny + 2) - 4;

  // Allocate arrays
  Ap = new int[(nx+2) * (ny+2) + 1];
  Ai = new int[nz];
  Ax = new double[nz];

  // Allocate forward and backward operators
  Abackward = new double[nz];

  int offset = 0;
  int apoffset = 0;
  int val;
  int nonzero;

  for (i = 0; i < nx+2; i++) {
    for(j = 0; j < ny+2; j++) {
      Ap[apoffset++]= offset;

      if (i == 0 || j==0 || i == nx + 1 || j == ny + 1){
	// boundary points
	Ai[offset] = IJ(i,j,nx + 2);
	Ax[offset++]           = 1.0;
      } else {  
	// interior points
	Ai[offset] = IJ(i-1, j, nx +2);
	Ax[offset++] = 1.0;

	Ai[offset] = IJ(i, j-1, nx +2);
	Ax[offset++] = 1.0;

	Ai[offset] = IJ(i, j, nx +2);
	Ax[offset++] = -4.0;
      
	Ai[offset] = IJ(i, j+1, nx +2);
	Ax[offset++] = 1.0;

	Ai[offset] = IJ(i+1, j, nx +2);
	Ax[offset++] = 1.0;
      }
    }
  }

  Ap[apoffset]= offset;
}

void LaplacianOp::set_lambda(double lambda){



  int offset =0 ;
  int i,j;
  for (i = 0; i < nx+2; i++) {
    for(j = 0; j < ny+2; j++) {

      if (i == 0 || j==0 || i == nx + 1 || j == ny + 1){
	// boundary points
	Abackward[offset++] = 1.0 - lambda/2.0 * 1.0;
      } else {  
	// interior points
	Abackward[offset++] = - lambda/2.0 * 1.0;

	Abackward[offset++] = - lambda/2.0 * 1.0;

	Abackward[offset++] = 1.0 + lambda/2.0 * 4.0;
      
	Abackward[offset++] = - lambda/2.0 * 1.0;

	Abackward[offset++] = - lambda/2.0 * 1.0;
      }
    }
  }

  // Perform LU  decomposition
  int n_row = n;
  int n_col = n_row;
  
  int status;
  status =  umfpack_di_symbolic(n_row, n_col, Ap, Ai,
				Abackward, &Symbolic, Control, Info);
  status = umfpack_di_numeric(Ap, Ai, Abackward,
			      Symbolic, &Numeric, Control, Info);


}

void LaplacianOp::apply_laplacian(double *y, double *x){
  int i,j;

  fill_boundary(PERIODIC_BC, x, nx, ny);
  for (i = 0; i < nx+2; i++) {
    for(j = 0; j < ny+2; j++) {
      if (i == 0 || j==0 || i == nx + 1 || j == ny + 1){
	// boundary points
	y[IJ(i,j, nx+2)] = x[IJ(i,j, nx+2)];
      } else {  
	// interior points
	y[IJ(i,j, nx+2)] = -4.0 * x[IJ(i,j,nx+2)] +
	  x[IJ(i-1,j,nx+2)] + x[IJ(i+1,j,nx+2)] +
	  x[IJ(i,j-1,nx+2)] + x[IJ(i,j+1,nx+2)];
      }
    }
  }
}

void fill_boundary(const int bc_type, double* u, int nx, int ny){
  int i;
  switch (bc_type){ 

  case PERIODIC_BC:

    for (i = 1; i < nx+1; i++) {
      u[IJ(i,0,nx+2)] = u[IJ(i,ny,nx+2)];
      u[IJ(i,ny+1,nx+2)] = u[IJ(i,1,nx+2)];
    }

    for (i = 1; i < ny+1; i++) {
      u[IJ(0,i,nx+2)] = u[IJ(nx, i,nx+2)];
      u[IJ(nx+1,i,nx+2)] = u[IJ(1,i,nx+2)];
    }

    //Don't forget corners
    u[IJ(0,0, nx+2)] = u[IJ(ny, nx, nx+2)];
    u[IJ(ny+1,nx+1, nx+2)] = u[IJ(1, 1, nx+2)];

    u[IJ(0,nx+1, nx+2)] = u[IJ(ny, 1, nx+2)];
    u[IJ(ny+1,0, nx+2)] = u[IJ(1, nx, nx+2)];
    break;
  }
}

void free_solvers(LaplacianOp & lapl){
  umfpack_di_free_numeric(&lapl.Numeric);
  umfpack_di_free_symbolic(&( lapl.Symbolic ));
  free(lapl.Ap);
  free(lapl.Ai);
  free(lapl.Ax);
  free(lapl.Abackward);
}




/*************************************************************
 *          Laplacian Solver
 *************************************************************/

void LaplacianOp::laplacian_solve(double * Ax, double*x, double *b){
  
  int status;
  fill_boundary(PERIODIC_BC, x, nx, ny);
  status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax,
			    x, b, Numeric, Control, Info);
}

void LaplacianOp::backward_solve(double* x,double*  work){
  laplacian_solve(Abackward, x, work);
}

/* @doc: test for building laplacian operator
 *
 * just runs code. doesn't do any tests.
 */
// int test_setup_laplacian()
// {
//   int nx = 10;
//   int ny = 10;
//   setup_laplacian(nx, ny);

//   // printmatrix(1, (nx+2) * (ny+2)+1, lapl.Ap);
//   // printmatrix(1, (nx+2) * (ny+2), lapl.Ai);

//   return 0;
// }


// int test_solve_laplace(int n){
//   int nx = n;
//   int ny = n;

//   // steady state
//   double * b, *x;
//   b  = new double[(nx+2)*(ny+2)];
//   x  = new double[(nx+2)*(ny+2)];
//   int i, j;

//   const double L = 1.0;
//   int k = 2;
//   int l = 4;

//   double dx =  L / nx;
//   double dy =  L / ny;

//   for (i = 1; i < nx +1; i++) {
//     for (j = 1; j < ny +1; j++) {
//       b[IJ(i,j,nx+2)] = sin(2 * PI / L *3 * (i-1) *dx) * sin(2*PI/L*(j-1)*dy);
//     }
//   }

 
//   setup_laplacian(nx, ny);
//   setup_umfpack
//   laplacian_solve(x,b);


//   // Output file
//   print_state("forcing.txt", nx, ny, b);
//   print_state("solution.txt", nx, ny, x);
  
//   free_solvers();
//   free(b);
//   free(x);
//   return 0;
// }


