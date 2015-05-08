#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <math.h>

#include "util.hpp"

#include "laplace.h"

#define PERIODIC_BC    1
#define DIRICHLET_BC   2


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
LaplacianOp::~LaplacianOp() {

  umfpack_di_free_symbolic(&Symbolic);
  umfpack_di_free_numeric(&Numeric);


  free(Ai);
  free(Ax);
  free(Ap);
  free(Abackward);

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

void LaplacianOp::laplacian_solve(double * Ax, double*x, double *b){
  
  int status;
  fill_boundary(PERIODIC_BC, x, nx, ny);
  status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax,
			    x, b, Numeric, Control, Info);
}

void LaplacianOp::backward_solve(double* x,double*  work){
  laplacian_solve(Abackward, x, work);
}



