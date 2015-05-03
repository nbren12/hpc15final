#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "umfpack.h"

#define PERIODIC_BC    1
#define DIRICHLET_BC   2
#define IJ(i, j, nx)  (i)*(nx) + j

typedef struct {
  int nx;
  int ny;
  double dx;
  int bc_type;
  // More grid parameters
} Grid;


typedef struct {
  double t; // Current time
  // Grid
  double * data;
  Grid * grid;
} State;


typedef struct {
  int * Ai;
  int * Ap;
  double * Ax;
} LaplacianOp;


Grid grid;
LaplacianOp lapl;
State state;

void setup_laplacian(int nx, int ny){
  int i,j;

  // number of non zero entries
  // -4 for the corners
  int nz = 5 * nx * ny + 2 * ( nx + 2) + 2 * ( ny + 2) - 4;

  int *Ap, *Ai;
  double *Ax;

  // Allocate arrays
  Ap = calloc((nx+2) * (ny+2) + 1, sizeof(int));
  Ai = calloc(nz, sizeof(int));
  Ax = calloc(nz, sizeof(double));

  int offset = 0;
  int val;
  int nonzero;

  for (i = 0; i < nx+2; i++) {
    for(j = 0; j < ny+2; j++) {
      Ap[IJ(i, j, nx + 2)]= offset;

      if (i == 0 || j==0 || i == nx + 1 || j == ny + 1){
	// boundary points
	Ai[offset] = IJ(i,j,nx);
	Ax[offset++] = 1.0;
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

  lapl.Ap = Ap;
  lapl.Ai = Ai;
  lapl.Ax = Ax;

}

/* @doc: test for building laplacian operator
 *
 * just runs code. doesn't do any tests.
 */
int test_build_laplacian_operator()
{
  int nx = 100;
  int ny = 100;
  setup_laplacian(nx, ny);

  return 0;
}


int test_solve_laplace(){
  int nx = 1000;
  int ny = 1000;

  // steady state
  double * b, *x;
  b = calloc((nx+2)*(ny+2), sizeof(double));  // right hand side
  x = calloc((nx+2)*(ny+2), sizeof(double));  // solution
  int i, j;

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      b[IJ(i,j,nx)] = 0.0;
    }

  }

  


  setup_laplacian(nx,ny);


  // UMFPack stuff
  int status;
  void *Symbolic, *Numeric;
  double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;

  int n_row = (nx+2)*(nx+2);
  int n_col = n_row;
  
  status =  umfpack_di_symbolic(n_row, n_col, lapl.Ap, lapl.Ai,
				lapl.Ax, &Symbolic, Control, Info);

  status = umfpack_di_numeric(lapl.Ap, lapl.Ai, lapl.Ax,
			      Symbolic, &Numeric, Control, Info);

  status = umfpack_di_solve(UMFPACK_A, lapl.Ap, lapl.Ai, lapl.Ax,
			    x, b, Numeric, Control, Info);
  return 0 ;
}

int main(int argc, char *argv[])
{
 
  test_build_laplacian_operator();
  test_solve_laplace();
  return 0;
}
