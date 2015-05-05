#include <cstdlib>
#include <assert.h>
#include <math.h>

#include "umfpack.h"
#include "util.hpp"

#define PERIODIC_BC    1
#define DIRICHLET_BC   2
#define IJ(i, j, nx)  (i)*(nx) + j
#define PI  3.141592653589793 

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
  Ap = new int[(nx+2) * (ny+2) + 1];
  Ai = new int[nz];
  Ax = new double[nz];

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

void fill_boundary(const int bc_type, double* u, int nx, int ny){
  int i;
  if (bc_type == PERIODIC_BC){

    for (i = 0; i < nx; i++) {
      u[IJ(i,0,nx+2)] = u[IJ(i,ny,nx+2)];
      u[IJ(i,ny+1,nx+2)] = u[IJ(i,1,nx+2)];
    }

    for (i = 0; i < ny; i++) {
      u[IJ(0,i,nx+2)] = u[IJ(nx, i,nx+2)];
      u[IJ(nx+1,i,nx+2)] = u[IJ(1,i,nx+2)];
    }
  }
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


int test_solve_laplace(int n){
  int nx = n;
  int ny = n;

  // steady state
  double * b, *x;
  b  = new double[(nx+2)*(ny+2)];
  x  = new double[(nx+2)*(ny+2)];
  int i, j;

  const double L = 1.0;
  int k = 2;
  int l = 4;

  double dx =  L / nx;
  double dy =  L / ny;

  for (i = 1; i < nx +1; i++) {
    for (j = 1; j < ny +1; j++) {
      b[IJ(i,j,nx+2)] = sin(2 * PI / L * (i-1) *dx) * sin(2*PI/L*(j-1)*dy);
    }
  }

  fill_boundary(PERIODIC_BC, x, nx, ny);
  


  setup_laplacian(nx,ny);


  // UMFPack stuff
  int status;
  void *Symbolic, *Numeric;
  double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL] ;

  int n_row = (nx+2)*(nx+2);
  int n_col = n_row;
  
  status =  umfpack_di_symbolic(n_row, n_col, lapl.Ap, lapl.Ai,
				lapl.Ax, &Symbolic, Control, Info);

  assert(status==0);
  status = umfpack_di_numeric(lapl.Ap, lapl.Ai, lapl.Ax,
			      Symbolic, &Numeric, Control, Info);

  status = umfpack_di_solve(UMFPACK_A, lapl.Ap, lapl.Ai, lapl.Ax,
			    x, b, Numeric, Control, Info);


  // TODO Output answer to file
  printmatrix(nx+2, 1, lapl.Ap);
  
  return 0;
}


