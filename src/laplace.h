#include "umfpack.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int * Ai;
  int * Ap;
  double * Ax;
  double * Abackward;
  double lambda;
  int nx;
  int ny;
  int nz;
  int n;

  int status;
  void *Symbolic, *Numeric;
  double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
} LaplacianOp;

  int test_solve_laplace(int n);
  int test_setup_laplacian();
  int test_apply_laplace();

  void setup_laplacian(int nx, int ny, LaplacianOp & lapl);
  void set_lambda_cn(double lambda, LaplacianOp & lapl);
  void apply_laplacian(double *y, double *x, LaplacianOp & lapl);
  void backward_solve(double* x,double*  work, LaplacianOp & lapl);
  void free_solvers(LaplacianOp & lapl);

  void fill_boundary(const int bc_type, double* u, int nx, int ny);
#ifdef __cplusplus
}
#endif
