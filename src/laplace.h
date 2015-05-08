#include "umfpack.h"
#ifdef __cplusplus
extern "C" {
#endif

struct LaplacianOp{
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
  
  LaplacianOp(int nx, int ny);
  ~LaplacianOp();
  void set_lambda(double lambda);
  void apply_laplacian(double *y, double *x);
  void backward_solve(double* x,double*  work);
  void laplacian_solve(double * Ax, double*x, double *b);
  
};


  void fill_boundary(const int bc_type, double* u, int nx, int ny);
#ifdef __cplusplus
}
#endif
