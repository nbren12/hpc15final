#ifdef __cplusplus
extern "C" {
#endif

  int test_solve_laplace(int n);
  int test_setup_laplacian();
  int test_apply_laplace();

  void setup_laplacian(int nx, int ny);
  void set_lambda_cn(double lambda);
  void laplacian_solve(double * Ax, double*x, double *b);
  void apply_laplacian(double *y, double *x);
  void backward_solve(double* x,double*  work);
  void free_solvers();

  void fill_boundary(const int bc_type, double* u, int nx, int ny);
#ifdef __cplusplus
}
#endif
