#ifdef __cplusplus
extern "C" {
#endif

  int test_solve_laplace(int n);
  int test_setup_laplacian();
  int test_apply_laplace();

  void setup_solvers(int nx, int ny);
  void laplacian_solve(double*x, double *b);
  void free_solvers();
  
#ifdef __cplusplus
}
#endif
