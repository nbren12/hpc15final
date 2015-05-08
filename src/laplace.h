#include<armadillo>

struct LaplaceOp {
  arma::sp_mat L, I, Afor, Aback;

  LaplaceOp(int n);
  void set_lambda(double);
};
