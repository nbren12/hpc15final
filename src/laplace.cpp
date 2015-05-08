#include "laplace.h"
#define IJ(i,j,n) (i)*(n)+(j)


using namespace arma;

LaplaceOp::LaplaceOp(int n) : L(sp_mat(n*n, n*n)), I(sp_mat(n*n, n*n)),
			      Afor(sp_mat(n*n, n*n)), Aback(sp_mat(n*n, n*n)) {
  double lambda = 1.0;
  I = speye(n*n, n*n);
  int i,j;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if ( i > 0 )
	L(IJ(i,j,n), IJ(i-1,j,n)) = 1.0;

      if ( j > 0 )
	L(IJ(i,j,n), IJ(i,j-1,n)) = 1.0;

      L(IJ(i,j,n), IJ(i,j,n)) =  -4.0;

      if ( j < n-1 )
	L(IJ(i,j,n), IJ(i,j+1,n)) = 1.0;

      if ( i < n-1 )
	L(IJ(i,j,n), IJ(i+1,j,n)) = 1.0;
    }
  }

}

void LaplaceOp::set_lambda(double lambda){
  
  Aback = I - lambda/2.0 * L;
  Afor  = I + lambda/2.0 * L;
}

void test_laplace_matrix(){
  LaplaceOp lapl(10);
  lapl.set_lambda(1.0);
  lapl.Afor.print();
}


int main(int argc, char *argv[])
{
  test_laplace_matrix();
  return 0;
}


