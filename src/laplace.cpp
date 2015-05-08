#include<iostream>
#include<armadillo>

#define IJ(i,j,n) (i)*(n)+(j)

#define TFORMAT "output/times.txt"
#define FORMAT "output/%07d.txt"

struct LaplaceOp {
  arma::sp_mat L, I, Afor, Aback;

  LaplaceOp(int n);
  void set_lambda(double);
};


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

void periodic_boundary(vec u, int n){
  int i,j;
  int m = n -2;
  
  for (i = 1; i < n-1; i++)
    u(IJ(0,i,n)) = u(IJ(m,i,n));

  for (i = 1; i < n-1; i++)
    u(IJ(n-1,i,n)) = u(IJ(1,i,n));

  for (i = 1; i < n-1; i++)
    u(IJ(i,0,n)) = u(IJ(i,m,n));

  for (i = 1; i < n-1; i++)
    u(IJ(i,n-1,n)) = u(IJ(i,1,n));

  u(IJ(0,0,n)) = 0.0;
  u(IJ(m,0,n)) = 0.0;
  u(IJ(0,m,n)) = 0.0;
  u(IJ(m,m,n)) = 0.0;
}


template<typename A> void savefile(int count, A & u){
  // Setup output
  char filename[100];
  sprintf(filename, FORMAT, count);
  cout << "Outputting " << filename << endl;
  u.save(filename, arma_ascii);
}

void evolve_heat_equation_2d(double *x, int n, double dx,
			     int nt, double dt, int output_interval){

  int i = 0;

  // Wrap the array with a armadillo object
  vec u(x, (n+2) * (n+2), false);

  // Working array
  vec work(x, (n+2) * (n+2));

  // Setup timestepping
  double lambda = dt / dx /dx;
  LaplaceOp lapl(n+2);
  lapl.set_lambda(lambda);

  // Output times
  int count= 0;
  vec times(nt+1);
  for (i = 0; i < nt+1; i++) {
    times(i) = i * dt;
  }
  times.save(TFORMAT, arma_ascii);

  // Do time stepping
  savefile(count++, u);   // Output initial condition
  for (i = 1; i < nt + 1; i++) {
    cout << "Time = " << i * dt<< endl;
    // Forward step
    periodic_boundary(u,n+2);
    work = lapl.Afor * u;

    // Backward step
    periodic_boundary(work,n+2);
    u = spsolve(lapl.Aback,work);
    if (i % output_interval == 0) savefile(count++, u);

  }
}
