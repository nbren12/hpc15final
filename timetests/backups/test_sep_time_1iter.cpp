#include <cmath>
#include <iostream>
#include <cstdio>

#include "util.hpp"
#include "laplace.h"
#define OUTPUT_FORMAT "output/%07d.txt"
#define IJK(i,k,ni,nj) (i)*(ni)*(nj)+(k)


using namespace std;
void save_IC(double *x_block, double *x0, int nx, int ny, int k){
  int i = 0;
  int j = 0;
  for (i = 1; i < nx +1; i++) {
    for (j = 1; j < ny +1; j++) {
      x_block[(nx+2)*(ny+2)*k+IJ(i,j,nx+2)] = x0[IJ(i,j,nx+2)];
    }
  }

}

void get_IC(double *x_block, double *x0, int nx, int ny, int k){
  int i = 0;
  int j = 0;
  for (i = 1; i < nx +1; i++) {
    for (j = 1; j < ny +1; j++) {
      x0[IJ(i,j,nx+2)] = x_block[(nx+2)*(ny+2)*k+IJ(i,j,nx+2)];
    }
  }

}

double norm(double *x_fine, double *x_exact, int nx, int ny){
  double sum = 0;
  double R = 0;
  int i = 0;
  int j = 0;
  for (i = 1; i < nx +1; i++) {
    for (j = 1; j < ny +1; j++) {
       R = x_fine[IJ(i,j,nx+2)] - x_exact[IJ(i,j,nx+2)];
       //cout << "NORMdiff: " <<  R  << endl;
       sum = sum + R*R;
    }
    //cout << "NORMdiffsum: " <<  sum << endl;
  }
  R = sqrt(sum)/((nx+2)*(ny+2));
  cout << "R: " <<  R << endl;
  return R;
}

int test_evolve_heat_equation_2d(){
  // Setup grid
  const double L = 1.0;
  const int nx = 100;
  const int ny = nx;

  const double dx = 1.0/nx;
  const double dy = dx;

  // Allocate arrays
  double *x0, *x_exact, *x_fine;
  x0  = new double[(nx+2)*(ny+2)];
  x_fine  = new double[(nx+2)*(ny+2)];
  x_exact = new double[(nx+2)*(ny+2)];
  int i, j;

  // Initialize x0
  int k = 2;
  int l = 4;

  for (i = 1; i < nx +1; i++) {
    for (j = 1; j < ny +1; j++) {
      x0[IJ(i,j,nx+2)] = sin(2 * PI / L *3 * (i-1) *dx) * sin(2*PI/L*(j-1)*dy);
      x_exact[IJ(i,j,nx+2)] = x0[IJ(i,j,nx+2)];
    }
  }

  // Setup time stepping
  double T = 2; 		//final integration time
  int nt_fine_block = 1;
  int nt_corse = 3;
  int nt = nt_corse*nt_fine_block;
  double dt = T/nt;
  double dt_corse = T/nt_corse;
  double dt_fine_block = T/nt_fine_block;

  int output_interval = 100;

  // Allocate arrays for each processor
  double *x_block;
  x_block = new double[(nx+2)*(ny+2)*(nt_corse)];


  //EXACT SOLUTION (full run on the fine grid)
  evolve_heat_equation_2d(x_exact, nx, dx, nt, dt, nt+1);

  cout << "Finished fine run!!" << endl;

 
  //ONE ITERATION OF TIME PARALLEL RUN:

  // Solve heat equation on corse grid
  for (i = 0; i < nt_corse ; i++) {
    evolve_heat_equation_2d(x0, nx, dx, nt_corse, dt_corse, nt+1);
    save_IC(x_block,x0,nx,ny,i); //save x0 after each time step to use as IC for time parallel thing
  }
    cout << "Finished corse run!!" << endl;


  // Solve heat equation on fine grid
  for (i = 0; i < nt_corse ; i++) {
    get_IC(x_block,x_fine,nx,ny,i);  //get initial condition for the fine run
    cout << "Past get and Finishedstarting fine run!!" << endl;

    evolve_heat_equation_2d(x_fine, nx, dx, nt_fine_block, dt_fine_block, nt+1);
    save_IC(x_block,x_fine,nx,ny,i);
  }

  // COMPUTE NORM
  double R;
  R = norm(x_fine,x_exact,nx,ny);
        cout << "R: " <<  R << endl;

    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

  cout << "The 2 norm between the exact and time parallel result is: " << R << endl;

    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

  free(x0);
  free(x_exact);
  free(x_fine);
  free(x_block);
  return 0;
}


int main(int argc, char *argv[])
{
  test_evolve_heat_equation_2d();
  return 0;
}



