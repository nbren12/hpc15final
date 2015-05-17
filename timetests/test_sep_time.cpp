#include <cmath>
#include <iostream>
#include <cstdio>

#include "util.hpp"
#include "laplace.h"
#define OUTPUT_FORMAT "output/%07d.txt"


using namespace std;

void save_IC(double *x_block, double *x, int nx, int ny, int k){
  int i = 0;
  int j = 0;
  for (i = 0; i < nx + 2; i++) {
    for (j = 0; j < ny + 2; j++) {
      x_block[(nx+2)*(ny+2)*k+IJ(i,j,nx+2)] = x[IJ(i,j,nx+2)];
    }
  }
}

void get_IC(double *x_block, double *x, int nx, int ny, int k){
  int i = 0;
  int j = 0;
  for (i = 0; i < nx + 2; i++) {
    for (j = 0; j < ny + 2; j++) {
      x[IJ(i,j,nx+2)] = x_block[(nx+2)*(ny+2)*k+IJ(i,j,nx+2)];
    }
  }
}

void update(double *x_old, double *x_new, int nx, int ny){
  int i = 0;
  int j = 0;
  for (i = 0; i < nx + 2; i++) {
    for (j = 0; j < ny + 2; j++) {
      x_old[IJ(i,j,nx+2)] = x_new[IJ(i,j,nx+2)];
    }
  }
}

void update_x_block(double *x_block, double *x_block_update, int nx, int ny, int p){
  int i = 0;
  int j = 0;
  int k = 0;
  for (k = 0; k < p ; k++){
    for (i = 0; i < nx + 2; i++) {
      for (j = 0; j < ny + 2; j++) {
         x_block[(nx+2)*(ny+2)*k+IJ(i,j,nx+2)] =  x_block_update[(nx+2)*(ny+2)*k+IJ(i,j,nx+2)];
      }
    }
  }
}

double norm_abs(double *x_end, double *x_exact, int nx, int ny){
  double sum = 0;
  double R = 0;
  int i = 0;
  int j = 0;
  for (i = 1; i < nx +1; i++) {
    for (j = 1; j < ny +1; j++) {
      R = x_end[IJ(i,j,nx+2)] - x_exact[IJ(i,j,nx+2)];
       //cout << "NORMdiff: " <<  R  << endl;
      //R = x_end[IJ(i,j,nx+2)];
       sum = sum + R*R;
    }
    //cout << "NORMdiffsum: " <<  sum << endl;
  }
  R = sqrt(sum)/((nx+2)*(ny+2));
  //cout << "R: " <<  R << endl;
  return R;
}
double norm_rel(double *x_end, double *x_exact, int nx, int ny){
  double sum = 0;
  double sumD = 0;
  double R = 0;
  double RD = 0;
  int i = 0;
  int j = 0;
  for (i = 1; i < nx +1; i++) {
    for (j = 1; j < ny +1; j++) {
      R = x_end[IJ(i,j,nx+2)] - x_exact[IJ(i,j,nx+2)];
       //cout << "NORMdiff: " <<  R  << endl;
      RD = x_exact[IJ(i,j,nx+2)];
      sum = sum + R*R;
      sumD = sumD + RD*RD;
    }
    //cout << "NORMdiffsum: " <<  sum << endl;
  }
  R = sqrt(sum/sumD)/((nx+2)*(ny+2));
  //cout << "R: " <<  R << endl;
  return R;
}

int test_time_parallel(){
  // Setup grid
  const double L = 1.0;
  const int nx = 20;
  const int ny = nx;

  const double dx = 1.0/nx;
  const double dy = dx;

  // Allocate arrays
  double *x0, *x_exact, *x_fine, *x_end;
  x0  = new double[(nx+2)*(ny+2)];
  x_fine  = new double[(nx+2)*(ny+2)];
  x_exact = new double[(nx+2)*(ny+2)];
  x_end = new double[(nx+2)*(ny+2)];

  int i, j, p;

  // Initialize x0
  int k = 2;
  int l = 4;

  for (i = 0; i < nx + 2; i++) {
    for (j = 0; j < ny + 2; j++) {
      x0[IJ(i,j,nx+2)] = sin(2 * PI / L *3 * (i-1) *dx) * sin(2*PI/L*(j-1)*dy);
      x_exact[IJ(i,j,nx+2)] = x0[IJ(i,j,nx+2)];
    }
  }

  // Setup time stepping
  double T = 10; 		//final integration time
  int nt_fine_block = 20;
  int nt_corse = 4096;
  int nt = nt_corse*nt_fine_block;
  double dt = T/nt;
  double dt_corse = T/nt_corse;
  double dt_fine_block = dt;

  int output_interval = 100;

  int num_iter = nt_corse;

  double TOL = 0.1;  
  double *tol_abs, *tol_rel;
  tol_abs  = new double[num_iter+2];
  tol_rel  = new double[num_iter+2];
  double R;

  // Allocate arrays for each processor
  double *x_block, *x_block_update;
  x_block = new double[(nx+2)*(ny+2)*(nt_corse)];
  x_block_update = new double[(nx+2)*(ny+2)*(nt_corse)];


  //EXACT SOLUTION (full run on the fine grid)
  evolve_heat_equation_2d(x_exact, nx, dx, nt, dt, nt+1);

  cout << "FINISHED EXACT SOLUTION RUN!!" << endl;

 
  //TIME PARALLEL RUN:

  //SOLVE HEAT EQUATION ON CORSE GRID
  for (i = 0; i < nt_corse ; i++) {
    save_IC(x_block,x0,nx,ny,i); //save x0 after each time step to use as IC for time parallel thing
    evolve_heat_equation_2d(x0, nx, dx, 1, dt_corse, nt+1);
    //save_IC(x_block,x0,nx,ny,i); //save x0 after each time step to use as IC for time parallel thing
  }
    cout << "FINISHED CORSE RUN!!" << endl;
    //R = norm_abs(x0,x_exact,nx,ny);
    //tol_abs[0] = R;
    R = norm_rel(x0,x_exact,nx,ny);
    tol_rel[0] = R;

    //initialize x_block update same way as x_block
    update_x_block(x_block_update, x_block, nx, ny, nt_corse);

    cout << "STARTING FINE RUN!!" << endl;

    
    // SOLVE HEAT EQUATION ON FINE GRID
    for (p = 0; p < num_iter ; p++) {
      cout << "FINE ITERATION NUMBER " << p+1 << endl;
      for (i = 0; i < nt_corse ; i++) {
	get_IC(x_block,x_fine,nx,ny,i);  //get initial condition for the fine run
	evolve_heat_equation_2d(x_fine, nx, dx, nt_fine_block, dt_fine_block, nt+1);
	if (i+1<nt_corse){
	  save_IC(x_block_update,x_fine,nx,ny,i+1);
	}
	else if(i + 1 == nt_corse){
	  update(x_end,x_fine,nx,ny);
	}     
      }
      update_x_block(x_block, x_block_update, nx, ny, nt_corse);
      // COMPUTE NORM
      //R = norm_abs(x_end,x_exact,nx,ny);
      //tol_abs[p+1] = R;
      R = norm_rel(x_end,x_exact,nx,ny);
      tol_rel[p+1] = R;
      if (tol_rel[p+1] < TOL){
	num_iter = p+1;
      }
    }
    //R = norm_abs(x_exact,x_exact,nx,ny);
    //tol_abs[num_iter+1] = R;
    R = norm_rel(x_exact,x_exact,nx,ny);
    tol_rel[num_iter+1] = R;
    
    /* cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << "ABSOLUTE TOLERANCE: " << endl;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    for (p = 0; p < num_iter+1 ; p++) {
      cout << tol_abs[p] << endl;
    }
    cout <<"COMPARE EXACT TO EXACT TO CHECK" << endl;
    cout << tol_abs[num_iter+1] << endl;*/

        
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << "RELATIVE TOLERANCE: " << endl;
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    for (p = 0; p < num_iter+1 ; p++) {
      cout << tol_rel[p] << endl;
    }
    cout <<"COMPARE EXACT TO EXACT TO CHECK" << endl;
    cout << tol_rel[num_iter+1] << endl;
    cout << "NUMBER OF ITERATIONS NECESSARY WAS: " << num_iter << endl;
    cout << "NUMBER OF CORSE BLOCKS WAS: " << nt_corse << endl;

  
  free(x0);
  free(x_exact);
  free(x_fine);
  free(x_block);
  free(x_block_update);
  free(x_end);
  free(tol_rel);
  free(tol_abs);

  return 0;
}


int main(int argc, char *argv[])
{
  test_time_parallel();
  return 0;
}



