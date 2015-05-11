#include <iostream>
#include <fstream>

#define IJ(i, j, nx)  (i)*(nx) + ( j )
#define PI  3.141592653589793 

using namespace std;

template<typename Ptr> void printmatrix(const char *fname,
					int nx, int ny, Ptr arr)
{

  ofstream myfile;
  myfile.open(fname);
  
  int i,j;
  for (i = 0; i < nx; i++) {
    for (j=0; j < ny; j++) {
      myfile << arr[i * nx + j] << " ";
    }
    myfile << endl;
  }
  myfile.close();
}

template<typename Ptr> void print_state(const char* fname,
					int nx, int ny, Ptr arr){
  ofstream myfile;
  myfile.open(fname);
  
  int i,j;
  for (i = 1; i < nx+1; i++) {
    for (j=1; j < ny+1; j++) {
      myfile << arr[i * ( nx +2 ) + j] << " ";
    }
    myfile << endl;
  }
  myfile.close();
}
