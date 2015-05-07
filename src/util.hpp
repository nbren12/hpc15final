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
