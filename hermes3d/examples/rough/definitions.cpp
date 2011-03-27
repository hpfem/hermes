#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <limits>
#include <cmath>
 
class roughSurface
{
public:
	std::string fileName, name;
	double xOff, yOff ; // x,y offset of the patch
	double dx, dy ; // step size along x and y direction
	int Dx, Dy; // number of divisions along x and y
	double minZ, maxZ, maxDeltaZ;
	double** Z ;
	void setFilneName(std::string fName) ;
	int loadSurface(void) ;
	double interpolate(double x, double y) ;
	void findMinMax() ;
private:
	void holdThreshold(double &max, double &min, double in) ; // expands the range of a min and a max variable according to new input
	void compare(double &a, double &b) ; // swap two numbers such that the former is greater
	void sort(double &a, double &b, double &c, double &d) ; // sort a, b, c, d in descending order
};


void roughSurface::compare(double &a, double &b)
{
  double u ;
  if(a < b) {u=a; a=b; b=u;}
}

void roughSurface::holdThreshold(double &max, double &min, double in)
{
  if(max < in) {max = in ;} ;
  if(min > in) {min = in ;} ;    
}

void roughSurface::sort(double &a, double &b, double &c, double &d)
{
  compare(a, b) ;
  compare(a, c) ;
  compare(a, d) ;
  compare(b, c) ;
  compare(b, d) ;
  compare(c, d) ;
}

void roughSurface::setFilneName(std::string fName)
{
   this->fileName = fName ;
}

int roughSurface::loadSurface(void)
{
  std::ifstream myReadFile;
  myReadFile.open(this->fileName.c_str());
  if (!myReadFile) 
  {
    error("\nCannot read input surface file %s.", this->fileName.c_str());
  }
  else
  {
     myReadFile >> name ;
     myReadFile >> xOff >> yOff ;
     myReadFile >> dx >> dy ;
     myReadFile >> Dx >> Dy ;
     std::cout << "\nSurface data statisitcs dx, dy, Dx, Dy: " << dx << " " << dy << " " << Dx << " " << Dy << std::endl ;
     Z = new double*[Dy+1] ;
     for(int j=0 ; j < Dy+1 ; j++) Z[j] = new double[Dx+1] ;
     for(int j=0 ; j < Dy+1 ; j++) for(int i=0 ; i < Dx+1 ; i++) myReadFile >> Z[j][i] ;
     myReadFile.close() ;
     std::cout << "\n" << "The (x,y) origin offset for the surface is: " << xOff << ", " << yOff ;
     std::cout << "\n" << "The (x,y) dimensions of the surface are: " << dx*Dx << ", " << dy*Dy ;
     std::cout << "\n" << "The original (x,y) steps size are: " << dx << ", " << dy ;
     return 0 ;
  }
  return 0;
}

double roughSurface::interpolate(double x, double y)
{
  x = x - xOff ;
  y = y - yOff ;
  int i00 = (int) floor(x/dx) ;
  int j00 = (int) floor(y/dy) ;
  int shiftX = 1;
  int shiftY = 1;
  if( i00 > Dx || j00 > Dy || x < 0. || y < 0. )
  {
    std::cout << "(x,y) out of range! " << x << " " << y << std::endl ;
    return std::numeric_limits<double>::quiet_NaN() ;
  }
  else
  {
    if(i00 == Dx) shiftX = -1;
    if(j00 == Dy) shiftY = -1;
  }
  double x00 = i00*dx ;
  double y00 = j00*dy ;
  double xi1 = (x-x00)/dx ;
  double xi2 = (y-y00)/dy ;
  return xi1*xi2*Z[j00+shiftY][i00+shiftX] + xi1*(1-xi2)*Z[j00][i00+shiftX] + (1-xi1)*xi2*Z[j00+shiftY][i00] + (1-xi1)*(1-xi2)*Z[j00][i00] ;
}

void roughSurface::findMinMax()
{
	minZ=+std::numeric_limits<double>::max() ;
	maxZ=-std::numeric_limits<double>::max() ; 
	double a, b, c, d ;
	for(int j=0 ; j<Dy  ; j++) for(int i=0 ; i<Dx  ; i++)// sweep across all cells to find min and max of z
	{
		a=Z[j][i] ;
		b=Z[j+1][i] ;
		c=Z[j][i+1] ;
		d=Z[j+1][i+1] ;
		sort(a, b, c, d) ;
		compare(maxZ, a) ;
		compare(d, minZ) ;
	}
	std::cout << "\nMaximum height: " << maxZ ;
	std::cout << "\nMinimum height: " << minZ ;
}



