#ifndef _MGEOM_SEARCH_H
#define _MGEOM_SEARCH_H
#include <vector>
#define MIN(x,y) ((x<y)?(x):(y))
#define MAX(x,y) ((x>y)?(x):(y))
#define TOL 1.e-06
#include "mCompiler.h"

template <class T>
class Brick
{
public:
  vector<T*> Objects;
  Brick(){}
  T* operator [] (int i)
    {
      if(i < 0 || i >= Objects.size())throw i;
      T *obj = Objects[i];
      return obj;
    }
  int size(){return Objects.size();}
};

template <class T>
class mGeomSearch
{
  Brick<T> *theTable;
  int Nx,Ny,Nz;
  double Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
  int getBrickId (double X, double Y, double Z)
    {
      int Ix = (int)((double)Nx * (X-Xmin) / (Xmax-Xmin)); 
      int Iy = (int)((double)Ny * (Y-Ymin) / (Ymax-Ymin)); 
      int Iz = (int)((double)Nz * (Z-Zmin) / (Zmax-Zmin)); 
      Ix = MIN(Ix,Nx-1);
      Iy = MIN(Iy,Ny-1);
      Iz = MIN(Iz,Nz-1);
      if(Ix<0)Ix=0;
      if(Iy<0)Iy=0;
      if(Iz<0)Iz=0;
      int index = Ix + Iy * Nx + Iz * Nx * Ny;
      return index;
    }
  Brick<T> * getBrick (int index)
    {
      if(index <0 || index >= Nx*Ny*Nz)throw index;
      Brick<T> *b = &(theTable[index]);
      return b;
    }
public:
  mGeomSearch (double x1,
	      double x2,
	      double y1, 
	      double y2, 
	      double z1, 
	      double z2,
	      int nx = 10, 
	      int ny = 10, 
	      int nz = 10) :	Nx(nx), Ny(ny), Nz(nz)
    {
      Xmin = x1-TOL; Xmax  = x2+TOL;
      Ymin = y1-TOL; Ymax  = y2+TOL;
      Zmin = z1-TOL; Zmax  = z2+TOL;
      
      theTable = new Brick<T> [Nx*Ny*Nz];
    }
  
  ~mGeomSearch ()
    {
      delete [] theTable;
    };
  Brick<T> *getBrick(double X, double Y, double Z)
    {
      return (getBrick(getBrickId(X,Y,Z)));
    }
  bool AddObject ( T *);
};

template <class T> 
bool mGeomSearch<T> :: AddObject ( T * obj )
{
  double  XminBox,XmaxBox,YminBox,YmaxBox,ZmaxBox,ZminBox;
  int     Ix1,Ix2,Iy1,Iy2,Iz1,Iz2;
  int     i,j,k,index;
  Brick<T>   *pBrick;
  
  /*Template Objects must overload getBox function*/
  obj->getBox (XminBox,XmaxBox,YminBox,YmaxBox,ZmaxBox,ZminBox);
  
  Ix1 = (int)( (double)Nx * (XminBox - Xmin) /( Xmax - Xmin )); 
  Ix2 = (int)( (double)Nx * (XmaxBox - Xmin) /( Xmax - Xmin )); 
  Iy1 = (int)( (double)Ny * (YminBox - Ymin) /( Ymax - Ymin ));  
  Iy2 = (int)( (double)Ny * (YmaxBox - Ymin) /( Ymax - Ymin ));  
  Iz1 = (int)( (double)Nz * (ZminBox - Zmin) /( Zmax - Zmin ));  
  Iz2 = (int)( (double)Nz * (ZmaxBox - Zmin) /( Zmax - Zmin )); 
  
  
  Ix1 = MAX(Ix1,0);
  Ix2 = MIN(Ix2,Nx-1);
  Iy1 = MAX(Iy1,0);
  Iy2 = MIN(Iy2,Ny-1);
  Iz1 = MAX(Iz1,0);
  Iz2 = MIN(Iz2,Nz-1);
 
  for(i=Ix1;i<=Ix2;i++){
    for(j=Iy1;j<=Iy2;j++){
      for(k=Iz1;k<=Iz2;k++){
	index = i + j * Nx + k * Nx * Ny;
	pBrick = getBrick(index);
	pBrick->Objects.push_back(obj);
      }
    }
  }
  return true;
}
#undef TOL
#endif
