// mPoint.h: interface for the mPoint class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _MPOINT_H_
#define _MPOINT_H_
#include <stdio.h>
class mPoint  
{
  double pos[3];
 public:
  mPoint(){}
  mPoint(double x, double y, double z = 0.0)
{	
	pos[0] = x;	pos[1] = y; pos[2] = z;
}
  mPoint& operator +=(const mPoint &);
  mPoint operator +(const mPoint &) const;
  mPoint operator -(const mPoint &) const;
  mPoint operator *(double) const;
  double & operator () (int i) { return pos[i];}
  double operator () (int i) const { return pos[i];}
  int operator < (const mPoint &other) const;
  int operator > (const mPoint &other) const;
  void print();
};

#endif
