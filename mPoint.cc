// mPoint.cpp: implementation of the mPoint class.
//
//////////////////////////////////////////////////////////////////////

#include "mPoint.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

mPoint mPoint::operator * (double other) const 
{
	return mPoint(pos[0]*other, pos[1]*other, pos[2]*other);
}

mPoint mPoint::operator + (const mPoint &other) const 
{
	return mPoint(pos[0]+other.pos[0], pos[1]+other.pos[1], pos[2]+other.pos[2]);
}

mPoint& mPoint::operator += (const mPoint &other) 
{
	pos[0]+=other.pos[0]; 
	pos[1]+=other.pos[1];
	pos[2]+=other.pos[2];
	return *this;
}

mPoint mPoint::operator - (const mPoint &other) const 
{
	return mPoint(pos[0]-other.pos[0], pos[1]-other.pos[1], pos[2]-other.pos[2]);
}


int mPoint :: operator< (const mPoint &other) const
{
  if(pos[0]<other.pos[0])return 1;
  if(pos[0]>other.pos[0])return 0;
  if(pos[1]<other.pos[1])return 1;
  if(pos[1]>other.pos[1])return 0;
  if(pos[2]<other.pos[2])return 1;
  if(pos[2]>other.pos[2])return 0;  
  return 0;
}

int mPoint :: operator> (const mPoint &other) const
{
	// has to be implemented
	return 0;
}

void mPoint :: print()
{
	printf("%e %e %e \n",pos[0],pos[1],pos[2]);
}