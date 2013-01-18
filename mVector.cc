// mVector.cpp: implementation of the mVector class.
//
//////////////////////////////////////////////////////////////////////
#include <math.h>
#include "mVector.h"
#include "mTensor2.h"
#include "mPoint.h"
#include "mException.h"
#include "Constants.h"
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

mVector::mVector(double x, double y, double z)
{
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
}


mVector::mVector(const mPoint &from, const mPoint &to)
{
	pos[0] = to(0) - from(0);
	pos[1] = to(1) - from(1);
	pos[2] = to(2) - from(2);
}

mVector mVector::operator - (const mVector &other)
{
	return mVector(pos[0]-other.pos[0],pos[1]-other.pos[1],pos[2]-other.pos[2]);
}

// cross product
mVector mVector::operator % (const mVector &other)
{
	return mVector(	pos[1]*other.pos[2]-pos[2]*other.pos[1],
					pos[2]*other.pos[0]-pos[0]*other.pos[2],
					pos[0]*other.pos[1]-pos[1]*other.pos[0]);
}

mVector mVector::operator + (const mVector &other)
{
	return mVector(pos[0]+other.pos[0],pos[1]+other.pos[1],pos[2]+other.pos[2]);
}

mVector mVector::operator * (const double &other)
{
	return mVector(pos[0]*other,pos[1]*other,pos[2]*other);
}

mVector mVector::operator / (const double &other)
{
	return mVector(pos[0]/other,pos[1]/other,pos[2]/other);
}

mVector & mVector::operator += (const mVector &other)
{
	pos[0]+=other.pos[0];pos[1]+=other.pos[1];pos[2]+=other.pos[2];
	return *this;
}

mVector & mVector::operator /= (const double &other)
{
  double d = 1./other;
  pos[0]*=d;
  pos[1]*=d;
  pos[2]*=d;
  return *this;
}

mVector & mVector::operator -= (const mVector &other)
{
	pos[0]-=other.pos[0];pos[1]-=other.pos[1];pos[2]-=other.pos[2];
	return *this;
}

void mVector::normalize()
{
	double n =sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	if(n!=0.0) n = 1./n;
	pos[0]*=n;pos[1]*=n;pos[2]*=n;
}

double mVector::L2norm() const
{
	return sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
}

double mVector::angleRad(const mVector &v) const
{
	// have to copy not to modify original vectors
	mVector y(v);
	mVector x(*this);
	x.normalize();
	y.normalize();
	double cosA = x * y;
	mVector cross = x % y;
	double sinA = ::sqrt(cross.pos[0]*cross.pos[0]+
						cross.pos[1]*cross.pos[1]+
						cross.pos[2]*cross.pos[2]);
	return atan2(sinA,cosA);
}

double mVector::angleDeg(const mVector &v) const
{
	return angleRad(v) * mRad2Deg;
}

mVector &mVector::operator *= (mTensor2 &other)
{
  mVector m(*this);
  pos[0] = pos[1] = pos[2] = 0.0;
  for(int i=0;i<3;i++)
    {
      for(int j=0;j<3;j++)
	{
	  pos[i] += other(i,j) * m.pos[j];
	}
    }
  return *this;
}

