// mVector.h: interface for the mVector class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _MVECTOR_H_
#define _MVECTOR_H_

class mPoint;
class mTensor2;

class mVector  
{
	double pos[3];
public:
	mVector() {}
	~mVector() {}
	mVector(double x , double y, double z);
	mVector(const mPoint &from, const mPoint &to);
	mVector operator- (const mVector &other);
	mVector operator+ (const mVector &other);
	mVector & operator += (const mVector &other);
	mVector & operator *= (const double &other)
	  {
	    pos[0]*=other;pos[1]*=other;pos[2]*=other;
	    return *this;
	  }
	mVector & operator /= (const double &other);
	mVector & operator -= (const mVector &other);
	double operator * (const mVector &other)
	{
	  return pos[0] * other.pos[0] + pos[1] * other.pos[1] + pos[2]*other.pos[2];
	}
	mVector operator % (const mVector &other);
	mVector operator * (const double &other);
	mVector & operator *= (mTensor2 &other);
	mVector operator / (const double &other);
	double & operator () (int i) {return pos[i];}
	double & operator [] (int i) {return pos[i];}
	double angleRad (const mVector &v) const;
	double angleDeg (const mVector &v) const;
	void normalize();
	double L2norm() const;
	friend class mTensor2;
};
/*
mVector operator - (const mPoint &p1, const mPoint &p2); 
*/
#endif 
