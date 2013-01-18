// mTensor2.h: interface for the mTensor2 class.
//
//////////////////////////////////////////////////////////////////////
#ifndef _MTensor2_H_
#define _MTensor2_H_

class mVector;
class mTensor2  
{
	double pos[3][3];
public:
	mTensor2(){}
	mTensor2(double);
	mTensor2(const mVector &v1, const mVector &v2, const mVector &v3);
	mTensor2 operator- (const mTensor2 &other);
	mTensor2 operator+ (const mTensor2 &other);
	mTensor2 & operator += (const mTensor2 &other);
	mTensor2 & operator -= (const mTensor2 &other);
	mVector operator * (const mVector &other);
	mTensor2 operator * (const mTensor2 &other);
	double & operator() ( int i, int j){return pos[i][j];}
};

mTensor2 operator / (const mVector &p1, const mVector &p2); 
#endif
