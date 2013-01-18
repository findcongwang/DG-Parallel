// mTensor2.cpp: implementation of the mTensor2 class.
//
//////////////////////////////////////////////////////////////////////

#include "mTensor2.h"
#include "mVector.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

mTensor2::mTensor2(double init)
{
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
	  pos[i][j] = init;
}


mVector mTensor2::operator * (const mVector &other)
{
  mVector m(0,0,0);
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
	  m.pos[i] += pos[i][j] * other.pos[j];
return m;
}

mTensor2 mTensor2::operator * (const mTensor2 &other)
{
  mTensor2 m(0);
  for(int i=0;i<3;i++)
      for(int j=0;j<3;j++)
			for(int l=0;l<3;l++)
				m.pos[i][j] += pos[i][l] * other.pos[l][j];
  return m;
}

