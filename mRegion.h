// mRegion.h: interface for the mRegion class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _MREGION_H_
#define _MREGION_H_

#include "mEntity.h"
#include "Edge.h"
#include "Face.h"

class mRegion  : public mEntity
{
public:
  virtual ~mRegion()
    {
      for(int i=0;i<4;i++)
	  if(i <= getLevel() && theAdjacencies[i])
		for(int j=0;j<theAdjacencies[i]->size();j++)
		  theAdjacencies[i]->get(j)->del(this);
    }

  int getLevel()const{return 3;};
 Edge *commonEdge (mRegion *r1, mRegion *r2) 
  {//not done yet;
	return 0;}
  Face *commonFace (mRegion *r)
{
  for(int i=0;i<size(2);i++)
    {
      mEntity *f = get(2,i);
      for(int j=0;j<r->size(2);j++)
	     if(f == r->get(2,j))return (Face*)f;    
    }
  return 0;
}
};

#endif 
