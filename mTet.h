// mFace.h: interface for the mFace class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _MTET_H_
#define _MTET_H_
#include "mRegion.h"
class Edge;
class Vertex;
class Face;

class mTet : public mRegion
{
public:
  mTet (Edge *e1, Edge *e2, Edge *e3, Edge *e4, Edge *e5, Edge *e6, gEntity *classification);
  mTet (Face *f1, Face *f2, Face *f3, Face *f4, gEntity *classification);
  mTet (vector<Vertex*>& v, int geomO, gEntity *classification);
  virtual mType getType() const{return TET;} 
  virtual int getNbTemplates (int what) const;
  virtual mEntity* getTemplate (int ith , int what , int with) const;
};

#endif 
