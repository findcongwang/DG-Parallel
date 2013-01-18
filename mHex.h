// mHex.h: interface for the mHex class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _MHEX_H_
#define _MHEX_H_
#include "mRegion.h"
class Edge;
class Vertex;
class Face;
class mHex : public mRegion
{
 public:
	 mHex (){}
  mHex (Vertex *v1, 
	Vertex *v2, 
	Vertex *v3, 
	Vertex *v4, 
	Vertex *v5, 
	Vertex *v6, 
	Vertex *v7, 
	Vertex *v8, 
	gEntity *classification);
  mHex (Edge *e1, Edge *e2, Edge *e3, Edge *e4, Edge *e5, Edge *e6, 
	Edge *e7, Edge *e8, Edge *e9, Edge *e10, Edge *e11, Edge *e12, 
	gEntity *classification);
  mHex (Face *f1, Face *f2, Face *f3, Face *f4, Face *f5, 
	Face *f6, gEntity *classification);
  virtual mType getType() const{return HEX;} 
  virtual int getNbTemplates (int what) const;
  virtual mEntity* getTemplate (int ith , int what , int with) const;
  virtual void print() const;
};

#endif 
