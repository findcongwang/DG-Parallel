// Face.h: interface for the Face class.
//
//////////////////////////////////////////////////////////////////////

#ifndef _FACE_H_
#define _FACE_H_
#include "mEntity.h"
class gEntity;
class Vertex;
class Edge;
// The Face is implemented generally
// The assumption is that the number
// of edges is equal to the number of
// vertices and that 2 consecutive edges
// share one vertex
// Special constructors for triangles & quads

//Grant: This assumption has been broken upon addition of higher geometric order.
class Face : public mEntity 
{
protected :
	mType theType;
public:
  virtual ~Face();
  Face (mDownwardAdjacencyContainer *lis, gEntity *classification);
  Face (Vertex *v1, Vertex *v2, Vertex *v3, gEntity *classification); //kept for splitting
  Face (Edge *v1, Edge *v2, Edge *v3, gEntity *classification);
  Face (Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4, gEntity *classification);
  Face (Edge *v1, Edge *v2, Edge *v3, Edge *v4, gEntity *classification);
  //Generalized constructor. Not yet implemented fully.
  Face(vector<Vertex*>& v,int Type,int geomO, gEntity *classification);
  Vertex *commonVertex (Face *f1, Face *f2);
  Edge *commonEdge (Face *f1);
  virtual int getLevel() const {return 2;} 
  virtual mType getType() const{return theType;} 
  virtual int getNbTemplates (int what) const;
  virtual mEntity* getTemplate (int ith , int what , int with) const;
  virtual void print() const;
};


#endif 
