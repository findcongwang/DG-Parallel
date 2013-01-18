// Vertex.h: interface for the Vertex class.

#ifndef _VERTEX_H_
#define _VERTEX_H_

#include <iostream>
#include "mPoint.h"
#include "mEntity.h"

class mIdGenerator;
class gEntity;

using namespace std;

class Vertex : public mEntity 
{
  mPoint p;
 public:
  Vertex(int theId, const mPoint & ,gEntity *classif);
  Vertex(mIdGenerator &theIdGenerator , const mPoint &, gEntity *classif);
  mPoint point() const {return p;}
  virtual int getLevel() const {return 0;} 
  virtual mType getType() const{return VERTEX;} 
  void deleteId(mIdGenerator &theIdGenerator);
  friend ostream& operator << ( ostream &o, const Vertex &v)
    {
      o << "vertex(" << v.getId() <<")={"<<v.point()(0) <<","<<v.point()(1)<<","<<v.point()(2)<<"};\n";
      return o;
    }
  unsigned long getRAND();
  void print() const;
};

#endif 
