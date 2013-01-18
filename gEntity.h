/*Class gEntity.h keeps track of specially marked parts of the
pysical domain. 
 Id == identification number of a physical entity. For example,
part of the boundary, where to impose boundary conditions or part 
of the domain we want to tag for some special purpose. The number 
comes from the meshing software.
level:
0 -- point
1 -- edge
2 -- 2D element
3 -- 3D element
*/

#ifndef _GENTITY_H_
#define _GENTITY_H_

//typedef CUT_CELL 1000
class gEntity  
{
 protected:
  int iD;
  int level;
 public:
  gEntity(int i, int l):iD(i),level(l){}
  int getId() const {return iD;}
  int getLevel() const {return level;}
};

class gEntityLessThanKey
{
  public :
    bool operator () (gEntity *g1, gEntity *g2) const
    {
      if(g1->getLevel() < g2->getLevel())return true;
      if(g1->getLevel() > g2->getLevel())return false;
      if(g1->getId() < g2->getId())return true;
      return false;
    }
};
#endif
