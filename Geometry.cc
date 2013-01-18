#include "Geometry.h"
#include "mPoint.h"
#include "DGCell.h"
#include "Vertex.h"
#include "mEntity.h" 
#include "gEntity.h"
#include <math.h>
#include <list>

using namespace std;

DGCell* Geometry::findReflPoint(DGCell *cell, double x,double y,double z,double &xf,double &yf,double &zf)
{
  double r =  sqrt(x*x+y*y);
  double R = getR(x,y);
  double n0 = x/r;
  double n1 = y/r;
  double n2 = z/r;
  double dist = R-r;
  
  xf = x + 2.*n0*dist;
  yf = y + 2.*n1*dist;
  zf = z + 2.*n2*dist;
  mPoint p(xf,yf,zf);
  double val[MaxNbEqn];
  list<mEntity* > allNeigh;
  list<mEntity*>::const_iterator it,it_neigh;
  mUpwardAdjacencyContainer::iter itt;

  if (cell->inCell(p,val)) return cell;
  else {
    int Size0 = cell->getMeshEntity()->size(0);
    for(int i=0;i<Size0;++i)
      {
	mEntity *ent = cell->getMeshEntity()->get(0,i);
	mUpwardAdjacencyContainer::iter begin = ent->begin(2);
	mUpwardAdjacencyContainer::iter end = ent->end(2);
	for (itt= begin;itt!=end ; ++itt)
	  {
		mEntity *m = (*itt);
		allNeigh.remove(m);
	    allNeigh.push_back(m);
	    DGCell *neighbor = (DGCell*)(*itt)->getCell();
	    if (neighbor->inCell(p,val)) return neighbor;
	  }
      }
    
    //didn't found in closest cells
    int n = (int) allNeigh.size();   
    for(it = allNeigh.begin();it!=allNeigh.end();++it)
      {
	DGCell *Cell = (DGCell*) (*it)->getCell();
	int Size0 = Cell->getMeshEntity()->size(0);
	
	for(int i=0;i<Size0;++i)
	  {
	    mEntity *ent = Cell->getMeshEntity()->get(0,i);
	    
	    mUpwardAdjacencyContainer::iter begin = ent->begin(2);
	    mUpwardAdjacencyContainer::iter end = ent->end(2);
	    for (itt= begin;itt!=end ; ++itt)
	      {
			mEntity *m = (*itt);
	//		m->print();
		DGCell *neighbor = (DGCell*)(*itt)->getCell();
		if (neighbor->inCell(p,val))  return neighbor;
		bool tag = 0;
for (it_neigh= allNeigh.begin();it_neigh!=allNeigh.end() ; ++it_neigh)	
if (*it_neigh ==m) tag =1;
		if (!tag) allNeigh.push_back(m);

	      }
	//	printf("\n");
	  }
	
      }
  }
  return 0;
}

double Geometry::getR(double x,double y)
{
  if(x*x+y*y>1.384*1.384) return 1.384;
  else return 1.;

//  return 0.5;
}
