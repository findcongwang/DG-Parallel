#include "FieldEvaluator.h"
#include "DGCell.h"
#include "FunctionSpace.h"
#include <stdio.h>
#include <iostream>
DGCellFieldEvaluator::DGCellFieldEvaluator()
{
}

void DGCellFieldEvaluator::addCell(DGCell* c)
{
  theCells.push_back(c);
}

void DGCellFieldEvaluator::eval (const mPoint &space, double time, double *val) const
{
  for(list<DGCell*>::const_iterator it = theCells.begin();it!=theCells.end();++it)
    {
      // double u,v,w;
      //if((*it)->inCell(space,val,u,v,w))
      if((*it)->inCell(space,val))
	{
	  return;
	}
    }
  cout << "error in DGCellfieldEvaluator\n";
}

int DGCellFieldEvaluator::order () const
{
  int ord = 0;
  for(list<DGCell*>::const_iterator it = theCells.begin();it!=theCells.end();++it)
    {
      int o = (*it)->getFunctionSpace()->order();
      ord = (o>ord)?o:ord;
    }
  return ord;
}

