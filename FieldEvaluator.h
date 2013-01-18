#ifndef _FIELDEVALUATOR_H_
#define _FIELDEVALUATOR_H_

#include<list>
#include "mCompiler.h"

class mPoint;
class DGCell;
class FieldEvaluator
{
  public :
    virtual void eval (const mPoint &space, double time, double *val) const = 0;
    virtual int order () const = 0;
};

class DGCellFieldEvaluator : public FieldEvaluator
{
    list<DGCell*> theCells;
  public :
    DGCellFieldEvaluator();
    void addCell(DGCell*);
    virtual void eval (const mPoint &space, double time, double *val) const;
    virtual int order () const;
};

#endif
