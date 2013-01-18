#ifndef _CUTCELL_H_
#define _CUTCELL_H_

#include <vector>
#include <list>
#include "DGCell.h"
#include "mPoint.h"
#include "mVector.h"
#include "mAttachableDataContainer.h"
#include "Constants.h"
#include "FunctionSpace.h"
#include "ConservationLaw.h"
#include "Geometry.h"

class mEntity;
class Mapping;
class FunctionSpace;
class FieldEvaluator;
class DGAnalysis;
class GaussIntegrator;

using namespace std;

class cutCellGaussPoint 
{
 public:
  bool isInsideDomain;
  double*  fctsAtReflPoint;
  DGCell* reflCell;
  mVector N;
  cutCellGaussPoint() {}
  ~cutCellGaussPoint(){}
  cutCellGaussPoint(const cutCellGaussPoint &);
};

class CutCell : public DGCell
{
protected:
  vector<cutCellGaussPoint>  cutCellGaussPoints;
  virtual void init();
  Geometry *theGeometry;
  public:
  CutCell();
  CutCell(ConservationLaw*, mEntity*, FunctionSpace*,FunctionSpace *, Mapping *,GaussIntegrator *,Geometry*);
  virtual ~CutCell();
  virtual void L2ProjCutCell();
  virtual void interpolateRefl(const double* fct, mVector &N, DGCell* reflCell, double *Q) const 
  { 
    reflCell->interpolate(fct,Q);
    double uright_T,uright_N;
    uright_N = -(N(0)*Q[1]+N(1)*Q[2]);
    uright_T =  (N(1)*Q[1]-N(0)*Q[2]);
    Q[1] = uright_N*N(0) + uright_T * N(1);
    Q[2] = uright_N*N(1) - uright_T * N(0);
  }
}; 


#endif











