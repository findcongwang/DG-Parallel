#ifndef _ScalarLAW_H_
#define _ScalarLAW_H_
#include "ConservationLaw.h"

class ScalarLaw : public ConservationLaw
{
  FieldEvaluator *velocityField;
  // FieldEvaluator *exactField;
 public :
  ScalarLaw (FieldEvaluator *, FieldEvaluator *, FieldEvaluator *);
  virtual int getNbFieldsOfInterest() const ;
  virtual void getNameAndSize(int i, int &n, char *name) const ;
  // number of fields
  virtual int getNbFields() const {return 1;}
  // ith flux 
  virtual void Fi  (const mPoint &position,double *, mVector *flux) const;
  // Right Hand Side
  virtual void RHS (double *field, const mPoint &position, double time, double *rhs) const;
  // numerical flux, give the normal of the interface and
  // left and right fluxes, output the numerical flux
  // if left or right are null, we apply boundary conditions.
  virtual void riemannSolver ( mVector &n , const mPoint &position, double *left, 
			       double *right, 
			       double *flux) const;
  virtual int getEdgeOrientation( mVector &n , const mPoint &position,double *) const;
  // boundary fluxes
  virtual void boundaryFlux (mVector &n, int BoundaryId, const mPoint &position, 
			     double *field, double *flux, double T) const;
  virtual double maximumEigenValue (const double *field, const mPoint &position) const;
  virtual double maximumEigenValue (mVector &n, const double *field1, const double* field2) const {return 0.0;}
  virtual double maxEigenValueCube(double du,double dv,double dw,double *theMean, const mPoint &p) const;
  virtual void computeJacobian(mVector& n, const double* Q,double** Jac) const {}
  void getIthFieldOfInterest(int i, const mPoint &, double *field, double *res, double time) const ;
  virtual bool isPhysical ( double *val ) const;
  virtual double jumpQuantity ( double *val ) const;
  virtual void normalVelocity(mVector &n, const mPoint &position, double T, double *,
							   double &normalVelocity) const;
};


#endif





