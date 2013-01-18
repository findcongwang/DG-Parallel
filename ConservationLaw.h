#ifndef _CONSERVATIONLAW_H_
#define _CONSERVATIONLAW_H_

// Base class for conservation laws
// Possible Inheritances : euler, burgers, linear advection ... 
// Possible riemannSolvers : Roe, Exact, ...

class mVector;
class mPoint;
class FieldEvaluator;

class ConservationLaw 
{
protected:
  FieldEvaluator *farField;
  FieldEvaluator *exactField;
public :
  virtual double computeCFL(const int order, const double *Mean, const mPoint &) const {return 0;};
  // returns the far field (boundary conditions + initial field)
  FieldEvaluator *getFarField(){return farField;}
  FieldEvaluator *getExactSolution(){return exactField;}
  // number of fields
  virtual int getNbFields() const = 0;
  //Fluxes as a matrix
  virtual void Fi  (const mPoint &position,double *, mVector *flux) const = 0;
  // Right Hand Side
  virtual void RHS (double *field, const mPoint &position, double time, double *rhs) const = 0;
  // numerical flux, give the normal of the interface and
  // left and right fluxes, output the numerical flux
  // if left or right are null, we apply boundary conditions.
  virtual void riemannSolver ( mVector &n , 
			       const mPoint &position,
			       double *left, 
			       double *right, 
			       double *flux) const = 0;
  // boundary fluxes
  virtual void boundaryFlux (mVector &n, int BoundaryId, const mPoint &position, 
			     double *field, double *flux, double) const = 0;
  virtual void boundary(mVector &n, mVector &N, int BoundaryId, const mPoint &position, 
			double *field, double *flux, double) const{}
  // maximum eigenvalue
  virtual double maximumEigenValue (const double *field, const mPoint &position) const = 0;
  virtual double maximumEigenValue (mVector &n, const double *field1, const double* field2) const=0;
  virtual double maxEigenValueCube(double du,double dv,double dw,double *theMean, const mPoint &p) const = 0;
  virtual void computeJacobian(mVector& n, const double* Q,double** Jac) const = 0;
  // number of fields of interest, for post processing
  virtual int getNbFieldsOfInterest() const = 0;
  virtual void getNameAndSize(int i, int &n, char *name) const = 0;
  virtual void getIthFieldOfInterest(int i, const mPoint &position, double *field, double *res, double T) const= 0;
  // a function saying if values are physical or not (too much overshoots)
  virtual bool isPhysical ( double *val ) const = 0;
  // a function giving the actual jump quantity relevant for error estimation
  virtual double jumpQuantity (double *val) const = 0;
   virtual int getEdgeOrientation( mVector &,const mPoint &, double *) const=0;
   virtual void normalVelocity(mVector &n, const mPoint &position, double T, double *,
							   double &normalVelocity) const=0;
   virtual void compute_left_eigenvectors(double *U, mVector &n, double *rev) const{}
   virtual void compute_right_eigenvectors(double *U, mVector &n, double *rev) const{}
   friend class DGBarthLimiterEuler;
};

#endif
