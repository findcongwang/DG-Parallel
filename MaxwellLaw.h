#ifndef _MAXWELLLAW_H_
#define _MAXWELLLAW_H_

#include <set>
#include "ConservationLaw.h"
#include "mPoint.h"
#include <math.h>
using namespace std; 

class Maxwell : public ConservationLaw
{
 protected:
  set<int> Walls;
  public :
    Maxwell (set<int> &Walls, FieldEvaluator *);
  virtual int getNbFieldsOfInterest() const ;
  virtual void getNameAndSize(int i, int &n, char *name) const ;
  virtual int getEdgeOrientation( mVector &,const mPoint &,double *) const {return 0;}
  virtual void normalVelocity(mVector &n, const mPoint &position, double T, double *,
			      double &normalVelocity) const {};
};

class Maxwell2d : public Maxwell
{
  public :
    Maxwell2d (set<int> &Walls, FieldEvaluator *);
  virtual double computeCFL(const int order, const double *Mean, const mPoint &) const ;
  // number of fields
  virtual int getNbFields() const {return 3;}
  //flux as matrix
  virtual void Fi  (const mPoint &position , double *, mVector *flux) const;
  // Right Hand Side
  virtual void RHS (double *field, const mPoint &position, double time,double *rhs) const ;
  // numerical flux, give the normal of the interface and
  // left and right fluxes, output the numerical flux
  // if left or right are null, we apply boundary conditions.
  virtual void riemannSolver (mVector &n ,
			      const mPoint &position, 
			      double *left, 
			      double *right,
			      double *flux) const;
  // boundary fluxes
  virtual void boundaryFlux (mVector &n, int BoundaryId, const mPoint &position, 
			     double *field, double *flux, double T) const;
  void boundary(mVector &n, mVector &N, int BoundaryId, const mPoint &position, 
		double *field, double *flux, double T) const;
  void boundarySource(mVector &n, mVector &N, int BoundaryId, const mPoint &position,     //"fake" BC by Dale
		double *field, double *flux, double T) const;
  virtual double maximumEigenValue (const double *field, const mPoint &position) const;
  virtual double maximumEigenValue (mVector &n, const double *field1, const double* field2) const;
  virtual double maxEigenValueCube(double du,double dv,double dw,double *theMean, const mPoint &p) const;
  virtual void computeJacobian(mVector& n, const double* Q,double** Jac) const;
  void getIthFieldOfInterest(int i, const mPoint &, double *field, double *res, double time) const ;
  virtual bool isPhysical ( double *Q ) const;
  virtual double jumpQuantity ( double *Q ) const {return Q[0];}
  virtual int getEdgeOrientation( mVector &n , const mPoint &position, double *val) const;
  virtual void compute_left_eigenvectors(double *U, mVector &n, double *rev) const;
  virtual void compute_right_eigenvectors(double *U, mVector &n, double *rev) const;
};

#endif



