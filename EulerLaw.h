#ifndef _EULERLAW_H_
#define _EULERLAW_H_

#include <set>
#include "ConservationLaw.h"
#include "mPoint.h"
#include <math.h>
using namespace std;

class Euler : public ConservationLaw
{
 protected:
  set<int> Walls;
  const double gamma;
  int  riemann( double rho1,                /* inputs */
		double  u1,
		double  v1,
		double  w1,
		double  p1,
		double rho2,
		double  u2,
		double  v2,
		double  w2,
		double  p2,
		double *rhoav,               /* outputs */
		double  *uav,
		double  *utav,
		double  *uttav,
		double  *pav ) const;
  public :
    Euler (set<int> &Walls, FieldEvaluator *, double GAMMA=1.4);
  virtual int getNbFieldsOfInterest() const ;
  virtual void getNameAndSize(int i, int &n, char *name) const ;
  virtual int getEdgeOrientation( mVector &,const mPoint &,double *) const {return 0;}
  virtual void normalVelocity(mVector &n, const mPoint &position, double T, double *,
			      double &normalVelocity) const {};
};

class Euler2d : public Euler
{
  public :
    Euler2d (set<int> &Walls, FieldEvaluator *, double GAMMA=1.4);
  virtual double computeCFL(const int order, const double *Mean, const mPoint &) const ;
  // number of fields
  virtual int getNbFields() const {return 4;}
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

class Euler2dRoe : public Euler2d
{  
 public:
  Euler2dRoe (set<int> &Walls, FieldEvaluator *, double GAMMA=1.4);
  virtual void riemannSolver (mVector &n , 
			       const mPoint &position,
			       double *left, 
			       double *right, 
			       double *flux) const;
};

class Euler3d : public Euler
{
 protected:
  public :
    Euler3d (set<int> &Walls, FieldEvaluator *, double GAMMA=1.4);
  // number of fields
  virtual int getNbFields() const {return 5;}
  // ith flux 
  virtual void Fi  (const mPoint &position, double *, mVector *flux) const;
  // Right Hand Side
  virtual void RHS (double *field, const mPoint &position, double time, double *rhs) const;
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
  virtual double maximumEigenValue (const double *field, const mPoint &position) const;
  virtual double maximumEigenValue (mVector &n, const double *field1, const double* field2) const;
  virtual void computeJacobian(mVector& n, const double* Q,double** Jac) const ;
  virtual double maxEigenValueCube(double du,double dv,double dw,double *theMean, const mPoint &p) const;
  void getIthFieldOfInterest(int i, const mPoint &, double *field, double *res, double time) const ;
  virtual bool isPhysical ( double *val ) const;
  virtual double jumpQuantity ( double *Q ) const {return Q[0];}
  
};

class EulervanLeer : public Euler2d
{  
 public:
  EulervanLeer (set<int> &Walls, FieldEvaluator *, double GAMMA=1.4);
  virtual void riemannSolver (mVector &n , 
			       const mPoint &position,
			       double *left, 
			       double *right, 
			       double *flux) const;
};

class Euler2dExactRiemann : public Euler2d
{
 public:
  Euler2dExactRiemann (set<int> &Walls, FieldEvaluator *, double GAMMA=1.4);
  virtual void riemannSolver (mVector &n , 
			       const mPoint &position,
			       double *left, 
			       double *right,
			       double *fluxr) const;
};

class Euler2dHLLC: public Euler2d
{
 public:
  Euler2dHLLC (set<int> &Walls, FieldEvaluator *, double GAMMA=1.4);
  virtual void riemannSolver (mVector &n , 
			      const mPoint &position,
			      double *left, 
			      double *right,
			      double *fluxr) const;
};

class Euler3dRoe : public Euler3d
{
 public:
  Euler3dRoe (set<int> &Walls, FieldEvaluator *, double GAMMA=1.4);
  virtual void riemannSolver (mVector &n , 
			      const mPoint &position,
			      double *left, 
			      double *right, 
			      double *flux) const;
};
#endif



