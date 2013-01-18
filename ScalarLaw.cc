#include "ScalarLaw.h"
#include "mVector.h"
#include "mPoint.h"
#include "FieldEvaluator.h"
#include <stdio.h>
#include <math.h>
#include <cstring>

/*
Scalar Law solves dt(u) + v grad (u) = 0
with v = (vx,vy,vz)
*/

/*
Constructor takes as parameter a FieldEvaluator
that you have to provide in the main
*/

ScalarLaw::ScalarLaw (FieldEvaluator *f, FieldEvaluator *v, FieldEvaluator *e)
{
  farField = f;
  velocityField = v;
  exactField = e;
}

/*Only one flux because nbFields() = 1 (one unknown in Q,
  Q is of size 1)
*/
void ScalarLaw::Fi ( const mPoint &position, double *Q, mVector *flux) const
{
  double v[3];
	double crap = 0.0;
       	velocityField->eval(position,crap,v);
	flux[0](0) = v[0] * Q[0];
	flux[0](1) = v[1] * Q[0];
	flux[0](2) = v[2] * Q[0];
}

void ScalarLaw::boundaryFlux (mVector &n, int BoundaryId, const mPoint &position, 
			    double *Q, double *flux, double T) const
{	
	double v[3];
  velocityField->eval(position,T,v);
  double normalVelocity = (n(0) * v[0] + n(1) * v[1] + n(2) * v[2]);


  //entering the domain v*n < 0 so that my flux (v.n)*u
  // uses u as the far field
  if(normalVelocity < 0)
    {
      double val[1];
      farField->eval(position,T,val);
      flux[0] = normalVelocity * val[0];
    }
  else 
    {
      flux[0] = normalVelocity*Q[0];
    }
}

void ScalarLaw::normalVelocity(mVector &n, const mPoint &position, double T, double *field, double &normalVelocity) const
{
	double v[3];
  velocityField->eval(position,T,v);
  normalVelocity = (n(0) * v[0] + n(1) * v[1] + n(2) * v[2]);
 // printf("normal %f\n", normalVelocity); 
}

void ScalarLaw::riemannSolver( mVector &n , 
			       const mPoint &position,
			       double *u0, 
			       double *u1, 
			       double *flux) const
{
  double v[3];
  double crap = 0;
  velocityField->eval(position,crap,v);
  double vn = (n(0) * v[0] + n(1) * v[1] + n(2) * v[2]);

  double u[1];

  if(vn > 0) u[0] = u0[0];
  else u[0] = u1[0];

  mVector f;
  Fi(position,u,&f);
  flux[0] = f * n;
}

double ScalarLaw::maximumEigenValue (const double *theMean, const mPoint &position) const
{
  double v[3];
  double dummy = 0.0;
  velocityField->eval(position,dummy,v);
  //return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  return fabs(v[0])+fabs(v[1])+fabs(v[2]);
}

double ScalarLaw::maxEigenValueCube(double du,double dv,double dw,double *theMean, const mPoint &p) const
{
  double v[3];
  double dummy = 0.0;
  velocityField->eval(p,dummy,v);
  //return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  return fabs(v[0])+fabs(v[1])+fabs(v[2]);
}

/********** P O S T   P R O **********/

#include <stdio.h>
int ScalarLaw::getNbFieldsOfInterest() const
{
  return 5;
}

void ScalarLaw::getNameAndSize(int i, int &n, char *name) const
{
  switch(i)
    {
    case 0:
      n = 1;
      strcpy(name,"error");
      break;
    case 1:
      n = 1;
      strcpy(name,"exact");
      break; 
    case 2:
      n = 1;
      strcpy(name,"numEr");
      break; 
	case 3:
      n = 1;
      strcpy(name,"function");
	  break;
  case 4:
      n = 1;
      strcpy(name,"limit");
      break; 
    }
} 

void ScalarLaw::getIthFieldOfInterest(int i, const mPoint &p, double *field, double *res, double time) const
{
  switch(i)
    { 
    case 0:
      double exact[1];
      exactField->eval(p,time,exact);
      res[0] = exact[0] - field[0];
      break;
    case 1:
      double ex[1];
      exactField->eval(p,time,ex);
      res[0] = ex[0];
      break;
    case 2:
      res[0] = field[0];
      break; 
	case 3:
      res[0] = field[0];
      break;
	case 4:
	  res[0] = field[0];
	  break;
    }
}

int ScalarLaw::getEdgeOrientation( mVector &n , 
				    const mPoint &position, double *field) const
{
  //might depend on time ... should be changed
  double v[3];
  double crap = 0;
  velocityField->eval(position,crap,v);
  double vn = (n(0) * v[0] + n(1) * v[1] + n(2) * v[2]);
  //printf("vn %f \n", vn);
  // printf("v %f %f %f \n", v[0],v[1],v[2]);
  if (vn > 0) return (-1);
  else return (1);
}


void ScalarLaw::RHS (double *field, const mPoint &position, double t, double *rhs) const
{
  rhs[0] =  0.0;
  //rhs[0] = 6.*(position(1)-2.);
   ///rhs[0]=4./ field[0]/2 ;
  //rhs[0]= (3.*position(0)*position(0) +2.*position(1)*position(1)+ 1. );
 // rhs[0]=1.0;
  // printf("rhs %f\n", rhs[0]);
  
}

bool ScalarLaw::isPhysical ( double *Q ) const
{
  //if(Q[0] < 0 && Q[0] > 1.0)return false;
  return true;
}

double ScalarLaw::jumpQuantity ( double *Q ) const
{
  return Q[0];
}













