#include "LinearSystem.h"
#include "mVector.h"
#include "mPoint.h"
#include "FieldEvaluator.h"
#include <stdio.h>
#include <math.h>
#include <cstring>

/*
LinearSystem solves dt(u) + (A1,A2) grad (u) = 0
*/

/*
Constructor takes as parameter a FieldEvaluator
that you have to provide in the main
*/

LinearSystem::LinearSystem (FieldEvaluator *f, FieldEvaluator *e, FieldEvaluator *a1, FieldEvaluator *a2)
{
  farField = f;
  A1 = a1;
  A2 = a2;
  exactField = e;
}

void LinearSystem::Fi ( const mPoint &position, double *Q, mVector *flux) const
{
  double a1[9];
  double a2[9];
  double crap = 0.0;
  A1->eval(position,crap,a1);
  A2->eval(position,crap,a2);
  
  flux[0](0) = a1[0]*Q[0]+a1[1]*Q[1]+a1[2]*Q[2];
  flux[0](1) = a2[0]*Q[0]+a2[1]*Q[1]+a2[2]*Q[2];
  flux[0](2) = 0.0;
  flux[1](0) = a1[3]*Q[0]+a1[4]*Q[1]+a1[5]*Q[2];
  flux[1](1) = a2[3]*Q[0]+a2[4]*Q[1]+a2[5]*Q[2];
  flux[1](2) = 0.0;
  flux[2](0) = a1[6]*Q[0]+a1[7]*Q[1]+a1[8]*Q[2];
  flux[2](1) = a2[6]*Q[0]+a2[7]*Q[1]+a2[8]*Q[2];
  flux[2](2) = 0.0;
}


void LinearSystem::Fi ( int i , const mPoint &position, double *Q, mVector &flux) const
{
  double a1[4];
  double a2[4];
  double crap = 0.0;
  A1->eval(position,crap,a1);
  A2->eval(position,crap,a2);
  
  switch(i)
    {
    case 0:
	flux(0) = a1[0]*Q[0]+a1[1]*Q[1];
	flux(1) = a1[2]*Q[0]+a1[3]*Q[1];
	flux(2) = 0.0;
      break;
	case 1:
	flux(0) = a2[0]*Q[0]+a2[1]*Q[1];
	flux(1) = a2[2]*Q[0]+a2[3]*Q[1];
	flux(2) = 0.0;
		break;
    }
}

void LinearSystem::boundaryFlux (mVector &n, int BoundaryId, const mPoint &position, 
			    double *Q, double *flux, double T) const
{
  double v[3];
  farField->eval(position,T,v);
  riemannSolver(n,position,Q,v,flux); 
}

void LinearSystem::normalVelocity(mVector &n, const mPoint &position, double T, double *field, double &normalVelocity) const
{
 printf("in not done %f\n", normalVelocity); 
}

void LinearSystem::riemannSolver( mVector &n , 
			       const mPoint &position,
			       double *u0, 
			       double *u1, 
			       double *flux) const
{
  double Q[3],diff[3];
  double a1[9],a2[9];
  double crap = 0.0;
  A1->eval(position,crap,a1);
  A2->eval(position,crap,a2);
	
  Q[0] = 0.5 * (u0[0]+u1[0]);
  Q[1] = 0.5 * (u0[1]+u1[1]);
  Q[2] = 0.5 * (u0[2]+u1[2]);
  diff[0] = 0.5 * (u1[0]-u0[0]);
  diff[1] = 0.5 * (u1[1]-u0[1]);
  diff[2] = 0.5 * (u1[2]-u0[2]);
  double lambda = maximumEigenValue(Q,position);
  lambda= 5.*(fabs(n(0)+n(1)*2.)>fabs(n(0)*2.+n(1))?fabs(n(0)+n(1)*2.):fabs(n(0)*2.+n(1)));
  double lambda_1 =fabs(n(0)*3+n(1)*2.-1);
  double lambda_2 =fabs(n(0)*3+n(1)*2.+1);
  //lambda= ( lambda_1>lambda_2 ? lambda_1:lambda_2);
  flux[0] = (a1[0]*Q[0]+a1[1]*Q[1]+a1[2]*Q[2])*n(0) + (a2[0]*Q[0]+a2[1]*Q[1]+a2[2]*Q[2])*n(1);
  flux[1] = (a1[3]*Q[0]+a1[4]*Q[1]+a1[5]*Q[2])*n(0) + (a2[3]*Q[0]+a2[4]*Q[1]+a2[5]*Q[2])*n(1);
  flux[2] = (a1[6]*Q[0]+a1[7]*Q[1]+a1[8]*Q[2])*n(0) + (a2[6]*Q[0]+a2[7]*Q[1]+a2[8]*Q[2])*n(1);
  flux[0] -= lambda*diff[0];
  flux[1] -= lambda*diff[1];
  flux[2] -= lambda*diff[2];
	
}

double LinearSystem::maximumEigenValue (const double *theMean, const mPoint &position) const
{
  double v[3];
  double dummy = 0.0;
  ///Field->eval(position,dummy,v);
  return 10.;
}

double LinearSystem::maxEigenValueCube(double du,double dv,double dw,double *theMean, const mPoint &p) const
{
return 10.;
}

/********** P O S T   P R O **********/

#include <stdio.h>
int LinearSystem::getNbFieldsOfInterest() const
{
  return 7;
}

void LinearSystem::getNameAndSize(int i, int &n, char *name) const
{
  switch(i)
    {
    case 0:
      n = 1;
      strcpy(name,"velocity");
      break;
    case 1:
      n = 1;
      strcpy(name,"error");
      break;
    case 2:
      n = 1;
      strcpy(name,"error2");
      break; 
    case 3:
      n = 1;
      strcpy(name,"exact1");
      break; 
	case 4:
      n = 1;
      strcpy(name,"exact2");
      break; 
	  case 5:
      n = 1;
      strcpy(name,"solution1");
      break; 
    case 6:
      n = 1;
      strcpy(name,"limit");
      break; 
    }
} 

void LinearSystem::getIthFieldOfInterest(int i, const mPoint &p, double *field, double *res, double time) const
{
  switch(i)
    {
    case 0:
      res[0] = field[0];res[1]=0;res[2]=0;
      break;
    case 1:
		{double exact[3];
      exactField->eval(p,time,exact);
      res[0] = 2.*(exact[0] - field[0])-(exact[1]-field[1]);}
      break;
    case 2:
		{
		double exact[3];
      exactField->eval(p,time,exact);
      res[0] = exact[1] - field[1];
		}
      break;
    case 3:
		{
	  double ex[3];
      exactField->eval(p,time,ex);
      res[0] = ex[0];}
      break; 
	case 4:
		{
	  double ex[3];
      exactField->eval(p,time,ex);
      res[0] = ex[2];
		}
      break; 
	case 5:
      res[0] = field[0];
      break; 
    case 6:
      res[0] = field[0];
      break;
    }
}

int LinearSystem::getEdgeOrientation( mVector &n , 
				    const mPoint &position, double *field) const
{
  // copied from ScalarLaw - wrong in this case
  //might depend on time ... should be changed
  double v[3];
  double crap = 0;
//  velocityField->eval(position,crap,v);
  double vn = (n(0) * v[0] + n(1) * v[1] + n(2) * v[2]);
   if (vn > 0) return (-1);
  else return (1);
}


void LinearSystem::RHS (double *field, const mPoint &position, double time, double *rhs) const
{
  double s=sin(0.5*(position(0)+position(1)));
  double c=cos(0.5*(position(0)+position(1)));
  rhs[0] =  0.5*(18.*s+21.*c);
  rhs[1] =  0.5*(21.*s+12.*c);
  rhs[0] =  0.5*(18.*field[0]+21.*field[1]);
  rhs[1] =  0.5*(21.*field[0]+12.*field[1]);
  rhs[2] =  0.0;
/*    rhs[0] = - 2.5*c +2.*(position(0) + position(1)) ;
  rhs[1] =  0.5*c -6.*position(0);
  rhs[2] =  0.5*c -4.*position(1);*/
//  rhs[0]=rhs[1]=rhs[2]=0;
}

bool LinearSystem::isPhysical ( double *Q ) const
{
  return true;
}

double LinearSystem::jumpQuantity ( double *Q ) const
{
  return Q[0];
}













