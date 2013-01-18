#include "DGAnalysis.h"
#include "FieldEvaluator.h"
#include "mPoint.h"
#include "Burgers.h"
#include "DGLimiter.h"
#include "DGSensor.h"
#include "mDGMesh.h"
#include "Constants.h"
#include <stdio.h>
#include <math.h>
#include <cstring>
#include "mVector.h"
#define TOL 5e-15


Burgers2D::Burgers2D(FieldEvaluator *f, FieldEvaluator *e)
{
  farField = f;
  exactField = e;
}

void Burgers2D::Fi ( const mPoint &position , double *Q, mVector *flux) const

{
	flux[0](0) = Q[0] * Q[0] * 0.5;
    flux[0](1) = Q[0] * Q[0] * 0.5;
	flux[0](2) = 0.0;
}

void Burgers2D::boundaryFlux (mVector &n, int BoundaryId, const mPoint &position, 

			    double *Q, double *flux, double T) const

{

  double normalVelocity = Q[0] *  Q[0] * 0.5 *(n(0) + n(1));

  if(normalVelocity < 0)

    {

      double val[1];

      farField->eval(position,T,val);

      flux[0] = val[0] *  val[0] * 0.5 * (n(0) +n (1));

    }

  else 

    flux[0] = normalVelocity;

}



void Burgers2D::riemannSolver( mVector &n , 

			     const mPoint &position,

			     double *u0, 

			     double *u1,
			     

			     double *flux) const

{
/*
  double vn = u0[0]* u0[0] * 0.5 * (n(0)+ n(1));
  double u[1];

  if(vn > 0)u[0] = u0[0];

  else u[0] = u1[0];

  flux[0] = u[0] * u[0] * 0.5 * ( n(0)+n(1));
 */
 
  double maxIV=0;
	flux[0] =0.5 * ( (u0[0] * u0[0] + u1[0] * u1[0] ) * (n(0) + n(1) ) *0.5 -
	  maxIV * (u1[0] - u0[0]) ); 
}



double Burgers2D::maximumEigenValue (const double *theMean, const mPoint &position) const

{

  return 1;

}

double Burgers2D::maxEigenValueCube(double du,double dv,double dw,double *theMean, const mPoint &p) const
{
return 1.;
}

void Burgers2D::normalVelocity(mVector &n, const mPoint &position, double T,

			     double *Q, double &normalVelocity) const 

{

  normalVelocity = Q[0] * (n(0) + n(1));

  //printf("normal vel %f \n", normalVelocity);

}



/********** P O S T   P R O **********/



#include <stdio.h>

int Burgers2D::getNbFieldsOfInterest() const

{

  return 4;

}



void Burgers2D::getNameAndSize(int i, int &n, char *name) const

{

  switch(i)

    {

    case 0:

      n = 1;

      strcpy(name,"function");

      break;

    case 1:

      n = 1;

      strcpy(name,"error");

      break;

    case 2:

      n = 1;

      strcpy(name,"exact");

      break; 

    case 3:

      n = 1;

      strcpy(name,"numEr");

      break; 

    }

} 



void Burgers2D::getIthFieldOfInterest(int i, const mPoint &p, double *field, double *res, double time) const

{

  switch(i)

    {

    case 0:

      res[0] = field[0];

      break;

    case 1:

      double exact[1];

      exactField->eval(p,time,exact);

      res[0] = exact[0] - field[0];

      break;

    case 2:

      double ex[1];

      exactField->eval(p,time,ex);

      res[0] = ex[0];

      break;

    case 3:

      res[0] = field[0];

      break; 

    }

}



int Burgers2D::getEdgeOrientation( mVector &n , 

				    const mPoint &position,double *field) const

{

  //might depend on time ... should be changed

  double vn = field[0] * field[0]/2 *(n(0) + n(1));

  //printf("vn %f \n", vn);

  // printf("n %f %f %f \n", n(0),n(1),n(2));

  //outflow edge

  if (vn > 0) return (-1);

  //inflow edge

  else return (1);

}





void Burgers2D::RHS (double *field, const mPoint &position, double t, double *rhs) const

{

  rhs[0]=0;    

}



bool Burgers2D::isPhysical ( double *Q ) const

{

    return true;

}



double Burgers2D::jumpQuantity ( double *Q ) const

{

  return Q[0];

}





/**************************************************************************/
/***************** Initial and Boundary Conditions*************************/

class Exact : public FieldEvaluator
{
public:
  void eval(const mPoint &p, double time, double *val) const
  {
    double alpha = 0.5;
  double beta = 0.5;
  if (p(1) < 0.0001) val[0] = alpha + beta * sin(mPi * (p(0)+p(1)));
else {
  double z;
  double x0, y0,t0;
  double a, b,c,d,w;
  double half_int, f_z, f_a, f_c, f_w;
  int i;
  
  t0 = beta * time;
  x0 = p(0) - alpha * time;
  y0 = p(1) -alpha * time;
  
  if (x0 < 0) x0 = -x0;
  if (y0 < 0) y0 = -y0;
  
  a = -1.5;
  b =  1.5;
  c = -1.5;
  d =  1.5;

  i = 0;
  while (1)
    {half_int = (b - a) / 2.0;
    z = a + half_int;
    f_z = z - sin(mPi * (x0 - z * t0));
	f_w = w - sin(mPi * (y0 - w * t0));
    f_a = a - sin(mPi * (x0 - a * t0));
	f_c = c - sin(mPi * (x0 - c * t0));
    
    if ((f_z == 0.0 && f_w==0) || (half_int < TOL))
      break;
    else if ((f_z * f_a) > 0)
      a = z;
    else
      b = z;
	if ((f_w * f_c) > 0)
      c = w;
    else
      d = w;
    }
  
  if ((p(0) - alpha * p(0)) < 0) z = -z;
  if ((p(1) - alpha * p(1)) < 0) w = -w;
  
  val[0] = alpha + beta * z;
}
}
  int order() const
  {
    return 4;
  }
};

class Ex : public FieldEvaluator
{
public:
  void eval(const mPoint &p, double time, double *val) const
  {
	  val[0]=1;
    //if (p(0)< 0.5 && p(1) <0.5 ) val[0] = 0.5;
	 //else if (p(1) > 0.5 && p(0) >0.5 ) val[0]=-0.5;
      //else val[0]=0.25;
	 // val[0] = sin(p(0)+p(1));
//val[0] =1;
  }
  int order() const
  {
    return 1;
  }
};


class In : public FieldEvaluator
{
public:
  void eval(const mPoint &p, double time, double *val) const
  {
    if (time <  0.0001) val[0] = 1 + 0.5 * sin(3.1415926*(p(0))) * sin(3.1415926*(p(1))) ;
	else val[0]  = 1;
	  //if (p(0)< 0.5 && p(1) <0.5 ) val[0] = 0.5;
	 //else if (p(1) > 0.5 && p(0) >0.5 ) val[0]=-0.5;
    //else val[0]=0.25;
  }
  int order() const
  {
    return 1;
  }
};

void Burger2D( mDGMesh *theMesh, int order) 
{
  Ex e;
  In in;
  DGBarthLimiter myLimiter;
  Burgers2D theLaw(&in,&in);
  DGLineSensor sensor1  ("t1.dat",mPoint(0,.1,0),mPoint(.3,0,0),100,2.);
  DGLineSensor sensor2  ("t2.dat",mPoint(0,.0,0),mPoint(.3,0,0),100,2.);

  DGAnalysis analysis(theMesh,&theLaw,0,order);
  //DGAnalysis analysis(theMesh,&theLaw,&myLimiter,order);
  DGGmshSensor sensor0  ("bur",0.1);
  analysis.addSensor(&sensor0);
  analysis.addSensor(&sensor1);
  analysis.addSensor(&sensor2);
  analysis.run();
}














