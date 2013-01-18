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



Burgers::Burgers(FieldEvaluator *f, FieldEvaluator *e)

{

  farField = f;
  exactField = e;

}

void Burgers::Fi ( const mPoint &position , double *Q, mVector *flux) const
{
	flux[0](0) = Q[0] * Q[0]/2;
    flux[0](1) = Q[0];
	flux[0](2) = 0.0;
}


void Burgers::boundaryFlux (mVector &n, int BoundaryId, const mPoint &position, 

			    double *Q, double *flux, double T) const

{


  double normalVelocity = Q[0] * (n(0) * Q[0]/2 + n(1));


  if(normalVelocity < 0)


    {


      double val[1];


      farField->eval(position,T,val);


      flux[0] = val[0] * (n(0) * val[0] * 0.5 +n(1));


    }


  else 


    flux[0] = normalVelocity;


}





void Burgers::riemannSolver( mVector &n , 


			     const mPoint &position,


			     double *u0, 


			     double *u1, 

			       


			     double *flux) const


{
/*

  double vn = u0[0]*(n(0) * u0[0] * 0.5 + n(1));
  double u[1];

  if(vn > 0) u[0] = u0[0];

  else u[0] = u1[0];

  flux[0] = u[0]*(u[0]*n(0)*0.5+n(1));

  */
/*
 double maxIV;
  flux[0] =0.5 * ( (u0[0] * u0[0] + u1[0] * u1[0] ) * n(0) * 0.5 + (u0[0] + u1[0]) * n(1) -

	  maxIV * (u1[0] - u0[0]) ); 
*/
	flux[0] = 0;
}





double Burgers::maximumEigenValue (const double *leftMean, const mPoint &position) const


{


return 1;


}

double Burgers::maxEigenValueCube(double du,double dv,double dw,double *theMean, const mPoint &p) const
{
	return 1.;
}



void Burgers::normalVelocity(mVector &n, const mPoint &position, double T,


			     double *Q, double &normalVelocity) const 


{


  normalVelocity = (n(0) * Q[0] + n(1));


  //printf("normal vel %f \n", normalVelocity);


}





/********** P O S T   P R O **********/





#include <stdio.h>


int Burgers::getNbFieldsOfInterest() const


{


  return 4;


}





void Burgers::getNameAndSize(int i, int &n, char *name) const


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





void Burgers::getIthFieldOfInterest(int i, const mPoint &p, double *field, double *res, double time) const


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





int Burgers::getEdgeOrientation( mVector &n , 


				    const mPoint &position,double *field) const


{
 double vn = field[0]*(n(0) * field[0] * 0.5 + n(1));
 if (vn > 0) return (-1);
 else return (1);


}








void Burgers::RHS (double *field, const mPoint &position, double time, double *rhs) const


{


  rhs[0]=0;    


}





bool Burgers::isPhysical ( double *Q ) const


{


    return true;


}





double Burgers::jumpQuantity ( double *Q ) const


{


  return Q[0];


}


/********************************************************************************/

/******Initial and Boundary Conditions*****/

/********************************************************************************/


class Exact : public FieldEvaluator

{

public:

  void eval(const mPoint &p, double time, double *val) const

  {

    double alpha = 0.5;

  double beta = 0.5;

  if (p(1) < 0.0001) val[0] = alpha + beta * sin(mPi * p(0));

else {

  double z;

  double x0, t0;

  double a, b;

  double half_int, f_z, f_a;

  int i;

  

  t0 = beta * p(1);

  x0 = p(0) - alpha * p(1);

  

  if (x0 < 0) x0 = -x0;

  

  a = -1.5;

  b =  1.5;

  

  i = 0;

  while (1)

    {half_int = (b - a) / 2.0;

    z = a + half_int;

    f_z = z - sin(mPi * (x0 - z * t0));

    f_a = a - sin(mPi * (x0 - a * t0));

    

    if ((f_z == 0.0) || (half_int < TOL))

      break;

    else if ((f_z * f_a) > 0)

      a = z;

    else

      b = z;

    }

  

  if ((p(0) - alpha * p(1)) < 0) z = -z;

  

  val[0] = alpha + beta * z;

}

}

  int order() const

  {

    return 4;

  }

};







/*

class Exact : public FieldEvaluator

{

public:

  void eval(const mPoint &p, double time, double *val) const

  {

     if (p(0)<=0) val[0] = 0;

    else val[0]=1;

    }

  int order() const

  {

    return 1;

  }

};

 */

class Init : public FieldEvaluator

{

public:

  void eval(const mPoint &p, double time, double *val) const

  {

	  if (time <0.01) val[0]=0; else

	  if (p(0) >= 0 ) val[0] = 1;

	   else val[0] = -1;

	   }

  int order() const

  {

    return 1;

  }

};

/*



class Exact : public FieldEvaluator

{

public:

  void eval(const mPoint &p, double time, double *val) const

  {

     if (p(0)<=p(1)) val[0] = 1;

	 else if (p(1)<=p(0) & p(0) <1 ) val[0]=(p(0)-1)/(p(1)-1);

    else val[0]=0;

      }

  int order() const

  {

    return 1;

  }

};

*/







void Burger( mDGMesh *theMesh, int order) 

{

  Exact e;

  Init in;

  DGBarthLimiter myLimiter;

  Burgers theLaw(&in,&in);

  DGLineSensor sensor1  ("t1.dat",mPoint(0,.1,0),mPoint(.3,0,0),100,2.);

  DGLineSensor sensor2  ("t2.dat",mPoint(0,.0,0),mPoint(.3,0,0),100,2.);



  //DGAnalysis analysis(theMesh,&theLaw,0,order);

  DGAnalysis analysis(theMesh,&theLaw,&myLimiter,order);

  DGGmshSensor sensor0  ("bur",0.1);

  analysis.addSensor(&sensor0);

  analysis.addSensor(&sensor1);

  analysis.addSensor(&sensor2);

  analysis.run();

}





































