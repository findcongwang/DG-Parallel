#include "DGAnalysis.h"
#include "FieldEvaluator.h"
#include "mPoint.h"
#include "ScalarLaw.h"
#include "DGLimiter.h"
#include "DGSensor.h"
#include "mDGMesh.h"
#include <stdio.h>
#include <math.h>
#include "mVector.h"
#include "Constants.h"

class SuliExact : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    
	 //val[0] =  space(0) +  space(1) ;
    val[0]= sqrt(  space(1) - 2*space(0) + 1);    
  }
  int order() const
  {
    return 2;
  }
};

// Initial Conditiong

class SuliBC : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    if (time <= 0.0) val[0] = 0.0;
    else 
      { 

	if (space(0) < 0.01 && space(1) <1. && space(1) >0.5) val[0] = 0;
	if (space(0) < 0.01 && space(1) >0 && space(1) < 0.5) val[0] = 1.;

	if (space(1) < 0.01 && space(0) <1. && space(0) >0.5) val[0] = 0;
	if (space(1) < 0.01 && space(0) >0 && space(0) < 0.5) val[0] = 1.;

	if (space(0) > 0.99) val[0] = pow(sin(mPi*space(1)),2);
	      }
  }
  int order() const
  {
    return 5;
  }
};

class SuliVel : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    val[0] =10.* pow(space(1),2) -12.*space(0) +1.;
    val[1] = 1.+space(1);
	//val[0] = space(0);
    //val[1] =  -space(1);
    val[2] = 0.0;
  }
  int order() const
  {
    return 1;
  }
};

void Suli (mDGMesh *theMesh, int order)
{

	SuliBC init;
  SuliVel v;
  SuliExact e;

  DGLineSensor sensor1  ("t1.dat",mPoint(0,1),mPoint(1,1),100,.1);
  DGPointSensor sensor2 ("t2.dat",mPoint(0.0,0.0),1.0);

  DGBarthLimiter myLimiter;

  ScalarLaw theLaw (&init,&v,&e);
  //DGAnalysis analysis(theMesh,&theLaw,&myLimiter,order);
  DGAnalysis analysis(theMesh,&theLaw,0,order);
  // analysis.addSensor(&sensor1);
  //analysis.addSensor(&sensor2);
  DGGmshSensor sensor0  ("squar",0.1);
  analysis.addSensor(&sensor0);
  analysis.run();
}





