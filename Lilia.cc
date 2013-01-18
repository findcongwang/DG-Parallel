#include "DGAnalysis.h"
#include "FieldEvaluator.h"
#include "mPoint.h"
#include "ScalarLaw.h"
#include "LinearSystem.h"
#include "DGLimiter.h"
#include "DGSensor.h"
#include "mDGMesh.h"
#include <stdio.h>
#include <math.h>
#include "mVector.h"
#include "Constants.h"

class ExactOne : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
	 //if (space(1)-time*.402<-.625*(space(0)-.666*time)-0.8&& space(0)<2.) val[0]=1.; else val[0]=0; 
	//  val[0] = pow(sin(mPi*(2.*space(0)+1.*space(1)-3.*time)),1);
    val[0] = 1;
	  //val[0]=pow(space(0)-2.*space(1),1);
	// val[0]=pow(space(0)+space(1),1);          //linear
	 //val[0]=pow(space(0)+space(1),2);		  //quadratic
	 // val[0]=pow(space(0)+space(1),3);	      //cubic
	  //val[0] = 1;
	  //val[0]=space(0)+time;
	  val[0] = sin(space(0)+space(1));
     //val[0] = pow(2.*space(0)-3.*space(1)+space(2),4);
     //val[0] = pow(2.*space(0)-3.*space(1)+space(2),2);
  }
  int order() const
  {
    return 2;
  }
};

// Initial Conditiong

class CaseOne : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
     // if (space(0)>0.1 && space(0)<0.6 && space(1)>-0.25 && space(1) <0.25) val[0] = 1.0;
	  //else 
	  val[0]=0.0;
/*	 else   
	 {
	 double r=sqrt(pow(space(0)+0.45,2) + space(1)*space(1));
	 if (r<=0.35) val[0] = 1-r/0.35; 
	 else val[0]=0.0;
	 }*/
    
  double r=sqrt(pow(space(0)+0.5,2) + space(1)*space(1));
    if (r<=0.25) val[0] = pow(cos(2*mPi*r),2); 
   // else val[0]=0.0;
    //val[0] = pow(space(0),2)*pow(space(1),3)+pow(space(0),3)*pow(space(1),2)*3. +space(1)*space(1)+space(0)*space(0)+space(1)+space(0);
  }
  int order() const
  {
    return 3;
  }
};

class VelocityOne : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
   val[0] =1.;
    val[1] =-1.;

    val[0] = 2.*mPi*space(1);
    val[1] = -2.*mPi*space(0);
    val[2] = 0.0;
  }
  int order() const
  {
    return 1;
  }
};

void ScalarPb (mDGMesh *theMesh, int order)
{

  CaseOne i;
  VelocityOne v;
  ExactOne e;

  DGLineSensor sensor1  ("t1.dat",mPoint(-1,0),mPoint(1,0),256,1.0);
  //DGLineSensor sensor2  ("t2.dat",mPoint(-1,-1),mPoint(1,1),256,.1);
  //DGPointSensor sensor2 ("t2.dat",mPoint(0.0,0.0),1.0);
  
  //DGVertexLimiter myLimiter;
  //DGBarthLimiter myLimiter;
  //DGSuperBee myLimiter;
   DGMomentLimiter myLimiter;
  ScalarLaw theLaw (&i,&v,&i);
  
  //DGAnalysis analysis(theMesh,&theLaw,&myLimiter,order);
   DGAnalysis analysis(theMesh,&theLaw,0,order);
   analysis.addSensor(&sensor1);
  // analysis.addSensor(&sensor2);
  DGGmshSensor sensor0  ("squar",0.1);
  analysis.addSensor(&sensor0);
   //analysis.addSensor(&sensor2);
  analysis.run();
}

/********************Linear System****************************/

class Matrix1 : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    val[0] =8.;
    val[1] = -9.;
	val[2] =0.;
    val[3] =6.;
    val[4] = -13.;
	val[5] = 0.;
	val[6] = 0.;
	val[7] = 0.;
	val[8] = 0.;
	  double u=-3.;
	  double rho=1;
	  double K=1.;
/*	 val[0] =u;
    val[1] = K;
	val[2] = 0;
    val[3] = 1./rho;
    val[4] = u;
	val[5] = 0.;
	val[6] = 0.;
	val[7] = 0.;
	val[8] = u;	*/
  }
  int order() const { return 1;}
};

class Matrix2 : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    val[0] =13.;
    val[1] = -9.;
	val[2] = 0;
    val[3] = 6.;
    val[4] = -8.;
	val[5] = 0.;
	val[6] = 0.;
	val[7] = 0.;
	val[8] = 0.;
	  double v=-2.;
	  double rho=1.;
	  double K=1.;
	/* val[0] =v;
    val[1] = 0;
	val[2] = K;
    val[3] = 0;
    val[4] = v;
	val[5] = 0.;
	val[6] = 1./rho;
	val[7] = 0.;
	val[8] = v;*/
  }
  int order() const { return 1;}
};

class ExactSolution : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    val[0]=sin(0.5*(space(0)+space(1)));
    val[1]=cos(0.5*(space(0)+space(1)));
	val[2]=0.0; 
	/*val[0]=sin(0.5*(space(0)+space(1)));
    val[1]=space(0)*space(0);
	val[2]=space(1)*space(1);*/
	//val[0]=val[1]=val[2]=1.;
  }
  int order() const {return 2;}
};

void LinSystem (mDGMesh *theMesh, int order)
{
  Matrix1 a1;
  Matrix2 a2;
  ExactSolution e;

  DGLineSensor sensor1  ("t1.dat",mPoint(0,1),mPoint(1,0),3000,5.);
  DGPointSensor sensor2 ("t2.dat",mPoint(0.0,0.0),1.0);
  
  // DGVertexLimiter myLimiter;
  // DGBarthLimiter myLimiter;
  //DGSuperBee myLimiter;

  LinearSystem theLaw (&e,&e,&a1,&a2);
  
   //DGAnalysis analysis(theMesh,&theLaw,&myLimiter,order);
     DGAnalysis analysis(theMesh,&theLaw,0,order);
  DGGmshSensor sensor0  ("squar",5.0);
  analysis.addSensor(&sensor0);
  analysis.addSensor(&sensor1);
  
  analysis.run();
}

