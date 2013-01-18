#include "DGAnalysis.h"
#include "FieldEvaluator.h"
#include "mPoint.h"
#include "EulerLaw.h"
#include "MaxwellLaw.h"
#include "ScalarLaw.h"
#include "DGLimiter.h"
#include "DGSensor.h"                          
#include "mDGMesh.h"

#include "Constants.h"
#include <stdio.h>
#include <math.h>

// PROPAGATION OF A CONTACT DISCONTINUITY

class ShockTube1 : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    double rho,pres,velx,vely;
    if(space(0) < 0.5)
      {
	rho = 2.0;
	pres = 1.0;
	velx = 1.0;
	vely = 0.0;
      }
    else
      {
	rho = 1.0;
	pres = 1.0;
	velx = 1.0;
	vely = 0.0;
      }
    val[0] = rho;
    val[1] = rho*velx;
    val[2] = rho*vely;
    val[3] = 0.5 * rho * (velx*velx+vely*vely) + pres / (0.4);
  }
  int order() const
  {
    return 0;
  }
};

// THREE DISCONTINUITIES

class ShockTube2 : public FieldEvaluator //SOD porblem
{
  int dim;
public :
  ShockTube2(int nbVolumes)
  {
    if(nbVolumes) dim = 3;
    else dim = 2;
  }
  void eval (const mPoint &space, double time, double *val) const
  {
    double rho,pres,velx,vely,velz;
    if(space(0) < 0.5)
      {
	rho = 1.0;
	pres = 1.0;
	velx = 0.0;
	vely = 0.0;
	velz = 0.0;
      }
    else
      {
	rho = 0.125;
	pres = .1;
	velx = 0.0;
	vely = 0.0;
	velz = 0.0;
      }
    val[0] = rho;
    val[1] = rho*velx;
    val[2] = rho*vely;
    val[3] = 0.5 * rho * (velx*velx+vely*vely+velz*velz) + pres / (0.4);
    if(dim == 3)val[4] = rho*velz;
  }
  int order() const
  {
    return 0;
  }
};

class ShockTube3 : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    double rho,pres,velx,vely;
    //rho = sin(mPi*(1.*space(0)-2.*space(1)-time))+2.;
	rho = pow((2.*space(1)-1.*space(0)),2)+20.;
	//rho = space(0)*space(0)+2.*space(1)+2.*space(0)+space(1)*space(1);
	//rho=1.4;
	pres = 1.0;
	velx = 2.0;
	vely = 1.0;
    
	val[0] = rho;
    val[1] = rho*velx;
    val[2] = rho*vely;
    val[3] = 0.5 * rho * (velx*velx+vely*vely) + pres / (0.4);
  }
  int order() const
  {
    return 0;
  }
};

class ShockTube4 : public FieldEvaluator
{
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    double rho,pres,velx,vely;
    if(space(1) > 0.999)
      {
	//rho = space(0);
	rho = 2.4739;
	pres = 2.2685;
	velx = 2.5876;
	vely = -0.5438;
      }
    else
      {rho=1.4;
	//rho = space(0);
	pres = 1.0;
	velx = 2.9;
	vely = 0.0;
      }
    val[0] = rho;
    val[1] = rho*velx;
    val[2] = rho*vely;
    val[3] = 0.5 * rho * (velx*velx+vely*vely) + pres / (0.4);
  }
  int order() const
  {
    return 0;
  }
};

class ShockTube5 : public FieldEvaluator
{   ///SHU-OSHER problem
public :
  void eval (const mPoint &space, double time, double *val) const
  {
    double rho,pres,velx,vely;
    if(space(0) > 0.)
      {
	rho = 1 +0.2 * sin(5*space(0));
	pres = 1.;
	velx = -3.549648;
	vely = 0.;
      }
    else
      {rho=3.857143;
	pres = 10.333333;
	velx = -0.920279;
	vely = 0.0;
      }
	if(space(0) >= 10.) rho=1.0;
    val[0] = rho;
    val[1] = rho*velx;
    val[2] = rho*vely;
    val[3] = 0.5 * rho * (velx*velx+vely*vely) + pres / (0.4);
  }
  int order() const
  {
    return 0;
  }
};

void ShockTube (mDGMesh *theMesh, int order)
{
  set<int> walls;
  //walls.insert(200); // put the walls in this set
  //walls.insert(33); // put the walls in this set
  //walls.insert(55); // put the walls in this set
  //walls.insert(77); // put the walls in this set
  //ShockTube2 init(2);
   ShockTube3 init;
  // ShockTube2 init(theMesh->size(3));

  //DGLineSensor sensor1  ("tube_sensor2.dat",mPoint(-10,0.00),mPoint(10.,0.00),500,0.2);
  //DGPointSensor sensor2 ("tube_sensor3.dat",mPoint(0.3,0.005),1.0);

  //DGMomentLimiter myLimiter;
  //DGMomentLimiterEuler myLimiter;

  ConservationLaw *theLaw;
  if(theMesh->size(3) !=0)
    {
      theLaw = new Euler3d (walls,&init);
    }
  else
    {
	  theLaw = new Euler2dExactRiemann (walls,&init); 
      //theLaw = new Euler2d (walls,&init);
    }

  //DGAnalysis analysis(theMesh,theLaw,&myLimiter,order);
 DGAnalysis analysis(theMesh,theLaw,0,order);
  DGGmshSensor sensor0  ("tube",.002);
  analysis.addSensor(&sensor0);
  //analysis.addSensor(&sensor1);
  //analysis.addSensor(&sensor2);

  analysis.run();

  delete theLaw;
}

// Out flow
class OutFlow1 : public FieldEvaluator
{
  double machNumber,angle;
public :
  OutFlow1(double m, double a) : machNumber (m), angle(a) {}
  void eval (const mPoint &space, double time, double *val) const
  {
    double rho = 1.0;
    double pres = 1.0;
    double v = sqrt(1.4) * machNumber;  
    double anglerad = 3.1415927*angle/180.;
    val[0] = rho;
    val[1] = rho*cos(anglerad)*v ;
    val[2] = rho*sin(anglerad)*v;
    val[3] = 0.5 * rho * (v*v) + pres / (0.4);
  }
  int order() const
  {
    return 0;
  }
};

// Rayleigh Taylor
class RT1 : public FieldEvaluator
{
  double rho1,rho2,length,height;
  int dim;
public :
  RT1(int N , double r1 = 1.0, double r2 = 2.0, double l = 0.25, double h = 1.0) 
    : rho1 (r1), rho2(r2), length(l), height(h) {dim = (N)?3:2;}
  void eval (const mPoint &space, double time, double *val) const
  {
    double rho ;
    double pres ;
    double vx , vy, vz;  
    if(space(1) < 0.0)
      {
	rho = rho1;
	pres = 2.-space(1);
      }
    else
      {
	rho = rho2;
	pres = 2-2.*space(1);
      }
    
    double taper=6;                           /* Order for sin() taper */
    double machnum= (0.1);                 /* Mach number of perturbation */
    double GAMMA = 1.4;
    /* Min speed of sound is at top where p=1, rho=2 */
    double epsilon_z= machnum * sqrt(GAMMA*1.0/2.0);
    double epsilon_xy= -epsilon_z*taper/16.0;    
    /* Vertical */
    double x = space(0);
    double z = 0.5 + space(1);
    double y = space(2);

    if(dim == 3)
      {
	vy = 1 * epsilon_z *
	  cos(x/0.25 * 2*mPi) * cos(y/0.25 * 2*mPi) *
	  pow(sin(z*mPi),taper);
	
	/* Crosswise */
	
	vx = epsilon_xy *
	  sin(x/0.25 * 2*mPi) * cos(x/0.25 * 2*mPi) *
	  cos(z*mPi)*pow(sin(z*mPi),taper-1);
	
	vz = epsilon_xy *
	  cos(x/0.25 * 2*mPi) * sin(x/0.25 * 2*mPi) *
	  cos(z*mPi)*pow(sin(z*mPi),taper-1);
      }
    else
      {
	
	vy = 1 * epsilon_z *
	  cos(x/0.25 * 2*mPi) * cos(y/0.25 * 2*mPi) *
	  pow(sin(z*mPi),taper);
	
	/* Crosswise */
	vx = epsilon_xy *
	  sin(x/0.25 * 2*mPi) * cos(x/0.25 * 2*mPi) *
	  cos(z*mPi)*pow(sin(z*mPi),taper-1);
	vz = 0.0;
      }

    val[0] = rho;
    val[1] = rho*vx ;
    val[2] = rho*vy;
    val[3] = 0.5 * rho * (vx*vx+vy*vy+vz*vz) + pres / (0.4);
    if(dim == 3)val[4] = rho*vz;
  }
  int order() const
  {
    return 1;
  }
};





void OutFlow (mDGMesh *theMesh, double mach, double angle, int wall, int order)
{
  set<int> walls;
  walls.insert(wall); // put the walls in this set

  OutFlow1 init (mach,angle);

  //DGLineSensor sensor1  ("tube_sensor2.dat",mPoint(0,0.005),mPoint(1.,0.005),100,0.01);
  //DGPointSensor sensor2 ("tube_sensor3.dat",mPoint(0.3,0.005),1.0);

  DGBarthLimiter myLimiter;

  Euler2d theLaw (walls,&init);
  DGAnalysis analysis(theMesh,&theLaw,&myLimiter,order);
  //  analysis.addSensor(&sensor1);
  DGGmshSensor sensor0  ("tube",0.1);
  analysis.addSensor(&sensor0);

  analysis.run();
}



void RayleighTaylor (mDGMesh *theMesh, int wall, int order) 
{ 
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  RT1 init(theMesh->size(3)); 
 
  DGLineSensor sensor1  ("rt_sensor1.dat",mPoint(0.125,-0.5),mPoint(0.125,0.5),500,0.05); 
  //DGPointSensor sensor2 ("tube_sensor3.dat",mPoint(0.3,0.005),1.0); 
 
  //DGPhysicalLimiter myLimiter; 
	DGBarthLimiter myLimiter; 
  ConservationLaw *theLaw; 
  if(theMesh->size(3) !=0) 
    { 
      theLaw = new Euler3d (walls,&init,true); 
    } 
  else 
    { 
      theLaw = new Euler2d (walls,&init,true); 
    } 
  DGAnalysis analysis(theMesh,theLaw,&myLimiter,order); 
  //  analysis.addSensor(&sensor1); 
  DGGmshSensor sensor0  ("rt",0.05); 
  analysis.addSensor(&sensor0); 
  analysis.run(); 
  delete theLaw; 
} 
 
 
class TravellingShock : public FieldEvaluator 
{ 
protected : 
  int dim; 
  double M1,gamma,p1,rho1,anglerad; 
  double p2,rho2,v2,vshock; 
  double initXShockPosition; 
public : 
  TravellingShock(int d, double m, double ang, double g)  
    : dim(d), M1 (m),  
     gamma(g),p1(1.0),rho1(1.4),initXShockPosition(1./6.) 
  { 
    anglerad = mPi*ang/180.; 
    rho2 = rho1 * ( (gamma+1.) *M1*M1 ) / ((gamma-1)*M1*M1 +2.0); 
    p2 = p1 * (2.* gamma * M1*M1 - (gamma-1.0)) / (gamma+1.0); 
    double a1 =  sqrt (gamma*p1/rho1); 
    vshock = M1*a1; 
    v2 = vshock*(1. - rho1/rho2); 
  } 
 
void shockSide (const mPoint &space, double time, 
                 double &rho, double &p, double &v) const  
  { 
    double vshockx = vshock / sin (anglerad); 
    double posshockX = initXShockPosition + time * vshockx; 
 
    // y = a * x + b 
    // y = 0 => x = posshockX(t) 
    // dy/dx = tan (alpha) => a = tan (alpha) 
    // b = - tan (alpha) * posshockX(t) 
    // y = tan (alpha) * (x-posshockX(t))   
 
    double a =  tan (anglerad); 
    double shockEquation = space(1) - a * (space(0) - posshockX); 
//double shockEquation = space(0) - posshockX; 
     
    if(shockEquation < 0.0) 
      // post shock 
      { 
        v = 0.0; 
        p = p1; 
        rho = rho1; 
      } 
    else 
      { 
        v = v2; 
        p = p2; 
        rho = rho2; 
      } 
  } 
 
  virtual void eval (const mPoint &space, double time, double *val) const 
  { 
    double p,v,rho; 
 
    shockSide(space, time, rho, p, v); 
     
    val[0] = rho; 
    val[1] = rho*sin(anglerad)*v ; 
    val[2] = -rho*cos(anglerad)*v; 
    val[3] = 0.5 * rho * (v*v) + p / (gamma-1.0); 
    if(dim == 3)val[4] = val[2]; 
  } 
  int order() const 
  { 
    return 1; 
  } 
}; 


void DoubleMach (mDGMesh *theMesh, int wall, int order) 
{ 
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  TravellingShock init(2,10.0,60,1.4); 
 
  //DGPointSensor sensor2 ("tube_sensor3.dat",mPoint(0.3,0.005),1.0); 
 
  //DGBarthLimiterEuler myLimiter; 
  //DGBarthLimiter myLimiter; 
  //DGMomentLimiter myLimiter;
  DGMomentLimiterEuler myLimiter;
  ConservationLaw *theLaw; 
  theLaw = new Euler2dExactRiemann (walls,&init); 
  //theLaw = new  Euler2dHLLC(walls,&init);
  // theLaw = new Euler2d(walls,&init); 
    
   DGAnalysis analysis(theMesh,theLaw,&myLimiter,order); 
  //  analysis.addSensor(&sensor1); 
  DGGmshSensor sensor0  ("DM",0.1); 
  analysis.addSensor(&sensor0); 
  analysis.run(); 
  delete theLaw; 
} 

void TalkExample (mDGMesh *theMesh, int wall, int order) 
{ 
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  TravellingShock init(2,10.0,-45,1.4); 
 
  //DGPointSensor sensor2 ("tube_sensor3.dat",mPoint(0.3,0.005),1.0); 
 
  //DGBarthLimiterEuler myLimiter; 
	DGBarthLimiter myLimiter; 
  //DGMomentLimiter myLimiter;
  //DGMomentLimiterEuler myLimiter;
  ConservationLaw *theLaw; 
  theLaw = new Euler2dExactRiemann (walls,&init); 
  //theLaw = new  Euler2dHLLC(walls,&init);
  // theLaw = new Euler2d(walls,&init); 
    
   DGAnalysis analysis(theMesh,theLaw,&myLimiter,order); 
  //  analysis.addSensor(&sensor1); 
  DGGmshSensor sensor0  ("DM",0.01); 
  analysis.addSensor(&sensor0); 
  analysis.run(); 
  delete theLaw; 
} 

class SuperSonicVortex : public FieldEvaluator
{
public:
virtual void eval (const mPoint &space, double time, double *val) const 
  { 
    //Exact solution
    double p,vel,rho,theta; 
    const double r_in   = 1.0;         // radius inner
    const double r_out  = 1.384;       // and outer
    const double M_in   = 2.25;        //Mach inner 
    const double gamma  = 1.4;
    const double gm1_inv    = 1./(gamma -1.0);
    const double rho_in = 1.;
    const double sos_in = 1.0;           //on the inner circle, assuming gamma = 1.4;
  
    double r_inv = 1./sqrt(space(0) * space(0) + space(1)*space(1));
   
    theta = atan(space(1)/space(0));
    vel = sos_in*M_in*r_inv;
    rho = pow(1+1.0125*(1.-r_inv*r_inv),2.5);
    p = (1.0/gamma)*pow(rho,gamma);

    val[0] = rho;
    val[1] =  vel * sin(theta) * rho; 
    val[2] = -vel * cos(theta) * rho; 
    val[3] =  0.5 * rho * ( vel*vel ) + p * gm1_inv; 
     } 

  int order() const { return 0; } 
}; 

 void SuperVortex (mDGMesh *theMesh,int wall, int order) 
{ 
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  SuperSonicVortex init;
  
//   DGBarthLimiterEuler myLimiter; 
  //	DGVertexLimiter myLimiter; 
  ConservationLaw *theLaw; 
  theLaw	= new Euler2dExactRiemann (walls,&init); 
  //theLaw = new EulervanLeer (walls,&init); 
  //theLaw = new Euler2d (walls,&init); 
 DGBarthLimiter myLimiter;   
  //DGAnalysis analysis(theMesh,theLaw,&myLimiter,order); 
  DGAnalysis analysis(theMesh,theLaw,0,order); 
  //  analysis.addSensor(&sensor1); 
  DGGmshSensor sensor0  ("SV",5.); 
 // DGLineSensor sensor1  ("t1.dat",mPoint(1,0),mPoint(1.384,0.),3000,1.);
  analysis.addSensor(&sensor0); 
  //analysis.addSensor(&sensor1); 
  analysis.run(); 
  delete theLaw; 
} 

class FlowAroundCylinder : public FieldEvaluator
{
public:
virtual void eval (const mPoint &space, double time, double *val) const 
  { 
    //Inflow Boundary Conditions
    double Mach = 0.38;
   // Mach = 0.0;
    double gamma = 1.4;
    double p = 1.;
    double rho = gamma;
    double vel = Mach;
    double alpha = 2.0/180.*mPi;
    alpha=0.;
    //rho=1.4*exp(space(0));
    //p=exp(-space(1));
    //u=2.;v=1.;
    val[0] = rho;
    val[1] = rho*vel*cos(alpha); 
    val[2] = rho*vel*sin(alpha);
    //val[1] =rho*u;val[2] = rho*v;
    val[3] =  0.5 * rho * ( vel*vel ) + p / (gamma-1.0); 
    //val[3] =  0.5 * rho * ( u*u+v*v ) + p / (gamma-1.0);
  }
  int order() const { return 0; } 
}; 

void Cylinder (mDGMesh *theMesh,int wall, int order) 
{ 
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  FlowAroundCylinder init;
  DGBarthLimiter myLimiter;
  ConservationLaw *theLaw;
  theLaw  = new Euler2d (walls,&init); 
  theLaw  = new Euler2dExactRiemann (walls,&init); 
  //theLaw  = new Euler2dHLLC(walls,&init); 
// DGAnalysis analysis(theMesh,theLaw,&myLimiter,order);   
  DGAnalysis analysis(theMesh,theLaw,0,order); 
  DGGmshSensor sensor0  ("cyl",10.); 
    //DGLineSensor sensor1  ("t1.dat",mPoint(0,1),mPoint(1,0),3000,5.);
  analysis.addSensor(&sensor0); 
  //analysis.addSensor(&sensor1); 
  analysis.run(); 
  delete theLaw; 
} 

class FlowAroundSphere : public FieldEvaluator
{
public:
virtual void eval (const mPoint &space, double time, double *val) const 
  { 
    //Inflow Boundary Conditions
    double Mach = .38;
    double gamma = 1.4;
    double p = 1.;
    double rho = gamma;
    double vel = Mach;
    double alpha = 2.0/180.*mPi;
    alpha=0.;
    val[0] = rho;
    val[1] = rho*vel*cos(alpha); 
    val[2] = rho*vel*sin(alpha);
    val[4] = 0; ///?????
    val[3] =  0.5 * rho * ( vel*vel ) + p / (gamma-1.0); 
    //val[3] =  0.5 * rho * ( u*u+v*v ) + p / (gamma-1.0);
  }
  int order() const { return 0; } 
}; 

void Sphere (mDGMesh *theMesh,int wall, int order) 
{ 
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  FlowAroundSphere init;
  DGBarthLimiter myLimiter;
  ConservationLaw *theLaw;
  theLaw  = new Euler3d (walls,&init); 
  //theLaw  = new Euler2dExactRiemann (walls,&init); 
// DGAnalysis analysis(theMesh,theLaw,&myLimiter,order);   
  DGAnalysis analysis(theMesh,theLaw,0,order); 
  DGGmshSensor sensor0  ("sphere",1.); 
    //DGLineSensor sensor1  ("t1.dat",mPoint(0,1),mPoint(1,0),3000,5.);
  analysis.addSensor(&sensor0); 
  //analysis.addSensor(&sensor1); 
  analysis.run(); 
  delete theLaw; 
} 

double getAlpha() {return 1.25/180.*mPi;}
double getChord() {return 1.0;}

class NACA : public FieldEvaluator
{
public:
virtual void eval (const mPoint &space, double time, double *val) const 
  { 
    //Inflow Boundary Conditions
    double Mach = 0.72;
    double gamma = 1.4;
    double p = 1.;
    double rho = gamma;
    double vel = Mach;
    double u = Mach;
    double v=0.;
	double w=0.;
    double alpha = 1.250/180.*mPi;
	alpha=0.0;
    val[0] = rho;
    val[1] = rho*vel*cos(alpha); 
    val[2] = rho*vel*sin(alpha);
    val[3] =  0.5 * rho * ( vel*vel ) + p / (gamma-1.0); 
  }
  int order() const { return 0; } 
}; 

void FlowAroundAirfoil (mDGMesh *theMesh,int wall, int order) 
{ 
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  NACA init;
  DGBarthLimiter myLimiter;
  //VertLimiter myLimiter;
  ConservationLaw *theLaw;
  //theLaw  = new Euler2d (walls,&init); 
  theLaw  = new Euler3d (walls,&init); 
  //theLaw  = new Euler2dExactRiemann (walls,&init); 
  //DGAnalysis analysis(theMesh,theLaw,&myLimiter,order);   
  DGAnalysis analysis(theMesh,theLaw,0,order); 
  DGGmshSensor sensor0  ("cyl",1.); 
  //DGLineSensor sensor1  ("t1.dat",mPoint(0,1),mPoint(1,0),3000,5.);
  analysis.addSensor(&sensor0); 
  //analysis.addSensor(&sensor1); 
  analysis.run(); 
  delete theLaw; 
} 


class Ringleb : public FieldEvaluator
{
public:
virtual void eval (const mPoint &space, double time, double *val) const 
  { 
    //Inflow Boundary Conditions
    double gamma = 1.4;
    double vel = 0.0;
    double k=0.0;
    double c = sqrt(1.-0.5*(gamma-1.)*vel*vel);
    double p=1./1.4;
    double rho = pow(c,1./(gamma-1.));
    double J = 1./c + 1./(3*c*c*c)+1./(5*c*c*c*c*c)-0.5*log((1+c)/(1-c));
    double theta = 2*mPi-asin(vel/k);
    double x = (0.5/(k*k)-1./(vel*vel))/rho+J/2;
    double y = sqrt(1.-vel*vel/(k*k))/(k*rho*vel);
    double u =cos(theta);
    double v=vel*sin(theta);
    
    val[0] = rho;
    val[1] = rho*u; 
    val[2] = 0.0; 
    val[3] =  0.5 * rho * ( u*u ) + p / (gamma-1.0); 
    
  }
 int order() const { return 0; } 
};

class Forward_Facing_Step : public FieldEvaluator
{
public:
  virtual void eval (const mPoint &space, double time, double *val) const 
  { 
    //Inflow initial conditions
    const double gamma  = 1.4;
    const double rho= gamma;
    const double p = 1.;
    const double vel = 3.; //in x direction; v_y = 0
    
    val[0] = rho;
    val[1] = vel * rho; 
    val[2] = 0.0; 
    val[3] = 0.5*rho*vel*vel + p / (gamma-1.); 
  } 
  
  int order() const { return 0; } 
}; 

void ForwardFacingStep (mDGMesh *theMesh,int wall, int order) 
{ 
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  Forward_Facing_Step init;
  
  //   DGBarthLimiterEuler myLimiter; 
  //	DGVertexLimiter myLimiter;
  DGMomentLimiterEuler myLimiter;
  ConservationLaw *theLaw; 
  theLaw	= new Euler2dExactRiemann (walls,&init); 
  //theLaw = new EulervanLeer (walls,&init); 
  //theLaw = new Euler2dHLLC (walls,&init);
  //theLaw = new Euler2d (walls,&init); 
  //DGBarthLimiter myLimiter;   
  DGAnalysis analysis(theMesh,theLaw,&myLimiter,order); 
  //DGAnalysis analysis(theMesh,theLaw,0,order); 
  //  analysis.addSensor(&sensor1); 
  DGGmshSensor sensor0  ("FS",.5); 
  // DGLineSensor sensor1  ("t1.dat",mPoint(1,0),mPoint(1.384,0.),3000,1.);
  analysis.addSensor(&sensor0); 
  //analysis.addSensor(&sensor1); 
  analysis.run(); 
  delete theLaw; 
} 

class Spherical_Blast : public FieldEvaluator
{
public:
  virtual void eval (const mPoint &space, double time, double *val) const 
  { 
    //Inflow initial conditions
    const double gamma  = 1.4;
    const double rho= 1.0;
    double p;
    if (space(0)*space(0)+space(1)*space(1) <0.01) p = 100;
    else p = 1.;
    const double vel = 0.;
    
    val[0] = rho;
    val[1] = vel * rho; 
    val[2] = 0.0; 
    val[3] = 0.5*rho*vel*vel + p / (gamma-1.); 
  } 
  
  int order() const { return 0; } 
}; 

void SphericalBlast (mDGMesh *theMesh,int wall, int order) 
{ 
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  Spherical_Blast init;
  
  //   DGBarthLimiterEuler myLimiter; 
  //	DGVertexLimiter myLimiter;
  DGMomentLimiterEuler myLimiter;
  ConservationLaw *theLaw; 
  theLaw	= new Euler2dExactRiemann (walls,&init); 
  //theLaw = new EulervanLeer (walls,&init); 
  //theLaw = new Euler2dHLLC (walls,&init);
  //theLaw = new Euler2d (walls,&init); 
  //DGBarthLimiter myLimiter;   
 // DGAnalysis analysis(theMesh,theLaw,&myLimiter,order); 
  DGAnalysis analysis(theMesh,theLaw,0,order); 
  //  analysis.addSensor(&sensor1); 
  DGGmshSensor sensor0  ("Sphere",.05); 
  // DGLineSensor sensor1  ("t1.dat",mPoint(1,0),mPoint(1.384,0.),3000,1.);
  analysis.addSensor(&sensor0); 
  //analysis.addSensor(&sensor1); 
  analysis.run(); 
  delete theLaw; 
} 


class Vortex_Advection : public FieldEvaluator
{
public:
  virtual void eval (const mPoint &space, double time, double *val) const 
  { 
    //Inflow initial conditions
    const double gamma  = 1.4;
    const double rho= 1.0;
    double p=1.0;
	double x_c = 0.25;
    double y_c = 0.5;
	double r=sqrt((space(0)-x_c)*(space(0)-x_c)+(space(1)-y_c)+(space(1)-y_c));
	double alpha=0.204;
	double eps = 0.3;
	double vel=0.0;
    val[0] = rho;
    val[1] = vel * rho; 
    val[2] = 0.0; 
    val[3] = 0.5*rho*vel*vel + p / (gamma-1.); 
  } 
  
  int order() const { return 0; } 
}; 


void VortexAdvection (mDGMesh *theMesh,int wall, int order) 
{ 
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  Spherical_Blast init;
  
  //   DGBarthLimiterEuler myLimiter; 
  //	DGVertexLimiter myLimiter;
  DGMomentLimiterEuler myLimiter;
  ConservationLaw *theLaw; 
  theLaw	= new Euler2dExactRiemann (walls,&init); 
  //theLaw = new EulervanLeer (walls,&init); 
  //theLaw = new Euler2dHLLC (walls,&init);
  //theLaw = new Euler2d (walls,&init); 
  //DGBarthLimiter myLimiter;   
 // DGAnalysis analysis(theMesh,theLaw,&myLimiter,order); 
  DGAnalysis analysis(theMesh,theLaw,0,order); 
  //  analysis.addSensor(&sensor1); 
  DGGmshSensor sensor0  ("Sphere",.05); 
  // DGLineSensor sensor1  ("t1.dat",mPoint(1,0),mPoint(1.384,0.),3000,1.);
  analysis.addSensor(&sensor0); 
  //analysis.addSensor(&sensor1); 
  analysis.run(); 
  delete theLaw; 
} 

class Riemann_2D : public FieldEvaluator
{
public:
  virtual void eval (const mPoint &space, double time, double *val) const 
  { 
    //Inflow initial conditions
	double p,rho,u,v;
	  if (space(0)<0.5 && space(1)>0.5)
	  {
		  p = 0.35; rho = 0.5065; u = 0.8939; v = 0;
	  }
	  if (space(0)<0.5 && space(1)<0.5)
	  {
		  p = 1.1; rho = 1.1; u = 0.8939; v = 0.8939;
	  }
	  if (space(0)>0.5 && space(1)<0.5)
	  {
		  p = 0.35; rho = 0.5065; u = 0.; v = 0.8939;
	  }
	  if (space(0)>0.5 && space(1)>0.5)
	  {
		  p = 1.1; rho = 1.1; u = 0.; v = 0;
	  }
    const double gamma  = 1.4;
    val[0] = rho;
    val[1] = u * rho; 
    val[2] = v*rho; 
    val[3] = 0.5*rho*(u*u+v*v) + p / (gamma-1.); 
  } 
  
  int order() const { return 0; } 
}; 

void Riemann2D (mDGMesh *theMesh,int wall, int order) 
{ 
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  Riemann_2D init;
  
  //   DGBarthLimiterEuler myLimiter; 
  //	DGVertexLimiter myLimiter;
  DGMomentLimiterEuler myLimiter;
  ConservationLaw *theLaw; 
  theLaw	= new Euler2dExactRiemann (walls,&init); 
  //theLaw = new EulervanLeer (walls,&init); 
  //theLaw = new Euler2dHLLC (walls,&init);
  //theLaw = new Euler2d (walls,&init); 
  //DGBarthLimiter myLimiter;   
  DGAnalysis analysis(theMesh,theLaw,&myLimiter,order); 
 // DGAnalysis analysis(theMesh,theLaw,0,order); 
  //  analysis.addSensor(&sensor1); 
  DGGmshSensor sensor0  ("Riemann",.05); 
  // DGLineSensor sensor1  ("t1.dat",mPoint(1,0),mPoint(1.384,0.),3000,1.);
  analysis.addSensor(&sensor0); 
  //analysis.addSensor(&sensor1); 
  analysis.run(); 
  delete theLaw; 
} 

class DaleInit : public FieldEvaluator
{
public:
  virtual void eval (const mPoint &space, double time, double *val) const 
  { 
	  /* Maxwell INIT */
	  val[0] = 0.0;
	  val[1] = 0.0;
	  val[2] = 0.0; 
  }
 int order() const { return 0; } 
}; 

void DaleExample(mDGMesh *theMesh, int wall, int order)
{ 
  
  set<int> walls; 
  walls.insert(wall); // put the walls in this set 
  DaleInit init;
  
  //DGBarthLimiter myLimiter;  // ???
  ConservationLaw *theLaw;
  //theLaw = new Maxwell2d(walls,&init); 
  theLaw = new Euler2d(walls,&init);

  DGAnalysis analysis(theMesh,theLaw,0,order); 
  DGGmshSensor sensor0  ("Dale",0.015);   // Sensor output will be Dale-sample.pos or something like that, at 0.1, 0.2, 0.3... seconds
  //SummationSensor sensor1  ("Dale_1st_order_movie",0.1);   // Sensor output will be Dale-sample.pos or something like that, at 0.1, 0.2, 0.3... seconds
  analysis.addSensor(&sensor0);   //dwc** -- Sensor info
 // analysis.addSensor(&sensor1);

  analysis.run(); 
  delete theLaw; 
} 

