#include "MaxwellLaw.h"
#include "mVector.h"
#include "mPoint.h"
#include "FieldEvaluator.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>

Maxwell::Maxwell (set<int> &w, FieldEvaluator *f)
  : Walls(w)
{
  farField = f;
  exactField = f;
}

Maxwell2d::Maxwell2d (set<int> &Walls, FieldEvaluator *f)
  : Maxwell(Walls,f)
{}

void Maxwell2d::Fi ( const mPoint &position , double *Q, mVector *flux) const
{
	//Q = [Hx,Hy,Ez]


	flux[0](0) = 0.0;
	flux[0](1) = Q[2]; //dHx/dt = -dEz/dy  *** signs flipped because we want flux on LHS
	flux[0](2) = 0.0;

	flux[1](0) = -1*Q[2];
	flux[1](1) = 0.0; //dHy/dt = dEz/dx
	flux[1](2) = 0.0;
	
	flux[2](0) = -1*Q[1];
	flux[2](1) = Q[0]; //dEz/dt = dHy/dx - dHx/dy
	flux[2](2) = 0.0;
}

void Maxwell2d::boundary (mVector &n, mVector &N, int BoundaryId, const mPoint 
			    &position,
			    double *Q, double *flux, double T) const
{
  //  -- reflective curved surfaces (20,000) only (DOMAIN boundary only)

	double nx = N(0);
	double ny = N(1);
	double Hdotn = Q[0]*nx+Q[1]*ny;

	double f[3];
	f[0] = -Q[0] + 2*nx*Hdotn;
	f[1] = -Q[1] + 2*ny*Hdotn;
	f[2] = Q[2];
	riemannSolver(n,position,Q,f,flux);

	/*double Z[3];
	Z[0] = 10*f[0];
	Z[1] = 10*f[1];
	Z[2] = 10*f[2];
	riemannSolver(n,position,Q,Z,flux);*/
}

void Maxwell2d::boundarySource (mVector &n, mVector &N, int BoundaryId, const mPoint 
			    &position,
			    double *Q, double *flux, double T) const
{
  //  -- 50000 Conditions = "fake" BC condition by Dale

	double nx = n(0);      //N normal isn't working here
	double ny = n(1);
	double cosValue = cos(3.14159265*T/0.04);

	double f[3];
	f[0] = 0;
	f[1] = Q[1];
	f[2] = Q[2];

	if (1)
	{

		// Aimed opposite the (outward) normal
		/*
		f[0] = -ny*cosValue;
		f[1] = nx*cosValue;
		f[2] = cosValue;
		*/

		
		// Simulated reflection from incoming plane wave (velocity in positive x) in perfect sync hitting boundary
		
		double fakeQ[3];
		fakeQ[0] = 0;
		fakeQ[1] = -cosValue;
		fakeQ[2] = cosValue;

		double Hdotn = fakeQ[0]*nx+fakeQ[1]*ny;

		f[0] = 0;//-fakeQ[0] + 2*nx*Hdotn;
		f[1] = 100;//-fakeQ[1] + 2*ny*Hdotn;
		f[2] = 1;//fakeQ[2];
		

	}

	riemannSolver(n,position,Q,f,flux);
}



void Maxwell2d::boundaryFlux (mVector &n, int BoundaryId, const mPoint 
			    &position,
			    double *Q, double *flux, double T) const
{
	// This is used only on DOMAIN BOUNDARY

	if (BoundaryId==40000) // Reflection on flat surface
	{
	double nx = n(0);
	double ny = n(1);
	double Hdotn = Q[0]*n(0)+Q[1]*n(1);

	double f[3];
	f[0] = -Q[0] + 2*nx*Hdotn;
	f[1] = -Q[1] + 2*ny*Hdotn;
	f[2] = Q[2];
    riemannSolver(n,position,Q,f,flux);
	}
	else  //Default Case --> use farfield data on boundary
	{
		double f[3];
		farField->eval(position,T,f);
		riemannSolver(n,position,Q,f,flux);
	} 
 
  
}

void Maxwell2d::riemannSolver( mVector &n ,
				const mPoint &position,
				double *u0,
				double *u1,
				double *flux) const
{

	double nx = n(0);
	double ny = n(1);
	double root2 = sqrt(2.0);
	
	/* --- This component has been condensed below */
	/*
	double wL[3],wR[3],w[3],uR[3];

	wL[0] = -ny*u0[0]/root2 + nx*u0[1]/root2 + u0[2]/root2;
	wL[1] = ny*u0[0] + nx*u0[1];
	wL[2] = ny*u0[0]/root2 - nx*u0[1]/root2 + u0[2]/root2;

	wR[0] = -ny*u1[0]/root2 + nx*u1[1]/root2 + u1[2]/root2;
	wR[1] = ny*u1[0] + nx*u1[1];
	wR[2] = ny*u1[0]/root2 - nx*u1[1]/root2 + u1[2]/root2;

	// w = wRiemann

	w[0] = wR[0];
	w[1] = 999990.0; //w[1] = (wR[1]+wL[1])/2.0;
	w[2] = wL[2];

	// uR = uRiemann

	uR[0] = -ny*w[0]/root2 + nx*w[1] + ny*w[2]/root2;
	uR[1] = nx*w[0]/root2 + ny*w[1] - nx*w[2]/root2;
	uR[2] = w[0]/root2 + 0.0 + w[2]/root2;
	*/
	
	
	//This is a simpler way of expressing what is done above
	double uR2[3];
	uR2[0] = (u0[0]+u1[0])*(nx*ny/2 + ny*ny/2) + (u0[1]+u1[1])*(nx*nx/2-nx*ny/2) + (u0[2]-u1[2])*ny/2;
	uR2[1] = (u0[0]+u1[0])*(-nx*ny/2 + ny*ny/2) + (u0[1]+u1[1])*(nx*nx/2 + nx*ny/2) + (u0[2]-u1[2])*(-nx)/2;
	uR2[2] = ny*(u0[0]-u1[0])/2 + nx*(u1[1]-u0[1])/2 + (u0[2]+u1[2])/2;
	

	/*
	// Central Difference Flux - I wouldn't use this if I was you! =)
	double uR[3];
	uR[0] = (u0[0]+u1[0])/2;
	uR[1] = (u0[1]+u1[1])/2;
	uR[2] = (u0[2]+u1[2])/2;
	*/
	
  mVector f[3];
  Fi(position,uR2,f);
  flux[0] = f[0] * n;
  flux[1] = f[1] * n;
  flux[2] = f[2] * n;

}

double Maxwell2d::maximumEigenValue(const double *theMean, const mPoint &point) const
{
  return 1.;  //used to compute time step things
}

double Maxwell2d::maximumEigenValue(mVector &n, const double *theMean, const double* theMean1) const
{
return 1.0;  //not used
}

double Maxwell2d::computeCFL(const int order,const double *theMean, const mPoint &point) const 
{
  return 1.;  // not called
}

double Maxwell2d::maxEigenValueCube(double du,double dv,double dw,double *theMean, const mPoint &point) const
{
return 1.0;  //not used
}

 void Maxwell2d::computeJacobian(mVector& n, const double* Q,double** Jac) const 
{
  double nx, ny;
  nx=n(0); ny=n(1); 
  
  // Define the matrix
  Jac[0][0]=0.0;
  Jac[0][1]=0.0;
  Jac[0][2]=ny;  // not called!!!
  Jac[1][0]=0.0;
  Jac[1][1]=0.0;
  Jac[1][2]=-nx;
  Jac[2][0]=ny;
  Jac[2][1]=-nx;
  Jac[2][2]=0.0;
}

void Maxwell2d::compute_left_eigenvectors(double *U, mVector &n, double *lev) const
{
lev[0] = n(1);
lev[1] = -n(0);  //these don't even get called?
lev[2] = 1.0;
}

void Maxwell2d::compute_right_eigenvectors(double *U, mVector &n, double *rev) const
{
rev[0] = -n(1);
rev[1] = n(0);  //these don't even get called?
rev[2] = 1.0;
}

int Maxwell2d::getEdgeOrientation(mVector &n,const mPoint &position,double *Q) 
const
{
return (1);  //not used
}
/********** P O S T   P R O **********/

#include <stdio.h>
int Maxwell::getNbFieldsOfInterest() const
{
  return 5;
}

void Maxwell::getNameAndSize(int i, int &n, char *name) const
{
  switch(i)
    {
    case 0:
      n = 1;
      strcpy(name,"Intensity_RefCond");
      break;
    case 1:
      n = 1;
      strcpy(name,"HxField_RefCond");
      break;
    case 2:
      n = 1;
      strcpy(name,"HyField_RefCond");
      break;
    case 3:
      n = 1;
      strcpy(name,"EField_RefCond");
      break;
    case 5:
      n = 1;
      strcpy(name,"HxSquared");
      break;
    case 4:
      n = 1;
      strcpy(name,"order");
      break;
    case 6:
      n = 1;
      strcpy(name,"limit");
      break;

    }
}

void Maxwell2d::getIthFieldOfInterest(int i, const mPoint &p, double *field, 
double *res, double T) const
{

double gamma=1.4;
  switch(i)
    {
    case 0:
      {
	res[0] = field[2]*field[2];
	break;
      }
    case 1:
      { 
	res[0] = field[0];
	break;
      }
    case 2:
      {
    res[0] = field[1];
	break;
	  }
    case 3:
      res[0] = field[2];
      break;
    case 5:
      res[0] = field[0]*field[0]+field[1]*field[1];
      break;
    case 4:
      {
	 res[0] = field[0];
      break;
      }
      break;
    case 6:
      res[0] = field[0];
      break;
    }
}


void Maxwell2d::RHS(double *field, const mPoint &position, double time,double *rhs) const
{
	double x = position(0);
	double y = position(1);
	double freq = 44;
    rhs[0] = 0.0;
    rhs[1] = 0.0;          //if Hy = -Ez ==> pulse *should* travel directly to the right.
	rhs[2] = 0.0;

	/*
	//Oscillator Source
	// using a source width of 0.01, and a frequency to determine effective wavelength in propagation
	// since speed is fixed to 1, using cos(pi*t/pulse_width) should give a pulse wavelength of 2*pulse_width (or a pulse width of pulse_width)
	
	// NOTICE - the pulses are now sinusoids, instead of the square-wave-like behavior of tanh - tanh
	
	double pulse_width = 0.03; // 0.08, 0.06, 0.04, 0.02

	if(x>=-0.12 && x<=-0.1 && (time/pulse_width)<=(1/2))
	{
	  rhs[2] = 75*cos(3.14159265*time/pulse_width);
	  rhs[1] = -rhs[2];
	}	
	*/

	// reflection.msh
	double pulse_width = 0.06;
	if(x>=-0.1 && x<=-0. && time/pulse_width <= 1)
	{
	  rhs[2] = 75* sin(3.14159265*time/pulse_width);
	  //rhs[1] = -rhs[2];
	}
	
	
	/* this one works nicely but will not be truly planar 
	rhs[2] = 1000*(tanh((x+0.5)/0.01)-tanh((x+0.51)/0.01))*(tanh((y+0.5)/0.01)-tanh((y-0.5)/0.01));
	if (time>0.0015)
	{rhs[2]=0;}
	rhs[1] = -rhs[2];
	*/
	
	
	rhs[3] = 0.0;
	rhs[4] = 0.0;

}
 
bool Maxwell2d::isPhysical ( double *Q ) const
{
  return true;
}


