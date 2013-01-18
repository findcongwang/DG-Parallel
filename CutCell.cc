#include "DGCell.h"
#include "ConservationLaw.h"
#include "FunctionSpace.h"
#include "Integrator.h"
#include "Mapping.h"
#include "mEntity.h"
#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
#include "FieldEvaluator.h"
#include "mMesh.h"
#include "CutCell.h"
#include "Geometry.h"
#include <stdio.h>
#include <math.h>

extern void freeMatrix(double **);

bool isInsideTheDomain(double x,double y,double z)
{ 
//if (x*x+y*y<0.2499999999) return 0;
if (x*x+y*y<1. ||x*x+y*y>1.384*1.384) return 0;
else return 1; 
}

cutCellGaussPoint::cutCellGaussPoint (const cutCellGaussPoint &other)
{
  isInsideDomain = other.isInsideDomain;
  fctsAtReflPoint = other.fctsAtReflPoint;
  reflCell = other.reflCell;
  N = other.N;
}

CutCell::CutCell (ConservationLaw*l, 
		mEntity*e, 
		FunctionSpace*f,
		FunctionSpace *er,
		Mapping *m,
		GaussIntegrator *I,Geometry* g)
  : DGCell(l,e,f,er,m,I)
{
	theGeometry = g;
//	init();
}

CutCell::~CutCell ()
{
  delete theFieldsCoefficients;
  delete [] theMean;
  delete [] limSlope;
  delete [] theRightHandSide;
  delete theFunctionSpace;
  delete theMapping;
  freeMatrix(theInvertMassMatrix);
}

void CutCell::init ()
{
  int i;
  
  cutCellGaussPoints.reserve(nbPtGauss);
  //Computing the integration points, weights, Jacobian.
  //x,y,z - physical coordinates
  double u,v,w,xr,yr,zr;
  for(i=0;i<nbPtGauss;++i)
    {
      cutCellGaussPoint pg;
	  volumeGaussPoint& PG = pt(i);
      //theGaussIntegrator->iPoint(theMeshEntity,i,order,u,v,w,weight); 
      pg.isInsideDomain = isInsideTheDomain(PG.x,PG.y,PG.z);
		
	  if (!pg.isInsideDomain) {// find the corresponding point inside the domain and what element it belongs to 
		  pg.reflCell = theGeometry->findReflPoint(this, PG.x,PG.y,PG.z,xr,yr,zr);
          pg.reflCell->getMapping()->invert(xr,yr,zr,u,v,w);
		  pg.fctsAtReflPoint = new double [fSize];
      pg.reflCell->theFunctionSpace->fcts(u,v,w,pg.fctsAtReflPoint);
	  }
      double r = sqrt(PG.x*PG.x+PG.y*PG.y);   //generalize!!!!!
      pg.N(0) = PG.x/r;
      pg.N(1) = PG.y/r;
      pg.N(2) = 0;
      cutCellGaussPoints.push_back(pg);
}
  L2ProjCutCell ();	
}


void CutCell::L2ProjCutCell ()
{
 double val[MaxNbEqn];
  int i, j, k;
  const int s = cSize * fSize;

  for( i=0;i<s;++i) theRightHandSide[i] = 0.0;
  
  for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint &PG=pt(i);
      cutCellGaussPoint  &pg = cutCellGaussPoints[i];

      if (!pg.isInsideDomain) interpolateRefl(pg.fctsAtReflPoint,pg.N,pg.reflCell, val);
      else interpolate(PG.fcts, val);
      
      for( j=0;j<cSize;++j)
	for( k=0;k<fSize;++k)
	  theRightHandSide[k+fSize*j] += val[j] * PG.fcts[k] * PG.JacTimesWeight;
	double tmp = (val[1]*val[1] + val[2]*val[2])/(val[0]*val[0]);
    }

  if (theFunctionSpace->isOrthogonal())
    {
      double inv_volume = 1./(2.*volume);
      for(i=0;i<cSize;++i)
	for(j=0;j<fSize;++j) 
	  theFieldsCoefficients->get(i,j) = theRightHandSide[j+fSize*i]*inv_volume;
    }
  else
    {
      for(int k = 0;k<cSize;++k)		
	for(int i = 0;i<fSize;++i)
	  {
	    double dxi = 0.0;
	    for(int j = 0;j<fSize;++j){
	      dxi += theInvertMassMatrix[i][j] * theRightHandSide[j+fSize*k];
	    }
	    theFieldsCoefficients->get(k,i) = dxi;
	  }
    }
 for( i=0;i<s;++i) theRightHandSide[i] = 0.0; 

for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint &PG=pt(i);
     
      interpolate(PG.fcts, val); 
     double vel = (val[1]*val[1] + val[2]*val[2])/(val[0]*val[0]);
	 double rho = val[0];
	 double p = 0.4*(val[3]-vel*.5*rho);
	// printf("vel = %17.16e, rho = %17.16e, p = %17.16se\n", vel, rho,p);
    }

}


