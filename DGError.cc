#include "DGCell.h"
#include "ConservationLaw.h"
#include "FunctionSpace.h"
#include "Integrator.h"
#include "Mapping.h"
#include "mEntity.h"
#include "gEntity.h"
#include "Edge.h"
#include "Face.h"
#include "mPoint.h"
#include "Vertex.h"
#include "FieldEvaluator.h"
#include "Constants.h"
#include "mMesh.h"
#include <stdio.h>
#include <math.h>

extern void invmat (double **, double **, int);
extern double ** allocateMatrix(int );
extern void freeMatrix( double **v);

void DGCell::allocateErrorMassMatrices ()
{
  freeMatrix(theErrorMassMatrix);
  int Size = 2 * theFunctionSpace->order() + 3 ;	
  //  printf("size = %d\n",fSize);

  theErrorMassMatrix = allocateMatrix(Size);
  freeMatrix(theErrorInvertMassMatrix);
  theErrorInvertMassMatrix = allocateMatrix(Size);
  //print("Massorder = %d\n",order);
  
  for(int j=0;j<Size;j++)
      for(int k=0;k<Size;k++)
	  {
		 theErrorInvertMassMatrix[j][k] = 0.0;
		 theErrorMassMatrix[j][k] = 0.0;
	  }
  for(int k=0;k<Size;k++)
  {
	  theErrorInvertMassMatrix[k][k] = 1.0;
	  theErrorMassMatrix[k][k] = 1.0;
  }
}

void DGCell::computeErrorMassMatrices ()
{
	int Size = 2 * theFunctionSpace->order() + 3;
  /* for(int i=0;i<Size ;i++)
	{
		for(int j=0;j<Size ;j++)
		printf(" %f ", theErrorMassMatrix[i][j] );

		printf("\n");
	}  
	printf("\n"); */ 
	invmat(theErrorMassMatrix,theErrorInvertMassMatrix,Size);
}

void DGCell::computeErrorVolumeContribution ()
{
/*
 // we compute int_{Cell} F(U +E) grad w dCell
  
  int order = 3 * theErrorFunctionSpace->order() + theMapping->order() - 1;
  vector<mVector> grads;
  vector<double> functs;
  int eSize = theErrorFunctionSpace->size();
  int eOrder = theErrorFunctionSpace->pSize();
  int fOrder = theFunctionSpace->pSize();
  int fSize = theFunctionSpace->size();
  int cSize = theConservationLaw->nbFields();
 
  grads.reserve(eSize);
  functs.reserve(eSize);
  double u,v,w,weight,val[256],rhs[256], err[256];
  mPoint p;

  theGaussIntegrator gauss(theMeshEntity);
  int i,j,k;
  for(i=0;i<gauss.nbIntegrationPoints(order);i++)
    {
      theGaussIntegrator->iPoint(i,order,u,v,w,weight);
      theErrorFunctionSpace->fcts(u,v,w,functs);
      double detJac = theErrorFunctionSpace->grads(u,v,w,theMapping,grads);
      interpolate( u,v,w, val );
      interpolateError( u,v,w, err );
      for(j=0;j<cSize; j++)
      val[j]+=err[j] ; 
      theMapping->eval(u,v,w,p(0),p(1),p(2));
      //printf("val,u,v,w,p(0),p(1),p(2)  %f %f %f %f %f %f %f\n",val[0], u,v,w,p(0),p(1),p(2));
      theConservationLaw->RHS(val,p,rhs);
      //printf("RHS in volume %d  %f\n", 2 , theErrorRightHandSide[2]);
      for(j=0;j<cSize;j++)
	{
	  mVector flux;
	  theConservationLaw->Fi (j,p,val,flux);
	  // printf("flux (%d) = %f %f %f\n",j,flux(0),flux(1),flux(2));
	  for( k=0;k<eOrder ;k++)

	      theErrorRightHandSide[k+(eOrder+fOrder)*j] += (flux * 
							     grads[fSize+k]  + rhs[j] *functs[fSize+k]) * detJac * weight;
	      
	
    
  for( k=0;k<eOrder ;k++)
    {
      //printf("volume , rhs,rhs[0]  %f %f %f %f \n",flux *grads[fSize+k]* detJac * weight,detJac * weight*rhs[0] *functs[fSize+k],rhs[0],functs[fSize+k]);
      //printf("grads %f %f \n",grads[fSize+k](0), grads[fSize+k](1));
    }
	}
    }
*/
}


void DGBoundaryCell::computeOrientation()
{ 
/*
  // if (!mright) return;
  DGCell *right = 0;
  DGCell *left = (DGCell*)mleft->getCell();
  if (mright) right = (DGCell*)mright->getCell();
  double u,v,w,weight,valleft[256],x,y,z;
  GaussIntegrator gauss(theBoundaryEntity);
  mVector n;
  int order= computeErrorOrder(left,right);  
   for(int i=0;i<gauss.nbIntegrationPoints(order);i++)
    {
      gauss.iPoint(i,order,u,v,w,weight);
      double detJac = theMapping->detJac(u,v,w);
      left  ->interpolate(uleft[i],vleft[i],wleft[i],valleft);
      //   right ->interpolate(urght[i],vrght[i],wrght[i],valrght);

      // compute the normalvector
     if(left->theMeshEntity->find(theBoundaryEntity))
	left->theMapping->normalVector(theBoundaryEntity,uleft[i],vleft[i],wleft[i],n); 
 

  theMapping->eval(u,v,w,x,y,z);
  int orient = left->theConservationLaw->getEdgeOrientation(n,mPoint(x,y,z),valleft);
 
  mEntity *f = left->theMeshEntity;
   for(int j=0;j<f->size(1);j++)
    {
      Edge *e = (Edge*) f->get(1,j);
      if (e == theBoundaryEntity) left->edgeOrientation[j] = orient;		        if (mright)
     {
       f = right->theMeshEntity;
       for(int j=0;j<f->size(1);j++)
	 {
	   Edge *e = (Edge*) f->get(1,j);
	   if (e == theBoundaryEntity) right->edgeOrientation[j] = !orient;

	 }
     }
    }
    }
*/
}

void DGBoundaryCell::computeBoundaryRestrictions(double T)
{
/*
  DGCell *right = 0;
  DGCell *left = (DGCell*)mleft->getCell();
  if(mright)right = (DGCell*)mright->getCell();
  
  if(!left)printf("arrhgh\n");
  if(mright && !right)printf("arrhgh\n");
  
   int order = computeErrorOrder(left,right);
   //printf("order %d", order);
  double u,v,w,weight,valleft[256],valrght[256],erleft[256],erright[256];
  double riemann[256];
  static vector<double>  fctleft;
  static vector<double>  fctrght;
  int cSize = left->theConservationLaw->nbFields();
  int eSizeLeft = left->theErrorFunctionSpace->size();
  int fSizeLeft = left->theFunctionSpace->size();
  int fOrderLeft = left->theFunctionSpace->pSize();
  int eOrderLeft = left->theErrorFunctionSpace->pSize();
  int fSizeRight = 0;
  int eSizeRight = 0;
  int eOrderRight =0;
  int fOrderRight =0;
  if(right)
    {
      fSizeRight = right->theFunctionSpace->size();
      eSizeRight = right->theErrorFunctionSpace->size();
      fOrderRight = right->theFunctionSpace->pSize();
      eOrderRight = right->theErrorFunctionSpace->pSize();
    }
  
  
  fctleft.reserve(eSizeLeft);  
  if(right)fctrght.reserve(eSizeRight);  
  mVector n;
  GaussIntegrator gauss(theBoundaryEntity);
  
  for(int i=0;i<gauss.nbIntegrationPoints(order);i++)
    {
      // get integration point
      gauss.iPoint(i,order,u,v,w,weight);
      
      // get det jacobian
      double detJac = theMapping->detJac(u,v,w);
      
      // compute the normalvector
      if(left->theMeshEntity->find(theBoundaryEntity))
	left->theMapping->normalVector(theBoundaryEntity,uleft_e[i],vleft_e[i],wleft_e[i],n);
      else
	{
	  if(!right)printf("aaargh\n");
	  right->theMapping->normalVector(theBoundaryEntity,urght_e[i],vrght_e[i],wrght_e[i],n);
	  n *= -1.0;
	}
      //printf("normal %f %f \n", n(0),n(1));
      // get shape functions on both sides
      left->theErrorFunctionSpace->fcts(uleft_e[i],vleft_e[i],wleft_e[i],fctleft);
      if(right)
	right->theErrorFunctionSpace->fcts(urght_e[i],vrght_e[i],wrght_e[i],fctrght);
      
      mPoint p;
      theMapping->eval(u,v,w,p(0),p(1),p(2));
      //       printf("u,v,w,p(0),p(1),p(2) %f %f %f %f %f %f\n", u,v,w,p(0),p(1),p(2));
      int j,k;
      left ->interpolate(uleft_e[i],vleft_e[i],wleft_e[i],valleft);
      int orient = left->theConservationLaw->getEdgeOrientation(n,p,valleft);
      double normalVelocity=0;
      
      if (orient == -1) 
	{
	  left->theConservationLaw->normalVelocity(n,p,T,valleft,normalVelocity);
	  //printf("normalVlocity %f \n", normalVelocity);
	  for(j=0;j<fOrderLeft;j++)
	    {
	      for(k=0;k<eOrderLeft;k++)
		left->theErrorMassMatrix [j+eOrderLeft][k] += detJac * 
		  normalVelocity  *  fctleft[eSizeLeft-eOrderLeft+k] * 
		  fctleft[eSizeLeft-eOrderLeft-fOrderLeft+j] * weight;
	      for(k=0;k<fOrderLeft;k++)
		left->theErrorMassMatrix [j+eOrderLeft][k+eOrderLeft] +=detJac * normalVelocity * 
		fctleft[eSizeLeft+k-fOrderLeft-eOrderLeft] * 
		fctleft[eSizeLeft-eOrderLeft-fOrderLeft+j] * weight;
	    }
	}
      
      else if (right)
	{
	  n *= -1;
	  right ->interpolate(uleft_e[i],vleft_e[i],wleft_e[i],valrght);
	  right->theConservationLaw->normalVelocity(n,p,T,valrght,normalVelocity);
	  for(j=0;j<fOrderRight;j++)
	    {
	      for(k=0;k<eOrderRight;k++)
		{
		  right->theErrorMassMatrix [j+eOrderRight][k] += detJac * normalVelocity * 
		    fctrght[eSizeRight-eOrderRight+k] * 
		    fctrght[eSizeRight-eOrderRight-fOrderRight+j] * weight;
		  
		}
	      for(k=0;k<fOrderRight;k++)
		right->theErrorMassMatrix [j+eOrderRight][k+eOrderRight] += detJac * normalVelocity * 
		  fctrght[eSizeRight+k-fOrderRight-eOrderRight] * 
		  fctrght[eSizeRight-eOrderRight-fOrderRight+j] * weight;
	    }
	}
    }
  //  for (int k=0; k<eOrderLeft+ fOrderLeft; k++)
  //  for (int j=0; j<eOrderLeft+ fOrderLeft;j++)
  //    printf("err mass %d %f \n", k,left->theErrorMassMatrix[k][j]);
*/
}


void DGBoundaryCell::computeErrorJump(double T)
{
/*
DGCell *right = 0;
DGCell *left = (DGCell*)mleft->getCell();
if(mright)right = (DGCell*)mright->getCell();

if(!left)printf("arrhgh\n");
if(mright && !right)printf("arrhgh\n");
int order = computeErrorOrder(left,right);
double u,v,w,weight,valleft[256],valrght[256],erleft[256],erright[256];
double riemann[256];
int cSize = left->theConservationLaw->nbFields();
 mVector n;
 double super,superR, exact[1],diff[1],exactNV,approxNV;
 super = 0;
 superR =0;
 mPoint p;
 int orient;

 GaussIntegrator gauss(theBoundaryEntity);
 
 for(int i=0;i<gauss.nbIntegrationPoints(order);i++)
   {
     // get integration point
     gauss.iPoint(i,order,u,v,w,weight);
     
     // get det jacobian
     double detJac = theMapping->detJac(u,v,w);
     
     // compute the normalvector
     if(left->theMeshEntity->find(theBoundaryEntity))
       left->theMapping->normalVector(theBoundaryEntity,uleft_e[i],vleft_e[i],wleft_e[i],n);
     else
       {
	 if(!right)printf("aaargh\n");
	 right->theMapping->normalVector(theBoundaryEntity,urght_e[i],vrght_e[i],wrght_e[i],n);
	 n *= -1.0;
       }
     // get solutions on both sides
     left ->interpolate(uleft_e[i],vleft_e[i],wleft_e[i],valleft);
     //printf("valleft init  %f \n",valleft[0]);
     if(right)
       right->interpolate(urght_e[i],vrght_e[i],wrght_e[i],valrght);
     //if(right) printf("valright %f \n",valrght[0]);

     //get error on both sides
     left ->interpolateError(uleft_e[i],vleft_e[i],wleft_e[i],erleft);
     if(right)
       right->interpolateError(urght_e[i],vrght_e[i],wrght_e[i],erright);
     
	 //mPoint p;
    theMapping->eval(u,v,w,p(0),p(1),p(2));
     int j,k;
	 orient = left->theConservationLaw->getEdgeOrientation(n,p,valleft);
      double normalVelocity=0;
     if (right && orient!=-1) 
		 right->jumpError+=erright[0]*detJac*weight;

if ((p(0) +p(1)<=1.001) && (p(0) +p(1)>=0.999)) 
{
	 left->theConservationLaw->exactField->eval(p,T,exact);
	 diff[0]= exact[0]-valleft[0];
	 //diff[0]= (exact[0] * exact[0] -valleft[0]*valleft[0]) *0.5 *n(0)+ 
		 (exact[0] - valleft[0]) * n(1) ;
	 //left->theConservationLaw->normalVelocity(n,p,T,exact,exactNV);
	 //left->theConservationLaw->normalVelocity(n,p,T,valleft,approxNV);
	 //if (orient==-1)
	 super+=diff[0]*detJac*weight;
		 if(right) 
		 {
			 diff[0]= exact[0] -valrght[0];
			 n*= -1;
			 //right->theConservationLaw->normalVelocity(n,p,T,exact,exactNV);
             superR+=diff[0]*detJac*weight;
		 }
 
} 
}
// if ((p(0) +p(1)<=1.001) && (p(0) +p(1)>=0.999) && orient ==-1) 
//	 printf("super %e %e p (%f, %f) \n", super,superR, p(0), p(1));
// if ((p(0) +p(1)<=1.001) && (p(0) +p(1)>=0.999) && orient !=-1) 
//	 printf("super %e %e p (%f, %f) \n", superR,super, p(0), p(1));
*/
 }

 double DGCell::getJumpError() const
 {
 return jumpError;
 }

void DGBoundaryCell::computeErrorBoundaryContributions(double T)
{
/*
DGCell *right = 0;
DGCell *left = (DGCell*)mleft->getCell();
if(mright)right = (DGCell*)mright->getCell();

if(!left)printf("arrhgh\n");
if(mright && !right)printf("arrhgh\n");
int order = computeErrorOrder(left,right);
//printf("order %d", order);
double u,v,w,weight,valleft[256],valrght[256],erleft[256],erright[256];
double riemann[256];
static vector<double>  fctleft;
static vector<double>  fctrght;
int cSize = left->theConservationLaw->nbFields();
int eSizeLeft = left->theErrorFunctionSpace->size();
int fSizeLeft = left->theFunctionSpace->size();
int fOrderLeft = left->theFunctionSpace->pSize();
int eOrderLeft = left->theErrorFunctionSpace->pSize();
int fSizeRight = 0;
int eSizeRight = 0;
int eOrderRight =0;
int fOrderRight =0;
if(right)
{
  fSizeRight = right->theFunctionSpace->size();
  eSizeRight = right->theErrorFunctionSpace->size();
  fOrderRight = right->theFunctionSpace->pSize();
  eOrderRight = right->theErrorFunctionSpace->pSize();
  
}
 
 
 fctleft.reserve(eSizeLeft);  
 if(right)fctrght.reserve(eSizeRight);  

 mVector n;
 
 GaussIntegrator gauss(theBoundaryEntity);
 
 for(int i=0;i<gauss.nbIntegrationPoints(order);i++)
   {
     // get integration point
     gauss.iPoint(i,order,u,v,w,weight);
     
     // get det jacobian
     double detJac = theMapping->detJac(u,v,w);
     
     // compute the normalvector
     if(left->theMeshEntity->find(theBoundaryEntity))
       left->theMapping->normalVector(theBoundaryEntity,uleft_e[i],vleft_e[i],wleft_e[i],n);
     else
       {
	 if(!right)printf("aaargh\n");
	 right->theMapping->normalVector(theBoundaryEntity,urght_e[i],vrght_e[i],wrght_e[i],n);
	 n *= -1.0;
       }
     //     printf("\n");
     //printf("normal vector %f %f  \n",n[0],n[1]);
     // get shape functions on both sides
     left->theErrorFunctionSpace->fcts(uleft_e[i],vleft_e[i],wleft_e[i],fctleft);
     if(right)
       right->theErrorFunctionSpace->fcts(urght_e[i],vrght_e[i],wrght_e[i],fctrght);
     
     // get solutions on both sides
     left ->interpolate(uleft_e[i],vleft_e[i],wleft_e[i],valleft);
     //printf("valleft init  %f \n",valleft[0]);
     if(right)
       right->interpolate(urght_e[i],vrght_e[i],wrght_e[i],valrght);
     //if(right) printf("valright %f \n",valrght[0]);

     //get error on both sides
     left ->interpolateError(uleft_e[i],vleft_e[i],wleft_e[i],erleft);
     //printf("erleft %f \n",erleft[0]);
     //printf("uleft[i],vleft[i],wleft[i] %f %f %f \n", uleft_e[i],vleft_e[i],wleft_e[i]);
     if(right)
       right->interpolateError(urght_e[i],vrght_e[i],wrght_e[i],erright);
     //if(right) printf("erright %f \n",erright[0]);
     mPoint p;
     theMapping->eval(u,v,w,p(0),p(1),p(2));
     int j,k;
     for(j=0;j<cSize; j++)
       {
	 valleft[j]+= erleft[j];
	 //printf("valleft %f \n",valleft[0]);
	 if (right)  valrght[j]+= erright[j];
	 //	 if (right) printf("valright %f \n",valrght[0]);
       }
     // the coupling : a riemann solver
     if(right)
	 {
		 computeErrorMaxIV(n, left->theErrorMean, right->theErrorMean);
		 double maxEIV = getErrorMaxIV();
       left->theConservationLaw->riemannSolver (n,p,valleft,valrght,riemann);
	 }
     // or boundary conditions
     else
       {
	 // compute parametric coordinates on both sides
	 left->theConservationLaw->boundaryFlux (n,theBoundaryEntity->getClassification()->getId(),
						 p,valleft,riemann,T);
       }
     // left and right contributions
     for(j=0;j<cSize;j++)
	 {
	 for(k=0;k<eOrderLeft;k++)
	   {
	     left->theErrorRightHandSide[k+(fOrderLeft+eOrderLeft)*j] -= 
	       (riemann[j]*fctleft[fSizeLeft + k]) * detJac * weight;
	     //     printf("riemann[j]*fctleft[fSizeLeft + k] * detJac * weight %f %f %f  \n", riemann[j],fctleft[fSizeLeft + k] ,riemann[j]*fctleft[fSizeLeft + k] * detJac * weight );  
	   }
	 	 if(right)
	   for(k=0;k<eOrderRight;k++)
	     right->theErrorRightHandSide[(eOrderRight+fOrderRight)*j+k] += 
	     (riemann[j]*fctrght[fSizeRight+k]) * detJac * weight;
       }
   }
*/
}
  

void DGCell::setType()
{
  if ( edgeOrientation[0] == 1 && edgeOrientation[1] == 0 && 
       edgeOrientation[2] ==0) {type = 1; return;}
  if ( edgeOrientation[0] == 0 && edgeOrientation[1] == 1 && 
       edgeOrientation[2] ==0) {type = 2; return;}
  if ( edgeOrientation[0] == 0 && edgeOrientation[1] == 0 && 
       edgeOrientation[2] ==0) {type = 3; return;}
  if ( edgeOrientation[1] == 1 && edgeOrientation[2] == 1)  
        {type = 4; return;}
  if ( edgeOrientation[0] == 1 && edgeOrientation[1] == 1 )  
        {type = 5; return;}     
  if ( edgeOrientation[0] == 1 && edgeOrientation[2] == 1 ) 
        {type = 6; return;}
}

void DGCell::interpolateError(double u, double v, double w, double *field)
{
int eSize = theErrorFunctionSpace->size();
int fSize = theFunctionSpace->size();
int fOrder = theFunctionSpace->pSize();
int eOrder = theErrorFunctionSpace->pSize();
double  *fct;
  fct = new double [eSize];  
  theErrorFunctionSpace->fcts(u,v,w,fct);
  int i, j ;	
  for(i=0;i<cSize;i++)
    {
      field[i] = 0;;
    }
  
  for(j=0; j<eOrder; j++)
    {
        for(i=0; i<cSize; i++)
	{
	  field[i] += (fct[j+fSize] * (theErrorCoefficients[i])[j]);
	}
    }

  for(j=0;j<fOrder;j++)
    {
      for(i=0;i<cSize;i++)
	{
	  field[i] += (fct[j+fSize-fOrder] * (theErrorCoefficients[i])[j+eOrder]); 
	}
    }
  delete [] fct;
}	

 int DGBoundaryCell::computeErrorOrder(DGCell *left, DGCell *right) const
  {
   if(right)
    return        3 * ((left->theErrorFunctionSpace->order()>right->theErrorFunctionSpace->order()) ?
 	   left->theErrorFunctionSpace->order()  :right->theErrorFunctionSpace->order()) + theMapping->detJacOrder() ;
    else
    return 3*(left->theErrorFunctionSpace->order())+ theMapping->detJacOrder() +4;
  }


double DGCell::L1exact(double time) 
{
  double errorL1 = 0;
  double exact[MaxNbEqn];
  double val[MaxNbEqn];
  double tmp=0.0;
  for(int i=0;i<nbPtGauss;i++)
    {
      volumeGaussPoint &pg = pt(i);
      interpolate(pg.fcts,val);
      theConservationLaw->getExactSolution()->eval(mPoint(pg.x,pg.y,pg.z),time,exact);
      errorL1+= fabs( exact[0] - val[0])  * pg.JacTimesWeight;
	  tmp+= pg.JacTimesWeight;
    }

  return errorL1;
}

double DGCell::LinfError(double time) 
{
  double error = 0;
  double exact[MaxNbEqn];
  double val[MaxNbEqn];
  for(int i=0;i<nbPtGauss;i++)
    {
      volumeGaussPoint &pg = pt(i);
      interpolate(pg.fcts,val);
      theConservationLaw->getExactSolution()->eval(mPoint(pg.x,pg.y,pg.z),time,exact);
      if (error<fabs( exact[0] - val[0])) error = fabs( exact[0] - val[0]);
    }
  return error;
}

double DGCell::L2exact()
{
  double gamma=1.4;
  double errorL2 = 0;
  double val[MaxNbEqn];
  double ent,p;
  double ent_init = -gamma*log(gamma);
  for(int i=0;i<nbPtGauss;i++)
    {
      volumeGaussPoint &pg = pt(i);
      interpolate(pg.fcts,val);
	  //printf("val %e %e %e %e \n",val[0],val[1],val[2],val[3]);
      p = (gamma-1.)*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
      ent = log(p/pow(val[0],gamma));
      errorL2+= (ent - ent_init)* (ent - ent_init)* pg.JacTimesWeight;
    }
  return errorL2;
}

double DGCell::L1exactPressure(double time)
{
  double gamma=1.4;
  double errorL1 = 0;
  double val[MaxNbEqn], val_ex[MaxNbEqn];
  double p,p_ex;
  for(int i=0;i<nbPtGauss;i++)
    {
      volumeGaussPoint &pg = pt(i);
      mPoint pos(pg.x,pg.y,pg.z);
      interpolate(pg.fcts,val);
      theConservationLaw->getExactSolution()->eval(pos,time,val_ex);
      p = (gamma-1.)*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
      p_ex = (gamma-1.)*(val_ex[3] - 0.5*(val_ex[1]*val_ex[1]+
					  val_ex[2]*val_ex[2])/val_ex[0]);
      errorL1+= fabs(p-p_ex)* pg.JacTimesWeight;
    }
  return errorL1;
}

double DGCell::entropyErrorL1_2D(double time)
{
  double gamma=1.4;
  double errorL1 = 0;
  double val[MaxNbEqn],val_ex[MaxNbEqn];
  double ent_error,p,p_ex;
  for(int i=0;i<nbPtGauss;i++)
    {
      volumeGaussPoint &pg = pt(i);
      interpolate(pg.fcts,val);
      mPoint pos(pg.x,pg.y,pg.z);
      theConservationLaw->getExactSolution()->eval(pos,time,val_ex);
      p_ex = (gamma-1.)*(val_ex[3] - 0.5*(val_ex[1]*val_ex[1]+
					  val_ex[2]*val_ex[2])/val_ex[0]);
      p = (gamma-1.)*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
      ent_error = (p/p_ex)/pow(val[0]/val_ex[0],gamma)-1.;
      errorL1+= fabs(ent_error)* pg.JacTimesWeight;
    }
  return errorL1;
}

double DGCell::entropyErrorL1_3D(double time)
{
  double gamma=1.4;
  double errorL1 = 0;
  double val[MaxNbEqn],val_ex[MaxNbEqn];
  double ent_error,p,p_ex;
  for(int i=0;i<nbPtGauss;i++)
    {
      volumeGaussPoint &pg = pt(i);
      interpolate(pg.fcts,val);
      mPoint pos(pg.x,pg.y,pg.z);
      theConservationLaw->getExactSolution()->eval(pos,time,val_ex);
      p_ex = (gamma-1.)*(val_ex[3] - 0.5*(val_ex[1]*val_ex[1]+
					  val_ex[2]*val_ex[2]+
					  val_ex[4]*val_ex[4])/val_ex[0]);
      p = (gamma-1.)*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2]
+val[4]*val[4])/val[0]);
      ent_error = (p/p_ex)/pow(val[0]/val_ex[0],gamma)-1.;
      errorL1+= fabs(ent_error)* pg.JacTimesWeight;
    }
  return errorL1;
}

/*
double DGCell::entropyErrorL2(double time)
{
  double gamma=1.4;
  double errorL2 = 0;
  double val[MaxNbEqn],val_ex[MaxNbEqn];
  double ent_error,p,p_ex;
  for(int i=0;i<nbPtGauss;i++)
    {
      volumeGaussPoint &pg = pt(i);
      interpolate(pg.fcts,val);
      mPoint pos(pg.x,pg.y,pg.z);
      theConservationLaw->getExactSolution()->eval(pos,time,val_ex);
      p_ex = (gamma-1.)*(val_ex[3] - 0.5*(val_ex[1]*val_ex[1]+
					  val_ex[2]*val_ex[2])/val_ex[0]);
      p = (gamma-1.)*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
      ent_error = (p/p_ex)/pow(val[0]/val_ex[0],gamma)-1.;
      errorL2+= ent_error*ent_error* pg.JacTimesWeight;
    }
  return errorL2;
}
*/
double DGCell::entropyErrorL2_2D(double time)
{
  double gamma = 1.4;
  double error = 0.;
  double val[MaxNbEqn],val_ex[MaxNbEqn];
  double ent_error,p,p_ex;
  
  int order = 2*(theFunctionSpace->order()+1);
  double u,v,w,weight,x,y,z;
  GaussIntegrator gauss;
  int Nb = gauss.nbIntegrationPoints(theMeshEntity,order);
  
  for(int i=0;i<Nb;i++)
    {
      gauss.iPoint(theMeshEntity,i,order,u,v,w,weight);
      double detJac = theMapping->detJac(u,v,w); 
      theMapping->eval(u,v,w,x,y,z);
           
      interpolate(u,v,w,val); 
      mPoint pos(x,y,z);
      getConservationLaw()->getExactSolution()->eval(pos,time,val_ex);
      
      p_ex = (gamma-1.)*(val_ex[3] - 0.5*(val_ex[1]*val_ex[1]+
					  val_ex[2]*val_ex[2])/val_ex[0]);
      p = (gamma-1.)*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
      ent_error = (p/p_ex)/pow(val[0]/val_ex[0],gamma)-1.;
      error+= ent_error*ent_error* detJac*weight;
    }
  return error;
}
double DGCell::entropyErrorL2_3D(double time)
{
  double gamma = 1.4;
  double error = 0.;
  double val[MaxNbEqn],val_ex[MaxNbEqn];
  double ent_error,p,p_ex;
  
  int order = 2*(theFunctionSpace->order()+1);
  double u,v,w,weight,x,y,z;
  GaussIntegrator gauss;
  int Nb = gauss.nbIntegrationPoints(theMeshEntity,order);
  
  for(int i=0;i<Nb;i++)
    {
      gauss.iPoint(theMeshEntity,i,order,u,v,w,weight);
      double detJac = theMapping->detJac(u,v,w); 
      theMapping->eval(u,v,w,x,y,z);
           
      interpolate(u,v,w,val); 
      mPoint pos(x,y,z);
      getConservationLaw()->getExactSolution()->eval(pos,time,val_ex);
      
      p_ex = (gamma-1.)*(val_ex[3] - 0.5*(val_ex[1]*val_ex[1]+
					  val_ex[2]*val_ex[2]+
					  val_ex[4]*val_ex[4])/val_ex[0]);
      p = (gamma-1.)*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2]+
				    val_ex[4]*val_ex[4])/val[0]);
      ent_error = (p/p_ex)/pow(val[0]/val_ex[0],gamma)-1.;
      error+= ent_error*ent_error* detJac*weight;
    }
  return error;
}

void DGBoundaryCell::computeEntropy(int &nbPoints,double *x,double *m) const
{
  double field[5];
  double gm1=0.4;
  double gamma = 1.4;

  nbPoints = nbPtGauss;
  for(int i=0;i<nbPtGauss;i++)
    {
      boundaryGaussPoint *pg = pt(i);
      left->interpolate(pg->fctleft,field);
      x[i] = pg->x;
      double pres = gm1*(field[3] - 
			 0.5*(field[1]*field[1]+field[2]*field[2])/field[0]); 
      double exact[4]; 
      mPoint p(0,0,0);
      left->getConservationLaw()->getExactSolution()->eval(p,0,exact); 
      double p_ex = gm1*(exact[3] -  
			 0.5*(exact[1]*exact[1]+exact[2]*exact[2])/exact[0]);
      m[i]= pres/p_ex/pow(field[0]/exact[0],gamma)-1.;
  
    }
}

double DGCell::pressureErrorL2(double time)
{
  double gamma = 1.4;
  double error = 0.;
  double val[MaxNbEqn],val_ex[MaxNbEqn];
  double pres_error,p,p_ex;
  
  int order = 2*(theFunctionSpace->order()+1);
  double u,v,w,weight,x,y,z;
  GaussIntegrator gauss;
  int Nb = gauss.nbIntegrationPoints(theMeshEntity,order);
  
  for(int i=0;i<Nb;i++)
    {
      gauss.iPoint(theMeshEntity,i,order,u,v,w,weight);
      double detJac = theMapping->detJac(u,v,w); 
      theMapping->eval(u,v,w,x,y,z);
           
      interpolate(u,v,w,val); 
      mPoint pos(x,y,z);
      getConservationLaw()->getExactSolution()->eval(pos,time,val_ex);
      
      p_ex = (gamma-1.)*(val_ex[3] - 0.5*(val_ex[1]*val_ex[1]+
					  val_ex[2]*val_ex[2])/val_ex[0]);
      p = (gamma-1.)*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
      pres_error = p-p_ex;
      error+= pres_error*pres_error* detJac*weight;
    }
  return error;
}

double DGCell::densityErrorL2(double time)
{
  double gamma = 1.4;
  double error = 0.;
  double val[MaxNbEqn],val_ex[MaxNbEqn];
    double dens_error;
  
  int order = 2*(theFunctionSpace->order()+1);
  double u,v,w,weight,x,y,z;
  GaussIntegrator gauss;
  int Nb = gauss.nbIntegrationPoints(theMeshEntity,order);

  for(int i=0;i<Nb;i++)
    {
      gauss.iPoint(theMeshEntity,i,order,u,v,w,weight);
      double detJac = theMapping->detJac(u,v,w); 
      theMapping->eval(u,v,w,x,y,z);
           
      interpolate(u,v,w,val); 
      mPoint pos(x,y,z);
      getConservationLaw()->getExactSolution()->eval(pos,time,val_ex);

      dens_error = val[0]-val_ex[0];
      error+= dens_error*dens_error* detJac*weight;
    }
  return error;
}

double DGCell::computeMass() 
{
  double val[MaxNbEqn];
  double mass =0.0;;
  for(int i=0;i<nbPtGauss;i++)
    {
      volumeGaussPoint &pg = pt(i);
      interpolate(pg.fcts,val);
      mass +=val[0]* pg.JacTimesWeight;
    }
  return mass;
}


double DGBoundaryCell::totalPressure(double &length)
{
  double gamma=1.4;
  double L2pressure = 0;
  length = 0;
  double val[MaxNbEqn],exact[MaxNbEqn];
  double p,p_total,p_total_inf,p_inf,time=0.;
  
  int order = 2*(left->theFunctionSpace->order()+1) + left->theMapping->detJacOrder() +2;
  double u,v,w,weight,x,y,z,uleft,vleft,wleft;
  GaussIntegrator gauss;
  int Nb = gauss.nbIntegrationPoints(theBoundaryEntity,order);
  for(int i=0;i<Nb;i++)
    {
      gauss.iPoint(theBoundaryEntity,i,order,u,v,w,weight);
      double detJac = theMapping->detJac(u,v,w); 
      theMapping->eval(u,v,w,x,y,z);
      left->theMapping->invert(x,y,z,uleft,vleft,wleft);
      
      left->interpolate(uleft,vleft,wleft,val); 
      mPoint pos(x,y,z);
      left->getConservationLaw()->getExactSolution()->eval(pos,time,exact);
      
      p = (gamma-1.)*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
      p_inf = (gamma-1.)*(exact[3] - 0.5*(exact[1]*exact[1]+exact[2]*exact[2])/exact[0]);
      p_total_inf = p_inf*pow(1.+0.5*(gamma-1.)*(exact[1]*exact[1]+exact[2]*exact[2])/(exact[0]*gamma*p_inf),gamma/(gamma-1.));
      p_total=p*pow(1.+0.5*(gamma-1.)*(val[1]*val[1]+val[2]*val[2])/(val[0]*gamma*p),gamma/(gamma-1.));
      L2pressure+= (p_total-p_total_inf)*(p_total-p_total_inf)* detJac*weight;
      length += detJac*weight;
    }
  return L2pressure;
}
/*
double DGBoundaryCell::totalPressure(double &length)
{
  double gamma=1.4;
  double L2pressure = 0;
  double time=0;
  length = 0;
  double val[MaxNbEqn];
  double p,p_inf,p_total,p_total_inf,exact[MaxNbEqn];
  length =0;
  for(int i=0;i<nbPtGauss;i++)
    {
      boundaryGaussPoint *pg = pt(i);
      left->interpolate(pg->fctleft,val);
      mPoint pos(pg->x,pg->y,pg->z);
      left->getConservationLaw()->getExactSolution()->eval(pos,time,exact);
      
      p = (gamma-1.)*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
      p_inf = (gamma-1.)*(exact[3] - 0.5*(exact[1]*exact[1]+exact[2]*exact[2])/exact[0]);
      p_total_inf = p_inf*pow(1.+0.5*(gamma-1.)*(exact[1]*exact[1]+exact[2]*exact[2])/(exact[0]*gamma*p_inf),gamma/(gamma-1.));
      p_total=p*pow(1.+0.5*(gamma-1.)*(val[1]*val[1]+val[2]*val[2])/(val[0]*gamma*p),gamma/(gamma-1.));
      L2pressure+= (p_total-p_total_inf)*(p_total-p_total_inf)*pg->JacTimesWeight;
      length += pg->JacTimesWeight;
     }
  return L2pressure;
  }
*/
void DGBoundaryCell::computeLiftDrag(double &lift, double &drag)
{
  double gamma=1.4;
  double val[MaxNbEqn];
  double p;
  lift = 0.; drag = 0.0;
  for(int i=0;i<nbPtGauss;i++)
    {
      boundaryGaussPoint *pg = pt(i);
      left->interpolate(pg->fctleft,val);
      p = (gamma-1.)*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
      lift+= p*n(1)* pg->JacTimesWeight;
      drag+= p*n(0)* pg->JacTimesWeight;
    }
}

void DGBoundaryCell::computePressureCoefficient(int &nbPoints,double *x,double *p) const
{
  mPoint point;
  double val[5];
  double gm1=0.4;
  //  double val[5],v_inf[5],p,p_inf;
  // theConservationLaw->getExactSolution()->eval(point,0.0,v_inf);
  // p_inf = (gamma-1.)*(v_inf[3] - 0.5*(v_inf[1]*v_inf[1]+v_inf[2]*v_inf[2])/v_inf[0]);
  // double denom = 0.5 * (v_inf[1]*v_inf[1] + v_inf[2]*v_inf[2])/v_inf[0];
  nbPoints = nbPtGauss;
  for(int i=0;i<nbPtGauss;i++)
    {
      boundaryGaussPoint *pg = pt(i);
      left->interpolate(pg->fctleft,val);
      p[i]= gm1*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
      x[i] = pg->x;
      // presCoef[i] = (p_inf-p)/denom;
    }

}

void DGBoundaryCell::computeMachNumber(int &nbPoints,double *x,double *m) const
{
  double field[5];
  double gm1=0.4;
  double gamma = 1.4;

  nbPoints = nbPtGauss;
  for(int i=0;i<nbPtGauss;i++)
    {
      boundaryGaussPoint *pg = pt(i);
      left->interpolate(pg->fctleft,field);
      x[i] = pg->x;
      double inv_rho = 1./field[0];
      double vsq = (field[1]*field[1]+field[2]*field[2])*inv_rho;
      double p = gm1 * (field[3] - 0.5*vsq);
      double csq =  gamma*p*inv_rho;
      m[i] = sqrt(vsq*inv_rho/csq);
    }
}

void DGBoundaryCell::computeTotalPressureAtAPoint(int &nbPoints,double *x, double *y,double *p_t) const
{
  mPoint point;
  double val[5],p;
  double gm1=0.4;
  double gamma =1.4;
  nbPoints = nbPtGauss;
  for(int i=0;i<nbPtGauss;i++)
    {
      boundaryGaussPoint *pg = pt(i);
      left->interpolate(pg->fctleft,val);
      p= gm1*(val[3] - 0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
      p_t[i]=p*pow(1.+0.5*(gamma-1.)*(val[1]*val[1]+val[2]*val[2])/(val[0]*gamma*p),gamma/(gamma-1.));  
      x[i] = pg->x;
      y[i] = pg->y;
    }
}

double DGCell::L1approx() 
{
  int order = theErrorFunctionSpace->order() + theMapping->detJacOrder() + 6;
  double u,v,w,weight,err[3];
  double error = 0;
  int Nb = theGaussIntegrator->nbIntegrationPoints(theMeshEntity,order);
  for(int i=0;i<Nb;i++)
    {
      theGaussIntegrator->iPoint(theMeshEntity,i,order,u,v,w,weight);
      double detJac = theMapping->detJac(u,v,w);
      interpolateError(u,v,w,err);
      error+= fabs( err[0])  * detJac * weight;
    }
  //printf("Cell error %f \n", errorL1);
   return error;
}


double DGBoundaryCell::BoundaryIntegral()
{
	double val[4],exact[4];
	double itt = 0;
	double n0,n1;
	for(int i=0;i<nbPtGauss;i++)
    {
      boundaryGaussPoint *pg = pt(i);
//		if(pg->x+pg->y>0.99 && pg->x+pg->y<1.01)
if(pg->y<0.001)
	  //if(pg->n(0)!=0 && pg->n(1)!=0)
	{
	  //if (pg->n(0)>0) {
            left->interpolate(pg->fctleft,val); n0=n(0);n1=n(1);
//}
//	 else {right->interpolate(pg->fctrght,val); n0=pg->n(0);n1=pg->n(1);}//}
     left->theConservationLaw->
		 getExactSolution()->
		 eval(mPoint(pg->x,pg->y,pg->z),1.0,exact);
	//printf("%e %e\n",pg->x,pg->y);
	 //itt+= (2.*(exact[0] - val[0])-(exact[1] - val[1]))*pg->JacTimesWeight;
//	  itt+= 0.5*((exact[0] - val[0])+(exact[1] - val[1])*n0 + (exact[2]-val[2])*n1)*pg->JacTimesWeight;
	  //itt+= ((exact[1] - val[1])*n1 - (exact[2]-val[2])*n0)*pg->JacTimesWeight;
	  //itt+= ((exact[0] - val[0]))*pg->JacTimesWeight;
//	itt+= (-(exact[0] - val[0])+3.*(exact[1]-val[1]))*pg->JacTimesWeight;
	itt+=(exact[1]/exact[0]-val[1]/val[0]);
	}
	}
	return itt;
}

double DGCell::advanceErrorInTime (double dt)
{
/*
  int n = 2 * theFunctionSpace->order() + 3 ;

  double resid = 0.0;

  for(int k = 0;k<theConservationLaw->nbFields();k++)
    {
      for(int i = 0;i<n;i++)
	{
		 double dx = 0.0;
	  for(int j = 0;j<n;j++)
	    {
	      dx += theErrorInvertMassMatrix[i][j] * theErrorRightHandSide[j+n*k];
	      ///printf("invert mass %d %d  %f\n",i,j, theErrorInvertMassMatrix[i][j]);
	      //  printf("Error RHS %d %f \n",j, theErrorRightHandSide[j+n*k]);
	    }
	  //	  printf("dx (%d,%d) = %12.5E resid = %12.5E\n",i,k,dx,resid);
	  resid += dx*dx;
	  (theErrorCoefficients[k])[i] += dt * dx;
	  //printf("(theErrorCoefficients[k])[i]  %f \n",(theErrorCoefficients[k])[i] );
	}
    }
*/
//  return resid;
return 0.0;
}

void DGCell::ZeroErrorRHS ()
{
  // RESET Error RIGHT HAND SIDE
  int fSize = 2 * theFunctionSpace->order() + 3 ;
  for(int j=0;j<theConservationLaw->getNbFields();j++)
      for(int k=0;k<fSize;k++)
		 theErrorRightHandSide[k+fSize*j] = 0.0;
}


void DGCell::L2ProjError (FieldEvaluator *f, vector<double> &proj)
{
  int order = f->order() + theErrorFunctionSpace->order() +4 + theMapping->detJacOrder() ;
  // printf("order in projector %d \n", order);
  double u,v,w,weight,val[256],exact[256],x,y,z;
  double *fct;
 
  int cSize = theConservationLaw->getNbFields();
  int fSize = theFunctionSpace->size();
  int eSize = theErrorFunctionSpace->size();
  int eOrder = theErrorFunctionSpace->pSize();
  int fOrder = theFunctionSpace->pSize();
  fct = new double [eSize];
  
   proj.reserve(cSize * (eOrder + fOrder));
  int i, j, k;
  for(i=0;i<cSize * (fOrder + eOrder);i++) proj[i] = 0.0;
  int Nb = theGaussIntegrator->nbIntegrationPoints(theMeshEntity,order);
  for(i=0;i<Nb;i++)
    {
      theGaussIntegrator->iPoint(theMeshEntity,i,order,u,v,w,weight);
      double detJac = theMapping->detJac(u,v,w);    
      theMapping->eval(u,v,w,x,y,z);
      f->eval(mPoint(x,y,z),0.0,exact);
      interpolate(u,v,w,val);
      //printf("gp %f %f %f %f %f %f\n",u,v,detJac,x,y,val[0]);
      theErrorFunctionSpace->fcts(u,v,w,fct);
      for(j=0;j<cSize;j++)
	{
	  for(k=0;k<eOrder;k++)
	    {
	      proj[k+(fOrder+eOrder)*j] += (exact[j]-val[j]) * fct[fSize+k] * detJac * weight;
	      //   printf(" exact[j]-val[j]) * fct[fSize+k] * detJac * weight %f %f %f %f \n",exact[j]-val[j], fct[fSize+k] , detJac , weight);
	    }
	}
           for(j=0;j<cSize;j++)
	{
	  for( k=0;k<fOrder;k++)
	    {
	      proj[k+eOrder+(eOrder + fOrder)*j] += 0;
	    }
	}
     //   for(j=0;j<cSize;j++)
//  	{
//  	  for( k=0;k<fOrder;k++)
//  	    {
//  	      proj[k+eOrder+(eOrder + fOrder)*j] += (exact[j]-val[j]) * fct[k + fSize-fOrder] * detJac * weight;
//  	      // printf(" exact[j]-val[j]) * fct[fSize+k] * detJac * weight %f %f %f %f \n",exact[j]-val[j], fct[fSize+k-fOrder] , detJac , weight);
//  	    }
//  	}
    } 
  //   for(int i=0;i<fOrder+eOrder;i++)printf(" %f ",proj[i]);
//    printf("\n");
  delete []fct;
}

void DGCell::L2ProjInitialError (FieldEvaluator *f)
{
  //  double val[256];
  for(int i=0;i<theConservationLaw->getNbFields();i++)
    {
      for(int j=0;j< (2 * theFunctionSpace->order() + 3);j++)
		  (theErrorCoefficients[i])[j] = 0;
    }
  L2ProjError(f,theErrorRightHandSide);
  advanceErrorInTime(1.0);
 
  //interpolateError(0,0,0,val);
   //printf("%f\n",val[0]);
}

void DGCell::ZeroErrorMassMatrices()
{
	int eSize = theErrorFunctionSpace->pSize();
	int fSize = theFunctionSpace->pSize();
	int j,k;
	  for( k=0;k<eSize + fSize;k++)
		for(j=0;j<eSize + fSize;j++)
	      theErrorMassMatrix [k][j] = 0.0;

	for(k=0;k<eSize;k++)
		theErrorMassMatrix [k][k] = 1.0;
}

void DGCell::computeErrorMean ()
{
  double u,v,w,weight,val[256],err[256],vol;
 
  int i,j; 
  int order = theErrorFunctionSpace->order() + theMapping->detJacOrder();
  for(i=0;i<theConservationLaw->getNbFields();i++)theMean[i] = 0.0;
  vol = 0.0;
  int Nb = theGaussIntegrator->nbIntegrationPoints(theMeshEntity,order);
    for(i=0;i<Nb;i++)
    {
      theGaussIntegrator->iPoint(theMeshEntity,i,order,u,v,w,weight);
      double detJac = theMapping->detJac(u,v,w);
      interpolate( u,v,w, val );
	  interpolateError(u,v,w,err);
      for(j=0;j<theConservationLaw->getNbFields();j++)
		  val[j]+=err[j];
      vol += detJac*weight;
      for(j=0;j<theConservationLaw->getNbFields();j++)theErrorMean[j] +=
                                weight*detJac*val[j];
    }
  for(j=0;j<theConservationLaw->getNbFields();j++)theErrorMean[j] /= vol;
}

void DGBoundaryCell::computeFullJump()
{
  if(!mright)return;

  double valleft[MaxNbEqn],valrght[MaxNbEqn],jump;
  boundaryGaussPoint *pg;
  double error = 0; 
  for(int i=0;i<nbPtGauss;i++)
    {
      pg = pt(i);
      left->interpolate(pg->fctleft,valleft);
      right ->interpolate(pg->fctrght,valrght);
      jump = left->theConservationLaw->jumpQuantity(valleft)
	- right->theConservationLaw->jumpQuantity(valrght); 
      error += (jump*jump*jump*jump) * pg->JacTimesWeight;
    }
  left->fullJump += error;
  right->fullJump += error;
 }

void DGBoundaryCell :: IndicatorGrad (double *denoms, double *numers) const
{
  double grad1,gradb;
  double ul,vl,wl;
  double ur,vr,wr;
  double xl,yl,zl,xr,yr,zr,dx;
  double vall[MaxNbEqn],valr[MaxNbEqn];
  //int nbf = GlobalParamManager::Instance()->getLaw()->nbFields();
  mVector *gradsl;
  gradsl = new mVector[cSize]; 

  left->getMapping()->COG(ul,vl,wl);
  if(right)
    right->getMapping()->COG(ur,vr,wr);

  left->interpolate(ul,vl,wl,vall);
  if(right)
    right->interpolate(ur,vr,wr,valr);

  left->getMapping()->eval(ul,vl,wl,xl,yl,zl);
  if(right)
    right->getMapping()->eval(ur,vr,wr,xr,yr,zr);

  if(!right)
    {
      left->dinterpolate(ul,vl,wl,gradsl);
      double  f  = vall[0];
      grad1 = gradsl[0] * gradsl[0] / (f*f);      
    }
  if(right)
    {
      dx = sqrt((xl-xr)*(xl-xr)+(yl-yr)*(yl-yr)+(zl-zr)*(zl-zr));
      double  f  = 0.5 *  (vall[0] + valr[0]);
      double df = vall[0] - valr[0];
      gradb = (1./f) * df/dx;
      gradb *= gradb;
    }
    
  delete[] gradsl;

  if(right)
    numers[0] = gradb;
  else
    numers[0] = grad1;
    
  denoms[0] = 1.0;
}
