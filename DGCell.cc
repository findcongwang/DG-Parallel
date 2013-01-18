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
#include "Geometry.h"
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <cstdlib>

extern void invmat (double **, double **, int);

double ** allocateMatrix(int n)
{
  int i;
  double ** v = new double *[n];
  double *help = new double[n*n];
  for (i=0; i!= n; ++i)
    v[i] = &help[i*n];
  return v;
}



void freeMatrix( double **v)
{
  if(v)
    {
      if(v[0]) delete [] v[0];
      delete [] v;
    }
}


DGCell::DGCell (ConservationLaw*l, 
		mEntity*e, 
		FunctionSpace*f,
		FunctionSpace *er,
		Mapping *m,
		GaussIntegrator *I)
  : theConservationLaw(l),theMeshEntity(e),theFunctionSpace(f),
    theMapping(m),theGaussIntegrator(I)
{
  //theGeometry = g;
  int i;
  fSize                 = f->size();
  fOrder                = f->order();
  cSize                 = l->getNbFields();
  theFieldsCoefficients = new DG_Dofs(2,cSize,fSize);
  theMean               = new double[cSize];
  limSlope              = new double[cSize];
  theRightHandSide      = new double [fSize*cSize];
  deltaUp0              = new double[cSize];
  Up0                   = new double[cSize];
  Rp0                   = new double[cSize];
  
  for(i=0;i<cSize;++i)
    for(int j=0;j<fSize;++j) theFieldsCoefficients->get(i,j) = 0.;
  limit = 0;
  order = 2 * fOrder + theMapping->detJacOrder();

  theInvertMassMatrix = 0;
  init();
  
  pMin = pMax = ((Vertex*)theMeshEntity->get(0,0))->point();
  int Size0 = theMeshEntity->size(0);
  for(i=1;i<Size0;++i)
    {
      mPoint p1 = ((Vertex*)theMeshEntity->get(0,i))->point();
      if (p1(0) < pMin(0)) pMin(0) = p1(0);
      if (p1(1) < pMin(1)) pMin(1) = p1(1);
      if (p1(2) < pMin(2)) pMin(2) = p1(2);
      if (p1(0) > pMax(0)) pMax(0) = p1(0);
      if (p1(1) > pMax(1)) pMax(1) = p1(1);
      if (p1(2) > pMax(2)) pMax(2) = p1(2);
    }

  computeCellSize();
  computeVolume();
  computePerimeter();
  ZeroError();
  computeMaxH();
}

DGCell::~DGCell ()
{
  delete theFieldsCoefficients;
  delete [] theMean;
  delete [] limSlope;
  delete [] theRightHandSide;
  delete theFunctionSpace;
  freeMatrix(theInvertMassMatrix);
}

void DGCell::cleanup()
{
  delete theMapping;
  theMapping = 0;
}

volumeGaussPoint::volumeGaussPoint (const volumeGaussPoint &other)
{
  JacTimesWeight = other.JacTimesWeight;
  x = other.x;
  y = other.y;
  z = other.z;
  int fSize=other.grads.size();
  grads.reserve(fSize);
  grads.resize(fSize);
  fcts = other.fcts;
  for(int i=0;i<fSize;++i)
    grads[i] = other.grads[i];
}

void DGCell::init ()
{
  int i,j,k;
  double** theMassMatrix;
  if (!theFunctionSpace->isOrthogonal())
    {
      freeMatrix(theInvertMassMatrix);
      theMassMatrix       = allocateMatrix(fSize);
      theInvertMassMatrix = allocateMatrix(fSize);
      for(j=0;j<fSize;++j)
	for(k=0;k<fSize;++k)
	  theMassMatrix[j][k] = 0.0;
    }
  
  nbPtGauss = theGaussIntegrator->nbIntegrationPoints(theMeshEntity,order);
  volumeGaussPoints.reserve(nbPtGauss);
  //volumeGaussPoints.resize( nbPtGauss );
  int VP = volumeGaussPoints.size();
  //Computing the integration points, weights, Jacobian.
  //x,y,z - physical coordinates
  double u,v,w,weight;
  for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint pg;
	  theGaussIntegrator->iPoint(theMeshEntity,i,order,u,v,w,weight);
	 // u=0;v=1;w=0;
      theMapping->eval(u,v,w,pg.x,pg.y,pg.z);
      pg.grads.reserve(fSize);
      pg.grads.resize(fSize);
      detJac = theFunctionSpace->grads(u,v,w,theMapping,pg.grads);
      pg.JacTimesWeight = detJac*weight;
      for(j=0;j<fSize;++j) pg.grads[j] *= (pg.JacTimesWeight);
      //values of basis functions at (u,v,w)
      //theFunctionSpace->fcts(u,v,w,pg.fcts);
      pg.fcts=theGaussIntegrator->iFct(theMeshEntity,theFunctionSpace,i,order); 
      volumeGaussPoints.push_back(pg);
		//volumeGaussPoints[ i ] = pg;

      ///Computing the Mass Matrix
      if (!theFunctionSpace->isOrthogonal())
	{
	  for( j=0;j<fSize;++j) {
	    for( k=0;k<fSize;++k){
	      theMassMatrix [j][k] += pg.fcts[k] * pg.fcts[j] * pg.JacTimesWeight;
	      //printf("%6.3f",theMassMatrix [j][k]);

		}
	  }
	}
    }
  if (!theFunctionSpace->isOrthogonal())
    {
      invmat(theMassMatrix,theInvertMassMatrix,fSize);
      freeMatrix(theMassMatrix);
    }
}

void DGCell::computeVolumeContribution (double t)
{
  // we compute int_{Cell} F grad w dCell
  int i,j,k;
  double rhs[MaxNbEqn];
  double *RHS;
  double val[MaxNbEqn];
  mVector fluxes[MaxNbEqn];
 
  for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint  &pg = pt(i);
      interpolate(pg.fcts, val);
      
      if(!theConservationLaw->isPhysical(val))
	{
	  printf("non physical field found in VOLUME integral at point (%f,%f,%f) \n",pg.x, pg.y, pg.z ); 
	  //  printf("%f %f %f %f\n",val[0],val[1],val[2],val[3]);
	  double p;
	  int dim = theMeshEntity->getLevel();
	  if (dim==2)
	    p = (0.4)*(val[3]-0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
	  else 
	    p = (0.4)*(val[3]-0.5*(val[1]*val[1]+val[2]*val[2]+val[4]*val[4])/val[0]);
	  printf("rho=%e p=%e \n",val[0],p);
	  
	 
	  for(j=0;j<cSize;++j)
	    for (k=1;k<fSize;++k)
	      theFieldsCoefficients->get(j,k)= 0.0;
	  interpolate(pg.fcts, val);
	  if (dim==2)
	    p = (0.4)*(val[3]-0.5*(val[1]*val[1]+val[2]*val[2])/val[0]);
	  else 
	    p = (0.4)*(val[3]-0.5*(val[1]*val[1]+val[2]*val[2]+val[4]*val[4])/val[0]);
	  printf("Corrected values rho=%e p=%e \n",val[0],p);
	  if (p<0) exit(0);
	  ZeroRHS();
	  computeVolumeContribution(t);
	  break;
	}
      else
	{
	  const mPoint p(pg.x,pg.y,pg.z);  
	  theConservationLaw->RHS(val,p,t,rhs);
	  theConservationLaw->Fi(p,val,fluxes);
	  RHS = theRightHandSide;
	  const double *fcts = pg.fcts;
	  const double *rhs_const =rhs;
	  double tmp=0;
	  for(j=0;j<cSize;++j)
	    for(k=0;k<fSize;++k)		    
		{
	      (*RHS++)+= fluxes[j] *pg.grads[k]  + rhs_const[j] * fcts[k]*pg.JacTimesWeight;
	 //    for(int m=0;m<3;++m) printf("%e\n",pg.grads[k](m)); printf("\n");
		  //tmp+= fluxes[j] *pg.grads[k]  + rhs_const[j] * fcts[k]*pg.JacTimesWeight; //TEST
		}
		//printf("\n volume contribution rhs = %e \n",tmp); //TEST
	}
    }
 // printf("\n in volume \n");
  //for(j=0;j<3;++j) printf("%e\n",theRightHandSide[j]);
    /*for(j=0;j<cSize;++j)
	    for (k=0;k<fSize;++k)
	      printf("%e\n",theFieldsCoefficients->get(j,k));*/
	//printf("\n"); 
}

void DGCell::reverseVolumeContribution (double t)
{
  // we compute int_{Cell} F grad w dCell
  int i,j,k;
  double rhs[MaxNbEqn];
  double *RHS;
  double val[MaxNbEqn];
  mVector fluxes[MaxNbEqn];
 
  for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint  &pg = pt(i);
      interpolate(pg.fcts, val);
      
      const mPoint p(pg.x,pg.y,pg.z);  
      theConservationLaw->RHS(val,p,t,rhs);
      theConservationLaw->Fi(p,val,fluxes);
      RHS = theRightHandSide;
      const double *fcts = pg.fcts;
      const double *rhs_const =rhs;
      double tmp=0;
      for(j=0;j<cSize;++j)
	for(k=0;k<fSize;++k)
	{
	  (*RHS++)-= fluxes[j] *pg.grads[k]  + rhs_const[j] * fcts[k]*pg.JacTimesWeight;
   
	}
    }
  //printf("\n");
  //for(j=0;j<9;++j) printf("%e\n",theRightHandSide[j]);
}

void DGCell::dinterpolate(double u, double v, double w, mVector *grads)
{

  //  if(theMeshEntity->isAdjacencyCreated(2))printf("oups ...\n");

  vector<mVector>  gr;
  gr.reserve(fSize);  
  theFunctionSpace->grads(u,v,w,theMapping,gr);

  mVector g = gr[0];
  for(int i=0;i<cSize;++i)
    {
      // grads[i] = g * (theFieldsCoefficients[i])[0];
      grads[i] = g * theFieldsCoefficients->get(i,0);
    }
  
  for(int j=1;fSize;++j)
    {
      g = gr[j];
      for(int i=0;i<cSize;++i)
	{
	  grads[i] += (g * theFieldsCoefficients->get(i,j));
	}
    }
}
    
double DGCell::getError() const
{
	return sqrt(error);
}


void DGCell::computeMaxH()
{
double h = 0.;
mPoint p0 =((Vertex*)theMeshEntity->get(0,0))->point();
  
  for(int i=1;i<theMeshEntity->size(0);++i)
    {
      mPoint p1 = ((Vertex*)theMeshEntity->get(0,i))->point();
      mVector v(p0,p1);
	  if (h * h < v*v) h = sqrt(v*v);
	  p0 = p1;
	 } 	  
  maxH = h;
}

void DGCell::computeMinH()
{
double h = 0.;
mPoint p0 =((Vertex*)theMeshEntity->get(0,0))->point();
  
  for(int i=1;i<theMeshEntity->size(0);++i)
    {
      mPoint p1 = ((Vertex*)theMeshEntity->get(0,i))->point();
      mVector v(p0,p1);
	  if (h > v*v) h = v*v;
	  p0 = p1;
	 } 
  minH = h;
}

double DGCell::getMaxH() const
{
	return maxH;
}

double DGCell::getMinH() const
{
	return minH;
}

/************ B O U N D A R Y   T E R M S *************/

void DGBoundaryCell::check()
{
  mEntity *b[2] = {0,0};
  int k = 0;
  int n = theBoundaryEntity->getLevel();
  
  if(!theBoundaryEntity->isAdjacencyCreated(n+1))
    {
      printf("weird face because no upward adj :");theBoundaryEntity->print();
      if(theBoundaryEntity->parent())
	{
	  printf("this face has parent : ");
	  theBoundaryEntity->parent()->print();
	}      
      if(theBoundaryEntity->isAdjacencyCreated(2))
	{
	  printf("this face has childeren and should not be there\n ");
	}      
    }
mUpwardAdjacencyContainer::iter it;
mUpwardAdjacencyContainer::iter end_np1 = theBoundaryEntity->end(n+1);

  for( it= theBoundaryEntity->begin(n+1); it !=  end_np1; ++it)
    {
      b[k++] = (*it);
    }

  if(mleft != b[0] || mright != b[1])
    { 
      mleft = b[0];
      mright = b[1];
      right = 0;
      left = (DGCell*)mleft->getCell();
      if(mright)right = (DGCell*)mright->getCell();
      init();
    }
}

DGBoundaryCell::DGBoundaryCell (mEntity *ent, Mapping *m,GaussIntegrator *i)
  : theBoundaryEntity(ent),mleft(0),mright(0),theMapping(m),
    theGaussIntegrator(i) 
{
  left=right=0;
  check();
}


DGBoundaryCell::~DGBoundaryCell ()
{
  delete theMapping;
  for(int i=0;i<nbPtGauss;++i) delete boundaryGaussPoints[i];
}

int DGBoundaryCell::computeOrder() const
{
  if(right) return 
      2 * ((left->fOrder>right->fOrder) ?
	  left->fOrder : right->fOrder) +  theMapping->detJacOrder()+1;
  else return 2*(left->fOrder) + theMapping->detJacOrder()+1;
}


boundaryGaussPoint::boundaryGaussPoint()
{}

boundaryGaussPoint::~boundaryGaussPoint()
{}

void DGBoundaryCell::init()
{
  int i,j,k,qleft,qright;
  int order = computeOrder();
  cSize = left->theConservationLaw->getNbFields();
  lSize = left->fSize;
  if (right) rSize = right->fSize; else rSize = 0;
  //theBoundaryEntity->print();
  if(boundaryGaussPoints.size())
    {
      for(i=0;i<boundaryGaussPoints.size();++i)delete boundaryGaussPoints[i];
      boundaryGaussPoints.clear();
    }

  nbPtGauss = theGaussIntegrator->nbIntegrationPoints(theBoundaryEntity,order);
  boundaryGaussPoints.reserve(nbPtGauss);
  double u,v,w,weight,uleft,vleft,wleft,urght,vrght,wrght;
   
  for(j=0;j<left->theMeshEntity->size(0);j++)
	{
	  Vertex* vV = (Vertex*) left->theMeshEntity->get(0,j);
	  int found = 0;
	  for(k=0;k<theBoundaryEntity->size(0);k++) 
	    {
	      Vertex* vB = (Vertex*) theBoundaryEntity->get(0,k);
	      if (vV==vB) {found = 1; break;}
	    }
	  if (found == 0) {qleft=j-1; if (qleft<0) qleft=3; break;} 
	}
  if (right) {
     for(j=0;j<right->theMeshEntity->size(0);j++)
	{
	  Vertex* vV = (Vertex*) right->theMeshEntity->get(0,j);
	  int found = 0;
	  for(k=0;k<theBoundaryEntity->size(0);k++) 
	    {
	      Vertex* vB = (Vertex*) theBoundaryEntity->get(0,k);
	      if (vV==vB) {found = 1; break;}
	    }
	  if (found == 0) {qright=j-1; if (qright<0) qright=3; break;} 
	}
  }

  for(i=0;i<nbPtGauss;++i)
    {
      boundaryGaussPoint *pg = new boundaryGaussPoint();
      theGaussIntegrator->iPoint(theBoundaryEntity,i,order,u,v,w,weight);
      // compute parametric coordinates on both sides
	 
      theMapping->eval(u,v,w,pg->x,pg->y,pg->z);
      
      if(!left->theMapping->invert(pg->x,pg->y,pg->z,uleft,vleft,wleft))
	{
	  printf("l : impossible to compute invert mapping in ");theBoundaryEntity->print();
	  for(int ip=0;ip<theBoundaryEntity->size(0);ip++)
	    theBoundaryEntity->get(0,ip)->print();
	}
      
      if(right)
	if(!right->theMapping->invert(pg->x,pg->y,pg->z,urght,vrght,wrght))
	  {
	    printf("r : impossible to compute invert mapping in ");theBoundaryEntity->print();
	    for(int ip=0;ip<theBoundaryEntity->size(0);ip++)
	      theBoundaryEntity->get(0,ip)->print();
	  }
      
	
      // get the normal vector
      
      if(left->theMeshEntity->find(theBoundaryEntity))
	left->theMapping->normalVector(theBoundaryEntity,uleft,vleft,wleft,pg->n);
      else
	{
	  //	  printf("computing normal in the complex way  :  ");
	  bool found = false;
	  for(j=0;j<mleft->size(theBoundaryEntity->getLevel());++j)
	    {
	      mEntity *initial = mleft->get(theBoundaryEntity->getLevel(),j); 
	      list<mEntity*> leaves;	  
	      initial->getLeaves(leaves);
	      for(list<mEntity*>::const_iterator it = leaves.begin();it != leaves.end();++it)
		{
		  if(*it == theBoundaryEntity)
		    {
		      left->theMapping->normalVector(initial,uleft,vleft,wleft,pg->n);
		     		      //printf("%f %f %f %f\n",n(0),n(1),pg->x,pg->y);
		      found = true;
		      break;
		    }
		}
	    }
	 
	  if(!found)
	    {
	      printf("unable to compute normal vector\n");
	    }
	}
	  n=pg->n; 
      // get det of jacobian
      detJac = theMapping->detJac(u,v,w);
      pg->JacTimesWeight = detJac * weight;
      // get shape functions on both sides
      pg->fctleft = new double[lSize];
      left->theFunctionSpace->fcts(uleft,vleft,wleft,pg->fctleft);
      //printf("Nb of Gauss points= %d, nb of functions= %d \n", nbPtGauss,lSize); 
    /*  printf("{");
      for (q=0; q<lSize; q++) 
	if (q<lSize-1) printf("%17.16e, ", pg->fctleft[q]);
	else printf("%17.16e}, \n", pg->fctleft[q]);*/
      
	 double* ff=0;
	 ff = left->theFunctionSpace->fcts(left->fOrder,qleft,i);
     //pg->fctleft = left->theFunctionSpace->fcts(left->fOrder,qleft,i);
	//for(j=0;j<rSize; j++) if (fabs (pg->fctleft[j]-ff[j])>0.000000001) printf (" left %d  %d %d %e %e \n",q, i, j, pg->fctleft[j],ff[j]);
	//for(q=0;q<2; q++) printf (" left %e %e \n",pg->fctleft[q],ff[q]);
	//printf("\n");
      if(right)
	{
		pg->fctrght = new double [rSize];
	 right->theFunctionSpace->fcts(urght,vrght,wrght,pg->fctrght); 
	
      
     //pg->fctrght = right->theFunctionSpace->fcts(right->fOrder,qright,i);
	//  ff = right->theFunctionSpace->fcts(right->fOrder,qright,i);
	 // for(j=0;j<rSize; j++) if (fabs (pg->fctrght[j]-ff[j])>0.000000001) printf (" right %d %d %d %e %e \n",q, i, j, pg->fctrght[j],ff[j]);
    /*for(q=0;q<2; q++) printf (" right %e %e \n",pg->fctright[q],ff[q]);
	printf("\n");*/
	}
      boundaryGaussPoints.push_back(pg);
    }
  computeSize();
  //printf("\n");
}

void DGBoundaryCell::setPeriodicBC()
{
  mEntity* ent=theBoundaryEntity->getAttachedEntity(TYPE_CELL);
  DGBoundaryCell *symm = (DGBoundaryCell *)ent->getCell();
  DGCell* another = symm->left; 
  rSize = another->fSize; 	
  if (theBoundaryEntity->getClassification()->getId()==610 || theBoundaryEntity->getClassification()->getId()==510)
    for(int i=0;i<nbPtGauss;++i)
      {
	double urght,vrght,wrght;
	boundaryGaussPoint *pg = pt(i);      
	mPoint p(pg->x,pg->y,pg->z);
	if (theBoundaryEntity->getClassification()->getId()==510)
	  another->theMapping->invert(-pg->x,pg->y,pg->z,urght,vrght,wrght);
	if (theBoundaryEntity->getClassification()->getId()==610)
	  another->theMapping->invert(pg->x,-pg->y,pg->z,urght,vrght,wrght);

	pg->fctrght = new double [rSize];
	another->theFunctionSpace->fcts(urght,vrght,wrght,pg->fctrght); 
	}
}

void DGBoundaryCell::computeError()
{
  if(!mright)return;
  double valleft[MaxNbEqn],valrght[256];
  double error = 0;
  boundaryGaussPoint *pg;
  for(int i=0;i<nbPtGauss;++i)
    {
	pg = pt(i);
      left  ->interpolate(pg->fctleft,valleft);
      right ->interpolate(pg->fctrght,valrght);
      double jump = left->theConservationLaw->jumpQuantity(valleft)
	- right->theConservationLaw->jumpQuantity(valrght); 
      error += (jump*jump) * pg->JacTimesWeight;
    }
  left->error += error;
  right->error += error;
}

void DGBoundaryCell::computeJump()
{
/*  if(!mright)return;
  DGCell *left = (DGCell*)mleft->getCell();
  DGCell *right = (DGCell*)mright->getCell();
  int order = computeOrder(left,right);
  double u,v,w,weight,valleft[MaxNbEqn],valrght[MaxNbEqn];
  double jump = 0;
  mVector n;
  
  for(int i=0;i< theGaussIntegrator.nbIntegrationPoints(theBoundaryEntity,order);++i)
    {
      theGaussIntegrator.iPoint(theBoundaryEntityi,order,u,v,w,weight);
      double detJac = theMapping->detJac(u,v,w);

      left  ->interpolate(uleft[i],vleft[i],wleft[i],valleft);
	 
	  right ->interpolate(urght[i],vrght[i],wrght[i],valrght);
      //printf("%e %e \n", valleft, valrght);
  if(left->theMeshEntity->find(theBoundaryEntity))
	left->theMapping->normalVector(theBoundaryEntity,uleft_e[i],vleft_e[i],wleft_e[i],n);
     else
	{
	  if(!right)printf("aaargh\n");
	  right->theMapping->normalVector(theBoundaryEntity,urght_e[i],vrght_e[i],wrght_e[i],n);
	  n *= -1.0;
	}
      mPoint p;
      theMapping->eval(u,v,w,p(0),p(1),p(2));
      int orient = left->theConservationLaw->getEdgeOrientation(n,p,valleft);
	  if (orient == 1) 
	  {
		  left->jump+= (valrght[0] - valleft[0]) * detJac * weight;
	  }
	  else if(right) right->jump +=(valrght[0] - valleft[0]) * detJac * weight;
	  //printf(" %e %e %e %e %e \n", valrght[0], valleft[0], detJac, weight,(valrght - valleft) * detJac * weight);
	  
  }
  */
}

int DGBoundaryCell::computeBoundaryContributions(double T)
{
  int i,j,k;
  int nonPhysical = 0;
  double riemann[MaxNbEqn];
  boundaryGaussPoint *pg;
  double valleft[MaxNbEqn],valrght[MaxNbEqn];
   
  for(i=0;i<nbPtGauss;++i)
    {
      pg = pt(i);      
      mPoint p(pg->x,pg->y,pg->z);
      left ->interpolate(pg->fctleft, valleft);
      if(!left->theConservationLaw->isPhysical(valleft))
	{
	  printf("non physical field found on an EDGE/FACE at point (%f,%f,%f)\n",pg->x,pg->y,pg->z); 
	  double p;
	  int dim = left->theMeshEntity->getLevel();
	  if (dim==2)
	    p = (0.4)*(valleft[3]-0.5*(valleft[1]*valleft[1]+valleft[2]*valleft[2])/valleft[0]);
	  else 
	    p = (0.4)*(valleft[3]-0.5*(valleft[1]*valleft[1]+valleft[2]*valleft[2]+valleft[4]*valleft[4])/valleft[0]);
	  printf("rho=%e p=%e \n",valleft[0],p);
	  printf("Boundary Id =  %d \n",theBoundaryEntity->getClassification()->getId());
 	  nonPhysical =1;
	}

      if (right)
	  {
	right ->interpolate(pg->fctrght, valrght);
	
	if(!right->theConservationLaw->isPhysical(valrght))
	  {
	    printf("non physical field found (right) on an edge at point (%f,%f,%f)\n",pg->x,pg->y,pg->z); 
	    double p;
	    int dim = left->theMeshEntity->getLevel();
	    if (dim==2)
	      p= (0.4)*(valrght[3]-0.5*(valrght[1]*valrght[1] +valrght[2]*valrght[2])/valrght[0]);
	    else
	      p= (0.4)*(valrght[3]-0.5*(valrght[1]*valrght[1] +valrght[2]*valrght[2]+valrght[4]*valrght[4])/valrght[0]);
	    printf("rho = %e p = %e \n",valrght[0],p);
	    printf("Id =  %d \n",theBoundaryEntity->getClassification()->getId());
	    if (nonPhysical==1) nonPhysical = 3;
	    else nonPhysical = 2;
	  }
	  }
  //Periodic boundary conditions
  if (theBoundaryEntity->getClassification()->getId()==610 || theBoundaryEntity->getClassification()->getId()==510)
    {
      mEntity* ent=theBoundaryEntity->getAttachedEntity(TYPE_CELL);
      DGBoundaryCell *symm = (DGBoundaryCell *)ent->getCell();
      DGCell* another = symm->left; 
	  another->interpolate(pg->fctrght, valrght);
    }	

      if(right)      //inner cell
	left->theConservationLaw->
	  riemannSolver (pg->n,p,valleft,valrght,riemann);
      // or boundary conditions
      else
	if(theBoundaryEntity->getClassification()->getId()==20000) 
	  left->theConservationLaw->boundary (pg->n,pg->N,theBoundaryEntity->getClassification()->getId(),  p, valleft,riemann,T);
//	else if (theBoundaryEntity->getClassification()->getId()==50000)
	//	left->theConservationLaw->boundarySource (pg->n,pg->N,theBoundaryEntity->getClassification()->getId(),  p, valleft,riemann,T);
	else if (theBoundaryEntity->getClassification()->getId()==510 || theBoundaryEntity->getClassification()->getId()==610) 
	  left->theConservationLaw->riemannSolver(pg->n,p,valleft,valrght,riemann);
	else left->theConservationLaw->boundaryFlux (pg->n,theBoundaryEntity->getClassification()->getId(),  p, valleft,riemann,T);
      
	//printf("%f %f %f %f %f \n",riemann[0],riemann[1],riemann[2],riemann[3],riemann[4]); 
      // left and right contributions
      double *RHS_LEFT = left->theRightHandSide;
      const double *RIEMANN = riemann;
      const double *FUNCTION_LEFT  = pg->fctleft;
	  double tmp=0;
      for(j=0;j<cSize;++j)
		  for(k=0;k<lSize;++k)
		  {
	  (*(RHS_LEFT++)) -= RIEMANN[j]*FUNCTION_LEFT[k]*pg->JacTimesWeight;
tmp-=RIEMANN[j]*FUNCTION_LEFT[k]*pg->JacTimesWeight;
		  }
	 // printf("\n RHS_LEFT from boundary = %e \n",tmp); //TEST
	 // printf("\n in boundary: \n"); for(k=0;k<lSize;++k) printf("%e \n", left->theRightHandSide[k]); printf("\n");
      if (right)
	{
	  double  *RHS_RIGHT = right->theRightHandSide;
	  const double *FUNCTION_RIGHT  = pg->fctrght;  
	  for(j=0;j<cSize;++j)
	    for(k=0;k<rSize;++k)
	      (*(RHS_RIGHT++)) += RIEMANN[j]*FUNCTION_RIGHT[k]*pg->JacTimesWeight;
	}
    }
//  printf("---------");
  return nonPhysical;
}

void DGBoundaryCell::reverseBoundaryContributions(double T)
{
  int i,j,k;
  double riemann[MaxNbEqn];
  boundaryGaussPoint *pg;
  double valleft[MaxNbEqn],valrght[MaxNbEqn];
  
  for(i=0;i<nbPtGauss;++i)
    {
      pg= pt(i);      
      left ->interpolate(pg->fctleft, valleft);
      if(right) right ->interpolate(pg->fctrght, valrght);
      
      //Periodic boundary conditions
      if (theBoundaryEntity->getClassification()->getId()==610 || theBoundaryEntity->getClassification()->getId()==510)
	{
	  mEntity* ent=theBoundaryEntity->getAttachedEntity(TYPE_CELL);
	  DGBoundaryCell *symm = (DGBoundaryCell *)ent->getCell();
	  DGCell* another = symm->left; 
	  another->interpolate(pg->fctrght, valrght);
	}
      
      mPoint p(pg->x,pg->y,pg->z);
      if(right)      //inner cell
	left->theConservationLaw->riemannSolver(pg->n,p,valleft,valrght,riemann);
      // or boundary conditions
      else
	if(theBoundaryEntity->getClassification()->getId()==20000) 
	  left->theConservationLaw->boundary (pg->n,pg->N,theBoundaryEntity->getClassification()->getId(),  p, valleft,riemann,T);
	else if (theBoundaryEntity->getClassification()->getId()==510 || theBoundaryEntity->getClassification()->getId()==610) 
	  left->theConservationLaw->riemannSolver(pg->n,p,valleft,valrght,riemann);
	else left->theConservationLaw->boundaryFlux (pg->n,theBoundaryEntity->getClassification()->getId(),  p, valleft,riemann,T);
      
      
      // left and right contributions
      double *RHS_LEFT = left->theRightHandSide;
      const double *RIEMANN = riemann;
      const double *FUNCTION_LEFT  = pg->fctleft;
      for(j=0;j<cSize;++j)
	for(k=0;k<lSize;++k)
	  (*(RHS_LEFT++)) += RIEMANN[j]*FUNCTION_LEFT[k]*pg->JacTimesWeight;
      
      if (right)
	{
	  double  *RHS_RIGHT = right->theRightHandSide;
	  const double *FUNCTION_RIGHT  = pg->fctrght;  
	  for(j=0;j<cSize;++j)
	    for(k=0;k<rSize;++k)
	      (*(RHS_RIGHT++)) -= RIEMANN[j]*FUNCTION_RIGHT[k]*pg->JacTimesWeight;
	}
    }
}

void DGCell::setToZero(double t)
{
  int i,j,k;
  int level = theMeshEntity->getLevel();
  int nbBound = theMeshEntity->size(level-1);
  mEntity* ent;
  for(i=0;i<nbBound;i++) 
    {
      ent = theMeshEntity->get(level-1,i);
      DGBoundaryCell* cell=(DGBoundaryCell *)ent->getCell();
      cell->reverseBoundaryContributions(t);
    }
  reverseVolumeContribution(t);

  for(j=0;j<cSize;++j)
    for(k=1;k<fSize;++k) 
      theFieldsCoefficients->get(j,k) = 0.0;

  computeVolumeContribution(t);

  for(i=0;i<nbBound;i++) 
    {
      ent = theMeshEntity->get(level-1,i);
      DGBoundaryCell* cell=(DGBoundaryCell *)ent->getCell();
      cell->computeBoundaryContributions(t);
    }
}

void DGBoundaryCell::computeSize() //higher order not yet coded here.
{
  switch(theBoundaryEntity->getType())
    {
  case mEntity::VERTEX:
      size = 0.0;
      break;
  case mEntity::EDGE:        /*Edge*/
      {
	Vertex* v1 = (Vertex* )theBoundaryEntity->get(0,0);
	Vertex* v2 = (Vertex* )theBoundaryEntity->get(0,1);
	mPoint p1 = v1->point();
	mPoint p2 = v2->point();
	size = sqrt((p1(0)-p2(0))*(p1(0)-p2(0))+(p1(1)-p2(1))*(p1(1)-p2(1))+(p1(2)-p2(2))*(p1(2)-p2(2)));
     	}
      break;
  case mEntity::TRI:        /*Triangle*/
      size = 0.5*detJac;
      //printf("Needs to be done");
      break;
  default: printf("no size for diemnsions greater than 2\n");
    }
}

 void DGBoundaryCell::normalToCircle(mPoint &p1, mPoint &p2,mPoint &p3)
{
  double x,y,r;
  double c1 = p1(0)*p1(0) + p1(1)*p1(1);
  double c2 = c1 - p2(0)*p2(0) - p2(1)*p2(1);
  double c3 = c1 - p3(0)*p3(0) - p3(1)*p3(1);
  double d1 = p1(0)-p3(0);
  double d2 = p1(0)-p2(0);
  double d3 = p1(1)-p3(1);
  double d4 = p1(1)-p2(1);
  double center_x = 0.5*(c2*d3-c3*d4)/(d2*d3-d1*d4);
  double center_y = 0.5*(c2*d1-c3*d2)/(d4*d1-d3*d2);
  for(int i=0;i<nbPtGauss;++i)
    {
      x = boundaryGaussPoints[i]->x - center_x;
      y = boundaryGaussPoints[i]->y - center_y;
      r = sqrt(x*x +y*y);
	  //printf("radius %e \n", r);
	   if (n(0)*x +n(1)*y>0)
	   {
      boundaryGaussPoints[i]->N(0) = x/r;  
      boundaryGaussPoints[i]->N(1) = y/r;  
      boundaryGaussPoints[i]->N(2) = 0.0;
	   } else 
	   {
      boundaryGaussPoints[i]->N(0) = -x/r;  
      boundaryGaussPoints[i]->N(1) = -y/r;  
      boundaryGaussPoints[i]->N(2) = 0.0;
	   }
    }
}

 void DGBoundaryCell::normalToCircle(double center_x, double center_y, double center_z)
{
  double x,y,z,r;
  //printf(" %e %e \n",center_x,center_y);
  for(int i=0;i<nbPtGauss;++i)
    {
      x = boundaryGaussPoints[i]->x - center_x;
      y = boundaryGaussPoints[i]->y - center_y;
      z = boundaryGaussPoints[i]->z - center_z;
      r = sqrt(x*x +y*y + z*z);
//	  printf("radius %e \n",r);
	  if (n(0)*x +n(1)*y + n(2)*z>0)
	  {
      boundaryGaussPoints[i]->N(0) = x/r;  
      boundaryGaussPoints[i]->N(1) = y/r;  
      boundaryGaussPoints[i]->N(2) = z/r;
	  } 
	  else
	  {
      boundaryGaussPoints[i]->N(0) = -x/r;  
      boundaryGaussPoints[i]->N(1) = -y/r;  
      boundaryGaussPoints[i]->N(2) = -z/r;
	  } 
    }
}


 void DGBoundaryCell::computeCurvedNormals(mPoint &vert1, mPoint &vert2, mPoint &vert3, mVector &norm1, mVector &norm2, mVector &norm3) {
	 const double x1=vert1(0),y1=vert1(1),z1=vert1(2);
	 const double x2=vert2(0),y2=vert2(1),z2=vert2(2);
	 const double x3=vert3(0),y3=vert3(1),z3=vert3(2);
	 double s,t,norm;
	 double ptx,pty,ptz;

	 for(int i = 0; i < nbPtGauss; ++i) {
		 ptx = boundaryGaussPoints[i]->x;
		 pty = boundaryGaussPoints[i]->y;
		 ptz = boundaryGaussPoints[i]->z;

		 //use equation of a plane r=r0+s*u+t*v
		 //	in this case, r0=(x1,y,z1) u=(x2-x1,y2-y1,z2-z1) v=(x3-x1,y3-y1,z3-z1)
		 s = -1*((pty-y1)*(x3-x1)-(ptx-x1)*(y3-y1))/((y3-y1)*(x2-x1)-(y2-y1)*(x3-x1));
		 t = ((pty-y1)*(x2-x1)-(ptx-x1)*(y2-y1))/((y3-y1)*(x2-x1)-(y2-y1)*(x3-x1));

		 norm = 1./sqrt(((1-s-t)*norm1(0)+s*norm2(0)+t*norm3(0))*((1-s-t)*norm1(0)+s*norm2(0)+t*norm3(0))+((1-s-t)*norm1(1)+s*norm2(1)+t*norm3(1))*((1-s-t)*norm1(1)+s*norm2(1)+t*norm3(1))+((1-s-t)*norm1(2)+s*norm2(2)+t*norm3(2))*((1-s-t)*norm1(2)+s*norm2(2)+t*norm3(2)));

		 boundaryGaussPoints[i]->N(0) = norm*(1-s-t)*norm1(0)+s*norm2(0)+t*norm3(0);
		 boundaryGaussPoints[i]->N(1) = norm*(1-s-t)*norm1(1)+s*norm2(1)+t*norm3(1);
		 boundaryGaussPoints[i]->N(2) = norm*(1-s-t)*norm1(2)+s*norm2(2)+t*norm3(2);
	 }
 }

double DGBoundaryCell::computeRadius(mPoint &p1, mPoint &p2,mPoint &p3)
{
  double x,y;
  double c1 = p1(0)*p1(0) + p1(1)*p1(1);
  double c2 = c1 - p2(0)*p2(0) - p2(1)*p2(1);
  double c3 = c1 - p3(0)*p3(0) - p3(1)*p3(1);
  double d1 = p1(0)-p3(0);
  double d2 = p1(0)-p2(0);
  double d3 = p1(1)-p3(1);
  double d4 = p1(1)-p2(1);
  double center_x = 0.5*(c2*d3-c3*d4)/(d2*d3-d1*d4);
  double center_y = 0.5*(c2*d1-c3*d2)/(d4*d1-d3*d2);
  x = p1(0)-center_x;
  y = p1(1)-center_y;
  return sqrt(x*x+y*y);
 }

// Calculates the determinant of an nxn square matrix a
// Defined recursively using expansion by minors
double DGBoundaryCell::determinant(double **a, int n)
{
       double det=0;
       double **m=NULL;
       
       if (n==1)
       {
            det=a[0][0];   
       }
       else if(n==2)
       {
            det=a[0][0]*a[1][1]-a[1][0]*a[0][1];  
       }
       else if(n>2)
       {
            for(int k=0;k<n;k++)
            {
                    // Allocate space for minor
                    m=(double**)malloc((n-1)*sizeof(double *));
                    for(int i=0;i<n-1;i++)
                    {
                             // Allocate space for each column in the minor
                             m[i]=(double*)malloc((n-1)*sizeof(double));
                    }
                    // Create the minor matrix
                    for(int i=1;i<n;i++)
                    {
                            int c=0;
                            for (int j=0;j<n;j++)
                            {
                                if (j!=k)
                                {
                                         m[i-1][c]=a[i][j];
                                         c++;
                                }
                            }
                    }
                    // Calculate the determinant recursively
                    det+=pow(-1.,k)*a[0][k]*determinant(m,n-1);
                    // Free space used for the minor
                    for(int i=0; i<n-1;i++)
                    {
                            free(m[i]);
                    }
                    free(m);
            }
       }
       return det;
}

// Returns the radius of the sphere formed by four points.  If the four points
// are coplaner then there are zero or infinitly many solutions so return 0.
// The centre is returned as p[4].
double DGBoundaryCell::computeSphere(mPoint p[4])
{
       double r=0;
       
       // Allocate space for the minors
       double **a=(double**)malloc(4*sizeof(double *));
       for (int i=0;i<4;i++)
       {
           a[i]=(double*)malloc(4*sizeof(double));
       }
       
       // Find determinant M11
       for (int i=0;i<4;i++)
       {
           a[i][0]=p[i](0);
           a[i][1]=p[i](1);
           a[i][2]=p[i](2);
           a[i][3]=1;
       }
       double m11=determinant(a,4);
       
       // Find determinant M12
       for (int i=0;i<4;i++)
       {
           a[i][0]=p[i](0)*p[i](0)+p[i](1)*p[i](1)+p[i](2)*p[i](2);
           a[i][1]=p[i](1);
           a[i][2]=p[i](2);
           a[i][3]=1;
       }
       double m12=determinant(a,4);
       
       // Find determinant M13
       for (int i=0;i<4;i++)
       {
           a[i][0]=p[i](0);
           a[i][1]=p[i](0)*p[i](0)+p[i](1)*p[i](1)+p[i](2)*p[i](2);
           a[i][2]=p[i](2);
           a[i][3]=1;
       }
       double m13=determinant(a,4);
       
       // Find determinant M14
       for (int i=0;i<4;i++)
       {
           a[i][0]=p[i](0);
           a[i][1]=p[i](1);
           a[i][2]=p[i](0)*p[i](0)+p[i](1)*p[i](1)+p[i](2)*p[i](2);
           a[i][3]=1;
       }
       double m14=determinant(a,4);
       
       // Find determinant M15
       for (int i=0;i<4;i++)
       {
           a[i][0]=p[i](0)*p[i](0)+p[i](1)*p[i](1)+p[i](2)*p[i](2);
           a[i][1]=p[i](0);
           a[i][2]=p[i](1);
           a[i][3]=p[i](2);
       }
       double m15=determinant(a,4);
      // Calculate the centre and the radius
       // If M11=0 then points define one or infinitly many spheres - take r=0
       if (m11>=0.00000001||m11<=-0.00000001)
       {
           p[4](0)=0.5*m12/m11;
           p[4](1)=0.5*m13/m11;
           p[4](2)=0.5*m14/m11;
           r=sqrt(p[4](0)*p[4](0)+p[4](1)*p[4](1)+p[4](2)*p[4](2)-m15/m11);
       }
       else
       {
           p[4](0)=0.0;
           p[4](1)=0.0;
           p[4](2)=0.0;
       }
       
       // Free space used for minors
       for (int i=0;i<4;i++)
       {
           free(a[i]);
       }
       free(a);
       
       return r;
} 

// Compute of maximum eigenvalue of the Jacobian of the Flux
// through the boundary
double DGBoundaryCell::computeMaxEVal()
{
	if ( right)return left->theConservationLaw->maximumEigenValue(n,left->theMean,right->theMean);
	else return left->theConservationLaw->maximumEigenValue(n,left->theMean,left->theMean);
}

// Compute the Jacabian Matrix of the Flux through the boundary
void DGBoundaryCell::computeJacobian(int Side, double* Up0, double** Jac)
{
    mVector nn(n);
    if (Side==0) nn*=(-1.);
	left->theConservationLaw->computeJacobian(nn,Up0,Jac);
}


// Computes change in Flux along boundary cell.
void DGBoundaryCell::computeDeltaF(int Side, double* Up0, double* deltaUp0,double* flux)
{       
   mPoint point;
   mVector flux1[MaxNbEqn];
   mVector flux2[MaxNbEqn];
   double Q1[MaxNbEqn];
   
   for (int i=0; i<cSize; i++) Q1[i]=Up0[i]+deltaUp0[i];
   left->theConservationLaw->Fi(point,Q1,flux1);
   left->theConservationLaw->Fi(point,Up0,flux2);
  
   // Calculate the DELTA FLUX
    for (int i=0; i<5; i++)
    {
        flux[i]=(flux1[i](0)-flux2[i](0))*n(0)+(flux1[i](1)-flux2[i](1))*n(1)+(flux1[i](2)-flux2[i](2))*n(2);
        if (Side==0) flux[i]=flux[i]*(-1);
    }
}

void DGCell::ZeroError ()
{
  error = 0;
}

void DGCell::ZeroJump ()
{
  jump = 0;
  jumpError=0;
}

/*********************** U T I L S *****************************/

// compute  of variables and their derivatives on the cell
void DGCell::computeMean ()
{
  
 if(theFunctionSpace->isOrthogonal())
 {
	 double firstBasisFunction = pt(0).fcts[0];
	for(int j=0;j<cSize;++j)theMean[j] = theFieldsCoefficients->get(j,0)* firstBasisFunction;
	}
 else 
 {
	 int i,j; 
  for(i=0;i<cSize;++i)theMean[i] = 0.0;
  double vol = 0.0;
  double val[MaxNbEqn];
    for(i=0;i<nbPtGauss;++i)
    {
		volumeGaussPoint  &pg = pt(i);
      interpolate(pg.fcts, val);
      vol += pg.JacTimesWeight;
      for(j=0;j<cSize;++j) theMean[j] +=pg.JacTimesWeight*val[j];
    }
  for(j=0;j<cSize;++j)theMean[j] /= vol;
 }
 }

void DGCell::adaptTimeStep (double CFLMAX, double &DT)
{
  // CFL = a DT / DX
  // CFL < CFLMAX -> DT < DX CFLMAX / a
  double a = computeMAXEV();
  double deltat = cellSize*CFLMAX / (a * (2.*fOrder + 1.0));

  if (DT > deltat)DT = deltat;
}

double DGCell::computeMAXEV()
{
   switch(theMeshEntity->getType())
	{
   case mEntity::TRI:        /*Triangle*/  //Note: should higher orders be the same?
		{
  volumeGaussPoint &pg = pt(0);
  mPoint p(pg.x,pg.y,pg.z);  
  double a = theConservationLaw->maximumEigenValue(theMean,p);
  double val[MaxNbEqn];
  double u,v,w,b;
  u=0;v=0;w=0;
  interpolate(u,v,w,val);
  b  = theConservationLaw->maximumEigenValue(val,p);
  a = (b >a ? b : a);
  u=1.;v=0;w=0;
  interpolate(u,v,w,val);
  b  = theConservationLaw->maximumEigenValue(val,p);
  a = (b >a ? b : a);
   u=0;v=1.;w=0;
  interpolate(u,v,w,val);
  b  = theConservationLaw->maximumEigenValue(val,p);
  a = (b >a ? b : a);
  return a;
		}
		break;
   case mEntity::QUAD:			/*Square*/
		{
			mPoint p((pMax(0)+pMin(0))*0.5,(pMax(1)+pMin(1))*0.5);
			return theConservationLaw->maxEigenValueCube(pMax(0)-pMin(0),pMax(1)-pMin(1),0.0,theMean,p);
		}
		break;
   case mEntity::TET:
		{
			mPoint p((pMax(0)+pMin(0))*0.5,(pMax(1)+pMin(1))*0.5,(pMax(2)+pMin(2))*0.5);
		return theConservationLaw->maximumEigenValue(theMean,p);
		}
		break;
	default:
		printf("computeMAXEV: this type has not been coded yet\n");
		return 0.0;
   }
}

void DGCell::computeCellSize()
{	
  switch(theMeshEntity->getType())
    {
  case mEntity::TRI:        /*Triangle*/
      {									 
		//Note: This is repeated for higher order. Ie: Just reduced to a linear triangle, since there isn't necessarily a
	    //inscribed circle for a higher order triangle.  However, perhaps there is a better solution???
	mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
	mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
	mPoint p3 = ((Vertex*)theMeshEntity->get(0,2))->point();
	mVector A(p1,p2);
	mVector B(p2,p3);
	mVector C(p3,p1);									//Note: Consider changing Heron's formula
	double a = sqrt(A*A);								//to more numerically stable version, perhaps?
	double b = sqrt(B*B);
	double c = sqrt(C*C);
	double halfP = 0.5 * (a+b+c); 
	cellSize = sqrt(halfP*(halfP-a)*(halfP-b)*(halfP-c))/halfP*2.; //Heron's Formula modified to give half
																	//the radius of inscribed circle
      }
      break;
  case mEntity::QUAD:         /*Square*/
      {
	mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
	mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
	mVector a(p1,p2);
	cellSize = a.L2norm(); 
      }
	  break;
  case mEntity::TET:        /*Tet*/
      {
	mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
	mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
	mPoint p3 = ((Vertex*)theMeshEntity->get(0,2))->point();
	mPoint p4 = ((Vertex*)theMeshEntity->get(0,3))->point();
	mVector A(p1,p2);
	mVector B(p2,p3);
	mVector C(p3,p1);
	mVector D(p1,p4);
	mVector E(p2,p4);
	mVector F(p3,p4);
	cellSize = A.L2norm();
	if (cellSize > B.L2norm() )cellSize = B.L2norm();
	if (cellSize > C.L2norm() )cellSize = C.L2norm();
	if (cellSize > D.L2norm() )cellSize = D.L2norm();
	if (cellSize > E.L2norm() )cellSize = E.L2norm();
	if (cellSize > F.L2norm() )cellSize = F.L2norm();
	break;
      }
    default : printf("size function for this case is not coded yet\n");
      cellSize=0.0;
    }
}

void DGCell::computeVolume()
{	
  switch(theMeshEntity->getType())
    {
  case mEntity::TRI:        /*Triangle*/
	  if (theMeshEntity->getgeomOrder() == 1)  
		  volume = 0.5*detJac;
	  else  //higher order. uses quadrature.
	  {
		  volume=0;
		  for (int k=0;k<nbPtGauss;k++) 
		  {
			volume+=volumeGaussPoints[k].JacTimesWeight;
		  }
	  }
      //volumeRatio = 1./pt(0).detJac;
      break;
  case mEntity::QUAD:         /*Square*/  
      {

	//volumeRatio = 4./volume;
	if (theMeshEntity->getgeomOrder() == 1) 
		{
		mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
		mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
		mPoint p3 = ((Vertex*)theMeshEntity->get(0,2))->point();
		mVector a(p1,p2);
		mVector b(p2,p3);
		volume = a.L2norm()*b.L2norm();
		}
	else  //higher order. uses quadrature.
		{
		  volume=0;
		  for (int k=0;k<nbPtGauss;k++) 
		  {
			volume+=volumeGaussPoints[k].JacTimesWeight;
		  }
		}
      }
    break;
  case mEntity::TET:			/*Tetrahedron*/
      {
		  if (theMeshEntity->getgeomOrder() == 1)
		  {
			volume = detJac/6.;
			// volumeRatio = 1./pt(0).detJac;
		  }
		  else
		  {
			volume=0;
			for (int k=0;k<nbPtGauss;k++) 
			{
			   volume+=volumeGaussPoints[k].JacTimesWeight;
			}
		  }
      }
      break;
    default : printf("volume function for this case is not coded yet\n");
    }
}

void DGCell::computePerimeter()
{								
								
  switch(theMeshEntity->getType())
    {
  case mEntity::TRI:        /*Triangle*/
      {
    if (theMeshEntity->getgeomOrder() == 1)
		{
			mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
			mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
			mPoint p3 = ((Vertex*)theMeshEntity->get(0,2))->point();
			mVector a(p1,p2);
			mVector b(p2,p3);
			mVector c(p3,p1);
			perimeter = a.L2norm()+b.L2norm()+c.L2norm();
		}
	else      //uses arclength integration for higher order elements
		{										
			perimeter=0;									
			int NbSides = theMeshEntity->getNbTemplates(1);
			for(int k=0;k<NbSides;k++)	//iterate through edges and sum up arclengths
			  perimeter+=arclength(theMeshEntity->get(1,k), theGaussIntegrator);
		}
	  }
    break;
  case mEntity::QUAD:        /*Quad*/       
      {
    if (theMeshEntity->getgeomOrder() == 1)
		{
		mPoint p1 = ((Vertex*)theMeshEntity->get(0,0))->point();
		mPoint p2 = ((Vertex*)theMeshEntity->get(0,1))->point();
		mPoint p3 = ((Vertex*)theMeshEntity->get(0,2))->point();
		mPoint p4 = ((Vertex*)theMeshEntity->get(0,3))->point();
		mVector a(p1,p2);
		mVector b(p2,p3);
		mVector c(p3,p4);
		mVector d(p4,p1);
		perimeter = a.L2norm()+b.L2norm()+c.L2norm()+d.L2norm();
		}
	else      //uses arclength integration for higher order elements
		{										
			perimeter=0;									
			int NbSides = theMeshEntity->getNbTemplates(1);
			for(int k=0;k<NbSides;k++)	//iterate through edges and sum up arclengths
			  perimeter+=arclength(theMeshEntity->get(1,k), theGaussIntegrator);
		}
	  }
      break;
    default : //printf("perimeter function for this case is not coded yet\n");
		perimeter=0.0;
    }
}
//Finds the arclength of an nth-order edge by use of quadrature.
double DGCell::arclength(mEntity* theEdge, GaussIntegrator* theGaussIntegrator) const
{
	double arclength=0;		
	double u,v,w,weight;        //used to represent an integration point.
	double dxdu,dydu,dzdu;	    //represent partial derivatives of edge mapping.

	vector<Vertex*> verts;		  //used to store vertices composing the edge
	theEdge->getVertices(verts);
	mVector grads;                //holds the gradients of geometrical shape functions
	MeshMapping m(theEdge);       //a mapping for the edge is set up.

	int nbPt = theGaussIntegrator->nbIntegrationPoints(theEdge,order);

	for(int i=0;i<nbPt;i++)         //iterates through the integration points and
	{								//adds to the quadrature sum which represents arclength.
		dxdu=0;
		dydu=0;
		dzdu=0;
		theGaussIntegrator->iPoint(theEdge,i,order,u,v,w,weight);
		for(int j=0;j<((theEdge->getgeomOrder())+1);j++)
		{
			m.GradGeomShapeFunction(j,u,v,w,grads);
			dxdu+=(verts[j]->point()(0))*grads[0];
			dydu+=(verts[j]->point()(1))*grads[0];
			dzdu+=(verts[j]->point()(2))*grads[0];
		}

		arclength+=weight*sqrt(dxdu*dxdu+dydu*dydu+dzdu*dzdu);
	}
	return arclength;
}

void DGCell::L2Proj (FieldEvaluator *f, double *proj)
{
  double val[MaxNbEqn];
  int i, j, k;
  const int s = cSize * fSize;

  for( i=0;i<s;++i) proj[i] = 0.0;
  
  for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint &pg=pt(i);
      f->eval(mPoint(pg.x,pg.y,pg.z),0.0,val);
     
      for( j=0;j<cSize;++j)
	for( k=0;k<fSize;++k)
	  proj[k+fSize*j] += val[j] * pg.fcts[k] * pg.JacTimesWeight;
    }
}

void DGCell::L2Proj (double* dU)
{
  double val[MaxNbEqn];
  int i, j, k;
  const int dof = cSize * fSize;

  for( i=0;i<dof;++i) theRightHandSide[i] = 0.0;
  
  for(i=0;i<nbPtGauss;++i)
    {
      volumeGaussPoint &pg=pt(i);
      interpolate(pg.fcts,val);
     for( j=0;j<cSize;++j)	val[j] += dU[j];
      for( j=0;j<cSize;++j)
	for( k=0;k<fSize;++k)
	  theRightHandSide[k+fSize*j] += val[j] * pg.fcts[k] * pg.JacTimesWeight;
    }

  if (theFunctionSpace->isOrthogonal())
    {
	 double inv_Jac = 1./detJac;
	 double* coeff = theFieldsCoefficients->get();
      for(i=0;i<dof;++i) coeff[i] = theRightHandSide[i]*inv_Jac;
    }
  else printf("Multigrid projection is not done for not orthogonal spaces\n");
}


void DGCell::L2ProjInitial (FieldEvaluator *f)
{
  int i;
  L2Proj(f,theRightHandSide);		 
  if (theFunctionSpace->isOrthogonal())
    {
	 double inv_Jac = 1./detJac;
	 int dof = cSize*fSize;
	 double* coeff = theFieldsCoefficients->get();
      for(i=0;i<dof;++i)
	  coeff[i] = theRightHandSide[i]*inv_Jac;
    }
  else
    {

		for(int k = 0;k<cSize;++k)		
	for(int i = 0;i<fSize;++i)
	  {
	    double dxi = 0.0;
		for(int j = 0;j<fSize;++j)
	      dxi += theInvertMassMatrix[i][j] * theRightHandSide[j+fSize*k];
	    
		theFieldsCoefficients->get(k,i) = dxi;
	  }
    }
}




void DGCell::getBox(double &XminBox,double &YminBox,double &ZminBox, 
		    double &XMaxBox,double &YMaxBox,double &ZMaxBox)
{
  XminBox = pMin(0);
  YminBox = pMin(1);
  ZminBox = pMin(2);
  XMaxBox = pMax(0);
  YMaxBox = pMax(1);
  ZMaxBox = pMax(2);
}

bool DGCell::inCell (const mPoint &p, double *val)
{
  /*
  printf("%f %f %f %f %f %f %f %f %f\n",p(0),p(1),p(2),pMin(0),pMin(1),pMin(2),
	 pMax(0),pMax(1),pMax(2));
  */
  if(p(0) > pMax(0) ||
     p(1) > pMax(1) ||
     p(2) > pMax(2) ||
     p(0) < pMin(0) ||
     p(1) < pMin(1) ||
     p(2) < pMin(2))return false;

  //  printf("in box...\n");
  
  double u,v,w;
  theMapping->invert(p(0),p(1),p(2),u,v,w);
  //  printf("%f %f %f\n",u,v,w);
  if(!theMapping->inReferenceElement(u,v,w))return false;
  //  printf("ok\n");
  interpolate(u,v,w,val);
  return true;
}

void DGCell::write(ostream &o)
{
  o << theFunctionSpace->order() << " ";
  o.precision(16);
  for(int i=0;i<cSize;++i)
    for(int j=0;j<fSize;++j)
      o << theFieldsCoefficients->get(i,j) << " ";
  
}

void DGCell::read(istream &is)
{
  for(int i=0;i<cSize;++i)
    {
      for(int j=0;j<fSize;++j) 
        is >> theFieldsCoefficients->get(i,j);
	}  
}

void DGCell::adaptOrder(int newOrder)
{
  int i,j;
  delete theFunctionSpace;
  switch(theMeshEntity->getType())
    {
    case mEntity::TRI  :
		if (theMeshEntity->getgeomOrder() == 1) theFunctionSpace = new OrthogonalTriangleFunctionSpace(newOrder);
		else                                    theFunctionSpace = new TriangleFunctionSpace(newOrder);
		break;
    case mEntity::QUAD : theFunctionSpace = new QuadFunctionSpace(newOrder); break;
    case mEntity::TET  :
		if (theMeshEntity->getgeomOrder() == 1) theFunctionSpace = new OrthogonalTetFunctionSpace(newOrder);
		else                                    theFunctionSpace = new TetFunctionSpace(newOrder);
		break;
    case mEntity::HEX  : theFunctionSpace = new HexFunctionSpace(newOrder); break;
    }
  int oldFSize = fSize;
  fSize                 = theFunctionSpace->size();
  fOrder                = theFunctionSpace->order();
  order = 2 * fOrder + theMapping->detJacOrder();

  DG_Dofs *oldFieldsCoefficients;
  oldFieldsCoefficients = new DG_Dofs(2,cSize,oldFSize);
  oldFieldsCoefficients->operator =(theFieldsCoefficients);
  delete theFieldsCoefficients;
  theFieldsCoefficients = new DG_Dofs(2,cSize,fSize);
  if (fSize>oldFSize)
    {
      for(i=0;i<cSize;++i)
	{
    for(j=0;j<oldFSize;++j)
	{
      theFieldsCoefficients->get(i,j) = oldFieldsCoefficients->get(i,j);
	  //printf("%e \n", oldFieldsCoefficients->get(i,j));
	}
    for(j=oldFSize;j<fSize;++j) 
      theFieldsCoefficients->get(i,j) = 0.;
	}
    }
else  for(i=0;i<cSize;++i)
    for(j=0;j<fSize;++j)
	{
      theFieldsCoefficients->get(i,j) = oldFieldsCoefficients->get(i,j);
	  //printf("%e \n", oldFieldsCoefficients->get(i,j));
	}
	
	delete [] theRightHandSide;
	delete oldFieldsCoefficients;
  theRightHandSide = new double [fSize*cSize];
  ZeroRHS();
  volumeGaussPoints.clear();	
  init();
  int Sizenm1 = theMeshEntity->size(theMeshEntity->getLevel()-1);
 for(j=0;j<Sizenm1;j++)
    {
      mEntity *b = theMeshEntity->get(1,j);
      DGBoundaryCell *bc = (DGBoundaryCell*)b->getCell();
      bc->init();
    }
}


/*-------- Coeffs Management ------------*/

DG_Dofs::DG_Dofs (short n, short NbUnknowns, short FunctionSpaceSize)
  : Nb (n), NU (NbUnknowns), FS(FunctionSpaceSize)
{
  theFieldsCoefficients = new double* [Nb];
  for(int i=0;i<Nb;++i) theFieldsCoefficients[i] = new double [NbUnknowns * FunctionSpaceSize];
}

DG_Dofs* DG_Dofs::operator= (const DG_Dofs *other)
{
  Nb=other->Nb;
  NU=other->NU;
  FS=other->FS;
  //theFieldsCoefficients = new double *[Nb];
  //for(int i=0;i<Nb;++i) theFieldsCoefficients[i] = new double [NU * FS];
  for(int i=0;i<NU;++i) 
	  for (int j=0;j<FS;++j)
	  {
		  theFieldsCoefficients[0][i*FS+j] = other->theFieldsCoefficients[0][i*FS+j];
		  //printf(" %d %e %e \n",FS,theFieldsCoefficients[0][i*FS+j],other->theFieldsCoefficients[0][i*FS+j]);
	  }
  return this;
}



DG_Dofs::~DG_Dofs ()
{
  for(int i=0;i<Nb;++i) delete [] theFieldsCoefficients[i];
  delete [] theFieldsCoefficients;
}

void DG_Dofs::copy()   //copies solution coefficients kept in the fist vector to the second vector in DG_Dofs
{
const int N = NU*FS;
for (int i=0;i<N;++i)
theFieldsCoefficients[1][i] = theFieldsCoefficients[0][i];  
}

void DG_Dofs::swap()   //copies solution coefficients kept in the fist vector to the second vector in DG_Dofs
{
  const int N = NU*FS;
  double tmp;
  for (int i=0;i<N;++i)
    {
      tmp = theFieldsCoefficients[0][i] ; 
      theFieldsCoefficients[0][i] = theFieldsCoefficients[1][i];
      theFieldsCoefficients[1][i] = tmp; 
    }
}

void DG_Dofs::copyBack() //copies solution coeff from the second vector to the first.
{
const int N = NU*FS;
for (int i=0;i<N;++i)
theFieldsCoefficients[0][i] = theFieldsCoefficients[1][i];  
}

void DG_Dofs::copyLimitedCoeffs(int k) //copies solution coeff from the second vector to the first after limiting
{
const int N = NU*FS;
 for (int q=0;q<NU;q++)
   for (int i=q*FS+k;i<FS*(q+1);++i)
     theFieldsCoefficients[0][i] = theFieldsCoefficients[1][i];  
}

void DG_Dofs:: eval (vector<double> &ff, double *field, double t) 
{
  //  double *b = 0;
  double *a,A,F;
 
  //double frac = (t-t0)/(t1-t0);
  //  printf("%f %f %f %f %d\n",t,t0,t1,frac,frac != 0.0);
  
  for(int i=0;i<NU;++i)
    {
      a = & theFieldsCoefficients[0][i*FS];
      //  if(frac != 0.0) 
      //	{
      //  b = &theFieldsCoefficients[1][i*FS];
//}
      field[i] = 0.0; 
      for(int j=0;j<FS;++j)
	{
	  //  if(b)field[i] += ff[j] * ((1.-frac) * a[j] + (frac)*b[j]);
	  // else 
	  F = ff[j]; A = a[j];
	  //field[i] += ff[j] * a[j];
	  field[i] += F * A;
	 	}
    }      
}

void DG_Dofs:: eval (vector<double> &ff, double *field)
{
  double *a, A,F;
  //  double *k = &ff.front();
  for(int i=0;i<NU;++i)
    {
      a = & theFieldsCoefficients[0][i*FS];
      field[i] = 0.0;
      for(int j=0;j<FS;++j)
	{
	  F=ff[j]; A=a[j];
	  //field[i] += ff[j] * a[j];
	  field[i] += F * A;
	  //	  printf("%f ",a[j]);
	}
      //      printf("\n");
    }      
}
