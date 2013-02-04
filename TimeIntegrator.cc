#include "TimeIntegrator.h" 
#include "DGLimiter.h"
#include "mDGMesh.h"
#include "mMesh.h"
#include "mEntity.h"
#include "DGCell.h"
#include "FunctionSpace.h"
#include "ConservationLaw.h"
#include "FieldEvaluator.h"
#include "Mapping.h"
#include "parallel.h"
#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <time.h>

extern double ** allocateMatrix(int n);
extern void invmat (double **, double **, int);
extern void freeMatrix( double **v);

extern unsigned int numCPU;
extern unsigned int numThreads;
extern vector<infoWrapper> threadInfo[4];

extern double timeComputeVolumeCalls;
extern int numComputeVolumeCalls;

extern double timeComputeBoundaryCalls;
extern int numComputeBoundaryCalls;

timespec timer1, timer2;

TimeIntegrator::TimeIntegrator(mDGMesh *m,DGLimiter *l,int NbFields, int DOF):theMesh(m),theLimiter(l),
								     cSize(NbFields),dof(DOF)
{
 if(theMesh->size(3)) n = 3;
  else if(theMesh->size(2)) n = 2;
  else n = 1;
 }

RungeKutta4::RungeKutta4(mDGMesh *m,DGLimiter *l,int NbFields,int DOF):TimeIntegrator(m,l,NbFields,DOF)
{
  soln = new double [dof];
  soln_begin = soln;  
  solnsum = new double [dof];
  solnsum_begin=solnsum;
}

RungeKuttaTVD2::RungeKuttaTVD2(mDGMesh *m,DGLimiter *l,int NbFields,int DOF):TimeIntegrator(m,l,NbFields,DOF)
{
  soln = new double [dof];
  soln_begin = soln;  
}

RungeKuttaTVD2adaptive::RungeKuttaTVD2adaptive(mDGMesh *m,DGLimiter *l,int NbFields,int DOF):TimeIntegrator(m,l,NbFields,DOF)
{
  soln = new double [dof];
  soln_begin = soln;  
}

RungeKuttaTVD2adaptive::~RungeKuttaTVD2adaptive()
{
  soln=soln_begin;
  delete [] soln; 
}

RungeKuttaTVD3::RungeKuttaTVD3(mDGMesh *m,DGLimiter *l,int NbFields,int DOF):TimeIntegrator(m,l,NbFields,DOF)
{
  soln = new double [dof];
  soln_begin = soln;  
}

ForwardEuler::ForwardEuler(mDGMesh *m,DGLimiter *l,int NbFields,int DOF):TimeIntegrator(m,l,NbFields,DOF)
{}

Multigrid::Multigrid(mDGMesh *m,DGLimiter *l,int NbFields,int DOF):TimeIntegrator(m,l,NbFields,DOF)
{
	DMatrix=allocateMatrix(cSize);
	Jac=allocateMatrix(cSize);   
	rk2 = new RungeKuttaTVD2(m,l,NbFields, DOF);
	rk3 = new RungeKuttaTVD3(m,l,NbFields, DOF);
	rk4 = new RungeKutta4(m,l,NbFields, DOF);
}

double RungeKutta4::advanceInTime(double t, double dt)
  /*  Classical Fourth Order Runge-Kutta integration.*/
{
  double  dt_over_detJac,dx,dof;
  int i;
  
  soln=soln_begin;
  solnsum=solnsum_begin;
  mMesh::iter it;
  mMesh::iter mesh_begin = theMesh->begin(n);
  mMesh::iter mesh_end=theMesh->end(n);
  
  assembleVolume(t);
  assembleBoundary(t);
  
 
  /********Limiting?************************/
  
  /***********  Solve for K1   ***************************/
        
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      dof = cell->fSize * cSize;
      const double *RHS = cell->theRightHandSide;
      double* coeff = cell->theFieldsCoefficients->get();
      
      dt = cell->getTimeStep();
      if(cell->details->theFunctionSpace->isOrthogonal())
	{
	  dt_over_detJac = dt/cell->getDetJac();
	  for(i = 0;i<dof;i++)
	    {
	      *soln++ = coeff[i]; /*Copy Initial values****/   
	      *solnsum = RHS[i]*dt_over_detJac;
	      coeff[i]+= 0.5*RHS[i]*dt_over_detJac;
	      solnsum++;
	    }
	}
      else
	{
	  int fSize = cell->fSize;
	  for(int k = 0;k<cSize;k++)
	    for(int i = 0;i<fSize;i++)
	      {
		dx = 0.0;
		*soln++ = cell->theFieldsCoefficients->get(k,i); /*Copy Initial values****/   
		for(int j = 0;j<fSize;j++)
		  dx += cell->details->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSize*k];
		*(solnsum++) = dx*dt;
		cell->theFieldsCoefficients ->get(k,i)+= dx*dt*0.5;
	      }
	}
      cell->computeMean();
      cell->ZeroRHS();
    }
  
  if (theLimiter) limit(t);
  
  /*  Solve for K2   */
  
  assembleVolume(t+dt*0.5);
  assembleBoundary(t+dt*0.5);
  
  solnsum = solnsum_begin;//Reset pointers
  soln = soln_begin;
  
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      dof = cell->fSize * cSize;
      const double *RHS = cell->theRightHandSide;
      double* coeff = cell->theFieldsCoefficients->get();
      dt = cell->getTimeStep();
      if(cell->details->theFunctionSpace->isOrthogonal())
	{
	  dt_over_detJac = dt /cell->getDetJac();
	   for(i = 0;i<dof;i++)
	     {
	       (*solnsum++) += 2.*dt_over_detJac*RHS[i];
	       coeff[i]      = (*soln++)+0.5*dt_over_detJac*RHS[i];
	    }
	}
      else
	{
	  int fSize = cell->fSize;
	  for(int k = 0;k<cSize;k++)
	    for(int i = 0;i<fSize;i++)
	    {
	      double dx = 0.0;
	      for(int j = 0;j<fSize;j++)
		dx += cell->details->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSize*k];
	      *(solnsum++) += 2.0*dx*dt;
	      cell->theFieldsCoefficients ->get(k,i)= *(soln++) + dx*dt*0.5;
	    }
	}
      cell->computeMean();
      cell->ZeroRHS();
    }
  
  if (theLimiter) limit(t);
  
  /*  Solve for K3   */
  assembleVolume(t+dt*0.5);
  assembleBoundary(t+dt*0.5);  
  
  solnsum = solnsum_begin;//Reset pointers
  soln = soln_begin;
  
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      dof = cell->fSize * cSize;
      const double *RHS = cell->theRightHandSide;
      double* coeff = cell->theFieldsCoefficients->get();
      dt = cell->getTimeStep();
      if(cell->details->theFunctionSpace->isOrthogonal())	
	{
	  dt_over_detJac = dt /cell->getDetJac();
	  for(i = 0;i<dof;i++)
	    {
	      (*solnsum++) += 2.*dt_over_detJac*RHS[i];
	      coeff[i]=(*soln++)+dt_over_detJac*RHS[i];
	    }
	}
      else {
	int fSize = cell->fSize;
	for(int k = 0;k<cSize;k++)
	  for(int i = 0;i<fSize;i++){
	    dx = 0.0;
	    for(int j = 0;j<fSize;j++)
	      dx += cell->details->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSize*k];
	    *(solnsum++) += 2.0*dx*dt;
	    cell->theFieldsCoefficients->get(k,i) =*(soln++) + dx*dt;
	  }
      }
      cell->computeMean();
      cell->ZeroRHS();
    } 
  
  if (theLimiter) limit(t);
  
  /*  Solve for K4   */ 
  assembleVolume(t+dt);
  assembleBoundary(t+dt);  
  
  solnsum = solnsum_begin;//Reset pointers
  soln = soln_begin;
  
  double resid=0.0;

  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      dof = cell->fSize * cSize;
      const double *RHS = cell->theRightHandSide;
      double* coeff = cell->theFieldsCoefficients->get();
      dt = cell->getTimeStep();
      
      if(cell->details->theFunctionSpace->isOrthogonal())
	{
	  dt_over_detJac = dt /cell->getDetJac();
	  for(i = 0;i<dof;i++)
	    {
	      (*solnsum) += dt_over_detJac*RHS[i];
	      coeff[i] = (*soln++) +(*solnsum)/6.0;
	      resid+=(*solnsum)*(*solnsum)/36.0;
	      solnsum++;
	    }
	}
      else {
	int fSize = cell->fSize;
	for(int k = 0;k<cSize;k++)
	  for(int i = 0;i<fSize;i++) {
	    dx = 0.0;
	    for(int j = 0;j<fSize;j++)
	      dx += cell->details->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSize*k];
	    *solnsum += dx*dt;    
	    cell->theFieldsCoefficients->get(k,i) = (*soln++) +(*solnsum)/6.0;
	    resid+=(*solnsum)*(*solnsum)/36.0;
	    solnsum++;
	  }
      }  
      cell->computeMean();
      cell->ZeroRHS();
    }
  
  if (theLimiter) limit(t);
  
  return resid;
}

/***************************************************************/

double RungeKuttaTVD2::advanceInTime(double t, double dt)
/*  TVD Second  Order Runge-Kutta integration.*/
{
  double  dt_over_detJac,dof;
  int i,j,k,fSizek;
  double dx;
  soln=soln_begin; 
  mMesh::iter it;
  const mMesh::iter mesh_begin = theMesh->begin(n);
  const mMesh::iter mesh_end=theMesh->end(n);
 //  computeL2ProjInCutCells();
  assembleVolume(t);
  assembleBoundary(t);  
  
  /**************  SOLVE for K1  **************/
  
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m    = (*it);
      DGCell  *cell = (DGCell*)m->getCell();
      dof = cell->fSize * 
		  cSize;
      const double *RHS = cell->theRightHandSide;
      double* coeff = cell->theFieldsCoefficients->get();
      /* Copy initial values*/
      
      for(i = 0;i<dof;i++) *(soln++) = coeff[i];
      
      dt = cell->getTimeStep();
      if(cell->details->theFunctionSpace->isOrthogonal())
	{
	  dt_over_detJac = dt/cell->getDetJac();
	  for(i = 0;i<dof;i++) coeff[i] += RHS[i]*dt_over_detJac;
	}
      else
	{
	  int fSize = cell->fSize;
	  for(k = 0;k<cSize;k++)
	    {
	      int fSizek = fSize*k;
	    for(i = 0;i<fSize;i++)
	      {
		dx = 0.0;
		for(j = 0;j<fSize;j++){
		  dx += cell->details->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSizek];
		}
		cell->theFieldsCoefficients ->get(k,i)+= dx*dt;
	      }
	    }
	}
      cell->computeMean();
      cell->ZeroRHS();
    }
  
  if (theLimiter) limit(t);
  /**************  SOLVE for K2   ************/
  //computeL2ProjInCutCells();
  assembleVolume(t+dt);
  assembleBoundary(t+dt);  
  
  /*  Reset c values */
  double resid=0.0,tmp;
  soln=soln_begin;
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      dof = cell->fSize*cSize; 
      const double *RHS = cell->theRightHandSide;
      double* coeff = cell->theFieldsCoefficients->get();
      dt = cell->getTimeStep(); 
      if(cell->details->theFunctionSpace->isOrthogonal())
	{
	  dt_over_detJac = dt /cell->getDetJac();
	  for(i = 0;i<dof;i++)
	    {
	      coeff[i] = (coeff[i]+RHS[i]*dt_over_detJac+(*soln))*0.5;
	      tmp = coeff[i]-(*(soln++));
	      resid+=tmp*tmp;
	    }
	}
      else
	{
	  int fSize = cell->fSize;
	  for(k = 0;k<cSize;k++)
	    {
	      fSizek = fSize*k;
	      for(i = 0;i<fSize;i++)
		{
		  dx = 0.0;
		  for(j = 0;j<fSize;j++){
		    dx += cell->details->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSizek];
		  }
		  cell->theFieldsCoefficients->get(k,i) = ( cell->theFieldsCoefficients->get(k,i)+
							    dx*dt +  (*soln))*0.5;
		  tmp = cell->theFieldsCoefficients->get(k,i)-(*(soln++));
		  resid+=tmp*tmp;
		}
	    }
	}
      cell->computeMean();
      cell->ZeroRHS();
    }   
  
  if (theLimiter) limit(t+dt);
  //computeL2ProjInCutCells();
  return resid;
}

/*********************************************************/

double RungeKuttaTVD2adaptive::firstStage(int from, int to, double t, double DT)
{
  static int iter=1;
  int i,k,b,fSize;
  double dt;
  double resid=0.0;
  
  dt= DT/pow(2.,from);
  mMesh::iter it,bin_begin,bin_end,bin_begin_nm1,bin_end_nm1;
  if (from) assembleVolume(theMesh->beginSmallOnInterface(from-1,n),theMesh->endSmallOnInterface(from-1,n),t+dt);
  for (b=from;b<to;b++)
    {
	  if (theLimiter) limit(theMesh->beginBin(b,n),theMesh->endBin(b,n),t+dt);
	  if (theLimiter) limit(theMesh->beginLargeOnInterface(b,n),theMesh->endLargeOnInterface(b,n),t+dt);
	  if (theLimiter) limit(theMesh->beginSmallOnInterface(b,n),theMesh->endSmallOnInterface(b,n),t+dt);
      //if (b>0) computeSolutionOnBoundary(theMesh->beginInterfaceBin(b-1,n-1),theMesh->endInterfaceBin(b-1,n-1),t+dt);
      assembleVolume(theMesh->beginBin(b,n),theMesh->endBin(b,n),t+dt);
      assembleVolume(theMesh->beginLargeOnInterface(b,n),theMesh->endLargeOnInterface(b,n),t+dt);
      assembleVolume(theMesh->beginSmallOnInterface(b,n),theMesh->endSmallOnInterface(b,n),t+dt);	
      //if (b>0) assembleVolume(theMesh->beginSmallOnInterface(b-1,n),theMesh->endSmallOnInterface(b-1,n),t+dt);
      assembleBoundary(theMesh->beginBin(b,n-1),theMesh->endBin(b,n-1),t+dt); 
      assembleBoundary(theMesh->beginInterfaceBin(b,n-1),theMesh->endInterfaceBin(b,n-1),t+dt);
	  dt*=0.5;
    }
  if (from) assembleBoundary(theMesh->beginInterfaceBin(from-1,n-1),theMesh->endInterfaceBin(from-1,n-1),t+dt);
  
  
  // copy solution coefficients at t_n to soln,
  // set poitners to places in soln where inner elements
  //of bins and small and large elements on interface were copied
  
  //1st stage on all inner elements in all bins
  dt = DT/pow(2.,from); // set initial dt 
  for(b=from;b<to;b++)
    {
      bin_begin = theMesh->beginBin(b,n);
      bin_end = theMesh->endBin(b,n);
      if (iter==1) binptr[b] = soln; else soln = binptr[b];
      
      for(it = bin_begin;it != bin_end;++it)
	{
	  mEntity *m    = (*it);
	  DGCell  *cell = (DGCell*)m->getCell();
	  fSize         = cell->fSize;
	  const double *RHS = cell->theRightHandSide;
	  
	  for(k = 0;k<cSize;k++)
	    for(i = 0;i<fSize;i++)
	      {
		*(soln++) = cell->theFieldsCoefficients->get(k,i); /*Copy initial values*/
		cell->theFieldsCoefficients ->get(k,i)+=(*RHS++)/cell->getDetJac()*dt;
	      }
	  cell->computeMean();
	  cell->ZeroRHS();
	}
      dt *= 0.5;//reduce dt for the next bin
    }
  
  //1st stage on large elements on interface
  dt = DT/pow(2.,from);
  for(b=from;b<to-1;b++)
    {
      bin_begin = theMesh->beginLargeOnInterface(b,n);
      bin_end = theMesh->endLargeOnInterface(b,n);
      if (iter==1) binptrLI[b] = soln; else soln=binptrLI[b];
      dt *=0.5; // we compute 2 half steps
      for(it = bin_begin;it != bin_end;++it)
	{
	  mEntity *m    = (*it);
	  DGCell  *cell = (DGCell*)m->getCell();
	  fSize         = cell->fSize;
	  const double *RHS = cell->theRightHandSide;
	  
	  for(k = 0;k<cSize;k++)  //compute and save values at half step
	    for(i = 0;i<fSize;i++)
	      {
		*(soln++) = cell->theFieldsCoefficients->get(k,i); /*Copy initial values*/
		cell->theFieldsCoefficients ->get(k,i)+= (*RHS++)*dt/cell->getDetJac();
	      }
	  cell->theFieldsCoefficients->copy();   //save values at half step
	  
	  RHS = cell->theRightHandSide;
	  for(k = 0;k<cSize;k++)  //recompute coef  to full step values 
	    for(i = 0;i<fSize;i++)
	      cell->theFieldsCoefficients ->get(k,i)+= (*RHS++)*dt/cell->getDetJac();  
	  cell->computeMean();
	  cell->ZeroRHS();
	}
    }
  
  dt = DT/pow(2.,from); 
  if (from)
    {
      bin_begin = theMesh->beginSmallOnInterface(from-1,n);
      bin_end = theMesh->endSmallOnInterface(from-1,n);
      if (iter !=1) soln=binptrSI[from-1];
      for(it = bin_begin;it != bin_end;++it)
	{
	  mEntity *m    = (*it);
	  DGCell  *cell = (DGCell*)m->getCell();
	  fSize         = cell->fSize;
	  const double *RHS = cell->theRightHandSide;
	  
	  for(k = 0;k<cSize;k++)
	    for(i = 0;i<fSize;i++)
	      {
		*(soln++) = cell->theFieldsCoefficients->get(k,i); /*Copy initial values*/
		cell->theFieldsCoefficients ->get(k,i)+= (*RHS++)*dt/cell->getDetJac();
	      }
	  cell->computeMean();
	  cell->ZeroRHS();
	}
    } 

  for(b=from;b<to-1;b++)     // big step on small elements on interface
    {
      bin_begin = theMesh->beginSmallOnInterface(b,n);
      bin_end = theMesh->endSmallOnInterface(b,n);
      if (iter ==1) binptrSI[b] = soln; else soln=binptrSI[b];
      dt *=0.5;  //again, compute 2 half steps
      for(it = bin_begin;it != bin_end;++it)
	{
	  mEntity *m    = (*it);
	  DGCell  *cell = (DGCell*)m->getCell();
	  fSize         = cell->fSize;
	  const double *RHS = cell->theRightHandSide;
	  
	  for(k = 0;k<cSize;k++)
	    for(i = 0;i<fSize;i++)
	      {
		*(soln++) = cell->theFieldsCoefficients->get(k,i); /*Copy initial values*/
		cell->theFieldsCoefficients ->get(k,i)+= (*RHS++)*dt/cell->getDetJac();
	      }
	  cell->theFieldsCoefficients->copy();	//save values at half step in copy
	  
	  //compute values with dt as on neighboring large elements
	  RHS = cell->theRightHandSide;
	  for(k = 0;k<cSize;k++)
	    for(i = 0;i<fSize;i++)
	      cell->theFieldsCoefficients ->get(k,i)+= (*RHS++)*dt/cell->getDetJac(); 
	  
	  cell->computeMean();
	  cell->ZeroRHS();
	}
    }
  
  //if (theLimiter) limit(t);
  iter++;
  return 0.0;
}

    /************SECOND STAGE***********************/
double RungeKuttaTVD2adaptive::secondStage(int from, int to, double t, double DT)
{
  int i,k,b,fSize;
  double dt;
  double resid=0.0,tmp;
  
  mMesh::iter it,bin_begin,bin_end,bin_begin_nm1,bin_end_nm1;
  const mMesh::iter mesh_begin = theMesh->begin(n);
  const mMesh::iter mesh_end=theMesh->end(n);
  
  
  //  Completing the RK time step with a local dt
  dt = DT/pow(2.,from);
  for(b=from; b<to;b++)   //loop over all bins
    { 
      if (b>from)
	{
	  bin_begin = theMesh->beginSmallOnInterface(b-1,n); // set coef of small cells on interface to half step
	  bin_end = theMesh->endSmallOnInterface(b-1,n);
	  for(it = bin_begin;it != bin_end;++it)
	    {
	      mEntity *m    = (*it);
	      DGCell  *cell = (DGCell*)m->getCell();
	      cell->ZeroRHS();
	      cell->theFieldsCoefficients->copyBack();
	    }
	  //if (1) return 0.;
	}
      if (theLimiter) limit(theMesh->beginBin(b,n),theMesh->endBin(b,n),t+dt);
	  if (theLimiter) limit(theMesh->beginLargeOnInterface(b,n),theMesh->endLargeOnInterface(b,n),t+dt);
	
	  bin_begin = theMesh->beginBin(b,n);
      bin_end = theMesh->endBin(b,n);
      bin_begin_nm1 = theMesh->beginBin(b,n-1);
      bin_end_nm1 = theMesh->endBin(b,n-1);
//      computeSolutionOnBoundary(theMesh->beginBin(b,n-1),theMesh->endBin(b,n-1),t+dt);
//      computeSolutionOnBoundary(theMesh->beginInterfaceBin(b,n-1),theMesh->endInterfaceBin(b,n-1),t+dt);
//      if (b>0) computeSolutionOnBoundary(theMesh->beginInterfaceBin(b-1,n-1),theMesh->endInterfaceBin(b-1,n-1),t+dt);
      assembleVolume(theMesh->beginBin(b,n),theMesh->endBin(b,n),t+dt);
      assembleVolume(theMesh->beginLargeOnInterface(b,n),theMesh->endLargeOnInterface(b,n),t+dt);
      if (b>0) assembleVolume(theMesh->beginSmallOnInterface(b-1,n),theMesh->endSmallOnInterface(b-1,n),t+dt);
      assembleBoundary(theMesh->beginBin(b,n-1),theMesh->endBin(b,n-1),t+dt); 
      assembleBoundary(theMesh->beginInterfaceBin(b,n-1),theMesh->endInterfaceBin(b,n-1),t+dt);
      if (b>0) assembleBoundary(theMesh->beginInterfaceBin(b-1,n-1),theMesh->endInterfaceBin(b-1,n-1),t+dt);
      soln = binptr[b];
      for(it = bin_begin;it != bin_end;++it)	//start with the largest inner cells
	{
	  mEntity *m = (*it);
	  DGCell *cell = (DGCell*)m->getCell();
	  fSize = cell->fSize;
	  const double *RHS = cell->theRightHandSide;
	  
	  for(k = 0;k<cSize;k++)
	    for(i = 0;i<fSize;i++)
	      {
		cell->theFieldsCoefficients->get(k,i) = ( cell->theFieldsCoefficients->get(k,i)
							  +(*RHS++)*dt/cell->getDetJac()+(*soln))*0.5;
		tmp = cell->theFieldsCoefficients->get(k,i)-(*soln);
		++soln;
		resid+=tmp*tmp;
	      }
	  cell->computeMean();
	  cell->ZeroRHS();
	}
      //	if (1) return 0.;
      //2nd stage on large interface elements
      bin_begin = theMesh->beginLargeOnInterface(b,n);
      bin_end = theMesh->endLargeOnInterface(b,n);
      soln = binptrLI[b];
      for(it = bin_begin;it != bin_end;++it)
	{
	  mEntity *m = (*it);
	  DGCell *cell = (DGCell*)m->getCell();
	  fSize = cell->fSize;
	  const double *RHS = cell->theRightHandSide;
	  for(k = 0;k<cSize;k++)
	    for(i = 0;i<fSize;i++)
	      {
		cell->theFieldsCoefficients->get(k,i) = ( cell->theFieldsCoefficients->get(k,i)
							  +(*RHS++)*dt/cell->getDetJac()+(*soln))*0.5;
		//printf("%e\n",cell->theFieldsCoefficients ->get(k,i));
		tmp = cell->theFieldsCoefficients->get(k,i)-(*soln);
		++soln;
		resid+=tmp*tmp;
	      }
	  cell->theFieldsCoefficients->swap();
	  //cell->computeMean();
	  //cell->ZeroRHS();      ???/?????????
	}
      // if (1) return 0.;
      
      if (b>0)
	{
	  soln = binptrSI[b-1]; 
	  bin_begin = theMesh->beginSmallOnInterface(b-1,n); // set coef of small cells on interface to half step
	  bin_end = theMesh->endSmallOnInterface(b-1,n);
	  for(it = bin_begin;it != bin_end;++it)
	    {
	      mEntity *m    = (*it);
	      DGCell  *cell = (DGCell*)m->getCell();
	      fSize = cell->fSize;
	      const double *RHS = cell->theRightHandSide;
	      
	      for(k = 0;k<cSize;k++)
		for(i = 0;i<fSize;i++)
		  {
		    cell->theFieldsCoefficients->get(k,i) = ( cell->theFieldsCoefficients->get(k,i)
							      +(*RHS++)*dt/cell->getDetJac()+(*soln))*0.5;
		    //printf("%e\n",cell->theFieldsCoefficients ->get(k,i));
		    tmp = cell->theFieldsCoefficients->get(k,i)-(*soln);
		    ++soln;
		    resid+=tmp*tmp;
		  }
	      cell->computeMean();
	      cell->ZeroRHS();    
	    }
	}
      dt *= 0.5;
    }
  //if (1) return 0.;
  return 0.0; 
}

void RungeKuttaTVD2adaptive::LIatHalfStep(int b, double DT)
	{
//COmpute and set coefficients in LI to y_n+1/2;
//compute and store in soln solution values on LI for the second stage on small elements
//Already computed y_n+1 on LI is stored in coeff_copy

 int i,k,fSize;
 double dt;
 mMesh::iter it,bin_begin,bin_end;
 dt=DT/pow(2.,b+1);
  
 bin_begin = theMesh->beginLargeOnInterface(b,n); 
 bin_end = theMesh->endLargeOnInterface(b,n); 
 soln=binptrLI[b];	
 for(it = bin_begin;it != bin_end;++it)
   {
     mEntity *m    = (*it);
     DGCell  *cell = (DGCell*)m->getCell();
       fSize         = cell->fSize;
       for(k = 0;k<cSize;k++)
	 for(i = 0;i<fSize;i++)
	   {
	     cell->theFieldsCoefficients ->get(k,i)=(*soln) + (cell->theFieldsCoefficients->get(k,i)-(*soln))+
	       0.25*(cell->theFieldsCoefficients->getCopy(k,i)-2.*cell->theFieldsCoefficients->get(k,i)+(*soln));
	     *soln = cell->theFieldsCoefficients ->get(k,i)+0.5*(cell->theFieldsCoefficients->getCopy(k,i)-(*soln));
	     ++soln;
	   }
   }
}

void RungeKuttaTVD2adaptive::LIatFirstStage(int b)
{
  int i,k,fSize;
  mMesh::iter it,bin_begin,bin_end;
  bin_begin = theMesh->beginLargeOnInterface(b,n); 
  bin_end = theMesh->endLargeOnInterface(b,n); 
  soln=binptrLI[b];	
  for(it = bin_begin;it != bin_end;++it)
    {
      mEntity *m    = (*it);
      DGCell  *cell = (DGCell*)m->getCell();
	  fSize         = cell->fSize;
      for(k = 0;k<cSize;k++)
 	for(i = 0;i<fSize;i++)
	  {
	    cell->theFieldsCoefficients ->get(k,i)=(*soln);
	    ++soln;
	  }
    }
}

void RungeKuttaTVD2adaptive::finalResetLI(int b)
{
  mMesh::iter it,bin_begin,bin_end;
  bin_begin = theMesh->beginLargeOnInterface(b,n); 
  bin_end = theMesh->endLargeOnInterface(b,n); 
  for(it = bin_begin;it != bin_end;++it)
    {
      mEntity *m    = (*it);
      DGCell  *cell = (DGCell*)m->getCell();
      cell->theFieldsCoefficients ->copyBack();
      cell->computeMean();
      cell->ZeroRHS();
    } 
}


double RungeKuttaTVD2adaptive::advanceInTime(double t, double DT)
  /*  TVD Second  Order Runge-Kutta integration with local 
      accurate time-stepping.*/
{
	int i,j;
 double resid=0.0;
 soln=soln_begin;
 int NbBins = theMesh->getNbBins();
 for (i=0; i<NbBins; i++) isFlat[i] = 1;// all element are at the same level

  /*mMesh::iter it;
  const mMesh::iter mesh_begin = theMesh->begin(n);
  const mMesh::iter mesh_end=theMesh->end(n);
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m    = (*it);
      DGCell  *cell = (DGCell*)m->getCell();
	  cell->ZeroRHS();
	  cell->computeMean();
  }*/
 firstStage(0,NbBins,t,DT);
 secondStage(0,NbBins,t,DT);
 for (i=1; i<NbBins; i++) isFlat[i] = 0;
 
 /*
 LIatHalfStep(NbBins-2, DT);
 firstStage(NbBins-1,NbBins,t,DT);
 LIatFirstStage(NbBins-2);
 secondStage(NbBins-1,NbBins,t,DT);
 finalResetLI(NbBins-2);
*/
 //return 0.0;
 i = NbBins-1;
 while (i>0)
 {
 for (i=NbBins-1; i>=1; i--) 
	 if (!isFlat[i])
	 {
		 LIatHalfStep(i-1, DT);
		 firstStage(i,NbBins,t,DT);
		 LIatFirstStage(i-1);
		 secondStage(i,NbBins,t,DT);
		finalResetLI(i-1);
		isFlat[i]=1; 
		if (i!=NbBins-1) for (j=i+1; j<NbBins; j++) isFlat[j] = 0;
		break;
	 }
 }
 return 0.0;
}


/***************************************************************/

double RungeKuttaTVD3::advanceInTime(double t, double dt)
/*  TVD Second  Order Runge-Kutta integration.*/
{
  double  dt_over_detJac,dof;
  int i,j,k,fSize,fSizek;
  double dx;
  soln=soln_begin; 
  mMesh::iter it;
  const mMesh::iter mesh_begin = theMesh->begin(n);
  const mMesh::iter mesh_end=theMesh->end(n);
  
//  computeL2ProjInCutCells();
  assembleVolume(t);
  assembleBoundary(t);  
  
  /**************  SOLVE for K1  **************/
  
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m    = (*it);
      DGCell  *cell = (DGCell*)m->getCell();
      fSize         = cell->fSize;
      const double *RHS = cell->theRightHandSide;
      dof = fSize * cSize;
      double* coeff = cell->theFieldsCoefficients->get();
      /* Copy initial values*/
      
      for(i = 0;i<dof;i++) *(soln++) = coeff[i];
      
      dt = cell->getTimeStep();
      if(cell->details->theFunctionSpace->isOrthogonal())
	{
	  dt_over_detJac = dt/cell->getDetJac();
	  for(i = 0;i<dof;i++) coeff[i] += RHS[i]*dt_over_detJac;
	}
      else
	for(k = 0;k<cSize;k++)
	  {
	    int fSizek = fSize*k;
	    for(i = 0;i<fSize;i++)
	      {
		dx = 0.0;
		for(j = 0;j<fSize;j++){
		  dx += cell->details->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSizek];
		}
		cell->theFieldsCoefficients ->get(k,i)+= dx*dt;
		
	      }
	  }
      
      cell->computeMean();
      cell->ZeroRHS();
    }
  
  if (theLimiter) limit(t);

 /**************  SOLVE for K2  **************/
//computeL2ProjInCutCells();
  assembleVolume(t+dt);
  assembleBoundary(t+dt);  
  
  /*  Reset c values */
  soln=soln_begin;

  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m    = (*it);
      DGCell  *cell = (DGCell*)m->getCell();
      fSize         = cell->fSize;
      dof = fSize * cSize;
      const double *RHS = cell->theRightHandSide;
      double* coeff = cell->theFieldsCoefficients->get();

      dt = cell->getTimeStep();
      if(cell->details->theFunctionSpace->isOrthogonal())
	{
	  dt_over_detJac = dt/cell->getDetJac();
	   for(i = 0;i<dof;i++)
	     coeff[i]=0.25*(coeff[i]+RHS[i]*dt_over_detJac+3.*(*soln++));
	}
      else
	for(k = 0;k<cSize;k++)
	  {
	    int fSizek = fSize*k;
	    for(i = 0;i<fSize;i++)
	      {
		dx = 0.0;
		for(j = 0;j<fSize;j++){
		  dx += cell->details->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSizek];
		}
		cell->theFieldsCoefficients ->get(k,i)= 0.25*(cell->theFieldsCoefficients ->get(k,i)+ dx*dt+3.*(*soln++));
	      }
	  }
      
      cell->computeMean();
      cell->ZeroRHS();
    }
  
  if (theLimiter) limit(t);

  /**************  SOLVE for K3   ************/
  //computeL2ProjInCutCells();
  assembleVolume(t+dt*0.5);
  assembleBoundary(t+dt*0.5);  
  
  /*  Reset c values */
  double resid=0.0,tmp;
  soln=soln_begin;
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      fSize = cell->fSize;
      dof = fSize * cSize;
      double* coeff = cell->theFieldsCoefficients->get();
      const double *RHS = cell->theRightHandSide;
      dt = cell->getTimeStep(); 

      if(cell->details->theFunctionSpace->isOrthogonal())
	{
	  dt_over_detJac = dt /cell->getDetJac();
	  for(i = 0;i<dof;i++)
	    {
	      coeff[i]=(2.0*(coeff[i]+RHS[i]*dt_over_detJac)+(*soln))/3.0;
	      tmp = coeff[i]-(*(soln++));
	      resid+=tmp*tmp;
	    }
	}
      else
	for(k = 0;k<cSize;k++)
	  {
	    fSizek = fSize*k;
	    for(i = 0;i<fSize;i++)
	      {
		dx = 0.0;
		for(j = 0;j<fSize;j++){
		  dx += cell->details->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSizek];
		}
		cell->theFieldsCoefficients->get(k,i) = (2.0*( cell->theFieldsCoefficients->get(k,i)+
							  dx*dt)+(*soln))/3.0;
		tmp = cell->theFieldsCoefficients->get(k,i)-(*(soln++));
		resid+=tmp*tmp;
	      }
	  }
      
      cell->computeMean();
      cell->ZeroRHS();
    }   
  
  if (theLimiter) limit(t+dt);
  //computeL2ProjInCutCells();
  return resid;
}

/***************************************************************/
double ForwardEuler::advanceInTime(double t, double dt)
{
  // double resid = 0.0;
  double dx,dt_over_detJac,dof; 
  int i,j,k;
  double resid = 0;
  mMesh::iter mesh_begin = theMesh->begin(n);
  mMesh::iter mesh_end=theMesh->end(n);

  assembleVolume(t);
  assembleBoundary(t);
  
  for(mMesh::iter it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      dof = cell->fSize * cSize;
      const double *RHS = cell->theRightHandSide;
      double* coeff = cell->theFieldsCoefficients->get();
      dt = cell->getTimeStep();
      
      if(cell->details->theFunctionSpace->isOrthogonal())
	{
	  dt_over_detJac = dt /cell->getDetJac();
	  for(i = 0;i<dof;i++)
	    {
	      dx = coeff[i];
	      coeff[i] += RHS[i]*dt_over_detJac;
	      dx -= coeff[i];
	      resid += dx*dx;
	    }
	}
      else
	{
	  int fSize = cell->fSize;
	  for(k = 0;k<cSize;k++)
	    {
	      for(i = 0;i<fSize;i++)
		{
		  dx = 0.0;
		  for(j = 0;j<fSize;j++)
		    dx += cell->details->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSize*k];
		  
		  cell->theFieldsCoefficients->get(k,i) += dx*dt;
		  resid += dx*dx;
		}
	    }
	}
      cell->computeMean();
      cell->ZeroRHS();
    }
  //computeL2ProjInCutCells();
  if (theLimiter) limit(t);
  
  return resid;
}

double Multigrid::advanceInTime(double t, double dt)
{
  int NbIter = 1;
  int Id,otherId,Id1,Id2, Side;
  double temp[MaxNbEqn],deltaF[MaxNbEqn];
  DGCell* ncell;
  // Perform a timestep at the highest order approximation
  double resid=0;
  double u,v,w,val[MaxNbEqn],fcts[45];
  resid = rk2->advanceInTime(t,dt);
  //resid=0;
  dt*=5000.; 
  
  // Assemble RHS
  //assembleVolume(t+dt);
  //assembleBoundary(t+dt);
  
  // Perform a timestep at lowest approximaton level
  mMesh::iter it;
  mMesh::iter mesh_begin = theMesh->begin(n);
  mMesh::iter mesh_end=theMesh->end(n);
  list<mEntity*>::const_iterator itt,allSubs_end;
  
  // Loop through each cell in the mesh
  for(it = mesh_begin;it != mesh_end;++it)
    {
      // Get the cell and id
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      Id = m->getId(); 
      
      // Transfer the solution to the lowest approximation level
      
      const int fSize = cell->fSize;
      const double *RHS = cell->theRightHandSide;
      //double* coeff = cell->theFieldsCoefficients->get();
      //cell->getMapping()->COG(u,v,w);
      //cell->details->theFunctionSpace->fcts(u,v,w,fcts);
      //cell->interpolate(fcts,val);
	  double  func = cell->pt(0).fcts[0];

      for (int i=0;i<cSize;i++)
	{
	  //cell->details->Up0[i]=val[i];
	  cell->details->Up0[i] = cell->theFieldsCoefficients->get(i,0)*func;
	  cell->details->deltaUp0[i]=0.0;
	  //cell->details->Rp0[i] = 0.;
	  //for (int j=0; j<fSize; j++) cell->details->Rp0[i]+=RHS[i*fSize+j]*fcts[j];
	  // cell->details->Rp0[i]= cell->details->Rp0[i]/(fcts[0]*cell->getDetJac())*cell->getVolume();
	  cell->details->Rp0[i] = RHS[i*fSize]/(func*cell->getDetJac())*cell->getVolume();
	 // printf("rhs %e \n", cell->details->Rp0[i]);
	}

      allSubs_end=cell->details->allSubs.end();
      
      // We now need to calculate the block diagonal matrix D
      
      for (int i=0; i<cSize; i++)
		for (int j=0; j<cSize; j++)
	       DMatrix[i][j]=0.0;
      for(int i=0; i<cSize; i++)
	    DMatrix[i][i]=cell->getVolume()/dt;

      // Loop through boundaries of the cell
      for(itt = cell->details->allSubs.begin();itt!=allSubs_end;++itt) 
	{    
          // Reset the variables
          Id1=0; Id2=0; Side=0;
	  
          // Get the boundary of the cell    
	  DGBoundaryCell *bcell = (DGBoundaryCell*)(*itt)->getCell();      
	  
          // Get the id of the cell on the left
          Id1 = bcell->mleft->getId();
          
          // Get the id of the cell on the right -- if such a cell exists
          if (bcell->mright) Id2 = bcell->mright->getId();
                
          // Figure out which id is the cell and which is the neighbour
	  if (Id==Id1) 
          {
	    // the cell is pointed to by "left" 
	    Side=1;
          }
          else 
	    {
	    //the cell is pointed to by "right"
	      Side=0;
	    }
	  
	  // If all three IDs are the same <==> two neighboring cells have the same ID and we have a problem
          if (Id1==Id==Id2) {printf("We have 2 neigh with the same Id!");exit(0);}
	  
	  // Now add the neighbours's contribution to the D matrix
	  // Add the maximum eigenvalue down the main diagonal
	  double maxEVal = bcell->computeMaxEVal();
	  double bcellSize = bcell->getSize();
	  for(int i=0; i<cSize; i++)
	    DMatrix[i][i]+=0.5*maxEVal*bcellSize;        
	  
	  // Now add the jacobian matrix
	  bcell->computeJacobian(Side, cell->details->Up0,Jac);
	  for(int i=0; i<cSize; i++)
	    for(int j=0; j<cSize; j++)
	      DMatrix[i][j]+=0.5*Jac[i][j]*bcellSize;
	}
      
      // Invert the diagonal matrix
      invmat(DMatrix,cell->details->DInv,cSize);
    }
  for (int q =0; q<NbIter; q++)
    {
		mMesh::iter mesh_begin = theMesh->begin(n);
  mMesh::iter mesh_end=theMesh->end(n);
      for(it = mesh_begin;it != mesh_end;++it)
	{ 
	  mEntity *m = (*it);
	  DGCell *cell = (DGCell*)m->getCell();
	  int Id = m->getId(); 
	  // PERFORM FORWARD SWEEP
	  // We loop through the boundaries of the cell again
	  allSubs_end=cell->details->allSubs.end();
	  for(itt = cell->details->allSubs.begin();itt!=allSubs_end;++itt) 
	    {   
	      // Reset the variables
	      Id2=0; Side=0; 
	      
	      // Get the boundary of the cell    
	      DGBoundaryCell *bcell = (DGBoundaryCell*)(*itt)->getCell();      
	      
	      // Get the id of the cell on the left
	      Id1 = bcell->mleft->getId();
	      
	      // Get the id of the cell on the right
	      if (bcell->mright) Id2 = bcell->mright->getId(); 
	      
	      // Figure out which id is the cell and which is the neighbour
	      if (Id==Id1) 
		{
		  otherId=Id2; //cell is left, neighbor is right 
		  if (Id2) ncell = (DGCell*)bcell->mright->getCell();
		  Side=1;
		}
	      else 
		{
		  otherId = Id1; //cell is right, neighbor is left
		  ncell = (DGCell*)bcell->mleft->getCell();
		  Side=0;
		}
	      
	      // compute contributions of  all of those neighbours with id's less than current ID
	      // Only act cell if there is a neighbor
	      if (Id2 && (otherId<Id))
		{
		  // Update details->deltaUp0
		  double maxEVal; 
		  
		  // Find the deltaF
		  bcell->computeDeltaF(Side,ncell->details->Up0,ncell->details->deltaUp0,deltaF);
		  
		  // Find the maxEVal             
		  maxEVal=bcell->computeMaxEVal();
		  
		  for (int i=0; i<cSize; i++)
		    {
		      cell->details->deltaUp0[i]+=0.5*(deltaF[i]-maxEVal*ncell->details->deltaUp0[i])*bcell->getSize();
		    }
		}
	    }
	  // Now finalize the value for details->deltaUp0*
	  double a[5]={0.0,0.0,0.0,0.0,0.0};
	  for (int i=0; i<cSize; i++)
	    {
	      for (int j=0; j<cSize; j++)
		{
		  a[i]+=cell->details->DInv[i][j]*(cell->details->Rp0[j]-cell->details->deltaUp0[j]);
		}
	    } 
	  for (int i=0; i<cSize; i++)
	    {
	      cell->details->deltaUp0[i]=a[i];
	      //printf("%e \n", a[i]);
	      //printf("%e \n", cell->details->Rp0[i]);
	    }
	}
      
      // Loop through each cell in the mesh - backwards
      mesh_end=--(theMesh->end(n));
	  mesh_begin=--(theMesh->begin(n));
      
      for(it = mesh_end;it != mesh_begin;--it)
	{
	  // Get the cell and id
	  mEntity *m = (*it);
	  DGCell *cell = (DGCell*)m->getCell();
	  Id = m->getId();  
	  
	  allSubs_end=cell->details->allSubs.end();
	  
	  // PERFORM BACKWARD SWEEP
	  // We loop through the boundaries of the cell
	  // Loop through boundaries of the cell
	  for(int i=0; i<cSize; i++)
	    {
	      temp[i]=0.0;
	    }
	  
	  for(itt = cell->details->allSubs.begin();itt!=allSubs_end;++itt) 
	    {    
	      // Reset the variables
	      Id2=0; Side=0; 
	      
	      // Get the boundary of the cell    
	      DGBoundaryCell *bcell = (DGBoundaryCell*)(*itt)->getCell();      
	      
	      Id1 = bcell->mleft->getId();
	      
	      // Get the id of the cell on the right
	      if (bcell->mright) Id2 = bcell->mright->getId(); 
	      
	      // Figure out which id is the cell and which is the neighbour
	      if (Id==Id1) 
		{
		  otherId=Id2; 
		  if(Id2) ncell = (DGCell*)bcell->mright->getCell();
		  Side=1;
		}
	      else 
		{
		  otherId = Id1;
		  ncell = (DGCell*)bcell->mleft->getCell();       
		  Side=0;
		}
	      // Now we iterate over all of those neighbours with id's greater than current ID
	      // Only act if have a neighbor
	      if(Id2 &&(otherId>Id) )
		{          
		  // Update details->deltaUp0
		  double maxEVal; 
		  
		  // Find the deltaF
		  bcell->computeDeltaF(Side,ncell->details->Up0,ncell->details->deltaUp0,deltaF);
		  
		  // Find the maxEVal             
		  maxEVal=bcell->computeMaxEVal();
		  
		  for (int i=0; i<cSize; i++)
		    {
		      temp[i]+=0.5*(deltaF[i]-maxEVal*ncell->details->deltaUp0[i])*bcell->getSize();
		    }
		}
	    }
	  
	  // Now finalize the value for details->deltaUp0
	  // and update the coefficents on the highest level
	  double  func = cell->pt(0).fcts[0];
	  for (int i=0; i<cSize; i++)
	    {
	      double a=0.;
	      for (int j=0; j<cSize; j++)
		{
		  a+=cell->details->DInv[i][j]*temp[j];
		}
	      cell->details->deltaUp0[i]-=a;
	      //  resid+=cell->details->deltaUp0[i]*cell->details->deltaUp0[i];
	    }  
	  
	  if(q==NbIter-1) 
	    {
		  cell->L2Proj(cell->details->deltaUp0);
	      cell->computeMean();
	      cell->ZeroRHS();
	    } 
	  
	}
    }
  //if (theLimiter) limit(t);
  return resid; 
}

/***************************************************************/

void TimeIntegrator::assembleVolume(double t)
{
	//timer start
	clock_gettime(CLOCK_MONOTONIC, &timer1);

	mMesh::iter it;
	int rc;
	void *retval;

	//fork all but last chunk
	for (int i = 0; i < numThreads-1; ++i){
		threadInfo[n][i]._t = t;
		rc = pthread_create(&threadInfo[n][i]._id, NULL,
							parallelVolume, &threadInfo[n][i]);
		if (rc){
			printf("ERROR: return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}

	//work on last chunk
	infoWrapper* info = &threadInfo[n][numThreads-1];
	for(it = info->_begin; it != info->_end; ++it)
    {
		mEntity *m = (*it);
		DGCell *cell = (DGCell*)m->getCell();
		cell->computeVolumeContribution(t);		
    }

	//wait for all threads to finish
	for (int i = 0; i < numThreads-1; ++i){
		pthread_join(threadInfo[n][i]._id, &retval);
	}

	//update timer infos
	clock_gettime(CLOCK_MONOTONIC, &timer2);
	timeComputeVolumeCalls += diff(timer1,timer2).tv_sec;
	timeComputeVolumeCalls += diff(timer1,timer2).tv_nsec * 0.000000001;
	numComputeVolumeCalls++;
}

void TimeIntegrator::assembleBoundary(double t)
{ 
	//timer start
	clock_gettime(CLOCK_MONOTONIC, &timer1);

    mMesh::iter it;
    int rc;
    void *retval;
    int notPhysical;

    //fork all but last chunk
    for (int i = 0; i < numThreads-1; ++i){
        threadInfo[n-1][i]._t = t;
        rc = pthread_create(&threadInfo[n-1][i]._id, NULL,
                            parallelBoundary, &threadInfo[n-1][i]);
        if (rc){
            printf("ERROR: return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    //work on last chunk
    infoWrapper* info = &threadInfo[n-1][numThreads-1];
    for(it = info->_begin; it != info->_end; ++it)
    {
        mEntity *m = (*it);
        DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
        notPhysical = cell->computeBoundaryContributions(info->_t);

        switch(notPhysical)
        {
        case 0:
            //nothing to do
            break;
        case 1:
            cell->getLeftCell()->setToZero(info->_t);
            break;
        case 2:
            cell->getRightCell()->setToZero(info->_t);
            break;
        default:
            cell->getLeftCell()->setToZero(info->_t);
            cell->getRightCell()->setToZero(info->_t);
            break;
        }
    }

    //wait for all threads to finish
    for (int i = 0; i < numThreads-1; ++i){
        pthread_join(threadInfo[n-1][i]._id, &retval);
    }

	//update timer infos
	clock_gettime(CLOCK_MONOTONIC, &timer2);
	timeComputeBoundaryCalls += diff(timer1,timer2).tv_sec;
	timeComputeBoundaryCalls += diff(timer1,timer2).tv_nsec * 0.000000001;
	numComputeBoundaryCalls++;
}

void TimeIntegrator::limit(double time)
{   
  mMesh::iter it;
  mMesh::iter mesh_begin = theMesh->begin(n);
  mMesh::iter mesh_end=theMesh->end(n);
   
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      if(cell->getFunctionSpace()->order()>0) theLimiter->limit(cell,time);
    }
 
 for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
       if(cell->getFunctionSpace()->order()>0) cell->theFieldsCoefficients->copyBack();
      //if(cell->getFunctionSpace()->order()>0) cell->theFieldsCoefficients->copyLimitedCoeffs(cell->limit);
    }
}

/***************************************************************/

void TimeIntegrator::assembleVolume(const mMesh::iter begin,const mMesh::iter end,double t)
{
	/*
	mMesh::iter it;

	vector<pthread_t> threads;
	int numThreads = 0;
	int batches = 0;
	int maxThreads = MAX_THREADS;
	int rc;
	void *retval;

	for(it = begin;it != end; ++it){
		mEntity *m = (*it);
		DGCell *cell = (DGCell*)m->getCell();

		threads.push_back(0);		//add one thread
		numThreads++;
		infoComputeVolume info(cell, t, &threads[threads.size()-1]);
		rc = pthread_create(&threads[threads.size()-1], NULL,
							parallelVolume, &info);
		if (rc){
			printf("ERROR: return code from pthread_create() is %d\n", rc);
			exit(-1);
		}

		if(numThreads >= maxThreads){
			numThreads = 0;
			for (int i = batches*maxThreads; i < threads.size(); ++i){
				pthread_join(threads[i], &retval);
			}	
			batches++;
		}
		//cell->computeVolumeContribution(t);		
	}

	for (int i = batches*maxThreads; i < threads.size(); ++i){
		pthread_join(threads[i], &retval);
	}*/
}

void TimeIntegrator::assembleBoundary(const mMesh::iter begin,const mMesh::iter end,double t)
{  
/* for(mMesh::iter it = begin;it != end;++it)
    {
      mEntity *m = (*it);
      DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
      cell->computeBoundaryContributions(t);
    } */
}

void TimeIntegrator::limit(const mMesh::iter begin, const mMesh::iter end, double time)
{   
  mMesh::iter it;
   
  for(it = begin;it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      if(cell->getFunctionSpace()->order()>0) theLimiter->limit(cell,time);
    }
 
 /*for(it = begin;it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
       if(cell->getFunctionSpace()->order()>0) cell->theFieldsCoefficients->copyBack();
      //if(cell->getFunctionSpace()->order()>0) cell->theFieldsCoefficients->copyLimitedCoeffs(cell->limit);
    }*/
}

/**************************************************************/

void TimeIntegrator::computeFullJump()
{
  mMesh::iter it;
  mMesh::iter mesh_begin = theMesh->begin(n);
	mMesh::iter mesh_end=theMesh->end(n);
  for( it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      cell->ZeroFullJump();
    }  
  for(it = theMesh->begin(n-1);it != theMesh->end(n-1);++it)
    {
      mEntity *m = (*it);
      DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
      cell->computeFullJump();
    }
 }


void TimeIntegrator::computeL2ProjInCutCells()
{
  mMesh::iter it;
  mMesh::iter mesh_begin = theMesh->begin(n);
	mMesh::iter mesh_end=theMesh->end(n);
  for( it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
     if (m->getClassification()->getId()==5000) 
	 {
	  cell->L2ProjCutCell ();
	  cell->computeMean();
	 }
    }  
 }

