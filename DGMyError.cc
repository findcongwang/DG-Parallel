#include "DGAnalysis.h"
#include "DGMyError.h"
//#include "Mapping.h"
#include "FunctionSpace.h"
#include "ConservationLaw.h"
#include "mDGMesh.h"
#include "mEntity.h"
#include "DGLimiter.h"
#include "DGCell.h"
#include "FieldEvaluator.h"
#include "mImportExport.h"
#include <stdio.h>
#include <math.h>
#include <list>
#include <time.h>
#include "mCompiler.h"

DGMyError::DGMyError(mDGMesh *m,ConservationLaw *c):theMesh(m),theLaw(c) 
{
	TACT = 0.0;
  TEND = 0.2;
  DT   = 0.1;
  HREF = PREF = 0;
  REF_SAMPLING_TIME = 1.e10;
  LEVEL_OF_REFINEMENT = 1;
  CFLMAX = 0.1;

  if(theMesh->size(3))n = 3;
  else if(theMesh->size(2))n = 2;
  else n = 1;
 mMesh::iter it;
// compute initial solution
   //compute last p+1 rows in Mass matrices

  for(it = theMesh->begin(n);it != theMesh->end(n);++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
       cell->ZeroErrorMassMatrices();
    }
  
    for(it = theMesh->begin(n-1);it != theMesh->end(n-1);++it)
    {
      mEntity *m = (*it);
      DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
      cell->computeBoundaryRestrictions(TACT);
    }
  for(it = theMesh->begin(n);it != theMesh->end(n);++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      cell->computeErrorMassMatrices();
    }

	for(it = theMesh->begin(n);it != theMesh->end(n);++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      cell->L2ProjInitialError(c->getFarField());
	  cell->ZeroErrorMassMatrices();
      cell->ZeroErrorRHS();
	  cell->computeErrorMean();
    }	
}

void DGMyError::computeMyError(double dt)
{
	 mMesh::iter it;
	  double residError = 0.0;
	  // Assemble last p+1 rows of the Error Mass Matrix	  
	  for(it = theMesh->begin(n-1);it != theMesh->end(n-1);++it)
	    {
	      mEntity *m = (*it);
	 DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
      	 //Compute Error boundary contributions
	 cell->computeErrorBoundaryContributions(TACT);
	 cell->computeBoundaryRestrictions(TACT);
	// cell->computeJump();
     //cell->computeErrorJump(TACT);
	    }

	  for(it = theMesh->begin(n);it != theMesh->end(n);++it)
	  {
	    mEntity *m = (*it);
	    DGCell *cell = (DGCell*)m->getCell();
	   cell->computeErrorVolumeContribution();
		cell->computeErrorMassMatrices();
	//	if(theLimiter && cell->getFunctionSpace()->order()>0)theLimiter->limit(cell);
		cell->ZeroJump();
	  }
	
	  for(it = theMesh->begin(n);it != theMesh->end(n);++it)
	    {
	      mEntity *m = (*it);
	      DGCell *cell = (DGCell*)m->getCell();
	      residError += cell->advanceErrorInTime(dt);
	      cell->ZeroErrorRHS();
		  cell->ZeroJump();
	      cell->ZeroErrorMassMatrices();
	    }

}
