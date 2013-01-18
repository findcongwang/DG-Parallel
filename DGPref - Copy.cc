#include "mEntity.h"
#include "mDGMesh.h"
#include "DGCell.h"
#include "DGAnalysis.h"
#include "FunctionSpace.h"
#include "mImportExport.h"
#include "FieldEvaluator.h"
#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
#include "mMesh.h"
#include "mPoint.h"
#include "mVector.h"
#include "Constants.h"
#include <math.h>
#include <stdio.h>

int modifyP(DGCell *theCell,
	    double MaxErr, double MinErr,
	    double globalError,
	    int &pMax, int &pMin, int LEVEL_OF_REFINEMENT)
{
  int initialOrder = theCell->getFunctionSpace()->order();
  int order;
  double error = theCell->getError();
  
  MinErr = 1.e-3;
  
  double rate0 = log10(MinErr) + .1 * (log10(MaxErr) - log10(MinErr)) / 7.;  
  double rate1 = log10(MinErr) + 1. * (log10(MaxErr) - log10(MinErr)) / 7.;  
  double rate2 = log10(MinErr) + 2. * (log10(MaxErr) - log10(MinErr)) / 7.;  
  double rate3 = log10(MinErr) + 3. * (log10(MaxErr) - log10(MinErr)) / 7.;  
  double rate4 = log10(MinErr) + 4. * (log10(MaxErr) - log10(MinErr)) / 7.;  
  double rate5 = log10(MinErr) + 5. * (log10(MaxErr) - log10(MinErr)) / 7.;  
  double rate6 = log10(MinErr) + 6. * (log10(MaxErr) - log10(MinErr)) / 7.;  
  
  int MAXIM = LEVEL_OF_REFINEMENT + 1;

  if(log10(error) > rate6) order = (MAXIM>3)?MAXIM:3;
  else if(log10(error) > rate5) order = (MAXIM>3)?MAXIM:3;
  else if(log10(error) > rate4) order = (MAXIM>3)?MAXIM:3;
  else if(log10(error) > rate3) order = 3;
  else if(log10(error) > rate2) order = 2;
  else if(log10(error) > rate1) order = 1;
  else if(log10(error) < rate0) order = 1;
  else order = 1;
  
  if(order+1 < initialOrder)order = initialOrder-1;
  
  if(order > pMax)pMax=order;
  if(order < pMin)pMin=order;

  return order;
}

void DGAnalysis::smoothP(mEntity *ent)
{
  DGBoundaryCell *theCell = (DGBoundaryCell*)ent->getCell();
  if(!theCell->mright)return;
  DGCell *left  = (DGCell*)theCell->mleft->getCell();
  DGCell *right = (DGCell*)theCell->mright->getCell();
  
  int oleft  = left ->getFunctionSpace()->order();
  int oright = right->getFunctionSpace()->order();

  if(oleft > oright+1)modifyOrder(theCell->mright,oleft - 1);
  else if(oright > oleft+1)modifyOrder(theCell->mleft,oright - 1);

}  


void DGAnalysis::adaptH ()
{
  static int IT = 1;
  printf("--refinement--\n");
  double errMin, errMax;
  computeError(errMax,errMin);
  int pMax,pMin;
  list<mEntity*> splitList;
  list<mEntity*> unSplitList;
  mMesh::iter it;
  mMesh::iter end = theMesh->end(n);
  for(it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      int order = modifyP(cell,errMax,errMin,1.0,pMin,pMax,LEVEL_OF_REFINEMENT);
      int level = theMesh->getRefinementLevel(m); 
      
	if(level < LEVEL_OF_REFINEMENT)
	{
	if (order > (level+1))
		if(1)
	{
	splitList.push_back(m);   //m->print();
	}
	}
      /***begin creating meshes for test problems - local time stepping***/
   /* if (IT==1) 
	{
	  double x=0.,y=0.;
	  for (int i=0;i<m->getNbTemplates(0); i++)
	    {
	      Vertex* v = (Vertex*) m->get(0,i);
	      x += v->point()(0);
		y += v->point()(1);
	    }
	  if (0.25*x>-0.5 &&0.25*x<0.5&&0.25*y>-0.5 &&0.25*y<0.5) splitList.push_back(m);//m->print();
	}
      if (IT==3) 
	{
	  double x=0.,y=0.;
	  for (int i=0;i<m->getNbTemplates(0); i++)
	    {
	      Vertex* v = (Vertex*) m->get(0,i);
	      x += v->point()(0);
	      y += v->point()(1);
	    }
	    if (0.25*x>-0.25 &&0.25*x<0.25&&0.25*y>-0.25 &&0.25*y<0.25) splitList.push_back(m);//m->print();
	}
	*/
      /****end creating meshes for test problems *****/
    }  
  
 
  mMesh::iter end_split = theMesh->endsplit(n);
  for(it = theMesh->beginsplit(n);it != end_split;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
	  int order=1;
      if(theMesh->getRefinementDepth(m) > 0)
	{
		//unSplitList.push_back(m);
	 for(int i=0;i<m->size(n);i++)
	    {
	      DGCell *cells = (DGCell*)m->get(n,i)->getCell();
	      int o = modifyP(cells,errMax,errMin,1.0,pMin,pMax,LEVEL_OF_REFINEMENT);
	      if(o > 1)order = o;
	    }	
	  if(order == 1)
	    {
	      // project the solution form the subcells to the cell
	     unsplitCell(m);
	     unSplitList.push_back(m);
	    }
	}
    }  


  /********SPLITTING*******/
  printf("splitting %d cells\n",splitList.size());  
  list<mEntity*>::const_iterator itt;
  list<mEntity*>::const_iterator split_end = splitList.end();
   for(itt= splitList.begin();itt!=split_end;++itt) 
  {
    theMesh->split(*itt);  //split mesh cell
    splitCell(*itt);       //split DGCell= create new cells and copy values
  }

   /********UNSPLITTING*****/
   printf("unsplitting %d cells\n",unSplitList.size());
   list<mEntity*>::const_iterator unSplit_end = unSplitList.end();
   for(itt= unSplitList.begin();itt!=unSplit_end;++itt) 
     {
    unsplitCell(*itt);       //split DGCell= create new cells and copy values
     }
   theMesh->unsplit(n,unSplitList);
   // theMesh->unsplit(n,(list<mEntity*>) 0);
   

   theMesh->setPeriodicBC(n);
   for(it = theMesh->begin(n-1);it != theMesh->end(n-1);++it)  
     {
      mEntity *ent = (*it);
      //ent->print();
      addBoundaryCell(ent);
     }

   mImportExport io;
   char name[256];
   sprintf(name,"refined-%d.msh",IT++);
   io.exportGmshFile(name,theMesh);
   computeNormalsToCurvedBoundaries();
}


void DGAnalysis::computeError (double &errMax, double &errMin)
{
  mMesh::iter it;
   mMesh::iter end = theMesh->end(n);
  for( it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      cell->ZeroError();
    }  
  mMesh::iter endm1 = theMesh->end(n-1);
  for(it = theMesh->begin(n-1);it != endm1;++it)
    {
      mEntity *m = (*it);
      DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
      cell->computeError();
    }
  int k = 0;
  for(it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      if(!k)
	{
	  errMin = errMax = cell->getError(); 
	  k++;
	}
      else
	{
	  double err = cell->getError();
	  errMin = (errMin < err)?errMin:err;
	  errMax = (errMax > err)?errMax:err;
	}
    }
}

void DGAnalysis::computeFullJump()
{
  mMesh::iter it;
  mMesh::iter end = theMesh->end(n);
  mMesh::iter endm1 = theMesh->end(n-1);
  for( it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      cell->ZeroFullJump();
    }  
  for(it = theMesh->begin(n-1);it != endm1;++it)
    {
      mEntity *m = (*it);
      DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
      cell->computeFullJump();
    }
 }

void DGAnalysis::ExactError(double time, double &error)
{
  error = 0.0;
  mMesh::iter end = theMesh->end(n);
  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      error+=cell->L1exact(time);
    }
}

void DGAnalysis::LinfError(double time, double &error)
{
  error = 0.0;
  double cell_err;
  mMesh::iter end = theMesh->end(n);
  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      cell_err=cell->LinfError(time);
	  if (error<cell_err) error=cell_err;
    }
}


void DGAnalysis::entropyErrorL1(double time, double &error)
{
  error = 0.0;
  mMesh::iter end = theMesh->end(n);
  if (n==2)
  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      error+=cell->entropyErrorL1_2D(time);
    }
  else
  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      error+=cell->entropyErrorL1_3D(time);
    }
}

void DGAnalysis::entropyErrorL2(double time, double &error)
{
  error = 0.0;
  mMesh::iter end = theMesh->end(n);
  if (n==2)
  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      error+=cell->entropyErrorL2_2D(time);
    }
  else
	  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      error+=cell->entropyErrorL2_3D(time);
    }
  error = sqrt(error);
}

void DGAnalysis::PressureError(double time,double &error)
{
  error = 0.0;
  mMesh::iter end = theMesh->end(n);
  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      error+=cell->L1exactPressure(time);
      // printf("in DGPref error %f \n",cell->L1exact() );
    }
}

void DGAnalysis::pressureErrorL2(double time, double &error)
{
  error = 0.0;
  mMesh::iter end = theMesh->end(n);
  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      error+=cell->pressureErrorL2(time);
    }
  error = sqrt(error);
}

void DGAnalysis::densityErrorL2(double time, double &error)
{
  error = 0.0;
  mMesh::iter end = theMesh->end(n);
  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      error+=cell->densityErrorL2(time);
    }
  error = sqrt(error);
}

void DGAnalysis::ApproxError(double &error)
{
  error=0.0;
  int i=0;
  mMesh::iter end = theMesh->end(n);
  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      error+=cell->L1approx();
    }
}


void DGAnalysis::computeTotalPressure(double &error)
{
  double p_t = 0.0;
  double length =0.0;
  double dl;
  mMesh::iter end = theMesh->end(n-1);

  for(mMesh::iter it = theMesh->begin(n-1);it != end;++it)
    if( (*it)->getClassification()->getId() ==20000)
      {	
	DGBoundaryCell *cell = (DGBoundaryCell*) (*it)->getCell();
	p_t+=cell->totalPressure(dl);
	length +=dl;
      }
  error = sqrt(p_t)/length;
}

void DGAnalysis::computeExactTotalPressure(double &error)
{
  double p_t = 0.0;
  double length =0.0;
  double dl;
  mMesh::iter end = theMesh->end(n-1);

  for(mMesh::iter it = theMesh->begin(n-1);it != end;++it)
    if( (*it)->getClassification()->getId() ==50000)
      {	
	DGBoundaryCell *cell = (DGBoundaryCell*) (*it)->getCell();
	p_t+=cell->totalPressure(dl);
	length +=dl;
      }
  error = sqrt(p_t)/length;
}

double getChord();
double getAlpha();

void DGAnalysis::computeLiftDragCoeffs(double &C_l,double &C_d)
{
  double chord = getChord();
  double alpha = getAlpha();
  double lift, drag;
  mMesh::iter end = theMesh->end(n-1);
  double C_x = 0; double C_y = 0;
  for(mMesh::iter it = theMesh->begin(n-1);it != end;++it)
    if( (*it)->getClassification()->getId() ==20000)
      {	
	DGBoundaryCell *cell = (DGBoundaryCell*) (*it)->getCell();
	cell->computeLiftDrag(lift,drag);
	C_y+=lift;
	C_x+=drag;
      }
  double val[MaxNbEqn];
  mPoint p;
  theLaw->getExactSolution()->eval(p,0.0,val);
  double denom = 0.5 * (val[1]*val[1] + val[2]*val[2])*chord/val[0];
  double c=cos(alpha);
  double s = sin(alpha);
  printf("C_X=%e C_Y=%e \n",C_x,C_y);
  C_l = (-C_x*s + C_y*c)/denom;
  C_d = (C_x*c + C_y*s)/denom;
}

void DGAnalysis::plotPressureCoefficient() const
{
  double gm1=0.4;
  FILE *f=fopen("Pres_coeff.dat","w");
  double val[MaxNbEqn];
  mPoint point;
  theLaw->getExactSolution()->eval(point,0.0,val);
  double denom = 0.5 * (val[1]*val[1] + val[2]*val[2])/val[0];
  double p_inf = gm1*(val[3]-0.5 * (val[1]*val[1] + val[2]*val[2])/val[0]);
  mMesh::iter end = theMesh->end(n-1);
  for(mMesh::iter it = theMesh->begin(n-1);it != end;++it)
    if( (*it)->getClassification()->getId() ==20000)
      {	
	DGBoundaryCell *cell = (DGBoundaryCell*) (*it)->getCell();
	int NbPts;
    double x[10],p[10];
	cell->computePressureCoefficient(NbPts,x,p);
	for (int i = 0;i<NbPts;i++)
	fprintf(f,"%e %e\n",x[i],(p_inf-p[i])/denom);
      }
  fclose(f);
}

void DGAnalysis::plotMachOnSurface() const
{
  int NbPts;
  double x[10],m[10];
  double gm1=0.4;
  FILE *f=fopen("Mach.dat","w");
  mMesh::iter end = theMesh->end(n-1);
  for(mMesh::iter it = theMesh->begin(n-1);it != end;++it)
    if( (*it)->getClassification()->getId() ==20000)
      {	
	DGBoundaryCell *cell = (DGBoundaryCell*) (*it)->getCell();
	cell->computeMachNumber(NbPts,x,m);
	for (int i = 0;i<NbPts;i++)
	fprintf(f,"%e %e\n",x[i],m[i]);
      }
  fclose(f);
}

void DGAnalysis::plotPressureLossCoefficient() const
{
  int NbPts;
  double x[10],y[10],p_t[10];
  double gm1=0.4;
  double gamma=1.4;
  FILE *f=fopen("Pres_coeff_loss.dat","w");
  double val[MaxNbEqn];
  mPoint point;
  theLaw->getExactSolution()->eval(point,0.0,val);

  double p_inf = gm1*(val[3]-0.5 * (val[1]*val[1] + val[2]*val[2])/val[0]);
  double p_t_inf=p_inf*pow(1.+0.5*(gamma-1.)*(val[1]*val[1]+val[2]*val[2])/(val[0]*gamma*p_inf),gamma/(gamma-1.));  
 
  mMesh::iter end = theMesh->end(n-1);
  for(mMesh::iter it = theMesh->begin(n-1);it != end;++it)
    if( (*it)->getClassification()->getId() ==20000)
      {	
	DGBoundaryCell *cell = (DGBoundaryCell*) (*it)->getCell();
	cell->computeTotalPressureAtAPoint(NbPts,x,y,p_t);

	for (int i = 0;i<NbPts;i++) 
	  if (y[i]>0.0)
	    fprintf(f,"%e %e\n",x[i],p_t[i]/p_t_inf);
      }
  fclose(f);
}

void DGAnalysis::plotEntropyOnSurface() const
{
  int NbPts;
  double x[10],m[10];
  FILE *f=fopen("Entropy.dat","w");
  mMesh::iter end = theMesh->end(n-1);
  for(mMesh::iter it = theMesh->begin(n-1);it != end;++it)
    if( (*it)->getClassification()->getId() ==20000)
      {	
	DGBoundaryCell *cell = (DGBoundaryCell*) (*it)->getCell();
	cell->computeEntropy(NbPts,x,m);
	for (int i = 0;i<NbPts;i++)
	fprintf(f,"%e %e\n",x[i],m[i]);
      }
  fclose(f);
}

void DGAnalysis::computeTotalMassLoss(double &totalMassLoss) const
{
  double totalVolume = 0.0;
  double totalMass = 0.0;
  mMesh::iter end = theMesh->end(n);
  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {	
      DGCell *cell = (DGCell*) (*it)->getCell();
      totalMass += cell->computeMass();
      totalVolume +=cell->getVolume();
    }
  
  double val[MaxNbEqn];
  mPoint point;
  theLaw->getExactSolution()->eval(point,0.0,val);
  totalMassLoss = (totalMass/(totalVolume*val[0])-1.)*100.;
}


void DGAnalysis::computeBoundaryIntegral(double &error)
{
  error=0.0;
  mMesh::iter end = theMesh->end(n-1);
  for(mMesh::iter it = theMesh->begin(n-1);it != end;++it)
	//  if( (*it)->getClassification()->getId() ==55555)
    {
      mEntity *m = (*it);
      DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
      error+=cell->BoundaryIntegral();
    }
}


void DGAnalysis::ComputeOrientation()
{
  int i=0;
  mMesh::iter it;
  mMesh::iter end = theMesh->end(n-1);
  for(it = theMesh->begin(n-1);it != end;++it)
    {
      mEntity *m = (*it);
      DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
      cell->computeOrientation();
      i++;
    }
  for (it = theMesh->begin(n); it != theMesh->end(n); ++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      for(int i=0;i<3;i++)
	  printf("edges %d ", cell->edgeOrientation[i]);
      printf("\n");
    }
}


void DGAnalysis::adaptP ()
{
  int k = 0;
  double errMin, errMax;
  computeError(errMax,errMin);
  int pMax = 1;
  int pMin = 1;
  mMesh::iter end = theMesh->end(n);
  for(mMesh::iter it = theMesh->begin(n);it != end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      int order = modifyP(cell,errMax,errMin,1.0,pMax,pMin,2);
      modifyOrder(m,order);
    }
  mMesh::iter endm1 = theMesh->end(n-1);
  for(int i=pMin;i<pMax;i++)
    {
      for(mMesh::iter it = theMesh->begin(n-1);it != endm1;++it)
	{
	  mEntity *m = (*it);
	  smoothP(m);
	}
    }
}

int DGAnalysis::splitCell (mEntity* theMeshEntity) // returns the number of additional dof's created by splitting
{
  DGCell *cell = (DGCell*)theMeshEntity->getCell();
  DGCellFieldEvaluator dgfe;
  dgfe.addCell(cell);
  
  int sn = theMeshEntity->size(n);    //uses getnbtemplates
  short o = cell->getFunctionSpace()->order();
  int dof;
  for(int i=0;i<sn;i++)
    {
      mEntity *sub = theMeshEntity->get(n,i);
      addCell(sub,o,dof);
	  
      DGCell *cellsub = (DGCell*)sub->getCell();
      cellsub->L2ProjInitial(&dgfe);
      cellsub->computeMean();
      cellsub->ZeroRHS(); 
      cellsub->computeMaxH();  
      cellsub->ZeroJump();
    }
  DOF+=dof*(sn-1);
  
return dof*(sn-1);
}

void DGAnalysis::unsplitCell (mEntity* theMeshEntity)
{
  DGCell *cell = (DGCell*)theMeshEntity->getCell();
  DGCellFieldEvaluator dgfe;
  int size = theMeshEntity->size(n);
  for(int i=0;i<size;i++)
    {
      mEntity *sub = theMeshEntity->get(n,i);
      DGCell *cellsub = (DGCell*)sub->getCell();
      dgfe.addCell(cellsub);
    }
  cell->L2ProjInitial(&dgfe);
  cell->computeMean();
  cell->ZeroRHS();
}

void DGAnalysis::modifyOrder(mEntity *e, int order)
{
//Old Version - modify!!!

  DGCell *cell2 = (DGCell*)e->getCell();
  DGCellFieldEvaluator dgfe;
  dgfe.addCell(cell2);
  int dummy;
  addCell(e,order,dummy);                          ///CORRECT THIS!!!!!
  DGCell *cell = (DGCell*)e->getCell();
  cell->L2ProjInitial(&dgfe);
  cell->computeMean();
  cell->ZeroRHS();
  delete cell2;
  for(int j=0;j<e->size(n-1);j++)
    {
      mEntity *b = e->get(n-1,j);
      DGBoundaryCell *bc = (DGBoundaryCell*)b->getCell();
      DGCell *c2 = 0;
      DGCell *c1 = (DGCell*)bc->mleft->getCell();
      if(bc->mright)c2 = (DGCell*)bc->mright->getCell();
//      bc->computeInvertMappings();
    }
}

void DGAnalysis::increaseOrder(int newOrder)
{
	mMesh::iter it;
	mMesh::iter end = theMesh->end(n);
  for( it= theMesh->begin(n);it != end ;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();;
	  cell->adaptOrder(newOrder);
	 }
  computeNormalsToCurvedBoundaries();
}

void DGAnalysis::write(ostream &mesh,ostream &data)
{
  // first write the mesh
  mImportExport io;
  io.exportDGFile(mesh,theMesh);

  // then write the contributors
  data.precision(16);
  data <<  TACT << "\n";
  data << ITER << "\n";
  data <<  n << " " << theMesh->size(n) << " " << theMesh->size(n-1) << "\n";
  mMesh::iter it;
	mMesh::iter end = theMesh->end(n);
  for( it= theMesh->begin(n);it != end ;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();;
      theMesh->writeMEntity(data,m);
      cell->write(data);
    }
  mMesh::iter endnm1 = theMesh->end(n-1);
  for(it = theMesh->begin(n-1);it != endnm1;++it)
    {
      mEntity *m = (*it);
      DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
      theMesh->writeMEntity(data,m);
    }  
}

void DGAnalysis::read(istream &is)
{
  int order;
  /* construct cells and recover data */
  int i=0,NbVolCells=0,NbBoundCells=0;
  is >> TACT;
  is >> ITER;
  is >> n >> NbVolCells >> NbBoundCells;
  
  for(i=0;i<NbVolCells;i++)
    {
      mEntity *m = theMesh->readMEntity(is);
      is >> order;
	  int dof;
      DGCell *cell = addCell(m,order,dof);
	  DOF+=dof;
      cell->read(is);
      cell->computeMean();
      cell->ZeroRHS();
	  if (timeSteppingMode ==3) cell->allocateInvJacMatrix();
    }
  /*for(i=0;i<NbBoundCells;i++)
    {
      mEntity *m = theMesh->readMEntity(is);
      DGBoundaryCell *cell = addBoundaryCell(m);
    }*/
 // par->read(is);
}

 
mEntity* DGAnalysis :: eval (const mPoint & p, double *val)
{
  double u,v,w;
  return eval(p,val,u,v,w);
}

mEntity* DGAnalysis :: eval (const mPoint & p, double *val, double &u, double &v, double &w)
{
	mMesh::iter mesh_begin = theMesh->begin(n),it;
	mMesh::iter mesh_end = theMesh->end(n);
	for(it=mesh_begin;it != mesh_end;++it)
    {
      mEntity *ent = (*it);
      DGCell *cell =  (DGCell*)ent->getCell();
		 if(cell->inCell(p,val))
			return ent;
	}
  return 0;
}

void DGAnalysis ::computeNormalsToCurvedBoundaries()
{
	mMesh::iter it;
	mMesh::iter mesh_begin_nm1   = theMesh->begin(n-1);
	mMesh::iter mesh_end_nm1     = theMesh->end(n-1);

	if (n==3)

		    //************USES CENTRES OF CELLS***********************
	{
		for(it = mesh_begin_nm1;it != mesh_end_nm1;++it)
		{
			mEntity *ent = (*it);
			if (ent->getClassification()->getId() == 20000)
			{

				// Create the cell object
				DGBoundaryCell *cell = (DGBoundaryCell*)ent->getCell();

				double x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3;
				int nbPoints, isOpp=0;
				Face *f,*theFace;
				Face *f1 = 0;    // faces of neighboring elements
				Face *f2 = 0;
				Face *f3 = 0;
				mVector n1,n2,n3;  
				Vertex *v0 = (Vertex*) ent->get(0,0);
				Vertex *v1 = (Vertex*) ent->get(0,1);
				Vertex *v2 = (Vertex*) ent->get(0,2);

				int v_Id [3];
				v_Id[0]=v0->getId();
				v_Id[1]=v1->getId();
				v_Id[2]=v2->getId();

				int s = v0->size(2);
				mUpwardAdjacencyContainer::iter itt;
				mUpwardAdjacencyContainer::iter begin = v0->begin(2);
				mUpwardAdjacencyContainer::iter end = v0->end(2);	      
				for (itt= begin;itt!=end ; ++itt)
				{
					f = (Face *) (*itt); //searching for the first neighbor
					if (f->getClassification()->getId() == 20000&& !(f->equal(ent))) 
					{ 
						f1=f;

						int V_Id [3];
						V_Id[0]=f->get(0,0)->getId();
						V_Id[1]=f->get(0,1)->getId();
						V_Id[2]=f->get(0,2)->getId();
						int counter=0;
						for (int i=0; i<3; i++) {
							for (int j=0; j<3; j++) {
								if(V_Id[i]==v_Id[j]) counter++; 
							}
						}

						if (counter==2)
						{
							DGBoundaryCell *f1_cell = (DGBoundaryCell*)f1->getCell();
							if(f1_cell->n(0)*cell->n(0)+f1_cell->n(1)*cell->n(1)+f1_cell->n(2)*cell->n(2)<0) {
								isOpp = 1;
							}
							else {
								mPoint p1 = ((Vertex*) f->get(0,0))->point();
								mPoint p2 = ((Vertex*) f->get(0,1))->point();
								mPoint p3 = ((Vertex*) f->get(0,2))->point();

								x1 = (p1(0)+p2(0)+p3(0))/3.;
								y1 = (p1(1)+p2(1)+p3(1))/3.;
								z1 = (p1(2)+p2(2)+p3(2))/3.;
								double r= sqrt(x1*x1+y1*y1+z1*z1); 
								break;
							}
						}
					}
				}	      

				s = v1->size(2);
				begin = v1->begin(2); end = v1->end(2);	      
				for (itt= begin;itt!=end ; ++itt)
				{
					f = (Face *) (*itt); //searching for the second neighbor
					if (f->getClassification()->getId() == 20000&& !(f->equal(ent)) && !(f->equal(f1))) 
					{ 
						f2 = f;

						int V_Id [3];
						V_Id[0]=f->get(0,0)->getId();
						V_Id[1]=f->get(0,1)->getId();
						V_Id[2]=f->get(0,2)->getId();
						int counter=0;
						for (int i=0; i<3; i++){
							for (int j=0; j<3; j++){
								if(V_Id[i]==v_Id[j]) counter++; 
							}
						}
						if (counter == 2) {
							DGBoundaryCell *f2_cell = (DGBoundaryCell*)f1->getCell();
							if(f2_cell->n(0)*cell->n(0)+f2_cell->n(1)*cell->n(1)+f2_cell->n(2)*cell->n(2)<0) {
								isOpp = 2;
							}
							else {
							mPoint p1 = ((Vertex*) f->get(0,0))->point();
							mPoint p2 = ((Vertex*) f->get(0,1))->point();
							mPoint p3 = ((Vertex*) f->get(0,2))->point();
							x2 = (p1(0)+p2(0)+p3(0))/3.;
							y2 = (p1(1)+p2(1)+p3(1))/3.;
							z2 = (p1(2)+p2(2)+p3(2))/3.;
							double r= sqrt(x2*x2+y2*y2+z2*z2); 
							break;
							}
						}
					}
				}

				s = v2->size(2);
				begin = v2->begin(2); end = v2->end(2);	      
				for (itt= begin;itt!=end ; ++itt)
				{
					f = (Face *) (*itt); //searching for the third neighbor
					if (f->getClassification()->getId() == 20000&& !(f->equal(ent)) && !(f->equal(f1))&& !(f->equal(f2))) 
					{ 
						f3 = f;

						int V_Id [3];
						V_Id[0]=f->get(0,0)->getId();
						V_Id[1]=f->get(0,1)->getId();
						V_Id[2]=f->get(0,2)->getId();
						int counter=0;
						for (int i=0; i<3; i++){
							for (int j=0; j<3; j++){
								if(V_Id[i]==v_Id[j]) counter++; 
							}
						}

						if (counter == 2) {
							DGBoundaryCell *f3_cell = (DGBoundaryCell*)f1->getCell();
							if(f3_cell->n(0)*cell->n(0)+f3_cell->n(1)*cell->n(1)+f3_cell->n(2)*cell->n(2)<0) {
								isOpp = 3;
							}
							else {
								mPoint p1 = ((Vertex*) f->get(0,0))->point();
								mPoint p2 = ((Vertex*) f->get(0,1))->point();
								mPoint p3 = ((Vertex*) f->get(0,2))->point();
								x3 = (p1(0)+p2(0)+p3(0))/3.;
								y3 = (p1(1)+p2(1)+p3(1))/3.;
								z3 = (p1(2)+p2(2)+p3(2))/3.;
								double r= sqrt(x3*x3+y3*y3+z3*z3); 
								break;
							}
						}
					}
				}


				if(isOpp == 1) {
					for(int a=0; a<3; a++) {
						for(int b=0; b<3; b++) {
							if(f2->get(0,a)->getId() == f3->get(0,b)->getId()) { //checking if vertex is shared
								for(itt=begin;itt!=end;++itt) { //loop through all neighbours
									f = (Face *) (*itt);
									if (f->getClassification()->getId() == 20000 && !(f->equal(ent)) && !(f->equal(f2)) && !(f->equal(f3))) {
										if(f->get(0,0)->getId() == f2->get(0,a)->getId() || f->get(0,1)->getId() == f2->get(0,a)->getId() || f->get(0,2)->getId() == f2->get(0,a)->getId()) { //does it have common vertex?
											mPoint p1 = ((Vertex*) f->get(0,0))->point();
											mPoint p2 = ((Vertex*) f->get(0,1))->point();
											mPoint p3 = ((Vertex*) f->get(0,2))->point();

											x1 = (p1(0)+p2(0)+p3(0))/3.;
											y1 = (p1(1)+p2(1)+p3(1))/3.;
											z1 = (p1(2)+p2(2)+p3(2))/3.;
											double r= sqrt(x1*x1+y1*y1+z1*z1); 
											break;
										}
									}
								}
								break;
							}
						}
					}
				}

				if(isOpp == 2) {
					for(int a=0; a<3; a++) {
						for(int b=0; b<3; b++) {
							if(f1->get(0,a)->getId() == f3->get(0,b)->getId()) { //checking if vertex is shared
								for(itt=begin;itt!=end;++itt) { //loop through all neighbours
									f = (Face *) (*itt);
									if (f->getClassification()->getId() == 20000 && !(f->equal(ent)) && !(f->equal(f1)) && !(f->equal(f3))) {
										if(f->get(0,0)->getId() == f1->get(0,a)->getId() || f->get(0,1)->getId() == f1->get(0,a)->getId() || f->get(0,2)->getId() == f1->get(0,a)->getId()) { //does it have common vertex?
											mPoint p1 = ((Vertex*) f->get(0,0))->point();
											mPoint p2 = ((Vertex*) f->get(0,1))->point();
											mPoint p3 = ((Vertex*) f->get(0,2))->point();

											x2 = (p1(0)+p2(0)+p3(0))/3.;
											y2 = (p1(1)+p2(1)+p3(1))/3.;
											z2 = (p1(2)+p2(2)+p3(2))/3.;
											double r= sqrt(x2*x2+y2*y2+z2*z2);
											break;
										}
									}
								}
								break;
							}
						}
					}
				}

				if(isOpp == 3) {
					for(int a=0; a<3; a++) {
						for(int b=0; b<3; b++) {
							if(f1->get(0,a)->getId() == f2->get(0,b)->getId()) { //checking if vertex is shared
								for(itt=begin;itt!=end;++itt) { //loop through all neighbours
									f = (Face *) (*itt);
									if (f->getClassification()->getId() == 20000 && !(f->equal(ent)) && !(f->equal(f1)) && !(f->equal(f2))) {
										if(f->get(0,0)->getId() == f1->get(0,a)->getId() || f->get(0,1)->getId() == f1->get(0,a)->getId() || f->get(0,2)->getId() == f1->get(0,a)->getId()) { //does it have common vertex?
											mPoint p1 = ((Vertex*) f->get(0,0))->point();
											mPoint p2 = ((Vertex*) f->get(0,1))->point();
											mPoint p3 = ((Vertex*) f->get(0,2))->point();

											x3 = (p1(0)+p2(0)+p3(0))/3.;
											y3 = (p1(1)+p2(1)+p3(1))/3.;
											z3 = (p1(2)+p2(2)+p3(2))/3.;
											double r= sqrt(x3*x3+y3*y3+z3*z3); 
											break;
										}
									}
								}
								break;
							}
						}
					}
				}
					//find shared vertex
					//find a cell that shares the vertex
					//use its centroid

				// Get the centre of the current face
				mPoint p1 = ((Vertex*) ent->get(0,0))->point();
				mPoint p2 = ((Vertex*) ent->get(0,1))->point();
				mPoint p3 = ((Vertex*) ent->get(0,2))->point();
				x0 = (p1(0)+p2(0)+p3(0))/3.;
				y0 = (p1(1)+p2(1)+p3(1))/3.;
				z0 = (p1(2)+p2(2)+p3(2))/3.;

				// Create new points storing the centres of the elements
				// The centre of the sphere containing the four points will
				// be returned as point p[4].
				mPoint p[5];
				p[0](0)=x0;
				p[0](1)=y0;
				p[0](2)=z0;
				p[1](0)=x1;
				p[1](1)=y1;
				p[1](2)=z1;
				p[2](0)=x2;
				p[2](1)=y2;
				p[2](2)=z2;
				p[3](0)=x3;
				p[3](1)=y3;
				p[3](2)=z3;

				// Compute the centre and radius of the sphere passing through the four points
				double radius = cell->computeSphere(p);

				// Compute the normal to the circle (sphere in this case)
				cell->normalToCircle(p[4](0),p[4](1),p[4](2));
	  }
	} 
	return;
	}
	}

/*void DGAnalysis ::computeNormalsToCurvedBoundaries()
{
  mMesh::iter it;
  mMesh::iter mesh_begin_nm1   = theMesh->begin(n-1);
  mMesh::iter mesh_end_nm1     = theMesh->end(n-1);

  if (n==3)  

	  //************USES AVERAGES***********************
  {
	  for(it = mesh_begin_nm1;it != mesh_end_nm1;++it)
	  {
		  mEntity *ent = (*it);
		  if (ent->getClassification()->getId() == 20000)
		  {
			  Face *f,*theFace; 
			  Vertex *v0 = (Vertex*) ent->get(0,0);
			  Vertex *v1 = (Vertex*) ent->get(0,1);
			  Vertex *v2 = (Vertex*) ent->get(0,2);

			  mVector norm0;
			  mVector norm1;
			  mVector norm2;

			  double dist0 = 0.0, dist1 = 0.0, dist2 = 0.0;
			  double size0, size1, size2;
			 
			  mUpwardAdjacencyContainer::iter itt0;
			  mUpwardAdjacencyContainer::iter begin0 = v0->begin(2);
			  mUpwardAdjacencyContainer::iter end0 = v0->end(2);	  

			  mUpwardAdjacencyContainer::iter itt1;
			  mUpwardAdjacencyContainer::iter begin1 = v1->begin(2);
			  mUpwardAdjacencyContainer::iter end1 = v1->end(2);

			  mUpwardAdjacencyContainer::iter itt2;
			  mUpwardAdjacencyContainer::iter begin2 = v2->begin(2);
			  mUpwardAdjacencyContainer::iter end2 = v2->end(2);

			  //Get the vertices
			  mPoint p0 = ((Vertex*) ent->get(0,0))->point();
			  mPoint p1 = ((Vertex*) ent->get(0,1))->point();
			  mPoint p2 = ((Vertex*) ent->get(0,2))->point();


			  for (itt0 = begin0;itt0 != end0 ; ++itt0)//find normal on first vertex
			  {
				  f = (Face *) (*itt0);
				  if(f->getClassification()->getId() == 20000) {
					  DGBoundaryCell *cell = (DGBoundaryCell*)f->getCell();
					  DGBoundaryCell *ent_cell = (DGBoundaryCell*)ent->getCell();
					  if(cell->n(0)*ent_cell->n(0)+cell->n(1)*ent_cell->n(1)+cell->n(2)*ent_cell->n(2)>0) {
						  //dist0 = 1./sqrt((p0(0)-cell->n(0))*(p0(0)-cell->n(0))+(p0(1)-cell->n(1))*(p0(1)-cell->n(1))+(p0(1)-cell->n(1))*(p0(1)-cell->n(1)));
						  //size0 = cell->size;
						  norm0 += cell->n;
					  }
				  }
			  }
			  norm0.normalize();

			  for (itt1 = begin1;itt1 != end1 ; ++itt1)//find normal on second vertex
			  {
				  f = (Face *) (*itt1);
				  if(f->getClassification()->getId() == 20000) {
					  DGBoundaryCell *cell = (DGBoundaryCell*)f->getCell();
					  DGBoundaryCell *ent_cell = (DGBoundaryCell*)ent->getCell();
					  if(cell->n(0)*ent_cell->n(0)+cell->n(1)*ent_cell->n(1)+cell->n(2)*ent_cell->n(2)>0) {
						  //dist1 = 1./sqrt((p1(0)-cell->n(0))*(p1(0)-cell->n(0))+(p1(1)-cell->n(1))*(p1(1)-cell->n(1))+(p1(1)-cell->n(1))*(p1(1)-cell->n(1)));
						  //size1 = cell->size;
						  norm1 += cell->n;
					  }
				  }
			  }
			  norm1.normalize();

			  for (itt2 = begin2;itt2 != end2 ; ++itt2)//find normal on third vertex
			  {
				  f = (Face *) (*itt2);
				  if(f->getClassification()->getId() == 20000) {
					  DGBoundaryCell *cell = (DGBoundaryCell*)f->getCell();
					  DGBoundaryCell *ent_cell = (DGBoundaryCell*)ent->getCell();
					  if(cell->n(0)*ent_cell->n(0)+cell->n(1)*ent_cell->n(1)+cell->n(2)*ent_cell->n(2)>0) {
						  //dist2 = 1./sqrt((p2(0)-cell->n(0))*(p2(0)-cell->n(0))+(p2(1)-cell->n(1))*(p2(1)-cell->n(1))+(p2(1)-cell->n(1))*(p2(1)-cell->n(1)));
						  //size2 = cell->size;
						  norm2 += cell->n;
					  }
				  }
			  }
			  norm2.normalize();
			  

			  // Create the cell object
			  DGBoundaryCell *cell = (DGBoundaryCell*)ent->getCell();

			  // Compute the normal to the cell
			  cell->computeCurvedNormals(p0,p1,p2,norm0,norm1,norm2);
	  }
	} 
	return;
  }
}*/


/*  if (n==3)
  
    //************USES CENTRES OF CELLS***********************
    {
      for(it = mesh_begin_nm1;it != mesh_end_nm1;++it)
	{
	  mEntity *ent = (*it);
	  if (ent->getClassification()->getId() == 20000)
	    {
	      double x0,x1,x2,x3,y0,y1,y2,y3,z0,z1,z2,z3;
	      int nbPoints;
	      Face *f;
	      Face *f1 = 0;    // faces of neighboring elements
	      Face *f2 = 0;
	      Face *f3 = 0;
	      mVector n1,n2,n3;  
	      Vertex *v0 = (Vertex*) ent->get(0,0);
		  int s = v0->size(2);
	      mUpwardAdjacencyContainer::iter itt;
	      mUpwardAdjacencyContainer::iter begin = v0->begin(2);
	      mUpwardAdjacencyContainer::iter end = v0->end(2);	      
	      for (itt= begin;itt!=end ; ++itt)
		{
		  f = (Face *) (*itt); //searching for the first neighbor
		  if (f->getClassification()->getId() == 20000&& !(f->equal(ent))) 
		    { 
		      f1 = f;
			  mPoint p1 = ((Vertex*) f->get(0,0))->point();
			  mPoint p2 = ((Vertex*) f->get(0,1))->point();
			  mPoint p3 = ((Vertex*) f->get(0,2))->point();
		      x1 = (p1(0)+p2(0)+p3(0))/3.;
			  y1 = (p1(1)+p2(1)+p3(1))/3.;
			  z1 = (p1(2)+p2(2)+p3(2))/3.;
			  double r= sqrt(x1*x1+y1*y1+z1*z1); 
			  break;
		    }
		}	      
		  Vertex *v1 = (Vertex*) ent->get(0,1);
		  s = v1->size(2);
	      begin = v1->begin(2); end = v1->end(2);	      
	      for (itt= begin;itt!=end ; ++itt)
		{
		  f = (Face *) (*itt); //searching for the second neighbor
		  if (f->getClassification()->getId() == 20000&& !(f->equal(ent)) && !(f->equal(f1))) 
		    { 
		      f2 = f;
			  mPoint p1 = ((Vertex*) f->get(0,0))->point();
			  mPoint p2 = ((Vertex*) f->get(0,1))->point();
			  mPoint p3 = ((Vertex*) f->get(0,2))->point();
		      x2 = (p1(0)+p2(0)+p3(0))/3.;
			  y2 = (p1(1)+p2(1)+p3(1))/3.;
			  z2 = (p1(2)+p2(2)+p3(2))/3.;
			  double r= sqrt(x2*x2+y2*y2+z2*z2); 
			  break;
		    }
		}
		  Vertex *v2 = (Vertex*) ent->get(0,2);
		  s = v2->size(2);
	      begin = v2->begin(2); end = v2->end(2);	      
	      for (itt= begin;itt!=end ; ++itt)
		{
		  f = (Face *) (*itt); //searching for the third neighbor
		  if (f->getClassification()->getId() == 20000&& !(f->equal(ent)) && !(f->equal(f1))&& !(f->equal(f2))) 
		    { 
		      f3 = f;
			  mPoint p1 = ((Vertex*) f->get(0,0))->point();
			  mPoint p2 = ((Vertex*) f->get(0,1))->point();
			  mPoint p3 = ((Vertex*) f->get(0,2))->point();
		      x3 = (p1(0)+p2(0)+p3(0))/3.;
			  y3 = (p1(1)+p2(1)+p3(1))/3.;
			  z3 = (p1(2)+p2(2)+p3(2))/3.;
			  double r= sqrt(x3*x3+y3*y3+z3*z3); 
			  break;
		    }
		}
	
        // Get the centre of the current face
		mPoint p1 = ((Vertex*) ent->get(0,0))->point();
    	mPoint p2 = ((Vertex*) ent->get(0,1))->point();
    	mPoint p3 = ((Vertex*) ent->get(0,2))->point();
        x0 = (p1(0)+p2(0)+p3(0))/3.;
    	y0 = (p1(1)+p2(1)+p3(1))/3.;
    	z0 = (p1(2)+p2(2)+p3(2))/3.;
    	
    	// Create new points storing the centres of the elements
    	// The centre of the sphere containing the four points will
    	// be returned as point p[4].
 	    mPoint p[5];
 	    p[0](0)=x0;
 	    p[0](1)=y0;
 	    p[0](2)=z0;
 	    p[1](0)=x1;
 	    p[1](1)=y1;
 	    p[1](2)=z1;
 	    p[2](0)=x2;
 	    p[2](1)=y2;
 	    p[2](2)=z2;
 	    p[3](0)=x3;
 	    p[3](1)=y3;
 	    p[3](2)=z3;
 	     
 	    // Create the cell object
 	    DGBoundaryCell *cell = (DGBoundaryCell*)ent->getCell();
 	    
 	    // Compute the centre and radius of the sphere passing through the four points
 	    double radius = cell->computeSphere(p);
 	        	
    	// Compute the normal to the circle (sphere in this case)
    	cell->normalToCircle(p[4](0),p[4](1),p[4](2));
	  }
	} 
    return;
    }
*/


/*  for(it = mesh_begin_nm1;it != mesh_end_nm1;++it)
    {
      mEntity *ent = (*it);
      if (ent->getClassification()->getId() == 20000)
	{
	  int i=0;
	  double xleft,yleft,xright,yright;
	  int nbPoints;
	  Edge *e;
	  Edge *eleft = 0;  //edges on the boundary to the left and right from this edge
	  Edge *eright = 0;   // left and right are just names, left is not necessary on the left 
	  mVector nleft,nright;  
	  Vertex *v0 = (Vertex*) ent->get(0,0);
	  mUpwardAdjacencyContainer::iter itt;
	  mUpwardAdjacencyContainer::iter begin = v0->begin(1);
	  mUpwardAdjacencyContainer::iter end = v0->end(1);
	  int s = v0->size(1);
	  for (itt= begin;itt!=end ; ++itt)
	    {
	      e = (Edge *) (*itt); //searching for the first neighbor
	      if (e->getClassification()->getId() == 20000&& !(e->equal(ent))) 
		{ 
		  eleft = e;
		  DGBoundaryCell *leftCell = (DGBoundaryCell*)e->getCell();
		  nbPoints = leftCell->getNbPtGauss();
		  nleft = leftCell->n;
		  xleft = 0.5*(leftCell->boundaryGaussPoints[0]->x+leftCell->boundaryGaussPoints[nbPoints-1]->x);
		  yleft = 0.5*(leftCell->boundaryGaussPoints[0]->y+leftCell->boundaryGaussPoints[nbPoints-1]->y);
		  break;
		}
	    }
	  Vertex *v1 = (Vertex*) ent->get(0,1);
	  begin = v1->begin(1);
	  end = v1->end(1);
	  int s1 = v1->size(1);
	  for (itt= begin;itt!=end ; ++itt)
	    { 
		   e = (Edge *) (*itt);
		   if (e->getClassification()->getId() == 20000&& 
		       !(e->equal(ent))) 
		{ 
		  eright = e;
		  DGBoundaryCell *rightCell = (DGBoundaryCell*)e->getCell();
		  nbPoints = rightCell->getNbPtGauss();
		  nright = rightCell->n;
		  xright = 0.5*(rightCell->boundaryGaussPoints[0]->x+rightCell->boundaryGaussPoints[nbPoints-1]->x);
		  yright = 0.5*(rightCell->boundaryGaussPoints[0]->y+rightCell->boundaryGaussPoints[nbPoints-1]->y);
		  break;
		}
	    }
	  DGBoundaryCell *cell = (DGBoundaryCell*)ent->getCell();
	 
	  if(cell->n * nleft <=0.0 || !eleft ) //corner or no left element
	    {
	      mPoint p0 = v0->point();  
	      mPoint p1 = v1->point();
	      mPoint p2;
	      if (eright->vertex(0) != v1) p2 = eright->vertex(0)->point();
	      else p2 = eright->vertex(1)->point();
	      cell->normalToCircle(p0,p1,p2);
	    }  
	  else if(cell->n * nright <=0.0|| !eright) //corner or no right relfecting element 
	    {
	      mPoint p0 = v1->point();  
	      mPoint p1 = v0->point();
	      mPoint p2;
	      if (eleft->vertex(0) != v1) p2 = eleft->vertex(0)->point();
	      else p2 = eleft->vertex(1)->point();
	      cell->normalToCircle(p0,p1,p2);
	    }  
	  else          //smooth curve
	    {
	      double r1,r2,r,x,y,x_c,y_c,segment,radius;
	      mPoint p0 = v0->point();    //vertices of the edge 
	      mPoint p1 = v1->point();
	      mPoint p2;
	      if (eright->vertex(0) != v1) p2 = eright->vertex(0)->point(); //the 3d point is a vertex of the right neigh. 
	      else p2 = eright->vertex(1)->point();
	      r1=cell->computeRadius(p0,p1,p2);
	      if (eleft->vertex(0) != v1) p2 = eleft->vertex(0)->point(); //now of the left neigh. 
	      else p2 = eleft->vertex(1)->point();
	      r2=cell->computeRadius(p0,p1,p2);
	      r=(r1+r2)*0.5; //average of two radii
		  x = 0.5*(p0(0)+p1(0));
		  y = 0.5*(p0(1)+p1(1));
	      segment=(p0(0)-x)*(p0(0)-x) + (p0(1)-y)*(p0(1)-y);
	      radius = sqrt(r*r-segment);
	      double x_c1=x+cell->n(0)*radius;
		  double y_c1=y+cell->n(1)*radius;
		  double x_c2=x-cell->n(0)*radius;
		  double y_c2=y-cell->n(1)*radius;
	      if ((p2(0)-x_c1)*(p2(0)-x_c1)+(p2(1)-y_c1)*(p2(1)-y_c1) < (p2(0)-x_c2)*(p2(0)-x_c2)+(p2(1)-y_c2)*(p2(1)-y_c2))
		  {x_c=x_c1; y_c=y_c1; }
		  else {x_c=x_c2; y_c=y_c2;}
	      cell->normalToCircle(x_c,y_c,0);
	    }  
	  //cell->normalToNACA();
	}
    }
}*/

void DGAnalysis::setPeriodicBC()
{
mMesh::iter it;
  mMesh::iter mesh_begin_nm1   = theMesh->begin(n-1);
  mMesh::iter mesh_end_nm1     = theMesh->end(n-1);
  for(it = mesh_begin_nm1;it != mesh_end_nm1;++it)
    {
      mEntity *ent = (*it);
      if (ent->getClassification()->getId() == 510 ||ent->getClassification()->getId() == 610)
	  {
      DGBoundaryCell* cell = (DGBoundaryCell*) ent->getCell();
      cell->setPeriodicBC();
	  }
  }
}
