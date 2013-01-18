#include "DGAnalysis.h"
#include "Mapping.h"
#include "FunctionSpace.h"
#include "ConservationLaw.h"
#include "mDGMesh.h"
#include "mEntity.h"
#include "DGLimiter.h"
#include "DGCell.h"
#include "FieldEvaluator.h"
#include "mImportExport.h"
#include "TimeIntegrator.h"
#include "Vertex.h"
#include "Edge.h"
#include "Integrator.h"
#include "CutCell.h"
#include "Geometry.h"
#include <stdio.h>
#include <math.h>
#include <list>
#include <time.h>
#include <cstring>
#include <cstdlib>

using namespace std;

static void recurGetAllsubs(mEntity *ent, list<mEntity*> &lis)
{
  int n = ent->getLevel();
  if(!ent->isAdjacencyCreated(n))lis.push_back(ent);
  else
    for(int i=0;i<ent->size(n);i++)recurGetAllsubs(ent->get(n,i),lis);
}

void DGAnalysis::parse()
{
  FILE *f = fopen(".dgrc","r");
  if(!f)return;
  char line[256],str[256];
  while(!feof(f))
    {
      fgets(line,255,f);
      if(line[0] != '#')
	{
	  sscanf(line,"%s",str);
	  if(!strcmp(str,"InitialTime"))sscanf(line,"%s %lf",str,&TACT);
	  if(!strcmp(str,"FinalTime"))sscanf(line,"%s %lf",str,&TEND);	                      
	  if(!strcmp(str,"TimeStepping"))
	    sscanf(line,"%s %d",str,&timeSteppingMode);	                      
	  if(!strcmp(str,"Recover")) {
	    recoverFile = new char[256];
	    sscanf(line,"%s %s",str,recoverFile);
	  }
	  if(!strcmp(str,"StepsToDump")) sscanf(line,"%s %d",str,&steps_to_dump);
	  if(!strcmp(str,"hRef"))
	    {
	      sscanf(line,"%s %lf %d",str,&REF_SAMPLING_TIME,&LEVEL_OF_REFINEMENT);	      
	      HREF = 1;
	    }
	  if(!strcmp(str,"pRef"))
	    {
	      sscanf(line,"%s %lf",str,&REF_SAMPLING_TIME);	      
	      PREF = 1;
	    }
	  if(!strcmp(str,"CFLMAX"))sscanf(line,"%s %lf",str,&CFLMAX);	      
	  if(!strcmp(str,"DT"))sscanf(line,"%s %lf",str,&DT);
	  if(!strcmp(str,"compute_Linf_error"))sscanf(line,"%s %d",str,&compute_Linf_error);
	  if(!strcmp(str,"compute_density_errorL1"))sscanf(line,"%s %d",str,&compute_density_errorL1);
	  if(!strcmp(str,"compute_density_errorL2"))sscanf(line,"%s %d",str,&compute_density_errorL2);
	  if(!strcmp(str,"compute_entropy_errorL1"))sscanf(line,"%s %d",str,&compute_entropy_errorL1);
	  if(!strcmp(str,"compute_entropy_errorL2"))sscanf(line,"%s %d",str,&compute_entropy_errorL2);
	  if(!strcmp(str,"compute_pressure_error"))sscanf(line,"%s %d",str,&compute_pressure_error);
	  if(!strcmp(str,"compute_pressure_errorL2"))sscanf(line,"%s %d",str,&compute_pressure_errorL2);
	  if(!strcmp(str,"compute_total_pressure"))sscanf(line,"%s %d",str,&compute_total_pressure);
	  if(!strcmp(str,"compute_total_mass_loss"))sscanf(line,"%s %d",str,&compute_total_mass_loss);
	  if(!strcmp(str,"compute_exact_total_pressure"))sscanf(line,"%s %d",str,&compute_exact_total_pressure);
	  if(!strcmp(str,"compute_lift_drag_coeffs"))sscanf(line,"%s %d",str,&compute_lift_drag_coeffs);
	  if(!strcmp(str,"plot_pressure_coefficient"))sscanf(line,"%s %d",str,&plot_pressure_coefficient);
	  if(!strcmp(str,"plot_mach_on_surface"))sscanf(line,"%s %d",str,&plot_mach_on_surface);
	  if(!strcmp(str,"plot_pressure_loss_coefficient"))sscanf(line,"%s %d",str,&plot_pressure_loss_coefficient);
	  if(!strcmp(str,"plot_entropy_on_surface"))sscanf(line,"%s %d",str,&plot_entropy_on_surface);
	}
    }
  fclose(f);
}

DGCell* DGAnalysis::addCell(mEntity *ent, int order, int& dof)
{
  FunctionSpace *fs;
  switch(ent->getType())
    {
    case mEntity::TRI  : 
		if (ent->getgeomOrder() == 1) fs = new OrthogonalTriangleFunctionSpace(order);
		else						  fs = new TriangleFunctionSpace(order);
		break;
    case mEntity::QUAD : fs = new QuadFunctionSpace(order); break;
    case mEntity::TET  : 
		if (ent->getgeomOrder() == 1) fs = new OrthogonalTetFunctionSpace(order);
		else						  fs = new TetFunctionSpace(order);
		break;
    case mEntity::HEX  : fs = new HexFunctionSpace(order); break;
    }
  
  //  er = new OrthogonalTriangleFunctionSpace(order +1);
  DGCell *cell;
  MeshMapping *theMapping;
  if (ent->getType()==mEntity::QUAD || ((ent->getgeomOrder()) > 1)  ) theMapping = new MeshMapping(ent);
  else theMapping = new ConstantMeshMapping(ent);     

  if (ent->getClassification()->getId()!=5000) 
    cell = new DGCell(theLaw, ent, fs ,0,theMapping,&theGaussIntegrator);
  else cell = new CutCell(theLaw, ent, fs ,0,theMapping,&theGaussIntegrator,theGeometry);
  ent->setCell(cell);
  dof = theLaw->getNbFields()*fs->size(); 
  return cell;
}

DGBoundaryCell* DGAnalysis::addBoundaryCell(mEntity *ent)
{
  DGBoundaryCell *cell = (DGBoundaryCell*)ent->getCell();
  if(cell) {cell->check(); return cell;}
  cell = new DGBoundaryCell(ent, new MeshMapping(ent),&theGaussIntegrator);
  ent->setCell(cell);
  return cell;
}


DGAnalysis::DGAnalysis(mDGMesh *m, ConservationLaw *c, DGLimiter *l, int order)
  : theMesh(m),theLaw(c),theLimiter(l)
{
  DOF=0;
  TACT = 0.0;
  TEND = 0.2;
  DT   = 0.001;
  ITER = 0;
  HREF = PREF = 0;
  REF_SAMPLING_TIME = 1.e10;
  LEVEL_OF_REFINEMENT = 1;
  CFLMAX = 0.5;
  recoverFile = 0;
  steps_to_dump = 1000;
  timeSteppingMode =_regular_;
  compute_density_errorL1 =0;
  compute_density_errorL2 = 0;
  compute_entropy_errorL1 = 0;
  compute_entropy_errorL2 = 0;
  compute_pressure_error  = 0;
  compute_pressure_errorL2= 0;
  compute_total_pressure  = 0;
  compute_total_mass_loss  = 0;
  compute_exact_total_pressure=0;
  compute_lift_drag_coeffs=0;
  plot_pressure_coefficient=0;
  plot_mach_on_surface     =0;
  plot_pressure_loss_coefficient = 0;
  plot_entropy_on_surface  =0;
  if(theMesh->size(3)) n = 3;
  else if(theMesh->size(2)) n = 2;
  else n = 1;
  
  mMesh::iter it;
  mMesh::iter mesh_begin       = theMesh->begin(n);
  mMesh::iter mesh_end         = theMesh->end(n);
  mMesh::iter mesh_begin_nm1   = theMesh->begin(n-1);
  mMesh::iter mesh_end_nm1     = theMesh->end(n-1);
   
  theGeometry = new Geometry;

  parse();
  
  if(recoverFile)
    {
      printf("recovering %s ...\n",recoverFile);      
      ifstream ifs;
      ifs.open(recoverFile);
      if(!ifs.is_open()) printf("Can't find the recovery file \n");
      else read(ifs);
      ifs.close();

	  for(it = mesh_begin_nm1;it != mesh_end_nm1;++it)
	{
	  mEntity *ent = (*it);
	  addBoundaryCell(ent);
	}
    }
  else
    {
      clock_t t1=clock();
      for(it = mesh_begin;it != mesh_end;++it)
	{
	  mEntity *ent = (*it);
	  int dof; // number of DoF on this element
	  DGCell *cell = addCell(ent,order,dof);
	  DOF+=dof;  // add to global DOF;
	  cell->L2ProjInitial(c->getFarField());
	  cell->computeMean();
	  cell->ZeroRHS(); 
	  if (timeSteppingMode ==3) cell->allocateInvJacMatrix();
	}
      clock_t t2=clock();
      printf("Volume cells created in %e sec \n", double (t2-t1)/CLOCKS_PER_SEC);
      t1=clock();
      for(it = mesh_begin_nm1;it != mesh_end_nm1;++it)
	{
	  mEntity *ent = (*it);
	  addBoundaryCell(ent);
	}
      t2=clock();
      printf("Boundary cells created in %e sec \n", double (t2-t1)/CLOCKS_PER_SEC);
    }
  
/*for(it = mesh_begin;it != mesh_end;++it) //cut cell
      {
	mEntity *m = (*it);
	if (m->getClassification()->getId()==5000)
	{
	DGCell *cell = (DGCell*)m->getCell();
	cell->init();
	cell->computeMean();
	}
	}
*/
  computeNormalsToCurvedBoundaries();
  if (timeSteppingMode==2) adaptH();
  if (timeSteppingMode==2) adaptH();
  setPeriodicBC();
  if (timeSteppingMode==5) sortCellsBySize();
  if (timeSteppingMode==5)
    theIntegrator = new RungeKuttaTVD2adaptive (theMesh,theLimiter,theLaw->getNbFields(),DOF);
  else  if (timeSteppingMode==3)
    theIntegrator = new Multigrid (theMesh,theLimiter,theLaw->getNbFields(),DOF);
  else
    {
      if (0 == order) theIntegrator = new ForwardEuler (theMesh,theLimiter,theLaw->getNbFields(),DOF);
      else if (order<2) theIntegrator = new RungeKuttaTVD2 (theMesh,theLimiter,theLaw->getNbFields(),DOF);
      else if (order<3) theIntegrator = new RungeKuttaTVD3 (theMesh,theLimiter,theLaw->getNbFields(),DOF);
      else theIntegrator = new RungeKutta4 (theMesh,theLimiter,theLaw->getNbFields(),DOF); 
    } 
  //if (theLimiter)
    for(it = mesh_begin;it != mesh_end;++it)
      {
	mEntity *m = (*it);
	DGCell *cell = (DGCell*)m->getCell();
	//cell->ZeroJump();
	int k;
	for(k=0;k<cell->theMeshEntity->size(n-1);k++)
	  recurGetAllsubs(cell->theMeshEntity->get(n-1,k),cell->allSubs);
	/*int size0=m->size(0);
	for(k=0;k<size0;k++)
	  {
	    mEntity *ent = cell->theMeshEntity->get(0,k);
	    //printf("size2= %d\n",ent->size(2));
	    mUpwardAdjacencyContainer::iter itt;
	    mUpwardAdjacencyContainer::iter begin = ent->begin(2);
	    mUpwardAdjacencyContainer::iter end = ent->end(2);
	    for (itt= begin;itt!=end ; ++itt)
	      {
		if((*itt)!=m) 
		  cell->allVert.push_back(*itt);
	      }
	    //	    list<mEntity*>::const_iterator lit;
	    //list<mEntity*>::const_iterator lit_end=cell->allVert.end();
	    //list<mEntity*>::const_iterator lit_begin=cell->allVert.begin();
	    //  ++begin;
	    //for(lit = vcell->allVert.begin();it!=;++it)
	    }*/
      }
    
  //if (theLimiter) computeFullJump();
  
  if(theLimiter) {
    for(it = mesh_begin;it != mesh_end;++it)
      {
	mEntity *m = (*it);
	DGCell *cell = (DGCell*)m->getCell();
	if(cell->getFunctionSpace()->order()>0) theLimiter->limit(cell,0);
      }
   
    for(mMesh::iter it = mesh_begin;it != mesh_end;++it)
      {
	mEntity *m = (*it);
	DGCell *cell = (DGCell*)m->getCell();
	if(cell->getFunctionSpace()->order()>0) cell->theFieldsCoefficients->copyBack();
	//if(cell->getFunctionSpace()->order()>0) cell->theFieldsCoefficients->copyLimitedCoeffs(cell->limit);
      } 
  }

   /*mEntity *M = *mesh_begin;
   DGCell *cell = (DGCell*)M->getCell();
   printf("Size of DGCell %d \n",sizeof(*cell));
   M = *mesh_begin_nm1;
   DGBoundaryCell *Bcell = (DGBoundaryCell*)M->getCell();
   printf("Size of DGBoundaryCell %d \n",sizeof(*Bcell));*/

  //FILE *myLog = fopen("DG.log","w");
  //fclose(myLog);  
}

DGAnalysis::~DGAnalysis () 
{
  delete theIntegrator;
}

void DGAnalysis::run()
{
  int P_ref_counter = 0;
  double dt;
  printf("dg running -- final time = %f...\n",TEND);
  double TOLD = 0.0;
  double cputime = 0.0;
  int sizen   = theMesh->size(n);
  int sizenm1 = theMesh->size(n-1);
  double error; ExactError(TACT,error); printf("Error in density L1 norm %e\n", error);
  while (TACT < TEND)
    {
      clock_t t1 = clock();
      if(HREF && TACT - TOLD > REF_SAMPLING_TIME)
	{
	  TOLD = TACT;
	  adaptH ();
	  // computeNormalsToCurvedBoundaries();
	  sizen   = theMesh->size(n);
	  sizenm1 = theMesh->size(n-1);
	}
    /*  if(PREF && TACT - TOLD > REF_SAMPLING_TIME)
	{
	  TOLD = TACT;
	  adaptP();
 	}*/
	   if(PREF && P_ref_counter>9)
	{
	  adaptP();
	  P_ref_counter=0;
 	}
	   else if (PREF) P_ref_counter++;
      //EVALUATE SENSORS
      evalSensors();
      //exit(0);
      // COMPUTE TIME STEP USING CFL CONDITION --------------
      dt = 500;//arbitrary number to start with; actual dt is computed below
      int N;
      mMesh::iter it;
      mMesh::iter mesh_end=theMesh->end(n);
      
      switch(timeSteppingMode)
	{
	case 0:
	  for(it = theMesh->begin(n);it != mesh_end;++it)
	    {
	      mEntity *m = (*it);
	      DGCell *cell = (DGCell*)m->getCell();
	      cell->adaptTimeStep(CFLMAX,dt);
	    }
	  N = (ceil)((TEND-TACT)/dt);
          dt = (double)(TEND-TACT)/(N);
	  for(it = theMesh->begin(n);it != mesh_end;++it)
	    {
	      DGCell *cell = (DGCell*)(*it)->getCell();
	      cell->setTimeStep(dt);
	    }
	  break;
	case 1:
	  dt=500000.;
	  for(it = theMesh->begin(n);it != mesh_end;++it)
	    {
	      mEntity *m = (*it);
	      DGCell *cell = (DGCell*)m->getCell();
	      cell->computeTimeStep(CFLMAX);
	      double tmp = cell->getTimeStep();
		  if (tmp<dt) dt=tmp;
	    }
	  break;
	case 2:
	  for(it = theMesh->begin(n);it != mesh_end;++it)
	    {
	      mEntity *m = (*it);
	      DGCell *cell = (DGCell*)m->getCell();
	      cell->adaptTimeStep(CFLMAX,dt);
	    } 
	  dt *=pow(2.,theMesh->getNbBins()-1); 
	  N = (ceil)((TEND-TACT)/dt);
          dt = (double)(TEND-TACT)/(N);
	  
	  for(it = theMesh->begin(n);it != mesh_end;++it)
	    {
	      DGCell *cell = (DGCell*)(*it)->getCell();
	      cell->setTimeStep(dt);
	    }
	  break;
	case 3:
	/* for(it = theMesh->begin(n);it != mesh_end;++it)
	    {
	      mEntity *m = (*it);
	      DGCell *cell = (DGCell*)m->getCell();
	      cell->adaptTimeStep(CFLMAX,dt);
	    }
	  N = (ceil)((TEND-TACT)/(dt));
          dt = (double)(TEND-TACT)/(N);
	  for(it = theMesh->begin(n);it != mesh_end;++it)
	    {
	      DGCell *cell = (DGCell*)(*it)->getCell();
	      cell->setTimeStep(dt);
	    }
	  */
	   dt=500000.;
	   
	  for(it = theMesh->begin(n);it != mesh_end;++it)
	    {
	      mEntity *m = (*it);
	      DGCell *cell = (DGCell*)m->getCell();
	      cell->computeTimeStep(CFLMAX);
	      double tmp = cell->getTimeStep();
		  if (tmp<dt) dt=tmp;
	    }
		
	  break;
	}
      
      // SOLVE AND ADVANCE IN TIME ---------------------------
      // THEN LIMIT ------------------------------------------
      double resid = 0.0;
      double residError = 0.0;
      
      resid = theIntegrator->advanceInTime(TACT, dt);
      // if(theError>0) theError->error(cell);
      
      //computeMyError();
      // evaluate all sensors that user has defined
      //  double em, eM; computeError(eM,em)
      //computeBoundaryIntegral(BI); printf("Boundary Integral %e\n",BI);
      
      //residError = 0;	
      
      if (compute_Linf_error) {
	double error;
		  LinfError(TACT+dt,error);
		  printf("Error in density Linf norm %e\n", error);
      }
      if (compute_density_errorL1) {
	double exactError;
		  ExactError(TACT+dt,exactError);
		  printf("Error in density L1 norm %e\n", exactError);
      }
      if (compute_density_errorL2) {
	double exactError;
	densityErrorL2(TACT+dt,exactError);
	printf("Error in density L2 norm %e\n", exactError);
	  }
      if (compute_entropy_errorL1) {
	double entropyError;
	entropyErrorL1(TACT+dt,entropyError); 
	printf(" Entropy error in L1 norm %e\n", entropyError);
      }
      if (compute_entropy_errorL2) {
	double entropyError;
	entropyErrorL2(TACT+dt,entropyError); 
	printf(" Entropy error in L2 norm %e\n", entropyError);
      }
      if (compute_pressure_error) { 
	double pressureError;
	PressureError(TACT+dt,pressureError); 
	printf(" Pressure error in L1 norm %e\n", pressureError);
      }
      if (compute_pressure_errorL2) {
	double pressureError;
	pressureErrorL2(TACT+dt,pressureError);
	printf(" Pressure error in L2 norm %e\n", pressureError);
      }
      if (compute_total_pressure) {
	double BP;
	computeTotalPressure(BP);
	printf(" Boundary Pressure %e\n", BP);
      }
      if (compute_total_mass_loss) {
	double totalMassLoss;
	computeTotalMassLoss(totalMassLoss);
	printf("Total Mass Loss % e\n",totalMassLoss); 
      }
      if (compute_exact_total_pressure) {
	double EBP;
	computeExactTotalPressure(EBP);
	printf(" Exact total Pressure %e\n", EBP);
      }
      if (compute_lift_drag_coeffs) { 
	double C_l,C_d;
	computeLiftDragCoeffs(C_l,C_d);
	printf("lift coeff = %e, drag coeff = %e \n", C_l,C_d);     
      }
      if (plot_pressure_coefficient) plotPressureCoefficient();
      if (plot_mach_on_surface) plotMachOnSurface();
      if (plot_pressure_loss_coefficient) plotPressureLossCoefficient();
      if (plot_entropy_on_surface) plotEntropyOnSurface();	  
      TACT += dt/1.;
      evalSensors();
      //exit(0);
      ITER++;
      clock_t t2 = clock();
      double cpu_t = (double)(t2-t1)/CLOCKS_PER_SEC;
      cputime += cpu_t;
      
      printf("%d %d iter %6d dt %12.5E T %12.5E ||rhs|| =  %12.5E ||Erhs|| =%12.5E cpu %12.5E cpu(ts) = %12.5E cpuleft = %12.5E\n",sizen,sizenm1,ITER,dt,TACT,sqrt(resid),sqrt(residError),cputime,cpu_t,cpu_t*(TEND-TACT)/dt);
      //     FILE *myLog = fopen("DG.log","a");
      //fprintf(myLog,"%d %d iter %6d dt %12.5E T %12.5E ||rhs|| =  %12.5E ||Erhs|| =%12.5E cpu %12.5E cpu(ts) = %12.5E cpuleft = %12.5E\n",sizen,sizenm1,ITER,dt,TACT,sqrt(resid),sqrt(residError),cputime,cpu_t,cpu_t*(TEND-TACT)/dt);
      //fclose(myLog);
      //      if (sqrt(resid)<1.0e-10) break;  
      if(!(ITER % steps_to_dump))   //WRITE DATA TO A FILE
	{
	  char recov[256],mesh[256];
	  sprintf(recov,"dg-recov-%d.dat",ITER);
	  sprintf(mesh,"mesh-recov-%d.msh",ITER);
	  ofstream data(recov);
	  ofstream m(mesh);
	  write(m,data);
	  // then write infos about the actual state
	  //par->write(ofs);
	  data.close();
	  m.close();
	}
      //ApproxError(approxError);
      //computeBoundaryIntegral();
    } 
  //WRITE SOLUTION AT FINAL TIME TO FILE
  char recov[256],mesh[256];
  sprintf(recov,"final-solution.dat");
  sprintf(mesh,"final-mesh.msh");
  ofstream data(recov);
  ofstream m(mesh);
  write(m,data);
  // then write infos about the actual state
  //par->write(ofs);
  data.close();
  m.close();
  
  ExactError(TACT,error);
  // evalSensors();
  printf(" Exact error %e \n", error); 
}

void DGAnalysis::sortCellsBySize()
{
  mMesh::iter it;
  mMesh::iter mesh_end=theMesh->end(n);
  mMesh::iter mesh_end_nm1=theMesh->end(n-1);
  int i;
  double cellSize;
  double maxCellSize = 0.0;
  double minCellSize = 1.0;
  
  for(it = theMesh->begin(n);it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      cellSize = cell->getSize();
      if (maxCellSize<cellSize) maxCellSize = cellSize;
      if (minCellSize>cellSize) minCellSize = cellSize;
      //printf("size %f\n",cellSize);
    }
  
  double ratio = maxCellSize/minCellSize;
  int NbBins = ceil(log(ratio)/log(2.));
  theMesh->setNbBins(NbBins);
  double* BB = new double [NbBins+1];
  BB[0] = maxCellSize;
  int* count = new int [NbBins];
  for (i=0; i<NbBins; i++)
    count[i]=0;
  
  //sort elements to bins according to their size
  if (NbBins >1)
    {
      for (i=1; i<=NbBins; i++)
	BB[i] = BB[i-1]*0.5;
      
      for(it = theMesh->begin(n);it != mesh_end;++it)
	{
	  mEntity *m = (*it);
	  DGCell *cell = (DGCell*)m->getCell();
	  cellSize = cell->getSize();
	  for (i=1; i<=NbBins; i++)
	    if (cellSize<=BB[i-1]*1.0001 &&cellSize>BB[i]*1.0001)
	      {
		theMesh->addToBin(i-1,m);
		cell->setBinNb(i-1);
		cell->limit=i-1;
		//printf("cell %f in %f %f \n",cellSize, BB[i], BB[i-1]);
		count[i-1]++;
		//	if (i>1) m->print();
	      }
	}
    }
  // sort boundaries into bins
  for(it = theMesh->begin(n-1);it != mesh_end_nm1;++it)
    {
      mEntity *m = (*it);
      //m->print();
      if (m->getId()==269)
	{
	  int s=1;
	}
      DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
      DGCell* leftCell = (DGCell* ) cell->getLeft()->getCell();
      mEntity *left = cell->getLeft();
      mEntity *right = cell->getRight();
      //cell->getLeft()->print();
      //if (right) right->print();
      //printf("\n");
      DGCell* rightCell;
      if (right) rightCell = (DGCell* ) right->getCell();
      int leftNb = leftCell->getBinNb();
      int rightNb=0;
      if (right) rightNb = rightCell->getBinNb();
      //		  if ((leftNb == rightNb || !right)&&leftNb==0) m->print();	
      if (leftNb == rightNb || !right) theMesh->addToBin(leftNb,m);			//inner edges
      else {
	int diff = (leftNb>rightNb ? leftNb-rightNb: rightNb-leftNb);     //interface edges
	if (diff == 1) {
	  int minNb = (leftNb<rightNb ? leftNb : rightNb); 
	  theMesh->addToInterfaceBin(minNb,m);
	  if (leftNb<rightNb) 
	    {
	      theMesh->addToLargeOnInterface(leftNb,left);
	      theMesh->addToSmallOnInterface(leftNb,right);
	    }
	  else
	    {
	      theMesh->addToLargeOnInterface(rightNb,right);
	      theMesh->addToSmallOnInterface(rightNb,left);
	    }
	  theMesh->delFromBin(leftNb,left);
	  theMesh->delFromBin(rightNb,right);
	  
	  //m->print();
	}
	else {printf("Too big a difference between elements\n"); exit(0);}
      } 
    }
  /*Checking if number of cells and edges in bins == cells & edges  in the mesh*/
  
  int sizeInMesh = theMesh->size(1);  //total number of edges
  int sizetotal = 0; //total number of edges in bins
  for (i=0;i<NbBins;i++)   //inner edges
    {
      int size2 = theMesh->sizeOfBin(i,1);
      sizetotal+=size2;
      printf("Number of edges in bin %d = %d\n",i,size2);
      //printf("count %d\n",count[i]);
    }
  for (i=0;i<NbBins-1;i++)   //interface edges
    {
      int size2 = theMesh->sizeOfInterfaceBin(i,1);
      sizetotal+=size2;
      printf("Number of edges in interface bin %d = %d\n",i,size2);
    }
  if (sizeInMesh==sizetotal) printf("Total number of edges %d\n",sizeInMesh); else
    {printf("Edges were lost in sorting!");exit(0);}
  
  sizeInMesh = theMesh->size(2);  //total number of cells
  sizetotal = 0; //total number of cells in bins
  for (i=0;i<NbBins;i++)   //inner edges
    {
      int size2 = theMesh->sizeOfBin(i,2);
      sizetotal+=size2;
      printf("Number of cells in bin %d = %d\n",i,size2);
      //printf("count %d\n",count[i]);
    }
  for (i=0;i<NbBins;i++)
    {
      int size2 = theMesh->sizeOfLargeOnInterface(i,2);
      sizetotal+=size2;
      printf("elements in Large on Interface %d %d\n",i, size2);
      size2 = theMesh->sizeOfSmallOnInterface(i,2);
      sizetotal+=size2;
      printf("elements in Small on Interface %d %d\n",i,size2);
    }
  if (sizeInMesh==sizetotal) printf("Total number of cells %d\n",sizeInMesh); else
    {printf("Cells were lost in sorting!");exit(0);}
  
  delete [] BB;
  delete [] count;
}













