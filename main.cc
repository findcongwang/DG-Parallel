#include "mImportExport.h"
#include "mDGMesh.h"
#include "mEntity.h"
#include "parallel.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <map>
#include <cstdlib>
#include <unistd.h>
#include <vector>

extern void ShockTube (mDGMesh *theMesh, int);
extern void OutFlow   (mDGMesh *theMesh, double mach, double angle, int wall, int);
extern void DoubleMach (mDGMesh *theMesh, int wall, int);
extern void ScalarPb  (mDGMesh *theMesh, int);
extern void Suli (mDGMesh *theMesh, int);
extern void RayleighTaylor  (mDGMesh *theMesh, int, int);
extern void Burger (mDGMesh *theMesh,int);
extern void Burger2D (mDGMesh *theMesh,int);
extern void SuperVortex (mDGMesh *theMesh, int wall, int);
extern void Cylinder (mDGMesh *theMesh, int wall, int);
extern void Sphere (mDGMesh *theMesh, int wall, int);
extern void FlowAroundAirfoil (mDGMesh *theMesh, int wall, int);
extern void LinSystem (mDGMesh *theMesh, int);
extern void ForwardFacingStep(mDGMesh *theMesh, int wall, int);
extern void SphericalBlast (mDGMesh *theMesh,int wall, int order);
extern void Riemann2D (mDGMesh *theMesh,int wall, int order);
extern void TalkExample (mDGMesh *theMesh,int wall, int order);
extern void DaleExample (mDGMesh *theMesh, int wall, int);

unsigned int numCPU;
unsigned int numThreads;
vector<infoWrapper> threadInfo[3];

int main(int argc, char *argv[])
{

  numCPU = sysconf( _SC_NPROCESSORS_ONLN );
  numThreads = numCPU;
  //printf("%d CPUs online\n", numCPU);

  mImportExport io;
  printf("reading the mesh ...\n");
  clock_t t1 = clock();
  mDGMesh *theMesh = new mDGMesh;
  io.import(argv[1],theMesh);
  clock_t t2 = clock();
  printf("Mesh read in %f seconds\n", double (t2-t1)/CLOCKS_PER_SEC);
  printf("creating connectivities ...\n");

  switch (atoi(argv[2]))
    {
    case 3:
      printf("mesh read. %f seconds needed\n",(double)(t2-t1)/CLOCKS_PER_SEC);
      t1 = clock();
      getchar();
      theMesh->createCellBoundaries(3,2);
      t2 = clock();
      printf("%d faces created %f seconds needed\n",theMesh->size(2),(double)(t2-t1)/CLOCKS_PER_SEC);
      getchar();
      t1 = clock();
      theMesh->createCellBoundaries(2,1);
      t2 = clock();
      printf("%d edges created %f seconds needed\n",theMesh->size(1),(double)(t2-t1)/CLOCKS_PER_SEC);
      getchar();
      t1 = clock();

	  int i;
      for(i=0;i<3;i++)
	{
	  int key = -1;
	  int NK = 0;
	  for(mMesh::iter it = theMesh->begin(0);it != theMesh->end(0); it++)
	    {
	      mEntity *e = *it;
	      if(key != e->getId())
		{
		  NK++;
		}
	    }
	  printf("PHI(%d) = %12.5E\n",i,(double)(theMesh->size(i))/(double)(NK));
	  NK = 0;
	}
      break;
    default:
      int dim = (theMesh->size(3) !=0)?3:2;
   
	  t1=clock();
      theMesh->createCellBoundaries(dim,dim-1); //creating edges
	  t2 = clock();
      printf("Edges created in %e seconds\n", double (t2-t1)/CLOCKS_PER_SEC);
      t1=clock();
	  theMesh->createConnections(dim-1,dim); //create poiters to cells from cell boundaries
	  //      theMesh->createConnections(0,dim); //connect faces to vertices  // needed for VertLimiter
	  // needed for reconstructing curved boundaries:
	  theMesh->createConnections(0,dim-1);
	  //if (dim==2) theMesh->createConnections(0,1);  //create pointers to edges from vertices   
     // if (dim==3) theMesh->createConnections(1,2);   // create pointers to faces from edges 
	  t2 = clock();
      printf("Connections created in %e seconds\n", double (t2-t1)/CLOCKS_PER_SEC);

       theMesh->setPeriodicBC(dim);
	  printf("%d cells \n", theMesh->size(dim));
      printf("%d boundaries ...\n",theMesh->size(dim-1));
      break;
    }

  printf("DG begins ...\n");
  
  switch (atoi(argv[2]))
    {
    case 1:
      ShockTube(theMesh,atoi(argv[3]));
      break;
    case 2:
      OutFlow (theMesh,3.0,0.0,200,atoi(argv[3]));
      break;
	case 9:
      LinSystem(theMesh,atoi(argv[3]));
      break;
    case 8:
      Suli(theMesh,atoi(argv[3]));
      break;
    case 4:
      ScalarPb (theMesh,atoi(argv[3]));
      break;
    case 5:
      RayleighTaylor (theMesh,200,atoi(argv[3]));
      break;
    case 6:
      Burger(theMesh,atoi(argv[3]));
      break;
    case 7:
      Burger2D(theMesh,atoi(argv[3]));
      break;
    case 12:
      DoubleMach(theMesh,200,atoi(argv[3]));
      break;
    case 13:
      SuperVortex(theMesh,200,atoi(argv[3]));
      break; 
    case 14:
      Cylinder(theMesh,200,atoi(argv[3]));
      break; 
	case 15:
      FlowAroundAirfoil(theMesh,200,atoi(argv[3]));
      break; 
    case 16:
      ForwardFacingStep(theMesh,200,atoi(argv[3]));
      break;
    case 17:
      SphericalBlast(theMesh,200,atoi(argv[3]));
      break;
	case 18:
      Riemann2D(theMesh,200,atoi(argv[3]));
      break;
	case 19:
      TalkExample(theMesh,200,atoi(argv[3]));
      break;
    case 20:
      Sphere(theMesh,200,atoi(argv[3]));
      break;
	case 88:
		DaleExample(theMesh,200,atoi(argv[3]));
		break;
    }
  return 0;
}

