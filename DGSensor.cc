#include "DGSensor.h"
#include "DGCell.h"
#include "mDGMesh.h"
#include "DGAnalysis.h"
#include "mEntity.h"
#include "ConservationLaw.h"
#include "Mapping.h"
#include "FunctionSpace.h"
#include "postProGmsh.h"
#include <stdio.h>
#include <iostream>
#include <cstring>
//#include <ostream.h>
#include "Constants.h"
using namespace std;

DGSensor::DGSensor(char *name, double DT)
  :  DTsampleFrequency(DT),TOLD(0.0)
{
#ifdef PARALLEL
  char np[256];
  sprintf(np,"%s-proc%d",name,Parutil.rank());
  fileName = new char [strlen(np)+1];
  strcpy(fileName,np);
#else
  fileName = new char [strlen(name)+1];
  strcpy(fileName,name);
#endif
}

DGSensor::~DGSensor()
{
  delete [] fileName;
}

void DGAnalysis::addSensor (DGSensor *s)
{
  allSensors.push_back(s);
}


DGLineSensor::DGLineSensor(char *name, const mPoint &p, const mPoint &q, int n, double dt)
  : DGSensor(name,dt),p1(p),p2(q),nbSamples(n)
{
  FILE *f = fopen (name,"w");
  if(f)fclose(f);
}

DGDataViseSensor::DGDataViseSensor(char *name, double dt)
  : DGSensor(name,dt),ITER(0)
{
}

DGPointSensor::DGPointSensor(char *name, const mPoint &p,  double dt)
  : DGSensor(name,dt),p1(p)
{
  FILE *f = fopen (name,"w");
  fprintf(f,"\n# evaluation of a point sensor\n");
  if(f)fclose(f);
}

DGTresholdLineSensor::DGTresholdLineSensor(char *name, int i, const mPoint &pp1,
					   const mPoint &pp2,  int s, double t)
  : DGSensor(name,0.0),ithValueOfInterest(i),p1(pp1),p2(pp2),nbSamples(s), treshold(t)
{
  FILE *f = fopen (fileName,"w");
  fprintf(f,"\n# evaluation of a treshold line sensor\n");
  if(f)fclose(f);
}

void DGAnalysis::evalSensors(bool doNow)
{
  for(list<DGSensor*>::const_iterator it = allSensors.begin(); it != allSensors.end();++it)
    {
      DGSensor *sens = *it;
      sens->eval(this,n,TACT,doNow);
    }
}


// only add the sensor to the geom search when
// needed !
void DGLineSensor::eval(DGAnalysis *theAnalysis, 
			int n, double T, bool doNow)
{
  if(T - TOLD <  DTsampleFrequency && T!=0 && !doNow) return;
 
  double val[MaxNbEqn], valofinterest[MaxNbEqn];
  TOLD = T;
//  TOLD = DTsampleFrequency*ITER;

  FILE *f = fopen (fileName,"w");
  fprintf(f,"\n# evaluation of a line sensor with %d samples at time %12.5E at  %d\n"
	  ,nbSamples,T);
  for(int i=0;i<nbSamples;i++)
    {
      mPoint pp = (p2-p1);
      mPoint p = p1 + (pp) * ((double)i/(double)(nbSamples-1));
      if(theAnalysis->eval(p,val))
	{	  
	  fprintf(f,"%d %12.5E %12.5E %12.5E %12.5E ",i+1,T,p(0),p(1),p(2));
	  for(int k=0;k<theAnalysis->getLaw()->getNbFieldsOfInterest();k++)
	    {
	      theAnalysis->getLaw()->getIthFieldOfInterest(k,p,val,valofinterest,T);
	      fprintf(f,"%12.5E ",valofinterest[0]);
	    }
	  fprintf(f,"\n");
	}
    }
  fclose(f);
}

void DGTresholdLineSensor::eval(DGAnalysis *theAnalysis,
				int n,  double T, bool doNow)
{
  double val1[256],val2[256];

  FILE *f = fopen (fileName,"a");

  if(!f)return;

  int k = 0;

  //  printf("eval %d %d\n",nbSamples,ithValueOfInterest);

  for(int i=0;i<nbSamples;i++)
    {
      mPoint pp = (p2-p1);
      mPoint p = p1 + (pp) * ((double)i/(double)(nbSamples-1));
      
      if(theAnalysis->eval(p,val1))
	{
	  //  if(k)printf("%f %f %f\n",val1[0],val2[0],treshold);
	  if(k++)
	    {
	      if((val1[ithValueOfInterest] > treshold && val2[ithValueOfInterest] < treshold)||
		 (val1[ithValueOfInterest] < treshold && val2[ithValueOfInterest] > treshold))
		{
		  {
		    fprintf(f,"%12.5E %12.5E %12.5E %12.5E\n",T,p(0),p(1),p(2));
		    fclose(f);
		    return;
		  }	     
		}
	    }
	  for(int l=0;l<theAnalysis->getLaw()->getNbFields();l++)
	    {
	      val2[l] = val1[l];
	    }
	}
    }
  fclose(f);
}

// only add the sensor to the geom search when
// needed !
void DGPointSensor::eval(DGAnalysis *theAnalysis,
			 int n, double T, bool doNow)
{
  double val[MaxNbEqn];
 
  FILE *f = fopen (fileName,"a");
  if(theAnalysis->eval(p1,val))
    {
      fprintf(f," %12.5E ",T);
      for(int k=0;k<theAnalysis->getLaw()->getNbFields();k++)
	{
	  fprintf(f,"%12.5E ",val[k]);
	}
      fprintf(f,"\n");
    }
  fclose(f);
}

/*

GMSH Output sensor

*/


DGGmshSensor::DGGmshSensor(char *name, double dt) 
  : DGSensor (name,dt), ITER(0)
{
   //triangle vertices
  u[0][0] = u[0][1] = 0.0; u[0][2] = 1.0;
  v[0][0] = v[0][2] = 0.0; v[0][1] = 1.0;
  //quad split in 4
  uu[0][0] = uu[0][1] = -1.0; uu[0][2] = 0.0;
  vv[0][0] = -1.0 ; vv[0][1] = 1.0; vv[0][2] = 0.0;
  uu[1][0] = -1.0 ; uu[1][1] =  1.0; uu[1][2] = 0.0;
  v[1][0] = v[1][1] =  1.0; v[1][2] = 0.0;
  uu[2][0] = uu[2][1] =  1.0; uu[2][2] = 0.0;
  vv[2][0] = 1.0; vv[2][1] = -1.0; vv[2][2] = 0.0;
  uu[3][0] = -1.0 ; uu[3][1] =  1.0; uu[3][2] = 0.0;
  vv[3][0] = vv[3][1] =  -1.0; vv[3][2] = 0.0;
  	
}

SummationSensor::SummationSensor(char *name, double dt) 
  : DGSensor (name,dt), ITER(0), COUNTER(0)
{}

extern void PascalgetIndices(int iFct, int &n, int &i);
class splitTri
{
  int level;
  mPoint t[3];
  int thisTri,iF,nF;
public :
  splitTri(int l, mPoint p1, mPoint p2, mPoint p3):level(l)
  {
    thisTri = 0;
    iF = nF = 0;
    t[0] = p1;t[1] = p2;t[2] = p3;
  }
  int nbTri (){return (level)*(level);}

  bool nextTri(mPoint &pt1, mPoint &pt2, mPoint &pt3)
  {
    int i=0,n=0;
    if(nF++ == level*level)return false;
    PascalgetIndices(iF,i,n);
    if(i == n || thisTri == 1)
      {
	iF++;
	thisTri = 0;
      }
    else
      {
	thisTri++;
      }

    double x1 = (double)i/(double)(level);
    double y1 = (double)n/(double)(level);
    double x2 = (double)(i+1)/(double)(level);
    double y2 = (double)(n+1)/(double)(level);
    if(thisTri == 0)
      {
	pt1 = t[1] + ((t[0]-t[1]) * x1) + ((t[2]-t[0]) * y1);
	pt2 = t[1] + ((t[0]-t[1]) * x2) + ((t[2]-t[0]) * y1);
	pt3 = t[1] + ((t[0]-t[1]) * x2) + ((t[2]-t[0]) * y2);
      }
    else
      {
	pt1 = t[1] + ((t[0]-t[1]) * x1) + ((t[2]-t[0]) * y1);
	pt2 = t[1] + ((t[0]-t[1]) * x1) + ((t[2]-t[0]) * y2);
	pt3 = t[1] + ((t[0]-t[1]) * x2) + ((t[2]-t[0]) * y2);
      }
    /*     printf(" tri %d %d %d in (%f %f) (%f %f) (%f %f) = (%f %f) (%f %f) (%f %f)\n",
	     iF,i,n,t[0](0),t[0](1),t[1](0),t[1](1),t[2](0),t[2](1),
	     pt1(0),pt1(1),pt2(0),pt2(1),pt3(0),pt3(1));
    */
     return true;
  }
};

void DGGmshSensor::eval(DGAnalysis *theAnalysis,
			int dim,  double T, bool doNow)
{
  if(T != 0.0 && T - TOLD <  DTsampleFrequency && !doNow) return;
  //TOLD = T;
  TOLD = DTsampleFrequency*ITER;

  ITER++;
  theAnalysis->exportGmsh(fileName, ITER);
  /*char name[256], field[256];
  for(int i=0;i<theAnalysis->getLaw()->getNbFieldsOfInterest();i++)
    {
      int k;
      theAnalysis->getLaw()->getNameAndSize(i, k, field);
      sprintf(name,"%s-sample%d",fileName,ITER);
      strcat(name,field);
      strcat(name,".pos");
      ofstream out(name);
      theAnalysis->exportGmsh(out,i);
      out.close();
      }*/
}

void DGDataViseSensor::eval(DGAnalysis *theAnalysis, int n, double T, bool doNow)
{
  if(T != 0.0 && T - TOLD <  DTsampleFrequency && !doNow) return;
  
  TOLD = T;
  ITER++;  
  char name[56];
  sprintf(name,"%s-sample%d.inf",fileName,ITER);
  
  FILE *f = fopen(name,"w");
  ofstream out(name);
  theAnalysis->exportDataVise(out);
  out.close();
}

void DGAnalysis::exportDataVise(ofstream &f)
{
/*
  double val[256],intr[256];
  int N = 0;
  mMesh::iter it;
  for( it= theMesh->begin(n);it != theMesh->end(n);++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      N += cell->nbPtGauss;
    }

  f << "0 -0.5 0\n";
  f << "0.25  0.5 0.25\n";
  f << N << " " << getLaw()->getNbFieldsOfInterest() << "\n";

  for(it = theMesh->begin(n);it != theMesh->end(n);++it)
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();
      for(int i=0;i<cell->nbPtGauss;i++)
	{
	  vptGaussInfo &pg = cell->pt(i);
	  cell->interpolate( pg.u,pg.v,pg.w, val,&pg.fcts);
	  f << pg.x << " " << pg.y << " " << pg.z << " " ;
	  for(int j=0;j<getLaw()->getNbFieldsOfInterest();j++)
	    {
	      getLaw()->getIthFieldOfInterest(j, mPoint(pg.x,pg.y,pg.z),val, intr,TACT);
	      f << intr[0] << " ";
	    }
	  f << "\n";
	}
    }*/
}

void DGAnalysis::exportGmsh(char *fileName, int iteration)
{
  int dim = n;
  if(dim == 3) {exportGmsh_Volume(fileName,iteration); return;}
  double u[4][3],v[4][3];
  
  int s = 0;
  int NN = 0;  //total number of triangles used in plotting
  int nbtris,nbspl;
  char field[256], name[256];
  mMesh::iter it;
  mMesh::iter mesh_begin  = theMesh->begin(2);
  mMesh::iter mesh_end = theMesh->end(2);

  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      if(dim == 2 || m->getClassification()->getLevel() == 2
	 || m->size(3) == 1)
	{
	  DGCell *cell;
	  DGBoundaryCell *bcell;
	  if(dim == 3)
	    {
	      bcell = (DGBoundaryCell*)m->getCell(); 
	      cell = (DGCell*)bcell->getLeft()->getCell();
	    }
	  else cell = (DGCell*)m->getCell();
	  if(m->getType() == mEntity::TRI)  nbtris = 1;      //one triangle if an element is a triangle
	  else if(m->getType() == mEntity::QUAD) nbtris = 4; //split quad into 4 triangles
	  
	  nbspl = cell->getFunctionSpace()->order();
	  if(nbspl == 0) nbspl = 1;
	  NN += nbtris * nbspl * nbspl;
	}
    }
  
  int NbVec=0;
  int NbFieldsOfInterest = getLaw()->getNbFieldsOfInterest();
  for (int q=0;q<NbFieldsOfInterest;q++)
    {
      int k;
      getLaw()->getNameAndSize(q, k, field);
      if (k==3) NbVec++;
    }
  
  PostProGmshView view(NN);
  PostProGmshView viewVector(NN*NbVec);
  PostProGmshShape * shape,*shapeVector;
  
  mPoint x1,x2,x3;
  double val1[MaxNbEqn],val2[MaxNbEqn],val3[MaxNbEqn];
  double intr1[3],intr2[3],intr3[3];
  
  
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      if(dim == 2 || m->getClassification()->getLevel() == 2
	 || m->size(3) == 1)
	{
	  DGCell *cell;
	  DGBoundaryCell *bcell;
	  if(dim == 3)
	    {
	      bcell = (DGBoundaryCell*)m->getCell(); 
	      cell = (DGCell*)bcell->getLeft()->getCell();
	    }
	  else cell = (DGCell*)m->getCell();
	  
	  
	  nbspl = cell->getFunctionSpace()->order();
	  //nbspl = nbspl*nbspl;
	  if(nbspl == 0) nbspl = 1;
	  
	  if(m->getType() == mEntity::TRI)                    //higher order code needed here?
	    {
	      nbtris = 1;
	      //triangle vertices
	      u[0][0] = u[0][1] = 0.0; u[0][2] = 1.0;
	      v[0][0] = v[0][2] = 0.0; v[0][1] = 1.0;  
	    }
	  else if(m->getType() == mEntity::QUAD)
	    {
	      nbtris = 4;
	      //quad split in 4
	      u[0][0] = u[0][1] = -1.0; u[0][2] = 0.0;
	      v[0][0] = -1.0 ; v[0][1] = 1.0; v[0][2] = 0.0;
	      u[1][0] = -1.0 ; u[1][1] =  1.0; u[1][2] = 0.0;
	      v[1][0] = v[1][1] =  1.0; v[1][2] = 0.0;
	      u[2][0] = u[2][1] =  1.0; u[2][2] = 0.0;
	      v[2][0] = 1.0; v[2][1] = -1.0; v[2][2] = 0.0;
	      u[3][0] = -1.0 ; u[3][1] =  1.0; u[3][2] = 0.0;
	      v[3][0] = v[3][1] =  -1.0; v[3][2] = 0.0; 
	    }
	  
	  mPoint p1,p2,p3;
	  
	  for(int j = 0;j<nbtris;j++)
	    {
	      splitTri tt(nbspl, mPoint (u[j][0],v[j][0],0.0),
			  mPoint (u[j][1],v[j][1],0.0),
			  mPoint (u[j][2],v[j][2],0.0));
	      
	      
	      while(tt.nextTri(p1,p2,p3))
		{
		  
		  if(dim == 3)
		    {
		      bcell->getMapping()->eval(p1(0),p1(1),0,x1(0),x1(1),x1(2));
		      cell->getMapping()->invert(x1(0),x1(1),x1(2),p1(0),p1(1),p1(2));
		      bcell->getMapping()->eval(p2(0),p2(1),0,x2(0),x2(1),x2(2));
		      cell->getMapping()->invert(x2(0),x2(1),x2(2),p2(0),p2(1),p2(2));
		      bcell->getMapping()->eval(p3(0),p3(1),0,x3(0),x3(1),x3(2));
		      cell->getMapping()->invert(x3(0),x3(1),x3(2),p3(0),p3(1),p3(2));
		    }
		  else
		    {
		      cell->getMapping()->eval(p1(0),p1(1),0.0,x1(0),x1(1),x1(2));
		      cell->getMapping()->eval(p2(0),p2(1),0.0,x2(0),x2(1),x2(2));
		      cell->getMapping()->eval(p3(0),p3(1),0.0,x3(0),x3(1),x3(2));
		    }
		  cell->interpolate(p1(0),p1(1),p1(2),val1);
		  cell->interpolate(p2(0),p2(1),p2(2),val2);
		  cell->interpolate(p3(0),p3(1),p3(2),val3);
		  
		  
		  view.NextInViewIsAScalarTriangle(NbFieldsOfInterest);
		  shape = view.getShape(s);
		  if (NbVec)
		  {
			  viewVector.NextInViewIsAVectorTriangle(NbVec);
			  shapeVector = viewVector.getShape(s);
		  }
		 s++; // s is a counter over all plotting triangles 
		  for (int q=0;q<NbFieldsOfInterest;q++)
		    {
			  int k;
			  getLaw()->getNameAndSize(q, k, field);
		      if (!strcmp(field,"limit")) 
			{
			  //cell->interpolateError(p1(0),p1(1),p1(2),val1);
			  //cell->interpolateError(p2(0),p2(1),p2(2),val2);
			  //cell->interpolateError(p3(0),p3(1),p3(2),val3);
			  val1[0] = val2[0] = val3[0] = cell->limitStatus();
			  }
			 if (!strcmp(field,"order")) 
			{	
			  val1[0] = val2[0] = val3[0] = (double ) cell->getFunctionSpace()->order();
			}
			
		      theLaw->getIthFieldOfInterest(q, x1, val1, intr1, TACT);
		      theLaw->getIthFieldOfInterest(q, x2, val2, intr2, TACT);
		      theLaw->getIthFieldOfInterest(q, x3, val3, intr3, TACT);

			

			  if (k==3)
			  {
		       shapeVector->setValue(0,intr1[0],intr1[1],intr1[2], 
						intr2[0],intr2[1],intr2[2], 
						intr3[0],intr3[1],intr3[2]);
			  shapeVector->setXYZ(x1(0),x1(1),x1(2),x2(0),x2(1),x2(2),x3(0),x3(1),x3(2));
			  }
			  else
			  {
			  shape->setValue(q,intr1[0],intr2[0],intr3[0]);
			  shape->setXYZ(x1(0),x1(1),x1(2),x2(0),x2(1),x2(2),x3(0),x3(1),x3(2));
			  }
		    }
		  
		}
	    }
	}
    }
  for(int i=0;i<NbFieldsOfInterest;i++)
    {
      int k;
      getLaw()->getNameAndSize(i, k, field);
      {
	sprintf(name,"%s-sample%d",fileName,iteration);
	strcat(name,field);
	strcat(name,".pos");
	ofstream out(name);
	out << "View \"Exported field\" {" << endl; 
	
	if (k==3) for(int j=0;j<viewVector.getCurrentShape();j++) viewVector.AllShapes[j]->print(out,0);
	else for(int j=0;j<view.getCurrentShape();j++) view.AllShapes[j]->print(out,i);
	out<< "};" << endl;
	out.close();
      }
    }
}


void SummationSensor::eval(DGAnalysis *theAnalysis,
			int dim,  double T, bool doNow)
{
  if (!ITER)  {theAnalysis->setupSummation();}
  else theAnalysis->doSummation();COUNTER++;

   if(T != 0.0 && T - TOLD <  DTsampleFrequency && !doNow) return;
  //TOLD = T;
  TOLD = DTsampleFrequency*ITER;
  ITER++;

  theAnalysis->exportSummation(fileName, ITER, COUNTER);
}

void DGAnalysis::setupSummation()
{
  int dim = n;
 
  double u[4][3],v[4][3];
  
  int s = 0;
  int NN = 0;  //total number of triangles used in plotting
  int nbtris,nbspl;
  mMesh::iter it;
  mMesh::iter mesh_begin  = theMesh->begin(2);
  mMesh::iter mesh_end = theMesh->end(2);

  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      if(dim == 2 || m->getClassification()->getLevel() == 2
	 || m->size(3) == 1)
	{
	  DGCell *cell;
	  DGBoundaryCell *bcell;
	  if(dim == 3)
	    {
	      bcell = (DGBoundaryCell*)m->getCell(); 
	      cell = (DGCell*)bcell->getLeft()->getCell();
	    }
	  else cell = (DGCell*)m->getCell();
	  if(m->getType() == mEntity::TRI)  nbtris = 1;      //one triangle if an element is a triangle
	  else if(m->getType() == mEntity::QUAD) nbtris = 4; //split quad into 4 triangles
	  
	  nbspl = cell->getFunctionSpace()->order();
	  if(nbspl == 0) nbspl = 1;
	  NN += nbtris * nbspl * nbspl;
	}
    }
  
  int NbFieldsOfInterest = 1;  
  
  view = new PostProGmshView(NN);
  PostProGmshShape * shape;
  
  mPoint x1,x2,x3;
  
  
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      if(dim == 2 || m->getClassification()->getLevel() == 2
	 || m->size(3) == 1)
	{
	  DGCell *cell;
	  DGBoundaryCell *bcell;
	  if(dim == 3)
	    {
	      bcell = (DGBoundaryCell*)m->getCell(); 
	      cell = (DGCell*)bcell->getLeft()->getCell();
	    }
	  else cell = (DGCell*)m->getCell();
	  
	  
	  nbspl = cell->getFunctionSpace()->order();
	  //nbspl = nbspl*nbspl;
	  if(nbspl == 0) nbspl = 1;
	  
	  if(m->getType() == mEntity::TRI)
	    {
	      nbtris = 1;
	      //triangle vertices
	      u[0][0] = u[0][1] = 0.0; u[0][2] = 1.0;
	      v[0][0] = v[0][2] = 0.0; v[0][1] = 1.0;  
	    }
	  else if(m->getType() == mEntity::QUAD)
	    {
	      nbtris = 4;
	      //quad split in 4
	      u[0][0] = u[0][1] = -1.0; u[0][2] = 0.0;
	      v[0][0] = -1.0 ; v[0][1] = 1.0; v[0][2] = 0.0;
	      u[1][0] = -1.0 ; u[1][1] =  1.0; u[1][2] = 0.0;
	      v[1][0] = v[1][1] =  1.0; v[1][2] = 0.0;
	      u[2][0] = u[2][1] =  1.0; u[2][2] = 0.0;
	      v[2][0] = 1.0; v[2][1] = -1.0; v[2][2] = 0.0;
	      u[3][0] = -1.0 ; u[3][1] =  1.0; u[3][2] = 0.0;
	      v[3][0] = v[3][1] =  -1.0; v[3][2] = 0.0; 
	    }
	  
	  mPoint p1,p2,p3;
	  
	  for(int j = 0;j<nbtris;j++)
	    {
	      splitTri tt(nbspl, mPoint (u[j][0],v[j][0],0.0),
			  mPoint (u[j][1],v[j][1],0.0),
			  mPoint (u[j][2],v[j][2],0.0));
	      
	      
	      while(tt.nextTri(p1,p2,p3))
		{
		  
		  if(dim == 3)
		    {
		      bcell->getMapping()->eval(p1(0),p1(1),0,x1(0),x1(1),x1(2));
		      cell->getMapping()->invert(x1(0),x1(1),x1(2),p1(0),p1(1),p1(2));
		      bcell->getMapping()->eval(p2(0),p2(1),0,x2(0),x2(1),x2(2));
		      cell->getMapping()->invert(x2(0),x2(1),x2(2),p2(0),p2(1),p2(2));
		      bcell->getMapping()->eval(p3(0),p3(1),0,x3(0),x3(1),x3(2));
		      cell->getMapping()->invert(x3(0),x3(1),x3(2),p3(0),p3(1),p3(2));
		    }
		  else
		    {
		      cell->getMapping()->eval(p1(0),p1(1),0.0,x1(0),x1(1),x1(2));
		      cell->getMapping()->eval(p2(0),p2(1),0.0,x2(0),x2(1),x2(2));
		      cell->getMapping()->eval(p3(0),p3(1),0.0,x3(0),x3(1),x3(2));
		    }

		  view->NextInViewIsAScalarTriangle(NbFieldsOfInterest);
		  shape = view->getShape(s);
		
		 s++; // s is a counter over all plotting triangles 
		  
			 
			  shape->setValue(0,0.0,0.0,0.0);
			  shape->setXYZ(x1(0),x1(1),x1(2),x2(0),x2(1),x2(2),x3(0),x3(1),x3(2));
		    }
		  
	  }
	}
    }
}

void DGAnalysis::doSummation()
{
  int dim = n;
  double u[4][3],v[4][3];
  
  int s = 0;
  int NN = 0;  //total number of triangles used in plotting
  int nbtris,nbspl;
  int NbFieldsOfInterest=1;
  char field[256];
  mMesh::iter it;
  mMesh::iter mesh_begin  = theMesh->begin(2);
  mMesh::iter mesh_end = theMesh->end(2);

  PostProGmshShape * shape;
  
  double val1[MaxNbEqn],val2[MaxNbEqn],val3[MaxNbEqn];
  double intr1[3],intr2[3],intr3[3];
  mPoint dummy;
  	  
	  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      if(dim == 2 || m->getClassification()->getLevel() == 2
	 || m->size(3) == 1)
	{
	  DGCell *cell;
	  DGBoundaryCell *bcell;
	  if(dim == 3)
	    {
	      bcell = (DGBoundaryCell*)m->getCell(); 
	      cell = (DGCell*)bcell->getLeft()->getCell();
	    }
	  else cell = (DGCell*)m->getCell();
	  
	  
	  nbspl = cell->getFunctionSpace()->order();
	  //nbspl = nbspl*nbspl;
	  if(nbspl == 0) nbspl = 1;
	  
	  if(m->getType() == mEntity::TRI)
	    {
	      nbtris = 1;
	      //triangle vertices
	      u[0][0] = u[0][1] = 0.0; u[0][2] = 1.0;
	      v[0][0] = v[0][2] = 0.0; v[0][1] = 1.0;  
	    }
	  else if(m->getType() == mEntity::QUAD)
	    {
	      nbtris = 4;
	      //quad split in 4
	      u[0][0] = u[0][1] = -1.0; u[0][2] = 0.0;
	      v[0][0] = -1.0 ; v[0][1] = 1.0; v[0][2] = 0.0;
	      u[1][0] = -1.0 ; u[1][1] =  1.0; u[1][2] = 0.0;
	      v[1][0] = v[1][1] =  1.0; v[1][2] = 0.0;
	      u[2][0] = u[2][1] =  1.0; u[2][2] = 0.0;
	      v[2][0] = 1.0; v[2][1] = -1.0; v[2][2] = 0.0;
	      u[3][0] = -1.0 ; u[3][1] =  1.0; u[3][2] = 0.0;
	      v[3][0] = v[3][1] =  -1.0; v[3][2] = 0.0; 
	    }
	  
	  mPoint p1,p2,p3;
	  for(int j = 0;j<nbtris;j++)
	    {
	      splitTri tt(nbspl, mPoint (u[j][0],v[j][0],0.0),
			  mPoint (u[j][1],v[j][1],0.0),
			  mPoint (u[j][2],v[j][2],0.0));
	      
	      
	      while(tt.nextTri(p1,p2,p3))
		{
		  
		  
		  cell->interpolate(p1(0),p1(1),p1(2),val1);
		  cell->interpolate(p2(0),p2(1),p2(2),val2);
		  cell->interpolate(p3(0),p3(1),p3(2),val3);
		  
		  shape = view->getShape(s);
		
		  s++; // s is a counter over all plotting triangles 
		  for (int q=0;q<NbFieldsOfInterest;q++)
		    {
			  int k;
			  getLaw()->getNameAndSize(q, k, field);
			  if(!strcmp(field,"Intensity_RefCond"))
			  {				 
		      theLaw->getIthFieldOfInterest(q, dummy, val1, intr1, TACT);
		      theLaw->getIthFieldOfInterest(q, dummy, val2, intr2, TACT);
		      theLaw->getIthFieldOfInterest(q, dummy, val3, intr3, TACT);			  
			 
			  shape->addValue(0,intr1[0],intr2[0],intr3[0]);
			  }
		    }
		  
		}
	    }
	}
    }
/* 
	sprintf(name,"%s-sample%d","Intensity_RefCond",iteration);
	strcat(name,"Intensity_RefCond");
	strcat(name,".pos");
	ofstream out(name);
	out << "View \"Exported field\" {" << endl; 
	
	for(int j=0;j<view->getCurrentShape();j++) view->AllShapes[j]->print(out,0);
	out<< "};" << endl;
	out.close();  
	*/
}


void DGAnalysis::exportSummation(char *fileName, int iteration, int COUNTER)
{
char name[256]; 
sprintf(name,"%s-sample%d",fileName,iteration);
	strcat(name,"Intensity_ave");
	strcat(name,".pos");
	ofstream out(name);
	out << "View \"Exported field\" {" << endl; 
	
	for(int j=0;j<view->getCurrentShape();j++) view->AllShapes[j]->print(out,0, COUNTER);
	out<< "};" << endl;
	out.close();  
}
void DGAnalysis::exportGmsh_Volume(char* fileName, int iteration)
{
  int dim = n;
  int s = 0;
  int NN = 0;
  int nbtets;
  char field[256],name[256];
 
  //counting tets
  mMesh::iter it;
  mMesh::iter mesh_begin=theMesh->begin(dim);
  mMesh::iter mesh_end=theMesh->end(dim);
  for( it = mesh_begin;it != mesh_end;++it)  
    {
      mEntity *m = (*it);
      DGCell *cell = (DGCell*)m->getCell();;
      int fOrder = cell->getFunctionSpace()->order();
      if (fOrder<=1) nbtets=1; else nbtets=8; 
      NN +=nbtets;
    }

  //Are there any vector fields?
  int NbVec=0;
  int NbFieldsOfInterest = getLaw()->getNbFieldsOfInterest();
  for (int q=0;q<NbFieldsOfInterest;q++)
    {
      int k;
      getLaw()->getNameAndSize(q, k, field);
      if (k==3) NbVec++;
    }
  
  PostProGmshView view(NN);   // will contain the scalar data from all cells for printing
  PostProGmshView viewVector(NN); // will contain the vector data from all cells for printing
  PostProGmshShape *shape,*shapeVector;
  
  int nbpts,k,i;
  mPoint x1,x2,x3,x4; // vertices of a tet
  double val[10][MaxNbEqn];// solution values at vertices
  double intr1[3],intr2[3],intr3[3],intr4[3]; // field of interest , i.e. pressure, at vertices
  
  double pts [10][3] = { {0,0,0}, {0,0,1}, {0,1,0}, {1,0,0}, {0,0,0.5}, {0,0.5,0}, {0.5,0,0}, {0,0.5,0.5} , {0.5,0,0.5}, {0.5,0.5,0} };
  int tet [1][4] = { {0,1,2,3}};
  int tet8 [8][4] = { {0,4,5,6}, {1,4,7,8}, {2,5,7,9}, {3,6,8,9}, {7,4,5,6}, {7,5,9,6}, {7,9,8,6}, {7,8,4,6} };
  double ptsXYZ[10][3];
  
  for(it = mesh_begin;it != mesh_end;++it)
    {
      mEntity *m = (*it);
      DGCell *cell;
      cell = (DGCell*)m->getCell();
      int fOrder = cell->getFunctionSpace()->order();
      if (fOrder<2) nbtets=1; else nbtets=8; 
      if (nbtets == 1) nbpts =4; else nbpts=10;	
     
      for (k=0;k<nbpts;k++)
	{
	  cell->getMapping()->eval(pts[k][0],pts[k][1],pts[k][2],ptsXYZ[k][0],ptsXYZ[k][1],ptsXYZ[k][2]);
	  cell->interpolate(pts[k][0],pts[k][1],pts[k][2],&val[k][0]);
	}
      
      for(int j = 0;j<nbtets;j++)
	{
	  shape = view.NextInViewIsAScalarTetrahedron(NbFieldsOfInterest);
	  if (NbVec)
	    shapeVector = viewVector.NextInViewIsAVectorTetrahedron(NbVec);
	  s++; // s is a counter over all tets 
	  for (int q=0;q<NbFieldsOfInterest;q++)
	    {
	      int k;
	      getLaw()->getNameAndSize(q, k, field);
	      
	      for (i=0; i<3;i++)
		{
		  if (nbtets==1) x1(i) =  ptsXYZ[tet[j][0]][i]; else x1(i) =  ptsXYZ[tet8[j][0]][i];
		  if (nbtets==1) x2(i) =  ptsXYZ[tet[j][1]][i]; else x2(i) =  ptsXYZ[tet8[j][1]][i];
		  if (nbtets==1) x3(i) =  ptsXYZ[tet[j][2]][i]; else x3(i) =  ptsXYZ[tet8[j][2]][i];
		  if (nbtets==1) x4(i) =  ptsXYZ[tet[j][3]][i]; else x4(i) =  ptsXYZ[tet8[j][3]][i];
		} 
	      if (!strcmp(field,"limit")) 
		{
		  double vall[MaxNbEqn];
		  vall[0]=cell->limitStatus();
		  theLaw->getIthFieldOfInterest(q, x1, vall, intr1, TACT);
		  theLaw->getIthFieldOfInterest(q, x2, vall, intr2, TACT);
		  theLaw->getIthFieldOfInterest(q, x3, vall, intr3, TACT);
		  theLaw->getIthFieldOfInterest(q, x4, vall, intr4, TACT);
		}
	      else if (nbtets==1) 
		{
		  theLaw->getIthFieldOfInterest(q, x1, &val[tet[j][0]][0], intr1, TACT);
		  theLaw->getIthFieldOfInterest(q, x2, &val[tet[j][1]][0], intr2, TACT);
		  theLaw->getIthFieldOfInterest(q, x3, &val[tet[j][2]][0], intr3, TACT);
		  theLaw->getIthFieldOfInterest(q, x4, &val[tet[j][3]][0], intr4, TACT);
		}
	      else 
		{
		  theLaw->getIthFieldOfInterest(q, x1, &val[tet8[j][0]][0], intr1, TACT);
		  theLaw->getIthFieldOfInterest(q, x2, &val[tet8[j][1]][0], intr2, TACT);
		  theLaw->getIthFieldOfInterest(q, x3, &val[tet8[j][2]][0], intr3, TACT);
		  theLaw->getIthFieldOfInterest(q, x4, &val[tet8[j][3]][0], intr4, TACT);
		}
	      
	      if (k==3)
		{
		  shapeVector->setValue(0,intr1[0],intr1[1],intr1[2], 
					intr2[0],intr2[1],intr2[2], 
					intr3[0],intr3[1],intr3[2],
					intr4[0],intr4[1],intr4[2]);
		  shapeVector->setXYZ(x1(0),x1(1),x1(2),x2(0),x2(1),x2(2),x3(0),x3(1),x3(2),x4(0),x4(1),x4(2));
		}
	      else
		{
		  shape->setValue(q,intr1[0],intr2[0],intr3[0],intr4[0]);
		  shape->setXYZ(x1(0),x1(1),x1(2),x2(0),x2(1),x2(2),x3(0),x3(1),x3(2),x4(0),x4(1),x4(2));
		}
	    }
	}
    }
  
  for(int i=0;i<NbFieldsOfInterest;i++)
    {
      int k;
      getLaw()->getNameAndSize(i, k, field);
      {
	sprintf(name,"%s-sample%d",fileName,iteration);
	strcat(name,field);
	strcat(name,".pos");
	ofstream out(name);
	out << "View \"Exported field\" {" << endl; 
	
	if (k==3) for(int j=0;j<viewVector.getCurrentShape();j++) viewVector.AllShapes[j]->print(out,0);
	else for(int j=0;j<view.getCurrentShape();j++) view.AllShapes[j]->print(out,i);
	out<< "};" << endl;
	out.close();
      }
    }
}

