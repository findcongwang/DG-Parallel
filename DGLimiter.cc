#include "mEntity.h"
#include "Mapping.h"
#include "Integrator.h"
#include "ConservationLaw.h"
#include "FunctionSpace.h"
#include "DGLimiter.h"
#include "DGCell.h"
#include "FieldEvaluator.h"
#include "Constants.h"
#include "gEntity.h"
#include <stdio.h>
#include <math.h>
#include <list>

using namespace std;

static void recurGetAllsubs(mEntity *ent, list<mEntity*> &lis)
{
  int n = ent->getLevel();
  if(!ent->isAdjacencyCreated(n))lis.push_back(ent);
  else
    {
      for(int i=0;i<ent->size(n);i++)recurGetAllsubs(ent->get(n,i),lis);
    }
}


/******************************************************************************/

void DGSuperBee::limit(DGCell *vcell, double time)
{
  double Max[MaxNbEqn],Min[MaxNbEqn];
  int n = vcell->theMeshEntity->getLevel();
  int kk = 0;
  int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  
  double a = 1./3.464101615137754;
  double I[2][2],II[2][2],III[2][2];
  double u,v,w,weight,val[MaxNbEqn];
  double uvol,vvol,wvol,x,y,z;
  double mean[4][5];
  
  I[0][0] = 0.5; I[0][1] = 0.5; I[1][0] = -a; I[1][1] = a;
  II[0][0] = 0.0; II[0][1] = -0.5; II[1][0] = 2*a; II[1][1] = a;
  III[0][0] = -0.5; III[0][1] = 0.0; III[1][0] = -a; III[1][1] = -2*a;
  
  for(int i=0;i<vcell->theMeshEntity->size(n-1);i++)
    {
      mEntity *bound = vcell->theMeshEntity->get(n-1,i);
      if(bound->size(n) == 2)
	{
	  DGCell *other;
	  if(bound->get(n,0) != vcell->theMeshEntity)
	    other = (DGCell*)bound->get(n,0)->getCell();
	  else
	    other = (DGCell*)bound->get(n,1)->getCell();
	  if(!kk)
	    for(int j=0;j<other->cSize;j++)
	      {
		Max[j] = Min[j] = other->theMean[j];
		mean[1][j] = other->theMean[j]; 
	      }
	  else
	    {
	      for(int j=0;j<other->cSize;j++)
		{
		  Max[j]  = (other->theMean[j] > Max[j])?other->theMean[j]:Max[j];
		  Min[j]  = (other->theMean[j] < Min[j])?other->theMean[j]:Min[j];
		  mean[kk+1][j] = other->theMean[j]; 
		}
	      // if ( Min[0]<0 || Min[3] <0 ) printf("NEGATIVE MEAN  \n");
	    }
	  kk++;
	}
    }
  
  list<mEntity *> allSubs;
  for(int k=0;k<vcell->theMeshEntity->size(n-1);k++)
  recurGetAllsubs(vcell->theMeshEntity->get(n-1,k),allSubs);
    for(list<mEntity*>::const_iterator it = allSubs.begin();it!=allSubs.end();++it)
    {
      //MODIFY!!!
      DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
            DGCell *right = 0;
      DGCell *left = (DGCell*)cell->mleft->getCell();
      if(cell->mright)right = (DGCell*)cell->mright->getCell();
      
      // printf(" left %e \n",left->theMean[0]);
      //if(cell->mright) printf(" right %e \n",right->theMean[0]);
      int order = cell->computeOrder();
      order =1;
      GaussIntegrator gauss;
      int Nb = gauss.nbIntegrationPoints(cell->theBoundaryEntity,order);
      
      for(int i=0;i<Nb;i++)
	{
	    gauss.iPoint(cell->theBoundaryEntity,i,order,u,v,w,weight);
	    cell->theMapping->eval(u,v,w,x,y,z);
	    vcell->theMapping->invert(x,y,z,uvol,vvol,wvol);
	    vcell->interpolate(uvol,vvol,wvol,val);
	    //  printf("u,v,w %f % f %f\n", u,v,w);
	    // printf("x,y,z %f % f %f\n",x,y,z );
	    // printf("uvol,vvol,wvol %f % f %f\n", uvol,vvol,wvol); 
	    // if values are not physical
	    
	    for(int j=0;j<vcell->cSize;j++)
	      {
		if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) && fabs( vcell->theMean[j]) > 1.e-13 )	    
		  {
		    
		    //    printf("MAX %f  MIn %f Mean %f val %f \n",Max[j], Min[j], vcell->theMean[j], val[j]);
		    
		    if(val[j] > Max[j] || val[j] < Min[j] )
		      {
			if (Max[j] > vcell->theMean[j] && Min[j] > vcell->theMean[j])
			  { 
			    for (int k=1;k<fSize;k++)
			      vcell->theFieldsCoefficients->get(j,k) = 0.0; return;
			  }
			if (Max[j]<vcell->theMean[j] && Min[j]<vcell->theMean[j]) 
			  {
			    for (int k=1;k<fSize;k++)
			      vcell->theFieldsCoefficients->get(j,k) = 0.0; return;
			  }
			
			double b[2];
			double diff[4];
			diff[1] = fabs(mean[1][j]- vcell->theMean[j]);
			diff[2] = fabs(mean[2][j]- vcell->theMean[j]);
			diff[3] = fabs(mean[3][j]- vcell->theMean[j]);
			double pos =0.0;
			double neg = 0.0;
			if (diff[1]>0.0) pos += diff[1];else neg += -diff[1];
			if (diff[2]>0.0) pos += diff[2];else neg += -diff[2];
			if (diff[3]>0.0) pos += diff[3];else neg += -diff[3];
	    		//printf(" means %e %e %e %e\n",mean[1][j],mean[2][j],mean[3][j],vcell->theMean[j]);
			//printf("uvol,vvol,wvol %f % f %f\n", uvol,vvol,wvol); 
			//printf("theMean= %e val=%e\n",vcell->theMean[j],val[j]);
			//printf("Original coeff %e %e %e \n",vcell->theFieldsCoefficients[0][0],vcell->theFieldsCoefficients[0][1],vcell->theFieldsCoefficients[0][2]);
			double phiplus, phiminus;
			phiplus = (neg>pos) ? 1.0 : neg/pos; 
			phiminus = (pos>neg) ? 1.0 : pos/neg; 

			if (diff[1]>0.0) diff[1] *= phiplus;
			else diff[1]*=phiminus;
			if (diff[2]>0.0) diff[2] *= phiplus;
			else diff[2]*=phiminus;
			if (diff[3]>0.0) diff[3] *= phiplus;
			else diff[3]*=phiminus;
			
			b[0]=diff[1];
			b[1]= diff[2];
			//printf("b %e %e \n",b[0],b[1]);
			//printf("CASE 1\n");
			vcell->theFieldsCoefficients->get(j,1) = I[0][0]*b[0]+I[0][1]*b[1];
			vcell->theFieldsCoefficients->get(j,2)=I[1][0]*b[0]+I[1][1]*b[1];
			vcell->interpolate(0.0,0.5,0.0,val);
			if (val[j]<Min[j] || val[j]>Max[j])
			  {
			    vcell->theFieldsCoefficients->get(j,1)=0.0;
			    vcell->theFieldsCoefficients->get(j,2)=0.0;
			    printf("set coef. to zero at 1!!!!!!!!\n");
			  }
			// printf("coeff %e %e %e \n",vcell->theFieldsCoefficients[0][0],vcell->theFieldsCoefficients[0][1],vcell->theFieldsCoefficients[0][2]);
			
			//	printf("j= %d\n",j);
			/* vcell->interpolate(0.5,0.0,-0.5,val);
		     printf("Updated val=%e",val[0]);
		     vcell->interpolate(0.5,0.5,-0.5,val);
		     printf("   val=%e",val[0]);
		     vcell->interpolate(0.0,0.5,-0.5,val);
		     printf("   val=%e",val[0]);
		     printf("\n \n");*/
		      }
		  }
	      }
	} 
    }  
}

/******************************************************************************/

void DGBarthLimiter::limit(DGCell *vcell, double time)
{
  int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  int cSize = vcell->cSize;
  int n = vcell->theMeshEntity->getLevel();
  vcell->limit=1;
  int i,kk=0;
  double phin,val[MaxNbEqn], Max[MaxNbEqn],Min[MaxNbEqn];
//  double u,v,w,uvol,vvol,wvol,x,y,z,weight;
  
  list<mEntity*>::const_iterator it,allSubs_end;
  double dx = vcell->cellSize;
  if (fabs(vcell->fullJump/(pow(dx,2*(fOrder+1))*vcell->getPerimeter()))>3. )
    vcell->limit = 1;
 
  if (vcell->limit) {     
    int nbFaces = vcell->theMeshEntity->size(n-1);
    for(i=0;i<nbFaces;i++) {
      mEntity *bound = vcell->theMeshEntity->get(n-1,i);
      if(bound->size(n) == 2)
	{
	  DGCell *other;
	  if(bound->get(n,0) != vcell->theMeshEntity)
	    other = (DGCell*)bound->get(n,0)->getCell();
	  else
	    other = (DGCell*)bound->get(n,1)->getCell();
	  if(!kk)
	    for(int j=0;j<cSize;j++)
	      Max[j] = Min[j] = other->theMean[j];
	  else
	    {
	      for(int j=0;j<cSize;j++)
		{
		  Max[j]  = (other->theMean[j] > Max[j])?other->theMean[j]:Max[j];
		  Min[j]  = (other->theMean[j] < Min[j])?other->theMean[j]:Min[j];
		}
	    }
	  kk++;
	}
    }    
    
    for(int j=0;j<vcell->cSize;j++) vcell->limSlope[j] = 1.0;
    
    allSubs_end=vcell->allSubs.end();
    boundaryGaussPoint *pg;
    for(it = vcell->allSubs.begin();it!=allSubs_end;++it) {
      DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
      DGCell *right = 0;
      DGCell *left = (DGCell*)cell->mleft->getCell();
      if(cell->mright)right = (DGCell*)cell->mright->getCell();  
      
      int Nb = cell->getNbPtGauss();
      
      for(int i=0;i<Nb;i++)
	{
	  pg = cell->pt(i);
	  vcell->interpolate(pg->fctleft, val);
	  for(int j=0;j<cSize;j++)
	    {
	      //if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )
	      //if(pg->x>0.5 )
		  if(true)
 {
		
		if(val[j] > Max[j]) {
		  phin=(Max[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		  vcell->limit=2;
		}
		else if(val[j] < Min[j]) {
		  phin=(Min[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		  vcell->limit=2;
		}
		else {phin = 1.0;vcell->limit=3;}
	      }
	      else {phin = 1.0; vcell->limit = 4;
	      }
	      if(phin<0.0) {phin = 0.0;vcell->limit=5;}
	      vcell->limSlope[j] = (vcell->limSlope[j]<phin) ? 
		vcell->limSlope[j]:phin;
	    }
	}
    }
    //printf("alha 0 %e \n",vcell->limSlope[0]);
    while(1)
      { 
	int fSize1 = vcell->theFunctionSpace->size(1);
	int fSizei = vcell->theFunctionSpace->size(fOrder);
	int fSizej = vcell->theFunctionSpace->size(fOrder-1);
	
	for(int j=0;j<cSize;j++)
	  {	
	    if(fOrder == 1)
	      for (int k=1;k<fSize1;k++)
		vcell->theFieldsCoefficients->get(j,k) *= vcell->limSlope[j];
	    else	
	      for (int k=fSizej;k<fSizei;k++)
		if(vcell->limSlope[j] < 0.9)vcell->theFieldsCoefficients->get(j,k)= 0.0;
	  }	 
	if(fOrder == 1) break;
	fOrder --;
      }
  }  
}  


//Limiting at the center point

/*
void DGBarthLimiter::limit(DGCell *vcell, double time)
{
  int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  int n = vcell->theMeshEntity->getLevel();
  vcell->limit=1;
  int i,k, kk=0;
  double u,v,w,weight,val[MaxNbEqn], Max[MaxNbEqn],Min[MaxNbEqn];
  double uvol,vvol,wvol,x,y,z,phin;
  
  list<mEntity *> allSubs;
  list<mEntity*>::const_iterator it,allSubs_end;
  double dx = vcell->getMaxH();
  if (fabs(vcell->fullJump/(pow(dx,2*(fOrder+1))*vcell->getPerimeter()))>1. )
    vcell->limit = 1;
  //printf("Jump %f, JUmp^3 %f, maxVal %f, dx %f\n",vcell->fullJump,pow(vcell->fullJump,3),vcell->maxVal,dx);
  
  if (vcell->limit)    
    {  
      int Size = vcell->theMeshEntity->size(n-1);
      for(i=0;i<Size;i++)
	{
	  mEntity *bound = vcell->theMeshEntity->get(n-1,i);
	  if(bound->size(n) == 2)
	    {
	      DGCell *other;
	      if(bound->get(n,0) != vcell->theMeshEntity)
		other = (DGCell*)bound->get(n,0)->getCell();
	      else
		other = (DGCell*)bound->get(n,1)->getCell();
	      if(!kk)
		for(int j=0;j<other->cSize;j++)
		  Max[j] = Min[j] = other->theMean[j];
	      else
		{
		  for(int j=0;j<other->cSize;j++)
		    {
		      Max[j]  = (other->theMean[j] > Max[j])?other->theMean[j]:Max[j];
		      Min[j]  = (other->theMean[j] < Min[j])?other->theMean[j]:Min[j];
		    }
		}
	      kk++;
	    }
	}    
      
      for(k=0;k<vcell->theMeshEntity->size(n-1);k++)
	recurGetAllsubs(vcell->theMeshEntity->get(n-1,k),allSubs);
      
      for(int j=0;j<vcell->cSize;j++) vcell->limSlope[j] = 1.0;
      allSubs_end=allSubs.end();
      
      for(it = allSubs.begin();it!=allSubs_end;++it)
	{if ((*it)->getClassification()->getId()!=20000) {
	DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
	int order =1;
	GaussIntegrator gauss;
	int Nb = gauss.nbIntegrationPoints(cell->theBoundaryEntity,order);
	for(int i=0;i<Nb;i++)
	    {
	      gauss.iPoint(cell->theBoundaryEntity,i,order,u,v,w,weight);
	      cell->theMapping->eval(u,v,w,x,y,z);
	      vcell->theMapping->invert(x,y,z,uvol,vvol,wvol);
	      if (x>0.5) {
		vcell->interpolate(uvol,vvol,wvol,val);
		for(int j=0;j<vcell->cSize;j++) {
		  //		  if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )
		  if(true)
		    {
		      if(val[j] > Max[j])
			{
			  //printf("j=%d x=%e y=%e val=%e max=%e\n",j,x,y,val[j],Max[j]);
			  phin=(Max[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
			  vcell->limit=2;
			}
		      else if(val[j] < Min[j])
			{
			  //printf("j=%d x=%e y=%e val=%e min = %e\n",j,x,y,val[j],Min[j]);
			  phin=(Min[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
			  vcell->limit=2;
			}
		      else {phin = 1.0;vcell->limit=3;}
		    }
		  else phin = 1.0;
		  if(phin<0.0) {phin = 0.0;vcell->limit=5;}
		  vcell->limSlope[j] = (vcell->limSlope[j]<phin) ? 
		    vcell->limSlope[j]:phin;
		}
	      }
	    }
	}
	}
      //printf("alha 0 %e \n",vcell->limSlope[0]);
      while(1)
	{ 
	  int fSize1 = vcell->theFunctionSpace->size(1);
	  int fSizei = vcell->theFunctionSpace->size(fOrder);
	  int fSizej = vcell->theFunctionSpace->size(fOrder-1);
	  
	  for(int j=0;j<vcell->cSize;j++)
	    {	
	      if(fOrder == 1)
		for (int k=1;k<fSize1;k++)
		  vcell->theFieldsCoefficients->get(j,k) *= vcell->limSlope[j];
	      else	
		for (int k=fSizej;k<fSizei;k++)
		  if(vcell->limSlope[j] < 0.9)vcell->theFieldsCoefficients->get(j,k)= 0.0;
	    }	 
	  if(fOrder == 1) break;
	  fOrder --;
	}
    }  
}
*/

/******************************************************************************/

void DGBarthLimiterEuler::limit(DGCell *vcell, double time)
{
  int fSize = vcell->fSize;
  int fOrder = vcell->fOrder;
  int cSize  = vcell->cSize;
  int n = vcell->theMeshEntity->getLevel();
  vcell->limit=5;
  int j,k,kk=0;
  double rhom,rhoM,um,uM,vm,vM,pm,pM;
//  double bExtrema[MaxNbEqn],Max[MaxNbEqn],Min[MaxNbEqn];
  double val[MaxNbEqn],alfa[MaxNbEqn];
  double const GAMMA=1.4;
  double const GAMMA_1=0.4;
list<mEntity *> allSubs;
  list<mEntity*>::const_iterator it;
    if (vcell->limit)    
    {  
      
      for(int i=0;i<vcell->theMeshEntity->size(n-1);i++)
	{
	  mEntity *bound = vcell->theMeshEntity->get(n-1,i);
	  if(bound->size(n) == 2)
	    {
	      DGCell *other;
	      if(bound->get(n,0) != vcell->theMeshEntity)
		other = (DGCell*)bound->get(n,0)->getCell();
	      else
		other = (DGCell*)bound->get(n,1)->getCell();
	      if(!kk)
		{
		  rhoM = rhom = other->theMean[0];
		  uM   = um   = other->theMean[1]/rhoM;
		  vM   = vm   = other->theMean[2]/rhoM;
		  pM   = pm   = (GAMMA_1)*(other->theMean[3] - (uM*uM + vM*vM) *0.5*rhoM);
		}
	      else
		{
		  double rhoO = other->theMean[0]; 
		  double invrhoO = 1./other->theMean[0];
		  double uO =other->theMean[1]*invrhoO;  
		  double vO =other->theMean[2]*invrhoO;
		  double pO   = (GAMMA_1)*(other->theMean[3] - (uO*uO + vO*vO) *0.5*rhoO);
		  
		  rhoM = ( rhoO > rhoM ) ? rhoO : rhoM ; 
		  rhom = ( rhoO < rhom ) ? rhoO : rhom ; 
		  uM   = ( uO   >  uM  ) ?   uO : uM ; 
		  um   = ( uO   <  um  ) ?   uO : um ; 
		  vM   = ( vO   >  vM  ) ?   vO : vM ; 		  
		  vm   = ( vO   <  vm  ) ?   vO : vm ; 		  
		  pM   = ( pO   >  pM  ) ?   pO : pM ; 		  
		  pm   = ( pO   <  pm  ) ?   pO : pm ; 		  
		}
	      kk++;
	    }
	}
       
      for(j=0;j<cSize;j++)  vcell->limSlope[j]= alfa[j] =1.0;
      
      for(k=0;k<vcell->theMeshEntity->size(n-1);k++)
	recurGetAllsubs(vcell->theMeshEntity->get(n-1,k),allSubs);      
     
      for(it = vcell->allSubs.begin();it!=vcell->allSubs.end();++it)
	{
	  DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
	  
	  for(int i=0;i<cell->nbPtGauss;i++)
	    {
	      boundaryGaussPoint *pg = cell->pt(i);
	      if(vcell==cell->left)vcell->interpolate(pg->fctleft,val);
	      else vcell->interpolate(pg->fctrght,val);
	      double u,v,p;
	      u = val[1]/val[0]; v=val[2]/val[0]; p = GAMMA_1*(val[3]-0.5*(u*u +v*v)*val[0]);
	      //if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )     
		{
		  if(val[0]>rhoM) alfa[0]=(rhoM-vcell->theMean[0])/(val[0]-vcell->theMean[0]);
		  if(val[0] < rhom) alfa[0] = (rhom-vcell->theMean[0])/(val[0]-vcell->theMean[0]);

		  if(u>uM) alfa[1]=(uM*val[0]-vcell->theMean[1])/(val[1]-vcell->theMean[1]);
		  if (u<um) alfa[1]=(um*val[0]-vcell->theMean[1])/(val[1]-vcell->theMean[1]);
		  
		  if(v>vM) alfa[2]=(vM*val[0]-vcell->theMean[2])/(val[2]-vcell->theMean[2]);
		  if(v<vm) alfa[2]=(vm*val[0]-vcell->theMean[2])/(val[2]-vcell->theMean[2]);	 

		  if(p>pM) alfa[3]=(pM/(GAMMA-1)+0.5*(u*u+v*v)*val[0]-vcell->theMean[3])/(val[3]-vcell->theMean[3]);  
		  if(p<pm) alfa[3]=(pm/(GAMMA-1)+0.5*(u*u+v*v)*val[0]-vcell->theMean[3])/(val[3]-vcell->theMean[3]);  
		  /*
		  {
		      phin = (Max[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      if (j==0) vcell->limit=0;
		      if (j==1) vcell->limit=1;
		      if (j==2) vcell->limit=2;
		      if (j==3) vcell->limit=3;
		      //printf("limiting %d\n",j);
		    }
		  else if(val[j] < Min[j])
		    {
		      phin = (Min[j]-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      if (j==0) vcell->limit=0;
		      if (j==1) vcell->limit=1;
		      if (j==2) vcell->limit=2;
		      if (j==3) vcell->limit=3;
		      //vcell->limit=3;
		      //printf("limiting %d\n",j);
		      //printf("val<Min val[%d]=%e Min=%e Mean=%e\n",j,val[j],Min[j],vcell->theMean[j]);
		    }
		    */
		}
	      for(int j=0;j<cSize;j++) 
		{
		  if (alfa[j]<0.0) alfa[j] = 0.0;
		  vcell->limSlope[j] = (vcell->limSlope[j]<alfa[j])?vcell->limSlope[j]:alfa[j];
		}
	    }
	}
	 // printf("alha 0 %e \n",vcell->limSlope[0]);
     // for(int j=0;j<cSize;j++) printf("vcell->limSlope[j] %e ",vcell->limSlope[j]);
	  //printf("\n");
      while(1)
	{ 
	  int fSize1 = vcell->theFunctionSpace->size(1);
	  int fSizei = vcell->theFunctionSpace->size(fOrder);
	  int fSizej = vcell->theFunctionSpace->size(fOrder-1);
	  
	  for(int j=0;j<cSize;j++)
	    {	
	      if(fOrder == 1)
		for (int k=1;k<fSize1;k++)
		  vcell->theFieldsCoefficients ->get(j,k)*= vcell->limSlope[j];
	      else	
		for (int k=fSizej;k<fSizei;k++)
		  if(vcell->limSlope[j] < 0.9)vcell->theFieldsCoefficients->get(j,k)= 0.0;
	    }	 
	  if(fOrder == 1) break;
	  fOrder --;
	}
    }
  
}  
/******************************************************************************/

void DGVertexLimiter::limit(DGCell *vcell, double time)
{
  int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  int n = vcell->theMeshEntity->getLevel();
  vcell->limit=5;
  int kk=0;
//s  double u,v,w,weight,val[MaxNbEqn], Max[MaxNbEqn],Min[MaxNbEqn];
//  double uvol,vvol,wvol,x,y,z;
  double phin;
  list<mEntity *> allSubs;
  list<mEntity*>::const_iterator it;
  //double dx = vcell->getMaxH();
  //  if (fabs(pow(vcell->fullJump,3)/vcell->maxVal/sqrt(pow(dx,3*(fOrder+1)))/dx)>1. )
  //vcell->limit = 1;
  //printf("Jump %f, JUmp^3 %f, maxVal %f, dx %f\n",vcell->fullJump,pow(vcell->fullJump,3),vcell->maxVal,dx);
  
if (vcell->limit)    
  {  
    double Max[MaxNbEqn],Min[MaxNbEqn];
    int n = vcell->theMeshEntity->getLevel();
    int k;
    int kk = 0;
    int fSize = vcell->fSize;
    int fOrder = vcell->theFunctionSpace->order();
    
    double val[MaxNbEqn];
    list<mEntity*>::const_iterator it;
    
    list<mEntity *> allVert;
    mEntity *ent;
    for(k=0;k<vcell->theMeshEntity->size(0);k++)
      {
	ent = vcell->theMeshEntity->get(0,k); 
	for(int i=0;i<ent->size(2);i++)
	  {
	    allVert.push_back(ent->get(n,i));
	    //ent->get(n,i)->print();
	  }
      }
    
    for(it = allVert.begin();it!=allVert.end();++it)
      {
	DGCell *cell = (DGCell*)(*it)->getCell();
	
	// if(cell->theMeshEntity->size() == 2)
	if(!kk)
	  {
	    for(int j=0;j<cell->cSize;j++)
	      Max[j] = Min[j] = cell->theMean[j];
	  }
	else
	  {
	    for(int j=0;j<cell->cSize;j++)
	      {
		Max[j]  = (cell->theMean[j] > Max[j])?cell->theMean[j]:Max[j];
		Min[j]  = (cell->theMean[j] < Min[j])?cell->theMean[j]:Min[j];
		// printf("Min= %e Max=%e Mean =%e\n",Min[j],Max[j],cell->theMean[j]);
		//(*it)->print();
	      }
	    //if ( Min[0]<0 || Min[3] <0 ) printf("NEGATIVE MEAN  \n");
	  }
	kk++;
      }
  //printf("Min= %e Max=%e\n",Min[0],Max[0]);
  
  for(int j=0;j<vcell->cSize;j++)
    vcell->limSlope[j] = 1.0;
  
  for(k=0;k<vcell->theMeshEntity->size(n-1);k++)
    recurGetAllsubs(vcell->theMeshEntity->get(n-1,k),allSubs);
  for(it = allSubs.begin();it!=allSubs.end();++it)
    {
      DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
      
      DGCell *right = 0;
      DGCell *left = (DGCell*)cell->mleft->getCell();
      if(cell->mright)right = (DGCell*)cell->mright->getCell();
      
      double V[3][3];
      V[0][0] = 0.0; V[0][1] = 0.0; V[0][2] = 0.0;
      V[1][0] = 1.0; V[1][1] = 0.0; V[1][2] = 0.0;
      V[2][0] = 0.0; V[2][1] = 1.0; V[2][2] = 0.0;
      for(int i=0;i<3;i++)
	for(int j=0;j<vcell->cSize;j++)
	  {
	    vcell->interpolate(V[i][0],V[i][1],V[i][2],val);
	    double x;
	    if(true)
	      //	    if(fabs(val[j] - vcell->theMean[j]) > 1.e-8 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )	    
	      {
		if(val[j] > Max[j])
		  {
		    x = Max[j];
		    phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		    vcell->limit=3;
		  }
		else if(val[j] < Min[j])
		  {
		    x = Min[j];
		    phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		    vcell->limit=3;
		  }
		else  phin = 1.0;
	      }
	    else       phin = 1.0;
	    if(phin<0.0) {phin = 0.0;vcell->limit=2;}
	    vcell->limSlope[j] = (vcell->limSlope[j]<phin)?vcell->limSlope[j]:phin;
	  }
    }   

  while(1)
    { 
      int fSize1 = vcell->theFunctionSpace->size(1);
      int fSizei = vcell->theFunctionSpace->size(fOrder);
      int fSizej = vcell->theFunctionSpace->size(fOrder-1);
      
      for(int j=0;j<vcell->cSize;j++)
	{	
	  if(fOrder == 1)
	    for (int k=1;k<fSize1;k++)
	      vcell->theFieldsCoefficients ->get(j,k)*= vcell->limSlope[j];
	  else	
	    for (int k=fSizej;k<fSizei;k++)
	      if(vcell->limSlope[j] < 0.9)vcell->theFieldsCoefficients ->get(j,k)= 0.0;
	}	 
      if(fOrder == 1) break;
      fOrder --;
    } 
  /*
///MODIFY!!!
    int order = cell->computeOrder();
    order=1;
      GaussIntegrator gauss(cell->theBoundaryEntity);
      for(int i=0;i<gauss.nbIntegrationPoints(order);i++)
	{
	  gauss.iPoint(i,order,u,v,w,weight);
	  cell->theMapping->eval(u,v,w,x,y,z);
	  vcell->theMapping->invert(x,y,z,uvol,vvol,wvol);
	  vcell->interpolate(uvol,vvol,wvol,val);
	  // printf("u=%f v=%f, x=%f y=%f,uvol=%f vvol=%f\n",u,v,x,y,uvol,vvol);
	  for(int j=0;j<vcell->cSize;j++)
	    {
	      double x;
	      if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )	    
		{
		  if(val[j] > Max[j])
		    {
		      x = Max[j];
		      phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      vcell->limit=3;
		      //    printf("val>Max val[%d]=%e Max=%e Mean=%e\n",j,val[j],Max[j],vcell->theMean[j]);
		    }
		  else if(val[j] < Min[j])
		    {
		      x = Min[j];
		      phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      vcell->limit=3;
		      //printf("val<Min val[%d]=%e Min=%e Mean=%e\n",j,val[j],Min[j],vcell->theMean[j]);
		    }
		  else 
		    {
		      phin = 1.0;
		    }
		}
	      else
		phin = 1.0;
	      if(phin<0.0) {phin = 0.0;vcell->limit=2;}
	      vcell->limSlope[j] = (vcell->limSlope[j]<phin)?vcell->limSlope[j]:phin;
	      }
	      }*/
  }
} 

/************************************************/
void VertLimiter::limit(DGCell *vcell, double time)
{
  const int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  const int cSize = vcell->cSize;
  const int n = vcell->theMeshEntity->getLevel();
  vcell->limit=1;
  int  kk=0;
  //double u,v,w,weight;
  //double uvol,vvol,wvol,x,y,z,phin;
  double phin,val[MaxNbEqn], Max[MaxNbEqn],Min[MaxNbEqn];
  list<mEntity*>::const_iterator it,itt;
  list<mEntity*>::const_iterator vertEnd=vcell->allVert.end(); 
  //double dx = vcell->getMaxH();
  //  if (fabs(pow(vcell->fullJump,3)/vcell->maxVal/sqrt(pow(dx,3*(fOrder+1)))/dx)>1. )
  //vcell->limit = 1;
  //printf("Jump %f, JUmp^3 %f, maxVal %f, dx %f\n",vcell->fullJump,pow(vcell->fullJump,3),vcell->maxVal,dx);
  
  if (vcell->limit) {  
    for(it = vcell->allVert.begin();it!=vertEnd;++it)
      {
	DGCell *cell = (DGCell*)(*it)->getCell();
	
	if(!kk)
	  {
	    for(int j=0;j<cSize;j++)
	      Max[j] = Min[j] = cell->theMean[j];
	  }
	else
	    for(int j=0;j<cSize;j++)
	      {
		Max[j]  = (cell->theMean[j] > Max[j])?cell->theMean[j]:Max[j];
		Min[j]  = (cell->theMean[j] < Min[j])?cell->theMean[j]:Min[j];
		//printf("Min= %e Max=%e Mean =%e\n",Min[j],Max[j],cell->theMean[j]);
		//(*it)->print();
	      }
	//if ( Min[0]<0 || Min[3] <0 ) printf("NEGATIVE MEAN  \n");
	
	kk++;
      }
    
    for(int j=0;j<cSize;j++)
      vcell->limSlope[j] = 1.0;
         boundaryGaussPoint *pg;
    list<mEntity*>::const_iterator endE = vcell->allSubs.end();
    for(it = vcell->allSubs.begin();it!=endE;++it)
      { if ((*it)->getClassification()->getId()!=20000) {
	DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();  
	//one-two poit switch starts here    
	/*int order=1;
	GaussIntegrator gauss;
	int Nb = gauss.nbIntegrationPoints( cell->theBoundaryEntity,order);
	for(int i=0;i<Nb;i++)
	  {
	    gauss.iPoint(cell->theBoundaryEntity,i,order,u,v,w,weight);
	    cell->theMapping->eval(u,v,w,x,y,z);
	    vcell->theMapping->invert(x,y,z,uvol,vvol,wvol);*/
	// and ends here
	int Nb = cell->getNbPtGauss();
	for(int i=0;i<Nb;i++)
	  {
	    pg = cell->pt(i);
	    vcell->interpolate(pg->fctleft, val);
	    //vcell->interpolate(uvol,vvol,wvol,val);
	    
	    for(int j=0;j<cSize;j++)
	      {
		double x;
		//if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )	
		if (1)//(pg->x>0.5)
		  {
		    if(val[j] > Max[j])
		      {
		      x = Max[j];
		      phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      vcell->limit=2;
		       if (j==0)  printf("val>Max val[%d]=%e Max=%e Mean=%e\n",j,val[j],Max[j],vcell->theMean[j]);
		      }
		  else if(val[j] < Min[j])
		    {
		      x = Min[j];
		      phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		      vcell->limit=2;
		     if (j==0) printf("val<Min val[%d]=%e Min=%e Mean=%e\n",j,val[j],Min[j],vcell->theMean[j]);
		    }
		    else 
		      {
			phin = 1.0;vcell->limit=3;
		      }
		  }
		else
		  phin = 1.0;
		if(phin<0.0) {phin = 0.0;vcell->limit=5;}
		vcell->limSlope[j] = (vcell->limSlope[j]<phin)?vcell->limSlope[j]:phin;
		if (vcell->limSlope[j]<1.) {
		  for(itt = vcell->allVert.begin();itt!=vertEnd;++itt){
		    DGCell *ccell = (DGCell*)(*itt)->getCell();
	if (j==0)	    printf(" %e ",ccell->theMean[j]);}
		 if (j==0) printf("val [%d] =%e MIN[j] =%e MAx[j] =%e vcell->limSlope[j]=%e \n",j,val[j],Min[j],Max[j],vcell->limSlope[j]);
		  if (j==0) printf(" %e %e\n",pg->x,pg->y);
		  printf("\n");
		  }
	      }
	  }
      }
      }
    while(1)
      { 
	int fSize1 = vcell->theFunctionSpace->size(1);
	int fSizei = vcell->theFunctionSpace->size(fOrder);
      int fSizej = vcell->theFunctionSpace->size(fOrder-1);
      
      for(int j=0;j<cSize;j++)
	{	
	  if(fOrder == 1)
	    for (int k=1;k<fSize1;k++)
	      vcell->theFieldsCoefficients ->get(j,k)*= vcell->limSlope[j];
	  else	
	    for (int k=fSizej;k<fSizei;k++)
	      if(vcell->limSlope[j] < 0.9)vcell->theFieldsCoefficients ->get(j,k)= 0.0;
	}	 
      if(fOrder == 1) break;
      fOrder --;
      } 
  }
} 


/**************************************************************/

void DGLimiter::computeMinMaxEdge(DGCell *vcell)
{
  double Max[MaxNbEqn],Min[MaxNbEqn];
  int n = vcell->theMeshEntity->getLevel();
  int k;
  int kk = 0;
  int fSize = vcell->fSize;
  int fOrder = vcell->theFunctionSpace->order();
  
  double u,v,w,weight,val[MaxNbEqn];
  double uvol,vvol,wvol,x,y,z,phin;
  
  

  for(int i=0;i<vcell->theMeshEntity->size(n-1);i++)
    {
      mEntity *bound = vcell->theMeshEntity->get(n-1,i);
      if(bound->size(n) == 2)
	{
	  DGCell *other;
	  if(bound->get(n,0) != vcell->theMeshEntity)
	    other = (DGCell*)bound->get(n,0)->getCell();
	  else
	    other = (DGCell*)bound->get(n,1)->getCell();
	  if(!kk)
	    {
	      for(int j=0;j<other->cSize;j++)
		Max[j] = Min[j] = other->theMean[j];
	    }
	  else
	    {
	      for(int j=0;j<other->cSize;j++)
		{
		  Max[j]  = (other->theMean[j] > Max[j])?other->theMean[j]:Max[j];
		  Min[j]  = (other->theMean[j] < Min[j])?other->theMean[j]:Min[j];
		  //  printf("Min= %e Max=%e\n",Min[j],Max[j]);
		}
	      // if ( Min[0]<0 || Min[3] <0 ) printf("NEGATIVE MEAN  \n");
	    }
	  kk++;
	}
    }
  
  //printf("Min= %e Max=%e\n",Min[0],Max[0]);

  list<mEntity *> allSubs;
  list<mEntity*>::const_iterator it;
  for(k=0;k<vcell->theMeshEntity->size(n-1);k++)
    recurGetAllsubs(vcell->theMeshEntity->get(n-1,k),allSubs);

  for(it = allSubs.begin();it!=allSubs.end();++it)
    {
      DGBoundaryCell *cell = (DGBoundaryCell*)(*it)->getCell();      
      ///DOn't need this!
      DGCell *right = 0;
      DGCell *left = (DGCell*)cell->mleft->getCell();
      if(cell->mright)right = (DGCell*)cell->mright->getCell();  
  for(int j=0;j<vcell->cSize;j++)
    vcell->limSlope[j] = 1.0;
  
  int order = cell->computeOrder();
  order=1;
  GaussIntegrator gauss;
  int Nb = gauss.nbIntegrationPoints(cell->theBoundaryEntity,order);
  for(int i=0;i<Nb;i++)
    {
      gauss.iPoint(cell->theBoundaryEntity,i,order,u,v,w,weight);
      cell->theMapping->eval(u,v,w,x,y,z);
      vcell->theMapping->invert(x,y,z,uvol,vvol,wvol);
      vcell->interpolate(uvol,vvol,wvol,val);
      // printf("u=%f v=%f, x=%f y=%f,uvol=%f vvol=%f\n",u,v,x,y,uvol,vvol);
      for(int j=0;j<vcell->cSize;j++)
	{
	  double x;
	  if(fabs(val[j] - vcell->theMean[j]) > 1.e-5 * fabs(val[j] + vcell->theMean[j]) &&fabs( vcell->theMean[j]) > 1.e-13 )	    
	    {
	      if(val[j] > Max[j])
		{
		  x = Max[j];
		  phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		  vcell->limit=3;
		  //    printf("val>Max val[%d]=%e Max=%e Mean=%e\n",j,val[j],Max[j],vcell->theMean[j]);
		}
	      else if(val[j] < Min[j])
		{
		  x = Min[j];
		  phin = (x-vcell->theMean[j])/(val[j]-vcell->theMean[j]);
		  vcell->limit=3;
		  //printf("val<Min val[%d]=%e Min=%e Mean=%e\n",j,val[j],Min[j],vcell->theMean[j]);
		}
	      else 
		{
		  phin = 1.0;
		}
	    }
	  else
	    phin = 1.0;
	  if(phin<0.0) {phin = 0.0;vcell->limit=2;}
	  vcell->limSlope[j] = (vcell->limSlope[j]<phin)?vcell->limSlope[j]:phin;
	}
    }
}
}

/*********************************************************************************/

double minmod(double,double,double,bool&);
double minmod(double,double,bool &);
#define absmin(a,b) (fabs(a)<fabs(b) ? a:b)
void DGMomentLimiter::limit(DGCell *vcell,double time)
{
  int d,k,i,q;
  int k_im1,km1_i,i_km1,im1_k;
  DGCell *left,*right,*top,*bottom;
  left=right=top=bottom=0;
  mPoint p1 = vcell->pMin;
  mPoint p2 = vcell->pMax;
  mPoint p=(p1+p2);
  p(0)/=2.;p(1)/=2.;p(2)/=2.;
  mPoint omin,omax;
  const int cSize = vcell->cSize;
  const int n = vcell->theMeshEntity->getLevel();
  const double small_number = 1.0e-10;
  for(i=0;i<vcell->theMeshEntity->size(n-1);i++)
    {
      mEntity *bound = vcell->theMeshEntity->get(n-1,i);
      DGCell *other=0;
      if(bound->get(n,0) != vcell->theMeshEntity)
	other = (DGCell*)bound->get(n,0)->getCell();
        else if (  bound->get(n,1)) other = (DGCell*)bound->get(n,1)->getCell();
      if (other)
	{
	  omin = other -> pMin;
	  omax = other -> pMax;
	  if (p(0) > omax(0)) left = other; 
	  if (p(0) < omin(0)) right = other; 
	  if (p(1) > omax(1)) bottom = other; 
	  if (p(1) < omin(1)) top = other; 
	}	 
    }
   
  double tmp1,tmp2;
  bool isLimited;
  bool stopLimiter = 0;
  
  vcell->theFieldsCoefficients->copy();

  for(q=0; q<cSize; q++) {
    d=(vcell->fOrder+1)*(vcell->fOrder+1)-1;
    vcell->limit = d;
    for (k=vcell->fOrder; k>0; k--)
      {
	if (stopLimiter) break;
	
	if (fabs(vcell->theFieldsCoefficients->getCopy(q,d))>small_number) 
	  {
	    if (top && bottom)
	      tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (top->theFieldsCoefficients->get(q,d-2)-vcell->theFieldsCoefficients->get(q,d-2))/sqrt(double (2*k+1)/(2*k-1)),
			    (vcell->theFieldsCoefficients->get(q,d-2)-bottom->theFieldsCoefficients->get(q,d-2))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
	    else if (top && !bottom) 
	      tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (top->theFieldsCoefficients->get(q,d-2)-vcell->theFieldsCoefficients->get(q,d-2))/sqrt(double (2*k+1)/(2*k-1)), isLimited);
	    else if (!top&&bottom)
	      tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (vcell->theFieldsCoefficients->get(q,d-2)-bottom->theFieldsCoefficients->get(q,d-2))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
	    
	    bool flag = isLimited;
	    
	    if (left && right)
	      tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (right->theFieldsCoefficients->get(q,d-1)-vcell->theFieldsCoefficients->get(q,d-1))/sqrt(double (2*k+1)/(2*k-1)),
			    (vcell->theFieldsCoefficients->get(q,d-1)-left->theFieldsCoefficients->get(q,d-1))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
	    else if (left && !right)
	      tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (vcell->theFieldsCoefficients->get(q,d-1)-left->theFieldsCoefficients->get(q,d-1))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
	    else if (!left && right)
	      tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
			    (right->theFieldsCoefficients->get(q,d-1)-vcell->theFieldsCoefficients->get(q,d-1))/sqrt(double (2*k+1)/(2*k-1)), isLimited);
	    
	    if (flag || isLimited) {
	      vcell->limit=d-1;
	      vcell->theFieldsCoefficients->getCopy(q,d) = absmin(tmp1,tmp2);
	      //vcell->getMeshEntity()->print();
	    } 
	    else {
	      stopLimiter=1;
	      break;
	    }
	  }
	//	else printf("A coeff smaller than small number\n");
	d--;
	
	for (i=k-1; i>=0; i--) 
	  {
	    if (fabs(vcell->theFieldsCoefficients->getCopy(q,d))>small_number)
	      {
		double tmp1,tmp2;
		if (i<k-1) i_km1 =(k-1)*(k-1)+2*i+1; else i_km1 = i*i+(k-1)*2; 
		if (i>0) im1_k = k*k+2*(i-1)+1;
		
		if (top && bottom)
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(top->theFieldsCoefficients->get(q,i_km1)-vcell->theFieldsCoefficients->get(q,i_km1))/sqrt(double (2*k+1)/(2*k-1)),
				(vcell->theFieldsCoefficients->get(q,i_km1)-bottom->theFieldsCoefficients->get(q,i_km1))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
		else if (top && !bottom) 
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(top->theFieldsCoefficients->get(q,i_km1)-vcell->theFieldsCoefficients->get(q,i_km1))/sqrt(double (2*k+1)/(2*k-1)), isLimited);
		else if (!top&&bottom)
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(vcell->theFieldsCoefficients->get(q,i_km1)-bottom->theFieldsCoefficients->get(q,i_km1))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
		
		bool flag = isLimited;
		
		if (left && right&&i)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(right->theFieldsCoefficients->get(q,im1_k)-vcell->theFieldsCoefficients->get(q,im1_k))/sqrt(double (2*i+1)/(2*i-1)),
				(vcell->theFieldsCoefficients->get(q,im1_k)-left->theFieldsCoefficients->get(q,im1_k))/sqrt(double (2*i+1)/(2*i-1)),isLimited);
		else if (left && !right&&i)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(vcell->theFieldsCoefficients->get(q,im1_k)-left->theFieldsCoefficients->get(q,im1_k))/sqrt(double (2*i+1)/(2*i-1)),isLimited);
		else if (!left && right&&i)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(right->theFieldsCoefficients->get(q,im1_k)-vcell->theFieldsCoefficients->get(q,im1_k))/sqrt(double (2*i+1)/(2*i-1)), isLimited);
		
		if (flag || isLimited) {
		  vcell->limit=d-1;
		  if (i) vcell->theFieldsCoefficients->getCopy(q,d) = absmin(tmp1,tmp2);
		  else vcell->theFieldsCoefficients->getCopy(q,d) = tmp2;
		} 
	      }
	    d--;
	    
	    if (fabs(vcell->theFieldsCoefficients->getCopy(q,d))>small_number)
	      {
		isLimited = 0;
		k_im1 = k*k+2*(i-1);
		km1_i = (k-1)*(k-1)+2*i;
		if (top && bottom &&i)
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(top->theFieldsCoefficients->get(q,k_im1)-vcell->theFieldsCoefficients->get(q,k_im1))/sqrt(double (2*i+1)/(2*i-1)),
				(vcell->theFieldsCoefficients->get(q,k_im1)-bottom->theFieldsCoefficients->get(q,k_im1))/sqrt(double (2*i+1)/(2*i-1)),isLimited);
		else if (top && !bottom && i) 
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(top->theFieldsCoefficients->get(q,k_im1)-vcell->theFieldsCoefficients->get(q,k_im1))/sqrt(double (2*i+1)/(2*i-1)), isLimited);
		else if (!top&&bottom && i)
		  tmp2 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(vcell->theFieldsCoefficients->get(q,k_im1)-bottom->theFieldsCoefficients->get(q,k_im1))/sqrt(double (2*i+1)/(2*i-1)),isLimited);
		
		bool flag = isLimited;
		
		if (left && right)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(right->theFieldsCoefficients->get(q,km1_i)-vcell->theFieldsCoefficients->get(q,km1_i))/sqrt(double (2*k+1)/(2*k-1)),
				(vcell->theFieldsCoefficients->get(q,km1_i)-left->theFieldsCoefficients->get(q,km1_i))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
		else if (left && !right)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(vcell->theFieldsCoefficients->get(q,km1_i)-left->theFieldsCoefficients->get(q,km1_i))/sqrt(double (2*k+1)/(2*k-1)),isLimited);
		else if (!left && right)
		  tmp1 = minmod(vcell->theFieldsCoefficients->getCopy(q,d),
				(right->theFieldsCoefficients->get(q,km1_i)-vcell->theFieldsCoefficients->get(q,km1_i))/sqrt(double (2*k+1)/(2*k-1)), isLimited);
		
		if (flag || isLimited) {
		  vcell->limit=d-1;
		  if (i) vcell->theFieldsCoefficients->getCopy(q,d) = absmin(tmp1,tmp2);
		  else vcell->theFieldsCoefficients->getCopy(q,d) = tmp1;
		} 
		else {
		  if (vcell->limit != d) 
		    {
		      stopLimiter=1;
		      break;
		    }
		}
	      }
	    d--;
	  }
      }
  }
}

inline double dotprod(int q,double *ev,double *x)
{	
  int cSize = 4;
  double sum =0.0;
  for(int s=0; s<4; s++) 
	  sum +=ev[q*cSize+s]*x[s];
  return sum;
}

inline double rightprod(int q,double *ev,double *x)
{
  int cSize = 4;
  double sum =0.0;
  for(int s=0; s<4; s++) sum +=ev[s*cSize+q]*x[s];
  return sum;
}

////*********************LIMITER FOR EULER EQUATIONS***********************************/

void DGMomentLimiterEuler::limit(DGCell *vcell,double time)
{
  const double small_number = 1.0e-05;
  int d,k,i,q;
  int k_im1,km1_i,i_km1,im1_k;
  DGCell *left,*right,*top,*bottom;
  left=right=top=bottom=0;
  mPoint p1 = vcell->pMin;
  mPoint p2 = vcell->pMax;
  mPoint p=(p1+p2);
  p(0)/=2.;p(1)/=2.;p(2)/=2.;
  mPoint omin,omax;
  const int cSize = vcell->cSize;
  const int n = vcell->theMeshEntity->getLevel();
  for(i=0;i<vcell->theMeshEntity->size(n-1);i++)
    {
      mEntity *bound = vcell->theMeshEntity->get(n-1,i);
      DGCell *other=0;
      if(bound->get(n,0) != vcell->theMeshEntity)
	other = (DGCell*)bound->get(n,0)->getCell();
      else if (  bound->get(n,1)) other = (DGCell*)bound->get(n,1)->getCell();
      if (other)
	{
	  omin = other -> pMin;
	  omax = other -> pMax;
	  
	  if (p(0) > omax(0)) if (bound->getClassification()->getId() != 510) left = other; else right =other;
	  if (p(0) < omin(0)) if (bound->getClassification()->getId() != 510) right = other; else left =other;
	  if (p(1) > omax(1)) if (bound->getClassification()->getId() != 610) bottom = other; else top =other;
	  if (p(1) < omin(1)) if (bound->getClassification()->getId() != 610) top = other; else bottom =other;
	}	 
    }
  
  double tmp[4];
  bool isLimited[4], wasLimited[4], continue_limiting[4], stopLimiter,previous_coeff_limited[4];

  stopLimiter = 0;
  for (q=0;q<4;q++) continue_limiting[q] = 1;
  
  vcell->theFieldsCoefficients->copy();
 
  double lev_x[16],lev_y[16],rev_x[16],rev_y[16];
  mVector pp1(1,0,0);
  vcell->getConservationLaw()->compute_left_eigenvectors(vcell->theMean,(mVector&) pp1,(double *)lev_x);
  mVector pp2(0,1,0);
  vcell->getConservationLaw()->compute_left_eigenvectors(vcell->theMean,pp2,(double*) lev_y);
   vcell->getConservationLaw()->compute_right_eigenvectors(vcell->theMean,(mVector&)pp1,(double*)rev_x);
  vcell->getConservationLaw()->compute_right_eigenvectors(vcell->theMean,pp2,(double*)rev_y);
  
  d=(vcell->fOrder+1)*(vcell->fOrder+1)-1; //the highest index in the coefficient
  vcell->limit = d; //for plotting - the highest coefficient that was not limited
  
  double slope[4],bottom_slope[4],top_slope[4],left_slope[4],right_slope[4];

  for (k=vcell->fOrder; k>0; k--) {// limiting (p,p), (p,p-1), (p-1,p) etc
    if (stopLimiter) break;
    
    for(q=0; q<cSize; q++) 
      {
	slope[q] =vcell->theFieldsCoefficients->getCopy(q,d);
	if (top) top_slope[q] = top->theFieldsCoefficients->get(q,d-2)-vcell->theFieldsCoefficients->get(q,d-2);
	if (bottom) bottom_slope[q] = vcell->theFieldsCoefficients->get(q,d-2)-bottom->theFieldsCoefficients->get(q,d-2);
	if (left) left_slope[q] = vcell->theFieldsCoefficients->get(q,d-1)-left->theFieldsCoefficients->get(q,d-1);
	if (right) right_slope[q] = right->theFieldsCoefficients->get(q,d-1)-vcell->theFieldsCoefficients->get(q,d-1);
      } 
    
    double ev;
    bool slope_limited=0;
    for(q=0; q<cSize; q++) {                 //limit coef (k,k) in vertical direction
      ev = dotprod(q,lev_y,slope);
      if (continue_limiting[q] &&(fabs(ev) > small_number))
	{
	  if (top && bottom)
	    tmp[q] = minmod(ev, dotprod(q,lev_y,top_slope)/sqrt(double (2*k+1)/(2*k-1)), 
			      dotprod(q,lev_y,bottom_slope)/sqrt(double (2*k+1)/(2*k-1)), wasLimited[q]);
	  else if (top && !bottom) 
	    tmp[q] = minmod(ev, dotprod(q,lev_y,top_slope)/sqrt(double (2*k+1)/(2*k-1)), wasLimited[q]);
	  else if (!top&&bottom)
	    tmp[q] = minmod(ev, dotprod(q,lev_y,bottom_slope)/sqrt(double (2*k+1)/(2*k-1)), wasLimited[q]);
	}
      else {
		  tmp[q] = ev;
		  wasLimited[q] = 0;
	  }
    }
    
    slope_limited = 0;
    for (q=0; q<cSize; q++) if (wasLimited[q]) slope_limited = 1;
    if (slope_limited) // if a component was limited, update the slope on (k,k)
      for(q=0; q<cSize; q++) 
	slope[q]=vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_y,tmp);
    
    for(q=0; q<cSize; q++) {           //limit coef (k,k) in vertical direction
      ev = dotprod(q,lev_x,slope);
      if (continue_limiting[q] &&(fabs(ev) > small_number))
	{
	  if (left && right)
	    tmp[q] = minmod(ev, dotprod(q,lev_x,right_slope)/sqrt(double (2*k+1)/(2*k-1)),
			    dotprod(q,lev_x,left_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	  else if (left && !right)
	    tmp[q] = minmod(ev,dotprod(q,lev_x,left_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	  else if (!left && right)
	    tmp[q] = minmod(ev,dotprod(q,lev_x,right_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	  
	  if (wasLimited[q] || isLimited[q]) {
	    if (q==0)
		  vcell->limit=d-1;
	  } 
	  else  continue_limiting[q] = 0;
	} 
      else  {
		  tmp[q] = ev; 
		  isLimited[q] = 0;
	  }
    }
    
    stopLimiter = 1;
    for(q=0; q<cSize; q++) if (continue_limiting[q]) stopLimiter = 0;  
    if (stopLimiter) break;
    
    slope_limited  =0;
    for (q=0; q<cSize; q++) if (isLimited[q]) slope_limited = 1;
    if (slope_limited) 
      for(q=0; q<cSize; q++) 
	vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_x,tmp);
    
    d--; // put the counter on next coef
  
	for (i=k-1; i>=0; i--) //limit (i,k)
      {
	if (i<k-1) i_km1 =(k-1)*(k-1)+2*i+1; else i_km1 = i*i+(k-1)*2; 
	if (i>0) im1_k = k*k+2*(i-1)+1; 
	for(q=0; q<cSize; q++) 
	  {
	    slope[q] =vcell->theFieldsCoefficients->getCopy(q,d);
	    if (top) top_slope[q] = top->theFieldsCoefficients->get(q,i_km1)-vcell->theFieldsCoefficients->get(q,i_km1);
	    if (bottom) bottom_slope[q] = vcell->theFieldsCoefficients->get(q,i_km1)-bottom->theFieldsCoefficients->get(q,i_km1);
	    if (left && i) left_slope[q] = vcell->theFieldsCoefficients->get(q,im1_k)-left->theFieldsCoefficients->get(q,im1_k);
	    if (right && i) right_slope[q] = right->theFieldsCoefficients->get(q,im1_k)-vcell->theFieldsCoefficients->get(q,im1_k);
	  } 
	if (i) {//limit (i,k) in x dir
	  for(q=0; q<cSize; q++) {
	    ev = dotprod(q,lev_x,slope);
	    if (continue_limiting[q] && (fabs(ev) > small_number)) 
	      {if (left && right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,right_slope)/sqrt(double (2*i+1)/(2*i-1)),
				dotprod(q,lev_x,left_slope)/sqrt(double (2*i+1)/(2*i-1)),wasLimited[q]);
	      else if (left && !right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,left_slope)/sqrt(double (2*i+1)/(2*i-1)),wasLimited[q]);
	      else if (!left && right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,right_slope)/sqrt(double (2*i+1)/(2*i-1)), wasLimited[q]);
	      }
	    else tmp[q] = ev;
	  }
	  slope_limited = 0;
	  for (q=0; q<cSize; q++) if (wasLimited[q]) slope_limited = 1;
	    for(q=0; q<cSize; q++) 
	      vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_x,tmp);
	}
	else for(q=0; q<cSize; q++) wasLimited[q] = 0;
	
	for(q=0; q<cSize; q++) {//limit (i,k) in y dir
	  ev = dotprod(q,lev_y,slope);
	  if (continue_limiting[q] && (fabs(ev) > small_number)) {
	    if (top && bottom)
	      tmp[q] = minmod(ev, dotprod(q,lev_y,top_slope)/sqrt(double (2*k+1)/(2*k-1)),
			      dotprod(q,lev_y,bottom_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	    else if (top && !bottom) 
	      tmp[q] = minmod(ev,dotprod(q,lev_y,top_slope)/sqrt(double (2*k+1)/(2*k-1)), isLimited[q]);
	    else if (!top&&bottom)
	      tmp[q] = minmod(ev,dotprod(q,lev_y,bottom_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]); 
	    
	    if (wasLimited[q] || isLimited[q]) 
	      {
		if (q==0) 
			vcell->limit=d-1;
			previous_coeff_limited[q]=1;
	      } 
	  else previous_coeff_limited[q] = 0;
	  }
	  else tmp[q] = ev;
	}
	
	slope_limited  =0;
	for (q=0; q<cSize; q++) if (isLimited[q]) slope_limited = 1;
	if (slope_limited) 
	  for(q=0; q<cSize; q++) 
	    vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_y,tmp);
	
	d--;
	
	k_im1 = k*k+2*(i-1);
	km1_i = (k-1)*(k-1)+2*i;
	
	for(q=0; q<cSize; q++)
	  {
	    slope[q] =vcell->theFieldsCoefficients->getCopy(q,d);
	    if (top && i) top_slope[q] = top->theFieldsCoefficients->get(q,k_im1)-vcell->theFieldsCoefficients->get(q,k_im1);
	    if (bottom&& i) bottom_slope[q] = vcell->theFieldsCoefficients->get(q,k_im1)-bottom->theFieldsCoefficients->get(q,k_im1);
	    if (left) left_slope[q] = vcell->theFieldsCoefficients->get(q,km1_i)-left->theFieldsCoefficients->get(q,km1_i);
	    if (right) right_slope[q] = right->theFieldsCoefficients->get(q,km1_i)-vcell->theFieldsCoefficients->get(q,km1_i);
	  } 
	
	if (i) { //limit (k,i) in y direction
	  for(q=0; q<cSize; q++) {
	    ev = dotprod(q,lev_y,slope);
	    if (continue_limiting[q] && (ev>small_number))
	      {
		if (top && bottom &&i)
		  tmp[q] = minmod(ev,dotprod(q,lev_y,top_slope)/sqrt(double (2*i+1)/(2*i-1)),
				  dotprod(q,lev_y,bottom_slope)/sqrt(double (2*i+1)/(2*i-1)),wasLimited[q]);
		else if (top && !bottom && i) 
		  tmp[q] = minmod(ev,dotprod(q,lev_y,top_slope)/sqrt(double (2*i+1)/(2*i-1)), wasLimited[q]);
		else if (!top&&bottom && i)
		tmp[q] = minmod(ev,dotprod(q,lev_y,bottom_slope)/sqrt(double (2*i+1)/(2*i-1)),wasLimited[q]);
	      }
	    else tmp[q] = ev;
	  }
	  slope_limited = 0;
	  for (q=0; q<cSize; q++) if (wasLimited[q]) slope_limited = 1;
	  if (slope_limited) 
	    for(q=0; q<cSize; q++) 
	      vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_y,tmp);
	  
	  if (slope_limited)
	    for(q=0; q<cSize; q++) 
	      slope[q] =vcell->theFieldsCoefficients->getCopy(q,d);
	}
	else for (q=0; q<cSize; q++) wasLimited[q] = 0;
	
	for(q=0; q<cSize; q++) { //limit (k,i) in x direction)
	  ev = dotprod(q,lev_x,slope);
	  if (continue_limiting[q] &&(fabs(ev) > small_number))
	    {
	      if (left && right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,right_slope)/sqrt(double (2*k+1)/(2*k-1)),
				dotprod(q,lev_x,left_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	      else if (left && !right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,left_slope)/sqrt(double (2*k+1)/(2*k-1)),isLimited[q]);
	      else if (!left && right)
		tmp[q] = minmod(ev,dotprod(q,lev_x,right_slope)/sqrt(double (2*k+1)/(2*k-1)), isLimited[q]);
	      
     	      if (wasLimited[q] || isLimited[q]) {
		if (q==0)
			vcell->limit=d-1;
	      } 
	      else {
		if (!previous_coeff_limited[q]) 
		  continue_limiting[q] = 0;
	      }
	    } 
	  else tmp[q] = ev;
	}
	
	stopLimiter = 1;
	for(q=0; q<cSize; q++) if (continue_limiting[q]) stopLimiter = 0;  
	if (stopLimiter) break;
	slope_limited  =0;
	for (q=0; q<cSize; q++) if (isLimited[q]) slope_limited = 1;
	if (slope_limited) 
	  for(q=0; q<cSize; q++) 
	    vcell->theFieldsCoefficients->getCopy(q,d) = dotprod(q,rev_x,tmp);
	
	d--;
      }
  }
}


#define SIGN(a)	  (a < 0.0 ? -1 : 1)
double minmod(double a,double b,double c, bool &isLimited)
{
if (SIGN(a)==SIGN(b) && SIGN(b)==SIGN(c))
{
double abs_a = fabs(a);
double abs_b = fabs(b);
double abs_c = fabs(c);
if (abs_a<abs_b && abs_a < abs_c) {isLimited=0;return a;}
else if (abs_b<abs_c) {isLimited=1;return b;}
else {isLimited=1; return c;}

}
else {isLimited = 1; return 0.0;}
}

double minmod(double a,double b,bool &isLimited)
{
if (SIGN(a)==SIGN(b))
{
if (fabs(a)<fabs(b)) {isLimited = 0; return a;}
else  {isLimited = 1;return b;}
}
else {isLimited = 1;return 0.0;}
}

