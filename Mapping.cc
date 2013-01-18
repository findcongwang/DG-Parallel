#include "Mapping.h"
#include "mEntity.h"
#include "mVector.h"
#include "mTensor2.h"
#include "mPoint.h"
#include "Vertex.h"
#include <math.h>
#include <stdio.h>
MeshMapping::MeshMapping (mEntity *e)
  :Mapping(),ent(e)
{
}

//Note: If wanted, this can be more efficient by storing the detjac order as part of the class
//      or could just have all the geomOrder cases hardcoded into this function via case structure.
int MeshMapping::detJacOrder() const //returns polynomial degree of the determinant of jacobian
{ 
	switch(ent->getType())			
	{
		
	case mEntity::QUAD  : return 1 + (geomOrder() - 1)*(geomOrder() - 1); //is this correct?
	case mEntity::TET   : return (geomOrder() - 1)*(geomOrder() - 1)*(geomOrder() - 1);
	case mEntity::TRI   : return (geomOrder() - 1)*(geomOrder() - 1); 
	case mEntity::EDGE  : return (geomOrder() - 1)*(geomOrder() - 1);//GRANT: Need to verify this is right 
	default             : return 0;
    
	//default             : return 0;
	}
}
							
int MeshMapping::geomOrder() const  
{									
								
	return (ent->getgeomOrder());
}

double MeshMapping::GeomShapeFunction (int iNod, double u, double v, double w) const
{ 
  switch(ent->getType())
    {
    case mEntity::VERTEX : return 0;
    case mEntity::EDGE   : return GeomShapeFunctionLine(iNod,u,v,w);
    case mEntity::TRI    : return GeomShapeFunctionTri(iNod,u,v,w);   
    case mEntity::QUAD   : return GeomShapeFunctionQuad(iNod,u,v,w);
    case mEntity::TET    : return GeomShapeFunctionTet(iNod,u,v,w);
    case mEntity::HEX    : return GeomShapeFunctionHex(iNod,u,v,w);
    }
  cout<<"unknown type \n";
  return 0;
}

void MeshMapping::GradGeomShapeFunction (int iNod, double u, double v, double w,mVector &grad) const
{
  switch(ent->getType())
    {
    case mEntity::VERTEX : grad[0] = grad[1] = grad[2] = 0.0;break;
	case mEntity::EDGE   : GradGeomShapeFunctionLine (iNod,u,v,w,grad);break;
    case mEntity::TRI    : GradGeomShapeFunctionTri  (iNod,u,v,w,grad);break;
    case mEntity::QUAD   : GradGeomShapeFunctionQuad (iNod,u,v,w,grad);break;
    case mEntity::TET    : GradGeomShapeFunctionTet (iNod,u,v,w,grad);break;
    case mEntity::HEX    : GradGeomShapeFunctionHex (iNod,u,v,w,grad);break;
    }
}

double MeshMapping::GeomShapeFunctionQuad (int iNod, double u, double v, double dum) const
{
	//v0,v1,v2,v3 are mapped to (-1,-1),(-1,1),(1,1),(-1,1), respectively
 if(geomOrder() == 1)
    {
      switch(iNod)
	{
	case 0:return 0.25 * (1.-u) * (1.-v); 
	case 1:return 0.25 * (1.+u) * (1.-v);
	case 2:return 0.25 * (1.+u) * (1.+v);
	case 3:return 0.25 * (1.-u) * (1.+v);
	default: return 0.;
	}
    }
 if(geomOrder() == 2)
    {
      switch(iNod)
	{
	case 0:return  0.25 * u * v * (1.-u) * (1.-v); 
	case 1:return -0.25 * u * v * (1.+u) * (1.-v);  //doublecheck these
	case 2:return  0.25 * u * v * (1.+u) * (1.+v);
	case 3:return -0.25 * u * v * (1.-u) * (1.+v);
	case 4:return -0.5 * v * (1.-u*u) * (1.-v); 
	case 5:return 0.5 * u * (1.+u) * (1.-v*v);
	case 6:return 0.5 * v * (1.-u*u) * (1.+v);
	case 7:return -0.5 * u * (1.-u) * (1-v*v);
	case 8:return (1.-u*u) * (1.-v*v);
	default: return 0.;
	}
    }
  cout << "implementation of a quad mesh mapping of higher order than 1 is not done yet\n";
  return 0.;
} 

double MeshMapping::GeomShapeFunctionTet (int iNod, double u, double v, double w) const
{
	double x = 1.-u-v-w;
  if(geomOrder() == 1)
    {
      switch(iNod)
	{
	case 0:return x;
	case 3:return u;    //GRANT  //why is it ordered in this way?
	case 2:return v;
	case 1:return w;
	}
    }
  if(geomOrder() == 2)
    {
      switch(iNod)
	{
	case 0:return x*(2*x-1);
	case 1:return u*(2*u-1);
	case 2:return v*(2*v-1);
	case 3:return w*(2*w-1);
	case 4:return 4*u*x;
	case 5:return 4*u*v;
	case 6:return 4*v*x;
	case 7:return 4*w*x;
	case 8:return 4*v*w;
	case 9:return 4*u*w;
	}
    }
  cout << "implementation of a tet mesh mapping of higher order than 2 is not done yet\n";
  return 0.;
} 

double MeshMapping::GeomShapeFunctionHex (int iNod, double u, double v, double w) const
{
  if(geomOrder() == 1)
    {
      switch(iNod)
	{
	case 0:return 0.125 * (1.-u) * (1.-v) * (1.-w);
	case 1:return 0.125 * (1.+u) * (1.-v) * (1.-w);
	case 2:return 0.125 * (1.+u) * (1.+v) * (1.-w);
	case 3:return 0.125 * (1.-u) * (1.+v) * (1.-w);
	case 4:return 0.125 * (1.-u) * (1.-v) * (1.+w);
	case 5:return 0.125 * (1.+u) * (1.-v) * (1.+w);
	case 6:return 0.125 * (1.+u) * (1.+v) * (1.+w);
	case 7:return 0.125 * (1.-u) * (1.+v) * (1.+w);
	}
    }
  cout << "implementation of a hex mesh mapping of higher order than 1 is not done yet\n";
  return 0.0;
} 


double MeshMapping::GeomShapeFunctionTri (int iNod, double u, double v, double dum) const
{
	//vertices v0,v1,v2 are mapped to (0,0),(1,0),(0,1) respectively
  double w = 1.-u-v;
  if(geomOrder() == 1)
    {
      switch(iNod)
	{
	case 0:return (w);
	case 1:return (u);
	case 2:return (v);
	}
    }
  else if(geomOrder() == 2)
//vertices v0,v1,v2 are mapped to (0,0),(1,0),(0,1) respectively
//vertices v3,v4,v5 are mapped to (0.5,0),(.5,.5),(0,.5) respectively
    {
      switch(iNod)
	{
	case 0:return (w * (2.0*w-1));
	case 1:return (u * (2.0*u-1));
	case 2:return (v * (2.0*v-1));
	case 3:return (4.0*u*w);
	case 4:return (4.0*u*v);
	case 5:return (4.0*v*w);
	}
    }
  else if(geomOrder() == 3)
    {
      switch(iNod)
	{
	case 0:return ((9./2.)*w*(w-(1./3.))*(w-(2./3.))); //primary points
	case 1:return ((9./2.)*u*(u-(1./3.))*(u-(2./3.)));
	case 2:return ((9./2.)*v*(v-(1./3.))*(v-(2./3.)));
	case 3:return ((27./2.)*u*w*(w-(1./3.)));					//secondary points
	case 4:return ((27./2.)*u*w*(u-(1./3.)));
	case 5:return ((27./2.)*u*v*(u-(1./3.)));
	case 6:return ((27./2.)*u*v*(v-(1./3.)));
	case 7:return ((27./2.)*v*w*(v-(1./3.)));
	case 8:return ((27./2.)*v*w*(w-(1./3.)));
	case 9:return (27.*u*v*w);                    //tertiary point
	}
    }
  else if(geomOrder() == 4)
    {
      switch(iNod)
	{
	case 0:return  ((32./3.)*w*(w-(1./4.))*(w-(2./4.))*(w-(3./4.))); //primary points
	case 1:return  ((32./3.)*u*(u-(1./4.))*(u-(2./4.))*(u-(3./4.)));
	case 2:return  ((32./3.)*v*(v-(1./4.))*(v-(2./4.))*(v-(3./4.)));
	case 3:return  ((128./3.)*u*w*(w-(1./4.))*(w-(2./4.)));			//secondary points		
	case 4:return  ((64.)*u*w*(u-(1./4.))*(w-(1./4.)));
	case 5:return  ((128./3.)*u*w*(u-(1./4.))*(u-(2./4.)));
	case 6:return  ((128./3.)*u*v*(u-(1./4.))*(u-(2./4.)));
	case 7:return  ((64.)*u*v*(u-(1./4.))*(v-(1./4.)));
	case 8:return  ((128./3.)*u*v*(v-(1./4.))*(v-(2./4.)));
	case 9:return  ((128./3.)*v*w*(v-(1./4.))*(v-(2./4.)));                 
	case 10:return ((64.)*v*w*(v-(1./4.))*(w-(1./4.)));
	case 11:return ((128./3.)*v*w*(w-(1./4.))*(w-(2./4.)));
	case 12:return (128.*u*v*w*(w-(1./4.)));              //tertiary points
	case 13:return (128.*u*v*w*(u-(1./4.)));				
	case 14:return (128.*u*v*w*(v-(1./4.)));
	}
    }
  cout << "implementation of a mesh mapping of higher order than 4 is not done yet\n";
  return 0.0;
} 


void MeshMapping::GradGeomShapeFunctionTri (int iNod, double u, double v, double dum , mVector &grad) const
{
	//grad = (d/du N_i,d/dv N_i,0), where N is a geometric shape function 
	// see Joe's notes chapter 5.4

  double w = 1.-u-v;
  grad[2] = 0.0;
  if(geomOrder() == 1)
    {
      switch(iNod)
	{
	case 0: grad[0]=-1;grad[1]=-1;break;
	case 1: grad[0]=1.;grad[1]=0.;break;
	case 2: grad[0]=0.;grad[1]=1.;break;
	}
    }
  else if(geomOrder() == 2)
    {
      switch(iNod)
	{
	case 0: grad[0]=grad[1]=-4.0*w+1.0;           break;     
	case 1: grad[0]=4.0*u-1.0; grad[1]=0.0;       break;       
	case 2: grad[0]=0.0;       grad[1]=4.0*v-1.0; break;
	case 3: grad[0]=4.0*(w-u); grad[1]=-4.0*u;    break;        
	case 4: grad[0]=4.0*v;     grad[1]=4.0*u;     break; 
	case 5: grad[0]=-4.0*v;    grad[1]=4.0*(w-v); break;
	}
    }
  else if(geomOrder() == 3)
    {
      switch(iNod)
	{
	case 0: grad[0]=grad[1]=-(27./2.)*w*w+9.*w-1.;                                                    break;     
	case 1: grad[0]=(27./2.)*u*u-9.*u+1.;                grad[1]=0.0;                                 break;       
	case 2: grad[0]=0.0;                                 grad[1]=(27./2.)*v*v-9.*v+1.;                break;
	case 3: grad[0]=(27./2.)*(w*w+(1./3.)*(u-w)-2.*u*w); grad[1]=(27./2.)*(-2.*u*w+(1./3.)*u);        break;        
	case 4: grad[0]=(27./2.)*(-u*u+(1./3.)*(u-w)+2.*u*w);grad[1]=(-27./2.)*u*(u-1./3.);       break; 
	case 5: grad[0]=27.*u*v-(9./2.)*v;                   grad[1]=(27./2.)*u*(u-1./3.);                break;
	case 6: grad[0]=(27./2.)*v*(v-1./3.);                grad[1]=27.*u*v-(9./2.)*u;                   break;
	case 7: grad[0]=-(27./2.)*v*(v-1./3.);               grad[1]=-(27./2.)*v*v+(9./2.)*(v-w)+27.*w*v; break;
	case 8: grad[0]=(9./2.)*v-27.*v*w;                   grad[1]=(27./2.)*w*w+(9./2.)*(v-w)-27.*v*w;  break;
	case 9: grad[0]=27.*(w-u)*v;                         grad[1]=27.*(w-v)*u;                         break;
	}
    }
  else if(geomOrder() == 4)
    {
      switch(iNod)
	{
	case 0:  grad[0]=grad[1]=-(128./3.)*w*w*w+48.*w*w-(44./3.)*w+1.;                                                               break;     
	case 1:  grad[0]=(128./3.)*u*u*u-48.*u*u+(44./3.)*u-1.;             grad[1]=0.0;                                               break;       
	case 2:  grad[0]=0.0;                                               grad[1]=(128./3.)*v*v*v-48.*v*v+(44./3.)*v-1.;             break;
	case 3:  grad[0]=(16./3.)*(-24.*u*w*w+12.*u*w-u+8.*w*w*w-6.*w*w+w); grad[1]=(16./3.)*(-24.*u*w*w+12.*u*w-u);                   break;        
	case 4:  grad[0]=(4.)*(-32.*u*u*w+4.*u*u+32.*u*w*w-u-4.*w*w+w);  grad[1]=(4.)*(-32.*u*u*w+4.*u*u+8.*u*w-u);              break; 
	case 5:  grad[0]=(16./3.)*(-8.*u*u*u+24.*u*u*w+6.*u*u-12.*u*w-u+w); grad[1]=(-128./3.)*u*(u-1./4.)*(u-2./4.);                  break;
	case 6:  grad[0]=(16./3.)*(24.*u*u*v-12.*u*v+v);                    grad[1]=(128./3.)*u*(u-1./4.)*(u-2./4.);                   break;
	case 7:  grad[0]=(4.)*(32.*u*v*v-8.*u*v-4.*v*v+v);               grad[1]=(4.)*(32.*u*u*v-4.*u*u-8.*u*v+u);               break;
	case 8:  grad[0]=(128./3.)*v*(v-1./4.)*(v-2./4.);                   grad[1]=(16./3.)*(24.*u*v*v-12.*u*v+u);                    break;
	case 9:  grad[0]=(-128./3.)*v*(v-1./4.)*(v-2./4.);                  grad[1]=(16./3.)*(-8.*v*v*v+24.*v*v*w+6.*v*v-12.*v*w-v+w); break;
	case 10: grad[0]=(4.)*(-32.*v*v*w+4.*v*v+8.*v*w-v);              grad[1]=(4.)*(-32.*v*v*w+4.*v*v+32.*v*w*w-v-4.*w*w+w);  break;     
	case 11: grad[0]=(16./3.)*(-24.*v*w*w+12.*v*w-v);                   grad[1]=(16./3.)*(-24.*v*w*w+12.*v*w-v+8.*w*w*w-6.*w*w+w); break;       
	case 12: grad[0]=32.*(-8.*u*v*w+u*v+4.*v*w*w-v*w);                  grad[1]=32.*(-8.*u*v*w+u*v+4.*u*w*w-u*w);                  break;
	case 13: grad[0]=32.*(-4.*u*u*v+8.*u*v*w+u*v-v*w);                  grad[1]=32.*(-4.*u*u*v+4.*u*u*w+u*v-u*w);                  break;        
	case 14: grad[0]=128.*(v*w-u*v)*(v-1./4.);                          grad[1]=128.*((u*w-u*v)*(v-1./4.)+u*v*w);                  break; 
	}
    }
} 

void MeshMapping::GradGeomShapeFunctionQuad (int iNod, double u, double v, double dum , mVector &grad) const
{
  grad[2] = 0.0;
  if(geomOrder() == 1)
    {
      switch(iNod)
	{
	case 0:grad[0] = -0.25 * (1.-v) ; grad[1]= -0.25 * (1.-u);break;
	case 1:grad[0] =  0.25 * (1.-v) ; grad[1]= -0.25 * (1.+u);break;
	case 2:grad[0] =  0.25 * (1.+v) ; grad[1]=  0.25 * (1.+u);break;
	case 3:grad[0] = -0.25 * (1.+v) ; grad[1]=  0.25 * (1.-u);break;
	}
    }
  if(geomOrder() == 2)
    {
      switch(iNod)
	{
	case 0:grad[0] = 0.25 * v * (1.-v) * (1.-2*u) ;   grad[1]= 0.25 * u * (1.-u) * (1.-2*v);   break;
	case 1:grad[0] =  -0.25 * v * (1.-v) * (1.+2*u) ; grad[1]= 0.25 * u * (1.+u) * (-1.+2*v);  break;
	case 2:grad[0] =  0.25 * v * (1.+v) * (1.+2*u) ;  grad[1]=  0.25 * u * (1.+u) * (1.+2*v);  break;
	case 3:grad[0] = 0.25 * v * (1.+v) * (-1.+2*u) ;  grad[1]=  -0.25 * u * (1.-u) * (1.+2*v); break;
	case 4:grad[0] = v * u * (1.-v) ;                 grad[1]= 0.5 * (1.-u*u) * (-1.+2*v);     break;
	case 5:grad[0] =  0.5 * (1.-v*v) * (1+2*u) ;      grad[1]= -u * v * (1.+u);                break;
	case 6:grad[0] = -u * v * (1.+v) ;                grad[1]=  0.5 * (1.-u*u) * (1.+2*v);     break;
	case 7:grad[0] = 0.5 * (1.-v*v) * (-1.+2*u) ;     grad[1]=  u * v * (1.-u);                break;
	case 8:grad[0] = -2*u * (1.-v*v) ;                grad[1]= (1.-u*u)*(-2*v);                break;
	}
    }
} 

void MeshMapping::GradGeomShapeFunctionTet (int iNod, double u, double v, double w , mVector &grad) const
{
  double x = 1.-u-v-w;
  if(geomOrder() == 1)
    {
      switch(iNod)
	{		
	case 0:grad[0] = -1. ; grad[1]= -1.; grad[2] = -1.;break;
	case 1:grad[0] =  0. ; grad[1]=  0.; grad[2] =  1.;break;
	case 2:grad[0] =  0. ; grad[1]=  1.; grad[2] =  0.;break;
	case 3:grad[0] =  1. ; grad[1]=  0.; grad[2] =  0.;break;
	}
    }
  if(geomOrder() == 2)
    {
      switch(iNod)
	{
	case 0:grad[0] = -4.*x+1.;  grad[1]= -4.*x+1.;  grad[2] = -4.*x+1.; break;
	case 1:grad[0] =  4.*u-1.;  grad[1]=  0.;       grad[2] =  0.;      break;
	case 2:grad[0] =  0. ;      grad[1]=  4.*v-1.;  grad[2] =  0.;      break;
	case 3:grad[0] =  0. ;      grad[1]=  0.;       grad[2] =  4.*w-1.; break;
	case 4:grad[0] =  4.*(x-u); grad[1]=  -4.*u;    grad[2] =  -4.*u;   break;
	case 5:grad[0] =  4.*v;     grad[1]=  4.*u;     grad[2] =  0.;      break;
	case 6:grad[0] =  -4.*v;    grad[1]=  4.*(x-v); grad[2] =  -4.*v;   break;
	case 7:grad[0] =  -4.*w;    grad[1]=  -4.*w;    grad[2] =  4.*(x-w);break;
	case 8:grad[0] =  0. ;      grad[1]=  4.*w;     grad[2] =  4.*v;    break;
	case 9:grad[0] =  4.*w;     grad[1]=  0.;       grad[2] =  4.*u;    break;
	}
    }
} 

void MeshMapping::GradGeomShapeFunctionHex (int iNod, double u, double v, double w , mVector &grad) const
{
  if(geomOrder() == 1)
    {
      switch(iNod)
	{
	case 0:grad[0] = -0.125*(1.-v)*(1.-w) ; grad[1]= -0.125*(1.-u)*(1.-w); grad[2] = -0.125*(1.-u)*(1.-v);break;
	case 1:grad[0] =  0.125*(1.-v)*(1.-w) ; grad[1]= -0.125*(1.+u)*(1.-w); grad[2] = -0.125*(1.+u)*(1.-v);break;
	case 2:grad[0] =  0.125*(1.+v)*(1.-w) ; grad[1]=  0.125*(1.+u)*(1.-w); grad[2] = -0.125*(1.+u)*(1.+v);break;
	case 3:grad[0] = -0.125*(1.+v)*(1.-w) ; grad[1]=  0.125*(1.-u)*(1.-w); grad[2] = -0.125*(1.-u)*(1.+v);break;
	case 4:grad[0] = -0.125*(1.-v)*(1.+w) ; grad[1]= -0.125*(1.-u)*(1.+w); grad[2] =  0.125*(1.-u)*(1.-v);break;
	case 5:grad[0] =  0.125*(1.-v)*(1.+w) ; grad[1]= -0.125*(1.+u)*(1.+w); grad[2] =  0.125*(1.+u)*(1.-v);break;
	case 6:grad[0] =  0.125*(1.+v)*(1.+w) ; grad[1]=  0.125*(1.+u)*(1.+w); grad[2] =  0.125*(1.+u)*(1.+v);break;
	case 7:grad[0] = -0.125*(1.+v)*(1.+w) ; grad[1]=  0.125*(1.-u)*(1.+w); grad[2] =  0.125*(1.-u)*(1.+v);break;
	}
    }
} 

//computational to physical
void MeshMapping::eval(double u, double v, double w, 
					   double &x, double &y, double &z) const
{
  x = 0.0, y = 0.0, z = 0.0;

  for(int i=0;i<ent->size(0);i++)
    {
      double f = GeomShapeFunction (i,u,v,w);
      mPoint p = ((Vertex*)(ent->get(0,i)))->point();
      x += p(0) * f;
      y += p(1) * f;
      z += p(2) * f;
    }
}
bool MeshMapping::invert(double xp, double yp, double zp, double &Upos, double &Vpos, double &Wpos) const
{
	//converts physical coordinates to computational
#define NR_PRECISION       1.e-12
#define NR_MAX_ITER        50
	
  double   x_est, y_est, z_est;
  double   u_new, v_new, w_new;
  double   Error = 1.0 ;
  int      iter = 0 ;
	
  // Upos = Vpos = Wpos = 0.0;
  COG(Upos,Vpos,Wpos);	
 // Upos=1./3.;
 // Vpos=1./3.;
 // Wpos=1./3.;
  mTensor2 InvJacMatrix;
	//Newton iterations for nonlinear mappings
// jac given by (4.5.4) in Joe's notes 
  while (Error > NR_PRECISION && iter < NR_MAX_ITER){
		
    iter++ ;

    jacInverse(Upos,Vpos,Wpos,InvJacMatrix); 
    eval(Upos,Vpos,Wpos,x_est,y_est,z_est);
      
    u_new = Upos + InvJacMatrix(0,0) * (xp-x_est) + InvJacMatrix(1,0) * (yp-y_est) +
      InvJacMatrix(2,0) * (zp-z_est) ;
    v_new = Vpos + InvJacMatrix(0,1) * (xp-x_est) + InvJacMatrix(1,1) * (yp-y_est) +
      InvJacMatrix(2,1) * (zp-z_est) ;
    w_new = Wpos + InvJacMatrix(0,2) * (xp-x_est) + InvJacMatrix(1,2) * (yp-y_est) +
      InvJacMatrix(2,2) * (zp-z_est) ;
    
    Error = (u_new - Upos) * (u_new - Upos) + 
      (v_new - Vpos) * (v_new - Vpos) + 
      (w_new - Wpos) * (w_new - Wpos) ;
    
    Upos = u_new;
    Vpos = v_new;
    Wpos = w_new;
  }
  if(Error > NR_PRECISION)
    {
      printf("impossible to find %f %f %f in ",xp,yp,zp);ent->print();
      for(int i=0;i<ent->size(0);i++)
	ent->get(0,i)->print();
   
      return false;
    }
  return true;
}

// certainly not the fastest way !!
double MeshMapping::detJac(double u, double v, double w) const 
{
  mTensor2 t;
  return jacInverse(u,v,w,t);
}

double MeshMapping::jacInverse(double u, double v, double w, mTensor2 &Invjac) const 
{

  double jac[3][3];
  jac[0][0] = jac[1][1] = jac[0][1] = jac[1][0] = 0.0;
  jac[1][2] = jac[2][1] = jac[0][2] = jac[2][0] = jac[2][2] = 0.0;
  int i;
  double DetJac;
  mVector du(0,0,0);
  switch(ent->getLevel())
    {
      // may take into account non x-y elements and curved elements, curved lines,...
    case 1:
      {
	for(i = 0;i<ent->size(0);i++)
	{
	  mPoint p = ((Vertex*)ent->get(0,i))->point();
	  
	  GradGeomShapeFunction  (i , u ,  v ,  w, du);
	  jac[0][0] += p(0) * du[0];
	  jac[0][1] += p(1) * du[0];
	  jac[0][2] += p(2) * du[0];
	}
	DetJac = sqrt ( jac[0][0]*jac[0][0] + jac[0][1]*jac[0][1]
			+ jac[0][2]*jac[0][2] );
	
	int vrai=0;
	for(int i=0;i<3;i++){
	  if(jac[0][i]==0){   // if this component of the normal vector is zero
	    vrai =1;
	    for(int j=0;j<3;j++){
	      if (j==i)
		jac[1][j] = 1; // then the component of the second normal must be one,
	      else
		jac[1][j] = 0; // and the other components of the second normal are zero.
	    }
	    continue;
	  }
	}                           
	// surface equation with normal vector n : n1*X + n2*Y + n3*Z = 0
	if(!vrai){ // looking for the second normal vector in the plane z = 0
	  double temp = sqrt ( jac[0][0]*jac[0][0] + jac[0][1]*jac[0][1] );
	  jac[1][0] = -jac[0][1]/temp;
	  jac[1][1] =  jac[0][0]/temp;
	  jac[1][2] =  0;
	}
	// The third normal vector
	jac[2][0] = ( jac[0][1]*jac[1][2] - jac[0][2]*jac[1][1] )/ DetJac;
	jac[2][1] = ( jac[0][2]*jac[1][0] - jac[0][0]*jac[1][2] )/ DetJac;
	jac[2][2] = ( jac[0][0]*jac[1][1] - jac[0][1]*jac[1][0] )/ DetJac; 
      }
      break;
    case 2:
      {
	//	printf("size = %d\n",ent->size(0));
	for(int i=0;i<ent->size(0);i++)
	{
	    GradGeomShapeFunction(i,u,v,w,du);
	    mPoint p = ((Vertex*)ent->get(0,i))->point();
	    //	    printf("%f %f %f %f %f %f\n",u,v,du[0],du[1],p(0),p(1));
	    jac[0][0] += p(0) *du[0];
	    jac[1][0] += p(0) *du[1];
	    jac[0][1] += p(1) *du[0];
	    jac[1][1] += p(1) *du[1];
	    jac[0][2] += p(2) *du[0];
	    jac[1][2] += p(2) *du[1];
	  }
	double d3 = jac[0][0]*jac[1][1] - jac[0][1] * jac[1][0];  
	double d2 = jac[0][2]*jac[1][0] - jac[0][0] * jac[1][2];  
	double d1 = jac[0][1]*jac[1][2] - jac[0][2] * jac[1][1];  
	
	DetJac = sqrt ( d1*d1 + d2*d2 + d3*d3 );
	
	jac[2][0] = d1/DetJac; 
	jac[2][1] = d2/DetJac; 
	jac[2][2] = d3/DetJac;
      }
      break;
    case 3:
      {
	for(int i=0;i<ent->size(0);i++)
	  {
	    GradGeomShapeFunction(i,u,v,w,du);
	    mPoint p = ((Vertex*)ent->get(0,i))->point();
	    jac[0][0] += p(0) *du[0];
	    jac[1][0] += p(0) *du[1];
	    jac[2][0] += p(0) *du[2];
	    jac[0][1] += p(1) *du[0];
	    jac[1][1] += p(1) *du[1];
	    jac[2][1] += p(1) *du[2];
	    jac[0][2] += p(2) *du[0];
	    jac[1][2] += p(2) *du[1];
	    jac[2][2] += p(2) *du[2];
	  }
	DetJac = (jac[0][0]*jac[1][1]*jac[2][2] + jac[0][2] *jac[1][0]*jac[2][1] +
		      jac[0][1]*jac[1][2]*jac[2][0] - jac[0][2] *jac[1][1]*jac[2][0] -
		      jac[0][0]*jac[1][2]*jac[2][1] - jac[0][1] *jac[1][0]*jac[2][2]);
      }
      break;
    }
  /*
      printf("\t|%5.2f %5.2f %5.2f|\n",jac[0][0],jac[0][1],jac[0][2]);
    printf("M=\t|%5.2f %5.2f %5.2f|\n",jac[1][0],jac[1][1],jac[1][2]);
    printf("\t|%5.2f %5.2f %5.2f|\n",jac[2][0],jac[2][1],jac[2][2]);
   */
    

  Invjac(0,0) = (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) / DetJac;
  Invjac(1,0) = -(jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]) / DetJac;
  Invjac(2,0) = (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]) / DetJac;
  
  Invjac(0,1) = -(jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1]) / DetJac;
  Invjac(1,1) = (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]) / DetJac;
  Invjac(2,1) = -(jac[0][0] * jac[2][1] - jac[0][1] * jac[2][0]) / DetJac;
  
  Invjac(0,2) = (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) / DetJac;
  Invjac(1,2) = -(jac[0][0] * jac[1][2] - jac[0][2] * jac[1][0]) / DetJac;
  Invjac(2,2) = (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]) / DetJac;

  /*printf("\t|%f %f %f|\n",Invjac(0,0),Invjac(0,1),Invjac(0,2));
  printf("\t|%f %f %5f|\n",Invjac(1,0),Invjac(1,1),Invjac(1,2));
  printf("\t|%f %f %f|\n\n",Invjac(2,0),Invjac(2,1),Invjac(2,2));*/
  
  return DetJac;
}

double MeshMapping::GeomShapeFunctionLine (int iNod, double u, double v, double w) const
{
  if(geomOrder() == 1)
    {
      switch(iNod)
	{
	case 0:return 0.5 * (1.+ u);
	case 1:return 0.5 * (1.- u);
	}
    }
  else if(geomOrder() == 2)
    {
      switch(iNod)
	{
	case 0:return (0.5 * (u*u - u));
	case 1:return (0.5 * (u*u + u));
	case 2:return (1. - u*u);
	}
	}
  else if(geomOrder() == 3)
    {
      switch(iNod)
	{
	case 0:return (-9./16.)*(u-1.)*(u+1./3.)*(u-1./3.);
	case 1:return (9./16.)*(u+1.)*(u+1./3.)*(u-1./3.);
	case 2:return (27./16.)*(u+1.)*(u-1.)*(u-1./3.);
	case 3:return (-27./16.)*(u+1.)*(u-1.)*(u+1./3.);
	}
    }
  else if(geomOrder() == 4)
    {
      switch(iNod)
	{
	case 0:return (2./3.)*(u+0.5)*u*(u-0.5)*(u-1.);
	case 1:return (2./3.)*(u+1.)*(u+0.5)*u*(u-0.5);
	case 2:return (-8./3.)*(u+1.)*u*(u-0.5)*(u-1.);
	case 3:return (4.)*(u+1.)*(u+0.5)*(u-0.5)*(u-1.);
	case 4:return (-8./3.)*(u+1.)*(u+0.5)*u*(u-1.);
	}
    }
   cout << "implementation of a mesh mapping of higher order than 4 is not done yet\n";
   return 0.0;
} 

void MeshMapping::GradGeomShapeFunctionLine (int iNod, double u, double v, double w, mVector &grad) const
{
  grad[1] = grad[2] = 0.0;
  if(geomOrder() == 1)
    {
      switch(iNod)
	{
	case 0:grad[0] = (+0.5);break;
	case 1:grad[0] = (-0.5);break;
	}
    }
  else if(geomOrder() == 2)
    {
      switch(iNod)
	{
	case 0:grad[0] = (u-0.5);break;
	case 1:grad[0] = (u+0.5);break;
	case 2:grad[0] = (-2.*u);break;
	}
    }
  else if(geomOrder() == 3)
    {
      switch(iNod)
	{
	case 0:grad[0] = (-9./16.) * (3.*u*u-2.*u-(1./9.));break;
	case 1:grad[0] = (9./16.) * (3.*u*u+2.*u-(1./9.));break;
	case 2:grad[0] = (27./16.) * (3.*u*u-(2./3.)*u-(1.));break;
	case 3:grad[0] = (-27./16.) * (3.*u*u+(2./3.)*u-(1.));break;
	}
    }
  else if(geomOrder() == 4)
    {
      switch(iNod)
	{
	/*case 0:grad[0] = (2./18.)*(16.*u*u*u-4.*u*u-2.*u+1.);break;
	case 1:grad[0] = (2./18.)*(16.*u*u*u+4.*u*u-2.*u-1.);break;
	case 2:grad[0] = (-8./9.)*(8.*u*u*u-3.*u*u-4.*u+1.); break;
	case 3:grad[0] = (4./3.)*(8.*u*u*u-5.*u);            break;
	case 4:grad[0] = (-8./9.)*(8.*u*u*u+3.*u*u-4.*u-1.); break;*/
	case 0:grad[0] = (1./6.)*(16.*u*u*u-12.*u*u-2.*u+1.);break;
	case 1:grad[0] = (1./6.)*(16.*u*u*u+12.*u*u-2.*u-1.);break;
	case 2:grad[0] = (-4./3.)*(8.*u*u*u-3.*u*u-4.*u+1.); break;
	case 3:grad[0] = (2.)*u*(8.*u*u-5.);                 break;
	case 4:grad[0] = (-4./3.)*(8.*u*u*u+3.*u*u-4.*u-1.); break;
	}
    }
  else cout << "implementation of a mesh mapping of higher order than 3 is not done yet\n";
} 
void MeshMapping::normalVector (mEntity *border , double u, double v, double w, mVector &n) const
{

  // In case of a mesh mapping, the normal to an entity which is on
  // a border = sum of grad geom shape functions which node are not
  // on the border entity

  n[0] = n[1] = n[2] = 0.0;
  
  mTensor2 jInv;
  double detJac = jacInverse(u,v,w,jInv);
  int size = ent->size(0);
  for(int i=0;i<size;i++)
    {
      Vertex *ver = (Vertex*)ent->get(0,i);
      if(!border->find(ver))
	{
		//cout<<i<<endl; 
	  mVector gr;
	  GradGeomShapeFunction(i,u,v,w,gr);
	  gr *= jInv;
	  n-= gr;
	}
    }
  n.normalize();
}

ConstantMeshMapping::ConstantMeshMapping(mEntity *m)
  : MeshMapping(m)
{
  det = MeshMapping::jacInverse(0,0,0,jac);
   if(det < 0)det = -det;
  
  //  printf("\t|%5.2f %5.2f %5.2f|\n",jac(0,0),jac(0,1),jac(0,2));
  //  printf("M=\t|%5.2f %5.2f %5.2f|\n",jac(1,0),jac(1,1),jac(1,2));
  //  printf("\t|%5.2f %5.2f %5.2f|\n",jac(2,0),jac(2,1),jac(2,2));
  
}

CylindricalCoordinatesMeshMapping::CylindricalCoordinatesMeshMapping(mEntity *m)
  : MeshMapping(m)
{  
}
double ConstantMeshMapping::detJac(double,double,double) const
{
  return det;
}

double ConstantMeshMapping::jacInverse(double,double,double,mTensor2 &j) const
{
  j = jac;
  return det;
}

double CylindricalCoordinatesMeshMapping::detJac(double u,double v,double w) const
{
  mTensor2 j;
  double r,z,teta;
  eval(u,v,w,r,z,teta);
  //  printf("r = %12.5E det = %12.5E\n",r,MeshMapping::jacInverse(u,v,w,j) );
  return (r * MeshMapping::jacInverse(u,v,w,j));
}

double CylindricalCoordinatesMeshMapping::jacInverse(double u,double v,double w,mTensor2 &j) const
{
  double r,z,teta;
  eval(u,v,w,r,z,teta);
  double det = MeshMapping::jacInverse (u,v,w,j);
  return (r * det);
}

int CylindricalCoordinatesMeshMapping::detJacOrder() const
{ 
  return 1; 
}

bool MeshMapping::inReferenceElement(double u, double v, double w) const
{
  switch(ent->getType())
    {
    case mEntity::TRI :
      if(u < 0 || v < 0.0 || 1.-u-v < 0.0)return false;
      break;
    case mEntity::QUAD :
      if(u < -1. || u > 1 || v < -1.0 || v > 1.0)return false;
      break;
	case mEntity::EDGE  :
      if(u < -1. || u > 1)return false;
      break;
    case mEntity::TET :
      if(u < 0 || v < 0.0 || w < 0.0 || 1.-u-v-w < 0.0)return false;
      break;
    case mEntity::HEX :
      if(u < -1. || u > 1 || v < -1.0 || v > 1.0 || w < -1.0 || w > 1.0)return false;
      break;
    case mEntity::WEDGE :
      if(u < 0 || v < 0.0 || 1.-u-v < 0.0 || w < -1.0 || w > 1.0)return false;
      break;
      
    }
  return true;
}


void MeshMapping::COG(double &u, double &v, double &w) const
{
  switch(ent->getType())
    {
    case mEntity::WEDGE :
	case mEntity::TRI :    
      u = v = 0.3333333333333;
      w = 0;
      break;
    case mEntity::HEX :
    case mEntity::QUAD :
      u = v =w =0.0;
      break;
	case mEntity::EDGE  :  
      u = v = w = 0.0;
      break;
    case mEntity::TET :
      u = v = w = 0.25;
      break;
    }
}

