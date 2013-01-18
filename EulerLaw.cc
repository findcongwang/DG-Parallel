#include "EulerLaw.h"
#include "mVector.h"
#include "mPoint.h"
#include "FieldEvaluator.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <cstring>

Euler::Euler (set<int> &w, FieldEvaluator *f, double Gamma)
  : Walls(w), gamma(Gamma)
{
  farField = f;
  exactField = f;
}

Euler2d::Euler2d (set<int> &Walls, FieldEvaluator *f, double Gamma)
  : Euler(Walls,f,Gamma)
{}

Euler2dRoe::Euler2dRoe (set<int> &Walls, FieldEvaluator *f, double Gamma)
  : Euler2d(Walls,f,Gamma)
{}

Euler2dHLLC::Euler2dHLLC (set<int> &Walls, FieldEvaluator *f, double Gamma)
  : Euler2d(Walls,f,Gamma)
{}

EulervanLeer::EulervanLeer (set<int> &Walls, FieldEvaluator *f, double Gamma)
  : Euler2d(Walls,f,Gamma)
{}

Euler2dExactRiemann::Euler2dExactRiemann (set<int> &Walls, FieldEvaluator 
*f, double Gamma)
  : Euler2d(Walls,f, Gamma)
{}

Euler3d::Euler3d (set<int> &Walls, FieldEvaluator *f, double Gamma)
  : Euler(Walls,f,Gamma)
{}

Euler3dRoe::Euler3dRoe (set<int> &Walls, FieldEvaluator *f, double Gamma)
  : Euler3d(Walls,f,Gamma)
{}

void Euler2d::Fi ( const mPoint &position , double *Q, mVector *flux) const
{
  double invrho = 1./Q[0];
  double Q1 = Q[1]; 
  double Q2 = Q[2];
  double q12 = Q1*Q2*invrho;
  double q11 = Q1*Q1*invrho;
  double q22 = Q2*Q2*invrho;
  double p = (gamma-1.)*(Q[3] - 0.5*(q11+q22));
  double qq = invrho*(Q[3]+p);
  flux[0](0) = Q1;
  flux[0](1) = Q2;
  flux[0](2) = 0.0;
  flux[1](0) = q11 + p;
  flux[1](1) = q12;
  flux[1](2) = 0.0;
  flux[2](0) = q12;
  flux[2](1) = q22 + p;
  flux[2](2) = 0.0;
  flux[3](0) = Q1*qq;
  flux[3](1) = Q2*qq;
  flux[3](2) = 0.0;
}

inline double computeP2d (const double* Q, double gamma)
{
  return (gamma -1.)*(Q[3] - 0.5*(Q[1]*Q[1]+Q[2]*Q[2])/Q[0]);
}

inline double computeP3d (const double* Q, double gamma=1.4)
{
  return (gamma -1.)*(Q[3] - 0.5*(Q[1]*Q[1]+Q[2]*Q[2]+Q[4]*Q[4])/Q[0]);
}

void Euler3d::Fi ( const mPoint &position , double *Q, mVector *flux) const
{
  double invrho = 1./Q[0];
  double q11 = Q[1]*Q[1]*invrho;
  double q12 = Q[1]*Q[2]*invrho;
  double q14 = Q[1]*Q[4]*invrho;
  double q22 = Q[2]*Q[2]*invrho;
  double q24 = Q[2]*Q[4]*invrho;
  double q44 = Q[4]*Q[4]*invrho;
  double p = computeP3d(Q);
  flux[0](0) = Q[1];
  flux[0](1) = Q[2];
  flux[0](2) = Q[4];
  flux[1](0) = q11+p;
  flux[1](1) = q12;
  flux[1](2) = q14;
  flux[2](0) = q12;
  flux[2](1) = q22+p;
  flux[2](2) = q24;
  flux[4](0) = q14;
  flux[4](1) = q24;
  flux[4](2) = q44+p;
  double qq = invrho*(Q[3]+p);
  flux[3](0) = Q[1]*qq;
  flux[3](1) = Q[2]*qq;
  flux[3](2) = Q[4]*qq;
}

void Euler2d::boundary (mVector &n, mVector &N, int BoundaryId, const mPoint 
			    &position,
			    double *Q, double *flux, double T) const
{
  if ( BoundaryId == 20000 ) 
    {
      double uright_n,uright_t,uright,vright;;
      double f[4];
      double gm1 = gamma-1.;
      double p = gm1*(Q[3] - 0.5*(Q[1]*Q[1]+Q[2]*Q[2])/Q[0]);
      
      //ALG I starts here
      /*      
      //set v dot N = 0 
      uright_t = (N(1)*Q[1]-N(0)*Q[2])/Q[0];
      uright = uright_t * N(1);
      vright = - uright_t * N(0);
      
      //compute new normal velicity and energy 
      uright_n = uright*n(0) + vright*n(1);
      double E = p/gm1+(uright*uright+vright*vright)*.5*Q[0];
      
      //compute flux
      flux[0] = Q[0]*uright_n;;
      flux[1] = Q[0]*uright*uright_n + n(0) * p;
      flux[2] = Q[0]*vright*uright_n +n(1) * p;
      flux[3] = (E+p)*uright_n; 
       */         
      //ALG I ends here
      
      //ALG II starts here
      
      double uright_N = -(N(0)*Q[1]+N(1)*Q[2])/Q[0];
      double uright_T = (N(1)*Q[1]-N(0)*Q[2])/Q[0];
      uright = uright_N*N(0) + uright_T * N(1);
      vright = uright_N*N(1) - uright_T * N(0);

      f[0] = Q[0];
      f[1] = uright*Q[0];
      f[2] = vright*Q[0];
      f[3] = Q[3];
      
      riemannSolver(n,position,Q,f,flux);
      //ALG II ends here
      
      //ALG III starts here
      /*double u_n = (n(0)*Q[1]+n(1)*Q[2])/Q[0];
      double u_t = (n(1)*Q[1]-n(0)*Q[2])/Q[0];
      double cosalpha = N(0)*n(0)+N(1)*n(1);
      double sinalpha = sqrt(1.-cosalpha*cosalpha);
      double tanalpha = sinalpha/cosalpha;
      if (fabs(cosalpha-1.)<0.0000000000001) tanalpha =0.0;
	  uright_t = u_t;

      if (u_n>0) uright_n = u_n-2.*(u_n-fabs(u_t)*tanalpha);
      if (u_n<0) uright_n = u_n-2.*(u_n+fabs(u_t)*tanalpha);	
      uright = uright_n*n(0)+uright_t*n(1);
      vright = uright_n*n(1)- uright_t * n(0);
      
      f[0] = Q[0];
      f[1] = uright*Q[0];
      f[2] = vright*Q[0];
      f[3] = Q[3];
      
      riemannSolver(n,position,Q,f,flux);*/
      
      //ALG III ends here

      return; 
    }
}


void Euler2d::boundaryFlux (mVector &n, int BoundaryId, const mPoint 
			    &position,
			    double *Q, double *flux, double T) const
{
 
  if (BoundaryId==30000)
    {
      const double n0 = n(0);
      const double n1 = n(1);
      double Qinf[4];
      farField->eval(position,T,Qinf);
      const double gm1 = gamma-1.;
      const double ogm1 = 1./gm1;
      const double inv_rho = 1./Q[0];
      const double inv_rho_inf = 1./Qinf[0];
      double p    = gm1*(Q[3] - 0.5*(Q[1]*Q[1]+Q[2]*Q[2])*inv_rho);
      double pinf = gm1*(Qinf[3] - 0.5*(Qinf[1]*Qinf[1]+Qinf[2]*Qinf[2])*inv_rho_inf);
      double vn     = (n0*Q[1]+n1*Q[2])*inv_rho;
      double vn_inf = (n0*Qinf[1]+n1*Qinf[2])*inv_rho_inf;
      double c     = sqrt(gamma*p*inv_rho);

      bool is_inflow  =(vn<0)?1:0;
      
      if (fabs(vn)>=c)           /*  supersonic */
	{                  
	  if (is_inflow) /*SPECIFY_FREESTREAM;*/
	    {
	      flux[0] = vn_inf  * Qinf[0]; 
	      flux[1] = vn_inf  * Qinf[1] + pinf * n0;
	      flux[2] = vn_inf  * Qinf[2] + pinf * n1;
	      flux[3] = vn_inf  * (pinf*ogm1 + 
		   0.5*(Qinf[1]*Qinf[1]+Qinf[2]*Qinf[2])*inv_rho_inf+pinf);
	    }
	  else      /*SIMPLE_EXTRAPOLATION;*/ 
	    {
	      flux[0] = vn * Q[0];
	      flux[1] = vn * Q[1] + p * n0;
	      flux[2] = vn * Q[2] + p * n1;
	      flux[3] = vn * (p*ogm1 + 
			      0.5*(Q[1]*Q[1]+Q[2]*Q[2])*inv_rho+p);
	    }
	}
      else 
	{             /* ...Subsonic Far field, use Riemann invariant BCs */
	  if (is_inflow)
	    { 
	      double Qb0,Qb1,Qb2,Qb3,vt;
	      double Jplus, Jminus, vn_bnd, c_bnd, p_bnd;
	      double c_inf = sqrt(gamma*pinf*inv_rho_inf);
	      Jplus   = -vn_inf + 2.*c_inf*ogm1; 
	      Jminus  = -vn  - c*2.*ogm1;
	      
	      vn_bnd  = -0.5 * (Jplus + Jminus); 
	      c_bnd   = 0.25 * (Jplus - Jminus)*gm1;
	      vt = (n1*Qinf[1]-n0*Qinf[2])*inv_rho_inf;
	      
	      Qb0 = Qinf[0]*pow(c_bnd/c_inf , 2.*ogm1);
	      Qb1 = Qb0*(n0*vn_bnd +n1*vt);
	      Qb2 = Qb0*(n1*vn_bnd -n0*vt);
	      p_bnd = Qb0*c_bnd*c_bnd/gamma; 
	      Qb3 = p_bnd*ogm1 + 0.5*(Qb1*Qb1+Qb2*Qb2)/Qb0;
	      
	      flux[0] = Qb0 * vn_bnd;
	      flux[1] = vn_bnd * Qb1 + p_bnd * n0;
	      flux[2] = vn_bnd * Qb2 + p_bnd * n1;
	      flux[3] = vn_bnd * (Qb3 + p_bnd);
	      
	    }
	  else 
	    {                               /*...is outflow */
	      double Qb0,Qb1,Qb2,Qb3,vt; 
	      double Jplus, Jminus, vn_bnd,c_bnd, p_bnd;
	      double c_inf = sqrt(gamma*pinf*inv_rho_inf);

	      Jminus  = vn_inf  - 2.*c_inf*ogm1; 
	      Jplus   = vn      + 2.*c*ogm1;
	      
	      vn_bnd  = 0.5 *(Jplus + Jminus); 
	      c_bnd   = 0.25 *(Jplus - Jminus)*gm1;
	      vt = (n1*Q[1]-n0*Q[2])*inv_rho;
	      
	      Qb0 = Q[0]*pow(c_bnd/c , 2.*ogm1);
	      Qb1 = Qb0*(n0*vn_bnd +n1*vt);
	      Qb2 = Qb0*(n1*vn_bnd-n0*vt);
	      p_bnd = Qb0*c_bnd*c_bnd/gamma; 
	      //p_bnd = p*pow(c_bnd/c , 2.*ogm1*gamma); 
	      Qb3 = p_bnd*ogm1 + 0.5*(Qb1*Qb1+Qb2*Qb2)/Qb0;
	     
	      flux[0] = Qb0 * vn_bnd;
	      flux[1] = vn_bnd * Qb1 + p_bnd * n0;
	      flux[2] = vn_bnd * Qb2 + p_bnd * n1;
	      flux[3] = vn_bnd * (Qb3 + p_bnd);
	   
	    }
	}
      return;
    }
  
  if (BoundaryId==40000)
    {
      double f[4];
      double vn = -(Q[1]*n(0) + Q[2]*n(1));   //moments, not velocity
      double vt = Q[1]*n(1) - Q[2]*n(0);
      f[0] = Q[0];
      f[1] = vn*n(0) + vt*n(1);
      f[2] = vn*n(1) - vt*n(0);
      f[3] = Q[3];
      riemannSolver(n,position,Q,f,flux);
      return;
    }
 
  double f[4];
  farField->eval(position,T,f);
  riemannSolver(n,position,Q,f,flux);
}

void Euler3d::boundary (mVector &n, mVector &N, int BoundaryId, const mPoint 
			    &position,
			    double *Q, double *flux, double T) const
{
  if ( BoundaryId == 20000 ) 
    {
      // New boundary conditions - Algorithm 2	  
      
      double N0=N(0);
      double N1=N(1);
      double N2=N(2);
      double scale;
      
      // Compute mutually orthonormal vectors tangental to physical surface.
      
      // Tangent Vector 1
      double T0, T1, T2;
      T0=N1+N2;
      T1=-N0;
      T2=-N0;
      scale = 1.0/sqrt(T0*T0+T1*T1+T2*T2);
      T0=T0*scale;
      T1=T1*scale;
      T2=T2*scale;
      
      // Tangent Vector 2
      double S0, S1, S2;
      S0=N0*(N2-N1);
      S1=N0*N0+N2*(N1+N2);
      S2=-N0*N0-N1*(N1+N2);
      scale = 1.0/sqrt(S0*S0+S1*S1+S2*S2);
      S0=S0*scale;
      S1=S1*scale;
      S2=S2*scale;
          
      // Compute the moments in the (N,S,T) reference frame, but negate the 
      // N component.     
      double VN, VT, VS;
      VN =-(Q[1]*N0+Q[2]*N1+Q[4]*N2); 
      VT =  Q[1]*T0+Q[2]*T1+Q[4]*T2;
      VS =  Q[1]*S0+Q[2]*S1+Q[4]*S2;
      
      // Rotate back to the (U,V,W) reference frame to get the ghost momemtes.
      // Remember that inverse is transpose - since othonormal vectors.
      double UG, VG, WG; 
      UG = VN*N0+VT*T0+VS*S0;
      VG = VN*N1+VT*T1+VS*S1;
      WG = VN*N2+VT*T2+VS*S2;
      
      // Create the ghost state
      double f[5];
      f[0] = Q[0];  //ghost state
      f[1] = UG;
      f[2] = VG;
      f[3] = Q[3];
      f[4] = WG;

      // Solve Riemann Problem
      riemannSolver(n,position,Q,f,flux);
      //printf("%e %e %e %e %e \n", flux[0],flux[1],flux[2],flux[3],flux[4]);
      return;
  
      // *** OLD BOUNDARY CONDITIONS ***
      /*
      double vn, vt,vs, scale;
      const double n0=n(0);
      const double n1=n(1);
      const double n2=n(2);
      
      double t0=n1+n2;  //first tangential vector
      double t1=-n0;
      double t2=-n0;
      scale = 1./sqrt(t0*t0+t1*t1+t2*t2);
      t0 = t0 * scale; t1 = t1 * scale; t2 = t2 * scale;

      double s0= n0*(n2-n1);  //second tangential vector
      double s1= n0*n0+n2*(n1+n2);
      double s2 = -n0*n0-n1*(n1+n2);
      scale = 1./sqrt(s0*s0+s1*s1+s2*s2);
      s0 = s0 * scale; s1 = s1 * scale; s2 = s2 * scale;
            
      vn =-(Q[1]*n0+Q[2]*n1+Q[4]*n2); //Moments, not velocities!
      vt =  Q[1]*t0+Q[2]*t1+Q[4]*t2;
      vs =  Q[1]*s0+Q[2]*s1+Q[4]*s2;
      
      double f[5];
      double ug, vg, wg; //ghost moments
      ug = vn*n0+vt*t0+vs*s0;
      vg = vn*n1+vt*t1+vs*s1;
      wg = vn*n2+vt*t2+vs*s2;
      
      f[0] = Q[0];  //ghost state
      f[1] = ug;
      f[2] = vg;
      f[3] = Q[3];
      f[4] = wg;
      riemannSolver(n,position,Q,f,flux);
	  //printf("%e %e %e %e %e \n", flux[0],flux[1],flux[2],flux[3],flux[4]);
      return;
      */

      double p = computeP3d(Q);  
      flux[0] = 0.0;
      flux[1] = n[0] * p;
      flux[2] = n[1] * p;
      flux[3] = 0.0;
      flux[4] = n[2] * p;
      //printf("%e %e %e %e %e \n", flux[0],flux[1],flux[2],flux[3],flux[4]);
      return;
    }
}


void Euler3d::boundaryFlux (mVector &n, int BoundaryId, const mPoint 
&position,
			    double *Q, double *flux, double T) const
{
  if (BoundaryId==30000) //needs moew work on tangential components
    {
      const double n0 = n(0);
      const double n1 = n(1);
      const double n2 = n(2);
      double Qinf[5];
      farField->eval(position,T,Qinf);
      const double gm1 = gamma-1.;
      const double ogm1 = 1./gm1;
      const double inv_rho = 1./Q[0];
      const double inv_rho_inf = 1./Qinf[0];
      double p    = gm1*(Q[3] - 0.5*(Q[1]*Q[1]+Q[2]*Q[2]+Q[4]*Q[4])*inv_rho);
      double pinf = gm1*(Qinf[3] - 0.5*(Qinf[1]*Qinf[1]+Qinf[2]*Qinf[2]+Qinf[4]*Qinf[4])*inv_rho_inf);
      double vn     = (n0*Q[1]+n1*Q[2]+n2*Q[4])*inv_rho;
      double vn_inf = (n0*Qinf[1]+n1*Qinf[2]+n2*Qinf[4])*inv_rho_inf;
      double c     = sqrt(gamma*p*inv_rho);
      
      bool is_inflow  =(vn<0)?1:0;
      
      if (fabs(vn)>=c)           /*  supersonic */
	{                  
	  if (is_inflow) /*SPECIFY_FREESTREAM;*/
	    {
			 flux[0] = vn_inf  * Qinf[0]; 
	      flux[1] = vn_inf  * Qinf[1] + pinf * n0;
	      flux[2] = vn_inf  * Qinf[2] + pinf * n1;
	      flux[3] = vn_inf  * (pinf*ogm1 + 
		   0.5*(Qinf[1]*Qinf[1]+Qinf[2]*Qinf[2])*inv_rho_inf+pinf);
	      flux[4] = vn_inf  * Qinf[4] + pinf * n2;
	    }
	  else      /*SIMPLE_EXTRAPOLATION;*/ 
	    {
	      flux[0] = vn * Q[0];
	      flux[1] = vn * Q[1] + p * n0;
	      flux[2] = vn * Q[2] + p * n1;
	      flux[3] = vn * (p*ogm1 + 
			      0.5*(Q[1]*Q[1]+Q[2]*Q[2])*inv_rho+p);
	      flux[4] = vn * Q[4] + p * n2;
	    } 
	}
      else 
	{             /* ...Subsonic Far field, use Riemann invariant BCs */
	  if (is_inflow)
	    { 
	      double Qb0,Qb1,Qb2,Qb3,Qb4,vt,vs;
	      double Jplus, Jminus, vn_bnd, c_bnd, p_bnd;
	      double c_inf = sqrt(gamma*pinf*inv_rho_inf);
	      Jplus   = -vn_inf + 2.*c_inf*ogm1; 
	      Jminus  = -vn  - c*2.*ogm1;
	      
	      vn_bnd  = -0.5 * (Jplus + Jminus); 
	      c_bnd   = 0.25 * (Jplus - Jminus)*gm1;
	      double t0=n1+n2;
	      double t1=-n0;
	      double t2=-n0;
	      double scale = 1./sqrt(t0*t0+t1*t1+t2*t2);
	      t0 = t0 * scale; t1 = t1 * scale; t2 = t2 * scale;
	      
	      double s0= n0*(n2-n1);
	      double s1= n0*n0+n2*(n1+n2);
	      double s2 = -n0*n0-n1*(n1+n2);
	      scale = 1./sqrt(s0*s0+s1*s1+s2*s2);
	      s0 = s0 * scale; s1 = s1 * scale; s2 = s2 * scale;
	      vt =  (Qinf[1]*t0+Qinf[2]*t1+Qinf[4]*t2)*inv_rho_inf;
	      vs =  (Qinf[1]*s0+Qinf[2]*s1+Qinf[4]*s2)*inv_rho_inf;
	      
	      Qb0 = Qinf[0]*pow(c_bnd/c_inf , 2.*ogm1);
	      Qb1 = Qb0*(n0*vn_bnd+t0*vt+s0*vs);
	      Qb2 = Qb0*(n1*vn_bnd+t1*vt+s1*vs);
	      Qb4 = Qb0*(n2*vn_bnd+t2*vt+s2*vs);
	      p_bnd = Qb0*c_bnd*c_bnd/gamma; 
	      Qb3 = p_bnd*ogm1 + 0.5*(Qb1*Qb1+Qb2*Qb2+Qb4*Qb4)/Qb0;
	      
	      flux[0] = Qb0 * vn_bnd;
	      flux[1] = vn_bnd * Qb1 + p_bnd * n0;
	      flux[2] = vn_bnd * Qb2 + p_bnd * n1;
	      flux[3] = vn_bnd * (Qb3 + p_bnd);
	      flux[4] = vn_bnd *Qb4 +p_bnd*n2;
	      
	    }
	  else 
	    {                               /*...is outflow */
	      double Qb0,Qb1,Qb2,Qb3,Qb4,vt,vs; 
	      double Jplus, Jminus, vn_bnd,c_bnd, p_bnd;
	      double c_inf = sqrt(gamma*pinf*inv_rho_inf);

	      Jminus  = vn_inf  - 2.*c_inf*ogm1; 
	      Jplus   = vn      + 2.*c*ogm1;
	      
	      vn_bnd  = 0.5 *(Jplus + Jminus); 
	      c_bnd   = 0.25 *(Jplus - Jminus)*gm1;
	      double t0=n1+n2;
	      double t1=-n0;
	      double t2=-n0;
	      double scale = 1./sqrt(t0*t0+t1*t1+t2*t2);
	      t0 = t0 * scale; t1 = t1 * scale; t2 = t2 * scale;
	      
	      double s0= n0*(n2-n1);
	      double s1= n0*n0+n2*(n1+n2);
	      double s2 = -n0*n0-n1*(n1+n2);
	      scale = 1./sqrt(s0*s0+s1*s1+s2*s2);
	      s0 = s0 * scale; s1 = s1 * scale; s2 = s2 * scale;
	      vt =  (Q[1]*t0+Q[2]*t1+Q[4]*t2)*inv_rho;
	      vs =  (Q[1]*s0+Q[2]*s1+Q[4]*s2)*inv_rho;
	      
	      Qb0 = Q[0]*pow(c_bnd/c , 2.*ogm1);
	      Qb1 = Qb0*(n0*vn_bnd+t0*vt+s0*vs);
	      Qb2 = Qb0*(n1*vn_bnd+t1*vt+s1*vs);
	      Qb4 = Qb0*(n2*vn_bnd+t2*vt+s2*vs);
	      p_bnd = Qb0*c_bnd*c_bnd/gamma; 
	      Qb3 = p_bnd*ogm1 + 0.5*(Qb1*Qb1+Qb2*Qb2+Qb4*Qb4)/Qb0;
	     
	      flux[0] = Qb0 * vn_bnd;
	      flux[1] = vn_bnd * Qb1 + p_bnd * n0;
	      flux[2] = vn_bnd * Qb2 + p_bnd * n1;
	      flux[3] = vn_bnd * (Qb3 + p_bnd);
	      flux[4] = vn_bnd *Qb4 +p_bnd*n2;
	    }
	}
      return;
    }

  if(BoundaryId == 40000)
    {
      double vn, vt,vs, scale;
      const double n0=n(0);
      const double n1=n(1);
      const double n2=n(2);
      
      double t0=n1+n2;
      double t1=-n0;
      double t2=-n0;
      scale = 1./sqrt(t0*t0+t1*t1+t2*t2);
      t0 = t0 * scale; t1 = t1 * scale; t2 = t2 * scale;

      double s0= n0*(n2-n1);
      double s1= n0*n0+n2*(n1+n2);
      double s2 = -n0*n0-n1*(n1+n2);
      scale = 1./sqrt(s0*s0+s1*s1+s2*s2);
      s0 = s0 * scale; s1 = s1 * scale; s2 = s2 * scale;
            
      vn =-(Q[1]*n0+Q[2]*n1+Q[4]*n2); //Moments, not velocities!
      vt =  Q[1]*t0+Q[2]*t1+Q[4]*t2;
      vs =  Q[1]*s0+Q[2]*s1+Q[4]*s2;
      
      double f[5];
      double ug, vg, wg; //ghost moments
      ug = vn*n0+vt*t0+vs*s0;
      vg = vn*n1+vt*t1+vs*s1;
      wg = vn*n2+vt*t2+vs*s2;
      
      f[0] = Q[0];  //ghost state
      f[1] = ug;
      f[2] = vg;
      f[3] = Q[3];
      f[4] = wg;
      riemannSolver(n,position,Q,f,flux);
            
      //printf("%e %e %e %e %e \n", flux[0],flux[1],flux[2],flux[3],flux[4],flux[5]);
      return;
      double p = computeP3d(Q);
      flux[0] = 0.0;
      flux[1] = n[0] * p;
      flux[2] = n[1] * p;
      flux[3] = 0.0;
      flux[4] = n[2] * p;
      //  printf("%e %e %e %e %e \n", flux[0],flux[1],flux[2],flux[3],flux[4],flux[5]);
      return;
    }
	
/*if(BoundaryId == 40000)
    {
	double p = computeP3d(Q);  
      flux[0] = 0.0;
      flux[1] = n[0] * p;
      flux[2] = n[1] * p;
      flux[3] = 0.0;
      flux[4] = n[2] * p;
      //printf("%e %e %e %e %e \n", flux[0],flux[1],flux[2],flux[3],flux[4]);
      return;
}*/
  double f[5];
  farField->eval(position,T,f);
  riemannSolver(n,position,Q,f,flux);
}


#define RIEM_TOL   ((double)1e-14)       /* Rel error to achieve in pstar */
#define RIEM_ITER  (12)                   /* Max number of iterations */

#define smallp ((double)1e-14)

#define max2(a,b)  ( (a)>(b) ? (a) : (b) )
#define max3(a,b,c) ( max2( a , max2(b,c) ) )

int Euler::riemann( double rhol,                /* inputs */
		     double  ulft,
		     double  vlft,
		     double  wlft,
		     double  plft,
		     double rhor,
		     double  urght,
		     double  vrght,
		     double  wrght,
		     double  prght,
		     double *rhoav,               /* outputs */
		     double  *uav,
		     double  *vav,
		     double  *wav,
		     double  *pav ) const
{
  
  const double rhol_inv= 1./ rhol; 
  const double rhor_inv= 1./ rhor; 
  const double gm1 = gamma-1.;
  int nonPhysical=0;
/*  const double rhol =1.0;
  const double rhol_inv = 1./1.;
  const double ulft =0.;
  const double vlft=0;
  const double plft=.01;
  const double rhor=1.; 
  const double rhor_inv =1./1.;
  const double  urght=0.;const double vrght=0;
  const double prght=100.;*/
  
  double clft,crght,tmp1;
  tmp1 = gamma*plft*rhol_inv;                 
  if(tmp1 < 0.0) 
  {
    printf("Nonphysical values found in the Riemann solver rhol = %f pl = %f  \n",rhol, plft);
    tmp1=fabs(tmp1);
    nonPhysical = 1;
    //  exit(0);
  }
  clft=sqrt(tmp1);
     
  tmp1 = gamma*prght*rhor_inv;
  if(tmp1 < 0.0)
  {
    printf("Nonphysical values found in the Riemann solver rhor = %f pr = %f  \n",rhor, prght);
    tmp1=fabs(tmp1);
    if (nonPhysical==1) nonPhysical = 3; else nonPhysical = 2;
    //exit(0);
  }
  crght=sqrt(tmp1);
 
  if (ulft>=clft && urght>=crght) //completely supersonic in left state
    {
      *rhoav = rhol;
      *uav  = ulft;
      *vav  = vlft;
      *wav  = wlft;
      *pav   = plft;
      return nonPhysical;
    }
  else if (urght<=-crght && ulft<=-clft) //comletely supersonic in right state
    {
      *rhoav = rhor;
      *uav  = urght;
      *vav  = vrght;
      *wav  = wrght;
      *pav   = prght;
      return nonPhysical;
    }
  else {
  int i;
  double tmp2,scrch1;
  double pstar, ustar, rhostar, cstar;
  double rhos, us, ps,cs;
  double prevpstar, relpstar;
  double flft,dflft,frght,dfrght;
 
  const double invgamma = 1./gamma;
  const double gp1  = gamma+1.;
  const double C1       = 0.5*(gp1)*invgamma;
  const double C3       = 0.5*gm1*invgamma; 
  const double C4       = 2./gm1;
  
  // Initial guess for pstar by two shock approximation; Toro p.128
  double udiff = urght - ulft;
  double p0 = 0.5*(plft+prght) -0.125*udiff*(rhol+rhor)*(clft+crght);
  p0 = max2(smallp,p0);
  double gl = sqrt(2./(rhol*(gm1*plft+gp1*p0)));
  double gr = sqrt(2./(rhor*(gm1*prght+gp1*p0)));
  pstar = (gl*plft+gr*prght-udiff)/(gl+gr);
  
  relpstar=RIEM_TOL;
  pstar=max2(RIEM_TOL,pstar);
  prevpstar=pstar;                           /* save intial value */
  
  for (i=0; i<RIEM_ITER && relpstar>=RIEM_TOL ; i++)
    {
      if (pstar>plft)
	{
	  tmp1 = 1./(gp1*pstar+gm1*plft);
	  tmp2 = sqrt(2.*tmp1*rhol_inv);
	  flft = (pstar-plft)*tmp2;
	  dflft = tmp2*(1.-0.5*gp1*(pstar-plft)*tmp1);
	}
      else
	{
	  tmp2 = pow(pstar/plft,C3);
	  flft = C4*clft*(tmp2-1.);
	  dflft = tmp2*clft*invgamma/pstar;
	}
      if (pstar>prght)	
	{
	  tmp1 = 1./(gp1*pstar+gm1*prght);
	  tmp2 = sqrt(2.*tmp1*rhor_inv);
	  frght = (pstar-prght)*tmp2;
	  dfrght = tmp2*(1.-0.5*gp1*(pstar-prght)*tmp1); 	  
	}
      else
	{
	  tmp2 = pow(pstar/prght,C3);
	  frght = C4*crght*(tmp2-1.);
	  dfrght = tmp2*crght*invgamma/pstar;
	}
      pstar -= (flft+frght+udiff)/(dflft+dfrght);
      pstar=max2(smallp,pstar);
      
      relpstar=fabs(pstar-prevpstar)/(.5*(prevpstar+pstar));// compute relative error
      prevpstar=pstar; 
    }
  
  if (relpstar>=RIEM_TOL)
    {
      cout << "Riemann solver failed to converge in " << RIEM_ITER << "iterations\n";
      cout << "Relative error is " << relpstar << "\n";
    }
  //Recompute F_l  and F_r with the final pstar to get ustar
  if (pstar>plft)
    flft = (pstar-plft)*sqrt(2.*rhol_inv/(gp1*pstar+gm1*plft));
  else
    flft = C4*clft*(pow(pstar/plft,C3)-1.);
  if (pstar>prght)	
    frght = (pstar-prght)*sqrt(2.*rhor_inv/(gp1*pstar+gm1*prght));
  else
    frght = C4*crght*(pow(pstar/prght,C3)-1.);
  
  ustar = (ulft+urght+frght-flft)*0.5;
  
  if (ustar>=0) //choose left state
    {
      ps=plft;
      us=ulft;  *vav= vlft; *wav = wlft;  
      rhos=rhol;
      cs=clft;
      scrch1 = 1.;
    }
  else         //choose right state
    {
      ps=prght;
      us=urght; *vav= vrght; *wav = wrght;  
      rhos=rhor;
      cs=crght;
      scrch1=-1.;
    }
  
  if ( pstar >= ps  ) //shock, choose star state
    {
      double C= gm1/gp1;
      tmp1 = pstar/ps;
      *rhoav = rhos*(tmp1+C)/(C*tmp1+1.);
      *uav   = ustar;
      *pav   = pstar;
    }// rarefaction, choose star or rarefaction state
  else //rarefaction
    { 
      rhostar = rhos*pow(pstar/ps,invgamma);
      cstar = sqrt(gamma*pstar/rhostar);
      
      if ((ustar-scrch1*cstar)*ustar<=0.0) 
	//in brackets speed of the rarefaction tail 
	{//choose star state
	  *rhoav = rhostar;
	  *uav   = ustar;
	  *pav   = pstar;
	}	
      else
	{//choose a state inside rarefaction wave
	  double C = 2./gp1;
	  tmp1 =C+scrch1*gm1*us/(cs*gp1); 
	  *rhoav = rhos*pow(tmp1,C4);
	  *uav   = C*(gm1*0.5*us+scrch1*cs); 
	  *pav   = ps*pow(tmp1,C4*gamma );
	}
    }
  if ((clft+crght)*2.<=gm1*(urght-ulft))
    {
      printf("VACUUM\n");
    }
  return nonPhysical;
  }
}

void Euler2d::riemannSolver( mVector &n ,
				const mPoint &position,
				double *u0,
				double *u1,
				double *flux) const
{
  double vn1 = n(0) * (u0[1] / u0[0]) + n(1) * (u0[2] / u0[0]);
  double vt1 = n(1) * (u0[1] / u0[0]) - n(0) * (u0[2] / u0[0]);

  double vn2 = n(0) * (u1[1] / u1[0]) + n(1) * (u1[2] / u1[0]);
  double vt2 = n(1) * (u1[1] / u1[0]) - n(0) * (u1[2] / u1[0]);

  double p1 = (gamma-1)*u0[3] - 
0.5*(gamma-1)*(u0[1]*u0[1]+u0[2]*u0[2])/u0[0];
  double p2 = (gamma-1)*u1[3] - 
0.5*(gamma-1)*(u1[1]*u1[1]+u1[2]*u1[2])/u1[0];

  double vnav,vtav,dummy,pav;
  double riemann_average[4];

  riemann (u0[0],vn1,vt1,0.0,p1,
	   u1[0],vn2,vt2,0.0,p2,
	   & riemann_average[0],
	   & vnav,
	   & vtav,
	   & dummy,
	   & pav);

  riemann_average[1] = riemann_average[0] * (n(0) * vnav + n(1) * vtav);
  riemann_average[2] = riemann_average[0] * (n(1) * vnav - n(0) * vtav);
  riemann_average[3] = (1./(gamma-1))*(pav + 
0.5*(gamma-1)*riemann_average[0]*(vnav*vnav+vtav*vtav));

  mVector f[4];
  Fi(position,riemann_average,f);
  flux[0] = f[0] * n;
  flux[1] = f[1] * n;
  flux[2] = f[2] * n;
  flux[3] = f[3] * n;
}

void Euler3d::riemannSolver( mVector &n ,
			     const mPoint &position,
			     double *Ql,
			     double *Qr,
			     double *flux) const
{
  static int count = 0;
  count++;
  if (count ==170)
    {
		double dnj=1.0;
    }
  
  const double n0=n(0);
  const double n1=n(1);
  const double n2=n(2);
  
  mVector t(n1+n2,-n0,-n0);  
  if (t.L2norm()<0.0000001) {t(0)=-n1;t(1)=n0+n2;t(2)=-n1;}
  if (t.L2norm()<0.0000001) {t(0)= -n2;t(1)=-n2;t(2)=n0+n1;}
  t.normalize();

  mVector s (n % t);
  s.normalize();
  
  double vnl = (Ql[1]*n0+Ql[2]*n1+Ql[4]*n2)/Ql[0];
  double vtl = (Ql[1]*t(0)+Ql[2]*t(1)+Ql[4]*t(2))/Ql[0];
  double vsl = (Ql[1]*s(0)+Ql[2]*s(1)+Ql[4]*s(2))/Ql[0];

  double vnr = (Qr[1]*n0+Qr[2]*n1+Qr[4]*n2)/Qr[0];
  double vtr = (Qr[1]*t(0)+Qr[2]*t(1)+Qr[4]*t(2))/Qr[0];
  double vsr = (Qr[1]*s(0)+Qr[2]*s(1)+Qr[4]*s(2))/Qr[0]; 
 
  double pl = computeP3d(Ql);
  double pr = computeP3d(Qr);

  double rhoav,vnav,vtav,vsav,pav;
  /*
  printf("left : %f %f %f %f %f\n",u0[0],u0[1],u0[2],u0[3],u0[4]);
  printf("left : %f %f %f %f %f\n",u1[0],u1[1],u1[2],u1[3],u1[4]);
  */

  riemann (Ql[0],vnl,vtl,vsl,pl,
	   Qr[0],vnr,vtr,vsr,pr,
	   & rhoav,
	   & vnav,
	   & vtav,
	   & vsav,
	   & pav);

  flux[0] = rhoav * vnav;
  flux[1] = flux[0] * (n0*vnav + t(0)*vtav + s(0)*vsav) + pav*n0;
  flux[2] = flux[0] * (n1*vnav + t(1)*vtav + s(1)*vsav) + pav*n1;
  flux[4] = flux[0] * (n2*vnav + t(2)*vtav + s(2)*vsav) + pav*n2;
  flux[3] = vnav * (pav/(gamma-1.)+pav+0.5*rhoav*(vnav*vnav+vtav*vtav+vsav*vsav));
  
 // printf("%f %f %f %f %f\n", flux[0],flux[1],flux[2],flux[3],flux[4]);
  
}

void Euler2dRoe::riemannSolver( mVector &n ,
			     const mPoint &position,
			     double *u0,
			     double *u1,
			     double *flux) const
{
  double rho_A, u_A, v_A, H_A;
  // density
  double rho_0 = u0[0], rho_1 = u1[0];
  double sqrt_rho_0 = sqrt(rho_0), sqrt_rho_1 = sqrt(rho_1);
  double inv_rho_0 = 1./rho_0, inv_rho_1 = 1./rho_1;
 
  rho_A = sqrt_rho_0 * sqrt_rho_1;
  // u velocity
  double inv_denom = 1./(sqrt_rho_0 + sqrt_rho_1);
  double u_0 = u0[1]*inv_rho_0, u_1 = u1[1]*inv_rho_1;
  u_A = (u_0 * sqrt_rho_0 + u_1 * sqrt_rho_1)*inv_denom;
  // v velocity
  double v_0 = u0[2]*inv_rho_0, v_1 = u1[2]*inv_rho_1;
  v_A = (v_0 * sqrt_rho_0 + v_1 * sqrt_rho_1)*inv_denom;
  double v_n_0 = u_0 * n(0) + v_0 * n(1);
  double v_n_1 = u_1 * n(0) + v_1 * n(1);
  double p_0 = (gamma-1.)*(u0[3] - 0.5*(u0[1]*u0[1]+u0[2]*u0[2])*inv_rho_0);
  double p_1 = (gamma-1.)*(u1[3] - 0.5*(u1[1]*u1[1]+u1[2]*u1[2])*inv_rho_1);

  double H_0 = (u0[3]+p_0)*inv_rho_0, H_1 = (u1[3]+p_1)*inv_rho_1;
  H_A = (H_0 * sqrt_rho_0 + H_1 * sqrt_rho_1)*inv_denom;

  // speed of sound
  double qq = H_A - 0.5*(u_A * u_A + v_A * v_A);

  double roe[4];
  roe[0] = roe[1] = roe[2] = roe[3] = 0.0;

  if(qq > 0.0)
    {
      double c_A = sqrt((gamma-1.0) * qq);
      // normal velocity
      double vn = u_A * n(0) + v_A * n(1);

      double evec[4][4];
      double eval[4];

      // eigen values
      eval[0] = fabs(vn);
      eval[1] = fabs(vn);
      eval[2] = fabs(vn+c_A);
      eval[3] = fabs(vn-c_A);
      // eigen vectors

      evec[0][0] = 1.;
      evec[0][1] = u_A;
      evec[0][2] = v_A;
      evec[0][3] = 0.5*(u_A*u_A+v_A*v_A);//H_A-qq

      evec[1][0] = 0.;
      evec[1][1] = n(1) * rho_A;
      evec[1][2] = -n(0) * rho_A;
      evec[1][3] = rho_A * (u_A * n(1) - v_A * n(0));

      double k = 0.5*rho_A/c_A;
      evec[2][0] = k;
      evec[2][1] = k*(u_A + c_A * n(0));
      evec[2][2] = k*(v_A + c_A * n(1));
      evec[2][3] = k*(H_A + c_A * vn);

      evec[3][0] = k;
      evec[3][1] = k*(u_A - c_A * n(0));
      evec[3][2] = k*(v_A - c_A * n(1));
      evec[3][3] = k*(H_A - c_A * vn);

      // characteristic variables

      double charv[4];
	  double inv_c_A = 1./c_A;
      double drho = rho_1-rho_0;
      double dp = p_1-p_0;
      double du = u_1-u_0;
      double dv = v_1-v_0;
      charv[0] = drho - dp*inv_c_A*inv_c_A;
      charv[1] = n(1) * du - n(0) * dv;
      charv[2] = n(0) * du + n(1) * dv + dp/rho_A * inv_c_A;
      charv[3] = -n(0)*du - n(1) * dv + dp/rho_A * inv_c_A;

      // finally, the flux

      for(int i=0;i<4;i++)
	    for(int j=0;j<4;j++)
	      roe[j] += eval[i] * charv[i] * evec[i][j];
    }
	else  printf("unphysical data's found on %f %f %f\n",position(0),position(1),position(2));
mVector f1[4],f2[4];
  Fi(position,u0,f1);
  Fi(position,u1,f2);
  // printf("roe %e %e %e %e \n",roe[0],roe[1],roe[2],roe[3]);
  //printf("normal vel in rieman solver %e \n",(roe[1]*n(0)+ roe[2]*n(1))/roe[0]);
  for(int i=0;i<4;i++)flux[i] = 0.5 * ((f1[i]+f2[i]) * n - roe[i]);
}

double Euler2d::maximumEigenValue(const double *theMean, const mPoint &point) const
{
  double rho_inv = 1./theMean[0];
  double v2 = (theMean[1]*theMean[1]+theMean[2]*theMean[2])*rho_inv;
  double p = (gamma-1.) * (theMean[3] - 0.5*v2);
  v2 *=rho_inv;
  double c = sqrt(gamma*p*rho_inv);
  return c + sqrt(v2);
}

double Euler2d::maximumEigenValue(mVector &n, const double *theMean, const double* theMean1) const
{
  double rho_inv = 1./theMean[0];
  double p = (gamma-1.) * (theMean[3] - 0.5*(theMean[1]*theMean[1]+theMean[2]*theMean[2])*rho_inv);
  double c = sqrt(gamma*p*rho_inv);
  double lambda1 =  c + fabs((theMean[1]*n(0)+theMean[2]*n(1))*rho_inv);

  rho_inv = 1./theMean1[0];
  p = (gamma-1.) * (theMean1[3] - 0.5*(theMean1[1]*theMean1[1]+theMean1[2]*theMean1[2])*rho_inv);
  c = sqrt(gamma*p*rho_inv);
  double lambda2 =  c + fabs((theMean1[1]*n(0)+theMean1[2]*n(1))*rho_inv);
  if (lambda1>lambda2) return lambda1; else return lambda2;
}

double Euler2d::computeCFL(const int order,const double *theMean, const mPoint &point) const 
{
  double rho_inv = 1./theMean[0];
  double v2 = (theMean[1]*theMean[1]+theMean[2]*theMean[2])*rho_inv;
  double p = (gamma-1.) * (theMean[3] - 0.5*v2);
  v2 *=rho_inv;
  double c = sqrt(gamma*p*rho_inv);
  // printf("eigenvalue %e \n",c + sqrt(v2));
  // max eigenvalue = vel sound + vel fluid
  return 1./((2*p+1)*(c + sqrt(v2)));
}

double Euler2d::maxEigenValueCube(double du,double dv,double dw,double *theMean, const mPoint &point) const
{
  double rho_inv = 1./theMean[0];
  double v2 = (theMean[1]*theMean[1]+theMean[2]*theMean[2])*rho_inv;
  double p = (gamma-1.) * (theMean[3] - 0.5*v2);
  double c = sqrt(gamma*p*rho_inv);
  return (fabs(theMean[1])+fabs(theMean[2]))*rho_inv + c*2.;
}

 void Euler2d::computeJacobian(mVector& n, const double* Q,double** Jac) const 
{
  double nx, ny;
  nx=n(0); ny=n(1); 
  
  // Define constants and variables
  const double invrho=1./Q[0];
  const double u=Q[1]*invrho;
  const double v=Q[2]*invrho;
  const double vn=nx*u+ny*v;
  const double e0=Q[3]*invrho;
  const double ek=0.5*(u*u+v*v);
  const double h0=e0*gamma-(gamma-1.)*ek;
  
  // Define the matrix
  Jac[0][0]=0;
  Jac[0][1]=nx;
  Jac[0][2]=ny;
  Jac[0][3]=0;
  Jac[1][0]=(gamma-1.)*ek*nx-u*vn;
  Jac[1][1]=vn-(gamma-2.)*u*nx;
  Jac[1][2]=u*ny-(gamma-1.)*v*nx;
  Jac[1][3]=(gamma-1.)*nx;
  Jac[2][0]=(gamma-1.)*ek*ny-v*vn;
  Jac[2][1]=v*nx-(gamma-1.)*u*ny;
  Jac[2][2]=vn-(gamma-2.)*v*ny;
  Jac[2][3]=(gamma-1.)*ny;
  Jac[3][0]=((gamma-1.)*ek-h0)*vn;
  Jac[3][1]=h0*nx-(gamma-1.)*u*vn;
  Jac[3][2]=h0*ny-(gamma-1.)*v*vn;
  Jac[3][3]=gamma*vn;
}

void Euler2d::compute_left_eigenvectors(double *U, mVector &n, double *lev) const
{
	const double gm1 = gamma-1.;
	const double rho_inv = 1./U[0];
	const double u = U[1]*rho_inv;
	const double v = U[2]*rho_inv;
	const double e_k = 0.5*(u*u+v*v);
	const double aa = gamma*gm1*(U[3]*rho_inv-e_k);
	const double a = sqrt(aa);
	const double aa_inv = 1./aa;
	const double vn = u*n(0)+v*n(1);
lev[0] = (gm1*e_k+a*vn)*0.5*aa_inv;
lev[1] = (-gm1*u-a*n(0))*0.5*aa_inv;
lev[2] = (-gm1*v-a*n(1))*0.5*aa_inv;
lev[3] =  gm1*0.5*aa_inv;
lev[4] = (aa-gm1*e_k)*aa_inv;
lev[5] = gm1*u*aa_inv;
lev[6] = gm1*v*aa_inv;
lev[7] = -gm1*aa_inv;
lev[8] = (gm1*e_k-a*vn)*0.5*aa_inv;
lev[9] = (-gm1*u+a*n(0))*0.5*aa_inv;
lev[10] = (-gm1*v+a*n(1))*0.5*aa_inv;
lev[11] = gm1*0.5*aa_inv;
lev[12] = v*n(0)-u*n(1);
lev[13] = n(1);
lev[14] = -n(0);
lev[15] = 0.0;
}

void Euler2d::compute_right_eigenvectors(double *U, mVector &n, double *rev) const
{
    const double n0 = n(0);
	const double n1 = n(1);
	const double u = U[1]/U[0];
	const double v = U[2]/U[0];
	const double e_k = 0.5*(u*u+v*v);
	const double a = sqrt(gamma*(gamma-1.)*(U[3]/U[0]-e_k));
	rev[0] = 1.;
rev[1] = 1.;
rev[2] = 1.;
rev[3] = 0.;
rev[4] = u-a*n0;
rev[5] = u;
rev[6] = u+a*n0;
rev[7] = n1;
rev[8] = v-a*n1;
rev[9] = v;
rev[10] = v+a*n1;
rev[11] = -n0;
const double vn = u*n0+v*n1;
const double h0 = a*a/(gamma-1.) + e_k;
rev[12] = h0-a*vn;
rev[13] = e_k;
rev[14] = h0+a*vn;
rev[15] = u*n1-v*n0;

  }

int Euler2d::getEdgeOrientation(mVector &n,const mPoint &position,double *Q) 
const
{
double p = (gamma-1)*Q[3] - 0.5*(gamma-1)*(Q[1]*Q[1]+Q[2]*Q[2])/Q[0];
   double flux[4];
      flux[1] = n[0] * p;
if (flux[1] > 0) return (-1);
else return (1);

}


double Euler3d::maximumEigenValue (const double *theMean,const mPoint &position) const
{
  double v2 = 
  (theMean[1]*theMean[1]+theMean[2]*theMean[2]+theMean[4]*theMean[4])/(theMean[0]*theMean[0]);
  double p = computeP3d(theMean);
  double c = sqrt(gamma*p/theMean[0]);
  // max eigenvalue = vel sound + vel fluid
  return c + sqrt(v2);
}

double Euler3d::maximumEigenValue(mVector &n, const double *theMean, const double* theMean1) const
{
  double rho_inv = 1./theMean[0];
  double p = (gamma-1.) * (theMean[3] - 0.5*(theMean[1]*theMean[1]+theMean[2]*theMean[2]+theMean[4]*theMean[4])*rho_inv);
  double c = sqrt(gamma*p*rho_inv);
  double lambda1 =  c + fabs((theMean[1]*n(0)+theMean[2]*n(1)+ theMean[4]*n(2))*rho_inv);

  rho_inv = 1./theMean1[0];
  p = (gamma-1.) * (theMean1[3] - 0.5*(theMean1[1]*theMean1[1]+theMean1[2]*theMean1[2]+theMean[4]*theMean[4])*rho_inv);
  c = sqrt(gamma*p*rho_inv);
  double lambda2 =  c + fabs((theMean1[1]*n(0)+theMean1[2]*n(1)+ theMean[4]*n(2))*rho_inv);
  if (lambda1>lambda2) return lambda1; else return lambda2;
}

double Euler3d::maxEigenValueCube(double du,double dv,double dw,double *theMean, const mPoint &point) const
{
  double rho_inv = 1./theMean[0];
  double v2 = (theMean[1]*theMean[1]+theMean[2]*theMean[2]+theMean[4]*theMean[4])*rho_inv;
  double p = (gamma-1.) * (theMean[3] - 0.5*v2);
  double c = sqrt(gamma*p*rho_inv);
  return (theMean[1]*du+theMean[2]*dv+theMean[4]*dw)*rho_inv + c*(du+dv+dw);
}

 void Euler3d::computeJacobian(mVector& n, const double* Q,double** Jac) const 
{
  double nx, ny, nz;
  nx=n(0); ny=n(1); nz=n(2);
  
  // Define constants and variables
  const double invrho=1./Q[0];
  const double u=Q[1]*invrho;
  const double v=Q[2]*invrho;
  const double w=Q[4]*invrho;
  const double vn=nx*u+ny*v+nz*w;
  const double e0=Q[3]*invrho;
  const double ek=0.5*(u*u+v*v+w*w);
  const double h0=e0*gamma-(gamma-1.)*ek;
  
  // Define the matrix
  Jac[0][0]=0;
  Jac[0][1]=nx;
  Jac[0][2]=ny;
  Jac[0][3]=0;
  Jac[0][4]=nz;     
  Jac[1][0]=(gamma-1.)*ek*nx-u*vn;
  Jac[1][1]=vn-(gamma-2.)*u*nx;
  Jac[1][2]=u*ny-(gamma-1.)*v*nx;
  Jac[1][3]=(gamma-1.)*nx;
  Jac[1][4]=u*nz-(gamma-1.)*w*nx;
  Jac[2][0]=(gamma-1.)*ek*ny-v*vn;
  Jac[2][1]=v*nx-(gamma-1.)*u*ny;
  Jac[2][2]=vn-(gamma-2.)*v*ny;
  Jac[2][3]=(gamma-1.)*ny;
  Jac[2][4]= v*nz-(gamma-1.)*w*ny;
  Jac[3][0]=((gamma-1.)*ek-h0)*vn;
  Jac[3][1]=h0*nx-(gamma-1.)*u*vn;
  Jac[3][2]=h0*ny-(gamma-1.)*v*vn;
  Jac[3][3]=gamma*vn;
  Jac[3][4]=h0*nz-(gamma-1.)*w*vn;
  Jac[4][0]=(gamma-1.)*ek*nz-w*vn;
  Jac[4][1]=w*nx-(gamma-1.)*u*nz;
  Jac[4][2]=w*ny-(gamma-1.)*v*nz;
  Jac[4][3]=(gamma-1.)*nz;
  Jac[4][4]=vn-(gamma-2.)*w*nz;
}

/********** P O S T   P R O **********/

#include <stdio.h>
int Euler::getNbFieldsOfInterest() const
{
  return 7;
}

void Euler::getNameAndSize(int i, int &n, char *name) const
{
  switch(i)
    {
    case 0:
      n = 3;
      strcpy(name,"velocity");
      break;
    case 1:
      n = 1;
      strcpy(name,"pressure");
      break;
    case 2:
      n = 1;
      strcpy(name,"mach");
      break;
    case 3:
      n = 1;
      strcpy(name,"density");
      break;
    case 5:
      n = 1;
      strcpy(name,"entropy");
      break;
    case 4:
      n = 1;
      strcpy(name,"error");
      break;
    case 6:
      n = 1;
      strcpy(name,"limit");
      break;

    }
}

void Euler2d::getIthFieldOfInterest(int i, const mPoint &p, double *field, 
double *res, double T) const
{
  switch(i)
    {
    case 0:
      {
	double inv_rho = 1./field[0];
	res[0] = field[1]*inv_rho;
	res[1] = field[2]*inv_rho;
	res[2] = 0.0;
	break;
      }
    case 1:
      { 
	double exact[4]; 
	farField->eval(p,T,exact);
	double p_ex = (gamma-1.)*(exact[3] -  
		        0.5*(exact[1]*exact[1]+exact[2]*exact[2])/exact[0]);
	double pfree = p_ex + 0.5*(exact[1]*exact[1]+exact[2]*exact[2])/exact[0];
	double denom = 0.5 * (exact[1]*exact[1] + exact[2]*exact[2])/exact[0];
	double p = (gamma-1.)*(field[3] - 
		        0.5*(field[1]*field[1]+field[2]*field[2])/field[0]);
	double pt = p + 0.5*(field[1]*field[1]+field[2]*field[2])/1.4; 
	double p_total=p*pow(1.+0.5*(gamma-1.)*(field[1]*field[1]+field[2]*field[2])/(field[0]*gamma*p),gamma/(gamma-1.));
	res[0]=(p-p_ex)/(0.5*(exact[1]*exact[1]+exact[2]*exact[2])/exact[0]); 
	res[0]=p;
//	res[0] = (p_ex-p)/denom; //pressure coefficient
    break;
      }
    case 2:
      {
	double inv_rho = 1./field[0];
	double vsq = (field[1]*field[1]+field[2]*field[2])*inv_rho;
	double p = (gamma-1) * (field[3] - 0.5*vsq);
	double csq =  gamma*p*inv_rho;
	res[0] = sqrt(vsq*inv_rho/csq);
	//printf("%f\n",sqrt(vsq/field[0]));
      }
      break;
    case 3:
      res[0] = field[0];
      break;
    case 4:
      double exact[4];
      farField->eval(p,T,exact);
      res[0] = (field[1]/field[0]-exact[1]/exact[0]); 
	  //res[0]= (field[1]-field[2])/field[0] - (exact[1]-exact[2])/exact[0];
	  res[0] = field[1]-exact[1];
      break;
    case 5:
      {
	double pres = (gamma-1.)*(field[3] - 
		        0.5*(field[1]*field[1]+field[2]*field[2])/field[0]); 
	double exact[4]; 
    farField->eval(p,T,exact); 
	double p_ex = (gamma-1.)*(exact[3] -  
		        0.5*(exact[1]*exact[1]+exact[2]*exact[2])/exact[0]);
	res[0]= pres/p_ex/pow(field[0]/exact[0],gamma)-1.;
      }
      break;
    case 6:
      res[0] = field[0];
      break;
    }
}


void Euler3d::getIthFieldOfInterest(int i, const mPoint &, double *field, 
double *res, double time) const
{
  switch(i)
    {
    case 0:
      {
	double inv_rho = 1./field[0];
	res[0] = field[1]*inv_rho;
	res[1] = field[2]*inv_rho;
	res[2] = field[4]*inv_rho;
      }
      break;
    case 1:
      res[0] =  computeP3d(field);
      break;
    case 2:
      {
		  double inv_rho = 1./field[0];
	double vsq = ((field[1]*field[1]+field[2]*field[2]+field[4]*field[4])*inv_rho);
	//double p = (gamma-1) * (field[3] - 0.5*vsq);
	double p = computeP3d(field);
	double csq =  gamma*p*inv_rho;
	res[0] = sqrt(vsq*inv_rho/csq);
	  }
      break;
    case 3:
      res[0] = field[0];
      break;
    case 4:
      res[0] = field[0];
      break;
    case 5:
      res[0] = field[0];
      break;
	case 6:
		res[0] = field[0];
      break;

    }
}

void Euler2d::RHS(double *field, const mPoint &position, double time,double *rhs) const
{
  rhs[0] = 0.0;
  rhs[1] = 0.0;
  rhs[2] = 0.0;
  rhs[3] = 0.0;
}
 
bool Euler2d::isPhysical ( double *Q ) const
{
  //check for negative density and pressure
  if(Q[0] < 0)return false;
  if((gamma-1.)*(Q[3] - 0.5*(Q[1]*Q[1]+Q[2]*Q[2])/Q[0])< 0.0)return false;
  return true;
}

void Euler3d::RHS (double *field, const mPoint &position, double time, double *rhs) const
{
  rhs[0] = 0.0;
  rhs[1] = 0.0;
  rhs[2] = 0.0;
  rhs[3] = 0.0;
  rhs[4] = 0.0;
}

bool Euler3d::isPhysical ( double *Q ) const
{
  if(Q[0] < 0) return false;
  double pres = computeP3d(Q);
  if(pres < 0) return false;
  return true;
}

/************************************************************************************/

void Euler2dExactRiemann::riemannSolver( mVector &n ,
				const mPoint &position,
				double *u0,
				double *u1,
				double *flux) const
{

  const double n0 = n[0] ,     n1 = n[1];
  const double ul = u0[1],     ur = u1[1];  //moments, not velocity!
  const double vl = u0[2],     vr = u1[2];
  const double rhol = u0[0],   rhor = u1[0];
  const double rhol_inv= 1./ rhol; 
  const double rhor_inv= 1./ rhor; 
  const double ulft = (n0 * ul + n1 * vl ) * rhol_inv;
  const double vlft = (n1 * ul - n0 * vl ) * rhol_inv;
  const double urght = (n0 * ur + n1* vr ) * rhor_inv;
  const double vrght = (n1 * ur - n0 * vr ) * rhor_inv;
  const double gm1 = gamma-1.;
  double plft = gm1*(u0[3] - 0.5*(ul*ul+vl*vl)*rhol_inv);
  double prght = gm1*(u1[3] - 0.5*(ur*ur+vr*vr)*rhor_inv);
  double vnav,vtav,pav,rhoav;

/*  const double rhol =1.0;
  const double rhol_inv = 1./1.;
  const double ulft =0.;
  const double vlft=0;
  const double plft=.01;
  const double rhor=1.; 
  const double rhor_inv =1./1.;
  const double  urght=0.;const double vrght=0;
  const double prght=100.;*/
  
  double clft,crght,tmp1;
  tmp1 = gamma*plft*rhol_inv;                 
  if(tmp1 < 0.0) 
  {
    printf("vel of sound negative...\n"); 
    printf("rho = %f p = %f at point (%f,%f) \n",rhol, plft,position(0),position(1));
    tmp1 = -tmp1;
	plft = -plft;
//	nonPhysical = 1;
    // exit(0);
  }
  clft=sqrt(tmp1);     
  tmp1 = gamma*prght*rhor_inv;
  if(tmp1 < 0.0)
  {
    printf("vel of sound negative...\n");
    printf("rho = %f p = %f at point (%f,%f)\n",rhor, prght,position(0),position(1));
    tmp1 = -tmp1;
	prght = -prght;
//	if (nonPhysical==1) nonPhysical = 3; else nonPhysical = 2;
    //exit(0);
  }
  crght=sqrt(tmp1);
 
  if (ulft>=clft && urght>=crght) //completely supersonic in left state
    {
      rhoav = rhol;
      vnav  = ulft;
      vtav  = vlft;
      pav   = plft;
    }
  else if (urght<=-crght && ulft<=-clft) //comletely supersonic in right state
    {
      rhoav = rhor;
      vnav  = urght;
      vtav  = vrght;
      pav   = prght;
    }
  else {
  int i;
  double tmp2,scrch1;
  double pstar, ustar, rhostar, cstar;
  double rhos, us, ps,cs;
  double prevpstar, relpstar;
  double flft,dflft,frght,dfrght;
 
  const double invgamma = 1./gamma;
  const double gp1  = gamma+1.;
  const double C1       = 0.5*(gp1)*invgamma;
  const double C3       = 0.5*gm1*invgamma; 
  const double C4       = 2./gm1;
  
  // Initial guess for pstar by two shock approximation; Toro p.128
  double udiff = urght - ulft;
  double p0 = 0.5*(plft+prght) -0.125*udiff*(rhol+rhor)*(clft+crght);
  p0 = max2(smallp,p0);
  double gl = sqrt(2./(rhol*(gm1*plft+gp1*p0)));
  double gr = sqrt(2./(rhor*(gm1*prght+gp1*p0)));
  pstar = (gl*plft+gr*prght-udiff)/(gl+gr);
  
  relpstar=RIEM_TOL;
  pstar=max2(RIEM_TOL,pstar);
  prevpstar=pstar;                           /* save intial value */
  
  for (i=0; i<RIEM_ITER && relpstar>=RIEM_TOL ; i++)
    {
      if (pstar>plft)
	{
	  tmp1 = 1./(gp1*pstar+gm1*plft);
	  tmp2 = sqrt(2.*tmp1*rhol_inv);
	  flft = (pstar-plft)*tmp2;
	  dflft = tmp2*(1.-0.5*gp1*(pstar-plft)*tmp1);
	}
      else
	{
	  tmp2 = pow(pstar/plft,C3);
	  flft = C4*clft*(tmp2-1.);
	  dflft = tmp2*clft*invgamma/pstar;
	}
      if (pstar>prght)	
	{
	  tmp1 = 1./(gp1*pstar+gm1*prght);
	  tmp2 = sqrt(2.*tmp1*rhor_inv);
	  frght = (pstar-prght)*tmp2;
	  dfrght = tmp2*(1.-0.5*gp1*(pstar-prght)*tmp1); 	  
	}
      else
	{
	  tmp2 = pow(pstar/prght,C3);
	  frght = C4*crght*(tmp2-1.);
	  dfrght = tmp2*crght*invgamma/pstar;
	}
      pstar -= (flft+frght+udiff)/(dflft+dfrght);
      pstar=max2(smallp,pstar);
      
      relpstar=fabs(pstar-prevpstar)/(.5*(prevpstar+pstar));// compute relative error
      prevpstar=pstar; 
    }
  
  if (relpstar>=RIEM_TOL)
    {
      cout << "Riemann solver failed to converge in " << RIEM_ITER << "iterations\n";
      cout << "Relative error is " << relpstar << "\n";
    }
  //Recompute F_l  and F_r with the final pstar to get ustar
  if (pstar>plft)
    flft = (pstar-plft)*sqrt(2.*rhol_inv/(gp1*pstar+gm1*plft));
  else
    flft = C4*clft*(pow(pstar/plft,C3)-1.);
  if (pstar>prght)	
    frght = (pstar-prght)*sqrt(2.*rhor_inv/(gp1*pstar+gm1*prght));
  else
    frght = C4*crght*(pow(pstar/prght,C3)-1.);
  
  ustar = (ulft+urght+frght-flft)*0.5;
  
  if (ustar>=0) //choose left state
    {
      ps=plft;
      us=ulft;  vtav= vlft;  
      rhos=rhol;
      cs=clft;
      scrch1 = 1.;
    }
  else         //choose right state
    {
      ps=prght;
      us=urght; vtav= vrght;  
      rhos=rhor;
      cs=crght;
      scrch1=-1.;
    }
  
  if ( pstar >= ps  ) //shock, choose star state
    {
      double C= gm1/gp1;
      tmp1 = pstar/ps;
      rhoav = rhos*(tmp1+C)/(C*tmp1+1.);
      vnav   = ustar;
      pav   = pstar;
    }// rarefaction, choose star or rarefaction state
  else //rarefaction
    { 
      rhostar = rhos*pow(pstar/ps,invgamma);
      cstar = sqrt(gamma*pstar/rhostar);
      
      if ((ustar-scrch1*cstar)*ustar<=0.0) 
	//in brackets speed of the rarefaction tail 
	{//choose star state
	  rhoav = rhostar;
	  vnav   = ustar;
	  pav   = pstar;
	}	
      else
	{//choose a state inside rarefaction wave
	  double C = 2./gp1;
	  tmp1 =C+scrch1*gm1*us/(cs*gp1); 
	  rhoav = rhos*pow(tmp1,C4);
	  vnav   = C*(gm1*0.5*us+scrch1*cs); 
	  pav   = ps*pow(tmp1,C4*gamma );
	}
    }
  if ((clft+crght)*2.<=gm1*(urght-ulft)) {
    printf("VACUUM\n");
    //exit(0);
  }
  }
 
  flux[0] = rhoav   * vnav;
  flux[1] = flux[0] * (n0 * vnav + n1 * vtav) + pav * n0;
  flux[2] = flux[0] * (n1 * vnav - n0 * vtav) + pav * n1;
  flux[3] = vnav    * (pav/gm1 + 0.5*rhoav*(vnav*vnav+vtav*vtav)+pav);
}

/**********************************************************/

#define SMALL_DENSITY  1.0e-12
#define SMALL_PRESS ((double)1e-14)
#define SQUARE(a) a*a

/*  -----------  prim_vanLeerFlux()----... use vanLeer 1982 ref, with the form
                                           Anderson et.al, AIAA Jol, V.24 No.9 
                                           Sep 1986
				        ...This version is coded assuming that
					   uL and uR are passed in already in
					   primative variables, and returns
					   the flux F in conservative variables
	       			   ---------------------------------- */
void EulervanLeer::riemannSolver( mVector &n , 
				const mPoint &position,
				double *u0, 
				double *u1, 
				double *flux) const
{
  double rhol   = u0[0];
  double rhor   = u1[0];     
  double invrhol, unl, utan0l, utan1l;
  double invrhor,unr, utan0r, utan1r;
  double pr,pl,cl,cr,odenom;
  double fL[4], fR[4];
  double fplus, fminus;
  const  double ogamma = 1./gamma;
  const double gm1 = gamma-1.;
   
  if (rhol < SMALL_DENSITY) 
    {
      printf("resetting density in vl from uL[RHO]=%g\n",u0[0]);
      rhol = rhor;
    }
  if (rhor < SMALL_DENSITY) 
    {
      printf("resetting density in vl from uR[RHO]=%g\n",u1[0]);
      rhor = rhol;
    }

  invrhol  = 1./rhol; invrhor  = 1./rhor;
  double vnl = (n(0) * u0[1] + n(1) * u0[2] ) * invrhol;
  double vtl = (n(1) * u0[1] - n(0) * u0[2] ) * invrhol;

  double vnr = (n(0) * u1[1] + n(1) * u1[2] ) * invrhor;
  double vtr = (n(1) * u1[1] - n(0) * u1[2] ) * invrhor;
  printf("vnl %e vnr%e\n",vnl,vnr);
  unl    = vnl;
  utan0l = vtl;
  utan1l = 0.0;//for now
     
  unr    = vnr;
  utan0r = vtr;
  utan1r = 0.0;

  pl= gm1*(u0[3] - 0.5*(u0[1]*u0[1]+u0[2]*u0[2])*invrhol);
  pr = gm1*(u1[3] - 0.5*(u1[1]*u1[1]+u1[2]*u1[2])*invrhor);
  pl  = max2(SMALL_PRESS,pl);
  pr  = max2(SMALL_PRESS,pr);/*..safe L/R pressures and sound speeds */
 
  cl = sqrt(gamma*pl*invrhol);
  cr = sqrt(gamma*pr*invrhor);

  odenom = 0.5/(gamma*gamma -1.);
 
  if ( unl >= cl ){           /* ...completely SuperSonic (in) left state */
    mVector f[4];
    Fi(position,u0,f);    
    fL[0] = f[0]*n;
    fL[1] = f[1]*n;
    fL[2] = f[2]*n;
    fL[3] = f[3]*n;
  }
  else if (unl > -cl){       /*  ...subsonic left state (either in or out) */
    const double f1l    = (gamma - n(0))*unl + 2.*cl*n(0);
    fplus       = rhol*0.25 * SQUARE(unl + cl) / cl;     
    fL[0] = fplus;
    fL[1] = fplus * f1l * ogamma;
    fL[2] = fplus * (utan0l);
    fL[3] = fplus *((utan0l*utan0l + utan1l*utan1l)*0.5 +
                                            f1l * f1l * odenom );
  }
  else{
    fL[0] = fL[1] = fL[2] = fL[3] = 0.0; //removed fL[4] 
  }
  
  if (unr <= -cr){           /* ...completely Supersonic (in) right state */
    mVector f[4];
    Fi(position,u1,f);    
    fR[0] = f[0]*n;
    fR[1] = f[1]*n;
    fR[2] = f[2]*n;
    fR[3] = f[3]*n;    
  }
  else if (unr < cr){       /* ...subsonic right state (either in or out) */
    const double f1r    = (gamma - n(0))*unr - 2.*cr*n(0);
    fminus =-rhor*0.25 * SQUARE(unr - cr) / cr;
    fR[0] = fminus;
    fR[1] = fminus * f1r * ogamma*n(0);
    fR[2] = fminus * (utan0r);
    fR[3] = fminus *((utan0r*utan0r + utan1r*utan1r)*0.5 +
                                                 f1r*f1r*odenom );
  }
  else{
    fR[0] = fR[1] = fR[2] = fR[3] = 0.0;// removed fR[4]= 
  }

             /* ...split fluxes contribute to make final Flux SEE NOTE.1 */
  flux[0] = (fL[0] + fR[0]);
  flux[1] = (fL[1] + fR[1]);
  flux[2] = (fL[2] + fR[2]);
  flux[3] = (fL[3] + fR[3]);
  return;
}
/* --------------------------------------------------------------------------
   NOTES:
   o  This routine has been converted to primative variables, so
      when uL and uR are indexed by thisMom and mom0/1 its referring to
      speeds (u,v,w) not momentum. the flux vector F IS indexed and returned
      in conservative variables, so its proper to index it by momentum.
   o  vanLeer splitting is on normal mach No...notice that there is some
      ambiguity if this is the normal mach no unl/cl or unr/cr ( or some
      middle state etc..). So here we're accumulating it in since
      van leer's transistion is smooth, just incase unl/cl would put us
      in a different case than unr/cr.
                                                                 -MJA/B
   -------------------------------------------------------------------------- */

void Euler2dHLLC::riemannSolver( mVector &n ,
                                const mPoint &position,
                                double *Ul,
                                double *Ur,
                                double *flux) const
{
 double fluxr[4];    
 double pl,pr,cl,cr,vnr,vnl,vtl,vtr,rhor,rhol;
 double ql, qr; 
 rhol=Ul[0];
 rhor=Ur[0];
 double rhoav = 0.5 * (rhol+rhor);
 
 const double n0=n[0];  
 const double n1=n[1];         
 
 //normal and tangential components of velocities
 vnl=(Ul[1]*n0+Ul[2]*n1)/rhol;
 vnr=(Ur[1]*n0+Ur[2]*n1)/rhor;
 vtl=(Ul[1]*n1-Ul[2]*n0)/rhol;
 vtr=(Ur[1]*n1-Ur[2]*n0)/rhor;

 pl=computeP2d(Ul,gamma);
 pr=computeP2d(Ur,gamma);
 cl=sqrt(gamma*pl/rhol);
 cr=sqrt(gamma*pr/rhor);
 double cav = 0.5 * (cl+cr);

 double pstar = 0.5*(pl+pr-(vnr-vnl)*rhoav*cav);
 double ustar = 0.5*(vnl+vnr-(pr-pl)/(rhoav*cav));
 if (pstar>pl) ql = sqrt(1.+0.5*(gamma+1.)/gamma*(pstar/pl-1.));
 else ql = 1.0;
 if (pstar>pr) qr = sqrt(1.+0.5*(gamma+1.)/gamma*(pstar/pr-1.));
 else qr = 1.0;
 double sl = vnl-cl*ql;
 double sr = vnr-cr*qr;
 
 if (sl>0.) 
   {
     flux[0] = rhol   * vnl;
     flux[1] = vnl*Ul[1] + pl * n0;
     flux[2] = vnl*Ul[2] + pl * n1;
     flux[3] = vnl*(pl + Ul[3]);
     return;
   }
 else 
   if (sr<0.0)
     {
       flux[0] = rhor   * vnr;
       flux[1] = vnr*Ur[1] + pr * n0;
       flux[2] = vnr*Ur[2] + pr * n1;
       flux[3] = vnr*(pr + Ur[3]);
       return;
     }
   else
     if (sl<0.0 && 0.0<=ustar) 
       {
	 double fact = rhol*(sl-vnl)/(sl-ustar);
	 flux[0] = rhol* vnl +sl*(fact-rhol);
	 flux[1] = vnl*Ul[1] + pl * n0+sl*(vnl*fact*n0+vtl*n1*rhol-Ul[1]);
	 flux[2] = vnl*Ul[2] + pl * n1+sl*(vnl*fact*n1-vtl*n0*rhol-Ul[2]);
	 flux[3] = vnl*(pl + Ul[3])+sl*(fact*Ul[3]/rhol+(ustar-vnl)*(ustar+pl/(rhol*(sl-vnl))));
	 return;
       }
     else
       if (ustar<0.0 && 0.0<=sr) 
	 {
	   double fact = rhor*(sr-vnr)/(sr-ustar);
	   flux[0] = rhol* vnl +sr*(fact-rhor);
	   flux[1] = vnr*Ur[1] + pr * n0+sr*(vnr*fact*n0+vtr*n1*rhor-Ur[1]);
	   flux[2] = vnr*Ur[2] + pr * n1+sr*(vnr*fact*n1-vtr*n0*rhor-Ur[2]);
	   flux[3] = vnr*(pr + Ur[3])+sr*(fact*Ur[3]/rhor+(ustar-vnr)*(ustar+pr/(rhor*(sr-vnr))));
	   return;
	 }
 
}


/* Got this from Peng*/
/*

void Euler2dHLLC::riemannSolver( mVector &n ,
                                const mPoint &position,
                                double *Ul,
                                double *Ur,
                                double *fluxl) const
{
 double fluxr[4];    
 double pl,pr,al,ar,ur,ul,vl,vr,rhor,rhol,Hl,Hr,Sl, Sr, Sstar;

rhol=Ul[0];
rhor=Ur[0];

double tan,cos,sin;

cos=n[0];  
sin=n[1];         


double project=cos*n[0]+sin*n[1];

//normal component of velocities
ul=(Ul[1]*cos+Ul[2]*sin)/rhol;
ur=(Ur[1]*cos+Ur[2]*sin)/rhor;
//tangential component of velocities
vl=(-Ul[1]*sin+Ul[2]*cos)/rhol;
vr=(-Ur[1]*sin+Ur[2]*cos)/rhor;
pl=computeP2d(Ul,gamma);
pr=computeP2d(Ur,gamma);
ar=sqrt(gamma*pr/Ur[0]);
al=sqrt(gamma*pl/Ul[0]);
Hl=(Ul[3]+pl)/rhol;   
Hr=(Ur[3]+pr)/rhor;

double utilda, atilda, Htilda;
utilda=(sqrt(rhol)*ul+sqrt(rhor)*ur)/(sqrt(rhol)+sqrt(rhor));
Htilda=(sqrt(rhol)*Hl+sqrt(rhor)*Hr)/(sqrt(rhol)+sqrt(rhor));
atilda=sqrt((gamma-1)*(Htilda-0.5*utilda*utilda));
Sl=utilda-atilda;
Sr=utilda+atilda;

Sstar=pr-pl+rhol*ul*(Sl-ul)-rhor*ur*(Sr-ur);
Sstar=Sstar/(rhol*(Sl-ul)-rhor*(Sr-ur));

double factl=(rhol*(Sl-ul))/(Sl-Sstar);
double factr=(rhor*(Sr-ur))/(Sr-Sstar);

double Ustarl[4],Ustarr[4];
Ustarl[0]=factl;
Ustarl[1]=factl*Sstar;  
Ustarl[2]=factl*vl;
Ustarl[3]=factl*(Ul[3]/rhol+(Sstar-ul)*(Sstar+pl/(rhol*(Sl-ul))));
Ustarr[0]=factr;      
Ustarr[1]=factr*Sstar;
Ustarr[2]=factr*vr;
Ustarr[3]=factr*(Ur[3]/rhor+(Sstar-ur)*(Sstar+pr/(rhor*(Sr-ur))));

double UstarGl[4],UstarGr[4];
UstarGl[0]=Ustarl[0];
UstarGl[1]=Ustarl[1]*cos-Ustarl[2]*sin;
UstarGl[2]=Ustarl[1]*sin+Ustarl[2]*cos;
UstarGl[3]=Ustarl[3];
UstarGr[0]=Ustarr[0];
UstarGr[1]=Ustarr[1]*cos-Ustarr[2]*sin; 
UstarGr[2]=Ustarr[1]*sin+Ustarr[2]*cos;
UstarGr[3]=Ustarr[3];

mVector ftmp1[5],ftmp2[5];

if (0.<=Sl)
{
 Fi(position,Ul,ftmp1);
 fluxl[0] = fluxr[0] = ftmp1[0] * n;
 fluxl[1] = fluxr[1] = ftmp1[1] * n;
 fluxl[2] = fluxr[2] = ftmp1[2] * n;
 fluxl[3] = fluxr[3] = ftmp1[3] * n;
}
else if (0.>=Sr)
{
 Fi(position,Ur,ftmp2);
 fluxl[0] = fluxr[0] = ftmp2[0] * n;   
 fluxl[1] = fluxr[1] = ftmp2[1] * n;   
 fluxl[2] = fluxr[2] = ftmp2[2] * n;
 fluxl[3] = fluxr[3] = ftmp2[3] * n;
}
else if(0.<=Sstar && 0.>=Sl)
{
Fi(position,Ul,ftmp1);
fluxl[0] = fluxr[0] = ftmp1[0]*n + Sl*(UstarGl[0]-Ul[0])*project;
fluxl[1] = fluxr[1] = ftmp1[1]*n + Sl*(UstarGl[1]-Ul[1])*project;
fluxl[2] = fluxr[2] = ftmp1[2]*n + Sl*(UstarGl[2]-Ul[2])*project;
fluxl[3] = fluxr[3] = ftmp1[3]*n + Sl*(UstarGl[3]-Ul[3])*project;
}
else if((0.<=Sr) && (0.>=Sstar))
{
 Fi(position,Ur,ftmp2);
fluxl[0] = fluxr[0] = ftmp2[0]*n + Sr*(UstarGr[0]-Ur[0])*project;
fluxl[1] = fluxr[1] = ftmp2[1]*n + Sr*(UstarGr[1]-Ur[1])*project;
fluxl[2] = fluxr[2] = ftmp2[2]*n + Sr*(UstarGr[2]-Ur[2])*project;
fluxl[3] = fluxr[3] = ftmp2[3]*n + Sr*(UstarGr[3]-Ur[3])*project;
}
}
*/

void Euler3dRoe::riemannSolver (mVector &n, 
				const mPoint &position,
				double *uLx, 
				double *uRx,
				double *flux) const
{
 int    l,k;
 const int NEQ = 5;
 double rho,u,v,w,halfrhoun,H,HL,HR,sqr_rhoL,sqr_rhoR,invz1,
        drho,du,dv,dw,dp,un,u2,c,c2,M2,rho2c,b[3],P[NEQ*NEQ],
        lamb[NEQ],dW[NEQ],dflux[NEQ],nx,ny,nz;
 double Gamma = gamma;
 double GM1 = gamma-1.0;
 double uL[5],uR[5];

 uL[0] = uLx[0]          ; uR[0] = uRx[0];
 uL[1] = uLx[1]/uL[0]    ; uR[1] = uRx[1]/uR[0];
 uL[2] = uLx[2]/uL[0]    ; uR[2] = uRx[2]/uR[0];
 uL[3] = uLx[4]/uL[0]    ; uR[3] = uRx[4]/uR[0];
 uL[4] = computeP3d (uLx); uR[4] = computeP3d (uRx);

 /* --- left contributions ---*/
 halfrhoun = 0.5*(uL[0]*(uL[1]*n(0)+uL[2]*n(1)+uL[3]*n(2)));
 H        = Gamma/GM1* uL[4]/uL[0]+0.5*(uL[1]*uL[1]+uL[2]*uL[2]+uL[3]*uL[3]);

 flux[0] = halfrhoun;
 flux[1] = halfrhoun*uL[1] + .5*uL[4]*n(0);
 flux[2] = halfrhoun*uL[2] + .5*uL[4]*n(1);
 flux[3] = halfrhoun*uL[3] + .5*uL[4]*n(2);
 flux[4] = halfrhoun*H; 

 /* --- right contributions ---*/
 halfrhoun = 0.5*(uR[0]*(uR[1]*n(0)+uR[2]*n(1)+uR[3]*n(2)));
 H        = Gamma/GM1* uR[4]/uR[0]+0.5*(uR[1]*uR[1]+uR[2]*uR[2]+uR[3]*uR[3]);

 flux[0] += halfrhoun;
 flux[1] += halfrhoun*uR[1] + .5*uR[4]*n(0);
 flux[2] += halfrhoun*uR[2] + .5*uR[4]*n(1);
 flux[3] += halfrhoun*uR[3] + .5*uR[4]*n(2);
 flux[4] += halfrhoun*H; 

 /* --- add dissipation ---*/
 nx=n(0); ny=n(1); nz=n(2);
 
 sqr_rhoL = sqrt(uL[0]);					
 sqr_rhoR = sqrt(uR[0]);					
 
 invz1  = 1./ (sqr_rhoL + sqr_rhoR);				  
 
 HL = Gamma/GM1 * uL[4]/uL[0] + 0.5 *				  
   (uL[1]*uL[1]+uL[2]*uL[2]+uL[3]*uL[3]);		  
 
 HR = Gamma/GM1 * uR[4]/uR[0] + 0.5 *				  
   (uR[1]*uR[1]+uR[2]*uR[2]+uR[3]*uR[3]);		  
 
 
 rho  =   sqr_rhoL * sqr_rhoR;					  
 u    = ( sqr_rhoL* uL[1] + sqr_rhoR * uR[1] ) * invz1;	  
 v    = ( sqr_rhoL* uL[2] + sqr_rhoR * uR[2] ) * invz1;	  
 w    = ( sqr_rhoL* uL[3] + sqr_rhoR * uR[3] ) * invz1;	  
 H    = ( sqr_rhoL* HL    + sqr_rhoR * HR    ) * invz1;	  
 
 drho   = uR[0] - uL[0];					  
 du     = uR[1] - uL[1];					  
 dv     = uR[2] - uL[2];					  
 dw     = uR[3] - uL[3];					  
 dp     = uR[4] - uL[4];					  
 
 un = u*nx + v*ny + w*nz;					  
 u2 = u*u + v*v + w*w;						  
 c2 = GM1 * ( H - 0.5 * u2);					  
 if(c2 >= 0.0)
   {
     c  = sqrt(c2);							  
     M2 = u2/c2;							  
     rho2c = rho * 0.5 / c;						  
     
     b[0] = 0.5 * u2 * nx + rho * ( v*nz - w*ny );			  
     b[1] = 0.5 * u2 * ny + rho * ( w*nx - u*nz );			  
     b[2] = 0.5 * u2 * nz + rho * ( u*ny - v*nx );			  
     
     P[ 0 * NEQ + 0] = nx;						  
     P[ 0 * NEQ + 1] = u * nx;					  
     P[ 0 * NEQ + 2] = v * nx + rho * nz;				  
     P[ 0 * NEQ + 3] = w * nx - rho * ny;				  
     P[ 0 * NEQ + 4] = b[0] ;					  
     
     P[ 1 * NEQ + 0] = ny;						  
     P[ 1 * NEQ + 1] = u * ny - rho * nz;				  
     P[ 1 * NEQ + 2] = v * ny ;					  
     P[ 1 * NEQ + 3] = w * ny + rho * nx;				  
     P[ 1 * NEQ + 4] = b[1] ;					  
     
     P[ 2 * NEQ + 0] = nz;						  
     P[ 2 * NEQ + 1] = u * nz + rho * ny;				  
     P[ 2 * NEQ + 2] = v * nz - rho * nx;				  
     P[ 2 * NEQ + 3] = w * nz ;					  
     P[ 2 * NEQ + 4] = b[2] ;					  
     
     P[ 3 * NEQ + 0] = rho2c;					  
     P[ 3 * NEQ + 1] = rho2c * (u + nx * c);				  
     P[ 3 * NEQ + 2] = rho2c * (v + ny * c);				  
     P[ 3 * NEQ + 3] = rho2c * (w + nz * c);				  
     P[ 3 * NEQ + 4] = rho2c * (H + c * un) ;			  
     
     P[ 4 * NEQ + 0] = rho2c;					  
     P[ 4 * NEQ + 1] = rho2c * (u - nx * c);				  
     P[ 4 * NEQ + 2] = rho2c * (v - ny * c);				  
     P[ 4 * NEQ + 3] = rho2c * (w - nz * c);				  
     P[ 4 * NEQ + 4] = rho2c * (H - c * un) ;			  
     
     lamb[0] = fabs(un) ;						  
     lamb[1] = fabs(un) ;						  
     lamb[2] = fabs(un) ;						  
     lamb[3] = fabs((un + c)) ;					  
     lamb[4] = fabs((un - c)) ;					  
     
     dW[0]= nx * drho           + nz * dv - ny * dw - nx/c2 * dp;	  
     dW[1]= ny * drho - nz * du           + nx * dw - ny/c2 * dp;	  
     dW[2]= nz * drho + ny * du - nx * dv           - nz/c2 * dp;	  
     dW[3]=             nx * du + ny * dv + nz * dw + 1./(rho*c) * dp; 
     dW[4]=           - nx * du - ny * dv - nz * dw + 1./(rho*c) * dp; 
     for (k=0;k<NEQ;k++) { dflux[k] = 0.; }				  
     
     for (k=0;k<NEQ;k++) {						  
       for (l=0;l<NEQ;l++) {						  
	 dflux[k] += lamb[l] * dW[l] * P[l*NEQ+k] ;			  
       }
     }						
			  
     /* --- subtract dissipation from flux ---*/
     for (k=0;k<NEQ;k++) { flux[k] -= dflux[k] * 0.5 ; }
   }
 // different conventions
 double temp = flux[3];flux[3]=flux[4];flux[4]=temp;
}












