#ifndef _TIMEINTEGRATOR_H_ 
#define _TIMEINTEGRATOR_H_
#include <stdio.h>
#include "mMesh.h"
class mDGMesh;
class DGLimiter;

class TimeIntegrator
{
  public:
  int n;// dimension of the mesh
  const int cSize;//number of equations
  int dof; //number of degrees of freedom
  double *soln,*soln_begin;
  double *solnsum,*solnsum_begin;
  mDGMesh *theMesh;
  DGLimiter *theLimiter;
  void assembleVolume(double t);
  void assembleBoundary(double t);
  void assembleVolume(const mMesh::iter begin,const mMesh::iter end,double t);
  void assembleBoundary(const mMesh::iter begin,const mMesh::iter end,double t);
  void limit(double time);
  void limit(const mMesh::iter begin, const mMesh::iter end, double time);
  void computeFullJump();
  void computeL2ProjInCutCells();
  TimeIntegrator(mDGMesh *,DGLimiter *, int cSize,int dof);
  virtual ~TimeIntegrator(){};
   virtual double advanceInTime(double t, double dt) = 0;
};

class ForwardEuler : public TimeIntegrator
{
 public:
  ForwardEuler(mDGMesh *m,DGLimiter *l,int NbEqn, int dof);
  virtual double advanceInTime(double t, double dt);
};

class RungeKuttaTVD2 : public TimeIntegrator
{
 public:
  RungeKuttaTVD2(mDGMesh *m,DGLimiter *l,int NbEqn, int dof);
  virtual double advanceInTime(double t, double dt);
};

class RungeKuttaTVD2adaptive : public TimeIntegrator
{
 protected:
  double* binptr[10];
  double* binptrLI[10];
  double* binptrSI[10];
  bool isFlat [10];
 public:
  RungeKuttaTVD2adaptive(mDGMesh *m,DGLimiter *l,int NbEqn, int dof);
  virtual ~RungeKuttaTVD2adaptive();
  virtual double advanceInTime(double t, double dt);
  
  double firstStage(int from, int to, double t, double DT);
  double secondStage(int from,int to, double t, double DT);
  void LIatHalfStep(int b, double DT);
  void LIatFirstStage(int b);
  void finalResetLI(int b);
  double createStaircase(int from, int to, double t, double DT);
  double lastCellsUp(int NbBins, double t, double DT);
};

class RungeKuttaTVD3 : public TimeIntegrator
{
 public:
  RungeKuttaTVD3(mDGMesh *m,DGLimiter *l,int NbEqn, int dof);
  virtual double advanceInTime(double t, double dt);
};

class RungeKutta4 : public TimeIntegrator
{
 public:
  RungeKutta4(mDGMesh *m,DGLimiter *l,int NbEqn, int dof);
  virtual double advanceInTime(double t, double dt);
};

class Multigrid : public TimeIntegrator
{
 protected:
  RungeKuttaTVD2* rk2;
  RungeKuttaTVD3* rk3;
  RungeKutta4* rk4;
  double** DMatrix;
  double** Jac;
 public:
  Multigrid(mDGMesh *m,DGLimiter *l,int NbEqn, int dof);
  virtual double advanceInTime(double t, double dt);
};
#endif
