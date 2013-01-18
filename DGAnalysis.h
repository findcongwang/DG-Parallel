#ifndef _DGANALYSIS_H_
#define _DGANALYSIS_H_
/*
This class is on the same level than FEAnalysis.
We want to implement something FAST for HP-ADAPTIVITY
We use the AOMD Mesh database.

Hypothesis :
  - Explicit time stepping this hypothesis
  is taken first for simplicity. we'll be
  able to go implicit without major modifications.
  We use Runge and Kutta.
  - Conservation laws : dt(u) + div (F(u)) = rhs
  this can be extended to viscous problems

What should be general :
  - use of mixed meshes, curved elements, non conforming
  - same code for 1-d, 2-d and 3-d !
  - use of any set of shape functions
  - use of any conservation law
*/
#include <iostream>
#include <fstream>
#include <list>
#include <set>
#include "Integrator.h"
class mPoint;
class mDGMesh;
class mEntity;
class ConservationLaw;
class DGCell;
class TimeIntegrator;
class DGBoundaryCell;
class DGLimiter;
class DGSensor;
class GaussIntegrator;
class Geometry;
class PostProGmshView;

using namespace std;
typedef enum TIME_STEPPING {_regular_,_local_,_localAccurate_} _TIME_STEPPING_;

template <class T> class mGeomSearch;

class lessThan
{
public:
	bool operator()(DGCell *one, DGCell *two)
	{
	if (one<two) return true; else return false;
	}
};

class DGAnalysis 
{
  unsigned long DOF;
  list<DGSensor*> allSensors;
  // mGeomSearch<DGCell*> *geomSearch;
  mDGMesh *theMesh;
  ConservationLaw *theLaw;
  DGLimiter *theLimiter;
  int compute_density_errorL2, compute_entropy_errorL2, compute_pressure_error,compute_pressure_errorL2;
  int compute_entropy_errorL1, compute_Linf_error;
  int compute_total_pressure, compute_total_mass_loss,compute_exact_total_pressure;
  int compute_lift_drag_coeffs, plot_pressure_coefficient,compute_density_errorL1;
  int plot_mach_on_surface, plot_pressure_loss_coefficient, plot_entropy_on_surface;
  double TEND, DT, TACT,REF_SAMPLING_TIME,CFLMAX;
  int ITER;
  int LEVEL_OF_REFINEMENT;
  _TIME_STEPPING_ timeSteppingMode;
  int n, HREF,PREF;
  int steps_to_dump; 
  char* recoverFile;
  TimeIntegrator *theIntegrator;
  GaussIntegrator theGaussIntegrator;
  Geometry *theGeometry;
  PostProGmshView *view;

  DGCell* addCell(mEntity *, int, int& dof);
  DGBoundaryCell* addBoundaryCell(mEntity *);
  int splitCell(mEntity*);
  void unsplitCell(mEntity*);	
  void adaptH ();
  void parse ();
  void modifyOrder(mEntity *, int order);
  void smoothP(mEntity *ent);
  void adaptP();
  void increaseOrder(int);
  void computeError(double &, double &);
  void computeFullJump();
  void ExactError(double, double &);
  void LinfError(double time, double &error);
  void entropyErrorL1(double time, double &);
  void entropyErrorL2(double time, double &);
  void densityErrorL2(double time,double &error);
  void PressureError(double time,double&);
  void pressureErrorL2(double time, double &error);
  void computeTotalPressure(double &error);
  void computeExactTotalPressure(double &error);
  void computeLiftDragCoeffs(double &,double &);
  void ApproxError(double &);
  void ComputeOrientation();
  void computeBoundaryIntegral(double &);
  void computeTotalMassLoss(double &totalMassLoss) const;
  void plotPressureCoefficient() const;
  void plotPressureLossCoefficient() const;
  void plotMachOnSurface() const;
  void plotEntropyOnSurface() const;
  
 public:
  DGAnalysis(mDGMesh *, ConservationLaw *, DGLimiter *, int);
  ~DGAnalysis();
  void run ();
  void addSensor  (DGSensor *);
  void evalSensors(bool doNow=0);
  void computeNormalsToCurvedBoundaries();
  void setPeriodicBC();
  void sortCellsBySize();
  virtual mEntity* eval ( const mPoint &, double *);
  virtual mEntity* eval ( const mPoint &, double * ,double &u, double &v, double &w);
  virtual ConservationLaw *getLaw(){return theLaw;}
  // export stuff...
  //virtual void exportMesh (char *,int);
  virtual void exportGmsh (char *file,int iteration);
  virtual void setupSummation();
  virtual void doSummation();
  virtual void exportSummation(char *fileName, int iteration, int);
  virtual void exportGmsh_Volume (char *file, int iteration);
  virtual void exportDataVise (ofstream &out);
  //virtual void exportSolutionOnVertices (ofstream &out);
  // write state of DGAnalysis to stream
  virtual void write(ostream &, ostream &);
  // read state of DGAnalysis from stream
  virtual void read(istream &);
};

#endif

