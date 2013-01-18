#ifndef _DGCELL_H_
#define _DGCELL_H_
#include <vector>
#include <list>
#include "mPoint.h"
#include "mVector.h"
#include "mAttachableDataContainer.h"
#include "Constants.h"
#include "FunctionSpace.h"
#include "ConservationLaw.h"

class mEntity;
class Mapping;
class FunctionSpace;
class FieldEvaluator;
class DGAnalysis;
class GaussIntegrator;
class Geometry;
using namespace std;
extern double ** allocateMatrix(int n);

class DG_Dofs 
{
  public :
    short NU,Nb,FS;
  ~DG_Dofs();
  double **theFieldsCoefficients; 
  DG_Dofs();
  DG_Dofs(short Nb, short NbUnknowns, short FunctionSpaceSize);
  DG_Dofs* operator= (const DG_Dofs*);
  double & get (int i, int j)  {return theFieldsCoefficients[0][i*FS+j];}
  double* get() {return theFieldsCoefficients[0];}
  double & getCopy (int i, int j)  {return theFieldsCoefficients[1][i*FS+j];}
  void copy();
  void swap();
  void copyBack();
  void copyLimitedCoeffs(int);
  void eval (vector<double> & sf, double *field, double t);
  void eval (vector<double> & sf, double *field);
};

class volumeGaussPoint
{
 public:
  double x,y,z;  //coordinates in physical space
  double JacTimesWeight; //precompute (Jac*weight)
  vector<mVector> grads;
  double*  fcts;
  volumeGaussPoint(){}
  ~volumeGaussPoint(){grads.clear();}
  volumeGaussPoint(const volumeGaussPoint &);
};

class DGCell:public mAttachableData
{
protected:
  ConservationLaw	*theConservationLaw;
  Geometry	     	*theGeometry;
  mEntity	        *theMeshEntity;
  double                *theRightHandSide;
  vector<double>	theErrorRightHandSide;
  vector<volumeGaussPoint>  volumeGaussPoints;
  DG_Dofs       	*theFieldsCoefficients;
  vector<double>	*theErrorCoefficients;
  FunctionSpace		*theFunctionSpace;  
  FunctionSpace		*theErrorFunctionSpace;
  Mapping	        *theMapping;
  GaussIntegrator       *theGaussIntegrator;
  double                **theInvertMassMatrix, *theMean,*theErrorMean,*limSlope,error, jump, fullJump,
						**theErrorMassMatrix,**theErrorInvertMassMatrix;
  double                errNum,errDen;
  mPoint                pMin,pMax; 
  short                 type; //???
  short			limit; 
  short			order;
  short			fSize;           // size of the Function Space
  short		      	fOrder;			//order of approximation
  short		       	cSize;            // number of equations
  short		       	nbPtGauss;
  short			binNb;
  double	       	cellSize;
  double	       	volume;
  double	       	perimeter;
  double	       	maxH, minH,jumpError, maxVal;
  double                detJac;
  double                dt;
  list<mEntity *>	allSubs;
  list<mEntity *>	allVert;

  void                  computeMassMatrices();
  void			allocateErrorMassMatrices();
  void          allocateInvJacMatrix(){DInv=allocateMatrix(cSize);}
  virtual void          init();
  double arclength(mEntity* theEdge, GaussIntegrator* theGaussIntegrator) const;

  public:
  
  double                *deltaUp0;
  double                **DInv;
  double                *Up0;
  double                *Rp0;
  
  short                 edgeOrientation[3];
  void                  setType();///???
  DGCell();
  DGCell(ConservationLaw*, mEntity*, FunctionSpace*,FunctionSpace *, Mapping *,GaussIntegrator *);
  virtual ~DGCell();
  void cleanup();
  volumeGaussPoint &  pt(int i) {return volumeGaussPoints[i];}
  Mapping *getMapping(){return theMapping;}
  FunctionSpace *getFunctionSpace(){return theFunctionSpace;}
  ConservationLaw *getConservationLaw(){return theConservationLaw;}
  mEntity* getMeshEntity(){return theMeshEntity;}
  list<mEntity*> getAllNeighbors() {return allVert;} 
  void computeVolumeContribution(double time);
  void reverseVolumeContribution (double t);
  void computeErrorVolumeContribution();
  void computeErrorMassMatrices();
  double advanceInTime(double dt);
  double advanceErrorInTime(double dt);
  void computeMean(); 
  void computeErrorMean();
  void setToZero(double t);
  void adaptTimeStep (double CFLMAX, double &DT);
  void adaptOrder(int order);
  void L2Proj        (FieldEvaluator *f, double *proj);
  void L2Proj (double* dU);
  void L2ProjError    (FieldEvaluator *f, vector<double> &proj);
  void L2ProjInitial (FieldEvaluator *f);
  void L2ProjInitialError (FieldEvaluator *f);
  virtual void L2ProjCutCell(){}
  void invertErrorMassMatrices (); 
  void ZeroRHS ()       // RESET RIGHT HAND SIDE
    {
      const int s = cSize*fSize;
      for(int j=0;j<s;++j)  theRightHandSide[j]= 0.0;
    }
  void ZeroErrorRHS ();
  void ZeroError (); 
  void ZeroJump ();
  void ZeroFullJump() {fullJump=0.0;}
  void ZeroErrorMassMatrices();
  inline void interpolateError (double u,double v,double w, double *field); 
  void dinterpolate (double u,  double v, double w, mVector *grads);  
  double L1exact(double time);
  double L2exact();
  double L1approx();
  double LinfError(double);
  double L1exactPressure(double time);
  double entropyErrorL1_2D(double time);
  double entropyErrorL1_3D(double time);
  double entropyErrorL2_2D(double time);
  double entropyErrorL2_3D(double time);
  double pressureErrorL2(double time);
  double densityErrorL2(double time);
  friend class DGBoundaryCell;
  friend class CutBoundaryCell;
  friend class CutCell;
  friend class DGAnalysis;
  friend class DGLimiter;
  friend class DGSuperBee;
  friend class DGBarthLimiter;
  friend class DGBarthLimiterEuler;
  friend class VertLimiter;
  friend class DGVertexLimiter;
  friend class DGVertexLimiterEuler;
  friend class DGPhysicalLimiter;
  friend class DGExpFilterLimiter;
  friend class DGMomentLimiter;
  friend class DGMomentLimiterEuler;
  friend class TimeIntegrator;
  friend class RungeKutta4;
  friend class RungeKuttaTVD2;
  friend class RungeKuttaTVD2adaptive;
  friend class RungeKuttaTVD3;
  friend class ForwardEuler;
  friend class Multigrid;
  void getBox(double &,double &,double &, double &,double &, double &);
  bool inCell (const mPoint &, double *val); 
  int limitStatus() const {return limit;}
  double getError() const; 
  void computeMaxH (); 
  void computeMinH (); 
  double getMaxH() const; 
  double getMinH() const;
  double getDetJac() const {return detJac;}
  void computeVolume();
  double getVolume() const {return volume;}
  void computePerimeter();
  double getPerimeter() const {return perimeter;}
  double getJumpError() const;
 GaussIntegrator*  getGaussIntegrator(){return theGaussIntegrator;}
  double computeMass();
  double limiter() const {return limSlope[0];}
  // write state of DGCell to stream
  void write(ostream &);
  // read state of DGCell from stream
  void read(istream &);
  void beginTimeStep(double);
  void computeCellSize();
  double getSize() {return cellSize;}
  void setBinNb(int b) {binNb=b;}
  short getBinNb() {return binNb;}
  double computeMAXEV();
  void computeTimeStep(double cfl)
    {
       dt=cellSize*cfl/((2.*fOrder+1.)*computeMAXEV());
    }
  double getTimeStep(){return dt;}
  void setTimeStep(double t){dt = t;}
  virtual void interpolateRefl(const double* fct, mVector &N, DGCell* reflCell, double *Q) const {
  double tmp=1.0;
  }
  inline void interpolate (const double *fct, double *field) const
    {
      const double *a =&(theFieldsCoefficients->theFieldsCoefficients[0])[0];
      for(int i=0;i<cSize;++i)
	{
	  double tmp = 0.0;
	  for(int j=0;j<fSize;++j) 
		  tmp +=fct[j]*a[i*fSize+j];
	  field[i]=tmp;
	}
  }
  inline void interpolate(const double* fct,double *val,DGCell* reflCell);
  inline void interpolate (double u, double v, double w,
			   double *field) ;
};

class boundaryGaussPoint
{
 protected:
  boundaryGaussPoint& operator = (const boundaryGaussPoint &){return *this;}
  boundaryGaussPoint(const boundaryGaussPoint &){}
 public:
  boundaryGaussPoint();
  ~boundaryGaussPoint();
  double x,y,z,JacTimesWeight;
  double *fctleft,*fctrght;
  mVector N; // normal to the circle, approximating a curved boundary
  mVector n; //normal to a curved edge
};

class DGBoundaryCell : public mAttachableData
{
protected:
  mEntity *mleft, *mright;
  DGCell *left, *right;
  mEntity *theBoundaryEntity;
  Mapping *theMapping;
  GaussIntegrator       *theGaussIntegrator;
  vector<boundaryGaussPoint*> boundaryGaussPoints;
  mVector n; //normal for a straight edge
  int cSize;       // number of equations
  int rSize;//size of the functional space on the right 
  int lSize;//size of the functional space on the left
  int nbPtGauss;// No of integrationpoints
  double size;   //size of the cell
  double detJac;
  
  public:
  void computeDeltaF(int Side, double* Up0, double* deltaUp0,double* flux);
  double computeMaxEVal();
  void computeJacobian(int Side, double* Up0, double** Jac);
  void check();
  DGBoundaryCell(mEntity *, Mapping*,GaussIntegrator *);
  virtual ~DGBoundaryCell();
  virtual int computeBoundaryContributions(double T);
  virtual void reverseBoundaryContributions(double t);
  void computeErrorBoundaryContributions(double T);
  // access to data
  mEntity *getEntity() const {return theBoundaryEntity;}
  Mapping *getMapping() const {return theMapping;}
  mEntity *getLeft() const {return mleft;}
  mEntity *getRight() const{return mright;}
  DGCell* getLeftCell() const {return left;}
  DGCell* getRightCell() const{return right;}
  void computeSize(); 
  double getSize() const {return size;} 
  boundaryGaussPoint * pt (int i) const {return boundaryGaussPoints[i];}
  int getNbPtGauss() const {return nbPtGauss;}
  virtual void init();
  friend class DGLimiter;
  friend class DGBarthLimiter;
  friend class DGBarthLimiterEuler;
  friend class VertLimiter;
  friend class DGVertexLimiter;
  friend class DGVertexLimiterEuler;
  friend class DGSuperBee;
  friend class DGPhysicalLimiter;
  friend class DGAnalysis;
  friend class TimeIntegrator;
  friend class RungeKutta4;
  friend class RungeKuttaTVD2;
  friend class ForwardEuler;
  friend class Multigrid;
  friend class DGCell;
  friend class CutCell;
  friend class CutBoundaryCell;
  int computeOrder() const;
  void computeCurvedNormals(mPoint &,mPoint &, mPoint &, mVector &, mVector &, mVector &);
  int computeErrorOrder(DGCell *left, DGCell *right) const;
  void computeError(); 
  void computeJump();
  void computeFullJump();
  void computeOrientation();
  void computeBoundaryRestrictions(double T); 
  void computeErrorJump(double T); 
  double BoundaryIntegral();
  double totalPressure(double &);
  void computeLiftDrag(double &,double &);
  void computePressureCoefficient(int &,double *x,double *val) const;
  void computeTotalPressureAtAPoint(int &nbPoints,double *x,double *y,double *p_t)const;
  void computeMachNumber(int &nbPoints,double *x,double *m) const;
  void computeEntropy(int &nbPoints,double *x,double *m) const;
  void normalToCircle(mPoint &, mPoint &,mPoint &);
  void normalToCircle(double center_x, double center_y, double center_z);
  double computeRadius(mPoint &, mPoint &,mPoint &);
  double computeSphere(mPoint p[4]);
  double determinant(double **a, int n);
  void setPeriodicBC();
  void IndicatorGrad (double *denoms, double *numers) const;
  // write state of DGBoundaryCell to stream
  void write(ostream &){}
  // read state of DGBoundaryCell from stream
  void read(istream &){}
};

/*****************  INLINING  ****************************/

inline void DGCell::interpolate(double u, double v, double w, double *field)
{
  double fct[MaxNbFunctions]; 
  theFunctionSpace->fcts(u,v,w,fct);
  int i,j;
  double *tmp=field;
  const double *a =&(theFieldsCoefficients->theFieldsCoefficients[0])[0];
  const double *k;
  --a;
  for(i=0;i<cSize;++i)
    {
      *tmp = 0.0;
      k=fct;
      for(j=0;j<fSize;++j) (*tmp) += (*(k++))*(*(++a));
      ++tmp;
    }
}

#endif











