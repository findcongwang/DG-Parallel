#ifndef H_ErrorEstimator
#define H_ErrorEstimator

#include "Field.h"
#include "MFace.h"
#include "InterpolationVisitor.h"

class pCompute
{
  public :
    static int p(Field<DofScalar > *f1, Field<DofArray<4> >*,
		 MFace *theFace,
		 double maxerr, double minerr,
		 double globalError, int &maximumP, int &minimumP);
    static int smoothP(MFace *theFace, Field<DofArray<4> >*);
};

template <class DofType>
class ErrorEstimator : public InterpolationVisitor<DofType>
{
public:
  ErrorEstimator(Field<DofScalar > *f1, 
		 Field<DofType > *f2, 
		 int pMax, 
		 int pMin,
		 int freq,
		 int totalNumberOfUnknownsPrescribed);
  void performErrorEstimation();
  virtual void visit (Interpolation<DofType > *interp);
protected:
  int whatToDoInVisit;
  Field<DofScalar > *theErrorField;
  Field<DofType >   *theField;
  int errorEstimatorFrequency;
  int totalNumberOfUnknownsPrescribed;
  int pMin,pMax;
  int iter;
  double theDenominator,MaxErr,MinErr;
  double globalError;
};

template <class DofType>
ErrorEstimator<DofType > ::ErrorEstimator(Field<DofScalar > *f, Field<DofType > *f2, 
					  int p1, int p2, int fq, int t)
  : theErrorField(f),theField(f2),pMax(p1),pMin(p2),errorEstimatorFrequency(fq),
    totalNumberOfUnknownsPrescribed(t),iter(0)
{
}

template <class DofType>
void ErrorEstimator<DofType > ::performErrorEstimation ()
{ 
  pMax = 0;
  pMin = 100;
  MaxErr = 0.0;
  MinErr = 1.e10;
  iter ++;
  if(iter == 1 || iter % errorEstimatorFrequency == 0)
    {
      globalError = 0.0;
      std::cout << "error estimation iter " << iter <<"\n";
      whatToDoInVisit = 0;
      theField->visitAllInterpolations(this);
      globalError = sqrt(globalError);
      std::cout << "p refinement ...\n";
      whatToDoInVisit = 1;
      theField->visitAllInterpolations(this);
      std::cout << "smoothing ...\n";
      whatToDoInVisit = 2;
      for(int i=pMin;i<=pMax;i++)
	theField->visitAllInterpolations(this);
    }
  whatToDoInVisit = 3;
  theField->visitAllInterpolations(this);
}

template <class DofType>
void ErrorEstimator<DofType >:: visit (Interpolation<DofType> *interp)
{
  MFace *theFace = (MFace*)interp->meshEnt();
  if(theFace->dim() != 2)return;
  DofGroup *dof = theFace->getDof(*theErrorField,0);
  double error = dof->cValue(0,0);
  switch (whatToDoInVisit)
    {
    case 0://performs error estimation
      // actually, the error estimation is done      
      globalError += error;
      if(error > MaxErr)MaxErr = error;
      if(error < MinErr)MinErr = error;
      break;
    case 1://adapt p
      theField->setInterpOrder(theFace,pCompute::p(theErrorField,theField,theFace,MaxErr,MinErr,globalError,pMax,pMin));
      break;
    case 2: //smooth
      theField->setInterpOrder(theFace,pCompute::smoothP(theFace,theField));
      break;      
    case 3: //cleanup
      dof->cSetValue(0,0.0,0);
      break;      
    }
}



#endif
