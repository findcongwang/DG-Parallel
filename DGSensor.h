#ifndef _DGSENSOR_H_
#define _DGSENSOR_H_

#include "mPoint.h"

class DGAnalysis;
class PostProGmshView;

class DGSensor
{
 protected:
  char *fileName;
  double DTsampleFrequency,TOLD;
 public:
  DGSensor(char *, double);
  virtual ~DGSensor();
  virtual void eval (DGAnalysis *, int, double, bool) = 0;
};

class DGLineSensor : public DGSensor
{
  mPoint p1,p2;
  int nbSamples;
 public :
  DGLineSensor (char *, const mPoint &, const mPoint &, int nbSamples ,double );
  virtual void eval (DGAnalysis *, int, double, bool);
};

class DGPointSensor : public DGSensor
{
  mPoint p1;
 public :
  DGPointSensor (char *, const mPoint &, double);
  virtual void eval (DGAnalysis *, int, double, bool);
};

class DGGmshSensor : public DGSensor
{
  int ITER;
  double u[1][3],v[1][3],w[1][3];
  double uu[4][3],vv[4][3],ww[4][3];
  public :
    DGGmshSensor(char *, double dt);
  virtual void eval (DGAnalysis *, int, double, bool);
};
class SummationSensor : public DGSensor
{
  int ITER, COUNTER;
  PostProGmshView* view;

  public :
  SummationSensor(char *, double dt);
  virtual void eval (DGAnalysis *, int, double, bool);
};
class DGDataViseSensor : public DGSensor
{
  int ITER;
  public :
    DGDataViseSensor(char *, double dt);
    virtual void eval (DGAnalysis *, int, double, bool);
};

class DGDXSensor : public DGSensor
{
  int ITER;
  public :
    DGDXSensor(char *, double dt);
    virtual void eval (DGAnalysis *, int, double, bool);
};

class DGTresholdLineSensor : public DGSensor
{
    int ithValueOfInterest;
    mPoint p1,p2;
    int nbSamples;
    double treshold;
  public :
    DGTresholdLineSensor (char *, 
			  int i,
			  const mPoint &, const mPoint &, 
			  int s, double t);
    virtual void eval (DGAnalysis *, int, double, bool);
};

#endif





























































