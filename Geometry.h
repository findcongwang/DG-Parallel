#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_
#include "DGCell.h"

using namespace std;

class DGCell;
class Geometry{
 //privite:
 public:
  Geometry(){}
  ~Geometry(){}
  DGCell* findReflPoint(DGCell*, double, double, double,double&,double&, double&);
  double getR(double x, double y);
};
#endif
