#ifndef _DGMYERROR_H_

#define _DGMYERROR_H_

class mDGMesh;
class ConservationLaw; 
class DGMyError
{

mDGMesh *theMesh;
ConservationLaw *theLaw;

double TEND, DT, TACT,REF_SAMPLING_TIME,CFLMAX;

  int LEVEL_OF_REFINEMENT;

  int n, HREF,PREF;

void computeMyError(double dt);

public:
	DGMyError(mDGMesh *,ConservationLaw *);
};

#endif
