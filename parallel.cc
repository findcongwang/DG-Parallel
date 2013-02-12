#include "parallel.h"
#include "TimeIntegrator.h"

extern unsigned int numCPU;
extern unsigned int numThreads;
extern vector<infoWrapper> threadInfo[4];
extern vector<infoWrapper> threadInfoBoundary[4];

/*mMesh::iter itCopy = it;
++itCopy;
mEntity *mNext = (*itCopy);
__builtin_prefetch(mNext->getCell(), 1, 3);*/

timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

infoWrapper::infoWrapper(mMesh::iter begin, mMesh::iter end, double t, pthread_t id)
	: _begin(begin), _end(end), _t(t), _id(id)
{
	_moreinfo = NULL;
}

void *parallelVolume(void *ptr){
	infoWrapper* info = (infoWrapper*)ptr;
	mMesh::iter it;
	//printf("Thread %ld started\n", *(info->idptr));
	for(it = info->_begin; it != info->_end; ++it)
    {
		mEntity *m = (*it);
		DGCell *cell = (DGCell*)m->getCell();
		cell->computeVolumeContribution(info->_t);		
    }
	//printf("Thread %ld done\n", *(info->idptr));
	pthread_exit(0);
}

void *parallelBoundary(void *ptr){
	infoWrapper* info = (infoWrapper*)ptr;
    mMesh::iter it;
    int notPhysical;

    for(it = info->_begin; it != info->_end; ++it)
    {
        mEntity *m = (*it);
        DGBoundaryCell *cell = (DGBoundaryCell*)m->getCell();
        notPhysical = cell->computeBoundaryContributions(info->_t);

	    switch(notPhysical)
	    {
		case 0:
			//nothing to do
			break;
		case 1:
			cell->getLeftCell()->setToZero(info->_t);
			break;
		case 2:
			cell->getRightCell()->setToZero(info->_t);
			break;
		default:
			cell->getLeftCell()->setToZero(info->_t);
			cell->getRightCell()->setToZero(info->_t);
			break;
	    }
    }

	pthread_exit(0);
}

void *parallelRungeKuttaTVD2K1(void *ptr)
{
	infoWrapper* info = (infoWrapper*)ptr;
	RungeKuttaTVD2 *rk = (RungeKuttaTVD2*)info->_moreinfo;
	mMesh::iter it;

	double dt, dt_over_detJac, dof;
    int i, j, k, fSizek;
    double dx;
    rk->soln = rk->soln_begin; 

	for(it = info->_begin; it != info->_end; ++it)
    {
		mEntity *m = (*it);
		DGCell *cell = (DGCell*)m->getCell();
		
		dof = cell->fSize * rk->cSize;
        const double *RHS = cell->theRightHandSide;
        double* coeff = cell->theFieldsCoefficients->get();
        /* Copy initial values*/

        for(i = 0; i < dof; i++) * (rk->soln++) = coeff[i];

        dt = cell->getTimeStep();
        if(cell->theFunctionSpace->isOrthogonal())
        {
            dt_over_detJac = dt/cell->getDetJac();
            for(i = 0;i<dof;i++) coeff[i] += RHS[i]*dt_over_detJac;
        }
        else
        {
            int fSize = cell->fSize;
            for(k = 0; k < rk->cSize; k++)
            {
                int fSizek = fSize*k;
                for(i = 0;i<fSize;i++)
                {
                    dx = 0.0;
                    for(j = 0;j<fSize;j++){
                        dx += cell->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSizek];
                    }
                    cell->theFieldsCoefficients ->get(k,i)+= dx*dt;
                }
            }
        }
        cell->computeMean();
        cell->ZeroRHS();
    }	
}

void *parallelRungeKuttaTVD2K2(void *ptr)
{
	infoWrapper* info = (infoWrapper*)ptr;
	RungeKuttaTVD2 *rk = (RungeKuttaTVD2*)info->_moreinfo;
	mMesh::iter it;

	double dt, dt_over_detJac, dof;
    int i, j, k, fSizek;
    double dx;
    
	/*  Reset c values */
    double resid = 0.0, tmp;
    rk->soln = rk->soln_begin; 

	for(it = info->_begin; it != info->_end; ++it)
    {
		mEntity *m = (*it);
		DGCell *cell = (DGCell*)m->getCell();
		
		dof = cell->fSize*rk->cSize; 
        const double *RHS = cell->theRightHandSide;
        double* coeff = cell->theFieldsCoefficients->get();
        dt = cell->getTimeStep(); 
        if(cell->theFunctionSpace->isOrthogonal())
        {
            dt_over_detJac = dt /cell->getDetJac();
            for(i = 0;i<dof;i++)
            {
                coeff[i] = (coeff[i]+RHS[i]*dt_over_detJac+(*rk->soln))*0.5;
                tmp = coeff[i]-(*(rk->soln++));
                resid+=tmp*tmp;
            }
        }
        else
        {
            int fSize = cell->fSize;
            for(k = 0;k<rk->cSize;k++)
            {
                fSizek = fSize*k;
                for(i = 0;i<fSize;i++)
                {
                    dx = 0.0;
                    for(j = 0;j<fSize;j++)
                    {
                        dx += cell->theInvertMassMatrix[i][j] * cell->theRightHandSide[j+fSizek];
                    }
                    cell->theFieldsCoefficients->get(k,i) = ( cell->theFieldsCoefficients->get(k,i)+
                                                            dx*dt +  (*rk->soln))*0.5;
                    tmp = cell->theFieldsCoefficients->get(k,i)-(*(rk->soln++));
                    resid+=tmp*tmp;
                }
            }
        }
        cell->computeMean();
        cell->ZeroRHS();
    }	
}