#include "parallel.h"

extern unsigned int numCPU;
extern unsigned int numThreads;
extern vector<infoWrapper> threadInfo[4];
extern vector<infoWrapper> threadInfoBoundary[4];

infoWrapper::infoWrapper(mMesh::iter begin, mMesh::iter end, double t, pthread_t id)
	: _begin(begin), _end(end), _t(t), _id(id){}

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