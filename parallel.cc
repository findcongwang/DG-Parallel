#include "parallel.h"

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

/*
void *parallelBoundary(void *ptr){
	infoComputeBoundary* info = (infoComputeBoundary*)ptr;

	printf("Thread %ld started\n", *(info->idptr));
	parallelBoundaryHelper(info->cell, info->t, info->checkPhysical, info->zeroList);
	printf("Thread %ld done\n", *(info->idptr));
	pthread_exit(0);
}
*/

/*
void *parallelBoundaryHelper(DGCell* cell, double t, bool checkPhysical, list<DGCell*> *zeroList){
	if(checkPhysical){
		int notPhysical;
		notPhysical = cell->computeBoundaryContributions(t);
		if (notPhysical){
			if (notPhysical==1) zeroList->push_back(cell->getLeftCell());
			else if (notPhysical==2) zeroList->push_back(cell->getRightCell());
			else{
				zeroList->push_back(cell->getLeftCell());
				zeroList->push_back(cell->getRightCell());
			}
		}
	}
	else{
		cell->computeBoundaryContributions(t);
	}
}
*/

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