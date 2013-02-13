#ifndef _PARALLEL_H_
#define _PARALLEL_H_

#include "DGCell.h"
#include "mDGMesh.h"
#include "mMesh.h"
#include <list>
#include <vector>
#include <pthread.h>
#include <time.h>

class infoWrapper{
public:
	pthread_t _id;
	mMesh::iter _begin;
	mMesh::iter _end;
	double _t;
	int _n;
	void *_moreinfo;
	infoWrapper(mMesh::iter begin, mMesh::iter end, double t, pthread_t id);
};

void *parallelVolume(void *info);
void *parallelBoundary(void *info);

//RungeKutta2
void *parallelRungeKuttaTVD2K1(void *info);
void *parallelRungeKuttaTVD2K2(void *info);

void *parallelLinfError(void *info);
void *parallelExportGmshP1(void *info);

timespec diff(timespec start, timespec end);

#endif