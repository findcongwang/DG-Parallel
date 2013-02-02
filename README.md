DG-Parallel
===========

TODO:
-----
+ Redo cell and boundary cell structures
	- isolate hot and cold fields to reduce structure size
	- group likely accessed together members together (within the structure)
	- arrange members in decreasing size
	- plan size to be cache aligned
	- use restricted const pointers to basis functions (reduce aliasing)
	- use prefetch functions to preload neighbours into memory
+ Design new data structures to store mesh
	
Runtime trials
--------------
1. with restrict const pointer to basis functions
2. with cache aligned structure size
3. with prefetching
