* ltcache.h
* declarations for the LoopTools cache
* this file is part of FeynHiggs
* last modified 21 Sep 12 th


	memindex cacheindex
	external cacheindex

	integer ncaches
	parameter (ncaches = 2)

	ComplexType cache(2,ncaches)
	common /ltcache/ cache
