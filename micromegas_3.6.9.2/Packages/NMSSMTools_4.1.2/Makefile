.PHONY: clean init

ifeq ($(MAKECMDGOALS),init)
init: 
#	$(MAKE) -C sources/micromegas
	$(MAKE) -C sources
#	$(MAKE) -C microRun  main=main.c

#microRun/lib/lpath.h:
#	echo \#define LPATH \"`pwd`\" > microRun/lib/lpath.h

else 

ifeq ($(MAKECMDGOALS),clean)
clean:
#	$(MAKE) -C sources/micromegas clean
	$(MAKE) -C sources clean
	$(MAKE) -C main clean
#	$(MAKE) -C microRun clean
#	-find . -name '*.DS_Store' -print | xargs rm -rf
#	-find . -name '*.dSYM' -print | xargs rm -rf
#	-find . -name '._*' -print | xargs rm -rf
else

#micrO = sources/micromegas

#ifeq ($(wildcard $(micrO)/CalcHEP_src/FlagsForMake),)
#$(error Use '[g]make init' for initialization)
#endif

all:
	$(MAKE) -C main

endif

endif
