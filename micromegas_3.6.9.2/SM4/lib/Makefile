.PHONY: all clean

# read compiler flags
ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif

# path to your project files.

MDIR = $(subst /lib*?,,$(CURDIR)*?)


all:
	echo External  function for the model are not needed. 

clean: 
	rm -f *.o *.a *.so *.dll *.exe so_located
