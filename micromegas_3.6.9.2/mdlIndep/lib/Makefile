.PHONY: all clean

# read compiler flags
ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif

# path to your project files.

MDIR = $(subst /lib*?,,$(CURDIR)*?)

aLibOBJ = $(subst .c,.o,$(wildcard *.c))

ifneq ($(aLibOBJ), )

all:aLib.a
aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f $(aLibOBJ)
else

all:	
	echo External  function for the model are not needed. 
endif

clean: 
	rm -f *.o *.a *.so *.dll *.exe so_located
