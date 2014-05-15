.PHONY: all clean

# read compiler flags
ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif

# path to your project files.

MDIR = $(subst /lib*?,,$(CURDIR)*?)
 
aLibOBJ = $(subst .c,.o,$(wildcard *.c))


all: aLib.a

aLib.a:aLib.a($(aLibOBJ))


clean: 
	rm -f *.o *.a *.so *.dll *.exe so_located
