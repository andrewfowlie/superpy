.PHONY: all clean

ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake 
endif


all:  aLib.a

aLibOBJ = $(subst .c,.o,$(wildcard *.c))


aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f *.o

clean: 
	rm -f *.o *.a 
