.PHONY: clean all

Tools = LoopTools-2.1

ifneq ($(MAKECMDGOALS),clean)
include ../../sources/micromegas/CalcHEP_src/FlagsForMake
endif

aLibOBJ = $(notdir $(subst .c,.o,$(wildcard *.c))) fort.o


all:  aLib.a  $(Tools)/libooptools.a nngg/lGamma.exe

$(Tools)/libooptools.a:
	make -C $(Tools)
	
nngg/lGamma.exe: $(Tools)/libooptools.a
	 make -C nngg

aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f $(aLibOBJ)

aLib.a(calchep_func.o):lpath.h

clean: 
	rm -f lpath.h *.o *.a  
	make -C nngg clean 
	make -C $(Tools) clean
