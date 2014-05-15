.PHONY: clean all

NMSSMTools = ../../Packages/NMSSMTools_4.1.2
CALCHEP =    ../../CalcHEP_src
all:  aLib.a $(NMSSMTools)/main/nmhdecay 

ifneq ($(MAKECMDGOALS),clean)
include $(CALCHEP)/FlagsForMake
endif

aLibOBJ = $(notdir $(subst .c,.o,$(wildcard *.c))) fort.o


aLib.a:aLib.a($(aLibOBJ))
	$(RANLIB) aLib.a
	rm -f $(aLibOBJ)

aLib.a(calchep_func.o):lpath.h

lpath.h:
	echo \#define LPATH \"`pwd`\"  >lpath.h
	echo \#define NMSSMTOOLS \"$(NMSSMTools)\" >> lpath.h


$(NMSSMTools)/main/nmhdecay:
	make -C $(NMSSMTools) init
	make -C $(NMSSMTools)/main


clean: 
	rm -f lpath.h *.o *.a  
	make -C nngg clean
