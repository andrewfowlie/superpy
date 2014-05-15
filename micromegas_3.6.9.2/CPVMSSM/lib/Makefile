.PHONY: all


ifneq ($(MAKECMDGOALS),clean)
include ../../CalcHEP_src/FlagsForMake
endif

CPSUPERH =  ../../Packages/CPsuperH2.3

all:  $(CPSUPERH)/cpsuperh.exe aLib.a 


OBJ = $(subst .c,.o,$(wildcard *.c)) 



aLib.a:  aLib.a($(OBJ))
	$(RANLIB) aLib.a
	rm -f $(OBJ)

aLib.a(cpHiggs.o):localpath.h cpHiggs.c

localpath.h:
	echo \#define CPsuperH   \"$(CURDIR)/$(CPSUPERH)/cpsuperh.exe\"  >  localpath.h
	        
$(CPSUPERH)/cpsuperh.exe:
	$(MAKE) -C $(CPSUPERH)

clean: 
	rm -f *.o *.a  localpath.h
	$(MAKE) -C nngg clean
  
