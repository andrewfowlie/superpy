  
.PHONY: all clean

ifneq ($(MAKECMDGOALS),clean)
include ../CalcHEP_src/FlagsForMake
endif

all:micromegas.a 

OBJ = $(notdir $(subst .c,.o,$(filter-out $(special),$(wildcard *.c)))) 

#ifneq ($(FC),)
#   override OBJ += faux.o
#endif

micromegas.a:micromegas.a($(OBJ)) 
	$(RANLIB) micromegas.a 
	rm -f *.o

clean:
	rm -f *.a  *.o 
