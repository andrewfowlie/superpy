
.PHONY: all libs clean

# read compiler flags
ifneq ($(MAKECMDGOALS),clean)
AllFlags = ../CalcHEP_src/FlagsForMake
ifeq (,$(wildcard $(AllFlags) )) 
$(error File $(AllFlags) is absent. Presumably you forgot to compile main code)
endif 
include ../CalcHEP_src/FlagsForMake 
endif 

cLib = $(CALCHEP)/lib
# files to compile
SSS = $(wildcard lib/*.a) ../sources/micromegas.a  $(cLib)/dynamic_me.a ../sources/micromegas.a \
 work/work_aux.a $(cLib)/sqme_aux.$(SO) $(cLib)/libSLHAplus.a   \
 ../CalcHEP_src/lib/num_c.a ../CalcHEP_src/lib/ntools.a   ../CalcHEP_src/lib/serv.a $(LX11) 


all: $(main) libs
ifeq (,$(main)) 
	@echo Main program is not specified. Use gmake main='<code of main program>'
else  
  ifeq (.cpp,$(suffix $(main)))
	$(CXX) $(CXXFLAGS) -o $(main:.cpp=) $(main)  $(SSS)  $(lDL)  -lm
  else  
    ifeq (.c,$(suffix $(main)))
	$(CC)  $(CFLAGS) -o $(main:.c=) $(main)  $(SSS) $(lDL)  -lm 
    else  
      ifeq (.F,$(suffix $(main)))
        ifeq ($(FC), )
	$(error  Fortran compiler was not detected.)
        else 
	$(FC) $(FFLAGS)  -ffixed-line-length-none   -o $(main:.F=) $(main) $(cLib)/faux.o  $(SSS)  $(lDL)  -lm
        endif
      else
	$(error Only .c, .cpp, and .F  main programs  are supported).
      endif
    endif
  endif  
endif


libs:
	$(MAKE) -C work
	$(MAKE) -C lib

clean: 
	$(MAKE) -C lib  clean
	$(MAKE) -C work clean 
	rm -fr HB.slha  calchep/tmp/* calchep/results/* lanhep/*.mdl
	../sources/cleanexe
