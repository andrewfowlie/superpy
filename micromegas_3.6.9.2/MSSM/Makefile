
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
 work/work_aux.a  $(wildcard lib/*.a)   $(cLib)/sqme_aux.$(SO) $(cLib)/libSLHAplus.a   \
 $(cLib)/num_c.a $(cLib)/ntools.a   $(cLib)/serv.a $(LX11) 


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
	$(FC) $(FFLAGS) -ffixed-line-length-none -o $(main:.F=) $(main) $(cLib)/faux.o  $(SSS)   $(lDL)  -lm
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
	rm -fr calchep/tmp/* calchep/results/*  lanhep/*.mdl 
	rm -f suspect2_lha.in suspect2_lha.out HB.in HB.slha nngg.out nngg.in Key.dat output.flha
	../sources/cleanexe
	$(MAKE) -C lib  clean
	$(MAKE) -C work clean 
