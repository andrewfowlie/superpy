.PHONY: clean

ifneq ($(MAKECMDGOALS),clean)
include ../../../CalcHEP_src/FlagsForMake
endif

micrO  = ../sources/micromegas

micrOr = ../microRun
cLib = $(CALCHEP)/lib

LIB =  ../sources/nmhdecay.a 
#   $(micrOr)/lib/aLib.a $(micrOr)/lib/nngg/nngg.a  $(micrO)/sources/micromegas.a $(micrOr)/work/work_aux.a \
#     $(cLib)/sqme_aux.so  $(cLib)/dynamic_me.a  $(cLib)/libSLHAplus.a $(cLib)/num_c.a  $(cLib)/ntools.a  $(cLib)/serv.a \
#      $(micrOr)/lib/LoopTools-2.1/libooptools.a



EXEC = nmhdecay nmspec

all: $(EXEC)

%.o: %.f
	$(FC) $(FFLAGS) -o $@ -c $<

nmhdecay:	nmhdecay.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ 
	
nmhdecay_grid:	nmhdecay_grid.o  $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ 

nmhdecay_rand:	nmhdecay_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmspec: 	nmspec.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ 

nmspec_grid: 	nmspec_grid.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmspec_rand: 	nmspec_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmgmsb: 	nmgmsb.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmgmsb_grid: 	nmgmsb_grid.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

nmgmsb_rand: 	nmgmsb_rand.o $(LIB)
	$(FC) $(FFLAGS) -o $@ $^ $(lDL) $(LX11)

clean:
	rm -f *.o $(EXEC)
