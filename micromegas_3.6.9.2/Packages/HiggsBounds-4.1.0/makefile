
# --- variables defined by configure ---
# If you want to edit these variables, change ./configure, not ./makefile


#----------------------------

# Put you favourite compiler and compiler options in here:

# e.g.

#  a) GNU compiler (http://gcc.gnu.org/wiki/GFortran): 
F90C = gfortran
F77C = gfortran
F90FLAGS = -fbounds-check -ffixed-line-length-none
# use the compiler flag -fopenmp for parallel processing

#  b) g95 compiler (www.g95.org):
#F90C = g95
#F77C = g95
#F90FLAGS = -fbounds-check -ffree-line-length-huge

#  c)  Portland Group Inc. compiler: 
#F90C = pgf90
#F77C = pgf90
#F90FLAGS = -C -Ktrap=fp -Mbounds

#  d) Intel compiler:
#F90C = ifort
#F77C = ifort
#F90FLAGS =  

#if using the NAG compiler, you need the compiler flag -DNAGf90Fortran
#because the subroutines flush,iargc,getarg need to use modules
#caution: the use of the NAG compiler is unsupported)

#(nb. the f77 compiler is currently only used for compiling the example programs)

#------------------------------

# If you wish to use HiggsBounds in conjunction with FeynHiggs or CPsuperH,
# make sure these paths indicate where you've stored these packages on your system,

# for example:

FHINCLUDE = -I/blah/blah/FeynHiggs-2.x.x/i586-linux/include
FHLIBS = -L/blah/blah/FeynHiggs-2.x.x/i586-linux/lib -lFH

CPSUPERHLIBS = -L/blah/blah/CPsuperH2 -lcpsuperh2 

# Note that care must be taken to ensure that consistent compilers are used
# If you get error messages similar to 
#   ... undefined reference to initialize_higgsbounds__...
# and the paths the libraries are correct, it's a good idea to check 
# that the same compilers are used for each package.
# (you may be able to work around this 
#   e.g. if FeynHiggs was compiled with g77 it may be possible to compile HiggsBounds
#  with gfortran and add -lg2c to the end of the variable FHLIBS)

#-----------------------------

HBLIBS = -L./ -lHB

EXE = HiggsBounds

#-----------------------------

#USECHISQMODS = $(CHISQMODS)
#ADDITIONALDEFINE = -DenableCHISQ

# --- end defs by configure ---


CHISQMODS = extra_bits_for_chisquared.MOD

MODS = usefulbits.mod store_pathname.mod \
 S95tables_type1.MOD S95tables_type2.MOD\
 S95tables_type3.MOD \
 interpolate.mod \
 theory_colliderSfunctions.MOD theory_XS_SM_functions.MOD theory_BRfunctions.MOD \
 S95tables.mod \
 PDGnumbering.mod string_manip.mod SLHA_manip.mod \
 extra_bits_for_SLHA.mod \
 extra_bits_for_web.MOD\
 $(USECHISQMODS) \
 input.MOD output.MOD channels.MOD theo_manip.mod \

tempMODS = $(MODS:.mod=.o)

OBJSbasic = $(tempMODS:.MOD=.o)

OBJScommandline = $(OBJSbasic) \
 HiggsBounds.o 

OBJSsubroutines =  $(OBJSbasic) \
 HiggsBounds_subroutines.o access_SM.o

.SUFFIXES: .exe .o .mod .f90 .F .F90 .MOD

#as advised in http://gcc.gnu.org/wiki/GfortranFAQ
%.o : %.mod

default: HiggsBounds

.f90.mod:
	$(F90C) $(F90FLAGS) $(ADDITIONALDEFINE) -c $< -o $*.o 

.f90.o:
	$(F90C) $(F90FLAGS) $(ADDITIONALDEFINE) -c $< -o $*.o 

.F90.MOD:
	$(F90C) $(F90FLAGS) $(ADDITIONALDEFINE) -c $< -o $*.o 

.F90.o:
	$(F90C) $(F90FLAGS) $(ADDITIONALDEFINE) -c $< -o $*.o 

.F.o:
	$(F77C) -c $< -o $*.o	

.mod.o:
	$(F90C) $(F90FLAGS) -c $*.f90 -o $*.o 

.MOD.o:
	$(F90C) $(F90FLAGS) -c $*.F90 -o $*.o 

HiggsBounds: HBwithSLHA libHB.a $(MODS) $(OBJScommandline) $(OBJSsubroutines)
	$(F90C) $(F90FLAGS) $(OBJScommandline) -o $(EXE) $(HBLIBS)
	$(F90C) $(F90FLAGS) AllAnalyses.F90 -o AllAnalyses $(HBLIBS)
	rm -f Expt_tables/S95_t1.binary 
	rm -f Expt_tables/S95_t2.binary 
	rm -f Theory_tables/BRSM.binary 
	touch Expt_tables/S95_t1.binary 
	touch Expt_tables/S95_t2.binary 
	touch Theory_tables/BRSM.binary 

libHB:  $(MODS) $(OBJSsubroutines)
	ar -rv libHB.a $(OBJSsubroutines)
	ranlib libHB.a	
	rm -f Expt_tables/S95_t1.binary 
	rm -f Expt_tables/S95_t2.binary 
	rm -f Theory_tables/BRSM.binary 	
	touch Expt_tables/S95_t1.binary 
	touch Expt_tables/S95_t2.binary 
	touch Theory_tables/BRSM.binary 

libHB.a: libHB

HBwithFH:  libHB.a
	$(F90C) $(F90FLAGS) $(FHINCLUDE) example_programs/HBwithFH.F -o example_programs/HBwithFH $(FHLIBS) $(HBLIBS)
	$(F90C) $(F90FLAGS) $(FHINCLUDE) example_programs/HBwithFH_dm.F -o example_programs/HBwithFH_dm $(FHLIBS) $(HBLIBS)	
	$(F90C) $(F90FLAGS) $(FHINCLUDE) example_programs/HBSLHAinputblocksfromFH_extras.F example_programs/HBSLHAinputblocksfromFH.F90 -o example_programs/HBSLHAinputblocksfromFH $(FHLIBS) $(HBLIBS)

HBwithCPsuperH:  libHB.a
	$(F77C) example_programs/HBwithCPsuperH.f -o example_programs/HBwithCPsuperH $(CPSUPERHLIBS) $(HBLIBS)
	@echo 'now run ./HBwithCPsuperH < HBwithCPsuperH.input in the folder example_programs'

HBweb: libHB.a HiggsBounds
	@echo 'check that WEBVERSION is defined'
	$(F90C) $(F90FLAGS) extract_SM_results_for_web.f90 -o extract_SM_results_for_web $(HBLIBS) 

HBwithSLHA: libHB.a
	$(F90C) $(F90FLAGS) example_programs/HBwithSLHA.F90 -o example_programs/HBwithSLHA $(HBLIBS)

HBchisq: libHB.a
	$(F90C) $(F90FLAGS) example_programs/HBchisq.F90 -o example_programs/HBchisq $(HBLIBS)
	$(F90C) $(F90FLAGS) example_programs/HBchisqwithSLHA.F90 -o example_programs/HBchisqwithSLHA $(HBLIBS)

clean:
	rm -f *.o *.mod *.MOD *.a
	rm -f store_pathname.f90

	rm -f *.o *.mod *.MOD *.a *~ 
	rm -f example_programs/*~ 
	rm -f store_pathname.f90
	rm -f example_programs/HBwithFH
	rm -f example_programs/HBwithCPsuperH
	rm -f example_programs/example-SM_vs_4thGen
	rm -f extract_SM_results_for_web
	rm -f HiggsBounds
#	rm -f example_data/*results.dat
#	rm -f example_data/*Key.dat
	rm -f cs-ratios_sigma-bg-Hb/Tevatron*~
	rm -f example_programs/example-4thGen-results.dat
	rm -f example_programs/example-SM-results.dat
	rm -f example_programs/Key.dat
	rm -f example_programs/HBwithCPsuperH_effC.f
	rm -f example_programs/HBwithFH_effC.F
	rm -f example_programs/debug_channels.txt
	rm -f example_programs/debug_predratio.txt
	rm -f README_old

