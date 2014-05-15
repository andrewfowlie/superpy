
# --- variables defined by configure ---
# If you want to edit these variables, change ./configure, not ./makefile

#----------------------------
# Put you favourite compiler and compiler options in here:

F90C = gfortran
F77C = gfortran
F90FLAGS = -fbounds-check -ffixed-line-length-none
# F90FLAGS = -fbounds-check -ffixed-line-length-none -fPIC
#F90FLAGS = -fbounds-check -Wall -W  

# Please specify the path to your HiggsBounds (version 4.1.0 or higher) distribution.
# Note that you have to build HiggsBounds successfully BEFORE you build HiggsSignals. 

HiggsBounds = ../HiggsBounds-4.1.0
HBLIBS =-L../HiggsBounds-4.1.0
HBINCLUDE =-I../HiggsBounds-4.1.0

#other useful gfortran flags:
# -pg           can then use gprof ./programname to look at how efficient program is
# -fopenmp      switches on parallel processing (see HiggsBounds.F90)
# -O3           sets a high level of optimisation 


#F90C = g95
#F77C = g95
#F90FLAGS = -fbounds-check -ffree-line-length-huge

#F90C = pgf90
#F77C = pgf90
#F90FLAGS = -C -Ktrap=fp -Mbounds

#F90C = f95
#F77C = f95
#F90FLAGS = -C 

#F90C = ifort
#F77C = ifort
#F90FLAGS = -C 

#if using the NAG compiler, you need the compiler flag -DNAGf90Fortran
#because the subroutines flush,iargc,getarg need to use modules
#caution: the use of the NAG compiler is unsupported

#------------------------------

# If you wish to use HiggsSignals in conjunction with FeynHiggs,
# make sure these paths indicate where you've stored these packages on your system,

# for example:

FHINCLUDE = -I/blah/blah/FeynHiggs-2.x.x/i586-linux/include
FHLIBS = -L/blah/blah/FeynHiggs-2.x.x/i586-linux/lib -lFH

# Note that care must be taken to ensure that consistent compilers are used
# If you get error messages similar to 
#   ... undefined reference to initialize_higgsbounds__...
# and the paths the libraries are correct, it's a good idea to check 
# that the same compilers are used for each package.
# (you may be able to work around this 
#   e.g. if FeynHiggs was compiled with g77 it may be possible to compile HiggsBounds
#  with gfortran and add -lg2c to the end of the variable FHLIBS)

#-----------------------------


HSLIBS = -L./ -lHS

EXE = HiggsSignals

#-----------------------------




MODS = usefulbits_HS.mod store_pathname_HS.mod datatables.mod combinatorics.mod \
	   numerics.mod io.mod expt_syst.mod pc_chisq.mod mc_chisq.mod all_chisq.mod

tempMODS = $(MODS:.mod=.o)

OBJSbasic = $(tempMODS:.MOD=.o)

OBJScommandline = $(OBJSbasic) \
 HiggsSignals.o

OBJSsubroutines =  $(OBJSbasic) \
 HiggsSignals_subroutines.o

.SUFFIXES: .exe .o .mod .f90 .F .F90 .MOD

#as advised in http://gcc.gnu.org/wiki/GfortranFAQ
%.o : %.mod

default: HiggsSignals

.f90.mod:
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) $(ADDITIONALDEFINE) -c $< -o $*.o -lHB

.f90.o:
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) $(ADDITIONALDEFINE) -c $< -o $*.o -lHB

.F90.MOD:
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) $(ADDITIONALDEFINE) -c $< -o $*.o -lHB

.F90.o:
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) $(ADDITIONALDEFINE) -c $< -o $*.o -lHB

.F.o:
	$(F77C) -c $< -o $*.o

.mod.o: $(HiggsBounds)/libHB.a
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) -c $*.f90 -o $*.o -lHB

.MOD.o: $(HiggsBounds)/libHB.a
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) -c $*.F90 -o $*.o -lHB

HiggsSignals: $(HiggsBounds)/libHB.a libHS.a $(MODS) $(OBJScommandline)
	$(F90C) $(F90FLAGS) $(HBINCLUDE) $(OBJScommandline) -o $(EXE) $(HSLIBS) -lHS $(HBLIBS) -lHB

libHS:  $(MODS) $(OBJSsubroutines)
	ar -rv libHS.a $(OBJSsubroutines)
	ranlib libHS.a	

libHS.a: libHS

$(HiggsBounds)/libHB.a:
	cd $(HiggsBounds); ./configure; $(MAKE)


HSexamples:  libHS.a
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HSeffC.f90 -o example_programs/HSeffC $(HSLIBS) -lHB
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HShadr.f90 -o example_programs/HShadr $(HSLIBS) -lHB
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HSscaleUncertainties.f90 -o example_programs/HSscaleUncertainties $(HSLIBS) -lHB
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HSwithToys.f90 -o example_programs/HSwithToys $(HSLIBS) -lHB
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HS_efficiencies.f90 -o example_programs/HS_efficiencies $(HSLIBS) -lHB	
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HSwithSLHA.f90 -o example_programs/HSwithSLHA $(HSLIBS) -lHB
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HBandHSwithSLHA.f90 -o example_programs/HBandHSwithSLHA $(HSLIBS) -lHB
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HS_mass.f90 -o example_programs/HS_mass $(HSLIBS) -lHB

HSeffC: libHS.a
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HSeffC.f90 -o example_programs/HSeffC $(HSLIBS) -lHB

HS_mass: libHS.a
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HS_mass.f90 -o example_programs/HS_mass $(HSLIBS) -lHB

HShadr: libHS.a
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HShadr.f90 -o example_programs/HShadr $(HSLIBS) -lHB

HSscaleUncertainties: libHS.a
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HSscaleUncertainties.f90 -o example_programs/HSscaleUncertainties $(HSLIBS) -lHB

HSwithToys: libHS.a
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HSwithToys.f90 -o example_programs/HSwithToys $(HSLIBS) -lHB

HS_efficiencies: libHS.a
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HS_efficiencies.f90 -o example_programs/HS_efficiencies $(HSLIBS) -lHB

HSwithSLHA: libHS.a
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HSwithSLHA.f90 -o example_programs/HSwithSLHA $(HSLIBS) -lHB

HBandHSwithSLHA: libHS.a
	$(F90C) $(F90FLAGS) $(HBLIBS) $(HBINCLUDE) example_programs/HBandHSwithSLHA.f90 -o example_programs/HBandHSwithSLHA $(HSLIBS) -lHB

HBandHSwithFH: libHS.a
	$(F77C) $(FHINCLUDE) $(HBLIBS) $(HBINCLUDE) example_programs/HBandHSwithFH.F -o example_programs/HBandHSwithFH $(FHLIBS) $(HSLIBS) -lHB
	
clean:
	rm -f *.o *.mod *.MOD *.a
	rm -f store_pathname_HS.f90
	rm -f *.o *.mod *.MOD *.a *~ 
	rm -f example_programs/*~
	rm -f example_data/*-fromHB
	rm -f example_data/*-fromHS
	rm -f example_data/*-fromHBandHS
	rm -f store_pathname_HS.f90	
	rm -f example_programs/*.txt	
	rm -f example_programs/tmp/*
	rm -f Expt_tables/mutables.binary
	rm -f Expt_tables/analyses.txt	
	rm -f temp*.txt
	rm -f HiggsSignals
