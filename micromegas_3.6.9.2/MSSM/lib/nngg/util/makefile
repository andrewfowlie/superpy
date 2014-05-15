

LIB = ../nngg.a


VEGAS_C = vegas/Vegas.c
VEGAS_F = vegas/vegas-f.c
VEGAS_H = vegas/decl.h
VEGAS_SRCS = vegas/util.c vegas/debug.c vegas/common.c \
  vegas/Sobol.c vegas/ChiSquare.c \
  vegas/Grid.c vegas/Integrate.c

SUAVE_C = suave/Suave.c
SUAVE_F = suave/suave-f.c
SUAVE_H = suave/decl.h
SUAVE_SRCS = suave/util.c suave/debug.c suave/common.c \
  suave/Sobol.c suave/ChiSquare.c suave/Grid.c \
  suave/Fluct.c suave/Sample.c suave/Integrate.c

DIVONNE_C = divonne/Divonne.c
DIVONNE_F = divonne/divonne-f.c
DIVONNE_H = divonne/decl.h
DIVONNE_SRCS = divonne/util.c divonne/debug.c divonne/common.c \
  divonne/KorobovCoeff.c divonne/Sobol.c divonne/ChiSquare.c \
  divonne/Rule.c divonne/Sample.c divonne/FindMinimum.c \
  divonne/Explore.c divonne/Split.c divonne/Integrate.c \

CUHRE_C = cuhre/Cuhre.c
CUHRE_F = cuhre/cuhre-f.c
CUHRE_H = cuhre/decl.h
CUHRE_SRCS = cuhre/util.c cuhre/debug.c cuhre/common.c \
  cuhre/ChiSquare.c cuhre/Rule.c cuhre/Integrate.c

#CUBA = \
#  $(VEGAS_C) $(VEGAS_F) \
#  $(SUAVE_C) $(SUAVE_F) \
#  $(DIVONNE_C) $(DIVONNE_F) \
#  $(CUHRE_C) $(CUHRE_F)

CUBA = Vegas.o vegas-f.o Suave.o suave-f.o Divonne.o divonne-f.o Cuhre.o cuhre-f.o


VPATH := $(VPATH) ./kin ./vegas ./suave ./diag ./univariate ./system

OBJS = \
  logfile.o \
  Delta.o \
  ThreeMom.o \
  Li2.o \
  VecSet.o \
  SInvariant.o \
  TInvariant.o \
  Pair.o \
  Eps.o \
  SxS.o \
  SeS.o \
  VxS.o \
  VeS.o \
  BxS.o \
  BeS.o \
  Jacobi.o \
  Eigen2x2.o \
  Eigensystem.o \
  SingularValues.o \
  PseudoEigensystem.o \
  ApplyEigen.o \
  ApplySVD.o \
  ReverseEigen.o \
  ReverseSVD.o \
  Patterson.o \
  Gauss.o \
  $(CUBA)


$(LIB): $(LIB)($(OBJS))

$(LIB)(Gauss.o Patterson.o): univariate/config.h

$(LIB)(Eigensystem.o SingularValues.o PseudoEigensystem.o \
  ApplyEigen.o ApplySVD.o): diag/config.h



$(LIB)(Vegas.o): $(VEGAS_C) $(VEGAS_H) $(VEGAS_SRCS)
$(LIB)(vegas-f.o): $(VEGAS_F) $(VEGAS_H)

$(LIB)(Suave.o): $(SUAVE_C) $(SUAVE_H) $(SUAVE_SRCS)
$(LIB)(suave-f.o): $(SUAVE_F) $(SUAVE_H)

$(LIB)(Divonne.o): $(DIVONNE_C) $(DIVONNE_H) $(DIVONNE_SRCS)
$(LIB)(divonne-f.o): $(DIVONNE_F) $(DIVONNE_H)

$(LIB)(Cuhre.o): $(CUHRE_C) $(CUHRE_H) $(CUHRE_SRCS)
$(LIB)(cuhre-f.o): $(CUHRE_F) $(CUHRE_H)


clean:
	$(RM) $(OBJS)

distclean: clean
	$(RM) $(LIB)

