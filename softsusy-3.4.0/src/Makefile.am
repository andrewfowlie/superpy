ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

lib_LTLIBRARIES = libsoft.la 

bin_PROGRAMS = $(top_builddir)/softsusy.x \
               $(top_builddir)/softsusy-nmssm.x \
               $(top_builddir)/softpoint.x \
               $(top_builddir)/rpvsoftsusy.x \
               $(top_builddir)/rpvneut.x 

__top_builddir__softsusy_x_SOURCES = main.cpp
__top_builddir__softsusy_x_LDADD = libsoft.la  

__top_builddir__softsusy_nmssm_x_SOURCES = main-nmssm.cpp
__top_builddir__softsusy_nmssm_x_LDADD = libsoft.la

__top_builddir__softpoint_x_SOURCES = softpoint.cpp 
__top_builddir__softpoint_x_LDADD = libsoft.la

__top_builddir__rpvsoftsusy_x_SOURCES = rpvmain.cpp
__top_builddir__rpvsoftsusy_x_LDADD = libsoft.la 

__top_builddir__rpvneut_x_SOURCES = rpvneutmain.cpp
__top_builddir__rpvneut_x_LDADD = libsoft.la 

pkginclude_HEADERS = \
def.h main.h numerics.h \
softpars.h softpars.cpp susy.h linalg.h \
physpars.h softpoint.h twoloophiggs.h \
nmssm2loop.h dilogwrap.h \
lowe.h mycomplex.h rge.h softsusy.h softsusy.cpp \
utils.h mssmUtils.h rpvsoft.h rpvsusypars.h \
tensor.h rpvmain.h flavoursoft.h \
xpr-base.h xpr-matrix.h xpr-vector.h rpvneut.h \
nmssmsusy.h nmssmsoftpars.h  \
nmssmsoftsusy.h nmssmUtils.h

libsoft_la_SOURCES = def.cpp linalg.cpp lowe.cpp numerics.cpp  rpvneut.cpp \
physpars.cpp rge.cpp nmssm2loop.f dilogwrap.f \
 nmssmsoftsusy.cpp nmssmsusy.cpp nmssmsoftpars.cpp nmssmUtils.cpp \
susy.cpp utils.cpp mssmUtils.cpp \
tensor.cpp rpvsusypars.cpp rpvsoft.cpp flavoursoft.cpp twoloophiggs.f 

libsoft_la_LDFLAGS = @FLIBS@

EXTRA_DIST = nmProcessSpec.f
