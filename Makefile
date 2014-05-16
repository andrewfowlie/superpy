#########################################################################
#                                                                       #
#     M a k e f i l e
#                                                                       #
#########################################################################


# SuperPy - compile the codes. You might need to tweak individual makefiles
# for your own system. You will need LAPACK, gcc compilers and python and python-dev.

default: all

all: multinest softsusy python micromegas superiso feynhiggs higgsbounds higgssignals susyh nmssm

multinest:
	cd MultiNest-master/build/; cmake .. && make
	export LD_LIBRARY_PATH=./MultiNest-master/lib/:$LD_LIBRARY_PATH

softsusy:
	cd ./softsusy-*; ./configure F77=gfortran # Neccesary on my Ubuntu install.
	make -C softsusy-*

feynhiggs:
	cd ./FeynHiggs-*; ./configure
	make -C FeynHiggs-*
	cd ./FeynHiggs-*; make install
	cd ./FeynHiggs-*; gfortran -o SuperPyFH -Ibuild example/SuperPy.F -Lbuild -lFH

micromegas:
	make -C micromegas_*
	make -C micromegas_*/MSSM main=main.c
	make -C micromegas_*/NMSSM main=main.c

superiso:
	make -C superiso_v*
	make -C superiso_v* slha.c

higgsbounds:
	cd ./HiggsBounds-*; ./configure
	make -C HiggsBounds*

higgssignals:
	cd ./HiggsSignals-*; ./configure
	make -C HiggsSignals-*
	# You may have to alter FHINCLUDE AND FHLIBS by hand in HiggsSignals-*/configure.
	# Paths for HiggsBounds *should* be picked up automatically.
	cd ./HiggsSignals-*; make superpy

susyh:
	make -C susyhit

nmssm:
	make -C NMSSMTools_* init
	make -C NMSSMTools_*

# Build the Python libraries. You might need to sudo these commands.
# Also, if you have your own machine, rather than a networked machine,
# you might want to install things globally, rather than locally, by
# removing the --user argument.

python:
	cd ./pyslha-*; python setup.py install --user
	cd ./PyMultiNest-master; python setup.py install --user

clean:
	make -C MultiNest-master/build clean
	make -C softsusy-* clean
	make -C FeynHiggs-* clean
	make -C micromegas_* clean
	make -C superiso_v* clean
	make -C HiggsBounds-* clean
	make -C HiggsSignals-* clean
	make -C susyhit clean
	make -C NMSSMTools_*
	-rm *~ *.pyc ./SuperPlot/*.pyc ./SuperPlot/*~ ./fastlim-1.0/*pyc *.out


