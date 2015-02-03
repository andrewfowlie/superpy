#all: sprog
all: sfile

FOPT = -O0 -g -fno-automatic -fbounds-check -Wall -fbacktrace
#FOPT = -O -fno-automatic -Wall
F77 = gfortran-4.8
COMP = $(F77) $(FOPT)

obj/rombint.o: lib/rombint.f
	$(COMP) -c lib/rombint.f
	mv rombint.o obj
obj/bsg_nl.o: lib/bsg_nl.f
	$(COMP) -c lib/bsg_nl.f
	mv bsg_nl.o obj
obj/ddg_fun.o: lib/ddg_fun.f
	$(COMP) -c lib/ddg_fun.f
	mv ddg_fun.o obj
obj/dd_gamma.o: lib/dd_gamma.f
	$(COMP) -c lib/dd_gamma.f
	mv dd_gamma.o obj
obj/dd_gluon.o: lib/dd_gluon.f
	$(COMP) -c lib/dd_gluon.f
	mv dd_gluon.o obj
obj/dd_vv.o: lib/dd_vv.f
	$(COMP) -c lib/dd_vv.f
	mv dd_vv.o obj
obj/dd_ll.o: lib/dd_ll.f
	$(COMP) -c lib/dd_ll.f
	mv dd_ll.o obj
obj/dd_mix.o: lib/dd_mix.f
	$(COMP) -c lib/dd_mix.f
	mv dd_mix.o obj
obj/uu_mix.o: lib/uu_mix.f
	$(COMP) -c lib/uu_mix.f
	mv uu_mix.o obj
#obj/sff_fun0.o: lib/sff_fun0.f
#	$(COMP) -c lib/sff_fun0.f
#	mv sff_fun0.o obj
#obj/sdd_vert0.o: lib/sdd_vert0.f
#	$(COMP) -c lib/sdd_vert0.f
#	mv sdd_vert0.o obj
#obj/pdd_vert0.o: lib/pdd_vert0.f
#	$(COMP) -c lib/pdd_vert0.f
#	mv pdd_vert0.o obj
obj/zdd_vert0.o: lib/zdd_vert0.f
	$(COMP) -c lib/zdd_vert0.f
	mv zdd_vert0.o obj
obj/d_self0.o: lib/d_self0.f
	$(COMP) -c lib/d_self0.f
	mv d_self0.o obj
#obj/puu_vert0.o: lib/puu_vert0.f
#	$(COMP) -c lib/puu_vert0.f
#	mv puu_vert0.o obj
#obj/suu_vert0.o: lib/suu_vert0.f
#	$(COMP) -c lib/suu_vert0.f
#	mv suu_vert0.o obj
#obj/u_self0.o: lib/u_self0.f
#	$(COMP) -c lib/u_self0.f
#	mv u_self0.o obj
obj/phen_4q.o: lib/phen_4q.f
	$(COMP) -c lib/phen_4q.f
	mv phen_4q.o obj
obj/phen_2q.o: lib/phen_2q.f
	$(COMP) -c lib/phen_2q.f
	mv phen_2q.o obj
obj/mh_diag.o: lib/mh_diag.f
	$(COMP) -c lib/mh_diag.f
	mv mh_diag.o obj
#
obj/cd_fun.o: lib/cd_fun.f
	$(COMP) -c lib/cd_fun.f
	mv cd_fun.o obj
obj/b_fun.o: lib/b_fun.f
	$(COMP) -c lib/b_fun.f
	mv b_fun.o obj
obj/db_fun.o: lib/db_fun.f
	$(COMP) -c lib/db_fun.f
	mv db_fun.o obj
obj/c_fun.o: lib/c_fun.f
	$(COMP) -c lib/c_fun.f
	mv c_fun.o obj
obj/c_fun_exp.o: lib/c_fun_exp.f
	$(COMP) -c lib/c_fun_exp.f
	mv c_fun_exp.o obj
obj/qcd_fun.o: lib/qcd_fun.f
	$(COMP) -c lib/qcd_fun.f
	mv qcd_fun.o obj
obj/eisch1.o: lib/eisch1.f
	$(COMP) -c lib/eisch1.f
	mv eisch1.o obj
obj/vg_def.o: lib/vg_def.f
	$(COMP) -c lib/vg_def.f
	mv vg_def.o obj
obj/vf_def.o: lib/vf_def.f
	$(COMP) -c lib/vf_def.f
	mv vf_def.o obj
obj/mh_init.o: lib/mh_init.f
	$(COMP) -c lib/mh_init.f
	mv mh_init.o obj
#
obj/cdm_q.o: lib/cdm_q.f
	$(COMP) -c lib/cdm_q.f
	mv cdm_q.o obj
obj/edm_q.o: lib/edm_q.f 
	$(COMP) -c lib/edm_q.f
	mv edm_q.o obj
obj/vegas.o: lib/vegas.f
	$(COMP) -c lib/vegas.f
	mv vegas.o obj
obj/cdm_g.o: lib/cdm_g.f
	$(COMP) -c lib/cdm_g.f
	mv cdm_g.o obj
obj/sflav_io.o: lib/sflav_io.f
	$(COMP) -c lib/sflav_io.f
	mv sflav_io.o obj
obj/yuk_ren.o: lib/yuk_ren.f
	$(COMP) -c lib/yuk_ren.f
	mv yuk_ren.o obj
obj/q_self0_dlim.o: lib/q_self0_dlim.f
	$(COMP) -c lib/q_self0_dlim.f
	mv q_self0_dlim.o obj
obj/l_self0_dlim.o: lib/l_self0_dlim.f
	$(COMP) -c lib/l_self0_dlim.f
	mv l_self0_dlim.o obj
obj/ll_gamma.o: lib/ll_gamma.f
	$(COMP) -c lib/ll_gamma.f
	mv ll_gamma.o obj
obj/phen_2l.o: lib/phen_2l.f
	$(COMP) -c lib/phen_2l.f
	mv phen_2l.o obj
obj/sflav_main.o: lib/sflav_main.f
	$(COMP) -c lib/sflav_main.f
	mv sflav_main.o obj

obj/suu_vert.o: lib/suu_vert.f
	$(COMP) -c lib/suu_vert.f
	mv suu_vert.o obj
obj/u_self.o: lib/u_self.f
	$(COMP) -c lib/u_self.f
	mv u_self.o obj
obj/vh_def.o: lib/vh_def.f
	$(COMP) -c lib/vh_def.f
	mv vh_def.o obj
obj/sff_fun.o: lib/sff_fun.f
	$(COMP) -c lib/sff_fun.f
	mv sff_fun.o obj

obj/uu_gluon.o: lib/uu_gluon.f
	$(COMP) -c lib/uu_gluon.f
	mv uu_gluon.o obj
obj/vff_fun.o: lib/vff_fun.f
	$(COMP) -c lib/vff_fun.f
	mv vff_fun.o obj


libfcnc.a: obj/bsg_nl.o obj/ddg_fun.o obj/dd_gamma.o obj/dd_gluon.o	\
	obj/dd_vv.o obj/dd_ll.o obj/dd_mix.o obj/uu_mix.o		\
	obj/zdd_vert0.o obj/d_self0.o obj/phen_4q.o obj/phen_2q.o	\
	obj/mh_diag.o obj/cd_fun.o obj/b_fun.o obj/db_fun.o		\
	obj/c_fun.o obj/qcd_fun.o obj/vg_def.o obj/vf_def.o		\
	obj/mh_init.o obj/eisch1.o obj/rombint.o obj/cdm_q.o		\
	obj/cdm_g.o obj/vegas.o obj/edm_q.o obj/sflav_io.o		\
	obj/yuk_ren.o obj/q_self0_dlim.o obj/sflav_main.o		\
	obj/l_self0_dlim.o obj/ll_gamma.o obj/phen_2l.o obj/u_self.o	\
	obj/sff_fun.o obj/suu_vert.o obj/vh_def.o obj/c_fun_exp.o	\
	obj/vff_fun.o obj/uu_gluon.o

	ar crs libfcnc.a obj/bsg_nl.o obj/ddg_fun.o obj/dd_gamma.o	\
	obj/dd_gluon.o obj/dd_vv.o obj/dd_ll.o obj/dd_mix.o		\
	obj/uu_mix.o obj/zdd_vert0.o obj/d_self0.o obj/phen_4q.o	\
	obj/phen_2q.o obj/mh_diag.o obj/cd_fun.o obj/b_fun.o		\
	obj/db_fun.o obj/c_fun.o obj/qcd_fun.o obj/vg_def.o		\
	obj/vf_def.o obj/mh_init.o obj/eisch1.o obj/rombint.o		\
	obj/cdm_q.o obj/cdm_g.o obj/vegas.o obj/edm_q.o			\
	obj/sflav_io.o obj/yuk_ren.o obj/q_self0_dlim.o			\
	obj/sflav_main.o obj/l_self0_dlim.o obj/ll_gamma.o		\
	obj/phen_2l.o obj/u_self.o obj/sff_fun.o obj/suu_vert.o		\
	obj/vh_def.o obj/c_fun_exp.o obj/vff_fun.o obj/uu_gluon.o


#
#	Useful but cause problems on BSD systems:
#	ar ts libfcnc.a
#

sscan:	susy_flavor_scan.f mssm_constraints.f libfcnc.a 
	$(COMP) -o sscan susy_flavor_scan.f  -L. -lfcnc
	@echo "Executing TEST..."
	./sscan

sfile:	susy_flavor_file.f libfcnc.a 
	$(COMP) -o sfile susy_flavor_file.f  -L. -lfcnc
	@echo "Executing SUSY_FLAVOR..."
	./sfile
#	cat susy_flavor.out

sprog:	susy_flavor_prog.f libfcnc.a 
	$(COMP) -o sprog susy_flavor_prog.f  -L. -lfcnc
	@echo "Executing SUSY_FLAVOR..."
	./sprog

clean:
	rm -rf *.o obj/*.o *.a *~ lib/*~ sprog sfile sscan

tar:	
	rm -rf obj/*
	tar zcvf susy_flavor.tgz susy_flavor*  obj lib/*.f Makefile
	mv susy_flavor.tgz $(HOME)

