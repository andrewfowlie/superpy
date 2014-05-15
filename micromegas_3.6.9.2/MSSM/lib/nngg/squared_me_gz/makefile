PREFIX = AZ

FFLAGS += -DPREFIX=$(PREFIX)

LIB = ../nngg.a

OBJS = \
  gzabbr_s.o \
  gzabbr_angle.o \
  gzabbr_hel1.o \
  gzabbr_hel2.o \
  gzvert_neu1.o \
  gzvert_neu1_quu1.o \
  gzvert_neu1_qud1.o \
  gzvert_lpd1_neu1.o \
  gzvert_cha1_neu11.o \
  gzvert_cha1_neu12.o \
  gzvert.o \
  gzvert_sqv1.o \
  gzvert_squ1.o \
  gzvert_sqe1.o \
  gzvert_sqd1.o \
  gzvert_sle1.o \
  gzvert_sld1.o \
  gzvert_quu1.o \
  gzvert_qud1.o \
  gzvert_lpd1.o \
  gzvert_cha1_cha2.o \
  gzvert_his1.o \
  gzvert_his1_sqv1.o \
  gzvert_his1_squ1.o \
  gzvert_his1_sqe1.o \
  gzvert_his1_sqd1.o \
  gzvert_his1_sle1.o \
  gzvert_his1_sld1.o \
  gzvert_his1_quu1.o \
  gzvert_his1_qud1.o \
  gzvert_his1_lpd1.o \
  gzvert_cha1_cha2_his1.o \
  gzbox_quu11.o \
  gzbox_quu12.o \
  gzbox_qud11.o \
  gzbox_qud12.o \
  gzbox_lpd11.o \
  gzbox_lpd12.o \
  gzbox_cha11.o \
  gzbox_cha12.o \
  gzbox_cha13.o \
  gzbox_cha14.o \
  gzbox_cha15.o \
  gzbox_cha16.o \
  gzbox_cha17.o \
  gzbox_cha1_cha21.o \
  gzbox_cha1_cha22.o \
  gzbox_cha1_cha23.o \
  gzbox_cha1_cha24.o \
  gzbox_cha1_cha25.o \
  gzSquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): vars.h renconst.h model.h util.h looptools.h

clean distclean:
	$(RM) $(OBJS) $(LIB)

