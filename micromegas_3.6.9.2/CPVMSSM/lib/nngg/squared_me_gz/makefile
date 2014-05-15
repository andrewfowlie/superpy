LIB = ../nngg.a

OBJS = \
  $(DIR)/gzabbr1_s.o \
  $(DIR)/gzabbr1_angle.o \
  $(DIR)/gzabbr1_hel1.o \
  $(DIR)/gzabbr1_hel2.o \
  $(DIR)/gzvert_neu5.o \
  $(DIR)/gzvert_neu5_quu5.o \
  $(DIR)/gzvert_neu5_qud5.o \
  $(DIR)/gzvert_lpd5_neu5.o \
  $(DIR)/gzvert_cha5_neu51.o \
  $(DIR)/gzvert_cha5_neu52.o \
  $(DIR)/gzvert.o \
  $(DIR)/gzvert_sqv5.o \
  $(DIR)/gzvert_squ5.o \
  $(DIR)/gzvert_sqe5.o \
  $(DIR)/gzvert_sqd5.o \
  $(DIR)/gzvert_sle5.o \
  $(DIR)/gzvert_sld5.o \
  $(DIR)/gzvert_quu5.o \
  $(DIR)/gzvert_qud5.o \
  $(DIR)/gzvert_lpd5.o \
  $(DIR)/gzvert_cha5_cha6.o \
  $(DIR)/gzvert_his5.o \
  $(DIR)/gzvert_his5_sqv5.o \
  $(DIR)/gzvert_his5_squ5.o \
  $(DIR)/gzvert_his5_sqe5.o \
  $(DIR)/gzvert_his5_sqd5.o \
  $(DIR)/gzvert_his5_sle5.o \
  $(DIR)/gzvert_his5_sld5.o \
  $(DIR)/gzvert_his5_quu5.o \
  $(DIR)/gzvert_his5_qud5.o \
  $(DIR)/gzvert_his5_lpd5.o \
  $(DIR)/gzvert_cha5_cha6_his5.o \
  $(DIR)/gzbox_quu51.o \
  $(DIR)/gzbox_quu52.o \
  $(DIR)/gzbox_quu53.o \
  $(DIR)/gzbox_qud51.o \
  $(DIR)/gzbox_qud52.o \
  $(DIR)/gzbox_qud53.o \
  $(DIR)/gzbox_lpd51.o \
  $(DIR)/gzbox_lpd52.o \
  $(DIR)/gzbox_cha51.o \
  $(DIR)/gzbox_cha52.o \
  $(DIR)/gzbox_cha5_cha61.o \
  $(DIR)/gzbox_cha5_cha62.o \
  $(DIR)/gzbox_cha5_cha63.o \
  $(DIR)/gzbox_cha5_cha64.o \
  $(DIR)/gzSquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): vars.h renconst.h model.h util.h looptools.h

clean distclean:
	$(RM) $(OBJS) $(LIB)

