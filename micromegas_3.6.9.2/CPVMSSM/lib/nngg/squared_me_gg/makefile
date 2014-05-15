LIB = ../nngg.a

OBJS = \
  $(DIR)/abbr1_s.o \
  $(DIR)/abbr1_angle.o \
  $(DIR)/abbr1_hel.o \
  $(DIR)/vert.o \
  $(DIR)/vert_sqv5.o \
  $(DIR)/vert_squ5.o \
  $(DIR)/vert_sqe5.o \
  $(DIR)/vert_sqd5.o \
  $(DIR)/vert_sle5.o \
  $(DIR)/vert_sld5.o \
  $(DIR)/vert_quu5.o \
  $(DIR)/vert_qud5.o \
  $(DIR)/vert_lpd5.o \
  $(DIR)/vert_cha5.o \
  $(DIR)/vert_his5.o \
  $(DIR)/vert_his5_sqv5.o \
  $(DIR)/vert_his5_squ5.o \
  $(DIR)/vert_his5_sqe5.o \
  $(DIR)/vert_his5_sqd5.o \
  $(DIR)/vert_his5_sle5.o \
  $(DIR)/vert_his5_sld5.o \
  $(DIR)/vert_his5_quu5.o \
  $(DIR)/vert_his5_qud5.o \
  $(DIR)/vert_his5_lpd5.o \
  $(DIR)/vert_cha5_his5.o \
  $(DIR)/box_quu51.o \
  $(DIR)/box_quu52.o \
  $(DIR)/box_qud51.o \
  $(DIR)/box_qud52.o \
  $(DIR)/box_lpd51.o \
  $(DIR)/box_lpd52.o \
  $(DIR)/box_cha51.o \
  $(DIR)/box_cha52.o \
  $(DIR)/box_cha53.o \
  $(DIR)/SquaredME.o


$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): vars.h renconst.h model.h util.h looptools.h

clean distclean:
	$(RM) $(OBJS) $(LIB)

