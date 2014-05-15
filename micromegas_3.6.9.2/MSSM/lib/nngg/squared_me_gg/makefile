LIB = ../nngg.a

OBJS = \
  abbr_s.o \
  abbr_angle.o \
  abbr_hel.o \
  vert.o \
  vert_sqv1.o \
  vert_squ1.o \
  vert_sqe1.o \
  vert_sqd1.o \
  vert_sle1.o \
  vert_sld1.o \
  vert_quu1.o \
  vert_qud1.o \
  vert_lpd1.o \
  vert_cha1.o \
  vert_his1.o \
  vert_his1_sqv1.o \
  vert_his1_squ1.o \
  vert_his1_sqe1.o \
  vert_his1_sqd1.o \
  vert_his1_sle1.o \
  vert_his1_sld1.o \
  vert_his1_quu1.o \
  vert_his1_qud1.o \
  vert_his1_lpd1.o \
  vert_cha1_his1.o \
  box_quu1.o \
  box_qud1.o \
  box_lpd1.o \
  box_cha11.o \
  box_cha12.o \
  box_cha13.o \
  box_cha14.o \
  box_cha15.o \
  box_cha16.o \
  SquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): vars.h 

#renconst.h model.h util.h looptools.h

clean distclean:
	$(RM) $(OBJS) $(LIB)

