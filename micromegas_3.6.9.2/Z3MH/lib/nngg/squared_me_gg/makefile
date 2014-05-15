LIB = ../nngg.a

OBJS = \
  $(DIR)/abbr1_s.o \
  $(DIR)/abbr1_angle.o \
  $(DIR)/abbr1_hel.o \
  $(DIR)/vert.o \
  $(DIR)/box.o \
  $(DIR)/SquaredME.o


$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): vars.h renconst.h model.h util.h looptools.h

clean distclean:
	$(RM) $(OBJS) $(LIB)

