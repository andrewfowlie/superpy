LIB = ../nngg.a

OBJS = \
  CalcRenConst.o \
  ZeroRenConst.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): renconst.h model.h util.h looptools.h

clean distclean:
	$(RM) $(OBJS) $(LIB)

