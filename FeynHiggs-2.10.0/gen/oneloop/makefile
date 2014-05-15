all: mfv nmfv
	./combine

SE = \
  se_h0h0.F dse_h0h0.F \
  se_HHHH.F dse_HHHH.F \
  se_A0A0.F dse_A0A0.F \
  se_h0HH.F dse_h0HH.F \
  se_h0A0.F dse_h0A0.F \
  se_HHA0.F dse_HHA0.F \
  se_HmHp.F dse_HmHp.F \
  se_G0G0.F dse_G0G0.F \
  se_h0G0.F dse_h0G0.F \
  se_HHG0.F dse_HHG0.F \
  se_A0G0.F dse_A0G0.F \
  se_HmGp.F dse_HmGp.F \
  se_GmGp.F dse_GmGp.F \
  se_WW.F   dse_WW.F \
  se_ZZ.F   dse_ZZ.F \
  se_A0Z.F \
  td_h0.F \
  td_HH.F \
  td_A0.F

.SUFFIXES:

mfv: $(addprefix f-mfv/, $(SE))

nmfv: $(addprefix f-nmfv/, $(SE))

%.F:
	./se $@


clean:
	rm -fr f-mfv f-nmfv f

