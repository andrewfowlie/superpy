default: frontend mma$(ML)

all: frontend mma1


lib: libFH.a CFeynHiggs.h CSLHA.h fcc

frontend: lib FeynHiggs$(EXE) table$(EXE)

mma0:

mma1: lib MFeynHiggs$(EXE)


MAIN = SetFlags.o RetrieveFlags.o \
  SetSMPara.o RetrieveSMPara.o GetSMPara.o \
  SetPara.o RetrievePara.o RetrieveOSPara.o GetPara.o GetTLPara.o \
  SetFV.o RetrieveFV.o GetNMFV.o SetSLHA.o \
  Para.o Sfermions.o SetDebug.o GetSelf.o HiggsCorr.o \
  Uncertainties.o Couplings.o SelectUZ.o Constraints.o Flavour.o \
  HiggsProd.o AddSelf.o CalcSelf.o Error.o DRbartoOS.o AlphaS.o \
  MDcorr.o Deltab2L.o \
  Output.o OutputSLHA.o
#MAIN += gghNLO.o
#MAIN += DRbar2OS_t.o

REC = RecordIndex.o ClearRecord.o \
  ReadRecord.o NextToken.o SLHARecord.o \
  LoopRecord.o SetRecord.o RetrieveRecord.o \
  LoadTable.o TableRecord.o

LT = cache.o ini.o auxB.o auxCD.o A0.o \
  B0.o B1.o B00.o B11.o DB0.o DB1.o DB00.o DB11.o \
  C0.o Cget.o Dget.o \
  A0del.o myG.o TF.o T134.o T134p.o Li2.o Dilog.o CD0.o D0.o

UTIL = Gauss.o ODESolve.o EigenSf.o SfEigensystem.o \
  HEigensystem.o SEigensystem.o TakagiFactor.o SVD.o LU.o

SLHA = SLHAClear.o SLHARead.o SLHAWrite.o \
  SLHAGetInfo.o SLHAPutInfo.o SLHACopyString.o \
  SLHANewDecay.o SLHAFindDecay.o SLHAAddDecay.o SLHAGetDecay.o \
  SLHADecayTable.o SLHAValid.o SLHAExist.o SLHAPDGName.o

ONELOOP = se_h0h0.o se_HHHH.o se_A0A0.o se_HmHp.o \
  se_h0HH.o se_h0A0.o se_HHA0.o \
  se_G0G0.o se_h0G0.o se_HHG0.o se_A0G0.o \
  se_GmGp.o se_HmGp.o \
  se_A0Z.o se_ZZ.o se_WW.o td_h0.o td_HH.o td_A0.o \
  dse_h0h0.o dse_HHHH.o dse_A0A0.o dse_HmHp.o \
  dse_h0HH.o dse_h0A0.o dse_HHA0.o \
  dse_G0G0.o dse_h0G0.o dse_HHG0.o dse_A0G0.o \
  dse_GmGp.o dse_HmGp.o \
  dse_ZZ.o dse_WW.o

TLPS = TLShifts.o TLasat.o TLasab.o TLatab.o TLatat.o \
  Pi11t.o Pi12t.o Pi22t.o Pi11b.o Pi12b.o Pi22b.o Pitt.o \
  Resum4H.o

TLHR = TLhr_abbr0.o TLhr_abbr1_1.o TLhr_abbr1_2.o TLhr_abbr1_3.o \
  TLhr_abbr1_4.o TLhr_abbr1_5.o TLhr_abbr1_6.o TLhr_abbr1_7.o \
  TLhr_abbr1_8.o TLhr_abbr1_9.o TLhr_abbr1_10.o TLhr_abbr1_11.o \
  TLhr_abbr1_12.o TLhr_abbr1_13.o TLhr_abbr1_14.o TLhr_abbr1_15.o \
  TLhr_abbr1_16.o TLhr_abbr1_17.o TLhr_abbr1_18.o \
  TLhr_abbr2_1.o TLhr_abbr2_2.o TLhr_abbr2_3.o \
  TLhr_abbr3_1.o TLhr_abbr3_2.o TLhr_abbr3_3.o TLhr_abbr3_4.o \
  TLhr_rcs.o TLhr_self1.o TLhr_self2.o TLhr_self3.o \
  TLhr.o MTrepara.o

TWOLOOP = $(TLPS) $(TLHR) InvTri.o Phi.o

GM2_2L_0 = gm2_2L_1.o gm2_2L_2.o

GM2_2L_FULL = gm2_2L_full1.o gm2_2L_full2.o gm2_2L_full3.o \
  gm2_2L_full4.o gm2_2L_full5.o gm2_2L_full6.o gm2_2L_full7.o \
  gm2_2L_full8.o gm2_2L_full9.o gm2_2L_full10.o gm2_2L_full11.o

GM2_2L_1 = $(GM2_2L_0) $(GM2_2L_FULL)

GM2 = gm2_1L_1.o $(GM2_2L_$(GM2FULL))

CONSTRAINTS = $(GM2) EWPO.o DeltaRho_1L.o DeltaRho_2L.o EDM.o \
  Bsgamma.o DeltaMs.o Bsll.o

Bsgamma.o DeltaMs.o Bsll.o: debug.h

POBS =
#POBS += drho.o

HGAGA = hgaga-SquaredME.o \
  hgaga-abbr1s.o hgaga-abbr1h.o hgaga-v1.o \
  hgaga-v1_Cha4.o hgaga-v1_Gen4.o \
  hgaga-v1_All4_Ind1_Ind2.o

HGAGASM = hgagaSM-SquaredME.o \
  hgagaSM-abbr1s.o hgagaSM-abbr1h.o \
  hgagaSM-v1.o hgagaSM-v1_Gen4.o

HGAZ = hgaZ-SquaredME.o hgaZ-CalcRenConst.o \
  hgaZ-abbr1s.o hgaZ-abbr1h.o \
  hgaZ-v1.o hgaZ-v1_Gen4.o \
  hgaZ-v1_All4_All5_Ind1_Ind2_Ind3.o \
  hgaZ-v1_Cha4_Cha5.o

HGAZSM = hgaZSM-SquaredME.o hgaZSM-CalcRenConst.o \
  hgaZSM-abbr1s.o hgaZSM-abbr1h.o \
  hgaZSM-v1.o hgaZSM-v1_Gen4.o

HGG = hgg-SquaredME.o \
  hgg-abbr1s.o hgg-abbr1h.o \
  hgg-v1_Gen4.o hgg-v1_All4_Ind1_Ind2.o

HGGSM = hggSM-SquaredME.o \
  hggSM-abbr1s.o hggSM-abbr1h.o \
  hggSM-v1_Gen4.o

HEE = hee-SquaredME.o hee-CalcRenConst.o \
  hee-abbr0s.o hee-abbr0h.o \
  hee-abbr1s.o hee-abbr1h.o \
  hee-v0.o hee-v1.o \
  hee-v1_All4_All5_Ind1_Ind2_Ind3.o \
  hee-v1_All4_All5_Ind1_Ind2_Ind3_Ind4.o \
  hee-v1_All4_All5_Ind1_Ind2_Neu4.o \
  hee-v1_All4_Ind1.o \
  hee-v1_All4_Neu4_Neu5.o \
  hee-v1_Cha4_Cha5.o \
  hee-v1_Cha4_Cha5_Gen4.o \
  hee-v1_Cha4_Gen4.o \
  hee-v1.o \
  hee-v1_Gen4.o \
  hee-v1_Neu4_Neu5.o

HEESM = heeSM-SquaredME.o heeSM-CalcRenConst.o \
  heeSM-abbr0h.o \
  heeSM-abbr1s.o heeSM-abbr1h.o \
  heeSM-v0.o heeSM-v1.o \

HUU = huu-SquaredME.o huu-CalcRenConst.o \
  huu-abbr0s.o huu-abbr0h.o \
  huu-abbr1s.o huu-abbr1h.o \
  huu-v0.o huu-v1.o \
  huu-v1_All4_All5_Cha4_Ind1_Ind2_Ind3_Ind4.o \
  huu-v1_All4_All5_Ind1_Ind2.o \
  huu-v1_All4_All5_Ind1_Ind2_Ind3.o \
  huu-v1_All4_All5_Ind1_Ind2_Ind3_Ind4.o \
  huu-v1_All4_All5_Ind1_Ind2_Neu4.o \
  huu-v1_All4_Cha4_Cha5_Ind1_Ind2.o \
  huu-v1_All4_Ind1.o \
  huu-v1_All4_Neu4_Neu5.o \
  huu-v1_Cha4_Cha5.o \
  huu-v1_Gen4.o \
  huu-v1_Neu4_Neu5.o

HUUSM = huuSM-SquaredME.o huuSM-CalcRenConst.o \
  huuSM-abbr0s.o huuSM-abbr0h.o \
  huuSM-abbr1s.o huuSM-abbr1h.o \
  huuSM-v0.o huuSM-v1.o huuSM-v1_Gen4.o

HDD = hdd-SquaredME.o hdd-CalcRenConst.o \
  hdd-abbr0s.o hdd-abbr0h.o \
  hdd-abbr1s.o hdd-abbr1h.o \
  hdd-v0.o hdd-v1.o \
  hdd-v1_All4_All5_Cha4_Ind1_Ind2_Ind3_Ind4.o \
  hdd-v1_All4_All5_Ind1_Ind2.o \
  hdd-v1_All4_All5_Ind1_Ind2_Ind3.o \
  hdd-v1_All4_All5_Ind1_Ind2_Ind3_Ind4.o \
  hdd-v1_All4_All5_Ind1_Ind2_Neu4.o \
  hdd-v1_All4_Cha4_Cha5_Ind1_Ind2.o \
  hdd-v1_All4_Ind1.o \
  hdd-v1_All4_Neu4_Neu5.o \
  hdd-v1_Cha4_Cha5.o \
  hdd-v1_Gen4.o \
  hdd-v1_Neu4_Neu5.o

HDDSM = hddSM-SquaredME.o hddSM-CalcRenConst.o \
  hddSM-abbr0s.o hddSM-abbr0h.o \
  hddSM-abbr1s.o hddSM-abbr1h.o \
  hddSM-v0.o hddSM-v1.o hddSM-v1_Gen4.o

fvHHH = hhh-CalcRenConst.o hhh-SquaredME.o \
  hhh-abbr0s.o hhh-abbr1s1.o hhh-abbr1s2.o \
  hhh-v0.o \
  hhh-v1_All4_All5_All6_Ind1_Ind2_Ind3_Ind4_Ind5_Ind6.o \
  hhh-v1_All4_All5_Ind1_Ind2_Ind3.o hhh-v1_All4_All5_Ind1_Ind2_Ind3_Ind4.o \
  hhh-v1_All4_Ind1.o \
  hhh-v1_Cha4_Cha5_Cha6.o hhh-v1_Cha4_Cha5.o hhh-v1.o hhh-v1_Gen4.o \
  hhh-v1_Gen4_Sfe4.o \
  hhh-v1_Gen4_Sfe4_Sfe5.o \
  hhh-v1_Gen4_Sfe4_Sfe5_Sfe6.o \
  hhh-v1_Neu4_Neu5.o \
  hhh-v1_Neu4_Neu5_Neu6.o

HHH = hhh-CalcRenConst.o hhh-SquaredME.o \
  hhh-abbr0s.o hhh-abbr1s1.o hhh-abbr1s2.o \
  hhh-v0.o hhh-v1.o hhh-v1_Cha4_Cha5.o \
  hhh-v1_Cha4_Cha5_Cha6.o hhh-v1_Gen4.o \
  hhh-v1_Gen4_Sfe4.o hhh-v1_Gen4_Sfe4_Sfe5.o \
  hhh-v1_Gen4_Sfe4_Sfe5_Sfe6.o hhh-v1_Neu4_Neu5.o \
  hhh-v1_Neu4_Neu5_Neu6.o

DECAYS = VecSet.o \
  $(HGAGA) $(HGAGASM) \
  $(HGAZ) $(HGAZSM) \
  $(HGG) $(HGGSM) \
  $(HEE) $(HEESM) \
  $(HUU) $(HUUSM) \
  $(HDD) $(HDDSM) \
  $(HHH)

MOST = $(MAIN) $(ONELOOP) $(TWOLOOP) $(DECAYS) $(CONSTRAINTS) $(POBS)

ALL = $(MOST) $(REC) $(LT) $(UTIL) $(SLHA)


libFH.a: $(ALL)
	ar cru libFH.a $?
	-ranlib libFH.a


$(ALL): externals.h types.h

$(MOST) LoopRecord.o ReadRecord.o TableRecord.o: FH.h const.h

$(REC) SetSMPara.o Para.o SetSLHA.o SetNMFV.o: FHRecord.h

ReadRecord.o LoopRecord.o SetRecord.o LoadTable.o TableRecord.o: debug.h

LoopRecord.o SetRecord.o: Defaults.h

SetSLHA.o: SLHA.h SLHADefs.h

$(SLHA): SLHADefs.h

SLHARead.o SLHAWrite.o SLHAPDGName.o: PDG.h

SLHARead.o: SLHAReadBlocks.h

SLHAWrite.o: SLHAWriteBlocks.h

SLHARead.o SLHAWrite.o SLHANewDecay.o SLHAFindDecay.o \
  SLHAAddDecay.o SLHAGetDecay.o SLHADecayTable.o: SLHAIntern.h

SLHARecord.o OutputSLHA.o: DeltaAlfaTop.h

SetFlags.o OutputSLHA.o: version.h

SetFlags.o RetrieveFlags.o SetSMPara.o RetrieveSMPara.o \
  RetrievePara.o SetSLHA.o SetCKM.o RetrieveCKM.o SetNMFV.o \
  RetrieveNMFV.o Para.o GetPara.o GetSelf.o CalcSelf.o \
  HiggsCorr.o Uncertainties.o Couplings.o SelectUZ.o \
  HiggsProd.o Constraints.o DRbartoOS.o Error.o: debug.h

CalcSelf.o: OneLoop.h

Couplings.o HiggsProd.o: FHCouplings.h

Couplings.o: Dilog.o NHiggsDecayFitsLO.f NHiggsDecayFitsNLO.f

HiggsProd.o: gghcTev.f gghkTev.f NHiggsProdFits-Tev.f NHiggsProdFits-LHC.f \
  gghcLHC.f gghkLHC.f NHiggsProdFits-LHC7.f NHiggsProdFits-LHC8.f \
  NHiggsProdFits-LHC14.f CHiggsProdFits-LHC8.f CHiggsProdFits-LHC14.f

Para.o DRbartoOS.o: looptools.h ltcache.h

$(LT): lt.h

Cget.o: ltcache.h

HEigensystem.o SEigensystem.o TakagiFactor.o SVD.o Gauss.o LU.o: Util.h

$(ONELOOP): looptools.h ltcache.h debug.h

$(TLHR): TLhr.h TLhrvars.h FH.h looptools.h ltcache.h debug.h

MTrepara.o: MTrepara.h

$(TLPS): TLps.h TLpsvars.h looptools.h ltcache.h

Constraints.o $(GM2): gm2.h

$(GM2_2L_0): gm2_2L.h

$(GM2_2L_FULL): gm2_2L_full.h

DeltaRho_1L.o $(GM2): looptools.h ltcache.h

dr.o drho.o: pobs.h

Bsgamma.o: bsgamma.h bsglu.h looptools.h ltcache.h

DeltaMs.o: deltaMs.h looptools.h ltcache.h

Bsll.o: bsll.h looptools.h ltcache.h

$(HGAGA): hgaga-vars.h

$(HGAGASM): hgagaSM-vars.h

$(HGAZ): hgaZ-vars.h hgaZ-renconst.h

$(HGAZSM): hgaZSM-vars.h hgaZSM-renconst.h

$(HGG): hgg-vars.h

$(HGGSM): hggSM-vars.h

$(HEE): hee-vars.h hee-renconst.h

$(HEESM): heeSM-vars.h heeSM-renconst.h

$(HUU): huu-vars.h huu-renconst.h

$(HUUSM): huuSM-vars.h huuSM-renconst.h

$(HDD): hdd-vars.h hdd-renconst.h

$(HDDSM): hddSM-vars.h hddSM-renconst.h

$(HHH): hhh-vars.h hhh-renconst.h

$(DECAYS) Couplings.o: Decay.h inline.h looptools.h ltcache.h debug.h

VecSet.o: Decay.h


FHRecord.h: FHRecord.h.in RecordIndices.h
	cat RecordIndices.h FHRecord.h.in > FHRecord.h

CFeynHiggs.h: CFeynHiggs.h.in RecordIndices.h ftypes.h
	sed "s:NOUNDERSCORE:$(NOUNDERSCORE):" ftypes.h RecordIndices.h CFeynHiggs.h.in > CFeynHiggs.h

CSLHA.h: CSLHA.h.in ftypes.h
	sed "s:NOUNDERSCORE:$(NOUNDERSCORE):" ftypes.h CSLHA.h.in > CSLHA.h


FE = FeynHiggs.o pager.o

FeynHiggs$(EXE): $(FE) libFH.a
	$(FC) $(FFLAGS) -o FeynHiggs$(EXE) $(FE) libFH.a

MFeynHiggs$(EXE): MFeynHiggs$(QUAD).tm CFeynHiggs.h CSLHA.h fcc fortranflush.o libFH.a
	NM="$(NM)" DLLTOOL="$(DLLTOOL)" \
	CC="./fcc" REALCC="$(CC) $(CFLAGS)" \
	CXX="./f++" REALCXX="$(CXX) $(CXXFLAGS)" \
	PATH="$$PATH:." \
	  mcc -o MFeynHiggs$(EXE) $(MCFLAGS) -DQUAD=0$(QUAD) MFeynHiggs$(QUAD).tm fortranflush.o libFH.a
	rm -f MFeynHiggs$(QUAD).tm.*

MFeynHiggs1.tm: MFeynHiggs.tm
	sed '/:Begin:/,/:End:/ s/Real/Real128/g' MFeynHiggs.tm > MFeynHiggs1.tm


FeynHiggs.o: version.h FHRecord.h SLHA.h


fcc: fcc.in makefile
	sed "s:^fldflags=.*:fldflags=\"$(LDFLAGS)\":" fcc.in > fcc
	chmod 755 fcc
	rm -f f++
	ln -s fcc f++

table$(EXE): table.c
	$(CC) $(CFLAGS) -o table$(EXE) table.c

