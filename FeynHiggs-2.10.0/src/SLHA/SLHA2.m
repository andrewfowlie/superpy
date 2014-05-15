Needs["GenSLHA`"]

Unprotect[D];
Remove[D]

MAXMSG = 15;
MSGLEN = 5;
CPLXSIZE = 16

desc = {
  ModSel[modse] -> {
    Model,
    GridPts,
    Content,
    RPV,
    CPV,
    FV
  },
  SMInputs[sminp] -> {
    invAlfaMZ,
    GF,
    AlfasMZ,
    MZ,
    Mf[{t,1,4},{g,3}],
      Mnu1 == Mf[1,1],
      Me == Mf[2,1],
      Mu == Mf[3,1],
      Md == Mf[4,1],
      Mnu2 == Mf[1,2],
      Mmu == Mf[2,2],
      Mc == Mf[3,2],
      Ms == Mf[4,2],
      Mnu3 == Mf[1,3],
      Mtau == Mf[2,3],
      Mt == Mf[3,3],
      Mb == Mf[4,3]
  },
  MinPar[minpa] -> {
    M0,
      Lambda == M0,
    M12,
      Mmess == M12,
      M32 == M12,
    TB,
    signMUE,
    A,
      N5 == A,
    cgrav
  },
  ExtPar[extpa] -> {
    Q,
    M1,
    M2,
    M3,
    Af[{t,2,4}],
      Atau == Af[2],
      At == Af[3],
      Ab == Af[4],
    MHu2,
    MHd2,
    MUE,
    MA02,
    TB,
    MA0,
    MHp,
    Struct[{q,5}][
      MSS[{g,3}]
    ],
      MSL[g] == MSS[g,1],
      MSE[g] == MSS[g,2],
      MSQ[g] == MSS[g,3],
      MSU[g] == MSS[g,4],
      MSD[g] == MSS[g,5],
    N5[{g,3}],
    lambda,
    kappa,
    Alambda,
    Akappa,
    lambdaS,
    xiF,
    xiS,
    MUEprime,
    mS2prime,
    mS2
  },
  QExtPar[qxtpa] -> {
    QM1,
    QM2,
    QM3,
    QAf[{t,2,4}],
      QAtau == QAf[2],
      QAt == QAf[3],
      QAb == QAf[4],
    QMHu2,
    QMHd2,
    QMUE,
    QMA02,
    QTB,
    QMSS[{q,5}],
      QMSL == QMSS[1],
      QMSE == QMSS[2],
      QMSQ == QMSS[3],
      QMSU == QMSS[4],
      QMSD == QMSS[5]
  },
  NMSSMRun[nmrun] -> {
    Q,
    lambda,
    kappa,
    Alambda,
    Akappa,
    lambdaS,
    xiF,
    xiS,
    MUEprime,
    mS2prime,
    mS2
  },
  Mass[mass] -> {
    Mf[{t,4},{g,3}],
    MSf[{s,2},{t,4},{g,3}],
    MZ,
    MW,
    Mh0,
    MHH,
    MA0,
    MHp,
      MH1 == Mh0,
      MH2 == MHH,
    MH3,
      MA1 == MA0,
    MA2,
    MNeu[{n,5}],
    MCha[{c,2}],
    MGl,
    MGrav
  },
  DMass[tmass] -> {
    Q,
    DeltaMh0,
    DeltaMHH,
    DeltaMA0,
    DeltaMHp
  },
  NMix[nmix] -> {
    ZNeu[{n1,4},{n2,4}]
  },
  UMix[umix] -> {
    UCha[{c1,2},{c2,2}]
  },
  VMix[vmix] -> {
    VCha[{c1,2},{c2,2}]
  },
  SfMix :> {
    Struct[{t,2,4}][
      USf[{s1,2},{s2,2}]
    ]
  },
  StauMix[staum] -> {
    USf[{s1,2},{s2,2}] -> SfMix**USf[s1,s2,2]
  },
  StopMix[stopm] -> {
    USf[{s1,2},{s2,2}] -> SfMix**USf[s1,s2,3]
  },
  SbotMix[sbotm] -> {
    USf[{s1,2},{s2,2}] -> SfMix**USf[s1,s2,4]
  },
  Alpha[alfa] -> {
    Alpha
  },
  DAlpha[talfa] -> {
    DeltaAlpha
  },
  HMix[hmix] -> {
    Q,
    MUE,
    TB,
    VEV,
    MA02
  },
  Gauge[gauge] -> {
    Q,
    g1,
    g2,
    g3
  },
  MSoft[msoft] -> {
    Q,
    M1,
    M2,
    M3,
    MHu2,
    MHd2,
    Struct[{q,5}][
      MSS[{g,3}]
    ],
      MSL[g] == MSS[g,1],
      MSE[g] == MSS[g,2],
      MSQ[g] == MSS[g,3],
      MSU[g] == MSS[g,4],
      MSD[g] == MSS[g,5]
  },
  Af :> {
    Struct[{t,2,4}][
      Q,
      Af[{g1,3},{g2,3}]
    ]
  },
  Ae[ae] -> {
    Q -> Af**Q[2],
    Af[{g1,3},{g2,3}] -> Af**Af[g1,g2,2],
    Atau -> Af[3,3]
  },
  Au[au] -> {
    Q -> Af**Q[3],
    Af[{g1,3},{g2,3}] -> Af**Af[g1,g2,3],
    At -> Af[3,3]
  },
  Ad[ad] -> {
    Q -> Af**Q[4],
    Af[{g1,3},{g2,3}] -> Af**Af[g1,g2,4],
    Ab -> Af[3,3]
  },
  Yf :> {
    Struct[{t,2,4}][
      Q,
      Yf[{g1,3},{g2,3}]
    ]
  },
  Ye[ye] -> {
    Q -> Yf**Q[2],
    Yf[{g1,3},{g2,3}] -> Yf**Yf[g1,g2,2],
    Ytau -> Yf[3,3]
  },
  Yu[yu] -> {
    Q -> Yf**Q[3],
    Yf[{g1,3},{g2,3}] -> Yf**Yf[g1,g2,3],
    Yt -> Yf[3,3]
  },
  Yd[yd] -> {
    Q -> Yf**Q[4],
    Yf[{g1,3},{g2,3}] -> Yf**Yf[g1,g2,4],
    Yb -> Yf[3,3]
  },
  RVLamLLEIn[lllei] -> {
    lamLLE[{i,3},{j,3},{k,3}]
  },
  RVLamLQDIn[llqdi] -> {
    lamLQD[{i,3},{j,3},{k,3}]
  },
  RVLamUDDIn[luddi] -> {
    lamUDD[{i,3},{j,3},{k,3}]
  },
  RVLamLLE[llle] -> {
    Q,
    lamLLE[{i,3},{j,3},{k,3}]
  },
  RVLamLQD[llqd] -> {
    Q,
    lamLQD[{i,3},{j,3},{k,3}]
  },
  RVLamUDD[ludd] -> {
    Q,
    lamUDD[{i,3},{j,3},{k,3}]
  },
  RVTLLEIn[tllei] -> {
    TLLE[{i,3},{j,3},{k,3}]
  },
  RVTLQDIn[tlqdi] -> {
    TLQD[{i,3},{j,3},{k,3}]
  },
  RVTUDDIn[tuddi] -> {
    TUDD[{i,3},{j,3},{k,3}]
  },
  RVTLLE[tlle] -> {
    Q,
    TLLE[{i,3},{j,3},{k,3}]
  },
  RVTLQD[tlqd] -> {
    Q,
    TLQD[{i,3},{j,3},{k,3}]
  },
  RVTUDD[tudd] -> {
    Q,
    TUDD[{i,3},{j,3},{k,3}]
  },
  RVKappaIn[rki] -> {
    kappa[{i,3}]
  },
  RVKappa[rk] -> {
    Q,
    kappa[{i,3}]
  },
  RVDIn[rdi] -> {
    D[{i,3}]
  },
  RVD[rd] -> {
    Q,
    D[{i,3}]
  },
  RVSnVEVIn[rvevi] -> {
    VEV[{i,3}]
  },
  RVSnVEV[rvev] -> {
    Q,
    VEV[{i,3}]
  },
  RVM2LH1In[rmlhi] -> {
    M2LH1[{i,3}]
  },
  RVM2LH1[rmlh] -> {
    Q,
    M2LH1[{i,3}]
  },
  RVNMix[rnmix] -> {
    ZNeu[{n1,7},{n2,7}]
  },
  RVUMix[rumix] -> {
    UCha[{c1,5},{c2,5}]
  },
  RVVMix[rvmix] -> {
    VCha[{c1,5},{c2,5}]
  },
  RVHMix[rhmix] -> {
    UH[{h1,5},{h2,5}]
  },
  RVAMix[ramix] -> {
    UA[{h1,5},{h2,5}]
  },
  RVLMix[rlmix] -> {
    CLep[{l1,8},{l2,8}]
  },
  VCKMIn[vckmi] -> {
    lambda,
    A,
    rhobar,
    etabar
  },
  VCKM[vckm] -> {
    Q,
    VCKM[{g1,3},{g2,3}]
  },
  UPMNSIn[umnsi] -> {
    theta12,
    theta23,
    theta13,
    delta13,
    alpha1,
    alpha2
  },
  UPMNS[umns] -> {
    Q,
    UPMNS[{g1,3},{g2,3}]
  },
  MSS2In :> {
    Struct[{q,5}][
      MSS2[{g1,3},{g2,3}]
    ]
  },
  MSL2In[msl2i] -> {
    MSL2[{g1,3},{g2,3}] -> MSS2In**MSS2[g1,g2,1]
  },
  MSE2In[mse2i] -> {
    MSE2[{g1,3},{g2,3}] -> MSS2In**MSS2[g1,g2,2]
  },
  MSQ2In[msq2i] -> {
    MSQ2[{g1,3},{g2,3}] -> MSS2In**MSS2[g1,g2,3]
  },
  MSU2In[msu2i] -> {
    MSU2[{g1,3},{g2,3}] -> MSS2In**MSS2[g1,g2,4]
  },
  MSD2In[msd2i] -> {
    MSD2[{g1,3},{g2,3}] -> MSS2In**MSS2[g1,g2,5]
  },
  MSS2 :> {
    Struct[{q,5}][
      Q,
      MSS2[{g1,3},{g2,3}]
    ]
  },
  MSL2[msl2] -> {
    Q -> MSS2**Q[1],
    MSL2[{g1,3},{g2,3}] -> MSS2**MSS2[g1,g2,1]
  },
  MSE2[mse2] -> {
    Q -> MSS2**Q[2],
    MSE2[{g1,3},{g2,3}] -> MSS2**MSS2[g1,g2,2]
  },
  MSQ2[msq2] -> {
    Q -> MSS2**Q[3],
    MSQ2[{g1,3},{g2,3}] -> MSS2**MSS2[g1,g2,3]
  },
  MSU2[msu2] -> {
    Q -> MSS2**Q[4],
    MSU2[{g1,3},{g2,3}] -> MSS2**MSS2[g1,g2,4]
  },
  MSD2[msd2] -> {
    Q -> MSS2**Q[5],
    MSD2[{g1,3},{g2,3}] -> MSS2**MSS2[g1,g2,5]
  },
  TfIn :> {
    Struct[{t,2,4}][
      Tf[{g1,3},{g2,3}]
    ]
  },
  TeIn[tei] -> {
    Tf[{g1,3},{g2,3}] -> TfIn**Tf[g1,g2,2]
  },
  TuIn[tui] -> {
    Tf[{g1,3},{g2,3}] -> TfIn**Tf[g1,g2,3]
  },
  TdIn[tdi] -> {
    Tf[{g1,3},{g2,3}] -> TfIn**Tf[g1,g2,4]
  },
  Tf :> {
    Struct[{t,2,4}][
      Q,
      Tf[{g1,3},{g2,3}]
    ]
  },
  Te[te] -> {
    Q -> Tf**Q[2],
    Tf[{g1,3},{g2,3}] -> Tf**Tf[g1,g2,2]
  },
  Tu[tu] -> {
    Q -> Tf**Q[3],
    Tf[{g1,3},{g2,3}] -> Tf**Tf[g1,g2,3]
  },
  Td[td] -> {
    Q -> Tf**Q[4],
    Tf[{g1,3},{g2,3}] -> Tf**Tf[g1,g2,4]
  },
  ASfMix :> {
    Struct[{t,4}][
      UASf[{s1,6},{s2,6}]
    ]
  },
  SnuMix[snmix] -> {
    UASf[{s1,6},{s2,6}] -> ASfMix**UASf[s1,s2,1]
  },
  SelMix[slmix] -> {
    UASf[{s1,6},{s2,6}] -> ASfMix**UASf[s1,s2,2]
  },
  USqMix[usmix] -> {
    UASf[{s1,6},{s2,6}] -> ASfMix**UASf[s1,s2,3]
  },
  DSqMix[tsmix] -> {
    UASf[{s1,6},{s2,6}] -> ASfMix**UASf[s1,s2,4]
  },
  SnSMix[ssmix] -> {
    US[{g1,3},{g2,3}]
  },
  SnAMix[samix] -> {
    UA[{g1,3},{g2,3}]
  },
  CVHMix[hcmix] -> {
    UH[{h1,4},{h2,4}]
  },
  NMNMix[nnmix] -> {
    ZNeu[{n1,5},{n2,5}]
  },
  NMHMix[nhmix] -> {
    UH[{h1,3},{h2,3}]
  },
  NMAMix[namix] -> {
    UA[{h1,3},{h2,3}]
  },
  PrecObs[prcob] -> {
    DeltaRho,
    MWMSSM,
    MWSM,
    SW2effMSSM,
    SW2effSM,
    gminus2mu,
    EDMeTh,
    EDMn,
    EDMHg,
    bsgammaMSSM,
    bsgammaSM,
    DeltaMsMSSM,
    DeltaMsSM,
    BsmumuMSSM,
    BsmumuSM
  },
  SPInfo[spinf] -> {
    NLines,
    Severity,
    Code[{n,MAXMSG}],
    Text[{i,MSGLEN},{n,MAXMSG}],
    Len == MSGLEN CPLXSIZE
  },
  DCInfo[tcinf] -> {
    NLines,
    Severity,
    Code[{n,MAXMSG}],
    Text[{i,MSGLEN},{n,MAXMSG}],
    Len == MSGLEN CPLXSIZE
  },
  Decays[tecys] -> {
    Data[{n,4096}]
  }
}


WriteString["SLHADefs.h", SLHADefs[desc]]

WriteString["SLHAReadBlocks.h",
  SLHANames[desc, "blockname", ""]]

WriteString["SLHAWriteBlocks.h",
  SLHANames[desc, "blockname", "IM"]]


