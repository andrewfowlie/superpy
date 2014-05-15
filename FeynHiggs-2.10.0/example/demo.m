(*
   demo.m
   demonstration program for calling FeynHiggs from Mathematica
   this file is part of FeynHiggs
   last modified 5 May 11 th
*)


Install["build/MFeynHiggs"]


mhiggs[i_, TB_] :=
Block[ {mssmpart, fieldren, tanbren, higgsmix, p2approx,
looplevel, runningMT, botResum, tlCplxApprox,
invAlfa, AlfasMZ, GF, MU, MD, MC, MS, MB, MW, MZ,
CKMlambda, CKMA, CKMrho, CKMeta,
scalefactor, MT, MA0, MHp,
MSL, MSE, MSQ, MSU, MSD, MUE, Af, M1, M2, M3,
Qtau, Qt, Qb},
  mssmpart = 4;
  fieldren = 0;
  tanbren = 0;
  higgsmix = 2;
  p2approx = 0;
  looplevel = 2;
  runningMT = 1;
  botResum = 1;
  tlCplxApprox = 3;
  FHSetFlags[mssmpart, fieldren, tanbren, higgsmix, p2approx,
    looplevel, runningMT, botResum, tlCplxApprox];

  invAlfa = -1;
  AlfasMZ = -1;
  GF = -1;
  ME = -1;
  MU = -1;
  MD = -1;
  MM = -1;
  MC = -1;
  MS = -1;
  ML = -1;
  MB = -1;
  MW = -1;
  MZ = -1;

  CKMlambda = -1;
  CKMA = -1;
  CKMrhobar = -1;
  CKMetabar = -1;

  FHSetSMPara[invAlfa, AlfasMZ, GF,
    ME, MU, MD, MM, MC, MS, ML, MB,
    MW, MZ,
    CKMlambda, CKMA, CKMrho, CKMeta];

  scalefactor = 1;
  MT = 170.9;

  MA0 = 250;
  MHp = -1;
  _MSL = _MSE = _MSQ = _MSU = _MSD = 1000;

  _Af = 2000;
  MUE = 200;

  M1 = 0;
  M2 = 500;
  M3 = 800;

  Qtau = Qt = Qb = 0;

  FHSetPara[scalefactor,
    MT, TB, MA0, MHp,
    MSL[3], MSE[3], MSQ[3], MSU[3], MSD[3],
    MSL[2], MSE[2], MSQ[2], MSU[2], MSD[2],
    MSL[1], MSE[1], MSQ[1], MSU[1], MSD[1],
    MUE,
    Af[2,3], Af[3,3], Af[4,3],
    Af[2,2], Af[3,2], Af[4,2],
    Af[2,1], Af[3,1], Af[4,1],
    M1, M2, M3,
    Qtau, Qt, Qb];

  mass[i, FHHiggsCorr[]]
]


mass[i_, rul_List] := (MHiggs /. rul)[[i]]

_mass = $Failed



Plot[ mhiggs[1, TB], {TB, 1.5, 12} ]

