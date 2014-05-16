#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=51;
static ModelPrtclsStr ModelPrtcls_[51]=
{
  {"A","A", 22, "0","0",2,1,0}
, {"Z","Z", 23, "MZ","wZ",2,1,0}
, {"W+","W-", 24, "MW","wW",2,1,3}
, {"G","G", 21, "0","0",2,8,0}
, {"ne","Ne", 12, "0","0",1,1,0}
, {"e","E", 11, "0","0",1,1,-3}
, {"nm","Nm", 14, "0","0",1,1,0}
, {"m","M", 13, "0","0",1,1,-3}
, {"nl","Nl", 16, "0","0",1,1,0}
, {"l","L", 15, "Ml","0",1,1,-3}
, {"u","U", 2, "Mq","0",1,3,2}
, {"d","D", 1, "Mq","0",1,3,-1}
, {"s","S", 3, "Mq","0",1,3,-1}
, {"c","C", 4, "Mc","0",1,3,2}
, {"t","T", 6, "Mt","wt",1,3,2}
, {"b","B", 5, "Mb","0",1,3,-1}
, {"h1","h1", 25, "Mh1","wh1",0,1,0}
, {"h2","h2", 35, "Mh2","wh2",0,1,0}
, {"h3","h3", 45, "Mh3","wh3",0,1,0}
, {"ha","ha", 36, "Mha","wha",0,1,0}
, {"hb","hb", 46, "Mhb","whb",0,1,0}
, {"H+","H-", 37, "MHc","wHc",0,1,3}
, {"~1+","~1-", 1000024, "MC1","wC1",1,1,3}
, {"~2+","~2-", 1000037, "MC2","wC2",1,1,3}
, {"~o1","~o1", 1000022, "MNE1","0",1,1,0}
, {"~o2","~o2", 1000023, "MNE2","wNE2",1,1,0}
, {"~o3","~o3", 1000025, "MNE3","wNE3",1,1,0}
, {"~o4","~o4", 1000035, "MNE4","wNE4",1,1,0}
, {"~o5","~o5", 1000045, "MNE5","wNE5",1,1,0}
, {"~g","~g", 1000021, "MSG","wSG",1,8,0}
, {"~eL","~EL", 1000011, "MSeL","wSeL",0,1,-3}
, {"~eR","~ER", 2000011, "MSeR","wSeR",0,1,-3}
, {"~mL","~ML", 1000013, "MSmL","wSmL",0,1,-3}
, {"~mR","~MR", 2000013, "MSmR","wSmR",0,1,-3}
, {"~l1","~L1", 1000015, "MSl1","wSl1",0,1,-3}
, {"~l2","~L2", 2000015, "MSl2","wSl2",0,1,-3}
, {"~ne","~Ne", 1000012, "MSne","wSne",0,1,0}
, {"~nm","~Nm", 1000014, "MSnm","wSnm",0,1,0}
, {"~nl","~Nl", 1000016, "MSnl","wSnl",0,1,0}
, {"~uL","~UL", 1000002, "MSuL","wSuL",0,3,2}
, {"~uR","~UR", 2000002, "MSuR","wSuR",0,3,2}
, {"~dL","~DL", 1000001, "MSdL","wSdL",0,3,-1}
, {"~dR","~DR", 2000001, "MSdR","wSdR",0,3,-1}
, {"~cL","~CL", 1000004, "MScL","wScL",0,3,2}
, {"~cR","~CR", 2000004, "MScR","wScR",0,3,2}
, {"~sL","~SL", 1000003, "MSsL","wSsL",0,3,-1}
, {"~sR","~SR", 2000003, "MSsR","wSsR",0,3,-1}
, {"~t1","~T1", 1000006, "MSt1","wSt1",0,3,2}
, {"~t2","~T2", 2000006, "MSt2","wSt2",0,3,2}
, {"~b1","~B1", 1000005, "MSb1","wSb1",0,3,-1}
, {"~b2","~B2", 2000005, "MSb2","wSb2",0,3,-1}
};
ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=160;
int nModelFunc=341;
static char*varNames_[501]={
 "EE","SW","MZ","alfSMZ","MbMb","McMc","Mtp","Ml","wt","Mq"
,"MG1","MG2","MG3","tb","mu","Lambda","Kappa","aLambda","aKappa","At"
,"Ab","Al","Ml2","Ml3","Mr2","Mr3","Mq2","Mq3","Mu2","Mu3"
,"Md2","Md3","MA","MP","Mh1","Mh2","Mh3","Mha","Mhb","MHc"
,"MNE1","MNE2","MNE3","MNE4","MNE5","MC1","MC2","MSuL","MSuR","MSdL"
,"MSdR","MScL","MScR","MSsL","MSsR","MSt1","MSt2","MSb1","MSb2","MSeL"
,"MSeR","MSmL","MSmR","MSl1","MSl2","MSne","MSnm","MSnl","Zn11","Zn12"
,"Zn13","Zn14","Zn15","Zn21","Zn22","Zn23","Zn24","Zn25","Zn31","Zn32"
,"Zn33","Zn34","Zn35","Zn41","Zn42","Zn43","Zn44","Zn45","Zn51","Zn52"
,"Zn53","Zn54","Zn55","Zu11","Zu12","Zu21","Zu22","Zv11","Zv12","Zv21"
,"Zv22","Zt11","Zt12","Zt21","Zt22","Zb11","Zb12","Zb21","Zb22","Zl11"
,"Zl12","Zl21","Zl22","Zh11","Zh12","Zh13","Zh21","Zh22","Zh23","Zh31"
,"Zh32","Zh33","Za11","Za12","Za13","Za21","Za22","Za23","la1","la2"
,"la3","la4","la5","la6","la7","la1s","la2s","la3s","la4s","la5s"
,"la6s","la7s","la8s","aa1","aa2","aa3","aa4","aa5","aa6","B1"
,"B2","X","dMb","Q","Au","Ad","Am","Maux","tB","MSG"
,"CW","C2W","MW","LamQCD","sb","cb","t2b","xvev","Pa12","Pa22"
,"Pa11","Pa21","Td3","hMM11","hMM12","hMM13","hMM22","hMM23","hMM33","idS"
,"Zh11_","Zh12_","Zh13_","Zh21_","Zh22_","Zh23_","Zh31_","Zh32_","Zh33_","MA11"
,"MA12","MA22","idA","Pa11_","Pa12_","Pa21_","Pa22_","Mt","Mb","Mc"
,"xH2","yH2","dMd","Td2","fiuu","fidd","ficc","fiss","Zuu11","Zuu12"
,"Zuu21","Zuu22","Zdd11","Zdd12","Zdd21","Zdd22","Zcc11","Zcc12","Zcc21","Zcc22"
,"Zss11","Zss12","Zss21","Zss22","StMM11","StMM12","StMM22","MtMM","AT","SbMM11"
,"SbMM12","SbMM22","MbMM","AB","dX3","dX2","dX1","NMM11","NMM12","NMM13"
,"NMM14","NMM15","NMM22","NMM23","NMM24","NMM25","NMM33","NMM34","NMM35","NMM44"
,"NMM45","NMM55","MG1I","MG2I","PI","aQCD","Mbp","Mcp","rhF_c","ihF_c"
,"rAF_c","iAF_c","rhF_b","ihF_b","rAF_b","iAF_b","rhF_t","ihF_t","rAF_t","iAF_t"
,"rhF_l","ihF_l","rAF_l","iAF_l","rhS_eL","ihS_eL","rhS_eR","ihS_eR","rhS_mL","ihS_mL"
,"rhS_mR","ihS_mR","rhS_l1","ihS_l1","rhS_l2","ihS_l2","rhS_uL","ihS_uL","rhS_uR","ihS_uR"
,"rhS_cL","ihS_cL","rhS_cR","ihS_cR","rhS_t1","ihS_t1","rhS_t2","ihS_t2","rhS_dL","ihS_dL"
,"rhS_dR","ihS_dR","rhS_sL","ihS_sL","rhS_sR","ihS_sR","rhS_b1","ihS_b1","rhS_b2","ihS_b2"
,"rhS_Hc","ihS_Hc","rhV_W","ihV_W","rhF_c1","ihF_c1","rAF_c1","iAF_c1","rhF_c2","ihF_c2"
,"rAF_c2","iAF_c2","McR","MbR","MtR","rhF1_c","ihF1_c","rAF1_c","iAF1_c","rhF1_b"
,"ihF1_b","rAF1_b","iAF1_b","rhF1_t","ihF1_t","rAF1_t","iAF1_t","tau2_uL","rhS1_uL","ihS1_uL"
,"tau2_uR","rhS1_uR","ihS1_uR","tau2_cL","rhS1_cL","ihS1_cL","tau2_cR","rhS1_cR","ihS1_cR","tau2_t1"
,"rhS1_t1","ihS1_t1","tau2_t2","rhS1_t2","ihS1_t2","tau2_dL","rhS1_dL","ihS1_dL","tau2_dR","rhS1_dR"
,"ihS1_dR","tau2_sL","rhS1_sL","ihS1_sL","tau2_sR","rhS1_sR","ihS1_sR","tau2_b1","rhS1_b1","ihS1_b1"
,"tau2_b2","rhS1_b2","ihS1_b2","ah1F_c","ah1F_b","ah1F_t","ah1F_l","ah2F_c","ah2F_b","ah2F_t"
,"ah2F_l","ah3F_c","ah3F_b","ah3F_t","ah3F_l","ah1S_eL","ah1S_eR","ah1S_mL","ah1S_mR","ah1S_uL"
,"ah1S_uR","ah1S_cL","ah1S_cR","ah1S_dL","ah1S_dR","ah1S_sL","ah1S_sR","ah1S_l1","ah1S_l2","ah1S_t1"
,"ah1S_t2","ah1S_b1","ah1S_b2","ah2S_eL","ah2S_eR","ah2S_mL","ah2S_mR","ah2S_uL","ah2S_uR","ah2S_cL"
,"ah2S_cR","ah2S_dL","ah2S_dR","ah2S_sL","ah2S_sR","ah2S_l1","ah2S_l2","ah2S_t1","ah2S_t2","ah2S_b1"
,"ah2S_b2","ah3S_eL","ah3S_eR","ah3S_mL","ah3S_mR","ah3S_uL","ah3S_uR","ah3S_cL","ah3S_cR","ah3S_dL"
,"ah3S_dR","ah3S_sL","ah3S_sR","ah3S_l1","ah3S_l2","ah3S_t1","ah3S_t2","ah3S_b1","ah3S_b2","ah1V_W"
,"B00000","B00001","B00002","B00003","B00004","ah1S_Hc","ah2V_W","B00005","B00006","B00007"
,"B00008","B00009","ah2S_Hc","ah3V_W","B00010","B00011","B00012","B00013","B00014","ah3S_Hc"
,"ah1F_c1","ah1F_c2","ah2F_c1","ah2F_c2","ah3F_c1","ah3F_c2","ahaF_c","ahaF_b","ahaF_t","ahaF_l"
,"ahbF_c","ahbF_b","ahbF_t","ahbF_l","ahaF_c1","ahaF_c2","ahbF_c1","ahbF_c2","Rqcd","lnTop"
,"Ctop","Cq","Csq","alphaE0","Qu","Qd","Fodd","LGGh1","LAAh1","LGGh2"
,"LAAh2","LGGh3","LAAh3","LGGha","LAAha","LGGhb","LAAhb","aSMhF_f","aSMhV_W","LGGSM"
,"LAASM"};
char**varNames=varNames_;
static REAL varValues_[501]={
   3.128530E-01,  4.820000E-01,  9.120000E+01,  1.184000E-01,  4.500000E+00,  1.300000E+00,  1.730700E+02,  1.777000E+00,  1.442000E+00,  5.000000E-02
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  2.000000E+02,  2.000000E+02,  2.020000E+02,  2.000000E+02,  1.000000E+03,  1.000000E+03,  1.000000E+03,  1.000000E+03
,  1.000000E+03,  1.000000E+03,  2.000000E+02,  1.000000E+03,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  1.000000E+02,  0.000000E+00,  0.000000E+00,  0.000000E+00,  1.000000E+00,  0.000000E+00,  0.000000E+00
};
REAL*varValues=varValues_;
int calcMainFunc(void)
{
   int i;
   static REAL * VV=NULL;
   static int iQ=-1;
   static int cErr=1;
   REAL *V=varValues;
   FError=0;
   if(VV && cErr==0)
   { for(i=0;i<nModelVars;i++) if(i!=iQ && VV[i]!=V[i]) break;
     if(i==nModelVars)      {if(iQ>=0 && VV[iQ]!=V[iQ]) goto FirstQ; else return 0;} 
   }
  cErr=1;
   V[160]=sqrt(1-pow(V[1],2));
   if(!isfinite(V[160]) || FError) return 160;
   V[161]=pow(V[160],2)-pow(V[1],2);
   if(!isfinite(V[161]) || FError) return 161;
   V[162]=V[2]*V[160];
   if(!isfinite(V[162]) || FError) return 162;
   V[163]=initQCD5(V[3],V[5],V[4],V[6]);
   if(!isfinite(V[163]) || FError) return 163;
   V[164]=V[158]/(sqrt(1+pow(V[158],2)));
   if(!isfinite(V[164]) || FError) return 164;
   V[165]=sqrt(1-pow(V[164],2));
   if(!isfinite(V[165]) || FError) return 165;
   V[166]=2*V[158]/(1-pow(V[158],2));
   if(!isfinite(V[166]) || FError) return 166;
   V[167]=V[14]/(V[15]);
   if(!isfinite(V[167]) || FError) return 167;
   V[168]=V[124];

   V[169]=V[127];

   V[170]=(V[122]*V[169]>0 ? V[169] : -V[169]);
   if(!isfinite(V[170]) || FError) return 170;
   V[171]=(V[125]*V[168]>0 ? V[168] : -V[168]);
   if(!isfinite(V[171]) || FError) return 171;
   V[172]=V[152]/(1+V[152]);
   if(!isfinite(V[172]) || FError) return 172;
   V[173]=4/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[165]*V[165]*V[128]-6/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[165]*V[133]*V[164]+1*V[149]/(V[165])*V[164]+1*V[147]/(V[165])*V[164]*V[167]+1*V[148]/(V[165])*V[164]*V[167]+1/(V[165])*V[138]*V[164]*pow(V[167],2)+1/(V[165])*V[139]*V[164]*pow(V[167],2)+1/(V[165])*V[140]*V[164]*pow(V[167],2)+2/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)/(V[165])*V[134]*pow(V[164],3);
   if(!isfinite(V[173]) || FError) return 173;
   V[174]=4/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[165]*V[130]*V[164]+4/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[165]*V[131]*V[164]+4/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[165]*V[132]*V[164]-6/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[165]*V[165]*V[133]-6/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[134]*pow(V[164],2)-1*V[138]*pow(V[167],2)-1*V[139]*pow(V[167],2)-1*V[140]*pow(V[167],2)-1*V[147]*V[167]-1*V[148]*V[167]-1*V[149];
   if(!isfinite(V[174]) || FError) return 174;
   V[175]=2/(V[0])*V[162]*V[1]*M_SQRT2*V[165]*V[135]*V[167]-2/(V[0])*V[162]*V[1]*M_SQRT2*V[138]*V[164]*V[167]-2/(V[0])*V[162]*V[1]*M_SQRT2*V[139]*V[164]*V[167]-2/(V[0])*V[162]*V[1]*M_SQRT2*V[140]*V[164]*V[167]+2/(V[0])*V[162]*V[1]*M_SQRT2*V[143]*V[165]-1/(V[0])*V[162]*V[1]*M_SQRT2*V[147]*V[164]-1/(V[0])*V[162]*V[1]*M_SQRT2*V[148]*V[164];
   if(!isfinite(V[175]) || FError) return 175;
   V[176]=4/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[129]*pow(V[164],2)-6/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[165]*V[134]*V[164]+1*V[149]*V[165]/(V[164])+1*V[147]*V[165]/(V[164])*V[167]+1*V[148]*V[165]/(V[164])*V[167]+2/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*pow(V[165],3)*V[133]/(V[164])+1*V[165]*V[138]/(V[164])*pow(V[167],2)+1*V[165]*V[139]/(V[164])*pow(V[167],2)+1*V[165]*V[140]/(V[164])*pow(V[167],2);
   if(!isfinite(V[176]) || FError) return 176;
   V[177]=2/(V[0])*V[162]*V[1]*M_SQRT2*V[136]*V[164]*V[167]-2/(V[0])*V[162]*V[1]*M_SQRT2*V[165]*V[138]*V[167]-2/(V[0])*V[162]*V[1]*M_SQRT2*V[165]*V[139]*V[167]-2/(V[0])*V[162]*V[1]*M_SQRT2*V[165]*V[140]*V[167]+2/(V[0])*V[162]*V[1]*M_SQRT2*V[144]*V[164]-1/(V[0])*V[162]*V[1]*M_SQRT2*V[147]*V[165]-1/(V[0])*V[162]*V[1]*M_SQRT2*V[148]*V[165];
   if(!isfinite(V[177]) || FError) return 177;
   V[178]=4*V[137]*pow(V[167],2)+8*V[141]*pow(V[167],2)+8*V[142]*pow(V[167],2)+3*V[145]*V[167]+3*V[146]*V[167]-1*V[151]/(V[167])-2/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[143]*V[165]*V[165]/(V[167])-2/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[144]*pow(V[164],2)/(V[167])+2/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[147]*V[165]*V[164]/(V[167])+2/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[148]*V[165]*V[164]/(V[167]);
   if(!isfinite(V[178]) || FError) return 178;
   V[179]=rDiagonal(3,V[173],V[174],V[175],V[176],V[177],V[178]);
   if(!isfinite(V[179]) || FError) return 179;
   V[180]=MixMatrix(V[179],1,1);
   if(!isfinite(V[180]) || FError) return 180;
   V[181]=MixMatrix(V[179],1,2);
   if(!isfinite(V[181]) || FError) return 181;
   V[182]=MixMatrix(V[179],1,3);
   if(!isfinite(V[182]) || FError) return 182;
   V[183]=MixMatrix(V[179],2,1);
   if(!isfinite(V[183]) || FError) return 183;
   V[184]=MixMatrix(V[179],2,2);
   if(!isfinite(V[184]) || FError) return 184;
   V[185]=MixMatrix(V[179],2,3);
   if(!isfinite(V[185]) || FError) return 185;
   V[186]=MixMatrix(V[179],3,1);
   if(!isfinite(V[186]) || FError) return 186;
   V[187]=MixMatrix(V[179],3,2);
   if(!isfinite(V[187]) || FError) return 187;
   V[188]=MixMatrix(V[179],3,3);
   if(!isfinite(V[188]) || FError) return 188;
   V[189]=2/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[165]*V[133]/(V[164])-4/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)*V[132]+2/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)/(V[165])*V[134]*V[164]+1*V[149]/(V[165])/(V[164])+1*V[147]/(V[165])/(V[164])*V[167]+1*V[148]/(V[165])/(V[164])*V[167]+1/(V[165])*V[138]/(V[164])*pow(V[167],2)+1/(V[165])*V[139]/(V[164])*pow(V[167],2)+1/(V[165])*V[140]/(V[164])*pow(V[167],2);
   if(!isfinite(V[189]) || FError) return 189;
   V[190]=(2*V[139]*V[167]-2*V[140]*V[167]+V[147]-V[148])*V[162]*V[1]*M_SQRT2/(V[0]);
   if(!isfinite(V[190]) || FError) return 190;
   V[191]=(4*V[165]*V[139]*V[164]+4*V[165]*V[140]*V[164]+(V[147]*V[165]*V[164]-V[144]*pow(V[164],2)-V[143]*V[165]*V[165]+V[148]*V[165]*V[164])/(V[167]))*2/(pow(V[0],2))*pow(V[162],2)*pow(V[1],2)-16*V[141]*pow(V[167],2)-1*V[145]*V[167]-9*V[146]*V[167]-4*V[150]-1*V[151]/(V[167])-4*V[142]*pow(V[167],2);
   if(!isfinite(V[191]) || FError) return 191;
   V[192]=rDiagonal(2,V[189],V[190],V[191]);
   if(!isfinite(V[192]) || FError) return 192;
   V[193]=MixMatrix(V[192],1,1);
   if(!isfinite(V[193]) || FError) return 193;
   V[194]=MixMatrix(V[192],1,2);
   if(!isfinite(V[194]) || FError) return 194;
   V[195]=MixMatrix(V[192],2,1);
   if(!isfinite(V[195]) || FError) return 195;
   V[196]=MixMatrix(V[192],2,2);
   if(!isfinite(V[196]) || FError) return 196;
 FirstQ:
 cErr=1;
   V[197]=MtEff(V[153]);
   if(!isfinite(V[197]) || FError) return 197;
   V[198]=MbEff(V[153]);
   if(!isfinite(V[198]) || FError) return 198;
   V[199]=McEff(V[153]);
   if(!isfinite(V[199]) || FError) return 199;
   V[200]=pow(V[53],2)/(pow(V[159],2));
   if(!isfinite(V[200]) || FError) return 200;
   V[201]=pow(V[54],2)/(pow(V[159],2));
   if(!isfinite(V[201]) || FError) return 201;
   V[202]=-3/(double)((2))*alphaQCD(V[159])/(3.1415)*V[158]*V[14]/(V[159])*(V[200]*log(V[200])/(1-V[200])-V[201]*log(V[201])/(1-V[201]))/(V[200]-V[201]);
   if(!isfinite(V[202]) || FError) return 202;
   V[203]=V[202]/(1+V[202]);
   if(!isfinite(V[203]) || FError) return 203;
   V[204]=atan(-2*V[9]*(V[154]-V[14]/(V[158]))/(pow(V[48],2)-pow(V[47],2)))/(2);
   if(!isfinite(V[204]) || FError) return 204;
   V[205]=atan(-2*V[9]*(1-V[203])*(V[155]-V[14]*V[158])/(pow(V[50],2)-pow(V[49],2)))/(2);
   if(!isfinite(V[205]) || FError) return 205;
   V[206]=atan(-2*V[199]*(V[154]-V[14]/(V[158]))/(pow(V[52],2)-pow(V[51],2)))/(2);
   if(!isfinite(V[206]) || FError) return 206;
   V[207]=atan(-2*V[9]*(1-V[203])*(V[155]-V[14]*V[158])/(pow(V[54],2)-pow(V[53],2)))/(2);
   if(!isfinite(V[207]) || FError) return 207;
   V[208]=cos(V[204]);
   if(!isfinite(V[208]) || FError) return 208;
   V[209]=sin(V[204]);
   if(!isfinite(V[209]) || FError) return 209;
   V[210]=-sin(V[204]);
   if(!isfinite(V[210]) || FError) return 210;
   V[211]=cos(V[204]);
   if(!isfinite(V[211]) || FError) return 211;
   V[212]=cos(V[205]);
   if(!isfinite(V[212]) || FError) return 212;
   V[213]=sin(V[205]);
   if(!isfinite(V[213]) || FError) return 213;
   V[214]=-sin(V[205]);
   if(!isfinite(V[214]) || FError) return 214;
   V[215]=cos(V[205]);
   if(!isfinite(V[215]) || FError) return 215;
   V[216]=cos(V[206]);
   if(!isfinite(V[216]) || FError) return 216;
   V[217]=sin(V[206]);
   if(!isfinite(V[217]) || FError) return 217;
   V[218]=-sin(V[206]);
   if(!isfinite(V[218]) || FError) return 218;
   V[219]=cos(V[206]);
   if(!isfinite(V[219]) || FError) return 219;
   V[220]=cos(V[207]);
   if(!isfinite(V[220]) || FError) return 220;
   V[221]=sin(V[207]);
   if(!isfinite(V[221]) || FError) return 221;
   V[222]=-sin(V[207]);
   if(!isfinite(V[222]) || FError) return 222;
   V[223]=cos(V[207]);
   if(!isfinite(V[223]) || FError) return 223;
   V[224]=V[101]*pow(V[55],2)*V[101]+V[103]*pow(V[56],2)*V[103];
   if(!isfinite(V[224]) || FError) return 224;
   V[225]=V[101]*pow(V[55],2)*V[102]+V[103]*pow(V[56],2)*V[104];
   if(!isfinite(V[225]) || FError) return 225;
   V[226]=V[102]*pow(V[55],2)*V[102]+V[104]*pow(V[56],2)*V[104];
   if(!isfinite(V[226]) || FError) return 226;
   V[227]=MtRun(sqrt(V[55]*V[56]));
   if(!isfinite(V[227]) || FError) return 227;
   V[228]=V[225]/(V[227])+V[15]*V[165]/(V[164])*V[167];
   if(!isfinite(V[228]) || FError) return 228;
   V[229]=V[105]*pow(V[57],2)*V[105]+V[107]*pow(V[58],2)*V[107];
   if(!isfinite(V[229]) || FError) return 229;
   V[230]=V[105]*pow(V[57],2)*V[106]+V[107]*pow(V[58],2)*V[108];
   if(!isfinite(V[230]) || FError) return 230;
   V[231]=V[106]*pow(V[57],2)*V[106]+V[108]*pow(V[58],2)*V[108];
   if(!isfinite(V[231]) || FError) return 231;
   V[232]=MbRun(sqrt(V[57]*V[58]));
   if(!isfinite(V[232]) || FError) return 232;
   V[233]=V[230]/(V[232])+V[15]*V[164]/(V[165])*V[167];
   if(!isfinite(V[233]) || FError) return 233;
   V[234]=(V[224]-V[229]+pow(V[232],2)-pow(V[227],2)+2*pow(V[162],2)*pow(V[164],2)-pow(V[162],2))/(pow(V[162],2)*pow(V[164],2));
   if(!isfinite(V[234]) || FError) return 234;
   V[235]=(pow(V[51],2)-pow(V[53],2)+pow(V[162],2)*(2*pow(V[164],2)*V[1]-2*pow(V[164],2)*pow(V[1],3)+pow(V[1],3)-V[1])/(pow(V[160],2)*V[1]))/(pow(V[162],2)*pow(V[164],2));
   if(!isfinite(V[235]) || FError) return 235;
   V[236]=(pow(V[47],2)-pow(V[49],2)+pow(V[162],2)*(2*pow(V[164],2)*V[1]-2*pow(V[164],2)*pow(V[1],3)+pow(V[1],3)-V[1])/(pow(V[160],2)*V[1]))/(pow(V[162],2)*pow(V[164],2));
   if(!isfinite(V[236]) || FError) return 236;
   V[237]=V[68]*V[40]*V[68]+V[73]*V[41]*V[73]+V[78]*V[42]*V[78]+V[83]*V[43]*V[83]+V[88]*V[44]*V[88];
   if(!isfinite(V[237]) || FError) return 237;
   V[238]=V[68]*V[40]*V[69]+V[73]*V[41]*V[74]+V[78]*V[42]*V[79]+V[83]*V[43]*V[84]+V[88]*V[44]*V[89];
   if(!isfinite(V[238]) || FError) return 238;
   V[239]=V[68]*V[40]*V[70]+V[73]*V[41]*V[75]+V[78]*V[42]*V[80]+V[83]*V[43]*V[85]+V[88]*V[44]*V[90];
   if(!isfinite(V[239]) || FError) return 239;
   V[240]=V[68]*V[40]*V[71]+V[73]*V[41]*V[76]+V[78]*V[42]*V[81]+V[83]*V[43]*V[86]+V[88]*V[44]*V[91];
   if(!isfinite(V[240]) || FError) return 240;
   V[241]=V[68]*V[40]*V[72]+V[73]*V[41]*V[77]+V[78]*V[42]*V[82]+V[83]*V[43]*V[87]+V[88]*V[44]*V[92];
   if(!isfinite(V[241]) || FError) return 241;
   V[242]=V[69]*V[40]*V[69]+V[74]*V[41]*V[74]+V[79]*V[42]*V[79]+V[84]*V[43]*V[84]+V[89]*V[44]*V[89];
   if(!isfinite(V[242]) || FError) return 242;
   V[243]=V[69]*V[40]*V[70]+V[74]*V[41]*V[75]+V[79]*V[42]*V[80]+V[84]*V[43]*V[85]+V[89]*V[44]*V[90];
   if(!isfinite(V[243]) || FError) return 243;
   V[244]=V[69]*V[40]*V[71]+V[74]*V[41]*V[76]+V[79]*V[42]*V[81]+V[84]*V[43]*V[86]+V[89]*V[44]*V[91];
   if(!isfinite(V[244]) || FError) return 244;
   V[245]=V[69]*V[40]*V[72]+V[74]*V[41]*V[77]+V[79]*V[42]*V[82]+V[84]*V[43]*V[87]+V[89]*V[44]*V[92];
   if(!isfinite(V[245]) || FError) return 245;
   V[246]=V[70]*V[40]*V[70]+V[75]*V[41]*V[75]+V[80]*V[42]*V[80]+V[85]*V[43]*V[85]+V[90]*V[44]*V[90];
   if(!isfinite(V[246]) || FError) return 246;
   V[247]=V[70]*V[40]*V[71]+V[75]*V[41]*V[76]+V[80]*V[42]*V[81]+V[85]*V[43]*V[86]+V[90]*V[44]*V[91];
   if(!isfinite(V[247]) || FError) return 247;
   V[248]=V[70]*V[40]*V[72]+V[75]*V[41]*V[77]+V[80]*V[42]*V[82]+V[85]*V[43]*V[87]+V[90]*V[44]*V[92];
   if(!isfinite(V[248]) || FError) return 248;
   V[249]=V[71]*V[40]*V[71]+V[76]*V[41]*V[76]+V[81]*V[42]*V[81]+V[86]*V[43]*V[86]+V[91]*V[44]*V[91];
   if(!isfinite(V[249]) || FError) return 249;
   V[250]=V[71]*V[40]*V[72]+V[76]*V[41]*V[77]+V[81]*V[42]*V[82]+V[86]*V[43]*V[87]+V[91]*V[44]*V[92];
   if(!isfinite(V[250]) || FError) return 250;
   V[251]=V[72]*V[40]*V[72]+V[77]*V[41]*V[77]+V[82]*V[42]*V[82]+V[87]*V[43]*V[87]+V[92]*V[44]*V[92];
   if(!isfinite(V[251]) || FError) return 251;
   V[252]=V[237];

   V[253]=V[242];

   V[254]=4*atan(1);
   if(!isfinite(V[254]) || FError) return 254;
   V[255]=alphaQCD(V[153])/(V[254]);
   if(!isfinite(V[255]) || FError) return 255;
   V[256]=V[4]*(1+4/(double)((3))*alphaQCD(V[4])/(V[254]));
   if(!isfinite(V[256]) || FError) return 256;
   V[257]=V[5]*(1+4/(double)((3))*alphaQCD(V[5])/(V[254]));
   if(!isfinite(V[257]) || FError) return 257;
   V[258]=creal(HggF(pow(V[153]/(2)/(V[257]),2)));
   if(!isfinite(V[258]) || FError) return 258;
   V[259]=cimag(HggF(pow(V[153]/(2)/(V[257]),2)));
   if(!isfinite(V[259]) || FError) return 259;
   V[260]=creal(HggA(pow(V[153]/(2)/(V[257]),2)));
   if(!isfinite(V[260]) || FError) return 260;
   V[261]=cimag(HggA(pow(V[153]/(2)/(V[257]),2)));
   if(!isfinite(V[261]) || FError) return 261;
   V[262]=creal(HggF(pow(V[153]/(2)/(V[256]),2)));
   if(!isfinite(V[262]) || FError) return 262;
   V[263]=cimag(HggF(pow(V[153]/(2)/(V[256]),2)));
   if(!isfinite(V[263]) || FError) return 263;
   V[264]=creal(HggA(pow(V[153]/(2)/(V[256]),2)));
   if(!isfinite(V[264]) || FError) return 264;
   V[265]=cimag(HggA(pow(V[153]/(2)/(V[256]),2)));
   if(!isfinite(V[265]) || FError) return 265;
   V[266]=creal(HggF(pow(V[153]/(2)/(V[6]),2)));
   if(!isfinite(V[266]) || FError) return 266;
   V[267]=cimag(HggF(pow(V[153]/(2)/(V[6]),2)));
   if(!isfinite(V[267]) || FError) return 267;
   V[268]=creal(HggA(pow(V[153]/(2)/(V[6]),2)));
   if(!isfinite(V[268]) || FError) return 268;
   V[269]=cimag(HggA(pow(V[153]/(2)/(V[6]),2)));
   if(!isfinite(V[269]) || FError) return 269;
   V[270]=creal(HggF(pow(V[153]/(2)/(V[7]),2)));
   if(!isfinite(V[270]) || FError) return 270;
   V[271]=cimag(HggF(pow(V[153]/(2)/(V[7]),2)));
   if(!isfinite(V[271]) || FError) return 271;
   V[272]=creal(HggA(pow(V[153]/(2)/(V[7]),2)));
   if(!isfinite(V[272]) || FError) return 272;
   V[273]=cimag(HggA(pow(V[153]/(2)/(V[7]),2)));
   if(!isfinite(V[273]) || FError) return 273;
   V[274]=creal(HggS(pow(V[153]/(2)/(V[59]),2)));
   if(!isfinite(V[274]) || FError) return 274;
   V[275]=cimag(HggS(pow(V[153]/(2)/(V[59]),2)));
   if(!isfinite(V[275]) || FError) return 275;
   V[276]=creal(HggS(pow(V[153]/(2)/(V[60]),2)));
   if(!isfinite(V[276]) || FError) return 276;
   V[277]=cimag(HggS(pow(V[153]/(2)/(V[60]),2)));
   if(!isfinite(V[277]) || FError) return 277;
   V[278]=creal(HggS(pow(V[153]/(2)/(V[61]),2)));
   if(!isfinite(V[278]) || FError) return 278;
   V[279]=cimag(HggS(pow(V[153]/(2)/(V[61]),2)));
   if(!isfinite(V[279]) || FError) return 279;
   V[280]=creal(HggS(pow(V[153]/(2)/(V[62]),2)));
   if(!isfinite(V[280]) || FError) return 280;
   V[281]=cimag(HggS(pow(V[153]/(2)/(V[62]),2)));
   if(!isfinite(V[281]) || FError) return 281;
   V[282]=creal(HggS(pow(V[153]/(2)/(V[63]),2)));
   if(!isfinite(V[282]) || FError) return 282;
   V[283]=cimag(HggS(pow(V[153]/(2)/(V[63]),2)));
   if(!isfinite(V[283]) || FError) return 283;
   V[284]=creal(HggS(pow(V[153]/(2)/(V[64]),2)));
   if(!isfinite(V[284]) || FError) return 284;
   V[285]=cimag(HggS(pow(V[153]/(2)/(V[64]),2)));
   if(!isfinite(V[285]) || FError) return 285;
   V[286]=creal(HggS(pow(V[153]/(2)/(V[47]),2)));
   if(!isfinite(V[286]) || FError) return 286;
   V[287]=cimag(HggS(pow(V[153]/(2)/(V[47]),2)));
   if(!isfinite(V[287]) || FError) return 287;
   V[288]=creal(HggS(pow(V[153]/(2)/(V[48]),2)));
   if(!isfinite(V[288]) || FError) return 288;
   V[289]=cimag(HggS(pow(V[153]/(2)/(V[48]),2)));
   if(!isfinite(V[289]) || FError) return 289;
   V[290]=creal(HggS(pow(V[153]/(2)/(V[51]),2)));
   if(!isfinite(V[290]) || FError) return 290;
   V[291]=cimag(HggS(pow(V[153]/(2)/(V[51]),2)));
   if(!isfinite(V[291]) || FError) return 291;
   V[292]=creal(HggS(pow(V[153]/(2)/(V[52]),2)));
   if(!isfinite(V[292]) || FError) return 292;
   V[293]=cimag(HggS(pow(V[153]/(2)/(V[52]),2)));
   if(!isfinite(V[293]) || FError) return 293;
   V[294]=creal(HggS(pow(V[153]/(2)/(V[55]),2)));
   if(!isfinite(V[294]) || FError) return 294;
   V[295]=cimag(HggS(pow(V[153]/(2)/(V[55]),2)));
   if(!isfinite(V[295]) || FError) return 295;
   V[296]=creal(HggS(pow(V[153]/(2)/(V[56]),2)));
   if(!isfinite(V[296]) || FError) return 296;
   V[297]=cimag(HggS(pow(V[153]/(2)/(V[56]),2)));
   if(!isfinite(V[297]) || FError) return 297;
   V[298]=creal(HggS(pow(V[153]/(2)/(V[49]),2)));
   if(!isfinite(V[298]) || FError) return 298;
   V[299]=cimag(HggS(pow(V[153]/(2)/(V[49]),2)));
   if(!isfinite(V[299]) || FError) return 299;
   V[300]=creal(HggS(pow(V[153]/(2)/(V[50]),2)));
   if(!isfinite(V[300]) || FError) return 300;
   V[301]=cimag(HggS(pow(V[153]/(2)/(V[50]),2)));
   if(!isfinite(V[301]) || FError) return 301;
   V[302]=creal(HggS(pow(V[153]/(2)/(V[53]),2)));
   if(!isfinite(V[302]) || FError) return 302;
   V[303]=cimag(HggS(pow(V[153]/(2)/(V[53]),2)));
   if(!isfinite(V[303]) || FError) return 303;
   V[304]=creal(HggS(pow(V[153]/(2)/(V[54]),2)));
   if(!isfinite(V[304]) || FError) return 304;
   V[305]=cimag(HggS(pow(V[153]/(2)/(V[54]),2)));
   if(!isfinite(V[305]) || FError) return 305;
   V[306]=creal(HggS(pow(V[153]/(2)/(V[57]),2)));
   if(!isfinite(V[306]) || FError) return 306;
   V[307]=cimag(HggS(pow(V[153]/(2)/(V[57]),2)));
   if(!isfinite(V[307]) || FError) return 307;
   V[308]=creal(HggS(pow(V[153]/(2)/(V[58]),2)));
   if(!isfinite(V[308]) || FError) return 308;
   V[309]=cimag(HggS(pow(V[153]/(2)/(V[58]),2)));
   if(!isfinite(V[309]) || FError) return 309;
   V[310]=creal(HggS(pow(V[153]/(2)/(V[39]),2)));
   if(!isfinite(V[310]) || FError) return 310;
   V[311]=cimag(HggS(pow(V[153]/(2)/(V[39]),2)));
   if(!isfinite(V[311]) || FError) return 311;
   V[312]=creal(HggV(pow(V[153]/(2)/(V[162]),2)));
   if(!isfinite(V[312]) || FError) return 312;
   V[313]=cimag(HggV(pow(V[153]/(2)/(V[162]),2)));
   if(!isfinite(V[313]) || FError) return 313;
   V[314]=creal(HggF(pow(V[153]/(2)/(V[45]),2)));
   if(!isfinite(V[314]) || FError) return 314;
   V[315]=cimag(HggF(pow(V[153]/(2)/(V[45]),2)));
   if(!isfinite(V[315]) || FError) return 315;
   V[316]=creal(HggA(pow(V[153]/(2)/(V[45]),2)));
   if(!isfinite(V[316]) || FError) return 316;
   V[317]=cimag(HggA(pow(V[153]/(2)/(V[45]),2)));
   if(!isfinite(V[317]) || FError) return 317;
   V[318]=creal(HggF(pow(V[153]/(2)/(V[46]),2)));
   if(!isfinite(V[318]) || FError) return 318;
   V[319]=cimag(HggF(pow(V[153]/(2)/(V[46]),2)));
   if(!isfinite(V[319]) || FError) return 319;
   V[320]=creal(HggA(pow(V[153]/(2)/(V[46]),2)));
   if(!isfinite(V[320]) || FError) return 320;
   V[321]=cimag(HggA(pow(V[153]/(2)/(V[46]),2)));
   if(!isfinite(V[321]) || FError) return 321;
   V[322]=V[199]*McRun(V[153]/(2))/(McRun(V[199]));
   if(!isfinite(V[322]) || FError) return 322;
   V[323]=V[198]*MbRun(V[153]/(2))/(MbRun(V[198]));
   if(!isfinite(V[323]) || FError) return 323;
   V[324]=V[6]*MtRun(V[153]/(2))/(MtRun(V[197]));
   if(!isfinite(V[324]) || FError) return 324;
   V[325]=creal(HggF(pow(V[153]/(2)/(V[322]),2))*(1+V[255]*Hgam1F(pow(V[153]/(2)/(V[322]),2))));
   if(!isfinite(V[325]) || FError) return 325;
   V[326]=cimag(HggF(pow(V[153]/(2)/(V[322]),2))*(1+V[255]*Hgam1F(pow(V[153]/(2)/(V[322]),2))));
   if(!isfinite(V[326]) || FError) return 326;
   V[327]=creal(HggA(pow(V[153]/(2)/(V[322]),2))*(1+V[255]*Hgam1A(pow(V[153]/(2)/(V[322]),2))));
   if(!isfinite(V[327]) || FError) return 327;
   V[328]=cimag(HggA(pow(V[153]/(2)/(V[322]),2))*(1+V[255]*Hgam1A(pow(V[153]/(2)/(V[322]),2))));
   if(!isfinite(V[328]) || FError) return 328;
   V[329]=creal(HggF(pow(V[153]/(2)/(V[323]),2))*(1+V[255]*Hgam1F(pow(V[153]/(2)/(V[323]),2))));
   if(!isfinite(V[329]) || FError) return 329;
   V[330]=cimag(HggF(pow(V[153]/(2)/(V[323]),2))*(1+V[255]*Hgam1F(pow(V[153]/(2)/(V[323]),2))));
   if(!isfinite(V[330]) || FError) return 330;
   V[331]=creal(HggA(pow(V[153]/(2)/(V[323]),2))*(1+V[255]*Hgam1A(pow(V[153]/(2)/(V[323]),2))));
   if(!isfinite(V[331]) || FError) return 331;
   V[332]=cimag(HggA(pow(V[153]/(2)/(V[323]),2))*(1+V[255]*Hgam1A(pow(V[153]/(2)/(V[323]),2))));
   if(!isfinite(V[332]) || FError) return 332;
   V[333]=creal(HggF(pow(V[153]/(2)/(V[324]),2))*(1+V[255]*Hgam1F(pow(V[153]/(2)/(V[324]),2))));
   if(!isfinite(V[333]) || FError) return 333;
   V[334]=cimag(HggF(pow(V[153]/(2)/(V[324]),2))*(1+V[255]*Hgam1F(pow(V[153]/(2)/(V[324]),2))));
   if(!isfinite(V[334]) || FError) return 334;
   V[335]=creal(HggA(pow(V[153]/(2)/(V[324]),2))*(1+V[255]*Hgam1A(pow(V[153]/(2)/(V[324]),2))));
   if(!isfinite(V[335]) || FError) return 335;
   V[336]=cimag(HggA(pow(V[153]/(2)/(V[324]),2))*(1+V[255]*Hgam1A(pow(V[153]/(2)/(V[324]),2))));
   if(!isfinite(V[336]) || FError) return 336;
   V[337]=pow(V[153]/(2*V[47]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[47])),6/(double)((23)))),2);
   if(!isfinite(V[337]) || FError) return 337;
   V[338]=creal(HggS(V[337])*(1+V[255]*Hgam1S(V[337])));
   if(!isfinite(V[338]) || FError) return 338;
   V[339]=cimag(HggS(V[337])*(1+V[255]*Hgam1S(V[337])));
   if(!isfinite(V[339]) || FError) return 339;
   V[340]=pow(V[153]/(2*V[48]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[48])),6/(double)((23)))),2);
   if(!isfinite(V[340]) || FError) return 340;
   V[341]=creal(HggS(V[340])*(1+V[255]*Hgam1S(V[340])));
   if(!isfinite(V[341]) || FError) return 341;
   V[342]=cimag(HggS(V[340])*(1+V[255]*Hgam1S(V[340])));
   if(!isfinite(V[342]) || FError) return 342;
   V[343]=pow(V[153]/(2*V[51]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[51])),6/(double)((23)))),2);
   if(!isfinite(V[343]) || FError) return 343;
   V[344]=creal(HggS(V[343])*(1+V[255]*Hgam1S(V[343])));
   if(!isfinite(V[344]) || FError) return 344;
   V[345]=cimag(HggS(V[343])*(1+V[255]*Hgam1S(V[343])));
   if(!isfinite(V[345]) || FError) return 345;
   V[346]=pow(V[153]/(2*V[52]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[52])),6/(double)((23)))),2);
   if(!isfinite(V[346]) || FError) return 346;
   V[347]=creal(HggS(V[346])*(1+V[255]*Hgam1S(V[346])));
   if(!isfinite(V[347]) || FError) return 347;
   V[348]=cimag(HggS(V[346])*(1+V[255]*Hgam1S(V[346])));
   if(!isfinite(V[348]) || FError) return 348;
   V[349]=pow(V[153]/(2*V[55]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[55])),6/(double)((23)))),2);
   if(!isfinite(V[349]) || FError) return 349;
   V[350]=creal(HggS(V[349])*(1+V[255]*Hgam1S(V[349])));
   if(!isfinite(V[350]) || FError) return 350;
   V[351]=cimag(HggS(V[349])*(1+V[255]*Hgam1S(V[349])));
   if(!isfinite(V[351]) || FError) return 351;
   V[352]=pow(V[153]/(2*V[56]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[56])),6/(double)((23)))),2);
   if(!isfinite(V[352]) || FError) return 352;
   V[353]=creal(HggS(V[352])*(1+V[255]*Hgam1S(V[352])));
   if(!isfinite(V[353]) || FError) return 353;
   V[354]=cimag(HggS(V[352])*(1+V[255]*Hgam1S(V[352])));
   if(!isfinite(V[354]) || FError) return 354;
   V[355]=pow(V[153]/(2*V[49]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[49])),6/(double)((23)))),2);
   if(!isfinite(V[355]) || FError) return 355;
   V[356]=creal(HggS(V[355])*(1+V[255]*Hgam1S(V[355])));
   if(!isfinite(V[356]) || FError) return 356;
   V[357]=cimag(HggS(V[355])*(1+V[255]*Hgam1S(V[355])));
   if(!isfinite(V[357]) || FError) return 357;
   V[358]=pow(V[153]/(2*V[50]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[50])),6/(double)((23)))),2);
   if(!isfinite(V[358]) || FError) return 358;
   V[359]=creal(HggS(V[358])*(1+V[255]*Hgam1S(V[358])));
   if(!isfinite(V[359]) || FError) return 359;
   V[360]=cimag(HggS(V[358])*(1+V[255]*Hgam1S(V[358])));
   if(!isfinite(V[360]) || FError) return 360;
   V[361]=pow(V[153]/(2*V[53]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[53])),6/(double)((23)))),2);
   if(!isfinite(V[361]) || FError) return 361;
   V[362]=creal(HggS(V[361])*(1+V[255]*Hgam1S(V[361])));
   if(!isfinite(V[362]) || FError) return 362;
   V[363]=cimag(HggS(V[361])*(1+V[255]*Hgam1S(V[361])));
   if(!isfinite(V[363]) || FError) return 363;
   V[364]=pow(V[153]/(2*V[54]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[54])),6/(double)((23)))),2);
   if(!isfinite(V[364]) || FError) return 364;
   V[365]=creal(HggS(V[364])*(1+V[255]*Hgam1S(V[364])));
   if(!isfinite(V[365]) || FError) return 365;
   V[366]=cimag(HggS(V[364])*(1+V[255]*Hgam1S(V[364])));
   if(!isfinite(V[366]) || FError) return 366;
   V[367]=pow(V[153]/(2*V[57]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[57])),6/(double)((23)))),2);
   if(!isfinite(V[367]) || FError) return 367;
   V[368]=creal(HggS(V[367])*(1+V[255]*Hgam1S(V[367])));
   if(!isfinite(V[368]) || FError) return 368;
   V[369]=cimag(HggS(V[367])*(1+V[255]*Hgam1S(V[367])));
   if(!isfinite(V[369]) || FError) return 369;
   V[370]=pow(V[153]/(2*V[58]*pow(alphaQCD(V[153]/(2))/(alphaQCD(V[58])),6/(double)((23)))),2);
   if(!isfinite(V[370]) || FError) return 370;
   V[371]=creal(HggS(V[370])*(1+V[255]*Hgam1S(V[370])));
   if(!isfinite(V[371]) || FError) return 371;
   V[372]=cimag(HggS(V[370])*(1+V[255]*Hgam1S(V[370])));
   if(!isfinite(V[372]) || FError) return 372;
   V[373]=-V[0]/(V[162])*V[199]/(V[1])*V[114]/(V[164])/(2)/(V[199]);
   if(!isfinite(V[373]) || FError) return 373;
   V[374]=-V[0]/(V[162])*V[198]/(V[1])/(V[165])/(V[164])*(V[164]*V[113]-V[164]*V[172]*V[113]+V[172]*V[114]*V[165])/(2)/(V[198]);
   if(!isfinite(V[374]) || FError) return 374;
   V[375]=-V[0]/(V[162])*V[197]/(V[1])*V[114]/(V[164])/(2)/(V[197]);
   if(!isfinite(V[375]) || FError) return 375;
   V[376]=-V[0]/(V[162])*V[7]/(V[1])*V[113]/(V[165])/(2)/(V[7]);
   if(!isfinite(V[376]) || FError) return 376;
   V[377]=-V[0]/(V[162])*V[199]/(V[1])*V[117]/(V[164])/(2)/(V[199]);
   if(!isfinite(V[377]) || FError) return 377;
   V[378]=-V[0]/(V[162])*V[198]/(V[1])/(V[165])/(V[164])*(V[164]*V[116]-V[164]*V[172]*V[116]+V[172]*V[117]*V[165])/(2)/(V[198]);
   if(!isfinite(V[378]) || FError) return 378;
   V[379]=-V[0]/(V[162])*V[197]/(V[1])*V[117]/(V[164])/(2)/(V[197]);
   if(!isfinite(V[379]) || FError) return 379;
   V[380]=-V[0]/(V[162])*V[7]/(V[1])*V[116]/(V[165])/(2)/(V[7]);
   if(!isfinite(V[380]) || FError) return 380;
   V[381]=-V[0]/(V[162])*V[199]/(V[1])*V[120]/(V[164])/(2)/(V[199]);
   if(!isfinite(V[381]) || FError) return 381;
   V[382]=-V[0]/(V[162])*V[198]/(V[1])/(V[165])/(V[164])*(V[164]*V[119]-V[164]*V[172]*V[119]+V[172]*V[120]*V[165])/(2)/(V[198]);
   if(!isfinite(V[382]) || FError) return 382;
   V[383]=-V[0]/(V[162])*V[197]/(V[1])*V[120]/(V[164])/(2)/(V[197]);
   if(!isfinite(V[383]) || FError) return 383;
   V[384]=-V[0]/(V[162])*V[7]/(V[1])*V[119]/(V[165])/(2)/(V[7]);
   if(!isfinite(V[384]) || FError) return 384;
   V[385]=1/(pow(V[160],2))*V[0]*V[162]/(V[1])*(V[161]*V[113]*V[165]-V[161]*V[164]*V[114])/(2)/(pow(V[59],2))/(2);
   if(!isfinite(V[385]) || FError) return 385;
   V[386]=1/(pow(V[160],2))*V[0]*V[162]*V[1]*(V[113]*V[165]-V[164]*V[114])/(pow(V[60],2))/(2);
   if(!isfinite(V[386]) || FError) return 386;
   V[387]=1/(pow(V[160],2))*V[0]*V[162]/(V[1])*(V[161]*V[113]*V[165]-V[161]*V[164]*V[114])/(2)/(pow(V[61],2))/(2);
   if(!isfinite(V[387]) || FError) return 387;
   V[388]=1/(pow(V[160],2))*V[0]*V[162]*V[1]*(V[113]*V[165]-V[164]*V[114])/(pow(V[62],2))/(2);
   if(!isfinite(V[388]) || FError) return 388;
   V[389]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[208],2)*V[236]+3*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[208],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[208],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[208],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[208],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[209],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[209],2)-6*pow(V[160],2)*V[0]*V[15]*V[9]*V[113]*V[208]*V[209]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[9]*M_SQRT2*V[115]*V[208]*V[209]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[9]*V[114]*V[208]*V[209]+6*pow(V[160],2)*V[0]*pow(V[9],2)*V[114])/(6)/(pow(V[47],2))/(2);
   if(!isfinite(V[389]) || FError) return 389;
   V[390]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[210],2)*V[236]+3*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[210],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[210],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[210],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[210],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[211],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[211],2)-6*pow(V[160],2)*V[0]*V[15]*V[9]*V[113]*V[210]*V[211]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[9]*M_SQRT2*V[115]*V[210]*V[211]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[9]*V[114]*V[210]*V[211]+6*pow(V[160],2)*V[0]*pow(V[9],2)*V[114])/(6)/(pow(V[48],2))/(2);
   if(!isfinite(V[390]) || FError) return 390;
   V[391]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[216],2)*V[114]*V[235]+3*V[164]*V[0]*pow(V[162],2)*pow(V[216],2)*V[113]*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[216],2)*V[113]*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[216],2)*V[114]+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[216],2)*V[114]+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[217],2)*V[113]*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[217],2)*V[114]-6*pow(V[160],2)*V[0]*V[15]*V[199]*V[216]*V[217]*V[113]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[199]*M_SQRT2*V[216]*V[217]*V[115]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[199]*V[216]*V[217]*V[114]+6*pow(V[160],2)*V[0]*pow(V[199],2)*V[114])/(6)/(pow(V[51],2))/(2);
   if(!isfinite(V[391]) || FError) return 391;
   V[392]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[218],2)*V[114]*V[235]+3*V[164]*V[0]*pow(V[162],2)*pow(V[218],2)*V[113]*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[218],2)*V[113]*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[218],2)*V[114]+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[218],2)*V[114]+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[219],2)*V[113]*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[219],2)*V[114]-6*pow(V[160],2)*V[0]*V[15]*V[199]*V[218]*V[219]*V[113]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[199]*M_SQRT2*V[218]*V[219]*V[115]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[199]*V[218]*V[219]*V[114]+6*pow(V[160],2)*V[0]*pow(V[199],2)*V[114])/(6)/(pow(V[52],2))/(2);
   if(!isfinite(V[392]) || FError) return 392;
   V[393]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[212],2)*V[114]*V[165]*V[236]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[212],2)*V[113]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[212],2)*V[113]-3*V[164]*V[0]*pow(V[162],2)*pow(V[212],2)*V[114]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[212],2)*V[114]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[213],2)*V[113]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[213],2)*V[114]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[212]*V[213]*V[114]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[212]*V[213]*V[115]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[212]*V[213]*V[113]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[113])/(6)/(pow(V[49],2))/(2);
   if(!isfinite(V[393]) || FError) return 393;
   V[394]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[214],2)*V[114]*V[165]*V[236]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[214],2)*V[113]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[214],2)*V[113]-3*V[164]*V[0]*pow(V[162],2)*pow(V[214],2)*V[114]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[214],2)*V[114]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[215],2)*V[113]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[215],2)*V[114]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[214]*V[215]*V[114]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[214]*V[215]*V[115]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[214]*V[215]*V[113]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[113])/(6)/(pow(V[50],2))/(2);
   if(!isfinite(V[394]) || FError) return 394;
   V[395]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[220],2)*V[165]*V[235]+3*pow(V[165],2)*V[0]*pow(V[162],2)*V[113]*pow(V[220],2)-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[113]*pow(V[220],2)-3*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[220],2)*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[220],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[113]*pow(V[221],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[221],2)*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[114]*V[220]*V[221]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[115]*V[220]*V[221]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[113]*V[220]*V[221]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[113])/(6)/(pow(V[53],2))/(2);
   if(!isfinite(V[395]) || FError) return 395;
   V[396]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[222],2)*V[165]*V[235]+3*pow(V[165],2)*V[0]*pow(V[162],2)*V[113]*pow(V[222],2)-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[113]*pow(V[222],2)-3*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[222],2)*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[222],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[113]*pow(V[223],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[223],2)*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[114]*V[222]*V[223]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[115]*V[222]*V[223]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[113]*V[222]*V[223]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[113])/(6)/(pow(V[54],2))/(2);
   if(!isfinite(V[396]) || FError) return 396;
   V[397]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(V[161]*pow(V[165],2)*V[0]*pow(V[162],2)*V[113]*pow(V[109],2)-V[161]*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[109],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[113]*pow(V[110],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[110],2)*V[165]+2*pow(V[160],2)*V[0]*V[15]*V[7]*V[114]*V[109]*V[110]*V[167]+2*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[7]*M_SQRT2*V[115]*V[109]*V[110]-2*pow(V[160],2)*V[21]*V[0]*V[7]*V[113]*V[109]*V[110]-2*pow(V[160],2)*V[0]*pow(V[7],2)*V[113])/(2)/(pow(V[63],2))/(2);
   if(!isfinite(V[397]) || FError) return 397;
   V[398]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(V[161]*pow(V[165],2)*V[0]*pow(V[162],2)*V[113]*pow(V[111],2)-V[161]*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[111],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[113]*pow(V[112],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[114]*pow(V[112],2)*V[165]+2*pow(V[160],2)*V[0]*V[15]*V[7]*V[114]*V[111]*V[112]*V[167]+2*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[7]*M_SQRT2*V[115]*V[111]*V[112]-2*pow(V[160],2)*V[21]*V[0]*V[7]*V[113]*V[111]*V[112]-2*pow(V[160],2)*V[0]*pow(V[7],2)*V[113])/(2)/(pow(V[64],2))/(2);
   if(!isfinite(V[398]) || FError) return 398;
   V[399]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[101],2)*V[234]+3*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[101],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[101],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[101],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[101],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[102],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[102],2)-6*pow(V[160],2)*V[0]*V[15]*V[227]*V[113]*V[101]*V[102]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[227]*M_SQRT2*V[115]*V[101]*V[102]*V[165]+6*pow(V[160],2)*V[228]*V[0]*V[227]*V[114]*V[101]*V[102]+6*pow(V[160],2)*V[0]*pow(V[227],2)*V[114])/(6)/(pow(V[55],2))/(2);
   if(!isfinite(V[399]) || FError) return 399;
   V[400]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[103],2)*V[234]+3*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[103],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[103],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[103],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[103],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[113]*pow(V[104],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[114]*pow(V[104],2)-6*pow(V[160],2)*V[0]*V[15]*V[227]*V[113]*V[103]*V[104]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[227]*M_SQRT2*V[115]*V[103]*V[104]*V[165]+6*pow(V[160],2)*V[228]*V[0]*V[227]*V[114]*V[103]*V[104]+6*pow(V[160],2)*V[0]*pow(V[227],2)*V[114])/(6)/(pow(V[56],2))/(2);
   if(!isfinite(V[400]) || FError) return 400;
   V[401]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[105],2)*V[114]*V[165]*V[234]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[105],2)*V[113]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[105],2)*V[113]-3*V[164]*V[0]*pow(V[162],2)*pow(V[105],2)*V[114]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[105],2)*V[114]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[106],2)*V[113]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[106],2)*V[114]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[232]*V[105]*V[106]*V[114]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[232]*M_SQRT2*V[105]*V[106]*V[115]-6*pow(V[160],2)*V[233]*V[0]*V[232]*V[105]*V[106]*V[113]-6*pow(V[160],2)*V[0]*pow(V[232],2)*V[113])/(6)/(pow(V[57],2))/(2);
   if(!isfinite(V[401]) || FError) return 401;
   V[402]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[107],2)*V[114]*V[165]*V[234]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[107],2)*V[113]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[107],2)*V[113]-3*V[164]*V[0]*pow(V[162],2)*pow(V[107],2)*V[114]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[107],2)*V[114]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[108],2)*V[113]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[108],2)*V[114]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[232]*V[107]*V[108]*V[114]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[232]*M_SQRT2*V[107]*V[108]*V[115]-6*pow(V[160],2)*V[233]*V[0]*V[232]*V[107]*V[108]*V[113]-6*pow(V[160],2)*V[0]*pow(V[232],2)*V[113])/(6)/(pow(V[58],2))/(2);
   if(!isfinite(V[402]) || FError) return 402;
   V[403]=1/(pow(V[160],2))*V[0]*V[162]/(V[1])*(V[161]*V[116]*V[165]-V[161]*V[164]*V[117])/(2)/(pow(V[59],2))/(2);
   if(!isfinite(V[403]) || FError) return 403;
   V[404]=1/(pow(V[160],2))*V[0]*V[162]*V[1]*(V[116]*V[165]-V[164]*V[117])/(pow(V[60],2))/(2);
   if(!isfinite(V[404]) || FError) return 404;
   V[405]=1/(pow(V[160],2))*V[0]*V[162]/(V[1])*(V[161]*V[116]*V[165]-V[161]*V[164]*V[117])/(2)/(pow(V[61],2))/(2);
   if(!isfinite(V[405]) || FError) return 405;
   V[406]=1/(pow(V[160],2))*V[0]*V[162]*V[1]*(V[116]*V[165]-V[164]*V[117])/(pow(V[62],2))/(2);
   if(!isfinite(V[406]) || FError) return 406;
   V[407]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[208],2)*V[236]+3*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[208],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[208],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[208],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[208],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[209],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[209],2)-6*pow(V[160],2)*V[0]*V[15]*V[9]*V[116]*V[208]*V[209]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[9]*M_SQRT2*V[118]*V[208]*V[209]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[9]*V[117]*V[208]*V[209]+6*pow(V[160],2)*V[0]*pow(V[9],2)*V[117])/(6)/(pow(V[47],2))/(2);
   if(!isfinite(V[407]) || FError) return 407;
   V[408]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[210],2)*V[236]+3*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[210],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[210],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[210],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[210],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[211],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[211],2)-6*pow(V[160],2)*V[0]*V[15]*V[9]*V[116]*V[210]*V[211]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[9]*M_SQRT2*V[118]*V[210]*V[211]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[9]*V[117]*V[210]*V[211]+6*pow(V[160],2)*V[0]*pow(V[9],2)*V[117])/(6)/(pow(V[48],2))/(2);
   if(!isfinite(V[408]) || FError) return 408;
   V[409]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[216],2)*V[117]*V[235]+3*V[164]*V[0]*pow(V[162],2)*pow(V[216],2)*V[116]*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[216],2)*V[116]*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[216],2)*V[117]+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[216],2)*V[117]+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[217],2)*V[116]*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[217],2)*V[117]-6*pow(V[160],2)*V[0]*V[15]*V[199]*V[216]*V[217]*V[116]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[199]*M_SQRT2*V[216]*V[217]*V[118]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[199]*V[216]*V[217]*V[117]+6*pow(V[160],2)*V[0]*pow(V[199],2)*V[117])/(6)/(pow(V[51],2))/(2);
   if(!isfinite(V[409]) || FError) return 409;
   V[410]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[218],2)*V[117]*V[235]+3*V[164]*V[0]*pow(V[162],2)*pow(V[218],2)*V[116]*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[218],2)*V[116]*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[218],2)*V[117]+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[218],2)*V[117]+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[219],2)*V[116]*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[219],2)*V[117]-6*pow(V[160],2)*V[0]*V[15]*V[199]*V[218]*V[219]*V[116]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[199]*M_SQRT2*V[218]*V[219]*V[118]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[199]*V[218]*V[219]*V[117]+6*pow(V[160],2)*V[0]*pow(V[199],2)*V[117])/(6)/(pow(V[52],2))/(2);
   if(!isfinite(V[410]) || FError) return 410;
   V[411]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[212],2)*V[117]*V[165]*V[236]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[212],2)*V[116]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[212],2)*V[116]-3*V[164]*V[0]*pow(V[162],2)*pow(V[212],2)*V[117]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[212],2)*V[117]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[213],2)*V[116]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[213],2)*V[117]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[212]*V[213]*V[117]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[212]*V[213]*V[118]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[212]*V[213]*V[116]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[116])/(6)/(pow(V[49],2))/(2);
   if(!isfinite(V[411]) || FError) return 411;
   V[412]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[214],2)*V[117]*V[165]*V[236]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[214],2)*V[116]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[214],2)*V[116]-3*V[164]*V[0]*pow(V[162],2)*pow(V[214],2)*V[117]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[214],2)*V[117]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[215],2)*V[116]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[215],2)*V[117]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[214]*V[215]*V[117]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[214]*V[215]*V[118]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[214]*V[215]*V[116]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[116])/(6)/(pow(V[50],2))/(2);
   if(!isfinite(V[412]) || FError) return 412;
   V[413]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[220],2)*V[165]*V[235]+3*pow(V[165],2)*V[0]*pow(V[162],2)*V[116]*pow(V[220],2)-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[116]*pow(V[220],2)-3*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[220],2)*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[220],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[116]*pow(V[221],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[221],2)*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[117]*V[220]*V[221]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[118]*V[220]*V[221]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[116]*V[220]*V[221]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[116])/(6)/(pow(V[53],2))/(2);
   if(!isfinite(V[413]) || FError) return 413;
   V[414]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[222],2)*V[165]*V[235]+3*pow(V[165],2)*V[0]*pow(V[162],2)*V[116]*pow(V[222],2)-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[116]*pow(V[222],2)-3*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[222],2)*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[222],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[116]*pow(V[223],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[223],2)*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[117]*V[222]*V[223]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[118]*V[222]*V[223]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[116]*V[222]*V[223]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[116])/(6)/(pow(V[54],2))/(2);
   if(!isfinite(V[414]) || FError) return 414;
   V[415]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(V[161]*pow(V[165],2)*V[0]*pow(V[162],2)*V[116]*pow(V[109],2)-V[161]*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[109],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[116]*pow(V[110],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[110],2)*V[165]+2*pow(V[160],2)*V[0]*V[15]*V[7]*V[117]*V[109]*V[110]*V[167]+2*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[7]*M_SQRT2*V[118]*V[109]*V[110]-2*pow(V[160],2)*V[21]*V[0]*V[7]*V[116]*V[109]*V[110]-2*pow(V[160],2)*V[0]*pow(V[7],2)*V[116])/(2)/(pow(V[63],2))/(2);
   if(!isfinite(V[415]) || FError) return 415;
   V[416]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(V[161]*pow(V[165],2)*V[0]*pow(V[162],2)*V[116]*pow(V[111],2)-V[161]*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[111],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[116]*pow(V[112],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[117]*pow(V[112],2)*V[165]+2*pow(V[160],2)*V[0]*V[15]*V[7]*V[117]*V[111]*V[112]*V[167]+2*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[7]*M_SQRT2*V[118]*V[111]*V[112]-2*pow(V[160],2)*V[21]*V[0]*V[7]*V[116]*V[111]*V[112]-2*pow(V[160],2)*V[0]*pow(V[7],2)*V[116])/(2)/(pow(V[64],2))/(2);
   if(!isfinite(V[416]) || FError) return 416;
   V[417]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[101],2)*V[234]+3*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[101],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[101],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[101],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[101],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[102],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[102],2)-6*pow(V[160],2)*V[0]*V[15]*V[227]*V[116]*V[101]*V[102]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[227]*M_SQRT2*V[118]*V[101]*V[102]*V[165]+6*pow(V[160],2)*V[228]*V[0]*V[227]*V[117]*V[101]*V[102]+6*pow(V[160],2)*V[0]*pow(V[227],2)*V[117])/(6)/(pow(V[55],2))/(2);
   if(!isfinite(V[417]) || FError) return 417;
   V[418]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[103],2)*V[234]+3*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[103],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[103],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[103],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[103],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[116]*pow(V[104],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[117]*pow(V[104],2)-6*pow(V[160],2)*V[0]*V[15]*V[227]*V[116]*V[103]*V[104]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[227]*M_SQRT2*V[118]*V[103]*V[104]*V[165]+6*pow(V[160],2)*V[228]*V[0]*V[227]*V[117]*V[103]*V[104]+6*pow(V[160],2)*V[0]*pow(V[227],2)*V[117])/(6)/(pow(V[56],2))/(2);
   if(!isfinite(V[418]) || FError) return 418;
   V[419]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[105],2)*V[117]*V[165]*V[234]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[105],2)*V[116]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[105],2)*V[116]-3*V[164]*V[0]*pow(V[162],2)*pow(V[105],2)*V[117]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[105],2)*V[117]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[106],2)*V[116]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[106],2)*V[117]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[232]*V[105]*V[106]*V[117]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[232]*M_SQRT2*V[105]*V[106]*V[118]-6*pow(V[160],2)*V[233]*V[0]*V[232]*V[105]*V[106]*V[116]-6*pow(V[160],2)*V[0]*pow(V[232],2)*V[116])/(6)/(pow(V[57],2))/(2);
   if(!isfinite(V[419]) || FError) return 419;
   V[420]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[107],2)*V[117]*V[165]*V[234]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[107],2)*V[116]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[107],2)*V[116]-3*V[164]*V[0]*pow(V[162],2)*pow(V[107],2)*V[117]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[107],2)*V[117]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[108],2)*V[116]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[108],2)*V[117]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[232]*V[107]*V[108]*V[117]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[232]*M_SQRT2*V[107]*V[108]*V[118]-6*pow(V[160],2)*V[233]*V[0]*V[232]*V[107]*V[108]*V[116]-6*pow(V[160],2)*V[0]*pow(V[232],2)*V[116])/(6)/(pow(V[58],2))/(2);
   if(!isfinite(V[420]) || FError) return 420;
   V[421]=1/(pow(V[160],2))*V[0]*V[162]/(V[1])*(V[161]*V[119]*V[165]-V[161]*V[164]*V[120])/(2)/(pow(V[59],2))/(2);
   if(!isfinite(V[421]) || FError) return 421;
   V[422]=1/(pow(V[160],2))*V[0]*V[162]*V[1]*(V[119]*V[165]-V[164]*V[120])/(pow(V[60],2))/(2);
   if(!isfinite(V[422]) || FError) return 422;
   V[423]=1/(pow(V[160],2))*V[0]*V[162]/(V[1])*(V[161]*V[119]*V[165]-V[161]*V[164]*V[120])/(2)/(pow(V[61],2))/(2);
   if(!isfinite(V[423]) || FError) return 423;
   V[424]=1/(pow(V[160],2))*V[0]*V[162]*V[1]*(V[119]*V[165]-V[164]*V[120])/(pow(V[62],2))/(2);
   if(!isfinite(V[424]) || FError) return 424;
   V[425]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[208],2)*V[236]+3*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[208],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[208],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[208],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[208],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[209],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[209],2)-6*pow(V[160],2)*V[0]*V[15]*V[9]*V[119]*V[208]*V[209]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[9]*M_SQRT2*V[121]*V[208]*V[209]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[9]*V[120]*V[208]*V[209]+6*pow(V[160],2)*V[0]*pow(V[9],2)*V[120])/(6)/(pow(V[47],2))/(2);
   if(!isfinite(V[425]) || FError) return 425;
   V[426]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[210],2)*V[236]+3*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[210],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[210],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[210],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[210],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[211],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[211],2)-6*pow(V[160],2)*V[0]*V[15]*V[9]*V[119]*V[210]*V[211]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[9]*M_SQRT2*V[121]*V[210]*V[211]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[9]*V[120]*V[210]*V[211]+6*pow(V[160],2)*V[0]*pow(V[9],2)*V[120])/(6)/(pow(V[48],2))/(2);
   if(!isfinite(V[426]) || FError) return 426;
   V[427]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[216],2)*V[120]*V[235]+3*V[164]*V[0]*pow(V[162],2)*pow(V[216],2)*V[119]*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[216],2)*V[119]*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[216],2)*V[120]+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[216],2)*V[120]+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[217],2)*V[119]*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[217],2)*V[120]-6*pow(V[160],2)*V[0]*V[15]*V[199]*V[216]*V[217]*V[119]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[199]*M_SQRT2*V[216]*V[217]*V[121]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[199]*V[216]*V[217]*V[120]+6*pow(V[160],2)*V[0]*pow(V[199],2)*V[120])/(6)/(pow(V[51],2))/(2);
   if(!isfinite(V[427]) || FError) return 427;
   V[428]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[218],2)*V[120]*V[235]+3*V[164]*V[0]*pow(V[162],2)*pow(V[218],2)*V[119]*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[218],2)*V[119]*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[218],2)*V[120]+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[218],2)*V[120]+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[219],2)*V[119]*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*pow(V[219],2)*V[120]-6*pow(V[160],2)*V[0]*V[15]*V[199]*V[218]*V[219]*V[119]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[199]*M_SQRT2*V[218]*V[219]*V[121]*V[165]+6*pow(V[160],2)*V[154]*V[0]*V[199]*V[218]*V[219]*V[120]+6*pow(V[160],2)*V[0]*pow(V[199],2)*V[120])/(6)/(pow(V[52],2))/(2);
   if(!isfinite(V[428]) || FError) return 428;
   V[429]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[212],2)*V[120]*V[165]*V[236]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[212],2)*V[119]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[212],2)*V[119]-3*V[164]*V[0]*pow(V[162],2)*pow(V[212],2)*V[120]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[212],2)*V[120]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[213],2)*V[119]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[213],2)*V[120]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[212]*V[213]*V[120]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[212]*V[213]*V[121]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[212]*V[213]*V[119]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[119])/(6)/(pow(V[49],2))/(2);
   if(!isfinite(V[429]) || FError) return 429;
   V[430]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[214],2)*V[120]*V[165]*V[236]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[214],2)*V[119]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[214],2)*V[119]-3*V[164]*V[0]*pow(V[162],2)*pow(V[214],2)*V[120]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[214],2)*V[120]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[215],2)*V[119]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[215],2)*V[120]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[214]*V[215]*V[120]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[214]*V[215]*V[121]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[214]*V[215]*V[119]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[119])/(6)/(pow(V[50],2))/(2);
   if(!isfinite(V[430]) || FError) return 430;
   V[431]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[220],2)*V[165]*V[235]+3*pow(V[165],2)*V[0]*pow(V[162],2)*V[119]*pow(V[220],2)-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[119]*pow(V[220],2)-3*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[220],2)*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[220],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[119]*pow(V[221],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[221],2)*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[120]*V[220]*V[221]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[121]*V[220]*V[221]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[119]*V[220]*V[221]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[119])/(6)/(pow(V[53],2))/(2);
   if(!isfinite(V[431]) || FError) return 431;
   V[432]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[222],2)*V[165]*V[235]+3*pow(V[165],2)*V[0]*pow(V[162],2)*V[119]*pow(V[222],2)-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[119]*pow(V[222],2)-3*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[222],2)*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[222],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[119]*pow(V[223],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[223],2)*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[9]*V[120]*V[222]*V[223]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[9]*M_SQRT2*V[121]*V[222]*V[223]-6*pow(V[160],2)*V[155]*V[0]*V[9]*V[119]*V[222]*V[223]-6*pow(V[160],2)*V[0]*pow(V[9],2)*V[119])/(6)/(pow(V[54],2))/(2);
   if(!isfinite(V[432]) || FError) return 432;
   V[433]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(V[161]*pow(V[165],2)*V[0]*pow(V[162],2)*V[119]*pow(V[109],2)-V[161]*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[109],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[119]*pow(V[110],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[110],2)*V[165]+2*pow(V[160],2)*V[0]*V[15]*V[7]*V[120]*V[109]*V[110]*V[167]+2*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[7]*M_SQRT2*V[121]*V[109]*V[110]-2*pow(V[160],2)*V[21]*V[0]*V[7]*V[119]*V[109]*V[110]-2*pow(V[160],2)*V[0]*pow(V[7],2)*V[119])/(2)/(pow(V[63],2))/(2);
   if(!isfinite(V[433]) || FError) return 433;
   V[434]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(V[161]*pow(V[165],2)*V[0]*pow(V[162],2)*V[119]*pow(V[111],2)-V[161]*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[111],2)*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*V[119]*pow(V[112],2)-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[120]*pow(V[112],2)*V[165]+2*pow(V[160],2)*V[0]*V[15]*V[7]*V[120]*V[111]*V[112]*V[167]+2*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[7]*M_SQRT2*V[121]*V[111]*V[112]-2*pow(V[160],2)*V[21]*V[0]*V[7]*V[119]*V[111]*V[112]-2*pow(V[160],2)*V[0]*pow(V[7],2)*V[119])/(2)/(pow(V[64],2))/(2);
   if(!isfinite(V[434]) || FError) return 434;
   V[435]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[101],2)*V[234]+3*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[101],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[101],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[101],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[101],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[102],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[102],2)-6*pow(V[160],2)*V[0]*V[15]*V[227]*V[119]*V[101]*V[102]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[227]*M_SQRT2*V[121]*V[101]*V[102]*V[165]+6*pow(V[160],2)*V[228]*V[0]*V[227]*V[120]*V[101]*V[102]+6*pow(V[160],2)*V[0]*pow(V[227],2)*V[120])/(6)/(pow(V[55],2))/(2);
   if(!isfinite(V[435]) || FError) return 435;
   V[436]=-1/(pow(V[160],2))/(V[162])/(V[1])/(V[164])*(3*pow(V[160],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[103],2)*V[234]+3*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[103],2)*V[165]-4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[103],2)*V[165]-3*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[103],2)+4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[103],2)+4*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*V[119]*pow(V[104],2)*V[165]-4*pow(V[1],2)*pow(V[164],2)*V[0]*pow(V[162],2)*V[120]*pow(V[104],2)-6*pow(V[160],2)*V[0]*V[15]*V[227]*V[119]*V[103]*V[104]*V[167]-6*pow(V[160],2)*V[1]*V[15]*V[162]*V[227]*M_SQRT2*V[121]*V[103]*V[104]*V[165]+6*pow(V[160],2)*V[228]*V[0]*V[227]*V[120]*V[103]*V[104]+6*pow(V[160],2)*V[0]*pow(V[227],2)*V[120])/(6)/(pow(V[56],2))/(2);
   if(!isfinite(V[436]) || FError) return 436;
   V[437]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[105],2)*V[120]*V[165]*V[234]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[105],2)*V[119]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[105],2)*V[119]-3*V[164]*V[0]*pow(V[162],2)*pow(V[105],2)*V[120]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[105],2)*V[120]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[106],2)*V[119]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[106],2)*V[120]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[232]*V[105]*V[106]*V[120]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[232]*M_SQRT2*V[105]*V[106]*V[121]-6*pow(V[160],2)*V[233]*V[0]*V[232]*V[105]*V[106]*V[119]-6*pow(V[160],2)*V[0]*pow(V[232],2)*V[119])/(6)/(pow(V[57],2))/(2);
   if(!isfinite(V[437]) || FError) return 437;
   V[438]=1/(pow(V[160],2))/(V[162])/(V[1])/(V[165])*(3*pow(V[160],2)*V[164]*V[0]*pow(V[162],2)*pow(V[107],2)*V[120]*V[165]*V[234]+3*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[107],2)*V[119]-2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[107],2)*V[119]-3*V[164]*V[0]*pow(V[162],2)*pow(V[107],2)*V[120]*V[165]+2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[107],2)*V[120]*V[165]+2*pow(V[1],2)*pow(V[165],2)*V[0]*pow(V[162],2)*pow(V[108],2)*V[119]-2*pow(V[1],2)*V[164]*V[0]*pow(V[162],2)*pow(V[108],2)*V[120]*V[165]+6*pow(V[160],2)*V[0]*V[15]*V[232]*V[107]*V[108]*V[120]*V[167]+6*pow(V[160],2)*V[1]*V[164]*V[15]*V[162]*V[232]*M_SQRT2*V[107]*V[108]*V[121]-6*pow(V[160],2)*V[233]*V[0]*V[232]*V[107]*V[108]*V[119]-6*pow(V[160],2)*V[0]*pow(V[232],2)*V[119])/(6)/(pow(V[58],2))/(2);
   if(!isfinite(V[438]) || FError) return 438;
   V[439]=V[0]*V[162]/(V[1])*(V[113]*V[165]+V[164]*V[114])/(pow(V[162],2))/(2);
   if(!isfinite(V[439]) || FError) return 439;
   V[440]=2*V[1]*pow(V[164],2)*V[162]*V[113]*V[165]*V[128]+2*V[1]*pow(V[165],2)*V[164]*V[162]*V[114]*V[129]+2*V[1]*pow(V[165],2)*V[162]*V[113]*V[165]*V[130]+2*V[1]*pow(V[164],3)*V[162]*V[114]*V[130]-2*V[1]*pow(V[165],2)*V[164]*V[162]*V[114]*V[131];
   if(!isfinite(V[440]) || FError) return 440;
   V[441]=V[440]-2*V[1]*pow(V[164],2)*V[162]*V[113]*V[165]*V[131]-2*V[1]*pow(V[165],2)*V[164]*V[162]*V[114]*V[132]-2*V[1]*pow(V[164],2)*V[162]*V[113]*V[165]*V[132]-6*V[1]*pow(V[164],3)*V[162]*V[113]*V[133]+4*V[1]*V[164]*V[162]*V[113]*V[133];
   if(!isfinite(V[441]) || FError) return 441;
   V[442]=V[441]-2*V[1]*pow(V[164],2)*V[162]*V[114]*V[165]*V[133]+6*V[1]*pow(V[164],2)*V[162]*V[114]*V[165]*V[134]-2*V[1]*V[162]*V[114]*V[165]*V[134]-2*V[1]*pow(V[165],2)*V[164]*V[162]*V[113]*V[134]+pow(V[164],2)*V[0]*M_SQRT2*V[115]*V[135]*V[167];
   if(!isfinite(V[442]) || FError) return 442;
   V[443]=V[442]+pow(V[165],2)*V[0]*M_SQRT2*V[115]*V[136]*V[167]+2*V[164]*V[0]*M_SQRT2*V[115]*V[165]*V[138]*V[167]+2*V[164]*V[0]*M_SQRT2*V[115]*V[165]*V[139]*V[167]+2*V[164]*V[0]*M_SQRT2*V[115]*V[165]*V[140]*V[167]+pow(V[164],2)*V[0]*M_SQRT2*V[115]*V[143];
   if(!isfinite(V[443]) || FError) return 443;
   V[444]=V[443]+pow(V[165],2)*V[0]*M_SQRT2*V[115]*V[144]+V[164]*V[0]*M_SQRT2*V[115]*V[147]*V[165]+V[164]*V[0]*M_SQRT2*V[115]*V[148]*V[165];
   if(!isfinite(V[444]) || FError) return 444;
   V[445]=-1/(V[0])*V[444]/(pow(V[39],2))/(2);
   if(!isfinite(V[445]) || FError) return 445;
   V[446]=V[0]*V[162]/(V[1])*(V[116]*V[165]+V[164]*V[117])/(pow(V[162],2))/(2);
   if(!isfinite(V[446]) || FError) return 446;
   V[447]=2*V[1]*pow(V[164],2)*V[162]*V[116]*V[165]*V[128]+2*V[1]*pow(V[165],2)*V[164]*V[162]*V[117]*V[129]+2*V[1]*pow(V[165],2)*V[162]*V[116]*V[165]*V[130]+2*V[1]*pow(V[164],3)*V[162]*V[117]*V[130]-2*V[1]*pow(V[165],2)*V[164]*V[162]*V[117]*V[131];
   if(!isfinite(V[447]) || FError) return 447;
   V[448]=V[447]-2*V[1]*pow(V[164],2)*V[162]*V[116]*V[165]*V[131]-2*V[1]*pow(V[165],2)*V[164]*V[162]*V[117]*V[132]-2*V[1]*pow(V[164],2)*V[162]*V[116]*V[165]*V[132]-6*V[1]*pow(V[164],3)*V[162]*V[116]*V[133]+4*V[1]*V[164]*V[162]*V[116]*V[133];
   if(!isfinite(V[448]) || FError) return 448;
   V[449]=V[448]-2*V[1]*pow(V[164],2)*V[162]*V[117]*V[165]*V[133]+6*V[1]*pow(V[164],2)*V[162]*V[117]*V[165]*V[134]-2*V[1]*V[162]*V[117]*V[165]*V[134]-2*V[1]*pow(V[165],2)*V[164]*V[162]*V[116]*V[134]+pow(V[164],2)*V[0]*M_SQRT2*V[118]*V[135]*V[167];
   if(!isfinite(V[449]) || FError) return 449;
   V[450]=V[449]+pow(V[165],2)*V[0]*M_SQRT2*V[118]*V[136]*V[167]+2*V[164]*V[0]*M_SQRT2*V[118]*V[165]*V[138]*V[167]+2*V[164]*V[0]*M_SQRT2*V[118]*V[165]*V[139]*V[167]+2*V[164]*V[0]*M_SQRT2*V[118]*V[165]*V[140]*V[167]+pow(V[164],2)*V[0]*M_SQRT2*V[118]*V[143];
   if(!isfinite(V[450]) || FError) return 450;
   V[451]=V[450]+pow(V[165],2)*V[0]*M_SQRT2*V[118]*V[144]+V[164]*V[0]*M_SQRT2*V[118]*V[147]*V[165]+V[164]*V[0]*M_SQRT2*V[118]*V[148]*V[165];
   if(!isfinite(V[451]) || FError) return 451;
   V[452]=-1/(V[0])*V[451]/(pow(V[39],2))/(2);
   if(!isfinite(V[452]) || FError) return 452;
   V[453]=V[0]*V[162]/(V[1])*(V[119]*V[165]+V[164]*V[120])/(pow(V[162],2))/(2);
   if(!isfinite(V[453]) || FError) return 453;
   V[454]=2*V[1]*pow(V[164],2)*V[162]*V[119]*V[165]*V[128]+2*V[1]*pow(V[165],2)*V[164]*V[162]*V[120]*V[129]+2*V[1]*pow(V[165],2)*V[162]*V[119]*V[165]*V[130]+2*V[1]*pow(V[164],3)*V[162]*V[120]*V[130]-2*V[1]*pow(V[165],2)*V[164]*V[162]*V[120]*V[131];
   if(!isfinite(V[454]) || FError) return 454;
   V[455]=V[454]-2*V[1]*pow(V[164],2)*V[162]*V[119]*V[165]*V[131]-2*V[1]*pow(V[165],2)*V[164]*V[162]*V[120]*V[132]-2*V[1]*pow(V[164],2)*V[162]*V[119]*V[165]*V[132]-6*V[1]*pow(V[164],3)*V[162]*V[119]*V[133]+4*V[1]*V[164]*V[162]*V[119]*V[133];
   if(!isfinite(V[455]) || FError) return 455;
   V[456]=V[455]-2*V[1]*pow(V[164],2)*V[162]*V[120]*V[165]*V[133]+6*V[1]*pow(V[164],2)*V[162]*V[120]*V[165]*V[134]-2*V[1]*V[162]*V[120]*V[165]*V[134]-2*V[1]*pow(V[165],2)*V[164]*V[162]*V[119]*V[134]+pow(V[164],2)*V[0]*M_SQRT2*V[121]*V[135]*V[167];
   if(!isfinite(V[456]) || FError) return 456;
   V[457]=V[456]+pow(V[165],2)*V[0]*M_SQRT2*V[121]*V[136]*V[167]+2*V[164]*V[0]*M_SQRT2*V[121]*V[165]*V[138]*V[167]+2*V[164]*V[0]*M_SQRT2*V[121]*V[165]*V[139]*V[167]+2*V[164]*V[0]*M_SQRT2*V[121]*V[165]*V[140]*V[167]+pow(V[164],2)*V[0]*M_SQRT2*V[121]*V[143];
   if(!isfinite(V[457]) || FError) return 457;
   V[458]=V[457]+pow(V[165],2)*V[0]*M_SQRT2*V[121]*V[144]+V[164]*V[0]*M_SQRT2*V[121]*V[147]*V[165]+V[164]*V[0]*M_SQRT2*V[121]*V[148]*V[165];
   if(!isfinite(V[458]) || FError) return 458;
   V[459]=-1/(V[0])*V[458]/(pow(V[39],2))/(2);
   if(!isfinite(V[459]) || FError) return 459;
   V[460]=1/(V[162])/(V[1])*M_SQRT2/(V[165])/(V[164])/(V[167])*(V[1]*V[164]*V[162]*V[247]*V[115]*V[94]*V[98]*V[165]-V[164]*V[0]*V[243]*V[113]*V[94]*V[97]*V[167]+V[0]*V[244]*V[114]*V[93]*V[98]*V[165]*V[167])/(2)/(V[45]);
   if(!isfinite(V[460]) || FError) return 460;
   V[461]=1/(V[162])/(V[1])*M_SQRT2/(V[165])/(V[164])/(V[167])*(V[1]*V[164]*V[162]*V[247]*V[115]*V[96]*V[100]*V[165]-V[164]*V[0]*V[243]*V[113]*V[96]*V[99]*V[167]+V[0]*V[244]*V[114]*V[95]*V[100]*V[165]*V[167])/(2)/(V[46]);
   if(!isfinite(V[461]) || FError) return 461;
   V[462]=1/(V[162])/(V[1])*M_SQRT2/(V[165])/(V[164])/(V[167])*(V[1]*V[164]*V[162]*V[247]*V[118]*V[94]*V[98]*V[165]-V[164]*V[0]*V[243]*V[116]*V[94]*V[97]*V[167]+V[0]*V[244]*V[117]*V[93]*V[98]*V[165]*V[167])/(2)/(V[45]);
   if(!isfinite(V[462]) || FError) return 462;
   V[463]=1/(V[162])/(V[1])*M_SQRT2/(V[165])/(V[164])/(V[167])*(V[1]*V[164]*V[162]*V[247]*V[118]*V[96]*V[100]*V[165]-V[164]*V[0]*V[243]*V[116]*V[96]*V[99]*V[167]+V[0]*V[244]*V[117]*V[95]*V[100]*V[165]*V[167])/(2)/(V[46]);
   if(!isfinite(V[463]) || FError) return 463;
   V[464]=1/(V[162])/(V[1])*M_SQRT2/(V[165])/(V[164])/(V[167])*(V[1]*V[164]*V[162]*V[247]*V[121]*V[94]*V[98]*V[165]-V[164]*V[0]*V[243]*V[119]*V[94]*V[97]*V[167]+V[0]*V[244]*V[120]*V[93]*V[98]*V[165]*V[167])/(2)/(V[45]);
   if(!isfinite(V[464]) || FError) return 464;
   V[465]=1/(V[162])/(V[1])*M_SQRT2/(V[165])/(V[164])/(V[167])*(V[1]*V[164]*V[162]*V[247]*V[121]*V[96]*V[100]*V[165]-V[164]*V[0]*V[243]*V[119]*V[96]*V[99]*V[167]+V[0]*V[244]*V[120]*V[95]*V[100]*V[165]*V[167])/(2)/(V[46]);
   if(!isfinite(V[465]) || FError) return 465;
   V[466]=V[0]/(V[162])*V[199]*V[170]/(V[1])/(V[158])/(2)/(V[199])/(2);
   if(!isfinite(V[466]) || FError) return 466;
   V[467]=-V[0]/(V[162])*V[198]*V[170]/(V[1])/(V[165])/(V[164])*(V[172]-pow(V[164],2))/(2)/(V[198])/(2);
   if(!isfinite(V[467]) || FError) return 467;
   V[468]=V[0]/(V[162])*V[197]*V[170]/(V[1])/(V[158])/(2)/(V[197])/(2);
   if(!isfinite(V[468]) || FError) return 468;
   V[469]=V[0]/(V[162])*V[7]*V[170]/(V[1])*V[158]/(2)/(V[7])/(2);
   if(!isfinite(V[469]) || FError) return 469;
   V[470]=V[0]/(V[162])*V[199]*V[171]/(V[1])/(V[158])/(2)/(V[199])/(2);
   if(!isfinite(V[470]) || FError) return 470;
   V[471]=-V[0]/(V[162])*V[198]*V[171]/(V[1])/(V[165])/(V[164])*(V[172]-pow(V[164],2))/(2)/(V[198])/(2);
   if(!isfinite(V[471]) || FError) return 471;
   V[472]=V[0]/(V[162])*V[197]*V[171]/(V[1])/(V[158])/(2)/(V[197])/(2);
   if(!isfinite(V[472]) || FError) return 472;
   V[473]=V[0]/(V[162])*V[7]*V[171]/(V[1])*V[158]/(2)/(V[7])/(2);
   if(!isfinite(V[473]) || FError) return 473;
   V[474]=1/(V[162])/(V[1])*M_SQRT2/(V[165])/(V[164])/(V[167])*(pow(V[165],2)*V[0]*V[244]*V[170]*V[93]*V[98]*V[167]-pow(V[164],2)*V[0]*V[243]*V[170]*V[94]*V[97]*V[167]-V[1]*V[164]*V[162]*V[247]*V[168]*V[94]*V[98]*V[165])/(2)/(V[45])/(2);
   if(!isfinite(V[474]) || FError) return 474;
   V[475]=1/(V[162])/(V[1])*M_SQRT2/(V[165])/(V[164])/(V[167])*(pow(V[165],2)*V[0]*V[244]*V[170]*V[95]*V[100]*V[167]-pow(V[164],2)*V[0]*V[243]*V[170]*V[96]*V[99]*V[167]-V[1]*V[164]*V[162]*V[247]*V[168]*V[96]*V[100]*V[165])/(2)/(V[46])/(2);
   if(!isfinite(V[475]) || FError) return 475;
   V[476]=1/(V[162])/(V[1])*M_SQRT2/(V[165])/(V[164])/(V[167])*(pow(V[165],2)*V[0]*V[244]*V[171]*V[93]*V[98]*V[167]-pow(V[164],2)*V[0]*V[243]*V[171]*V[94]*V[97]*V[167]-V[1]*V[164]*V[162]*V[247]*V[169]*V[94]*V[98]*V[165])/(2)/(V[45])/(2);
   if(!isfinite(V[476]) || FError) return 476;
   V[477]=1/(V[162])/(V[1])*M_SQRT2/(V[165])/(V[164])/(V[167])*(pow(V[165],2)*V[0]*V[244]*V[171]*V[95]*V[100]*V[167]-pow(V[164],2)*V[0]*V[243]*V[171]*V[96]*V[99]*V[167]-V[1]*V[164]*V[162]*V[247]*V[169]*V[96]*V[100]*V[165])/(2)/(V[46])/(2);
   if(!isfinite(V[477]) || FError) return 477;
   V[478]=1+V[255]*(149/(double)((12))+V[255]*(68.6482-V[255]*212.447));
   if(!isfinite(V[478]) || FError) return 478;
   V[479]=2*log(V[6]/(V[153]));
   if(!isfinite(V[479]) || FError) return 479;
   V[480]=1+V[255]*(11/(double)((4))+V[255]*(6.1537-2.8542*V[479]+V[255]*(10.999-17.93*V[479]+5.47*pow(V[479],2))));
   if(!isfinite(V[480]) || FError) return 480;
   V[481]=1+11/(double)((4))*V[255];
   if(!isfinite(V[481]) || FError) return 481;
   V[482]=1+9/(double)((2))*V[255];
   if(!isfinite(V[482]) || FError) return 482;
   V[483]=1/(137.036);
   if(!isfinite(V[483]) || FError) return 483;
   V[484]=2/(double)((3));
   if(!isfinite(V[484]) || FError) return 484;
   V[485]=-1/(double)((3));
   if(!isfinite(V[485]) || FError) return 485;
   V[486]=1+V[255]*(221/(double)((12))+V[255]*(171.5-5*V[479]));
   if(!isfinite(V[486]) || FError) return 486;
   V[487]=-V[255]/(8)*1/(double)((2))*sqrt(V[478])*cabs(((V[262]+I*V[263])*V[374]+(V[258]+I*V[259])*V[373])*V[481]+(V[266]+I*V[267])*V[375]*V[480]+((V[286]+I*V[287])*V[389]+(V[288]+I*V[289])*V[390]+(V[298]+I*V[299])*V[393]+(V[300]+I*V[301])*V[394]+(V[302]+I*V[303])*V[395]+(V[304]+I*V[305])*V[396]+(V[290]+I*V[291])*V[391]+(V[292]+I*V[293])*V[392]+(V[306]+I*V[307])*V[401]+(V[308]+I*V[309])*V[402]+(V[294]+I*V[295])*V[399]+(V[296]+I*V[297])*V[400])*V[482]);
   if(!isfinite(V[487]) || FError) return 487;
   V[488]=-V[483]/(8*V[254])*cabs(3*pow(V[485],2)*(V[329]+I*V[330])*V[374]+3*pow(V[484],2)*((V[333]+I*V[334])*V[375]+(V[325]+I*V[326])*V[373])+3*pow(V[484],2)*((V[338]+I*V[339])*V[389]+(V[341]+I*V[342])*V[390]+(V[344]+I*V[345])*V[391]+(V[347]+I*V[348])*V[392]+(V[350]+I*V[351])*V[399]+(V[353]+I*V[354])*V[400])+3*pow(V[485],2)*((V[356]+I*V[357])*V[393]+(V[359]+I*V[360])*V[394]+(V[362]+I*V[363])*V[395]+(V[365]+I*V[366])*V[396]+(V[368]+I*V[369])*V[401]+(V[371]+I*V[372])*V[402])+(V[274]+I*V[275])*V[385]+(V[276]+I*V[277])*V[386]+(V[278]+I*V[279])*V[387]+(V[280]+I*V[281])*V[388]+(V[282]+I*V[283])*V[397]+(V[284]+I*V[285])*V[398]+(V[270]+I*V[271])*V[376]+(V[314]+I*V[315])*V[460]+(V[318]+I*V[319])*V[461]-(V[312]+I*V[313])*V[439]+(V[310]+I*V[311])*V[445]);
   if(!isfinite(V[488]) || FError) return 488;
   V[489]=-V[255]/(8)*1/(double)((2))*sqrt(V[478])*cabs(((V[262]+I*V[263])*V[378]+(V[258]+I*V[259])*V[377])*V[481]+(V[266]+I*V[267])*V[379]*V[480]+((V[286]+I*V[287])*V[407]+(V[288]+I*V[289])*V[408]+(V[298]+I*V[299])*V[411]+(V[300]+I*V[301])*V[412]+(V[302]+I*V[303])*V[413]+(V[304]+I*V[305])*V[414]+(V[290]+I*V[291])*V[409]+(V[292]+I*V[293])*V[410]+(V[306]+I*V[307])*V[419]+(V[308]+I*V[309])*V[420]+(V[294]+I*V[295])*V[417]+(V[296]+I*V[297])*V[418])*V[482]);
   if(!isfinite(V[489]) || FError) return 489;
   V[490]=-V[483]/(8*V[254])*cabs(3*pow(V[485],2)*(V[329]+I*V[330])*V[378]+3*pow(V[484],2)*((V[333]+I*V[334])*V[379]+(V[325]+I*V[326])*V[377])+3*pow(V[484],2)*((V[338]+I*V[339])*V[407]+(V[341]+I*V[342])*V[408]+(V[344]+I*V[345])*V[409]+(V[347]+I*V[348])*V[410]+(V[350]+I*V[351])*V[417]+(V[353]+I*V[354])*V[418])+3*pow(V[485],2)*((V[356]+I*V[357])*V[411]+(V[359]+I*V[360])*V[412]+(V[362]+I*V[363])*V[413]+(V[365]+I*V[366])*V[414]+(V[368]+I*V[369])*V[419]+(V[371]+I*V[372])*V[420])+(V[274]+I*V[275])*V[403]+(V[276]+I*V[277])*V[404]+(V[278]+I*V[279])*V[405]+(V[280]+I*V[281])*V[406]+(V[282]+I*V[283])*V[415]+(V[284]+I*V[285])*V[416]+(V[270]+I*V[271])*V[380]+(V[314]+I*V[315])*V[462]+(V[318]+I*V[319])*V[463]-(V[312]+I*V[313])*V[446]+(V[310]+I*V[311])*V[452]);
   if(!isfinite(V[490]) || FError) return 490;
   V[491]=-V[255]/(8)*1/(double)((2))*sqrt(V[478])*cabs(((V[262]+I*V[263])*V[382]+(V[258]+I*V[259])*V[381])*V[481]+(V[266]+I*V[267])*V[383]*V[480]+((V[286]+I*V[287])*V[425]+(V[288]+I*V[289])*V[426]+(V[298]+I*V[299])*V[429]+(V[300]+I*V[301])*V[430]+(V[302]+I*V[303])*V[431]+(V[304]+I*V[305])*V[432]+(V[290]+I*V[291])*V[427]+(V[292]+I*V[293])*V[428]+(V[306]+I*V[307])*V[437]+(V[308]+I*V[309])*V[438]+(V[294]+I*V[295])*V[435]+(V[296]+I*V[297])*V[436])*V[482]);
   if(!isfinite(V[491]) || FError) return 491;
   V[492]=-V[483]/(8*V[254])*cabs(3*pow(V[485],2)*(V[329]+I*V[330])*V[382]+3*pow(V[484],2)*((V[333]+I*V[334])*V[383]+(V[325]+I*V[326])*V[381])+3*pow(V[484],2)*((V[338]+I*V[339])*V[425]+(V[341]+I*V[342])*V[426]+(V[344]+I*V[345])*V[427]+(V[347]+I*V[348])*V[428]+(V[350]+I*V[351])*V[435]+(V[353]+I*V[354])*V[436])+3*pow(V[485],2)*((V[356]+I*V[357])*V[429]+(V[359]+I*V[360])*V[430]+(V[362]+I*V[363])*V[431]+(V[365]+I*V[366])*V[432]+(V[368]+I*V[369])*V[437]+(V[371]+I*V[372])*V[438])+(V[274]+I*V[275])*V[421]+(V[276]+I*V[277])*V[422]+(V[278]+I*V[279])*V[423]+(V[280]+I*V[281])*V[424]+(V[282]+I*V[283])*V[433]+(V[284]+I*V[285])*V[434]+(V[270]+I*V[271])*V[384]+(V[314]+I*V[315])*V[464]+(V[318]+I*V[319])*V[465]-(V[312]+I*V[313])*V[453]+(V[310]+I*V[311])*V[459]);
   if(!isfinite(V[492]) || FError) return 492;
   V[493]=-V[255]/(8)*1/(double)((2))*sqrt(V[486])*cabs((V[260]+I*V[261])*V[466]+(V[264]+I*V[265])*V[467]+(V[268]+I*V[269])*V[468]);
   if(!isfinite(V[493]) || FError) return 493;
   V[494]=-V[483]/(8*V[254])*cabs(3*pow(V[485],2)*(V[331]+I*V[332])*V[467]+3*pow(V[484],2)*((V[335]+I*V[336])*V[468]+(V[327]+I*V[328])*V[466])+(V[272]+I*V[273])*V[469]+(V[316]+I*V[317])*V[474]+(V[320]+I*V[321])*V[475]);
   if(!isfinite(V[494]) || FError) return 494;
   V[495]=-V[255]/(8)*1/(double)((2))*sqrt(V[486])*cabs((V[260]+I*V[261])*V[470]+(V[264]+I*V[265])*V[471]+(V[268]+I*V[269])*V[472]);
   if(!isfinite(V[495]) || FError) return 495;
   V[496]=-V[483]/(8*V[254])*cabs(3*pow(V[485],2)*(V[331]+I*V[332])*V[471]+3*pow(V[484],2)*((V[335]+I*V[336])*V[472]+(V[327]+I*V[328])*V[470])+(V[272]+I*V[273])*V[473]+(V[316]+I*V[317])*V[476]+(V[320]+I*V[321])*V[477]);
   if(!isfinite(V[496]) || FError) return 496;
   V[497]=-V[0]/(V[162])/(V[1])/(2);
   if(!isfinite(V[497]) || FError) return 497;
   V[498]=V[0]/(V[1])/(V[162])/(2);
   if(!isfinite(V[498]) || FError) return 498;
   V[499]=-V[255]/(8)*1/(double)((2))*sqrt(V[478])*fabs(V[497])*cabs((V[262]+I*V[263]+V[258]+I*V[259])*V[481]+(V[266]+I*V[267])*V[480]);
   if(!isfinite(V[499]) || FError) return 499;
   V[500]=-V[483]/(8*V[254])*cabs(3*pow(V[485],2)*(V[329]+I*V[330])*V[497]+3*pow(V[484],2)*((V[333]+I*V[334])*V[497]+(V[325]+I*V[326])*V[497])+(V[270]+I*V[271])*V[497]-(V[312]+I*V[313])*V[498]);
   if(!isfinite(V[500]) || FError) return 500;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   return 0;
}
