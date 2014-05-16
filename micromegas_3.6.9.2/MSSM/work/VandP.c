#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../CalcHEP_src/include/extern.h"
#include "../../CalcHEP_src/include/VandP.h"
#include "autoprot.h"
extern int  FError;
/*  Special model functions  */

int nModelParticles=48;
static ModelPrtclsStr ModelPrtcls_[48]=
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
, {"h","h", 25, "Mh","wh",0,1,0}
, {"H","H", 35, "MHH","wHh",0,1,0}
, {"H3","H3", 36, "MH3","wH3",0,1,0}
, {"H+","H-", 37, "MHc","wHc",0,1,3}
, {"~1+","~1-", 1000024, "MC1","wC1",1,1,3}
, {"~2+","~2-", 1000037, "MC2","wC2",1,1,3}
, {"~o1","~o1", 1000022, "MNE1","0",1,1,0}
, {"~o2","~o2", 1000023, "MNE2","wNE2",1,1,0}
, {"~o3","~o3", 1000025, "MNE3","wNE3",1,1,0}
, {"~o4","~o4", 1000035, "MNE4","wNE4",1,1,0}
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
int nModelVars=39;
int nModelFunc=358;
static char*varNames_[397]={
 "alfSMZ","MZ","MW","Mtp","MbMb","McMc","Q","EE","Ml","Mq"
,"Au","Ad","tb","wt","MH3","mu","MG1","MG2","MG3","Ml1"
,"Ml2","Ml3","Mr1","Mr2","Mr3","Mq1","Mq2","Mq3","Mu1","Mu2"
,"Mu3","Md1","Md2","Md3","At","Ab","Al","Am","Maux","CW"
,"SW","S2W","C2W","GF","Mh","MHH","MHc","MNE1","MNE2","MNE3"
,"MNE4","MC1","MC2","MSG","MSne","MSnm","MSnl","MSeL","MSeR","MSmL"
,"MSmR","MSl1","MSl2","MSdL","MSdR","MSuL","MSuR","MSsL","MSsR","MScL"
,"MScR","MSb1","MSb2","MSt1","MSt2","QSUSY","Zn11","Zn12","Zn13","Zn14"
,"Zn21","Zn22","Zn23","Zn24","Zn31","Zn32","Zn33","Zn34","Zn41","Zn42"
,"Zn43","Zn44","Zu11","Zu12","Zu21","Zu22","Zv11","Zv12","Zv21","Zv22"
,"Zt11","Zt12","Zt21","Zt22","Zb11","Zb12","Zb21","Zb22","Zl11","Zl12"
,"Zl21","Zl22","alpha","tB","dMb","dMl","calcL","la1","la3","la6"
,"ca","sa","sb","cb","c2b","Zh11","Zh12","Zh21","Zh22","Td3"
,"Lqcd","Mb","Mt","Mc","PI","Mbp","Mcp","Mhh11","Mhh12","Mhh22"
,"la5","la4","la7","la2","dMd","Td2","fiuu","fidd","ficc","fiss"
,"Zuu11","Zuu12","Zuu21","Zuu22","Zdd11","Zdd12","Zdd21","Zdd22","Zcc11","Zcc12"
,"Zcc21","Zcc22","Zss11","Zss12","Zss21","Zss22","MtMM","MbMM","StMM11","StMM12"
,"StMM22","SbMM11","SbMM12","SbMM22","AT","AB","dX3","dX2","dX1","Tl3"
,"NMM11","NMM12","NMM13","NMM14","NMM22","NMM23","NMM24","NMM33","NMM34","NMM44"
,"MG1I","MG2I","aQCD","rhF_c","ihF_c","rAF_c","iAF_c","rhF_b","ihF_b","rAF_b"
,"iAF_b","rhF_t","ihF_t","rAF_t","iAF_t","rhF_l","ihF_l","rAF_l","iAF_l","rhS_eL"
,"ihS_eL","rhS_eR","ihS_eR","rhS_mL","ihS_mL","rhS_mR","ihS_mR","rhS_l1","ihS_l1","rhS_l2"
,"ihS_l2","rhS_uL","ihS_uL","rhS_uR","ihS_uR","rhS_cL","ihS_cL","rhS_cR","ihS_cR","rhS_t1"
,"ihS_t1","rhS_t2","ihS_t2","rhS_dL","ihS_dL","rhS_dR","ihS_dR","rhS_sL","ihS_sL","rhS_sR"
,"ihS_sR","rhS_b1","ihS_b1","rhS_b2","ihS_b2","rhS_Hc","ihS_Hc","rhV_W","ihV_W","rhF_c1"
,"ihF_c1","rAF_c1","iAF_c1","rhF_c2","ihF_c2","rAF_c2","iAF_c2","McR","MbR","MtR"
,"rhF1_c","ihF1_c","rAF1_c","iAF1_c","rhF1_b","ihF1_b","rAF1_b","iAF1_b","rhF1_t","ihF1_t"
,"rAF1_t","iAF1_t","tau2_uL","rhS1_uL","ihS1_uL","tau2_uR","rhS1_uR","ihS1_uR","tau2_cL","rhS1_cL"
,"ihS1_cL","tau2_cR","rhS1_cR","ihS1_cR","tau2_t1","rhS1_t1","ihS1_t1","tau2_t2","rhS1_t2","ihS1_t2"
,"tau2_dL","rhS1_dL","ihS1_dL","tau2_dR","rhS1_dR","ihS1_dR","tau2_sL","rhS1_sL","ihS1_sL","tau2_sR"
,"rhS1_sR","ihS1_sR","tau2_b1","rhS1_b1","ihS1_b1","tau2_b2","rhS1_b2","ihS1_b2","ahF_c","ahF_b"
,"ahF_t","ahF_l","aHF_c","aHF_b","aHF_t","aHF_l","ahS_eL","ahS_eR","ahS_mL","ahS_mR"
,"ahS_uL","ahS_uR","ahS_cL","ahS_cR","ahS_dL","ahS_dR","ahS_sL","ahS_sR","ahS_l1","ahS_l2"
,"ahS_t1","ahS_t2","ahS_b1","ahS_b2","aHS_eL","aHS_eR","aHS_mL","aHS_mR","aHS_uL","aHS_uR"
,"aHS_cL","aHS_cR","aHS_dL","aHS_dR","aHS_sL","aHS_sR","aHS_l1","aHS_l2","aHS_t1","aHS_t2"
,"aHS_b1","aHS_b2","ahV_W","B00000","B00001","B00002","ahS_Hc","aHV_W","B00003","B00004"
,"B00005","aHS_Hc","ahF_c1","ahF_c2","aHF_c1","aHF_c2","aAF_c","aAF_b","aAF_t","aAF_l"
,"aAF_c1","aAF_c2","Rqcd","lnTop","Ctop","Cq","Csq","alphaE0","Qu","Qd"
,"Fodd","ahF_b0","aHF_b0","aAF_b0","ahF_l0","aHF_l0","aAF_l0","LGGh","LAAh","LGGH"
,"LAAH","LGGH3","LAAH3","aSMhF_f","aSMhV_W","LGGSM","LAASM"};
char**varNames=varNames_;
static REAL varValues_[397]={
   1.184000E-01,  9.130000E+01,  8.020000E+01,  1.730700E+02,  4.230000E+00,  1.270000E+00,  1.000000E+02,  3.123000E-01,  1.777000E+00,  5.000000E-02
,  0.000000E+00,  0.000000E+00,  1.000000E+01,  1.442000E+00,  1.000000E+03,  3.500000E+02,  2.000000E+02,  4.000000E+02,  8.000000E+02,  5.000000E+02
,  5.000000E+02,  5.000000E+02,  2.000000E+02,  2.000000E+02,  2.000000E+02,  1.000000E+03,  1.000000E+03,  1.000000E+03,  3.000000E+02,  3.000000E+02
,  3.000000E+02,  3.000000E+02,  3.000000E+02,  3.000000E+02, -1.000000E+03,  0.000000E+00,  0.000000E+00,  0.000000E+00,  1.000000E+00};
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
   V[39]=V[2]/(V[1]);
   if(!isfinite(V[39]) || FError) return 39;
   V[40]=sqrt(1-pow(V[39],2));
   if(!isfinite(V[40]) || FError) return 40;
   V[41]=2*V[40]*V[39];
   if(!isfinite(V[41]) || FError) return 41;
   V[42]=pow(V[39],2)-pow(V[40],2);
   if(!isfinite(V[42]) || FError) return 42;
   V[43]=pow(V[7],2)/(pow(2*V[40]*V[2],2))/(M_SQRT2);
   if(!isfinite(V[43]) || FError) return 43;
   V[44]=slhaVal("MASS",V[1],1,25);
   if(!isfinite(V[44]) || FError) return 44;
   V[45]=slhaVal("MASS",V[1],1,35);
   if(!isfinite(V[45]) || FError) return 45;
   V[46]=slhaVal("MASS",V[1],1,37);
   if(!isfinite(V[46]) || FError) return 46;
   V[47]=slhaVal("MASS",V[1],1,1000022);
   if(!isfinite(V[47]) || FError) return 47;
   V[48]=slhaVal("MASS",V[1],1,1000023);
   if(!isfinite(V[48]) || FError) return 48;
   V[49]=slhaVal("MASS",V[1],1,1000025);
   if(!isfinite(V[49]) || FError) return 49;
   V[50]=slhaVal("MASS",V[1],1,1000035);
   if(!isfinite(V[50]) || FError) return 50;
   V[51]=slhaVal("MASS",V[1],1,1000024);
   if(!isfinite(V[51]) || FError) return 51;
   V[52]=slhaVal("MASS",V[1],1,1000037);
   if(!isfinite(V[52]) || FError) return 52;
   V[53]=slhaVal("MASS",V[1],1,1000021);
   if(!isfinite(V[53]) || FError) return 53;
   V[54]=slhaVal("MASS",V[1],1,1000012);
   if(!isfinite(V[54]) || FError) return 54;
   V[55]=slhaVal("MASS",V[1],1,1000014);
   if(!isfinite(V[55]) || FError) return 55;
   V[56]=slhaVal("MASS",V[1],1,1000016);
   if(!isfinite(V[56]) || FError) return 56;
   V[57]=slhaVal("MASS",V[1],1,1000011);
   if(!isfinite(V[57]) || FError) return 57;
   V[58]=slhaVal("MASS",V[1],1,2000011);
   if(!isfinite(V[58]) || FError) return 58;
   V[59]=slhaVal("MASS",V[1],1,1000013);
   if(!isfinite(V[59]) || FError) return 59;
   V[60]=slhaVal("MASS",V[1],1,2000013);
   if(!isfinite(V[60]) || FError) return 60;
   V[61]=slhaVal("MASS",V[1],1,1000015);
   if(!isfinite(V[61]) || FError) return 61;
   V[62]=slhaVal("MASS",V[1],1,2000015);
   if(!isfinite(V[62]) || FError) return 62;
   V[63]=slhaVal("MASS",V[1],1,1000001);
   if(!isfinite(V[63]) || FError) return 63;
   V[64]=slhaVal("MASS",V[1],1,2000001);
   if(!isfinite(V[64]) || FError) return 64;
   V[65]=slhaVal("MASS",V[1],1,1000002);
   if(!isfinite(V[65]) || FError) return 65;
   V[66]=slhaVal("MASS",V[1],1,2000002);
   if(!isfinite(V[66]) || FError) return 66;
   V[67]=slhaVal("MASS",V[1],1,1000003);
   if(!isfinite(V[67]) || FError) return 67;
   V[68]=slhaVal("MASS",V[1],1,2000003);
   if(!isfinite(V[68]) || FError) return 68;
   V[69]=slhaVal("MASS",V[1],1,1000004);
   if(!isfinite(V[69]) || FError) return 69;
   V[70]=slhaVal("MASS",V[1],1,2000004);
   if(!isfinite(V[70]) || FError) return 70;
   V[71]=slhaVal("MASS",V[1],1,1000005);
   if(!isfinite(V[71]) || FError) return 71;
   V[72]=slhaVal("MASS",V[1],1,2000005);
   if(!isfinite(V[72]) || FError) return 72;
   V[73]=slhaVal("MASS",V[1],1,1000006);
   if(!isfinite(V[73]) || FError) return 73;
   V[74]=slhaVal("MASS",V[1],1,2000006);
   if(!isfinite(V[74]) || FError) return 74;
   V[75]=sqrt(V[73]*V[74]);
   if(!isfinite(V[75]) || FError) return 75;
   V[76]=slhaVal("NMIX",V[75],2,1,1);
   if(!isfinite(V[76]) || FError) return 76;
   V[77]=slhaVal("NMIX",V[75],2,1,2);
   if(!isfinite(V[77]) || FError) return 77;
   V[78]=slhaVal("NMIX",V[75],2,1,3);
   if(!isfinite(V[78]) || FError) return 78;
   V[79]=slhaVal("NMIX",V[75],2,1,4);
   if(!isfinite(V[79]) || FError) return 79;
   V[80]=slhaVal("NMIX",V[75],2,2,1);
   if(!isfinite(V[80]) || FError) return 80;
   V[81]=slhaVal("NMIX",V[75],2,2,2);
   if(!isfinite(V[81]) || FError) return 81;
   V[82]=slhaVal("NMIX",V[75],2,2,3);
   if(!isfinite(V[82]) || FError) return 82;
   V[83]=slhaVal("NMIX",V[75],2,2,4);
   if(!isfinite(V[83]) || FError) return 83;
   V[84]=slhaVal("NMIX",V[75],2,3,1);
   if(!isfinite(V[84]) || FError) return 84;
   V[85]=slhaVal("NMIX",V[75],2,3,2);
   if(!isfinite(V[85]) || FError) return 85;
   V[86]=slhaVal("NMIX",V[75],2,3,3);
   if(!isfinite(V[86]) || FError) return 86;
   V[87]=slhaVal("NMIX",V[75],2,3,4);
   if(!isfinite(V[87]) || FError) return 87;
   V[88]=slhaVal("NMIX",V[75],2,4,1);
   if(!isfinite(V[88]) || FError) return 88;
   V[89]=slhaVal("NMIX",V[75],2,4,2);
   if(!isfinite(V[89]) || FError) return 89;
   V[90]=slhaVal("NMIX",V[75],2,4,3);
   if(!isfinite(V[90]) || FError) return 90;
   V[91]=slhaVal("NMIX",V[75],2,4,4);
   if(!isfinite(V[91]) || FError) return 91;
   V[92]=slhaVal("UMIX",V[75],2,1,1);
   if(!isfinite(V[92]) || FError) return 92;
   V[93]=slhaVal("UMIX",V[75],2,1,2);
   if(!isfinite(V[93]) || FError) return 93;
   V[94]=slhaVal("UMIX",V[75],2,2,1);
   if(!isfinite(V[94]) || FError) return 94;
   V[95]=slhaVal("UMIX",V[75],2,2,2);
   if(!isfinite(V[95]) || FError) return 95;
   V[96]=slhaVal("VMIX",V[75],2,1,1);
   if(!isfinite(V[96]) || FError) return 96;
   V[97]=slhaVal("VMIX",V[75],2,1,2);
   if(!isfinite(V[97]) || FError) return 97;
   V[98]=slhaVal("VMIX",V[75],2,2,1);
   if(!isfinite(V[98]) || FError) return 98;
   V[99]=slhaVal("VMIX",V[75],2,2,2);
   if(!isfinite(V[99]) || FError) return 99;
   V[100]=slhaVal("STOPMIX",V[75],2,1,1);
   if(!isfinite(V[100]) || FError) return 100;
   V[101]=slhaVal("STOPMIX",V[75],2,1,2);
   if(!isfinite(V[101]) || FError) return 101;
   V[102]=slhaVal("STOPMIX",V[75],2,2,1);
   if(!isfinite(V[102]) || FError) return 102;
   V[103]=slhaVal("STOPMIX",V[75],2,2,2);
   if(!isfinite(V[103]) || FError) return 103;
   V[104]=slhaVal("SBOTMIX",V[75],2,1,1);
   if(!isfinite(V[104]) || FError) return 104;
   V[105]=slhaVal("SBOTMIX",V[75],2,1,2);
   if(!isfinite(V[105]) || FError) return 105;
   V[106]=slhaVal("SBOTMIX",V[75],2,2,1);
   if(!isfinite(V[106]) || FError) return 106;
   V[107]=slhaVal("SBOTMIX",V[75],2,2,2);
   if(!isfinite(V[107]) || FError) return 107;
   V[108]=slhaVal("STAUMIX",V[75],2,1,1);
   if(!isfinite(V[108]) || FError) return 108;
   V[109]=slhaVal("STAUMIX",V[75],2,1,2);
   if(!isfinite(V[109]) || FError) return 109;
   V[110]=slhaVal("STAUMIX",V[75],2,2,1);
   if(!isfinite(V[110]) || FError) return 110;
   V[111]=slhaVal("STAUMIX",V[75],2,2,2);
   if(!isfinite(V[111]) || FError) return 111;
   V[112]=slhaVal("ALPHA",V[75],0);
   if(!isfinite(V[112]) || FError) return 112;
   V[113]=slhaVal("HMIX",0,1,2);
   if(!isfinite(V[113]) || FError) return 113;
   V[114]=deltaMb();
   if(!isfinite(V[114]) || FError) return 114;
   V[115]=deltaMl();
   if(!isfinite(V[115]) || FError) return 115;
   V[116]=calcLambdas();
   if(!isfinite(V[116]) || FError) return 116;
   V[117]=Lambda1();
   if(!isfinite(V[117]) || FError) return 117;
   V[118]=Lambda3();
   if(!isfinite(V[118]) || FError) return 118;
   V[119]=Lambda6();
   if(!isfinite(V[119]) || FError) return 119;
   V[120]=cos(V[112]);
   if(!isfinite(V[120]) || FError) return 120;
   V[121]=sin(V[112]);
   if(!isfinite(V[121]) || FError) return 121;
   V[122]=V[113]/(sqrt(1+pow(V[113],2)));
   if(!isfinite(V[122]) || FError) return 122;
   V[123]=sqrt(1-pow(V[122],2));
   if(!isfinite(V[123]) || FError) return 123;
   V[124]=pow(V[123],2)-pow(V[122],2);
   if(!isfinite(V[124]) || FError) return 124;
   V[125]=V[120];

   V[126]=V[121];

   V[127]=-V[121];
   if(!isfinite(V[127]) || FError) return 127;
   V[128]=V[120];

   V[129]=V[114]/(1+V[114]);
   if(!isfinite(V[129]) || FError) return 129;
   V[130]=initQCD5(V[0],V[5],V[4],V[3]);
   if(!isfinite(V[130]) || FError) return 130;
 FirstQ:
 cErr=1;
   V[131]=MbEff(V[6]);
   if(!isfinite(V[131]) || FError) return 131;
   V[132]=MtEff(V[6]);
   if(!isfinite(V[132]) || FError) return 132;
   V[133]=McEff(V[6]);
   if(!isfinite(V[133]) || FError) return 133;
   V[134]=acos(-1);
   if(!isfinite(V[134]) || FError) return 134;
   V[135]=V[4]*(1+4/(double)((3))*alphaQCD(V[4])/(V[134]));
   if(!isfinite(V[135]) || FError) return 135;
   V[136]=V[5]*(1+4/(double)((3))*alphaQCD(V[5])/(V[134]));
   if(!isfinite(V[136]) || FError) return 136;
   V[137]=pow(V[44],2)*pow(V[126],2)+pow(V[45],2)*pow(V[125],2);
   if(!isfinite(V[137]) || FError) return 137;
   V[138]=pow(V[45],2)*V[125]*V[126]+pow(V[44],2)*V[128]*V[127];
   if(!isfinite(V[138]) || FError) return 138;
   V[139]=pow(V[44],2)*pow(V[128],2)+pow(V[45],2)*pow(V[127],2);
   if(!isfinite(V[139]) || FError) return 139;
   V[140]=2*V[119]*V[123]/(V[122])-V[117]*pow(V[123],2)/(pow(V[122],2))-pow(V[7],2)*(pow(V[14]*V[122],2)-V[137])/(pow(2*V[2]*V[40]*V[122],2));
   if(!isfinite(V[140]) || FError) return 140;
   V[141]=V[140]-1/(double)((2))*pow(V[7],2)*(pow(V[46],2)-pow(V[14],2))/(pow(V[2]*V[40],2));
   if(!isfinite(V[141]) || FError) return 141;
   V[142]=V[123]/(V[122])*(V[118]+V[141]-V[123]/(V[122])*V[119]-1/(double)((4))*pow(V[7],2)*(V[138]/(V[123])/(V[122])+pow(V[14],2))/(pow(V[2]*V[40],2)));
   if(!isfinite(V[142]) || FError) return 142;
   V[143]=2*V[142]*V[123]/(V[122])-V[140]*pow(V[123],2)/(pow(V[122],2))-pow(V[7],2)*(pow(V[14]*V[123],2)-V[139])/(pow(2*V[2]*V[40]*V[122],2));
   if(!isfinite(V[143]) || FError) return 143;
   V[144]=deltaMd();
   if(!isfinite(V[144]) || FError) return 144;
   V[145]=V[144]/(1+V[144]);
   if(!isfinite(V[145]) || FError) return 145;
   V[146]=atan(-2*V[9]*(V[10]-V[15]/(V[113]))/(pow(V[66],2)-pow(V[65],2)))/(2);
   if(!isfinite(V[146]) || FError) return 146;
   V[147]=atan(-2*V[9]*(1-V[145])*(V[11]-V[15]*V[113])/(pow(V[64],2)-pow(V[63],2)))/(2);
   if(!isfinite(V[147]) || FError) return 147;
   V[148]=atan(-2*V[133]*(V[10]-V[15]/(V[113]))/(pow(V[70],2)-pow(V[69],2)))/(2);
   if(!isfinite(V[148]) || FError) return 148;
   V[149]=atan(-2*V[9]*(1-V[145])*(V[11]-V[15]*V[113])/(pow(V[68],2)-pow(V[67],2)))/(2);
   if(!isfinite(V[149]) || FError) return 149;
   V[150]=cos(V[146]);
   if(!isfinite(V[150]) || FError) return 150;
   V[151]=sin(V[146]);
   if(!isfinite(V[151]) || FError) return 151;
   V[152]=-sin(V[146]);
   if(!isfinite(V[152]) || FError) return 152;
   V[153]=cos(V[146]);
   if(!isfinite(V[153]) || FError) return 153;
   V[154]=cos(V[147]);
   if(!isfinite(V[154]) || FError) return 154;
   V[155]=sin(V[147]);
   if(!isfinite(V[155]) || FError) return 155;
   V[156]=-sin(V[147]);
   if(!isfinite(V[156]) || FError) return 156;
   V[157]=cos(V[147]);
   if(!isfinite(V[157]) || FError) return 157;
   V[158]=cos(V[148]);
   if(!isfinite(V[158]) || FError) return 158;
   V[159]=sin(V[148]);
   if(!isfinite(V[159]) || FError) return 159;
   V[160]=-sin(V[148]);
   if(!isfinite(V[160]) || FError) return 160;
   V[161]=cos(V[148]);
   if(!isfinite(V[161]) || FError) return 161;
   V[162]=cos(V[149]);
   if(!isfinite(V[162]) || FError) return 162;
   V[163]=sin(V[149]);
   if(!isfinite(V[163]) || FError) return 163;
   V[164]=-sin(V[149]);
   if(!isfinite(V[164]) || FError) return 164;
   V[165]=cos(V[149]);
   if(!isfinite(V[165]) || FError) return 165;
   V[166]=MtRun(sqrt(V[73]*V[74]));
   if(!isfinite(V[166]) || FError) return 166;
   V[167]=MbRun(sqrt(V[71]*V[72]));
   if(!isfinite(V[167]) || FError) return 167;
   V[168]=V[100]*pow(V[73],2)*V[100]+V[102]*pow(V[74],2)*V[102];
   if(!isfinite(V[168]) || FError) return 168;
   V[169]=V[100]*pow(V[73],2)*V[101]+V[102]*pow(V[74],2)*V[103];
   if(!isfinite(V[169]) || FError) return 169;
   V[170]=V[101]*pow(V[73],2)*V[101]+V[103]*pow(V[74],2)*V[103];
   if(!isfinite(V[170]) || FError) return 170;
   V[171]=V[104]*pow(V[71],2)*V[104]+V[106]*pow(V[72],2)*V[106];
   if(!isfinite(V[171]) || FError) return 171;
   V[172]=V[104]*pow(V[71],2)*V[105]+V[106]*pow(V[72],2)*V[107];
   if(!isfinite(V[172]) || FError) return 172;
   V[173]=V[105]*pow(V[71],2)*V[105]+V[107]*pow(V[72],2)*V[107];
   if(!isfinite(V[173]) || FError) return 173;
   V[174]=(V[169]+V[15]*V[123]/(V[122])*V[166])/(V[166]);
   if(!isfinite(V[174]) || FError) return 174;
   V[175]=(V[172]+V[15]*V[122]/(V[123])*V[167])/(V[167]);
   if(!isfinite(V[175]) || FError) return 175;
   V[176]=-(V[168]-V[171]+pow(V[167],2)-pow(V[166],2)-V[124]*pow(V[2],2))/(pow(V[2],2)*pow(V[122],2));
   if(!isfinite(V[176]) || FError) return 176;
   V[177]=-(pow(V[69],2)-pow(V[67],2)-V[124]*pow(V[2],2))/(pow(V[2],2)*pow(V[122],2));
   if(!isfinite(V[177]) || FError) return 177;
   V[178]=-(pow(V[65],2)-pow(V[63],2)-V[124]*pow(V[2],2))/(pow(V[2],2)*pow(V[122],2));
   if(!isfinite(V[178]) || FError) return 178;
   V[179]=V[115]/(1+V[115]);
   if(!isfinite(V[179]) || FError) return 179;
   V[180]=V[76]*V[47]*V[76]+V[80]*V[48]*V[80]+V[84]*V[49]*V[84]+V[88]*V[50]*V[88];
   if(!isfinite(V[180]) || FError) return 180;
   V[181]=V[76]*V[47]*V[77]+V[80]*V[48]*V[81]+V[84]*V[49]*V[85]+V[88]*V[50]*V[89];
   if(!isfinite(V[181]) || FError) return 181;
   V[182]=V[76]*V[47]*V[78]+V[80]*V[48]*V[82]+V[84]*V[49]*V[86]+V[88]*V[50]*V[90];
   if(!isfinite(V[182]) || FError) return 182;
   V[183]=V[76]*V[47]*V[79]+V[80]*V[48]*V[83]+V[84]*V[49]*V[87]+V[88]*V[50]*V[91];
   if(!isfinite(V[183]) || FError) return 183;
   V[184]=V[77]*V[47]*V[77]+V[81]*V[48]*V[81]+V[85]*V[49]*V[85]+V[89]*V[50]*V[89];
   if(!isfinite(V[184]) || FError) return 184;
   V[185]=V[77]*V[47]*V[78]+V[81]*V[48]*V[82]+V[85]*V[49]*V[86]+V[89]*V[50]*V[90];
   if(!isfinite(V[185]) || FError) return 185;
   V[186]=V[77]*V[47]*V[79]+V[81]*V[48]*V[83]+V[85]*V[49]*V[87]+V[89]*V[50]*V[91];
   if(!isfinite(V[186]) || FError) return 186;
   V[187]=V[78]*V[47]*V[78]+V[82]*V[48]*V[82]+V[86]*V[49]*V[86]+V[90]*V[50]*V[90];
   if(!isfinite(V[187]) || FError) return 187;
   V[188]=V[78]*V[47]*V[79]+V[82]*V[48]*V[83]+V[86]*V[49]*V[87]+V[90]*V[50]*V[91];
   if(!isfinite(V[188]) || FError) return 188;
   V[189]=V[79]*V[47]*V[79]+V[83]*V[48]*V[83]+V[87]*V[49]*V[87]+V[91]*V[50]*V[91];
   if(!isfinite(V[189]) || FError) return 189;
   V[190]=V[180];

   V[191]=V[184];

   V[192]=alphaQCD(V[6])/(V[134]);
   if(!isfinite(V[192]) || FError) return 192;
   V[193]=creal(HggF(pow(V[6]/(2)/(V[136]),2)));
   if(!isfinite(V[193]) || FError) return 193;
   V[194]=cimag(HggF(pow(V[6]/(2)/(V[136]),2)));
   if(!isfinite(V[194]) || FError) return 194;
   V[195]=creal(HggA(pow(V[6]/(2)/(V[136]),2)));
   if(!isfinite(V[195]) || FError) return 195;
   V[196]=cimag(HggA(pow(V[6]/(2)/(V[136]),2)));
   if(!isfinite(V[196]) || FError) return 196;
   V[197]=creal(HggF(pow(V[6]/(2)/(V[135]),2)));
   if(!isfinite(V[197]) || FError) return 197;
   V[198]=cimag(HggF(pow(V[6]/(2)/(V[135]),2)));
   if(!isfinite(V[198]) || FError) return 198;
   V[199]=creal(HggA(pow(V[6]/(2)/(V[135]),2)));
   if(!isfinite(V[199]) || FError) return 199;
   V[200]=cimag(HggA(pow(V[6]/(2)/(V[135]),2)));
   if(!isfinite(V[200]) || FError) return 200;
   V[201]=creal(HggF(pow(V[6]/(2)/(V[3]),2)));
   if(!isfinite(V[201]) || FError) return 201;
   V[202]=cimag(HggF(pow(V[6]/(2)/(V[3]),2)));
   if(!isfinite(V[202]) || FError) return 202;
   V[203]=creal(HggA(pow(V[6]/(2)/(V[3]),2)));
   if(!isfinite(V[203]) || FError) return 203;
   V[204]=cimag(HggA(pow(V[6]/(2)/(V[3]),2)));
   if(!isfinite(V[204]) || FError) return 204;
   V[205]=creal(HggF(pow(V[6]/(2)/(V[8]),2)));
   if(!isfinite(V[205]) || FError) return 205;
   V[206]=cimag(HggF(pow(V[6]/(2)/(V[8]),2)));
   if(!isfinite(V[206]) || FError) return 206;
   V[207]=creal(HggA(pow(V[6]/(2)/(V[8]),2)));
   if(!isfinite(V[207]) || FError) return 207;
   V[208]=cimag(HggA(pow(V[6]/(2)/(V[8]),2)));
   if(!isfinite(V[208]) || FError) return 208;
   V[209]=creal(HggS(pow(V[6]/(2)/(V[57]),2)));
   if(!isfinite(V[209]) || FError) return 209;
   V[210]=cimag(HggS(pow(V[6]/(2)/(V[57]),2)));
   if(!isfinite(V[210]) || FError) return 210;
   V[211]=creal(HggS(pow(V[6]/(2)/(V[58]),2)));
   if(!isfinite(V[211]) || FError) return 211;
   V[212]=cimag(HggS(pow(V[6]/(2)/(V[58]),2)));
   if(!isfinite(V[212]) || FError) return 212;
   V[213]=creal(HggS(pow(V[6]/(2)/(V[59]),2)));
   if(!isfinite(V[213]) || FError) return 213;
   V[214]=cimag(HggS(pow(V[6]/(2)/(V[59]),2)));
   if(!isfinite(V[214]) || FError) return 214;
   V[215]=creal(HggS(pow(V[6]/(2)/(V[60]),2)));
   if(!isfinite(V[215]) || FError) return 215;
   V[216]=cimag(HggS(pow(V[6]/(2)/(V[60]),2)));
   if(!isfinite(V[216]) || FError) return 216;
   V[217]=creal(HggS(pow(V[6]/(2)/(V[61]),2)));
   if(!isfinite(V[217]) || FError) return 217;
   V[218]=cimag(HggS(pow(V[6]/(2)/(V[61]),2)));
   if(!isfinite(V[218]) || FError) return 218;
   V[219]=creal(HggS(pow(V[6]/(2)/(V[62]),2)));
   if(!isfinite(V[219]) || FError) return 219;
   V[220]=cimag(HggS(pow(V[6]/(2)/(V[62]),2)));
   if(!isfinite(V[220]) || FError) return 220;
   V[221]=creal(HggS(pow(V[6]/(2)/(V[65]),2)));
   if(!isfinite(V[221]) || FError) return 221;
   V[222]=cimag(HggS(pow(V[6]/(2)/(V[65]),2)));
   if(!isfinite(V[222]) || FError) return 222;
   V[223]=creal(HggS(pow(V[6]/(2)/(V[66]),2)));
   if(!isfinite(V[223]) || FError) return 223;
   V[224]=cimag(HggS(pow(V[6]/(2)/(V[66]),2)));
   if(!isfinite(V[224]) || FError) return 224;
   V[225]=creal(HggS(pow(V[6]/(2)/(V[69]),2)));
   if(!isfinite(V[225]) || FError) return 225;
   V[226]=cimag(HggS(pow(V[6]/(2)/(V[69]),2)));
   if(!isfinite(V[226]) || FError) return 226;
   V[227]=creal(HggS(pow(V[6]/(2)/(V[70]),2)));
   if(!isfinite(V[227]) || FError) return 227;
   V[228]=cimag(HggS(pow(V[6]/(2)/(V[70]),2)));
   if(!isfinite(V[228]) || FError) return 228;
   V[229]=creal(HggS(pow(V[6]/(2)/(V[73]),2)));
   if(!isfinite(V[229]) || FError) return 229;
   V[230]=cimag(HggS(pow(V[6]/(2)/(V[73]),2)));
   if(!isfinite(V[230]) || FError) return 230;
   V[231]=creal(HggS(pow(V[6]/(2)/(V[74]),2)));
   if(!isfinite(V[231]) || FError) return 231;
   V[232]=cimag(HggS(pow(V[6]/(2)/(V[74]),2)));
   if(!isfinite(V[232]) || FError) return 232;
   V[233]=creal(HggS(pow(V[6]/(2)/(V[63]),2)));
   if(!isfinite(V[233]) || FError) return 233;
   V[234]=cimag(HggS(pow(V[6]/(2)/(V[63]),2)));
   if(!isfinite(V[234]) || FError) return 234;
   V[235]=creal(HggS(pow(V[6]/(2)/(V[64]),2)));
   if(!isfinite(V[235]) || FError) return 235;
   V[236]=cimag(HggS(pow(V[6]/(2)/(V[64]),2)));
   if(!isfinite(V[236]) || FError) return 236;
   V[237]=creal(HggS(pow(V[6]/(2)/(V[67]),2)));
   if(!isfinite(V[237]) || FError) return 237;
   V[238]=cimag(HggS(pow(V[6]/(2)/(V[67]),2)));
   if(!isfinite(V[238]) || FError) return 238;
   V[239]=creal(HggS(pow(V[6]/(2)/(V[68]),2)));
   if(!isfinite(V[239]) || FError) return 239;
   V[240]=cimag(HggS(pow(V[6]/(2)/(V[68]),2)));
   if(!isfinite(V[240]) || FError) return 240;
   V[241]=creal(HggS(pow(V[6]/(2)/(V[71]),2)));
   if(!isfinite(V[241]) || FError) return 241;
   V[242]=cimag(HggS(pow(V[6]/(2)/(V[71]),2)));
   if(!isfinite(V[242]) || FError) return 242;
   V[243]=creal(HggS(pow(V[6]/(2)/(V[72]),2)));
   if(!isfinite(V[243]) || FError) return 243;
   V[244]=cimag(HggS(pow(V[6]/(2)/(V[72]),2)));
   if(!isfinite(V[244]) || FError) return 244;
   V[245]=creal(HggS(pow(V[6]/(2)/(V[46]),2)));
   if(!isfinite(V[245]) || FError) return 245;
   V[246]=cimag(HggS(pow(V[6]/(2)/(V[46]),2)));
   if(!isfinite(V[246]) || FError) return 246;
   V[247]=creal(HggV(pow(V[6]/(2)/(V[2]),2)));
   if(!isfinite(V[247]) || FError) return 247;
   V[248]=cimag(HggV(pow(V[6]/(2)/(V[2]),2)));
   if(!isfinite(V[248]) || FError) return 248;
   V[249]=creal(HggF(pow(V[6]/(2)/(V[51]),2)));
   if(!isfinite(V[249]) || FError) return 249;
   V[250]=cimag(HggF(pow(V[6]/(2)/(V[51]),2)));
   if(!isfinite(V[250]) || FError) return 250;
   V[251]=creal(HggA(pow(V[6]/(2)/(V[51]),2)));
   if(!isfinite(V[251]) || FError) return 251;
   V[252]=cimag(HggA(pow(V[6]/(2)/(V[51]),2)));
   if(!isfinite(V[252]) || FError) return 252;
   V[253]=creal(HggF(pow(V[6]/(2)/(V[52]),2)));
   if(!isfinite(V[253]) || FError) return 253;
   V[254]=cimag(HggF(pow(V[6]/(2)/(V[52]),2)));
   if(!isfinite(V[254]) || FError) return 254;
   V[255]=creal(HggA(pow(V[6]/(2)/(V[52]),2)));
   if(!isfinite(V[255]) || FError) return 255;
   V[256]=cimag(HggA(pow(V[6]/(2)/(V[52]),2)));
   if(!isfinite(V[256]) || FError) return 256;
   V[257]=V[136]*McRun(V[6]/(2))/(McRun(V[136]));
   if(!isfinite(V[257]) || FError) return 257;
   V[258]=V[135]*MbRun(V[6]/(2))/(MbRun(V[135]));
   if(!isfinite(V[258]) || FError) return 258;
   V[259]=V[3]*MtRun(V[6]/(2))/(MtRun(V[3]));
   if(!isfinite(V[259]) || FError) return 259;
   V[260]=creal(HggF(pow(V[6]/(2)/(V[257]),2))*(1+V[192]*Hgam1F(pow(V[6]/(2)/(V[257]),2))));
   if(!isfinite(V[260]) || FError) return 260;
   V[261]=cimag(HggF(pow(V[6]/(2)/(V[257]),2))*(1+V[192]*Hgam1F(pow(V[6]/(2)/(V[257]),2))));
   if(!isfinite(V[261]) || FError) return 261;
   V[262]=creal(HggA(pow(V[6]/(2)/(V[257]),2))*(1+V[192]*Hgam1A(pow(V[6]/(2)/(V[257]),2))));
   if(!isfinite(V[262]) || FError) return 262;
   V[263]=cimag(HggA(pow(V[6]/(2)/(V[257]),2))*(1+V[192]*Hgam1A(pow(V[6]/(2)/(V[257]),2))));
   if(!isfinite(V[263]) || FError) return 263;
   V[264]=creal(HggF(pow(V[6]/(2)/(V[258]),2))*(1+V[192]*Hgam1F(pow(V[6]/(2)/(V[258]),2))));
   if(!isfinite(V[264]) || FError) return 264;
   V[265]=cimag(HggF(pow(V[6]/(2)/(V[258]),2))*(1+V[192]*Hgam1F(pow(V[6]/(2)/(V[258]),2))));
   if(!isfinite(V[265]) || FError) return 265;
   V[266]=creal(HggA(pow(V[6]/(2)/(V[258]),2))*(1+V[192]*Hgam1A(pow(V[6]/(2)/(V[258]),2))));
   if(!isfinite(V[266]) || FError) return 266;
   V[267]=cimag(HggA(pow(V[6]/(2)/(V[258]),2))*(1+V[192]*Hgam1A(pow(V[6]/(2)/(V[258]),2))));
   if(!isfinite(V[267]) || FError) return 267;
   V[268]=creal(HggF(pow(V[6]/(2)/(V[259]),2))*(1+V[192]*Hgam1F(pow(V[6]/(2)/(V[259]),2))));
   if(!isfinite(V[268]) || FError) return 268;
   V[269]=cimag(HggF(pow(V[6]/(2)/(V[259]),2))*(1+V[192]*Hgam1F(pow(V[6]/(2)/(V[259]),2))));
   if(!isfinite(V[269]) || FError) return 269;
   V[270]=creal(HggA(pow(V[6]/(2)/(V[259]),2))*(1+V[192]*Hgam1A(pow(V[6]/(2)/(V[259]),2))));
   if(!isfinite(V[270]) || FError) return 270;
   V[271]=cimag(HggA(pow(V[6]/(2)/(V[259]),2))*(1+V[192]*Hgam1A(pow(V[6]/(2)/(V[259]),2))));
   if(!isfinite(V[271]) || FError) return 271;
   V[272]=pow(V[6]/(2*V[65]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[65])),6/(double)((23)))),2);
   if(!isfinite(V[272]) || FError) return 272;
   V[273]=creal(HggS(V[272])*(1+V[192]*Hgam1S(V[272])));
   if(!isfinite(V[273]) || FError) return 273;
   V[274]=cimag(HggS(V[272])*(1+V[192]*Hgam1S(V[272])));
   if(!isfinite(V[274]) || FError) return 274;
   V[275]=pow(V[6]/(2*V[66]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[66])),6/(double)((23)))),2);
   if(!isfinite(V[275]) || FError) return 275;
   V[276]=creal(HggS(V[275])*(1+V[192]*Hgam1S(V[275])));
   if(!isfinite(V[276]) || FError) return 276;
   V[277]=cimag(HggS(V[275])*(1+V[192]*Hgam1S(V[275])));
   if(!isfinite(V[277]) || FError) return 277;
   V[278]=pow(V[6]/(2*V[69]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[69])),6/(double)((23)))),2);
   if(!isfinite(V[278]) || FError) return 278;
   V[279]=creal(HggS(V[278])*(1+V[192]*Hgam1S(V[278])));
   if(!isfinite(V[279]) || FError) return 279;
   V[280]=cimag(HggS(V[278])*(1+V[192]*Hgam1S(V[278])));
   if(!isfinite(V[280]) || FError) return 280;
   V[281]=pow(V[6]/(2*V[70]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[70])),6/(double)((23)))),2);
   if(!isfinite(V[281]) || FError) return 281;
   V[282]=creal(HggS(V[281])*(1+V[192]*Hgam1S(V[281])));
   if(!isfinite(V[282]) || FError) return 282;
   V[283]=cimag(HggS(V[281])*(1+V[192]*Hgam1S(V[281])));
   if(!isfinite(V[283]) || FError) return 283;
   V[284]=pow(V[6]/(2*V[73]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[73])),6/(double)((23)))),2);
   if(!isfinite(V[284]) || FError) return 284;
   V[285]=creal(HggS(V[284])*(1+V[192]*Hgam1S(V[284])));
   if(!isfinite(V[285]) || FError) return 285;
   V[286]=cimag(HggS(V[284])*(1+V[192]*Hgam1S(V[284])));
   if(!isfinite(V[286]) || FError) return 286;
   V[287]=pow(V[6]/(2*V[74]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[74])),6/(double)((23)))),2);
   if(!isfinite(V[287]) || FError) return 287;
   V[288]=creal(HggS(V[287])*(1+V[192]*Hgam1S(V[287])));
   if(!isfinite(V[288]) || FError) return 288;
   V[289]=cimag(HggS(V[287])*(1+V[192]*Hgam1S(V[287])));
   if(!isfinite(V[289]) || FError) return 289;
   V[290]=pow(V[6]/(2*V[63]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[63])),6/(double)((23)))),2);
   if(!isfinite(V[290]) || FError) return 290;
   V[291]=creal(HggS(V[290])*(1+V[192]*Hgam1S(V[290])));
   if(!isfinite(V[291]) || FError) return 291;
   V[292]=cimag(HggS(V[290])*(1+V[192]*Hgam1S(V[290])));
   if(!isfinite(V[292]) || FError) return 292;
   V[293]=pow(V[6]/(2*V[64]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[64])),6/(double)((23)))),2);
   if(!isfinite(V[293]) || FError) return 293;
   V[294]=creal(HggS(V[293])*(1+V[192]*Hgam1S(V[293])));
   if(!isfinite(V[294]) || FError) return 294;
   V[295]=cimag(HggS(V[293])*(1+V[192]*Hgam1S(V[293])));
   if(!isfinite(V[295]) || FError) return 295;
   V[296]=pow(V[6]/(2*V[67]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[67])),6/(double)((23)))),2);
   if(!isfinite(V[296]) || FError) return 296;
   V[297]=creal(HggS(V[296])*(1+V[192]*Hgam1S(V[296])));
   if(!isfinite(V[297]) || FError) return 297;
   V[298]=cimag(HggS(V[296])*(1+V[192]*Hgam1S(V[296])));
   if(!isfinite(V[298]) || FError) return 298;
   V[299]=pow(V[6]/(2*V[68]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[68])),6/(double)((23)))),2);
   if(!isfinite(V[299]) || FError) return 299;
   V[300]=creal(HggS(V[299])*(1+V[192]*Hgam1S(V[299])));
   if(!isfinite(V[300]) || FError) return 300;
   V[301]=cimag(HggS(V[299])*(1+V[192]*Hgam1S(V[299])));
   if(!isfinite(V[301]) || FError) return 301;
   V[302]=pow(V[6]/(2*V[71]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[71])),6/(double)((23)))),2);
   if(!isfinite(V[302]) || FError) return 302;
   V[303]=creal(HggS(V[302])*(1+V[192]*Hgam1S(V[302])));
   if(!isfinite(V[303]) || FError) return 303;
   V[304]=cimag(HggS(V[302])*(1+V[192]*Hgam1S(V[302])));
   if(!isfinite(V[304]) || FError) return 304;
   V[305]=pow(V[6]/(2*V[72]*pow(alphaQCD(V[6]/(2))/(alphaQCD(V[72])),6/(double)((23)))),2);
   if(!isfinite(V[305]) || FError) return 305;
   V[306]=creal(HggS(V[305])*(1+V[192]*Hgam1S(V[305])));
   if(!isfinite(V[306]) || FError) return 306;
   V[307]=cimag(HggS(V[305])*(1+V[192]*Hgam1S(V[305])));
   if(!isfinite(V[307]) || FError) return 307;
   V[308]=-V[7]/(V[2])*V[133]/(V[40])*V[128]/(V[122])/(2)/(V[133]);
   if(!isfinite(V[308]) || FError) return 308;
   V[309]=-V[7]/(V[2])*V[131]/(V[40])/(V[123])/(V[122])*(V[122]*V[127]-V[122]*V[129]*V[127]+V[129]*V[128]*V[123])/(2)/(V[131]);
   if(!isfinite(V[309]) || FError) return 309;
   V[310]=-V[7]/(V[2])*V[132]/(V[40])*V[128]/(V[122])/(2)/(V[132]);
   if(!isfinite(V[310]) || FError) return 310;
   V[311]=-V[7]/(V[2])*V[8]/(V[40])/(V[123])/(V[122])*(V[122]*V[127]-V[122]*V[179]*V[127]+V[179]*V[128]*V[123])/(2)/(V[8]);
   if(!isfinite(V[311]) || FError) return 311;
   V[312]=-V[7]/(V[2])*V[133]/(V[40])*V[126]/(V[122])/(2)/(V[133]);
   if(!isfinite(V[312]) || FError) return 312;
   V[313]=-V[7]/(V[2])*V[131]/(V[40])/(V[123])/(V[122])*(V[122]*V[125]-V[122]*V[129]*V[125]+V[129]*V[126]*V[123])/(2)/(V[131]);
   if(!isfinite(V[313]) || FError) return 313;
   V[314]=-V[7]/(V[2])*V[132]/(V[40])*V[126]/(V[122])/(2)/(V[132]);
   if(!isfinite(V[314]) || FError) return 314;
   V[315]=-V[7]/(V[2])*V[8]/(V[40])/(V[123])/(V[122])*(V[122]*V[125]-V[122]*V[179]*V[125]+V[179]*V[126]*V[123])/(2)/(V[8]);
   if(!isfinite(V[315]) || FError) return 315;
   V[316]=1/(pow(V[39],2))*V[7]*V[2]/(V[40])*(V[42]*V[127]*V[123]-V[42]*V[122]*V[128])/(2)/(pow(V[57],2))/(2);
   if(!isfinite(V[316]) || FError) return 316;
   V[317]=1/(pow(V[39],2))*V[7]*V[2]*V[40]*(V[127]*V[123]-V[122]*V[128])/(pow(V[58],2))/(2);
   if(!isfinite(V[317]) || FError) return 317;
   V[318]=1/(pow(V[39],2))*V[7]*V[2]/(V[40])*(V[42]*V[127]*V[123]-V[42]*V[122]*V[128])/(2)/(pow(V[59],2))/(2);
   if(!isfinite(V[318]) || FError) return 318;
   V[319]=1/(pow(V[39],2))*V[7]*V[2]*V[40]*(V[127]*V[123]-V[122]*V[128])/(pow(V[60],2))/(2);
   if(!isfinite(V[319]) || FError) return 319;
   V[320]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[150],2)*V[178]+4*pow(V[40],2)*V[122]*pow(V[2],2)*V[127]*pow(V[150],2)*V[123]-3*V[122]*pow(V[2],2)*V[127]*pow(V[150],2)*V[123]+3*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[150],2)-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[150],2)-4*pow(V[40],2)*V[122]*pow(V[2],2)*V[127]*pow(V[151],2)*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[151],2)+6*pow(V[39],2)*V[9]*V[127]*V[150]*V[151]*V[15]-6*pow(V[39],2)*V[10]*V[9]*V[128]*V[150]*V[151]-6*pow(V[39],2)*pow(V[9],2)*V[128])/(6)/(pow(V[65],2))/(2);
   if(!isfinite(V[320]) || FError) return 320;
   V[321]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[152],2)*V[178]+4*pow(V[40],2)*V[122]*pow(V[2],2)*V[127]*pow(V[152],2)*V[123]-3*V[122]*pow(V[2],2)*V[127]*pow(V[152],2)*V[123]+3*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[152],2)-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[152],2)-4*pow(V[40],2)*V[122]*pow(V[2],2)*V[127]*pow(V[153],2)*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[153],2)+6*pow(V[39],2)*V[9]*V[127]*V[152]*V[153]*V[15]-6*pow(V[39],2)*V[10]*V[9]*V[128]*V[152]*V[153]-6*pow(V[39],2)*pow(V[9],2)*V[128])/(6)/(pow(V[66],2))/(2);
   if(!isfinite(V[321]) || FError) return 321;
   V[322]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*pow(V[158],2)*V[128]*V[177]+4*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[158],2)*V[127]*V[123]-3*V[122]*pow(V[2],2)*pow(V[158],2)*V[127]*V[123]+3*pow(V[122],2)*pow(V[2],2)*pow(V[158],2)*V[128]-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*pow(V[158],2)*V[128]-4*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[159],2)*V[127]*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*pow(V[159],2)*V[128]+6*pow(V[39],2)*V[133]*V[158]*V[159]*V[127]*V[15]-6*pow(V[39],2)*V[10]*V[133]*V[158]*V[159]*V[128]-6*pow(V[39],2)*pow(V[133],2)*V[128])/(6)/(pow(V[69],2))/(2);
   if(!isfinite(V[322]) || FError) return 322;
   V[323]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*pow(V[160],2)*V[128]*V[177]+4*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[160],2)*V[127]*V[123]-3*V[122]*pow(V[2],2)*pow(V[160],2)*V[127]*V[123]+3*pow(V[122],2)*pow(V[2],2)*pow(V[160],2)*V[128]-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*pow(V[160],2)*V[128]-4*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[161],2)*V[127]*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*pow(V[161],2)*V[128]+6*pow(V[39],2)*V[133]*V[160]*V[161]*V[127]*V[15]-6*pow(V[39],2)*V[10]*V[133]*V[160]*V[161]*V[128]-6*pow(V[39],2)*pow(V[133],2)*V[128])/(6)/(pow(V[70],2))/(2);
   if(!isfinite(V[323]) || FError) return 323;
   V[324]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*pow(V[154],2)*V[128]*V[123]*V[178]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[154],2)*V[127]-3*pow(V[123],2)*pow(V[2],2)*pow(V[154],2)*V[127]+3*V[122]*pow(V[2],2)*pow(V[154],2)*V[128]*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[154],2)*V[128]*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[155],2)*V[127]+2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[155],2)*V[128]*V[123]-6*pow(V[39],2)*V[9]*V[154]*V[155]*V[128]*V[15]+6*pow(V[39],2)*V[11]*V[9]*V[154]*V[155]*V[127]+6*pow(V[39],2)*pow(V[9],2)*V[127])/(6)/(pow(V[63],2))/(2);
   if(!isfinite(V[324]) || FError) return 324;
   V[325]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*pow(V[156],2)*V[128]*V[123]*V[178]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[156],2)*V[127]-3*pow(V[123],2)*pow(V[2],2)*pow(V[156],2)*V[127]+3*V[122]*pow(V[2],2)*pow(V[156],2)*V[128]*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[156],2)*V[128]*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[157],2)*V[127]+2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[157],2)*V[128]*V[123]-6*pow(V[39],2)*V[9]*V[156]*V[157]*V[128]*V[15]+6*pow(V[39],2)*V[11]*V[9]*V[156]*V[157]*V[127]+6*pow(V[39],2)*pow(V[9],2)*V[127])/(6)/(pow(V[64],2))/(2);
   if(!isfinite(V[325]) || FError) return 325;
   V[326]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*V[128]*pow(V[162],2)*V[123]*V[177]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[127]*pow(V[162],2)-3*pow(V[123],2)*pow(V[2],2)*V[127]*pow(V[162],2)+3*V[122]*pow(V[2],2)*V[128]*pow(V[162],2)*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*V[128]*pow(V[162],2)*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[127]*pow(V[163],2)+2*pow(V[40],2)*V[122]*pow(V[2],2)*V[128]*pow(V[163],2)*V[123]-6*pow(V[39],2)*V[9]*V[128]*V[162]*V[163]*V[15]+6*pow(V[39],2)*V[11]*V[9]*V[127]*V[162]*V[163]+6*pow(V[39],2)*pow(V[9],2)*V[127])/(6)/(pow(V[67],2))/(2);
   if(!isfinite(V[326]) || FError) return 326;
   V[327]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*V[128]*pow(V[164],2)*V[123]*V[177]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[127]*pow(V[164],2)-3*pow(V[123],2)*pow(V[2],2)*V[127]*pow(V[164],2)+3*V[122]*pow(V[2],2)*V[128]*pow(V[164],2)*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*V[128]*pow(V[164],2)*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[127]*pow(V[165],2)+2*pow(V[40],2)*V[122]*pow(V[2],2)*V[128]*pow(V[165],2)*V[123]-6*pow(V[39],2)*V[9]*V[128]*V[164]*V[165]*V[15]+6*pow(V[39],2)*V[11]*V[9]*V[127]*V[164]*V[165]+6*pow(V[39],2)*pow(V[9],2)*V[127])/(6)/(pow(V[68],2))/(2);
   if(!isfinite(V[327]) || FError) return 327;
   V[328]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(V[42]*pow(V[123],2)*pow(V[2],2)*V[127]*pow(V[108],2)-V[42]*V[122]*pow(V[2],2)*V[128]*pow(V[108],2)*V[123]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[127]*pow(V[109],2)-2*pow(V[40],2)*V[122]*pow(V[2],2)*V[128]*pow(V[109],2)*V[123]+2*pow(V[39],2)*V[8]*V[128]*V[108]*V[109]*V[15]-2*pow(V[39],2)*V[36]*V[8]*V[127]*V[108]*V[109]-2*pow(V[39],2)*pow(V[8],2)*V[127])/(2)/(pow(V[61],2))/(2);
   if(!isfinite(V[328]) || FError) return 328;
   V[329]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(V[42]*pow(V[123],2)*pow(V[2],2)*V[127]*pow(V[110],2)-V[42]*V[122]*pow(V[2],2)*V[128]*pow(V[110],2)*V[123]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[127]*pow(V[111],2)-2*pow(V[40],2)*V[122]*pow(V[2],2)*V[128]*pow(V[111],2)*V[123]+2*pow(V[39],2)*V[8]*V[128]*V[110]*V[111]*V[15]-2*pow(V[39],2)*V[36]*V[8]*V[127]*V[110]*V[111]-2*pow(V[39],2)*pow(V[8],2)*V[127])/(2)/(pow(V[62],2))/(2);
   if(!isfinite(V[329]) || FError) return 329;
   V[330]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[100],2)*V[176]+4*pow(V[40],2)*V[122]*pow(V[2],2)*V[127]*pow(V[100],2)*V[123]-3*V[122]*pow(V[2],2)*V[127]*pow(V[100],2)*V[123]+3*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[100],2)-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[100],2)-4*pow(V[40],2)*V[122]*pow(V[2],2)*V[127]*pow(V[101],2)*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[101],2)+6*pow(V[39],2)*V[166]*V[127]*V[100]*V[101]*V[15]-6*pow(V[39],2)*V[174]*V[166]*V[128]*V[100]*V[101]-6*pow(V[39],2)*pow(V[166],2)*V[128])/(6)/(pow(V[73],2))/(2);
   if(!isfinite(V[330]) || FError) return 330;
   V[331]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[102],2)*V[176]+4*pow(V[40],2)*V[122]*pow(V[2],2)*V[127]*pow(V[102],2)*V[123]-3*V[122]*pow(V[2],2)*V[127]*pow(V[102],2)*V[123]+3*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[102],2)-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[102],2)-4*pow(V[40],2)*V[122]*pow(V[2],2)*V[127]*pow(V[103],2)*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[128]*pow(V[103],2)+6*pow(V[39],2)*V[166]*V[127]*V[102]*V[103]*V[15]-6*pow(V[39],2)*V[174]*V[166]*V[128]*V[102]*V[103]-6*pow(V[39],2)*pow(V[166],2)*V[128])/(6)/(pow(V[74],2))/(2);
   if(!isfinite(V[331]) || FError) return 331;
   V[332]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*pow(V[104],2)*V[128]*V[123]*V[176]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[104],2)*V[127]-3*pow(V[123],2)*pow(V[2],2)*pow(V[104],2)*V[127]+3*V[122]*pow(V[2],2)*pow(V[104],2)*V[128]*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[104],2)*V[128]*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[105],2)*V[127]+2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[105],2)*V[128]*V[123]-6*pow(V[39],2)*V[167]*V[104]*V[105]*V[128]*V[15]+6*pow(V[39],2)*V[175]*V[167]*V[104]*V[105]*V[127]+6*pow(V[39],2)*pow(V[167],2)*V[127])/(6)/(pow(V[71],2))/(2);
   if(!isfinite(V[332]) || FError) return 332;
   V[333]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*pow(V[106],2)*V[128]*V[123]*V[176]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[106],2)*V[127]-3*pow(V[123],2)*pow(V[2],2)*pow(V[106],2)*V[127]+3*V[122]*pow(V[2],2)*pow(V[106],2)*V[128]*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[106],2)*V[128]*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[107],2)*V[127]+2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[107],2)*V[128]*V[123]-6*pow(V[39],2)*V[167]*V[106]*V[107]*V[128]*V[15]+6*pow(V[39],2)*V[175]*V[167]*V[106]*V[107]*V[127]+6*pow(V[39],2)*pow(V[167],2)*V[127])/(6)/(pow(V[72],2))/(2);
   if(!isfinite(V[333]) || FError) return 333;
   V[334]=1/(pow(V[39],2))*V[7]*V[2]/(V[40])*(V[42]*V[125]*V[123]-V[42]*V[122]*V[126])/(2)/(pow(V[57],2))/(2);
   if(!isfinite(V[334]) || FError) return 334;
   V[335]=1/(pow(V[39],2))*V[7]*V[2]*V[40]*(V[125]*V[123]-V[122]*V[126])/(pow(V[58],2))/(2);
   if(!isfinite(V[335]) || FError) return 335;
   V[336]=1/(pow(V[39],2))*V[7]*V[2]/(V[40])*(V[42]*V[125]*V[123]-V[42]*V[122]*V[126])/(2)/(pow(V[59],2))/(2);
   if(!isfinite(V[336]) || FError) return 336;
   V[337]=1/(pow(V[39],2))*V[7]*V[2]*V[40]*(V[125]*V[123]-V[122]*V[126])/(pow(V[60],2))/(2);
   if(!isfinite(V[337]) || FError) return 337;
   V[338]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[150],2)*V[178]+4*pow(V[40],2)*V[122]*pow(V[2],2)*V[125]*pow(V[150],2)*V[123]-3*V[122]*pow(V[2],2)*V[125]*pow(V[150],2)*V[123]+3*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[150],2)-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[150],2)-4*pow(V[40],2)*V[122]*pow(V[2],2)*V[125]*pow(V[151],2)*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[151],2)+6*pow(V[39],2)*V[9]*V[125]*V[150]*V[151]*V[15]-6*pow(V[39],2)*V[10]*V[9]*V[126]*V[150]*V[151]-6*pow(V[39],2)*pow(V[9],2)*V[126])/(6)/(pow(V[65],2))/(2);
   if(!isfinite(V[338]) || FError) return 338;
   V[339]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[152],2)*V[178]+4*pow(V[40],2)*V[122]*pow(V[2],2)*V[125]*pow(V[152],2)*V[123]-3*V[122]*pow(V[2],2)*V[125]*pow(V[152],2)*V[123]+3*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[152],2)-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[152],2)-4*pow(V[40],2)*V[122]*pow(V[2],2)*V[125]*pow(V[153],2)*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[153],2)+6*pow(V[39],2)*V[9]*V[125]*V[152]*V[153]*V[15]-6*pow(V[39],2)*V[10]*V[9]*V[126]*V[152]*V[153]-6*pow(V[39],2)*pow(V[9],2)*V[126])/(6)/(pow(V[66],2))/(2);
   if(!isfinite(V[339]) || FError) return 339;
   V[340]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*pow(V[158],2)*V[126]*V[177]+4*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[158],2)*V[125]*V[123]-3*V[122]*pow(V[2],2)*pow(V[158],2)*V[125]*V[123]+3*pow(V[122],2)*pow(V[2],2)*pow(V[158],2)*V[126]-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*pow(V[158],2)*V[126]-4*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[159],2)*V[125]*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*pow(V[159],2)*V[126]+6*pow(V[39],2)*V[133]*V[158]*V[159]*V[125]*V[15]-6*pow(V[39],2)*V[10]*V[133]*V[158]*V[159]*V[126]-6*pow(V[39],2)*pow(V[133],2)*V[126])/(6)/(pow(V[69],2))/(2);
   if(!isfinite(V[340]) || FError) return 340;
   V[341]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*pow(V[160],2)*V[126]*V[177]+4*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[160],2)*V[125]*V[123]-3*V[122]*pow(V[2],2)*pow(V[160],2)*V[125]*V[123]+3*pow(V[122],2)*pow(V[2],2)*pow(V[160],2)*V[126]-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*pow(V[160],2)*V[126]-4*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[161],2)*V[125]*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*pow(V[161],2)*V[126]+6*pow(V[39],2)*V[133]*V[160]*V[161]*V[125]*V[15]-6*pow(V[39],2)*V[10]*V[133]*V[160]*V[161]*V[126]-6*pow(V[39],2)*pow(V[133],2)*V[126])/(6)/(pow(V[70],2))/(2);
   if(!isfinite(V[341]) || FError) return 341;
   V[342]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*pow(V[154],2)*V[126]*V[123]*V[178]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[154],2)*V[125]-3*pow(V[123],2)*pow(V[2],2)*pow(V[154],2)*V[125]+3*V[122]*pow(V[2],2)*pow(V[154],2)*V[126]*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[154],2)*V[126]*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[155],2)*V[125]+2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[155],2)*V[126]*V[123]-6*pow(V[39],2)*V[9]*V[154]*V[155]*V[126]*V[15]+6*pow(V[39],2)*V[11]*V[9]*V[154]*V[155]*V[125]+6*pow(V[39],2)*pow(V[9],2)*V[125])/(6)/(pow(V[63],2))/(2);
   if(!isfinite(V[342]) || FError) return 342;
   V[343]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*pow(V[156],2)*V[126]*V[123]*V[178]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[156],2)*V[125]-3*pow(V[123],2)*pow(V[2],2)*pow(V[156],2)*V[125]+3*V[122]*pow(V[2],2)*pow(V[156],2)*V[126]*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[156],2)*V[126]*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[157],2)*V[125]+2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[157],2)*V[126]*V[123]-6*pow(V[39],2)*V[9]*V[156]*V[157]*V[126]*V[15]+6*pow(V[39],2)*V[11]*V[9]*V[156]*V[157]*V[125]+6*pow(V[39],2)*pow(V[9],2)*V[125])/(6)/(pow(V[64],2))/(2);
   if(!isfinite(V[343]) || FError) return 343;
   V[344]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*V[126]*pow(V[162],2)*V[123]*V[177]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[125]*pow(V[162],2)-3*pow(V[123],2)*pow(V[2],2)*V[125]*pow(V[162],2)+3*V[122]*pow(V[2],2)*V[126]*pow(V[162],2)*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*V[126]*pow(V[162],2)*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[125]*pow(V[163],2)+2*pow(V[40],2)*V[122]*pow(V[2],2)*V[126]*pow(V[163],2)*V[123]-6*pow(V[39],2)*V[9]*V[126]*V[162]*V[163]*V[15]+6*pow(V[39],2)*V[11]*V[9]*V[125]*V[162]*V[163]+6*pow(V[39],2)*pow(V[9],2)*V[125])/(6)/(pow(V[67],2))/(2);
   if(!isfinite(V[344]) || FError) return 344;
   V[345]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*V[126]*pow(V[164],2)*V[123]*V[177]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[125]*pow(V[164],2)-3*pow(V[123],2)*pow(V[2],2)*V[125]*pow(V[164],2)+3*V[122]*pow(V[2],2)*V[126]*pow(V[164],2)*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*V[126]*pow(V[164],2)*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[125]*pow(V[165],2)+2*pow(V[40],2)*V[122]*pow(V[2],2)*V[126]*pow(V[165],2)*V[123]-6*pow(V[39],2)*V[9]*V[126]*V[164]*V[165]*V[15]+6*pow(V[39],2)*V[11]*V[9]*V[125]*V[164]*V[165]+6*pow(V[39],2)*pow(V[9],2)*V[125])/(6)/(pow(V[68],2))/(2);
   if(!isfinite(V[345]) || FError) return 345;
   V[346]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(V[42]*pow(V[123],2)*pow(V[2],2)*V[125]*pow(V[108],2)-V[42]*V[122]*pow(V[2],2)*V[126]*pow(V[108],2)*V[123]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[125]*pow(V[109],2)-2*pow(V[40],2)*V[122]*pow(V[2],2)*V[126]*pow(V[109],2)*V[123]+2*pow(V[39],2)*V[8]*V[126]*V[108]*V[109]*V[15]-2*pow(V[39],2)*V[36]*V[8]*V[125]*V[108]*V[109]-2*pow(V[39],2)*pow(V[8],2)*V[125])/(2)/(pow(V[61],2))/(2);
   if(!isfinite(V[346]) || FError) return 346;
   V[347]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(V[42]*pow(V[123],2)*pow(V[2],2)*V[125]*pow(V[110],2)-V[42]*V[122]*pow(V[2],2)*V[126]*pow(V[110],2)*V[123]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*V[125]*pow(V[111],2)-2*pow(V[40],2)*V[122]*pow(V[2],2)*V[126]*pow(V[111],2)*V[123]+2*pow(V[39],2)*V[8]*V[126]*V[110]*V[111]*V[15]-2*pow(V[39],2)*V[36]*V[8]*V[125]*V[110]*V[111]-2*pow(V[39],2)*pow(V[8],2)*V[125])/(2)/(pow(V[62],2))/(2);
   if(!isfinite(V[347]) || FError) return 347;
   V[348]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[100],2)*V[176]+4*pow(V[40],2)*V[122]*pow(V[2],2)*V[125]*pow(V[100],2)*V[123]-3*V[122]*pow(V[2],2)*V[125]*pow(V[100],2)*V[123]+3*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[100],2)-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[100],2)-4*pow(V[40],2)*V[122]*pow(V[2],2)*V[125]*pow(V[101],2)*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[101],2)+6*pow(V[39],2)*V[166]*V[125]*V[100]*V[101]*V[15]-6*pow(V[39],2)*V[174]*V[166]*V[126]*V[100]*V[101]-6*pow(V[39],2)*pow(V[166],2)*V[126])/(6)/(pow(V[73],2))/(2);
   if(!isfinite(V[348]) || FError) return 348;
   V[349]=1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[122])*(3*pow(V[39],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[102],2)*V[176]+4*pow(V[40],2)*V[122]*pow(V[2],2)*V[125]*pow(V[102],2)*V[123]-3*V[122]*pow(V[2],2)*V[125]*pow(V[102],2)*V[123]+3*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[102],2)-4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[102],2)-4*pow(V[40],2)*V[122]*pow(V[2],2)*V[125]*pow(V[103],2)*V[123]+4*pow(V[40],2)*pow(V[122],2)*pow(V[2],2)*V[126]*pow(V[103],2)+6*pow(V[39],2)*V[166]*V[125]*V[102]*V[103]*V[15]-6*pow(V[39],2)*V[174]*V[166]*V[126]*V[102]*V[103]-6*pow(V[39],2)*pow(V[166],2)*V[126])/(6)/(pow(V[74],2))/(2);
   if(!isfinite(V[349]) || FError) return 349;
   V[350]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*pow(V[104],2)*V[126]*V[123]*V[176]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[104],2)*V[125]-3*pow(V[123],2)*pow(V[2],2)*pow(V[104],2)*V[125]+3*V[122]*pow(V[2],2)*pow(V[104],2)*V[126]*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[104],2)*V[126]*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[105],2)*V[125]+2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[105],2)*V[126]*V[123]-6*pow(V[39],2)*V[167]*V[104]*V[105]*V[126]*V[15]+6*pow(V[39],2)*V[175]*V[167]*V[104]*V[105]*V[125]+6*pow(V[39],2)*pow(V[167],2)*V[125])/(6)/(pow(V[71],2))/(2);
   if(!isfinite(V[350]) || FError) return 350;
   V[351]=-1/(pow(V[39],2))*V[7]/(V[2])/(V[40])/(V[123])*(3*pow(V[39],2)*V[122]*pow(V[2],2)*pow(V[106],2)*V[126]*V[123]*V[176]+2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[106],2)*V[125]-3*pow(V[123],2)*pow(V[2],2)*pow(V[106],2)*V[125]+3*V[122]*pow(V[2],2)*pow(V[106],2)*V[126]*V[123]-2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[106],2)*V[126]*V[123]-2*pow(V[40],2)*pow(V[123],2)*pow(V[2],2)*pow(V[107],2)*V[125]+2*pow(V[40],2)*V[122]*pow(V[2],2)*pow(V[107],2)*V[126]*V[123]-6*pow(V[39],2)*V[167]*V[106]*V[107]*V[126]*V[15]+6*pow(V[39],2)*V[175]*V[167]*V[106]*V[107]*V[125]+6*pow(V[39],2)*pow(V[167],2)*V[125])/(6)/(pow(V[72],2))/(2);
   if(!isfinite(V[351]) || FError) return 351;
   V[352]=V[7]*V[2]/(V[40])*(V[127]*V[123]+V[122]*V[128])/(pow(V[2],2))/(2);
   if(!isfinite(V[352]) || FError) return 352;
   V[353]=pow(V[122],2)*V[127]*V[123]*V[117]+pow(V[123],2)*V[122]*V[128]*V[143]+pow(V[123],2)*V[127]*V[123]*V[118]+pow(V[122],3)*V[128]*V[118]-pow(V[123],2)*V[122]*V[128]*V[141];
   if(!isfinite(V[353]) || FError) return 353;
   V[354]=V[353]-pow(V[122],2)*V[127]*V[123]*V[141]-pow(V[123],2)*V[122]*V[128]*V[140]-pow(V[122],2)*V[127]*V[123]*V[140]-3*pow(V[122],3)*V[127]*V[119]+2*V[122]*V[127]*V[119];
   if(!isfinite(V[354]) || FError) return 354;
   V[355]=V[354]-pow(V[122],2)*V[128]*V[123]*V[119]+3*pow(V[122],2)*V[128]*V[123]*V[142]-V[128]*V[123]*V[142]-pow(V[123],2)*V[122]*V[127]*V[142];
   if(!isfinite(V[355]) || FError) return 355;
   V[356]=-2/(V[7])*V[2]*V[40]*V[355]/(pow(V[46],2))/(2);
   if(!isfinite(V[356]) || FError) return 356;
   V[357]=V[7]*V[2]/(V[40])*(V[125]*V[123]+V[122]*V[126])/(pow(V[2],2))/(2);
   if(!isfinite(V[357]) || FError) return 357;
   V[358]=pow(V[122],2)*V[125]*V[123]*V[117]+pow(V[123],2)*V[122]*V[126]*V[143]+pow(V[123],2)*V[125]*V[123]*V[118]+pow(V[122],3)*V[126]*V[118]-pow(V[123],2)*V[122]*V[126]*V[141];
   if(!isfinite(V[358]) || FError) return 358;
   V[359]=V[358]-pow(V[122],2)*V[125]*V[123]*V[141]-pow(V[123],2)*V[122]*V[126]*V[140]-pow(V[122],2)*V[125]*V[123]*V[140]-3*pow(V[122],3)*V[125]*V[119]+2*V[122]*V[125]*V[119];
   if(!isfinite(V[359]) || FError) return 359;
   V[360]=V[359]-pow(V[122],2)*V[126]*V[123]*V[119]+3*pow(V[122],2)*V[126]*V[123]*V[142]-V[126]*V[123]*V[142]-pow(V[123],2)*V[122]*V[125]*V[142];
   if(!isfinite(V[360]) || FError) return 360;
   V[361]=-2/(V[7])*V[2]*V[40]*V[360]/(pow(V[46],2))/(2);
   if(!isfinite(V[361]) || FError) return 361;
   V[362]=-V[7]/(V[2])/(V[40])*M_SQRT2/(V[123])/(V[122])*(V[122]*V[185]*V[127]*V[93]*V[96]-V[186]*V[128]*V[92]*V[97]*V[123])/(2)/(V[51]);
   if(!isfinite(V[362]) || FError) return 362;
   V[363]=-V[7]/(V[2])/(V[40])*M_SQRT2/(V[123])/(V[122])*(V[122]*V[185]*V[127]*V[95]*V[98]-V[186]*V[128]*V[94]*V[99]*V[123])/(2)/(V[52]);
   if(!isfinite(V[363]) || FError) return 363;
   V[364]=-V[7]/(V[2])/(V[40])*M_SQRT2/(V[123])/(V[122])*(V[122]*V[185]*V[125]*V[93]*V[96]-V[186]*V[126]*V[92]*V[97]*V[123])/(2)/(V[51]);
   if(!isfinite(V[364]) || FError) return 364;
   V[365]=-V[7]/(V[2])/(V[40])*M_SQRT2/(V[123])/(V[122])*(V[122]*V[185]*V[125]*V[95]*V[98]-V[186]*V[126]*V[94]*V[99]*V[123])/(2)/(V[52]);
   if(!isfinite(V[365]) || FError) return 365;
   V[366]=V[7]/(V[2])*V[133]/(V[40])/(V[113])/(2)/(V[133])/(2);
   if(!isfinite(V[366]) || FError) return 366;
   V[367]=-V[7]/(V[2])*V[131]/(V[40])/(V[123])/(V[122])*(V[129]-pow(V[122],2))/(2)/(V[131])/(2);
   if(!isfinite(V[367]) || FError) return 367;
   V[368]=V[7]/(V[2])*V[132]/(V[40])/(V[113])/(2)/(V[132])/(2);
   if(!isfinite(V[368]) || FError) return 368;
   V[369]=-V[7]/(V[2])*V[8]/(V[40])/(V[123])/(V[122])*(V[179]-pow(V[122],2))/(2)/(V[8])/(2);
   if(!isfinite(V[369]) || FError) return 369;
   V[370]=V[7]/(V[2])/(V[40])*M_SQRT2/(V[123])/(V[122])*(pow(V[123],2)*V[186]*V[92]*V[97]-pow(V[122],2)*V[185]*V[93]*V[96])/(2)/(V[51])/(2);
   if(!isfinite(V[370]) || FError) return 370;
   V[371]=V[7]/(V[2])/(V[40])*M_SQRT2/(V[123])/(V[122])*(pow(V[123],2)*V[186]*V[94]*V[99]-pow(V[122],2)*V[185]*V[95]*V[98])/(2)/(V[52])/(2);
   if(!isfinite(V[371]) || FError) return 371;
   V[372]=1+V[192]*(149/(double)((12))+V[192]*(68.6482-V[192]*212.447));
   if(!isfinite(V[372]) || FError) return 372;
   V[373]=2*log(V[3]/(V[6]));
   if(!isfinite(V[373]) || FError) return 373;
   V[374]=1+V[192]*(11/(double)((4))+V[192]*(6.1537-2.8542*V[373]+V[192]*(10.999-17.93*V[373]+5.47*pow(V[373],2))));
   if(!isfinite(V[374]) || FError) return 374;
   V[375]=1+11/(double)((4))*V[192];
   if(!isfinite(V[375]) || FError) return 375;
   V[376]=1+9/(double)((2))*V[192];
   if(!isfinite(V[376]) || FError) return 376;
   V[377]=1/(137.036);
   if(!isfinite(V[377]) || FError) return 377;
   V[378]=2/(double)((3));
   if(!isfinite(V[378]) || FError) return 378;
   V[379]=-1/(double)((3));
   if(!isfinite(V[379]) || FError) return 379;
   V[380]=1+V[192]*(221/(double)((12))+V[192]*(171.5-5*V[373]));
   if(!isfinite(V[380]) || FError) return 380;
   V[381]=V[309]*(1+V[114])/(1-V[114]*V[120]/(V[121])/(V[113]));
   if(!isfinite(V[381]) || FError) return 381;
   V[382]=V[313]*(1+V[114])/(1+V[114]*V[121]/(V[120])/(V[113]));
   if(!isfinite(V[382]) || FError) return 382;
   V[383]=V[367]*(1+V[114])/(1-V[114]/(pow(V[113],2)));
   if(!isfinite(V[383]) || FError) return 383;
   V[384]=V[311]*(1+V[115])/(1-V[115]*V[120]/(V[121])/(V[113]));
   if(!isfinite(V[384]) || FError) return 384;
   V[385]=V[315]*(1+V[115])/(1+V[115]*V[121]/(V[120])/(V[113]));
   if(!isfinite(V[385]) || FError) return 385;
   V[386]=V[369]*(1+V[115])/(1-V[115]/(pow(V[113],2)));
   if(!isfinite(V[386]) || FError) return 386;
   V[387]=-V[192]/(8)*1/(double)((2))*sqrt(V[372])*cabs(((V[197]+I*V[198])*V[381]+(V[193]+I*V[194])*V[308])*V[375]+(V[201]+I*V[202])*V[310]*V[374]+((V[221]+I*V[222])*V[320]+(V[223]+I*V[224])*V[321]+(V[233]+I*V[234])*V[324]+(V[235]+I*V[236])*V[325]+(V[237]+I*V[238])*V[326]+(V[239]+I*V[240])*V[327]+(V[225]+I*V[226])*V[322]+(V[227]+I*V[228])*V[323]+(V[241]+I*V[242])*V[332]+(V[243]+I*V[244])*V[333]+(V[229]+I*V[230])*V[330]+(V[231]+I*V[232])*V[331])*V[376]);
   if(!isfinite(V[387]) || FError) return 387;
   V[388]=-V[377]/(8*V[134])*cabs(3*pow(V[379],2)*(V[264]+I*V[265])*V[381]+3*pow(V[378],2)*((V[268]+I*V[269])*V[310]+(V[260]+I*V[261])*V[308])+3*pow(V[378],2)*((V[273]+I*V[274])*V[320]+(V[276]+I*V[277])*V[321]+(V[279]+I*V[280])*V[322]+(V[282]+I*V[283])*V[323]+(V[285]+I*V[286])*V[330]+(V[288]+I*V[289])*V[331])+3*pow(V[379],2)*((V[291]+I*V[292])*V[324]+(V[294]+I*V[295])*V[325]+(V[297]+I*V[298])*V[326]+(V[300]+I*V[301])*V[327]+(V[303]+I*V[304])*V[332]+(V[306]+I*V[307])*V[333])+(V[209]+I*V[210])*V[316]+(V[211]+I*V[212])*V[317]+(V[213]+I*V[214])*V[318]+(V[215]+I*V[216])*V[319]+(V[217]+I*V[218])*V[328]+(V[219]+I*V[220])*V[329]+(V[205]+I*V[206])*V[384]+(V[249]+I*V[250])*V[362]+(V[253]+I*V[254])*V[363]-(V[247]+I*V[248])*V[352]+(V[245]+I*V[246])*V[356]);
   if(!isfinite(V[388]) || FError) return 388;
   V[389]=-V[192]/(8)*1/(double)((2))*sqrt(V[372])*cabs(((V[197]+I*V[198])*V[382]+(V[193]+I*V[194])*V[312])*V[375]+(V[201]+I*V[202])*V[314]*V[374]+((V[221]+I*V[222])*V[338]+(V[223]+I*V[224])*V[339]+(V[233]+I*V[234])*V[342]+(V[235]+I*V[236])*V[343]+(V[237]+I*V[238])*V[344]+(V[239]+I*V[240])*V[345]+(V[225]+I*V[226])*V[340]+(V[227]+I*V[228])*V[341]+(V[241]+I*V[242])*V[350]+(V[243]+I*V[244])*V[351]+(V[229]+I*V[230])*V[348]+(V[231]+I*V[232])*V[349])*V[376]);
   if(!isfinite(V[389]) || FError) return 389;
   V[390]=-V[377]/(8*V[134])*cabs(3*pow(V[379],2)*(V[264]+I*V[265])*V[382]+3*pow(V[378],2)*((V[268]+I*V[269])*V[314]+(V[260]+I*V[261])*V[312])+3*pow(V[378],2)*((V[273]+I*V[274])*V[338]+(V[276]+I*V[277])*V[339]+(V[279]+I*V[280])*V[340]+(V[282]+I*V[283])*V[341]+(V[285]+I*V[286])*V[348]+(V[288]+I*V[289])*V[349])+3*pow(V[379],2)*((V[291]+I*V[292])*V[342]+(V[294]+I*V[295])*V[343]+(V[297]+I*V[298])*V[344]+(V[300]+I*V[301])*V[345]+(V[303]+I*V[304])*V[350]+(V[306]+I*V[307])*V[351])+(V[209]+I*V[210])*V[334]+(V[211]+I*V[212])*V[335]+(V[213]+I*V[214])*V[336]+(V[215]+I*V[216])*V[337]+(V[217]+I*V[218])*V[346]+(V[219]+I*V[220])*V[347]+(V[205]+I*V[206])*V[385]+(V[249]+I*V[250])*V[364]+(V[253]+I*V[254])*V[365]-(V[247]+I*V[248])*V[357]+(V[245]+I*V[246])*V[361]);
   if(!isfinite(V[390]) || FError) return 390;
   V[391]=-V[192]/(8)*1/(double)((2))*sqrt(V[380])*cabs((V[195]+I*V[196])*V[366]+(V[199]+I*V[200])*V[383]+(V[203]+I*V[204])*V[368]);
   if(!isfinite(V[391]) || FError) return 391;
   V[392]=-V[377]/(8*V[134])*cabs(3*pow(V[379],2)*(V[266]+I*V[267])*V[383]+3*pow(V[378],2)*((V[270]+I*V[271])*V[368]+(V[262]+I*V[263])*V[366])+(V[207]+I*V[208])*V[386]+(V[251]+I*V[252])*V[370]+(V[255]+I*V[256])*V[371]);
   if(!isfinite(V[392]) || FError) return 392;
   V[393]=-V[7]/(V[2])/(V[40])/(2);
   if(!isfinite(V[393]) || FError) return 393;
   V[394]=V[7]/(V[40])/(V[2])/(2);
   if(!isfinite(V[394]) || FError) return 394;
   V[395]=-V[192]/(8)*1/(double)((2))*sqrt(V[372])*fabs(V[393])*cabs((V[197]+I*V[198]+V[193]+I*V[194])*V[375]+(V[201]+I*V[202])*V[374]);
   if(!isfinite(V[395]) || FError) return 395;
   V[396]=-V[377]/(8*V[134])*cabs(3*pow(V[379],2)*(V[264]+I*V[265])*V[393]+3*pow(V[378],2)*((V[268]+I*V[269])*V[393]+(V[260]+I*V[261])*V[393])+(V[205]+I*V[206])*V[393]-(V[247]+I*V[248])*V[394]);
   if(!isfinite(V[396]) || FError) return 396;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   return 0;
}
