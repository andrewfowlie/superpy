#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../sources/micromegas/CalcHEP_src/include/extern.h"
#include "../../sources/micromegas/CalcHEP_src/include/VandP.h"
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
, {"~o1","~o1", 1000022, "MNE1","wNE1",1,1,0}
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
int nModelVars=158;
int nModelFunc=216;
static char*varNames_[374]={
 "EE","SW","MZ","alfSMZ","MbMb","McMc","Mtp","Ml","wt","Mq"
,"MG1","MG2","MG3","tb","mu","Lambda","Kappa","aLambda","aKappa","At"
,"Ab","Al","Ml2","Ml3","Mr2","Mr3","Mq2","Mq3","Mu2","Mu3"
,"Md2","Md3","Mh1","Mh2","Mh3","Mha","Mhb","MHc","MNE1","MNE2"
,"MNE3","MNE4","MNE5","MC1","MC2","MSuL","MSuR","MSdL","MSdR","MScL"
,"MScR","MSsL","MSsR","MSt1","MSt2","MSb1","MSb2","MSeL","MSeR","MSmL"
,"MSmR","MSl1","MSl2","MSne","MSnm","MSnl","Zn11","Zn12","Zn13","Zn14"
,"Zn15","Zn21","Zn22","Zn23","Zn24","Zn25","Zn31","Zn32","Zn33","Zn34"
,"Zn35","Zn41","Zn42","Zn43","Zn44","Zn45","Zn51","Zn52","Zn53","Zn54"
,"Zn55","Zu11","Zu12","Zu21","Zu22","Zv11","Zv12","Zv21","Zv22","Zt11"
,"Zt12","Zt21","Zt22","Zb11","Zb12","Zb21","Zb22","Zl11","Zl12","Zl21"
,"Zl22","Zh11","Zh12","Zh13","Zh21","Zh22","Zh23","Zh31","Zh32","Zh33"
,"Za11","Za12","Za13","Za21","Za22","Za23","la1","la2","la3","la4"
,"la5","la6","la7","la1s","la2s","la3s","la4s","la5s","la6s","la7s"
,"la8s","aa1","aa2","aa3","aa4","aa5","aa6","B1","B2","X"
,"dMb","Q","Au","Ad","Maux","tB","MSG","wNE1","CW","C2W"
,"MW","LamQCD","sb","cb","t2b","xvev","Pa12","Pa22","Pa11","Pa21"
,"Td3","hMM11","hMM12","hMM13","hMM22","hMM23","hMM33","idS","Mh1_","Mh2_"
,"Mh3_","Zh11_","Zh12_","Zh13_","Zh21_","Zh22_","Zh23_","Zh31_","Zh32_","Zh33_"
,"MA11","MA12","MA22","idA","Pa11_","Pa12_","Pa21_","Pa22_","Mha_","Mhb_"
,"Mt","Mb","Mc","xH2","yH2","dMd","Td2","fiuu","fidd","ficc"
,"fiss","Zuu11","Zuu12","Zuu21","Zuu22","Zdd11","Zdd12","Zdd21","Zdd22","Zcc11"
,"Zcc12","Zcc21","Zcc22","Zss11","Zss12","Zss21","Zss22","StMM11","StMM12","StMM22"
,"MtMM","AT","SbMM11","SbMM12","SbMM22","MbMM","AB","dX3","dX2","dX1"
,"NMM11","NMM12","NMM13","NMM14","NMM15","NMM22","NMM23","NMM24","NMM25","NMM33"
,"NMM34","NMM35","NMM44","NMM45","NMM55","MG1I","MG2I","PI","aQCD","rhF_c"
,"ihF_c","rAF_c","iAF_c","rhF_b","ihF_b","rAF_b","iAF_b","rhF_t","ihF_t","rAF_t"
,"iAF_t","rhF_l","ihF_l","rAF_l","iAF_l","rhS_eL","ihS_eL","rhS_eR","ihS_eR","rhS_mL"
,"ihS_mL","rhS_mR","ihS_mR","rhS_l1","ihS_l1","rhS_l2","ihS_l2","rhS_uL","ihS_uL","rhS_uR"
,"ihS_uR","rhS_cL","ihS_cL","rhS_cR","ihS_cR","rhS_t1","ihS_t1","rhS_t2","ihS_t2","rhS_dL"
,"ihS_dL","rhS_dR","ihS_dR","rhS_sL","ihS_sL","rhS_sR","ihS_sR","rhS_b1","ihS_b1","rhS_b2"
,"ihS_b2","rhS_Hc","ihS_Hc","rhV_W","ihV_W","rhF_c1","ihF_c1","rAF_c1","iAF_c1","rhF_c2"
,"ihF_c2","rAF_c2","iAF_c2","McR","MbR","MtR","rhF1_c","ihF1_c","rAF1_c","iAF1_c"
,"rhF1_b","ihF1_b","rAF1_b","iAF1_b","rhF1_t","ihF1_t","rAF1_t","iAF1_t","tau2_uL","rhS1_uL"
,"ihS1_uL","tau2_uR","rhS1_uR","ihS1_uR","tau2_cL","rhS1_cL","ihS1_cL","tau2_cR","rhS1_cR","ihS1_cR"
,"tau2_t1","rhS1_t1","ihS1_t1","tau2_t2","rhS1_t2","ihS1_t2","tau2_dL","rhS1_dL","ihS1_dL","tau2_dR"
,"rhS1_dR","ihS1_dR","tau2_sL","rhS1_sL","ihS1_sL","tau2_sR","rhS1_sR","ihS1_sR","tau2_b1","rhS1_b1"
,"ihS1_b1","tau2_b2","rhS1_b2","ihS1_b2"};
char**varNames=varNames_;
static REAL varValues_[374]={
   3.128530E-01,  4.820000E-01,  9.120000E+01,  1.184000E-01,  4.500000E+00,  1.300000E+00,  1.745000E+02,  1.777000E+00,  1.442000E+00,  5.000000E-02
,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
,  0.000000E+00,  0.000000E+00,  2.000000E+02,  2.000000E+02,  2.020000E+02,  2.000000E+02,  1.000000E+03,  1.000000E+03,  1.000000E+03,  1.000000E+03
,  1.000000E+03,  1.000000E+03,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00
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
,  0.000000E+00,  1.000000E+02,  0.000000E+00,  0.000000E+00,  1.000000E+00,  0.000000E+00,  0.000000E+00,  0.000000E+00};
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
   V[158]=sqrt(1-pow(V[1],2));
   if(!finite(V[158]) || FError) return 158;
   V[159]=pow(V[158],2)-pow(V[1],2);
   if(!finite(V[159]) || FError) return 159;
   V[160]=V[2]*V[158];
   if(!finite(V[160]) || FError) return 160;
   V[161]=initQCD5(V[3],V[5],V[4],V[6]);
   if(!finite(V[161]) || FError) return 161;
   V[162]=V[155]/(sqrt(1+pow(V[155],2)));
   if(!finite(V[162]) || FError) return 162;
   V[163]=sqrt(1-pow(V[162],2));
   if(!finite(V[163]) || FError) return 163;
   V[164]=2*V[155]/(1-pow(V[155],2));
   if(!finite(V[164]) || FError) return 164;
   V[165]=V[14]/(V[15]);
   if(!finite(V[165]) || FError) return 165;
   V[166]=V[122];

   V[167]=V[125];

   V[168]=(V[120]*V[167]>0 ? V[167] : -V[167]);
   if(!finite(V[168]) || FError) return 168;
   V[169]=(V[123]*V[166]>0 ? V[166] : -V[166]);
   if(!finite(V[169]) || FError) return 169;
   V[170]=V[150]/(1+V[150]);
   if(!finite(V[170]) || FError) return 170;
   V[171]=4/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[163]*V[163]*V[126]-6/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[163]*V[131]*V[162]+1*V[147]/(V[163])*V[162]+1*V[145]/(V[163])*V[162]*V[165]+1*V[146]/(V[163])*V[162]*V[165]+1/(V[163])*V[136]*V[162]*pow(V[165],2)+1/(V[163])*V[137]*V[162]*pow(V[165],2)+1/(V[163])*V[138]*V[162]*pow(V[165],2)+2/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)/(V[163])*V[132]*pow(V[162],3);
   if(!finite(V[171]) || FError) return 171;
   V[172]=4/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[163]*V[128]*V[162]+4/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[163]*V[129]*V[162]+4/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[163]*V[130]*V[162]-6/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[163]*V[163]*V[131]-6/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[132]*pow(V[162],2)-1*V[136]*pow(V[165],2)-1*V[137]*pow(V[165],2)-1*V[138]*pow(V[165],2)-1*V[145]*V[165]-1*V[146]*V[165]-1*V[147];
   if(!finite(V[172]) || FError) return 172;
   V[173]=2/(V[0])*V[160]*V[1]*M_SQRT2*V[163]*V[133]*V[165]-2/(V[0])*V[160]*V[1]*M_SQRT2*V[136]*V[162]*V[165]-2/(V[0])*V[160]*V[1]*M_SQRT2*V[137]*V[162]*V[165]-2/(V[0])*V[160]*V[1]*M_SQRT2*V[138]*V[162]*V[165]+2/(V[0])*V[160]*V[1]*M_SQRT2*V[141]*V[163]-1/(V[0])*V[160]*V[1]*M_SQRT2*V[145]*V[162]-1/(V[0])*V[160]*V[1]*M_SQRT2*V[146]*V[162];
   if(!finite(V[173]) || FError) return 173;
   V[174]=4/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[127]*pow(V[162],2)-6/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[163]*V[132]*V[162]+1*V[147]*V[163]/(V[162])+1*V[145]*V[163]/(V[162])*V[165]+1*V[146]*V[163]/(V[162])*V[165]+2/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*pow(V[163],3)*V[131]/(V[162])+1*V[163]*V[136]/(V[162])*pow(V[165],2)+1*V[163]*V[137]/(V[162])*pow(V[165],2)+1*V[163]*V[138]/(V[162])*pow(V[165],2);
   if(!finite(V[174]) || FError) return 174;
   V[175]=2/(V[0])*V[160]*V[1]*M_SQRT2*V[134]*V[162]*V[165]-2/(V[0])*V[160]*V[1]*M_SQRT2*V[163]*V[136]*V[165]-2/(V[0])*V[160]*V[1]*M_SQRT2*V[163]*V[137]*V[165]-2/(V[0])*V[160]*V[1]*M_SQRT2*V[163]*V[138]*V[165]+2/(V[0])*V[160]*V[1]*M_SQRT2*V[142]*V[162]-1/(V[0])*V[160]*V[1]*M_SQRT2*V[145]*V[163]-1/(V[0])*V[160]*V[1]*M_SQRT2*V[146]*V[163];
   if(!finite(V[175]) || FError) return 175;
   V[176]=4*V[135]*pow(V[165],2)+8*V[139]*pow(V[165],2)+8*V[140]*pow(V[165],2)+3*V[143]*V[165]+3*V[144]*V[165]-1*V[149]/(V[165])-2/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[141]*V[163]*V[163]/(V[165])-2/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[142]*pow(V[162],2)/(V[165])+2/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[145]*V[163]*V[162]/(V[165])+2/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[146]*V[163]*V[162]/(V[165]);
   if(!finite(V[176]) || FError) return 176;
   V[177]=rDiagonal(3,V[171],V[172],V[173],V[174],V[175],V[176]);
   if(!finite(V[177]) || FError) return 177;
   V[178]=MassArray(V[177],1)/(sqrt(fabs(MassArray(V[177],1))));
   if(!finite(V[178]) || FError) return 178;
   V[179]=MassArray(V[177],2)/(sqrt(fabs(MassArray(V[177],2))));
   if(!finite(V[179]) || FError) return 179;
   V[180]=MassArray(V[177],3)/(sqrt(fabs(MassArray(V[177],3))));
   if(!finite(V[180]) || FError) return 180;
   V[181]=MixMatrix(V[177],1,1);
   if(!finite(V[181]) || FError) return 181;
   V[182]=MixMatrix(V[177],1,2);
   if(!finite(V[182]) || FError) return 182;
   V[183]=MixMatrix(V[177],1,3);
   if(!finite(V[183]) || FError) return 183;
   V[184]=MixMatrix(V[177],2,1);
   if(!finite(V[184]) || FError) return 184;
   V[185]=MixMatrix(V[177],2,2);
   if(!finite(V[185]) || FError) return 185;
   V[186]=MixMatrix(V[177],2,3);
   if(!finite(V[186]) || FError) return 186;
   V[187]=MixMatrix(V[177],3,1);
   if(!finite(V[187]) || FError) return 187;
   V[188]=MixMatrix(V[177],3,2);
   if(!finite(V[188]) || FError) return 188;
   V[189]=MixMatrix(V[177],3,3);
   if(!finite(V[189]) || FError) return 189;
   V[190]=2/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[163]*V[131]/(V[162])-4/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)*V[130]+2/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)/(V[163])*V[132]*V[162]+1*V[147]/(V[163])/(V[162])+1*V[145]/(V[163])/(V[162])*V[165]+1*V[146]/(V[163])/(V[162])*V[165]+1/(V[163])*V[136]/(V[162])*pow(V[165],2)+1/(V[163])*V[137]/(V[162])*pow(V[165],2)+1/(V[163])*V[138]/(V[162])*pow(V[165],2);
   if(!finite(V[190]) || FError) return 190;
   V[191]=(2*V[137]*V[165]-2*V[138]*V[165]+V[145]-V[146])*V[160]*V[1]*M_SQRT2/(V[0]);
   if(!finite(V[191]) || FError) return 191;
   V[192]=(4*V[163]*V[137]*V[162]+4*V[163]*V[138]*V[162]+(V[145]*V[163]*V[162]-V[142]*pow(V[162],2)-V[141]*V[163]*V[163]+V[146]*V[163]*V[162])/(V[165]))*2/(pow(V[0],2))*pow(V[160],2)*pow(V[1],2)-16*V[139]*pow(V[165],2)-1*V[143]*V[165]-9*V[144]*V[165]-4*V[148]-1*V[149]/(V[165])-4*V[140]*pow(V[165],2);
   if(!finite(V[192]) || FError) return 192;
   V[193]=rDiagonal(2,V[190],V[191],V[192]);
   if(!finite(V[193]) || FError) return 193;
   V[194]=MixMatrix(V[193],1,1);
   if(!finite(V[194]) || FError) return 194;
   V[195]=MixMatrix(V[193],1,2);
   if(!finite(V[195]) || FError) return 195;
   V[196]=MixMatrix(V[193],2,1);
   if(!finite(V[196]) || FError) return 196;
   V[197]=MixMatrix(V[193],2,2);
   if(!finite(V[197]) || FError) return 197;
   V[198]=MassArray(V[193],1)/(sqrt(fabs(MassArray(V[193],1))));
   if(!finite(V[198]) || FError) return 198;
   V[199]=MassArray(V[193],2)/(sqrt(fabs(MassArray(V[193],2))));
   if(!finite(V[199]) || FError) return 199;
 FirstQ:
 cErr=1;
   V[200]=MtEff(V[151]);
   if(!finite(V[200]) || FError) return 200;
   V[201]=MbEff(V[151]);
   if(!finite(V[201]) || FError) return 201;
   V[202]=McEff(V[151]);
   if(!finite(V[202]) || FError) return 202;
   V[203]=pow(V[51],2)/(pow(V[156],2));
   if(!finite(V[203]) || FError) return 203;
   V[204]=pow(V[52],2)/(pow(V[156],2));
   if(!finite(V[204]) || FError) return 204;
   V[205]=-3/(double)((2))*alphaQCD(V[156])/(3.1415)*V[155]*V[14]/(V[156])*(V[203]*log(V[203])/(1-V[203])-V[204]*log(V[204])/(1-V[204]))/(V[203]-V[204]);
   if(!finite(V[205]) || FError) return 205;
   V[206]=V[205]/(1+V[205]);
   if(!finite(V[206]) || FError) return 206;
   V[207]=atan(-2*V[9]*(V[152]-V[14]/(V[155]))/(pow(V[46],2)-pow(V[45],2)))/(2);
   if(!finite(V[207]) || FError) return 207;
   V[208]=atan(-2*V[9]*(1-V[206])*(V[153]-V[14]*V[155])/(pow(V[48],2)-pow(V[47],2)))/(2);
   if(!finite(V[208]) || FError) return 208;
   V[209]=atan(-2*V[202]*(V[152]-V[14]/(V[155]))/(pow(V[50],2)-pow(V[49],2)))/(2);
   if(!finite(V[209]) || FError) return 209;
   V[210]=atan(-2*V[9]*(1-V[206])*(V[153]-V[14]*V[155])/(pow(V[52],2)-pow(V[51],2)))/(2);
   if(!finite(V[210]) || FError) return 210;
   V[211]=cos(V[207]);
   if(!finite(V[211]) || FError) return 211;
   V[212]=sin(V[207]);
   if(!finite(V[212]) || FError) return 212;
   V[213]=-sin(V[207]);
   if(!finite(V[213]) || FError) return 213;
   V[214]=cos(V[207]);
   if(!finite(V[214]) || FError) return 214;
   V[215]=cos(V[208]);
   if(!finite(V[215]) || FError) return 215;
   V[216]=sin(V[208]);
   if(!finite(V[216]) || FError) return 216;
   V[217]=-sin(V[208]);
   if(!finite(V[217]) || FError) return 217;
   V[218]=cos(V[208]);
   if(!finite(V[218]) || FError) return 218;
   V[219]=cos(V[209]);
   if(!finite(V[219]) || FError) return 219;
   V[220]=sin(V[209]);
   if(!finite(V[220]) || FError) return 220;
   V[221]=-sin(V[209]);
   if(!finite(V[221]) || FError) return 221;
   V[222]=cos(V[209]);
   if(!finite(V[222]) || FError) return 222;
   V[223]=cos(V[210]);
   if(!finite(V[223]) || FError) return 223;
   V[224]=sin(V[210]);
   if(!finite(V[224]) || FError) return 224;
   V[225]=-sin(V[210]);
   if(!finite(V[225]) || FError) return 225;
   V[226]=cos(V[210]);
   if(!finite(V[226]) || FError) return 226;
   V[227]=V[99]*pow(V[53],2)*V[99]+V[101]*pow(V[54],2)*V[101];
   if(!finite(V[227]) || FError) return 227;
   V[228]=V[99]*pow(V[53],2)*V[100]+V[101]*pow(V[54],2)*V[102];
   if(!finite(V[228]) || FError) return 228;
   V[229]=V[100]*pow(V[53],2)*V[100]+V[102]*pow(V[54],2)*V[102];
   if(!finite(V[229]) || FError) return 229;
   V[230]=MtRun(sqrt(V[53]*V[54]));
   if(!finite(V[230]) || FError) return 230;
   V[231]=V[228]/(V[230])+V[15]*V[163]/(V[162])*V[165];
   if(!finite(V[231]) || FError) return 231;
   V[232]=V[103]*pow(V[55],2)*V[103]+V[105]*pow(V[56],2)*V[105];
   if(!finite(V[232]) || FError) return 232;
   V[233]=V[103]*pow(V[55],2)*V[104]+V[105]*pow(V[56],2)*V[106];
   if(!finite(V[233]) || FError) return 233;
   V[234]=V[104]*pow(V[55],2)*V[104]+V[106]*pow(V[56],2)*V[106];
   if(!finite(V[234]) || FError) return 234;
   V[235]=MbRun(sqrt(V[55]*V[56]));
   if(!finite(V[235]) || FError) return 235;
   V[236]=V[233]/(V[235])+V[15]*V[162]/(V[163])*V[165];
   if(!finite(V[236]) || FError) return 236;
   V[237]=(V[227]-V[232]+pow(V[235],2)-pow(V[230],2)+2*pow(V[160],2)*pow(V[162],2)-pow(V[160],2))/(pow(V[160],2)*pow(V[162],2));
   if(!finite(V[237]) || FError) return 237;
   V[238]=(pow(V[49],2)-pow(V[51],2)+pow(V[160],2)*(2*pow(V[162],2)*V[1]-2*pow(V[162],2)*pow(V[1],3)+pow(V[1],3)-V[1])/(pow(V[158],2)*V[1]))/(pow(V[160],2)*pow(V[162],2));
   if(!finite(V[238]) || FError) return 238;
   V[239]=(pow(V[45],2)-pow(V[47],2)+pow(V[160],2)*(2*pow(V[162],2)*V[1]-2*pow(V[162],2)*pow(V[1],3)+pow(V[1],3)-V[1])/(pow(V[158],2)*V[1]))/(pow(V[160],2)*pow(V[162],2));
   if(!finite(V[239]) || FError) return 239;
   V[240]=V[66]*V[38]*V[66]+V[71]*V[39]*V[71]+V[76]*V[40]*V[76]+V[81]*V[41]*V[81]+V[86]*V[42]*V[86];
   if(!finite(V[240]) || FError) return 240;
   V[241]=V[66]*V[38]*V[67]+V[71]*V[39]*V[72]+V[76]*V[40]*V[77]+V[81]*V[41]*V[82]+V[86]*V[42]*V[87];
   if(!finite(V[241]) || FError) return 241;
   V[242]=V[66]*V[38]*V[68]+V[71]*V[39]*V[73]+V[76]*V[40]*V[78]+V[81]*V[41]*V[83]+V[86]*V[42]*V[88];
   if(!finite(V[242]) || FError) return 242;
   V[243]=V[66]*V[38]*V[69]+V[71]*V[39]*V[74]+V[76]*V[40]*V[79]+V[81]*V[41]*V[84]+V[86]*V[42]*V[89];
   if(!finite(V[243]) || FError) return 243;
   V[244]=V[66]*V[38]*V[70]+V[71]*V[39]*V[75]+V[76]*V[40]*V[80]+V[81]*V[41]*V[85]+V[86]*V[42]*V[90];
   if(!finite(V[244]) || FError) return 244;
   V[245]=V[67]*V[38]*V[67]+V[72]*V[39]*V[72]+V[77]*V[40]*V[77]+V[82]*V[41]*V[82]+V[87]*V[42]*V[87];
   if(!finite(V[245]) || FError) return 245;
   V[246]=V[67]*V[38]*V[68]+V[72]*V[39]*V[73]+V[77]*V[40]*V[78]+V[82]*V[41]*V[83]+V[87]*V[42]*V[88];
   if(!finite(V[246]) || FError) return 246;
   V[247]=V[67]*V[38]*V[69]+V[72]*V[39]*V[74]+V[77]*V[40]*V[79]+V[82]*V[41]*V[84]+V[87]*V[42]*V[89];
   if(!finite(V[247]) || FError) return 247;
   V[248]=V[67]*V[38]*V[70]+V[72]*V[39]*V[75]+V[77]*V[40]*V[80]+V[82]*V[41]*V[85]+V[87]*V[42]*V[90];
   if(!finite(V[248]) || FError) return 248;
   V[249]=V[68]*V[38]*V[68]+V[73]*V[39]*V[73]+V[78]*V[40]*V[78]+V[83]*V[41]*V[83]+V[88]*V[42]*V[88];
   if(!finite(V[249]) || FError) return 249;
   V[250]=V[68]*V[38]*V[69]+V[73]*V[39]*V[74]+V[78]*V[40]*V[79]+V[83]*V[41]*V[84]+V[88]*V[42]*V[89];
   if(!finite(V[250]) || FError) return 250;
   V[251]=V[68]*V[38]*V[70]+V[73]*V[39]*V[75]+V[78]*V[40]*V[80]+V[83]*V[41]*V[85]+V[88]*V[42]*V[90];
   if(!finite(V[251]) || FError) return 251;
   V[252]=V[69]*V[38]*V[69]+V[74]*V[39]*V[74]+V[79]*V[40]*V[79]+V[84]*V[41]*V[84]+V[89]*V[42]*V[89];
   if(!finite(V[252]) || FError) return 252;
   V[253]=V[69]*V[38]*V[70]+V[74]*V[39]*V[75]+V[79]*V[40]*V[80]+V[84]*V[41]*V[85]+V[89]*V[42]*V[90];
   if(!finite(V[253]) || FError) return 253;
   V[254]=V[70]*V[38]*V[70]+V[75]*V[39]*V[75]+V[80]*V[40]*V[80]+V[85]*V[41]*V[85]+V[90]*V[42]*V[90];
   if(!finite(V[254]) || FError) return 254;
   V[255]=V[240];

   V[256]=V[245];

   V[257]=4*atan(1);
   if(!finite(V[257]) || FError) return 257;
   V[258]=alphaQCD(V[151])/(V[257]);
   if(!finite(V[258]) || FError) return 258;
   V[259]=creal(HggF(pow(V[151]/(2)/(V[202]),2)));
   if(!finite(V[259]) || FError) return 259;
   V[260]=cimag(HggF(pow(V[151]/(2)/(V[202]),2)));
   if(!finite(V[260]) || FError) return 260;
   V[261]=creal(HggA(pow(V[151]/(2)/(V[202]),2)));
   if(!finite(V[261]) || FError) return 261;
   V[262]=cimag(HggA(pow(V[151]/(2)/(V[202]),2)));
   if(!finite(V[262]) || FError) return 262;
   V[263]=creal(HggF(pow(V[151]/(2)/(V[201]),2)));
   if(!finite(V[263]) || FError) return 263;
   V[264]=cimag(HggF(pow(V[151]/(2)/(V[201]),2)));
   if(!finite(V[264]) || FError) return 264;
   V[265]=creal(HggA(pow(V[151]/(2)/(V[201]),2)));
   if(!finite(V[265]) || FError) return 265;
   V[266]=cimag(HggA(pow(V[151]/(2)/(V[201]),2)));
   if(!finite(V[266]) || FError) return 266;
   V[267]=creal(HggF(pow(V[151]/(2)/(V[200]),2)));
   if(!finite(V[267]) || FError) return 267;
   V[268]=cimag(HggF(pow(V[151]/(2)/(V[200]),2)));
   if(!finite(V[268]) || FError) return 268;
   V[269]=creal(HggA(pow(V[151]/(2)/(V[200]),2)));
   if(!finite(V[269]) || FError) return 269;
   V[270]=cimag(HggA(pow(V[151]/(2)/(V[200]),2)));
   if(!finite(V[270]) || FError) return 270;
   V[271]=creal(HggF(pow(V[151]/(2)/(V[7]),2)));
   if(!finite(V[271]) || FError) return 271;
   V[272]=cimag(HggF(pow(V[151]/(2)/(V[7]),2)));
   if(!finite(V[272]) || FError) return 272;
   V[273]=creal(HggA(pow(V[151]/(2)/(V[7]),2)));
   if(!finite(V[273]) || FError) return 273;
   V[274]=cimag(HggA(pow(V[151]/(2)/(V[7]),2)));
   if(!finite(V[274]) || FError) return 274;
   V[275]=creal(HggS(pow(V[151]/(2)/(V[57]),2)));
   if(!finite(V[275]) || FError) return 275;
   V[276]=cimag(HggS(pow(V[151]/(2)/(V[57]),2)));
   if(!finite(V[276]) || FError) return 276;
   V[277]=creal(HggS(pow(V[151]/(2)/(V[58]),2)));
   if(!finite(V[277]) || FError) return 277;
   V[278]=cimag(HggS(pow(V[151]/(2)/(V[58]),2)));
   if(!finite(V[278]) || FError) return 278;
   V[279]=creal(HggS(pow(V[151]/(2)/(V[59]),2)));
   if(!finite(V[279]) || FError) return 279;
   V[280]=cimag(HggS(pow(V[151]/(2)/(V[59]),2)));
   if(!finite(V[280]) || FError) return 280;
   V[281]=creal(HggS(pow(V[151]/(2)/(V[60]),2)));
   if(!finite(V[281]) || FError) return 281;
   V[282]=cimag(HggS(pow(V[151]/(2)/(V[60]),2)));
   if(!finite(V[282]) || FError) return 282;
   V[283]=creal(HggS(pow(V[151]/(2)/(V[61]),2)));
   if(!finite(V[283]) || FError) return 283;
   V[284]=cimag(HggS(pow(V[151]/(2)/(V[61]),2)));
   if(!finite(V[284]) || FError) return 284;
   V[285]=creal(HggS(pow(V[151]/(2)/(V[62]),2)));
   if(!finite(V[285]) || FError) return 285;
   V[286]=cimag(HggS(pow(V[151]/(2)/(V[62]),2)));
   if(!finite(V[286]) || FError) return 286;
   V[287]=creal(HggS(pow(V[151]/(2)/(V[45]),2)));
   if(!finite(V[287]) || FError) return 287;
   V[288]=cimag(HggS(pow(V[151]/(2)/(V[45]),2)));
   if(!finite(V[288]) || FError) return 288;
   V[289]=creal(HggS(pow(V[151]/(2)/(V[46]),2)));
   if(!finite(V[289]) || FError) return 289;
   V[290]=cimag(HggS(pow(V[151]/(2)/(V[46]),2)));
   if(!finite(V[290]) || FError) return 290;
   V[291]=creal(HggS(pow(V[151]/(2)/(V[49]),2)));
   if(!finite(V[291]) || FError) return 291;
   V[292]=cimag(HggS(pow(V[151]/(2)/(V[49]),2)));
   if(!finite(V[292]) || FError) return 292;
   V[293]=creal(HggS(pow(V[151]/(2)/(V[50]),2)));
   if(!finite(V[293]) || FError) return 293;
   V[294]=cimag(HggS(pow(V[151]/(2)/(V[50]),2)));
   if(!finite(V[294]) || FError) return 294;
   V[295]=creal(HggS(pow(V[151]/(2)/(V[53]),2)));
   if(!finite(V[295]) || FError) return 295;
   V[296]=cimag(HggS(pow(V[151]/(2)/(V[53]),2)));
   if(!finite(V[296]) || FError) return 296;
   V[297]=creal(HggS(pow(V[151]/(2)/(V[54]),2)));
   if(!finite(V[297]) || FError) return 297;
   V[298]=cimag(HggS(pow(V[151]/(2)/(V[54]),2)));
   if(!finite(V[298]) || FError) return 298;
   V[299]=creal(HggS(pow(V[151]/(2)/(V[47]),2)));
   if(!finite(V[299]) || FError) return 299;
   V[300]=cimag(HggS(pow(V[151]/(2)/(V[47]),2)));
   if(!finite(V[300]) || FError) return 300;
   V[301]=creal(HggS(pow(V[151]/(2)/(V[48]),2)));
   if(!finite(V[301]) || FError) return 301;
   V[302]=cimag(HggS(pow(V[151]/(2)/(V[48]),2)));
   if(!finite(V[302]) || FError) return 302;
   V[303]=creal(HggS(pow(V[151]/(2)/(V[51]),2)));
   if(!finite(V[303]) || FError) return 303;
   V[304]=cimag(HggS(pow(V[151]/(2)/(V[51]),2)));
   if(!finite(V[304]) || FError) return 304;
   V[305]=creal(HggS(pow(V[151]/(2)/(V[52]),2)));
   if(!finite(V[305]) || FError) return 305;
   V[306]=cimag(HggS(pow(V[151]/(2)/(V[52]),2)));
   if(!finite(V[306]) || FError) return 306;
   V[307]=creal(HggS(pow(V[151]/(2)/(V[55]),2)));
   if(!finite(V[307]) || FError) return 307;
   V[308]=cimag(HggS(pow(V[151]/(2)/(V[55]),2)));
   if(!finite(V[308]) || FError) return 308;
   V[309]=creal(HggS(pow(V[151]/(2)/(V[56]),2)));
   if(!finite(V[309]) || FError) return 309;
   V[310]=cimag(HggS(pow(V[151]/(2)/(V[56]),2)));
   if(!finite(V[310]) || FError) return 310;
   V[311]=creal(HggS(pow(V[151]/(2)/(V[37]),2)));
   if(!finite(V[311]) || FError) return 311;
   V[312]=cimag(HggS(pow(V[151]/(2)/(V[37]),2)));
   if(!finite(V[312]) || FError) return 312;
   V[313]=creal(HggV(pow(V[151]/(2)/(V[160]),2)));
   if(!finite(V[313]) || FError) return 313;
   V[314]=cimag(HggV(pow(V[151]/(2)/(V[160]),2)));
   if(!finite(V[314]) || FError) return 314;
   V[315]=creal(HggF(pow(V[151]/(2)/(V[43]),2)));
   if(!finite(V[315]) || FError) return 315;
   V[316]=cimag(HggF(pow(V[151]/(2)/(V[43]),2)));
   if(!finite(V[316]) || FError) return 316;
   V[317]=creal(HggA(pow(V[151]/(2)/(V[43]),2)));
   if(!finite(V[317]) || FError) return 317;
   V[318]=cimag(HggA(pow(V[151]/(2)/(V[43]),2)));
   if(!finite(V[318]) || FError) return 318;
   V[319]=creal(HggF(pow(V[151]/(2)/(V[44]),2)));
   if(!finite(V[319]) || FError) return 319;
   V[320]=cimag(HggF(pow(V[151]/(2)/(V[44]),2)));
   if(!finite(V[320]) || FError) return 320;
   V[321]=creal(HggA(pow(V[151]/(2)/(V[44]),2)));
   if(!finite(V[321]) || FError) return 321;
   V[322]=cimag(HggA(pow(V[151]/(2)/(V[44]),2)));
   if(!finite(V[322]) || FError) return 322;
   V[323]=V[202]*McRun(V[151]/(2))/(McRun(V[202]));
   if(!finite(V[323]) || FError) return 323;
   V[324]=V[201]*MbRun(V[151]/(2))/(MbRun(V[201]));
   if(!finite(V[324]) || FError) return 324;
   V[325]=V[6]*MtRun(V[151]/(2))/(MtRun(V[200]));
   if(!finite(V[325]) || FError) return 325;
   V[326]=creal(HggF(pow(V[151]/(2)/(V[323]),2))*(1+V[258]*Hgam1F(pow(V[151]/(2)/(V[323]),2))));
   if(!finite(V[326]) || FError) return 326;
   V[327]=cimag(HggF(pow(V[151]/(2)/(V[323]),2))*(1+V[258]*Hgam1F(pow(V[151]/(2)/(V[323]),2))));
   if(!finite(V[327]) || FError) return 327;
   V[328]=creal(HggA(pow(V[151]/(2)/(V[323]),2))*(1+V[258]*Hgam1A(pow(V[151]/(2)/(V[323]),2))));
   if(!finite(V[328]) || FError) return 328;
   V[329]=cimag(HggA(pow(V[151]/(2)/(V[323]),2))*(1+V[258]*Hgam1A(pow(V[151]/(2)/(V[323]),2))));
   if(!finite(V[329]) || FError) return 329;
   V[330]=creal(HggF(pow(V[151]/(2)/(V[324]),2))*(1+V[258]*Hgam1F(pow(V[151]/(2)/(V[324]),2))));
   if(!finite(V[330]) || FError) return 330;
   V[331]=cimag(HggF(pow(V[151]/(2)/(V[324]),2))*(1+V[258]*Hgam1F(pow(V[151]/(2)/(V[324]),2))));
   if(!finite(V[331]) || FError) return 331;
   V[332]=creal(HggA(pow(V[151]/(2)/(V[324]),2))*(1+V[258]*Hgam1A(pow(V[151]/(2)/(V[324]),2))));
   if(!finite(V[332]) || FError) return 332;
   V[333]=cimag(HggA(pow(V[151]/(2)/(V[324]),2))*(1+V[258]*Hgam1A(pow(V[151]/(2)/(V[324]),2))));
   if(!finite(V[333]) || FError) return 333;
   V[334]=creal(HggF(pow(V[151]/(2)/(V[325]),2))*(1+V[258]*Hgam1F(pow(V[151]/(2)/(V[325]),2))));
   if(!finite(V[334]) || FError) return 334;
   V[335]=cimag(HggF(pow(V[151]/(2)/(V[325]),2))*(1+V[258]*Hgam1F(pow(V[151]/(2)/(V[325]),2))));
   if(!finite(V[335]) || FError) return 335;
   V[336]=creal(HggA(pow(V[151]/(2)/(V[325]),2))*(1+V[258]*Hgam1A(pow(V[151]/(2)/(V[325]),2))));
   if(!finite(V[336]) || FError) return 336;
   V[337]=cimag(HggA(pow(V[151]/(2)/(V[325]),2))*(1+V[258]*Hgam1A(pow(V[151]/(2)/(V[325]),2))));
   if(!finite(V[337]) || FError) return 337;
   V[338]=pow(V[151]/(2*V[45]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[45])),6/(double)((23)))),2);
   if(!finite(V[338]) || FError) return 338;
   V[339]=creal(HggS(V[338])*(1+V[258]*Hgam1S(V[338])));
   if(!finite(V[339]) || FError) return 339;
   V[340]=cimag(HggS(V[338])*(1+V[258]*Hgam1S(V[338])));
   if(!finite(V[340]) || FError) return 340;
   V[341]=pow(V[151]/(2*V[46]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[46])),6/(double)((23)))),2);
   if(!finite(V[341]) || FError) return 341;
   V[342]=creal(HggS(V[341])*(1+V[258]*Hgam1S(V[341])));
   if(!finite(V[342]) || FError) return 342;
   V[343]=cimag(HggS(V[341])*(1+V[258]*Hgam1S(V[341])));
   if(!finite(V[343]) || FError) return 343;
   V[344]=pow(V[151]/(2*V[49]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[49])),6/(double)((23)))),2);
   if(!finite(V[344]) || FError) return 344;
   V[345]=creal(HggS(V[344])*(1+V[258]*Hgam1S(V[344])));
   if(!finite(V[345]) || FError) return 345;
   V[346]=cimag(HggS(V[344])*(1+V[258]*Hgam1S(V[344])));
   if(!finite(V[346]) || FError) return 346;
   V[347]=pow(V[151]/(2*V[50]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[50])),6/(double)((23)))),2);
   if(!finite(V[347]) || FError) return 347;
   V[348]=creal(HggS(V[347])*(1+V[258]*Hgam1S(V[347])));
   if(!finite(V[348]) || FError) return 348;
   V[349]=cimag(HggS(V[347])*(1+V[258]*Hgam1S(V[347])));
   if(!finite(V[349]) || FError) return 349;
   V[350]=pow(V[151]/(2*V[53]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[53])),6/(double)((23)))),2);
   if(!finite(V[350]) || FError) return 350;
   V[351]=creal(HggS(V[350])*(1+V[258]*Hgam1S(V[350])));
   if(!finite(V[351]) || FError) return 351;
   V[352]=cimag(HggS(V[350])*(1+V[258]*Hgam1S(V[350])));
   if(!finite(V[352]) || FError) return 352;
   V[353]=pow(V[151]/(2*V[54]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[54])),6/(double)((23)))),2);
   if(!finite(V[353]) || FError) return 353;
   V[354]=creal(HggS(V[353])*(1+V[258]*Hgam1S(V[353])));
   if(!finite(V[354]) || FError) return 354;
   V[355]=cimag(HggS(V[353])*(1+V[258]*Hgam1S(V[353])));
   if(!finite(V[355]) || FError) return 355;
   V[356]=pow(V[151]/(2*V[47]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[47])),6/(double)((23)))),2);
   if(!finite(V[356]) || FError) return 356;
   V[357]=creal(HggS(V[356])*(1+V[258]*Hgam1S(V[356])));
   if(!finite(V[357]) || FError) return 357;
   V[358]=cimag(HggS(V[356])*(1+V[258]*Hgam1S(V[356])));
   if(!finite(V[358]) || FError) return 358;
   V[359]=pow(V[151]/(2*V[48]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[48])),6/(double)((23)))),2);
   if(!finite(V[359]) || FError) return 359;
   V[360]=creal(HggS(V[359])*(1+V[258]*Hgam1S(V[359])));
   if(!finite(V[360]) || FError) return 360;
   V[361]=cimag(HggS(V[359])*(1+V[258]*Hgam1S(V[359])));
   if(!finite(V[361]) || FError) return 361;
   V[362]=pow(V[151]/(2*V[51]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[51])),6/(double)((23)))),2);
   if(!finite(V[362]) || FError) return 362;
   V[363]=creal(HggS(V[362])*(1+V[258]*Hgam1S(V[362])));
   if(!finite(V[363]) || FError) return 363;
   V[364]=cimag(HggS(V[362])*(1+V[258]*Hgam1S(V[362])));
   if(!finite(V[364]) || FError) return 364;
   V[365]=pow(V[151]/(2*V[52]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[52])),6/(double)((23)))),2);
   if(!finite(V[365]) || FError) return 365;
   V[366]=creal(HggS(V[365])*(1+V[258]*Hgam1S(V[365])));
   if(!finite(V[366]) || FError) return 366;
   V[367]=cimag(HggS(V[365])*(1+V[258]*Hgam1S(V[365])));
   if(!finite(V[367]) || FError) return 367;
   V[368]=pow(V[151]/(2*V[55]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[55])),6/(double)((23)))),2);
   if(!finite(V[368]) || FError) return 368;
   V[369]=creal(HggS(V[368])*(1+V[258]*Hgam1S(V[368])));
   if(!finite(V[369]) || FError) return 369;
   V[370]=cimag(HggS(V[368])*(1+V[258]*Hgam1S(V[368])));
   if(!finite(V[370]) || FError) return 370;
   V[371]=pow(V[151]/(2*V[56]*pow(alphaQCD(V[151]/(2))/(alphaQCD(V[56])),6/(double)((23)))),2);
   if(!finite(V[371]) || FError) return 371;
   V[372]=creal(HggS(V[371])*(1+V[258]*Hgam1S(V[371])));
   if(!finite(V[372]) || FError) return 372;
   V[373]=cimag(HggS(V[371])*(1+V[258]*Hgam1S(V[371])));
   if(!finite(V[373]) || FError) return 373;
   if(VV==NULL) 
   {  VV=malloc(sizeof(REAL)*nModelVars);
      for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0) iQ=i;
   }
   for(i=0;i<nModelVars;i++) VV[i]=V[i];
   cErr=0;
   return 0;
}
