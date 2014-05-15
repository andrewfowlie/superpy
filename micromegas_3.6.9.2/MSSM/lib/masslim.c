#include<math.h>
#include"../../sources/micromegas.h"
#include "pmodel.h"


int masslimits_(void)
{
  double tb,c2b,mc1,msne,msnm,msnl,mstau1,mn1,msel,mser,mse1,msmul,msmur,msmu1,Mh,Mhh;
  double mst1,msb1,msul,msur,msdl,msdr;
  int retcode=0;
  int limCharg=0,limSne=0,limSnm=0,limSnl=0,limSeR=0,limSmR=0 ,limSl1=0,limh=0,limSb1=0,limSt1=0,limSq=0;

  if(findVal("MNE1",&mn1))  return -1; mn1=fabs(mn1);  
  if(findVal("MC1",&mc1))   return -1; mc1=fabs(mc1);
  if(findVal("tB",&tb))     return -1; c2b=(1-tb*tb)/(1+tb*tb);
  if(findVal("MSeL",&msel)) return -1; msel=fabs(msel);
  if(findVal("MSeR",&mser)) return -1; mser=fabs(mser); 
 /*      mse1=(msel<mser)? msel:mser; */
mse1=mser;
  if(findVal("MSne",&msne)) return -1; msne=fabs(msne);
  if(findVal("MSnm",&msnm)) return -1; msnm=fabs(msnm); 
  if(findVal("MSnl",&msnl)) return -1; msnl=fabs(msnl);
  if(findVal("MSmL",&msmul))return -1; msmul=fabs(msmul);
  if(findVal("MSmR",&msmur))return -1; msmur=fabs(msmur); 
/*      msmu1=(msmul<msmur)? msmul:msmur;*/
  msmu1=msmur;
  if(findVal("MSl1",&mstau1))return -1;mstau1=fabs(mstau1);
  if(findVal("MSt1",&mst1)) return -1; mst1=fabs(mst1); 
  if(findVal("MSb1",&msb1)) return -1; msb1=fabs(msb1); 
  if(findVal("MSuL",&msul)) return -1; msul=fabs(msul); 
  if(findVal("MSuR",&msur)) return -1; msur=fabs(msur); 
  if(findVal("MSdL",&msdl)) return -1; msdl=fabs(msdl); 
  if(findVal("MSdR",&msdr)) return -1; msdr=fabs(msdr); 

  Mh=findValW("Mh");
  Mhh=findValW("MHH");
  if(Mh>128) limh=1;
  if(Mh<123&&Mhh>128) limh=1;
  if(Mh<123&&Mhh<123) limh=1;

  if(mc1-mn1>60.)
  {if((msne>43. &&msne<55. &&mc1<89.5)   ||
       (msne>55. &&msne<65. &&mc1<80.)   ||
       (msne>65. &&msne<75. &&mc1<74.4)  ||
       (msne>75. &&msne<85. &&mc1<73.1)  ||
       (msne>85. &&msne<95. &&mc1<73.6)  ||
       (msne>95. &&msne<105.&&mc1<75.1)  ||
       (msne>105.&&msne<115.&&mc1<76.9)  ||
       (msne>115.&&msne<125.&&mc1<80.1)  ||
       (msne>125.&&msne<135.&&mc1<83.14) ||
       (msne>135.&&msne<145.&&mc1<85.37) ||
       (msne>145.&&msne<155.&&mc1<88.4)  ||
       (msne>155.&&msne<165.&&mc1<91.04) ||
       (msne>165.&&msne<175.&&mc1<92.4)  ||
       (msne>175.&&msne<185.&&mc1<94.7)  ||
       (msne>185.&&msne<195.&&mc1<95.5)  ||
       (msne>195.&&msne<205.&&mc1<96.4)  ||
       (msne>205.&&msne<275.&&mc1<99.4)  ||
       (msne>275.&&msne<325.&&mc1<100.5) ||
       (msne>325.&&msne<375.&&mc1<101)   ||
       (msne>375.&&msne<425.&&mc1<101.2) ||
       (msne>425.&&mc1<101.5)
      )  limCharg=1;
  }
  
  if(mc1-mn1<60.0)
   {/*if(msne<45.) return 2;*/
    if((msne>43. &&msne<55. &&mc1<101.)  ||
       (msne>55. &&msne<65. &&mc1<98.)   ||
       (msne>65. &&msne<75. &&mc1<95.)   ||
       (msne>75. &&msne<85. &&mc1<91.)   ||
       (msne>85. &&msne<95. &&mc1<90.)   ||
       (msne>95. &&msne<105.&&mc1<92.)   ||
       (msne>105.&&msne<115.&&mc1<93.)   ||
       (msne>115.&&msne<125.&&mc1<96.)   ||
       (msne>125.&&msne<135.&&mc1<97.)   ||
       (msne>135.&&msne<145.&&mc1<99.)   ||
       (msne>145.&&msne<155.&&mc1<100.)  ||
       (msne>155.&&msne<165.&&mc1<101.)  ||
       (msne>165.&&msne<175.&&mc1<101.)  ||
       (msne>175.&&msne<225.&&mc1<101.)  ||
       (msne>225.&&msne<275.&&mc1<102.)  ||
       (msne>275.&&msne<325.&&mc1<102.5) ||
       (msne>325.&&msne<375.&&mc1<102.9) ||
       (msne>375.&&mc1<103.) 
      ) limCharg=1; 
  }
	 
	 		
  if(msne<43.0)           limSne=1;
  if(msnm<43.0)           limSnm=1;
  if(msnl<43.0)           limSnl=1; 
  if(mse1<99.4 &&mn1>40.) limSeR=1;
  if(mse1<100.5&&mn1<40.) limSeR=1;
  if(msmu1<95.0)          limSmR=1;
  if(msul<40.0)           limSq=1;
  if(msur<40.0)           limSq=1;
  if(msdl<40.0)           limSq=1;
  if(msdr<40.0)           limSq=1;
  if(mst1<63.0)           limSt1=1;
  if(msb1<89 && msb1-mn1>8.0)   limSb1=1;
  
 

/* These limits are valid only if Mstau-mn1> 3GeV*/	 

  if((mstau1 < 80.5&& mn1<10.)           ||
     (mstau1 < 82. && mn1>10.&& mn1<20.) ||
     (mstau1 < 84. && mn1>20.&& mn1<30.) ||
     (mstau1 < 85. && mn1>30.&& mn1<40.) ||
     (mstau1 < 87.1&& mn1>40.&& mn1<50.) ||
     (mstau1 < 88. && mn1>50.&& mn1<60.) ||
     (mstau1 < 85. && mn1>60.&& mn1<75.) 
    ) limSl1=1;	


  if(limCharg){retcode+=1;   printf("WARNING: Chargino below LEP limit \n");}
  if(limSne)  {retcode+=2;   printf("WARNING: Sneutrino-e below LEP limit \n");}
  if(limSnm)  {retcode+=4;   printf("WARNING: Sneutrino-mu below LEP limit \n");}
  if(limSnl)  {retcode+=8;   printf("WARNING: Sneutrino-tau below LEP limit \n");} 
  if(limSeR)  {retcode+=16;  printf("WARNING: Selectron-R below LEP limit \n");} 
  if(limSmR)  {retcode+=32;  printf("WARNING: Smuon-R below LEP limit \n");}
  if(limSl1)  {retcode+=64;  printf("WARNING: Stau1 below LEP limit \n");}
   if(limh)   {retcode+=128; printf("WARNING: there is no Higgs particle with a mass between 123-128GeV\n");}
  if(limSt1)  {retcode+=256; printf("WARNING: Stop1 below LEP limit \n");}
  if(limSb1)  {retcode+=512; printf("WARNING: Sbot1 below LEP limit \n");}
  if(limSq)   {retcode+=1024;printf("WARNING: Squarks below LEP limit \n");}



  return retcode;
}
