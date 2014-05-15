#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include<string.h>
#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include "pmodel.h"
#include "localpath.h"


double cpHiggs(double EE,  double alfaSMZ,double Mtp,   double McMt, double MbMt, double Ml,   double SW, double tb,  double MHc,
               double aMu, double fiMu,   double aM1,   double fiM1, double aM2,  double fiM2,
               double  aM3,double  fiM3,   double Mq3,  double Mu3,  double Md3,  double Ml3,  double  Mr3, 
               double aAt, double fiAt,    double aAb,  double fiAb, double  aAl, double fiAl,
               double aAu, double fiAu,    double aAd,  double fiAd, double  aAe, double fiAe, double  aAm, double fiAm)
{ FILE *f;
  char * command;
  int err;  
 f=fopen("CPsuperH.in","w");  
 
 fprintf(f,"%E              ! SMPARA( 1) = 1/AEM(MZ)  \n",4*M_PI/EE/EE );
 fprintf(f,"%E             ! SMPARA( 2) = AS(MZ)  \n",alfaSMZ);
 fprintf(f,"91.187D0            ! SMPARA( 3) = MZ in GeV \n");
 fprintf(f,"%E             ! SMPARA( 4) = sin^2 Theta_W \n",SW*SW);
 fprintf(f,"0.5D-3              ! SMPARA( 5) = m_e in GeV \n");
 fprintf(f,"0.1065D0            ! SMPARA( 6) = m_mu in GeV \n");
 fprintf(f,"%E             ! SMPARA( 7) = m_tau in GeV \n",Ml);
 fprintf(f,"6.D-3               ! SMPARA( 8) = m_d (m_t) in GeV \n");
 fprintf(f,"0.115D0             ! SMPARA( 9) = m_s (m_t) in GeV \n");
 fprintf(f,"%E             ! SMPARA(10) = m_b (m_t) in GeV \n",MbMt);
 fprintf(f,"3.D-3               ! SMPARA(11) = m_u (m_t) in GeV \n");
 fprintf(f,"%E              ! SMPARA(12) = m_c (m_t) in GeV \n",McMt);
 fprintf(f,"%E             ! SMPARA(13) = m_t^POLE  in GeV \n",Mtp);
 fprintf(f,"2.118D0             ! SMPARA(14) = Gam_W  in GeV \n" );
 fprintf(f,"2.4952D0            ! SMPARA(15) = Gam_Z  in GeV \n");
/*new*/ 
 fprintf(f,"0.2272D0            ! SMPARA(16) = lambda_CKM\n");
 fprintf(f,"0.8180D0            ! SMPARA(17) = A_CKM\n");
 fprintf(f,"0.2210D0            ! SMPARA(18) = rho^bar_CKM\n");
 fprintf(f,"0.3400D0            ! SMPARA(19) = eta^bar_CKM\n");
  fprintf(f,"125.                ! SMPARA(20) = M_HSM\n");
/* end of new*/

 fprintf(f,"%E      ! SSPARA( 1) = tan beta \n",             tb      );
 fprintf(f,"%E      ! SSPARA( 2) = m_H^pm^POLE in GeV \n",   MHc     );
 fprintf(f,"%E      ! SSPARA( 3) = |mu| in GeV \n",          aMu     );
 fprintf(f,"%E      ! SSPARA( 4) = Phi_mu in Degree \n",     fiMu    );
 fprintf(f,"%E      ! SSPARA( 5) = |M_1| in GeV \n",         aM1     );
 fprintf(f,"%E      ! SSPARA( 6) = Phi_1 in Degree \n",      fiM1    );
 fprintf(f,"%E      ! SSPARA( 7) = |M_2| in GeV \n",         aM2     );
 fprintf(f,"%E      ! SSPARA( 8) = Phi_2 in Degree \n",      fiM2    );
 fprintf(f,"%E      ! SSPARA( 9) = |M_3| in GeV \n",         aM3     );
 fprintf(f,"%E      ! SSPARA(10) = Phi_3 in Degree \n",      fiM3    );
 fprintf(f,"%E      ! SSPARA(11) = m_Q3 in GeV \n",          Mq3     );
 fprintf(f,"%E      ! SSPARA(12) = m_U3 in GeV \n",          Mu3     );
 fprintf(f,"%E      ! SSPARA(13) = m_D3 in GeV \n",          Md3     );
 fprintf(f,"%E      ! SSPARA(14) = m_L3 in GeV \n",          Ml3     );
 fprintf(f,"%E      ! SSPARA(15) = m_E3 in GeV \n",          Mr3     );
 fprintf(f,"%E      ! SSPARA(16) = |A_t| in GeV \n",         aAt     );
 fprintf(f,"%E      ! SSPARA(17) = Phi_{A_t} in Degree \n",  fiAt    );
 fprintf(f,"%E      ! SSPARA(18) = |A_b| in GeV \n",         aAb     );
 fprintf(f,"%E      ! SSPARA(19) = Phi_{A_b} in Degree \n",  fiAb    );
 fprintf(f,"%E      ! SSPARA(20) = |A_tau| in GeV \n",       aAl     );
 fprintf(f,"%E      ! SSPARA(21) = Phi_{A_tau} in Degree \n",fiAl    );
/* new */
fprintf(f,"%E              ! SSPARA(22) = Hierarchy factor between first 2 and third generations M_Q\n",findValW("Mq2")/Mq3);
fprintf(f,"%E              ! SSPARA(23) = Hierarchy factor between first 2 and third generations M_U\n",findValW("Mu2")/Mu3);
fprintf(f,"%E              ! SSPARA(24) = Hierarchy factor between first 2 and third generations M_D\n",findValW("Md2")/Md3);
fprintf(f,"%E              ! SSPARA(25) = Hierarchy factor between first 2 and third generations M_L\n",findValW("Ml2")/Ml3);
fprintf(f,"%E              ! SSPARA(26) = Hierarchy factor between first 2 and third generations M_E\n",findValW("Mr2")/Mr3);
/*===*/

fprintf(f,"%E               ! SSPARA(27) = |A_e| in GeV\n",aAe);
fprintf(f,"%E               ! SSPARA(28) = Phi_{A_e} in Degree\n",fiAe);
fprintf(f,"%E               ! SSPARA(29) = |A_mu| in GeV\n",aAm);
fprintf(f,"%E               ! SSPARA(30) = Phi_{A_mu} in Degree\n",fiAm);
fprintf(f,"%E               ! SSPARA(31) = |A_u| in GeV\n",aAu);
fprintf(f,"0.               ! SSPARA(32) = Phi_{A_u} in Degree\n");
fprintf(f,"%E               ! SSPARA(33) = |A_c| in GeV\n",aAu);
fprintf(f,"0.               ! SSPARA(34) = Phi_{A_c} in Degree\n",fiAu);
fprintf(f,"%E               ! SSPARA(35) = |A_d| in GeV\n",aAd);
fprintf(f,"0.               ! SSPARA(36) = Phi_{A_d} in Degree\n");
fprintf(f,"%E               ! SSPARA(37) = |A_s| in GeV\n",aAd);
fprintf(f,"%E               ! SSPARA(38) = Phi_{A_s} in Degree\n",fiAd);


 fprintf(f,"0                   ! IFLAG_H(1)  if 1, print input parameters \n");
 fprintf(f,"1                   ! IFLAG_H(2)  if 1, print Higgs sector \n");
 fprintf(f,"1                   ! IFLAG_H(3)  if 1, print sfermion sector \n");
 fprintf(f,"1                   ! IFLAG_H(4)  if 1, print -ino sector \n");
 fprintf(f,"0                   ! IFLAG_H(5)  print couplings : 1=H_1 2=H_2 3=H_3 4=H^pm 5=Higgs self 6=ALL \n");
 fprintf(f,"5                   ! IFLAG_H(6)  print decay widths and brs : 1=H_1 2=H_2 3=H_3 4=H^pm 5=ALL \n");
 fprintf(f,"0                   ! IFLAG_H(10) if 0, include the threshold corrections  \n");
 fprintf(f,"0                   ! IFLAG_H(11) use pole mass (0) or eff. pot. mass (1) \n");
/* new */


 fprintf(f,"5                   ! IFLAG_H(12) 5 or 0 for full improvement\n");
 fprintf(f,"0                   ! IFLAG_H(13) 1 Not to include the off-diagonal absorptive parts\n");
 fprintf(f,"1                   ! IFLAG_H(14) 1 to print FILLDHPG results\n");
/* fprintf(f,"1                   ! IFLAG_H(15) 1 to print HiggsEDM results\n");*/
 fprintf(f,"1                   ! IFLAG_H(16) 1 to print FILLBOBS results\n");
 fprintf(f,"1                   ! IFLAG_H(17) 1 to print b -> s gamma details\n");
 fprintf(f,"1                   ! IFLAG_H(18) 1 to print EDM results\n");
 fprintf(f,"1                   ! IFLAG_H(19) 1 to print fllmuon results\n");
 fprintf(f,"1                   ! IFLAG_H(20) 1 to print SLHA2\n");
/*===*/
 fclose(f);

 command=malloc(strlen(CPsuperH)+100);
 sprintf(command,"%s < CPsuperH.in > CPsuperH.out",CPsuperH);
 err=System(command);
 free(command);

 err=slhaRead("cpsuperh2_slha.out",4);
 return err;
}
