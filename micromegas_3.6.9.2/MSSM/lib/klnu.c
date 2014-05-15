
 /* New physics contribution to  Ratio R_l23 in MSSM with main tanbeta NLO contributions
 This is the ratio of 2body/3body Kaon decay as defined in M. Antonelli et al (FlaviaNet Working group) arXiv:0801.1817
*/

#include "../../sources/micromegas.h"
#include "pmodel.h"
#include "pmodel_aux.h"

#define mK  0.494
#define  mdms  1./20.2 /* md/ms from PDG PRD86 010001 (2012)*/
#define  mDs  1.968  /* from PDG PRD86 010001 (2012)*/
#define msc  0.08  /* mc/ms */
#define TauDs  7.596e+11 /*  in GeV-1 from PDG : 500+/-7 10^-15sec */
#define fDs  0.2486  /*  +/- 0.003 from lattice, 0910.2928 */
#define Vcs  0.98 /* +/- 0.01 +/- 0.1 PDG from B-Klnu*/
#define GF  1.166e-5

double Rl23_(void)
{
	double rl23,mHp,xHp;
	struct read_param_tag param;
	double eps_d=0;
	
	if(read_prm(&param))
	{
		puts("klnu: can not read parameters.");
		return 0.0;
	}
	
/*	dump_prm(&param);
*/
       
     mHp=param.Mhc;
	 xHp=mK*mK/mHp/mHp;
     rl23=(1-xHp*(1-mdms)*param.tb*param.tb/(1+deltaMd()));
      	
	return rl23;
}

double rl23_(void) { return Rl23_();}

double dtaunu_(double * dmunu)
{
	double taunu,lnu,lmu,mHp,xHp;
	struct read_param_tag param;
    double Mm=0.1057;
		
	if(read_prm(&param))
	{
		puts("dtaunu: can not read parameters.");
		return 0.0;
	}
	
    mHp=param.Mhc;
    xHp=mDs*mDs/mHp/mHp;
    lnu=GF*GF/8./M_PI*Vcs*Vcs*fDs*fDs*TauDs*mDs
    *(1+xHp/(1+msc)*(1-msc*param.tb*param.tb/(1+deltaMc())))
    *(1+xHp/(1+msc)*(1-msc*param.tb*param.tb/(1+deltaMc())));
    
    lmu=lnu*Mm*Mm*(1-Mm*Mm/mDs/mDs)*(1-Mm*Mm/mDs/mDs);
    lnu*=(1-findValW("Ml")*findValW("Ml")/mDs/mDs)*(1-findValW("Ml")*findValW("Ml")/mDs/mDs);
    
  
 /*   printf("constant =%.3e tau =%.3e %.3e rat1%.3e rat2%.3e \n",GF*GF/8./M_PI*Vcs*Vcs*fDs*fDs*TauDs*mDs,TauDs,findValW("Ml"),
          (1-findValW("Ml")*findValW("Ml")/mDs/mDs), 
           (1+xHp/(1+msc)*(1-msc*param.tb*param.tb/(1+deltaMc()))));*/
    *dmunu=lmu;
	return lnu*findValW("Ml")*findValW("Ml");
}
