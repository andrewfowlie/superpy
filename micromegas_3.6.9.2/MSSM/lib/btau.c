/* New physics contribution to  B->tau nu in MSSM with main
 tanbeta NLO contributions
*/

#include "../../sources/micromegas.h"
#include "pmodel.h"
#include "pmodel_aux.h"


				

/*double btaunu_(double *deltam)
*/
double btaunu_(void)
{
	double rtau,mBu,mHp,xHp;
	struct read_param_tag param;
	double eps_b=0, eps_bp=0, eps_tps=0;
	
	if(read_prm(&param))
	{
		puts("btaunu: can not read parameters.");
		return 0.0;
	}
	
/*	dump_prm(&param);
*/
       
         mHp=param.Mhc;
	 mBu=5.279;
	 xHp=mBu*mBu/mHp/mHp;
	     

/*Large tan beta */
/* Calculates the eps_b eps_b' ..*/	
	calc_eps(&param, &eps_b, &eps_bp, &eps_tps);
		
         rtau=(1-xHp*param.tb*param.tb/(1+eps_b*param.tb))*(1-xHp*param.tb*param.tb/(1+eps_b*param.tb));
      	
	return rtau;
}



