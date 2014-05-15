#include "src/include.h"

/*--------------------------------------------------------*/
/* Calculation of the observables in the Standard Model   */
/*--------------------------------------------------------*/

int main()
{ 
	struct parameters param;
		
	Init_param(&param);
	slha_adjust(&param);
	param.SM=1;
	
	printf("\n");
	
	printf("SuperIso v3.3 - F. Mahmoudi\n\n");
	printf("Standard Model predictions\n\n");
	printf("Observable\t\t\tValue\n\n");

	double C0b[11],C0spec[11],C1b[11],C1spec[11],C0w[11],C1w[11],C2w[11],C2b[11],Cpb[11];
	double complex CQpb[3],CQ0b[3],CQ1b[3];
	CQ0b[1]=CQ0b[2]=CQ1b[1]=CQ1b[2]=CQpb[1]=CQpb[2]=0.;
	double obs[Nobs_BKsll+1];

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_1S/2.;
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	printf("BR(b->s gamma)\t\t\t%.3e\n",bsgamma(C0b,C1b,C2b,mu_b,mu_W,&param));
	
	double lambda_h=0.5;
	double mu_spec=sqrt(lambda_h*param.mass_b);		
	C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0spec,C1spec,mu_spec,&param);
	printf("delta0(B->K* gamma)\t\t%.3e\n\n",delta0(C0b,C0spec,C1b,C1spec,&param,mu_b,mu_spec,lambda_h));
	
	mu_b=param.mass_b_pole;
	C_calculator_base2(C0w,C1w,mu_W,C0b,C1b,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);
	printf("BR(Bs->mu mu)\t\t\t%.3e\n",Bsmumu(C0b,C1b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b));
	printf("BR(Bs->mu mu)_untag\t\t%.3e\n",Bsmumu_untag(C0b,C1b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b));
	printf("BR(Bd->mu mu)\t\t\t%.3e\n\n",Bdmumu(C0b,C1b,CQ0b,CQ1b,&param,mu_b));

	mu_b=param.mass_b_pole;
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);
	
	printf("BR(B->K* mu mu)_low\t\t%.3e\n",BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b));
	printf("AFB(B->K* mu mu)_low\t\t%.3e\n",obs[1]);
	printf("FL(B->K* mu mu)_low\t\t%.3e\n",obs[2]);
	printf("AT1(B->K* mu mu)_low\t\t%.3e\n",obs[4]);
	printf("AT2(B->K* mu mu)_low\t\t%.3e\n",obs[5]);
	printf("AT3(B->K* mu mu)_low\t\t%.3e\n",obs[6]);
	printf("AT4(B->K* mu mu)_low\t\t%.3e\n",obs[7]);
	printf("AT5(B->K* mu mu)_low\t\t%.3e\n",obs[8]);
	printf("AI(B->K* mu mu)_low\t\t%.3e\n\n",AI_BKstarmumu_lowq2(C0b,C1b,C2b,&param,mu_b));
	
	printf("BR(B->K* mu mu)_high\t\t%.3e\n",BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b));
	printf("AFB(B->K* mu mu)_high\t\t%.3e\n",obs[1]);
	printf("FL(B->K* mu mu)_high\t\t%.3e\n",obs[2]);
	printf("HT1(B->K* mu mu)_high\t\t%.3e\n",obs[9]);
	printf("HT2(B->K* mu mu)_high\t\t%.3e\n",obs[10]);
	printf("HT3(B->K* mu mu)_high\t\t%.3e\n",obs[11]);
	printf("AI(B->K* mu mu)_high\t\t%.3e\n\n",AI_BKstarmumu_highq2(C0b,C1b,C2b,&param,mu_b));
		
	double smin=pow(2.*param.mass_mu,2.);
	double smax=pow(param.m_Bd-param.m_Kstar,2.)*0.999; 
	BRBKstarmumu(smin,smax,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
	printf("q0^2(AFB(B->K* mu mu))\t\t%.3e\n",obs[0]);
	printf("q0^2(AI(B->K* mu mu))\t\t%.3e\n\n",AI_BKstarmumu_zero(C0b,C1b,C2b,&param,mu_b));
	
	mu_W=120.;
	mu_b=5.;
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	printf("BR(B->Xs mu mu)_low\t\t%.3e\n",BRBXsmumu_lowq2(C0b,C1b,C2b,CQ0b,CQ1b,&param,mu_b));
	printf("BR(B->Xs mu mu)_high\t\t%.3e\n",BRBXsmumu_highq2(C0b,C1b,C2b,CQ0b,CQ1b,&param,mu_b));
	printf("q0^2(AFB(B->Xs mu mu)\t\t%.3e\n",A_BXsmumu_zero(C0b,C1b,C2b,CQ0b,CQ1b,&param,mu_b));
	printf("BR(B->Xs tau tau)_high\t\t%.3e\n\n",BRBXstautau_highq2(C0b,C1b,C2b,CQ0b,CQ1b,&param,mu_b));
	
	printf("BR(B->tau nu)\t\t\t%.3e\n",Btaunu(&param));
      	printf("R(B->tau nu)\t\t\t%.3e\n",RBtaunu(&param));
      	printf("BR(B->D tau nu)\t\t\t%.3e\n",BDtaunu(&param));
      	printf("BR(B->D tau nu)/BR(B->D e nu)\t%.3e\n",BDtaunu_BDenu(&param));
     	printf("BR(Ds->tau nu)\t\t\t%.3e\n",Dstaunu(&param));
     	printf("BR(Ds->mu nu)\t\t\t%.3e\n",Dsmunu(&param));
     	printf("BR(D->mu nu)\t\t\t%.3e\n",Dmunu(&param));
      	printf("BR(K->mu nu)/BR(pi->mu nu)\t%.3e\n",Kmunu_pimunu(&param));
     	printf("Rmu23(K->mu nu)\t\t\t%.3e\n\n",Rmu23(&param));

	return 1;
}
