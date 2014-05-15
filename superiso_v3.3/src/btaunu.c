#include "include.h"

double Btaunu(struct parameters* param)
/* computes the branching ratio of B-> tau nu */
{
	double Vub=cabs(param->Vub);

	if(param->SM==1)	return param->life_B/hbar*param->m_B/8./pi*pow(param->Gfermi*Vub*param->mass_tau_pole*param->f_B*(1.-param->mass_tau_pole*param->mass_tau_pole/param->m_B/param->m_B),2.);

	if(param->THDM_model==0) return param->life_B/hbar*param->m_B/8./pi*pow(param->Gfermi*Vub*param->mass_tau_pole*param->f_B*(1.-param->mass_tau_pole*param->mass_tau_pole/param->m_B/param->m_B)	
	*(1.-param->m_B*param->m_B/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_0(param)*param->tan_beta)),2.);

	else return param->life_B/hbar*param->m_B/8./pi*pow(param->Gfermi*Vub*param->mass_tau_pole*param->f_B*(1.-param->mass_tau_pole*param->mass_tau_pole/param->m_B/param->m_B)	
	*(1.-param->m_B*param->m_B/param->mass_H/param->mass_H*param->lambda_l[3][3]*param->lambda_d[3][3]),2.);
}

/*--------------------------------------------------------------------*/

double Btaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> tau nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Btaunu(&param);
}

/*--------------------------------------------------------------------*/

double RBtaunu(struct parameters* param)
/* computes the ratio of BR(B-> tau nu)_MSSM/BR(B-> tau nu)_SM */
{

	if(param->SM==1) return 1.;
	
	if(param->THDM_model==0) return pow(1.-param->m_B*param->m_B/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_0(param)*param->tan_beta),2.);

	else return pow(1.-param->m_B*param->m_B/param->mass_H/param->mass_H*param->lambda_l[3][3]*param->lambda_d[3][3],2.);
}

/*--------------------------------------------------------------------*/

double RBtaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> tau nu)_MSSM/BR(B-> tau nu)_SM */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return RBtaunu(&param);
}
