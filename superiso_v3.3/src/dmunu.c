#include "include.h"

double Dmunu(struct parameters* param)
/* computes the branching ratio of D -> mu nu */
{
	double Vcd=cabs(param->Vcd);
	
	if(param->SM==1) return 
param->m_D/8./pi*pow(param->Gfermi*Vcd*param->mass_mu*param->f_D*(1.-param->mass_mu*param->mass_mu/param->m_D/param->m_D),2.)*param->life_D/hbar;
	
	double mc=running_mass(param->mass_c,param->mass_c,param->m_D,param->mass_top_pole,param->mass_b,param);
	
	if(param->THDM_model>0) return param->m_D/8./pi*pow(param->Gfermi*Vcd*param->mass_mu*param->f_D*(1.-param->mass_mu*param->mass_mu/param->m_D/param->m_D)*(1.+param->m_D*param->m_D/param->mass_H/param->mass_H*(mc*param->lambda_u[2][2]-param->mass_d*param->lambda_d[1][1])*param->lambda_l[2][2]/(param->mass_d+mc)),2.)*param->life_D/hbar;

	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);
	double epsilon0=-2./3.*alphas_MSOFT/pi*param->mu_Q/param->mass_gluino*H2(param->MqL2_Q*param->MqL2_Q/param->mass_gluino/param->mass_gluino,param->McR_Q*param->McR_Q/param->mass_gluino/param->mass_gluino);
	
	return param->m_D/8./pi*pow(param->Gfermi*Vcd*param->mass_mu*param->f_D*(1.-param->mass_mu*param->mass_mu/param->m_D/param->m_D)*(1.+param->m_D*param->m_D/param->mass_H/param->mass_H*(mc-param->mass_d*param->tan_beta*param->tan_beta/(1.+epsilon0*param->tan_beta))/(param->mass_d+mc)),2.)*param->life_D/hbar;
}

/*--------------------------------------------------------------------*/

double Dmunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(D -> mu nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return Dmunu(&param);
}
