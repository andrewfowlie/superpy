#include "include.h"


double GBDlnu(double w)
{
	double G1=1.03;
	double rho2=1.17;
	double zw=(sqrt(w+1.)-sqrt(2.))/(sqrt(w+1.)+sqrt(2.));
	
	return G1*(1.-8.*rho2*zw+(51.*rho2-10.)*zw*zw-(252.*rho2-84.)*zw*zw*zw);
}

/*--------------------------------------------------------------------*/

double tBDlnu(double w, double m_B, double m_D)
{
	return m_B*m_B+m_D*m_D-2.*w*m_D*m_B;
}

/*--------------------------------------------------------------------*/

double rhoV(double w, double ml, double m_B, double m_D)
{
	return 4.*pow(1.+m_D/m_B,2.)*pow(m_D/m_B,3.)*pow(w*w-1.,1.5)*pow(1.-ml*ml/tBDlnu(w,m_B,m_D),2.)*(1.+ml*ml/2./tBDlnu(w,m_B,m_D))*pow(GBDlnu(w),2.);
}

/*--------------------------------------------------------------------*/

double rhoS(double w, double ml, double m_B, double m_D)
{
	double Deltaw=0.46;
	
	return 1.5*m_B*m_B/tBDlnu(w,m_B,m_D)/(1.+ml*ml/2./tBDlnu(w,m_B,m_D))*(1.+w)/(1.-w)*Deltaw*Deltaw;
}

/*--------------------------------------------------------------------*/

double dGammaBDlnu_dw(double w, double ml, struct parameters* param)
{
	double Vcb=cabs(param->Vcb);

	if(param->SM==1) return param->Gfermi*param->Gfermi*Vcb*Vcb*pow(param->m_B,5.)/192./pow(pi,3.)*rhoV(w,ml,param->m_B,param->m_D0)*(1.-ml*ml/param->m_B/param->m_B*rhoS(w,ml,param->m_B,param->m_D0));
	
	double mc=running_mass(param->mass_c,param->mass_c,param->m_B,param->mass_top_pole,param->mass_b,param);
	double mb=running_mass(param->mass_b,param->mass_b,param->m_B,param->mass_top_pole,param->mass_b,param);

	if(param->THDM_model==0) return param->Gfermi*param->Gfermi*Vcb*Vcb*pow(param->m_B,5.)/192./pow(pi,3.)*rhoV(w,ml,param->m_B,param->m_D0)*	(1.-ml*ml/param->m_B/param->m_B*pow(1.-tBDlnu(w,param->m_B,param->m_D0)/(mb-mc)*mb/param->mass_H/param->mass_H*param->tan_beta*param->tan_beta/(1.+epsilon_0(param)*param->tan_beta),2.)*rhoS(w,ml,param->m_B,param->m_D0));

	else 
	{	
		double lambdal;
		if(fabs(1.-ml/param->mass_e)<1.e-2) lambdal=param->lambda_l[1][1];
		else if(fabs(1.-ml/param->mass_mu)<1.e-2) lambdal=param->lambda_l[2][2];
		else lambdal=param->lambda_l[3][3];

return param->Gfermi*param->Gfermi*Vcb*Vcb*pow(param->m_B,5.)/192./pow(pi,3.)*rhoV(w,ml,param->m_B,param->m_D0)*	(1.-ml*ml/param->m_B/param->m_B*pow(1.-tBDlnu(w,param->m_B,param->m_D0)/(mb-mc)/param->mass_H/param->mass_H*(param->lambda_d[3][3]*mb-param->lambda_u[2][2]*mc)*lambdal,2.)*rhoS(w,ml,param->m_B,param->m_D0));
	}

}

/*--------------------------------------------------------------------*/

double GammaBDlnu(double ml, struct parameters* param)
{
	int ie;
	int nmax=100.;
	double Gamma=0.;
	double w;
	double wmin=1.;
	double wmax=(1.+param->m_D0*param->m_D0/param->m_B/param->m_B-ml*ml/param->m_B/param->m_B)/2./(param->m_D0/param->m_B);
	
	for(ie=1;ie<=nmax;ie++)
	{
		w=wmin+(wmax-wmin)*ie/nmax;
		Gamma+=dGammaBDlnu_dw(w,ml,param);
	}
	Gamma*=(wmax-wmin)/nmax;

	return Gamma;
}

/*--------------------------------------------------------------------*/

double BDtaunu(struct parameters* param)
/* computes the branching ratio of B-> D0 tau nu */
{
	return param->life_B/hbar*GammaBDlnu(param->mass_tau_pole,param);
}

/*--------------------------------------------------------------------*/

double BDtaunu_BDenu(struct parameters* param)
/* computes the ratio BR(B-> D0 tau nu)/BR(B-> D0 e nu) */
{

	return GammaBDlnu(param->mass_tau_pole,param)/GammaBDlnu(param->mass_e,param);
	
}

/*--------------------------------------------------------------------*/

double BDtaunu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> D0 tau nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return BDtaunu(&param);
}

/*--------------------------------------------------------------------*/

double BDtaunu_BDenu_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B-> D0 tau nu)/BR(B-> D0 e nu) */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return BDtaunu_BDenu(&param);
}

