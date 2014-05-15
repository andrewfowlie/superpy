#include "include.h"


double F1N(double x)
{
	if(x==1.) return 1.;
	return 2./pow(1.-x,4.)*(1.-6.*x+3.*x*x+2.*pow(x,3.)-6.*x*x*log(x));
}

/*--------------------------------------------------------------------*/

double F2N(double x)
{
	if(x==1.) return 1.;
	return 3./pow(1.-x,3.)*(1.-x*x+2.*x*log(x));
}

/*--------------------------------------------------------------------*/

double F1C(double x)
{
	if(x==1.) return 1.;
	return 2./pow(1.-x,4.)*(2.+3.*x-6.*x*x+pow(x,3.)+6.*x*log(x));
}

/*--------------------------------------------------------------------*/

double F2C(double x)
{
	if(x==1.) return 1.;
	return -3./2./pow(1.-x,3.)*(3.-4.*x+x*x+2.*log(x));
}

/*--------------------------------------------------------------------*/

double fPS(double x)
{
	if(x<0.25)
	{
		double y=sqrt(1.-4.*x);
		return 2.*x/y*(Li2(1.-(1.-y)/2./x)-Li2(1.-(1.+y)/2./x));
	}
	else if(x<2.5) return -1.365013496e-1-1.858455623*log(1.+x)-5.996763746e-1*log(1.+x)*log(1.+x)+4.390843985e-1*sqrt(x)*log(1.+x)-1.444359743e-1*x*log(1.+x)+3.852425143*sqrt(x);
	else if(x<100.)  return 4.304425955e-1+6.766323794e-2*log(1.+x)-1.584446296e-1*log(1.+x)*log(1.+x)-2.787080541e-1*sqrt(x)*log(1.+x)+1.557845370e-3*x*log(1.+x)+2.139180566*sqrt(x);
	else if(x<10000.)  return 2.025445594+9.960866255e-1*log(x)+1.122896720e-4*sqrt(x);
	else return 2.000835136+9.9992369e-1*log(x)+2.327105016e-7*sqrt(x);
}

/*--------------------------------------------------------------------*/

double fS(double x)
{
	return (2.*x-1.)*fPS(x)-2.*x*(2.+log(x));
}

/*--------------------------------------------------------------------*/

double fft(double x)
{
	return x/2.*(2.+log(x)-fPS(x));
}

/*--------------------------------------------------------------------*/

double muonf(double z)
{
	if(z<1.e-8) return 0.;
	
	double int1=0.;

	int ie;
	
	int ne=5000;
	double dx=0.5/ne;
	
	double x=0.5;
	int1+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	
	for(ie=1;ie<=ne-2;ie++)
	{	
		x-=dx;
		int1+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z);
	}	
	x=dx;
	int1+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	
	int1*=dx;
	
	double int2=0.;

	x=dx;
	dx=dx/ne;
	int2+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	
	for(ie=1;ie<=ne-1;ie++)
	{	
		x-=dx;
		int2+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z);
	}	
	x=dx;
	int2+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	
	x=1.e-10;
	int2+=(1.-2.*x*(1.-x))/(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	int2*=dx;
		
	return int1+int2;
}

/*--------------------------------------------------------------------*/

double muong(double z)
{
	if(z<1.e-8) return 0.;
	
	double int1=0.;

	int ie;
	
	int ne=5000;
	double dx=0.5/ne;
	
	double x=0.5;
	int1+=1./(x*(1.-x)-z)*log(x*(1.-x)/z)/2.;
	
	for(ie=1;ie<=ne-2;ie++)
	{	
		x-=dx;
		int1+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	}	
	x=dx;
	int1+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	
	int1*=dx;
	
	double int2=0.;

	x=dx;
	dx=dx/ne;
	int2+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	
	for(ie=1;ie<=ne-1;ie++)
	{	
		x-=dx;
		int2+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	}	
	x=dx;
	int2+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	
	x=1.e-10;
	int2+=1./(x*(1.-x)-z)*log(x*(1.-x)/z);
	int2*=dx;
		
	return int1+int2;
}

/*--------------------------------------------------------------------*/

double muonI1(double a)
{
	double int1=0.;

	if(a>=20.)
	{ 
		int1=(3.+a*a)*log(a*a)/pow(a,4.)-7./6./a/a;
		return int1;
	}
	
	int ie;
	
	int ne=500;
	double dx=1./ne;
	
	double x=0.;
	
	for(ie=1;ie<=ne-1;ie++)
	{	
		x=ie*dx;
		int1+=x*x*(2.-x)/(x*x+a*a*(1.-x));
	}	
	
	int1+=1./2.;
	return int1*dx;
}

/*--------------------------------------------------------------------*/

double muonI2(double a)
{
	double int1=0.;

	if(a>=20.)
	{ 
		int1=(5.+a*a)*log(a*a)/pow(a,4.)-11./6./a/a;
		
		return int1;
	}
	
	int ie;
	
	int ne=500;
	double dx=1./ne;
	
	double x=0.;
	
	for(ie=1;ie<=ne-1;ie++)
	{	
		x=ie*dx;
		int1+=x*x*x/(x*x+a*a*(1.-x));
	}
	
	int1+=0.5;
	return int1*dx;
}

/*--------------------------------------------------------------------*/

double muonI3(double a)
{
	double int1=0.;	
	
	if(a>=20.)
	{ 
		int1=-(1.+2.*a*a)/12./pow(a,4.);
		return int1;
	}
	
	int ie;
	
	int ne=500;
	double dx=1./ne;
	
	double x=0.;
	
	for(ie=1;ie<=ne-1;ie++)
	{	
		x=ie*dx;
		int1+=x*(x-1.)/(x-1.+a*a);
	}	
	
	return int1*dx;
}

/*--------------------------------------------------------------------*/

double muon_gm2(struct parameters* param)
/* computes the muon anomalous magnetic moment a_mu */
{
	if(param->SM==1) return 0.;

#ifdef SM_ChargedHiggs
	return 0.;
#endif	
	int ie,me,ke,je;

	if(param->THDM_model>0) 
	{
		double gmuon_2l=0.;
		double mass_f[9];
		double T3_f[9];
		double charge_f[9];
		double ncol_f[9];
		double sba=sin(atan(param->tan_beta)-param->alpha);
		double cba=cos(atan(param->tan_beta)-param->alpha);
		double v=1./sqrt(sqrt(2.)*param->Gfermi);
		
		double rho_f[9];
		double kappa_f[9];
		
		mass_f[0]=param->mass_u;
		
		for(ie=0;ie<=2;ie++) 
		{
			T3_f[ie]=0.5;
			charge_f[ie]=2./3.;
			ncol_f[ie]=3.;
		}
		
		mass_f[3]=param->mass_d;
		mass_f[4]=param->mass_s;
		
		for(ie=3;ie<=5;ie++) 
		{
			T3_f[ie]=-0.5;
			charge_f[ie]=-1./3.;
			ncol_f[ie]=3.;
		}
		
		mass_f[6]=param->mass_e;
		mass_f[7]=param->mass_mu;
		mass_f[8]=param->mass_tau_pole;
		
		for(ie=6;ie<=8;ie++) 
		{
			T3_f[ie]=-0.5;
			charge_f[ie]=-1.;
			ncol_f[ie]=1.;
		}
		
		/* h0 case */
		
		mass_f[1]=running_mass(param->mass_c,param->mass_c,param->mass_h0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[2]=running_mass(param->mtmt,param->mtmt,param->mass_h0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[5]=running_mass(param->mass_b,param->mass_b,param->mass_h0,param->mass_top_pole,param->mass_b_pole,param);
		
		for(ie=0;ie<=8;ie++) kappa_f[ie]=sqrt(2.)*mass_f[ie]/v;
		
		rho_f[0]=kappa_f[0]*param->lambda_u[1][1];
		rho_f[1]=kappa_f[1]*param->lambda_u[2][2];
		rho_f[2]=kappa_f[2]*param->lambda_u[3][3];
		rho_f[3]=kappa_f[3]*param->lambda_d[1][1];
		rho_f[4]=kappa_f[4]*param->lambda_d[2][2];
		rho_f[5]=kappa_f[5]*param->lambda_d[3][3];
		rho_f[6]=kappa_f[6]*param->lambda_l[1][1];
		rho_f[7]=kappa_f[7]*param->lambda_l[2][2];
		rho_f[8]=kappa_f[8]*param->lambda_l[3][3];
		
		for(ie=0;ie<=8;ie++) gmuon_2l+=
		ncol_f[ie]*charge_f[ie]*charge_f[ie]*mass_f[7]*mass_f[ie]/param->mass_h0/param->mass_h0
		*(-(kappa_f[ie]*sba+rho_f[ie]*cba)*(kappa_f[7]*sba+rho_f[7]*cba)/2.*muonf(pow(mass_f[ie]/param->mass_h0,2.)));
		
		/* H0 case */
		
		mass_f[1]=running_mass(param->mass_c,param->mass_c,param->mass_H0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[2]=running_mass(param->mtmt,param->mtmt,param->mass_H0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[5]=running_mass(param->mass_b,param->mass_b,param->mass_H0,param->mass_top_pole,param->mass_b_pole,param);
		
		for(ie=0;ie<=8;ie++) kappa_f[ie]=sqrt(2.)*mass_f[ie]/v;
		
		rho_f[0]=kappa_f[0]*param->lambda_u[1][1];
		rho_f[1]=kappa_f[1]*param->lambda_u[2][2];
		rho_f[2]=kappa_f[2]*param->lambda_u[3][3];
		rho_f[3]=kappa_f[3]*param->lambda_d[1][1];
		rho_f[4]=kappa_f[4]*param->lambda_d[2][2];
		rho_f[5]=kappa_f[5]*param->lambda_d[3][3];
		rho_f[6]=kappa_f[6]*param->lambda_l[1][1];
		rho_f[7]=kappa_f[7]*param->lambda_l[2][2];
		rho_f[8]=kappa_f[8]*param->lambda_l[3][3];
		
		for(ie=0;ie<=8;ie++) gmuon_2l+=
		ncol_f[ie]*charge_f[ie]*charge_f[ie]*mass_f[7]*mass_f[ie]/param->mass_H0/param->mass_H0
		*(-(kappa_f[ie]*cba-rho_f[ie]*sba)*(kappa_f[7]*cba-rho_f[7]*sba)/2.*muonf(pow(mass_f[ie]/param->mass_H0,2.)));
		
		
		/* A0 case */
		
		mass_f[1]=running_mass(param->mass_c,param->mass_c,param->mass_A0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[2]=running_mass(param->mtmt,param->mtmt,param->mass_A0,param->mass_top_pole,param->mass_b_pole,param);
		mass_f[5]=running_mass(param->mass_b,param->mass_b,param->mass_A0,param->mass_top_pole,param->mass_b_pole,param);
		
		for(ie=0;ie<=8;ie++) kappa_f[ie]=sqrt(2.)*mass_f[ie]/v;
		
		rho_f[0]=kappa_f[0]*param->lambda_u[1][1];
		rho_f[1]=kappa_f[1]*param->lambda_u[2][2];
		rho_f[2]=kappa_f[2]*param->lambda_u[3][3];
		rho_f[3]=kappa_f[3]*param->lambda_d[1][1];
		rho_f[4]=kappa_f[4]*param->lambda_d[2][2];
		rho_f[5]=kappa_f[5]*param->lambda_d[3][3];
		rho_f[6]=kappa_f[6]*param->lambda_l[1][1];
		rho_f[7]=kappa_f[7]*param->lambda_l[2][2];
		rho_f[8]=kappa_f[8]*param->lambda_l[3][3];
		
		for(ie=0;ie<=8;ie++) gmuon_2l+=
		ncol_f[ie]*charge_f[ie]*charge_f[ie]*mass_f[7]*mass_f[ie]/param->mass_A0/param->mass_A0
		*(-2.*T3_f[ie]*rho_f[ie]*rho_f[7]/2.*muong(pow(mass_f[ie]/param->mass_A0,2.)));
		
		return gmuon_2l/param->inv_alpha_em/4./pi/pi/pi;
	
	}

	/* SUSY */

	double mass_smu[3],mass_charg[3];
	double n_L[6][3],n_R[6][3],xim[6][3],xk[3],c_L[3],c_R[3],X[3][3];	
	int nb_neut;
	if(param->mass_neut[5]==0.) nb_neut=4; else nb_neut=5;


	mass_charg[1]=param->mass_cha1;
	mass_charg[2]=param->mass_cha2;


	double M2smu11=pow(param->mass_mul,2.);
	double M2smu22=pow(param->mass_mur,2.);
	double M2smu12=(param->A_mu-param->mu_Q*param->tan_beta)*param->mass_mu;

	mass_smu[1]=sqrt((M2smu11+M2smu22-sqrt(pow(M2smu11-M2smu22,2.)+4.*M2smu12*M2smu12))/2.);
	mass_smu[2]=sqrt((M2smu11+M2smu22+sqrt(pow(M2smu11-M2smu22,2.)+4.*M2smu12*M2smu12))/2.);

	X[1][1]=cos(atan2(-2.*M2smu12,-M2smu11+M2smu22)/2.);
	X[1][2]=sin(atan2(-2.*M2smu12,-M2smu11+M2smu22)/2.);
	X[2][1]=-X[1][2];
	X[2][2]=X[1][1];
	
	double ymu=param->g2*param->mass_mu/sqrt(2.)/param->mass_W/cos(atan(param->tan_beta));
	
	for(ie=1;ie<=nb_neut;ie++) for(me=1;me<=2;me++)
	{
		n_L[ie][me]=1./sqrt(2.)*(param->gp*param->neut_mix[ie][1]+param->g2*param->neut_mix[ie][2])*X[me][1]
		-ymu*param->neut_mix[ie][3]*X[me][2];
		n_R[ie][me]=sqrt(2.)*param->gp*param->neut_mix[ie][1]*X[me][2]+ymu*param->neut_mix[ie][3]*X[me][1];
		xim[ie][me]=pow(param->mass_neut[ie]/mass_smu[me],2.);
	}
	
	for(ke=1;ke<=2;ke++)
	{
		c_L[ke]=-param->g2*param->charg_Vmix[ke][1];
		c_R[ke]=ymu*param->charg_Umix[ke][2];
		xk[ke]=pow(mass_charg[ke]/param->mass_numl,2.);
	}
		
	double a_neut=0.;
	for(ie=1;ie<=nb_neut;ie++) for(me=1;me<=2;me++)
	{
		a_neut+=
		-param->mass_mu/12./pow(mass_smu[me],2.)*(pow(n_L[ie][me],2.)+pow(n_R[ie][me],2.))*F1N(xim[ie][me])
		+param->mass_neut[ie]/3./pow(mass_smu[me],2.)*n_L[ie][me]*n_R[ie][me]*F2N(xim[ie][me]);
	}
	a_neut*=param->mass_mu/16./pi/pi;
	
	/* --- */
		
	double a_charg=0.;
	for(ke=1;ke<=2;ke++)
	{
		a_charg+=
		param->mass_mu/12./pow(param->mass_numl,2.)*(pow(c_L[ke],2.)+pow(c_R[ke],2.))*F1C(xk[ke])
		+2.*mass_charg[ke]/3./pow(param->mass_numl,2.)*c_L[ke]*c_R[ke]*F2C(xk[ke]);
	}
	a_charg*=param->mass_mu/16./pi/pi;

	/* --- */

	double G_mu=param->g2/4./sqrt(2.)/param->mass_W/param->mass_W;
	
	double c=G_mu*param->mass_mu*param->mass_mu/4./sqrt(2.)/pi/pi;
	
	double S12=-sin(param->alpha);	
	double S22=cos(param->alpha);
	double P12=1.;
	double mA=param->mass_A0;
	double a_HiggsNMSSM=0.;
	
	if((param->mass_A02!=0.)&&(param->mass_H03!=0.))
	{
		mA=param->mass_A02;
		S12=param->H0_mix[1][2];
		S22=param->H0_mix[2][2];
		double S23=param->H0_mix[2][3];
		double P11=param->A0_mix[1][1];
		P12=param->A0_mix[1][2];
		
		double a_H03=c*S23*S23*(1.+param->tan_beta*param->tan_beta)*muonI1(param->mass_H03*param->mass_H03/param->mass_mu/param->mass_mu);
		
		double a_a1=-c*P11*P11*param->tan_beta*param->tan_beta*muonI2(param->mass_A0*param->mass_A0/param->mass_mu/param->mass_mu);
		
		a_HiggsNMSSM+=a_H03+a_a1;
	}
		
	double a_h0=c*S12*S12*(1.+param->tan_beta*param->tan_beta)*muonI1(param->mass_h0*param->mass_h0/param->mass_mu/param->mass_mu);
		
	double a_H0=c*S22*S22*(1.+param->tan_beta*param->tan_beta)*muonI1(param->mass_H0*param->mass_H0/param->mass_mu/param->mass_mu);
	
	double a_A0=-c*P12*P12*param->tan_beta*param->tan_beta*muonI2(mA*mA/param->mass_mu/param->mass_mu);
		
	double a_Hpm=c*param->tan_beta*param->tan_beta*muonI3(param->mass_H*param->mass_H/param->mass_mu/param->mass_mu);

	/* --- */

	double a_SUSYQED=(a_neut+a_charg+a_h0+a_H0+a_A0+a_Hpm+a_HiggsNMSSM)*(1.-4./param->inv_alpha_em/pi*log(param->MSOFT_Q/param->mass_mu));
	
	/* --- */

	double mass_stop[3],mass_sbot[3],mass_stau[3],mass_H[4],mass_A[3];
	double sw2=pow(sin(atan(param->gp/param->g2)),2.);

	mass_H[1]=param->mass_h0;
	mass_H[2]=param->mass_H0;
	mass_H[3]=param->mass_H03;
	
	mass_A[1]=param->mass_A0;
	mass_A[2]=param->mass_A02;

	mass_stop[1]=param->mass_t1;
	mass_stop[2]=param->mass_t2;

	mass_sbot[1]=param->mass_b1;
	mass_sbot[2]=param->mass_b2;

	mass_stau[1]=param->mass_tau1;
	mass_stau[2]=param->mass_tau2;
	
	double a_charH=0.;
	double a_sfH=0.;
	double a_bosEW=0.;
	double c_Lh=0.;
	
	double beta=atan(param->tan_beta);
	
	if((param->mass_A02==0.)&&(param->mass_H03==0.))
	{
		double lmu[4],lcharg[4][3],lst[3][3],lsb[3][3],lstau[3][3];
		double sa=sin(param->alpha);
		double ca=cos(param->alpha);
		double sb=sin(beta);
		double cb=cos(beta);
	
		lmu[1]=-sa/cb;
		lmu[2]=ca/cb;
		lmu[3]=param->tan_beta;
	
		for(ie=1;ie<=2;ie++)
		{
			lcharg[1][ie] = sqrt(2.)*param->mass_W/mass_charg[ie]*(param->charg_Umix[ie][1]*param->charg_Vmix[ie][2]*ca - param->charg_Umix[ie][2]*param->charg_Vmix[ie][1]*sa);

			lcharg[2][ie] = sqrt(2.)*param->mass_W/mass_charg[ie]*(param->charg_Umix[ie][1]*param->charg_Vmix[ie][2]*sa + param->charg_Umix[ie][2]*param->charg_Vmix[ie][1]*ca);
		
			lcharg[3][ie] = -sqrt(2.)*param->mass_W/mass_charg[ie]*(param->charg_Umix[ie][1]*param->charg_Vmix[ie][2]*cb + param->charg_Umix[ie][2]*param->charg_Vmix[ie][1]*sb);
		
			lst[1][ie] = 2.*param->mtmt/mass_stop[ie]/mass_stop[ie]/sb*(param->mu_Q*sa+param->A_t*ca)*param->stop_mix[ie][1]*param->stop_mix[ie][2];
		
			lst[2][ie] = 2.*param->mtmt/mass_stop[ie]/mass_stop[ie]/sb*(-param->mu_Q*ca+param->A_t*sa)*param->stop_mix[ie][1]*param->stop_mix[ie][2];

			lsb[1][ie] = 2.*param->mass_b/mass_sbot[ie]/mass_sbot[ie]/cb*(-param->mu_Q*ca-param->A_b*sa)*param->sbot_mix[ie][1]*param->sbot_mix[ie][2];
		
			lsb[2][ie] = 2.*param->mass_b/mass_sbot[ie]/mass_sbot[ie]/cb*(-param->mu_Q*sa+param->A_b*ca)*param->sbot_mix[ie][1]*param->sbot_mix[ie][2];
		
			lstau[1][ie] = 2.*param->mass_tau_pole/mass_stau[ie]/mass_stau[ie]/cb*(-param->mu_Q*ca-param->A_tau*sa)*param->stau_mix[ie][1]*param->stau_mix[ie][2];
		
			lstau[2][ie] = 2.*param->mass_tau_pole/mass_stau[ie]/mass_stau[ie]/cb*(-param->mu_Q*sa+param->A_tau*ca)*param->stau_mix[ie][1]*param->stau_mix[ie][2];
		}

		for(ke=1;ke<=2;ke++)
		{
			for(je=1;je<=2;je++) a_charH+=lmu[je]*lcharg[je][ke] *fS(mass_charg[ke]*mass_charg[ke]/mass_H[je]/mass_H[je]);
		
			a_charH += lmu[3]*lcharg[3][ke]*fPS(mass_charg[ke]*mass_charg[ke]/param->mass_A0/param->mass_A0);	
		}
		
		double const1=param->mass_mu*param->mass_mu/8./pi/pi/param->inv_alpha_em/param->inv_alpha_em/param->mass_W/param->mass_W/sw2;
		a_charH*=const1;
	
	/* --- */

		for(ie=1;ie<=2;ie++) for(ke=1;ke<=2;ke++)
		{
			a_sfH+=4./3.*lmu[ke]*lst[ke][ie]*fft(mass_stop[ie]*mass_stop[ie]/mass_H[ke]/mass_H[ke]);
			a_sfH+=1./3.*lmu[ke]*lsb[ke][ie]*fft(mass_sbot[ie]*mass_sbot[ie]/mass_H[ke]/mass_H[ke]);
			a_sfH+=lmu[ke]*lstau[ke][ie]*fft(mass_stau[ie]*mass_stau[ie]/mass_H[ke]/mass_H[ke]);
		}
		a_sfH*=const1;
		
		c_Lh=cos(2.*beta)*param->mass_Z*param->mass_Z/cb*(ca*cos(param->alpha+beta)/param->mass_H0/param->mass_H0+sa*sin(param->alpha+beta)/param->mass_h0/param->mass_h0);
		
	}
	else
	{
		double lmuh[4],lmua[3],lchargh[4][3],lcharga[3][3],lsth[4][3],lsbh[4][3],lstauh[4][3];
		double cb=cos(beta);
		double sb=sin(beta);
		double s=param->lambdaSNMSSM/param->lambdaNMSSM;
		double vu=sqrt(sb*sb/sqrt(2.)/param->Gfermi);
		double vd=vu/param->tan_beta;
	
		for(ie=1;ie<=3;ie++) lmuh[ie]=param->H0_mix[ie][2]/cb;
		for(ie=1;ie<=2;ie++) lmua[ie]=param->A0_mix[ie][2]*param->tan_beta;
		
	
		for(ke=1;ke<=2;ke++) for(ie=1;ie<=3;ie++)
		{
			lchargh[ie][ke] = sqrt(2.)*param->mass_W/mass_charg[ke]/param->g2*
			(param->lambdaNMSSM*param->charg_Umix[ke][2]*param->charg_Vmix[ke][2]*param->H0_mix[ie][3]
			+param->g2*(param->charg_Umix[ke][1]*param->charg_Vmix[ke][2]*param->H0_mix[ie][1]
			+param->charg_Umix[ke][2]*param->charg_Vmix[ke][1]*param->H0_mix[ie][2]));
			
		
			lsth[ie][ke] = 2.*sqrt(2.)*param->mass_W/param->g2/mass_stop[ke]/mass_stop[ke]*(
				      param->yut[3]*(param->A_t*param->H0_mix[ie][1]-param->lambdaNMSSM*(s*param->H0_mix[ie][2]+vd*param->H0_mix[ie][3]))*param->stop_mix[ke][1]*param->stop_mix[ke][2]
				      + (param->yut[3]*param->yut[3]*vu*param->H0_mix[ie][1]-param->gp*param->gp/3.*(vu*param->H0_mix[ie][1]-vd*param->H0_mix[ie][2]))*param->stop_mix[ke][2]*param->stop_mix[ke][2]
				      +(param->yut[3]*param->yut[3]*vu*param->H0_mix[ie][1]-(3.*param->g2*param->g2-param->gp*param->gp)/12.*(vu*param->H0_mix[ie][1]-vd*param->H0_mix[ie][2]))*param->stop_mix[ke][1]*param->stop_mix[ke][1]);

			
			 				
			lsbh[ie][ke] = 2.*sqrt(2.)*param->mass_W/param->g2/mass_sbot[ke]/mass_sbot[ke]*(
				      param->yub[3]*(param->A_b*param->H0_mix[ie][2]-param->lambdaNMSSM*(s*param->H0_mix[ie][1]+vu*param->H0_mix[ie][3]))*param->sbot_mix[ke][1]*param->sbot_mix[ke][2]
				      + (param->yub[3]*param->yub[3]*vd*param->H0_mix[ie][2]+param->gp*param->gp/6.*(vu*param->H0_mix[ie][1]-vd*param->H0_mix[ie][2]))*param->sbot_mix[ke][2]*param->sbot_mix[ke][2]
				      +(param->yub[3]*param->yub[3]*vd*param->H0_mix[ie][2]+(3.*param->g2*param->g2+param->gp*param->gp)/12.*(vu*param->H0_mix[ie][1]-vd*param->H0_mix[ie][2]))*param->sbot_mix[ke][1]*param->sbot_mix[ke][1]);
			
			
			lstauh[ie][ke] = 2.*sqrt(2.)*param->mass_W/param->g2/mass_stau[ke]/mass_stau[ke]*(
				      param->yutau[3]*(param->A_tau*param->H0_mix[ie][2]-param->lambdaNMSSM*(s*param->H0_mix[ie][1]+vu*param->H0_mix[ie][3]))*param->stau_mix[ke][1]*param->stau_mix[ke][2]
				      + (param->yutau[3]*param->yutau[3]*vd*param->H0_mix[ie][2]+param->gp*param->gp/2.*(vu*param->H0_mix[ie][1]-vd*param->H0_mix[ie][2]))*param->stau_mix[ke][2]*param->stau_mix[ke][2]
				      +(param->yutau[3]*param->yutau[3]*vd*param->H0_mix[ie][2]+(param->g2*param->g2-param->gp*param->gp)/4.*(vu*param->H0_mix[ie][1]-vd*param->H0_mix[ie][2]))*param->stau_mix[ke][1]*param->stau_mix[ke][1]);
			
		}
		
		for(ke=1;ke<=2;ke++) for(ie=1;ie<=2;ie++) lcharga[ie][ke] =
			sqrt(2.)*param->mass_W/mass_charg[ke]/param->g2*
			(param->lambdaNMSSM*param->charg_Umix[ke][2]*param->charg_Vmix[ke][2]*param->A0_mix[ie][2]
			-param->g2*(param->charg_Umix[ke][1]*param->charg_Vmix[ke][2]*cos(beta)
			+param->charg_Umix[ke][2]*param->charg_Vmix[ke][1]*sin(beta))*param->A0_mix[ie][1]);

	
	/* --- */

		for(ke=1;ke<=2;ke++)
		{
			for(je=1;je<=3;je++) a_charH+=lmuh[je]*lchargh[je][ke] *fS(mass_charg[ke]*mass_charg[ke]/mass_H[je]/mass_H[je]);
		
			for(je=1;je<=2;je++) a_charH += lmua[je]*lcharga[je][ke]*fPS(mass_charg[ke]*mass_charg[ke]/mass_A[je]/mass_A[je]);	
		}
		
		double const1=param->mass_mu*param->mass_mu/8./pi/pi/param->inv_alpha_em/param->inv_alpha_em/param->mass_W/param->mass_W/sw2;
		a_charH*=const1;
	
	/* --- */


		for(ie=1;ie<=2;ie++) for(ke=1;ke<=3;ke++)
		{
			a_sfH+=4./3.*lmuh[ke]*lsth[ke][ie]*fft(mass_stop[ie]*mass_stop[ie]/mass_H[ke]/mass_H[ke]);
			a_sfH+=1./3.*lmuh[ke]*lsbh[ke][ie]*fft(mass_sbot[ie]*mass_sbot[ie]/mass_H[ke]/mass_H[ke]);
			a_sfH+=lmuh[ke]*lstauh[ke][ie]*fft(mass_stau[ie]*mass_stau[ie]/mass_H[ke]/mass_H[ke]);
		}
		a_sfH*=const1;
		
	/* --- */
	
		for(ie=1;ie<=3;ie++) c_Lh+=param->H0_mix[ie][2]*(param->H0_mix[ie][2]-param->tan_beta*param->H0_mix[ie][1])/mass_H[ie]/mass_H[ie];

		c_Lh*=cos(2.*beta)*param->mass_Z*param->mass_Z;
	
	}
	
	/* --- */
	
	double cbos_L=(98.+9.*c_Lh+23.*pow(1.-4.*sw2,2.))/30.;
	double cbos_L_SM=(107.+23.*pow(1.-4.*sw2,2.))/30.;
	double cbos_0=0;
	double cbos_0_SM=0;
	 a_bosEW=5.*param->Gfermi*param->mass_mu*param->mass_mu/24./sqrt(2.)/pi/pi/pi/param->inv_alpha_em*((cbos_L-cbos_L_SM)*log(param->mass_mu*param->mass_mu/param->mass_W/param->mass_W)+cbos_0-cbos_0_SM); 

	return a_SUSYQED+a_charH+a_sfH+a_bosEW;
}

/*--------------------------------------------------------------------*/

double muon_gm2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating a_mu */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	return muon_gm2(&param);
}
