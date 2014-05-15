#include "include.h"

/*----------------------------------------------------------------------*/

double complex h_bkll(double q2, double mq, double mu)
{
	if(mq==0.) return 4./9.*(2./3.+I*pi+log(mu*mu/q2));
	
	double z=4.*mq*mq/q2;
	
	if(z>1.) return -4./9.*(log(mq*mq/mu/mu)-2./3.-z)
	-4./9.*(2.+z)*sqrt(z-1.)*atan(1./sqrt(z-1.));
	
	else return -4./9.*(log(mq*mq/mu/mu)-2./3.-z)
	-4./9.*(2.+z)*sqrt(1.-z)*(log((1.+sqrt(1.-z))/sqrt(z))-I*pi/2.);
}

/*----------------------------------------------------------------------*/

double phi_Kstar(double u, double a1, double a2)
{
	double x=2.*u-1.;
	double C1=3.*x;
	double C2=-1.5+15./2.*x*x;

	return 6.*u*(1.-u)*(1.+a1*C1+a2*C2);
}

/*----------------------------------------------------------------------*/

double complex B0_bkll(double s, double mq)
{
	double epsilon=1.e-10;
	return -2.*csqrt(4.*(mq*mq-I*epsilon)/s-1.)*catan(1./csqrt(4.*(mq*mq-I*epsilon)/s-1.));
}

/*----------------------------------------------------------------------*/

double complex L1_bkll(double complex x)
{
	return clog((x-1.)/x)*clog(1.-x)-pi*pi/6.+CLi2(x/(x-1.));
}

/*----------------------------------------------------------------------*/

double complex I1_bkll(double u, double mq, double q2, struct parameters* param)
{
	if(mq==0.) return 1.;
	
	double epsilon=1.e-10;

	double complex xp=0.5+csqrt(0.25-(mq*mq-I*epsilon)/((1.-u)*param->m_Bd*param->m_Bd+u*q2));
	double complex xm=0.5-csqrt(0.25-(mq*mq-I*epsilon)/((1.-u)*param->m_Bd*param->m_Bd+u*q2));
	double complex yp=0.5+csqrt(0.25-(mq*mq-I*epsilon)/q2);
	double complex ym=0.5-csqrt(0.25-(mq*mq-I*epsilon)/q2);

	return 1.+2.*mq*mq/(1.-u)/(param->m_Bd*param->m_Bd-q2)*(L1_bkll(xp)+L1_bkll(xm)-L1_bkll(yp)-L1_bkll(ym));
}


/*----------------------------------------------------------------------*/

double complex tperp_bkll(double u, double mq, double q2, double E_Kstar, struct parameters* param)
{
	return 2.*param->m_Bd/(1.-u)/E_Kstar*I1_bkll(u,mq,q2,param)+q2/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(B0_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,mq)-B0_bkll(q2,mq));
}

/*----------------------------------------------------------------------*/

double complex tpar_bkll(double u, double mq, double q2, double E_Kstar, struct parameters* param)
{
	return 2.*param->m_Bd/(1.-u)/E_Kstar*I1_bkll(u,mq,q2,param)+((1.-u)*param->m_Bd*param->m_Bd+u*q2)/(1.-u)/(1.-u)/E_Kstar/E_Kstar*(B0_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,mq)-B0_bkll(q2,mq));
}

/*----------------------------------------------------------------------*/

double dGamma_BKstarmumu_dq2(double q2, double obs[][3], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	double mc=mc_pole(param);

	int ie,je;
	double shat=q2/param->m_Bd/param->m_Bd;
	
	double beta_l=sqrt(1.-4.*param->mass_mu*param->mass_mu/q2);

	double alpha_em=1./133.;

	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	
	double mu_f=sqrt(mu_b*0.5);
	
	double alphas_muf=alphas_running(mu_f,param->mass_top_pole,param->mass_b_pole,param);
	double eta=alphas_muf/alphas_running(1.,param->mass_top_pole,param->mass_b_pole,param);

	double Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_Kstar=(param->m_Bd*param->m_Bd+param->m_Kstar*param->m_Kstar-q2)/2./param->m_Bd;

	int nf=5;
	double f_K_perp=param->f_K_perp;
	f_K_perp*=pow(eta,4./3./(11.-2./3.*nf));

	double f_K_par=param->f_K_par;

	double V=0.923/(1.-q2/5.32/5.32)-0.511/(1.-q2/49.4);
	double A1=0.290/(1.-q2/40.38);
	double A2=-0.084/(1.-q2/52.)+0.342/(1.-q2/52.)/(1.-q2/52.);
	double xi_par=(param->m_Bd+param->m_Kstar)/2./E_Kstar*A1-(param->m_Bd-param->m_Kstar)/param->m_Bd*A2;
	double xi_perp=param->m_Bd/(param->m_Bd+param->m_Kstar)*V;
	
	double complex ALperp=0.;
	double complex ARperp=0.;
	double complex ALpar=0.;
	double complex ARpar=0.;
	double complex AL0=0.;
	double complex AR0=0.;
	double complex At=0.;
	double complex AS=0.;

	if(q2<7.)
	{
	 	double C7eff=Cmub[7];
		double C8eff=Cmub[8];
		double C9=Cmub[9];
		double C10=Cmub[10];
	
		double C7effp=Cpb[7];
		double C9p=Cpb[9];
		double C10p=Cpb[10];
	
		double C1bar=Cmub[1]/2.;
		double C2bar=Cmub[2]-Cmub[1]/6.;
		double C3bar=Cmub[3]-Cmub[4]/6.+16.*Cmub[5]-8./3.*Cmub[6];
		double C4bar=Cmub[4]/2.+8.*Cmub[6];
		double C5bar=Cmub[3]-Cmub[4]/6.+4.*Cmub[5]-2./3.*Cmub[6];
		double C6bar=Cmub[4]/2.+2.*Cmub[6];
	
		double CQ1=CQ0b[1]+alphas_mub/4./pi*CQ1b[1];
		double CQ2=CQ0b[2]+alphas_mub/4./pi*CQ1b[2];
		double CQ1p=CQpb[1];
		double CQ2p=CQpb[2];

		double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
		double mb=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */

		double complex Y=4./3.*Cmub[3]+64./9.*Cmub[5]+64./27.*Cmub[6]
		+h_bkll(q2,mc,mu_b)*(4./3.*Cmub[1]+Cmub[2]+6.*Cmub[3]+60.*Cmub[5])
		+h_bkll(q2,param->mass_b_pole,mu_b)*(-7./2.*Cmub[3]-2./3.*Cmub[4]-38.*Cmub[5]-32./3.*Cmub[6])
		+h_bkll(q2,0.,mu_b)*(-1./2.*Cmub[3]-2./3.*Cmub[4]-8.*Cmub[5]-32./3.*Cmub[6]);

		double complex Cperpp0=C7eff+C7effp+q2/2./mb/param->m_Bd*Y;
		double complex Cperpm0=C7eff-C7effp+q2/2./mb/param->m_Bd*Y;
		double complex Cparm0=-(C7eff-C7effp)-param->m_Bd/2./mb*Y;

		double DeltaM=6.*log(mb/mu_b)-4.*(1.-mu_f/mb);
		double L=-(mb*mb-q2)/q2*log(1.-q2/mb/mb);
	
		double mchat=mc/param->m_Bd;
		double z=mchat*mchat;	

		double complex Cperppf=(C7eff+C7effp)*(2.*log(mb/mu_b)-L+DeltaM);
		double complex Cperpmf=(C7eff-C7effp)*(2.*log(mb/mu_b)-L+DeltaM);
		double complex Cparmf=-(C7eff-C7effp)*(2.*log(mb/mu_b)+2.*L+DeltaM);

		double complex Cperpnf=(-C2bar*F27_bsll(shat,z,log(mu_b/mb))-C8eff*F87_bsll(shat,log(mu_b/mb))
-q2/2./mb/param->m_Bd*(C2bar*F29_bsll(shat,z,log(mu_b/mb))+2.*C1bar*(F19_bsll(shat,z,log(mu_b/mb))+1./6.*F29_bsll(shat,z,log(mu_b/mb)))+C8eff*F89_bsll(shat)))/4.*3.;
		double complex Cparnf=(C2bar*F27_bsll(shat,z,log(mu_b/mb))+C8eff*F87_bsll(shat,log(mu_b/mb))
+param->m_Bd/2./mb*(C2bar*F29_bsll(shat,z,log(mu_b/mb))+2.*C1bar*(F19_bsll(shat,z,log(mu_b/mb))+1./6.*F29_bsll(shat,z,log(mu_b/mb)))+C8eff*F89_bsll(shat)))/4.*3.;

		double complex Cperpp1=Cperppf+Cperpnf;
		double complex Cperpm1=Cperpmf+Cperpnf;
		double complex Cparm1=Cparmf+Cparnf;
	
		double complex Cperpp=Cperpp0+alphas_mub*4./3./4./pi*Cperpp1;
		double complex Cperpm=Cperpm0+alphas_mub*4./3./4./pi*Cperpm1;
		double complex Cparm=Cparm0+alphas_mub*4./3./4./pi*Cparm1;

		double Xi_perp=1.;
		double Xi_par=param->m_Kstar/E_Kstar;
	
		double eu=2./3.;
		double ed=-1./3.;
		double eq=-1./3.;
	
		double a1perp=param->a1perp;
		double a2perp=param->a2perp;
		double a1par=param->a1par;
		double a2par=param->a2par;

		a1perp*=pow(eta,4./(11.-2./3.*nf));
		a2perp*=pow(eta,4./3.*(1.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

		a1par*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
		a2par*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

		double u;
	
		double complex int_perpp,int_perpm,int_parm;
		double complex Tperpp0,Tperppf,Tperppnf,Tperpp;
		double complex Tperpm0,Tperpmf,Tperpmnf,Tperpm;
		double complex Tparp0,Tparpf,Tparpnf,Tparp;
		double complex Tparm0,Tparmf,Tparmnf,Tparm;

		double lambda_Bp=param->lambda_Bp;
		lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_b,2.))*(1.-2.*1.4);

		double omega0=2.*(param->m_Bd-mb)/3.;
		double complex lambda_Bm=1./(exp(-q2/param->m_Bd/omega0)/omega0*(-Ei(q2/param->m_Bd/omega0)+I*pi));

		int n1=50;
		int_perpp=int_perpm=0.;
		for(ie=1;ie<=n1-1;ie++)
		{
			u=(double)ie/n1;
		
			Tperpp0=0.;
			Tperppf=(C7eff+C7effp)*2.*param->m_Bd/E_Kstar/(1.-u);

			Tperpm0=0.;
			Tperpmf=(C7eff-C7effp)*2.*param->m_Bd/E_Kstar/(1.-u);

			Tperppnf=-4.*ed*C8eff/(u+(1.-u)*q2/param->m_Bd/param->m_Bd)
			+param->m_Bd/2./mb*(eu*tperp_bkll(u,mc,q2,E_Kstar,param)*(C2bar+C4bar-C6bar)
			+ed*tperp_bkll(u,mb,q2,E_Kstar,param)*(C3bar+C4bar-C6bar-4.*mb/param->m_Bd*C5bar)
			+ed*tperp_bkll(u,param->mass_d,q2,E_Kstar,param)*C3bar);
		
			Tperpmnf=Tperppnf;
		
			Tperpp=Tperpp0+alphas_muf*4./3./4./pi*(Tperppf+Tperppnf);
			Tperpm=Tperpp0+alphas_muf*4./3./4./pi*(Tperpmf+Tperpmnf);
	
			int_perpp+=phi_Kstar(u,a1perp,a2perp)*Tperpp/n1;
			int_perpm+=phi_Kstar(u,a1perp,a2perp)*Tperpm/n1;
		}
		int_perpp*=1./lambda_Bp;	
		int_perpm*=1./lambda_Bp;	
	
		double complex Tauperpp=xi_perp*Cperpp+pi*pi/3.*param->f_B*f_K_perp/param->m_Bd*Xi_perp*int_perpp;
		double complex Tauperpm=xi_perp*Cperpm+pi*pi/3.*param->f_B*f_K_perp/param->m_Bd*Xi_perp*int_perpm;

		int n2=50;
		int_parm=0.;
		for(ie=1;ie<=n2-1;ie++)
		{
			u=(double)ie/n2;
		
			Tparp0=0.;
		
			Tparpf=(C7eff-C7effp)*4.*param->m_Bd/E_Kstar/(1.-u);

			Tparpnf=param->m_Bd/mb*(eu*tpar_bkll(u,mc,q2,E_Kstar,param)*(C2bar+C4bar-C6bar)
			+ed*tpar_bkll(u,mb,q2,E_Kstar,param)*(C3bar+C4bar-C6bar)
			+ed*tpar_bkll(u,param->mass_d,q2,E_Kstar,param)*C3bar);
	
			Tparp=Tparp0+alphas_muf*4./3./4./pi*(Tparpf+Tparpnf);
	
	
			Tparm0=-eq*4.*param->m_Bd/mb*(C3bar+3.*C4bar);
		
			Tparmf=0.;

			Tparmnf=eq*(8.*C8eff/((1.-u)+u*q2/param->m_Bd/param->m_Bd)
			+6.*param->m_Bd/mb*(h_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,mc,mu_b)*(C2bar+C4bar+C6bar)
			+h_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,param->mass_b_pole,mu_b)*(C3bar+C4bar+C6bar)
			+h_bkll((1.-u)*param->m_Bd*param->m_Bd+u*q2,0.,mu_b)*(C3bar+3.*C4bar+3.*C6bar)
			-8./27.*(C3bar-C5bar-15.*C6bar)));
	
			Tparm=Tparm0+alphas_muf*4./3./4./pi*(Tparmf+Tparmnf);
	
			int_parm+=(phi_Kstar(u,a1par,a2par)*(Tparp/lambda_Bp+Tparm/lambda_Bm))/n1;
		}

		double complex Tauparm=xi_par*Cparm+pi*pi/3.*param->f_B*f_K_par/param->m_Bd*Xi_par*int_parm;

		double lambda=pow(param->m_Bd,4.)+pow(param->m_Kstar,4.)+q2*q2-2.*(param->m_Bd*param->m_Bd*param->m_Kstar*param->m_Kstar+param->m_Kstar*param->m_Kstar*q2+param->m_Bd*param->m_Bd*q2);
	
		double N=sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)/param->m_Bd*pow(cabs(param->Vtb*conj(param->Vts)),2.)*shat*sqrt(lambda)*beta_l);
	
		ALperp=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)-(C10+C10p))*V/(param->m_Bd+param->m_Kstar)+2.*mb/q2*Tauperpp);
	
		ARperp=N*sqrt(2.)*sqrt(lambda)*(((C9+C9p)+(C10+C10p))*V/(param->m_Bd+param->m_Kstar)+2.*mb/q2*Tauperpp);
	
		ALpar=-N*sqrt(2.)*(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(((C9-C9p)-(C10-C10p))*A1/(param->m_Bd-param->m_Kstar)+4.*mb/param->m_Bd*E_Kstar/q2*Tauperpm);
	
		ARpar=-N*sqrt(2.)*(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(((C9-C9p)+(C10-C10p))*A1/(param->m_Bd-param->m_Kstar)+4.*mb/param->m_Bd*E_Kstar/q2*Tauperpm);
	
		AL0=-N/2./param->m_Kstar/sqrt(q2)*(((C9-C9p)-(C10-C10p))*((param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar-q2)*(param->m_Bd+param->m_Kstar)*A1-lambda*A2/(param->m_Bd+param->m_Kstar))+2.*mb*(2.*E_Kstar/param->m_Bd*(param->m_Bd*param->m_Bd+3.*param->m_Kstar*param->m_Kstar-q2)*Tauperpm-lambda/(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(Tauperpm+Tauparm)));
	
		AR0=-N/2./param->m_Kstar/sqrt(q2)*(((C9-C9p)+(C10-C10p))*((param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar-q2)*(param->m_Bd+param->m_Kstar)*A1-lambda*A2/(param->m_Bd+param->m_Kstar))+2.*mb*(2.*E_Kstar/param->m_Bd*(param->m_Bd*param->m_Bd+3.*param->m_Kstar*param->m_Kstar-q2)*Tauperpm-lambda/(param->m_Bd*param->m_Bd-param->m_Kstar*param->m_Kstar)*(Tauperpm+Tauparm)));

		At=N/sqrt(q2)*sqrt(lambda)*(2.*(C10-C10p)+q2/param->mass_mu*(CQ2-CQ2p)/param->mass_b_pole)*E_Kstar/param->m_Kstar*xi_par;
	
		AS=-2.*N*sqrt(lambda)*(CQ1-CQ1p)/param->mass_b_pole*E_Kstar/param->m_Kstar*xi_par;
		
	}
	else
	{
		double mb=running_mass(param->mass_b,param->mass_b,mu_b,param->mass_top_pole,param->mass_b,param);
		
		double z=4.*mb*mb/q2;
		
		double complex x1=0.5+0.5*I*csqrt(z-1.);
		double complex x2=0.5-0.5*I*csqrt(z-1.);
		double complex x3=0.5+0.5*I/csqrt(z-1.);
		double complex x4=0.5-0.5*I/csqrt(z-1.);
		
		double complex A=
		-104./243.*log(mb*mb/mu_b/mu_b)+4.*shat/27./(1.-shat)*(Li2(shat)+log(shat)*log(1.-shat))
		+1./729./(1.-shat)/(1.-shat)*(6.*shat*(29.-47.*shat)*log(shat)+785.-1600.*shat+833.*shat*shat+6.*pi*I*(20.-49.*shat+47.*shat*shat))
		-2./243./pow(1.-shat,3.)*(2.*csqrt(z-1.)*(-4.+9.*shat-15.*shat*shat+4.*shat*shat*shat)*(pi/2.-catan(cabs(z-1.)))+9.*shat*shat*shat*log(shat)*log(shat)+18.*pi*I*shat*(1.-2.*shat)*log(shat))
		+2.*shat/243./pow(1.-shat,4.)*(36.*cpow(pi/2.-catan(csqrt(z-1.)),2.)+pi*pi*(-4.+9.*shat-9.*shat*shat+3.*shat*shat*shat));
		
		double complex B=
		8./243./shat*((4.-34.*shat-17.*pi*I*shat)*log(mb*mb/mu_b/mu_b)+8.*shat*pow(log(mb*mb/mu_b/mu_b),2.)+17.*shat*log(shat)*log(mb*mb/mu_b/mu_b))
		+(2.+shat)*csqrt(z-1.)/729./shat*(-48.*log(mb*mb/mu_b/mu_b)*(pi/2.-catan(csqrt(z-1.)))-18.*pi*clog(z-1.)+3.*I*clog(z-1.)*clog(z-1.)
		-24.*I*CLi2(-x2/x1)-5.*pi*pi*I+6.*I*(-9.*clog(x1)*clog(x1)+clog(x2)*clog(x2)-2.*clog(x4)*clog(x4)+6.*clog(x1)*clog(x2)-4.*clog(x1)*clog(x3)+8.*clog(x1)*clog(x4))
		-12.*pi*(2.*clog(x1)+clog(x3)+clog(x4)))
		-2./243./shat/(1.-shat)*(4.*shat*(-8.+17.*shat)*(Li2(shat)+log(shat)*log(1.-shat))
		+3.*(2.+shat)*(3.-shat)*clog(x2/x1)*clog(x2/x1)+12.*pi*(-6.-shat+shat*shat)*(pi/2.-catan(csqrt(z-1.))))
		+2./2187./shat/(1.-shat)/(1.-shat)*(-18.*shat*(120.-211.*shat+73.*shat*shat)*log(shat)-288.-8.*shat+934.*shat*shat-692.*shat*shat*shat+18.*pi*I*shat*(82.-173.*shat+73.*shat*shat))
		-4./243./shat/pow(1.-shat,3.)*(-2.*csqrt(z-1.)*(4.-3.*shat-18.*shat*shat+16.*shat*shat*shat-5.*pow(shat,4.))*(pi/2.-catan(csqrt(z-1.)))-9.*shat*shat*shat*log(shat)*log(shat)+2.*pi*I*shat*(8.-33.*shat+51.*shat*shat-17.*shat*shat*shat)*log(shat))
		+2./729./shat/pow(1.-shat,4.)*(72.*(3.-8.*shat+2.*shat*shat)*cpow(pi/2.-catan(csqrt(z-1.)),2.)-pi*pi*(54.-53.*shat-286.*shat*shat+612.*pow(shat,3.)-446.*pow(shat,4.)+113.*pow(shat,5.)));
		
		double complex C=-16./81.*log(q2/mu_b/mu_b)+428./243.-64./27.*zeta3+16./81.*pi*I;
		
		double kappa=1.-2.*alphas_mub/3./pi*log(mu_b/mb);
		
		double complex C9eff=Cmub[9]
		+h_bkll(q2,0.,mu_b)*(4./3.*Cmub[1]+Cmub[2]+11./2.*Cmub[3]-2./3.*Cmub[4]+52.*Cmub[5]-32./3.*Cmub[6])
		-1./2.*h_bkll(q2,mb,mu_b)*(7.*Cmub[3]+4./3.*Cmub[4]+76.*Cmub[5]+64./3.*Cmub[6])
	+4./3.*(Cmub[3]+16./3.*Cmub[5]+16./9.*Cmub[6])
	+alphas_mub/4./pi*(Cmub[1]*(B+4.*C)-3.*Cmub[2]*(2.*B-C)-Cmub[8]*F89_bsll(shat))
	+8.*mc*mc/q2*(4./9.*Cmub[1]+1./3.*Cmub[2]+2.*Cmub[3]+20.*Cmub[5]);
		double complex C7eff=Cmub[7]
		+alphas_mub/4./pi*((Cmub[1]-6.*Cmub[2])*A-Cmub[8]*F87_bsll(shat,log(mu_b/mb)));
		
		double lambda_hat=1.+shat*shat+pow(param->m_Kstar/param->m_Bd,4.)-2.*(shat+shat*param->m_Kstar*param->m_Kstar/param->m_Bd/param->m_Bd+param->m_Kstar*param->m_Kstar/param->m_Bd/param->m_Bd);
		
		double N=sqrt(param->Gfermi*param->Gfermi*alpha_em*alpha_em/3./1024./pow(pi,5.)*param->m_Bd*pow(cabs(param->Vtb*conj(param->Vts)),2.)*shat*sqrt(lambda_hat));
		
		double f_perp=N*param->m_Bd*sqrt(2.*lambda_hat)/(1.+param->m_Kstar/param->m_Bd)*V;
		
		double f_par=N*param->m_Bd*sqrt(2.)*(1.+param->m_Kstar/param->m_Bd)*A1;
		
		double f_0=N*param->m_Bd*((1.-shat-param->m_Kstar*param->m_Kstar/param->m_Bd/param->m_Bd)*pow(1.+param->m_Kstar/param->m_Bd,2.)*A1-lambda_hat*A2)/(2.*param->m_Kstar/param->m_Bd*(1.+param->m_Kstar/param->m_Bd)*sqrt(shat));
		
		ALperp=I*((C9eff-Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_perp;		
		ARperp=I*((C9eff+Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_perp;
		
		ALpar=-I*((C9eff-Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_par;		
		ARpar=-I*((C9eff+Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_par;
		
		AL0=-I*((C9eff-Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_0;		
		AR0=-I*((C9eff+Cmub[10])+kappa*2.*mb/param->m_Bd/shat*C7eff)*f_0;

		AS=At=0.;
	}
			
	double A02=AL0*conj(AL0)+AR0*conj(AR0);
	double Apar2=ALpar*conj(ALpar)+ARpar*conj(ARpar);
	double Aperp2=ALperp*conj(ALperp)+ARperp*conj(ARperp);
	
	
	double J1s=0.25*(2.+beta_l*beta_l)*(Aperp2 + Apar2) + 4.*param->mass_mu*param->mass_mu/q2*creal(ALperp*conj(ARperp)+ALpar*conj(ARpar));

	double J1c=A02 + 4.*param->mass_mu*param->mass_mu/q2*(At*conj(At)+2.*creal(AL0*conj(AR0)))+beta_l*beta_l*AS*conj(AS);

	double J2s=0.25*beta_l*beta_l*(Aperp2+Apar2);

	double J2c=-beta_l*beta_l*A02;
	
	double J3=0.5*beta_l*beta_l*(Aperp2-Apar2);
	
	double J4=1./sqrt(2.)*beta_l*beta_l*creal(AL0*conj(ALpar)+AR0*conj(ARpar));
	
	double J5=sqrt(2.)*beta_l*(creal(AL0*conj(ALperp)-AR0*conj(ARperp))-param->mass_mu/sqrt(q2)*(creal(ALpar*conj(AS)+ARpar*conj(AS))));
	
	double J6s=2.*beta_l*creal(ALpar*conj(ALperp)-ARpar*conj(ARperp));

	double J6c=4.*beta_l*param->mass_mu/sqrt(q2)*creal(AL0*conj(AS)+AR0*conj(AS));

	double J7=sqrt(2.)*beta_l*(cimag(AL0*conj(ALpar)-AR0*conj(ARpar))+param->mass_mu/sqrt(q2)*(cimag(ALperp*conj(AS)+ARperp*conj(AS))));

	double J8=1./sqrt(2.)*beta_l*beta_l*cimag(AL0*conj(ALperp)+AR0*conj(ARperp));
	
	double J9=beta_l*beta_l*cimag(ALpar*conj(ALperp)+ARpar*conj(ARperp));
	
	double dGamma_BKstarmumu_dq2=3./4.*(2.*J1s+J1c-(2.*J2s+J2c)/3.);

	double AFB[3],FL[3],FT[3],AT1[3],AT2[3],AT3[3],AT4[3],AT5[3],HT1[3],HT2[3],HT3[3],alpha_Kstar[3],AIm[3],P2[3],P3[3],P6[3];

	AFB[0]=-3./8.*(2.*J6s+J6c)/dGamma_BKstarmumu_dq2;
	AFB[1]=-3./8.*(2.*J6s+J6c);
	AFB[2]=dGamma_BKstarmumu_dq2;

	FL[0]=(3.*J1c-J2c)/4./dGamma_BKstarmumu_dq2;
	FL[1]=(3.*J1c-J2c)/4.;
	FL[2]=dGamma_BKstarmumu_dq2;

	FT[0]=4.*J2s/dGamma_BKstarmumu_dq2;
	FT[1]=4.*J2s;
	FT[2]=dGamma_BKstarmumu_dq2;

	AT1[0]=-2.*creal(ALpar*conj(ALperp)+ARpar*conj(ARperp))/(Apar2+Aperp2);
	AT1[1]=-2.*creal(ALpar*conj(ALperp)+ARpar*conj(ARperp));
	AT1[2]=Apar2+Aperp2;

	AT2[0]=J3/2./J2s;
	AT2[1]=J3;
	AT2[2]=2.*J2s;

	AT3[0]=sqrt(fabs((4.*J4*J4+beta_l*beta_l*J7*J7)/(-2.*J2c*(2.*J2s+J3))));
	AT3[1]=sqrt(fabs(4.*J4*J4+beta_l*beta_l*J7*J7));
	AT3[2]=sqrt(fabs(-2.*J2c*(2.*J2s+J3)));
	
	AT4[0]=sqrt((beta_l*beta_l*J5*J5+4.*J8*J8)/(4.*J4*J4+beta_l*beta_l*J7*J7));
	AT4[1]=sqrt(beta_l*beta_l*J5*J5+4.*J8*J8);
	AT4[2]=sqrt(4.*J4*J4+beta_l*beta_l*J7*J7);
	
	AT5[0]=cabs(ALperp*conj(ARpar)+ALpar*conj(ARperp))/(Apar2+Aperp2);
	AT5[1]=cabs(ALperp*conj(ARpar)+ALpar*conj(ARperp));
	AT5[2]=Apar2+Aperp2;
	
	HT1[0]=sqrt(2.)*J4/sqrt(fabs(-J2c*(2.*J2s-J3)));
	HT1[1]=sqrt(2.)*J4;
	HT1[2]=sqrt(fabs(-J2c*(2.*J2s-J3)));
	
	HT2[0]=beta_l*J5/sqrt(fabs(-2.*J2c*(2.*J2s+J3)));
	HT2[1]=beta_l*J5;
	HT2[2]=sqrt(fabs(-2.*J2c*(2.*J2s+J3)));
	
	HT3[0]=beta_l*J6s/2./sqrt(fabs(4.*J2s*J2s-J3*J3));
	HT3[1]=beta_l*J6s;
	HT3[2]=2.*sqrt(fabs(4.*J2s*J2s-J3*J3));

	alpha_Kstar[0]=-(2.*J2s+J2c)/2./J2s;
	alpha_Kstar[1]=-(2.*J2s+J2c);
	alpha_Kstar[2]=2.*J2s;

	AIm[0]=J9/dGamma_BKstarmumu_dq2;
	AIm[1]=J9;
	AIm[2]=dGamma_BKstarmumu_dq2;

	P2[0]=beta_l*J6s/8./J2s;
	P2[1]=beta_l*J6s;
	P2[2]=8.*J2s;

	P3[0]=-J9/4./J2s;
	P3[1]=-J9;
	P3[2]=4.*J2s;

	P6[0]=-beta_l*J7/sqrt(fabs(-2.*J2c*(2.*J2s+J3)));
	P6[1]=-beta_l*J7;
	P6[2]=sqrt(fabs(-2.*J2c*(2.*J2s+J3)));
	
	for(je=0;je<=Nobs_BKsll;je++) for(ie=0;ie<=2;ie++) obs[je][ie]=0.;
	
	for(ie=0;ie<=2;ie++)
	{
		obs[1][ie]=AFB[ie];
		obs[2][ie]=FL[ie];
		obs[3][ie]=FT[ie];
		obs[4][ie]=AT1[ie];
		obs[5][ie]=AT2[ie];
		obs[6][ie]=AT3[ie];
		obs[7][ie]=AT4[ie];
		obs[8][ie]=AT5[ie];
		obs[9][ie]=HT1[ie];
		obs[10][ie]=HT2[ie];
		obs[11][ie]=HT3[ie];
		obs[12][ie]=alpha_Kstar[ie];
		obs[13][ie]=AIm[ie];
		obs[14][ie]=P2[ie];
		obs[15][ie]=P3[ie];
		obs[16][ie]=P6[ie];
	}

	return dGamma_BKstarmumu_dq2;
}

/*----------------------------------------------------------------------*/

double dGamma_BKstarmumu_dq2_calculator(double q2, double obs[][3], char name[])
/* "container" function scanning the SLHA file "name" and calculating dGamma/dq2(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	return dGamma_BKstarmumu_dq2(q2,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarmumu(double smin, double smax, double obs[], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	int ie,je;
	int nmax=100.;
	double Gamma=0.;
	double s;
	
	double obs_num[Nobs_BKsll+1],obs_den[Nobs_BKsll+1];
	for(je=0;je<=Nobs_BKsll;je++) obs_num[je]=obs_den[je]=0.;

	obs[0]=0.; /* zero AFB */
	obs[1]=0.; /* integrated AFB */
	obs[2]=0.; /* integrated FL */
	obs[3]=0.; /* integrated FT */
	obs[4]=0.; /* integrated AT1 */
	obs[5]=0.; /* integrated AT2 */
	obs[6]=0.; /* integrated AT3 */
	obs[7]=0.; /* integrated AT4 */
	obs[8]=0.; /* integrated AT5 */
	obs[9]=0.; /* integrated HT1 */
	obs[10]=0.; /* integrated HT2 */
	obs[11]=0.; /* integrated HT3 */
	obs[12]=0.; /* integrated alpha */
	obs[13]=0.; /* integrated AIm */
	obs[14]=0.; /* integrated P2 */
	obs[15]=0.; /* integrated P3 */
	obs[16]=0.; /* integrated P6 */
	
	double dobs[Nobs_BKsll+1][3],dAFBtmp;
	double s0m,s0p,s0;
		
	dAFBtmp=0.;
	s0=0.;
	s0p=1.;
	
	for(ie=1;ie<=nmax;ie++)
	{
		dAFBtmp=dobs[1][0];
		s0m=s0p;
		s=smin+(smax-smin)*ie/nmax;
		s0p=s;	
		Gamma+=dGamma_BKstarmumu_dq2(s,dobs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
		
		for(je=1;je<=Nobs_BKsll;je++) 
		{
			obs_num[je]+=dobs[je][1];
			obs_den[je]+=dobs[je][2];
		}
		
		if((dAFBtmp/dobs[1][0]<0.)&&(ie>1)) s0=(dobs[1][0]*s0m-dAFBtmp*s0p)/(dobs[1][0]-dAFBtmp);
	}
	Gamma*=(smax-smin)/nmax;
	obs[0]=s0;
	for(je=1;je<=Nobs_BKsll;je++) obs[je]=obs_num[je]/obs_den[je];

	return param->life_Bd/hbar*Gamma;
}

/*----------------------------------------------------------------------*/

double BRBKstarmumu_lowq2(double obs[], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKstarmumu(1.,6.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarmumu_highq2(double obs[], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b)
{
	return BRBKstarmumu(14.18,16.,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarmumu_lowq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRBKstarmumu_highq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double A_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[1];
}

/*----------------------------------------------------------------------*/

double A_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;

	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[1];
}

/*----------------------------------------------------------------------*/

double FL_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[2];
}

/*----------------------------------------------------------------------*/

double FL_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;

	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[2];
}

/*----------------------------------------------------------------------*/

double FT_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[3];
}

/*----------------------------------------------------------------------*/

double FT_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[3];
}

/*----------------------------------------------------------------------*/

double AT1_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[4];
}

/*----------------------------------------------------------------------*/

double AT1_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[4];
}

/*----------------------------------------------------------------------*/

double AT2_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[5];
}

/*----------------------------------------------------------------------*/

double AT2_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[5];
}

/*----------------------------------------------------------------------*/

double AT3_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[6];
}

/*----------------------------------------------------------------------*/

double AT3_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[6];
}

/*----------------------------------------------------------------------*/

double AT4_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[7];
}

/*----------------------------------------------------------------------*/

double AT4_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[7];
}

/*----------------------------------------------------------------------*/

double AT5_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[8];
}

/*----------------------------------------------------------------------*/

double AT5_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[8];
}

/*----------------------------------------------------------------------*/

double HT1_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[9];
}

/*----------------------------------------------------------------------*/

double HT1_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[9];
}

/*----------------------------------------------------------------------*/

double HT2_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[10];
}

/*----------------------------------------------------------------------*/

double HT2_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[10];
}

/*----------------------------------------------------------------------*/

double HT3_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[11];
}

/*----------------------------------------------------------------------*/

double HT3_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[11];
}

/*----------------------------------------------------------------------*/

double alpha_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[12];
}

/*----------------------------------------------------------------------*/

double alpha_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[12];
}

/*----------------------------------------------------------------------*/

double AIm_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[13];
}

/*----------------------------------------------------------------------*/

double AIm_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[13];
}

/*----------------------------------------------------------------------*/

double P2_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[14];
}

/*----------------------------------------------------------------------*/

double P2_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[14];
}
/*----------------------------------------------------------------------*/

double P3_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[15];
}

/*----------------------------------------------------------------------*/

double P3_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[15];
}
/*----------------------------------------------------------------------*/

double P6_BKstarmumu_lowq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[16];
}

/*----------------------------------------------------------------------*/

double P6_BKstarmumu_highq2_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[16];
}

/*----------------------------------------------------------------------*/

double BRobs_BKstarmumu_lowq2_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) and all the other observables */

	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarmumu_lowq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double BRobs_BKstarmumu_highq2_calculator(char name[], double obs[])
{
/* "container" function scanning the SLHA file "name" and calculating BR(B->Kstar mu+ mu-) and all the other observables */
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	return BRBKstarmumu_highq2(obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double A_BKstarmumu_zero_calculator(char name[])
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	struct parameters param;
	double obs[Nobs_BKsll+1];
	
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;

	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b,&param);

	double smin=pow(2.*param.mass_mu,2.);
	double smax=pow(param.m_Bd-param.m_Kstar,2.)*0.999; 
	
	BRBKstarmumu(smin,smax,obs,C0b,C1b,C2b,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b);
		
	return obs[0];
}

/*----------------------------------------------------------------------*/

double dAI_BKstarmumu_dq2(double q2, double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b)
{
	double mc=mc_pole(param);

	int ie;
	double shat=q2/param->m_B/param->m_B;
	
	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);

	double mu_f=sqrt(mu_b*0.5);
	double alphas_muf=alphas_running(mu_f,param->mass_top_pole,param->mass_b_pole,param);
	double eta=alphas_muf/alphas_running(1.,param->mass_top_pole,param->mass_b_pole,param);

	double alphas_mbpole=alphas_running(param->mass_b_pole,param->mass_top_pole,param->mass_b_pole,param);
	double mb=param->mass_b_pole-4.*alphas_mbpole*mu_f/3./pi; /* mb(PS)_muf */

	double Cmub[11];
	for(ie=1;ie<=10;ie++) Cmub[ie]=C0b[ie]+alphas_mub/4./pi*C1b[ie]+pow(alphas_mub/4./pi,2.)*C2b[ie];
	
	double E_Kstar=(param->m_B*param->m_B+param->m_Kstar*param->m_Kstar-q2)/2./param->m_B;

	int nf=5;
	double f_K_perp=param->f_K_perp;
	f_K_perp*=pow(eta,4./3./(11.-2./3.*nf));

	double f_K_par=param->f_K_par;

	double V=0.923/(1.-q2/5.32/5.32)-0.511/(1.-q2/49.4);
	double A1=0.290/(1.-q2/40.38);
	double A2=-0.084/(1.-q2/52.)+0.342/(1.-q2/52.)/(1.-q2/52.);
	double xi_par=(param->m_B+param->m_Kstar)/2./E_Kstar*A1-(param->m_B-param->m_Kstar)/param->m_B*A2;
	double xi_perp=param->m_B/(param->m_B+param->m_Kstar)*V;

	double C7eff=Cmub[7];
	double C8eff=Cmub[8];
	double C9=Cmub[9];
	double C10=Cmub[10];
		
	double C1bar=Cmub[1]/2.;
	double C2bar=Cmub[2]-Cmub[1]/6.;
	double C3bar=Cmub[3]-Cmub[4]/6.+16.*Cmub[5]-8./3.*Cmub[6];
	double C4bar=Cmub[4]/2.+8.*Cmub[6];
	double C5bar=Cmub[3]-Cmub[4]/6.+4.*Cmub[5]-2./3.*Cmub[6];
	double C6bar=Cmub[4]/2.+2.*Cmub[6];
	
	double complex Y=4./3.*Cmub[3]+64./9.*Cmub[5]+64./27.*Cmub[6]
	+h_bkll(q2,mc,mu_b)*(4./3.*Cmub[1]+Cmub[2]+6.*Cmub[3]+60.*Cmub[5])
	+h_bkll(q2,mb,mu_b)*(-7./2.*Cmub[3]-2./3.*Cmub[4]-38.*Cmub[5]-32./3.*Cmub[6])
	+h_bkll(q2,0.,mu_b)*(-1./2.*Cmub[3]-2./3.*Cmub[4]-8.*Cmub[5]-32./3.*Cmub[6]);

	double complex C90perp=C9+Y+2.*mb*param->m_B/q2*C7eff;
	double complex C90par=C9+Y+2.*mb/param->m_B*C7eff;

	double a1perp=param->a1perp;
	double a2perp=param->a2perp;
	double a1par=param->a1par;
	double a2par=param->a2par;

	a1perp*=pow(eta,4./(11.-2./3.*nf));
	a2perp*=pow(eta,4./3.*(1.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

	a1par*=pow(eta,4./3.*(1.-1./3.+2.)/(11.-2./3.*nf));
	a2par*=pow(eta,4./3.*(1.-1./6.+4.*(1./2.+1./3.))/(11.-2./3.*nf));

	double u;

	double lambda_Bp=param->lambda_Bp;
	lambda_Bp /= 1.+alphas_muf/3./pi*log(pow(mu_b,2.))*(1.-2.*1.4);

	double omega0=2.*(param->m_Bd-mb)/3.;
	double complex lambda_Bm=1./(exp(-q2/param->m_Bd/omega0)/omega0*(-Ei(q2/param->m_Bd/omega0)+I*pi));

	double complex integ1,integ2,integ3;
	integ1=integ2=integ3=0.;
	double complex Fperp=0.;
	double complex Xperp=0.;
	double complex Fpar=0.;
	double complex FV;
	double x;
	
	double zeta3A=param->zeta3A;
	double zeta3V=param->zeta3V;
	double wA10=param->wA10;
	double deltatp=param->deltatp;
	double deltatm=param->deltatm;
	
	int n1=50;
	for(ie=1;ie<=n1-1;ie++)
	{
		u=(double)ie/n1;
		x=(1.-u)*param->m_B*param->m_B+u*q2;
		FV=3./4.*(h_bkll(x,mc,mu_b)*(C2bar+C4bar+C6bar)+h_bkll(x,mb,mu_b)*(C3bar+C4bar+C6bar)+h_bkll(x,0.,mu_b)*(C3bar+3.*C4bar+3.*C6bar)-8./27.*(C3bar-C5bar-15.*C6bar));
		
		integ1+=phi_Kstar(u,a1perp,a2perp)/((1.-u)+u*shat)*FV/n1;
		
		integ2+=((3./4.*(1.+pow(2.*u-1.,2.))+a1par*3./2.*pow(2.*u-1.,3.)+(3./7.*a2par+5.*zeta3A)*(3.*pow(2.*u-1.,2.)-1.)+(9./122.*a2par+105./16.*zeta3V-15./64.*zeta3A*wA10)*(3.-30.*pow(2.*u-1.,2.)+35.*pow(2.*u-1.,4.))+3.*deltatp+3.*deltatm*(2.*u-1.))-1./4.*(6.*(1.-2.*u)*(1.+a1par*(2.*u-1.)+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(5.*pow(2.*u-1.,2.)-1.))+6.*u*(1.-u)*(2.*a1par*u+(a2par/4.+5./3.*zeta3A*(1.-3./16.*wA10)+35./4.*zeta3V)*(20.*u*(2.*u-1.)))+18.*deltatp*(1.-2.*u)-12.*deltatm))*FV/n1;
	
		integ3+=phi_Kstar(u,a1par,a2par)*FV/n1;
		
		Fperp+=phi_Kstar(u,a1perp,a2perp)/((1.-u)+u*shat)/3./n1;
		Xperp+=(u<=1.-0.5/param->m_B)*phi_Kstar(u,a1perp,a2perp)/pow((1.-u)+u*shat,2.)/3.*(0.5/param->m_B)/n1;
		Fpar+=2.*phi_Kstar(u,a1par,a2par)/((1.-u)+u*shat)/n1;
	}
	double rho=0.;
	double phi=0.;
	Xperp=Fperp+(1.+rho*(cos(phi)+I*sin(phi)))*Xperp;
	
	double complex lambda_u=conj(param->Vus)*param->Vub;
	double complex lambda_t=conj(param->Vts)*param->Vtb;

	double complex K1upar_a=-lambda_u/lambda_t*(C1bar/3.+C2bar)+(C4bar+C3bar/3.);
	double complex K1dpar_a=C4bar+C3bar/3.;
	
	double complex K1perp_a=-(C6bar+C5bar/3.)*Fperp;
	
	double complex K2uperp_a=-lambda_u/lambda_t*(C1bar/3.+C2bar)+(C4bar+C3bar/3.);
	double complex K2dperp_a=C4bar+C3bar/3.;
	
	double complex K1perp_b=C8eff*mb/param->m_B*4./9.*alphas_mub/4./pi*Xperp;
	double complex K2perp_b=0.;
	double complex K1par_b=-C8eff*mb/param->m_B*4./9.*alphas_mub/4./pi*Fpar;
			
	double complex K1perp_c=4./9.*alphas_mub/4./pi*2./3.*integ1;
	double complex K2perp_c=-4./9.*alphas_mub/4./pi*1./2.*integ2;
	
	double complex K1par_c=-4./9.*alphas_mub/4./pi*2.*integ3;
	
	double complex K1perp=K1perp_a+K1perp_b+K1perp_c;
	double complex K2uperp=K2uperp_a+K2perp_c;
	double complex K2dperp=K2dperp_a+K2perp_c;
	double complex K1upar=K1upar_a+K1par_b+K1par_c;
	double complex K1dpar=K1dpar_a+K1par_b+K1par_c;
	
	double eu=2./3.;
	double ed=-1./3.;
	
	double complex bd_perp=24.*pi*pi*param->m_B*param->f_B*ed/q2/xi_perp/C90perp*(f_K_perp/param->m_B*K1perp+f_K_par*param->m_Kstar/6./lambda_Bp/param->m_B*K2dperp/(1.-q2/param->m_B/param->m_B));
	
	double complex bu_perp=24.*pi*pi*param->m_B*param->f_B*eu/q2/xi_perp/C90perp*(f_K_perp/param->m_B*K1perp+f_K_par*param->m_Kstar/6./lambda_Bp/param->m_B*K2uperp/(1.-q2/param->m_B/param->m_B));
	
	double complex bd_par=24.*pi*pi*param->f_B*ed*param->m_Kstar/param->m_B/E_Kstar/xi_par/C90par*(f_K_par/3./lambda_Bm*K1dpar);
	
	double complex bu_par=24.*pi*pi*param->f_B*eu*param->m_Kstar/param->m_B/E_Kstar/xi_par/C90par*(f_K_par/3./lambda_Bm*K1upar);	
	
	double dAI_dq2=creal(bd_perp-bu_perp)*pow(cabs(C90perp),2.)/(pow(cabs(C90perp),2.)+C10*C10)*(1.+0.25/q2*pow(E_Kstar*param->m_B/param->m_Kstar*xi_par/xi_perp*cabs(C90par/C90perp),2.)*creal(bd_par-bu_par)/creal(bd_perp-bu_perp))/(1.+0.25/q2*pow(E_Kstar*param->m_B/param->m_Kstar*xi_par/xi_perp,2.)*(pow(cabs(C90par),2.)+C10*C10)/(pow(cabs(C90perp),2.)+C10*C10));
		
	return dAI_dq2;
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu(double smin, double smax, double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b)
{
	int ie,je;
	int nmax=100.;
	double AI=0.;
	double s;
		
	for(ie=1;ie<=nmax;ie++)
	{
		s=smin+(smax-smin)*ie/nmax;
		AI+=dAI_BKstarmumu_dq2(s,C0b,C1b,C2b,param,mu_b);
	}
	AI*=(smax-smin)/nmax;

	return AI;
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_lowq2(double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b)
{
	return AI_BKstarmumu(1.,6.,C0b,C1b,C2b,param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_highq2(double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b)
{
	return AI_BKstarmumu(14.18,16.,C0b,C1b,C2b,param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_zero(double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b)
{
	double smin=pow(2.*param->mass_mu,2.);
	double smax=pow(param->m_B-param->m_Kstar,2.)*0.999; 

	double stemp;

	while(fabs(1.-smin/smax)>1.e-4)
	{
		stemp=(smax+smin)/2.;
		
		if(dAI_BKstarmumu_dq2(stemp,C0b,C1b,C2b,param,mu_b)>0.) smin=stemp; else smax=stemp;
	}
	return stemp;
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_lowq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating AI(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);

	return AI_BKstarmumu_lowq2(C0b,C1b,C2b,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_highq2_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating AI(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);

	return AI_BKstarmumu_highq2(C0b,C1b,C2b,&param,mu_b);
}

/*----------------------------------------------------------------------*/

double AI_BKstarmumu_zero_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the zero of AI(B->Kstar mu+ mu-) */
{
	double C0b[11],C1b[11],C2b[11],C0w[11],C1w[11],C2w[11];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;

	double mu_W=2.*param.mass_W;
	double mu_b=param.mass_b_pole;
				
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);

	return AI_BKstarmumu_zero(C0b,C1b,C2b,&param,mu_b);
}

/*----------------------------------------------------------------------*/
