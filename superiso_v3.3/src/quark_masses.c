#include "include.h"


double running_mass(double quark_mass, double Qinit, double Qfin,  double mtop, double mbot, struct parameters* param)
/* computes the running quark mass at the energy Qfin, from a given running quark mass quark_mass at energy Qinit */
{

	double alphas_Qinit,alphas_Qfin,running_mass;
	double beta0,beta1,beta2,gamma0,gamma1,gamma2;
	int nf;

	alphas_Qinit=alphas_running(Qinit,mtop,mbot,param);
	
	double R_Qinit,R_Qfin;
		
	if(Qinit <= mbot) 
	/* 4 active flavors at Qinit */
	{
		nf=4;
	
		beta0=11.-2./3.*nf;
		beta1=51.-19./3.*nf;
		beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
		gamma0=2.;
		gamma1=101./12.-5./18.*nf;
		gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

		R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));

		if(Qfin <= mbot) 
		{
		/* if Qinit and Qfin are in the same range, just calculate R_Qfin */
		/* 4 active flavors at Qinit and Qfin */
			nf=4;

			alphas_Qfin=alphas_running(Qfin,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
		}
		else if(Qfin <= mtop)
		{
		/* if Qinit and Qfin are NOT in the same range, evolve first the running mass towards the range limit(s), then from the limit(s) towards Qfin */
		/* 4 active flavors at Qinit, 5 active flavors at Qfin */
			nf=4;

			alphas_Qfin=alphas_running(mbot,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;


			alphas_Qinit=alphas_Qfin;

			nf=5;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=5;

			alphas_Qfin=alphas_running(Qfin,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}
		else
		/* 4 active flavors at Qinit, 6 active flavors at Qfin */
		{
			nf=4;

			alphas_Qfin=alphas_running(mbot,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;


			alphas_Qinit=alphas_Qfin;

			nf=5;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=5;

			alphas_Qfin=alphas_running(mtop,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
	
	
			alphas_Qinit=alphas_Qfin;

			nf=6;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=6;
	
			alphas_Qfin=alphas_running(Qfin,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}
	}
	else if(Qinit <= mtop) 
	/* 5 active flavors at Qinit */
	{
		nf=5;
	
		beta0=11.-2./3.*nf;
		beta1=51.-19./3.*nf;
		beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
		gamma0=2.;
		gamma1=101./12.-5./18.*nf;
		gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

		R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));

		if(Qfin <= mbot) 
		/* 5 active flavors at Qinit, 4 active flavors at Qfin */
		{
			nf=5;

			alphas_Qfin=alphas_running(mbot,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;

			alphas_Qinit=alphas_Qfin;

			nf=4;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=4;

			alphas_Qfin=alphas_running(Qfin,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}
		else if(Qfin <= mtop)
		/* 5 active flavors at Qinit and Qfin */
		{
			nf=5;

			alphas_Qfin=alphas_running(Qfin,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
		}
		else
		/* 5 active flavors at Qinit, 6 active flavors at Qfin */
		{
			nf=5;

			alphas_Qfin=alphas_running(mtop,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
	
			alphas_Qinit=alphas_Qfin;


			nf=6;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=6;
	
			alphas_Qfin=alphas_running(Qfin,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}	
	}
	else
	/* 6 active flavors at Qinit */
	{
		nf=6;
	
		beta0=11.-2./3.*nf;
		beta1=51.-19./3.*nf;
		beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
		gamma0=2.;
		gamma1=101./12.-5./18.*nf;
		gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

		R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));

		if(Qfin >= mtop) 
		/* 6 active flavors at Qinit and Qfin */
		{
			nf=6;

			alphas_Qfin=alphas_running(Qfin,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
		}
		else if(Qfin >= mbot)
		/* 6 active flavors at Qinit, 5 active flavors at Qfin */
		{
			nf=6;

			alphas_Qfin=alphas_running(mtop,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
	
			alphas_Qinit=alphas_Qfin;


			nf=5;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=5;

			alphas_Qfin=alphas_running(Qfin,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}
		else
		/* 6 active flavors at Qinit, 4 active flavors at Qfin */
		{
			nf=6;

			alphas_Qfin=alphas_running(mtop,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*quark_mass;
	
			alphas_Qinit=alphas_Qfin;


			nf=5;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			alphas_Qfin=alphas_running(mbot,mtop,mbot,param);

			nf=5;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
	
			alphas_Qinit=alphas_Qfin;


			nf=4;
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);

			R_Qinit = pow(beta0/2./pi*alphas_Qinit,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qinit/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qinit/pi,2.));
	
	
			nf=4;

			alphas_Qfin=alphas_running(Qfin,mtop,mbot,param);
	
			beta0=11.-2./3.*nf;
			beta1=51.-19./3.*nf;
			beta2=2857.-5033./9.*nf+325./27.*nf*nf;
	
			gamma0=2.;
			gamma1=101./12.-5./18.*nf;
			gamma2=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf-140./81.*nf*nf);
	
			R_Qfin = pow(beta0/2./pi*alphas_Qfin,2.*gamma0/beta0)*(1.+(2.*gamma1/beta0-beta1*gamma0/beta0/beta0)*alphas_Qfin/pi+1./2.*(pow(2.*gamma1/beta0-beta1*gamma0/beta0/beta0,2.)+2.*gamma2/beta0-beta1*gamma1/beta0/beta0-beta2*gamma0/16./beta0/beta0+beta1*beta1*gamma0/2./beta0/beta0/beta0)*pow(alphas_Qfin/pi,2.));

			running_mass=R_Qfin/R_Qinit*running_mass;
		}
	}
	
	return running_mass;
}

/*--------------------------------------------------------------------*/

double mb_pole(struct parameters* param)
/* computes the b pole mass */
{
	
	double alphas_mb=alphas_running(param->mass_b,param->mass_top_pole,param->mass_b,param);	

 	return param->mass_b*(1.+alphas_mb/pi*(4./3.+alphas_mb/pi*((13.4434-1.0414*4.+1.0414*4./3.*((param->mass_u+param->mass_d+param->mass_s+param->mass_c)/param->mass_b))
	/* +alphas_mb/pi*(190.595-4.*(26.655-4.*0.6527)) */
	)));
}

/*--------------------------------------------------------------------*/

double mc_pole(struct parameters* param)
/* computes the c pole mass */{
	
	double alphas_mc=alphas_running(param->mass_c,param->mass_top_pole,param->mass_b_pole,param);	

 	return param->mass_c*(1+alphas_mc/pi*(4./3.+alphas_mc/pi*((13.4434-1.0414*3.+1.0414*4./3.*((param->mass_u+param->mass_d+param->mass_s)/param->mass_c)))));
}

/*--------------------------------------------------------------------*/

double mb_1S(struct parameters* param)
/* computes the 1S b mass */
{
	double mb=param->mass_b_pole;
	double mu=param->mass_b/2.;
	double as=alphas_running(mu,param->mass_top_pole,param->mass_b_pole,param);
	double L=log(mu/(4./3.*as*mb));
	double beta0=11.-8./3.;
	double beta1=62.-32./3.;
	double a1=31./3.-40./9.;
	double a2=(4343./162.+4.*pi*pi-pow(pi,4.)/4.+22./3.*zeta3)*9.-(1798./81.+56./3.*zeta3)*6.
	-(55./3.-16.*zeta3)*8./3.+1600./81.;
	
	return mb*(1.-2./9.*pow(as,2.)
	/* -2./9.*pow(as,3.)/pi*(beta0*(L+1.)+a1/2.) */
	)
	-2./9.*pow(as,2.)
	/* -2./9.*pow(as,4.)/pi/pi*(beta0*beta0*(3./4.*L*L+L+zeta3/2.+pi*pi/24.+1./4.)
	+beta0*a1/2.*(3./2.*L+1.)+beta1/4.*(L+1.)+a1*a1/16.+a2/8.+(3.-1./36.)*4./3.*pi*pi) */
	;
}

/*--------------------------------------------------------------------*/

double mt_mt(struct parameters* param)
/* computes the top mass mt(mt) */
{
	
	double alphas_mtop=alphas_running(param->mass_top_pole,param->mass_top_pole,param->mass_b,param); 

 	return param->mass_top_pole*(1.-4./3.*alphas_mtop/pi);
}
