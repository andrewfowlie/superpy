#include "include.h"


double alphas_running(double Q, double mtop, double mbot, struct parameters* param)
/* computes the QCD coupling constant alphas at the energy Q, using the matching scales between the flavors mtop and mbot */
{
	double beta0,beta1,beta2,alphas_running,Lambda4,Lambda5,Lambda6,Lambda_min,Lambda_max,Lambda_moy,alphas_min,alphas_max,alphas_moy;
	int nf;

	double MZ=param->mass_Z;
	double alphas_MZ=param->alphas_MZ;

	nf=5;
	beta0 = 11.-2./3.*nf;
	beta1=51.-19./3.*nf;
	beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

	if(param->Lambda5==-1.) return -1.;

	if(param->Lambda5==0.)
	{
		Lambda_min=1.e-3;
		Lambda_max=1.;
		alphas_min=0.;

		while((fabs(1.-alphas_min/alphas_MZ)>=1.e-4)&&(fabs(1.-Lambda_min/Lambda_max)>1.e-5))
		{
			alphas_min=4.*pi/beta0/log(pow(MZ/Lambda_min,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(MZ/Lambda_min,2.)))/log(pow(MZ/Lambda_min,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(MZ/Lambda_min,2.)),2.)*(pow(log(log(pow(MZ/Lambda_min,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			alphas_max=4.*pi/beta0/log(pow(MZ/Lambda_max,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(MZ/Lambda_max,2.)))/log(pow(MZ/Lambda_max,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(MZ/Lambda_max,2.)),2.)*(pow(log(log(pow(MZ/Lambda_max,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			Lambda_moy=(Lambda_min+Lambda_max)/2.;

			alphas_moy=4.*pi/beta0/log(pow(MZ/Lambda_moy,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(MZ/Lambda_moy,2.)))/log(pow(MZ/Lambda_moy,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(MZ/Lambda_moy,2.)),2.)*(pow(log(log(pow(MZ/Lambda_moy,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			if((alphas_MZ>=alphas_min)&&(alphas_MZ<=alphas_moy))
				Lambda_max=Lambda_moy;
			else Lambda_min=Lambda_moy;
		}

		Lambda5=Lambda_min;
		param->Lambda5=Lambda5;
		
		if(fabs(1.-Lambda_min/Lambda_max)<=1.e-5)
		{
			param->Lambda5=-1.;
			return -1.;
		}
	}
	else Lambda5=param->Lambda5;

	if((Q<=mtop)&&(Q>=mbot))
	/* 5 active flavors */
	{
		alphas_running=4.*pi/beta0/log(pow(Q/Lambda5,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda5,2.)))/log(pow(Q/Lambda5,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda5,2.)),2.)*(pow(log(log(pow(Q/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
		return alphas_running;
	}
	else if((Q>mtop))
	/* 6 active flavors */
	{
		alphas_running=4.*pi/beta0/log(pow(mtop/Lambda5,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda5,2.)))/log(pow(mtop/Lambda5,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda5,2.)),2.)*(pow(log(log(pow(mtop/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

		nf=6;
		beta0 = 11.-2./3.*nf;
		beta1=51.-19./3.*nf;
		beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

		Lambda_min=1.e-3;
		Lambda_max=1.;
		alphas_min=0.;

		while((fabs(1.-alphas_min/alphas_running)>=1.e-4)&&(fabs(1.-Lambda_min/Lambda_max)>1.e-5))
		{
			alphas_min=4.*pi/beta0/log(pow(mtop/Lambda_min,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda_min,2.)))/log(pow(mtop/Lambda_min,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda_min,2.)),2.)*(pow(log(log(pow(mtop/Lambda_min,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			alphas_max=4.*pi/beta0/log(pow(mtop/Lambda_max,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda_max,2.)))/log(pow(mtop/Lambda_max,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda_max,2.)),2.)*(pow(log(log(pow(mtop/Lambda_max,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			Lambda_moy=(Lambda_min+Lambda_max)/2.;
			alphas_moy=4.*pi/beta0/log(pow(mtop/Lambda_moy,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda_moy,2.)))/log(pow(mtop/Lambda_moy,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda_moy,2.)),2.)*(pow(log(log(pow(mtop/Lambda_moy,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			if((alphas_running>=alphas_min)&&(alphas_running<=alphas_moy))
				Lambda_max=Lambda_moy;
			else Lambda_min=Lambda_moy;
		}

		Lambda6=Lambda_min;
		
		if(fabs(1.-Lambda_min/Lambda_max)<=1.e-5)
		{
			param->Lambda5=-1.;
			return -1.;
		}
	alphas_running=4.*pi/beta0/log(pow(Q/Lambda6,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda6,2.)))/log(pow(Q/Lambda6,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda6,2.)),2.)*(pow(log(log(pow(Q/Lambda6,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

		return alphas_running;
	}
	else
	/* 4 active flavors */
	{
		alphas_running=4.*pi/beta0/log(pow(mbot/Lambda5,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/Lambda5,2.)))/log(pow(mbot/Lambda5,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/Lambda5,2.)),2.)*(pow(log(log(pow(mbot/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

		nf=4;
		beta0 = 11.-2./3.*nf;
		beta1=51.-19./3.*nf;
		beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

		Lambda_min=1.e-3;
		Lambda_max=1.;
		alphas_min=0.;

		while((fabs(1.-alphas_min/alphas_running)>=1.e-4)&&(fabs(1.-Lambda_min/Lambda_max)>1.e-5))
		{
			alphas_min=4.*pi/beta0/log(pow(mbot/Lambda_min,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/Lambda_min,2.)))/log(pow(mbot/Lambda_min,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/Lambda_min,2.)),2.)*(pow(log(log(pow(mbot/Lambda_min,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			alphas_max=4.*pi/beta0/log(pow(mbot/Lambda_max,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/Lambda_max,2.)))/log(pow(mbot/Lambda_max,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/Lambda_max,2.)),2.)*(pow(log(log(pow(mbot/Lambda_max,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			Lambda_moy=(Lambda_min+Lambda_max)/2.;
			alphas_moy=4.*pi/beta0/log(pow(mbot/Lambda_moy,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mbot/Lambda_moy,2.)))/log(pow(mbot/Lambda_moy,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mbot/Lambda_moy,2.)),2.)*(pow(log(log(pow(mbot/Lambda_moy,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

			if((alphas_running>=alphas_min)&&(alphas_running<=alphas_moy))
				Lambda_max=Lambda_moy;
			else Lambda_min=Lambda_moy;
		}
		Lambda4=Lambda_min;
		
		if(fabs(1.-Lambda_min/Lambda_max)<=1.e-5)
		{
			param->Lambda5=-1.;
			return -1.;
		}
	alphas_running=4.*pi/beta0/log(pow(Q/Lambda4,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda4,2.)))/log(pow(Q/Lambda4,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda4,2.)),2.)*(pow(log(log(pow(Q/Lambda4,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

		return alphas_running;
	}
}
