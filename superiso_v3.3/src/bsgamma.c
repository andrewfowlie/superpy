#include "include.h"


double phi77(double delta)
{
	return -2./3.*log(delta)*log(delta)-7./3.*log(delta)-31./9.+10./3.*delta+1./3.*delta*delta
	-2./9.*delta*delta*delta+1./3.*delta*(delta-4.)*log(delta);
}

/*---------------------------------------------------------------------*/

double phi78(double delta)
{
	return 8./9.*(Li2(1.-delta)-pi*pi/6.-delta*log(delta)+9.*delta/4.-delta*delta/4.+delta*delta*delta/12.);
}

/*---------------------------------------------------------------------*/

double phi88(double delta, double b)
{
	return (4.*Li2(1.-delta)-2.*pi*pi/3.+8.*log(1.-delta)-delta*(2.+delta)*log(delta)+
	+7.*delta+3.*delta*delta-2.*delta*delta*delta/3.-2.*(2.*delta+delta*delta+4.*log(1.-delta))*log(b))/27.;
}

/*---------------------------------------------------------------------*/

double G1(double t)
{
/* G1(t)=|G(t)/t+1/2|^2 */

	if (t==0.) return 0.;
	else if (t<4.) return pow(-2.*pow(atan(sqrt(t/(4.-t))),2.)/t+0.5,2.);
	else return pow((-pi*pi/2.+2.*pow(log((sqrt(t)+sqrt(t-4.))/2.),2.))/t+0.5,2.)+pow(2.*pi*log((sqrt(t)+sqrt(t-4.))/2.)/t,2.);
}

/*---------------------------------------------------------------------*/

double phi22(double delta, double z)
{
	if (z==0.) return 0.;

	int ie;

	int n1=500;
	int n2=500;

	double t=0.;
	double int1=(1.-z*t)*G1(t);
	
	for(ie=1;ie<=n1;ie++)
	{
		t+=(1.-delta)/z/n1;
		int1 +=(1.-z*t)*G1(t);
	}
	int1 *= delta*(1.-delta)/z/n1;
	
	t=(1.-delta)/z;
	double int2=pow(1.-z*t,2.)*G1(t);

	for(ie=1;ie<=n2;ie++)
	{
		t+=delta/z/n2;
		int2 +=pow(1.-z*t,2.)*G1(t);
	}
	int2 *= delta/z/n2;

	return 16./27.*z*(int1+int2);
}

/*---------------------------------------------------------------------*/

double phi11(double delta, double z)
{
	return phi22(delta,z)/36.;
}

/*---------------------------------------------------------------------*/

double phi12(double delta, double z)
{
	return -phi22(delta,z)/3.;
}

/*---------------------------------------------------------------------*/

double G2(double t)
{
/* G2(t)=Re(G(t)+t/2) */
	
	if (t<4.) return -2.*pow(atan(sqrt(t/(4.-t))),2.)+t/2.;
	else return -pi*pi/2.+2.*pow(log((sqrt(t)+sqrt(t-4.))/2.),2.)+t/2;
}

/*---------------------------------------------------------------------*/

double phi27(double delta, double z)
{
	if (z==0.) return 0.;
		
	int ie;
	int n1=500;
	int n2=500;

	double t=0.;
	double int1=G2(t);
	
	for(ie=1;ie<=n1;ie++)
	{
		t+=(1.-delta)/z/n1;
		int1 +=G2(t);
	}
	int1 *= delta*(1.-delta)/z/n1;
	
	t=(1.-delta)/z;
	double int2=pow(1.-z*t,2.)*G2(t);

	for(ie=1;ie<=n2;ie++)
	{
		t+=delta/z/n2;
		int2 +=(1.-z*t)*G2(t);
	}
	int2 *= delta/z/n2;

	return -8./9.*z*z*(int1+int2);
}

/*---------------------------------------------------------------------*/

double phi17(double delta, double z)
{
	return -phi27(delta,z)/6.;
}

/*---------------------------------------------------------------------*/

double phi18(double delta, double z)
{
	return phi27(delta,z)/18.;
}

/*---------------------------------------------------------------------*/

double phi28(double delta, double z)
{
	return -phi27(delta,z)/3.;
}

/*---------------------------------------------------------------------*/

double phi47(double delta)
{
	return -1./54.*delta*(1.-delta+delta*delta/3.) + phi27(delta,1.)/12.;
}

/*---------------------------------------------------------------------*/

double phi48(double delta)
{
	return -phi47(delta)/3.;
}

/*---------------------------------------------------------------------*/

double F2_nf(double z)
{
	double F2_nf;
	
	if (z==0.) return 0.;
	F2_nf=-log(1.-z)*log(1.-z)/(1.-z)/2.-13./36.*log(1.-z)/(1.-z)+(-pi*pi/18.+85./72.)/(1.-z)
	+(z*z-3.)/6./(z-1.)*Li2(1.-z)+(z*z-3.)/6./(z-1.)*log(1.-z)*log(z)-(1.+z)*log(1.-z)*log(1.-z)/4.-(6.*z*z-25.*z-1.)*log(1.-z)/36.+log(1.-z)/z/2.-(1.+z)*pi*pi/36.+(-49.+38.*z*z-55.*z)/72.;

	return F2_nf;
}

/*---------------------------------------------------------------------*/

double phi77_2beta(double delta,double mu, struct parameters* param)
{
	double Lb=log(mu*mu/param->mass_b_1S/param->mass_b_1S);
	int nf=5;
	double beta0=11.-2./3.*nf;

	int ie;
	int n1=500;

	double t=0.;
	double int1=F2_nf(t);
	
	for(ie=1;ie<=n1;ie++)
	{
		t+=(1.-delta)/n1;
		int1 +=F2_nf(t);
	}
	int1 *=(1.-delta)/n1;

	return beta0*(phi77(delta)*Lb+4.*int1);
}

/*---------------------------------------------------------------------*/

double F2_a(double z)
{
	double F2_a;
	if(z==0.) return 0.;

	F2_a=0.5*pow(log(1.-z),3.)/(1.-z)+21./8.*pow(log(1.-z),2.)/(1.-z)
	+(-pi*pi/6.+271./48.)*log(1.-z)/(1.-z)
	+(425./96.-pi*pi/6.-zeta3/2.)/(1.-z)
	+(4.*z-4.*z*z+1.+z*z*z)/2./(z-1.)*(Li3(z/(2.-z))-Li3(-z/(2.-z))-2.*Li3(1./(2.-z))+zeta3/4.)
	+((z*z*z-2.*z*z+2.*z-3.)/2./(z-1.)*log(1.-z)-(-140.*pow(z,4.)+219.*pow(z,3.)-124.*z*z+28.*z+27.*pow(z,5.)+9.*pow(z,6.)+pow(z,8.)-6.*pow(z,7.)-6.)/12./z/pow(z-1.,3.))*Li2(z-1.)
-2.*pow(z-1.,2.)*Li3(z-1.)+((2.*pow(z,3.)-9.*z*z-2.*z+11.)/4./(z-1.)*log(1.-z)-(-27.*z*z+8.*pow(z,6.)-9.+21.*z-3.*z*z*z+64.*pow(z,4.)-46.*pow(z,5.))/12./z/pow(z-1.,3.))*Li2(1.-z)
-(-17.*z*z+4.*z+4.*z*z*z+11.)/4./(z-1.)*Li3(1.-z)-(2.*z*z*z+13.-9.*z*z)/4./(z-1.)*Li3(z)+(4.*z-4.*z*z+1.+z*z*z)/6./(z-1.)*pow(log(2.-z),3.) 
+(-(4.*z-4.*z*z+1.+z*z*z)/2./(z-1.)*pow(log(1.-z),2.)-(-140.*pow(z,4.)+219.*pow(z,3.)-124.*z*z+28.*z+27.*pow(z,5.)+9.*pow(z,6.)+pow(z,8.)-6.*pow(z,7.)-6.)/12./z/pow(z-1.,3.)*log(1.-z) - (4.*z-4.*z*z+1.+z*z*z)/(z-1.)*pi*pi/12.)*log(2.-z)
+(z*z*z-2.*z*z+2.*z+1.)/4./z*pow(log(1.-z),3.)+(pow(z,5.)-3.*pow(z,4.)+5.*pow(z,3.)+7.*z*z+5.*z-9.)/24./z*pow(log(1.-z),2.)
+(-(z*z+8.*z-11.)/8./(z-1.)*pow(log(1.-z),2.)-(-27.*z*z+8.*pow(z,6.)-9.+21.*z-3.*z*z*z+64.*pow(z,4.)-46.*pow(z,5.))/12./z/pow(z-1.,3.)*log(1.-z))*log(z)
+((-z*z+z-3.)*pi*pi/12.-(4.*pow(z,5.)+151.*z+2.*pow(z,4.)-48.*z*z-41.*z*z*z-36.)/48./z*(z-1.))*log(1.-z)
-(z-2.)*(pow(z,4.)-z*z*z-11.*z*z+13.*z+3.)/z*pi*pi/72. + (z*z*z-11.*z*z-2.*z+18.)/4./(z-1.)*zeta3-(8.*pow(z,4.)-244.*z*z*z+175.*z*z+598.*z-569.)/96./(z-1.);
	return F2_a;
}

/*---------------------------------------------------------------------*/

double F2_na(double z)
{
	double F2_na;
	if(z==0.) return 0.;

	F2_na=11./8.*pow(log(1.-z),2.)/(1.-z)+(pi*pi/12.+95./144.)*log(1.-z)/(1.-z)+(zeta3/4.-905./288.+17.*pi*pi/72.)/(1.-z)
	-(4.*z-4.*z*z+1.+z*z*z)/4./(z-1.)*(Li3(z/(2.-z))-Li3(-z/(2.-z))-2.*Li3(1./(2.-z))+zeta3/4.)+pow(z-1.,2.)*Li3(z-1.)
	+(-(z*z*z-2.*z*z+2.*z-3.)/4./(z-1.)*log(1.-z)+(-140.*pow(z,4.)+219.*z*z*z-124.*z*z+28.*z+27.*pow(z,5.)+9.*pow(z,6.)+pow(z,8.)-6.*pow(z,7.)-6.)/24./z/pow(z-1.,3.))*Li2(z-1.)
	+(z*(3.-z)/4.*log(1.-z)+(1.+z)*(2.*pow(z,4.)-29.*pow(z,3.)+73.*z*z-57.*z+15.)/24./pow(z-1.,3.))*Li2(1.-z)+(4.*z-4.*z*z+1.+z*z*z)/4./(z-1.)*Li3(z)
	+(z-3.)*z/2.*Li3(1.-z)-(4.*z-4.*z*z+1.+z*z*z)/12./(z-1.)*pow(log(2.-z),3.)+((4.*z-4.*z*z+1.+z*z*z)/4./(z-1.)*pow(log(1.-z),2.)
	+(-140.*pow(z,4.)+219.*pow(z,3.)-124.*z*z+28.*z+27.*pow(z,5.)+9.*pow(z,6.)+pow(z,8.)-6.*pow(z,7.)-6.)/24./z/pow(z-1.,3.)*log(1.-z)+(4.*z-4.*z*z+1.+z*z*z)*pi*pi/24./(z-1.))*log(2.-z)
	+(1.+z)*(2.*pow(z,4.)-29.*z*z*z+73.*z*z-57.*z+15.)/24./pow(z-1.,3.)*log(1.-z)*log(z)-pow(z-1.,2.)/8.*pow(log(1.-z),3.)-(z+2.)*(z*z*z-5.*z*z+9.*z-35.)/48.*pow(log(1.-z),2.)
	+((z*z-z+3.)*pi*pi/24.+(6.*pow(z,5.)+72.-392.*pow(z,3.)+51.*pow(z,4.)+219.*z*z+92.*z)/144./z/(z-1.))*log(1.-z)
	+(pow(z,5.)-3.*pow(z,4.)-3.*pow(z,3.)+34.*z*z-24.*z+3.)/z*pi*pi/144.
	-(z*z*z-10.*z*z+6.*z+7.)/8./(z-1.)*zeta3+(12.*pow(z,4.)-754.*pow(z,3.)+1191.*z*z+264.*z-761.)/288./(z-1.);	
	
	return F2_na;
}

/*---------------------------------------------------------------------*/

double phi77_2rem(double delta, struct parameters* param)
{
	int ie;
	int n1=500;
	double alphas_upsilon=alphas_running(param->mass_b_1S,param->mass_top_pole,param->mass_b_pole,param);

	double t=0.;
	double int1=F2_a(t);
	double int2=F2_na(t);
	double int3=F2_nf(t);
	
	for(ie=1;ie<=n1;ie++)
	{
		t+=(1.-delta)/n1;
		int1 +=F2_a(t);
		int2 +=F2_na(t);
		int3 +=F2_nf(t);
	}
	int1 *=(1.-delta)/n1;
	int2 *=(1.-delta)/n1;
	int3 *=(1.-delta)/n1;

	return -4.*(16./9.*int1+4.*int2+29./3.*int3)-8.*pi*alphas_upsilon/27./delta*(2.*delta*log(delta)*log(delta)+(4.+7.*delta-2.*delta*delta+pow(delta,3.))*log(delta)+7.-8./3.*delta-7.*delta*delta+4.*pow(delta,3.)-4./3.*pow(delta,4.));
}

/*---------------------------------------------------------------------*/

double Re_a(double z)
{
	if(z==1.) return 4.0859;
	if (z==0.) return 0.;
		
	double Lz=log(z);

	return 16./9.*((5./2.-pi*pi/3.-3.*zeta3+(5./2.-3./4.*pi*pi)*Lz+Lz*Lz/4.+Lz*Lz*Lz/12.)*z
	+(7./4.+2./3.*pi*pi-0.5*pi*pi*Lz-Lz*Lz/4.+Lz*Lz*Lz/12.)*z*z+(-7./6.-pi*pi/4.+2.*Lz-3./4.*Lz*Lz)*z*z*z
	+(457./216.-5./18.*pi*pi-Lz/72.-5./6.*Lz*Lz)*pow(z,4.)+(35101./8640.-35./72.*pi*pi-185./144.*Lz-35./24.*Lz*Lz)*pow(z,5.)
	+(67801./8000.-21./20.*pi*pi-3303./800.*Lz-63./20.*Lz*Lz)*pow(z,6.));
}

/*---------------------------------------------------------------------*/

double Re_b(double z)
{

	if(z==1.) return 0.0316;
	if (z==0.) return 0.;
		
		double Lz=log(z);

	return -8./9.*((-3.+pi*pi/6.-Lz)*z-2./3.*pi*pi*pow(z,1.5)+(0.5+pi*pi-2.*Lz-0.5*Lz*Lz)*z*z
	+(-25./12.-pi*pi/9.-19./18.*Lz+2.*Lz*Lz)*z*z*z+(-1376./225.+137./30.*Lz+2.*Lz*Lz+2./3.*pi*pi)*pow(z,4.)
	+(-131317./11760.+887./84.*Lz+5.*Lz*Lz+5./3.*pi*pi)*pow(z,5.)+(-2807617./97200.+16597./540.*Lz+14.*Lz*Lz+14./3.*pi*pi)*pow(z,6.));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

double bsgamma(double C0[], double C1[], double C2[], double mu, double mu_W, struct parameters* param)
/* computes the inclusive branching ratio of B -> Xs gamma */
/* C0, C1, C2: respectively LO, NLO and NNLO contributions of the Wilson coefficients at scale mu=O(mb) */
{
	int ie,je;	
	double alpha_em=1./137.036;
	double BR_BXcenu_exp=0.1061;
	double Cbr=0.580;
	double E0=1.6;
	
	double alphas_mu=alphas_running(mu,param->mass_top_pole,param->mass_b_pole,param);


	double P0 = C0[7]*C0[7];
	double P1_1 = 2.*C0[7]*C1[7];
	double P2_1 = C1[7]*C1[7] + 2.*C0[7]*C2[7];
	

	double K1[9][9], phi1[9][9], gamma0eff[9][9];
	for (ie=1;ie<=8;ie++) for (je=1;je<=8;je++) K1[ie][je]=phi1[ie][je]=gamma0eff[ie][je]=0.;
	
	gamma0eff[1][7]=-208./243.; 
	gamma0eff[2][7]=416./81.;
	gamma0eff[3][7]=-176./81.;
	gamma0eff[4][7]=-152./243.;
	gamma0eff[5][7]=-6272./81.;
	gamma0eff[6][7]=4624./243.;
	gamma0eff[7][7]=32./3.;
	gamma0eff[8][7]=-32./9.;
	
	
	double z=(param->mass_c/param->mass_b_1S)*(param->mass_c/param->mass_b_1S);
		
	double delta=1.-2.*E0/param->mass_b_1S;

	phi1[2][2]=phi22(delta,z);
	phi1[1][1]=phi11(delta,z);
	phi1[1][2]=phi12(delta,z);
	phi1[2][7]=phi27(delta,z);
	phi1[1][7]=phi17(delta,z);
	phi1[1][8]=phi18(delta,z);
	phi1[2][8]=phi28(delta,z);
	phi1[4][7]=phi47(delta);
	phi1[4][8]=phi48(delta);	
	phi1[7][7]=phi77(delta);
	phi1[7][8]=phi78(delta);
	phi1[8][8]=phi88(delta,param->mass_b_1S/param->mass_s);
	
	double Lb=log(mu*mu/param->mass_b_1S/param->mass_b_1S);
	
	double Re_r1[7];
	for (ie=1;ie<=6;ie++) Re_r1[ie]=0.;
	Re_r1[2]=-1666./243.+2.*(Re_a(z)+Re_b(z));
	Re_r1[1]=-Re_r1[2]/6.;
	double Xb = -0.1684;
	Re_r1[3]=2392./243.+8.*pi/3./sqrt(3.)+32./9.*Xb-Re_a(1.)+2.*Re_b(1.);
	Re_r1[4]=-761./729.-4.*pi/9./sqrt(3.)-16./27.*Xb+Re_a(1.)/6.+5./3.*Re_b(1.)+2.*Re_b(z);
	Re_r1[5]=56680./243.+32.*pi/3./sqrt(3.)+128./9.*Xb-16.*Re_a(1.)+32.*Re_b(1.);
	Re_r1[6]=5710./729.-16.*pi/9./sqrt(3.)-64./27.*Xb-10./3.*Re_a(1.)+44./3.*Re_b(1.)+12.*Re_a(z)+20.*Re_b(z);
	
	
	for (ie=1;ie<=8;ie++) for (je=1;je<=8;je++) K1[ie][je]=2.*phi1[ie][je];
	for (ie=1;ie<=8;ie++) K1[ie][ie]=4.*phi1[ie][ie];
	for (ie=1;ie<=6;ie++) K1[ie][7]=Re_r1[ie]-0.5*gamma0eff[ie][7]*Lb+2.*phi1[ie][7];
	K1[7][7]=-182./9.+8./9.*pi*pi-gamma0eff[7][7]*Lb+4.*phi1[7][7];
	K1[7][8]=44./9.-8./27.*pi*pi-0.5*gamma0eff[8][7]*Lb+2.*phi1[7][8];
	
	for(je=1;je<=8;je++) for(ie=1;ie<=8;ie++) K1[ie][je]=K1[je][ie];
		
	double P1_2=0.;
	double P2_3=0.;
	for (ie=1;ie<=8;ie++) for (je=1;je<=8;je++)
	{ 
		P1_2+=C0[ie]*C0[je]*K1[ie][je];
		P2_3+=2.*C0[ie]*C1[je]*K1[ie][je];
	}
	
	

	double Lz=log(z);

	double Re_r2=
	67454./6561.-124./729.*pi*pi-4./1215.*(11280.-1520.*pi*pi-171.*pow(pi,4.)-5760.*zeta3+6840.*Lz
	-1440.*pi*pi*Lz-2520.*zeta3*Lz+120.*Lz*Lz+100.*pow(Lz,3.)-30.*pow(Lz,4.))*z
	-64./243.*pi*pi*(43.-12.*log(2.)-3.*Lz)*pow(z,1.5)-2./1215.*(11475.-380.*pi*pi+96.*pow(pi,4.)+7200.*zeta3
	-1110.*Lz-1560.*pi*pi*Lz+1440.*zeta3*Lz+990.*Lz*Lz+260.*pow(Lz,3.)-60.*pow(Lz,4.))*z*z
	+2240./243.*pi*pi*pow(z,2.5)-2./2187.*(62471.-2424.*pi*pi-33264.*zeta3-19494.*Lz-504.*pi*pi*Lz
	-5184.*Lz*Lz+2160.*pow(Lz,3.))*pow(z,3.)-2464./6075.*pi*pi*pow(z,3.5)+(-15103841./546750.+7912./3645.*pi*pi+2368./81.*zeta3
	+147038./6075.*Lz+352./243.*pi*pi*Lz+88./243.*Lz*Lz-512./243.*pow(Lz,3.))*pow(z,4.);

	int nf=5;
	double beta0=11.-2./3.*nf;

	double K2_beta[9][9], K2_rem[9][9],phi2_beta[9][9];
	for (ie=1;ie<=8;ie++) for (je=1;je<=8;je++) K2_beta[ie][je]=K2_rem[ie][je]=phi2_beta[ie][je]=0.;
	
	phi2_beta[7][7]=phi77_2beta(delta,mu,param);
	
	for (ie=1;ie<=8;ie++) for (je=1;je<=8;je++) K2_beta[ie][je]=2.*phi2_beta[ie][je]; 
	
	for (ie=1;ie<=8;ie++) K2_beta[ie][ie]=4.*phi2_beta[ie][ie]; 
	
	K2_beta[2][7]=beta0*(-1.5*Re_r2+2.*(Re_a(z)+Re_b(z)-290./81.)*Lb-100./81.*Lb*Lb)+2.*phi2_beta[2][7]; 
	K2_beta[1][7]=-K2_beta[2][7]/6.; 		
	K2_beta[7][7]=
	beta0*(-3803./54.-46./27.*pi*pi+80./3.*zeta3+(8./9.*pi*pi-98./3.)*Lb-16./3.*Lb*Lb)+4.*phi2_beta[7][7];
	K2_beta[7][8]= beta0*(1256./81.-64./81.*pi*pi-32./9.*zeta3+(188./27.-8./27.*pi*pi)*Lb+8./9.*Lb*Lb)+2.*phi2_beta[7][8]; 

	for(je=1;je<=8;je++) for(ie=1;ie<=8;ie++) K2_beta[ie][je]=K2_beta[je][ie];
   
	double P2_2_beta=0.;
	for (ie=1;ie<=8;ie++) for (je=1;je<=8;je++) P2_2_beta+=C0[ie]*C0[je]*K2_beta[ie][je];
	
	

	double x1,x2,x3,x4,x5;
	
	x1=C0[2]*C0[2]+C0[1]*C0[1]/36.-C0[1]*C0[2]/3.;
	
	x2=9./8.*(-2.*4736./729.*C0[2]*C0[7]+(1./3.*4736./729.+140./27.)*C0[1]*C0[7]
	-32./9.*C0[7]*C0[7]+32./27.*C0[7]*C0[8]);
	

	x5=0.;
	double z0=1.e10;
	
	double Lz0=log(z0);
	double LD=Lb-Lz0;
	double Lc=0.;
	
	double phi77_2=phi77_2rem(delta,param);
	double alphas_upsilon=alphas_running(param->mass_b_1S,param->mass_top_pole,param->mass_b_pole,param);
	
	double K22rem=pow(218./243.-208./81.*LD,2.);
	double K11rem=K22rem/36.;
	double K12rem=-K22rem/6.;
	double K27rem=(218./243.-208./81.*LD)*K1[7][7]+(127./324.-35./27.*LD)*K1[7][8]+2./3.*(1.-LD)*(K1[4][7]-beta0*(26./81.-4./27.*Lb))
	-4736./729.*LD*LD+1150./729.*LD-1617980./19683.+20060./243.*zeta3+1664./81.*Lc;
	double K28rem=(218./243.-208./81.*LD)*K1[7][8]+(127./324.-35./27.*LD)*K1[8][8]+2./3.*(1.-LD)*K1[4][8];
	double K17rem=-1./6.*K27rem+(5./16.-3./4.*LD)*K1[7][8]-1237./729.+232./27.*zeta3+70./27.*LD*LD-20./27.*LD;
	double K18rem=-1./6.*K28rem+(5./16.-3./4.*LD)*K1[8][8];
	double K77rem=(K1[7][7]-4.*phi77(delta)+2./3.*Lz)*K1[7][7]-32./9.*LD*LD+224./27.*LD-628487./729.-628./405.*pow(pi,4.)
	+31823./729.*pi*pi+428./27.*pi*pi*log(2.)+26590./81.*zeta3-160./3.*Lb*Lb-2720./9.*Lb+256./27.*pi*pi*Lb
	+512./27.*pi*alphas_upsilon+4.*phi77_2;
	
	double K78rem=(-50./3.+8./3.*pi*pi-2./3.*LD)*K1[7][8]+16./27.*LD*LD-112./81.*LD+364./243.;
	double K88rem=(-50./3.+8./3.*pi*pi-2./3.*LD)*K1[8][8];
	
	double P22rem_z0 = C0[2]*C0[2]*K22rem + C0[1]*C0[1]*K11rem + 2.*C0[1]*C0[2]*K12rem + 2.* C0[2]*C0[7]*K27rem + 2.*C0[2]*C0[8]*K28rem + 2.*C0[1]*C0[7]*K17rem + 2.*C0[1]*C0[8]*K18rem + C0[7]*C0[7]*K77rem + 2.* C0[7]*C0[8]*K78rem + C0[8]*C0[8]*K88rem;

	double Re_r12_z0 = -1666./243.+2.*(4./3.*Lz0+34./9.-4./81.*Lz0+8./81.); 
	double Im_r1_2_z0 = -80./81.*pi+2.*(4./9.+4./81.)*pi;
	double Re_r2_z0 = 8./9.*Lz0*Lz0+112./243.*Lz0+27650./6561.;
	double z0_dRe_r1_2_dz0=2.*(4./3.-4./81.);


	double z1=1.e20;
	
	double Lz1=log(z1);
	LD=Lb-Lz1;
	
	K22rem=pow(218./243.-208./81.*LD,2.);
	K11rem=K22rem/36.;
	K12rem=-K22rem/6.;
	K27rem=(218./243.-208./81.*LD)*K1[7][7]+(127./324.-35./27.*LD)*K1[7][8]+2./3.*(1.-LD)*(K1[4][7]-beta0*(26./81.-4./27.*Lb))
	-4736./729.*LD*LD+1150./729.*LD-1617980./19683.+20060./243.*zeta3+1664./81.*Lc;
	K28rem=(218./243.-208./81.*LD)*K1[7][8]+(127./324.-35./27.*LD)*K1[8][8]+2./3.*(1.-LD)*K1[4][8];
	K17rem=-1./6.*K27rem+(5./16.-3./4.*LD)*K1[7][8]-1237./729.+232./27.*zeta3+70./27.*LD*LD-20./27.*LD;
	K18rem=-1./6.*K28rem+(5./16.-3./4.*LD)*K1[8][8];
	K77rem=(K1[7][7]-4.*phi77(delta)+2./3.*Lz)*K1[7][7]-32./9.*LD*LD+224./27.*LD-628487./729.-628./405.*pow(pi,4.)
	+31823./729.*pi*pi+428./27.*pi*pi*log(2.)+26590./81.*zeta3-160./3.*Lb*Lb-2720./9.*Lb+256./27.*pi*pi*Lb
	+512./27.*pi*alphas_upsilon+4.*phi77_2;
	K78rem=(-50./3.+8./3.*pi*pi-2./3.*LD)*K1[7][8]+16./27.*LD*LD-112./81.*LD+364./243.;
	K88rem=(-50./3.+8./3.*pi*pi-2./3.*LD)*K1[8][8];
	
	double P22rem_z1 = C0[2]*C0[2]*K22rem + C0[1]*C0[1]*K11rem + 2.*C0[1]*C0[2]*K12rem + 2.* C0[2]*C0[7]*K27rem + 2.*C0[2]*C0[8]*K28rem + 2.*C0[1]*C0[7]*K17rem + 2.*C0[1]*C0[8]*K18rem + C0[7]*C0[7]*K77rem + 2.* C0[7]*C0[8]*K78rem + C0[8]*C0[8]*K88rem;

	double Re_r12_z1 = -1666./243.+2.*(4./3.*Lz1+34./9.-4./81.*Lz1+8./81.);
	double Im_r1_2_z1 = -80./81.*pi+2.*(4./9.+4./81.)*pi;
	double Re_r2_z1 = 8./9.*Lz1*Lz1+112./243.*Lz1+27650./6561.;
	double z1_dRe_r1_2_dz1=2.*(4./3.-4./81.);

	double a= P22rem_z0- (x1*(Re_r12_z0*Re_r12_z0-(1666./243.)*(1666./243.)+Im_r1_2_z0*Im_r1_2_z0-(80./81.*pi*80./81.*pi))
	+x2*(Re_r2_z0-(67454./6561.-124./729.*pi*pi))
	+x5);
	double b= P22rem_z1- (x1*(Re_r12_z1*Re_r12_z1-(1666./243.)*(1666./243.)+Im_r1_2_z1*Im_r1_2_z1-(80./81.*pi*80./81.*pi))
	+x2*(Re_r2_z1-(67454./6561.-124./729.*pi*pi))
	+x5);
	
	double a3=(Re_r12_z0-(-1666./243.));
	double b3=(Re_r12_z1-(-1666./243.));
	
	double a4=z0_dRe_r1_2_dz0;
	double b4=z1_dRe_r1_2_dz1;
	
	x3=(a*b4-a4*b)/(a3*b4-b3*a4);
	x4=(a*b3-a3*b)/(a4*b3-b4*a3);

#ifdef DEBUG
	printf("---------------\n");
	printf("Branching Ratio\n");
	printf("---------------\n");
	printf("cas (a) : x1=%f\t x2=%f\t x3=%f\t x4=%f\t x5=%f\n",x1,x2,x3,x4,x5);
#endif


	double Im_r1_2 = -80./81.*pi+2.*pi*(16./9.*((4.-pi*pi/3.+Lz+Lz*Lz)*z/2.+(0.5-pi*pi/6.-Lz-0.5*Lz*Lz)*z*z+pow(z,3.)+5./9.*pow(z,4.))-8./9.*(-z+(1.-2.*Lz)*z*z+(-10./9.+4./3.*Lz)*pow(z,3.)+pow(z,4.)));
	
	double dRe_r1_2_dz=2.*(Re_a(z+z*1.e-4)-Re_a(z-z*1.e-4)+Re_b(z+z*1.e-4)-Re_b(z-z*1.e-4))/(z*2.e-4);

	double P2_2_rem1=x1*(Re_r1[2]*Re_r1[2]-(1666./243.)*(1666./243.)+Im_r1_2*Im_r1_2-(80./81.*pi*80./81.*pi))
	+x2*(Re_r2-(67454./6561.-124./729.*pi*pi))
	+x3*(Re_r1[2]-(-1666./243.))
	+x4*z*dRe_r1_2_dz
	+x5;

	
/* --------------------------------------------- */

	double K1z0[9][9], phi1z0[9][9];
	for (ie=1;ie<=8;ie++) for (je=1;je<=8;je++) K1z0[ie][je]=phi1z0[ie][je]=0.;
		
	phi1z0[2][2]=phi22(delta,0.);
	phi1z0[1][1]=phi11(delta,0.);
	phi1z0[2][7]=phi27(delta,0.);
	phi1z0[1][7]=phi17(delta,0.);
	phi1z0[1][8]=phi18(delta,0.);
	phi1z0[2][8]=phi28(delta,0.);
	phi1z0[4][7]=phi1[4][7];
	phi1z0[4][8]=phi1[4][8];	
	phi1z0[7][7]=phi1[7][7];
	phi1z0[7][8]=phi1[7][8];
	phi1z0[8][8]=phi1[8][8];
	
	double Re_r1z0[7];
	for (ie=1;ie<=6;ie++) Re_r1z0[ie]=0.;
	Re_r1z0[2]=-1666./243.+2.*(Re_a(0.)+Re_b(0.));
	Re_r1z0[1]=-1./6.*Re_r1z0[2];
	Re_r1z0[3]=2392./243.+8.*pi/3./sqrt(3.)+32./9.*Xb-Re_a(1.)+2.*Re_b(1.);
	Re_r1z0[4]=-761./729.-4.*pi/9./sqrt(3.)-16./27.*Xb+Re_a(1.)/6.+5./3.*Re_b(1.)+2.*Re_b(0.);
	Re_r1z0[5]=56680./243.+32.*pi/3./sqrt(3.)+128./9.*Xb-16.*Re_a(1.)+32.*Re_b(1.);
	Re_r1z0[6]=5710./729.-16.*pi/9./sqrt(3.)-64./27.*Xb-10./3.*Re_a(1.)+44./3.*Re_b(1.)+12.*Re_a(0.)+20.*Re_b(0.);
	
	
	for (ie=1;ie<=8;ie++) for (je=1;je<=8;je++) K1z0[ie][je]=2.*phi1z0[ie][je];
	for (ie=1;ie<=8;ie++) K1z0[ie][ie]=2.*2.*phi1z0[ie][ie];
	for (ie=1;ie<=6;ie++) K1z0[ie][7]=Re_r1z0[ie]-0.5*gamma0eff[ie][7]*Lb+2.*phi1z0[ie][7];
	K1z0[7][7]=-182./9.+8./9.*pi*pi-gamma0eff[7][7]*Lb+4.*phi1z0[7][7];
	K1z0[7][8]=44./9.-8./27.*pi*pi-0.5*gamma0eff[8][7]*Lb+2.*phi1z0[7][8];
	
	for(je=1;je<=8;je++) for(ie=1;ie<=8;ie++) K1z0[ie][je]=K1z0[je][ie];
		
	double P2_3z0=0.;
	for (ie=1;ie<=8;ie++) for (je=1;je<=8;je++) P2_3z0+=2.*C0[ie]*C1[je]*K1z0[ie][je];
	

	x5=-P2_3z0-P2_1;
	

	LD=Lb-Lz0;
	
	K22rem=pow(218./243.-208./81.*LD,2.);
	K11rem=K22rem/36.;
	K12rem=-K22rem/6.;
	K27rem=(218./243.-208./81.*LD)*K1[7][7]+(127./324.-35./27.*LD)*K1[7][8]+2./3.*(1.-LD)*(K1[4][7]-beta0*(26./81.-4./27.*Lb))
	-4736./729.*LD*LD+1150./729.*LD-1617980./19683.+20060./243.*zeta3+1664./81.*Lc;
	K28rem=(218./243.-208./81.*LD)*K1[7][8]+(127./324.-35./27.*LD)*K1[8][8]+2./3.*(1.-LD)*K1[4][8];
	K17rem=-1./6.*K27rem+(5./16.-3./4.*LD)*K1[7][8]-1237./729.+232./27.*zeta3+70./27.*LD*LD-20./27.*LD;
	K18rem=-1./6.*K28rem+(5./16.-3./4.*LD)*K1[8][8];
	K77rem=(K1[7][7]-4.*phi77(delta)+2./3.*Lz)*K1[7][7]-32./9.*LD*LD+224./27.*LD-628487./729.-628./405.*pow(pi,4.)
	+31823./729.*pi*pi+428./27.*pi*pi*log(2.)+26590./81.*zeta3-160./3.*Lb*Lb-2720./9.*Lb+256./27.*pi*pi*Lb
	+512./27.*pi*alphas_upsilon+4.*phi77_2;
	
	K78rem=(-50./3.+8./3.*pi*pi-2./3.*LD)*K1[7][8]+16./27.*LD*LD-112./81.*LD+364./243.;
	K88rem=(-50./3.+8./3.*pi*pi-2./3.*LD)*K1[8][8];
	
	P22rem_z0 = C0[2]*C0[2]*K22rem + C0[1]*C0[1]*K11rem + 2.*C0[1]*C0[2]*K12rem + 2.* C0[2]*C0[7]*K27rem + 2.*C0[2]*C0[8]*K28rem + 2.*C0[1]*C0[7]*K17rem + 2.*C0[1]*C0[8]*K18rem + C0[7]*C0[7]*K77rem + 2.* C0[7]*C0[8]*K78rem + C0[8]*C0[8]*K88rem;

	Re_r12_z0 = -1666./243.+2.*(4./3.*Lz0+34./9.-4./81.*Lz0+8./81.);
	Im_r1_2_z0 = -80./81.*pi+2.*(4./9+4./81.)*pi;
	Re_r2_z0 = 8./9.*Lz0*Lz0+112./243.*Lz0+27650./6561.;
	z0_dRe_r1_2_dz0=2.*(4./3.-4./81.);


	LD=Lb-Lz1;
	
	K22rem=pow(218./243.-208./81.*LD,2.);
	K11rem=K22rem/36.;
	K12rem=-K22rem/6.;
	K27rem=(218./243.-208./81.*LD)*K1[7][7]+(127./324.-35./27.*LD)*K1[7][8]+2./3.*(1.-LD)*(K1[4][7]-beta0*(26./81.-4./27.*Lb))
	-4736./729.*LD*LD+1150./729.*LD-1617980./19683.+20060./243.*zeta3+1664./81.*Lc;
	K28rem=(218./243.-208./81.*LD)*K1[7][8]+(127./324.-35./27.*LD)*K1[8][8]+2./3.*(1.-LD)*K1[4][8];
	K17rem=-1./6.*K27rem+(5./16.-3./4.*LD)*K1[7][8]-1237./729.+232./27.*zeta3+70./27.*LD*LD-20./27.*LD;
	K18rem=-1./6.*K28rem+(5./16.-3./4.*LD)*K1[8][8];
	K77rem=(K1[7][7]-4.*phi77(delta)+2./3.*Lz)*K1[7][7]-32./9.*LD*LD+224./27.*LD-628487./729.-628./405.*pow(pi,4.)
	+31823./729.*pi*pi+428./27.*pi*pi*log(2.)+26590./81.*zeta3-160./3.*Lb*Lb-2720./9.*Lb+256./27.*pi*pi*Lb
	+512./27.*pi*alphas_upsilon+4.*phi77_2;
	K78rem=(-50./3.+8./3.*pi*pi-2./3.*LD)*K1[7][8]+16./27.*LD*LD-112./81.*LD+364./243.;
	K88rem=(-50./3.+8./3.*pi*pi-2./3.*LD)*K1[8][8];
	
	P22rem_z1 = C0[2]*C0[2]*K22rem + C0[1]*C0[1]*K11rem + 2.*C0[1]*C0[2]*K12rem + 2.* C0[2]*C0[7]*K27rem + 2.*C0[2]*C0[8]*K28rem + 2.*C0[1]*C0[7]*K17rem + 2.*C0[1]*C0[8]*K18rem + C0[7]*C0[7]*K77rem + 2.* C0[7]*C0[8]*K78rem + C0[8]*C0[8]*K88rem;

	Re_r12_z1 = -1666./243.+2.*(4./3.*Lz1+34./9.-4./81.*Lz1+8./81.);
	Im_r1_2_z1 = -80./81.*pi+2.*(4./9.+4./81.)*pi;
	Re_r2_z1 = 8./9.*Lz1*Lz1+112./243.*Lz1+27650./6561.;
	z1_dRe_r1_2_dz1=2.*(4./3.-4./81.);


	a= P22rem_z0- (x1*(Re_r12_z0*Re_r12_z0-(1666./243.)*(1666./243.)+Im_r1_2_z0*Im_r1_2_z0-(80./81.*pi*80./81.*pi))
	+x2*(Re_r2_z0-(67454./6561.-124./729.*pi*pi))
	+x5);
	b= P22rem_z1- (x1*(Re_r12_z1*Re_r12_z1-(1666./243.)*(1666./243.)+Im_r1_2_z1*Im_r1_2_z1-(80./81.*pi*80./81.*pi))
	+x2*(Re_r2_z1-(67454./6561.-124./729.*pi*pi))
	+x5);
	
	a3=(Re_r12_z0-(-1666./243.));
	b3=(Re_r12_z1-(-1666./243.));
	
	a4=z0_dRe_r1_2_dz0;
	b4=z1_dRe_r1_2_dz1;
	
	x3=(a*b4-a4*b)/(a3*b4-b3*a4);
	x4=(a*b3-a3*b)/(a4*b3-b4*a3);

#ifdef DEBUG
	printf("cas (b) : x1=%f\t x2=%f\t x3=%f\t x4=%f\t x5=%f\n\n",x1, x2, x3, x4, x5);
#endif


	double P2_2_rem2=x1*(Re_r1[2]*Re_r1[2]-(1666./243.)*(1666./243.)+Im_r1_2*Im_r1_2-(80./81.*pi*80./81.*pi))
	+x2*(Re_r2-(67454./6561.-124./729.*pi*pi))
	+x3*(Re_r1[2]-(-1666./243.))
	+x4*z*dRe_r1_2_dz
	+x5;
	
	
	double P2_2_rem=(P2_2_rem1+P2_2_rem2)/2.;
	double P2_2=P2_2_beta+P2_2_rem;

	if(fabs((alphas_mu/4./pi*alphas_mu/4./pi *(P2_1 + P2_2 + P2_3))/(P0 + alphas_mu/4./pi * (P1_1 + P1_2))) > 0.4) P2_1=P2_2=P2_3=0.;
	

	double P_E0 = P0 + alphas_mu/4./pi * (P1_1 + P1_2) + alphas_mu/4./pi*alphas_mu/4./pi *(P2_1 + P2_2 + P2_3);


	
	double Cem[9],C0w[9],C1w[9],C2w[9];
	CW_calculator(C0w,C1w,C2w,mu_W,param);

	double eta_mu=alphas_running(mu_W,param->mass_top_pole,param->mass_b_pole,param)/alphas_mu;

	double r=param->mass_b/param->mass_b_1S;
	
	double Kc0=1.4107*pow(eta_mu,14./23.)-0.8380*pow(eta_mu,16./23.)
	-0.4286*pow(eta_mu,6./23.)-0.0714*pow(eta_mu,-12./23.)
	-0.6494*pow(eta_mu,0.4086)-0.0380*pow(eta_mu,-0.4230)
	-0.0185*pow(eta_mu,-0.8994)-0.0057*pow(eta_mu,0.1456);

	double Kt0=pow(eta_mu,4./23.)*(C0w[7]+23./36.)-8./3.*(pow(eta_mu,4./23.)-pow(eta_mu,2./23.))*(C0w[8]+1./3.);

	double N_E0= -1./18.*(Kc0+r*Kt0)*(pow(eta_mu,6./23.)+pow(eta_mu,-12./23.))*param->lambda2/param->mass_c/param->mass_c;
	
	/* if((P2_1==0.)&&(P2_2==0.)&&(P2_3==0.)) N_E0=0.; */



	Cem[2]= -190./8073.*pow(eta_mu,-35./23.)-359./3105. *pow(eta_mu,-17./23.)+4276./121095. *pow(eta_mu,-12./23.) +350531./1009125. *pow(eta_mu,-9./23.)+2./4347. *pow(eta_mu,-7./23.)-5956./15525.*pow(eta_mu,6./23.) +38380./169533. *pow(eta_mu,14./23.)-748./8625. *pow(eta_mu,16./23.);
	
	Cem[8]= -32./575.*pow(eta_mu,-9./23.)+32./1449.*pow(eta_mu,-7./23.)+640./1449.*pow(eta_mu,14./23.)-704./1725.*pow(eta_mu,16./23.);
			
	Cem[7]= (32./75.*pow(eta_mu,-9./23.)-40./69.*pow(eta_mu,-7/23.)+88./575.*pow(eta_mu,16./23.))*C0w[7]+Cem[8]*C0w[8]+Cem[2];


	double P_em=alpha_em*(2.*Cem[7]*C0[7]-C0[7]*C0[7]*2.*alphas_mu*log(param->mass_W/mu)/pi)/alphas_mu;

#ifdef DEBUG
	for (ie=1;ie<=8;ie++) printf("C0[%d]=%f\t C1[%d]=%f\n",ie,C0[ie],ie,C1[ie]);
	printf("C2[7]=%f\n\n",C2[7]);
#endif

	double BRinc=BR_BXcenu_exp*pow(cabs(conj(param->Vts)*param->Vtb/param->Vcb),2.)*6.*alpha_em/pi/Cbr*(P_E0+P_em+N_E0);

	return BRinc;
}

/*---------------------------------------------------------------------*/

double bsgamma_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calculating the inclusive branching ratio of b -> s gamma */
{
	double C0w[11],C1w[11],C2w[11],C0b[11],C1b[11],C2b[11];
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return 0.;
	
	double mu_W=2.*param.mass_W;
		
	double mu_b=param.mass_b_1S/2.;

	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b,C1b,C2b,mu_b,&param);
	return bsgamma(C0b,C1b,C2b,mu_b,mu_W,&param);
}
