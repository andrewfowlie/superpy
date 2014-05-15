#include "include.h"

/*--------------------------------------------------------------------*/

double Bplus(double x, double y)
{
	return y/(x-y)*(log(y)/(y-1.)-log(x)/(x-1.));
}

/*---------------------------------------------------------------------*/

double D2(double x, double y)
{
	if(fabs(x-y)<1.e-5)
	{
		if(fabs(x-1.)<1.e-5) return -0.5;	
		return (1.-x+log(x))/(1.-x)/(1.-x);
	}

	return (D3(x)-D3(y))/(x-y);
}

/*---------------------------------------------------------------------*/

double D3(double x)
{
	if(fabs(x)<1.e-5) return 0.;
	if(fabs(x-1.)<1.e-5) return -1.;
	
	return x*log(x)/(1.-x);
}

/*--------------------------------------------------------------------*/

double h10(double x)
{
	if(fabs(1.-x)<1.e-5) return h10(0.9999);
	
	return (3.*x*x-2.*x)/3./pow(x-1.,4.)*log(x) 
	-(8.*x*x+5.*x-7.)/18./pow(x-1.,3.);
}

/*--------------------------------------------------------------------*/

double h20(double x)
{
	if(fabs(1.-x)<1.e-5) return h20(0.9999);
	
	return (-6.*x*x+4.*x)/3./pow(x-1.,3.)*log(x) 
	+(7.*x-5.)/3./pow(x-1.,2.);
}

/*--------------------------------------------------------------------*/

double h30(double x)
{
	if(fabs(1.-x)<1.e-5) return h30(0.9999);
	
	return (-6.*x*x*x+9.*x*x-2)/9./pow(x-1.,4.)*log(x) 
	+(52.*x*x-101.*x+43.)/54./pow(x-1.,3.);
}

/*----------------------------------------------------------------------*/

double h40(double x)
{
	if(fabs(1.-x)<1.e-5) return h40(0.9999);
	
	return -log(x)/3./pow(x-1.,4.)+(2.*x*x-7.*x+11.)/18./pow(x-1.,3.);
}

/*----------------------------------------------------------------------*/

double h50(double x)
{
	if(fabs(1.-x)<1.e-5) return h50(0.9999);
	
	return -x/pow(x-1.,4.)*log(x)+(-x*x+5.*x+2.)/6./pow(x-1.,3.);
}

/*----------------------------------------------------------------------*/

double h60(double x)
{
	if(fabs(1.-x)<1.e-5) return h60(0.9999);

	return 2.*x/pow(x-1.,3.)*log(x)-(x+1.)/pow(x-1.,2.);
}

/*----------------------------------------------------------------------*/

double f20(double x)
{
	if(fabs(1.-x)<1.e-5) return f20(0.9999);

	return -x/(x-1.)*(1.-1./(x-1.)*log(x));
}

/*----------------------------------------------------------------------*/

double f30(double x, double y)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f30(0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f30(0.9999,y);
	if(fabs(1.-y)<1.e-5) return f30(x,0.9999);

	if(fabs(1.-x/y)<1.e-5) return f30(y*0.9998,y);

	return x*log(x)/(x-1.)/(x-y)+y*log(y)/(y-1.)/(y-x);
}

/*----------------------------------------------------------------------*/

double f40(double x, double y)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f40(0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f40(0.9999,y);
	if(fabs(1.-y)<1.e-5) return f40(x,0.9999);

	if(fabs(1.-x/y)<1.e-5) return f40(y*0.9998,y);
	
	return x*x*log(x)/(x-1.)/(x-y)+y*y*log(y)/(y-1.)/(y-x);
}

/*----------------------------------------------------------------------*/

double f50(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f50(0.9996,0.9998,1.0002);

	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f50(0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f50(0.9998,y,1.0002);
	if((fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f50(x,0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f50(0.9999,y,z);
	if(fabs(1.-y)<1.e-5) return f50(x,0.9999,z);
	if(fabs(1.-z)<1.e-5) return f50(x,y,0.9999);
 
	if(fabs(1.-x/y)<1.e-5) return f50(y*0.9998,y,z);
	if(fabs(1.-x/y)<1.e-5) return f50(y*0.9998,y,z);
	if(fabs(1.-y/z)<1.e-5) return f50(x,z*0.9998,z);
	if(fabs(1.-x/z)<1.e-5) return f50(x,y,x*0.9998);

	return x*x*log(x)/(x-1.)/(x-y)/(x-z)+y*y*log(y)/(y-1.)/(y-x)/(y-z)
	+z*z*log(z)/(z-1.)/(z-x)/(z-y);
}

/*----------------------------------------------------------------------*/

double f60(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f60(0.9996,0.9998,1.0002);

	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f60(0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f60(0.9998,y,1.0002);
	if((fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f60(x,0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f60(0.9999,y,z);
	if(fabs(1.-y)<1.e-5) return f60(x,0.9999,z);
	if(fabs(1.-z)<1.e-5) return f60(x,y,0.9999);

	if(fabs(1.-x/y)<1.e-5) return f60(y*0.9998,y,z);
	if(fabs(1.-y/z)<1.e-5) return f60(x,z*0.9998,z);
	if(fabs(1.-x/z)<1.e-5) return f60(x,y,x*0.9998); 
	
	return x*log(x)/(x-1.)/(x-y)/(x-z)+y*log(y)/(y-1.)/(y-x)/(y-z)
	+z*log(z)/(z-1.)/(z-x)/(z-y);
}

/*----------------------------------------------------------------------*/

double f70(double x, double y)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f70(0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f70(0.9999,y);
	if(fabs(1.-y)<1.e-5) return f70(x,0.9999);

	if(fabs(1.-x/y)<1.e-5) return f70(y*0.9998,y);
	
	return x*log(x)/(x-1.)/(x-y)+x*log(y)/(y-1.)/(y-x);
}

/*----------------------------------------------------------------------*/

double f80(double x)
{
	if(fabs(x-1.)<1.e-5) return 1.;
	
	return x*log(x)/(x-1.);
}

/*----------------------------------------------------------------------*/

double f90(double w, double x, double y, double z)
{
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(0.9996,0.9998,1.0002,1.0004);

	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f90(0.9996,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(w,0.9996,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(0.9996,x,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(0.9996,0.9998,y,1.0002);

	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)) return f90(0.9998,1.0002,y,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f90(0.9998,x,1.0002,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(0.9998,x,y,1.0002);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f90(w,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(w,0.9998,y,1.0002);
	if((fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f90(w,x,0.9998,1.0002);
	
	if(fabs(1.-w)<1.e-5) return f90(0.9999,x,y,z);
	if(fabs(1.-x)<1.e-5) return f90(w,0.9999,y,z);
	if(fabs(1.-y)<1.e-5) return f90(w,x,0.9999,z);
	if(fabs(1.-z)<1.e-5) return f90(w,x,y,0.9999);

	if(fabs(1.-w/x)<1.e-5) return f90(x*0.9998,x,y,z);
	if(fabs(1.-w/y)<1.e-5) return f90(y*0.9998,x,y,z);
	if(fabs(1.-w/z)<1.e-5) return f90(z*0.9998,x,y,z);
	if(fabs(1.-x/y)<1.e-5) return f90(w,y*0.9998,y,z);
	if(fabs(1.-y/z)<1.e-5) return f90(w,x,z*0.9998,z);
	if(fabs(1.-x/z)<1.e-5) return f90(w,x,y,x*0.9998); 


	return w*w*log(w)/(w-1.)/(w-x)/(w-y)/(w-z) +x*x*log(x)/(x-1.)/(x-w)/(x-y)/(x-z) 
	+y*y*log(y)/(y-1.)/(y-x)/(y-w)/(y-z) +z*z*log(z)/(z-1.)/(z-x)/(z-y)/(z-w);
}

/*----------------------------------------------------------------------*/

double f100(double w, double x, double y, double z)
{
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(0.9996,0.9998,1.0002,1.0004);

	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f100(0.9996,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(w,0.9996,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(0.9996,x,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(0.9996,0.9998,y,1.0002);

	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)) return f100(0.9998,1.0002,y,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f100(0.9998,x,1.0002,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(0.9998,x,y,1.0002);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f100(w,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(w,0.9998,y,1.0002);
	if((fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return f100(w,x,0.9998,1.0002);

	if(fabs(1.-w)<1.e-5) return f100(0.9999,x,y,z);
	if(fabs(1.-x)<1.e-5) return f100(w,0.9999,y,z);
	if(fabs(1.-y)<1.e-5) return f100(w,x,0.9999,z);
	if(fabs(1.-z)<1.e-5) return f100(w,x,y,0.9999);

	if(fabs(1.-w/x)<1.e-5) return f100(x*0.9998,x,y,z);
	if(fabs(1.-w/y)<1.e-5) return f100(y*0.9998,x,y,z);
	if(fabs(1.-w/z)<1.e-5) return f100(z*0.9998,x,y,z);
	if(fabs(1.-x/y)<1.e-5) return f100(w,y*0.9998,y,z);
	if(fabs(1.-y/z)<1.e-5) return f100(w,x,z*0.9998,z);
	if(fabs(1.-x/z)<1.e-5) return f100(w,x,y,x*0.9998); 
	
	return w*log(w)/(w-1.)/(w-x)/(w-y)/(w-z) +x*log(x)/(x-1.)/(x-w)/(x-y)/(x-z) 
	+y*log(y)/(y-1.)/(y-x)/(y-w)/(y-z) +z*log(z)/(z-1.)/(z-x)/(z-y)/(z-w);
}

/*----------------------------------------------------------------------*/

double f110(double x, double y)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return f110(0.9998,1.0002);

	if(fabs(1.-x)<1.e-5) return f110(0.9999,y);
	if(fabs(1.-y)<1.e-5) return f110(x,0.9999);

	if(fabs(1.-x/y)<1.e-5) return f110(y*0.9998,y);
	
	return x*log(x)/(x-y)+x*log(y)/(y-x);
}

/*----------------------------------------------------------------------*/

double h11(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h11(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h11(0.99,y);
	if(fabs(1.-y)<1.e-3) return h11(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h11(y*0.98,y);

	return ((-48.*x*x*x-104.*x*x+64.*x)*Li2(1.-1./x)
	+(-378.*x*x*x-1566.*x*x+850.*x+86.)/9./(x-1.)*log(x)
	+(2060.*x*x*x+3798.*x*x-2664.*x-170.)/27.
	+((12.*x*x*x-124.*x*x+64.*x)/(x-1.)*log(x)+(-56.*x*x*x+258.*x*x+24.*x-82.)/3.)*y)/9./pow(x-1.,4.);
}

/*----------------------------------------------------------------------*/

double h21(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h21(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h21(0.99,y);
	if(fabs(1.-y)<1.e-3) return h21(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h21(y*0.98,y);

	return ((224.*x*x-96.*x)*Li2(1.-1./x)
	+(-24.*x*x*x+352.*x*x-128.*x-32.)*log(x)/(x-1.)
	+(-340.*x*x+132.*x+40.)
	+((-24.*x*x*x+176.*x*x-80.*x)*log(x)/(x-1.)+(-28.*x*x-108.*x+64.))*y)/9./pow(x-1.,3.);
}

/*----------------------------------------------------------------------*/

double h31(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h31(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h31(0.99,y);
	if(fabs(1.-y)<1.e-3) return h31(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h31(y*0.98,y);

	return (32.*x*x*x+120.*x*x-384.*x+128.)*Li2(1.-1./x)/81./pow(x-1.,4.)
	+(-108.*x*x*x*x+1058.*x*x*x-898.*x*x-1098.*x+710.)/81./pow(x-1.,5.)*log(x)
	+(-304.*x*x*x-13686.*x*x+29076.*x-12062.)/729./pow(x-1.,4.)
	+((540.*x*x*x-972.*x*x+232.*x+56.)/81./pow(x-1.,5.)*log(x)
	+(-664.*x*x*x+54.*x*x+1944.*x-902.)/243./pow(x-1.,4.))*y;
}

/*----------------------------------------------------------------------*/

double h41(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h41(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h41(0.99,y);
	if(fabs(1.-y)<1.e-3) return h41(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h41(y*0.98,y);

	return (-562.*x*x*x+1101.*x*x-420.*x+101.)/54./pow(x-1.,4.)*Li2(1.-1./x)
	+(-562.*x*x*x+1604.*x*x-799.*x+429.)/54./pow(x-1.,5.)*log(x)
	+(17470.*x*x*x-47217.*x*x+31098.*x-13447.)/972./pow(x-1.,4.)
	+((89.*x+55.)*log(x)/27./pow(x-1.,5.)+(38.*x*x*x-135.*x*x+54.*x-821.)/162./pow(x-1.,4.))*y;
}

/*----------------------------------------------------------------------*/

double h51(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h51(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h51(0.99,y);
	if(fabs(1.-y)<1.e-3) return h51(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h51(y*0.98,y);

	return ((9.*x*x*x+46.*x*x+49.*x)*Li2(1.-1./x)/2.
	+(81.*x*x*x+594.*x*x+1270.*x+71.)*log(x)/18./(x-1.)
	+(-923.*x*x*x-3042.*x*x-6921.*x-1210.)/108.
	+((10.*x*x+38.*x)/(x-1.)*log(x)+(-7.*x*x*x+30.*x*x-141.*x-26.)/3.)*y)/3./pow(x-1.,4.);
}

/*----------------------------------------------------------------------*/

double h61(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h61(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h61(0.99,y);
	if(fabs(1.-y)<1.e-3) return h61(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h61(y*0.98,y);

	return ((-32.*x*x-24.*x)*Li2(1.-1./x)
	+(-52.*x*x-109.*x-7.)*log(x)/(x-1.)
	+(95.*x*x+180.*x+61.)/2.
	+((-20.*x*x-52.*x)/(x-1.)*log(x)+(-2.*x*x+60.*x+14.))*y)/3./pow(x-1.,3.);
}

/*----------------------------------------------------------------------*/

double h71(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return h71(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return h71(0.99,y);
	if(fabs(1.-y)<1.e-3) return h71(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return h71(y*0.98,y);
	
	return (-20.*x*x*x+60.*x*x-60.*x-20.)*Li2(1.-1./x)/27./pow(x-1.,4.)
	+(-60.*x*x+240.*x+4.)/81./pow(x-1.,4.)*log(x)
	+(132.*x*x-382.*x+186.)/81./pow(x-1.,3.)
	+(20.*log(x)/27./pow(x-1.,4.)+(-20.*x*x+70.*x-110.)/81./pow(x-1.,3.))*y;
}

/*----------------------------------------------------------------------*/

double f31(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f31(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f31(0.99,y);
	if(fabs(1.-y)<1.e-3) return f31(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f31(y*0.98,y);
	
	return -28.*y/3./(y-1.)/(x-y)+2.*x*(11.*x+3.*y)/3./(x-1.)/(x-y)/(x-y)*log(x)
	+2.*y*(25.*x-11.*x*y-11.*y-3.*y*y)/3./(y-1.)/(y-1.)/(x-y)/(x-y)*log(y)
	+4.*(1.+y)/(x-1.)/(y-1)*Li2(1.-1./y)+4.*(x+y)/(x-1.)/(x-y)*Li2(1.-x/y);
}

/*----------------------------------------------------------------------*/

double f41(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f41(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f41(0.99,y);
	if(fabs(1.-y)<1.e-3) return f41(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f41(y*0.98,y);
	
	return (59.*x*(1.-y)-y*(59.-3.*y))/6./(y-1.)/(x-y)
	+4.*x*(7.*x*x-3.*x*y+3.*y*y)/3./(x-1.)/(x-y)/(x-y)*log(x)+2.*log(y)*log(y)
	+4.*y*y*(18.*x-11.*x*y-11.*y+4.*y*y)/3./(y-1.)/(y-1.)/(x-y)/(x-y)*log(y)
	+4.*(1.+y*y)/(x-1.)/(y-1)*Li2(1.-1./y)+4.*(x*x+y*y)/(x-1.)/(x-y)*Li2(1.-x/y);
}

/*----------------------------------------------------------------------*/

double f51(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f51(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f51(0.99,y);
	if(fabs(1.-y)<1.e-3) return f51(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f51(y*0.98,y);

	return (-83.-27.*x*(y-1.)+27.*y)/6./(x-1.)/(y-1.)
	
	-4.*x*(1.+x*(12.+y)-y-6.*x*x)/3./(x-1.)/(x-1.)/(x-y)*log(x)
	+2.*(1.+6.*x*x*(y-1.)-3.*x*x*x*(y-1.)+x*(3.*y-4.))/3./(x-1.)/(x-1.)/(y-1.)/(x-y)*log(x)*log(x)
	-4.*y*(3.*x*x*(y-1.)+x*y*(3.-2.*y)+y*y*(y-2.))/3./(x-1.)/(y-1)/(x-y)/(x-y)*Li2(1.-x/y)
	-4.*(1.-3.*x-x*x*(3.-6.*y)-x*x*x)/3./(x-1.)/(y-1)/(x-y)*Li2(1.-1./x)
	
	-4.*y*(1.+y*(12.+x)-x-6.*y*y)/3./(y-1.)/(y-1.)/(y-x)*log(y)
	+2.*(1.+6.*y*y*(x-1.)-3.*y*y*y*(x-1.)+y*(3.*x-4.))/3./(y-1.)/(y-1.)/(x-1.)/(y-x)*log(y)*log(y)
	-4.*x*(3.*y*y*(x-1.)+y*x*(3.-2.*x)+x*x*(x-2.))/3./(y-1.)/(x-1)/(y-x)/(y-x)*Li2(1.-y/x)
	-4.*(1.-3.*y-y*y*(3.-6.*x)-y*y*y)/3./(y-1.)/(x-1)/(y-x)*Li2(1.-1./y)
	+4.*log(x)*(f40(x,y)+(f40(1.0001*x,y)-f40(0.9999*x,y))/0.0002+(f40(x,1.0001*y)-f40(x,0.9999*y))/0.0002);
}

/*----------------------------------------------------------------------*/

double f81(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f81(0.96,0.98,1.02);

	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f81(0.98,1.02,z);
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f81(0.98,y,1.02);
	if((fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f81(x,0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f81(0.99,y,z);
	if(fabs(1.-y)<1.e-3) return f81(x,0.99,z);
	if(fabs(1.-z)<1.e-3) return f81(x,y,0.99);

	if(fabs(1.-x/y)<5.e-3) return f81(y*0.98,y,z);
	if(fabs(1.-y/z)<5.e-3) return f81(x,z*0.98,z);
	if(fabs(1.-x/z)<5.e-3) return f81(x,y,x*0.98);

	return -28.*y*y/3./(y-1.)/(x-y)/(y-z)
	+4.*x*(7.*x*x-3.*x*y+3.*y*y)/3./(x-1.)/(x-y)/(x-y)/(x-z)*log(x)
	+4.*z*(7.*z*z-3.*z*y+3.*y*y)/3./(z-1.)/(z-y)/(z-y)/(z-x)*log(z)
	-4.*y*y*(x*(4.*y*y+18.*z-11.*y*(1.+z))+y*(3.*y*y-11.*z+4.*y*(1.+z)))/3./(y-1.)/(y-1.)/(x-y)/(x-y)/(y-z)/(y-z)*log(y)
	-4.*(1.+y*y)/(x-1.)/(y-1.)/(z-1.)*Li2(1.-1./y)
	+4.*(x*x+y*y)/(x-1.)/(x-y)/(x-z)*Li2(1.-x/y)
	+4.*(z*z+y*y)/(z-1.)/(z-y)/(z-x)*Li2(1.-z/y);
}

/*----------------------------------------------------------------------*/

double f91(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f91(0.96,0.98,1.02);

	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f91(0.98,1.02,z);
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f91(0.98,y,1.02);
	if((fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f91(x,0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f91(0.99,y,z);
	if(fabs(1.-y)<1.e-3) return f91(x,0.99,z);
	if(fabs(1.-z)<1.e-3) return f91(x,y,0.99);

	if(fabs(1.-x/y)<5.e-3) return f91(y*0.98,y,z);
	if(fabs(1.-y/z)<5.e-3) return f91(x,z*0.98,z);
	if(fabs(1.-x/z)<5.e-3) return f91(x,y,x*0.98);
		
	return -28.*y/3./(y-1.)/(x-y)/(y-z)
	+2.*x*(11.*x+3.*y)/3./(x-1.)/(x-y)/(x-y)/(x-z)*log(x)
	+2.*z*(11.*z+3.*y)/3./(z-1.)/(z-y)/(z-y)/(z-x)*log(z)
	+2.*y*(x*(3.*y*y-25.*z+11.*y*(1.+x))+y*(-17.*y*y+11.*z+3.*y*(1.+z)))/3./(y-1.)/(y-1.)/(x-y)/(x-y)/(y-z)/(y-z)*log(y)
	-4.*(1.+y)/(x-1.)/(y-1.)/(z-1.)*Li2(1.-1./y)
	+4.*(x+y)/(x-1.)/(x-y)/(x-z)*Li2(1.-x/y)
	+4.*(z+y)/(z-1.)/(z-y)/(z-x)*Li2(1.-z/y);
}

/*----------------------------------------------------------------------*/

double f111(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f111(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f111(0.99,y);
	if(fabs(1.-y)<1.e-3) return f111(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f111(y*0.98,y);
	
	return 4.*x*(8.*y+(x-1.)*(x-y)*pi*pi)/3./y/(x-1.)/(x-y)
	-8.*x*(x*x-7.*y+3.*x*(1.+y))/3./(x-y)/(x-y)/(x-1.)/(x-1.)*log(x)
	-8.*x*(3.*x-7.*y)/3./(x-y)/(x-y)/(y-1.)*log(y)
	-8.*x/(y-1.)*Li2(1.-1./x)
	+8.*x/y/(y-1.)*Li2(1.-y/x);
}

/*----------------------------------------------------------------------*/

double f121(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f121(0.96,0.98,1.02);

	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f121(0.98,1.02,z);
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f121(0.98,y,1.02);
	if((fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f121(x,0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f121(0.99,y,z);
	if(fabs(1.-y)<1.e-3) return f121(x,0.99,z);
	if(fabs(1.-z)<1.e-3) return f121(x,y,0.99);

	if(fabs(1.-x/y)<5.e-3) return f121(y*0.98,y,z);
	if(fabs(1.-y/z)<5.e-3) return f121(x,z*0.98,z);
	if(fabs(1.-x/z)<5.e-3) return f121(x,y,x*0.98);
	
	return -28.*y*y/3./(x-y)/(y-1.)/(y-z)
	+4.*x*x*(6.*x+y)/3./(x-1.)/(x-y)/(x-y)/(x-z)*log(x)
	+4.*z*z*(6.*z+y)/3./(z-1.)/(z-y)/(z-y)/(z-x)*log(z)
	-4.*y*y*(x*(6.*y*y+20.*z-13.*y*(1.+z))+y*(y*y-13.*z+6.*y*(1.+z)))/3./(x-y)/(x-y)/(y-1.)/(y-1.)/(y-z)/(y-z)*log(y);
	
}

/*----------------------------------------------------------------------*/

double f131(double x, double y, double z)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f131(0.96,0.98,1.02);

	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f131(0.98,1.02,z);
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f131(0.98,y,1.02);
	if((fabs(1.-y)<1.e-3)&&(fabs(1.-z)<1.e-3)) return f131(x,0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f131(0.99,y,z);
	if(fabs(1.-y)<1.e-3) return f131(x,0.99,z);
	if(fabs(1.-z)<1.e-3) return f131(x,y,0.99);

	if(fabs(1.-x/y)<5.e-3) return f131(y*0.98,y,z);
	if(fabs(1.-y/z)<5.e-3) return f131(x,z*0.98,z);
	if(fabs(1.-x/z)<5.e-3) return f131(x,y,x*0.98);
	
	return -28.*y/3./(x-y)/(y-1.)/(y-z)
	+4.*x*(6.*x+y)/3./(x-1.)/(x-y)/(x-y)/(x-z)*log(x)
	+4.*z*(6.*z+y)/3./(z-1.)/(z-y)/(z-y)/(z-x)*log(z)
	+4.*y*(x*(y*y-13.*z+6.*y*(1.+z))+y*(y-8.*y*y+6.*z+y*z))/3./(x-y)/(x-y)/(y-1.)/(y-1.)/(y-z)/(y-z)*log(y);
}

/*----------------------------------------------------------------------*/

double f141(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f141(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f141(0.99,y);
	if(fabs(1.-y)<1.e-3) return f141(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f141(y*0.98,y);
	
	return 32.*x*x/3./(x-1.)/(x-y)
	-8.*x*x*(7.*x*(1.+y)-11.*y-3.*x*x)/3./(x-1.)/(x-1.)/(x-y)/(x-y)*log(x)
	-8.*x*y*(3.*x-7.*y)/3./(x-y)/(x-y)/(y-1.)*log(y)
	-8.*x/(y-1.)*Li2(1.-1./x)
	+8.*x/(y-1.)*Li2(1.-y/x);
}

/*----------------------------------------------------------------------*/

double f151(double x)
{
	if(fabs(1.-x)<1.e-3) return f151(0.99);

	return (1.-3.*x)/(x-1.)+2.*x/(x-1.)/(x-1.)*log(x)+2.*x/(x-1.)*Li2(1.-1./x);
}


/*----------------------------------------------------------------------*/

double f161(double x)
{
	if(fabs(1.-x)<1.e-3) return f161(0.99);

	return 28./3./(x-1.)-4.*x*(13.-6.*x)/3./(x-1.)/(x-1.)*log(x);
}

/*----------------------------------------------------------------------*/

double f171(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f171(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f171(0.99,y);
	if(fabs(1.-y)<1.e-3) return f171(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f171(y*0.98,y);
	
	return -28./3./(x-1.)/(y-1.)
	+4.*y*(10.-3.*y)/3./(x-y)/(y-1.)/(y-1.)*log(y)
	-4.*y/(x-y)/(y-1.)/(y-1.)*log(y)*log(y)
	+(4.*(13.*x-6.*x*x-3.*y-7.*x*y+3.*x*x*y)/3./(x-1.)/(x-1.)/(x-y)/(y-1.)+4.*y*log(y)/(x-y)/(y-1.)/(y-1.))*log(x);
}

/*----------------------------------------------------------------------*/

double f181(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f181(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f181(0.99,y);
	if(fabs(1.-y)<1.e-3) return f181(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f181(y*0.98,y);
	
	return -28.*y/3./(x-y)/(y-1.)
	+4.*x*(6.*x+y)/3./(x-1.)/(x-y)/(x-y)*log(x)
	-4.*y*(y*(6.+y)-x*(13.-6.*y))/3./(x-y)/(x-y)/(y-1.)/(y-1.)*log(y);
}

/*----------------------------------------------------------------------*/

double f191(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return f191(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return f191(0.99,y);
	if(fabs(1.-y)<1.e-3) return f191(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return f191(y*0.98,y);
	
	return -28.*(x*(y-1.)+y)/3./(x-y)/(y-1.)
	+4.*x*x*(6.*x+y)/3./(x-1.)/(x-y)/(x-y)*log(x)
	+4.*y*y*(x*(20.-13.*y)-y*(13.-6.*y))/3./(x-y)/(x-y)/(y-1.)/(y-1.)*log(y);
}

/*----------------------------------------------------------------------*/

double q11(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q11(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q11(0.99,y);
	if(fabs(1.-y)<1.e-3) return q11(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q11(y*0.98,y);

	return 4./3./(x-y)*(x*x*log(x)/pow(x-1.,4.)-y*y*log(y)/pow(y-1.,4.)) +(4.*x*x*y*y+10.*x*y*y-2.*y*y+10.*x*x*y-44.*x*y+10.*y-2.*x*x+10.*x+4.)/9./pow(x-1.,3.)/pow(y-1.,3.);
}

/*----------------------------------------------------------------------*/

double q21(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q21(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q21(0.99,y);
	if(fabs(1.-y)<1.e-3) return q21(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q21(y*0.98,y);
	
	return 4./3./(x-y)*(x*log(x)/pow(x-1.,4.)-y*log(y)/pow(y-1.,4.)) +(-2.*x*x*y*y+10.*x*y*y+4.*y*y+10.*x*x*y-20.*x*y-14.*y+4.*x*x-14.*x+22.)/9./pow(x-1.,3.)/pow(y-1.,3.);
}

/*----------------------------------------------------------------------*/

double q31(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q31(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q31(0.99,y);
	if(fabs(1.-y)<1.e-3) return q31(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q31(y*0.98,y);

	return 8./3./(x-y)*(-x*x*log(x)/pow(x-1.,3.)+y*y*log(y)/pow(y-1.,3.)) +(-12.*x*y+4.*y+4.*x+4.)/3./pow(x-1.,2.)/pow(y-1.,2.);
}

/*----------------------------------------------------------------------*/

double q41(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q41(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q41(0.99,y);
	if(fabs(1.-y)<1.e-3) return q41(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q41(y*0.98,y);

	return 8./3./(x-y)*(-x*log(x)/pow(x-1.,3.)+y*log(y)/pow(y-1.,3.)) +(-4.*x*y-4.*y-4.*x+12.)/3./pow(x-1.,2.)/pow(y-1.,2.);
}

/*----------------------------------------------------------------------*/

double q51(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q51(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q51(0.99,y);
	if(fabs(1.-y)<1.e-3) return q51(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q51(y*0.98,y);
	
	return 4./27./(x-y)*((6.*x*x*x-9.*x*x+2.)*log(x)/pow(x-1.,4.)-(6.*y*y*y-9.*y*y+2.)*log(y)/pow(y-1.,4.))
	+(104.*x*x*y*y-202.*x*y*y+86.*y*y-202.*x*x*y+380.*x*y-154.*y+86.*x*x-154.*x+56.)/81./pow(x-1.,3.)/pow(y-1.,3.);
}

/*----------------------------------------------------------------------*/

double q61(double x, double y)
{
	if((fabs(1.-x)<1.e-3)&&(fabs(1.-y)<1.e-3)) return q61(0.98,1.02);

	if(fabs(1.-x)<1.e-3) return q61(0.99,y);
	if(fabs(1.-y)<1.e-3) return q61(x,0.99);

	if(fabs(1.-x/y)<5.e-3) return q61(y*0.98,y);

	return 4./9./(x-y)*(log(x)/pow(x-1.,4.)-log(y)/pow(y-1.,4.))
	+(4.*x*x*y*y-14.*x*y*y+22.*y*y-14.*x*x*y+52.*x*y-62.*y+22.*x*x-62.*x+52.)/27./pow(x-1.,3.)/pow(y-1.,3.);
}

/*----------------------------------------------------------------------*/

double A0t(double x)
{
	return (-46.+159.*x-153.*x*x+22.*x*x*x)/36./pow(1.-x,3.)+x*x*(-3.*x+2.)/2./pow(1.-x,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double A1t(double x, double l)
{
	return (32.*pow(x,4.)+244.*pow(x,3.)-160.*x*x+16.*x)/9./pow(1.-x,4.)*Li2(1.-1./x)
	+(-774.*x*x*x*x-2826.*x*x*x+1994.*x*x-130.*x+8.)/81./pow(1.-x,5.)*log(x)
	+(-94.*x*x*x*x-18665.*x*x*x+20682.*x*x-9113.*x+2006.)/243./pow(1.-x,4.)
	+((-12.*x*x*x*x-92.*x*x*x+56.*x*x)/3./pow(1.-x,5.)*log(x)+(-68.*x*x*x*x-202.*x*x*x-804.*x*x+794.*x-152.)/27./pow(1.-x,4.))*l;
}

/*--------------------------------------------------------------------*/

double B0t(double x)
{
	return x/4./pow(1.-x,2.)*log(x)+1./4./(1.-x);
}

/*----------------------------------------------------------------------*/

double C0t(double x)
{
	return (3.*x+2.)*x/8./pow(1.-x,2.)*log(x)+(-x+6.)*x/8./(1.-x);
}

/*----------------------------------------------------------------------*/

double D0t(double x)
{
	return (-3.*pow(x,4.)+30.*pow(x,3.)-54.*x*x+32.*x-8.)/18./pow(1.-x,4.)*log(x)+(-47.*pow(x,3.)+237.*x*x-312.*x+104.)/108./pow(1.-x,3.);
}

/*----------------------------------------------------------------------*/

double B1t(double x, double l)
{
	return -2.*x/(1.-x)/(1.-x)*Li2(1.-1./x) +(-x+17.)*x/3./pow(1.-x,3.)*log(x) +(13.*x+3)/3./(1.-x)/(1.-x) +((2.*x+2)*x/pow(1.-x,3.)*log(x)+4.*x/(1.-x)/(1.-x))*l;
}

/*----------------------------------------------------------------------*/

double C1t(double x, double l)
{
	return (-x*x-4.)*x/(1.-x)/(1.-x)*Li2(1.-1./x) +(3.*x*x+14.*x+23.)*x/3./pow(1.-x,3.)*log(x) +(4.*x*x+7.*x+29.)*x/3./(1.-x)/(1.-x) +((8.*x+2.)*x/pow(1.-x,3.)*log(x)+(x*x+x+8.)*x/(1.-x)/(1.-x))*l;
}

/*----------------------------------------------------------------------*/

double D1t(double x, double l)
{
	return (380.*pow(x,4.)-1352.*pow(x,3.)+1656.*x*x-784.*x+256.)/81./pow(1.-x,4.)*Li2(1.-1./x) +(304.*pow(x,4.)+1716.*pow(x,3.)-4644.*x*x+2768.*x-720.)/81./pow(1.-x,5.)*log(x) +(-6175.*pow(x,4.)+41608.*pow(x,3.)-66723.*x*x+33106.*x-7000.)/729./pow(1.-x,4.) +((648.*pow(x,4.)-720.*pow(x,3.)-232.*x*x-160.*x+32.)/81./pow(1.-x,5.)*log(x)
	+(-352.*pow(x,4.)+4912.*pow(x,3.)-8280.*x*x+3304.*x-880.)/243./pow(1.-x,4.))*l;	
}

/*----------------------------------------------------------------------*/

double F0t(double x)
{
	return  (5.*x*x*x-9.*x*x+30.*x-8.)/12./pow(1.-x,3.)+3.*x*x/2./pow(1.-x,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double F1t(double x,double l)
{
	return (4.*pow(x,4.)-40.*pow(x,3.)-41.*x*x-x)/3./pow(1.-x,4.)*Li2(1.-1./x)
	+(-144.*x*x*x*x+3177.*x*x*x+3661.*x*x+250.*x-32.)/108./pow(1.-x,5.)*log(x)
	+(-247.*x*x*x*x+11890.*x*x*x+31779.*x*x-2966.*x+1016.)/648./pow(1.-x,4.)
	+((17.*x*x*x+31.*x*x)/pow(1.-x,5.)*log(x)+(-35.*x*x*x*x+170.*x*x*x+447.*x*x+338.*x-56.)/18./pow(1.-x,4.))*l;
}

/*----------------------------------------------------------------------*/

double E0t(double x)
{
	return (-9.*x*x+16.*x-4.)/6./pow(1.-x,4.)*log(x)+(-7.*x*x*x-21.*x*x+42.*x+4.)/36./pow(1.-x,3.);
}

/*----------------------------------------------------------------------*/

double G1t(double x, double l)
{
	return (10.*pow(x,4.)-100.*pow(x,3.)+30.*x*x+160.*x-40.)/27./pow(1.-x,4.)*Li2(1.-1./x)
	+(30.*pow(x,3.)-42.*x*x-332.*x+68.)/81./pow(1.-x,4.)*log(x)
	+(-6.*pow(x,3.)-293.*x*x+161.*x+42.)/81./pow(1.-x,3.)
	+((90.*x*x-160.*x+40.)/27./pow(1.-x,4.)*log(x)+(35.*pow(x,3.)+105.*x*x-210.*x-20.)/81./pow(1.-x,3.))*l;
}

/*----------------------------------------------------------------------*/

double E1t(double x, double l)
{
	return (515.*pow(x,4.)-614.*pow(x,3.)-81.*x*x-190.*x+40.)/54./pow(1.-x,4.)*Li2(1.-1./x)
	+(-1030.*pow(x,4.)+435.*pow(x,3.)+1373.*x*x+1950.*x-424.)/108./pow(1.-x,5.)*log(x)
	+(-29467.*pow(x,4.)+45604.*pow(x,3.)-30237.*x*x+66532.*x-10960.)/1944./pow(1.-x,4.)
	+((-1125.*pow(x,3.)+1685.*x*x+380.*x-76.)/54./pow(1.-x,5.)*log(x)
	+(133.*pow(x,4.)-2758.*pow(x,3.)-2061.*x*x+11522.*x-1652.)/324./pow(1.-x,4.))*l;
}

/*----------------------------------------------------------------------*/

double T(double x)
{
	return -(16.*x+8.)*sqrt(4.*x-1.)*Cl2(2.*asin(0.5/sqrt(x)))+(16.*x+20./3.)*log(x)+32.*x+112./9.;
}

/*----------------------------------------------------------------------*/

double F7_1(double x)
{
	return x*(7.-5.*x-8.*x*x)/24./pow(x-1.,3.)+x*x*(3.*x-2.)/4./pow(x-1.,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double F7_2(double x)
{
	return x*(3.-5.*x)/12./pow(x-1.,2.)+x*(3.*x-2.)/6./pow(x-1.,3.)*log(x);
}

/*----------------------------------------------------------------------*/

double F8_1(double x)
{
	return x*(2.+5.*x-x*x)/8./pow(x-1.,3.)-3.*x*x/4./pow(x-1.,4.)*log(x);
}

/*----------------------------------------------------------------------*/

double F8_2(double x)
{
	return x*(3.-x)/4./pow(x-1.,2.)-x/2./pow(x-1.,3.)*log(x);
}

/*----------------------------------------------------------------------*/

double H2(double x, double y)
{
	return D2(x,y);
}

/*----------------------------------------------------------------------*/

double B(double m1, double m2, double Q)
{
	double x=pow(m2/m1,2.);

	if(fabs(x-1.)<1.e-5) return -0.5*log(m2*m2/Q/Q);

	return 0.5*(0.5+1./(1.-x)+log(x)/pow(1.-x,2.)-log(m2*m2/Q/Q));
}

/*----------------------------------------------------------------------*/

double G7H(double x, double lu, double ld)
{
	return lu*(ld*4.*x*(4.*(-3.+7.*x-2.*x*x)*Li2(1.-1./x)+(8.-14.*x-3.*x*x)*log(x)*log(x)/(x-1.)
	+2.*(-3.-x+12.*x*x-2.*x*x*x)*log(x)/(x-1.)+3.*(7.-13.*x+2.*x*x))/9./pow((x-1.),3.)
	+lu*2.*x*(x*(18.-37.*x+8.*x*x)*Li2(1.-1./x)+x*(-14.+23.*x+3.*x*x)*log(x)*log(x)/(x-1.)+
	(-50.+251.*x-174.*x*x-192.*x*x*x+21.*x*x*x*x)*log(x)/9./(x-1.)+(797.-5436.*x+7569.*x*x-1202.*x*x*x)/108.)/pow((x-1.),4.)/9.);
}

/*----------------------------------------------------------------------*/
	
double Delta7H(double x, double lu, double ld)
{
	return 2.*x/9./pow((x-1.),4.)*lu*
	(ld*((x-1.)*(21.-47.*x+8.*x*x)+2.*(-8.+14.*x+3.*x*x)*log(x))
	+lu*((-31.-18.*x+135.*x*x-14.*x*x*x)/6.+x*(14.-23.*x-3.*x*x)*log(x)/(x-1.)));
}

/*----------------------------------------------------------------------*/

double G8H(double x, double lu, double ld)
{
	return lu*(ld*x*(0.5*(-36.+25.*x-17.*x*x)*Li2(1.-1./x)+(19.+17.*x)*log(x)*log(x)/(x-1.)
	+0.25*(-3.-187.*x+12.*x*x-14.*x*x*x)*log(x)/(x-1.)+3.*(143.-44.*x+29.*x*x)/8.)/3./pow((x-1.),3.)
	+ lu*x*(x*(30.-17.*x+13.*x*x)*Li2(1.-1./x)-x*(31.+17.*x)*log(x)*log(x)/(x-1.)+
	(-226.+817.*x+1353.*x*x+318.*x*x*x+42.*x*x*x*x)*log(x)/36./(x-1.)+(1130.-18153.*x+7650.*x*x-4451.*x*x*x)/216.)/pow((x-1.),4.)/6.);
}

/*----------------------------------------------------------------------*/

double Delta8H(double x, double lu, double ld)
{
	return (x/6./pow((x-1.),4.))*lu*(ld*((x-1.)*(81.-16.*x+7.*x*x)-2.*(19.+17.*x)*log(x))
	+lu*((-38.-261.*x+18.*x*x-7.*x*x*x)/6.+x*(31.+17.*x)*log(x)/(x-1.)));	
}

/*----------------------------------------------------------------------*/

double EH(double x, double lu)
{
	return lu*lu*x*((x-1.)*(16.-29.*x+7.*x*x)+6.*(3.*x-2.)*log(x))/36./pow((x-1.),4.);
}

/*----------------------------------------------------------------------*/

double G4H(double x, double lu)
{
	return lu*lu*((515.*x*x*x-906.*x*x+99.*x+182.)*x*Li2(1.-1./x)/54./pow(x-1.,4.)
	+(1030.*x*x*x-2763.*x*x-15.*x+980.)*x*log(x)/108./pow(x-1.,5.)
	+(-29467.*x*x*x+68142.*x*x-6717.*x-18134.)*x/1944./pow(x-1.,4.));
}

/*----------------------------------------------------------------------*/

double Delta4H(double x, double lu)
{
	return -lu*lu*((-375.*x*x-95.*x+182.)*x*log(x)/54./pow(x-1.,5.)
	+(133.*x*x*x-108.*x*x+4023.*x-2320.)*x/324./pow(x-1.,4.));
}

/*----------------------------------------------------------------------*/

double G3H(double x, double lu)
{
	return lu*lu*((10.*x*x*x+30.*x-20.)*x*Li2(1.-1./x)/27./pow(x-1.,4.)
	+(30.*x*x-66.*x-56.)*x*log(x)/81./pow(x-1.,4.)
	+(6.*x*x-187.*x+213.)*x/81./pow(x-1.,3.));
}

/*----------------------------------------------------------------------*/

double Delta3H(double x, double lu)
{
	return -lu*lu*((-30.*x+20.)*x*log(x)/27./pow(x-1.,4.)
	+(-35.*x*x+145.*x-80.)*x/81./pow(x-1.,3.));
}

/*----------------------------------------------------------------------*/

double C9llH0(double x, double y, double lu)
{
	return x/y/8.*lu*lu*(-log(y)/(y-1.)+1.)*y*y/(y-1.);
}

/*----------------------------------------------------------------------*/

double D9H0(double x, double lu)
{
	return lu*lu*((-3.*x*x*x+6.*x-4.)*x/18./pow(x-1.,4.)*log(x)+(47.*x*x-79.*x+38.)*x/108./pow(x-1.,3.));	
}

/*----------------------------------------------------------------------*/

double C9llH1(double x, double y, double lu, double L)
{
	if(fabs(y-1.)<1.e-5) return C9llH1(x,0.9999,lu,L);
	
 	return x/y/8.*lu*lu*((-8.*y*y*y+16.*y*y)/pow(y-1.,2.)*Li2(1.-1./y)
	+(-24.*y*y*y+88.*y*y)/3./pow(y-1.,3.)*log(y)
	+(32.*y*y*y-96.*y*y)/3./(y-1.)/(y-1.)
	+(16.*y*y*log(y)/pow(y-1.,3.)+(8.*y*y*y-24.*y*y)/(y-1.)/(y-1.))*L);
}

/*----------------------------------------------------------------------*/

double D9H1(double x, double lu, double L)
{
	if(fabs(x-1.)<1.e-5) return D9H1(0.9999,lu,L);
	
	return lu*lu*((380.*x*x*x-528.*x*x+72.*x+128.)*x/81.*pow(1.-x,4.)*Li2(1.-1./x)
	+(596.*x*x*x-672.*x*x+64.*x+204.)*x*log(x)/81./pow(x-1.,5.)
	+(-6175.*x*x*x+9138.*x*x-3927.*x-764.)*x/729./pow(x-1.,4.)
	+((432.*x*x*x-456.*x*x+40.*x+128.)*x*log(x)/81./pow(x-1.,5.)
	+(-352.*x*x*x-972.*x*x+1944.*x-1052.)*x/243./pow(x-1.,4.))*L);	
}

/*----------------------------------------------------------------------*/

double C7t2mt(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.4) return 12.06+12.93*z+3.013*z*log(z)+96.71*z*z+52.73*z*z*log(z)+147.9*pow(z,3.)
	+187.7*pow(z,3.)*log(z)-144.9*pow(z,4.)+236.1*pow(z,4.)*log(z);
	
	else return 11.74+0.3642*w+0.1155*w*w-0.003145*pow(w,3.)-0.03263*pow(w,4.)-0.03528*pow(w,5.)
	-0.03076*pow(w,6.)-0.02504*pow(w,7.)-0.01985*pow(w,8.);
}

/*----------------------------------------------------------------------*/

double C7c2MW(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.4) return 1.525-0.1165*z+0.01975*z*log(z)+0.06283*z*z+0.005349*z*z*log(z)
	+0.01005*pow(z*log(z),2.)-0.04202*pow(z,3.)+0.01535*pow(z,3.)*log(z)-0.00329*z*pow(z*log(z),2.)
	+0.002372*pow(z,4.)-0.0007910*pow(z,4.)*log(z);
		
	else return 1.432+0.06709*w+0.01257*w*w+0.004710*pow(w,3.)+0.002373*pow(w,4.)
	+0.001406*pow(w,5.)+0.0009216*pow(w,6.)+0.00064730*pow(w,7.)+0.0004779*pow(w,8.);
}

/*----------------------------------------------------------------------*/

double C8t2mt(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.35) return -0.8954-7.043*z-98.34*z*z-46.21*z*z*log(z)-127.1*pow(z,3.)
	-181.6*pow(z,3.)*log(z)+535.8*pow(z,4.)-76.76*pow(z,4.)*log(z);
	
	else return -0.6141-0.8975*w-0.03492*w*w+0.06791*pow(w,3.)+0.07966*pow(w,4.)
	+0.07226*pow(w,5.)+0.06132*pow(w,6.)+0.05096*pow(w,7.)+0.04216*pow(w,8.);
}

/*----------------------------------------------------------------------*/

double C8c2MW(double x)
{
	double z=1./x;
	double w=1.-z;
	double y=sqrt(z);
	
	if(y<0.35) return -1.870+0.1010*z-0.1218*z*log(z)+0.1045*z*z-0.03748*z*z*log(z)
	+0.01151*pow(z*log(z),2.)-0.01023*pow(z,3.)+0.004342*pow(z,3.)*log(z)+0.0003031*z*pow(z*log(z),2.)
	-0.001537*pow(z,4.)+0.0007532*pow(z,4.)*log(z);
	
	else return -1.676-0.1179*w-0.02926*w*w-0.01297*pow(w,3.)-0.007296*pow(w,4.)
	-0.004672*pow(w,5.)-0.003248*pow(w,6.)-0.002389*pow(w,7.)-0.001831*pow(w,8.);
}

/*----------------------------------------------------------------------*/

double epsilon_0(struct parameters* param)
{
	if(param->SM==1) return 0.;

#ifdef SM_ChargedHiggs
	return 0.;
#endif

 	double sw2=pow(sin(atan(param->gp/param->g2)),2.);	
	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);

	return 2./3.*alphas_MSOFT/pi*((param->A_b/param->tan_beta-param->mu_Q)/param->mass_gluino*
H2(param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino)	
	-0.5*(B(param->mass_gluino,param->mass_b1, param->MSOFT_Q)+B(param->mass_gluino,param->mass_b2, param->MSOFT_Q))/param->tan_beta )
	 +1./param->inv_alpha_em/sw2/4./pi*(param->mu_Q*param->M2_Q)*( +param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1,param->mu_Q*param->mu_Q/param->mass_b1/param->mass_b1)/param->mass_b1/param->mass_b1/2. +param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2,param->mu_Q*param->mu_Q/param->mass_b2/param->mass_b2)/param->mass_b2/param->mass_b2/2.);
}

/*----------------------------------------------------------------------*/

double epsilon_2(struct parameters* param)
{
	if(param->SM==1) return 0.;

#ifdef SM_ChargedHiggs
	return 0.;
#endif

 	double sw2=pow(sin(atan(param->gp/param->g2)),2.);	
	return param->yut[3]*param->yut[3]/16./pi/pi*(param->mu_Q/param->tan_beta-param->A_t)*(param->charg_Umix[1][2]*param->charg_Vmix[1][2]/param->mass_cha1 *H2(param->mass_t1*param->mass_t1/param->mass_cha1/param->mass_cha1,param->mass_t2*param->mass_t2/param->mass_cha1/param->mass_cha1) +param->charg_Umix[2][2]*param->charg_Vmix[2][2]/param->mass_cha2*H2(param->mass_t1*param->mass_t1/param->mass_cha2/param->mass_cha2,param->mass_t2*param->mass_t2/param->mass_cha2/param->mass_cha2))	+1./param->inv_alpha_em/sw2/4./pi*(param->mu_Q*param->M2_Q)*(param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1,param->mu_Q*param->mu_Q/param->mass_t1/param->mass_t1)/param->mass_t1/param->mass_t1 +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2,param->mu_Q*param->mu_Q/param->mass_t2/param->mass_t2)/param->mass_t2/param->mass_t2);
}

/*----------------------------------------------------------------------*/

double epsilon_b(struct parameters* param)
{
	return epsilon_0(param)+epsilon_2(param);
}

/*----------------------------------------------------------------------*/

double epsilon_bp(struct parameters* param)
{
	if(param->SM==1) return 0.;

#ifdef SM_ChargedHiggs
	return 0.;
#endif
 	double sw2=pow(sin(atan(param->gp/param->g2)),2.);
	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);
	int ie;
	int nb_neut;
	if(param->mass_neut[5]==0.) nb_neut=4; else nb_neut=5;

	double epsilonbp=2./3.*alphas_MSOFT/pi*(param->A_b/param->tan_beta-param->mu_Q)/param->mass_gluino*( param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_b2*param->mass_b2/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_b1*param->mass_b1/param->mass_gluino/param->mass_gluino));	
	for(ie=1;ie<=nb_neut;ie++)	epsilonbp+=param->yut[3]*param->yut[3]/16./pi/pi*param->neut_mix[ie][4]*param->neut_mix[ie][3]*(param->A_t-param->mu_Q/param->tan_beta)/param->mass_neut[ie]*( param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]));
	epsilonbp+=1./param->inv_alpha_em/sw2/4./pi*(param->mu_Q*param->M2_Q)*( (param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_t1/param->mass_t1,param->mu_Q*param->mu_Q/param->mass_t1/param->mass_t1)/param->mass_t1/param->mass_t1 +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_t2/param->mass_t2,param->mu_Q*param->mu_Q/param->mass_t2/param->mass_t2)/param->mass_t2/param->mass_t2)/2. +(param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->M2_Q*param->M2_Q/param->mass_b1/param->mass_b1,param->mu_Q*param->mu_Q/param->mass_b1/param->mass_b1)/param->mass_b1/param->mass_b1 +param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->M2_Q*param->M2_Q/param->mass_b2/param->mass_b2,param->mu_Q*param->mu_Q/param->mass_b2/param->mass_b2)/param->mass_b2/param->mass_b2));

	return epsilonbp;
}

/*----------------------------------------------------------------------*/

double epsilon_0p(struct parameters* param)
{
	if(param->SM==1) return 0.;

#ifdef SM_ChargedHiggs
	return 0.;
#endif

	double alphas_MSOFT=alphas_running(param->MSOFT_Q,param->mass_top_pole,param->mass_b_pole,param);
	int ie;
	int nb_neut;
	if(param->mass_neut[5]==0.) nb_neut=4; else nb_neut=5;
	
	double epsilon0p=-2./3.*alphas_MSOFT/pi*(param->mu_Q+param->A_t/param->tan_beta)/param->mass_gluino*( param->stop_mix[1][1]*param->stop_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_gluino/param->mass_gluino,param->mass_stl*param->mass_stl/param->mass_gluino/param->mass_gluino) +param->stop_mix[1][2]*param->stop_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_gluino/param->mass_gluino,param->mass_stl*param->mass_stl/param->mass_gluino/param->mass_gluino));			
	for(ie=1;ie<=nb_neut;ie++)	epsilon0p+=param->yub[3]*param->yub[3]/16./pi/pi*param->neut_mix[ie][4]*param->neut_mix[ie][3]*(param->mu_Q/param->tan_beta)/param->mass_neut[ie]*( param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][1]*param->stop_mix[1][1]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t1*param->mass_t1/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][1]*param->sbot_mix[1][1]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b2*param->mass_b2/param->mass_neut[ie]/param->mass_neut[ie]) +param->stop_mix[1][2]*param->stop_mix[1][2]*param->sbot_mix[1][2]*param->sbot_mix[1][2]*H2(param->mass_t2*param->mass_t2/param->mass_neut[ie]/param->mass_neut[ie],param->mass_b1*param->mass_b1/param->mass_neut[ie]/param->mass_neut[ie]));

	return epsilon0p;
}

/*----------------------------------------------------------------------*/

double epsilon_1p(struct parameters* param)
{
	if(param->SM==1) return 0.;

#ifdef SM_ChargedHiggs
	return 0.;
#endif

	return 1./16./pi/pi*(param->yub[3]*param->yub[3]*param->A_b/param->mu_Q*H2(pow(param->MqL3_Q/param->mu_Q,2.),pow(param->MbR_Q/param->mu_Q,2.))
	-param->g2*param->g2*param->M2_Q/param->mu_Q*H2(pow(param->MqL3_Q/param->mu_Q,2.),pow(param->M2_Q/param->mu_Q,2.)));
}

/*----------------------------------------------------------------------*/

void CW_calculator(double C0w[], double C1w[], double C2w[], double mu_W, struct parameters* param)
/* calculates the LO (C0w), NLO (C1w) and NNLO (C2w) contributions to the Wilson coefficients at scale mu_W, using the parameters of the structure param */
{
	int ie;
	for(ie=1;ie<=10;ie++) C0w[ie]=C1w[ie]=C2w[ie]=0.;
	
	double mass_top_muW=running_mass(param->mtmt,param->mtmt,mu_W,param->mass_top_pole,param->mass_b,param);
	double mass_b_muW=running_mass(param->mass_b,param->mass_b,mu_W,param->mass_top_pole,param->mass_b,param);

	double epsilonbp,epsilon0p,epsilon0,epsilon2,epsilon1p,epsilonb;
	if(param->THDM_model==0)
	{	
		epsilonbp=epsilon_bp(param);
		epsilon0p=epsilon_0p(param);
		epsilon0=epsilon_0(param);
		epsilon2=epsilon_2(param);
		epsilon1p=epsilon_1p(param);
		epsilonb=epsilon0+epsilon2;	
	}

/*----------------------------------------------------------------------*/
/* STANDARD MODEL */ 
/*----------------------------------------------------------------------*/

	double L=log(mu_W*mu_W/param->mass_W/param->mass_W);
 	double sw2=pow(sin(atan(param->gp/param->g2)),2.);

	//printf("sw2=%.5e\n",sw2);

	double xt= pow(mass_top_muW/param->mass_W,2.);
	double yt= pow(mass_top_muW/param->mass_H,2.);

	/* LO */
	
	double C2SM_0 = 1.;
	double C7SM_0 = -0.5*A0t(xt)-23./36.;
	double C8SM_0 = -0.5*F0t(xt)-1./3.;
	double C9SM_0 = (1.-4.*sw2)/sw2*C0t(xt)-B0t(xt)/sw2-D0t(xt) +1./4./sw2+38./27.-4./9.*L;
	double C10SM_0 = (B0t(xt)-C0t(xt))/sw2-1./4./sw2;
	
	//printf("C10SM=%.5e\n",C10SM_0);

	/* epsilon corrections */
	
	double C7SMeps_0,C8SMeps_0;	
	if((param->THDM_model==0)&&(param->SM==0)) 
	{
		C7SMeps_0= (epsilonb-epsilonbp)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F7_2(xt);
		C8SMeps_0= (epsilonb-epsilonbp)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F8_2(xt);
	}

	/* NLO */
	
	double C1SM_1 = 15.+6.*L;
	double C4SM_1 = E0t(xt)-7./9.+2./3.*L;
	double C7SM_1 = -0.5*A1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))+713./243.+4./81.*L-4./9.*C4SM_1;
	double C8SM_1 = -0.5*F1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))+91./324.-4./27.*L-C4SM_1/6.;
	double C9SM_1 = (1.-4.*sw2)/sw2*C1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))-B1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))/sw2-D1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW)) +1./sw2+524./729.-128./243.*pi*pi-16./3.*L-128./81.*L*L;
	double C10SM_1 = (B1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))-C1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW)))/sw2-1./sw2;

	/* NNLO */
	
	double C1SM_2 = -T(xt)+7987./72.+17.*pi*pi/3.+475./6.*L+17.*L*L;
	double C2SM_2 = 127./18.+4./3.*pi*pi+46./3.*L+4.*L*L;
	double C3SM_2 = G1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))-680./243.-20./81.*pi*pi-68./81.*L-20./27.*L*L;
	double C4SM_2 = E1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))+950./243.+10./81.*pi*pi+124./27.*L+10./27.*L*L;
	double C5SM_2 = -G1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))/10.+2./15.*E0t(xt)+68./243.+2./81.*pi*pi+14./81.*L+2./27.*L*L;
	double C6SM_2 = -3./16.*G1t(xt,log(mu_W*mu_W/mass_top_muW/mass_top_muW))+E0t(xt)/4.+85./162.+5./108.*pi*pi+35./108.*L+5./36.*L*L;

	double xtW=pow(running_mass(param->mtmt,param->mtmt,param->mass_W,param->mass_top_pole,param->mass_b,param)/param->mass_W,2.);
	double xtt=pow(param->mtmt/param->mass_W,2.);
	
	double C7SM_2 = (C7t2mt(xtt)+log(mu_W*mu_W/mass_top_muW/mass_top_muW)*((-592.*pow(xt,5.)-22.*pow(xt,4.)+12814.*pow(xt,3.)-6376.*xt*xt+512.*xt)/27./pow(xt-1.,5.)*Li2(1.-1./xt)
	+(-26838.*pow(xt,5.)+25938.*pow(xt,4.)+627367.*pow(xt,3.)-331956.*xt*xt+16989.*xt-460.)/729./pow(xt-1.,6.)*log(xt)
	+(34400.*pow(xt,5.)+276644.*pow(xt,4.)-2668324.*pow(xt,3.)+1694437.*xt*xt-323354.*xt+53077.)/2187./pow(xt-1.,5.)
	+log(mu_W*mu_W/mass_top_muW/mass_top_muW)*((-63.*pow(xt,5.)+532.*pow(xt,4.)+2089.*pow(xt,3.)-1118.*xt*xt)/9./pow(xt-1.,6.)*log(xt)
	+(1186.*pow(xt,5.)-2705.*pow(xt,4.)-24791.*pow(xt,3.)-16099.*xt*xt+19229.*xt-2740.)/162./pow(xt-1.,5.))) )
	-(C7c2MW(xtW)+13763./2187.*log(mu_W*mu_W/param->mass_W/param->mass_W)+814./729.*pow(log(mu_W*mu_W/param->mass_W/param->mass_W),2.));

	double C8SM_2 = (C8t2mt(xtt)+log(mu_W*mu_W/mass_top_muW/mass_top_muW)*((-148.*pow(xt,5.)+1052.*pow(xt,4.)-4811.*pow(xt,3.)-3520.*xt*xt-61.*xt)/18./pow(xt-1.,5.)*Li2(1.-1./xt)
	+(-15984.*pow(xt,5.)+152379.*pow(xt,4.)-1358060.*pow(xt,3.)-1201653.*xt*xt-74190.*xt+9188.)/1944./pow(xt-1.,6.)*log(xt)
	+(109669.*pow(xt,5.)-1112675.*pow(xt,4.)+6239377.*pow(xt,3.)+8967623.*xt*xt+768722.*xt-42796.)/11664./pow(xt-1.,5.)
	+log(mu_W*mu_W/mass_top_muW/mass_top_muW)*((-139.*pow(xt,4.)-2938.*pow(xt,3.)-2683.*xt*xt)/12./pow(xt-1.,6.)*log(xt)
	+(1295.*pow(xt,5.)-7009.*pow(xt,4.)+29495.*pow(xt,3.)+64513.*xt*xt+17458.*xt-2072.)/216./pow(xt-1.,5.))) )
	-(C8c2MW(xtW)+16607./5832.*log(mu_W*mu_W/param->mass_W/param->mass_W)+397./486.*pow(log(mu_W*mu_W/param->mass_W/param->mass_W),2.));
	
/*----------------------------------------------------------------------*/
/* CHARGED HIGGS */
/*----------------------------------------------------------------------*/
	
	double C7Heps_0,C8Heps_0,C7Heps2_0,C8Heps2_0;
	double lu,ld;

	/* epsilon corrections */
	if((param->THDM_model==0)&&(param->SM==0)) 
	{
		C7Heps_0=(-epsilon0p-epsilonb)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F7_2(yt);
		C8Heps_0=(-epsilon0p-epsilonb)/(1.+epsilonb*param->tan_beta)*param->tan_beta*F8_2(yt);

		C7Heps2_0=0.;
		C8Heps2_0=0.;

		if((param->mass_A02==0.)&&(param->mass_H03==0.))
		{
			C7Heps2_0=-epsilon2*epsilon1p*pow(param->tan_beta,2.)/(1.+epsilonb*param->tan_beta)/(1.+epsilon0*param->tan_beta)*F7_2(yt);
		C7Heps2_0+=epsilon2/pow(1.+epsilonb*param->tan_beta,2.)*(1.+pow(param->tan_beta,2.))/(1.+epsilon0*param->tan_beta)/72.		*((cos(param->alpha)+sin(param->alpha)*param->tan_beta)*(-sin(param->alpha)+epsilonb*cos(param->alpha))*pow(mass_b_muW/param->mass_h0,2.)
		+(sin(param->alpha)-cos(param->alpha)*param->tan_beta)*(cos(param->alpha)+epsilonb*sin(param->alpha))*pow(mass_b_muW/param->mass_H0,2.)			+(-cos(atan(param->tan_beta))-sin(atan(param->tan_beta))*param->tan_beta)*(sin(atan(param->tan_beta))-epsilonb*cos(atan(param->tan_beta)))*pow(mass_b_muW/param->mass_A0,2.));
	
		C8Heps2_0=-epsilon2*epsilon1p*pow(param->tan_beta,2.)/(1.+epsilonb*param->tan_beta)/(1.+epsilon0*param->tan_beta)*F8_2(yt);
		C8Heps2_0+=epsilon2/pow(1.+epsilonb*param->tan_beta,2.)*(1.+pow(param->tan_beta,2.))/(1.+epsilon0*param->tan_beta)/72.		*((cos(param->alpha)+sin(param->alpha)*param->tan_beta)*(-sin(param->alpha)+epsilonb*cos(param->alpha))*pow(mass_b_muW/param->mass_h0,2.)
		+(sin(param->alpha)-cos(param->alpha)*param->tan_beta)*(cos(param->alpha)+epsilonb*sin(param->alpha))*pow(mass_b_muW/param->mass_H0,2.)			+(-cos(atan(param->tan_beta))-sin(atan(param->tan_beta))*param->tan_beta)*(sin(atan(param->tan_beta))-epsilonb*cos(atan(param->tan_beta)))*pow(mass_b_muW/param->mass_A0,2.));
		}
		else
				{		C7Heps2_0=-epsilon2*epsilon1p*pow(param->tan_beta,2.)/(1.+epsilonb*param->tan_beta)/(1.+epsilon0*param->tan_beta)*F7_2(yt);
		C7Heps2_0+=epsilon2/pow(1.+epsilonb*param->tan_beta,2.)*(1.+pow(param->tan_beta,2.))/(1.+epsilon0*param->tan_beta)/72.	*((param->H0_mix[1][1]+param->H0_mix[1][2]*param->tan_beta)*(-param->H0_mix[1][2]+epsilonb*param->H0_mix[1][1])*pow(mass_b_muW/param->mass_h0,2.)
		+(param->H0_mix[2][1]+param->H0_mix[2][2]*param->tan_beta)*(-param->H0_mix[2][2]+epsilonb*param->H0_mix[2][1])*pow(mass_b_muW/param->mass_H0,2.)
		+(param->H0_mix[3][1]+param->H0_mix[3][2]*param->tan_beta)*(-param->H0_mix[3][2]+epsilonb*param->H0_mix[3][1])*pow(mass_b_muW/param->mass_H03,2.)

+(param->A0_mix[1][1]+param->A0_mix[1][2]*param->tan_beta)*(-param->A0_mix[1][2]+epsilonb*param->A0_mix[1][1])*pow(mass_b_muW/param->mass_A0,2.)
		+(param->A0_mix[2][1]+param->A0_mix[2][2]*param->tan_beta)*(-param->A0_mix[2][2]+epsilonb*param->A0_mix[2][1])*pow(mass_b_muW/param->mass_A02,2.));
		C8Heps2_0=-epsilon2*epsilon1p*pow(param->tan_beta,2.)/(1.+epsilonb*param->tan_beta)/(1.+epsilon0*param->tan_beta)*F8_2(yt);
		C8Heps2_0+=-3.*epsilon2/pow(1.+epsilonb*param->tan_beta,2.)*(1.+pow(param->tan_beta,2.))/(1.+epsilon0*param->tan_beta)/72.
		*((param->H0_mix[1][1]+param->H0_mix[1][2]*param->tan_beta)*(-param->H0_mix[1][2]+epsilonb*param->H0_mix[1][1])*pow(mass_b_muW/param->mass_h0,2.)
		+(param->H0_mix[2][1]+param->H0_mix[2][2]*param->tan_beta)*(-param->H0_mix[2][2]+epsilonb*param->H0_mix[2][1])*pow(mass_b_muW/param->mass_H0,2.)
		+(param->H0_mix[3][1]+param->H0_mix[3][2]*param->tan_beta)*(-param->H0_mix[3][2]+epsilonb*param->H0_mix[3][1])*pow(mass_b_muW/param->mass_H03,2.)		

+(param->A0_mix[1][1]+param->A0_mix[1][2]*param->tan_beta)*(-param->A0_mix[1][2]+epsilonb*param->A0_mix[1][1])*pow(mass_b_muW/param->mass_A0,2.)
		+(param->A0_mix[2][1]+param->A0_mix[2][2]*param->tan_beta)*(-param->A0_mix[2][2]+epsilonb*param->A0_mix[2][1])*pow(mass_b_muW/param->mass_A02,2.));
		}

		lu=1./param->tan_beta;
		ld=-param->tan_beta;
	}
	else
	{
		lu=param->lambda_u[3][3];
		ld=param->lambda_d[3][3];
	}

	/* LO */
	
	double C7H_0=1./3.*lu*lu*F7_1(yt) - lu*ld*F7_2(yt);
	double C8H_0=1./3.*lu*lu*F8_1(yt) - lu*ld*F8_2(yt);

	double C9H_0=(1.-4.*sw2)/sw2*C9llH0(xt,yt,lu)-D9H0(yt,lu);
	double C10H_0=-C9llH0(xt,yt,lu)/sw2;

	//printf("C10H=%.5e\n",C10H_0);

	/* NLO */

	double C4H_1=EH(yt,lu);
	
 	double C7H_1= G7H(yt,lu,ld)+Delta7H(yt,lu,ld)*log(pow(mu_W/param->mass_H,2.))-4./9.*C4H_1;
	double C8H_1= G8H(yt,lu,ld)+Delta8H(yt,lu,ld)*log(pow(mu_W/param->mass_H,2.))-1./6.*C4H_1;
	double C9H_1=(1.-4.*sw2)/sw2*C9llH1(xt,yt,lu,log(pow(mu_W/param->mass_H,2.)))-D9H1(yt,lu,log(pow(mu_W/param->mass_H,2.)));
	double C10H_1=-C9llH1(xt,yt,lu,log(pow(mu_W/param->mass_H,2.)))/sw2;


	/* NNLO */	

	double C3H_2=G3H(yt,lu)+Delta3H(yt,lu)*log(pow(mu_W/param->mass_H,2.));
	double C4H_2=G4H(yt,lu)+Delta4H(yt,lu)*log(pow(mu_W/param->mass_H,2.));
	double C5H_2=-C3H_2/10.+2./15.*C4H_1;
	double C6H_2=-3./16.*C3H_2+1./4.*C4H_1;
	

/*----------------------------------------------------------------------*/
/* CHARGINOS */	
/*----------------------------------------------------------------------*/	
	double C4charg_1,C4charg_2;
	double C3charg_2,C5charg_2,C6charg_2;
	double C7charg_0,C8charg_0,C7_chargeps_0,C8_chargeps_0,C7charg_1,C8charg_1;
	double C9charg_0,C9charg_1,C10charg_0,C10charg_1;
	double C7four_1,C8four_1,C9four_1,C10four_1,C4four_2;
	double C1squark_2;
	
	double Gamma_UL[7][4],Gamma_UR[7][4],Gamma_NL[4][4],Gamma_NR[4][4];
	double Gamma_U[7][7],I_LR[7][7],P_U[7][7];
	double X_UL[3][7][4],X_UR[3][7][4],X_NL[3][4][4],X_NR[3][4][4];
	double MU[4],MD[4],ME[4],VCKM[4][4],Mch[3],MsqU[7],MsqD[7],Msn[4],mintmp;
	double kappa,ag,aY,cosb,sinb,st,ct,alphas_mg;
	int ae,be,ce,de,ee,fe,ge,je,ke;
	
	if((param->THDM_model==0)&&(param->SM==0)) 
	{
		alphas_mg=alphas_running(param->mass_gluino,param->mass_top_pole,param->mass_b_pole,param);
		ag=1.-7./12./pi*alphas_mg;
		aY=1.+alphas_mg/4./pi;
		
		kappa=1./(param->g2*param->g2*param->Vtb*param->Vts);
		
		VCKM[1][1]=param->Vud;
		VCKM[1][2]=param->Vus;
		VCKM[1][3]=-(param->Vts*param->Vtb+param->Vcs*param->Vcb)/param->Vus; /* Vub from unitarity */
		VCKM[2][1]=param->Vcd;
		VCKM[2][2]=param->Vcs;
		VCKM[2][3]=param->Vcb;
		VCKM[3][1]=param->Vtd;
		VCKM[3][2]=param->Vts;
		VCKM[3][3]=param->Vtb;
		
		sinb=sin(atan(param->tan_beta));
		cosb=cos(atan(param->tan_beta));
		ct=param->stop_mix[2][2];
		st=param->stop_mix[1][2];
		
		MU[1]=param->mass_u;
		MU[2]=param->mass_c;
		MU[3]=mass_top_muW;

		MD[1]=param->mass_u;
		MD[2]=param->mass_s;
		MD[3]=mass_b_muW;

		ME[1]=param->mass_e;
		ME[2]=param->mass_mu;
		ME[3]=param->mass_tau;

		Mch[1]=param->mass_cha1;
		Mch[2]=param->mass_cha2;
		
		MsqU[1]=param->mass_upl;
		MsqU[2]=param->mass_chl;
		MsqU[3]=param->mass_t1;
		MsqU[4]=param->mass_upr;
		MsqU[5]=param->mass_chr;
		MsqU[6]=param->mass_t2;
		
		Msn[1]=param->mass_nuel;
		Msn[2]=param->mass_numl;
		Msn[3]=param->mass_nutl;		
		if(param->sU_mix[1][1]*param->sU_mix[1][2]*param->sU_mix[1][3]*param->sU_mix[1][4]*param->sU_mix[1][5]*param->sU_mix[1][6]!=0.)
		{
		
			for(je=1;je<=5;je++)
			{
				for(ie=je+1;ie<=6;ie++) if (MsqU[je]>MsqU[ie])
				{
					mintmp=MsqU[je];
					MsqU[je]=MsqU[ie];
					MsqU[ie]=mintmp;
				}
			}
		
		
			for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UL[ae][ie]=param->sU_mix[ae][ie];
			for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UR[ae][ie]=param->sU_mix[ae][ie+3];
		}
		else
		{
			Gamma_UL[1][1]=1.;
			Gamma_UL[2][1]=0.;
			Gamma_UL[3][1]=0.;
			Gamma_UL[4][1]=0.;
			Gamma_UL[5][1]=0.;
			Gamma_UL[6][1]=0.;
			Gamma_UL[1][2]=0.;
			Gamma_UL[2][2]=1.;
			Gamma_UL[3][2]=0.;
			Gamma_UL[4][2]=0.;
			Gamma_UL[5][2]=0.;
			Gamma_UL[6][2]=0.;
			Gamma_UL[1][3]=0.;
			Gamma_UL[2][3]=0.;
			Gamma_UL[3][3]=ct;
			Gamma_UL[4][3]=0.;
			Gamma_UL[5][3]=0.;
			Gamma_UL[6][3]=-st;
		
			Gamma_UR[1][1]=0.;
			Gamma_UR[2][1]=0.;
			Gamma_UR[3][1]=0.;
			Gamma_UR[4][1]=1.;
			Gamma_UR[5][1]=0.;
			Gamma_UR[6][1]=0.;
			Gamma_UR[1][2]=0.;
			Gamma_UR[2][2]=0.;
			Gamma_UR[3][2]=0.;
			Gamma_UR[4][2]=0.;
			Gamma_UR[5][2]=1.;
			Gamma_UR[6][2]=0.;
			Gamma_UR[1][3]=0.;
			Gamma_UR[2][3]=0.;
			Gamma_UR[3][3]=st;
			Gamma_UR[4][3]=0.;
			Gamma_UR[5][3]=0.;
			Gamma_UR[6][3]=ct;
		}


		for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++)
		{
			Gamma_U[ae][ie]=Gamma_UL[ae][ie];
			Gamma_U[ae][ie+3]=Gamma_UR[ae][ie];
		}
		
		for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) I_LR[ae][be]=0.;
		for(ae=1;ae<=3;ae++) I_LR[ae][ae]=1.;
		for(ae=4;ae<=6;ae++) I_LR[ae][ae]=-1.;
		
		for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++) for(de=1;de<=6;de++) P_U[ae][be]=Gamma_U[ae][ce]*I_LR[ce][de]*Gamma_U[be][de];
		
		for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++) if(ae==be) Gamma_NL[ae][be]=Gamma_NR[ae][be]=1.; else Gamma_NL[ae][be]=Gamma_NR[ae][be]=0.;
				
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
		{
			X_UL[ie][ae][be]=0.;
			for(ce=1;ce<=3;ce++) X_UL[ie][ae][be]+=-param->g2*(ag*param->charg_Vmix[ie][1]*Gamma_UL[ae][ce]
			-aY*param->charg_Vmix[ie][2]*Gamma_UR[ae][ce]*MU[ce]/(sqrt(2.)*param->mass_W*sinb))*VCKM[ce][be];
		
			X_UR[ie][ae][be]=0.;
			for(ce=1;ce<=3;ce++) X_UR[ie][ae][be]+=param->g2*aY*param->charg_Umix[ie][2]*Gamma_UL[ae][ce]*VCKM[ce][be]*MD[be]/(sqrt(2)*param->mass_W*cosb);
		}
	
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++)
		{
			X_NL[ie][ae][be]=-param->g2*param->charg_Vmix[ie][1]*Gamma_NL[ae][be];
			X_NR[ie][ae][be]=param->g2*param->charg_Umix[ie][2]*Gamma_NL[ae][be]*ME[be]/(sqrt(2.)*param->mass_W*cosb);
		}

		/* LO */

		C7charg_0=0.;
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C7charg_0+=pow(param->mass_W/Mch[ie],2.)*(X_UL[ie][ae][2]*X_UL[ie][ae][3]*h10(pow(MsqU[ae]/Mch[ie],2.)) + Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ae][3]*h20(pow(MsqU[ae]/Mch[ie],2.)));		
		C7charg_0*=-0.5*kappa; 
		

		C8charg_0=0.;
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C8charg_0+=pow(param->mass_W/Mch[ie],2.)*(X_UL[ie][ae][2]*X_UL[ie][ae][3]*h50(pow(MsqU[ae]/Mch[ie],2.)) + Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ae][3]*h60(pow(MsqU[ae]/Mch[ie],2.)));		
		C8charg_0*=-0.5*kappa; 
		

		C7_chargeps_0=0.;		
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C7_chargeps_0+= -epsilonb/(1.+epsilonb*param->tan_beta)*param->tan_beta*pow(param->mass_W/Mch[ie],2.)*(Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ae][3]*h20(pow(MsqU[ae]/Mch[ie],2.)));				
		C7_chargeps_0*=-0.5*kappa; 
		

		C8_chargeps_0=0.;		
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C8_chargeps_0+= -epsilonb/(1.+epsilonb*param->tan_beta)*param->tan_beta*pow(param->mass_W/Mch[ie],2.)*(Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ae][3]*h60(pow(MsqU[ae]/Mch[ie],2.)));		
		C8_chargeps_0*=-0.5*kappa; 

		double B0c1=0.;
		double B0c2=0.;
		for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
		{ 	B0c1+=X_UL[je][ae][2]*X_UL[ie][ae][3]/Mch[ie]/Mch[ie]*(0.5*X_NL[ie][be][2]*X_NL[je][be][2]*f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));	 B0c2+=X_UL[je][ae][2]*X_UL[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NR[ie][be][2]*X_NR[je][be][2]*fabs(Mch[je]/Mch[ie])*f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
		}
		
		double B90c=-(B0c1-B0c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2;
		
		double B100c=(B0c1+B0c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2;
		
		
		double C90c=0.;
		for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) C90c+=X_UL[je][ae][2]*X_UL[ie][ae][3]*(2.*fabs(Mch[je]/Mch[ie])*f30(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))*param->charg_Umix[je][1]*param->charg_Umix[ie][1] -f40(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))*param->charg_Vmix[je][1]*param->charg_Vmix[ie][1]);
		
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=3;ce++) C90c+=X_UL[ie][be][2]*X_UL[ie][ae][3]*f40(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))*Gamma_UL[be][ce]*Gamma_UL[ae][ce];
		C90c*=-kappa/8.;	
		
		
		double D90c=0.;
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) D90c+=pow(param->mass_W/Mch[ie],2.)*X_UL[ie][ae][2]*X_UL[ie][ae][3]*h30(pow(MsqU[ae]/Mch[ie],2.));
		D90c*=kappa;

		C9charg_0=(1.-4.*sw2)/sw2*C90c-B90c/sw2-D90c;
		C10charg_0=(B100c-C90c)/sw2;

		//printf("C10charg=%.5e\n",C10charg_0);


		/* NLO */
		
		C4charg_1=0.;		
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C4charg_1+= pow(param->mass_W/Mch[ie],2.)*(X_UL[ie][ae][2]*X_UL[ie][ae][3]*h40(pow(MsqU[ae]/Mch[ie],2.)));	
		C4charg_1*=kappa; 


		C7charg_1=0.;
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C7charg_1+=pow(param->mass_W/Mch[ie],2.)*(X_UL[ie][ae][2]*X_UL[ie][ae][3]*h11(pow(MsqU[ae]/Mch[ie],2.),log(pow(mu_W/MsqU[ae],2.))) + Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ae][3]*h21(pow(MsqU[ae]/Mch[ie],2.),log(pow(mu_W/MsqU[ae],2.))));		
		C7charg_1*=-0.5*kappa; 
		

		C8charg_1=0.;
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C8charg_1+=pow(param->mass_W/Mch[ie],2.)*(X_UL[ie][ae][2]*X_UL[ie][ae][3]*h51(pow(MsqU[ae]/Mch[ie],2.),log(pow(mu_W/MsqU[ae],2.))) + Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ae][3]*h61(pow(MsqU[ae]/Mch[ie],2.),log(pow(mu_W/MsqU[ae],2.))));	
		C8charg_1*=-0.5*kappa; 


		double B1c1=0.;
		double B1c2=0.;
		
		for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
		{
		B1c1+=X_UL[je][ae][2]*X_UL[ie][ae][3]/Mch[ie]/Mch[ie]*(0.5*X_NL[ie][be][2]*X_NL[je][be][2]*(
		f81(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.))
		+4.*(f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.))
			+(f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*1.0001,pow(Msn[be]/Mch[ie],2.))-f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*0.9999,pow(Msn[be]/Mch[ie],2.)))/0.0002
		)*log(pow(mu_W/MsqU[ae],2.))));
		B1c2+=X_UL[je][ae][2]*X_UL[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NR[ie][be][2]*X_NR[je][be][2]*fabs(Mch[je]/Mch[ie])*(
		f91(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.))
		+4.*(f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.))						+(f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*1.0001,pow(Msn[be]/Mch[ie],2.))-f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*0.9999,pow(Msn[be]/Mch[ie],2.)))/0.0002
		)*log(pow(mu_W/MsqU[ae],2.))));
		}
		
		double B91c=-(B1c1-B1c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2;
		
		double B101c=(B1c1+B1c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2;		
	
	
		double C91c=0.;
		for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) C91c+=X_UL[je][ae][2]*X_UL[ie][ae][3]*(
		(2.*fabs(Mch[je]/Mch[ie])*(f31(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))
			+4.*(f30(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))+(f30(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*1.0001)-f30(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*0.9999))/0.0002)*log(pow(mu_W/MsqU[ae],2.)))*param->charg_Umix[je][1]*param->charg_Umix[ie][1])
		
		 -(f41(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))
			+4.*(f40(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))+(f40(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*1.0001)-f40(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*0.9999))/0.0002)*log(pow(mu_W/MsqU[ae],2.)))*param->charg_Vmix[je][1]*param->charg_Vmix[ie][1]);
		
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=3;ce++) C91c+=X_UL[ie][be][2]*X_UL[ie][ae][3]*(f51(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))
			+4.*(f40(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))+(f40(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.)*1.0001)-f40(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.)*0.9999))/0.0002+(f40(pow(MsqU[ae]/Mch[ie],2.)*1.0001,pow(MsqU[be]/Mch[ie],2.))-f40(pow(MsqU[ae]/Mch[ie],2.)*0.9999,pow(MsqU[be]/Mch[ie],2.)))/0.0002)*log(pow(mu_W/MsqU[ae],2.)))*Gamma_UL[be][ce]*Gamma_UL[ae][ce];
		C91c*=-kappa/8.;	
		
		
		double D91c=0.;
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) D91c+=pow(param->mass_W/Mch[ie],2.)*X_UL[ie][ae][2]*X_UL[ie][ae][3]*h31(pow(MsqU[ae]/Mch[ie],2.),log(pow(mu_W/MsqU[ae],2.)));
		D91c*=kappa;

		C9charg_1=(1.-4.*sw2)/sw2*C91c-B91c/sw2-D91c;
		C10charg_1=(B101c-C91c)/sw2;

		C7four_1=0.;
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++) C7four_1+=pow(param->mass_W/Mch[ie],2.)*P_U[ae][be]*MsqU[be]/Mch[ie]*P_U[be][ce]*(1.+log(pow(mu_W/MsqU[be],2.)))*
		(
		X_UL[ie][ae][2]*X_UL[ie][ce][3]*(-q11(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.))+2./3.*q21(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.)))
		+Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ce][3]*(-q31(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.))+2./3.*q41(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.)))
		);
		C7four_1*=-0.5*kappa;


		C8four_1=0.;
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++) C8four_1+=pow(param->mass_W/Mch[ie],2.)*P_U[ae][be]*MsqU[be]/Mch[ie]*P_U[be][ce]*(1.+log(pow(mu_W/MsqU[be],2.)))*
		(
		X_UL[ie][ae][2]*X_UL[ie][ce][3]*q21(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.))
		+Mch[ie]/mass_b_muW*X_UL[ie][ae][2]*X_UR[ie][ce][3]*q41(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.))
		);
		C8four_1*=-0.5*kappa;


		double B1f1=0.;
		double B1f2=0.;
		for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(fe=1;fe<=3;fe++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++)
		{ B1f1+=pow(param->mass_W/Mch[ie],2.)*P_U[ae][be]*pow(MsqU[be]/Mch[ie],2.)*P_U[be][ce]*(1.+log(pow(mu_W/MsqU[be],2.)))
		*X_UL[je][ae][2]*X_UL[ie][ce][3]*(
		0.5*f90(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.),pow(Msn[fe]/Mch[ie],2.))*X_NL[ie][fe][2]*X_NL[je][fe][2]);
		
		B1f2+=pow(param->mass_W/Mch[ie],2.)*P_U[ae][be]*pow(MsqU[be]/Mch[ie],2.)*P_U[be][ce]*(1.+log(pow(mu_W/MsqU[be],2.)))
		*X_UL[je][ae][2]*X_UL[ie][ce][3]*(fabs(Mch[je]/Mch[ie])*f100(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.),pow(Msn[fe]/Mch[ie],2.))*X_NR[ie][fe][2]*X_NR[je][fe][2]
		);
		}
		double B91f=(B1f1-B1f2)*2./3.*kappa/param->g2/param->g2;
				
		double B101f=-(B1f1+B1f2)*2./3.*kappa/param->g2/param->g2;
		
	
		double C91f=0.;
		for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(de=1;de<=6;de++) for(ke=1;ke<=6;ke++) C91f+=P_U[de][ke]*pow(MsqU[ke]/Mch[ie],2.)*P_U[ke][ae]*(1.+log(pow(mu_W/MsqU[ke],2.)))*X_UL[je][de][2]*X_UL[ie][ae][3]*(
2.*fabs(Mch[je]/Mch[ie])*f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[de]/Mch[ie],2.))*param->charg_Umix[je][1]*param->charg_Umix[ie][1]	-f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[de]/Mch[ie],2.))*param->charg_Vmix[je][1]*param->charg_Vmix[ie][1]);
		
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++)  for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++) for(fe=1;fe<=3;fe++) for(ke=1;ke<=6;ke++) C91f+=P_U[be][ke]*pow(MsqU[ke]/Mch[ie],2.)*P_U[ke][ae]*(1.+log(pow(mu_W/MsqU[ke],2.)))*X_UL[ie][ce][2]*X_UL[ie][ae][3]*	f50(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.))*Gamma_UL[ce][fe]*Gamma_UL[be][fe];
		
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(ce=1;ce<=6;ce++) for(de=1;de<=6;de++) for(fe=1;fe<=3;fe++) for(ke=1;ke<=6;ke++) C91f+=P_U[de][ke]*pow(MsqU[ke]/Mch[ie],2.)*P_U[ke][ce]*(1.+log(pow(mu_W/MsqU[ke],2.)))*X_UL[ie][de][2]*X_UL[ie][ae][3]*
f50(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.),pow(MsqU[de]/Mch[ie],2.))*Gamma_UL[ce][fe]*Gamma_UL[ae][fe];
		
		C91f*=kappa/6.;			
		
		
		double D91f=0.;		
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++) D91f+= pow(param->mass_W/Mch[ie],2.)*P_U[ae][be]*pow(MsqU[be]/Mch[ie],2.)*P_U[be][ce]*(1.+log(pow(mu_W/MsqU[be],2.)))*X_UL[ie][ae][2]*X_UL[ie][ce][3]*q51(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.));
		D91f*=kappa; 

		C9four_1=(1.-4.*sw2)/sw2*C91f-B91f/sw2-D91f;
		C10four_1=(B101f-C91f)/sw2;


		/* NNLO */

		C3charg_2=0.;
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C3charg_2+= pow(param->mass_W/Mch[ie],2.)*X_UL[ie][ae][2]*X_UL[ie][ae][3]*h71(pow(MsqU[ae]/Mch[ie],2.),log(pow(mu_W/MsqU[ae],2.)));	
		C3charg_2*=kappa; 
		
		
		C4charg_2=0.;		
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C4charg_2+= pow(param->mass_W/Mch[ie],2.)*X_UL[ie][ae][2]*X_UL[ie][ae][3]*h41(pow(MsqU[ae]/Mch[ie],2.),log(pow(mu_W/MsqU[ae],2.)));	
		C4charg_2*=kappa; 

		C5charg_2=-C3charg_2/10.+2./15.*C4charg_1;
		C6charg_2=-3./16.*C3charg_2+1./4.*C4charg_1;


		C4four_2=0.;		
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++) C4four_2+= pow(param->mass_W/Mch[ie],2.)*P_U[ae][be]*MsqU[be]/Mch[ie]*P_U[be][ce]*(1.+log(pow(mu_W/MsqU[be],2.)))*X_UL[ie][ae][2]*X_UL[ie][ce][3]*q61(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.));
		C4four_2*=kappa;

	
		
/*----------------------------------------------------------------------*/
/* SQUARKS */	
/*----------------------------------------------------------------------*/	
		
		MsqD[1]=param->mass_dnl;
		MsqD[2]=param->mass_stl;
		MsqD[3]=param->mass_b1;
		MsqD[4]=param->mass_dnr;
		MsqD[5]=param->mass_str;
		MsqD[6]=param->mass_b2;
		
		int test=1;
		for(ae=1;ae<=6;ae++) test=test&&(fabs(MsqU[ae])>param->mass_W/2.)&&(fabs(MsqD[ae])>param->mass_W/2.);
		
		if(test)
		{		
			C1squark_2=-208./3.;
			double xsqa;
			for(ae=1;ae<=6;ae++) 
			{
				xsqa=pow(MsqU[ae]/param->mass_W,2.);
				C1squark_2+=-2.*pow(4.*xsqa-1.,1.5)*Cl2(2.*asin(0.5/sqrt(xsqa))) +8.*(xsqa-1./3.)*log(xsqa)+16.*xsqa;
			
				xsqa=pow(MsqD[ae]/param->mass_W,2.);
				C1squark_2+=-2.*pow(4.*xsqa-1.,1.5)*Cl2(2.*asin(0.5/sqrt(xsqa))) +8.*(xsqa-1./3.)*log(xsqa)+16.*xsqa;
			}
		}
		else C1squark_2=0.;
	}

/*----------------------------------------------------------------------*/
/* C TOTAL */	
/*----------------------------------------------------------------------*/

	if(param->SM==1) 
	{
		C7SMeps_0=C7H_0=C7Heps_0=C7Heps2_0=C7charg_0=C7_chargeps_0=C7H_1=C7charg_1=C7four_1=0.;
		C8SMeps_0=C8H_0=C8Heps_0=C8Heps2_0=C8charg_0=C8_chargeps_0=C8H_1=C8charg_1=C8four_1=0.;
		C3H_2=C3charg_2=0.;
		C4H_1=C4charg_1=C4H_2=C4charg_2=C4four_2=0.;
		C5H_2=C5charg_2=0.;
		C6H_2=C6charg_2=0.;
		C9H_0=C9H_1=C9charg_0=C9charg_1=C9four_1=0.;
		C10H_0=C10H_1=C10charg_0=C10charg_1=C10four_1=0.;
		C1squark_2=0.;
	}

#ifdef SM_ChargedHiggs
	C7SMeps_0=C7Heps_0=C7Heps2_0=C7charg_0=C7_chargeps_0=C7charg_1=C7four_1=0.;
	C8SMeps_0=C8Heps_0=C8Heps2_0=C8charg_0=C8_chargeps_0=C8charg_1=C8four_1=0.;
	C3H_2=C3charg_2=0.;
	C4H_1=C4charg_1=C4H_2=C4charg_2=C4four_2=0.;
	C5H_2=C5charg_2=0.;
	C6H_2=C6charg_2=0.;
	C9H_0=C9H_1=C9charg_0=C9charg_1=C9four_1=0.;
	C10H_0=C10H_1=C10charg_0=C10charg_1=C10four_1=0.;
	C1squark_2=0.;
#endif

	double alphas_mu=alphas_running(mu_W,param->mass_top_pole,param->mass_b_pole,param);

	if(cabs(C7H_1)*alphas_mu/4./pi>cabs(C0w[7])) C7H_1*=cabs(C0w[7])/cabs(C7H_1)*4.*pi/alphas_mu;
	if(cabs(C7charg_1)*alphas_mu/4./pi>cabs(C0w[7])) C7charg_1*=cabs(C0w[7])/cabs(C7charg_1)*4.*pi/alphas_mu;
	if(cabs(C7four_1)*alphas_mu/4./pi>cabs(C0w[7])) C7four_1*=cabs(C0w[7])/cabs(C7four_1)*4.*pi/alphas_mu;

	if(cabs(C8H_1)*alphas_mu/4./pi>cabs(C0w[8])) C8H_1*=cabs(C0w[8])/cabs(C8H_1)*4.*pi/alphas_mu;
	if(cabs(C8charg_1)*alphas_mu/4./pi>cabs(C0w[8])) C8charg_1*=cabs(C0w[8])/cabs(C8charg_1)*4.*pi/alphas_mu;
	if(cabs(C8four_1)*alphas_mu/4./pi>cabs(C0w[8])) C8four_1*=cabs(C0w[8])/cabs(C8four_1)*4.*pi/alphas_mu;   

	if(cabs(C9H_1)*alphas_mu/4./pi>cabs(C0w[9])) C9H_1*=cabs(C0w[9])/cabs(C9H_1)*4.*pi/alphas_mu;
	if(cabs(C9charg_1)*alphas_mu/4./pi>cabs(C0w[9])) C9charg_1*=cabs(C0w[9])/cabs(C9charg_1)*4.*pi/alphas_mu;
	if(cabs(C9four_1)*alphas_mu/4./pi>cabs(C0w[9])) C9four_1*=cabs(C0w[9])/cabs(C9four_1)*4.*pi/alphas_mu;

	if(cabs(C10H_1)*alphas_mu/4./pi>cabs(C0w[10])) C10H_1*=cabs(C0w[10])/cabs(C10H_1)*4.*pi/alphas_mu;
	if(cabs(C10charg_1)*alphas_mu/4./pi>cabs(C0w[10])) C10charg_1*=cabs(C0w[10])/cabs(C10charg_1)*4.*pi/alphas_mu;
	if(cabs(C10four_1)*alphas_mu/4./pi>cabs(C0w[10])) C10four_1*=cabs(C0w[10])/cabs(C10four_1)*4.*pi/alphas_mu;

	if(param->THDM_model==0)
	{
		C0w[2]= C2SM_0;
		C0w[7]= C7SM_0+C7SMeps_0+C7H_0+C7Heps_0+C7Heps2_0+C7charg_0+C7_chargeps_0;
		C0w[8]= C8SM_0+C8SMeps_0+C8H_0+C8Heps_0+C8Heps2_0+C8charg_0+C8_chargeps_0;
		C0w[9]= C9SM_0+C9H_0+C9charg_0;
		C0w[10]= C10SM_0+C10H_0+C10charg_0;
		C1w[1]= C1SM_1;
		C1w[4]= C4SM_1+C4H_1+C4charg_1;
		C1w[7]= C7SM_1+C7H_1+C7charg_1+C7four_1;
		C1w[8]= C8SM_1+C8H_1+C8charg_1+C8four_1;
		C1w[9]= C9SM_1+C9H_1+C9charg_1+C9four_1;
		C1w[10]= C10SM_1+C10H_1+C10charg_1+C10four_1;
		C2w[1]= C1SM_2+C1squark_2;
		C2w[2]= C2SM_2;
		C2w[3]= C3SM_2+C3H_2+C3charg_2;
		C2w[4]= C4SM_2+C4H_2+C4charg_2+C4four_2;
		C2w[5]= C5SM_2+C5H_2+C5charg_2;
		C2w[6]= C6SM_2+C6H_2+C6charg_2;
		C2w[7]= C7SM_2;
		C2w[8]= C8SM_2;
	}
	else
	{
		C0w[2]= C2SM_0;
		C0w[7]= C7SM_0+C7H_0;
		C0w[8]= C8SM_0+C8H_0;
		C0w[9]= C9SM_0+C9H_0;
		C0w[10]= C10SM_0+C10H_0;
		C1w[1]= C1SM_1;
		C1w[4]= C4SM_1+C4H_1;	
		C1w[7]= C7SM_1+C7H_1;
		C1w[8]= C8SM_1+C8H_1;
		C1w[9]= C9SM_1+C9H_1;
		C1w[10]= C10SM_1+C10H_1;
		C2w[1]= C1SM_2;
		C2w[2]= C2SM_2;
		C2w[3]= C3SM_2+C3H_2;
		C2w[4]= C4SM_2+C4H_2;
		C2w[5]= C5SM_2+C5H_2;
		C2w[6]= C6SM_2+C6H_2;
		C2w[7]= C7SM_2;
		C2w[8]= C8SM_2;
	}

	return;
}

/*-----------------------------------------------------------------------------*/

void C_calculator_base1(double C0w[], double C1w[], double C2w[], double mu_W, double C0b[], double C1b[], double C2b[], double mu, struct parameters* param)
/* calculates the LO (C0b), NLO (C1b) and NNLO (C2b) contributions to the Wilson coefficients at scale mu, using the LO (C0w), NLO (C1w) and NNLO (C2w) contributions to the Wilson coefficients at scale mu_W and the parameters of the structure param, in the standard operator basis */
{
	int ie;

	for(ie=1;ie<=10;ie++) C0b[ie]=C1b[ie]=C2b[ie]=0.;

	double alphas_muW=alphas_running(mu_W,param->mass_top_pole,param->mass_b_pole,param);
	double alphas_mu=alphas_running(mu,param->mass_top_pole,param->mass_b_pole,param);	
	double eta_mu=alphas_muW/alphas_mu;

	double C0w7= C0w[7]-1./3.*C0w[3]-4./9.*C0w[4]-20./3.*C0w[5]-80./9.*C0w[6]; 
	double C1w7= C1w[7]-1./3.*C1w[3]-4./9.*C1w[4]-20./3.*C1w[5]-80./9.*C1w[6]; 
	double C2w7= C2w[7]-1./3.*C2w[3]-4./9.*C2w[4]-20./3.*C2w[5]-80./9.*C2w[6]; 

	double C0w8= C0w[8]+C0w[3]-1./6.*C0w[4]+20.*C0w[5]-10./3.*C0w[6]; 
	double C1w8= C1w[8]+C1w[3]-1./6.*C1w[4]+20.*C1w[5]-10./3.*C1w[6]; 
	double C2w8= C2w[8]+C2w[3]-1./6.*C2w[4]+20.*C2w[5]-10./3.*C2w[6]; 

	double m00[10][10][10]={ 
	/*m001[10][10]*/{
		{0 , 0 , 0.333333 , 0.666667 , 0 , 0 , 0 , 0 , 0 , 0}, 
		{0 , 0 , 1. , -1. , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0},
	},
	    
	/*m002[10][10]*/{
		{0 , 0 , 0.222222 , -0.222222 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0.666667 , 0.333333 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0}
	},
	    
	/*m003[10][10]*/{
		{0 , 0 , 0.010582 , 0.0246914 , -0.0129383 , -0.049685 , 0.00917753 , 0.0181725 , 0 , 0 },
		{0 , 0 , 0.031746 , -0.037037 , -0.065915 , 0.0594717 , -0.021794 , 0.0335283 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.28966 , -0.20057 , 0.175618 , 1.31461 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.193297 , 0.157873 , 0.142818 , -0.107394 , 0 , 0 },
		{0 , 0 , 0 , 0 , -5.6115 , -4.14773 , 0.844839 , 8.91439 , 0 , 0 },
		{0 , 0 , 0 , 0 , -2.58268 , 2.01866 , 0.0268102 , 0.537203 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0}
	},
	/*m004[10][10]*/{
		{0 , 0 , 0.015873 , -0.0740741 , 0.0046498 , 0.0144314 , 0.056241 , -0.0171212 , 0 , 0 },
		{0 , 0 , 0.047619 , 0.111111 , 0.0236887 , -0.0172741 , -0.133556 , -0.0315887 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.104099 , 0.0582573 , 1.07621 , -1.23856 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.0694676 , -0.0458555 , 0.875207 , 0.101181 , 0 , 0 },
		{0 , 0 , 0 , 0 , 2.01667 , 1.20474 , 5.17728 , -8.39869 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.928167 , -0.586338 , 0.164296 , -0.506125 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0}
	},
	/*m005[10][10]*/{
		{0 , 0 , -0.0026455 , -0.00617284 , 0.00183608 , 0.00831955 , -0.000414444 , -0.000922849 , 0 , 0 },
		{0 , 0 , -0.00793651 , 0.00925926 , 0.00935402 , -0.00995829 , 0.000984185 , -0.00170266 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.0411058 , 0.0335846 , -0.00793065 , -0.0667598 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.0274309 , -0.0264352 , -0.00644947 , 0.00545376 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.79633 , 0.69452 , -0.0381518 , -0.452698 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.366508 , -0.338017 , -0.00121071 , -0.0272807 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0}
	},
	/*m006[10][10]*/{
		{0 , 0 , -0.00396825 , 0.0185185 , 0.00212341 , -0.013584 , -0.00433812 , 0.0012484 , 0 , 0 },
		{0 , 0 , -0.0119048 , -0.0277778 , 0.0108178 , 0.0162596 , 0.0103018 , 0.0023033 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.0475384 , -0.0548362 , -0.0830125 , 0.0903103 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.0317235 , 0.0431627 , -0.0675086 , -0.00737765 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.920947 , -1.134 , -0.399346 , 0.612394 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.423863 , 0.551905 , -0.0126729 , 0.0369043 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0}
	},
	/*m007[10][10]*/{
		{0.578391 , -0.392139 , -0.142857 , 0.047619 , -0.127461 , 0.0317363 , 0.00781139 , -0.00310131 , 0 , 0 },
		{2.29959 , -1.08798 , -0.428571 , -0.0714286 , -0.649357 , -0.0379876 , -0.0185498 , -0.00572193 , 0 , 0 },
		{8.078 , -5.27767 , 0 , 0 , -2.85357 , 0.128114 , 0.149476 , -0.224351 , 0 , 0 },
		{5.70644 , -3.84123 , 0 , 0 , -1.90425 , -0.100841 , 0.121559 , 0.0183278 , 0 , 0 },
		{202.901 , -149.467 , 0 , 0 , -55.2813 , 2.64936 , 0.719079 , -1.52133 , 0 , 0 },
		{86.4618 , -59.6604 , 0 , 0 , -25.443 , -1.28942 , 0.0228193 , -0.0916788 , 0 , 0 },
		{0 , 1. , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{2.66667 , -2.66667 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0}
	},
	/*m008[10][10]*/{
		{0.216897 , 0 , 0 , 0 , -0.179308 , -0.072956 , 0.0240377 , 0.0113296 , 0 , 0 },
		{0.862347 , 0 , 0 , 0 , -0.913494 , 0.0873264 , -0.0570826 , 0.0209031 , 0 , 0 },
		{3.02925 , 0 , 0 , 0 , -4.01431 , -0.294511 , 0.459976 , 0.819591 , 0 , 0 },
		{2.13991 , 0 , 0 , 0 , -2.67884 , 0.231816 , 0.374068 , -0.0669542 , 0 , 0 },
		{76.0879 , 0 , 0 , 0 , -77.7679 , -6.0904 , 2.2128 , 5.55764 , 0 , 0 },
		{32.4232 , 0 , 0 , 0 , -35.7924 , 2.96414 , 0.0702211 , 0.334917 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{1. , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0}
	},
	/*m009[10][10]*/{
		{0 , 0 , -0.0328407 , -0.040404 , 0.00212648 , -0.0288998 , -0.0173508 , -0.000992832 , 0.118362 , 0 },
		{0 , 0 , -0.0985222 , 0.0606061 , 0.0108335 , 0.0345923 , 0.0412032 , -0.00183178 , -0.0468811 , 0 },
		{0 , 0 , 0 , 0 , 0.0476072 , -0.116664 , -0.332019 , -0.0718224 , 0.472898 , 0 },
		{0 , 0 , 0 , 0 , 0.0317694 , 0.0918285 , -0.270009 , 0.00586734 , 0.140544 , 0 },
		{0 , 0 , 0 , 0 , 0.922279 , -2.41257 , -1.59723 , -0.487028 , 3.57455 , 0 },
		{0 , 0 , 0 , 0 , 0.424476 , 1.17418 , -0.0506868 , -0.0293495 , -1.51862 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1. , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0}
	},
	/*m0010[10][10]*/{
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1.}
	}
	};

/*----------------------------------------*/

	double m10[10][10][10]={
	/*m101[10][10]*/{
		{0 , 0 , -1.98687 , 0.730099 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , -2.96062 , -4.09515 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m102[10][10]*/{
		{0 , 0 , -1.32458 , -0.243366 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , -1.97375 , 1.36505 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m103[10][10]*/{
		{0 , 0 , -0.0630753 , 0.0270407 , 0.0680822 , -0.190991 , -0.0678562 , -0.208007 , 0 , 0 },
		{0 , 0 , -0.0939879 , -0.151672 , -0.232688 , 0.228774 , 0.145536 , -0.476004 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.633194 , -0.979258 , -5.02512 , -0.724774 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.0769351 , -0.696071 , 0.363374 , -2.53039 , 0 , 0 },
		{0 , 0 , 0 , 0 , -28.9856 , -60.1823 , -16.284 , -100.862 , 0 , 0 },
		{0 , 0 , 0 , 0 , 11.5767 , -3.23026 , 23.2477 , -35.8962 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m104[10][10]*/{
		{0 , 0 , -0.094613 , -0.0811221 , -0.0244675 , 0.055475 , -0.415831 , 0.195974 , 0 , 0 },
		{0 , 0 , -0.140982 , 0.455016 , 0.0836237 , -0.0664493 , 0.89186 , 0.448468 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.227558 , 0.284434 , -30.7946 , 0.682845 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.0276491 , 0.20218 , 2.2268 , 2.384 , 0 , 0 },
		{0 , 0 , 0 , 0 , 10.4169 , 17.4805 , -99.7905 , 95.0275 , 0 , 0 },
		{0 , 0 , 0 , 0 , -4.16046 , 0.938257 , 142.465 , 33.8196 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m105[10][10]*/{
		{0 , 0 , 0.0157688 , -0.00676017 , -0.00966156 , 0.0319807 , 0.00306429 , 0.0105632 , 0 , 0 },
		{0 , 0 , 0.023497 , 0.037918 , 0.0330208 , -0.0383072 , -0.00657219 , 0.0241729 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.0898568 , 0.163973 , 0.226928 , 0.0368061 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.0109179 , 0.116554 , -0.0164094 , 0.1285 , 0 , 0 },
		{0 , 0 , 0 , 0 , 4.11336 , 10.0773 , 0.735364 , 5.12208 , 0 , 0 },
		{0 , 0 , 0 , 0 , -1.64286 , 0.540894 , -1.04983 , 1.82291 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m106[10][10]*/{
		{0 , 0 , 0.0236532 , 0.0202805 , -0.0111735 , -0.0522173 , 0.0320749 , -0.0142895 , 0 , 0 },
		{0 , 0 , 0.0352454 , -0.113754 , 0.0381882 , 0.0625471 , -0.0687931 , -0.0327002 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.103918 , -0.267731 , 2.37532 , -0.04979 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.0126264 , -0.190307 , -0.171763 , -0.173831 , 0 , 0 },
		{0 , 0 , 0 , 0 , 4.75706 , -16.4539 , 7.69728 , -6.92897 , 0 , 0 },
		{0 , 0 , 0 , 0 , -1.89995 , -0.883159 , -10.9889 , -2.46597 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m107[10][10]*/{
		{0.00214516 , -1.44984 , 0.851517 , 0.0521499 , 0.670706 , 0.121996 , -0.0577553 , 0.0354984 , 0 , 0 },
		{9.93719 , -7.48777 , 1.26884 , -0.292511 , -2.29231 , -0.146129 , 0.123872 , 0.0812348 , 0 , 0 },
		{-5.07366 , 7.68709 , 0 , 0 , 6.23786 , 0.625502 , -4.2771 , 0.12369 , 0 , 0 },
		{-8.68399 , 8.5586 , 0 , 0 , 0.75792 , 0.444616 , 0.309283 , 0.431835 , 0 , 0 },
		{2421.27 , -1982.64 , 0 , 0 , -285.55 , 38.4415 , -13.86 , 17.2132 , 0 , 0 },
		{-639.781 , 576.021 , 0 , 0 , 114.047 , 2.06333 , 19.7871 , 6.12604 , 0 , 0 },
		{0 , 7.81516 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{17.9842 , -18.7604 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m108[10][10]*/{
		{0.000804435 , 0 , 0 , 0 , 0.943528 , -0.280446 , -0.177728 , -0.129681 , 0 , 0 },
		{3.72645 , 0 , 0 , 0 , -3.22474 , 0.335924 , 0.381186 , -0.296763 , 0 , 0 },
		{-1.90262 , 0 , 0 , 0 , 8.77523 , -1.43791 , -13.1618 , -0.451857 , 0 , 0 },
		{-3.25649 , 0 , 0 , 0 , 1.06622 , -1.02209 , 0.951746 , -1.57756 , 0 , 0 },
		{907.978 , 0 , 0 , 0 , -401.702 , -88.3699 , -42.651 , -62.8823 , 0 , 0 },
		{-239.918 , 0 , 0 , 0 , 160.438 , -4.74322 , 60.8902 , -22.3794 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{6.74407 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m109[10][10]*/{
		{0 , 0 , 0.195751 , -0.0442484 , -0.0111896 , -0.111092 , 0.128287 , 0.0113642 , -0.359616 , 0 },
		{0 , 0 , 0.291686 , 0.248191 , 0.0382434 , 0.133069 , -0.275146 , 0.026006 , -0.879366 , 0 },
		{0 , 0 , 0 , 0 , -0.104069 , -0.569596 , 9.50039 , 0.0395972 , -0.485592 , 0 },
		{0 , 0 , 0 , 0 , -0.0126447 , -0.404877 , -0.686986 , 0.138245 , 0.417243 , 0 },
		{0 , 0 , 0 , 0 , 4.76394 , -35.0057 , 30.7862 , 5.51051 , 62.3651 , 0 },
		{0 , 0 , 0 , 0 , -1.90269 , -1.87892 , -43.9516 , 1.96115 , 54.4557 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m1010[10][10]*/{
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	}};
	
/*------------------*/

	double m11[10][10][10]={
	/*m111[10][10]*/{
		{0 , 0 , 1.98687 , -0.730099 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 5.96062 , 1.09515 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m112[10][10]*/{
		{0 , 0 , 0.657915 , 0.910033 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 1.97375 , -1.36505 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m113[10][10]*/{
		{0 , 0 , -0.180311 , -1.08877 , 0.322029 , 1.39533 , 0.108473 , -0.12195 , 0 , 0 },
		{0 , 0 , -0.540933 , 1.63315 , 1.64059 , -1.67018 , -0.257592 , -0.224998 , 0 , 0 },
		{0 , 0 , 0 , 0 , 7.20951 , 5.63273 , 2.0757 , -8.82197 , 0 , 0 },
		{0 , 0 , 0 , 0 , 4.81108 , -4.43364 , 1.68803 , 0.720687 , 0 , 0 },
		{0 , 0 , 0 , 0 , 139.668 , 116.483 , 9.9855 , -59.8218 , 0 , 0 },
		{0 , 0 , 0 , 0 , 64.2815 , -56.6914 , 0.316881 , -3.60501 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m114[10][10]*/{
		{0 , 0 , 0.740116 , -1.35098 , -0.821069 , 0.596066 , 0.767049 , 0.4334 , 0 , 0 },
		{0 , 0 , 2.22035 , 2.02647 , -4.18298 , -0.713476 , -1.82152 , 0.799626 , 0 , 0 },
		{0 , 0 , 0 , 0 , -18.3819 , 2.40622 , 14.678 , 31.3526 , 0 , 0 },
		{0 , 0 , 0 , 0 , -12.2667 , -1.89398 , 11.9366 , -2.56126 , 0 , 0 },
		{0 , 0 , 0 , 0 , -356.107 , 49.7599 , 70.6109 , 212.602 , 0 , 0 },
		{0 , 0 , 0 , 0 , -163.897 , -24.2177 , 2.24077 , 12.8119 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m115[10][10]*/{
		{0 , 0 , 0.0133317 , 0.124044 , -0.0327528 , -0.157615 , -0.00846373 , 0.0165001 , 0 , 0 },
		{0 , 0 , 0.0399951 , -0.186066 , -0.166861 , 0.188661 , 0.0200989 , 0.0304428 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.733261 , -0.636264 , -0.161959 , 1.19363 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.489323 , 0.500817 , -0.13171 , -0.0975107 , 0 , 0 },
		{0 , 0 , 0 , 0 , -14.2052 , -13.1577 , -0.779131 , 8.09404 , 0 , 0 },
		{0 , 0 , 0 , 0 , -6.53792 , 6.40376 , -0.024725 , 0.487766 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m116[10][10]*/{
		{0 , 0 , -0.0871453 , 0.127868 , 0.0823762 , -0.024661 , -0.0620538 , -0.0347122 , 0 , 0 },
		{0 , 0 , -0.261436 , -0.191801 , 0.41967 , 0.0295185 , 0.14736 , -0.0640442 , 0 , 0 },
		{0 , 0 , 0 , 0 , 1.84422 , -0.0995522 , -1.18744 , -2.51111 , 0 , 0 },
		{0 , 0 , 0 , 0 , 1.23069 , 0.0783596 , -0.965664 , 0.205138 , 0 , 0 },
		{0 , 0 , 0 , 0 , 35.7275 , -2.05871 , -5.71238 , -17.0279 , 0 , 0 },
		{0 , 0 , 0 , 0 , 16.4435 , 1.00196 , -0.181277 , -1.02614 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m117[10][10]*/{
		{-4.35186 , 3.06463 , 1.51692 , -0.501267 , 0.393366 , -0.624545 , 0.226779 , 0.0495597 , 0 , 0 },
		{-17.3023 , 8.50271 , 4.55076 , 0.7519 , 2.00402 , 0.747564 , -0.538534 , 0.0914379 , 0 , 0 },
		{-60.7794 , 41.2459 , 0 , 0 , 8.80659 , -2.52118 , 4.33955 , 3.58519 , 0 , 0 },
		{-42.9356 , 30.0198 , 0 , 0 , 5.87685 , 1.98447 , 3.52907 , -0.292883 , 0 , 0 },
		{-1526.64 , 1168.11 , 0 , 0 , 170.607 , -52.1373 , 20.8762 , 24.3112 , 0 , 0 },
		{-650.544 , 466.256 , 0 , 0 , 78.5215 , 25.3748 , 0.662487 , 1.46505 , 0 , 0 },
		{0 , -7.81516 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{-20.0642 , 20.8404 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m118[10][10]*/{
		{-1.46277 , 0 , 0.468732 , 2.65964 , 0.644814 , -3.07871 , 0.0599697 , 0.35185 , 0 , 0 },
		{-5.81573 , 0 , 1.40619 , -3.98945 , 3.28504 , 3.68514 , -0.142411 , 0.649165 , 0 , 0 },
		{-20.4295 , 0 , 0 , 0 , 14.436 , -12.4282 , 1.14756 , 25.4531 , 0 , 0 },
		{-14.4317 , 0 , 0 , 0 , 9.63346 , 9.78253 , 0.933232 , -2.07932 , 0 , 0 },
		{-513.142 , 0 , 0 , 0 , 279.663 , -257.012 , 5.52053 , 172.598 , 0 , 0 },
		{-218.664 , 0 , 0 , 0 , 128.714 , 125.086 , 0.175189 , 10.4012 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{-6.74407 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m119[10][10]*/{
		{0 , 0 , 0.291812 , 0.0483636 , -0.0330782 , -0.0269467 , 0.0199899 , -0.109398 , 0 , 0 },
		{0 , 0 , 0.875436 , -0.0725455 , -0.168519 , 0.0322545 , -0.0474702 , -0.201839 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.740548 , -0.108779 , 0.382519 , -7.91392 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.494185 , 0.0856224 , 0.311077 , 0.646506 , 0 , 0 },
		{0 , 0 , 0 , 0 , -14.3464 , -2.24952 , 1.84017 , -53.6643 , 0 , 0 },
		{0 , 0 , 0 , 0 , -6.60289 , 1.09482 , 0.0583962 , -3.23394 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m1110[10][10]*/{
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	}};
	    
/*-----------------------------*/

	double m20[10][10][10]={ 
	/*m201[10][10]*/{
		{0 , 0 , -6.98112 , 15.6429 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , -13.4086 , -52.1664 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m202[10][10]*/{
		{0 , 0 , -4.65408 , -5.2143 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , -8.93903 , 17.3888 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m203[10][10]*/{
		{0 , 0 , -0.221623 , 0.579366 , 1.17316 , -1.38853 , 0.135066 , -1.02839 , 0 , 0 },
		{0 , 0 , -0.425668 , -1.93209 , 2.92926 , 3.32871 , 2.69087 , -0.856585 , 0 , 0 },
		{0 , 0 , 0 , 0 , 15.6719 , -10.1006 , -0.0487765 , -11.775 , 0 , 0 },
		{0 , 0 , 0 , 0 , 15.9377 , 1.94304 , 14.6555 , 15.6499 , 0 , 0 },
		{0 , 0 , 0 , 0 , 39.256 , -955.766 , -18.5938 , -1469.65 , 0 , 0 },
		{0 , 0 , 0 , 0 , 324.798 , 59.1356 , 65.9104 , 599.472 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m204[10][10]*/{
		{0 , 0 , -0.332434 , -1.7381 , -0.421611 , 0.40331 , 0.827699 , 0.968895 , 0 , 0 },
		{0 , 0 , -0.638502 , 5.79627 , -1.05272 , -0.966852 , 16.49 , 0.807031 , 0 , 0 },
		{0 , 0 , 0 , 0 , -5.63218 , 2.93381 , -0.298908 , 11.0938 , 0 , 0 },
		{0 , 0 , 0 , 0 , -5.72772 , -0.564372 , 89.8106 , -14.7446 , 0 , 0 },
		{0 , 0 , 0 , 0 , -14.1079 , 277.61 , -113.945 , 1384.63 , 0 , 0 },
		{0 , 0 , 0 , 0 , -116.727 , -17.1764 , 403.907 , -564.792 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m205[10][10]*/{
		{0 , 0 , 0.0554057 , -0.144842 , -0.166483 , 0.232504 , -0.00609938 , 0.0522244 , 0 , 0 },
		{0 , 0 , 0.106417 , 0.483023 , -0.415692 , -0.557379 , -0.121516 , 0.0434998 , 0 , 0 },
		{0 , 0 , 0 , 0 , -2.224 , 1.69131 , 0.00220268 , 0.59797 , 0 , 0 },
		{0 , 0 , 0 , 0 , -2.26172 , -0.325354 , -0.661821 , -0.794747 , 0 , 0 },
		{0 , 0 , 0 , 0 , -5.57083 , 160.039 , 0.839672 , 74.6328 , 0 , 0 },
		{0 , 0 , 0 , 0 , -46.0922 , -9.90202 , -2.97642 , -30.4429 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m206[10][10]*/{
		{0 , 0 , 0.0831086 , 0.434525 , -0.192536 , -0.379626 , -0.0638441 , -0.0706474 , 0 , 0 },
		{0 , 0 , 0.159626 , -1.44907 , -0.480743 , 0.910075 , -1.27195 , -0.058845 , 0 , 0 },
		{0 , 0 , 0 , 0 , -2.57203 , -2.76153 , 0.0230561 , -0.808912 , 0 , 0 },
		{0 , 0 , 0 , 0 , -2.61566 , 0.531229 , -6.92748 , 1.07511 , 0 , 0 },
		{0 , 0 , 0 , 0 , -6.4426 , -261.308 , 8.7891 , -100.961 , 0 , 0 },
		{0 , 0 , 0 , 0 , -53.3051 , 16.1678 , -31.1551 , 41.182 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m207[10][10]*/{
		{-91.2189 , 81.3165 , 2.99191 , 1.11735 , 11.5573 , 0.886925 , 0.11496 , 0.175504 , 0 , 0 },
		{-212.414 , 167.658 , 5.74652 , -3.72617 , 28.8574 , -2.12622 , 2.29032 , 0.146185 , 0 , 0 },
		{-527.503 , 450.594 , 0 , 0 , 154.39 , 6.45179 , -0.0415157 , 2.00952 , 0 , 0 },
		{-749.114 , 579.46 , 0 , 0 , 157.009 , -1.24112 , 12.4739 , -2.67081 , 0 , 0 },
		{7764.36 , -5063.44 , 0 , 0 , 386.727 , 610.496 , -15.826 , 250.809 , 0 , 0 },
		{-22612.5 , 19501.2 , 0 , 0 , 3199.72 , -37.7729 , 56.0992 , -102.306 , 0 , 0 },
		{0 , 44.4252 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{15.4051 , -18.7662 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m208[10][10]*/{
		{-34.2071 , 0 , 0 , 0 , 16.2584 , -2.03888 , 0.353763 , -0.641144 , 0 , 0 },
		{-79.6551 , 0 , 0 , 0 , 40.5956 , 4.88778 , 7.04792 , -0.534035 , 0 , 0 },
		{-197.814 , 0 , 0 , 0 , 217.191 , -14.8315 , -0.127755 , -7.3411 , 0 , 0 },
		{-280.918 , 0 , 0 , 0 , 220.875 , 2.8531 , 38.3855 , 9.75688 , 0 , 0 },
		{2911.64 , 0 , 0 , 0 , 544.035 , -1403.42 , -48.7008 , -916.245 , 0 , 0 },
		{-8479.69 , 0 , 0 , 0 , 4501.27 , 86.833 , 172.632 , 373.738 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{5.77691 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m209[10][10]*/{
		{0 , 0 , 0.687795 , -0.948054 , -0.192814 , -0.807654 , -0.255352 , 0.0561848 , -0.64363 , 0 },
		{0 , 0 , 1.32104 , 3.1616 , -0.481439 , 1.93618 , -5.08731 , 0.0467986 , -13.5825 , 0 },
		{0 , 0 , 0 , 0 , -2.57575 , -5.87514 , 0.0922157 , 0.643316 , 7.77564 , 0 },
		{0 , 0 , 0 , 0 , -2.61944 , 1.13019 , -27.7073 , -0.855015 , 16.0333 , 0 },
		{0 , 0 , 0 , 0 , -6.45192 , -555.931 , 35.1531 , 80.2925 , 102.043 , 0 },
		{0 , 0 , 0 , 0 , -53.3822 , 34.3969 , -124.609 , -32.7515 , -98.8845 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m2010[10][10]*/{
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	}};
	
/*------------------------*/

	double m21[10][10][10]={ 
	/*m211[10][10]*/{
		{0 , 0 , -11.843 , -0.799566 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , -17.6471 , 4.48479 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m212[10][10]*/{
		{0 , 0 , -3.92158 , 0.996621 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , -5.8435 , -5.59008 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m213[10][10]*/{
		{0 , 0 , 1.07476 , -1.19236 , -1.69453 , 5.36372 , -0.80202 , 1.39587 , 0 , 0 },
		{0 , 0 , 1.6015 , 6.68801 , 5.79148 , -6.42479 , 1.72015 , 3.19432 , 0 , 0 },
		{0 , 0 , 0 , 0 , -15.7599 , 27.5011 , -59.394 , 4.86374 , 0 , 0 },
		{0 , 0 , 0 , 0 , -1.91488 , 19.5482 , 4.29486 , 16.9807 , 0 , 0 },
		{0 , 0 , 0 , 0 , 721.438 , 1690.14 , -192.468 , 676.858 , 0 , 0 },
		{0 , 0 , 0 , 0 , -288.139 , 90.7174 , 274.774 , 240.889 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m214[10][10]*/{
		{0 , 0 , -4.41155 , -1.47952 , 4.32051 , 2.2913 , -5.67136 , -4.96081 , 0 , 0 },
		{0 , 0 , -6.5736 , 8.29867 , -14.7664 , -2.74457 , 12.1637 , -11.3524 , 0 , 0 },
		{0 , 0 , 0 , 0 , 40.1826 , 11.7481 , -419.995 , -17.2853 , 0 , 0 },
		{0 , 0 , 0 , 0 , 4.88232 , 8.35068 , 30.3704 , -60.3479 , 0 , 0 },
		{0 , 0 , 0 , 0 , -1839.43 , 722.001 , -1361. , -2405.5 , 0 , 0 },
		{0 , 0 , 0 , 0 , 734.661 , 38.7531 , 1943.02 , -856.099 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m215[10][10]*/{
		{0 , 0 , -0.0794652 , 0.135847 , 0.172347 , -0.605877 , 0.0625786 , -0.188865 , 0 , 0 },
		{0 , 0 , -0.11841 , -0.761969 , -0.589038 , 0.725734 , -0.134216 , -0.4322 , 0 , 0 },
		{0 , 0 , 0 , 0 , 1.6029 , -3.10648 , 4.63429 , -0.658076 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.194757 , -2.20813 , -0.335112 , -2.29752 , 0 , 0 },
		{0 , 0 , 0 , 0 , -73.3757 , -190.915 , 15.0175 , -91.5805 , 0 , 0 },
		{0 , 0 , 0 , 0 , 29.3059 , -10.2473 , -21.4396 , -32.5928 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m216[10][10]*/{
		{0 , 0 , 0.51944 , 0.140034 , -0.433468 , -0.0947977 , 0.458809 , 0.397325 , 0 , 0 },
		{0 , 0 , 0.774012 , -0.785455 , 1.48148 , 0.113551 , -0.984038 , 0.909241 , 0 , 0 },
		{0 , 0 , 0 , 0 , -4.03144 , -0.486051 , 33.9773 , 1.38443 , 0 , 0 },
		{0 , 0 , 0 , 0 , -0.489833 , -0.345492 , -2.45695 , 4.83342 , 0 , 0 },
		{0 , 0 , 0 , 0 , 184.547 , -29.8713 , 110.104 , 192.663 , 0 , 0 },
		{0 , 0 , 0 , 0 , -73.7071 , -1.60333 , -157.189 , 68.5673 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m217[10][10]*/{
		{-0.0161403 , 11.3308 , -9.04178 , -0.548961 , -2.06991 , -2.40077 , -1.67674 , -0.567273 , 0 , 0 },
		{-74.7681 , 58.5182 , -13.4731 , 3.07914 , 7.07444 , 2.8757 , 3.59622 , -1.29815 , 0 , 0 },
		{38.1746 , -60.0759 , 0 , 0 , -19.2511 , -12.3094 , -124.172 , -1.97659 , 0 , 0 },
		{65.3389 , -66.8869 , 0 , 0 , -2.33907 , -8.74966 , 8.97905 , -6.90083 , 0 , 0 },
		{-18217.8 , 15494.7 , 0 , 0 , 881.254 , -756.496 , -402.382 , -275.071 , 0 , 0 },
		{4813.75 , -4501.7 , 0 , 0 , -351.969 , -40.6047 , 574.456 , -97.8957 , 0 , 0 },
		{0 , -61.0768 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{-135.314 , 146.616 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m218[10][10]*/{
		{-0.00542516 , 0 , -2.79393 , 2.9127 , -3.39304 , -11.8347 , -0.4434 , -4.02737 , 0 , 0 },
		{-25.1314 , 0 , -4.16321 , -16.3374 , 11.5966 , 14.1759 , 0.950989 , -9.21626 , 0 , 0 },
		{12.8314 , 0 , 0 , 0 , -31.5568 , -60.6794 , -32.8362 , -14.0329 , 0 , 0 },
		{21.962 , 0 , 0 , 0 , -3.83425 , -43.1317 , 2.37443 , -48.9926 , 0 , 0 },
		{-6123.46 , 0 , 0 , 0 , 1444.57 , -3729.17 , -106.406 , -1952.87 , 0 , 0 },
		{1618.02 , 0 , 0 , 0 , -576.955 , -200.162 , 151.91 , -695.012 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{-45.4824 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m219[10][10]*/{
		{0 , 0 , -1.73938 , 0.0529653 , 0.174059 , -0.103584 , -0.1478 , 1.25219 , 0 , 0 },
		{0 , 0 , -2.59183 , -0.297084 , -0.594891 , 0.124075 , 0.316996 , 2.86553 , 0 , 0 },
		{0 , 0 , 0 , 0 , 1.61883 , -0.531101 , -10.9454 , 4.36311 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0.196693 , -0.377514 , 0.791476 , 15.2328 , 0 , 0 },
		{0 , 0 , 0 , 0 , -74.1049 , -32.6399 , -35.4688 , 607.188 , 0 , 0 },
		{0 , 0 , 0 , 0 , 29.5971 , -1.75193 , 50.6366 , 216.094 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m2110[10][10]*/{
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	}};
	
/*--------------------------*/

	double m22[10][10][10]={ 
	/*m221[10][10]*/{
		{0 , 0 , 18.8241 , -14.8433 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 56.4723 , 22.265 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m222[10][10]*/{
		{0 , 0 , 4.92751 , 7.86582 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 14.7825 , -11.7987 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m223[10][10]*/{
		{0 , 0 , 0.663515 , -12.8257 , -4.84529 , 10.7966 , 1.68798 , 1.12844 , 0 , 0 },
		{0 , 0 , 1.99055 , 19.2386 , -24.6846 , -12.9233 , -4.00845 , 2.08197 , 0 , 0 },
		{0 , 0 , 0 , 0 , -108.475 , 43.5841 , 32.3005 , 81.6323 , 0 , 0 },
		{0 , 0 , 0 , 0 , -72.3881 , -34.3059 , 26.2678 , -6.66873 , 0 , 0 },
		{0 , 0 , 0 , 0 , -2101.46 , 901.305 , 155.387 , 553.549 , 0 , 0 },
		{0 , 0 , 0 , 0 , -967.189 , -438.657 , 4.93107 , 33.3582 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m224[10][10]*/{
		{0 , 0 , 2.70469 , -28.4843 , -2.29684 , 29.6401 , 5.93931 , 2.70068 , 0 , 0 },
		{0 , 0 , 8.11407 , 42.7264 , -11.7014 , -35.4784 , -14.1041 , 4.98277 , 0 , 0 },
		{0 , 0 , 0 , 0 , -51.4211 , 119.652 , 113.652 , 195.37 , 0 , 0 },
		{0 , 0 , 0 , 0 , -34.3145 , -94.1806 , 92.4259 , -15.9602 , 0 , 0 },
		{0 , 0 , 0 , 0 , -996.165 , 2474.37 , 546.745 , 1324.8 , 0 , 0 },
		{0 , 0 , 0 , 0 , -458.482 , -1204.25 , 17.3504 , 79.8358 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m225[10][10]*/{
		{0 , 0 , -0.122007 , 0.839185 , 0.541041 , -0.515314 , -0.120171 , -0.142009 , 0 , 0 },
		{0 , 0 , -0.36602 , -1.25878 , 2.75636 , 0.616818 , 0.285371 , -0.262006 , 0 , 0 },
		{0 , 0 , 0 , 0 , 12.1127 , -2.08024 , -2.29954 , -10.273 , 0 , 0 },
		{0 , 0 , 0 , 0 , 8.08309 , 1.6374 , -1.87007 , 0.839227 , 0 , 0 },
		{0 , 0 , 0 , 0 , 234.656 , -43.0187 , -11.0624 , -69.6614 , 0 , 0 },
		{0 , 0 , 0 , 0 , 107.999 , 20.9368 , -0.351054 , -4.19797 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m226[10][10]*/{
		{0 , 0 , -0.774763 , 2.37177 , 0.576243 , -2.08569 , -0.655564 , -0.230319 , 0 , 0 },
		{0 , 0 , -2.32429 , -3.55766 , 2.9357 , 2.49652 , 1.55677 , -0.424939 , 0 , 0 },
		{0 , 0 , 0 , 0 , 12.9008 , -8.41959 , -12.5446 , -16.6615 , 0 , 0 },
		{0 , 0 , 0 , 0 , 8.60901 , 6.62723 , -10.2017 , 1.36111 , 0 , 0 },
		{0 , 0 , 0 , 0 , 249.923 , -174.115 , -60.3481 , -112.981 , 0 , 0 },
		{0 , 0 , 0 , 0 , 115.026 , 84.74 , -1.91509 , -6.80853 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m227[10][10]*/{
		{7.90885 , -6.52974 , 7.73725 , -8.85138 , -3.90022 , -3.36511 , 3.63242 , 1.41727 , 0 , 0 },
		{31.4443 , -18.1165 , 23.2117 , 13.2771 , -19.8699 , 4.02795 , -8.62594 , 2.61487 , 0 , 0 },
		{110.458 , -87.8817 , 0 , 0 , -87.3171 , -13.5844 , 69.5086 , 102.527 , 0 , 0 },
		{78.0291 , -63.9626 , 0 , 0 , -58.2688 , 10.6925 , 56.5267 , -8.37564 , 0 , 0 },
		{2774.44 , -2488.86 , 0 , 0 , -1691.57 , -280.921 , 334.383 , 695.234 , 0 , 0 },
		{1182.27 , -993.442 , 0 , 0 , -778.539 , 136.722 , 10.6114 , 41.8965 , 0 , 0 },
		{0 , 16.6516 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{36.4636 , -44.4043 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m228[10][10]*/{
		{8.612 , 0 , -8.57734 , 9.92405 , 15.2774 , 28.49 , -10.4305 , -3.43545 , 0 , 0 },
		{34.2399 , 0 , -25.732 , -14.8861 , 77.8316 , -34.1018 , 24.7694 , -6.33842 , 0 , 0 },
		{120.278 , 0 , 0 , 0 , 342.027 , 115.009 , -199.594 , -248.523 , 0 , 0 },
		{84.9664 , 0 , 0 , 0 , 228.243 , -90.5262 , -162.317 , 20.3024 , 0 , 0 },
		{3021.11 , 0 , 0 , 0 , 6625.99 , 2378.36 , -960.183 , -1685.24 , 0 , 0 },
		{1287.38 , 0 , 0 , 0 , 3049.59 , -1157.53 , -30.4705 , -101.556 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{39.7055 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m229[10][10]*/{
		{0 , 0 , 4.15307 , -0.462666 , -0.340353 , -1.03257 , 0.0809144 , 0.216668 , 0 , 0 },
		{0 , 0 , 12.4592 , 0.693998 , -1.73395 , 1.23595 , -0.192148 , 0.399754 , 0 , 0 },
		{0 , 0 , 0 , 0 , -7.61975 , -4.16829 , 1.54835 , 15.674 , 0 , 0 },
		{0 , 0 , 0 , 0 , -5.08484 , 3.28095 , 1.25917 , -1.28044 , 0 , 0 },
		{0 , 0 , 0 , 0 , -147.615 , -86.199 , 7.44859 , 106.285 , 0 , 0 },
		{0 , 0 , 0 , 0 , -67.9394 , 41.9523 , 0.236374 , 6.405 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	},
	    
	/*m2210[10][10]*/{
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
		{0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }
	}};

/*--------------------------*/

	int je,ke,le;
	double ai[10]={14./23., 16./23., 6./23., -12./23., 0.408619, -0.422989, -0.899395, 0.145649, -1., -1.}; 
	double U0[10][10],U1[10][10],U2[10][10];


	for (ie=0; ie<10; ie++) for (je=0; je<10; je++) U0[ie][je]=U1[ie][je]=U2[ie][je]=0.;

	for (ke=0; ke<10; ke++) for (le=0; le<10; le++) for (ie=0; ie<10; ie++)
	{
		U0[ke][le]=U0[ke][le] + m00[ke][le][ie]*pow(eta_mu,ai[ie]);
		U1[ke][le]=U1[ke][le] + m10[ke][le][ie]*pow(eta_mu,ai[ie])+m11[ke][le][ie]*pow(eta_mu,ai[ie]-1.);
		U2[ke][le]=U2[ke][le] + m20[ke][le][ie]*pow(eta_mu,ai[ie])+m21[ke][le][ie]*pow(eta_mu,ai[ie]-1.)+m22[ke][le][ie]*pow(eta_mu,ai[ie]-2.);
	}

	for (ie=1; ie<=8; ie++) for (je=1; je<=8; je++)
	{
		if (je<=6)
		{
			C0b[ie]=C0b[ie]+U0[ie-1][je-1]*C0w[je];
			C1b[ie]=C1b[ie]+eta_mu*(U0[ie-1][je-1]*C1w[je]+U1[ie-1][je-1]*C0w[je]);
			C2b[ie]=C2b[ie]+eta_mu*eta_mu*(U0[ie-1][je-1]*C2w[je]+U1[ie-1][je-1]*C1w[je]+U2[ie-1][je-1]*C0w[je]);
		}
		if (je==7)
		{
			C0b[ie]=C0b[ie]+U0[ie-1][je-1]*C0w7;
			C1b[ie]=C1b[ie]+eta_mu*(U0[ie-1][je-1]*C1w7+U1[ie-1][je-1]*C0w7);
			C2b[ie]=C2b[ie]+eta_mu*eta_mu*(U0[ie-1][je-1]*C2w7+U1[ie-1][je-1]*C1w7+U2[ie-1][je-1]*C0w7);
		}
		if (je==8)
		{
			C0b[ie]=C0b[ie]+U0[ie-1][je-1]*C0w8;
			C1b[ie]=C1b[ie]+eta_mu*(U0[ie-1][je-1]*C1w8+U1[ie-1][je-1]*C0w8);
			C2b[ie]=C2b[ie]+eta_mu*eta_mu*(U0[ie-1][je-1]*C2w8+U1[ie-1][je-1]*C1w8+U2[ie-1][je-1]*C0w8);
		}
	}

	for (je=1; je<=8; je++)
	{
		C0b[9]=C0b[9]+4.*pi/alphas_mu*(U0[9-1][je-1]*C0w[je]);
		C1b[9]=C1b[9]+4.*pi/alphas_mu*(eta_mu*(U0[9-1][je-1]*C1w[je]+U1[9-1][je-1]*C0w[je]));
		C2b[9]=C2b[9]+4.*pi/alphas_mu*(eta_mu*eta_mu*(U0[9-1][je-1]*C2w[je]+U1[9-1][je-1]*C1w[je]+U2[9-1][je-1]*C0w[je]));
	}
	C1b[9]=C1b[9]+4.*pi/alphas_mu*(eta_mu*(U0[9-1][9-1]*C0w[je]));
	C2b[9]=C2b[9]+4.*pi/alphas_mu*(eta_mu*eta_mu*(U0[9-1][9-1]*C1w[je]+U1[9-1][9-1]*C0w[je]));

	C0b[10]=C0w[10];
	C1b[10]=eta_mu*C1w[10];

	return;	
}

/*-----------------------------------------------------------------------------*/

void C_calculator_base2(double C0w[], double C1w[], double mu_W, double C0b[], double C1b[], double mu, struct parameters* param)
/* calculates the LO (C0b) and NLO (C1b) contributions to the Wilson coefficients at scale mu, using the LO (C0w) and NLO (C1w) contributions to the Wilson coefficients at scale mu_W and the parameters of the structure param, in the traditional operator basis */
{
	int ie;
	for(ie=1;ie<=10;ie++) C0b[ie]=C1b[ie]=0.;

	double alphas_muW=alphas_running(mu_W,param->mass_top_pole,param->mass_b_pole,param);
	double alphas_mu=alphas_running(mu,param->mass_top_pole,param->mass_b_pole,param);	
	double eta_mu=alphas_muW/alphas_mu;
	
	double C0w7= C0w[7]-1./3.*C0w[5]-C0w[6]; 
	double C1w7= C1w[7]-1./3.*C1w[5]-C1w[6]; 

	double C0w8= C0w[8]+C0w[5]; 
	double C1w8= C1w[8]+C1w[5]; 

 	C0b[1]= (1./2.*pow(eta_mu,6./23.) -1./2.*pow(eta_mu,-12./23.))*C0w[2];
	C0b[2]= (1./2.*pow(eta_mu,6./23.) +1./2.*pow(eta_mu,-12./23.))*C0w[2];
	C0b[3]= (-1./14.*pow(eta_mu,6./23.) +1./6.*pow(eta_mu,-12./23.) +0.0509*pow(eta_mu,0.4086) -0.1403*pow(eta_mu,-0.4230) -0.01126*pow(eta_mu,-0.8994) +0.0054*pow(eta_mu,0.1456))*C0w[2];
	C0b[4]= (-1./14.*pow(eta_mu,6./23.) -1./6.*pow(eta_mu,-12./23.) +0.0984*pow(eta_mu,0.4086) +0.1214*pow(eta_mu,-0.4230) +0.0156*pow(eta_mu,-0.8994) +0.0026*pow(eta_mu,0.1456))*C0w[2];
	C0b[5]= (-0.0397*pow(eta_mu,0.4086) +0.0117*pow(eta_mu,-0.4230) -0.0025*pow(eta_mu,-0.8994) +0.0304*pow(eta_mu,0.1456))*C0w[2];
	C0b[6]= (0.0335*pow(eta_mu,0.4086) +0.0239*pow(eta_mu,-0.4230) -0.0462*pow(eta_mu,-0.8994) -0.0112*pow(eta_mu,0.1456))*C0w[2];
 
	C0b[7]= pow(eta_mu,16./23.)*C0w[7] + 8./3.*(pow(eta_mu,14./23.)-pow(eta_mu,16./23.))*C0w[8] + C0w[2] * (2.2996*pow(eta_mu,14./23.) -1.0880*pow(eta_mu,16./23.) -3./7.*pow(eta_mu,6./23.) -1./14.*pow(eta_mu,-12./23.) -0.6494*pow(eta_mu,0.4086) -0.0380*pow(eta_mu,-0.4230) -0.0185*pow(eta_mu,-0.8994) -0.0057*pow(eta_mu,0.1456));

	C0b[8]= pow(eta_mu,14./23.)*C0w[8] + C0w[2] * (0.8623*pow(eta_mu,14./23.) -0.9135*pow(eta_mu,0.4086) +0.0873*pow(eta_mu,-0.4230) -0.0571*pow(eta_mu,-0.8994) +0.0209*pow(eta_mu,0.1456));
	
	C0b[9]=C0w[9]+4.*pi/alphas_muW*(-4./33.*(1.-pow(eta_mu,11./23.))+8./87.*(1.-pow(eta_mu,29./23.)))*C0w[2];
	C0b[10]=C0w[10];

	C1b[1]= (C0w[2] *0.8136+1.0197*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)+(C0w[2] *0.7142+2.9524*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.);
	
	C1b[2]= (C0w[2] *0.8136+1.0197*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)-(C0w[2] *0.7142+2.9524*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.);
	
	C1b[3]= (-0.0766*C0w[2]-0.1457*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)+(-0.1455*C0w[2]-0.9841*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.)
	        +(0.1494*eta_mu*C1w[4]-0.8848*C0w[2]+0.2303*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(-0.3726*eta_mu*C1w[4]+0.4137*C0w[2]+1.4672*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(0.0738*eta_mu*C1w[4]-0.0114*C0w[2]+0.0971*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(-0.0173*eta_mu*C1w[4]+0.1722*C0w[2]-0.0213*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
	
	C1b[4]= (-0.2353*C0w[2]-0.1457*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)+(-0.0397*C0w[2]+0.9841*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.)
	        +(0.2885*eta_mu*C1w[4]+0.4920*C0w[2]+0.4447*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(0.3224*eta_mu*C1w[4]-0.2758*C0w[2]-1.2696*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(-0.1025*eta_mu*C1w[4]+0.0019*C0w[2]-0.1349*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(-0.0084*eta_mu*C1w[4]-0.1449*C0w[2]-0.0104*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
	
	C1b[5]= 0.0397*C0w[2]*pow(eta_mu,6./23.)+0.0926*C0w[2]*pow(eta_mu,-12./23.)
	        +(-0.1163*eta_mu*C1w[4]+0.7342*C0w[2]-0.1792*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(0.0310*eta_mu*C1w[4]-0.1262*C0w[2]-0.1221*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(0.0162*eta_mu*C1w[4]-0.1209*C0w[2]+0.0213*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(-0.0975*eta_mu*C1w[4]-0.1085*C0w[2]-0.1197*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
	
	C1b[6]= -0.1191*C0w[2]*pow(eta_mu,6./23.)-0.2778*C0w[2]*pow(eta_mu,-12./23.)
	        +(0.0982*eta_mu*C1w[4]-0.5544*C0w[2]+0.1513*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(0.0634*eta_mu*C1w[4]+0.1915*C0w[2]-0.2497*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(0.3026*eta_mu*C1w[4]-0.2744*C0w[2]+0.3983*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(0.0358*eta_mu*C1w[4]+0.3568*C0w[2]+0.0440*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);

	C1b[7]= pow(eta_mu,39./23.)*C1w7 +8./3.*(pow(eta_mu,37./23.)-pow(eta_mu,39./23.))*C1w8
		+(297664./14283.*pow(eta_mu,16./23.)-7164416./357075.*pow(eta_mu,14./23.)+256868./14283.*pow(eta_mu,37./23.)-6698884./357075.*pow(eta_mu,39./23.))*C0w[8]
		+37208./4761.*(pow(eta_mu,39./23.)-pow(eta_mu,16./23.))*C0w[7]
		+(4661194./816831.*eta_mu*C1w[4]-17.3023*C0w[2]+14.8088*eta_mu*C1w[1]/15.)*pow(eta_mu,14./23.)
		+(-8516./2217.*eta_mu*C1w[4]+8.5027*C0w[2]-10.8090*eta_mu*C1w[1]/15.)*pow(eta_mu,16./23.)
		+ (4.5508*C0w[2]-0.8740*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)
		+ (0.7519*C0w[2]+0.4218*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.)
		+ (-1.9043*eta_mu*C1w[4]+2.0040*C0w[2]-2.9347*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+ (-0.1008*eta_mu*C1w[4]+0.7476*C0w[2]+0.3971*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+ (0.1216*eta_mu*C1w[4]-0.5385*C0w[2]+0.1600*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+ (0.0183*eta_mu*C1w[4]+0.0914*C0w[2]+0.0225*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
	C1b[9]=eta_mu*(C1w[9]+4.*pi/alphas_muW*(-4./33.*(1.-pow(eta_mu,11./23.))+8./87.*(1.-pow(eta_mu,29./23.)))*C1w[2]);
	C1b[10]=eta_mu*C1w[10];
		
	return;	
}

/*----------------------------------------------------------------------*/

void Cprime_calculator(double Cpb[], double complex CQpb[], double mu_W, double mu, struct parameters* param)
{
	int ie;
	for(ie=1;ie<=10;ie++) Cpb[ie]=0.;
	for(ie=1;ie<=2;ie++) CQpb[ie]=0.;
	
	double alphas_muW=alphas_running(mu_W,param->mass_top_pole,param->mass_b_pole,param);
	double alphas_mu=alphas_running(mu,param->mass_top_pole,param->mass_b_pole,param);	
	double eta_mu=alphas_muW/alphas_mu;

	double mass_c_muW=running_mass(param->mass_c,param->mass_c,mu_W,param->mass_top_pole,param->mass_b_pole,param);
	double mass_b_muW=running_mass(param->mass_b,param->mass_b,mu_W,param->mass_top_pole,param->mass_b,param);
	double mass_top_muW=running_mass(param->mtmt,param->mtmt,mu_W,param->mass_top_pole,param->mass_b,param);
	
	/* Wilson coefficients C7 and C8 prime in the SM */ 
	double xt= pow(mass_top_muW/param->mass_W,2.);
	double C7pSM = param->mass_s/mass_b_muW*(-0.5*A0t(xt)-23./36.);
	Cpb[7]=pow(eta_mu,16./23.)*C7pSM;
	double C8pSM = param->mass_s/mass_b_muW*(-0.5*F0t(xt)-1./3.);
	Cpb[8]=pow(eta_mu,14./23.)*C8pSM;
	
	if((param->SM==1)||(param->THDM_model>0)||(param->mass_A02!=0.)||(param->mass_H03!=0.)) return;


 	double sw2=pow(sin(atan(param->gp/param->g2)),2.);
	double MU[4];
		
	MU[1]=param->mass_u;
	MU[2]=mass_c_muW;
	MU[3]=mass_top_muW;

	double epsfac;
	if(param->THDM_model>0) epsfac=1.;
	else epsfac=pow((1.+epsilon_b(param)*param->tan_beta),2.);
	
	double yt= pow(mass_top_muW/param->mass_H,2.);
	double yb= pow(mass_top_muW,2.)/param->mass_c/mass_b_muW;
	double z= pow(param->mass_H/param->mass_W,2.);

	double Gamma_UL[7][4],Gamma_UR[7][4],Gamma_NL[4][4],Gamma_NR[4][4];
	double Gamma_U[7][7], G_aimn[7][3][4][4];
	double X_UL[3][7][4],X_UR[3][7][4],X_NL[3][4][4],X_NR[3][4][4];
	double MD[4],ME[4],VCKM[4][4],Mch[3],MsqU[7],MsqD[7],Msn[4],mintmp;
	double kappa,ag,aY,cosb,sinb,st,ct,alphas_mg;
	double a0a,a0b,a0c,a0Q1,a0Q2,a1,Dp,Dm;
	int ae,be,ce,de,ee,fe,ge,je,ke,me,ne;
	
	VCKM[1][1]=param->Vud;
	VCKM[1][2]=param->Vus;
	VCKM[1][3]=-(param->Vts*param->Vtb+param->Vcs*param->Vcb)/param->Vus; /* Vub from unitarity */
	VCKM[2][1]=param->Vcd;
	VCKM[2][2]=param->Vcs;
	VCKM[2][3]=param->Vcb;
	VCKM[3][1]=param->Vtd;
	VCKM[3][2]=param->Vts;
	VCKM[3][3]=param->Vtb;
	
	sinb=sin(atan(param->tan_beta));
	cosb=cos(atan(param->tan_beta));
	ct=param->stop_mix[2][2];
	st=param->stop_mix[1][2];
	
	MD[1]=param->mass_u;
	MD[2]=param->mass_s;
	MD[3]=mass_b_muW;

	ME[1]=param->mass_e;
	ME[2]=param->mass_mu;
	ME[3]=param->mass_tau;

	Mch[1]=param->mass_cha1;
	Mch[2]=param->mass_cha2;
	
	MsqU[1]=param->mass_upl;
	MsqU[2]=param->mass_chl;
	MsqU[3]=param->mass_t1;
	MsqU[4]=param->mass_upr;
	MsqU[5]=param->mass_chr;
	MsqU[6]=param->mass_t2;
	
	Msn[1]=param->mass_nuel;
	Msn[2]=param->mass_numl;
	Msn[3]=param->mass_nutl;

	
	if(param->THDM_model==0)
	{	alphas_mg=alphas_running(param->mass_gluino,param->mass_top_pole,param->mass_b_pole,param);
		ag=1.-7./12./pi*alphas_mg;
		aY=1.+alphas_mg/4./pi;
		
		kappa=1./(param->g2*param->g2*param->Vtb*param->Vts);
					if(param->sU_mix[1][1]*param->sU_mix[1][2]*param->sU_mix[1][3]*param->sU_mix[1][4]*param->sU_mix[1][5]*param->sU_mix[1][6]!=0.)
		{
		
			for(je=1;je<=5;je++)
			{
				for(ie=je+1;ie<=6;ie++) if (MsqU[je]>MsqU[ie])
				{
					mintmp=MsqU[je];
					MsqU[je]=MsqU[ie];
					MsqU[ie]=mintmp;
				}
			}
		
		
			for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UL[ae][ie]=param->sU_mix[ae][ie];
			for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UR[ae][ie]=param->sU_mix[ae][ie+3];
		}
		else
		{
			Gamma_UL[1][1]=1.;
			Gamma_UL[2][1]=0.;
			Gamma_UL[3][1]=0.;
			Gamma_UL[4][1]=0.;
			Gamma_UL[5][1]=0.;
			Gamma_UL[6][1]=0.;
			Gamma_UL[1][2]=0.;
			Gamma_UL[2][2]=1.;
			Gamma_UL[3][2]=0.;
			Gamma_UL[4][2]=0.;
			Gamma_UL[5][2]=0.;
			Gamma_UL[6][2]=0.;
			Gamma_UL[1][3]=0.;
			Gamma_UL[2][3]=0.;
			Gamma_UL[3][3]=ct;
			Gamma_UL[4][3]=0.;
			Gamma_UL[5][3]=0.;
			Gamma_UL[6][3]=-st;
		
			Gamma_UR[1][1]=0.;
			Gamma_UR[2][1]=0.;
			Gamma_UR[3][1]=0.;
			Gamma_UR[4][1]=1.;
			Gamma_UR[5][1]=0.;
			Gamma_UR[6][1]=0.;
			Gamma_UR[1][2]=0.;
			Gamma_UR[2][2]=0.;
			Gamma_UR[3][2]=0.;
			Gamma_UR[4][2]=0.;
			Gamma_UR[5][2]=1.;
			Gamma_UR[6][2]=0.;
			Gamma_UR[1][3]=0.;
			Gamma_UR[2][3]=0.;
			Gamma_UR[3][3]=st;
			Gamma_UR[4][3]=0.;
			Gamma_UR[5][3]=0.;
			Gamma_UR[6][3]=ct;
		}

		for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++)
		{
			Gamma_U[ae][ie]=Gamma_UL[ae][ie];
			Gamma_U[ae][ie+3]=Gamma_UR[ae][ie];
		}
						
		for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++) if(ae==be) Gamma_NL[ae][be]=Gamma_NR[ae][be]=1.; else Gamma_NL[ae][be]=Gamma_NR[ae][be]=0.;
				
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
		{
			X_UL[ie][ae][be]=0.;
			for(ce=1;ce<=3;ce++) X_UL[ie][ae][be]+=-param->g2*(ag*param->charg_Vmix[ie][1]*Gamma_UL[ae][ce]
			-aY*param->charg_Vmix[ie][2]*Gamma_UR[ae][ce]*MU[ce]/(sqrt(2.)*param->mass_W*sinb))*VCKM[ce][be];
		
			X_UR[ie][ae][be]=0.;
			for(ce=1;ce<=3;ce++) X_UR[ie][ae][be]+=param->g2*aY*param->charg_Umix[ie][2]*Gamma_UL[ae][ce]*VCKM[ce][be]*MD[be]/(sqrt(2)*param->mass_W*cosb);
		}
	
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++)
		{
			X_NL[ie][ae][be]=-param->g2*param->charg_Vmix[ie][1]*Gamma_NL[ae][be];
			X_NR[ie][ae][be]=param->g2*param->charg_Umix[ie][2]*Gamma_NL[ae][be]*ME[be]/(sqrt(2.)*param->mass_W*cosb);
		}
		
		for(ae=1;ae<=6;ae++) for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++)
		{
			G_aimn[ae][ie][me][ne]=0.5/sqrt(2.)*(sqrt(2.)*param->mass_W*param->charg_Vmix[ie][1]*Gamma_UL[ae][ne]*ag-MU[ne]*param->charg_Vmix[ie][2]*Gamma_UR[ae][ne]*aY)*(VCKM[me][3]*VCKM[ne][2]/VCKM[3][3]/VCKM[3][2]);
		}

	}
	
	/* Wilson coefficients C7 and C8 prime */ 
	double ld=-param->tan_beta;
	double C7pH=param->mass_s*mass_b_muW/mass_top_muW/mass_top_muW*1./3.*ld*ld*F7_1(yt);
	double C8pH=param->mass_s*mass_b_muW/mass_top_muW/mass_top_muW*1./3.*ld*ld*F8_1(yt);

	double C7pcharg=0.;
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C7pcharg+=pow(param->mass_W/Mch[ie],2.)*(X_UR[ie][ae][2]*X_UR[ie][ae][3]*h10(pow(MsqU[ae]/Mch[ie],2.)) + Mch[ie]/mass_b_muW*X_UR[ie][ae][2]*X_UL[ie][ae][3]*h20(pow(MsqU[ae]/Mch[ie],2.)));		
	C7pcharg*=-0.5*kappa; 
		
	double C8pcharg=0.;
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C8pcharg+=pow(param->mass_W/Mch[ie],2.)*(X_UR[ie][ae][2]*X_UR[ie][ae][3]*h50(pow(MsqU[ae]/Mch[ie],2.)) + Mch[ie]/mass_b_muW*X_UR[ie][ae][2]*X_UL[ie][ae][3]*h60(pow(MsqU[ae]/Mch[ie],2.)));		
	C8pcharg*=-0.5*kappa; 

	Cpb[7]+=pow(eta_mu,16./23.)*(C7pH+C7pcharg);
	Cpb[8]+=pow(eta_mu,14./23.)*(C8pH+C8pcharg);
	/* printf("C7p=%.5e\n",Cpb[7]); */
	/* printf("C8p=%.5e\n",Cpb[8]); */
	
	
	/* Wilson coefficients C9 and C10 prime */ 	
	double C10pH = -mass_b_muW*param->mass_s*(param->tan_beta*param->tan_beta/8./param->mass_W/param->mass_W
	+pow(param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_H,2.))*f20(yt)/sw2;
	
	double C9pH =(4.*sw2-1.)*C10pH - param->mass_s*mass_b_muW/mass_top_muW/mass_top_muW*D9H0(yt,ld);
	
	
	double B10pc=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++) B10pc+=-X_UR[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(0.5*X_NR[ie][be][2]*X_NR[je][be][2]*f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)) +X_NL[ie][be][2]*X_NL[je][be][2]*fabs(Mch[je]/Mch[ie])*f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
	B10pc*=kappa*param->mass_W*param->mass_W/2./param->g2/param->g2;
		
	double C9pc=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) C9pc+=X_UR[je][ae][2]*X_UR[ie][ae][3]*(2.*fabs(Mch[je]/Mch[ie])*f30(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))*param->charg_Vmix[je][1]*param->charg_Vmix[ie][1] -f40(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))*param->charg_Umix[je][1]*param->charg_Umix[ie][1]);
		
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=3;ce++) C9pc+=-X_UR[ie][be][2]*X_UR[ie][ae][3]*f40(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))*Gamma_UR[be][ce]*Gamma_UR[ae][ce];
	C9pc*=-kappa/8.;
		
	double C10pcharg=(B10pc-C9pc)/sw2;
	
	double B9pc=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++) B9pc+=X_UR[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(0.5*X_NR[ie][be][2]*X_NR[je][be][2]*f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)) -X_NL[ie][be][2]*X_NL[je][be][2]*fabs(Mch[je]/Mch[ie])*f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
	B9pc*=kappa*param->mass_W*param->mass_W/2./param->g2/param->g2;
	
	double D9pc=0.;
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) D9pc+=pow(param->mass_W/Mch[ie],2.)*X_UR[ie][ae][2]*X_UR[ie][ae][3]*h30(pow(MsqU[ae]/Mch[ie],2.));
	D9pc*=kappa;

	double C9pcharg=(1.-4.*sw2)/sw2*C9pc-B9pc/sw2-D9pc;
	
	Cpb[9]+=C9pH+C9pcharg;
	Cpb[10]+=C10pH+C10pcharg;	
	/* printf("C9p=%.5e\n",Cpb[9]); */
	/* printf("C10p=%.5e\n",Cpb[10]); */


	/* Wilson coefficients CQ1 and CQ2 prime */ 
	double NQ1pH=-param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*xt*f30(xt,z);
	
	double BQ1pH=param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*f70(xt,z);
	
	double complex CQ1pH=(NQ1pH+BQ1pH)*param->mass_s/sw2;
	
	double complex CQ2pH=CQ1pH;
	
	
	double BQ1pc1=0.;
	double BQ1pc2=0.;
	double BQ2pc1=0.;
	double BQ2pc2=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
	{ 	BQ1pc1+=X_UR[je][ae][2]*X_UL[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NL[ie][be][2]*X_NR[je][be][2]*f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.))
	); 	BQ1pc2+=X_UR[je][ae][2]*X_UL[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NR[ie][be][2]*X_NL[je][be][2]*fabs(Mch[je]/Mch[ie])*f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
	}
	
	double BQ1pc=(BQ1pc1+BQ1pc2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;
	double BQ2pc=(BQ1pc1-BQ1pc2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;

	
	double NQ1pc=0.;
	double NQ2pc=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++)
		{ Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]+Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]-Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
		}
		a0a=-(fabs(Mch[ie]/Mch[je])*f30(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))*param->charg_Umix[ie][2]*param->charg_Vmix[je][1])*kron(ae,be);

a0b=-(f40(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))*param->charg_Umix[je][2]*param->charg_Vmix[ie][1])*kron(ae,be);
		a0c=1./param->mass_W*f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))*kron(ie,je);
		a0Q1=a0a+a0b+Dp*a0c;
		a0Q2=-a0a+a0b+Dm*a0c;
		a1=Mch[ie]/sqrt(2.)/param->mass_W*f80(pow(MsqU[ae]/Mch[ie],2.))*kron(ie,je)*kron(ae,be);
		
		NQ1pc+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q1+a1*param->tan_beta);
		NQ2pc+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q2+a1*param->tan_beta);

	}
	NQ1pc*=param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*param->mass_s/sw2;	NQ2pc*=param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*param->mass_s/sw2;
	
	double complex CQ1pcharg=NQ1pc+BQ1pc;
	CQpb[1]=CQ1pH+CQ1pcharg;
	CQpb[1]/=epsfac;
	

	double complex CQ2pcharg=NQ2pc+BQ2pc;
	CQpb[2]=CQ2pH+CQ2pcharg;
	CQpb[2]/=epsfac;

	int nf=5;
	double beta0 = 11.-2./3.*nf;
	CQpb[1]*=pow(eta_mu,-4./beta0);
	CQpb[2]*=pow(eta_mu,-4./beta0);
	/* printf("CQ1p=%.5e\n",CQpb[1]);
	printf("CQ2p=%.5e\n",CQpb[2]); */

	return;
}


/*-----------------------------------------------------------------------------*/

	void CQ_calculator(double complex CQ0b[], double complex CQ1b[], double mu_W, double mu, struct parameters* param)
{
	int ie;
	for(ie=1;ie<=2;ie++) CQ0b[ie]=CQ1b[ie]=0.;
	
	if(param->SM==1) return;

 	double sw2=pow(sin(atan(param->gp/param->g2)),2.);
	double MU[4];
	
	double mass_c_muW=running_mass(param->mass_c,param->mass_c,mu_W,param->mass_top_pole,param->mass_b_pole,param);
	double mass_b_muW=running_mass(param->mass_b,param->mass_b,mu_W,param->mass_top_pole,param->mass_b,param);
	double mass_top_muW=running_mass(param->mtmt,param->mtmt,mu_W,param->mass_top_pole,param->mass_b,param);
	
	MU[1]=param->mass_u;
	MU[2]=mass_c_muW;
	MU[3]=mass_top_muW;

	double alphas_muW=alphas_running(mu_W,param->mass_top_pole,param->mass_b_pole,param);
	double alphas_mu=alphas_running(mu,param->mass_top_pole,param->mass_b_pole,param);	
	double eta_mu=alphas_muW/alphas_mu;

	int nf=5;
	double beta0 = 11.-2./3.*nf;
	
	double complex CQ1H_0,CQ2H_0;

	if(param->THDM_model==0) 
	{	
		param->lambda_u[3][3]=1./param->tan_beta;
		param->lambda_d[2][2]=-param->tan_beta;
		param->lambda_d[3][3]=-param->tan_beta;
		param->lambda_l[2][2]=-param->tan_beta;
	}

	/* Wilson coefficients CQ1 et CQ2 in 2HDM */ 
 
	double complex CQ1box=-param->mass_mu*(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3])*param->lambda_l[2][2]*pow(1./2./param->mass_W,2.)/sw2*Bplus(param->mass_H*param->mass_H/param->mass_W/param->mass_W,MU[3]*MU[3]/param->mass_W/param->mass_W)/4.;

	double complex CQ2box=-CQ1box;
		
	double Pplus=-D2(param->mass_H*param->mass_H/param->mass_W/param->mass_W,MU[3]*MU[3]/param->mass_W/param->mass_W)*MU[3]*MU[3]/param->mass_W/param->mass_W;
	
	double complex CQ1peng1=-param->mass_mu*(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3])*param->lambda_l[2][2]/4./sw2*Pplus*(sin(param->alpha)*sin(param->alpha)/param->mass_h0/param->mass_h0+cos(param->alpha)*cos(param->alpha)/param->mass_H0/param->mass_H0)/4.;
	
	double complex CQ2peng1=param->mass_mu*(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3])*param->lambda_l[2][2]/4./sw2*Pplus/param->mass_A0/param->mass_A0/4.;
		
	double complex CQ1peng2=param->mass_mu*(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3])*param->lambda_l[2][2]/4./sw2*Pplus*(sin(param->alpha)*sin(param->alpha)/param->mass_h0/param->mass_h0*(param->mass_H*param->mass_H-param->mass_h0*param->mass_h0)/param->mass_W/param->mass_W+cos(param->alpha)*cos(param->alpha)/param->mass_H0/param->mass_H0*(param->mass_H*param->mass_H-param->mass_H0*param->mass_H0)/param->mass_W/param->mass_W)/4.;
	
	double complex CQ2peng2=-param->mass_mu*(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3])*param->lambda_l[2][2]/4./sw2*Pplus/param->mass_A0/param->mass_A0*(param->mass_H*param->mass_H-param->mass_A0*param->mass_A0)/param->mass_W/param->mass_W/4.;
 
	
	double complex CQ1self=-param->mass_mu*param->lambda_d[3][3]*param->lambda_l[2][2]/4./sw2*(param->mass_H*param->mass_H/param->mass_W/param->mass_W+((param->lambda_d[3][3]+param->mass_s/mass_b_muW*param->lambda_d[2][2])*param->lambda_u[3][3]))*Pplus*(sin(param->alpha)*sin(param->alpha)/param->mass_h0/param->mass_h0+cos(param->alpha)*cos(param->alpha)/param->mass_H0/param->mass_H0)/4.;
	
	double complex CQ2self=param->mass_mu*param->lambda_d[3][3]*param->lambda_l[2][2]/4./sw2*(param->mass_H*param->mass_H/param->mass_W/param->mass_W+((param->lambda_d[3][3]-param->mass_s/mass_b_muW*param->lambda_d[2][2])*param->lambda_u[3][3]))*Pplus/param->mass_A0/param->mass_A0/4.;
		
	
	CQ1box *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/param->mass_s*param->lambda_u[3][3]);
	
	CQ2box *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3]);
	
	CQ1peng1 *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3]);
	
	CQ2peng1 *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3]);

	CQ1peng2 *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3]);
	
	CQ2peng2 *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3]);
	
	CQ1self *= 1.+param->mass_s/mass_b_muW*param->lambda_d[2][2]/param->lambda_d[3][3];

	CQ2self *= 1.+param->mass_s/mass_b_muW*param->lambda_d[2][2]/param->lambda_d[3][3];
	
	CQ1H_0=(CQ1box+CQ1peng1+CQ1peng2+CQ1self);
	CQ2H_0=(CQ2box+CQ2peng1+CQ2peng2+CQ2self);
		
	if(param->THDM_model>0)
	{
		CQ0b[1]=(CQ1H_0)*mass_b_muW/sw2;
		CQ0b[2]=(CQ2H_0)*mass_b_muW/sw2;
	
		CQ0b[1]*=pow(eta_mu,-4./beta0);
		CQ0b[2]*=pow(eta_mu,-4./beta0);

		return;
	}
	
	double epsfac=pow((1.+epsilon_b(param)*param->tan_beta),2.);
	
	double xt= pow(mass_top_muW/param->mass_W,2.);
	double yt= pow(mass_top_muW/param->mass_H,2.);
	double yb= pow(mass_top_muW,2.)/param->mass_c/mass_b_muW;
	double z= pow(param->mass_H/param->mass_W,2.);

	double Gamma_UL[7][4],Gamma_UR[7][4],Gamma_NL[4][4],Gamma_NR[4][4];
	double Gamma_U[7][7], G_aimn[7][3][4][4],I_LR[7][7],P_U[7][7];
	double X_UL[3][7][4],X_UR[3][7][4],X_NL[3][4][4],X_NR[3][4][4];
	double MD[4],ME[4],VCKM[4][4],Mch[3],MsqU[7],MsqD[7],Msn[4],mintmp;
	double kappa,ag,aY,cosb,sinb,st,ct,alphas_mg,temp;
	double a0,a0Q1,a0Q2,a1,a0a,a0b,a0c,a0p,a2p,Dp,Dm,Dpm;
	int ae,be,ce,de,ee,fe,ge,je,ke,me,ne;

		
	VCKM[1][1]=param->Vud;
	VCKM[1][2]=param->Vus;
	VCKM[1][3]=-(param->Vts*param->Vtb+param->Vcs*param->Vcb)/param->Vus; /* Vub from unitarity */
	VCKM[2][1]=param->Vcd;
	VCKM[2][2]=param->Vcs;
	VCKM[2][3]=param->Vcb;
	VCKM[3][1]=param->Vtd;
	VCKM[3][2]=param->Vts;
	VCKM[3][3]=param->Vtb;
		
	sinb=sin(atan(param->tan_beta));
	cosb=cos(atan(param->tan_beta));
	ct=param->stop_mix[2][2];
	st=param->stop_mix[1][2];
		
	MD[1]=param->mass_u;
	MD[2]=param->mass_s;
	MD[3]=mass_b_muW;

	ME[1]=param->mass_e;
	ME[2]=param->mass_mu;
	ME[3]=param->mass_tau;

	Mch[1]=param->mass_cha1;
	Mch[2]=param->mass_cha2;
		
	MsqU[1]=param->mass_upl;
	MsqU[2]=param->mass_chl;
	MsqU[3]=param->mass_t1;
	MsqU[4]=param->mass_upr;
	MsqU[5]=param->mass_chr;
	MsqU[6]=param->mass_t2;
		
	Msn[1]=param->mass_nuel;
	Msn[2]=param->mass_numl;
	Msn[3]=param->mass_nutl;

	if((param->mass_A02==0.)&&(param->mass_H03==0.))
	{
		alphas_mg=alphas_running(param->mass_gluino,param->mass_top_pole,param->mass_b_pole,param);
		ag=1.-7./12./pi*alphas_mg;
		aY=1.+alphas_mg/4./pi;
	
		kappa=1./(param->g2*param->g2*param->Vtb*param->Vts);
		
	if(param->sU_mix[1][1]*param->sU_mix[1][2]*param->sU_mix[1][3]*param->sU_mix[1][4]*param->sU_mix[1][5]*param->sU_mix[1][6]!=0.)
		{
		
			for(je=1;je<=5;je++)
			{
				for(ie=je+1;ie<=6;ie++) if (MsqU[je]>MsqU[ie])
				{
					mintmp=MsqU[je];
					MsqU[je]=MsqU[ie];
					MsqU[ie]=mintmp;
				}
			}
		
		
			for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UL[ae][ie]=param->sU_mix[ae][ie];
			for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UR[ae][ie]=param->sU_mix[ae][ie+3];
		}
		else
		{
			Gamma_UL[1][1]=1.;
			Gamma_UL[2][1]=0.;
			Gamma_UL[3][1]=0.;
			Gamma_UL[4][1]=0.;
			Gamma_UL[5][1]=0.;
			Gamma_UL[6][1]=0.;
			Gamma_UL[1][2]=0.;
			Gamma_UL[2][2]=1.;
			Gamma_UL[3][2]=0.;
			Gamma_UL[4][2]=0.;
			Gamma_UL[5][2]=0.;
			Gamma_UL[6][2]=0.;
			Gamma_UL[1][3]=0.;
			Gamma_UL[2][3]=0.;
			Gamma_UL[3][3]=ct;
			Gamma_UL[4][3]=0.;
			Gamma_UL[5][3]=0.;
			Gamma_UL[6][3]=-st;
		
			Gamma_UR[1][1]=0.;
			Gamma_UR[2][1]=0.;
			Gamma_UR[3][1]=0.;
			Gamma_UR[4][1]=1.;
			Gamma_UR[5][1]=0.;
			Gamma_UR[6][1]=0.;
			Gamma_UR[1][2]=0.;
			Gamma_UR[2][2]=0.;
			Gamma_UR[3][2]=0.;
			Gamma_UR[4][2]=0.;
			Gamma_UR[5][2]=1.;
			Gamma_UR[6][2]=0.;
			Gamma_UR[1][3]=0.;
			Gamma_UR[2][3]=0.;
			Gamma_UR[3][3]=st;
			Gamma_UR[4][3]=0.;
			Gamma_UR[5][3]=0.;
			Gamma_UR[6][3]=ct;
		}

		for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++)
		{
			Gamma_U[ae][ie]=Gamma_UL[ae][ie];
			Gamma_U[ae][ie+3]=Gamma_UR[ae][ie];
		}
		
		for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) I_LR[ae][be]=0.;
		for(ae=1;ae<=3;ae++) I_LR[ae][ae]=1.;
		for(ae=4;ae<=6;ae++) I_LR[ae][ae]=-1.;
		
		for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++) for(de=1;de<=6;de++) P_U[ae][be]=Gamma_U[ae][ce]*I_LR[ce][de]*Gamma_U[be][de];
		
		for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++) if(ae==be) Gamma_NL[ae][be]=Gamma_NR[ae][be]=1.; else Gamma_NL[ae][be]=Gamma_NR[ae][be]=0.;
				
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
		{
			X_UL[ie][ae][be]=0.;
			for(ce=1;ce<=3;ce++) X_UL[ie][ae][be]+=-param->g2*(ag*param->charg_Vmix[ie][1]*Gamma_UL[ae][ce]-aY*param->charg_Vmix[ie][2]*Gamma_UR[ae][ce]*MU[ce]/(sqrt(2.)*param->mass_W*sinb))*VCKM[ce][be];
		
			X_UR[ie][ae][be]=0.;
			for(ce=1;ce<=3;ce++) X_UR[ie][ae][be]+=param->g2*aY*param->charg_Umix[ie][2]*Gamma_UL[ae][ce]*VCKM[ce][be]*MD[be]/(sqrt(2)*param->mass_W*cosb);
		}
	
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++)
		{
			X_NL[ie][ae][be]=-param->g2*param->charg_Vmix[ie][1]*Gamma_NL[ae][be];
			X_NR[ie][ae][be]=param->g2*param->charg_Umix[ie][2]*Gamma_NL[ae][be]*ME[be]/(sqrt(2.)*param->mass_W*cosb);
		}
		
		for(ae=1;ae<=6;ae++) for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++)
		{
			G_aimn[ae][ie][me][ne]=0.5/sqrt(2.)*(sqrt(2.)*param->mass_W*param->charg_Vmix[ie][1]*Gamma_UL[ae][ne]*ag-MU[ne]*param->charg_Vmix[ie][2]*Gamma_UR[ae][ne]*aY)*(VCKM[me][3]*VCKM[ne][2]/VCKM[3][3]/VCKM[3][2]);
		}
	
	/* LO */
	
	double NQ10H=-param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*xt*f30(xt,z);
	double BQ10H=param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*f70(xt,z);
	CQ1H_0=(NQ10H+BQ10H)*mass_b_muW/sw2; 
	CQ2H_0=-CQ1H_0;
	
	
	double BQ10c1=0.;
	double BQ10c2=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
	{ BQ10c1+=X_UL[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NR[ie][be][2]*X_NL[je][be][2]*f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
BQ10c2+=X_UL[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NL[ie][be][2]*X_NR[je][be][2]*fabs(Mch[je]/Mch[ie])*f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
	}
	double BQ10c=(BQ10c1+BQ10c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;
	double BQ20c=-(BQ10c1-BQ10c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;
	
	
	double NQ10c=0.;
	double NQ20c=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++) 
		{
			Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]+Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
			Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]-Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
		}
			a0a=-(fabs(Mch[ie]/Mch[je])*f30(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))*param->charg_Umix[ie][2]*param->charg_Vmix[je][1])*kron(ae,be);
a0b=-(f40(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))*param->charg_Umix[je][2]*param->charg_Vmix[ie][1])*kron(ae,be);
		a0c=1./param->mass_W*f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))*kron(ie,je);
		a0Q1=a0a+a0b+Dp*a0c;
		a0Q2=-a0a+a0b+Dm*a0c;
		
		a1=Mch[ie]/sqrt(2.)/param->mass_W*f80(pow(MsqU[ae]/Mch[ie],2.))*kron(ie,je)*kron(ae,be);
		NQ10c+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q1+a1*param->tan_beta);
		NQ20c+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q2+a1*param->tan_beta);
	}
	NQ10c*=param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;
        NQ20c*=-param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;


	double complex CQ1charg_0=NQ10c+BQ10c;
	double complex CQ2charg_0=NQ20c+BQ20c;
	
	CQ0b[1]=CQ1H_0+CQ1charg_0;
	CQ0b[1]/=epsfac;
		
	CQ0b[2]=CQ2H_0+CQ2charg_0;
	CQ0b[2]/=epsfac;
	
	/* NLO - Charged Higgs */
	
	double NQ11H=-param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*(f141(xt,z)+8.*xt*(f30(xt,z)+xt*(f30(xt*1.0001,z)-f30(xt*0.9999,z))/0.0002)*log(mu_W*mu_W/mass_top_muW/mass_top_muW));
	double BQ11H=param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*(f111(xt,z)+8.*(f70(xt*1.0001,z)-f70(xt*0.9999,z))/0.0002*log(mu_W*mu_W/mass_top_muW/mass_top_muW));
	double CQ1H_1=(NQ11H+BQ11H)*mass_b_muW/sw2;
	double CQ2H_1=-CQ1H_1;

	/* NLO - charginos */
	
	double BQ11c1=0.;
	double BQ11c2=0.;
	
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++) BQ11c1+=X_UL[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NR[ie][be][2]*X_NL[je][be][2]*(f121(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.))+4.*(f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*1.0001,pow(Msn[be]/Mch[ie],2.))-f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*0.9999,pow(Msn[be]/Mch[ie],2.)))/0.0002*log(pow(mass_top_muW/MsqU[ae],2.))));
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++) BQ11c2+=X_UL[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NL[ie][be][2]*X_NR[je][be][2]*fabs(Mch[je]/Mch[ie])*(f131(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.))+4.*(f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*1.0001,pow(Msn[be]/Mch[ie],2.))-f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*0.9999,pow(Msn[be]/Mch[ie],2.)))/0.0002*log(pow(mass_top_muW/MsqU[ae],2.))));
	
	double BQ11c=(BQ11c1+BQ11c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;
	double BQ21c=-(BQ11c1-BQ11c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;
	
	
	double NQ11c=0.;
	double NQ21c=0.;

	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++)
		{ 	Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]+Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
	Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]-Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
		}
			a0a=-(fabs(Mch[ie]/Mch[je])*(f181(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))+4.*(f30(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.)*1.0001)-f30(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.)*0.9999))/0.0002*log(pow(mass_top_muW/MsqU[ae],2.)))*param->charg_Umix[ie][2]*param->charg_Vmix[je][1])*kron(ae,be);
					a0b=-((f191(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))+4.*(f40(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.)*1.0001)-f40(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.)*0.9999))/0.0002*log(pow(mass_top_muW/MsqU[ae],2.)))*param->charg_Umix[je][2]*param->charg_Vmix[ie][1])*kron(ae,be);
					a0c=1./param->mass_W*(f171(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))+4.*(f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))+(f30(pow(MsqU[ae]/Mch[ie],2.)*1.0001,pow(MsqU[be]/Mch[ie],2.))-f30(pow(MsqU[ae]/Mch[ie],2.)*0.9999,pow(MsqU[be]/Mch[ie],2.)))/0.0002+(f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.)*1.0001)-f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.)*0.9999))/0.0002)*log(pow(mass_top_muW/MsqU[ae],2.)))*kron(ie,je);
	
		a0Q1=a0a+a0b+Dp*a0c;
		a0Q2=-a0a+a0b+Dm*a0c;
			a0p=4.*G_aimn[ae][ie][me][ne]/param->mass_W/(VCKM[me][3]*VCKM[ne][2]/VCKM[3][3]/VCKM[3][2])/param->charg_Umix[je][2]*f151(pow(MsqU[ae]/Mch[ie],2.))*kron(ie,je)*kron(ae,be)*kron(me,ne);
				a1=Mch[ie]/sqrt(2.)/param->mass_W*(f161(pow(MsqU[ae]/Mch[ie],2.))+4.*(f80(pow(MsqU[ae]/Mch[ie],2.)*1.0001)-f80(pow(MsqU[ae]/Mch[ie],2.)*0.9999))/0.0002*log(pow(mass_top_muW/MsqU[ae],2.)))*kron(ie,je)*kron(ae,be);
		a2p=Gamma_UL[be][me]*(VCKM[me][3]*VCKM[ne][2]/VCKM[3][3]/VCKM[3][2])*param->charg_Umix[je][2]/2./param->mass_W*f151(pow(MsqU[ae]/Mch[ie],2.))*kron(ie,je)*kron(ae,be)*kron(me,ne);
		
		NQ11c+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q1+a1*param->tan_beta)
		+G_aimn[ae][ie][me][ne]*param->charg_Umix[je][2]*a0p
		+Gamma_UL[be][me]*param->charg_Umix[je][2]*a2p*pow(param->mass_s*param->tan_beta,2.);	
		NQ21c+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q2+a1*param->tan_beta)
		+G_aimn[ae][ie][me][ne]*param->charg_Umix[je][2]*a0p
		+Gamma_UL[be][me]*param->charg_Umix[je][2]*a2p*pow(param->mass_s*param->tan_beta,2.);
	}
		NQ11c*=param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;
	NQ21c*=-param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;
	
	double complex CQ1charg_1=NQ11c+BQ11c;
	
	if(cabs(CQ1charg_1)*alphas_mu/4./pi>cabs(CQ0b[1])) CQ1charg_1*=cabs(CQ0b[1])/cabs(CQ1charg_1)*4.*pi/alphas_mu;
	if(cabs(CQ1H_1)*alphas_mu/4./pi>cabs(CQ0b[1])) CQ1H_1*=cabs(CQ0b[1])/cabs(CQ1H_1)*4.*pi/alphas_mu;
	CQ1b[1]=CQ1H_1+CQ1charg_1;
	
	CQ1b[1]/=epsfac;
	
	double complex CQ2charg_1=NQ21c+BQ21c;
	
	if(cabs(CQ2charg_1)*alphas_mu/4./pi>cabs(CQ0b[2])) CQ2charg_1*=cabs(CQ0b[2])/cabs(CQ2charg_1)*4.*pi/alphas_mu;
	if(cabs(CQ2H_1)*alphas_mu/4./pi>cabs(CQ0b[2])) CQ2H_1*=cabs(CQ0b[2])/cabs(CQ2H_1)*4.*pi/alphas_mu;
	CQ1b[2]=CQ2H_1+CQ2charg_1;
		
	CQ1b[2]/=epsfac;
	
	
	/* Wilson coefficient CQ1 */ 
	/* NLO  - four points */
	double BQ11f1=0.;
	double BQ11f2=0.;
	
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(fe=1;fe<=3;fe++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++)
	{ BQ11f1+=-X_UL[je][be][2]*X_UR[ie][ae][3]*pow(param->mass_W/Mch[ie],2.)*P_U[ae][ce]*MsqU[ce]/Mch[ie]*P_U[ce][be]*(1.+log(pow(mass_top_muW/MsqU[ce],2.)))	*(f90(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.),pow(Msn[fe]/Mch[ie],2.))*X_NR[ie][fe][2]*X_NL[je][fe][2]);
		BQ11f2+=-X_UL[je][be][2]*X_UR[ie][ae][3]*pow(param->mass_W/Mch[ie],2.)*P_U[ae][ce]*MsqU[ce]/Mch[ie]*P_U[ce][be]*(1.+log(pow(mass_top_muW/MsqU[ce],2.)))	*(fabs(Mch[je]/Mch[ie])*f100(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.),pow(Msn[fe]/Mch[ie],2.))*X_NL[ie][fe][2]*X_NR[je][fe][2]);
	}
		
	double BQ11f=(BQ11f1+BQ11f2)*2./3.*kappa/param->g2/param->g2/sw2;
	double BQ21f=-(BQ11f1-BQ11f2)*2./3.*kappa/param->g2/param->g2/sw2;
	
	double NQ11f=0.;
	double NQ21f=0.;

	for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) for(ae=1;ae<=6;ae++) for(de=1;de<=6;de++) for(ke=1;ke<=6;ke++) NQ11f+=G_aimn[ae][ie][me][ne]*Gamma_UL[de][me]*param->charg_Umix[ie][2]*P_U[ae][ke]*MsqU[ke]/Mch[ie]*P_U[ke][de]*(1.+log(pow(mass_top_muW/MsqU[ke],2.)))*param->tan_beta*Mch[ie]/sqrt(2.)*f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[de]/Mch[ie],2.));

	NQ21f=NQ11f;

	for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++) for(ke=1;ke<=6;ke++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++) 
		{		Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[be][fe]*Gamma_UL[ce][fe]+Gamma_UL[be][fe]*Gamma_UR[ce][fe]); 
Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[be][fe]*Gamma_UL[ce][fe]-Gamma_UL[be][fe]*Gamma_UR[ce][fe]); 
		}
		temp=G_aimn[ae][ie][me][ne]*Gamma_UL[ce][me]*param->charg_Umix[ie][2]*P_U[ae][ke]*MsqU[ke]/Mch[ie]*P_U[ke][be]*(1.+log(pow(mass_top_muW/MsqU[ke],2.)))*f60(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.));	
		NQ11f+=Dp*temp;
		NQ21f+=Dm*temp;
	}
	
	for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) for(ae=1;ae<=6;ae++) for(ce=1;ce<=6;ce++)  for(de=1;de<=6;de++) for(ke=1;ke<=6;ke++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++)
		{ Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[ce][fe]+Gamma_UL[ae][fe]*Gamma_UR[ce][fe]); Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[ce][fe]-Gamma_UL[ae][fe]*Gamma_UR[ce][fe]); 
		}
	
	temp=G_aimn[ae][ie][me][ne]*Gamma_UL[de][me]*param->charg_Umix[ie][2]*P_U[ce][ke]*MsqU[ke]/Mch[ie]*P_U[ke][de]*(1.+log(pow(mass_top_muW/MsqU[ke],2.)))*f60(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.),pow(MsqU[de]/Mch[ie],2.));	
	NQ11f+=temp*Dp;
	NQ21f+=temp*Dm;
	}
	
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) for(ae=1;ae<=6;ae++) for(de=1;de<=6;de++) for(ke=1;ke<=6;ke++)
	{
		temp=-G_aimn[ae][ie][me][ne]*Gamma_UL[de][me]*param->charg_Umix[je][2]*P_U[ae][ke]*MsqU[ke]/Mch[je]*P_U[ke][de]*(1.+log(pow(mass_top_muW/MsqU[ke],2.)))*param->mass_W*(fabs(Mch[ie]/Mch[je])*f60(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.),pow(MsqU[de]/Mch[je],2.))*param->charg_Umix[ie][2]*param->charg_Vmix[je][1]+f50(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.),pow(MsqU[de]/Mch[je],2.))*param->charg_Umix[je][2]*param->charg_Vmix[ie][1]); 
	NQ11f+=temp;
	NQ21f+=-temp;
	}

	
	for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) for(ae=1;ae<=6;ae++) for(de=1;de<=6;de++)for(ee=1;ee<=6;ee++) for(ge=1;ge<=6;ge++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++)
		{ Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ee][fe]*Gamma_UL[ge][fe]+Gamma_UL[ee][fe]*Gamma_UR[ge][fe]); 
Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ee][fe]*Gamma_UL[ge][fe]-Gamma_UL[ee][fe]*Gamma_UR[ge][fe]); 
		}
	
	temp=	-G_aimn[ae][ie][me][ne]*Gamma_UL[de][me]*param->charg_Umix[ie][2]*P_U[ae][ee]*(1.+log(pow(mass_top_muW/MsqU[ge],2.))-f110(pow(MsqU[ee]/Mch[ie],2.),pow(MsqU[ge]/Mch[ie],2.)))*P_U[ge][de]*f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[de]/Mch[ie],2.));	
	NQ11f+=Dp*temp;
	NQ21f+=Dm*temp;
	}
		NQ11f*=-4./3.*param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;

NQ21f*=4./3.*param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;


	double complex CQ1four_1=NQ11f+BQ11f;

	if(cabs(CQ1four_1)*alphas_mu/4./pi>cabs(CQ0b[1])) CQ1four_1*=cabs(CQ0b[1])/cabs(CQ1four_1)*4.*pi/alphas_mu;
		
	CQ1b[1]+=CQ1four_1;

	double complex CQ2four_1=NQ21f+BQ21f;
	
	if(cabs(CQ2four_1)*alphas_mu/4./pi>cabs(CQ0b[2])) CQ2four_1*=cabs(CQ0b[2])/cabs(CQ2four_1)*4.*pi/alphas_mu;
	CQ1b[2]+=CQ2four_1;

	CQ0b[1]*=pow(eta_mu,-4./beta0);
	CQ0b[2]*=pow(eta_mu,-4./beta0);
	CQ1b[1]*=pow(eta_mu,-4./beta0)*eta_mu;
	CQ1b[2]*=pow(eta_mu,-4./beta0)*eta_mu;
	}
	
/* NMSSM */

	if((param->mass_A02!=0.)||(param->mass_H03!=0.))
	{
		double s=param->lambdaSNMSSM/param->lambdaNMSSM;
		double v=sqrt(1./sqrt(2.)/param->Gfermi);
		
		double v_deltam_s=v/s*(sqrt(2.)*param->AlambdaNMSSM-2.*param->kappaNMSSM*s)/(sqrt(2.)*param->AlambdaNMSSM+param->kappaNMSSM*s);
		
		CQ0b[1]=0.;
		CQ0b[2]=0.;
		CQ1b[1]=0.;
		CQ1b[2]=0.;
		
		double mH0[4],mA0[3],mstop[3];
		
		mstop[0]=param->mass_upr;
		mstop[1]=param->mass_t1;
		mstop[2]=param->mass_t2;
		
		mH0[1]=param->mass_h0;
		mH0[2]=param->mass_H0;
		mH0[3]=param->mass_H03;
		mA0[1]=param->mass_A0;
		mA0[2]=param->mass_A02;
		
		double Ralj[3][3][3],Qalj[4][3][3],G1[4][4][3][3];
		double T1[3][4][4],T2[4][4][4];
		
		double TU[4][4];
		for(ie=1;ie<=3;ie++) TU[ie][1]=TU[1][ie]=0.;
		TU[1][1]=1.;
		for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) TU[ie+1][je+1]=param->stop_mix[ie][je];

		
		double vu=sqrt(pow(sin(atan(param->tan_beta)),2.)/sqrt(2.)/param->Gfermi);
		double vd=vu/param->tan_beta;
		int le;

		for(ae=1;ae<=2;ae++) for(le=1;le<=2;le++) for(je=1;je<=2;je++) Ralj[ae][le][je]=-param->g2/sqrt(2.)*(param->A0_mix[ae][1]*param->charg_Umix[2][le]*param->charg_Vmix[2][je]+param->A0_mix[ae][2]*param->charg_Umix[1][le]*param->charg_Vmix[2][je])-param->lambdaNMSSM/sqrt(2.)*param->A0_mix[ae][3]*param->charg_Umix[2][le]*param->charg_Vmix[2][je];
		
		for(ae=1;ae<=3;ae++) for(le=1;le<=2;le++) for(je=1;je<=2;je++) Qalj[ae][le][je]=param->g2/sqrt(2.)*(param->H0_mix[ae][1]*param->charg_Umix[2][le]*param->charg_Vmix[2][je]+param->H0_mix[ae][2]*param->charg_Umix[1][le]*param->charg_Vmix[2][je])-param->lambdaNMSSM/sqrt(2.)*param->H0_mix[ae][3]*param->charg_Umix[2][le]*param->charg_Vmix[2][je];
		
		for(ie=1;ie<=3;ie++) for(ke=1;ke<=3;ke++) for(je=1;je<=2;je++) for(le=1;le<=2;le++) G1[ie][ke][je][le]=(TU[ie][2]*TU[ke][2]-kron(ie,1)*kron(ke,1))*param->charg_Vmix[1][le]*param->charg_Umix[2][je]-mass_top_muW/sqrt(2.)/sin(atan(param->tan_beta))/param->mass_W*TU[ie][3]*TU[ke][2]*param->charg_Vmix[2][le]*param->charg_Umix[2][je];
		
		for(ae=1;ae<=2;ae++) for(ie=1;ie<=3;ie++) for(ke=1;ke<=3;ke++) T1[ae][ie][ke]=(TU[ie][3]*TU[ke][2]-TU[ie][2]*TU[ke][3])*((param->lambdaNMSSM/sqrt(2.)*(vd*param->A0_mix[ae][3]+s*param->A0_mix[ae][1]))-param->A_u*param->A0_mix[ae][2]);
		
		for(ae=1;ae<=3;ae++) for(ie=1;ie<=3;ie++) for(ke=1;ke<=3;ke++) T2[ae][ie][ke]=-mass_top_muW/2./param->mass_W*(2.*mass_top_muW*param->H0_mix[ae][2]*(TU[ie][2]*TU[ke][2]+TU[ie][3]*TU[ke][3])	+((param->lambdaNMSSM/sqrt(2.)*(vd*param->H0_mix[ae][3]+s*param->H0_mix[ae][1]))+param->A_u*param->H0_mix[ae][2])*(TU[ie][3]*TU[ke][2]+TU[ie][2]*TU[ke][3]))	+param->mass_Z/2./sqrt(1.-sw2)*(1.-4./3.*sw2)*param->H0_mix[ae][2]*(TU[ie][1]*TU[ke][1]+TU[ie][2]*TU[ke][2])
	+2./3.*param->mass_W*sw2/(1.-sw2)*param->H0_mix[ae][2]*TU[ie][3]*TU[ke][3];
			
		double complex CQ1H=0.;
		for(ae=1;ae<=3;ae++) CQ1H+=(param->mass_H*param->mass_H/param->mass_W/param->mass_W*param->H0_mix[ae][1]*param->H0_mix[ae][1]*f30(param->mass_H*param->mass_H/mass_top_muW/mass_top_muW,param->mass_W*param->mass_W/mass_top_muW/mass_top_muW)	+mass_top_muW*mass_top_muW*mH0[ae]*mH0[ae]/param->mass_W/param->mass_W/param->mass_H/param->mass_H*f30(mass_top_muW*mass_top_muW/param->mass_H/param->mass_H,mass_top_muW*mass_top_muW/param->mass_W/param->mass_W))
		/mH0[ae]/mH0[ae];
		
		CQ1H*=-param->mass_mu/4.*param->tan_beta*param->tan_beta;
		
		double complex CQ2H=0.;
		for(ae=1;ae<=2;ae++) CQ2H+=((param->mass_H*param->mass_H/param->mass_W/param->mass_W*param->A0_mix[ae][1]*param->A0_mix[ae][1]+kron(ae,2)*param->A0_mix[ae][1])*f30(param->mass_H*param->mass_H/mass_top_muW/mass_top_muW,param->mass_W*param->mass_W/mass_top_muW/mass_top_muW)	+mass_top_muW*mass_top_muW*mA0[ae]*mA0[ae]/param->mass_W/param->mass_W/param->mass_H/param->mass_H*f30(mass_top_muW*mass_top_muW/param->mass_H/param->mass_H,mass_top_muW*mass_top_muW/param->mass_W/param->mass_W))/mA0[ae]/mA0[ae];
		
		CQ2H*=param->mass_mu/4.*param->tan_beta*param->tan_beta;
		
		double complex CAH=-I*param->lambdaNMSSM*param->AlambdaNMSSM/param->g2/param->mass_W*param->tan_beta*f30(param->mass_H*param->mass_H/mass_top_muW/mass_top_muW,param->mass_W*param->mass_W/mass_top_muW/mass_top_muW);
		
		double complex CQ1c=0.;
		
		for(ae=1;ae<=3;ae++) for(ie=1;ie<=3;ie++) for(ke=1;ke<=3;ke++) for(je=1;je<=2;je++) for(le=1;le<=2;le++) CQ1c+=G1[ie][ke][je][le]/mH0[ae]/mH0[ae]*( sqrt(2.)*param->H0_mix[ae][1]*param->H0_mix[ae][1]*Mch[je]/param->mass_W/cos(atan(param->tan_beta))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/Mch[je],2.))
			-2.*sqrt(2.)*param->H0_mix[ae][1]/param->g2*kron(ie,ke)*(Qalj[ae][le][je]*f40(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.))+Mch[je]/Mch[le]*Qalj[ae][je][le]*f30(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.)))		+2.*sqrt(2.)*param->H0_mix[ae][1]*T2[ae][ie][ke]*Mch[je]/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow(Mch[je]/mstop[ke-1],2.))
			+mH0[ae]*mH0[ae]/Mch[je]/Mch[je]*kron(ie,ke)*(param->charg_Umix[2][je]*param->charg_Vmix[1][le]*f50(pow(mstop[ie-1]/Mch[je],2.),pow(Mch[le]/Mch[je],2.),pow(param->mass_nutl/Mch[le],2.))
		
		-Mch[le]/Mch[je]*param->charg_Umix[2][le]*param->charg_Vmix[1][je]* f60(pow(mstop[ie-1]/Mch[je],2.),pow(Mch[le]/Mch[je],2.),pow(param->mass_nutl/Mch[le],2.))));
			
		CQ1c*=param->mass_mu/4.*param->tan_beta*param->tan_beta;
		
		double complex CQ2c=0.;
		
		for(ae=1;ae<=2;ae++) for(ie=1;ie<=3;ie++) for(ke=1;ke<=3;ke++) for(je=1;je<=2;je++) for(le=1;le<=2;le++) CQ2c+=G1[ie][ke][je][le]/mA0[ae]/mA0[ae]*(sqrt(2.)*param->A0_mix[ae][1]*param->A0_mix[ae][1]*Mch[je]/param->mass_W/cos(atan(param->tan_beta))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/Mch[je],2.))
			-2.*sqrt(2.)*param->A0_mix[ae][1]/param->g2*kron(ie,ke)*(-Ralj[ae][le][je]*f40(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.))+Mch[je]/Mch[le]*Ralj[ae][je][le]*f30(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.)))			-sqrt(2.)*param->A0_mix[ae][1]*T1[ae][ie][ke]*mass_top_muW*Mch[je]/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow(Mch[je]/mstop[ke-1],2.))
				+mA0[ae]*mA0[ae]/Mch[je]/Mch[je]*kron(ie,ke)*(param->charg_Umix[2][je]*param->charg_Vmix[1][le]*f50(pow(mstop[ie-1]/Mch[je],2.),pow(Mch[le]/Mch[je],2.),pow(param->mass_nutl/Mch[le],2.))
			-Mch[le]/Mch[je]*param->charg_Umix[2][le]*param->charg_Vmix[1][je]*f60(pow(mstop[ie-1]/Mch[je],2.),pow(Mch[le]/Mch[je],2.),pow(param->mass_nutl/Mch[le],2.))));
			
		CQ2c*=-param->mass_mu/4.*param->tan_beta*param->tan_beta;		
				
		double complex CAc=0.;
		for(ie=1;ie<=3;ie++) for(je=1;je<=2;je++) for(le=1;le<=2;le++) CAc+=I*param->tan_beta/sqrt(2.)*G1[ie][ie][je][le]*(v_deltam_s*kron(le,je)*fabs(Mch[je]/param->mass_W)*f80(pow(mstop[ie-1]/Mch[je],2.))-(Ralj[1][je][le]*fabs(Mch[je]/Mch[le])*f30(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.))-Ralj[1][le][je]*f40(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.))));
				
		CQ0b[1]=(CQ1H+CQ1c)*mass_b_muW/sw2/epsfac;
		CQ0b[2]=(CQ2H+CQ2c)*mass_b_muW/sw2/epsfac;

		double complex CA=CAH+CAc;

		if(param->mass_A0>mu_W) CQ0b[2]+=-v_deltam_s/2.*mass_b_muW/sw2*param->mass_mu*CA/param->mass_A0/param->mass_A0;

		CQ0b[1]*=pow(eta_mu,-4./beta0);
		CQ0b[2]*=pow(eta_mu,-4./beta0);
		CQ1b[1]*=pow(eta_mu,-4./beta0)*eta_mu;
		CQ1b[2]*=pow(eta_mu,-4./beta0)*eta_mu;
		
		if((param->mass_A0>param->mass_b_pole)&&(param->mass_A0<mu_W))
		{	
			double alphas_Ma1=alphas_running(param->mass_A0,param->mass_top_pole,param->mass_b_pole,param);	
			double eta_a1=alphas_Ma1/alphas_mu;
			double mass_b_ma1=running_mass(param->mass_b,param->mass_b,param->mass_A0,param->mass_top_pole,param->mass_b,param);
			CQ0b[2]+=-v_deltam_s/2.*mass_b_ma1/sw2*param->mass_mu*CA/param->mass_A0/param->mass_A0*pow(eta_a1,-4./beta0);
		}
		
		if(param->mass_A0<param->mass_b_pole)
		{	
			double width_A0=1.e-6;	
CQ0b[2]+=v_deltam_s/2.*param->mass_b/sw2*param->mass_mu*CA/(param->m_Bs*param->m_Bs-param->mass_A0*param->mass_A0+I*param->mass_A0*width_A0);
		}
	}
			
	return;
}
