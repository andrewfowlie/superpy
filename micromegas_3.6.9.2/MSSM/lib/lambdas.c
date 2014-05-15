/* lambdas.F -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

//#include "f2c.h"
#include <math.h>
#include <stdlib.h>

/* Common Block Declarations */

    extern double MbRun(double);
    extern double MtRun(double);
    extern double alphaQCD(double);
    extern double initQCD(double, double, double, double);
    extern double slhaVal(char*, double, int, ...);
    extern double slhaRead(char*fname,int mode); 
#define max(x,y)  (x)>(y) ? (x):(y)

struct {
    double gf, alph, amtau, ammuon, amz, amw;
} param_hdec__;

#define param_hdec__1 param_hdec__

struct {
    double lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7;
} hself_hdec__;

#define hself_hdec__1 hself_hdec__

/* Table of constant values */

static int c__1 = 1;
static int c__2 = 2;
static int c__6 = 6;
static int c__36 = 36;
static int c_b17 = 1000021;
static int c_b19 = 1000037;
static int c__3 = 3;
static double c_b24 = 1.4;
static int c__5 = 5;
static int c_b28 = 1000006;
static int c_b30 = 2000006;
static int c__43 = 43;
static int c__46 = 46;
static int c__49 = 49;
static double c_b83 = 4.;
static double c_b90 = .25;
static double c_b98 = .5;

#include <stdio.h>

double calcLambdas(void)
{  
   FILE*f;
    /* System generated locals */
    double d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static double v, q0;
    extern /* Subroutine */ int subh1_hdec__(double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *);
    static double ca, ad, cb, cf, ma, qb, sa, sb, au, pi, ft, xb, mu, qt, 
	    xt, ad0, au0, dmb, mch, bet, rmb, amu, mhp, err, hmp, rmt, amh12, 
	    amh22, albp, qcdl, amdr, tanb, alsp, altp, amsq, mglu, amur, zero,
	     mtop, amdl0, mchi0, susy, amdr0, amul0, amur0;
    extern int iargc_(void);
    static double tanba;
    
/*<       GF=         1.16639E-05 >*/
    param_hdec__1.gf = 1.16639e-5f;
/*<       ALPH=       0.00729735254 >*/
    param_hdec__1.alph = .00729735254f;
/*<       AMTAU=      1.777 >*/
    param_hdec__1.amtau = 1.777f;
/*<       AMMUON=     0.105658367 >*/
    param_hdec__1.ammuon = .105658367f;
/*<       AMZ=        91.1884003 >*/
    param_hdec__1.amz = 91.1884003f;
/*<       AMW=        79.9467442 >*/
    param_hdec__1.amw = 79.9467442f;
/*<       PI=4*DATAN(1D0) >*/
    pi = atan(1.) * 4;
/*<       V=1.D0/DSQRT(DSQRT(2.D0)*GF) >*/
    v = 1. / sqrt(sqrt(2.) * param_hdec__1.gf);
/*<       zero=0 >*/
    zero = 0.;
/*<       MTOP=     slhaVal1('SMINPUTS',zero,6) >*/
    mtop = slhaVal("SMINPUTS", zero, 1,c__6);
/*      MA=       slhaVal1('EXTPAR',zero,26) */
/*<       MA=       slhaVal1('MASS',zero,36) >*/
    ma = slhaVal("MASS", zero, 1,c__36);
/*<       MGLU=     slhaVal1('MASS',zero,1000021)   >*/
    mglu = slhaVal("MASS", zero,1, c_b17);
/*<       MCHI0=    slhaVal1('MASS',zero,1000037) >*/
    mchi0 = slhaVal("MASS", zero,1, c_b19);
/*<       MGLU=     slhaVal1('MASS',zero,1000021)   >*/
    mglu = slhaVal("MASS", zero,1, c_b17);
/*<        >*/
    d__1 = slhaVal("SMINPUTS", zero,1,c__3);
    d__2 = slhaVal("SMINPUTS", zero,1,c__5);
    qcdl = initQCD(d__1, c_b24, d__2, mtop);
/*<        >*/
    susy = sqrt(slhaVal("MASS", zero,1, c_b28) * slhaVal("MASS", zero,1,c_b30));
/*<       AMSQ= slhaVal1('MSOFT',susy,43) >*/
    amsq = slhaVal("MSOFT", susy,1, c__43);
/*<       AMUR= slhaVal1('MSOFT',susy,46) >*/
    amur = slhaVal("MSOFT", susy,1, c__46);
/*<       AMDR= slhaVal1('MSOFT',susy,49) >*/
    amdr = slhaVal("MSOFT", susy,1, c__49);
/*     susy = DSQRT(2*AMSQ**2+AMUR**2+AMDR**2)/2 */
/*<       TANB=     slhaVal1('HMIX',susy,2) >*/
    tanb = slhaVal("HMIX", susy,1, c__2);
/*<       MU=       slhaVal1('HMIX',susy,1)  >*/
    mu = slhaVal("HMIX", susy,1, c__1);
/*<       AD=slhaVal2('AD',susy,3,3) >*/
    ad = slhaVal("AD", susy,2, c__3, c__3);
/*<       AU=slhaVal2('AU',susy,3,3) >*/
    au = slhaVal("AU", susy,2, c__3, c__3);
/*<       ft=1.017363287 >*/
    ft = 1.017363287f;
/*<         BET=DATAN(TANB) >*/
    bet = atan(tanb);
/*<         SB = DSIN(BET) >*/
    sb = sin(bet);
/*<         CB = DCOS(BET) >*/
    cb = cos(bet);
/*  ============ TRANSFORMATION OF INPUT FOR SUBH ========== */
/*<       CF = 4/3.D0 >*/
    cf = 1.3333333333333333;
/*<       Q0 = susy >*/
    q0 = susy;
/*<       ALSP = alphaQCD(Q0)/PI >*/
    alsp = alphaQCD(q0) / pi;
/*<       ALTP = MtRun(Q0)**2/2/PI/V**2/SB**2 / PI *ft >*/
    d__1 = MtRun(q0);
    d__2 = v;
    d__3 = sb;
    altp = d__1 * d__1 / 2 / pi / (d__2 * d__2) / (d__3 * d__3) / pi * ft;
/*<       ALBP = MbRun(Q0)**2/2/PI/V**2/CB**2 / PI >*/
    d__1 = MbRun(q0);
    d__2 = v;
    d__3 = cb;
    albp = d__1 * d__1 / 2 / pi / (d__2 * d__2) / (d__3 * d__3) / pi;
/*<       RMT = MtRun(MTOP)*ft >*/
    rmt = MtRun(mtop) * ft;
/*<       RMB = MbRun(MTOP) >*/
    rmb = MbRun(mtop);
/*<       QT = DSQRT(DMAX1(AMSQ**2+RMT**2,AMUR**2+RMT**2)) >*/
/* Computing MAX */
    d__3 = amsq;
    d__4 = rmt;
    d__5 = amur;
    d__6 = rmt;
    d__1 = d__3 * d__3 + d__4 * d__4, d__2 = d__5 * d__5 + d__6 * d__6;
    qt = sqrt((max(d__1,d__2)));
/*<       QB = DSQRT(DMAX1(AMSQ**2+RMB**2,AMDR**2+RMB**2)) >*/
/* Computing MAX */
    d__3 = amsq;
    d__4 = rmb;
    d__5 = amdr;
    d__6 = rmb;
    d__1 = d__3 * d__3 + d__4 * d__4, d__2 = d__5 * d__5 + d__6 * d__6;
    qb = sqrt((max(d__1,d__2)));
/*<       AMH12 = MA**2*SB**2 - AMZ**2/2*(CB**2-SB**2) - AMU**2 >*/
    d__1 = ma;
    d__2 = sb;
    d__3 = param_hdec__1.amz;
    d__4 = cb;
    d__5 = sb;
    d__6 = amu;
    amh12 = d__1 * d__1 * (d__2 * d__2) - d__3 * d__3 / 2 * (d__4 * d__4 - 
	    d__5 * d__5) - d__6 * d__6;
/*<       AMH22 = MA**2*CB**2 + AMZ**2/2*(CB**2-SB**2) - AMU**2 >*/
    d__1 = ma;
    d__2 = cb;
    d__3 = param_hdec__1.amz;
    d__4 = cb;
    d__5 = sb;
    d__6 = amu;
    amh22 = d__1 * d__1 * (d__2 * d__2) + d__3 * d__3 / 2 * (d__4 * d__4 - 
	    d__5 * d__5) - d__6 * d__6;
/*<       XB = AMSQ**2 + AMDR**2 + AMH12 + AD**2 >*/
    d__1 = amsq;
    d__2 = amdr;
    d__3 = ad;
    xb = d__1 * d__1 + d__2 * d__2 + amh12 + d__3 * d__3;
/*<       XT = AMSQ**2 + AMUR**2 + AMH22 + AU**2 >*/
    d__1 = amsq;
    d__2 = amur;
    d__3 = au;
    xt = d__1 * d__1 + d__2 * d__2 + amh22 + d__3 * d__3;
/*<        >*/
    d__1 = qb;
    d__2 = q0;
    ad0 = ad + (cf * alsp * mglu + altp * 3 / 2 * au + albp / 4 * ad) * log(
	    d__1 * d__1 / (d__2 * d__2));
/*<        >*/
    d__1 = amsq;
    d__2 = mglu;
    d__3 = qb;
    d__4 = q0;
    amdl0 = sqrt(d__1 * d__1 + (-cf * alsp * (d__2 * d__2) + (altp * xt + 
	    albp * xb) / 4) * log(d__3 * d__3 / (d__4 * d__4)));
/*<        >*/
    d__1 = amdr;
    d__2 = mglu;
    d__3 = qb;
    d__4 = q0;
    amdr0 = sqrt(d__1 * d__1 + (-cf * alsp * (d__2 * d__2) + albp * xb / 4) * 
	    log(d__3 * d__3 / (d__4 * d__4)));
/*<        >*/
    d__1 = qt;
    d__2 = q0;
    au0 = au + (cf * alsp * mglu + altp / 4 * au + albp * 3 / 2 * ad) * log(
	    d__1 * d__1 / (d__2 * d__2));
/*<        >*/
    d__1 = amsq;
    d__2 = mglu;
    d__3 = qt;
    d__4 = q0;
    amul0 = sqrt(d__1 * d__1 + (-cf * alsp * (d__2 * d__2) + (altp * xt + 
	    albp * xb) / 4) * log(d__3 * d__3 / (d__4 * d__4)));
/*<        >*/
    d__1 = amur;
    d__2 = mglu;
    d__3 = qt;
    d__4 = q0;
    amur0 = sqrt(d__1 * d__1 + (-cf * alsp * (d__2 * d__2) + altp * xt / 4) * 
	    log(d__3 * d__3 / (d__4 * d__4)));
/*<        >*/
    subh1_hdec__(&ma, &tanb, &amul0, &amdl0, &amur0, &amdr0, &mtop, &au0, &
	    ad0, &mu, &mchi0, &mhp, &hmp, &mch, &sa, &ca, &tanba, &mglu, &dmb)
	    ;
    return 0;
} 

double Lambda1(void){ return hself_hdec__1.lambda1;}
double Lambda2(void){ return hself_hdec__1.lambda2;}
double Lambda3(void){ return hself_hdec__1.lambda3;}
double Lambda4(void){ return hself_hdec__1.lambda4;}
double Lambda5(void){ return hself_hdec__1.lambda5;}
double Lambda6(void){ return hself_hdec__1.lambda6;}
double Lambda7(void){ return hself_hdec__1.lambda7;}

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*       Program based on the work by M. Carena, M. Quiros */
/*       and C.E.M. Wagner, "Effective potential methods and */
/*       the Higgs mass spectrum in the MSSM", Nucl. Phys. */
/*       B461 (1996) 407. */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*<        >*/
/* Subroutine */ int subh1_hdec__(double *ma, double *tanb, 
	double *mul, double *mdl, double *mur, double *md, 
	double *mtop, double *au, double *ad, double *mu, 
	double *mchi0, double *mhp, double *hmp, double *mch, 
	double *sa, double *ca, double *tanba, double *mglu, 
	double *deltamb)
{
    /* System generated locals */
    double d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12, d__13, d__14;

    /* Builtin functions */
    double atan(double), sqrt(double), log(double), pow_dd(
	    double *, double *), asin(double), sin(double), 
	    cos(double);

    /* Local variables */
    static double deltal3p4, dlambdap2, cos2alpha, sin2alpha, v;
    extern /* Subroutine */ int gfun_hdec__(double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *);
    static double g1, g2, g3, m2[4]	/* was [2][2] */, hd, mb, td, pi, cw, 
	    hu, mw, tp, mz, tu, sw, m2p[4]	/* was [2][2] */, tdl, tdp, 
	    tpd, tul, mch2, mh2p, hm2p, mchi, cosb, sinb, tdpd, tglu, dlam1, 
	    dlam2, dlam3, dlam4, dlam5, dlam6, dlam7, cos2b, trm2p, alpha, 
	    tchar, sqbma, tanbt;
    static double rmtop;
    static double alpha1, alpha2, alpha3, detm2p, alpha3z, mssusy, 
	    deltal12, dlambda1, dlambda2, dlambda3, dlambda4, deltam112, 
	    deltam122, deltam222;

/*<       IMPLICIT REAL*8(A-H,L,M,O-Z) >*/
/*<       real*8 MbRun,MtRun >*/
/*<       DIMENSION VH(2,2),M2(2,2),M2P(2,2) >*/
/*<       COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW >*/
/*<        >*/
/*<       MCHI = MCHI0 >*/
    mchi = *mchi0;
/*<       TANBA = TANB >*/
    *tanba = *tanb;
/*<       TANBT = TANB >*/
    tanbt = *tanb;
/*<       PI = 4*DATAN(1D0) >*/
    pi = atan(1.) * 4;
/*<       MZ = AMZ >*/
    mz = param_hdec__1.amz;
/*<       MW = AMW >*/
    mw = param_hdec__1.amw;
/*<       V  = 1/DSQRT(2*DSQRT(2D0)*GF) >*/
    v = 1 / sqrt(sqrt(2.) * 2 * param_hdec__1.gf);
/*<       CW = AMW**2/AMZ**2 >*/
    d__1 = param_hdec__1.amw;
    d__2 = param_hdec__1.amz;
    cw = d__1 * d__1 / (d__2 * d__2);
/*<       SW = 1-CW >*/
    sw = 1 - cw;
/*<       ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI >*/
    d__1 = param_hdec__1.amw * 2 / v / sqrt(2.);
    alpha2 = d__1 * d__1 / 4 / pi;
/*<       ALPHA1  = ALPHA2*SW/CW >*/
    alpha1 = alpha2 * sw / cw;
/*<       ALPHA3Z = alphaQCD(AMZ) >*/
    alpha3z = alphaQCD(param_hdec__1.amz);
/*<       ALPHA3  = alphaQCD(MTOP) >*/
    alpha3 = alphaQCD(*mtop);
/*<       MB      = MbRun(MTOP) >*/
    mb = MbRun(*mtop);
/*<       RMTOP   = MtRun(MTOP) >*/
    rmtop = MtRun(*mtop);
/*<       TUL = LOG((MUL**2+MTOP**2)/MTOP**2) >*/
    d__1 = *mul;
    d__2 = *mtop;
    d__3 = *mtop;
    tul = log((d__1 * d__1 + d__2 * d__2) / (d__3 * d__3));
/*<       TDL = LOG((MDL**2+MTOP**2)/MTOP**2) >*/
    d__1 = *mdl;
    d__2 = *mtop;
    d__3 = *mtop;
    tdl = log((d__1 * d__1 + d__2 * d__2) / (d__3 * d__3));
/*<       TU = LOG((MUR**2 + MTOP**2)/MTOP**2) >*/
    d__1 = *mur;
    d__2 = *mtop;
    d__3 = *mtop;
    tu = log((d__1 * d__1 + d__2 * d__2) / (d__3 * d__3));
/*<       TD = LOG((MD**2 + MTOP**2)/MTOP**2) >*/
    d__1 = *md;
    d__2 = *mtop;
    d__3 = *mtop;
    td = log((d__1 * d__1 + d__2 * d__2) / (d__3 * d__3));
/*<       SINB = TANB/DSQRT(1.D0 + TANB**2) >*/
    d__1 = *tanb;
    sinb = *tanb / sqrt(d__1 * d__1 + 1.);
/*<       COSB = SINB/TANB >*/
    cosb = sinb / *tanb;
/*<        >*/
    if (*ma > *mtop) {
	d__1 = pi;
	d__2 = rmtop;
	d__3 = v;
	d__4 = sinb;
	d__5 = mb;
	d__6 = v;
	d__7 = cosb;
	d__8 = *ma;
	d__9 = *mtop;
	*tanba = *tanb * (1. - .09375 / (d__1 * d__1) * (d__2 * d__2 / (d__3 *
		 d__3) / (d__4 * d__4) - d__5 * d__5 / (d__6 * d__6) / (d__7 *
		 d__7)) * log(d__8 * d__8 / (d__9 * d__9)));
    }
/*<       IF(MA.LT.MTOP.OR.MA.EQ.MTOP) TANBT = TANBA >*/
    if (*ma < *mtop || *ma == *mtop) {
	tanbt = *tanba;
    }
/*<       SINB = TANBT/DSQRT(1.D0 + TANBT**2) >*/
    d__1 = tanbt;
    sinb = tanbt / sqrt(d__1 * d__1 + 1.);
/*<       COSB = 1.D0/DSQRT(1.D0 + TANBT**2) >*/
    d__1 = tanbt;
    cosb = 1. / sqrt(d__1 * d__1 + 1.);
/*<       COS2B = (TANBT**2 - 1.D0)/(TANBT**2 + 1.D0) >*/
    d__1 = tanbt;
    d__2 = tanbt;
    cos2b = (d__1 * d__1 - 1.) / (d__2 * d__2 + 1.);
/*<       G1 = DSQRT(ALPHA1*4.D0*PI) >*/
    g1 = sqrt(alpha1 * 4. * pi);
/*<       G2 = DSQRT(ALPHA2*4.D0*PI) >*/
    g2 = sqrt(alpha2 * 4. * pi);
/*<       G3 = DSQRT(ALPHA3*4.D0*PI) >*/
    g3 = sqrt(alpha3 * 4. * pi);
/*<       HU = RMTOP/V/SINB >*/
    hu = rmtop / v / sinb;
/*<       HD =  MB/V/COSB >*/
    hd = mb / v / cosb;

/*<       IF(MUL.GT.MUR) TP = TUL - TU >*/
    if (*mul > *mur) {
	tp = tul - tu;
    }
/*<       IF(MUL.LT.MUR.OR.MUL.EQ.MUR) TP = TU - TUL >*/
    if (*mul < *mur || *mul == *mur) {
	tp = tu - tul;
    }
/*<       IF(MUL.GT.MUR) TDP = TU >*/
    if (*mul > *mur) {
	tdp = tu;
    }
/*<       IF(MUL.LT.MUR.OR.MUL.EQ.MUR) TDP = TUL >*/
    if (*mul < *mur || *mul == *mur) {
	tdp = tul;
    }
/*<       IF(MDL.GT.MD) TPD = TDL - TD >*/
    if (*mdl > *md) {
	tpd = tdl - td;
    }
/*<       IF(MDL.LT.MD.OR.MDL.EQ.MD) TPD = TD - TDL >*/
    if (*mdl < *md || *mdl == *md) {
	tpd = td - tdl;
    }
/*<       IF(MDL.GT.MD) TDPD = TD >*/
    if (*mdl > *md) {
	tdpd = td;
    }
/*<       IF(MDL.LT.MD.OR.MDL.EQ.MD) TDPD = TDL >*/
    if (*mdl < *md || *mdl == *md) {
	tdpd = tdl;
    }
/*<       IF(MDL.GT.MD) DLAMBDA1 = 6./96./PI**2*G1**2*HD**2*TPD >*/
    if (*mdl > *md) {
	d__1 = pi;
	d__2 = g1;
	d__3 = hd;
	dlambda1 = .0625f / (d__1 * d__1) * (d__2 * d__2) * (d__3 * d__3) * 
		tpd;
    }
/*<        >*/
    if (*mdl < *md || *mdl == *md) {
	d__1 = pi;
	d__2 = hd;
	d__3 = g1;
	d__4 = g2;
	dlambda1 = .09375f / (d__1 * d__1) * (d__2 * d__2) * (d__3 * d__3 / 
		3.f + d__4 * d__4) * tpd;
    }
/*<       IF(MUL.GT.MUR) DLAMBDA2 =12./96./PI**2*G1**2*HU**2*TP >*/
    if (*mul > *mur) {
	d__1 = pi;
	d__2 = g1;
	d__3 = hu;
	dlambda2 = .125f / (d__1 * d__1) * (d__2 * d__2) * (d__3 * d__3) * tp;
    }
/*<        >*/
    if (*mul < *mur || *mul == *mur) {
	d__1 = pi;
	d__2 = hu;
	d__3 = g1;
	d__4 = g2;
	dlambda2 = .09375f / (d__1 * d__1) * (d__2 * d__2) * (-(d__3 * d__3) /
		 3.f + d__4 * d__4) * tp;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*  dlambdap1 and dlambdap2 are the new log corrections due to */
/*  the presence of the gluino mass. They are in general very small, */
/*  and only present if there is a hierarchy of masses between the */
/*  two stops. */


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*<         dlambdap2 = 0 >*/
    dlambdap2 = 0.;
/*<         tglu = log(mglu**2/mtop**2) >*/
    d__1 = *mglu;
    d__2 = *mtop;
    tglu = log(d__1 * d__1 / (d__2 * d__2));
/*<         if(mglu.lt.mur.or.mglu.lt.mul) then >*/
    if (*mglu < *mur || *mglu < *mul) {
/*<         if(mul.gt.mur.and.mglu.gt.mur) then >*/
	if (*mul > *mur && *mglu > *mur) {
/*<         dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tul**2-tglu**2) >*/
	    d__2 = pi;
	    d__1 = d__2 * d__2 * 16.f;
/* Computing 4th power */
	    d__3 = hu, d__3 *= d__3;
	    d__4 = tul;
	    d__5 = tglu;
	    dlambdap2 = -4.f / (d__1 * d__1) * (d__3 * d__3) * (d__4 * d__4 - 
		    d__5 * d__5);
/*<         endif >*/
	}
/*<         if(mul.gt.mur.and.mglu.lt.mur) then >*/
	if (*mul > *mur && *mglu < *mur) {
/*<         dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tul**2-tu**2) >*/
	    d__2 = pi;
	    d__1 = d__2 * d__2 * 16.f;
/* Computing 4th power */
	    d__3 = hu, d__3 *= d__3;
	    d__4 = tul;
	    d__5 = tu;
	    dlambdap2 = -4.f / (d__1 * d__1) * (d__3 * d__3) * (d__4 * d__4 - 
		    d__5 * d__5);
/*<         endif >*/
	}
/*<         if(mul.gt.mur.and.mglu.eq.mur) then >*/
	if (*mul > *mur && *mglu == *mur) {
/*<         dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tul**2-tu**2) >*/
	    d__2 = pi;
	    d__1 = d__2 * d__2 * 16.f;
/* Computing 4th power */
	    d__3 = hu, d__3 *= d__3;
	    d__4 = tul;
	    d__5 = tu;
	    dlambdap2 = -4.f / (d__1 * d__1) * (d__3 * d__3) * (d__4 * d__4 - 
		    d__5 * d__5);
/*<         endif >*/
	}
/*<         if(mur.gt.mul.and.mglu.gt.mul) then >*/
	if (*mur > *mul && *mglu > *mul) {
/*<         dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tglu**2) >*/
	    d__2 = pi;
	    d__1 = d__2 * d__2 * 16.f;
/* Computing 4th power */
	    d__3 = hu, d__3 *= d__3;
	    d__4 = tu;
	    d__5 = tglu;
	    dlambdap2 = -4.f / (d__1 * d__1) * (d__3 * d__3) * (d__4 * d__4 - 
		    d__5 * d__5);
/*<         endif >*/
	}
/*<         if(mur.gt.mul.and.mglu.lt.mul) then >*/
	if (*mur > *mul && *mglu < *mul) {
/*<         dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tul**2) >*/
	    d__2 = pi;
	    d__1 = d__2 * d__2 * 16.f;
/* Computing 4th power */
	    d__3 = hu, d__3 *= d__3;
	    d__4 = tu;
	    d__5 = tul;
	    dlambdap2 = -4.f / (d__1 * d__1) * (d__3 * d__3) * (d__4 * d__4 - 
		    d__5 * d__5);
/*<         endif >*/
	}
/*<         if(mur.gt.mul.and.mglu.eq.mul) then >*/
	if (*mur > *mul && *mglu == *mul) {
/*<         dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tul**2) >*/
	    d__2 = pi;
	    d__1 = d__2 * d__2 * 16.f;
/* Computing 4th power */
	    d__3 = hu, d__3 *= d__3;
	    d__4 = tu;
	    d__5 = tul;
	    dlambdap2 = -4.f / (d__1 * d__1) * (d__3 * d__3) * (d__4 * d__4 - 
		    d__5 * d__5);
/*<         endif >*/
	}
/*<         endif >*/
    }
/*<       DLAMBDA3 = 0. >*/
    dlambda3 = 0.f;
/*<       DLAMBDA4 = 0. >*/
    dlambda4 = 0.f;
/*<       IF(MDL.GT.MD) DLAMBDA3 = -1./32./PI**2*G1**2*HD**2*TPD >*/
    if (*mdl > *md) {
	d__1 = pi;
	d__2 = g1;
	d__3 = hd;
	dlambda3 = -.03125f / (d__1 * d__1) * (d__2 * d__2) * (d__3 * d__3) * 
		tpd;
    }
/*<        >*/
    if (*mdl < *md || *mdl == *md) {
	d__1 = pi;
	d__2 = hd;
	d__3 = g2;
	d__4 = g1;
	dlambda3 = .046875f / (d__1 * d__1) * (d__2 * d__2) * (d__3 * d__3 - 
		d__4 * d__4 / 3.f) * tpd;
    }
/*<        >*/
    if (*mul > *mur) {
	d__1 = pi;
	d__2 = g1;
	d__3 = hu;
	dlambda3 -= .0625f / (d__1 * d__1) * (d__2 * d__2) * (d__3 * d__3) * 
		tp;
    }
/*<        >*/
    if (*mul < *mur || *mul == *mur) {
	d__1 = pi;
	d__2 = hu;
	d__3 = g2;
	d__4 = g1;
	dlambda3 += .046875f / (d__1 * d__1) * (d__2 * d__2) * (d__3 * d__3 + 
		d__4 * d__4 / 3.f) * tp;
    }
/*<       IF(MUL.LT.MUR) DLAMBDA4 = -3./32./PI**2*G2**2*HU**2*TP >*/
    if (*mul < *mur) {
	d__1 = pi;
	d__2 = g2;
	d__3 = hu;
	dlambda4 = -.09375f / (d__1 * d__1) * (d__2 * d__2) * (d__3 * d__3) * 
		tp;
    }
/*<        >*/
    if (*mdl < *md) {
	d__1 = pi;
	d__2 = g2;
	d__3 = hd;
	dlambda4 -= .09375f / (d__1 * d__1) * (d__2 * d__2) * (d__3 * d__3) * 
		tpd;
    }

/*<        >*/
    d__1 = g1;
    d__2 = g2;
    d__3 = hd;
    d__4 = pi;
    d__5 = pi;
    d__6 = hd;
    d__7 = hu;
    d__8 = g3;
    d__9 = pi;
    d__10 = pi;
    d__11 = hd;
    d__12 = hu;
    d__13 = g3;
    d__14 = pi;
    hself_hdec__1.lambda1 = (d__1 * d__1 + d__2 * d__2) / 4.f * (1.f - d__3 * 
	    d__3 * 3.f * (tpd + tdpd) / 8.f / (d__4 * d__4)) + pow(hd,
	    c_b83) * 3.f / 16.f / (d__5 * d__5) * tpd * ((d__6 * d__6 * 3.f / 
	    2.f + d__7 * d__7 / 2.f - d__8 * d__8 * 8.f) * (tpd + tdpd * 2.f) 
	    / 16.f / (d__9 * d__9) + 1.f) + pow(hd, c_b83) * 3.f / 8.f / 
	    (d__10 * d__10) * tdpd * ((d__11 * d__11 * 3.f / 2.f + d__12 * 
	    d__12 / 2.f - d__13 * d__13 * 8.f) * tdpd / 16.f / (d__14 * d__14)
	     + 1.f) + dlambda1;

/*<        >*/
    d__1 = g1;
    d__2 = g2;
    d__3 = hu;
    d__4 = pi;
    d__5 = pi;
    d__6 = hu;
    d__7 = hd;
    d__8 = g3;
    d__9 = pi;
    d__10 = pi;
    d__11 = hu;
    d__12 = hd;
    d__13 = g3;
    d__14 = pi;
    hself_hdec__1.lambda2 = (d__1 * d__1 + d__2 * d__2) / 4.f * (1.f - d__3 * 
	    d__3 * 3.f * (tp + tdp) / 8.f / (d__4 * d__4)) + pow(hu, 
	    c_b83) * 3.f / 16.f / (d__5 * d__5) * tp * ((d__6 * d__6 * 3.f / 
	    2.f + d__7 * d__7 / 2.f - d__8 * d__8 * 8.f) * (tp + tdp * 2.f) / 
	    16.f / (d__9 * d__9) + 1.f) + pow(hu, c_b83) * 3.f / 8.f / (
	    d__10 * d__10) * tdp * ((d__11 * d__11 * 3.f / 2.f + d__12 * 
	    d__12 / 2.f - d__13 * d__13 * 8.f) * tdp / 16.f / (d__14 * d__14) 
	    + 1.f) + dlambda2 + dlambdap2;

/*<        >*/
    d__1 = g2;
    d__2 = g1;
    d__3 = hu;
    d__4 = pi;
    d__5 = hd;
    d__6 = pi;
    hself_hdec__1.lambda3 = (d__1 * d__1 - d__2 * d__2) / 4.f * (1.f - d__3 * 
	    d__3 * 3.f * (tp + tdp) / 16.f / (d__4 * d__4) - d__5 * d__5 * 
	    3.f * (tpd + tdpd) / 16.f / (d__6 * d__6)) + dlambda3;

/*<        >*/
    d__1 = g2;
    d__2 = hu;
    d__3 = pi;
    d__4 = hd;
    d__5 = pi;
    hself_hdec__1.lambda4 = -(d__1 * d__1) / 2.f * (1.f - d__2 * d__2 * 3.f * 
	    (tp + tdp) / 16.f / (d__3 * d__3) - d__4 * d__4 * 3.f * (tpd + 
	    tdpd) / 16.f / (d__5 * d__5)) + dlambda4;

/*< 	LAMBDA5 = 0. >*/
    hself_hdec__1.lambda5 = 0.f;
/*< 	LAMBDA6 = 0. >*/
    hself_hdec__1.lambda6 = 0.f;
/*< 	LAMBDA7 = 0. >*/
    hself_hdec__1.lambda7 = 0.f;

/*     THIS IS THE CONTRIBUTION FROM LIGHT CHARGINOS/NEUTRALINOS */
/*     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*<   	 MSSUSY=DSQRT(0.5D0*(MUL**2+MUR**2)+MTOP**2) >*/
    d__1 = *mul;
    d__2 = *mur;
    d__3 = *mtop;
    mssusy = sqrt((d__1 * d__1 + d__2 * d__2) * .5 + d__3 * d__3);
/*< 	IF(MCHI.GT.MSSUSY)GOTO 3790 >*/
    if (mchi > mssusy) {
	goto L3790;
    }
/*< 	IF(MCHI.LT.MTOP) MCHI=MTOP >*/
    if (mchi < *mtop) {
	mchi = *mtop;
    }
/*< 	TCHAR=LOG(MSSUSY**2/MCHI**2) >*/
    d__1 = mssusy;
    d__2 = mchi;
    tchar = log(d__1 * d__1 / (d__2 * d__2));
/*< 	DELTAL12=(9./64./PI**2*G2**4+5./192./PI**2*G1**4)*TCHAR >*/
    d__1 = pi;
/* Computing 4th power */
    d__2 = g2, d__2 *= d__2;
    d__3 = pi;
/* Computing 4th power */
    d__4 = g1, d__4 *= d__4;
    deltal12 = (.140625f / (d__1 * d__1) * (d__2 * d__2) + 
	    .026041666666666668f / (d__3 * d__3) * (d__4 * d__4)) * tchar;
/*< 	DELTA >*/
    d__1 = pi;
/* Computing 4th power */
    d__2 = g2, d__2 *= d__2;
    d__3 = pi;
/* Computing 4th power */
    d__4 = g1, d__4 *= d__4;
    d__5 = pi;
    d__6 = g1;
    d__7 = g2;
    deltal3p4 = (.046875f / (d__1 * d__1) * (d__2 * d__2) + 
	    .036458333333333336f / (d__3 * d__3) * (d__4 * d__4) + .125f / (
	    d__5 * d__5) * (d__6 * d__6) * (d__7 * d__7)) * tchar;
/*< 	DELTAM112=2.*DELTAL12*V**2*COSB**2 >*/
    d__1 = v;
    d__2 = cosb;
    deltam112 = deltal12 * 2.f * (d__1 * d__1) * (d__2 * d__2);
/*< 	DELTAM222=2.*DELTAL12*V**2*SINB**2 >*/
    d__1 = v;
    d__2 = sinb;
    deltam222 = deltal12 * 2.f * (d__1 * d__1) * (d__2 * d__2);
/*< 	DELTAM122=2.*DELTAL3P4*V**2*SINB*COSB >*/
    d__1 = v;
    deltam122 = deltal3p4 * 2.f * (d__1 * d__1) * sinb * cosb;
/* --EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I */
/*<         DLAM1 = DELTAM112/2.D0/V**2/COSB**2 >*/
    d__1 = v;
    d__2 = cosb;
    dlam1 = deltam112 / 2. / (d__1 * d__1) / (d__2 * d__2);
/*<         DLAM2 = DELTAM222/2.D0/V**2/SINB**2 >*/
    d__1 = v;
    d__2 = sinb;
    dlam2 = deltam222 / 2. / (d__1 * d__1) / (d__2 * d__2);
/*<        >*/
    d__1 = v;
    d__2 = g1;
    d__3 = g2;
    d__4 = g1;
    d__5 = g2;
    dlam3 = deltam122 / 2. / (d__1 * d__1) / sinb / cosb * (d__2 * d__2 - 
	    d__3 * d__3) / (d__4 * d__4 + d__5 * d__5);
/*<        >*/
    d__1 = v;
    d__2 = g2;
    d__3 = g1;
    d__4 = g2;
    dlam4 = deltam122 / 2. / (d__1 * d__1) / sinb / cosb * (d__2 * d__2 * 2) /
	     (d__3 * d__3 + d__4 * d__4);
/*<         LAMBDA1 = LAMBDA1+DLAM1 >*/
    if(isfinite(dlam1)) hself_hdec__1.lambda1 += dlam1;
/*<         LAMBDA2 = LAMBDA2+DLAM2 >*/
    if(isfinite(dlam2)) hself_hdec__1.lambda2 += dlam2;
/*<         LAMBDA3 = LAMBDA3+DLAM3 >*/
    if(isfinite(dlam3)) hself_hdec__1.lambda3 += dlam3;
/*<         LAMBDA4 = LAMBDA4+DLAM4 >*/
    if(isfinite(dlam4)) hself_hdec__1.lambda4 += dlam4;
/* --END OF EXTENSION */
/*<  3790	CONTINUE >*/
L3790:
/* CCCCCCCCCCCCCC    END OF CHARGINOS AND NEUTRALINOS  CCCCCCCCCCCC */
/* --EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I */
/*<        >*/
    gfun_hdec__(ma, tanba, mul, mdl, mur, md, mtop, au, ad, mu, mglu, &dlam1, 
	    &dlam2, &dlam3, &dlam4, &dlam5, &dlam6, &dlam7, deltamb);
/*<       LAMBDA1 = LAMBDA1+DLAM1 >*/
    if(isfinite(dlam1)) hself_hdec__1.lambda1 += dlam1;
/*<       LAMBDA2 = LAMBDA2+DLAM2 >*/
    if(isfinite(dlam2)) hself_hdec__1.lambda2 += dlam2;
/*<       LAMBDA3 = LAMBDA3+DLAM3 >*/
    if(isfinite(dlam3)) hself_hdec__1.lambda3 += dlam3;
/*<       LAMBDA4 = LAMBDA4+DLAM4 >*/
    if(isfinite(dlam4)) hself_hdec__1.lambda4 += dlam4;
/*<       LAMBDA5 = LAMBDA5+DLAM5 >*/
    if(isfinite(dlam5)) hself_hdec__1.lambda5 += dlam5;
/*<       LAMBDA6 = LAMBDA6+DLAM6 >*/
    if(isfinite(dlam6)) hself_hdec__1.lambda6 += dlam6;
/*<       LAMBDA7 = LAMBDA7+DLAM7 >*/
    if(isfinite(dlam7)) hself_hdec__1.lambda7 += dlam7;
/*<        >*/
    d__1 = v;
    d__2 = cosb;
    d__3 = sinb;
    d__4 = *ma;
    d__5 = sinb;
    m2[0] = d__1 * d__1 * 2.f * (hself_hdec__1.lambda1 * (d__2 * d__2) + 
	    hself_hdec__1.lambda6 * 2.f * cosb * sinb + hself_hdec__1.lambda5 
	    * (d__3 * d__3)) + d__4 * d__4 * (d__5 * d__5);
/*<        >*/
    d__1 = v;
    d__2 = cosb;
    d__3 = sinb;
    d__4 = *ma;
    d__5 = cosb;
    m2[3] = d__1 * d__1 * 2.f * (hself_hdec__1.lambda5 * (d__2 * d__2) + 
	    hself_hdec__1.lambda7 * 2.f * cosb * sinb + hself_hdec__1.lambda2 
	    * (d__3 * d__3)) + d__4 * d__4 * (d__5 * d__5);
/*<        >*/
    d__1 = v;
    d__2 = cosb;
    d__3 = sinb;
    d__4 = *ma;
    m2[2] = d__1 * d__1 * 2.f * (hself_hdec__1.lambda6 * (d__2 * d__2) + (
	    hself_hdec__1.lambda3 + hself_hdec__1.lambda4) * cosb * sinb + 
	    hself_hdec__1.lambda7 * (d__3 * d__3)) - d__4 * d__4 * sinb * 
	    cosb;
/*<       M2(2,1) = M2(1,2) >*/
    m2[1] = m2[2];
/*<       M2P(1,1) = M2(1,1) >*/
    m2p[0] = m2[0];
/*<       M2P(2,2) = M2(2,2) >*/
    m2p[3] = m2[3];
/*<       M2P(1,2) = M2(1,2) >*/
    m2p[2] = m2[2];
/*<       M2P(2,1) = M2(2,1) >*/
    m2p[1] = m2[1];
/* --END OF EXTENSION */
/*<       TRM2P  = M2P(1,1) + M2P(2,2) >*/
    trm2p = m2p[0] + m2p[3];
/*<       DETM2P = M2P(1,1)*M2P(2,2) - M2P(1,2)*M2P(2,1) >*/
    detm2p = m2p[0] * m2p[3] - m2p[2] * m2p[1];
/*<       MH2P = (TRM2P - DSQRT(TRM2P**2 - 4.D0* DETM2P))/2.D0 >*/
    d__1 = trm2p;
    mh2p = (trm2p - sqrt(d__1 * d__1 - detm2p * 4.)) / 2.;
/*<       HM2P = (TRM2P + DSQRT(TRM2P**2 - 4.D0* DETM2P))/2.D0 >*/
    d__1 = trm2p;
    hm2p = (trm2p + sqrt(d__1 * d__1 - detm2p * 4.)) / 2.;
/* !!!!!!!!!!!!!!!!!!! */
/*<       MCH2=MA**2+(LAMBDA5-LAMBDA4)*V**2 >*/
    d__1 = *ma;
    d__2 = v;
    mch2 = d__1 * d__1 + (hself_hdec__1.lambda5 - hself_hdec__1.lambda4) * (
	    d__2 * d__2);
/* !!!!!!!!!!!!!!!!!!! */
/*<       MCH=DSQRT(MCH2) >*/
    *mch = sqrt(mch2);
/*<       HMP = DSQRT(HM2P)  >*/
    *hmp = sqrt(hm2p);
/*<       IF(MH2P.LT.0.)GOTO 5555 >*/
    if (mh2p < 0.f) {
	goto L5555;
    }
/*<       MHP = DSQRT(MH2P)  >*/
    *mhp = sqrt(mh2p);

/*<       SIN2ALPHA = 2.*M2P(1,2)/DSQRT(TRM2P**2-4.D0*DETM2P) >*/
    d__1 = trm2p;
    sin2alpha = m2p[2] * 2.f / sqrt(d__1 * d__1 - detm2p * 4.);
/*<       COS2ALPHA = (M2P(1,1)-M2P(2,2))/DSQRT(TRM2P**2-4.D0*DETM2P) >*/
    d__1 = trm2p;
    cos2alpha = (m2p[0] - m2p[3]) / sqrt(d__1 * d__1 - detm2p * 4.);
/*<       IF(COS2ALPHA.GT.0.) ALPHA = DASIN(SIN2ALPHA)/2.D0 >*/
    if (cos2alpha > 0.f) {
	alpha = asin(sin2alpha) / 2.;
    }
/*<       IF(COS2ALPHA.LT.0.) ALPHA = -PI/2.D0-DASIN(SIN2ALPHA)/2.D0 >*/
    if (cos2alpha < 0.f) {
	alpha = -pi / 2. - asin(sin2alpha) / 2.;
    }
/*<       SA = DSIN(ALPHA) >*/
    *sa = sin(alpha);
/*<       CA = DCOS(ALPHA)   >*/
    *ca = cos(alpha);
/*<       SQBMA = (SINB*CA - COSB*SA)**2 >*/
    d__1 = sinb * *ca - cosb * *sa;
    sqbma = d__1 * d__1;
/*< 5555  RETURN >*/
L5555:
    return 0;
/*<       END >*/
} /* subh1_hdec__ */

/* CCCCCCCCCCCCCCCCCCCCCCC NON DEGENERATE STOP/SBOTTOM EFFECTS CCCCCCCCC */

/*<        >*/
/* Subroutine */ int gfun_hdec__(double *ma, double *tanb, double 
	*mul, double *mdl, double *mur, double *md, double *
	mtop, double *at, double *ab, double *mu, double *
	mglu, double *dlam1, double *dlam2, double *dlam3, 
	double *dlam4, double *dlam5, double *dlam6, double *
	dlam7, double *deltamb)
{
    /* System generated locals */
    double d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12, d__13, d__14, d__15, d__16;

    /* Builtin functions */
    double sqrt(double), atan(double), log(double), pow_dd(
	    double *, double *);

    /* Local variables */
    static double v, g1, g2, g3, hb, g32, al[4]	/* was [2][2] */, mb, 
	    tb;
    extern /* Subroutine */ int delmb_hdec__(double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *);
    static double td, pi, cw, ht, vi, mw, tq, mz, tt, tu, sw, bb2, f1b, 
	    f2b, h1b, h2b, al1, al2, h1i, md2, h2i, bt2, f1t, f2t, h1t, h2t, 
	    vh1[4]	/* was [2][2] */, vh2[4]	/* was [2][2] */, msb,
	     tqd, mst, stw, mdl2, mul2, mur2, facb, fact, cosb, sinb, alst, 
	    htst, mbot2, mbot4, sbot1, sbot2, bt2st, mtop2, mtop4, stop1, 
	    stop2, tanba, cosba, cosbb, sinba, sinbb, cosbt, sbot12, sbot22, 
	    sinbt;
    static double stop12, stop22, dlam1b, rmtop, dlam2b, dlam3b, dlam4b;
    static double alpha1, alpha2, alpha3, dlam5b, dlam6b, dlam7b, dlam1t, 
	    dlam2t, dlam3t, dlam4t, dlam5t, dlam6t, dlam7t, tanbsb, factor, 
	    tanbst, msusyb, alpha3z, msusyt, deltamt;

/*<         IMPLICIT REAL*8 (A-H,L,M,O-Z) >*/
/*<         real*8 MbRun,MtRun >*/
/*<        >*/
/*<         COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW >*/
/*<         G(X,Y) = 2.D0 - (X+Y)/(X-Y)*DLOG(X/Y) >*/
/*<         IF(DABS(MU).LT.0.000001) MU = 0.000001 >*/
    if (abs(*mu) < 1e-6f) {
	*mu = 1e-6f;
    }
/*<         MUL2  = MUL**2 >*/
    d__1 = *mul;
    mul2 = d__1 * d__1;
/*<         MDL2  = MDL**2 >*/
    d__1 = *mdl;
    mdl2 = d__1 * d__1;
/*<         MUR2  = MUR**2 >*/
    d__1 = *mur;
    mur2 = d__1 * d__1;
/*<         MD2   = MD**2 >*/
    d__1 = *md;
    md2 = d__1 * d__1;
/*<         TANBA = TANB >*/
    tanba = *tanb;
/*<         SINBA = TANBA/DSQRT(TANBA**2+1.D0) >*/
    d__1 = tanba;
    sinba = tanba / sqrt(d__1 * d__1 + 1.);
/*<         COSBA = SINBA/TANBA         >*/
    cosba = sinba / tanba;
/*<         SINB = TANB/DSQRT(TANB**2+1.D0) >*/
    d__1 = *tanb;
    sinb = *tanb / sqrt(d__1 * d__1 + 1.);
/*<         COSB = SINB/TANB >*/
    cosb = sinb / *tanb;
/*<        MB = MbRun(MTOP) >*/
    mb = MbRun(*mtop);
/*<       PI = 4*DATAN(1D0) >*/
    pi = atan(1.) * 4;
/*<       MZ = AMZ >*/
    mz = param_hdec__1.amz;
/*<       MW = AMW >*/
    mw = param_hdec__1.amw;
/*<       V  = 1/DSQRT(2*DSQRT(2D0)*GF) >*/
    v = 1 / sqrt(sqrt(2.) * 2 * param_hdec__1.gf);
/*<       CW = AMW**2/AMZ**2 >*/
    d__1 = param_hdec__1.amw;
    d__2 = param_hdec__1.amz;
    cw = d__1 * d__1 / (d__2 * d__2);
/*<       SW = 1-CW >*/
    sw = 1 - cw;
/*<       ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI >*/
    d__1 = param_hdec__1.amw * 2 / v / sqrt(2.);
    alpha2 = d__1 * d__1 / 4 / pi;
/*<       ALPHA1  = ALPHA2*SW/CW >*/
    alpha1 = alpha2 * sw / cw;
/*<       ALPHA3Z = alphaQCD(AMZ) >*/
    alpha3z = alphaQCD(param_hdec__1.amz);
/*<       ALPHA3  = alphaQCD(MTOP) >*/
    alpha3 = alphaQCD(*mtop);
/*<       G1 = DSQRT(ALPHA1*4.*PI) >*/
    g1 = sqrt(alpha1 * 4.f * pi);
/*<       G2 = DSQRT(ALPHA2*4.*PI) >*/
    g2 = sqrt(alpha2 * 4.f * pi);
/*<       G3 = DSQRT(ALPHA3*4.*PI) >*/
    g3 = sqrt(alpha3 * 4.f * pi);
/*<         IF(MUL.GT.MUR) MST = MUL >*/
    if (*mul > *mur) {
	mst = *mul;
    }
/*<         IF(MUR.GT.MUL.OR.MUR.EQ.MUL) MST = MUR >*/
    if (*mur > *mul || *mur == *mul) {
	mst = *mur;
    }
/*<         MSUSYT = DSQRT(MST**2  + MTOP**2) >*/
    d__1 = mst;
    d__2 = *mtop;
    msusyt = sqrt(d__1 * d__1 + d__2 * d__2);
/*< 	IF(MDL.GT.MD) MSB = MDL >*/
    if (*mdl > *md) {
	msb = *mdl;
    }
/*< 	IF(MD.GT.MDL.OR.MD.EQ.MDL) MSB = MD >*/
    if (*md > *mdl || *md == *mdl) {
	msb = *md;
    }
/*< 	MSUSYB = DSQRT(MSB**2 + MB**2) >*/
    d__1 = msb;
    d__2 = mb;
    msusyb = sqrt(d__1 * d__1 + d__2 * d__2);
/*< 	TT = LOG(MSUSYT**2/MTOP**2) >*/
    d__1 = msusyt;
    d__2 = *mtop;
    tt = log(d__1 * d__1 / (d__2 * d__2));
/*< 	TB = LOG(MSUSYB**2/MTOP**2) >*/
    d__1 = msusyb;
    d__2 = *mtop;
    tb = log(d__1 * d__1 / (d__2 * d__2));
/*<         RMTOP   = MtRun(MTOP) *1.017363287 >*/
    rmtop = MtRun(*mtop) * 1.017363287f;
/*<         HT = RMTOP/V/SINB >*/
    ht = rmtop / v / sinb;
/*<         HTST = RMTOP/V >*/
    htst = rmtop / v;
/*<         HB =  MB/V/COSB >*/
    hb = mb / v / cosb;
/*<         G32 = ALPHA3*4.*PI >*/
    g32 = alpha3 * 4.f * pi;
/*<         BT2 = -(8.*G32 - 9.*HT**2/2. - HB**2/2.)/(4.*PI)**2 >*/
    d__1 = ht;
    d__2 = hb;
    d__3 = pi * 4.f;
    bt2 = -(g32 * 8.f - d__1 * d__1 * 9.f / 2.f - d__2 * d__2 / 2.f) / (d__3 *
	     d__3);
/*< 	BB2 = -(8.*G32 - 9.*HB**2/2. - HT**2/2.)/(4.*PI)**2 >*/
    d__1 = hb;
    d__2 = ht;
    d__3 = pi * 4.f;
    bb2 = -(g32 * 8.f - d__1 * d__1 * 9.f / 2.f - d__2 * d__2 / 2.f) / (d__3 *
	     d__3);
/*<         AL2 = 3./8./PI**2*HT**2 >*/
    d__1 = pi;
    d__2 = ht;
    al2 = .375f / (d__1 * d__1) * (d__2 * d__2);
/*<         BT2ST = -(8.*G32 - 9.*HTST**2/2.)/(4.*PI)**2 >*/
    d__1 = htst;
    d__2 = pi * 4.f;
    bt2st = -(g32 * 8.f - d__1 * d__1 * 9.f / 2.f) / (d__2 * d__2);
/*<         ALST = 3./8./PI**2*HTST**2 >*/
    d__1 = pi;
    d__2 = htst;
    alst = .375f / (d__1 * d__1) * (d__2 * d__2);
/*<         AL1 = 3./8./PI**2*HB**2 >*/
    d__1 = pi;
    d__2 = hb;
    al1 = .375f / (d__1 * d__1) * (d__2 * d__2);
/*<         AL(1,1) = AL1 >*/
    al[0] = al1;
/*<         AL(1,2) = (AL2+AL1)/2. >*/
    al[2] = (al2 + al1) / 2.f;
/*<         AL(2,1) = (AL2+AL1)/2. >*/
    al[1] = (al2 + al1) / 2.f;
/*<         AL(2,2) = AL2 >*/
    al[3] = al2;
/*< 	IF(MA.GT.MTOP) THEN >*/
    if (*ma > *mtop) {
/*<         VI = V*(1. + 3./32./PI**2*HTST**2*LOG(MTOP**2/MA**2)) >*/
	d__1 = pi;
	d__2 = htst;
	d__3 = *mtop;
	d__4 = *ma;
	vi = v * (.09375f / (d__1 * d__1) * (d__2 * d__2) * log(d__3 * d__3 / 
		(d__4 * d__4)) + 1.f);
/*<         H1I = VI*COSBA >*/
	h1i = vi * cosba;
/*<         H2I = VI*SINBA >*/
	h2i = vi * sinba;
/*<         H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYT**2))**.25 >*/
	d__2 = pi;
	d__3 = hb;
	d__4 = *ma;
	d__5 = msusyt;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h1t = h1i * pow(d__1, c_b90);
/*<         H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYT**2))**.25 >*/
	d__2 = pi;
	d__3 = ht;
	d__4 = *ma;
	d__5 = msusyt;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h2t = h2i * pow(d__1, c_b90);
/*<         H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYB**2))**.25 >*/
	d__2 = pi;
	d__3 = hb;
	d__4 = *ma;
	d__5 = msusyb;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h1b = h1i * pow(d__1, c_b90);
/*<         H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYB**2))**.25 >*/
	d__2 = pi;
	d__3 = ht;
	d__4 = *ma;
	d__5 = msusyb;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h2b = h2i * pow(d__1, c_b90);
/*< 	ELSE >*/
    } else {
/*< 	VI =  V >*/
	vi = v;
/*< 	H1I = VI*COSB >*/
	h1i = vi * cosb;
/*< 	H2I = VI*SINB >*/
	h2i = vi * sinb;
/*<         H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYT**2))**.25 >*/
	d__2 = pi;
	d__3 = hb;
	d__4 = *mtop;
	d__5 = msusyt;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h1t = h1i * pow(d__1, c_b90);
/*<         H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYT**2))**.25 >*/
	d__2 = pi;
	d__3 = ht;
	d__4 = *mtop;
	d__5 = msusyt;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h2t = h2i * pow(d__1, c_b90);
/*<         H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYB**2))**.25 >*/
	d__2 = pi;
	d__3 = hb;
	d__4 = *mtop;
	d__5 = msusyb;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h1b = h1i * pow(d__1, c_b90);
/*<         H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYB**2))**.25 >*/
	d__2 = pi;
	d__3 = ht;
	d__4 = *mtop;
	d__5 = msusyb;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h2b = h2i * pow(d__1, c_b90);
/*< 	END IF >*/
    }
/*<         TANBST = H2T/H1T >*/
    tanbst = h2t / h1t;
/*<         SINBT = TANBST/(1.+TANBST**2)**.5 >*/
    d__2 = tanbst;
    d__1 = d__2 * d__2 + 1.f;
    sinbt = tanbst / pow(d__1, c_b98);
/*<         COSBT = SINBT/TANBST >*/
    cosbt = sinbt / tanbst;
/*<         TANBSB = H2B/H1B >*/
    tanbsb = h2b / h1b;
/*<         SINBB = TANBSB/(1.+TANBSB**2)**.5 >*/
    d__2 = tanbsb;
    d__1 = d__2 * d__2 + 1.f;
    sinbb = tanbsb / pow(d__1, c_b98);
/*<         COSBB = SINBB/TANBSB >*/
    cosbb = sinbb / tanbsb;
/*<        >*/
    delmb_hdec__(ma, tanb, mul, mdl, mur, md, at, ab, mu, mglu, mtop, &
	    deltamt, deltamb, &stop12, &stop22, &sbot12, &sbot22);
/*<         IF(STOP22.LT.0.) GOTO 4237 >*/
    if (stop22 < 0.f) {
	goto L4237;
    }
/*<         IF(SBOT22.LT.0.) GOTO 4237 >*/
    if (sbot22 < 0.f) {
	goto L4237;
    }
/*<         STOP1 = STOP12**.5 >*/
    stop1 = pow(stop12, c_b98);
/*<         STOP2 = STOP22**.5 >*/
    stop2 = pow(stop22, c_b98);
/*<         SBOT1 = SBOT12**.5 >*/
    sbot1 = pow(sbot12, c_b98);
/*<         SBOT2 = SBOT22**.5 >*/
    sbot2 = pow(sbot22, c_b98);
/*<         mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt) >*/
    mtop4 = pow(rmtop, c_b83) * (bt2 * 2.f * tt + 1.f - al2 * tt - 
	    deltamt * 4.f);
/*     * /(1.+deltamt)**4. */
/*<        >*/
    d__1 = *deltamb + 1.f;
    mbot4 = pow(mb, c_b83) * (bb2 * 2.f * tb + 1.f - al1 * tb) / pow(
	    d__1, c_b83);
/*<         MTOP2 = DSQRT(MTOP4) >*/
    mtop2 = sqrt(mtop4);
/*<         MBOT2 = DSQRT(MBOT4) >*/
    mbot2 = sqrt(fabs(mbot4));       //  fabs by AP
/*<         mb = mb/(1+deltamb) >*/
    mb /= *deltamb + 1;
/*<         VH1(1,1) = 1./TANBST >*/
    vh1[0] = 1.f / tanbst;
/*<         VH1(2,1) = -1. >*/
    vh1[1] = -1.f;
/*<         VH1(1,2) = -1. >*/
    vh1[2] = -1.f;
/*<         VH1(2,2) = TANBST >*/
    vh1[3] = tanbst;
/*<         VH2(1,1) = TANBST >*/
    vh2[0] = tanbst;
/*<         VH2(1,2) = -1. >*/
    vh2[2] = -1.f;
/*<         VH2(2,1) = -1. >*/
    vh2[1] = -1.f;
/*<         VH2(2,2) = 1./TANBST >*/
    vh2[3] = 1.f / tanbst;
/* CCCCCCCCCCCCCCCCCCCCCCCCCCC  D-terms CCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/*< 	STW=SW >*/
    stw = sw;
/*< 	F1T=( >*/
    f1t = (mul2 - mur2) / (stop12 - stop22) * (.5f - stw * 
	    1.3333333333333333f) * log(stop1 / stop2) + (.5f - stw * 
	    .66666666666666663f) * log(stop1 * stop2 / (mul2 + mtop2)) + stw *
	     .66666666666666663f * log(stop1 * stop2 / (mur2 + mtop2));
/*< 	F1B=( >*/
    f1b = (mdl2 - md2) / (sbot12 - sbot22) * (stw * .66666666666666663f - .5f)
	     * log(sbot1 / sbot2) + (stw * .33333333333333331f - .5f) * log(
	    sbot1 * sbot2 / (mdl2 + mbot2)) - stw * .33333333333333331f * log(
	    sbot1 * sbot2 / (md2 + mbot2));
/*< 	F2T=1 >*/
    f2t = 1 / (stop12 - stop22) * (log(stop12 / stop22) * -.5f + (stw * 
	    1.3333333333333333f - .5f) * (mul2 - mur2) / (stop12 - stop22) * (
	    2. - (stop12 + stop22) / (stop12 - stop22) * log(stop12 / stop22))
	    );
/*< 	F2B=1 >*/
    f2b = 1 / (sbot12 - sbot22) * (log(sbot12 / sbot22) * .5f + (stw * 
	    -.66666666666666663f + .5f) * (mdl2 - md2) / (sbot12 - sbot22) * (
	    2. - (sbot12 + sbot22) / (sbot12 - sbot22) * log(sbot12 / sbot22))
	    );
/* ************************************************************* */

/* --EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I */

/* TRAFOS APPROXIMATE -> EXACT: */

/* (i)  1/M_{SUSY}^2 -> LOG(M1^2/M2^2) / (M1^2-M2^2) */

/* (ii) 1/M_{SUSY}^4 -> -6 G(M1^2,M2^2) / (M1^2-M2^2)^2 */

/* Then use results of Phys. Lett. B355 (1995) 209 in order to */
/* obtain the results for lambda_1 - lambda_7 according to */
/* Nucl. Phys. B461 (1996) 407. Perform a full evolution from */
/* M_SUSY -> m_t for lambdas (anomalous dimensions, v_i). */

/* - ht^2*hb^2 terms neglected in lambda_3,4 (according to */
/*   Nucl. Phys. B461 (1996) 407) */

/* ************************************************************* */
/*<        >*/
/* Computing 4th power */
    d__1 = sinbt, d__1 *= d__1;
    d__3 = *mu;
    d__4 = stop1;
    d__5 = stop2;
    d__2 = d__3 * d__3 / (d__4 * d__4 - d__5 * d__5);
    d__6 = mz;
    d__7 = *mu;
    d__8 = tanbst;
    d__9 = cosbt;
    dlam1t = mtop4 / (d__1 * d__1) * (d__2 * d__2) * (2. - (stop12 + stop22) /
	     (stop12 - stop22) * log(stop12 / stop22)) - d__6 * d__6 * mtop2 *
	     (d__7 * d__7) / (d__8 * d__8) * f2t / (d__9 * d__9);
/*<        >*/
/* Computing 4th power */
    d__1 = cosbb, d__1 *= d__1;
    d__2 = sbot1;
    d__3 = sbot2;
    d__4 = *ab;
    d__5 = sbot1;
    d__6 = sbot2;
    d__7 = sbot1;
    d__8 = sbot2;
/* Computing 4th power */
    d__9 = cosbb, d__9 *= d__9;
    d__11 = *ab;
    d__12 = sbot1;
    d__13 = sbot2;
    d__10 = d__11 * d__11 / (d__12 * d__12 - d__13 * d__13);
    d__14 = mz;
    d__15 = *ab;
    d__16 = cosbb;
    dlam1b = mbot4 / (d__1 * d__1) * (log(d__2 * d__2 * (d__3 * d__3) / (mdl2 
	    + mbot2) / (md2 + mbot2)) + d__4 * d__4 * 2 / (d__5 * d__5 - d__6 
	    * d__6) * log(d__7 * d__7 / (d__8 * d__8))) + mbot4 / (d__9 * 
	    d__9) * (d__10 * d__10) * (2. - (sbot12 + sbot22) / (sbot12 - 
	    sbot22) * log(sbot12 / sbot22)) + d__14 * d__14 * (mbot2 * 2 * 
	    f1b - mbot2 * (d__15 * d__15) * f2b) / (d__16 * d__16);
/*<        >*/
/* Computing 4th power */
    d__1 = sinbt, d__1 *= d__1;
    d__2 = stop1;
    d__3 = stop2;
    d__4 = *at;
    d__5 = stop1;
    d__6 = stop2;
    d__7 = stop1;
    d__8 = stop2;
/* Computing 4th power */
    d__9 = sinbt, d__9 *= d__9;
    d__11 = *at;
    d__12 = stop1;
    d__13 = stop2;
    d__10 = d__11 * d__11 / (d__12 * d__12 - d__13 * d__13);
    d__14 = mz;
    d__15 = *at;
    d__16 = sinbt;
    dlam2t = mtop4 / (d__1 * d__1) * (log(d__2 * d__2 * (d__3 * d__3) / (mul2 
	    + mtop2) / (mur2 + mtop2)) + d__4 * d__4 * 2 / (d__5 * d__5 - 
	    d__6 * d__6) * log(d__7 * d__7 / (d__8 * d__8))) + mtop4 / (d__9 *
	     d__9) * (d__10 * d__10) * (2. - (stop12 + stop22) / (stop12 - 
	    stop22) * log(stop12 / stop22)) + d__14 * d__14 * (mtop2 * -2 * 
	    f1t + mtop2 * (d__15 * d__15) * f2t) / (d__16 * d__16);
/*<        >*/
/* Computing 4th power */
    d__1 = cosbb, d__1 *= d__1;
/* Computing 4th power */
    d__2 = *mu, d__2 *= d__2;
    d__4 = sbot1;
    d__5 = sbot2;
    d__3 = d__4 * d__4 - d__5 * d__5;
    d__6 = mz;
    d__7 = *mu;
    d__8 = tanbsb;
    d__9 = sinbb;
    dlam2b = mbot4 / (d__1 * d__1) * (d__2 * d__2) / (d__3 * d__3) * (2. - (
	    sbot12 + sbot22) / (sbot12 - sbot22) * log(sbot12 / sbot22)) + 
	    d__6 * d__6 * mbot2 * (d__7 * d__7) * (d__8 * d__8) * f2b / (d__9 
	    * d__9);
/*<        >*/
/* Computing 4th power */
    d__1 = sinbt, d__1 *= d__1;
    d__2 = *mu;
    d__3 = stop1;
    d__4 = stop2;
    d__5 = stop1;
    d__6 = stop2;
    d__7 = *at;
    d__8 = stop1;
    d__9 = stop2;
    d__10 = mz;
    d__11 = *at;
    d__12 = *mu;
    dlam3t = mtop4 / (d__1 * d__1) * (d__2 * d__2) / (d__3 * d__3 - d__4 * 
	    d__4) * (log(d__5 * d__5 / (d__6 * d__6)) / 2. + d__7 * d__7 / (
	    d__8 * d__8 - d__9 * d__9) * (2. - (stop12 + stop22) / (stop12 - 
	    stop22) * log(stop12 / stop22))) + d__10 * d__10 * (mtop2 / 
	    tanbst * f1t - mtop2 * (d__11 * d__11 - d__12 * d__12) / tanbst / 
	    2.f * f2t) / sinbt / cosbt / 2;
/*    *  + MTOP2*MBOT2/(SINBT**2*COSBB**2)*( */
/*    *    LOG(STOP1**2*STOP2**2/(MQ2+MTOP2)/(MUR2+MTOP2)) */
/*    *  + LOG(SBOT1**2*SBOT2**2/(MQ2+MBOT2)/(MD2+MBOT2)) */
/*    *  + ((AT+AB)**2/2-MU**2)*( */
/*    *      1.D0/(STOP1**2-SBOT1**2)*LOG(STOP1**2/SBOT1**2) */
/*    *    + 1.D0/(STOP2**2-SBOT2**2)*LOG(STOP2**2/SBOT2**2)) */
/*    *  - (MU**2-AT*AB)**2*( */
/*    *    - 1.D0/(STOP1**2-SBOT1**2)**2*G(STOP12,SBOT12) */
/*    *    - 1.D0/(STOP2**2-SBOT2**2)**2*G(STOP22,SBOT22))) */
/*<        >*/
/* Computing 4th power */
    d__1 = cosbb, d__1 *= d__1;
    d__2 = *mu;
    d__3 = sbot1;
    d__4 = sbot2;
    d__5 = sbot1;
    d__6 = sbot2;
    d__7 = *ab;
    d__8 = sbot1;
    d__9 = sbot2;
    d__10 = mz;
    d__11 = *ab;
    d__12 = *mu;
    dlam3b = mbot4 / (d__1 * d__1) * (d__2 * d__2) / (d__3 * d__3 - d__4 * 
	    d__4) * (log(d__5 * d__5 / (d__6 * d__6)) / 2. + d__7 * d__7 / (
	    d__8 * d__8 - d__9 * d__9) * (2. - (sbot12 + sbot22) / (sbot12 - 
	    sbot22) * log(sbot12 / sbot22))) + d__10 * d__10 * (-mbot2 * 
	    tanbsb * f1b + mbot2 * (d__11 * d__11 - d__12 * d__12) * tanbsb / 
	    2.f * f2b) / sinbb / cosbb / 2;
/*<        >*/
/* Computing 4th power */
    d__1 = sinbt, d__1 *= d__1;
    d__2 = *mu;
    d__3 = stop1;
    d__4 = stop2;
    d__5 = stop1;
    d__6 = stop2;
    d__7 = *at;
    d__8 = stop1;
    d__9 = stop2;
    d__10 = mz;
    d__11 = *at;
    d__12 = *mu;
    dlam4t = mtop4 / (d__1 * d__1) * (d__2 * d__2) / (d__3 * d__3 - d__4 * 
	    d__4) * (log(d__5 * d__5 / (d__6 * d__6)) / 2. + d__7 * d__7 / (
	    d__8 * d__8 - d__9 * d__9) * (2. - (stop12 + stop22) / (stop12 - 
	    stop22) * log(stop12 / stop22))) + d__10 * d__10 * (mtop2 / 
	    tanbst * f1t - mtop2 * (d__11 * d__11 - d__12 * d__12) / tanbst / 
	    2.f * f2t) / sinbt / cosbt / 2;
/*    *  - MTOP2*MBOT2/(SINBT**2*COSBB**2)*( */
/*    *    LOG(STOP1**2*STOP2**2/(MQ2+MTOP2)/(MUR2+MTOP2)) */
/*    *  + LOG(SBOT1**2*SBOT2**2/(MQ2+MBOT2)/(MD2+MBOT2)) */
/*    *  + ((AT+AB)**2/2-MU**2)*( */
/*    *      1.D0/(STOP1**2-SBOT1**2)*LOG(STOP1**2/SBOT1**2) */
/*    *    + 1.D0/(STOP2**2-SBOT2**2)*LOG(STOP2**2/SBOT2**2)) */
/*    *  - (MU**2-AT*AB)**2*( */
/*    *    - 1.D0/(STOP1**2-SBOT1**2)**2*G(STOP12,SBOT12) */
/*    *    - 1.D0/(STOP2**2-SBOT2**2)**2*G(STOP22,SBOT22))) */
/*<        >*/
/* Computing 4th power */
    d__1 = cosbb, d__1 *= d__1;
    d__2 = *mu;
    d__3 = sbot1;
    d__4 = sbot2;
    d__5 = sbot1;
    d__6 = sbot2;
    d__7 = *ab;
    d__8 = sbot1;
    d__9 = sbot2;
    d__10 = mz;
    d__11 = *ab;
    d__12 = *mu;
    dlam4b = mbot4 / (d__1 * d__1) * (d__2 * d__2) / (d__3 * d__3 - d__4 * 
	    d__4) * (log(d__5 * d__5 / (d__6 * d__6)) / 2. + d__7 * d__7 / (
	    d__8 * d__8 - d__9 * d__9) * (2. - (sbot12 + sbot22) / (sbot12 - 
	    sbot22) * log(sbot12 / sbot22))) + d__10 * d__10 * (-mbot2 * 
	    tanbsb * f1b + mbot2 * (d__11 * d__11 - d__12 * d__12) * tanbsb / 
	    2.f * f2b) / sinbb / cosbb / 2;
/*<        >*/
/* Computing 4th power */
    d__1 = sinbt, d__1 *= d__1;
    d__2 = *mu;
    d__3 = *at;
    d__5 = stop1;
    d__6 = stop2;
    d__4 = d__5 * d__5 - d__6 * d__6;
    dlam5t = mtop4 / (d__1 * d__1) * (d__2 * d__2 * (d__3 * d__3)) / (d__4 * 
	    d__4) * (2. - (stop12 + stop22) / (stop12 - stop22) * log(stop12 /
	     stop22));
/*<        >*/
/* Computing 4th power */
    d__1 = cosbb, d__1 *= d__1;
    d__2 = *mu;
    d__3 = *ab;
    d__5 = sbot1;
    d__6 = sbot2;
    d__4 = d__5 * d__5 - d__6 * d__6;
    dlam5b = mbot4 / (d__1 * d__1) * (d__2 * d__2 * (d__3 * d__3)) / (d__4 * 
	    d__4) * (2. - (sbot12 + sbot22) / (sbot12 - sbot22) * log(sbot12 /
	     sbot22));
/*<        >*/
/* Computing 4th power */
    d__1 = sinbt, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = *mu;
    d__4 = stop1;
    d__5 = stop2;
    d__3 = d__4 * d__4 - d__5 * d__5;
    d__6 = mz;
    dlam6t = mtop4 / (d__1 * d__1) * (-(d__2 * (d__2 * d__2)) * *at) / (d__3 *
	     d__3) * (2. - (stop12 + stop22) / (stop12 - stop22) * log(stop12 
	    / stop22)) + d__6 * d__6 * mtop2 * *mu * *at / tanbst * f2t / (
	    sinbt * 2 * cosbt);
/*<        >*/
/* Computing 4th power */
    d__1 = cosbb, d__1 *= d__1;
    d__2 = sbot1;
    d__3 = sbot2;
    d__4 = sbot1;
    d__5 = sbot2;
    d__6 = *ab;
    d__8 = sbot1;
    d__9 = sbot2;
    d__7 = d__8 * d__8 - d__9 * d__9;
    d__10 = mz;
    dlam6b = mbot4 / (d__1 * d__1) * *mu * *ab * (-1. / (d__2 * d__2 - d__3 * 
	    d__3) * log(d__4 * d__4 / (d__5 * d__5)) - d__6 * d__6 / (d__7 * 
	    d__7) * (2. - (sbot12 + sbot22) / (sbot12 - sbot22) * log(sbot12 /
	     sbot22))) - d__10 * d__10 * (-mbot2 * *ab * *mu * tanbsb * f2b) /
	     (sinbb * 2 * cosbb);
/*<        >*/
/* Computing 4th power */
    d__1 = sinbt, d__1 *= d__1;
    d__2 = stop1;
    d__3 = stop2;
    d__4 = stop1;
    d__5 = stop2;
    d__6 = *at;
    d__8 = stop1;
    d__9 = stop2;
    d__7 = d__8 * d__8 - d__9 * d__9;
    d__10 = mz;
    dlam7t = mtop4 / (d__1 * d__1) * *mu * *at * (-1. / (d__2 * d__2 - d__3 * 
	    d__3) * log(d__4 * d__4 / (d__5 * d__5)) - d__6 * d__6 / (d__7 * 
	    d__7) * (2. - (stop12 + stop22) / (stop12 - stop22) * log(stop12 /
	     stop22))) - d__10 * d__10 * mtop2 * *at * *mu / tanbst * f2t / (
	    sinbt * 2 * cosbt);
/*<        >*/
/* Computing 4th power */
    d__1 = cosbb, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = *mu;
    d__4 = sbot1;
    d__5 = sbot2;
    d__3 = d__4 * d__4 - d__5 * d__5;
    d__6 = mz;
    dlam7b = mbot4 / (d__1 * d__1) * (-(d__2 * (d__2 * d__2)) * *ab) / (d__3 *
	     d__3) * (2. - (sbot12 + sbot22) / (sbot12 - sbot22) * log(sbot12 
	    / sbot22)) - d__6 * d__6 * mbot2 * *mu * *ab * tanbsb * f2b / (
	    sinbb * 2 * cosbb);
/*<        TQ = LOG((MUL2 + MTOP2)/MTOP2) >*/
    tq = log((mul2 + mtop2) / mtop2);
/*<        TU = LOG((MUR2+MTOP2)/MTOP2) >*/
    tu = log((mur2 + mtop2) / mtop2);
/*<        TQD = LOG((MDL2 + MB**2)/MB**2) >*/
    d__1 = mb;
    d__2 = mb;
    tqd = log((mdl2 + d__1 * d__1) / (d__2 * d__2));
/*<        TD = LOG((MD2+MB**2)/MB**2) >*/
    d__1 = mb;
    d__2 = mb;
    td = log((md2 + d__1 * d__1) / (d__2 * d__2));
/*<         FACT = 3.D0/(16.D0*PI**2*(H1T**2+H2T**2)**2) >*/
    d__1 = pi;
    d__3 = h1t;
    d__4 = h2t;
    d__2 = d__3 * d__3 + d__4 * d__4;
    fact = 3. / (d__1 * d__1 * 16. * (d__2 * d__2));
/*<         FACB = 3.D0/(16.D0*PI**2*(H1B**2+H2B**2)**2) >*/
    d__1 = pi;
    d__3 = h1b;
    d__4 = h2b;
    d__2 = d__3 * d__3 + d__4 * d__4;
    facb = 3. / (d__1 * d__1 * 16. * (d__2 * d__2));
/*<         DLAM1 = FACT*DLAM1T*(1.-AL1*TT) + FACB*DLAM1B*(1.-AL1*TB) >*/
    *dlam1 = fact * dlam1t * (1.f - al1 * tt) + facb * dlam1b * (1.f - al1 * 
	    tb);
/*<         DLAM2 = FACT*DLAM2T*(1.-AL2*TT) + FACB*DLAM2B*(1.-AL2*TB) >*/
    *dlam2 = fact * dlam2t * (1.f - al2 * tt) + facb * dlam2b * (1.f - al2 * 
	    tb);
/*<        >*/
    *dlam3 = fact * dlam3t * (1.f - (al1 + al2) / 2 * tt) + facb * dlam3b * (
	    1.f - (al1 + al2) / 2 * tb);
/*<        >*/
    *dlam4 = fact * dlam4t * (1.f - (al1 + al2) / 2 * tt) + facb * dlam4b * (
	    1.f - (al1 + al2) / 2 * tb);
/*<        >*/
    *dlam5 = fact * dlam5t * (1.f - (al1 + al2) / 2 * tt) + facb * dlam5b * (
	    1.f - (al1 + al2) / 2 * tb);
/*<        >*/
    *dlam6 = fact * dlam6t * (1.f - (al1 * 3 + al2) / 4 * tt) + facb * dlam6b 
	    * (1.f - (al1 * 3 + al2) / 4 * tb);
/*<        >*/
    *dlam7 = fact * dlam7t * (1.f - (al1 + al2 * 3) / 4 * tt) + facb * dlam7b 
	    * (1.f - (al1 + al2 * 3) / 4 * tb);
/*<         FACTOR = 1.D0 >*/
    factor = 1.;
/*<         DLAM1 = DLAM1 * FACTOR >*/
    *dlam1 *= factor;
/*<         DLAM2 = DLAM2 * FACTOR >*/
    *dlam2 *= factor;
/*<         DLAM3 = DLAM3 * FACTOR >*/
    *dlam3 *= factor;
/*<         DLAM4 = DLAM4 * FACTOR >*/
    *dlam4 *= factor;
/*<         DLAM5 = DLAM5 * FACTOR >*/
    *dlam5 *= factor;
/*<         DLAM6 = DLAM6 * FACTOR >*/
    *dlam6 *= factor;
/*<         DLAM7 = DLAM7 * FACTOR >*/
    *dlam7 *= factor;
/* --END OF EXTENSION */
/*<         GOTO 4236 >*/
    goto L4236;
/*<  4237   CONTINUE >*/
L4237:
/*<         DLAM1 = -1.D+15 >*/
    *dlam1 = 0;
/*<         DLAM2 = -1.D+15 >*/
    *dlam2 = 0;
/*<         DLAM3 = -1.D+15 >*/
    *dlam3 = 0;
/*<         DLAM4 = -1.D+15 >*/
    *dlam4 = 0;
/*<         DLAM5 = -1.D+15 >*/
    *dlam5 = 0;
/*<         DLAM6 = -1.D+15 >*/
    *dlam6 = 0;
/*<         DLAM7 = -1.D+15 >*/
    *dlam7 = 0;
/*< 4236    RETURN >*/
L4236:
    return 0;
/*<         END >*/
} /* gfun_hdec__ */

/*<        >*/
/* Subroutine */ int delmb_hdec__(double *ma, double *tanb, 
	double *mul, double *mdl, double *mur, double *md, 
	double *at, double *ab, double *mu, double *mglu, 
	double *mtop, double *deltamt, double *deltamb, 
	double *stop12, double *stop22, double *sbot12, 
	double *sbot22)
{
    /* System generated locals */
    double d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11;

    /* Builtin functions */
    double sqrt(double), atan(double), log(double), pow_dd(
	    double *, double *);

    /* Local variables */
    static double v, g1, g2, g3, hb, g32, mb, tb, pi, cw, ht, vi, mw, mz, 
	    tt, sw, bb2, h1b, h2b, al1, al2, h1i, md2, h2i, bt2, h1t, h2t, 
	    msb, mst, mdl2, mul2, mur2, cosb, sinb, alst, htst, mbot2, mbot4, 
	    sbot1, sbot2, bt2st, mtop2, mtop4, stop1, stop2, tanba, cosba, 
	    cosbb, sinba, sinbb, cosbt, sinbt;
    static double rmtop;
    static double alpha1, alpha2, alpha3;
    extern double t_hdec__(double *, double *, double *);
    static double tanbsb, tanbst, msusyb, alpha3z, msusyt;

/*<         IMPLICIT REAL*8 (A-H,L,M,O-Z) >*/
/*<         real*8 MbRun,MtRun >*/
/*<         COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW >*/
/*<         IF(DABS(MU).LT.0.000001) MU = 0.000001 >*/
    if (abs(*mu) < 1e-6f) {
	*mu = 1e-6f;
    }
/*<         MUL2  = MUL**2 >*/
    d__1 = *mul;
    mul2 = d__1 * d__1;
/*<         MDL2  = MDL**2 >*/
    d__1 = *mdl;
    mdl2 = d__1 * d__1;
/*<         MUR2  = MUR**2 >*/
    d__1 = *mur;
    mur2 = d__1 * d__1;
/*<         MD2   = MD**2 >*/
    d__1 = *md;
    md2 = d__1 * d__1;
/*<         TANBA = TANB >*/
    tanba = *tanb;
/*<         SINBA = TANBA/DSQRT(TANBA**2+1.D0) >*/
    d__1 = tanba;
    sinba = tanba / sqrt(d__1 * d__1 + 1.);
/*<         COSBA = SINBA/TANBA         >*/
    cosba = sinba / tanba;
/*<         SINB = TANB/DSQRT(TANB**2+1.D0) >*/
    d__1 = *tanb;
    sinb = *tanb / sqrt(d__1 * d__1 + 1.);
/*<         COSB = SINB/TANB >*/
    cosb = sinb / *tanb;
/*<       RMTOP = MtRun(MTOP) >*/
    rmtop = MtRun(*mtop);
/*<       MB = MbRun(MTOP) >*/
    mb = MbRun(*mtop);
/*<       PI = 4*DATAN(1D0) >*/
    pi = atan(1.) * 4;
/*<       MZ = AMZ >*/
    mz = param_hdec__1.amz;
/*<       MW = AMW >*/
    mw = param_hdec__1.amw;
/*<       V  = 1/DSQRT(2*DSQRT(2D0)*GF) >*/
    v = 1 / sqrt(sqrt(2.) * 2 * param_hdec__1.gf);
/*<       CW = AMW**2/AMZ**2 >*/
    d__1 = param_hdec__1.amw;
    d__2 = param_hdec__1.amz;
    cw = d__1 * d__1 / (d__2 * d__2);
/*<       SW = 1-CW >*/
    sw = 1 - cw;
/*<       ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI >*/
    d__1 = param_hdec__1.amw * 2 / v / sqrt(2.);
    alpha2 = d__1 * d__1 / 4 / pi;
/*<       ALPHA1  = ALPHA2*SW/CW >*/
    alpha1 = alpha2 * sw / cw;
/*<       ALPHA3Z = alphaQCD(AMZ) >*/
    alpha3z = alphaQCD(param_hdec__1.amz);
/*<       ALPHA3  = alphaQCD(MTOP) >*/
    alpha3 = alphaQCD(*mtop);
/*<       G1 = DSQRT(ALPHA1*4.*PI) >*/
    g1 = sqrt(alpha1 * 4.f * pi);
/*<       G2 = DSQRT(ALPHA2*4.*PI) >*/
    g2 = sqrt(alpha2 * 4.f * pi);
/*<       G3 = DSQRT(ALPHA3*4.*PI) >*/
    g3 = sqrt(alpha3 * 4.f * pi);
/*<         IF(MUL.GT.MUR) MST = MUL >*/
    if (*mul > *mur) {
	mst = *mul;
    }
/*<         IF(MUR.GT.MUL.OR.MUR.EQ.MUL) MST = MUR >*/
    if (*mur > *mul || *mur == *mul) {
	mst = *mur;
    }
/*<         MSUSYT = DSQRT(MST**2  + MTOP**2) >*/
    d__1 = mst;
    d__2 = *mtop;
    msusyt = sqrt(d__1 * d__1 + d__2 * d__2);
/*< 	IF(MDL.GT.MD) MSB = MDL >*/
    if (*mdl > *md) {
	msb = *mdl;
    }
/*< 	IF(MD.GT.MDL.OR.MD.EQ.MDL) MSB = MD >*/
    if (*md > *mdl || *md == *mdl) {
	msb = *md;
    }
/*< 	MSUSYB = DSQRT(MSB**2 + MB**2) >*/
    d__1 = msb;
    d__2 = mb;
    msusyb = sqrt(d__1 * d__1 + d__2 * d__2);
/*< 	TT = LOG(MSUSYT**2/MTOP**2) >*/
    d__1 = msusyt;
    d__2 = *mtop;
    tt = log(d__1 * d__1 / (d__2 * d__2));
/*< 	TB = LOG(MSUSYB**2/MTOP**2) >*/
    d__1 = msusyb;
    d__2 = *mtop;
    tb = log(d__1 * d__1 / (d__2 * d__2));
/*<         HT = RMTOP/V/SINB >*/
    ht = rmtop / v / sinb;
/*<         HTST = RMTOP/V >*/
    htst = rmtop / v;
/*<         HB =  MB/V/COSB >*/
    hb = mb / v / cosb;
/*<         G32 = ALPHA3*4.*PI >*/
    g32 = alpha3 * 4.f * pi;
/*<         BT2 = -(8.*G32 - 9.*HT**2/2. - HB**2/2.)/(4.*PI)**2 >*/
    d__1 = ht;
    d__2 = hb;
    d__3 = pi * 4.f;
    bt2 = -(g32 * 8.f - d__1 * d__1 * 9.f / 2.f - d__2 * d__2 / 2.f) / (d__3 *
	     d__3);
/*< 	BB2 = -(8.*G32 - 9.*HB**2/2. - HT**2/2.)/(4.*PI)**2 >*/
    d__1 = hb;
    d__2 = ht;
    d__3 = pi * 4.f;
    bb2 = -(g32 * 8.f - d__1 * d__1 * 9.f / 2.f - d__2 * d__2 / 2.f) / (d__3 *
	     d__3);
/*<         AL2 = 3./8./PI**2*HT**2 >*/
    d__1 = pi;
    d__2 = ht;
    al2 = .375f / (d__1 * d__1) * (d__2 * d__2);
/*<         BT2ST = -(8.*G32 - 9.*HTST**2/2.)/(4.*PI)**2 >*/
    d__1 = htst;
    d__2 = pi * 4.f;
    bt2st = -(g32 * 8.f - d__1 * d__1 * 9.f / 2.f) / (d__2 * d__2);
/*<         ALST = 3./8./PI**2*HTST**2 >*/
    d__1 = pi;
    d__2 = htst;
    alst = .375f / (d__1 * d__1) * (d__2 * d__2);
/*<         AL1 = 3./8./PI**2*HB**2 >*/
    d__1 = pi;
    d__2 = hb;
    al1 = .375f / (d__1 * d__1) * (d__2 * d__2);
/*<         IF(MA.GT.MTOP) THEN >*/
    if (*ma > *mtop) {
/*<         VI = V*(1. + 3./32./PI**2*HTST**2*LOG(MTOP**2/MA**2)) >*/
	d__1 = pi;
	d__2 = htst;
	d__3 = *mtop;
	d__4 = *ma;
	vi = v * (.09375f / (d__1 * d__1) * (d__2 * d__2) * log(d__3 * d__3 / 
		(d__4 * d__4)) + 1.f);
/*<         H1I = VI*COSBA >*/
	h1i = vi * cosba;
/*<         H2I = VI*SINBA >*/
	h2i = vi * sinba;
/*<         H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYT**2))**.25 >*/
	d__2 = pi;
	d__3 = hb;
	d__4 = *ma;
	d__5 = msusyt;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h1t = h1i * pow(d__1, c_b90);
/*<         H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYT**2))**.25 >*/
	d__2 = pi;
	d__3 = ht;
	d__4 = *ma;
	d__5 = msusyt;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h2t = h2i * pow(d__1, c_b90);
/*<         H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYB**2))**.25 >*/
	d__2 = pi;
	d__3 = hb;
	d__4 = *ma;
	d__5 = msusyb;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h1b = h1i * pow(d__1, c_b90);
/*<         H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYB**2))**.25 >*/
	d__2 = pi;
	d__3 = ht;
	d__4 = *ma;
	d__5 = msusyb;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h2b = h2i * pow(d__1, c_b90);
/*<         ELSE >*/
    } else {
/*<         VI =  V >*/
	vi = v;
/*<         H1I = VI*COSB >*/
	h1i = vi * cosb;
/*<         H2I = VI*SINB >*/
	h2i = vi * sinb;
/*<         H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYT**2))**.25 >*/
	d__2 = pi;
	d__3 = hb;
	d__4 = *mtop;
	d__5 = msusyt;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h1t = h1i * pow(d__1, c_b90);
/*<         H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYT**2))**.25 >*/
	d__2 = pi;
	d__3 = ht;
	d__4 = *mtop;
	d__5 = msusyt;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h2t = h2i * pow(d__1, c_b90);
/*<         H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYB**2))**.25 >*/
	d__2 = pi;
	d__3 = hb;
	d__4 = *mtop;
	d__5 = msusyb;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h1b = h1i * pow(d__1, c_b90);
/*<         H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYB**2))**.25 >*/
	d__2 = pi;
	d__3 = ht;
	d__4 = *mtop;
	d__5 = msusyb;
	d__1 = .375f / (d__2 * d__2) * (d__3 * d__3) * log(d__4 * d__4 / (
		d__5 * d__5)) + 1.f;
	h2b = h2i * pow(d__1, c_b90);
/*<         END IF >*/
    }
/*<         TANBST = H2T/H1T >*/
    tanbst = h2t / h1t;
/*<         SINBT = TANBST/(1.+TANBST**2)**.5 >*/
    d__2 = tanbst;
    d__1 = d__2 * d__2 + 1.f;
    sinbt = tanbst / pow(d__1, c_b98);
/*<         COSBT = SINBT/TANBST >*/
    cosbt = sinbt / tanbst;
/*<         TANBSB = H2B/H1B >*/
    tanbsb = h2b / h1b;
/*<         SINBB = TANBSB/(1.+TANBSB**2)**.5 >*/
    d__2 = tanbsb;
    d__1 = d__2 * d__2 + 1.f;
    sinbb = tanbsb / pow(d__1, c_b98);
/*<         COSBB = SINBB/TANBSB >*/
    cosbb = sinbb / tanbsb;
/*<         deltamt = 0 >*/
    *deltamt = 0.;
/*<         deltamb = 0 >*/
    *deltamb = 0.;
/*<         mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt) >*/
    mtop4 = pow(rmtop, c_b83) * (bt2 * 2.f * tt + 1.f - al2 * tt - *
	    deltamt * 4.f);
/*     * /(1.+deltamt)**4. */
/*<        >*/
    d__1 = *deltamb + 1.f;
    mbot4 = pow(mb, c_b83) * (bb2 * 2.f * tb + 1.f - al1 * tb) / pow(
	    d__1, c_b83);
/*<         MTOP2 = DSQRT(MTOP4) >*/
    mtop2 = sqrt(mtop4);
/*< 	MBOT2 = DSQRT(MBOT4) >*/
    mbot2 = sqrt(sqrt(mbot4));  // fabs  by AP
/*<        >*/
    d__1 = g2;
    d__2 = g1;
    d__3 = h1t;
    d__4 = h2t;
    d__7 = g2;
    d__8 = g1;
    d__9 = h1t;
    d__10 = h2t;
    d__6 = (d__7 * d__7 - d__8 * d__8 * 5.f / 3.f) / 4.f * (d__9 * d__9 - 
	    d__10 * d__10) + mul2 - mur2;
    d__11 = *at - *mu / tanbst;
    d__5 = d__6 * d__6 * .25f + mtop2 * (d__11 * d__11);
    *stop12 = (mul2 + mur2) * .5f + mtop2 + (d__1 * d__1 + d__2 * d__2) * 
	    .125f * (d__3 * d__3 - d__4 * d__4) + pow(d__5, c_b98);
/*<        >*/
    d__1 = g2;
    d__2 = g1;
    d__3 = h1t;
    d__4 = h2t;
    d__7 = g2;
    d__8 = g1;
    d__9 = h1t;
    d__10 = h2t;
    d__6 = (d__7 * d__7 - d__8 * d__8 * 5.f / 3.f) / 4.f * (d__9 * d__9 - 
	    d__10 * d__10) + mul2 - mur2;
    d__11 = *at - *mu / tanbst;
    d__5 = d__6 * d__6 * .25f + mtop2 * (d__11 * d__11);
    *stop22 = (mul2 + mur2) * .5f + mtop2 + (d__1 * d__1 + d__2 * d__2) * 
	    .125f * (d__3 * d__3 - d__4 * d__4) - pow(d__5, c_b98);
/*<         IF(STOP22.LT.0.) GOTO 4237 >*/
    if (*stop22 < 0.f) {
	goto L4237;
    }
/*<        >*/
    d__1 = g2;
    d__2 = g1;
    d__3 = h1b;
    d__4 = h2b;
    d__7 = g1;
    d__8 = g2;
    d__9 = h1b;
    d__10 = h2b;
    d__6 = (d__7 * d__7 / 3.f - d__8 * d__8) / 4.f * (d__9 * d__9 - d__10 * 
	    d__10) + mdl2 - md2;
    d__11 = *ab - *mu * tanbsb;
    d__5 = d__6 * d__6 * .25f + mbot2 * (d__11 * d__11);
    *sbot12 = (mdl2 + md2) * .5f - (d__1 * d__1 + d__2 * d__2) * .125f * (
	    d__3 * d__3 - d__4 * d__4) + pow(d__5, c_b98);
/*<        >*/
    d__1 = g2;
    d__2 = g1;
    d__3 = h1b;
    d__4 = h2b;
    d__7 = g1;
    d__8 = g2;
    d__9 = h1b;
    d__10 = h2b;
    d__6 = (d__7 * d__7 / 3.f - d__8 * d__8) / 4.f * (d__9 * d__9 - d__10 * 
	    d__10) + mdl2 - md2;
    d__11 = *ab - *mu * tanbsb;
    d__5 = d__6 * d__6 * .25f + mbot2 * (d__11 * d__11);
    *sbot22 = (mdl2 + md2) * .5f - (d__1 * d__1 + d__2 * d__2) * .125f * (
	    d__3 * d__3 - d__4 * d__4) - pow(d__5, c_b98);
/*<         IF(SBOT22.LT.0.) GOTO 4237 >*/
    if (*sbot22 < 0.f) {
	goto L4237;
    }
/*<         STOP1 = STOP12**.5 >*/
    stop1 = pow(*stop12, c_b98);
/*<         STOP2 = STOP22**.5 >*/
    stop2 = pow(*stop22, c_b98);
/*<         SBOT1 = SBOT12**.5 >*/
    sbot1 = pow(*sbot12, c_b98);
/*<         SBOT2 = SBOT22**.5 >*/
    sbot2 = pow(*sbot22, c_b98);
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*     Here is the definition of deltamb and deltamt, which */
/*     are the vertex corrections to the bottom and top quark */
/*     mass, keeping the dominant QCD and top Yukawa coupling */
/*     induced corrections. */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*<        >*/
    d__1 = ht;
    d__2 = pi * 4.f;
    *deltamb = alpha3 * -2 / 3.f / pi * *mglu * (*ab - *mu * *tanb) * 
	    t_hdec__(&sbot1, &sbot2, mglu) + d__1 * d__1 / (d__2 * d__2) * (*
	    at - *mu / *tanb) * *mu * *tanb * t_hdec__(&stop1, &stop2, mu);
/*<        >*/
    *deltamt = alpha3 * -2.f / 3.f / pi * (*at - *mu / *tanb) * *mglu * 
	    t_hdec__(&stop1, &stop2, mglu);
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*   Here the new values of the top and bottom quark masses at */
/*   the scale MS are defined, to be used in the effective */
/*   potential approximation. They are just the old ones, but */
/*   including the finite corrections deltamt and deltamb. */
/*   The deltamb corrections can become large and are resummed */
/*   to all orders, as suggested in the two recent works by M. Carena, */
/*   S. Mrenna and C.E.M. Wagner, as well as in the work by M. Carena, */
/*   D. Garcia, U. Nierste and C.E.M. Wagner, to appear. The top */
/*   quark mass corrections are small and are kept in the perturbative */
/*   formulation. The function T(X,Y,Z) is necessary for the calculation. */
/*   the entries are masses and NOT their squares ! */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*<         mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt) >*/
    mtop4 = pow(rmtop, c_b83) * (bt2 * 2.f * tt + 1.f - al2 * tt - *
	    deltamt * 4.f);
/*     * /(1.+deltamt)**4. */
/*<        >*/
    d__1 = *deltamb + 1.f;
    mbot4 = pow(mb, c_b83) * (bb2 * 2.f * tb + 1.f - al1 * tb) / pow(
	    d__1, c_b83);
/*<         MTOP2 = DSQRT(MTOP4) >*/
    mtop2 = sqrt(mtop4);
/*< 	MBOT2 = DSQRT(MBOT4) >*/
    mbot2 = sqrt(fabs(mbot4));   // fabs  by AP
/*<        >*/
    d__1 = g2;
    d__2 = g1;
    d__3 = h1t;
    d__4 = h2t;
    d__7 = g2;
    d__8 = g1;
    d__9 = h1t;
    d__10 = h2t;
    d__6 = (d__7 * d__7 - d__8 * d__8 * 5.f / 3.f) / 4.f * (d__9 * d__9 - 
	    d__10 * d__10) + mul2 - mur2;
    d__11 = *at - *mu / tanbst;
    d__5 = d__6 * d__6 * .25f + mtop2 * (d__11 * d__11);
    *stop12 = (mul2 + mur2) * .5f + mtop2 + (d__1 * d__1 + d__2 * d__2) * 
	    .125f * (d__3 * d__3 - d__4 * d__4) + pow(d__5, c_b98);
/*<        >*/
    d__1 = g2;
    d__2 = g1;
    d__3 = h1t;
    d__4 = h2t;
    d__7 = g2;
    d__8 = g1;
    d__9 = h1t;
    d__10 = h2t;
    d__6 = (d__7 * d__7 - d__8 * d__8 * 5.f / 3.f) / 4.f * (d__9 * d__9 - 
	    d__10 * d__10) + mul2 - mur2;
    d__11 = *at - *mu / tanbst;
    d__5 = d__6 * d__6 * .25f + mtop2 * (d__11 * d__11);
    *stop22 = (mul2 + mur2) * .5f + mtop2 + (d__1 * d__1 + d__2 * d__2) * 
	    .125f * (d__3 * d__3 - d__4 * d__4) - pow(d__5, c_b98);
/*<         IF(STOP22.LT.0.) GOTO 4237 >*/
    if (*stop22 < 0.f) {
	goto L4237;
    }
/*<        >*/
    d__1 = g2;
    d__2 = g1;
    d__3 = h1b;
    d__4 = h2b;
    d__7 = g1;
    d__8 = g2;
    d__9 = h1b;
    d__10 = h2b;
    d__6 = (d__7 * d__7 / 3.f - d__8 * d__8) / 4.f * (d__9 * d__9 - d__10 * 
	    d__10) + mdl2 - md2;
    d__11 = *ab - *mu * tanbsb;
    d__5 = d__6 * d__6 * .25f + mbot2 * (d__11 * d__11);
    *sbot12 = (mdl2 + md2) * .5f - (d__1 * d__1 + d__2 * d__2) * .125f * (
	    d__3 * d__3 - d__4 * d__4) + pow(d__5, c_b98);
/*<        >*/
    d__1 = g2;
    d__2 = g1;
    d__3 = h1b;
    d__4 = h2b;
    d__7 = g1;
    d__8 = g2;
    d__9 = h1b;
    d__10 = h2b;
    d__6 = (d__7 * d__7 / 3.f - d__8 * d__8) / 4.f * (d__9 * d__9 - d__10 * 
	    d__10) + mdl2 - md2;
    d__11 = *ab - *mu * tanbsb;
    d__5 = d__6 * d__6 * .25f + mbot2 * (d__11 * d__11);
    *sbot22 = (mdl2 + md2) * .5f - (d__1 * d__1 + d__2 * d__2) * .125f * (
	    d__3 * d__3 - d__4 * d__4) - pow(d__5, c_b98);
/*< 4237    RETURN >*/
L4237:
    return 0;
/*<         END >*/
} /* delmb_hdec__ */

/*<       FUNCTION T_HDEC(X,Y,Z) >*/
double t_hdec__(double *x, double *y, double *z__)
{
    /* System generated locals */
    double ret_val, d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, 
	    d__10, d__11, d__12, d__13, d__14, d__15, d__16, d__17, d__18;

    /* Builtin functions */
    double log(double);

/*<       implicit real*8(a-h,l,m,o-z) >*/
/*<       if(x.eq.y) x = x - 0.00001 >*/
    if (*x == *y) {
	*x += -1e-5f;
    }
/*<       if(x.eq.z) x = x - 0.00002 >*/
    if (*x == *z__) {
	*x += -2e-5f;
    }
/*<       if(y.eq.z) y = y - 0.00003 >*/
    if (*y == *z__) {
	*y += -3e-5f;
    }
/*       write(*,*) 'xyz',x,y,z */
/*<        >*/
    d__1 = *x;
    d__2 = *y;
    d__3 = *x;
    d__4 = *y;
    d__5 = *x;
    d__6 = *z__;
    d__7 = *z__;
    d__8 = *x;
    d__9 = *y;
    d__10 = *z__;
    d__11 = *y;
    d__12 = *z__;
    d__13 = *x;
    d__14 = *y;
    d__15 = *y;
    d__16 = *z__;
    d__17 = *x;
    d__18 = *z__;
    ret_val = (d__1 * d__1 * (d__2 * d__2) * log(d__3 * d__3 / (d__4 * d__4)) 
	    + d__5 * d__5 * (d__6 * d__6) * log(d__7 * d__7 / (d__8 * d__8)) 
	    + d__9 * d__9 * (d__10 * d__10) * log(d__11 * d__11 / (d__12 * 
	    d__12))) / ((d__13 * d__13 - d__14 * d__14) * (d__15 * d__15 - 
	    d__16 * d__16) * (d__17 * d__17 - d__18 * d__18));
/*<       return >*/
    return ret_val;
/*<       end >*/
} /* t_hdec__ */

