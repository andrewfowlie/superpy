#include <math.h>
#include <stdio.h>

#define abs(x) fabs(x)	
		
static double mysqrt(double x)
{
	if(x<0)
		return -sqrt(-x);
	else
		return sqrt(x);
}

static int dia_(int *, double *, double *, double *);

static double neuz[32][32], neur[32];
	
extern double neut1_(double *MZ, double *MW, double *SW, double *MG1, 
	double *MG2, double *vv, double *tb, double *hk, double *hl, double *xv);

double neu5_(double *m1, double *m2, double *sw, double *mw, double *ee, 
	double *sb, double *cb, double *hk, double *hl, double *xvev)
{
	double neum[32][32];
	int n=5;
	int i, j, f=1;
	double vev=sqrt(2.0)*(*mw)*(*sw)/(*ee), mz=(*mw)/sqrt(1.0-(*sw)*(*sw)), 
	tb=(*sb)/(*cb);

	neut1_(&mz,mw,sw,m1,m2,&vev,&tb,hk,hl,xvev);
	return 1.0;
	
	neum[0][0]=(*m1);
	neum[1][1]=(*m2);
	neum[2][2]=neum[3][3]=0.0;
	neum[4][4]=2.0*(*hk)*(*xvev);
	
	neum[0][1]=neum[1][0]=0.0;
	neum[0][2]=neum[2][0]=-(*mw)/sqrt(1.0-(*sw)*(*sw))*(*sw)*(*cb);
	neum[0][3]=neum[3][0]= (*mw)/sqrt(1.0-(*sw)*(*sw))*(*sw)*(*sb);
	neum[0][4]=neum[4][0]=0.0;
	
	neum[1][2]=neum[2][1]= (*mw)*(*cb);
	neum[1][3]=neum[3][1]=-(*mw)*(*sb);
	neum[1][4]=neum[4][1]=0;
				
	neum[2][3]=neum[3][2]=-(*hl)*(*xvev);
	neum[2][4]=neum[4][2]=-vev*(*hl)*(*sb);
	
	neum[3][4]=neum[4][3]=-vev*(*hl)*(*cb);

	
	if(dia_(&n,(double *)neum,(double *)neur, (double *)neuz))
		return 0;
	
	while(f)
	{
		f=0;
		for(i=0;i<4;i++)
			if(fabs(neur[i])>fabs(neur[i+1]))
			{
				double tmp;
				f=1;
				tmp=neur[i];
				neur[i]=neur[i+1];
				neur[i+1]=tmp;
				for(j=0;j<5;j++)
				{
					tmp=neuz[i][j];
					neuz[i][j]=neuz[i+1][j];
					neuz[i+1][j]=tmp;
				}
			}
	}
	
	return 1.0;
}


double neudiag5_(double *di, double *dj, double *tk)
{
	int i,j;
	
	if((*tk)<0.5)
		return 0.0;
	i=floor((*di)+0.5);
	j=floor((*dj)+0.5);
	
	if(i==0)
		return neur[j-1];
	else
		return neuz[j-1][i-1];
}

static double hisz[32][32], hisr[32];


double higs1_(double *swf, double *mwf, double *eef, double *sbf, 
	double *cbf, double *hkf, double *hlf, double *xvevf, double *hksf, 
	double *hlsf, double *dlh2f)
{
	double sw; double mw; double ee; double sb; 
	double cb; double hk; double hl; double xvev; double hks; double hls;
	double dlh2;
	double neum[32][32];
	int n=3;
	int i, j, f=1;
	double cw;
	double vev;
	
	sw= *swf;
	mw= *mwf;
	ee= *eef;
	sb= *sbf;
	cb= *cbf;
	hk= *hkf;
	hl= *hlf;
	xvev= *xvevf;
	hks= *hksf;
	hls= *hlsf;
	dlh2=0.0;
	
	cw=sqrt(1.0-sw*sw);
	vev=sqrt(2.0)*mw*sw/ee;
	
	neum[0][0]= 1.0/cw/cw*mw*mw*cb*cb+hl*hls/cb*sb*xvev+hk*hl/cb*sb*xvev*xvev;
	neum[1][1]= 1.0/cw/cw*mw*mw*sb*sb+hl*hls/sb*cb*xvev+hk*hl/sb*cb*xvev*xvev
			+2.0*vev*vev*sb*sb*dlh2;
	neum[2][2]=4*hk*hk*xvev*xvev+vev*vev*hl*hls/xvev*cb*sb+hk*hks*xvev;
	
	neum[0][1]=neum[1][0]= -1.0/cw/cw*mw*mw*cb*sb+2.0*vev*vev*cb*sb*hl*hl
			-xvev*hl*(hls+xvev*hk);
	neum[0][2]=neum[2][0]=vev*(2.0*hl*hl*xvev*cb-2*hk*hl*xvev*sb-hl*hls*sb);
	
	neum[1][2]=neum[2][1]=vev*(2.0*hl*hl*xvev*sb-2*hk*hl*xvev*cb-hl*hls*cb);

	
	if(dia_(&n,(double *)neum,(double *)hisr, (double *)hisz))
	{
		puts("CP-even higgs mass matrix diagonalization failed!");
		return 0;
	}

	return 1.0;
}

double higs2_(double *di, double *dj, double *tk)
{
	int i,j;
	
	if((*tk)<0.5)
		return 0.0;
	i=floor((*di)+0.5);
	j=floor((*dj)+0.5);
	
	if(i==0)
		return mysqrt(hisr[j-1]);
	else
		return hisz[j-1][i-1];
}

static double hipz[32][32], hipr[32];

double higp1_(double *swf, double *mwf, double *eef, double *sbf, 
	double *cbf, double *hkf, double *hlf, double *xvevf, double *hksf, 
	double *hlsf)
{
    double sw; double mw; double ee; double sb; 
	double cb; double hk; double hl; double xvev; double hks; double hls;
	double neum[32][32];
	int n=3;
	int i, j, f=1;
	double cw;
	double vev;
	
	sw= *swf;
	mw= *mwf;
	ee= *eef;
	sb= *sbf;
	cb= *cbf;
	hk= *hkf;
	hl= *hlf;
	xvev= *xvevf;
	hks= *hksf;
	hls= *hlsf;
	
	cw=sqrt(1.0-sw*sw);
	vev=sqrt(2.0)*mw*sw/ee;
	
	neum[0][0]= hl*xvev*sb/cb*(hls+hk*xvev);
	neum[1][1]= hl*xvev*cb/sb*(hls+hk*xvev);
	neum[2][2]=4*vev*vev*hk*hl*cb*sb-3*hk*hks*xvev+vev*vev*hl*hls/xvev*cb*sb;
	
	neum[0][1]=neum[1][0]= hl*xvev*(hls+hk*xvev);
	neum[0][2]=neum[2][0]= hl*vev*sb*(hls-2*hk*xvev);
	
	neum[1][2]=neum[2][1]= hl*vev*cb*(hls-2*hk*xvev);

	
	if(dia_(&n,(double *)neum,(double *)hipr, (double *)hipz))
		return 0;

	while(f)
	{
		f=0;
		for(i=0;i<2;i++)
			if(fabs(hipr[i])<0.1 || (fabs(hipr[i])>0.1 && 
					fabs(hipr[i+1])>0.1 && hipr[i]>hipr[i+1]))
			{
				double tmp;
				f=1;
				tmp=hipr[i];
				hipr[i]=hipr[i+1];
				hipr[i+1]=tmp;
				for(j=0;j<3;j++)
				{
					tmp=hipz[i][j];
					hipz[i][j]=hipz[i+1][j];
					hipz[i+1][j]=tmp;
				}
			}
	}
	return 1.0;
}

double higp2_(double *di, double *dj, double *tk)
{
	int i,j;
	
	if((*tk)<0.5)
		return 0.0;
	i=floor((*di)+0.5);
	j=floor((*dj)+0.5);
	if(i==0)
		return mysqrt(hipr[j-1]);
	else
		return hipz[j-1][i-1];
}

typedef double doublereal;
typedef int integer;


/* Table of constant values */

static integer c__32 = 32;
static /* Subroutine */ int eisch1_();

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*    This is just a CERNLIB E202 general complex matrix diagonalization c*/
/*    procedure                                                          c*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
static doublereal d_sign(x, y)
doublereal *x, *y;
{
    /* System generated locals */
    doublereal ret_val;

    if (*y > (float)0.) {
	ret_val = *x;
    } else {
	ret_val = -(*x);
    }
    return ret_val;
} /* sign_ */

/* Subroutine */ static int dia_(n, a, am, z__)
integer *n;
doublereal *a, *am, *z__;
{
    static integer ierr, i__, i, j;
    static doublereal w[1024];
    static doublereal ac[1024]	/* was [32][32] */, zc[1024]	/* was [32][
	    32] */;

    /* Parameter adjustments */
    z__ -= 33;
    --am;
    a -= 33;

    /* Function Body */
    for (i__ = 1; i__ <= 32; ++i__) {
	for (j = 1; j <= 32; ++j) {
	    ac[i__ + (j << 5) - 33] = (float)0.;
/* L1: */
	    zc[i__ + (j << 5) - 33] = (float)0.;
	}
    }
    eisch1_(&c__32, n, &a[33], ac, &am[1], &z__[33], zc, &ierr, (doublereal *)w);
    
/*	for(j=1;j<=(*n);j++)
		printf("%f ",mysqrt(am[j]));
	puts("");
	for(i=0;i<(*n);i++)
		{
		printf("c %dx ",i);
		for(j=0;j<(*n);j++)
			printf(" %9.2e ",zc[i*(*n)+j]);
		printf("\n");
		}
	for(j=0;j<(*n);j++)
		printf("%f ",zc[j]);
	puts("");*/
	
	
	return ierr;
} /* dia_ */


static /* Subroutine */ int htridi_(), htribk_(), tql2_();

/* Subroutine */ static int eisch1_(nm, n, ar, ai, wr, zr, zi, ierr, work)
integer *nm, *n;
doublereal *ar, *ai, *wr, *zr, *zi;
integer *ierr;
doublereal *work;
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, zr_dim1, zr_offset, 
	    zi_dim1, zi_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;

    /* Parameter adjustments */
    zi_dim1 = *nm;
    zi_offset = zi_dim1 + 1;
    zi -= zi_offset;
    zr_dim1 = *nm;
    zr_offset = zr_dim1 + 1;
    zr -= zr_offset;
    ai_dim1 = *nm;
    ai_offset = ai_dim1 + 1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = ar_dim1 + 1;
    ar -= ar_offset;
    --wr;
    --work;

    /* Function Body */
    htridi_(nm, n, &ar[ar_offset], &ai[ai_offset], &wr[1], &zi[zi_offset], &
	    zi[zi_offset], &work[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L50: */
	    zr[i__ + j * zr_dim1] = 0.;
	}
/* L100: */
	zr[i__ + i__ * zr_dim1] = 1.;
    }
    tql2_(nm, n, &wr[1], &zi[zi_offset], &zr[zr_offset], ierr);
    if (*ierr != 0) {
	return 0;
    }
    htribk_(nm, n, &ar[ar_offset], &ai[ai_offset], &work[1], n, &zr[zr_offset]
	    , &zi[zi_offset]);
    return 0;
} /* eisch1_ */

/* Subroutine */ static int htridi_(nm, n, ar, ai, d__, e, e2, tau)
integer *nm, *n;
doublereal *ar, *ai, *d__, *e, *e2, *tau;
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal f, g, h__;
    static integer i__, j, k, l;
    static doublereal scale, fi, gi, hh;
    static integer ii;
    static doublereal si;
    static integer jp1;

    /* Parameter adjustments */
    tau -= 3;
    --e2;
    --e;
    --d__;
    ai_dim1 = *nm;
    ai_offset = ai_dim1 + 1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = ar_dim1 + 1;
    ar -= ar_offset;

    /* Function Body */
    tau[(*n << 1) + 1] = 1.;
    tau[(*n << 1) + 2] = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L100: */
	d__[i__] = ar[i__ + i__ * ar_dim1];
    }
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *n + 1 - ii;
	l = i__ - 1;
	h__ = 0.;
	scale = 0.;
	if (l < 1) {
	    goto L130;
	}
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
/* L120: */
	    scale = scale + (d__1 = ar[i__ + k * ar_dim1], abs(d__1)) + (d__2 
		    = ai[i__ + k * ai_dim1], abs(d__2));
	}
	if (scale != 0.) {
	    goto L140;
	}
	tau[(l << 1) + 1] = 1.;
	tau[(l << 1) + 2] = 0.;
L130:
	e[i__] = 0.;
	e2[i__] = 0.;
	goto L290;
L140:
	i__2 = l;
	for (k = 1; k <= i__2; ++k) {
	    ar[i__ + k * ar_dim1] /= scale;
	    ai[i__ + k * ai_dim1] /= scale;
	    h__ = h__ + ar[i__ + k * ar_dim1] * ar[i__ + k * ar_dim1] + ai[
		    i__ + k * ai_dim1] * ai[i__ + k * ai_dim1];
/* L150: */
	}
	e2[i__] = scale * scale * h__;
	g = sqrt(h__);
	e[i__] = scale * g;
/* Computing 2nd power */
	d__1 = ar[i__ + l * ar_dim1];
/* Computing 2nd power */
	d__2 = ai[i__ + l * ai_dim1];
	f = sqrt(d__1 * d__1 + d__2 * d__2);
	if (f == 0.) {
	    goto L160;
	}
	tau[(l << 1) + 1] = (ai[i__ + l * ai_dim1] * tau[(i__ << 1) + 2] - ar[
		i__ + l * ar_dim1] * tau[(i__ << 1) + 1]) / f;
	si = (ar[i__ + l * ar_dim1] * tau[(i__ << 1) + 2] + ai[i__ + l * 
		ai_dim1] * tau[(i__ << 1) + 1]) / f;
	h__ += f * g;
	g = g / f + 1.;
	ar[i__ + l * ar_dim1] = g * ar[i__ + l * ar_dim1];
	ai[i__ + l * ai_dim1] = g * ai[i__ + l * ai_dim1];
	if (l == 1) {
	    goto L270;
	}
	goto L170;
L160:
	tau[(l << 1) + 1] = -tau[(i__ << 1) + 1];
	si = tau[(i__ << 1) + 2];
	ar[i__ + l * ar_dim1] = g;
L170:
	f = 0.;
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    g = 0.;
	    gi = 0.;
	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
		g = g + ar[j + k * ar_dim1] * ar[i__ + k * ar_dim1] + ai[j + 
			k * ai_dim1] * ai[i__ + k * ai_dim1];
		gi = gi - ar[j + k * ar_dim1] * ai[i__ + k * ai_dim1] + ai[j 
			+ k * ai_dim1] * ar[i__ + k * ar_dim1];
/* L180: */
	    }
	    jp1 = j + 1;
	    if (l < jp1) {
		goto L220;
	    }
	    i__3 = l;
	    for (k = jp1; k <= i__3; ++k) {
		g = g + ar[k + j * ar_dim1] * ar[i__ + k * ar_dim1] - ai[k + 
			j * ai_dim1] * ai[i__ + k * ai_dim1];
		gi = gi - ar[k + j * ar_dim1] * ai[i__ + k * ai_dim1] - ai[k 
			+ j * ai_dim1] * ar[i__ + k * ar_dim1];
/* L200: */
	    }
L220:
	    e[j] = g / h__;
	    tau[(j << 1) + 2] = gi / h__;
	    f = f + e[j] * ar[i__ + j * ar_dim1] - tau[(j << 1) + 2] * ai[i__ 
		    + j * ai_dim1];
/* L240: */
	}
	hh = f / (h__ + h__);
	i__2 = l;
	for (j = 1; j <= i__2; ++j) {
	    f = ar[i__ + j * ar_dim1];
	    g = e[j] - hh * f;
	    e[j] = g;
	    fi = -ai[i__ + j * ai_dim1];
	    gi = tau[(j << 1) + 2] - hh * fi;
	    tau[(j << 1) + 2] = -gi;
	    i__3 = j;
	    for (k = 1; k <= i__3; ++k) {
		ar[j + k * ar_dim1] = ar[j + k * ar_dim1] - f * e[k] - g * ar[
			i__ + k * ar_dim1] + fi * tau[(k << 1) + 2] + gi * ai[
			i__ + k * ai_dim1];
		ai[j + k * ai_dim1] = ai[j + k * ai_dim1] - f * tau[(k << 1) 
			+ 2] - g * ai[i__ + k * ai_dim1] - fi * e[k] - gi * 
			ar[i__ + k * ar_dim1];
/* L260: */
	    }
	}
L270:
	i__3 = l;
	for (k = 1; k <= i__3; ++k) {
	    ar[i__ + k * ar_dim1] = scale * ar[i__ + k * ar_dim1];
	    ai[i__ + k * ai_dim1] = scale * ai[i__ + k * ai_dim1];
/* L280: */
	}
	tau[(l << 1) + 2] = -si;
L290:
	hh = d__[i__];
	d__[i__] = ar[i__ + i__ * ar_dim1];
	ar[i__ + i__ * ar_dim1] = hh;
	ai[i__ + i__ * ai_dim1] = scale * scale * h__;
/* L300: */
    }
    return 0;
} /* htridi_ */

/* Subroutine */ static int htribk_(nm, n, ar, ai, tau, m, zr, zi)
integer *nm, *n;
doublereal *ar, *ai, *tau;
integer *m;
doublereal *zr, *zi;
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, zr_dim1, zr_offset, 
	    zi_dim1, zi_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal h__;
    static integer i__, j, k, l;
    static doublereal s, si;

    /* Parameter adjustments */
    tau -= 3;
    ai_dim1 = *nm;
    ai_offset = ai_dim1 + 1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = ar_dim1 + 1;
    ar -= ar_offset;
    zi_dim1 = *nm;
    zi_offset = zi_dim1 + 1;
    zi -= zi_offset;
    zr_dim1 = *nm;
    zr_offset = zr_dim1 + 1;
    zr -= zr_offset;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    zi[k + j * zi_dim1] = -zr[k + j * zr_dim1] * tau[(k << 1) + 2];
	    zr[k + j * zr_dim1] *= tau[(k << 1) + 1];
/* L50: */
	}
    }
    if (*n == 1) {
	goto L200;
    }
    i__2 = *n;
    for (i__ = 2; i__ <= i__2; ++i__) {
	l = i__ - 1;
	h__ = ai[i__ + i__ * ai_dim1];
	if (h__ == 0.) {
	    goto L140;
	}
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    s = 0.;
	    si = 0.;
	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
		s = s + ar[i__ + k * ar_dim1] * zr[k + j * zr_dim1] - ai[i__ 
			+ k * ai_dim1] * zi[k + j * zi_dim1];
		si = si + ar[i__ + k * ar_dim1] * zi[k + j * zi_dim1] + ai[
			i__ + k * ai_dim1] * zr[k + j * zr_dim1];
/* L110: */
	    }
	    s /= h__;
	    si /= h__;
	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
		zr[k + j * zr_dim1] = zr[k + j * zr_dim1] - s * ar[i__ + k * 
			ar_dim1] - si * ai[i__ + k * ai_dim1];
		zi[k + j * zi_dim1] = zi[k + j * zi_dim1] - si * ar[i__ + k * 
			ar_dim1] + s * ai[i__ + k * ai_dim1];
/* L120: */
	    }
/* L130: */
	}
L140:
	;
    }
L200:
    return 0;
} /* htribk_ */

/* Subroutine */ static int tql2_(nm, n, d__, e, z__, ierr)
integer *nm, *n;
doublereal *d__, *e, *z__;
integer *ierr;
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(), d_sign();

    /* Local variables */
    static doublereal b, c__, f, g, h__;
    static integer i__, j, k, l, m;
    static doublereal p, r__, s;
    static integer ii;
    static doublereal machep;
    static integer mml;

    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    --e;
    --d__;

    /* Function Body */
    machep = 1.1920928955078125e-7;
    *ierr = 0;
    if (*n == 1) {
	goto L1001;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L100: */
	e[i__ - 1] = e[i__];
    }
    f = 0.;
    b = 0.;
    e[*n] = 0.;
    i__1 = *n;
    for (l = 1; l <= i__1; ++l) {
	j = 0;
	h__ = machep * ((d__1 = d__[l], abs(d__1)) + (d__2 = e[l], abs(d__2)))
		;
	if (b < h__) {
	    b = h__;
	}
	i__2 = *n;
	for (m = l; m <= i__2; ++m) {
	    if ((d__1 = e[m], abs(d__1)) <= b) {
		goto L120;
	    }
/* L110: */
	}
L120:
	if (m == l) {
	    goto L220;
	}
L130:
	if (j == 30) {
	    goto L1000;
	}
	++j;
	p = (d__[l + 1] - d__[l]) / (e[l] * 2.);
	r__ = sqrt(p * p + 1.);
	h__ = d__[l] - e[l] / (p + d_sign(&r__, &p));
	i__2 = *n;
	for (i__ = l; i__ <= i__2; ++i__) {
/* L140: */
	    d__[i__] -= h__;
	}
	f += h__;
	p = d__[m];
	c__ = 1.;
	s = 0.;
	mml = m - l;
	i__2 = mml;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = m - ii;
	    g = c__ * e[i__];
	    h__ = c__ * p;
	    if (abs(p) < (d__1 = e[i__], abs(d__1))) {
		goto L150;
	    }
	    c__ = e[i__] / p;
	    r__ = sqrt(c__ * c__ + 1.);
	    e[i__ + 1] = s * p * r__;
	    s = c__ / r__;
	    c__ = 1. / r__;
	    goto L160;
L150:
	    c__ = p / e[i__];
	    r__ = sqrt(c__ * c__ + 1.);
	    e[i__ + 1] = s * e[i__] * r__;
	    s = 1. / r__;
	    c__ *= s;
L160:
	    p = c__ * d__[i__] - s * g;
	    d__[i__ + 1] = h__ + s * (c__ * g + s * d__[i__]);
	    i__3 = *n;
	    for (k = 1; k <= i__3; ++k) {
		h__ = z__[k + (i__ + 1) * z_dim1];
		z__[k + (i__ + 1) * z_dim1] = s * z__[k + i__ * z_dim1] + c__ 
			* h__;
		z__[k + i__ * z_dim1] = c__ * z__[k + i__ * z_dim1] - s * h__;
/* L180: */
	    }
/* L200: */
	}
	e[l] = s * p;
	d__[l] = c__ * p;
	if ((d__1 = e[l], abs(d__1)) > b) {
	    goto L130;
	}
L220:
	d__[l] += f;
/* L240: */
    }
    i__1 = *n;
    for (ii = 2; ii <= i__1; ++ii) {
	i__ = ii - 1;
	k = i__;
	p = d__[i__];
	i__2 = *n;
	for (j = ii; j <= i__2; ++j) {
	    if (d__[j] >= p) {
		goto L260;
	    }
	    k = j;
	    p = d__[j];
L260:
	    ;
	}
	if (k == i__) {
	    goto L300;
	}
	d__[k] = d__[i__];
	d__[i__] = p;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    p = z__[j + i__ * z_dim1];
	    z__[j + i__ * z_dim1] = z__[j + k * z_dim1];
	    z__[j + k * z_dim1] = p;
/* L280: */
	}
L300:
	;
    }
    goto L1001;
L1000:
    *ierr = l;
L1001:
    return 0;
} /* tql2_ */


