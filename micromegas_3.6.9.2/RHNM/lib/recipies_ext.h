#ifndef __RECIPIES_EXT_H__
#define __RECIPIES_EXT_H__

extern 	void  adi(double **a, double **b, double **c, double **d, double **e,
							double **f, double **g, double **u, int jmax, int k, double alpha, double beta, double eps);
extern 	void amebsa(float **p,float y[],int ndim, float pb[],float *yb,float ftol,
	    					float (*funk)(float[]),int *iter,float temptr);
extern 	void amoeba(float** p, float y[], int ndim, float ftol, float (*funk)(float*), int* nfunk);
extern 	float amotsa(float **p,float y[],float psum[],int ndim,float pb[],
	     					float *yb,float (*funk)(float[]),int ihi,float *yhi,float fac);
extern 	void  anneal(float *x, float *y, int *iorder, int ncity);
extern 	void  avevar(float *data, int n, float *ave, float *svar);
extern 	void  balanc(float **a, int n);
extern 	void  bcucof(float *y, float *y1, float *y2, float *y12, float d1,
		float d2, float **c);
extern 	void  bcuint(float *y, float *y1, float *y2, float *y12, float x1l, 
		float x1u, float x2l, float x2u, float x1, float x2, 
		float *ansy, float *ansy1, float *ansy2);
extern 	float bessi(int n, float x);
extern 	float bessi0(float x);
extern 	float bessi1(float x);
extern 	float bessj(int n, float x);
extern 	float bessj0(float x);
extern 	float bessj1(float x);
extern 	float bessk(int n, float x);
extern 	float bessk0(float x);
extern 	float bessk1(float x);
extern 	float bessy(int n, float x);
extern 	float bessy0(float x);
extern 	float bessy1(float x);
extern 	float beta(float z, float w);
extern 	float betacf(float a, float b, float x);
extern 	float betai(float a, float b, float x);
extern 	float bico(int n, int k);
extern 	void  bksub(int ne, int nb, int jf, int k1, int k2, float ***c);
extern 	float bnldev(float pp, int n, int *idum);
extern 	float brent(float ax, float bx, float cx, float (*f)(float), float tol,
		float *xmin);
extern 	void  bsstep(float *y, float *dydx, int nv, float *xx, float htry,
		float eps, float *yscal, float *hdid, float *hnext, 
		void (*derivs)(float,float *,float *));
extern 	void  caldat(long julian, int *mm, int *id, int *iyyy);
extern 	float cel(float qqc, float pp, float aa, float bb);
extern 	void  chder(float a, float b, float *c, float *cder, int n);
extern 	float chebev(float a, float b, float *c, int m, float x);
extern 	void  chebft(float a, float b, float *c, int n, float (*func)(float));
extern 	void  chebpc(float *c, float *d, int n);
extern 	void  chint(float a, float b, float *c, float *cint, int n);
extern 	void  chsone(float *bins, float *ebins, int nbins, int knstrn, 
		float *df, float *chsq, float *prob);
extern 	void  chstwo(float *bins1, float *bins2, int nbins, int knstrn, 
		float *df, float *chsq, float *prob);
extern 	void  cntab1(int **nn, int n1, int nj, float *chisq, float *df, 
		float *prob, float *cramrv, float *ccc);
extern 	void  cntab2(int **nn, int ni, int nj, float *h, float *hx, float *hy, 
		float *hygx, float *hxgy, float *uygx, float *uxgy,
		float *uxy);
extern 	void  convlv(float *data, int n, float *respns, int m, int isign, 
		float *ans);
extern 	void  correl(float *data1, float *data2, int n, float *ans);
extern 	void  cosft(float *y, int n, int isign);
extern 	void  covsrt(float **covar, int ma, int *lista, int mfit);
extern 	void  crank(int n, float *w, float *s);
extern 	float dbrent(float ax, float bx, float cx, float (*f)(float),
		float (*df)(float), float tol, float *xmin);
extern 	void  ddpoly(float *c, int nc, float x, float *pd, int nd);
extern 	void  des(immense inp, immense key, int *newkey, int isw, immense *out);
extern 	unsigned long getbit(immense source, int bitno, int nbits);
extern 	void  ks(immense key, int n, great *kn);
extern 	void  cyfun(unsigned long ir, great k, unsigned long *iout);
extern 	float df1dim(float x);
extern 	void  dfpmin(float *p, int n, float ftol, int *iter, float *fret, 
		float (*func)(float *), void (*dfunc)(float *,float *));
extern 	void  difeq(int k, int k1, int k2, int jsf, int is1, int isf, 
		int *indexv, int ne, float **s, float **y);
extern 	void  dlinmin(float *p, float *xi, int n, float *fret,
		float (*func)(float *), void (*dfunc)(float *,float *));
extern 	void  eclass(int *nf, int n, int *lista, int *listb, int m);
extern 	void  eclazz(int *nf, int n, int (*equiv)(int,int));
extern 	void  eigsrt(float *d, float **v, int n);
extern 	float el2(float x, float qqc, float aa, float bb);
extern 	void  elmhes(float **a, int n);
#ifdef _WINDOWS
extern 	float erf(float x);
extern 	float erfc(float x);
#endif
extern 	float erfcc(float x);
extern 	void  eulsum(float *sum, float term, int jterm, float *wksp);
extern 	float evlmem(float fdt, float *cof, int m, float pm);
extern 	float expdev(int *idum);
extern 	float f1dim(float x);
extern 	float factln(int n);
extern 	float factrl(int n);
extern 	void  fgauss(float x, float *a, float *y, float *dyda, int na);
extern 	void  fit(float *x, float *y, int ndata, float *sig, int mwt, float *a, 
		float *b, float *siga, float *sigb, float *chi2, float *q);
extern 	void  fixrts(float *d, int npoles);
extern 	void  fleg(float x, float *pl, int nl);
extern 	void  flmoon(int n, int nph, long *jd, float *frac);
extern 	void  four1(float *data, int nn, int isign);
extern 	void  fourn(float *data, int *nn, int ndim, int isign);
extern 	void  fpoly(float x, float *p, int np);
extern 	void  frprmn(float *p, int n, float ftol, int *iter, float *fret, 
		float (*func)(float *), void (*dfunc)(float *,float *));
extern 	void  ftest(float *data1, int n1, float *data2, int n2, float *f, 
		float *prob);
extern 	float gamdev(int ia, int *idum);
extern 	float gammln(float xx);
extern 	float gammp(float a, float x);
extern 	float gammq(float a, float x);
extern 	float gasdev(int *idum);
extern 	void  gauleg(double x1, double x2, double *x, double *w, int n);
extern 	void  gaussj(float **a, int n, float **b, int m);
extern 	void  gcf(float *gammcf, float a, float x, float *gln);
extern 	float golden(float ax, float bx, float cx, float (*f)(float), float tol, 
		float *xmin);
extern 	void  gser(float *gamser, float a, float x, float *gln);
extern 	void  hqr(float **a, int n, float *wr, float *wi);
extern 	void  hunt(float *xx, int n, float x, int *jlo);
extern 	void  indexx(int n, float *arrin, int *indx);
extern 	int   irbit1(unsigned long int *iseed);
extern 	int   irbit2(unsigned long int *iseed);
extern 	void  jacobi(float **a, int n, float *d, float **v, int *nrot);
extern 	long  julday(int mm, int id, int iyyy);
extern 	void  kendl1(float *data1, float *data2, int n, float *tau, float *z,
		float *prob);
extern 	void  kendl2(float **tab, int i, int j, float *tau, float *z,
		float *prob);
extern 	void  ksone(float *data, int n, float (*func)(float), float *d,
		float *prob);
extern 	void  kstwo(float *data1, int n1, float *data2, int n2, float *d,
		float *prob);
extern 	void  laguer(fcomplex *a, int m, fcomplex *x, float eps, int polish);
extern 	void  lfit(float *x, float *y, float *sig, int ndata, float *a, int ma, 
		int *lista, int mfit, float **covar, float *chisq,
		void (*funcs)(float,float *,int));
extern 	void  linmin(float *p, float *xi, int n, float *fret, float (*func)(float*));
extern 	void  locate(float *xx, int n, float x, int *j);
extern 	void  lubksb(float **a, int n, int *indx, float *b);
extern 	void  ludcmp(float **a, int n, int *indx, float *d);
extern 	void  mdian1(float *x, int n, float *xmed);
extern 	void  mdian2(float *x, int n, float *xmed);
extern 	void  medfit(float *x, float *y, int ndata, float *a, float *b,
		float *abdev);
extern 	void  memcof(float *data, int n, int m, float *pm, float *cof);
extern 	float midexp(float (*funk)(float), float aa, float bb, int n);
extern 	float midinf(float (*funk)(float), float aa, float bb, int n);
extern 	float midpnt(float (*func)(float), float a, float b, int n);
extern 	float midsql(float (*funk)(float), float aa, float bb, int n);
extern 	float midsqu(float (*funk)(float), float aa, float bb, int n);
extern 	void  mmid(float *y, float *dydx, int nvar, float xs, float htot,
		int nstep, float *yout,
		void (*derivs)(float,float *,float *));
extern 	void  mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, 
		float *fc, float (*func)(float));
extern 	void  mnewt(int ntrial, float *x, int n, float tolx, float tolf);
extern 	void  moment(float *data, int n, float *ave, float *adev, float *sdev, 
		float *svar, float *skew, float *curt);
extern 	void  mprove(float **a, float **alud, int n, int *indx, float *b, 
		float *x);
extern 	void  mrqcof(float *x, float *y, float *sig, int ndata, float *a, int ma, 
		int *lista, int mfit, float **alpha, float *beta, float
		*chisq, void (*funcs)(float,float *,float *,float *,int));
extern 	void  mrqmin(float *x, float *y, float *sig, int ndata, float *a,
		int ma, int *lista, int mfit, float **covar, float **alpha, 
		float *chisq, void (*funcs)(float,float *,float *,float *,
		int),float *alamda);
extern 	void  odeint(float *ystart, int nvar, float x1, float x2, float eps,
		float h1, float hmin, int *nok, int *nbad,
		void (*derivs)(float,float *,float *),
		void  (*rkqc)(float *,float *,int,float *,float,float,float
		*,float *,float *,void (*)(float,float *,float *)));
extern 	void  pcshft(float a, float b, float *d, int n);
extern 	void  pearsn(float *x, float *y, int n, float *r, float *prob, float *z);
extern 	void  piksr2(int n, float *arr, float *brr);
extern 	void  piksrt(int n, float *arr);
extern 	void  pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
		float ***c, float **s);
extern 	float plgndr(int l, int m, float x);
extern 	float poidev(float xm, int *idum);
extern 	void  polcoe(float *x, float *y, int n, float *cof);
extern 	void  polcof(float *xa, float *ya, int n, float *cof);
extern 	void  poldiv(float *u, int n, float *v, int nv, float *q, float *r);
extern 	void  polin2(float *x1a, float *x2a, float **ya, int m, int n, float x1, 
		float x2, float *y, float *dy);
extern 	void  polint(float *xa, float *ya, int n, float x, float *y, float *dy);
extern 	void  polintd(double *xa, double *ya, int n, double x, double *y, double *dy);
extern 	void  powell(float *p, float **xi, int n, float ftol, int *iter,
		float *fret, float (*func)(float *));
extern 	void  predic(float *data, int ndata, float *d, int npoles, 
		float *future, int nfut);
extern 	float probks(float alam);
extern 	void  pzextr(int iest, float xest, float *yest, float *yz, float *dy,
							int nv, int nuse);
extern 	void  qcksrt(int n, float *arr);
extern 	float qgaus(float (*func)(float), float a, float b);
extern 	float qromb(float (*func)(float), float a, float b);
extern 	double qrombd(double (*func)(double), double a, double b);
extern 	float qromo(float (*func)(), float a, float b, 
							float (*choose)(float (*func)(), float,float,int));
extern 	void  qroot(float *p, int n, float *b, float *c, float eps);
extern 	float qsimp(float (*func)(float), float a, float b);
extern 	float qtrap(float (*func)(float), float a, float b);
extern 	float quad3d(float (*func)(float,float,float), float x1, float x2);
extern 	float ran0(int *idum);
extern 	float ran1(int *idum);
extern 	float ran2(long *idum);
extern 	float ran3(int *idum);
extern 	float ran4(int *idum);
extern 	void  rank(int n, int *indx, int *irank);
extern 	void  ratint(float *xa, float *ya, int n, float x, float *y, float *dy);
extern 	void  realft(float *data, int n, int isign);
extern 	void  red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
		int ic1, int jc1, int jcf, int kc, float ***c, float **s);
extern 	void  rk4(float *y, float *dydx, int n, float x, float h, float *yout,
		void (*derivs)(float,float *,float *));
extern 	void  rkdumb(float *vstart, int nvar, float x1, float x2, int nstep, 
		void (*derivs)(float,float *,float *));
extern 	void  rkqc(float *y, float *dydx, int n, float *x, float htry, 
		float eps, float *yscal, float *hdid, float *hnext, 
		void (*derivs)(float,float *,float *));
extern 	float rofunc(float b);
extern 	float rtbis(float (*func)(float), float x1, float x2, float xacc);
extern 	float rtflsp(float (*func)(float), float x1, float x2, float xacc);
extern 	float rtnewt(void (*funcd)(float,float *,float *), float x1, float x2,
		float xacc);
extern 	float rtsafe(void (*funcd)(float,float *,float *), float x1, float x2,
		float xacc);
extern 	float rtsec(float (*func)(float), float x1, float x2, float xacc);
extern 	void  rzextr(int iest, float xest, float *yest, float *yz, float *dy,
		int nv, int nuse);
extern 	void  scrsho(float (*fx)(float));
extern 	void  shell(int n, float *arr);
extern 	void  shoot(int nvar, float *v, float *delv, int n2, float x1, float x2,
		float eps, float h1, float hmin, float *f, float *dv);
extern 	void  shootf(int nvar, float *v1, float *v2, float *delv1, float *delv2,
		int n1, int n2, float x1, float x2, float xf, float eps, 
		float h1, float hmin, float *f, float *dv1, float *dv2);
extern 	void  simp1(float **a, int mm, int *ll, int nll, int iabf, int *kp,
		float *bmax);
extern 	void  simp2(float **a, int n, int *l2, int nl2, int *ip, int kp,
		float *q1);
extern 	void  simp3(float **a, int i1, int k1, int ip, int kp);
extern 	void  simplx(float **a, int m, int n, int m1, int m2, int m3, 
		int *icase, int *izrov, int *iposv);
extern 	void  sinft(float *y, int n);
extern 	void  smooft(float *y, int n, float pts);
extern 	void  sncndn(float uu, float emmc, float *sn, float *cn, float *dn);
extern 	void  solvde(int itmax, float conv, float slowc, float *scalv,
		int *indexv, int ne, int nb, int m, float **y, float ***c,
		float **s);
extern 	void  sor(double **a, double **b, double **c, double **d, double **e, 
		double **f, double **u, int jmax, double rjac);
extern 	void  sort(int n, float *ra);
extern 	void  sort2(int n, float *ra, float *rb);
extern 	void  sort3(int n, float *ra, float *rb, float *rc);
extern 	void  sparse(float *b, int n, float *x, float *rsq);
extern 	void  spctrm(FILE *fp, float *p, int m, int k, int ovrlap);
extern 	void  spear(float *data1, float *data2, int n, float *d, float *zd,
		float *probd, float *rs, float *probrs);
extern 	void  splie2(float *x1a, float *x2a, float **ya, int m, int n,
		float **y2a);
extern 	void  splin2(float *x1a, float *x2a, float **ya, float **y2a, int m,
		int n, float x1, float x2, float *y);
extern 	void  spline(float *x, float *y, int n, float yp1, float ypn, float *y2);
extern 	void  splint(float *xa, float *ya, float *y2a, int n, float x, float *y);
extern 	void  svbksb(float **u, float *w, float **v, int m, int n, float *b,
		float *x);
extern 	void  svdcmp(float **a, int m, int n, float *w, float **v);
extern 	void  svdfit(float *x, float *y, float *sig, int ndata, float *a, 
		int ma, float **u, float **v, float *w, float *chisq,
		void (*funcs)(float,float *,int));
extern 	void  svdvar(float **v, int ma, float *w, float **cvm);
extern 	void  toeplz(float *r, float *x, float *y, int n);
extern 	void  tptest(float *data1, float *data2, int n, float *t, float *prob);
extern 	void  tqli(float *d, float *e, int n, float **z);
extern 	float trapzd(float (*func)(float), float a, float b, int n);
extern 	double trapzdou(double (*func)(double), double a, double b, int n);
extern 	void  tred2(float **a, int n, float *d, float *e);
extern 	void  tridag(float *a, float *b, float *c, float *r, float *u, int n);
extern 	void  ttest(float *data1, int n1, float *data2, int n2, float *t,
		float *prob);
extern 	void  tutest(float *data1, int n1, float *data2, int n2, float *t,
		float *prob);
extern 	void  twofft(float *data1, float *data2, float *fft1, float *fft2,
		int n);
extern 	void  vander(float *x, float *w, float *q, int n);
extern 	int   zbrac(float (*func)(float), float *x1, float *x2);
extern 	void  zbrak(float (*fx)(float), float x1, float x2, int n, float *xb1,
		float *xb2, int *nb);
extern 	float zbrent(float (*func)(float), float x1, float x2, float tol);
extern 	void  zroots(fcomplex *a, int m, fcomplex *roots, int polish);
#endif /*__RECIPIES_EXT_H__*/
