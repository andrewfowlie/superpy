#ifndef  __PMODEL__
#define  __PMODEL__

#ifdef __cplusplus
extern "C" {
#endif

extern double simps( double (*func)(double),double a,double b, double  eps);
extern void bessjy(double  x, double xnu, double *rj, double *ry);
extern int readVarRHNM(char * fname);
extern double MLZP(double, double);
extern double g6f(double c1);
extern double g8f(double c1);
extern double gztot(double c1, double c2, double c3, double c4);
extern double g1f(double c1,double c2);
extern double SAVE(double ee, double sw, double Mtop, double v, double g10, 
double mzp, double rc);
extern double sigmap(double);
extern double sigman(double);
extern double sigmanh(double,double, double, double);
extern double sigmaph(double,double, double, double);
extern double sigmage(double,double, double, double);
extern double sigmaxe(double,double, double, double);
extern double sigmahe(double);
extern double zinv(double, double);
extern int rdVarRHNM(char * fname);
#ifdef __cplusplus
}
#endif


#endif
