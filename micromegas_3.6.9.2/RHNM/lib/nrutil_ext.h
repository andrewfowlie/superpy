#ifndef __NRUTIL_EXT_H__
#define __NRUTIL_EXT_H__

extern  void nrerror(char *error_text);
extern  float *vector(int nl,int nh);
extern  int *ivector(int nl,int nh);
extern  double *dvector(int nl,int nh);
extern  float **matrix(int nrl,int nrh,int ncl,int nch);
extern  double **dmatrix(int nrl,int nrh,int ncl,int nch);
extern  int **imatrix(int nrl,int nrh,int ncl,int nch);
extern  unsigned char **umatrix(int nrl,int nrh,int ncl,int nch);
extern  float **submatrix(float **a,int oldrl,int oldrh,int oldcl,int oldch,int newrl,int newcl);
extern  float ***tab3d(int nrl,int nrh,int ncl,int nch,int nzl,int nzh);
extern  int ***itab3d(int nrl,int nrh,int ncl,int nch,int nzl,int nzh);
extern  double ***dtab3d(int nrl,int nrh,int ncl,int nch,int nzl,int nzh);
extern  void free_vector(float *v,int nl,int nh);
extern  void free_ivector(int *v,int nl,int nh);
extern  void free_dvector(double *v,int nl,int nh);
extern  void free_matrix(float **m,int nrl,int nrh,int ncl,int nch);
extern  void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
extern  void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);
extern  void free_umatrix(unsigned char **m,int nrl,int nrh,int ncl,int nch);
extern  void free_tab3d(float ***m,int nrl,int nrh,int ncl,int nch,int nzl,int nzh);
extern  void free_itab3d(int ***m,int nrl,int nrh,int ncl,int nch,int nzl,int nzh);
extern  void free_dtab3d(double ***m,int nrl,int nrh,int ncl,int nch,int nzl,int nzh);
extern  void free_submatrix(float **b,int nrl,int nrh,int ncl,int nch);
extern  float **convert_matrix(float *a,int nrl,int nrh,int ncl,int nch);
extern  void free_convert_matrix(float **b,int nrl,int nrh,int ncl,int nch);

#endif /* __NRUTIL_EXT_H__ */
