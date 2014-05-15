/*********************************************************/
/* File: cweight.h - header file for module of colour    */
/*                   weigh calculation.                  */
/* Author: A.Kryukov (kryukov@theory.npi.msu.su)         */
/*-------------------------------------------------------*/
/* Release:  25/05/98  Add new external function cwtarg0 */
/*********************************************************/
#ifndef __CWEIGHT_
#define __CWEIGHT_
  extern void  cwtarg(vcsect * g);
  extern void c_basis_coef(vampl*g,int pow,int nc,int*chains,long*num,long*den);
  extern int NcInfLimit;

  extern int  generateColorWeights(csdiagram*csdiagr,
                          int cBasisPower,int nC,int*cChains,
                          long * cCoefN,long * cCoefD);

  extern int infCbases(int np,         /* number of particles */
                       int * cweight,  /* array of particle color weights */ 
                       int *nc,        /* number of color chains */
                       int *pow,       /* power of basis */
                       int ** chains   /* returns array   which descibes   
                                          (*pow) basis elements, 
                                          each of them contains  (*nc)  chains,
                                          each of them is a couple of 
                                          particle numbers 
                                        */ 
                       );

#endif
