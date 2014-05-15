#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <complex.h>
#include <string.h>
#include <strings.h>

/*--------------------------------------------------------------------*/

/*#define DEBUG*/
/*#define SM_ChargedHiggs*/

/*--------------------------------------------------------------------*/

#define pi    3.1415926535897932385
#define zeta3 1.2020569031595942855
#define hbar  6.58211889e-25 /* in GeV.s */

#define Nobs_BKsll 16

/*--------------------------------------------------------------------*/

typedef struct parameters
/* structure containing all the scanned parameters from the SLHA file */
{
	int SM;
	int model; /* CMSSM=1, GMSB=2, AMSB=3 */
	int generator; /* ISAJET=1, SOFTSUSY=3, SPHENO=4, SUSPECT=5, NMSSMTOOLS=6 */
	double Q; /* Qmax ; default = M_EWSB = sqrt(m_stop1*mstop2) */
	double m0,m12,tan_beta,sign_mu,A0; /* CMSSM parameters */
	double Lambda,Mmess,N5,cgrav,m32; /* AMSB, GMSB parameters */
	double mass_Z,mass_W,mass_b,mass_top_pole,mass_tau_pole; /* SM parameters */
	double inv_alpha_em,alphas_MZ,Gfermi,GAUGE_Q; /* SM parameters */
	double charg_Umix[3][3],charg_Vmix[3][3],stop_mix[3][3],sbot_mix[3][3],stau_mix[3][3],neut_mix[6][6],mass_neut[6],alpha; /* mass mixing matrices */
	double Min,M1_Min,M2_Min,M3_Min,At_Min,Ab_Min,Atau_Min,M2H1_Min,M2H2_Min,mu_Min,M2A_Min,tb_Min,mA_Min; /* optional input parameters at scale Min */
	double MeL_Min,MmuL_Min,MtauL_Min,MeR_Min,MmuR_Min,MtauR_Min; /* optional input parameters at scale Min */
	double MqL1_Min,MqL2_Min,MqL3_Min,MuR_Min,McR_Min,MtR_Min,MdR_Min,MsR_Min,MbR_Min; /* optional input parameters at scale Min */
	double N51,N52,N53,M2H1_Q,M2H2_Q; /* optional input parameters (N51...3: GMSB) */
	double mass_d,mass_u,mass_s,mass_c,mass_t,mass_e,mass_nue,mass_mu,mass_num,mass_tau,mass_nut; /* SM masses */
	double mass_gluon,mass_photon,mass_Z0; /* SM masses */
	double mass_h0,mass_H0,mass_A0,mass_H,mass_dnl,mass_upl,mass_stl,mass_chl,mass_b1,mass_t1; /* Higgs & superparticle masses */
	double mass_el,mass_nuel,mass_mul,mass_numl,mass_tau1,mass_nutl,mass_gluino,mass_cha1,mass_cha2; /* superparticle masses */
	double mass_dnr,mass_upr,mass_str,mass_chr,mass_b2,mass_t2,mass_er,mass_mur,mass_tau2; /* superparticle masses */
	double mass_nuer,mass_numr,mass_nutr,mass_graviton,mass_gravitino; /* superparticle masses */
	double gp,g2,gp_Q,g2_Q,g3_Q,YU_Q,yut[4],YD_Q,yub[4],YE_Q,yutau[4]; /* couplings */
	double HMIX_Q,mu_Q,tanb_GUT,Higgs_VEV,mA2_Q,MSOFT_Q,M1_Q,M2_Q,M3_Q; /* parameters at scale Q */
	double MeL_Q,MmuL_Q,MtauL_Q,MeR_Q,MmuR_Q,MtauR_Q,MqL1_Q,MqL2_Q,MqL3_Q,MuR_Q,McR_Q,MtR_Q,MdR_Q,MsR_Q,MbR_Q; /* masses at scale Q */
	double AU_Q,A_u,A_c,A_t,AD_Q,A_d,A_s,A_b,AE_Q,A_e,A_mu,A_tau; /* trilinear couplings */
	
	/* SLHA2 */
	int NMSSM,RV,CPV,FV;
	double mass_nutau2,mass_e2,mass_nue2,mass_mu2,mass_numu2,mass_d2,mass_u2,mass_s2,mass_c2;
	double CKM_lambda,CKM_A,CKM_rhobar,CKM_etabar;
	double PMNS_theta12,PMNS_theta23,PMNS_theta13,PMNS_delta13,PMNS_alpha1,PMNS_alpha2;
	double lambdaNMSSM_Min,kappaNMSSM_Min,AlambdaNMSSM_Min,AkappaNMSSM_Min,lambdaSNMSSM_Min,xiFNMSSM_Min,xiSNMSSM_Min,mupNMSSM_Min,mSp2NMSSM_Min,mS2NMSSM_Min,mass_H03,mass_A02,NMSSMRUN_Q,lambdaNMSSM,kappaNMSSM,AlambdaNMSSM,AkappaNMSSM,lambdaSNMSSM,xiFNMSSM,xiSNMSSM,mupNMSSM,mSp2NMSSM,mS2NMSSM; /* NMSSM parameters */
	double PMNSU_Q,CKM_Q,IMCKM_Q,MSE2_Q,MSU2_Q,MSD2_Q,MSL2_Q,MSQ2_Q,TU_Q,TD_Q,TE_Q;
	double CKM[4][4],IMCKM[4][4]; /* CKM matrix */
	double H0_mix[4][4],A0_mix[4][4]; /* Higgs mixing matrices */
	double sU_mix[7][7],sD_mix[7][7],sE_mix[7][7], sNU_mix[4][4]; /* mixing matrices */
	double sCKM_msq2[4][4],sCKM_msl2[4][4],sCKM_msd2[4][4],sCKM_msu2[4][4],sCKM_mse2[4][4]; /* super CKM matrices */
	double PMNS_U[4][4]; /* PMNS mixing matrices */
	double TU[4][4],TD[4][4],TE[4][4]; /* trilinear couplings */
	
	/* non-SLHA*/
	double mass_b_1S,mass_b_pole,mtmt;
	double Lambda5; /* Lambda QCD */
	
	/* Flavour constants */
	double f_B,f_Bs,f_Ds,f_D,fK_fpi;
	double f_K_par,f_K_perp;
	double m_B,m_Bs,m_Bd,m_pi,m_Ds,m_K,m_Kstar,m_D0,m_D;
	double life_pi,life_K,life_B,life_Bs,life_Bd,life_D,life_Ds;
	double a1par,a2par,a1perp,a2perp;
	double zeta3A,zeta3V,wA10,deltatp,deltatm;
	double lambda_Bp,rho1,lambda2;
	
	/* CKM matrix */
	double complex Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb;
	
	/* 2HDM */
	int THDM_model;
	double lambda_u[4][4],lambda_d[4][4],lambda_l[4][4];

	/* NMSSMTools */
	int NMSSMcoll,NMSSMtheory,NMSSMups1S,NMSSMetab1S;
	
	/* Decay widths */
	double width_Z,width_W;
}
parameters;

/*--------------------------------------------------------------------*/
/* Prototypes */

/* isajet.c */
int isajet_cmssm(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, char name[]);
int isajet_gmsb(double Lambda, double Mmess, double tanb, int N5, double cGrav, double sgnmu, double mtop, char name[]);
int isajet_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, char name[]);
int isajet_nuhm(double m0, double m12, double tanb, double A0, double mu, double mA, double mtop, char name[]);
int isajet_mmamsb(double alpha, double m32, double tanb, double sgnmu, double mtop, char name[]);
int isajet_hcamsb(double alpha, double m32, double tanb, double sgnmu, double mtop, char name[]);

/* softsusy.c */
int softsusy_cmssm(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_gmsb(double Lambda, double Mmess, double tanb, int N5, double cGrav, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_nuhm(double m0, double m12, double tanb, double A0, double mu, double mA, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_mssm(double m1, double m2, double m3, double tanb, double mA, double at, double ab, double atau, double mu, double mer, double mel, double mstaul, double mstaur, double mql, double mq3l, double mqur, double mqtr, double mqdr, double mqbr, double Q, double mtop, double mbot, double alphas_mz, char name[]);
int softsusy_slhain(char name_in[], char name_out[]);

/* leshouches.c */ 
int Les_Houches_Reader(char name[], struct parameters* param);
void Init_param(struct parameters* param);
void slha_adjust(struct parameters* param);
int test_slha(char name[]);

/* alphas.c */
double alphas_running(double Q, double mtop, double mbot, struct parameters* param);

/* quark_masses.c */
double running_mass(double quark_mass, double Qinit, double Qfin,  double mtop, double mbot, struct parameters* param);
double mb_pole(struct parameters* param);
double mc_pole(struct parameters* param);
double mb_1S(struct parameters* param);
double mt_mt(struct parameters* param);

/* general.c */
double max(double x, double y);
double min(double x, double y);
double Ei1(double x);
double Ei2(double x);
double Ei3(double x);
double Ei(double x);
double complex polylog(int n, int m, double x);
double complex cd(double x, double y);
double complex hpl_base1(int i, double complex x);
double complex hpl_base2(int i1, int i2, double complex x);
double complex hpl_base3(int i1, int i2, int i3, double complex x);
double complex hpl_base4(int i1, int i2, int i3, int i4, double complex x);
double complex hpl1(int i, double complex x);
double complex hpl2(int i1, int i2, double complex x);
double complex hpl3(int i1, int i2, int i3, double complex x);
double complex hpl4(int i1, int i2, int i3, int i4, double complex x);
double Li2(double x);
double Li3(double x);
double Li4(double x);
double complex CLi2(double complex x);
double complex CLi3(double complex x);
double complex CLi4(double complex x);
double Cl2(double x);
double Cl3(double x);
double I0(double x);
double I1(double x);
double K0(double x);
double K1(double x);
double K2(double x);
double K3(double x);
double K4(double x);
double Lbessel(double x);
double Mbessel(double x);
double Nbessel(double x);
double expcor(double x);
double K0exp(double x,double z);
double K1exp(double x,double z);
double K2exp(double x,double z);
double kron(int x, int y);
int test_integer(char name[]);
int test_file(char *name);

/* wilson.c */
double Bplus(double x, double y);
double D2(double x, double y);
double D3(double x);
double h10(double x);
double h20(double x);
double h30(double x);
double h40(double x);
double h50(double x);
double h60(double x);
double f20(double x);
double f30(double x,double y);
double f40(double x,double y);
double f50(double x,double y,double z);
double f60(double x,double y,double z);
double f70(double x, double y);
double f80(double x);
double f90(double w,double x,double y,double z);
double f100(double w,double x,double y,double z);
double f110(double x,double y);
double h11(double x,double y);
double h21(double x,double y);
double h31(double x,double y);
double h41(double x,double y);
double h51(double x,double y);
double h61(double x,double y);
double h71(double x,double y);
double f31(double x, double y);
double f41(double x,double y);
double f51(double x,double y);
double f61(double x,double y);
double f71(double x,double y,double z);
double f81(double x,double y,double z);
double f91(double x,double y,double z);
double f111(double x, double y);
double f121(double x, double y, double z);
double f131(double x, double y, double z);
double f141(double x, double y);
double f151(double x);
double f161(double x);
double f171(double x, double y);
double f181(double x, double y);
double f191(double x, double y);
double q11(double x,double y);
double q21(double x,double y);
double q31(double x,double y);
double q41(double x,double y);
double q51(double x,double y);
double q61(double x,double y);
double A0t(double x);
double A1t(double x, double l);
double B0t(double x);
double C0t(double x);
double D0t(double x);
double B1t(double x, double l);
double C1t(double x, double l);
double D1t(double x, double l);
double F0t(double x);
double F1t(double x,double l);
double E0t(double x);
double G1t(double x, double l);
double E1t(double x, double l);
double T(double x);
double F7_1(double x);
double F7_2(double x);
double F8_1(double x);
double F8_2(double x);
double H2(double x, double y);
double B(double m1, double m2, double Q);
double G7H(double x, double lu, double ld);
double Delta7H(double x, double lu, double ld);
double G8H(double x, double lu, double ld);
double Delta8H(double x, double lu, double ld);
double EH(double x, double lu);
double G4H(double x, double lu);
double Delta4H(double x, double lu);
double G3H(double x, double lu);
double Delta3H(double x, double lu);
double C9llH0(double x, double y, double lu);
double D9H0(double x, double lu);
double C9llH1(double x, double y, double lu, double L);
double D9H1(double x, double lu, double L);
double C7t2mt(double x);
double C7c2MW(double x);
double C8t2mt(double x);
double C8c2MW(double x);
double epsilon_0(struct parameters* param);
double epsilon_2(struct parameters* param);
double epsilon_b(struct parameters* param);
double epsilon_bp(struct parameters* param);
double epsilon_0p(struct parameters* param);
double epsilon_1p(struct parameters* param);
void CW_calculator(double C0w[], double C1w[], double C2w[], double mu_W, struct parameters* param); 
void C_calculator_base1(double C0w[], double C1w[], double C2w[], double mu_W, double C0b[], double C1b[], double C2b[], double mu, struct parameters* param); 
void C_calculator_base2(double C0w[], double C1w[], double mu_W, double C0b[], double C1b[], double mu, struct parameters* param); 
void Cprime_calculator(double Cpb[], double complex CQpb[], double mu_W, double mu, struct parameters* param);
void CQ_calculator(double complex CQ0b[], double complex CQ1b[], double mu_W, double mu, struct parameters* param);

/* bsgamma.c */
double phi77(double delta);
double phi78(double delta);
double phi88(double delta, double b);
double G1(double t);
double phi22(double delta, double z);
double phi11(double delta, double z);
double phi12(double delta, double z);
double G2(double t);
double phi27(double delta, double z);
double phi17(double delta, double z);
double phi18(double delta, double z);
double phi28(double delta, double z);
double phi47(double delta);
double phi48(double delta);
double F2_nf(double z);
double F2_a(double z);
double F2_na(double z);
double phi77_2beta(double delta,double mu, struct parameters* param);
double phi77_2rem(double delta, struct parameters* param);
double Re_a(double z);
double Re_b(double z);
double bsgamma(double C0[], double C1[], double C2[], double mu, double mu_W, struct parameters* param);
double bsgamma_calculator(char name[]);

/* isospin.c */
double F_perp(double a1_perp, double a2_perp);
double X_perp1(double a1_perp, double a2_perp);
double X_perp2(double a1_perp, double a2_perp);
double complex G(double s, double x);
double complex G_perp(double s, double a1_perp, double a2_perp);
double complex H_perp(double s, double a1_par, double a2_par, double zeta3A, double zeta3V, double wA10, double deltatp, double deltatm);
double H8_perp(double a1_perp, double a2_perp);
double complex h(double u,double s);
double complex H2_perp(double s, double a1_perp, double a2_perp);
double delta0(double C0[],double C0_spec[],double C1[],double C1_spec[],struct parameters* param,double mub,double muspec, double lambda_h);
double delta0_calculator(char name[]);

/* excluded_masses.c */
int excluded_Higgs_mass_calculator(char name[]);
int excluded_Higgs_masses(struct parameters* param);
int excluded_SUSY_mass_calculator(char name[]);
int excluded_SUSY_masses(struct parameters* param);
int excluded_mass_calculator(char name[]);
int excluded_masses(struct parameters* param);
int charged_LSP_calculator(char name[]);
int charged_LSP(struct parameters* param);

/* gmuon.c */
double F1N(double x);
double F2N(double x);
double F1C(double x);
double F2C(double x);
double fPS(double x);
double fS(double x);
double fft(double x);
double muonI1(double a);
double muonI2(double a);
double muonI3(double a);
double muonf(double z);
double muong(double z);
double muon_gm2_calculator(char name[]);
double muon_gm2(struct parameters* param);

/* btaunu.c */
double Btaunu(struct parameters* param);
double Btaunu_calculator(char name[]);
double RBtaunu(struct parameters* param);
double RBtaunu_calculator(char name[]);

/* bdtaunu.c */
double GBDlnu(double w);
double tBDlnu(double w, double m_B, double m_D);
double rhoV(double w, double ml, double m_B, double m_D);
double rhoS(double w, double ml, double m_B, double m_D);
double dGammaBDlnu_dw(double w, double ml, struct parameters* param);
double GammaBDlnu(double ml, struct parameters* param);
double BDtaunu_BDenu(struct parameters* param);
double BDtaunu_BDenu_calculator(char name[]);
double BDtaunu(struct parameters* param);
double BDtaunu_calculator(char name[]);

/* bsmumu.c */
double Bsmumu(double C0b[], double C1b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[],struct parameters* param, double mu_b);
double Bsmumu_untag(double C0b[], double C1b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[],struct parameters* param, double mu_b);
double Bdmumu(double C0b[], double C1b[], double complex CQ0b[], double complex CQ1b[],struct parameters* param, double mu_b);
double Bsmumu_calculator(char name[]);
double Bsmumu_untag_calculator(char name[]);
double Bdmumu_calculator(char name[]);

/* kmunu.c */
double Kmunu_pimunu(struct parameters* param);
double Kmunu_pimunu_calculator(char name[]);
double Rmu23(struct parameters* param);
double Rmu23_calculator(char name[]);

/* dslnu.c */
double Dstaunu(struct parameters* param);
double Dstaunu_calculator(char name[]);
double Dsmunu(struct parameters* param);
double Dsmunu_calculator(char name[]);

/* dmunu.c */
double Dmunu(struct parameters* param);
double Dmunu_calculator(char name[]);

/* 2hdmc.c */
int thdmc_types(double l1, double l2, double l3, double l4, double l5, double l6, double l7, double m12_2, double tanb, int type, char name[]);

/* spheno.c */
int spheno_cmssm(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int spheno_gmsb(double Lambda, double Mmess, double tanb, int N5, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int spheno_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);

/* suspect.c */
int suspect_cmssm(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int suspect_gmsb(double Lambda, double Mmess, double tanb, int N5, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int suspect_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int suspect_nuhm(double m0, double m12, double tanb, double A0, double mu, double mA, double mtop, double mbot, double alphas_mz, char name[]);
int suspect_mssm(double m1, double m2, double m3, double tanb, double mA, double at, double ab, double atau, double mu, double mer, double mel, double mstaul, double mstaur, double mql, double mq3l, double mqur, double mqtr, double mqdr, double mqbr, double Q, double mtop, double mbot, double alphas_mz, char name[]);

/* higgsbounds.c */
double higgsbounds(char name[], struct parameters* param);
double higgsbounds_calculator(char name[]);

/* flha.c */
void flha_generator(char name[], char name_output[]);

/* nmssmtools.c */
int nmssmtools_cnmssm(double m0, double m12, double tanb, double A0, double lambda, double AK, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int nmssmtools_nnuhm(double m0, double m12, double tanb, double A0, double MHDGUT, double MHUGUT, double lambda, double AK, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int nmssmtools_ngmsb(double Lambda, double Mmess, double tanb, int N5, double lambda, double AL, double Del_h, double sgnmu, double mtop, double mbot, double alphas_mz, char name[]);
int NMSSM_collider_excluded(char name[]);
int NMSSM_theory_excluded(char name[]);
int NMSSM_upsilon_excluded(char name[]);
int NMSSM_etab_excluded(char name[]);

/* bsll.c */
double complex g_bsll(double z, double s);
double Rcchad(double s, struct parameters* param);
double complex g_bsll_parametrized(double z, double s, struct parameters* param);
double glambda_bsll(double z);
double grho_bsll(double z);
double f_bsll(double x);
double h_bsll(double z);
double kappa_bsll(double x, double alfas);
double sigma_bsll(double s);
double sigma7_bsll(double s, double L);
double sigma9_bsll(double s);
double f7_bsll(double s);
double f9_bsll(double s);
double complex Gm1_bsll(double t);
double complex G0_bsll(double t);
double complex Di23_bsll(double s, double w, double z);
double complex Di27_bsll(double s, double w, double z);
double tau77_bsll(double s);
double tau99_bsll(double s);
double tau79_bsll(double s);
double tau710_bsll(double s);
double tau910_bsll(double s);
double tau22_bsll(double w, double s, double z);
double integ_tau22(double s, double z);
double tau78_bsll(double s);
double tau88_bsll(double s);
double tau89_bsll(double s);
double complex tau27_bsll(double w, double s, double z);
double complex integ_tau27(double s, double z);
double complex tau28_bsll(double w, double s, double z);
double complex integ_tau28(double s, double z);
double complex tau29_bsll(double w, double s, double z);
double complex integ_tau29(double s, double z);
double complex tau810_bsll(double s);
double complex dtau210_bsll(double w, double s, double z);
double complex tau210_bsll(double s, double z);
double complex F_bsll(double r);
double complex F87_bsll(double s, double L);
double F89_bsll(double s);
double w77em(double s, double L);
double w79em(double s, double L);
double w710em(double s, double L);
double w99em(double s, double L);
double w910em(double s, double L);
double w1010em(double s, double L);
double w22em(double s, double L, double mub);
double complex w27em(double s, double L, double mub);
double complex w29em(double s, double L, double mub);
double complex w210em(double s, double L, double a, double mub);
double dBR_BXsmumu_dshat(double shat,  double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], struct parameters* param, double mu_b);
double dBR_BXsmumu_dshat_calculator(double shat, char name[]);
double BRBXsmumu_lowq2(double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], struct parameters* param, double mu_b);
double BRBXsmumu_highq2(double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], struct parameters* param, double mu_b);
double BRBXsmumu_lowq2_calculator(char name[]);
double BRBXsmumu_highq2_calculator(char name[]);
double A_BXsmumu(double shat, double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], struct parameters* param, double mu_b);
double A_BXsmumu_lowq2(double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], struct parameters* param, double mu_b);
double A_BXsmumu_highq2(double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], struct parameters* param, double mu_b);
double A_BXsmumu_lowq2_calculator(char name[]);
double A_BXsmumu_highq2_calculator(char name[]);
double A_BXsmumu_zero(double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], struct parameters* param, double mu_b);
double A_BXsmumu_zero_calculator(char name[]);
double dBR_BXstautau_dshat(double shat,  double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], struct parameters* param, double mu_b);
double dBR_BXstautau_dshat_calculator(double shat, char name[]);
double BRBXstautau_highq2(double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], struct parameters* param, double mu_b);
double BRBXstautau_highq2_calculator(char name[]);
double A_BXstautau(double shat, double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], struct parameters* param, double mu_b);
double A_BXstautau_highq2(double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], struct parameters* param, double mu_b);
double A_BXstautau_highq2_calculator(char name[]);

/* bsll_extra.c */
double complex F17_bsll(double s, double z, double L);
double complex F27_bsll(double s, double z, double L);
double complex F19_bsll(double s, double z, double L);
double complex F29_bsll(double s, double z, double L);
double F_17re(double muh, double z, double sh, int maxpow);
double F_17im(double muh, double z, double sh, int maxpow);
double F_19re(double muh, double z, double sh, int maxpow);
double F_19im(double muh, double z, double sh, int maxpow);
double F_27re(double muh, double z, double sh, int maxpow);
double F_27im(double muh, double z, double sh, int maxpow);
double F_29re(double muh, double z, double sh, int maxpow);
double F_29im(double muh, double z, double sh, int maxpow);
double DeltaF_29re(double muh, double z, double sh, int maxpow);
double DeltaF_29im(double muh, double z, double sh, int maxpow);
double DeltaF_19re(double muh, double z, double sh, int maxpow);
double DeltaF_19im(double muh, double z, double sh, int maxpow);

/* bkll.c */
double complex h_bkll(double q2, double mq, double mu);
double phi_Kstar(double u, double a1, double a2);
double complex B0_bkll(double s, double mq);
double complex L1_bkll(double complex x);
double complex I1_bkll(double u, double mq, double q2, struct parameters* param);
double complex tperp_bkll(double u, double mq, double q2, double E_Kstar, struct parameters* param);
double complex tpar_bkll(double u, double mq, double q2, double E_Kstar, struct parameters* param);
double dGamma_BKstarmumu_dq2(double q2, double obs[][3], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double dGamma_BKstarmumu_dq2_calculator(double q2, double obs[][3], char name[]);
double BRBKstarmumu(double smin, double smax, double obs[], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKstarmumu_lowq2(double obs[], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKstarmumu_highq2(double obs[], double C0b[], double C1b[], double C2b[], double complex CQ0b[], double complex CQ1b[], double Cpb[], double complex CQpb[], struct parameters* param, double mu_b);
double BRBKstarmumu_lowq2_calculator(char name[]);
double BRBKstarmumu_highq2_calculator(char name[]);
double A_BKstarmumu_lowq2_calculator(char name[]);
double A_BKstarmumu_highq2_calculator(char name[]);
double FL_BKstarmumu_lowq2_calculator(char name[]);
double FL_BKstarmumu_highq2_calculator(char name[]);
double FT_BKstarmumu_lowq2_calculator(char name[]);
double FT_BKstarmumu_highq2_calculator(char name[]);
double AT1_BKstarmumu_lowq2_calculator(char name[]);
double AT1_BKstarmumu_highq2_calculator(char name[]);
double AT2_BKstarmumu_lowq2_calculator(char name[]);
double AT2_BKstarmumu_highq2_calculator(char name[]);
double AT3_BKstarmumu_lowq2_calculator(char name[]);
double AT3_BKstarmumu_highq2_calculator(char name[]);
double AT4_BKstarmumu_lowq2_calculator(char name[]);
double AT4_BKstarmumu_highq2_calculator(char name[]);
double AT5_BKstarmumu_lowq2_calculator(char name[]);
double AT5_BKstarmumu_highq2_calculator(char name[]);
double HT1_BKstarmumu_lowq2_calculator(char name[]);
double HT1_BKstarmumu_highq2_calculator(char name[]);
double HT2_BKstarmumu_lowq2_calculator(char name[]);
double HT2_BKstarmumu_highq2_calculator(char name[]);
double HT3_BKstarmumu_lowq2_calculator(char name[]);
double HT3_BKstarmumu_highq2_calculator(char name[]);
double alpha_BKstarmumu_lowq2_calculator(char name[]);
double alpha_BKstarmumu_highq2_calculator(char name[]);
double AIm_BKstarmumu_lowq2_calculator(char name[]);
double AIm_BKstarmumu_highq2_calculator(char name[]);
double P2_BKstarmumu_lowq2_calculator(char name[]);
double P2_BKstarmumu_highq2_calculator(char name[]);
double P3_BKstarmumu_lowq2_calculator(char name[]);
double P3_BKstarmumu_highq2_calculator(char name[]);
double P6_BKstarmumu_lowq2_calculator(char name[]);
double P6_BKstarmumu_highq2_calculator(char name[]);
double BRobs_BKstarmumu_lowq2_calculator(char name[],double obs[]);
double BRobs_BKstarmumu_highq2_calculator(char name[],double obs[]);
double A_BKstarmumu_zero_calculator(char name[]);
double dAI_BKstarmumu_dq2(double q2, double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b);
double AI_BKstarmumu(double smin, double smax, double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b);
double AI_BKstarmumu_lowq2(double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b);
double AI_BKstarmumu_highq2(double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b);
double AI_BKstarmumu_zero(double C0b[], double C1b[], double C2b[], struct parameters* param, double mu_b);
double AI_BKstarmumu_lowq2_calculator(char name[]);
double AI_BKstarmumu_highq2_calculator(char name[]);
double AI_BKstarmumu_zero_calculator(char name[]);
