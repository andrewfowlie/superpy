#include "include.h"
#include "softsusy.h"

int softsusy_cmssm(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, double mbot, double alphas_mz, char name[])
/* generates a SLHA file for a CMSSM parameter space point using SOFTSUSY */
{
	FILE *tmp;
	char tmp_char[200];

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"%s.ss",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"Block MODSEL		     # Select model\n");
	fprintf(tmp,"    1    1		     # cmssm\n");
	fprintf(tmp,"Block SMINPUTS		     # Standard Model inputs\n");
	fprintf(tmp,"    1	1.279340000e+02	     # alpha^(-1) SM MSbar(MZ)\n");
	fprintf(tmp,"    2      1.166370000e-05	     # G_Fermi\n");
	fprintf(tmp,"    3      %.10e	     # alphas(MZ) SM MSbar\n",alphas_mz);
	fprintf(tmp,"    4      9.118760000e+01	     # MZ(pole)\n");
	fprintf(tmp,"    5	%.10e	     # mb(mb) SM MSbar\n",mbot);
	fprintf(tmp,"    6      %.10e	     # mtop(pole)\n",mtop);
	fprintf(tmp,"    7	1.777000000e+00	     # mtau(pole)\n");
	fprintf(tmp,"Block MINPAR		     # Input parameters\n");
	fprintf(tmp,"    1      %.10e	     # m0\n",m0);
	fprintf(tmp,"    2      %.10e	     # m12\n",m12);
	fprintf(tmp,"    3      %.10e	     # tanb\n",tanb);
	fprintf(tmp,"    4      %d	     # sign(mu)\n",(int)sgnmu);
	fprintf(tmp,"    5      %.10e	     # A0\n",A0);
	fclose(tmp);

	sprintf(tmp_char,"%s leshouches < %s.ss > %s", SOFTSUSY, name, name);
	system(tmp_char);

	sprintf(tmp_char,"rm -f %s.ss",name);
	system(tmp_char);
	return 1;
}

/*--------------------------------------------------------------------*/

int softsusy_gmsb(double Lambda, double Mmess, double tanb, int N5, double cGrav, double sgnmu, double mtop, double mbot, double alphas_mz, char name[])
/* generates a SLHA file for a GMSB parameter space point using SOFTSUSY */
{
	FILE *tmp;
	char tmp_char[200];

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	if(Lambda>Mmess) 
	{
		tmp=fopen(name,"w");
		fclose(tmp);
		return 0;
	}

	sprintf(tmp_char,"%s.ss",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"Block MODSEL		     # Select model\n");
	fprintf(tmp,"    1    2		     # gmsb\n");
	fprintf(tmp,"Block SMINPUTS		     # Standard Model inputs\n");
	fprintf(tmp,"    1	1.279340000e+02	     # alpha^(-1) SM MSbar(MZ)\n");
	fprintf(tmp,"    2      1.166370000e-05	     # G_Fermi\n");
	fprintf(tmp,"    3      %.10e	     # alphas(MZ) SM MSbar\n",alphas_mz);
	fprintf(tmp,"    4      9.118760000e+01	     # MZ(pole)\n");
	fprintf(tmp,"    5	%.10e	     # mb(mb) SM MSbar\n",mbot);
	fprintf(tmp,"    6      %.10e	     # mtop(pole)\n",mtop);
	fprintf(tmp,"    7	1.777000000e+00	     # mtau(pole)\n");
	fprintf(tmp,"Block MINPAR		     # Input parameters\n");
	fprintf(tmp,"    1      %.10e	     # scale\n",Lambda);
	fprintf(tmp,"    2      %.10e	     # Mmess\n",Mmess);
	fprintf(tmp,"    3      %.10e	     # tanb\n",tanb);
	fprintf(tmp,"    4      %d	     # sign(mu)\n",(int)sgnmu);
	fprintf(tmp,"    5      %d	     # N5\n",N5);
	fprintf(tmp,"    6      %.10e	     # cGrav\n",cGrav);
	fclose(tmp);

	sprintf(tmp_char,"%s leshouches < %s.ss > %s", SOFTSUSY, name, name);
	system(tmp_char);

	sprintf(tmp_char,"rm -f %s.ss",name);
	system(tmp_char);
	return 1;
}

/*--------------------------------------------------------------------*/

int softsusy_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, double mbot, double alphas_mz, char name[])
/* generates a SLHA file for an AMSB parameter space point using SOFTSUSY */
{
	FILE *tmp;
	char tmp_char[200];

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"%s.ss",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"Block MODSEL		     # Select model\n");
	fprintf(tmp,"    1    3		     # amsb\n");
	fprintf(tmp,"Block SMINPUTS		     # Standard Model inputs\n");
	fprintf(tmp,"    1	1.279340000e+02	     # alpha^(-1) SM MSbar(MZ)\n");
	fprintf(tmp,"    2      1.166370000e-05	     # G_Fermi\n");
	fprintf(tmp,"    3      %.10e	     # alphas(MZ) SM MSbar\n",alphas_mz);
	fprintf(tmp,"    4      9.118760000e+01	     # MZ(pole)\n");
	fprintf(tmp,"    5	%.10e	     # mb(mb) SM MSbar\n",mbot);
	fprintf(tmp,"    6      %.10e	     # mtop(pole)\n",mtop);
	fprintf(tmp,"    7	1.777000000e+00	     # mtau(pole)\n");
	fprintf(tmp,"Block MINPAR		     # Input parameters\n");
	fprintf(tmp,"    1      %.10e	     # m0\n",m0);
	fprintf(tmp,"    2      %.10e	     # m12\n",m32);
	fprintf(tmp,"    3      %.10e	     # tanb\n",tanb);
	fprintf(tmp,"    4      %d	     # sign(mu)\n",(int)sgnmu);
	fclose(tmp);

	sprintf(tmp_char,"%s leshouches < %s.ss > %s", SOFTSUSY, name, name);
	system(tmp_char);

	sprintf(tmp_char,"rm -f %s.ss",name);
	system(tmp_char);
	return 1;
}

/*--------------------------------------------------------------------*/

int softsusy_nuhm(double m0, double m12, double tanb, double A0, double mu, double mA, double mtop, double mbot, double alphas_mz, char name[])
/* generates a SLHA file for NUHM parameter space point using SOFTSUSY */
{
	FILE *tmp;
	char tmp_char[200];

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"%s.ss",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"Block MODSEL		     # Select model\n");
	fprintf(tmp,"    1    1		     # cmssm\n");
	fprintf(tmp,"Block SMINPUTS		     # Standard Model inputs\n");
	fprintf(tmp,"    1	1.279340000e+02	     # alpha^(-1) SM MSbar(MZ)\n");
	fprintf(tmp,"    2      1.166370000e-05	     # G_Fermi\n");
	fprintf(tmp,"    3      %.10e	     # alphas(MZ) SM MSbar\n",alphas_mz);
	fprintf(tmp,"    4      9.118760000e+01	     # MZ(pole)\n");
	fprintf(tmp,"    5	%.10e	     # mb(mb) SM MSbar\n",mbot);
	fprintf(tmp,"    6      %.10e	     # mtop(pole)\n",mtop);
	fprintf(tmp,"    7	1.777000000e+00	     # mtau(pole)\n");
	fprintf(tmp,"Block MINPAR		     # Input parameters\n");
	fprintf(tmp,"    1      %.10e	     # m0\n",m0);
	fprintf(tmp,"    2      %.10e	     # m12\n",m12);
	fprintf(tmp,"    3      %.10e	     # tanb\n",tanb);
	fprintf(tmp,"    5      %.10e	     # A0\n",A0);
	fprintf(tmp,"Block EXTPAR		     # External Input parameters\n");
	fprintf(tmp,"    23     %.10e	     # mu\n",mu);
	fprintf(tmp,"    26     %.10e	     # mA pole\n",mA);
	fclose(tmp);

	sprintf(tmp_char,"%s leshouches < %s.ss > %s", SOFTSUSY, name, name);
	system(tmp_char);

	sprintf(tmp_char,"rm -f %s.ss",name);
	system(tmp_char);
	return 1;
}

/*--------------------------------------------------------------------*/

int softsusy_mssm(double m1, double m2, double m3, double tanb, double mA, double at, double ab, double atau, double mu, double mer, double mel, double mstaul, double mstaur, double mql, double mq3l, double mqur, double mqtr, double mqdr, double mqbr, double Q, double mtop, double mbot, double alphas_mz, char name[])
/* generates a SLHA file for a pMSSM parameter space point using SOFTSUSY */
{
	FILE *tmp;
	char tmp_char[200];

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"%s.ss",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"Block MODSEL		     # Select model\n");
	fprintf(tmp,"    1    0		     #MSSM\n");
	fprintf(tmp,"Block SMINPUTS		     # Standard Model inputs\n");
	fprintf(tmp,"    1	1.279340000e+02	     # alpha^(-1) SM MSbar(MZ)\n");
	fprintf(tmp,"    2      1.166370000e-05	     # G_Fermi\n");
	fprintf(tmp,"    3      %.10e	     # alphas(MZ) SM MSbar\n",alphas_mz);
	fprintf(tmp,"    4      9.118760000e+01	     # MZ(pole)\n");
	fprintf(tmp,"    5	%.10e	     # mb(mb) SM MSbar\n",mbot);
	fprintf(tmp,"    6      %.10e	     # mtop(pole)\n",mtop);
	fprintf(tmp,"    7	1.777000000e+00	     # mtau(pole)\n");
	fprintf(tmp,"Block MINPAR		     # SUSY input parameters\n");
	fprintf(tmp,"    3      %.10e	     # tanb\n",tanb);
	fprintf(tmp,"Block EXTPAR		     # Input parameters\n");
	fprintf(tmp,"    0      %.10e	     # Q\n",Q);
	fprintf(tmp,"    1      %.10e	     # M1\n",m1);
	fprintf(tmp,"    2      %.10e	     # M2\n",m2);
	fprintf(tmp,"    3      %.10e	     # M3\n",m3);
	fprintf(tmp,"   11      %.10e	     # At\n",at);
	fprintf(tmp,"   12      %.10e	     # Ab\n",ab);
	fprintf(tmp,"   13      %.10e	     # Atau\n",atau);
	fprintf(tmp,"   23      %.10e	     # Mu\n",mu);
	fprintf(tmp,"   26      %.10e	     # MA\n",mA);
	fprintf(tmp,"   31      %.10e	     # MeL\n",mel);
	fprintf(tmp,"   32      %.10e	     # MmuL\n",mel);
	fprintf(tmp,"   33      %.10e	     # MstauL\n",mstaul);
	fprintf(tmp,"   34      %.10e	     # MeR\n",mer);
	fprintf(tmp,"   35      %.10e	     # MmuR\n",mer);
	fprintf(tmp,"   36      %.10e	     # MstauR\n",mstaur);
	fprintf(tmp,"   41      %.10e	     # Mq1L\n",mql);
	fprintf(tmp,"   42      %.10e	     # Mq2L\n",mql);
	fprintf(tmp,"   43      %.10e	     # Mq3L\n",mq3l);
	fprintf(tmp,"   44      %.10e	     # MquR\n",mqur);
	fprintf(tmp,"   45      %.10e	     # MqcR\n",mqur);
	fprintf(tmp,"   46      %.10e	     # MqtR\n",mqtr);
	fprintf(tmp,"   47      %.10e	     # MqdR\n",mqdr);
	fprintf(tmp,"   48      %.10e	     # MqsR\n",mqdr);
	fprintf(tmp,"   49      %.10e	     # MqbR\n",mqbr);

	fclose(tmp);

	sprintf(tmp_char,"%s leshouches < %s.ss > %s", SOFTSUSY, name, name);
	system(tmp_char);

	sprintf(tmp_char,"rm -f %s.ss",name);
	system(tmp_char);
	return 1;
}

/*--------------------------------------------------------------------*/

int softsusy_slhain(char name_in[], char name_out[])
/* generates a SLHA file from a SLHA input file using SOFTSUSY */
{
	char tmp_char[500];
	
	sprintf(tmp_char,"%s leshouches < %s > %s", SOFTSUSY, name_in, name_out);
	system(tmp_char);

	return 1;
}

