#include "include.h"
#include "spheno.h"

int spheno_cmssm(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, double mbot, double alphas_mz, char name[])
/* generates a SLHA file for a CMSSM parameter space point using SPHENO */
{
	FILE *tmp;
	char tmp_char[300],namedir[300];
	char *curdir;
	curdir=getcwd(NULL, 500);
	
	sprintf(namedir,"%s.sptmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);

	tmp=fopen("LesHouches.in","w");
	
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

	sprintf(tmp_char,"%s > spheno.out", SPHENO);
	system(tmp_char);

	sprintf(tmp_char,"mv -f SPheno.spc %s/%s",curdir,name);
	system(tmp_char);
	
	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}

/*--------------------------------------------------------------------*/

int spheno_gmsb(double Lambda, double Mmess, double tanb, int N5, double sgnmu, double mtop, double mbot, double alphas_mz, char name[])
/* generates a SLHA file for a GMSB parameter space point using SPHENO */
{
	FILE *tmp;
	char tmp_char[300],namedir[300];
	char *curdir;
	curdir=getcwd(NULL, 500);
	
	sprintf(namedir,"%s.sptmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	if(Lambda>Mmess) 
	{
		tmp=fopen(name,"w");
		fclose(tmp);
		return 0;
	}
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);

	tmp=fopen("LesHouches.in","w");
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
	fclose(tmp);

	sprintf(tmp_char,"%s > spheno.out", SPHENO);
	system(tmp_char);

	sprintf(tmp_char,"mv -f SPheno.spc %s/%s",curdir,name);
	system(tmp_char);
	
	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}

/*--------------------------------------------------------------------*/

int spheno_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, double mbot, double alphas_mz, char name[])
/* generates a SLHA file for an AMSB parameter space point using SPHENO */
{
	FILE *tmp;
	char tmp_char[300],namedir[300];
	char *curdir;
	curdir=getcwd(NULL, 500);
	
	sprintf(namedir,"%s.sptmp",name);
	
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);

	tmp=fopen("LesHouches.in","w");
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
	fprintf(tmp,"    2      %.10e	     # m32\n",m32);
	fprintf(tmp,"    3      %.10e	     # tanb\n",tanb);
	fprintf(tmp,"    4      %d	     # sign(mu)\n",(int)sgnmu);
	fclose(tmp);

	sprintf(tmp_char,"%s > spheno.out", SPHENO);
	system(tmp_char);
	
	sprintf(tmp_char,"mv -f SPheno.spc %s/%s",curdir,name);
	system(tmp_char);
	
	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}
