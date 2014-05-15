#include "include.h"
#include "nmssmtools.h"

int nmssmtools_cnmssm(double m0, double m12, double tanb, double A0, double lambda, double AK, double sgnmu, double mtop, double mbot, double alphas_mz, char name[])
/* generates a SLHA file for a CNMSSM parameter space point using NMSSMtools */
{
	FILE *tmp,*tmp2;
	char tmp_char[500],namedir[300];
	char *curdir;
	int dummy;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.nmtmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);

	sprintf(tmp_char,"mkdir -p %s/%s",namedir,namedir);
 	system(tmp_char);

	chdir(namedir);

	sprintf(tmp_char,"ln -s %s/../EXPCON .",NMSSMTools);
	system(tmp_char);

	chdir(namedir);

	tmp=fopen("inp","w");
	fprintf(tmp,"BLOCK MODSEL		     # Select model\n");
	fprintf(tmp,"    1    1		     # CNMSSM\n");
	fprintf(tmp,"    3    1		     # NMSSM PARTICLE CONTENT\n");
	fprintf(tmp,"    9    0		# FLAG FOR MICROMEGAS (0=NO, 1=YES)\n");	
	fprintf(tmp,"    10    0		# 1 POINT\n");	
	fprintf(tmp,"BLOCK SMINPUTS		     # Standard Model inputs\n");
	fprintf(tmp,"    1	1.279340000e+02	     # alpha^(-1) SM MSbar(MZ)\n");
	fprintf(tmp,"    2      1.166370000e-05	     # G_Fermi\n");
	fprintf(tmp,"    3      %.10e	     # alpha_s(MZ) SM MSbar\n",alphas_mz);
	fprintf(tmp,"    4      9.118760000e+01	     # MZ(pole)\n");
	fprintf(tmp,"    5	%.10e	     # mb(mb) SM MSbar\n",mbot);
	fprintf(tmp,"    6      %.10e	     # mtop(pole)\n",mtop);
	fprintf(tmp,"    7	1.777000000e+00	     # mtau(pole)\n");
	fprintf(tmp,"BLOCK MINPAR		     # Input parameters\n");
	fprintf(tmp,"    1      %.10e	     # m0\n",m0);
	fprintf(tmp,"    2      %.10e	     # m12\n",m12);
	fprintf(tmp,"    3      %.10e	     # tanb\n",tanb);
	fprintf(tmp,"    4      %d	     # sign(mu)\n",(int)sgnmu);
	fprintf(tmp,"    5      %.10e	     # A0\n",A0);
	fprintf(tmp,"BLOCK EXTPAR		     \n");
	fprintf(tmp,"    61	%.10e		# L\n",lambda);
	fprintf(tmp,"    64	%.10e		# AK\n",AK);
	fclose(tmp);

	sprintf(tmp_char,"%s/nmspec > nmssmtools.out",NMSSMTools);
	system(tmp_char);
	
	system("rm -f inp nmssmtools.out");
	
	tmp=fopen("spectr","r");
	sprintf(tmp_char,"%s/%s",curdir,name);
	tmp2=fopen(tmp_char,"w");
	while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
	fclose(tmp);
	fclose(tmp2);
	
	tmp=fopen("decay","r");
	sprintf(tmp_char,"%s/%s",curdir,name);
	tmp2=fopen(tmp_char,"a");
	while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
	fclose(tmp);
	fclose(tmp2);
	
	chdir(curdir);
	
	sprintf(tmp_char,"rm -rf %s.nmtmp",name);
 	system(tmp_char);
			
	return 1;
}

/*--------------------------------------------------------------------*/

int nmssmtools_nnuhm(double m0, double m12, double tanb, double A0, double MHDGUT, double MHUGUT, double lambda, double AK, double sgnmu, double mtop, double mbot, double alphas_mz, char name[])
/* generates a SLHA file for a NNUHM parameter space point using NMSSMtools */
{
	FILE *tmp,*tmp2;
	char tmp_char[500],namedir[300];
	char *curdir;
	int dummy;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.nmtmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);

	sprintf(tmp_char,"mkdir -p %s/%s",namedir,namedir);
 	system(tmp_char);

	chdir(namedir);

	sprintf(tmp_char,"ln -s %s/../EXPCON .",NMSSMTools);
	system(tmp_char);

	chdir(namedir);

	tmp=fopen("inp","w");
	fprintf(tmp,"BLOCK MODSEL		     # Select model\n");
	fprintf(tmp,"    1    1		     # CNMSSM\n");
	fprintf(tmp,"    3    1		     # NMSSM PARTICLE CONTENT\n");
	fprintf(tmp,"    9    0		# FLAG FOR MICROMEGAS (0=NO, 1=YES)\n");	
	fprintf(tmp,"    10    0		# 1 POINT\n");	
	fprintf(tmp,"BLOCK SMINPUTS		     # Standard Model inputs\n");
	fprintf(tmp,"    1	1.279340000e+02	     # alpha^(-1) SM MSbar(MZ)\n");
	fprintf(tmp,"    2      1.166370000e-05	     # G_Fermi\n");
	fprintf(tmp,"    3      %.10e	     # alpha_s(MZ) SM MSbar\n",alphas_mz);
	fprintf(tmp,"    4      9.118760000e+01	     # MZ(pole)\n");
	fprintf(tmp,"    5	%.10e	     # mb(mb) SM MSbar\n",mbot);
	fprintf(tmp,"    6      %.10e	     # mtop(pole)\n",mtop);
	fprintf(tmp,"    7	1.777000000e+00	     # mtau(pole)\n");
	fprintf(tmp,"BLOCK MINPAR		     # Input parameters\n");
	fprintf(tmp,"    1      %.10e	     # m0\n",m0);
	fprintf(tmp,"    2      %.10e	     # m12\n",m12);
	fprintf(tmp,"    3      %.10e	     # tanb\n",tanb);
	fprintf(tmp,"    4      %d	     # sign(mu)\n",(int)sgnmu);
	fprintf(tmp,"    5      %.10e	     # A0\n",A0);
	fprintf(tmp,"BLOCK EXTPAR		     \n");
	fprintf(tmp,"    61	%.10e		# L\n",lambda);
	fprintf(tmp,"    64	%.10e		# AK\n",AK);
	fprintf(tmp,"    21	%.10e		# MHDGUT\n",MHDGUT);
	fprintf(tmp,"    22	%.10e		# MHUGUT\n",MHUGUT);
	fclose(tmp);

	sprintf(tmp_char,"%s/nmspec > nmssmtools.out",NMSSMTools);
	system(tmp_char);
	
	system("rm -f inp nmssmtools.out");
	
	tmp=fopen("spectr","r");
	sprintf(tmp_char,"%s/%s",curdir,name);
	tmp2=fopen(tmp_char,"w");
	while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
	fclose(tmp);
	fclose(tmp2);
	
	tmp=fopen("decay","r");
	sprintf(tmp_char,"%s/%s",curdir,name);
	tmp2=fopen(tmp_char,"a");
	while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
	fclose(tmp);
	fclose(tmp2);
	
	chdir(curdir);
	
	sprintf(tmp_char,"rm -rf %s.nmtmp",name);
 	system(tmp_char);
			
	return 1;
}

/*--------------------------------------------------------------------*/

int nmssmtools_ngmsb(double Lambda, double Mmess, double tanb, int N5, double lambda, double AL, double Del_h, double sgnmu, double mtop, double mbot, double alphas_mz, char name[])
/* generates a SLHA file for a NGMSB parameter space point using NMSSMtools */
{
	FILE *tmp,*tmp2;
	char tmp_char[500],namedir[300];
	char *curdir;
	int dummy;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.nmtmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);

	sprintf(tmp_char,"mkdir -p %s/%s",namedir,namedir);
 	system(tmp_char);

	chdir(namedir);

	sprintf(tmp_char,"ln -s %s/../EXPCON .",NMSSMTools);
	system(tmp_char);

	chdir(namedir);

	tmp=fopen("inp","w");
	fprintf(tmp,"BLOCK MODSEL		     # Select model\n");
	fprintf(tmp,"    1    2		     # NGMSB\n");
	fprintf(tmp,"    3    1		     # NMSSM PARTICLE CONTENT\n");
	fprintf(tmp,"    9    0		# FLAG FOR MICROMEGAS (0=NO, 1=YES)\n");	
	fprintf(tmp,"    10    0		# 1 POINT\n");	
	fprintf(tmp,"BLOCK SMINPUTS		     # Standard Model inputs\n");
	fprintf(tmp,"    1	1.279340000e+02	     # alpha^(-1) SM MSbar(MZ)\n");
	fprintf(tmp,"    2      1.166370000e-05	     # G_Fermi\n");
	fprintf(tmp,"    3      %.10e	     # alpha_s(MZ) SM MSbar\n",alphas_mz);
	fprintf(tmp,"    4      9.118760000e+01	     # MZ(pole)\n");
	fprintf(tmp,"    5	%.10e	     # mb(mb) SM MSbar\n",mbot);
	fprintf(tmp,"    6      %.10e	     # mtop(pole)\n",mtop);
	fprintf(tmp,"    7	1.777000000e+00	     # mtau(pole)\n");
	fprintf(tmp,"BLOCK MINPAR		     # Input parameters\n");
	fprintf(tmp,"    1      %.10e	     # Lambda\n",Lambda);
	fprintf(tmp,"    2      %.10e	     # Mmess\n",Mmess);
	fprintf(tmp,"    3      %.10e	     # tanb\n",tanb);
	fprintf(tmp,"    4      %d	     # sign(mu)\n",(int)sgnmu);
	fprintf(tmp,"    5      %d	     # N5\n",N5);
	fprintf(tmp,"BLOCK EXTPAR		     \n");
	fprintf(tmp,"    61	%.10e		# L\n",lambda);
	fprintf(tmp,"    63	%.10e		# AL\n",AL);
	fprintf(tmp,"    66      0.D0		# XiF at M_mess in GeV**2\n");
	fprintf(tmp,"    67	0.D0 		# XiS at M_mess in GeV**3 (If MS is not an input)\n");
	fprintf(tmp,"    68	0.D0		# Mu' at M_mess in GeV\n");
	fprintf(tmp,"    69	0.D0		# m_s'**2 at M_mess in GeV**2\n");
	fprintf(tmp,"    71	%.10e		# Del_h, Shift in m_H at M_mess\n",Del_h);

	fclose(tmp);

	sprintf(tmp_char,"%s/nmgmsb > nmssmtools.out",NMSSMTools);
	system(tmp_char);
	
	system("rm -f inp nmssmtools.out");
	
	tmp=fopen("spectr","r");
	sprintf(tmp_char,"%s/%s",curdir,name);
	tmp2=fopen(tmp_char,"w");
	while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
	fclose(tmp);
	fclose(tmp2);
	
	tmp=fopen("decay","r");
	sprintf(tmp_char,"%s/%s",curdir,name);
	tmp2=fopen(tmp_char,"a");
	while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
	fclose(tmp);
	fclose(tmp2);
	
	chdir(curdir);
	
	sprintf(tmp_char,"rm -rf %s.nmtmp",name);
 	system(tmp_char);
			
	return 1;
}

/*--------------------------------------------------------------------*/

int NMSSM_collider_excluded(char name[])
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return param.NMSSMcoll;
}

/*--------------------------------------------------------------------*/

int NMSSM_theory_excluded(char name[])
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return param.NMSSMtheory;
}

/*--------------------------------------------------------------------*/

int NMSSM_upsilon_excluded(char name[])
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return param.NMSSMups1S;
}

/*--------------------------------------------------------------------*/

int NMSSM_etab_excluded(char name[])
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return param.NMSSMetab1S;
}
