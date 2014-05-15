#include "include.h" 
#include "isajet.h"

/*--------------------------------------------------------------------*/

int isajet_cmssm(double m0, double m12, double tanb, double A0, double sgnmu, double mtop, char name[])
/* generates a SLHA file for a CMSSM parameter space point using ISAJET */
{
	FILE *tmp,*tmp2;
	char tmp_char[300],namedir[300];
	char *curdir;
	int dummy;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.istmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"%s.is1",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"%s.is2\n",name);
	fprintf(tmp,"%s\n",name);
	fprintf(tmp,"/\n");
	fprintf(tmp,"1\n");
	fprintf(tmp,"%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n",m0,m12,A0,tanb,sgnmu,mtop);
	fprintf(tmp,"0\n");
	fprintf(tmp,"/\n");
	fclose(tmp);

	sprintf(tmp_char,"%s < %s.is1 > %s.is3",ISAJET,name,name);
	system(tmp_char);

	sprintf(tmp_char,"mv %s %s/",name,curdir);
	system(tmp_char);

	if(test_file("ISALHD.out"))
	{
		tmp=fopen("ISALHD.out","r");
		sprintf(tmp_char,"%s/%s",curdir,name);
		tmp2=fopen(tmp_char,"a");
		while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
		fclose(tmp);
		fclose(tmp2);
	}

	sprintf(tmp_char,"rm -f %s.is1 %s.is2 %s.is3 ISALHD.out",name,name,name);
	system(tmp_char);
	
	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}

/*--------------------------------------------------------------------*/

int isajet_gmsb(double Lambda, double Mmess, double tanb, int N5, double cGrav, double sgnmu, double mtop, char name[])
/* generates a SLHA file for a GMSB parameter space point using ISAJET */
{
	FILE *tmp,*tmp2;
	char tmp_char[300],namedir[300];
	char *curdir;
	int dummy;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.istmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"%s.is1",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"%s.is2\n",name);
	fprintf(tmp,"%s\n",name);
	fprintf(tmp,"/\n");
	fprintf(tmp,"2\n");
	fprintf(tmp,"%.5e,%.5e,%d,%.5e,%.5e,%.5e,%.5e\n", Lambda, Mmess, N5, tanb, sgnmu, mtop, cGrav);
	fprintf(tmp,"0\n");
	fprintf(tmp,"/\n");
	fclose(tmp);

	sprintf(tmp_char,"%s < %s.is1 > %s.is3",ISAJET,name,name);
	system(tmp_char);

	sprintf(tmp_char,"mv %s %s/",name,curdir);
	system(tmp_char);

	if(test_file("ISALHD.out"))
	{
		tmp=fopen("ISALHD.out","r");
		sprintf(tmp_char,"%s/%s",curdir,name);
		tmp2=fopen(tmp_char,"a");
		while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
		fclose(tmp);
		fclose(tmp2);
	}

	sprintf(tmp_char,"rm -f %s.is1 %s.is2 %s.is3 ISALHD.out",name,name,name);
	system(tmp_char);
	
	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}

/*--------------------------------------------------------------------*/

int isajet_amsb(double m0, double m32, double tanb, double sgnmu, double mtop, char name[])
/* generates a SLHA file for a AMSB parameter space point using ISAJET */
{
	FILE *tmp,*tmp2;
	char tmp_char[300],namedir[300];
	char *curdir;
	int dummy;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.istmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"%s.is1",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"%s.is2\n",name);
	fprintf(tmp,"%s\n",name);
	fprintf(tmp,"/\n");
	fprintf(tmp,"7\n");
	fprintf(tmp,"%.5e,%.5e,%.5e,%.5e,%.5e\n", m0, m32, tanb, sgnmu, mtop);
	fprintf(tmp,"0\n");
	fprintf(tmp,"/\n");
	fclose(tmp);

	sprintf(tmp_char,"%s < %s.is1 > %s.is3",ISAJET,name,name);
	system(tmp_char);

	sprintf(tmp_char,"mv %s %s/",name,curdir);
	system(tmp_char);

	if(test_file("ISALHD.out"))
	{
		tmp=fopen("ISALHD.out","r");
		sprintf(tmp_char,"%s/%s",curdir,name);
		tmp2=fopen(tmp_char,"a");
		while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
		fclose(tmp);
		fclose(tmp2);
	}

	sprintf(tmp_char,"rm -f %s.is1 %s.is2 %s.is3 ISALHD.out",name,name,name);
	system(tmp_char);
	
	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}
/*--------------------------------------------------------------------*/

int isajet_nuhm(double m0, double m12, double tanb, double A0, double mu, double mA, double mtop, char name[])
/* generates a SLHA file for a NUHM parameter space point using ISAJET */
{
	FILE *tmp,*tmp2;
	char tmp_char[300],namedir[300];
	char *curdir;
	int dummy;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.istmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"%s.is1",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"%s.is2\n",name);
	fprintf(tmp,"%s\n",name);
	fprintf(tmp,"/\n");
	fprintf(tmp,"3\n");
	fprintf(tmp,"%.5e,%.5e,%.5e,%.5e,%.5e,%.5e\n",m0,m12,A0,tanb,mu/fabs(mu),mtop);
	fprintf(tmp,"8\n");
	fprintf(tmp,"%.5e,%.5e\n",mu,mA);
	fprintf(tmp,"0\n");
	fprintf(tmp,"0\n");
	fprintf(tmp,"/\n");
	fclose(tmp);

	sprintf(tmp_char,"%s < %s.is1 > %s.is3",ISAJET,name,name);
	system(tmp_char);

	sprintf(tmp_char,"mv %s %s/",name,curdir);
	system(tmp_char);

	if(test_file("ISALHD.out"))
	{
		tmp=fopen("ISALHD.out","r");
		sprintf(tmp_char,"%s/%s",curdir,name);
		tmp2=fopen(tmp_char,"a");
		while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
		fclose(tmp);
		fclose(tmp2);
	}

	sprintf(tmp_char,"rm -f %s.is1 %s.is2 %s.is3 ISALHD.out",name,name,name);
	system(tmp_char);
	
	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}

/*--------------------------------------------------------------------*/

int isajet_mmamsb(double alpha, double m32, double tanb, double sgnmu, double mtop, char name[])
/* generates a SLHA file for a Mixed-Moduli AMSB parameter space point using ISAJET */
{
	FILE *tmp,*tmp2;
	char tmp_char[300],namedir[300];
	char *curdir;
	int dummy;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.istmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"%s.is1",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"%s.is2\n",name);
	fprintf(tmp,"%s\n",name);
	fprintf(tmp,"/\n");
	fprintf(tmp,"9\n");
	fprintf(tmp,"%.5e,%.5e,%.5e,%.5e,%.5e\n", alpha, m32, tanb, sgnmu, mtop);
	fprintf(tmp,"/\n");
	fprintf(tmp,"/\n");
	fprintf(tmp,"0\n");
	fprintf(tmp,"/\n");
	fclose(tmp);

	sprintf(tmp_char,"%s < %s.is1 > %s.is3",ISAJET,name,name);
	system(tmp_char);

	sprintf(tmp_char,"mv %s %s/",name,curdir);
	system(tmp_char);

	if(test_file("ISALHD.out"))
	{
		tmp=fopen("ISALHD.out","r");
		sprintf(tmp_char,"%s/%s",curdir,name);
		tmp2=fopen(tmp_char,"a");
		while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
		fclose(tmp);
		fclose(tmp2);
	}

	sprintf(tmp_char,"rm -f %s.is1 %s.is2 %s.is3 ISALHD.out",name,name,name);
	system(tmp_char);
	
	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}

/*--------------------------------------------------------------------*/

int isajet_hcamsb(double alpha, double m32, double tanb, double sgnmu, double mtop, char name[])
/* generates a SLHA file for a Hypercharged AMSB parameter space point using ISAJET */
{
	FILE *tmp,*tmp2;
	char tmp_char[300],namedir[300];
	char *curdir;
	int dummy;
	curdir=getcwd(NULL, 500);

	sprintf(namedir,"%s.istmp",name);

	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namedir);
 	system(tmp_char);
	
	chdir(namedir);

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);

	sprintf(tmp_char,"%s.is1",name);

	tmp=fopen(tmp_char,"w");
	fprintf(tmp,"%s.is2\n",name);
	fprintf(tmp,"%s\n",name);
	fprintf(tmp,"/\n");
	fprintf(tmp,"10\n");
	fprintf(tmp,"%.5e,%.5e,%.5e,%.5e,%.5e\n", alpha, m32, tanb, sgnmu, mtop);
	fprintf(tmp,"0\n");
	fprintf(tmp,"/\n");
	fclose(tmp);

	sprintf(tmp_char,"%s < %s.is1 > %s.is3",ISAJET,name,name);
	system(tmp_char);

	sprintf(tmp_char,"mv %s %s/",name,curdir);
	system(tmp_char);

	if(test_file("ISALHD.out"))
	{
		tmp=fopen("ISALHD.out","r");
		sprintf(tmp_char,"%s/%s",curdir,name);
		tmp2=fopen(tmp_char,"a");
		while((dummy=getc(tmp))!=EOF) putc(dummy,tmp2);
		fclose(tmp);
		fclose(tmp2);
	}

	sprintf(tmp_char,"rm -f %s.is1 %s.is2 %s.is3 ISALHD.out",name,name,name);
	system(tmp_char);
	
	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namedir);
 	system(tmp_char);

	return 1;
}
