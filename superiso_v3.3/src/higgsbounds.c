#include "include.h"
#include "higgsbounds.h"


double higgsbounds(char name[], struct parameters* param)
/* function calling HiggsBounds using the SLHA file "name" and the "param" structure already read in leshouches.c */
{
	int ie;
	char dum[5],tmp_char[500],namemod[200],namemod2[200],dummy[500];
	FILE *input,*output;
	char *curdir;
	double hb;
	curdir=getcwd(NULL, 500);

	sprintf(namemod,"%s.hbtmp",name);
	sprintf(namemod2,"%s/tmp.slha",namemod);
	
	sprintf(tmp_char,"rm -rf %s",namemod);
 	system(tmp_char);
	
	sprintf(tmp_char,"mkdir -p %s",namemod);
 	system(tmp_char);
	
	input=fopen(name,"r");
	output=fopen(namemod2,"w");
		
	while(EOF != fscanf(input,"%c",dummy))
	{
			if(!strncasecmp("\n",dummy,1)) 
			{
				fprintf(output,"%c",dummy[0]);
				if(EOF != fscanf(input,"%c",dummy))
	
			if(!strncasecmp("d",dummy,1)) 
			{
				dum[0]=dummy[0];
				fscanf(input,"%c",dummy);
				
			if(!strncasecmp("e",dummy,1)) 
			{
				dum[1]=dummy[0];
				fscanf(input,"%c",dummy);
			
			if(!strncasecmp("c",dummy,1)) 
			{
				dum[2]=dummy[0];
				fscanf(input,"%c",dummy);
			
			if(!strncasecmp("a",dummy,1)) 
			{
				dum[3]=dummy[0];
				fscanf(input,"%c",dummy);
							
			if(!strncasecmp("y",dummy,1)) 
			{
				dum[4]=dummy[0];
				
				if(!strncasecmp("decay",dum,5)) 
				{
					while(EOF != fscanf(input,"%c",dummy))
					if(!strncasecmp("\n",dummy,1))
					{
						fscanf(input,"%c",dummy);
		
					if(!strncasecmp("b",dummy,1))
					{
						dum[0]=dummy[0];
						fscanf(input,"%c",dummy);
						
					if(!strncasecmp("l",dummy,1))
					{
						dum[1]=dummy[0];
						fscanf(input,"%c",dummy);
	
					if(!strncasecmp("o",dummy,1))
					{
						dum[2]=dummy[0];
						fscanf(input,"%c",dummy);
						
					if(!strncasecmp("c",dummy,1))
					{
						dum[3]=dummy[0];
						fscanf(input,"%c",dummy);
	
					if(!strncasecmp("k",dummy,1))
					{
						dum[4]=dummy[0];
					
					if(!strncasecmp("block",dum,5))
					{
						fprintf(output,"%s",dum);
						break;
					}
					}
					}
					}
					}
					}
					}
				}
				else fprintf(output,"%s",dum); 
			}
			} else fprintf(output,"%c%c%c%c",dum[0],dum[1],dum[2],dummy[0]);
			} else fprintf(output,"%c%c%c",dum[0],dum[1],dummy[0]);
			} else fprintf(output,"%c%c",dum[0],dummy[0]);
			} else fprintf(output,"%c",dummy[0]);
			} else fprintf(output,"%c",dummy[0]);
		}
	
		fclose(input);
		fclose(output);		
	
	chdir(namemod);
	
	sprintf(namemod2,"%s > %s.tmp",HBwithFH,name);
	system(namemod2);
	
	sprintf(namemod2,"%s.tmp",name);
	
	output=fopen(namemod2,"r");
	
	while((EOF != fscanf(output,"%s",dummy))&&(strcasecmp("of",dummy)));
	
	if(!strcasecmp("of",dummy)) 
	{
		fscanf(output,"%s",dummy);
		hb=atof(dummy);
	}
	else hb=-1.;
	
	fclose(output);
	
	chdir(curdir);
	sprintf(tmp_char,"rm -rf %s",namemod);
 	system(tmp_char);
	
	return hb;
}


double higgsbounds_calculator(char name[])
/* "container" function scanning the SLHA file "name" and calling HiggsBounds */
{
	struct parameters param;

	Init_param(&param);

	if(!Les_Houches_Reader(name,&param)) return -1.;
	return higgsbounds(name, &param);
}



