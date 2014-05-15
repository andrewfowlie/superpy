#include "include.h" 
#include "2hdmc.h"


int thdmc_types(double l1, double l2, double l3, double l4, double l5, double l6, double l7, double m12_2, double tanb, int type, char name[])
/* generates a LHA file for a 2HDM parameter space point using 2HDMC */
{
	FILE *tmp;
	char tmp_char[500];
	char *dir;

	sprintf(tmp_char,"rm -f %s",name);
	system(tmp_char);
	
	dir=getcwd(NULL, 500);

	chdir(THDMC);

	sprintf(tmp_char,"./CalcGen %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %d %s/%s > 2hdmc_%s",l1,l2,l3,l4,l5,l6,l7,m12_2,tanb,type,dir,name,name);
	system(tmp_char);

	sprintf(tmp_char,"rm -f 2hdmc_%s",name);
	system(tmp_char);

	chdir(dir);
	free(dir);

	return 1;
}

