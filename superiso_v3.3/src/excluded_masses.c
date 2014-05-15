#include "include.h"

int excluded_Higgs_masses(struct parameters* param)
/* tests whether the SUSY point is excluded by the collider contraints */
/* if excluded, return 1, otherwise 0 */
{
	int excluded=0;

	if(param->SM==1) return excluded;

	if(param->mass_h0!=0.) excluded=(excluded||(fabs(param->mass_h0)<111.)); /* Higgs */
	if(param->mass_H!=0.) excluded=(excluded||(fabs(param->mass_H)<79.3)); /* charged Higgs */
	if(param->mass_A0!=0.) excluded=(excluded||(fabs(param->mass_A0)<93.4)); /* CP-odd Higgs */

	return excluded;
}

/*--------------------------------------------------------------------*/

int excluded_SUSY_masses(struct parameters* param)
/* tests whether the SUSY point is excluded by the collider contraints */
/* if excluded, return 1, otherwise 0 */
{
	int excluded=0;

	if(param->SM==1) return excluded;

#ifdef SM_ChargedHiggs	
	return excluded;
#endif

	if(param->mass_neut[1]!=0.) excluded=(excluded||(fabs(param->mass_neut[1])<46.)); /* neutralino 1 */
	if(param->mass_neut[2]!=0.) if(param->tan_beta<40.) excluded=(excluded||(fabs(param->mass_neut[2])<62.4)); /* neutralino 2 */
	if(param->mass_neut[3]!=0.) if(param->tan_beta<40.) excluded=(excluded||(fabs(param->mass_neut[3])<99.9)); /* neutralino 3 */
	if(param->mass_neut[4]!=0.) if(param->tan_beta<40.) excluded=(excluded||(fabs(param->mass_neut[4])<116.)); /* neutralino 4 */
	if(param->mass_cha1!=0.) if((param->tan_beta<40.)&&(fabs(param->mass_cha1)-fabs(param->mass_neut[1])>5.)) excluded=(excluded||(fabs(param->mass_cha1)<94.)); /* chargino */
	if(param->mass_er!=0.) excluded=(excluded||(fabs(param->mass_er)<73.)); /* selectron R */
	if(param->mass_el!=0.) excluded=(excluded||(fabs(param->mass_el)<107.)); /* selectron L */
	if(param->mass_mur!=0.) if((param->tan_beta<40.)&&(fabs(param->mass_mur)-fabs(param->mass_neut[1])>10.)) excluded=(excluded||(fabs(param->mass_mur)<94.)); /* smuon R*/
	if(param->mass_mul!=0.) if((param->tan_beta<40.)&&(fabs(param->mass_mur)-fabs(param->mass_neut[1])>10.)) excluded=(excluded||(fabs(param->mass_mul)<94.)); /* smuon L*/
	if(param->mass_nuel!=0.) if((param->tan_beta<40.)&&(fabs(param->mass_er)-fabs(param->mass_neut[1])>10.)) excluded=(excluded||(fabs(param->mass_nuel)<94.)); /* sneutrino */
	if(param->mass_numl!=0.) if((param->tan_beta<40.)&&(fabs(param->mass_mur)-fabs(param->mass_neut[1])>10.)) excluded=(excluded||(fabs(param->mass_numl)<94.)); /* sneutrino */
	if(param->mass_nutl!=0.) if((param->tan_beta<40.)&&(fabs(param->mass_tau1)-fabs(param->mass_neut[1])>10.)) excluded=(excluded||(fabs(param->mass_nutl)<94.)); /* sneutrino */
	if(param->mass_dnr!=0.)  if(fabs(param->mass_dnr)-fabs(param->mass_neut[1])>10.) excluded=(excluded||(fabs(param->mass_dnr)<100.)); /* squark d_R */
	if(param->mass_dnl!=0.)  if(fabs(param->mass_dnl)-fabs(param->mass_neut[1])>10.) excluded=(excluded||(fabs(param->mass_dnl)<100.)); /* squark d_L */
	if(param->mass_upr!=0.)  if(fabs(param->mass_upr)-fabs(param->mass_neut[1])>10.) excluded=(excluded||(fabs(param->mass_upr)<100.)); /* squark u_R */
	if(param->mass_upl!=0.)  if(fabs(param->mass_upl)-fabs(param->mass_neut[1])>10.) excluded=(excluded||(fabs(param->mass_upl)<100.)); /* squark u_L */
	if(param->mass_t1!=0.) if(fabs(param->mass_t1)-fabs(param->mass_neut[1])>10.) excluded=(excluded||(fabs(param->mass_t1)<95.7)); /* stop */	
	if(param->mass_gluino!=0.) excluded=(excluded||(fabs(param->mass_gluino)<195.)); /* gluino */
	if(param->mass_b1!=0.)  if(fabs(param->mass_b1)-fabs(param->mass_neut[1])>5.) excluded=(excluded||(fabs(param->mass_b1)<100.)); /* sbottom */
	if(param->mass_tau1!=0.) if(fabs(param->mass_tau1)-fabs(param->mass_neut[1])>15.) excluded=(excluded||(fabs(param->mass_tau1)<81.9)); /* stau */
	
	return excluded;
}

/*--------------------------------------------------------------------*/

int excluded_masses(struct parameters* param)
{
	return excluded_Higgs_masses(param)||excluded_SUSY_masses(param);
}

/*--------------------------------------------------------------------*/

int excluded_Higgs_mass_calculator(char name[])
/* "container" function scanning the SLHA file "name" and checking if the SUSY point is excluded by the Higgs mass contraints */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return excluded_Higgs_masses(&param);
}
/*--------------------------------------------------------------------*/

int excluded_mass_calculator(char name[])
/* "container" function scanning the SLHA file "name" and checking if the SUSY point is excluded by the mass contraints */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return excluded_masses(&param);
}

/*--------------------------------------------------------------------*/

int excluded_SUSY_mass_calculator(char name[])
/* "container" function scanning the SLHA file "name" and checking if the SUSY point is excluded by the SUSY mass contraints */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return excluded_SUSY_masses(&param);
}

/*--------------------------------------------------------------------*/

int charged_LSP(struct parameters* param)
/* tests whether the SUSY point corresponds to a charged LSP (NLSP if the LSP is a gravitino) */
/* if the LSP is charged, return 1, otherwise 0 */
{
	double neutral_mass;
	int charged_LSP=0;
	int test_sneutrinoL=0;
	int test_gluino=0;
	
#ifdef SMONLY	
	return charged_LSP;
#endif

#ifdef SM_ChargedHiggs	
	return charged_LSP;
#endif

	neutral_mass=fabs(param->mass_neut[1]);
	
	if(param->mass_nuel!=0.) 
	{
		neutral_mass=min(fabs(param->mass_nuel),neutral_mass);
		test_sneutrinoL=((test_sneutrinoL)||(neutral_mass==fabs(param->mass_nuel)));
	}
	if(param->mass_numl!=0.) 
	{
		neutral_mass=min(fabs(param->mass_numl),neutral_mass);
		test_sneutrinoL=((test_sneutrinoL)||(neutral_mass==fabs(param->mass_numl)));
	}
	if(param->mass_nutl!=0.) 
	{
		neutral_mass=min(fabs(param->mass_nutl),neutral_mass);
		test_sneutrinoL=((test_sneutrinoL)||(neutral_mass==fabs(param->mass_nutl)));
	}
	if(param->mass_gluino!=0.)
	{
		neutral_mass=min(fabs(param->mass_gluino),neutral_mass);
		test_gluino=(neutral_mass==fabs(param->mass_gluino));
	}
	
	if(param->mass_el!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_el)<neutral_mass));
	if(param->mass_er!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_er)<neutral_mass));
	if(param->mass_mul!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_mul)<neutral_mass));
	if(param->mass_mur!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_mur)<neutral_mass));
	if(param->mass_tau1!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_tau1)<neutral_mass));
	if(param->mass_upl!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_upl)<neutral_mass));
	if(param->mass_upr!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_upr)<neutral_mass));
	if(param->mass_chl!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_chl)<neutral_mass));
	if(param->mass_chr!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_chr)<neutral_mass));
	if(param->mass_t1!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_t1)<neutral_mass));
	if(param->mass_dnl!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_dnl)<neutral_mass));
	if(param->mass_dnr!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_dnr)<neutral_mass));
	if(param->mass_stl!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_stl)<neutral_mass));
	if(param->mass_str!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_str)<neutral_mass));
	if(param->mass_b1!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_b1)<neutral_mass));
	if(param->mass_cha1!=0.) charged_LSP=(charged_LSP||(fabs(param->mass_cha1)<neutral_mass));
	
	if(!charged_LSP)
	{
		if(test_gluino) return 3;
		if(test_sneutrinoL) return 2;
	}
	return charged_LSP;
}

/*--------------------------------------------------------------------*/

int charged_LSP_calculator(char name[])
/* "container" function scanning the SLHA file "name" and checking if the SUSY point corresponds to a charged LSP */
{
	struct parameters param;
		
	Init_param(&param);
	
	if(!Les_Houches_Reader(name,&param)) return -1;

	return charged_LSP(&param);
}

