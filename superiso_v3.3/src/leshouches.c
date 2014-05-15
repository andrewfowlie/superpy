#include "include.h"

void Init_param(struct parameters* param)
/* initializes the "param" structure by initializing the parameters with either 0 or a value from the PDG2010 */
{
	int ie,je;
	
	param->SM=0;
	param->model=-3; /* this parameter is used to test whether the scan of the SLHA file succeeds. */
	param->generator=0;
	param->Q=0.;
	param->m0=0.;
	param->m12=0.;
	param->tan_beta=0.;
	param->sign_mu=0.;
	param->A0=0.;
	param->mass_W=0.;
	param->Lambda=0.;
	param->Mmess=0.;
	param->N5=0.;
	param->cgrav=0.;
	param->m32=0.;
	param->mass_Z=0.;
	param->mass_b=0.;
	param->mass_top_pole=0.;
	param->mass_tau_pole=0.;
	param->inv_alpha_em=0.;
	param->alphas_MZ=0.;
	param->alpha=0.;
	param->Gfermi=0.;
	param->GAUGE_Q=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++)
	{
		param->charg_Umix[ie][je]=0.;
		param->charg_Vmix[ie][je]=0.;
		param->stop_mix[ie][je]=0.;
		param->sbot_mix[ie][je]=0.;
		param->stau_mix[ie][je]=0.;
	}
	for(ie=1;ie<=5;ie++) for(je=1;je<=5;je++) param->neut_mix[ie][je]=0.;
	for(ie=1;ie<=5;ie++) param->mass_neut[ie]=0.;
	for(ie=1;ie<=3;ie++) param->yub[ie]=param->yut[ie]=param->yutau[ie]=0.;
		
	param->Min=0.;
	param->M1_Min=0.;
	param->M2_Min=0.;
	param->M3_Min=0.;
	param->At_Min=0.;
	param->Ab_Min=0.;
	param->Atau_Min=0.;
	param->M2H1_Min=0.;
	param->M2H2_Min=0.;
	param->mu_Min=0.;
	param->M2A_Min=0.;
	param->tb_Min=0.;
	param->mA_Min=0.;
	param->MeL_Min=0.;
	param->MmuL_Min=0.;
	param->MtauL_Min=0.;
	param->MeR_Min=0.;
	param->MmuR_Min=0.;
	param->MtauR_Min=0.;
	param->MqL1_Min=0.;
	param->MqL2_Min=0.;
	param->MqL3_Min=0.;
	param->MuR_Min=0.;
	param->McR_Min=0.;
	param->MtR_Min=0.;
	param->MdR_Min=0.;
	param->MsR_Min=0.;
	param->MbR_Min=0.;
	param->N51=0.;
	param->N52=0.;
	param->N53=0.;
	param->M2H1_Q=0.;
	param->M2H2_Q=0.;
	param->mass_h0=0.;
	param->mass_H0=0.;
	param->mass_A0=0.;
	param->mass_H=0.;
	param->mass_dnl=0.;
	param->mass_upl=0.;
	param->mass_stl=0.;
	param->mass_chl=0.;
	param->mass_b1=0.;
	param->mass_t1=0.;
	param->mass_el=0.;
	param->mass_nuel=0.;
	param->mass_mul=0.;
	param->mass_numl=0.;
	param->mass_tau1=0.;
	param->mass_nutl=0.;
	param->mass_gluino=0.;
	param->mass_cha1=0.;
	param->mass_cha2=0.;
	param->mass_dnr=0.;
	param->mass_upr=0.;
	param->mass_str=0.;
	param->mass_chr=0.;
	param->mass_b2=0.;
	param->mass_t2=0.;
	param->mass_er=0.;
	param->mass_mur=0.;
	param->mass_tau2=0.;
	param->gp=0.;
	param->g2=0.;
	param->gp_Q=0.;
	param->g2_Q=0.;
	param->g3_Q=0.;
	param->YU_Q=0.;
	param->YD_Q=0.;
	param->YE_Q=0.;
	param->HMIX_Q=0.;
	param->mu_Q=0.;
	param->tanb_GUT=0.;
	param->Higgs_VEV=0.;
	param->mA2_Q=0.;
	param->MSOFT_Q=0.;
	param->M1_Q=0.;
	param->M2_Q=0.;
	param->M3_Q=0.;
	param->MeL_Q=0.;
	param->MmuL_Q=0.;
	param->MtauL_Q=0.;
	param->MeR_Q=0.;
	param->MmuR_Q=0.;
	param->MtauR_Q=0.;
	param->MqL1_Q=0.;
	param->MqL2_Q=0.;
	param->MqL3_Q=0.;
	param->MuR_Q=0.;
	param->McR_Q=0.;
	param->MtR_Q=0.;
	param->MdR_Q=0.;
	param->MsR_Q=0.;
	param->MbR_Q=0.;
	param->AU_Q=0.;
	param->A_u=0.;
	param->A_c=0.;
	param->A_t=0.;
	param->AD_Q=0.;
	param->A_d=0.;
	param->A_s=0.;
	param->A_b=0.;
	param->AE_Q=0.;
	param->A_e=0.;
	param->A_mu=0.;
	param->A_tau=0.;
	param->mass_graviton=0.;
	param->mass_gravitino=0.;
	param->mass_nuer=0.;
	param->mass_numr=0.;
	param->mass_nutr=0.;
	param->mass_t=0.;
	param->mass_tau=0.;
	param->mass_gluon=0.;
	param->mass_nue=0.;
	param->mass_num=0.;
	param->mass_nut=0.;
	param->mass_photon=0.;
	param->mass_Z0=0.;

	/* SLHA2 */
	param->NMSSM=0;
	param->RV=0;
	param->CPV=0;
	param->FV=0;
	param->mass_nutau2=0.;
	param->mass_e2=0.;
	param->mass_nue2=0.;
	param->mass_mu2=0.;
	param->mass_numu2=0.;
	param->mass_d2=0.;
	param->mass_u2=0.;
	param->mass_s2=0.;
	param->mass_c2=0.;
	param->CKM_lambda=0.;
	param->CKM_A=0.;
	param->CKM_rhobar=0.;
	param->CKM_etabar=0.;
	param->PMNS_theta12=0.;
	param->PMNS_theta23=0.;
	param->PMNS_theta13=0.;
	param->PMNS_delta13=0.;
	param->PMNS_alpha1=0.;
	param->PMNS_alpha2=0.;
	param->lambdaNMSSM_Min=0.;
	param->kappaNMSSM_Min=0.;
	param->AlambdaNMSSM_Min=0.;
	param->AkappaNMSSM_Min=0.;
	param->lambdaSNMSSM_Min=0.;
	param->xiFNMSSM_Min=0.;
	param->xiSNMSSM_Min=0.;
	param->mupNMSSM_Min=0.;
	param->mSp2NMSSM_Min=0.;
	param->mS2NMSSM_Min=0.;
	param->mass_H03=0.;
	param->mass_A02=0.;
	param->NMSSMRUN_Q=0.;
	param->lambdaNMSSM=0.;
	param->kappaNMSSM=0.;
	param->AlambdaNMSSM=0.;
	param->AkappaNMSSM=0.;
	param->lambdaSNMSSM=0.;
	param->xiFNMSSM=0.;
	param->xiSNMSSM=0.;
	param->mupNMSSM=0.;
	param->mSp2NMSSM=0.;
	param->mS2NMSSM=0.;
	param->PMNSU_Q=0.;
	param->CKM_Q=0.;
	param->IMCKM_Q=0.;
	param->MSE2_Q=0.;
	param->MSU2_Q=0.;
	param->MSD2_Q=0.;
	param->MSL2_Q=0.;
	param->MSQ2_Q=0.;
	param->TU_Q=0.;
	param->TD_Q=0.;
	param->TE_Q=0.;
	param->NMSSMcoll=0;
	param->NMSSMtheory=0;
	param->NMSSMups1S=0;
	param->NMSSMetab1S=0;
	
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		param->sU_mix[ie][je]=0.;
		param->sD_mix[ie][je]=0.;
		param->sE_mix[ie][je]=0.;
	}
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++)
	{
		param->H0_mix[ie][je]=0.;
		param->A0_mix[ie][je]=0.;
		param->sNU_mix[ie][je]=0.;
		param->sCKM_msq2[ie][je]=0.;
		param->sCKM_msl2[ie][je]=0.;
		param->sCKM_msd2[ie][je]=0.;
		param->sCKM_msu2[ie][je]=0.;
		param->sCKM_mse2[ie][je]=0.;
		param->PMNS_U[ie][je]=0.;
		param->CKM[ie][je]=0.;
		param->IMCKM[ie][je]=0.;
		param->TU[ie][je]=0.;
		param->TD[ie][je]=0.;
		param->TE[ie][je]=0.;
	}
	
	/* non-SLHA*/
	param->mass_b_1S=0.;
	param->mass_b_pole=0.;
	param->mtmt=0.;
	param->Lambda5=0.;

	/* 2HDM */
	param->THDM_model=0;
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++)
	{
		param->lambda_u[ie][je]=0.;
		param->lambda_d[ie][je]=0.;
		param->lambda_l[ie][je]=0.;
	}
	
	/* Flavour physics */
	param->f_B=0.194;
	param->f_Bs=0.234;
	param->f_Ds=0.248;
	param->f_D=0.207;
	param->fK_fpi=1.193;
	param->f_K_par=0.220;
	param->f_K_perp=0.185;
	param->m_B=5.27917;
	param->m_Bs=5.3663;
	param->m_Bd=5.27950;
	param->m_pi=0.1396;
	param->m_K=0.4937;
	param->m_Kstar=0.89594;
	param->m_D0=1.86483;
	param->m_D=1.86960;
	param->m_Ds=1.96847;
	param->life_pi=2.6033e-8;
	param->life_K=1.2380e-8;  
	param->life_B=1.641e-12;
	param->life_Bs=1.472e-12;
	param->life_Bd=1.519e-12;
	param->life_D=1.040e-12;
	param->life_Ds=5.e-13;
	param->a1perp=0.10;
	param->a2perp=0.13;
	param->a1par=0.10;
	param->a2par=0.09;
	param->zeta3A=0.032;
	param->zeta3V=0.013;
	param->wA10=-2.1;
	param->deltatp=0.16;
	param->deltatm=-0.16;
	param->lambda_Bp=0.46;
	param->rho1=0.06;
	param->lambda2=0.12;
	
	/* CKM matrix */
	param->Vud=0.97428;
	param->Vus=0.2253;
	param->Vub=0.00347;
	param->Vcd=-0.2252;
	param->Vcs=0.97345;
	param->Vcb=0.0410;
	param->Vtd=0.00862;
	param->Vts=-0.0403;
	param->Vtb=0.999152;
	param->CKM[1][1]=param->Vud;
	param->CKM[1][2]=param->Vus;
	param->CKM[1][3]=param->Vub;
	param->CKM[2][1]=param->Vcd;
	param->CKM[2][2]=param->Vcs;
	param->CKM[2][3]=param->Vcb;
	param->CKM[3][1]=param->Vtd;
	param->CKM[3][2]=param->Vts;
	param->CKM[3][3]=param->Vtb;
	
	/* masses and coupling from PDG 2011 */
	param->mass_u = 2.49e-3;
	param->mass_d = 5.05e-3;
	param->mass_s = 0.101;
	param->mass_c = 1.29;
	param->mass_b = 4.19;
	param->mass_top_pole = 172.9;
	
	param->mass_e = 0.511e-3;
	param->mass_mu= 0.105658;
	param->mass_tau_pole=1.7768;
	param->mass_tau=param->mass_tau_pole;
	
	param->mass_Z=91.1876;
	param->alphas_MZ=0.1184;
	param->mass_W=80.399;

	param->gp=param->gp_Q=3.57458e-1;
	param->g2=param->g2_Q=6.51908e-1;	
	param->inv_alpha_em=1.27916e2;
	param->Gfermi=1.16637000e-5;

	param->width_Z=2.495;
	param->width_W=2.085;

	return;
}

/*--------------------------------------------------------------------*/

int Les_Houches_Reader(char name[], struct parameters* param)
/* reads the SLHA file "name" and puts all the values into the "param" structure */
{
	FILE *lecture;
	char dummy[500],dummy2[500],nda[500],id1[500],id2[500];
	double version,Qtemp;
	int ie,je;
	
	if(!test_file(name))
	{
		param->model=-4;
		return 0;
	}
	
	lecture = fopen(name,"r");
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"MODSEL"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block")))
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 0: fscanf(lecture,"%d",&param->THDM_model); break;
					case 1:	fscanf(lecture,"%d",&param->model); break;
					case 3:	fscanf(lecture,"%d",&param->NMSSM); break;
					case 4:	fscanf(lecture,"%d",&param->RV); break;
					case 5:	fscanf(lecture,"%d",&param->CPV); break;
					case 6:	fscanf(lecture,"%d",&param->FV); break;
					case 12: fscanf(lecture,"%lf",&param->Q); break;
				}	
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
	}
	fclose(lecture);		

	if(param->NMSSM != 0) param->model=param->NMSSM; 
	if(param->RV != 0) param->model=-2;
	if(param->CPV != 0) param->model=-2;
	if(param->THDM_model !=0) param->model=param->THDM_model;
	
	if(param->model<0) return 0;	
	
	lecture = fopen(name,"r");
	while(EOF != fscanf(lecture,"%s",dummy))
	{
		if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
		if(!strcasecmp(dummy,"SPINFO"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: 	fscanf(lecture,"%s",dummy); 
							if(!strncasecmp(dummy,"ISA",3)) param->generator=1; 
							if(!strncasecmp(dummy,"SOFTSUSY",8)) param->generator=3; 
							if(!strncasecmp(dummy,"SPHENO",6)) param->generator=4; 
							if(!strncasecmp(dummy,"SUSPECT",7)) param->generator=5; 
							if(!strncasecmp(dummy,"NMSSMTOOLS",10)) param->generator=6; 
							if(!strncasecmp(dummy,"2HDMC",5)) param->generator=10; 
							break;
					case 2: if(param->generator==1) 
						{
							fscanf(lecture,"%lf",&version); 
							if(version>=7.80) param->generator=2;
						}
						break;
						
					case 3: if(param->generator==6)
						{
							if(EOF!=fscanf(lecture,"%s",dummy)) sprintf(dummy2,"%s",dummy);									while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) sprintf(dummy2,"%s%s",dummy2,dummy);
							if(!strcasecmp("# Chargino too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Neutralinos too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Charged Higgs too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strncasecmp("# Excluded by ee",dummy2,16)) {param->NMSSMcoll=1; break;}
							if(!strncasecmp("# Excluded by stop",dummy2,18)) {param->NMSSMcoll=1; break;}
							if(!strncasecmp("# Excluded by sbottom",dummy2,21)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Squark/gluino too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Selectron/smuon too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Stau too light",dummy2)) {param->NMSSMcoll=1; break;}
							if(!strcasecmp("# Landau Pole below MGUT",dummy2)) {param->NMSSMtheory=1; break;}
							if(!strcasecmp("# Unphysical global minimum",dummy2)) {param->NMSSMtheory=1; break;}
							if(!strcasecmp("# Higgs soft masses >> Msusy",dummy2)) {param->NMSSMtheory=1; break;}
							if(!strncasecmp("# Excluded by Upsilon",dummy2,21)) {param->NMSSMups1S=1; break;}
							if(!strncasecmp("# Excluded etab(1S)",dummy2,19)) {param->NMSSMetab1S=1; break;}
						}
						break;
						
					case 4: if(param->generator==3)
						{
							if(EOF!=fscanf(lecture,"%s",dummy)) sprintf(dummy2,"%s",dummy);											while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) sprintf(dummy2,"%s%c",dummy2,dummy[0]);
							if(strcasecmp("Point invalid: stau LSP",dummy2)) {param->model=-1; fclose(lecture); return 0;}
						}
						else
						{
							param->model=-1; fclose(lecture); return 0;
						}
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"SMINPUTS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: fscanf(lecture,"%lf",&param->inv_alpha_em); break;
					case 2: fscanf(lecture,"%lf",&param->Gfermi); break;
					case 3: fscanf(lecture,"%lf",&param->alphas_MZ); break;
					case 4: fscanf(lecture,"%lf",&param->mass_Z); break;
					case 5: fscanf(lecture,"%lf",&param->mass_b); break;
					case 6: fscanf(lecture,"%lf",&param->mass_top_pole); break;
					case 7: fscanf(lecture,"%lf",&param->mass_tau_pole); break;
					case 8: fscanf(lecture,"%lf",&param->mass_nutau2); break;
					case 11: fscanf(lecture,"%lf",&param->mass_e2); break;
					case 12: fscanf(lecture,"%lf",&param->mass_nue2); break;
					case 13: fscanf(lecture,"%lf",&param->mass_mu2); break;
					case 14: fscanf(lecture,"%lf",&param->mass_numu2); break;
					case 21: fscanf(lecture,"%lf",&param->mass_d2); break;
					case 22: fscanf(lecture,"%lf",&param->mass_u2); break;
					case 23: fscanf(lecture,"%lf",&param->mass_s2); break;
					case 24: fscanf(lecture,"%lf",&param->mass_c2); break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"VCKMIN"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: fscanf(lecture,"%lf",&param->CKM_lambda); break;
					case 2: fscanf(lecture,"%lf",&param->CKM_A); break;
					case 3: fscanf(lecture,"%lf",&param->CKM_rhobar); break;
					case 4: fscanf(lecture,"%lf",&param->CKM_etabar); break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"UPMNSIN"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: fscanf(lecture,"%lf",&param->PMNS_theta12); break;
					case 2: fscanf(lecture,"%lf",&param->PMNS_theta23); break;
					case 3: fscanf(lecture,"%lf",&param->PMNS_theta13); break;
					case 4: fscanf(lecture,"%lf",&param->PMNS_delta13); break;
					case 5: fscanf(lecture,"%lf",&param->PMNS_alpha1); break;
					case 6: fscanf(lecture,"%lf",&param->PMNS_alpha2); break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MINPAR"))
		{
			switch(param->model)
			{
				case 1:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
							if(test_integer(dummy)) switch(atoi(dummy))
							{
								case 1: fscanf(lecture,"%lf",&param->m0); break;
								case 2: fscanf(lecture,"%lf",&param->m12); break;
								case 3: fscanf(lecture,"%lf",&param->tan_beta); break;
								case 4: fscanf(lecture,"%lf",&param->sign_mu); break;
								case 5: fscanf(lecture,"%lf",&param->A0); break;
							}
					}
					break;
				}
	
				case 2:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 1: fscanf(lecture,"%lf",&param->Lambda); break;
							case 2: fscanf(lecture,"%lf",&param->Mmess); break;
							case 3: fscanf(lecture,"%lf",&param->tan_beta); break;
							case 4: fscanf(lecture,"%lf",&param->sign_mu); break;
							case 5: fscanf(lecture,"%lf",&param->N5); break;
							case 6: fscanf(lecture,"%lf",&param->cgrav); break;
						}
					}
					break;
				}
	
				case 3:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 1: fscanf(lecture,"%lf",&param->m32); break;
							case 2: fscanf(lecture,"%lf",&param->m0); break;
							case 3: fscanf(lecture,"%lf",&param->tan_beta); break;
							case 4: fscanf(lecture,"%lf",&param->sign_mu); break;
						}
					}
					break;
				}

				default:
				{
					while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
					{
						if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 3: fscanf(lecture,"%lf",&param->tan_beta); break;
						}
					}
					break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"EXTPAR"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 0: fscanf(lecture,"%lf",&param->Min); break;
					case 1: fscanf(lecture,"%lf",&param->M1_Min); break;
					case 2: fscanf(lecture,"%lf",&param->M2_Min); break;
					case 3: fscanf(lecture,"%lf",&param->M3_Min); break;	
					case 11: fscanf(lecture,"%lf",&param->At_Min); break;
					case 12: fscanf(lecture,"%lf",&param->Ab_Min); break;
					case 13: fscanf(lecture,"%lf",&param->Atau_Min); break;
					case 21: fscanf(lecture,"%lf",&param->M2H1_Min); break;
					case 22: fscanf(lecture,"%lf",&param->M2H2_Min); break;
					case 23: fscanf(lecture,"%lf",&param->mu_Min); break;
					case 24: fscanf(lecture,"%lf",&param->M2A_Min); break;
					case 25: fscanf(lecture,"%lf",&param->tb_Min); break;
					case 26: fscanf(lecture,"%lf",&param->mA_Min); break;
					case 31: fscanf(lecture,"%lf",&param->MeL_Min); break;
					case 32: fscanf(lecture,"%lf",&param->MmuL_Min); break;
					case 33: fscanf(lecture,"%lf",&param->MtauL_Min); break;
					case 34: fscanf(lecture,"%lf",&param->MeR_Min); break;
					case 35: fscanf(lecture,"%lf",&param->MmuR_Min); break;
					case 36: fscanf(lecture,"%lf",&param->MtauR_Min); break;
					case 41: fscanf(lecture,"%lf",&param->MqL1_Min); break;
					case 42: fscanf(lecture,"%lf",&param->MqL2_Min); break;
					case 43: fscanf(lecture,"%lf",&param->MqL3_Min); break;
					case 44: fscanf(lecture,"%lf",&param->MuR_Min); break;
					case 45: fscanf(lecture,"%lf",&param->McR_Min); break;
					case 46: fscanf(lecture,"%lf",&param->MtR_Min); break;
					case 47: fscanf(lecture,"%lf",&param->MdR_Min); break;
					case 48: fscanf(lecture,"%lf",&param->MsR_Min); break;
					case 49: fscanf(lecture,"%lf",&param->MbR_Min); break;
					case 51: fscanf(lecture,"%lf",&param->N51); break;
					case 52: fscanf(lecture,"%lf",&param->N52); break;
					case 53: fscanf(lecture,"%lf",&param->N53); break;
					case 61: fscanf(lecture,"%lf",&param->lambdaNMSSM_Min); break;
					case 62: fscanf(lecture,"%lf",&param->kappaNMSSM_Min); break;
					case 63: fscanf(lecture,"%lf",&param->AlambdaNMSSM_Min); break;
					case 64: fscanf(lecture,"%lf",&param->AkappaNMSSM_Min); break;
					case 65: fscanf(lecture,"%lf",&param->lambdaSNMSSM_Min); break;
					case 66: fscanf(lecture,"%lf",&param->xiFNMSSM_Min); break;
					case 67: fscanf(lecture,"%lf",&param->xiSNMSSM_Min); break;
					case 68: fscanf(lecture,"%lf",&param->mupNMSSM_Min); break;
					case 69: fscanf(lecture,"%lf",&param->mSp2NMSSM_Min); break;
					case 70: fscanf(lecture,"%lf",&param->mS2NMSSM_Min); break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MASS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: fscanf(lecture,"%lf",&param->mass_d); if(isnan(param->mass_d)) {param->model=-3; return 0;} break;
					case 2: fscanf(lecture,"%lf",&param->mass_u); if(isnan(param->mass_s)) {param->model=-3; return 0;} break;
					case 3: fscanf(lecture,"%lf",&param->mass_s); if(isnan(param->mass_s)) {param->model=-3; return 0;} break;
					case 4: fscanf(lecture,"%lf",&param->mass_c); if(isnan(param->mass_c)) {param->model=-3; return 0;} break;
					case 6: fscanf(lecture,"%lf",&param->mass_t); if(isnan(param->mass_t)) {param->model=-3; return 0;} break;
					case 11: fscanf(lecture,"%lf",&param->mass_e); if(isnan(param->mass_e)) {param->model=-3; return 0;} break;
					case 12: fscanf(lecture,"%lf",&param->mass_nue); if(isnan(param->mass_nue)) {param->model=-3; return 0;} break;
					case 13: fscanf(lecture,"%lf",&param->mass_mu); if(isnan(param->mass_mu)) {param->model=-3; return 0;} break;
					case 14: fscanf(lecture,"%lf",&param->mass_num); if(isnan(param->mass_num)) {param->model=-3; return 0;} break;
					case 15: fscanf(lecture,"%lf",&param->mass_tau); if(isnan(param->mass_tau)) {param->model=-3; return 0;} break;
					case 16: fscanf(lecture,"%lf",&param->mass_nut); if(isnan(param->mass_nut)) {param->model=-3; return 0;} break;
					case 21: fscanf(lecture,"%lf",&param->mass_gluon); break;
					case 22: fscanf(lecture,"%lf",&param->mass_photon); break;
					case 23: fscanf(lecture,"%lf",&param->mass_Z0); if(isnan(param->mass_Z0)) {param->model=-3; return 0;} break;
					case 24: fscanf(lecture,"%lf",&param->mass_W); if(isnan(param->mass_W)) {param->model=-3; return 0;} break;
					case 25: fscanf(lecture,"%lf",&param->mass_h0); if(isnan(param->mass_h0)) {param->model=-3; return 0;} break;
					case 35: fscanf(lecture,"%lf",&param->mass_H0); if(isnan(param->mass_H0)) {param->model=-3; return 0;} break;
					case 36: fscanf(lecture,"%lf",&param->mass_A0); if(isnan(param->mass_A0)) {param->model=-3; return 0;} break;
					case 37: fscanf(lecture,"%lf",&param->mass_H); if(isnan(param->mass_H)) {param->model=-3; return 0;} break;
					case 39: fscanf(lecture,"%lf",&param->mass_graviton); break;
					case 45: fscanf(lecture,"%lf",&param->mass_H03); if(isnan(param->mass_H03)) {param->model=-3; return 0;} break;
					case 46: fscanf(lecture,"%lf",&param->mass_A02); if(isnan(param->mass_A02)) {param->model=-3; return 0;} break;
					case 1000001: fscanf(lecture,"%lf",&param->mass_dnl); if(isnan(param->mass_dnl)) {param->model=-3; return 0;} break;
					case 1000002: fscanf(lecture,"%lf",&param->mass_upl); if(isnan(param->mass_upl)) {param->model=-3; return 0;} break;
					case 1000003: fscanf(lecture,"%lf",&param->mass_stl); if(isnan(param->mass_stl)) {param->model=-3; return 0;} break;
					case 1000004: fscanf(lecture,"%lf",&param->mass_chl); if(isnan(param->mass_chl)) {param->model=-3; return 0;} break;
					case 1000005: fscanf(lecture,"%lf",&param->mass_b1); if(isnan(param->mass_b1)) {param->model=-3; return 0;} break;
					case 1000006: fscanf(lecture,"%lf",&param->mass_t1); if(isnan(param->mass_t1)) {param->model=-3; return 0;} break;
					case 1000011: fscanf(lecture,"%lf",&param->mass_el); if(isnan(param->mass_el)) {param->model=-3; return 0;} break;
					case 1000012: fscanf(lecture,"%lf",&param->mass_nuel); if(isnan(param->mass_nuel)) {param->model=-3; return 0;} break;
					case 1000013: fscanf(lecture,"%lf",&param->mass_mul); if(isnan(param->mass_mul)) {param->model=-3; return 0;} break;
					case 1000014: fscanf(lecture,"%lf",&param->mass_numl); if(isnan(param->mass_numl)) {param->model=-3; return 0;} break;
					case 1000015: fscanf(lecture,"%lf",&param->mass_tau1); if(isnan(param->mass_tau1)) {param->model=-3; return 0;} break;
					case 1000016: fscanf(lecture,"%lf",&param->mass_nutl); if(isnan(param->mass_nutl)) {param->model=-3; return 0;} break;
					case 1000021: fscanf(lecture,"%lf",&param->mass_gluino); if(isnan(param->mass_gluino)) {param->model=-3; return 0;} break;
					case 1000022: fscanf(lecture,"%lf",&param->mass_neut[1]); if(isnan(param->mass_neut[1])) {param->model=-3; return 0;} break;
					case 1000023: fscanf(lecture,"%lf",&param->mass_neut[2]); if(isnan(param->mass_neut[2])) {param->model=-3; return 0;} break;
					case 1000024: fscanf(lecture,"%lf",&param->mass_cha1); if(isnan(param->mass_cha1)) {param->model=-3; return 0;} break;
					case 1000025: fscanf(lecture,"%lf",&param->mass_neut[3]); if(isnan(param->mass_neut[3])) {param->model=-3; return 0;} break;
					case 1000035: fscanf(lecture,"%lf",&param->mass_neut[4]); if(isnan(param->mass_neut[4])) {param->model=-3; return 0;} break;
					case 1000037: fscanf(lecture,"%lf",&param->mass_cha2); if(isnan(param->mass_cha2)) {param->model=-3; return 0;} break;
					case 1000039: fscanf(lecture,"%lf",&param->mass_gravitino); break;
					case 1000045: fscanf(lecture,"%lf",&param->mass_neut[5]); if(isnan(param->mass_neut[5])) {param->model=-3; return 0;} break;
					case 2000001: fscanf(lecture,"%lf",&param->mass_dnr); if(isnan(param->mass_dnr)) {param->model=-3; return 0;} break;
					case 2000002: fscanf(lecture,"%lf",&param->mass_upr); if(isnan(param->mass_upr)) {param->model=-3; return 0;} break;
					case 2000003: fscanf(lecture,"%lf",&param->mass_str); if(isnan(param->mass_str)) {param->model=-3; return 0;} break;
					case 2000004: fscanf(lecture,"%lf",&param->mass_chr); if(isnan(param->mass_chr)) {param->model=-3; return 0;} break;
					case 2000005: fscanf(lecture,"%lf",&param->mass_b2); if(isnan(param->mass_b2)) {param->model=-3; return 0;} break;
					case 2000006: fscanf(lecture,"%lf",&param->mass_t2); if(isnan(param->mass_t2)) {param->model=-3; return 0;} break;
					case 2000011: fscanf(lecture,"%lf",&param->mass_er); if(isnan(param->mass_er)) {param->model=-3; return 0;} break;
					case 2000012: fscanf(lecture,"%lf",&param->mass_nuer); break;
					case 2000013: fscanf(lecture,"%lf",&param->mass_mur); if(isnan(param->mass_mur)) {param->model=-3; return 0;} break;
					case 2000014: fscanf(lecture,"%lf",&param->mass_numr); break;
					case 2000015: fscanf(lecture,"%lf",&param->mass_tau2); if(isnan(param->mass_tau2)) {param->model=-3; return 0;} break;
					case 2000016: fscanf(lecture,"%lf",&param->mass_nutr); break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"ALPHA"))
		{
			while((EOF != fscanf(lecture,"%s",dummy))&&(strcasecmp(dummy,"Block"))) 
			{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(atoi(dummy) != 0) param->alpha = atof(dummy);
			}		
		}
		else if(!strcasecmp(dummy,"STOPMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->stop_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"SBOTMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sbot_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"STAUMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->stau_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"NMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->neut_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"NMNMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->neut_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"UMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));

				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->charg_Umix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"VMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->charg_Vmix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"GAUGE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->GAUGE_Q!=0.)&&(Qtemp>param->GAUGE_Q)) break;
					param->GAUGE_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: fscanf(lecture,"%lf",&param->gp_Q); break;
					case 2: fscanf(lecture,"%lf",&param->g2_Q); break;
					case 3: fscanf(lecture,"%lf",&param->g3_Q); break;	
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"YU"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->YU_Q!=0.)&&(Qtemp>param->YU_Q)) break;
					param->YU_Q=Qtemp;
				}
				else 				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						if(ie==atoi(dummy)) fscanf(lecture,"%lf",&param->yut[ie]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"YD"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->YD_Q!=0.)&&(Qtemp>param->YD_Q)) break;
					param->YD_Q=Qtemp;
				}
				else
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						if(ie==atoi(dummy)) fscanf(lecture,"%lf",&param->yub[ie]);
					}
				}			
			}		
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"YE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->YE_Q!=0.)&&(Qtemp>param->YE_Q)) break;
					param->YE_Q=Qtemp;
				}
				else
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						if(ie==atoi(dummy)) fscanf(lecture,"%lf",&param->yutau[ie]);
					}
				}			
			}		
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"HMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->HMIX_Q!=0.)&&(Qtemp>param->HMIX_Q)) break;
					param->HMIX_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: fscanf(lecture,"%lf",&param->mu_Q); break;
					case 2: fscanf(lecture,"%lf",&param->tanb_GUT); break;
					case 3: fscanf(lecture,"%lf",&param->Higgs_VEV); break;	
					case 4: fscanf(lecture,"%lf",&param->mA2_Q); break;
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"NMHMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->H0_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"NMAMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->A0_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSOFT"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSOFT_Q!=0.)&&(Qtemp>param->MSOFT_Q)) break;
					param->MSOFT_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: fscanf(lecture,"%lf",&param->M1_Q); break;
					case 2: fscanf(lecture,"%lf",&param->M2_Q); break;
					case 3: fscanf(lecture,"%lf",&param->M3_Q); break;	
					case 21: fscanf(lecture,"%lf",&param->M2H1_Q); break;
					case 22: fscanf(lecture,"%lf",&param->M2H2_Q); break;
					case 31: fscanf(lecture,"%lf",&param->MeL_Q); break;
					case 32: fscanf(lecture,"%lf",&param->MmuL_Q); break;
					case 33: fscanf(lecture,"%lf",&param->MtauL_Q); break;
					case 34: fscanf(lecture,"%lf",&param->MeR_Q); break;
					case 35: fscanf(lecture,"%lf",&param->MmuR_Q); break;
					case 36: fscanf(lecture,"%lf",&param->MtauR_Q); break;
					case 41: fscanf(lecture,"%lf",&param->MqL1_Q); break;
					case 42: fscanf(lecture,"%lf",&param->MqL2_Q); break;
					case 43: fscanf(lecture,"%lf",&param->MqL3_Q); break;
					case 44: fscanf(lecture,"%lf",&param->MuR_Q); break;
					case 45: fscanf(lecture,"%lf",&param->McR_Q); break;
					case 46: fscanf(lecture,"%lf",&param->MtR_Q); break;
					case 47: fscanf(lecture,"%lf",&param->MdR_Q); break;
					case 48: fscanf(lecture,"%lf",&param->MsR_Q); break;
					case 49: fscanf(lecture,"%lf",&param->MbR_Q); break;
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"AU"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->AU_Q!=0.)&&(Qtemp>param->AU_Q)) break;
					param->AU_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 1: fscanf(lecture,"%lf",&param->A_u); break;	
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 2: fscanf(lecture,"%lf",&param->A_c); break;	
						}
						break;
					}
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 3: fscanf(lecture,"%lf",&param->A_t); break;	
						}
						break;
					}
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"AD"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->AD_Q!=0.)&&(Qtemp>param->AD_Q)) break;
					param->AD_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 1: fscanf(lecture,"%lf",&param->A_d); break;	
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 2: fscanf(lecture,"%lf",&param->A_s); break;	
						}
						break;
					}
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 3: fscanf(lecture,"%lf",&param->A_b); break;	
						}
						break;
					}
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"AE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->AE_Q!=0.)&&(Qtemp>param->AE_Q)) break;
					param->AE_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 1: fscanf(lecture,"%lf",&param->A_e); break;	
						}
						break;
					}
					case 2: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 2: fscanf(lecture,"%lf",&param->A_mu); break;	
						}
						break;
					}
					case 3: 
					{
						fscanf(lecture,"%s",dummy);
						if(test_integer(dummy)) switch(atoi(dummy))
						{
							case 3: fscanf(lecture,"%lf",&param->A_tau); break;	
						}
						break;
					}
				}
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}		
		else if(!strcasecmp(dummy,"NMSSMRUN"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->NMSSMRUN_Q!=0.)&&(Qtemp>param->NMSSMRUN_Q)) break;
					param->NMSSMRUN_Q=Qtemp;
				}
				else if(test_integer(dummy)) switch(atoi(dummy))
				{
					case 1: fscanf(lecture,"%lf",&param->lambdaNMSSM); break;
					case 2: fscanf(lecture,"%lf",&param->kappaNMSSM); break;
					case 3: fscanf(lecture,"%lf",&param->AlambdaNMSSM); break;
					case 4: fscanf(lecture,"%lf",&param->AkappaNMSSM); break;
					case 5: fscanf(lecture,"%lf",&param->lambdaSNMSSM); break;
					case 6: fscanf(lecture,"%lf",&param->xiFNMSSM); break;
					case 7: fscanf(lecture,"%lf",&param->xiSNMSSM); break;
					case 8: fscanf(lecture,"%lf",&param->mupNMSSM); break;
					case 9: fscanf(lecture,"%lf",&param->mSp2NMSSM); break;
					case 10: fscanf(lecture,"%lf",&param->mS2NMSSM); break;
				}
			}
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"USQMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sU_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"DSQMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sD_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"SELMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sE_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"SELMIX"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sNU_mix[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSQ2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSQ2_Q!=0.)&&(Qtemp>param->MSQ2_Q)) break;
					param->MSQ2_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sCKM_msq2[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSL2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSL2_Q!=0.)&&(Qtemp>param->MSL2_Q)) break;
					param->MSL2_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sCKM_msl2[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSD2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSD2_Q!=0.)&&(Qtemp>param->MSD2_Q)) break;
					param->MSD2_Q=Qtemp;
				}
				if(!strcasecmp(dummy,"Q=")) fscanf(lecture,"%lf",&param->MSD2_Q);
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sCKM_msd2[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSU2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSU2_Q!=0.)&&(Qtemp>param->MSU2_Q)) break;
					param->MSU2_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sCKM_msu2[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"MSE2"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->MSE2_Q!=0.)&&(Qtemp>param->MSE2_Q)) break;
					param->MSE2_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->sCKM_mse2[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"VCKM"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->CKM_Q!=0.)&&(Qtemp>param->CKM_Q)) break;
					param->CKM_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->CKM[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"IMVCKM"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->IMCKM_Q!=0.)&&(Qtemp>param->IMCKM_Q)) break;
					param->IMCKM_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->IMCKM[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"UPMNS"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->PMNSU_Q!=0.)&&(Qtemp>param->PMNSU_Q)) break;
					param->PMNSU_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->PMNS_U[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"TU"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->TU_Q!=0.)&&(Qtemp>param->TU_Q)) break;
					param->TU_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->TU[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"TD"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->TD_Q!=0.)&&(Qtemp>param->TD_Q)) break;
					param->TD_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->TD[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"TE"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(!strncasecmp(dummy,"Q=MGUT=",7)) break;
				if(!strncasecmp(dummy,"Q=",2))
				{
					fscanf(lecture,"%lf",&Qtemp);
					if((param->TE_Q!=0.)&&(Qtemp>param->TE_Q)) break;
					param->TE_Q=Qtemp;
				}
				else 

				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->TE[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"UCOUPL"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->lambda_u[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"DCOUPL"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->lambda_d[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}
		else if(!strcasecmp(dummy,"LCOUPL"))
		{
			while((EOF != fscanf(lecture,"%s",dummy) && strcasecmp(dummy,"Block") && strcasecmp(dummy,"Decay")))
	 		{
				if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
				if(test_integer(dummy))
				{
					ie=atoi(dummy);
					fscanf(lecture,"%s",dummy);
					if(test_integer(dummy))
					{
						je=atoi(dummy);	
						fscanf(lecture,"%lf",&param->lambda_l[ie][je]);
					}
				}			
			}	
			if(!strcasecmp(dummy,"Decay")) while((!fseek(lecture,-1,SEEK_CUR))&&(EOF != fscanf(lecture,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(lecture,-1,SEEK_CUR);
		}		
	}
	fclose(lecture);

	if(param->model<0) return 0;	

	if(alphas_running(param->mass_b/2.,param->mass_top_pole, param->mass_b,param)<0.) 
	{	
		param->model=-3;
		return 0;
	}
	
	if(alphas_running(2.*param->mass_top_pole,param->mass_top_pole, param->mass_b,param)<0.)
	{	
		param->model=-3;
		return 0;
	}
	
	slha_adjust(param);
	
	return 1;	
}

/*--------------------------------------------------------------------*/

void slha_adjust(struct parameters* param)
{
	double dum;
	double mass[7],invmat[7][7];
	int ie,je;
	int iemax=0;

	if(param->mass_Z0==0.) param->mass_Z0=param->mass_Z;

	if(param->MSOFT_Q==0.) param->MSOFT_Q=max(param->TU_Q,max(param->TD_Q,max(param->TE_Q,max(param->PMNSU_Q,max(param->CKM_Q,max(param->MSE2_Q,max(param->MSU2_Q,max(param->MSD2_Q,max(param->MSL2_Q,max(param->MSQ2_Q,max(param->NMSSMRUN_Q,max(param->YU_Q,max(param->YD_Q,max(param->YE_Q,max(param->HMIX_Q,max(param->GAUGE_Q,max(param->AU_Q,max(param->AD_Q,param->AE_Q))))))))))))))))));
	
	if(param->Q==0.) param->Q=sqrt(param->mass_t1*param->mass_t2);
	
	if(param->MSOFT_Q==0.) param->TU_Q=param->TD_Q=param->TE_Q=param->PMNSU_Q=param->CKM_Q=param->MSE2_Q=param->MSU2_Q=param->MSD2_Q=param->MSL2_Q=param->MSQ2_Q=param->MSOFT_Q=param->NMSSMRUN_Q=param->YU_Q=param->YD_Q=param->YE_Q=param->HMIX_Q=param->GAUGE_Q=param->AU_Q=param->AD_Q=param->AE_Q=param->Q;
			
	if(param->tan_beta==0.) param->tan_beta=param->tanb_GUT;
	
	if((param->tan_beta*param->MSOFT_Q)==0.) param->model=-3;
	
	if(param->MeL_Q==0.) param->MeL_Q=param->mass_el;
	if(param->MmuL_Q==0.) param->MmuL_Q=param->mass_mul;
	if(param->MtauL_Q==0.) param->MtauL_Q=param->mass_tau1;
	if(param->MeR_Q==0.) param->MeR_Q=param->mass_er;
	if(param->MmuR_Q==0.) param->MmuR_Q=param->mass_mur;
	if(param->MtauR_Q==0.) param->MtauR_Q=param->mass_tau2;
	if(param->MqL1_Q==0.) param->MqL1_Q=param->mass_dnl;
	if(param->MqL2_Q==0.) param->MqL2_Q=param->mass_stl;
	if(param->MqL3_Q==0.) param->MqL3_Q=param->mass_b1;
	if(param->MuR_Q==0.) param->MuR_Q=param->mass_upr;
	if(param->McR_Q==0.) param->McR_Q=param->mass_chr;
	if(param->MtR_Q==0.) param->MtR_Q=param->mass_t2;
	if(param->MdR_Q==0.) param->MdR_Q=param->mass_dnr;
	if(param->MsR_Q==0.) param->MsR_Q=param->mass_str;
	if(param->MbR_Q==0.) param->MbR_Q=param->mass_b2;
		
	if(param->A_tau==0.) param->A_tau=param->TE[3][3];
	if(param->A_b==0.) param->A_b=param->TD[3][3];
	if(param->A_t==0.) param->A_t=param->TU[3][3];
	
 	if(param->stop_mix[1][1]==0.) 
	{
		mass[1]=param->mass_upl;
		mass[2]=param->mass_chl;
		mass[3]=param->mass_t1;
		mass[4]=param->mass_upr;
		mass[5]=param->mass_chr;
		mass[6]=param->mass_t2;

		for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) invmat[ie][je]=param->sU_mix[je][ie];
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[1][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[1][ie]);
		}
		param->mass_upl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[2][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[2][ie]);
		}
		param->mass_chl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[3][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[3][ie]);
		}
		param->mass_t2=mass[iemax];		
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[4][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[4][ie]);
		}
		param->mass_upr=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[5][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[5][ie]);
		}
		param->mass_chr=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[6][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[6][ie]);
		}
		param->mass_t1=mass[iemax];

		param->stop_mix[1][1]=param->sU_mix[1][3];
		param->stop_mix[1][2]=param->sU_mix[1][6];
	}	
			
	dum=atan(param->stop_mix[1][2]/param->stop_mix[1][1]);
	
	if(param->generator==1) dum=atan(param->stop_mix[2][1]/param->stop_mix[1][1]);
	
	param->stop_mix[1][1]=cos(dum);
	param->stop_mix[2][1]=-sin(dum);
	param->stop_mix[1][2]=sin(dum);
	param->stop_mix[2][2]=cos(dum);
		
	if(param->sbot_mix[1][1]==0.)
	{
		mass[1]=param->mass_dnl;
		mass[2]=param->mass_stl;
		mass[3]=param->mass_b1;
		mass[4]=param->mass_dnr;
		mass[5]=param->mass_str;
		mass[6]=param->mass_b2;

		for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) invmat[ie][je]=param->sD_mix[je][ie];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[1][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[1][ie]);
		}
		param->mass_dnl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[2][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[2][ie]);
		}
		param->mass_stl=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[3][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[3][ie]);
		}
		param->mass_b1=mass[iemax];		

		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[4][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[4][ie]);
		}
		param->mass_dnr=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[5][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[5][ie]);
		}
		param->mass_str=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[6][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[6][ie]);
		}
		param->mass_b2=mass[iemax];

		param->sbot_mix[1][1]=param->sD_mix[1][3];
		param->sbot_mix[1][2]=param->sD_mix[1][6];
	}
		
	dum=atan(param->sbot_mix[1][2]/param->sbot_mix[1][1]);

	if(param->generator==1) dum=atan(param->sbot_mix[2][1]/param->sbot_mix[1][1]);

	param->sbot_mix[1][1]=cos(dum);
	param->sbot_mix[2][1]=-sin(dum);
	param->sbot_mix[1][2]=sin(dum);
	param->sbot_mix[2][2]=cos(dum);
		
	if(param->stau_mix[1][1]==0.)
	{
		mass[1]=param->mass_el;
		mass[2]=param->mass_mul;
		mass[3]=param->mass_tau1;
		mass[4]=param->mass_er;
		mass[5]=param->mass_mur;
		mass[6]=param->mass_tau2;

		for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) invmat[ie][je]=param->sE_mix[je][ie];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[1][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[1][ie]);
		}
		param->mass_el=mass[iemax];

		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[2][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[2][ie]);
		}
		param->mass_mul=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[3][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[3][ie]);
		}
		param->mass_tau2=mass[iemax];		
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[4][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[4][ie]);
		}
		param->mass_er=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[5][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[5][ie]);
		}
		param->mass_mur=mass[iemax];
		
		dum=0.;
		for(ie=1;ie<=6;ie++) if(fabs(invmat[6][ie])>dum) 
		{
			iemax=ie;
			dum=fabs(invmat[6][ie]);
		}
		param->mass_tau1=mass[iemax];

		param->stau_mix[1][1]=param->sE_mix[1][3];
		param->stau_mix[1][2]=param->sE_mix[1][6];
	}
		
	dum=atan(param->stau_mix[1][2]/param->stau_mix[1][1]);

	if(param->generator==1) dum=atan(param->stau_mix[2][1]/param->stau_mix[1][1]);

	param->stau_mix[1][1]=cos(dum);
	param->stau_mix[2][1]=-sin(dum);
	param->stau_mix[1][2]=sin(dum);
	param->stau_mix[2][2]=cos(dum);

	if(param->neut_mix[1][1]>0.) for(ie=1;ie<=5;ie++) for(je=1;je<=5;je++) param->neut_mix[ie][je]=-param->neut_mix[ie][je];
	
	param->mass_b_pole=mb_pole(param);
 	param->mass_b_1S=mb_1S(param);
	
	param->mtmt=mt_mt(param);
	
	if(param->CKM_lambda*param->CKM_A*param->CKM_rhobar!=0.)
	{
		double s12,s13,s23,c12,c13,c23;
		double complex expid;
	
		s12=param->CKM_lambda;
		s23=param->CKM_A*param->CKM_lambda*param->CKM_lambda;
		s13=cabs(param->CKM_A*pow(param->CKM_lambda,3.)*(param->CKM_rhobar+I*param->CKM_etabar)*sqrt(1.-pow(param->CKM_A*param->CKM_lambda*param->CKM_lambda,2.))/sqrt(1.-param->CKM_lambda*param->CKM_lambda)/(1.-pow(param->CKM_A*param->CKM_lambda*param->CKM_lambda,2.)*(param->CKM_rhobar+I*param->CKM_etabar)));
		expid=(param->CKM_A*pow(param->CKM_lambda,3.)*(param->CKM_rhobar+I*param->CKM_etabar)*sqrt(1.-pow(param->CKM_A*param->CKM_lambda*param->CKM_lambda,2.))/sqrt(1.-param->CKM_lambda*param->CKM_lambda)/(1.-pow(param->CKM_A*param->CKM_lambda*param->CKM_lambda,2.)*(param->CKM_rhobar+I*param->CKM_etabar)))/s13;
	
		c12=sqrt(1.-s12*s12);
		c13=sqrt(1.-s13*s13);
		c23=sqrt(1.-s23*s23);

		param->CKM[1][1]=c12*c13;
		param->CKM[1][2]=s12*c13;
		param->CKM[1][3]=creal(s13/expid);
		param->IMCKM[1][3]=cimag(s13/expid);
		param->CKM[2][1]=creal(-s12*c23-c12*s23*s13*expid);
		param->IMCKM[2][1]=cimag(-s12*c23-c12*s23*s13*expid);
		param->CKM[2][2]=creal(c12*c23-s12*s23*s13*expid);
		param->IMCKM[2][2]=cimag(c12*c23-s12*s23*s13*expid);
		param->CKM[2][3]=s23*c13;
		param->CKM[3][1]=creal(s12*s23-c12*c23*s13*expid);
		param->IMCKM[3][1]=cimag(s12*s23-c12*c23*s13*expid);
		param->CKM[3][2]=creal(-c12*s23-s12*c23*s13*expid);
		param->IMCKM[3][2]=cimag(-c12*s23-s12*c23*s13*expid);
		param->CKM[3][3]=c23*c13;
	}

	param->Vud=param->CKM[1][1]+I*param->IMCKM[1][1];
	param->Vus=param->CKM[1][2]+I*param->IMCKM[1][2];
	param->Vub=param->CKM[1][3]+I*param->IMCKM[1][3];
	param->Vcd=param->CKM[2][1]+I*param->IMCKM[2][1];
	param->Vcs=param->CKM[2][2]+I*param->IMCKM[2][2];
	param->Vcb=param->CKM[2][3]+I*param->IMCKM[2][3];
	param->Vtd=param->CKM[3][1]+I*param->IMCKM[3][1];
	param->Vts=param->CKM[3][2]+I*param->IMCKM[3][2];
	param->Vtb=param->CKM[3][3]+I*param->IMCKM[3][3];
	
	return;
}

/*--------------------------------------------------------------------*/

int test_slha(char name[])
/* "container" function scanning the SLHA file "name" and checking if it is valid. A negative value indicates a problem with the SLHA file. */
{
	struct parameters param;
		
	Init_param(&param);
	if(Les_Houches_Reader(name,&param))
	{
		if(param.FV!=0) return 2;
		if(param.NMSSM!=0) return 3;
		if(param.THDM_model!=0) return 10;
		return 1;
	}

	return param.model;
}

