#include "include.h"

void flha_generator(char name[], char name_output[])
{
	struct parameters param;
	FILE *slha,*output;
	char dummy[100];
	int test;
	
	output=fopen(name_output,"w");
	
	Init_param(&param);
	Les_Houches_Reader(name,&param);

	fprintf(output,"# SuperIso output in Flavour Les Houches Accord format\n");
	fprintf(output,"Block FCINFO  # Program information\n");
	fprintf(output,"     1     SUPERISO         # flavour calculator\n");
	fprintf(output,"     2     3.3              # version number\n");
	
	fprintf(output,"Block MODSEL  # Model selection\n");

	if(test_file(name))
	{
		slha=fopen(name,"r");
		test=1;
		while((EOF != fscanf(slha,"%s",dummy))&&(test))
		{
			if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(slha,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
			if(!strcasecmp(dummy,"MODSEL"))
			{
				while((EOF != fscanf(slha,"%s",dummy))&&(strcasecmp(dummy,"Block")&&strcasecmp(dummy,"Decay"))) 
				{
					if(!strncasecmp("#",dummy,1)) while ((EOF!=fscanf(slha,"%c",dummy))&&(strncasecmp("\n",dummy,1)));
					if((atoi(dummy)*((atoi(dummy)-atof(dummy))==0.)))
					{
						fprintf(output,"     %s",dummy);
						while ((EOF!=fscanf(slha,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fprintf(output,"%c",dummy[0]);
						fprintf(output,"\n");
						test=0;
					}
				}
				if(!strcasecmp(dummy,"Decay")) while((!fseek(slha,-1,SEEK_CUR))&&(EOF != fscanf(slha,"%c",dummy))&&(strncasecmp("\n",dummy,1))) fseek(slha,-1,SEEK_CUR);
			}
		}
		fclose(slha);
	}
	else fprintf(output,"     1     0     # Unknown\n");
	if(test) fprintf(output,"     1     0     # Unknown\n");

	fprintf(output,"Block SMINPUTS  # Standard Model inputs\n");
	fprintf(output,"     1     %.8e   # alpha_em^(-1)\n",param.inv_alpha_em);
	fprintf(output,"     2     %.8e   # G_Fermi\n",param.Gfermi);
	fprintf(output,"     3     %.8e   # alpha_s(M_Z)\n",param.alphas_MZ);
	fprintf(output,"     4     %.8e   # m_{Z}(pole)\n",param.mass_Z);
	fprintf(output,"     5     %.8e   # m_{b}(m_{b})\n",param.mass_b);
	fprintf(output,"     6     %.8e   # m_{top}(pole)\n",param.mass_top_pole);
	fprintf(output,"     7     %.8e   # m_{tau}(pole)\n",param.mass_tau_pole);


	fprintf(output,"Block FMASS  # Mass spectrum in GeV\n");
	fprintf(output,"#PDG_code  mass        scheme scale particle\n");
	fprintf(output,"     5     %.8e   3   0   # b (1S)\n",param.mass_b_1S);
	fprintf(output,"   211     %.8e   0   0   # pi+\n",param.m_pi);
	fprintf(output,"   313     %.8e   0   0   # K*\n",param.m_Kstar);
	fprintf(output,"   321     %.8e   0   0   # K+\n",param.m_K);
	fprintf(output,"   421     %.8e   0   0   # D0\n",param.m_D);
	fprintf(output,"   431     %.8e   0   0   # D_s+\n",param.m_Ds);
	fprintf(output,"   521     %.8e   0   0   # B+\n",param.m_B);
	fprintf(output,"   531     %.8e   0   0   # B_s\n",param.m_Bs);
	

	fprintf(output,"Block FLIFE  # Lifetime in sec\n");
	fprintf(output,"#PDG_code  lifetime         particle\n");
	fprintf(output,"   211     %.8e   # pi+\n",param.life_pi);
	fprintf(output,"   321     %.8e   # K+\n",param.life_K);
	fprintf(output,"   431     %.8e   # D_s+\n",param.life_Ds);
	fprintf(output,"   521     %.8e   # B+\n",param.life_B);
	fprintf(output,"   531     %.8e   # B_s\n",param.life_Bs);

	fprintf(output,"Block FCONST  # Decay constant in GeV\n");
	fprintf(output,"#PDG_code number decay_constant scheme scale particle\n");
	fprintf(output,"   431     1   %.8e     0      0   # D_s+\n",param.f_Ds);
	fprintf(output,"   521     1   %.8e     0      0   # B+\n",param.f_B);
	fprintf(output,"   531     1   %.8e     0      0   # B_s\n",param.f_Bs);

	fprintf(output,"Block FCONSTRATIO  # Ratio of decay constants\n");
	fprintf(output,"#PDG_code1 code2  nb1 nb2 ratio            scheme scale comment\n");
	fprintf(output,"   321     211    1   1   %.8e     0      0   # f_K/f_pi\n",param.fK_fpi);

	double C0w[11],C1w[11],C2w[11],C0b1[11],C1b1[11],C2b1[11],C0b2[11],C1b2[11],C0spec[11],C1spec[11],C0b3[11],C1b3[11],Cpb[11];
	double complex CQ0b[3],CQ1b[3],CQpb[3];
	double mu_W,mu_b1,mu_b2,mu_b3,lambda_h,mu_spec;		

	mu_W=2.*param.mass_W;
	
	mu_b1=param.mass_b_1S/2.;
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b1,C1b1,C2b1,mu_b1,&param);
	
	mu_b2=param.mass_b_1S/2.;
	lambda_h=0.5;
	mu_spec=sqrt(lambda_h*param.mass_b);
	C_calculator_base2(C0w,C1w,mu_W,C0b2,C1b2,mu_b2,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0spec,C1spec,mu_spec,&param);
	
	mu_b3=param.mass_b_pole;			
	C_calculator_base2(C0w,C1w,mu_W,C0b3,C1b3,mu_b3,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b3,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b3,&param);

	fprintf(output,"Block FOBS  # Flavour observables\n");
	fprintf(output,"# ParentPDG type  value       q   NDA   ID1   ID2  ID3 ... comment\n");
	fprintf(output,"    5    1   %.8e   0     2     3    22        # BR(b->s gamma)\n",bsgamma(C0b1,C1b1,C2b1,mu_b1,mu_W,&param));
	fprintf(output,"  521    4   %.8e   0     2   313    22        # Delta0(B->K* gamma)\n",delta0(C0b2,C0spec,C1b2,C1spec,&param,mu_b2,mu_spec,lambda_h));
	fprintf(output,"  531    1   %.8e   0     2    13   -13        # BR(B_s->mu+ mu-)\n",Bsmumu(C0b3,C1b3,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b3));
	fprintf(output,"  521    1   %.8e   0     2   -15    16        # BR(B_u->tau nu)\n",Btaunu(&param));
	fprintf(output,"  521    2   %.8e   0     2   -15    16        # R(B_u->tau nu)\n",RBtaunu(&param));
	fprintf(output,"  431    1   %.8e   0     2   -15    16        # BR(D_s->tau nu)\n",Dstaunu(&param));
	fprintf(output,"  431    1   %.8e   0     2   -13    14        # BR(D_s->mu nu)\n",Dsmunu(&param));
	fprintf(output,"  521    1   %.8e   0     3   421   -15    16  # BR(B+->D0 tau nu)\n",BDtaunu(&param));
	fprintf(output,"  521   11   %.8e   0     3   421   -15    16  # BR(B+->D0 tau nu)/BR(B+-> D0 e nu)\n",BDtaunu_BDenu(&param));
	fprintf(output,"  321   11   %.8e   0     2   -13    14        # BR(K->mu nu)/BR(pi->mu nu)\n",Kmunu_pimunu(&param));
	fprintf(output,"  321   12   %.8e   0     2   -13    14        # R_mu23\n",Rmu23(&param));

	param.SM=1;
	
	mu_W=2.*param.mass_W;
	mu_b1=param.mass_b_1S/2.;
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base1(C0w,C1w,C2w,mu_W,C0b1,C1b1,C2b1,mu_b1,&param);
	
	mu_b2=param.mass_b_1S/2.;
	lambda_h=0.5;
	mu_spec=sqrt(lambda_h*param.mass_b);
	CW_calculator(C0w,C1w,C2w,mu_W,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0b2,C1b2,mu_b2,&param);
	C_calculator_base2(C0w,C1w,mu_W,C0spec,C1spec,mu_spec,&param);
	
	mu_b3=param.mass_b_pole;		
	C_calculator_base2(C0w,C1w,mu_W,C0b3,C1b3,mu_b3,&param);
	CQ_calculator(CQ0b,CQ1b,mu_W,mu_b3,&param);
	Cprime_calculator(Cpb,CQpb,mu_W,mu_b3,&param);
	
	fprintf(output,"Block FOBSSM  # SM predictions for flavour observables\n");
	fprintf(output,"# ParentPDG type  value       q   NDA   ID1   ID2  ID3 ... comment\n");
	fprintf(output,"    5    1   %.8e   0     2     3    22        # BR(b->s gamma)\n",bsgamma(C0b1,C1b1,C2b1,mu_b1,mu_W,&param));
	fprintf(output,"  521    4   %.8e   0     2   313    22        # Delta0(B->K* gamma)\n",delta0(C0b2,C0spec,C1b2,C1spec,&param,mu_b2,mu_spec,lambda_h));
	fprintf(output,"  531    1   %.8e   0     2    13   -13        # BR(B_s->mu+ mu-)\n",Bsmumu(C0b3,C1b3,CQ0b,CQ1b,Cpb,CQpb,&param,mu_b3));
	fprintf(output,"  521    1   %.8e   0     2   -15    16        # BR(B_u->tau nu)\n",Btaunu(&param));
	fprintf(output,"  521    2   %.8e   0     2   -15    16        # R(B_u->tau nu)\n",RBtaunu(&param));
	fprintf(output,"  431    1   %.8e   0     2   -15    16        # BR(D_s->tau nu)\n",Dstaunu(&param));
	fprintf(output,"  431    1   %.8e   0     2   -13    14        # BR(D_s->mu nu)\n",Dsmunu(&param));
	fprintf(output,"  521    1   %.8e   0     3   421   -15    16  # BR(B+->D0 tau nu)\n",BDtaunu(&param));
	fprintf(output,"  521   11   %.8e   0     3   421   -15    16  # BR(B+->D0 tau nu)/BR(B+-> D0 e nu)\n",BDtaunu_BDenu(&param));
	fprintf(output,"  321   11   %.8e   0     2   -13    14        # BR(K->mu nu)/BR(pi->mu nu)\n",Kmunu_pimunu(&param));
	fprintf(output,"  321   12   %.8e   0     2   -13    14        # R_mu23\n",Rmu23(&param));

	fclose(output);
	
	return;

}

