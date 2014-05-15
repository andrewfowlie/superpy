* const.h
* some constants
* this file is part of FeynHiggs
* last modified 30 May 12 th


	RealType pi, degree, sqrt2, hbar_c2, GeV_cm, bogus
	parameter (pi = 3.1415926535897932384626433832795029D0)
	parameter (degree = pi/180D0)
	parameter (sqrt2 = 1.41421356237309504880168872421D0)

	parameter (hbar_c2 = 3.8937966D8)
*	  = hbar c^2 in picobarn

	parameter (GeV_cm = 1.98D-14)
*	  = GeV^-1 in cm

	parameter (bogus = -1D123)
*	  some weird number likely to noticeably distort the final result;
*	  used for initializing arrays to check that all components
*	  have been calculated

	ComplexType cI
	parameter (cI = (0D0, 1D0))

	RealType Qe, Qu, Qd
	parameter (Qe = -1, Qu = 2/3D0, Qd = -1/3D0)

	RealType Alfa0, DeltaAlfaLept, DeltaAlfaHad5, DeltaAlfa
	parameter (Alfa0 = 1/137.0359895D0)
	parameter (DeltaAlfaLept = .031497687D0)
	parameter (DeltaAlfaHad5 = .027547D0)
	parameter (DeltaAlfa = DeltaAlfaLept + DeltaAlfaHad5)

	RealType C_F, C_A, T_F
	parameter (C_F = 4/3D0, C_A = 3, T_F = 1/2D0)

	RealType invAlfaMZ_default, GF_default, AlfasMZ_default
	RealType MZ_default, MW_default
	RealType ME_default, MM_default, ML_default
	RealType MU_default, MC_default
	RealType MD_default, MS_default, MB_default
	RealType CKMlambda_default, CKMA_default
	RealType CKMrhobar_default, CKMetabar_default

	parameter (invAlfaMZ_default = (1 - DeltaAlfa)/Alfa0)
	parameter (GF_default = 1.16639D-5)
	parameter (AlfasMZ_default = .118D0)

	parameter (MZ_default = 91.1875D0)
	parameter (MW_default = 80.392D0)

	parameter (ME_default = .510998902D-3)
	parameter (MM_default = 105.658357D-3)
	parameter (ML_default = 1777.03D-3)
	parameter (MU_default = 3D-3)
	parameter (MC_default = 1.286D0)
	parameter (MD_default = 6D-3)
	parameter (MS_default = 95D-3)
	parameter (MB_default = 4.2D0)

	parameter (CKMlambda_default = .2253D0)
	parameter (CKMA_default = .808D0)
	parameter (CKMrhobar_default = .132D0)
	parameter (CKMetabar_default = .341D0)

	RealType GammaW, GammaZ
	parameter (GammaW = 2.118D0)
	parameter (GammaZ = 2.4952D0)

	RealType uncomputable
	parameter (uncomputable = -2)

	integer valid
	parameter (valid = 4711)
