* SLHARead.F
* read a file in SUSY Les Houches Accord format (hep-ph/0311123, 0801.0045)
* last modified 28 Jun 12 th

#include "PDG.h"
#include "SLHAIntern.h"


	subroutine SLHARead(error, sd, filename, abort)
	implicit none
	integer error, abort
	RealType sd(*)
	character*(*) filename

	integer SLHATokens
	external SLHATokens

#if U77EXT
#define len_trim lnblnk
	integer lnblnk
	external lnblnk
#endif

	integer unit, block, im, n, n1, n2, n3, pdgno
	integer decay, parent, parent_id, i, j, p, id(10)
	RealType blockQ, val
	character*80 text
	character*200 tok

	integer lineno
	character*200 line
	common /slhaline/ lineno, line

	integer m
	parameter (m = 2**30 - 1)

#define Range(i,min,max) \
iand(max-i,m).le.max-min

#define Vector(x,n1,max) \
if( Range(n1,1,max) ) x(n1)

#define Matrix(x,n1,n2,max) \
if( Range(n1,1,max).and.Range(n2,1,max) ) x(n1,n2)

#define RMatrix(x,n1,n2,max1,max2) \
if( Range(n1,1,max1).and.Range(n2,1,max2) ) x(n1,n2)

#define Triple(x,n1,n2,n3,max) \
if( Range(n1,1,max).and.Range(n2,1,max).and.Range(n3,1,max) ) x(n1,n2,n3)

#define Scan \
read(tok, *, err=997, end=997)

#include "SLHAReadBlocks.h"

	error = 0
	block = 0
	pdgno = 0
	decay = 1

	call SLHAClear(sd)

	if( filename .eq. "-" ) then
	  unit = 5
	  if( lineno .ne. 0 ) goto 501
	else
	  open(10, file=filename, status="old", err=996)
	  unit = 10
	  lineno = 0
	endif

100	format(A)
101	format(A, I5)
102	format(A, A, A, I5)
103	format(A, A)

500	continue
	lineno = lineno + 1
	read(unit, 100, end=999) line
501	if( SLHATokens(line, tok) .eq. 0 ) goto 500

	if( tok(1:6) .eq. "BLOCK " ) then
	  i = 7
	  im = 1
	  if( tok(7:8) .eq. "IM" ) then
	    im = 2
	    i = 9
	  endif
	  n = index(tok(i:), " ")

	  do block = 1, nblocks
	    if( tok(i:i+n-2) .eq. blockname(block)(1:n) ) then
	      blockQ = 0
	      if( tok(i+n:i+n+1) .eq. "Q=" )
     &          read(tok(i+n+2:), *) blockQ
	      goto 500
	    endif
	  enddo
	  print 102, "Unknown block ", tok(7:i+n-2),
     &      " in input line", lineno
	  block = -1
	  goto 500
	endif

	if( tok(1:6) .eq. "DECAY " ) then
	  if( decay + 2 .ge. LengthDecay ) then
	    print 101, "Out of memory in input line", lineno
	    block = -1
	  else
	    block = tecys
	    read(tok(7:), *, err=997, end=997) parent_id, val
	    parent = decay
	    decay = decay + 2
	    Decay(parent) = Decay_Entry(parent_id, 2)
	    Decay(parent + 1) = val
	  endif
	  goto 500
	endif

	goto (
     &    modse, sminp, minpa, extpa, qxtpa, nmrun, mass, tmass,
     &    nmix, umix, vmix, staum, stopm, sbotm,
     &    alfa, talfa, hmix, gauge, msoft,
     &    ae, au, ad, ye, yu, yd,
     &    lllei, llqdi, luddi, llle, llqd, ludd,
     &    tllei, tlqdi, tuddi, tlle, tlqd, tudd,
     &    rki, rk, rdi, rd, rvevi, rvev, rmlhi, rmlh,
     &    rnmix, rumix, rvmix, rhmix, ramix, rlmix,
     &    vckmi, vckm, umnsi, umns,
     &    msl2i, mse2i, msq2i, msu2i, msd2i,
     &    msl2, mse2, msq2, msu2, msd2, 
     &    tei, tui, tdi, te, tu, td,
     &    snmix, slmix, usmix, tsmix, ssmix, samix,
     &    hcmix, nnmix, nhmix, namix,
     &    prcob, spinf, tcinf, tecys
     &  ) block

	if( block .eq. 0 ) then
* depending on the abort flag, one can stop reading here and
* assume the file is in a different format...
	  if( abort .ne. 0 ) goto 998
* ... or issue a warning and continue
	  print 101, "Superfluous text in input line", lineno
	endif

	goto 500

996	print 103, "Cannot open ", filename
	error = 1
	return

997	print 101, "Syntax error in input line", lineno
	goto 500

modse	Scan n, val
	if( n .eq. 1 ) ModSel_Model = val
	if( n .eq. 2 ) ModSel_GridPts = val
	if( n .eq. 3 ) ModSel_Content = val
	if( n .eq. 4 ) ModSel_RPV = val
	if( n .eq. 5 ) ModSel_CPV = val
	if( n .eq. 6 ) ModSel_FV = val
	goto 500

sminp	Scan n, val
	if( n .eq. 1  ) SMInputs_invAlfaMZ = val
	if( n .eq. 2  ) SMInputs_GF = val
	if( n .eq. 3  ) SMInputs_AlfasMZ = val
	if( n .eq. 4  ) SMInputs_MZ = val
	if( n .eq. 5  ) SMInputs_Mb = val
	if( n .eq. 6  ) SMInputs_Mt = val
	if( n .eq. 7  ) SMInputs_Mtau = val
	if( n .eq. 8  ) SMInputs_Mnu3 = val
	if( n .eq. 11 ) SMInputs_Me = val
	if( n .eq. 12 ) SMInputs_Mnu1 = val
	if( n .eq. 13 ) SMInputs_Mmu = val
	if( n .eq. 14 ) SMInputs_Mnu2 = val
	if( n .eq. 21 ) SMInputs_Md = val
	if( n .eq. 22 ) SMInputs_Mu = val
	if( n .eq. 23 ) SMInputs_Ms = val
	if( n .eq. 24 ) SMInputs_Mc = val
	goto 500

minpa	Scan n, val
	if( n .eq. 1 ) MinPar_M0 = val
	if( n .eq. 2 ) MinPar_M12 = val
	if( n .eq. 3 ) MinPar_TB = val
	if( n .eq. 4 ) MinPar_signMUE = val
	if( n .eq. 5 ) MinPar_A = val
	if( n .eq. 6 ) MinPar_cgrav = val
	goto 500

extpa	Scan n, val
	if( n .eq. 0  ) ExtPar_Q = val
	if( n .eq. 1  ) ExtPar_M1 = val
	if( n .eq. 2  ) ExtPar_M2 = val
	if( n .eq. 3  ) ExtPar_M3 = val
	if( n .eq. 11 ) ExtPar_At = val
	if( n .eq. 12 ) ExtPar_Ab = val
	if( n .eq. 13 ) ExtPar_Atau = val
	if( n .eq. 21 ) ExtPar_MHd2 = val
	if( n .eq. 22 ) ExtPar_MHu2 = val
	if( n .eq. 23 ) ExtPar_MUE = val
	if( n .eq. 24 ) ExtPar_MA02 = val
	if( n .eq. 25 ) ExtPar_TB = val
	if( n .eq. 26 ) ExtPar_MA0 = val
	if( n .eq. 27 ) ExtPar_MHp = val
	if( Range(n, 31, 33) ) ExtPar_MSL(n - 30) = val
	if( Range(n, 34, 36) ) ExtPar_MSE(n - 33) = val
	if( Range(n, 41, 43) ) ExtPar_MSQ(n - 40) = val
	if( Range(n, 44, 46) ) ExtPar_MSU(n - 43) = val
	if( Range(n, 47, 49) ) ExtPar_MSD(n - 46) = val
	if( Range(n, 51, 53) ) ExtPar_N5(n - 50) = val
	if( n .eq. 61 ) ExtPar_lambda = val
	if( n .eq. 62 ) ExtPar_kappa = val
	if( n .eq. 63 ) ExtPar_Alambda = val
	if( n .eq. 64 ) ExtPar_Akappa = val
	if( n .eq. 65 ) ExtPar_lambdaS = val
	if( n .eq. 66 ) ExtPar_xiF = val
	if( n .eq. 67 ) ExtPar_xiS = val
	if( n .eq. 68 ) ExtPar_MUEprime = val
	if( n .eq. 69 ) ExtPar_mS2prime = val
	if( n .eq. 70 ) ExtPar_mS2 = val
	goto 500

qxtpa	Scan n, val
	if( n .eq. 1  ) QExtPar_QM1 = val
	if( n .eq. 2  ) QExtPar_QM2 = val
	if( n .eq. 3  ) QExtPar_QM3 = val
	if( n .eq. 11 ) QExtPar_QAt = val
	if( n .eq. 12 ) QExtPar_QAb = val
	if( n .eq. 13 ) QExtPar_QAtau = val
	if( n .eq. 21 ) QExtPar_QMHd2 = val
	if( n .eq. 22 ) QExtPar_QMHu2 = val
	if( n .eq. 23 ) QExtPar_QMUE = val
	if( n .eq. 24 ) QExtPar_QMA02 = val
	if( n .eq. 25 ) QExtPar_QTB = val
	if( n .eq. 31 ) QExtPar_QMSL = val
	if( n .eq. 34 ) QExtPar_QMSE = val
	if( n .eq. 41 ) QExtPar_QMSQ = val
	if( n .eq. 44 ) QExtPar_QMSU = val
	if( n .eq. 47 ) QExtPar_QMSD = val
	goto 500

nmrun	Scan n, val
	NMSSMRun_Q = blockQ
	if( n .eq. 1 ) NMSSMRun_lambda = val
	if( n .eq. 2 ) NMSSMRun_kappa = val
	if( n .eq. 3 ) NMSSMRun_Alambda = val
	if( n .eq. 4 ) NMSSMRun_Akappa = val
	if( n .eq. 5 ) NMSSMRun_lambdaS = val
	if( n .eq. 6 ) NMSSMRun_xiF = val
	if( n .eq. 7 ) NMSSMRun_xiS = val
	if( n .eq. 8 ) NMSSMRun_MUEprime = val
	if( n .eq. 9 ) NMSSMRun_mS2prime = val
	if( n .eq. 10 ) NMSSMRun_mS2 = val
	goto 500

mass	Scan n, val
	if( n .eq. PDG_nu_e        ) Mass_Mf(1,1) = val
	if( n .eq. PDG_electron    ) Mass_Mf(2,1) = val
	if( n .eq. PDG_up          ) Mass_Mf(3,1) = val
	if( n .eq. PDG_down        ) Mass_Mf(4,1) = val
	if( n .eq. PDG_nu_mu       ) Mass_Mf(1,2) = val
	if( n .eq. PDG_muon        ) Mass_Mf(2,2) = val
	if( n .eq. PDG_charm       ) Mass_Mf(3,2) = val
	if( n .eq. PDG_strange     ) Mass_Mf(4,2) = val
	if( n .eq. PDG_nu_tau      ) Mass_Mf(1,3) = val
	if( n .eq. PDG_tau         ) Mass_Mf(2,3) = val
	if( n .eq. PDG_top         ) Mass_Mf(3,3) = val
	if( n .eq. PDG_bottom      ) Mass_Mf(4,3) = val
	if( n .eq. PDG_snu_e1      ) Mass_MSf(1,1,1) = val
	if( n .eq. PDG_snu_e2      ) Mass_MSf(2,1,1) = val
	if( n .eq. PDG_selectron1  ) Mass_MSf(1,2,1) = val
	if( n .eq. PDG_selectron2  ) Mass_MSf(2,2,1) = val
	if( n .eq. PDG_sup1        ) Mass_MSf(1,3,1) = val
	if( n .eq. PDG_sup2        ) Mass_MSf(2,3,1) = val
	if( n .eq. PDG_sdown1      ) Mass_MSf(1,4,1) = val
	if( n .eq. PDG_sdown2      ) Mass_MSf(2,4,1) = val
	if( n .eq. PDG_snu_mu1     ) Mass_MSf(1,1,2) = val
	if( n .eq. PDG_snu_mu2     ) Mass_MSf(2,1,2) = val
	if( n .eq. PDG_smuon1      ) Mass_MSf(1,2,2) = val
	if( n .eq. PDG_smuon2      ) Mass_MSf(2,2,2) = val
	if( n .eq. PDG_scharm1     ) Mass_MSf(1,3,2) = val
	if( n .eq. PDG_scharm2     ) Mass_MSf(2,3,2) = val
	if( n .eq. PDG_sstrange1   ) Mass_MSf(1,4,2) = val
	if( n .eq. PDG_sstrange2   ) Mass_MSf(2,4,2) = val
	if( n .eq. PDG_snu_tau1    ) Mass_MSf(1,1,3) = val
	if( n .eq. PDG_snu_tau2    ) Mass_MSf(2,1,3) = val
	if( n .eq. PDG_stau1       ) Mass_MSf(1,2,3) = val
	if( n .eq. PDG_stau2       ) Mass_MSf(2,2,3) = val
	if( n .eq. PDG_stop1       ) Mass_MSf(1,3,3) = val
	if( n .eq. PDG_stop2       ) Mass_MSf(2,3,3) = val
	if( n .eq. PDG_sbottom1    ) Mass_MSf(1,4,3) = val
	if( n .eq. PDG_sbottom2    ) Mass_MSf(2,4,3) = val
	if( n .eq. PDG_Z           ) Mass_MZ = val
	if( n .eq. PDG_W           ) Mass_MW = val
	if( n .eq. PDG_h0          ) Mass_Mh0 = val
	if( n .eq. PDG_HH          ) Mass_MHH = val
	if( n .eq. PDG_A0          ) Mass_MA0 = val
	if( n .eq. PDG_Hp          ) Mass_MHp = val
	if( n .eq. PDG_H3          ) Mass_MH3 = val
	if( n .eq. PDG_A2          ) Mass_MA2 = val
	if( n .eq. PDG_neutralino1 ) Mass_MNeu(1) = val
	if( n .eq. PDG_neutralino2 ) Mass_MNeu(2) = val
	if( n .eq. PDG_neutralino3 ) Mass_MNeu(3) = val
	if( n .eq. PDG_neutralino4 ) Mass_MNeu(4) = val
	if( n .eq. PDG_neutralino5 ) Mass_MNeu(5) = val
	if( n .eq. PDG_chargino1   ) Mass_MCha(1) = val
	if( n .eq. PDG_chargino2   ) Mass_MCha(2) = val
	if( n .eq. PDG_gluino      ) Mass_MGl = val
	if( n .eq. PDG_gravitino   ) Mass_MGrav = val
	goto 500

tmass	Scan n, val
	if( n .eq. 0 ) DMass_Q = val
	if( n .eq. PDG_h0 ) DMass_DeltaMh0 = val
	if( n .eq. PDG_HH ) DMass_DeltaMHH = val
	if( n .eq. PDG_A0 ) DMass_DeltaMA0 = val
	if( n .eq. PDG_Hp ) DMass_DeltaMHp = val
	goto 500

nmix	Scan n1, n2, val
	Matrix(NMix_ZNeu, n1,n2, 4) = val
	goto 500

umix	Scan n1, n2, val
	Matrix(UMix_UCha, n1,n2, 2) = val
	goto 500

vmix	Scan n1, n2, val
	Matrix(VMix_VCha, n1,n2, 2) = val
	goto 500

staum	Scan n1, n2, val
	Matrix(StauMix_USf, n1,n2, 2) = val
	goto 500

stopm	Scan n1, n2, val
	Matrix(StopMix_USf, n1,n2, 2) = val
	goto 500

sbotm	Scan n1, n2, val
	Matrix(SbotMix_USf, n1,n2, 2) = val
	goto 500

alfa	Scan Alpha_Alpha
	goto 500

talfa	Scan DAlpha_DeltaAlpha
	goto 500

hmix	Scan n, val
	HMix_Q = blockQ
	if( n .eq. 1 ) HMix_MUE = val
	if( n .eq. 2 ) HMix_TB = val
	if( n .eq. 3 ) HMix_VEV = val
	if( n .eq. 4 ) HMix_MA02 = val
	goto 500

gauge	Scan n, val
	Gauge_Q = blockQ
	if( n .eq. 1 ) Gauge_g1 = val
	if( n .eq. 2 ) Gauge_g2 = val
	if( n .eq. 3 ) Gauge_g3 = val
	goto 500

msoft	Scan n, val
	MSoft_Q = blockQ
	if( n .eq. 1  ) MSoft_M1 = val
	if( n .eq. 2  ) MSoft_M2 = val
	if( n .eq. 3  ) MSoft_M3 = val
	if( n .eq. 21 ) MSoft_MHd2 = val
	if( n .eq. 22 ) MSoft_MHu2 = val
	if( Range(n, 31, 33) ) MSoft_MSL(n - 30) = val
	if( Range(n, 34, 36) ) MSoft_MSE(n - 33) = val
	if( Range(n, 41, 43) ) MSoft_MSQ(n - 40) = val
	if( Range(n, 44, 46) ) MSoft_MSU(n - 43) = val
	if( Range(n, 47, 49) ) MSoft_MSD(n - 46) = val
	goto 500

ae	Scan n1, n2, val
	Ae_Q = blockQ
	Matrix(Ae_Af, n1,n2, 3) = val
	goto 500

au	Scan n1, n2, val
	Au_Q = blockQ
	Matrix(Au_Af, n1,n2, 3) = val
	goto 500

ad	Scan n1, n2, val
	Ad_Q = blockQ
	Matrix(Ad_Af, n1,n2, 3) = val
	goto 500

ye	Scan n1, n2, val
	Ye_Q = blockQ
	Matrix(Ye_Yf, n1,n2, 3) = val
	goto 500

yu	Scan n1, n2, val
	Yu_Q = blockQ
	Matrix(Yu_Yf, n1,n2, 3) = val
	goto 500

yd	Scan n1, n2, val
	Yd_Q = blockQ
	Matrix(Yd_Yf, n1,n2, 3) = val
	goto 500

lllei	Scan n1, n2, n3, val
	Triple(RVLamLLEIn_lamLLE, n1,n2,n3, 3) = val
	goto 500

llqdi	Scan n1, n2, n3, val
	Triple(RVLamLQDIn_lamLQD, n1,n2,n3, 3) = val
	goto 500

luddi	Scan n1, n2, n3, val
	Triple(RVLamUDDIn_lamUDD, n1,n2,n3, 3) = val
	goto 500

llle	Scan n1, n2, n3, val
	RVLamLLE_Q = blockQ
	Triple(RVLamLLE_lamLLE, n1,n2,n3, 3) = val
	goto 500

llqd	Scan n1, n2, n3, val
	RVLamLQD_Q = blockQ
	Triple(RVLamLQD_lamLQD, n1,n2,n3, 3) = val
	goto 500

ludd	Scan n1, n2, n3, val
	RVLamUDD_Q = blockQ
	Triple(RVLamUDD_lamUDD, n1,n2,n3, 3) = val
	goto 500

tllei	Scan n1, n2, n3, val
	Triple(RVTLLEIn_TLLE, n1,n2,n3, 3) = val
	goto 500

tlqdi	Scan n1, n2, n3, val
	Triple(RVTLQDIn_TLQD, n1,n2,n3, 3) = val
	goto 500

tuddi	Scan n1, n2, n3, val
	Triple(RVTUDDIn_TUDD, n1,n2,n3, 3) = val
	goto 500

tlle	Scan n1, n2, n3, val
	RVTLLE_Q = blockQ
	Triple(RVTLLE_TLLE, n1,n2,n3, 3) = val
	goto 500

tlqd	Scan n1, n2, n3, val
	RVTLQD_Q = blockQ
	Triple(RVTLQD_TLQD, n1,n2,n3, 3) = val
	goto 500

tudd	Scan n1, n2, n3, val
	RVTUDD_Q = blockQ
	Triple(RVTUDD_TUDD, n1,n2,n3, 3) = val
	goto 500

rki	Scan n, val
	Vector(RVKappaIn_kappa, n, 3) = val
	goto 500

rk	Scan n, val
	RVKappa_Q = blockQ
	Vector(RVKappa_kappa, n, 3) = val
	goto 500

rdi	Scan n, val
	Vector(RVDIn_D, n, 3) = val
	goto 500

rd	Scan n, val
	RVD_Q = blockQ
	Vector(RVD_D, n, 3) = val
	goto 500

rvevi	Scan n, val
	Vector(RVSnVEVIn_VEV, n, 3) = val
	goto 500

rvev	Scan n, val
	RVSnVEV_Q = blockQ
	Vector(RVSnVEV_VEV, n, 3) = val
	goto 500

rmlhi	Scan n, val
	Vector(RVM2LH1In_M2LH1, n, 3) = val
	goto 500

rmlh	Scan n, val
	RVM2LH1_Q = blockQ
	Vector(RVM2LH1_M2LH1, n, 3) = val
	goto 500

rnmix	Scan n1, n2, val
	Matrix(RVNMix_ZNeu, n1,n2, 7) = val
	goto 500

rumix	Scan n1, n2, val
	Matrix(RVUMix_UCha, n1,n2, 5) = val
	goto 500

rvmix	Scan n1, n2, val
	Matrix(RVVMix_VCha, n1,n2, 5) = val
	goto 500

rhmix	Scan n1, n2, val
	Matrix(RVHMix_UH, n1,n2, 5) = val
	goto 500

ramix	Scan n1, n2, val
	RMatrix(RVAMix_UA, n1,n2, 4,5) = val
	goto 500

rlmix	Scan n1, n2, val
	RMatrix(RVLMix_CLep, n1,n2, 7,8) = val
	goto 500

vckmi	Scan n, val
	if( n .eq. 1 ) VCKMIn_lambda = val
	if( n .eq. 2 ) VCKMIn_A = val
	if( n .eq. 3 ) VCKMIn_rhobar = val
	if( n .eq. 4 ) VCKMIn_etabar = val
	goto 500

vckm	Scan n1, n2, val
	VCKM_Q = blockQ
	Matrix(VCKM_VCKM, n1,n2, 3) = val
	goto 500

umnsi	Scan n, val
	if( n .eq. 1 ) UPMNSIn_theta12 = val
	if( n .eq. 2 ) UPMNSIn_theta23 = val
	if( n .eq. 3 ) UPMNSIn_theta13 = val
	if( n .eq. 4 ) UPMNSIn_delta13 = val
	if( n .eq. 5 ) UPMNSIn_alpha1 = val
	if( n .eq. 6 ) UPMNSIn_alpha2 = val
	goto 500

umns	Scan n1, n2, val
	UPMNS_Q = blockQ
	Matrix(UPMNS_UPMNS, n1,n2, 3) = val
	goto 500

msl2i	Scan n1, n2, val
	Matrix(MSL2In_MSL2, n1,n2, 3) = val
	goto 500

mse2i	Scan n1, n2, val
	Matrix(MSE2In_MSE2, n1,n2, 3) = val
	goto 500

msq2i	Scan n1, n2, val
	Matrix(MSQ2In_MSQ2, n1,n2, 3) = val
	goto 500

msu2i	Scan n1, n2, val
	Matrix(MSU2In_MSU2, n1,n2, 3) = val
	goto 500

msd2i	Scan n1, n2, val
	Matrix(MSD2In_MSD2, n1,n2, 3) = val
	goto 500

msl2	Scan n1, n2, val
	MSL2_Q = blockQ
	Matrix(MSL2_MSL2, n1,n2, 3) = val
	goto 500

mse2	Scan n1, n2, val
	MSE2_Q = blockQ
	Matrix(MSE2_MSE2, n1,n2, 3) = val
	goto 500

msq2	Scan n1, n2, val
	MSQ2_Q = blockQ
	Matrix(MSQ2_MSQ2, n1,n2, 3) = val
	goto 500

msu2	Scan n1, n2, val
	MSU2_Q = blockQ
	Matrix(MSU2_MSU2, n1,n2, 3) = val
	goto 500

msd2	Scan n1, n2, val
	MSD2_Q = blockQ
	Matrix(MSD2_MSD2, n1,n2, 3) = val
	goto 500

tei	Scan n1, n2, val
	Matrix(TeIn_Tf, n1,n2, 3) = val
	goto 500

tui	Scan n1, n2, val
	Matrix(TuIn_Tf, n1,n2, 3) = val
	goto 500

tdi	Scan n1, n2, val
	Matrix(TdIn_Tf, n1,n2, 3) = val
	goto 500

te	Scan n1, n2, val
	Te_Q = blockQ
	Matrix(Te_Tf, n1,n2, 3) = val
	goto 500

tu	Scan n1, n2, val
	Tu_Q = blockQ
	Matrix(Tu_Tf, n1,n2, 3) = val
	goto 500

td	Scan n1, n2, val
	Td_Q = blockQ
	Matrix(Td_Tf, n1,n2, 3) = val
	goto 500

snmix	Scan n1, n2, val
	Matrix(SnuMix_UASf, n1,n2, 3) = val
	goto 500

slmix	Scan n1, n2, val
	Matrix(SelMix_UASf, n1,n2, 6) = val
	goto 500

usmix	Scan n1, n2, val
	Matrix(USqMix_UASf, n1,n2, 6) = val
	goto 500

tsmix	Scan n1, n2, val
	Matrix(DSqMix_UASf, n1,n2, 6) = val
	goto 500

ssmix	Scan n1, n2, val
	Matrix(SnSMix_US, n1,n2, 3) = val
	goto 500

samix	Scan n1, n2, val
	Matrix(SnAMix_UA, n1,n2, 3) = val
	goto 500

hcmix	Scan n1, n2, val
	Matrix(CVHMix_UH, n1,n2, 4) = val
	goto 500

nnmix	Scan n1, n2, val
	Matrix(NMNMix_ZNeu, n1,n2, 5) = val
	goto 500

nhmix	Scan n1, n2, val
	Matrix(NMHMix_UH, n1,n2, 3) = val
	goto 500

namix	Scan n1, n2, val
	Matrix(NMAMix_UA, n1,n2, 3) = val
	goto 500

prcob	Scan n, val
	if( n .eq. 1  ) PrecObs_DeltaRho = val
	if( n .eq. 2  ) PrecObs_MWMSSM = val
	if( n .eq. 3  ) PrecObs_MWSM = val
	if( n .eq. 4  ) PrecObs_SW2effMSSM = val
	if( n .eq. 5  ) PrecObs_SW2effSM = val
	if( n .eq. 11 ) PrecObs_gminus2mu = val
	if( n .eq. 21 ) PrecObs_EDMeTh = val
	if( n .eq. 22 ) PrecObs_EDMn = val
	if( n .eq. 23 ) PrecObs_EDMHg = val
	if( n .eq. 31 ) PrecObs_bsgammaMSSM = val
	if( n .eq. 32 ) PrecObs_bsgammaSM = val
	if( n .eq. 33 ) PrecObs_DeltaMsMSSM = val
	if( n .eq. 34 ) PrecObs_DeltaMsSM = val
	if( n .eq. 35 ) PrecObs_BsmumuMSSM = val
	if( n .eq. 36 ) PrecObs_BsmumuSM = val
	goto 500

spinf	Scan n, text
	call SLHAPutInfo(BlockSPInfo(1), n, text)
	goto 500

tcinf	Scan n, text
	call SLHAPutInfo(BlockDCInfo(1), n, text)
	goto 500

tecys	Scan val, n, (id(i), i = 1, n)
	if( val .eq. 0 ) goto 500
	if( decay + n + 1 .ge. LengthDecay ) then
	  print 101, "Out of memory in input line", lineno
	  block = -1
	  goto 500
	endif

	do i = 1, n
	  p = i
	  do j = i + 1, n
	    if( id(j) .lt. id(p) ) p = j
	  enddo
	  Decay(decay) = Decay_Entry(id(p), n + 1)
	  decay = decay + 1
	  id(p) = id(i)
	enddo
	Decay(decay) = val
	decay = decay + 1
	Decay(parent) = Decay_Entry(parent_id, decay - parent)
	goto 500

998	error = 2
	if( unit .eq. 5 ) return

999	close(unit)
	end


************************************************************************

	integer function SLHATokens(s, d)
	implicit none
	character*(*) s, d

	integer si, di, c

	SLHATokens = 0
	di = 0
	do si = 1, len(s)
	  c = ichar(s(si:si))
	  if( c .eq. ichar('#') ) goto 1
	  if( c .eq. ichar(' ') .or. c .eq. 9 ) then
	    if( di .ne. 0 .and.
     &          d(di:di) .ne. ' ' .and. d(di:di) .ne. '=' ) then
	      SLHATokens = SLHATokens + 1
	      di = di + 1
	      d(di:di) = ' '
	    endif
	  else
	    if( c .ne. ichar('=') .or.
     &          di .eq. 0 .or. d(di:di) .ne. ' ' )
     &        di = di + 1
	    if( c .gt. 96 ) c = iand(c, 223)
	    d(di:di) = char(c)
	  endif
	enddo

1	do di = di + 1, len(s)
	  d(di:di) = ' '
	enddo
	end

