* SLHAWrite.F
* write a file in SUSY Les Houches Accord format (hep-ph/0311123, 0801.0045)
* last modified 28 Jun 12 th

#include "PDG.h"
#include "SLHAIntern.h"


	subroutine SLHAWrite(error, sd, filename)
	implicit none
	integer error
	RealType sd(*)
	character*(*) filename

#if U77EXT
#define len_trim lnblnk
	integer lnblnk
	external lnblnk
#endif

	integer SLHAExist
	external SLHAExist

	integer unit, im, i, j, k, n, decay, next, id
	integer namelhs, namelen, idslen
	character*(SPInfo_Len) text
	character*(11*PDGLen) name
	character*(10*PDGDigits) ids

#include "SLHAWriteBlocks.h"

	if( filename .eq. "-" ) then
	  unit = 6
	else
	  open(10, file=filename, status="unknown", err=997)
	  unit = 10
	endif

	error = 0
	im = 1

100	format("BLOCK ",A,:," Q=",E16.8)
101	format(I6,3x,A)
201	format(I10,3x,1P,E16.8,0P,"   # ",A)
202	format(I10,3x,I16,"   # ",A)
203	format(I10,3x,1P,E16.8,0P,"   # ",A,"(",I1,")")
204	format(I6,I4,3x,1P,E16.8,0P,"   # ",A,"(",I1,",",I1,")")
205	format(I6,I4,I4,3x,1P,E16.8,0P,"   # ",A,"(",I1,",",I1,",",I1,")")
206	format(13x,1P,E16.8,0P,"   # ",A)
300	format("DECAY",I10,3x,1P,E16.8,0P,"   # Gamma(",A,")")
301	format(I PDGDigits)
302	format(3x,1P,E16.8,0P,I4,A,"   # BR(",A,")")

#define Block(i) \
write(unit,100) blockname(i)(5-2*im:len_trim(blockname(i)))

#define Text(n,s) \
write(unit,101) int(n), s(1:len_trim(s))

#define Var(i,val,name) \
if( val .ne. invalid ) write(unit,201) i, val, name

#define IVar(i,val,name) \
if( val .ne. invalid ) write(unit,202) i, int(val), name

#define Vector(i,val,name) \
if( val(i) .ne. invalid ) write(unit,203) i, val(i), name, i

#define Matrix(i,j,val,name) \
if( val(i,j) .ne. invalid ) write(unit,204) i,j, val(i,j), name, i,j

#define Triple(i,j,k,val,name) \
if( val(i,j,k) .ne. invalid ) write(unit,205) i,j,k, val(i,j,k), name, i,j,k

#define Number(val,name) \
if( val .ne. invalid ) write(unit,206) val, name

	n = int(SPInfo_NLines)
	if( n .gt. 0 ) then
	  Block(spinf)
	  do i = 1, n
	    call SLHAGetInfo(BlockSPInfo(1), i, text)
	    Text(SPInfo_Code(i), text)
	  enddo
	endif

	n = int(DCInfo_NLines)
	if( n .gt. 0 ) then
	  Block(tcinf)
	  do i = 1, n
	    call SLHAGetInfo(BlockDCInfo(1), i, text)
	    Text(DCInfo_Code(i), text)
	  enddo
	endif

	do im = 1, SLHAExist(BlockModSel(1), LengthModSel)
	  Block(modse)
	  IVar(1, ModSel_Model, "Model")
	  IVar(2, ModSel_GridPts, "GridPts")
	  IVar(3, ModSel_Content, "Content")
	  IVar(4, ModSel_RPV, "RPV")
	  IVar(5, ModSel_CPV, "CPV")
	  IVar(6, ModSel_FV, "FV")
	enddo

	do im = 1, SLHAExist(BlockSMInputs(1), LengthSMInputs)
	  Block(sminp)
	  Var(1, SMInputs_invAlfaMZ, "invAlfaMZ")
	  Var(2, SMInputs_GF, "GF")
	  Var(3, SMInputs_AlfasMZ, "AlfasMZ")
	  Var(4, SMInputs_MZ, "MZ")
	  Var(5, SMInputs_Mb, "Mb")
	  Var(6, SMInputs_Mt, "Mt")
	  Var(7, SMInputs_Mtau, "Mtau")
	  Var(8, SMInputs_Mnu3, "Mnu3")
	  Var(11, SMInputs_Me, "Me")
	  Var(12, SMInputs_Mnu1, "Mnu1")
	  Var(13, SMInputs_Mmu, "Mmu")
	  Var(14, SMInputs_Mnu2, "Mnu2")
	  Var(21, SMInputs_Md, "Md")
	  Var(22, SMInputs_Mu, "Mu")
	  Var(23, SMInputs_Ms, "Ms")
	  Var(24, SMInputs_Mc, "Mc")
	enddo

	do im = 1, SLHAExist(BlockMinPar(1), LengthMinPar)
	  Block(minpa)
	  Var(1, MinPar_M0, "M0")
	  Var(2, MinPar_M12, "M12")
	  Var(3, MinPar_TB, "TB")
	  Var(4, MinPar_signMUE, "signMUE")
	  Var(5, MinPar_A, "A")
	  Var(6, MinPar_cgrav, "cgrav")
	enddo

	do im = 1, SLHAExist(BlockExtPar(1), LengthExtPar)
	  Block(extpa)
	  Var(0, ExtPar_Q, "Q")
	  Var(1, ExtPar_M1, "M1")
	  Var(2, ExtPar_M2, "M2")
	  Var(3, ExtPar_M3, "M3")
	  Var(11, ExtPar_At, "At")
	  Var(12, ExtPar_Ab, "Ab")
	  Var(13, ExtPar_Atau, "Atau")
	  Var(21, ExtPar_MHd2, "MHd2")
	  Var(22, ExtPar_MHu2, "MHu2")
	  Var(23, ExtPar_MUE, "MUE")
	  Var(24, ExtPar_MA02, "MA02")
	  Var(25, ExtPar_TB, "TB")
	  Var(26, ExtPar_MA0, "MA0")
	  Var(27, ExtPar_MHp, "MHp")
	  Var(31, ExtPar_MSL(1), "MSL(1)")
	  Var(32, ExtPar_MSL(2), "MSL(2)")
	  Var(33, ExtPar_MSL(3), "MSL(3)")
	  Var(34, ExtPar_MSE(1), "MSE(1)")
	  Var(35, ExtPar_MSE(2), "MSE(2)")
	  Var(36, ExtPar_MSE(3), "MSE(3)")
	  Var(41, ExtPar_MSQ(1), "MSQ(1)")
	  Var(42, ExtPar_MSQ(2), "MSQ(2)")
	  Var(43, ExtPar_MSQ(3), "MSQ(3)")
	  Var(44, ExtPar_MSU(1), "MSU(1)")
	  Var(45, ExtPar_MSU(2), "MSU(2)")
	  Var(46, ExtPar_MSU(3), "MSU(3)")
	  Var(47, ExtPar_MSD(1), "MSD(1)")
	  Var(48, ExtPar_MSD(2), "MSD(2)")
	  Var(49, ExtPar_MSD(3), "MSD(3)")
	  Var(51, ExtPar_N5(1), "N5(1)")
	  Var(52, ExtPar_N5(2), "N5(2)")
	  Var(53, ExtPar_N5(3), "N5(3)")
	  Var(61, ExtPar_lambda, "lambda")
	  Var(62, ExtPar_kappa, "kappa")
	  Var(63, ExtPar_Alambda, "Alambda")
	  Var(64, ExtPar_Akappa, "Akappa")
	  Var(65, ExtPar_lambdaS, "lambdaS")
	  Var(66, ExtPar_xiF, "xiF")
	  Var(67, ExtPar_xiS, "xiS")
	  Var(68, ExtPar_MUEprime, "MUEprime")
	  Var(69, ExtPar_mS2prime, "mS2prime")
	  Var(70, ExtPar_mS2, "mS2")
	enddo

	do im = 1, SLHAExist(BlockQExtPar(1), LengthQExtPar)
	  Block(qxtpa)
	  Var( 1, QExtPar_QM1, "QM1")
	  Var( 2, QExtPar_QM2, "QM2")
	  Var( 3, QExtPar_QM3, "QM3")
	  Var(11, QExtPar_QAt, "QAt")
	  Var(12, QExtPar_QAb, "QAb")
	  Var(13, QExtPar_QAtau, "QAtau")
	  Var(21, QExtPar_QMHd2, "QMHd2")
	  Var(22, QExtPar_QMHu2, "QMHu2")
	  Var(23, QExtPar_QMUE, "QMUE")
	  Var(24, QExtPar_QMA02, "QMA02")
	  Var(25, QExtPar_QTB, "QTB")
	  Var(31, QExtPar_QMSL, "QMSL")
	  Var(34, QExtPar_QMSE, "QMSE")
	  Var(41, QExtPar_QMSQ, "QMSQ")
	  Var(44, QExtPar_QMSU, "QMSU")
	  Var(47, QExtPar_QMSD, "QMSD")
	enddo

	do im = 1, SLHAExist(BlockNMSSMRun(1), LengthNMSSMRun)
	  Block(nmrun)
	  Var(1, NMSSMRun_lambda, "lambda")
	  Var(2, NMSSMRun_kappa, "kappa")
	  Var(3, NMSSMRun_Alambda, "Alambda")
	  Var(4, NMSSMRun_Akappa, "Akappa")
	  Var(5, NMSSMRun_lambdaS, "lambdaS")
	  Var(6, NMSSMRun_xiF, "xiF")
	  Var(7, NMSSMRun_xiS, "xiS")
	  Var(8, NMSSMRun_MUEprime, "MUEprime")
	  Var(9, NMSSMRun_mS2prime, "mS2prime")
	  Var(10, NMSSMRun_mS2, "mS2")
	enddo

	do im = 1, SLHAExist(BlockMass(1), LengthMass)
	  Block(mass)
	  Var(PDG_nu_e, Mass_Mf(1,1), "Mf(1,1)")
	  Var(PDG_electron, Mass_Mf(2,1), "Mf(2,1)")
	  Var(PDG_up, Mass_Mf(3,1), "Mf(3,1)")
	  Var(PDG_down, Mass_Mf(4,1), "Mf(4,1)")
	  Var(PDG_nu_mu, Mass_Mf(1,2), "Mf(1,2)")
	  Var(PDG_muon, Mass_Mf(2,2), "Mf(2,2)")
	  Var(PDG_charm, Mass_Mf(3,2), "Mf(3,2)")
	  Var(PDG_strange, Mass_Mf(4,2), "Mf(4,2)")
	  Var(PDG_nu_tau, Mass_Mf(1,3), "Mf(1,3)")
	  Var(PDG_tau, Mass_Mf(2,3), "Mf(2,3)")
	  Var(PDG_top, Mass_Mf(3,3), "Mf(3,3)")
	  Var(PDG_bottom, Mass_Mf(4,3), "Mf(4,3)")
	  Var(PDG_snu_e1, Mass_MSf(1,1,1), "MSf(1,1,1)")
	  Var(PDG_snu_e2, Mass_MSf(2,1,1), "MSf(2,1,1)")
	  Var(PDG_selectron1, Mass_MSf(1,2,1), "MSf(1,2,1)")
	  Var(PDG_selectron2, Mass_MSf(2,2,1), "MSf(2,2,1)")
	  Var(PDG_sup1, Mass_MSf(1,3,1), "MSf(1,3,1)")
	  Var(PDG_sup2, Mass_MSf(2,3,1), "MSf(2,3,1)")
	  Var(PDG_sdown1, Mass_MSf(1,4,1), "MSf(1,4,1)")
	  Var(PDG_sdown2, Mass_MSf(2,4,1), "MSf(2,4,1)")
	  Var(PDG_snu_mu1, Mass_MSf(1,1,2), "MSf(1,1,2)")
	  Var(PDG_snu_mu2, Mass_MSf(2,1,2), "MSf(2,1,2)")
	  Var(PDG_smuon1, Mass_MSf(1,2,2), "MSf(1,2,2)")
	  Var(PDG_smuon2, Mass_MSf(2,2,2), "MSf(2,2,2)")
	  Var(PDG_scharm1, Mass_MSf(1,3,2), "MSf(1,3,2)")
	  Var(PDG_scharm2, Mass_MSf(2,3,2), "MSf(2,3,2)")
	  Var(PDG_sstrange1, Mass_MSf(1,4,2), "MSf(1,4,2)")
	  Var(PDG_sstrange2, Mass_MSf(2,4,2), "MSf(2,4,2)")
	  Var(PDG_snu_tau1, Mass_MSf(1,1,3), "MSf(1,1,3)")
	  Var(PDG_snu_tau2, Mass_MSf(2,1,3), "MSf(2,1,3)")
	  Var(PDG_stau1, Mass_MSf(1,2,3), "MSf(1,2,3)")
	  Var(PDG_stau2, Mass_MSf(2,2,3), "MSf(2,2,3)")
	  Var(PDG_stop1, Mass_MSf(1,3,3), "MSf(1,3,3)")
	  Var(PDG_stop2, Mass_MSf(2,3,3), "MSf(2,3,3)")
	  Var(PDG_sbottom1, Mass_MSf(1,4,3), "MSf(1,4,3)")
	  Var(PDG_sbottom2, Mass_MSf(2,4,3), "MSf(2,4,3)")
	  Var(PDG_Z, Mass_MZ, "MZ")
	  Var(PDG_W, Mass_MW, "MW")
	  Var(PDG_h0, Mass_Mh0, "Mh0")
	  Var(PDG_HH, Mass_MHH, "MHH")
	  Var(PDG_A0, Mass_MA0, "MA0")
	  Var(PDG_Hp, Mass_MHp, "MHp")
	  Var(PDG_neutralino1, Mass_MNeu(1), "MNeu(1)")
	  Var(PDG_neutralino2, Mass_MNeu(2), "MNeu(2)")
	  Var(PDG_neutralino3, Mass_MNeu(3), "MNeu(3)")
	  Var(PDG_neutralino4, Mass_MNeu(4), "MNeu(4)")
	  Var(PDG_neutralino5, Mass_MNeu(5), "MNeu(5)")
	  Var(PDG_chargino1, Mass_MCha(1), "MCha(1)")
	  Var(PDG_chargino2, Mass_MCha(2), "MCha(2)")
	  Var(PDG_gluino, Mass_MGl, "MGl")
	  Var(PDG_gravitino, Mass_MGrav, "MGrav")
	enddo

	do im = 1, SLHAExist(BlockDMass(1), LengthDMass)
	  Block(tmass)
	  Var(0, DMass_Q, "Q")
	  Var(PDG_h0, DMass_DeltaMh0, "Delta Mh0")
	  Var(PDG_HH, DMass_DeltaMHH, "Delta MHH")
	  Var(PDG_A0, DMass_DeltaMA0, "Delta MA0")
	  Var(PDG_Hp, DMass_DeltaMHp, "Delta MHp")
	enddo

	do im = 1, SLHAExist(BlockNMix(1), LengthNMix)
	  Block(nmix)
	  do i = 1, 4
	    do j = 1, 4
	      Matrix(i,j, NMix_ZNeu, "ZNeu")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockUMix(1), LengthUMix)
	  Block(umix)
	  do i = 1, 2
	    do j = 1, 2
	      Matrix(i,j, UMix_UCha, "UCha")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockVMix(1), LengthVMix)
	  Block(vmix)
	  do i = 1, 2
	    do j = 1, 2
	      Matrix(i,j, VMix_VCha, "VCha")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockStauMix(1), LengthStauMix)
	  Block(staum)
	  do i = 1, 2
	    do j = 1, 2
	      Matrix(i,j, StauMix_USf, "USf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockStopMix(1), LengthStopMix)
	  Block(stopm)
	  do i = 1, 2
	    do j = 1, 2
	      Matrix(i,j, StopMix_USf, "USf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockSbotMix(1), LengthSbotMix)
	  Block(sbotm)
	  do i = 1, 2
	    do j = 1, 2
	      Matrix(i,j, SbotMix_USf, "USf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockAlpha(1), LengthAlpha)
	  Block(alfa)
	  Number(Alpha_Alpha, "Alpha")
	enddo

	do im = 1, SLHAExist(BlockDAlpha(1), LengthDAlpha)
	  Block(talfa)
	  Number(DAlpha_DeltaAlpha, "Delta Alpha")
	enddo

	do im = 1, SLHAExist(BlockHMix(1), LengthHMix)
	  Block(hmix), HMix_Q
	  Var(1, HMix_MUE, "MUE")
	  Var(2, HMix_TB, "TB")
	  Var(3, HMix_VEV, "VEV")
	  Var(4, HMix_MA02, "MA02")
	enddo

	do im = 1, SLHAExist(BlockGauge(1), LengthGauge)
	  Block(gauge), Gauge_Q
	  Var(1, Gauge_g1, "g1")
	  Var(2, Gauge_g2, "g2")
	  Var(3, Gauge_g3, "g3")
	enddo

	do im = 1, SLHAExist(BlockMSoft(1), LengthMSoft)
	  Block(msoft), MSoft_Q
	  Var(1, MSoft_M1, "M1")
	  Var(2, MSoft_M2, "M2")
	  Var(3, MSoft_M3, "M3")
	  Var(21, MSoft_MHd2, "MHd2")
	  Var(22, MSoft_MHu2, "MHu2")
	  Var(31, MSoft_MSL(1), "MSL(1)")
	  Var(32, MSoft_MSL(2), "MSL(2)")
	  Var(33, MSoft_MSL(3), "MSL(3)")
	  Var(34, MSoft_MSE(1), "MSE(1)")
	  Var(35, MSoft_MSE(2), "MSE(2)")
	  Var(36, MSoft_MSE(3), "MSE(3)")
	  Var(41, MSoft_MSQ(1), "MSQ(1)")
	  Var(42, MSoft_MSQ(2), "MSQ(2)")
	  Var(43, MSoft_MSQ(3), "MSQ(3)")
	  Var(44, MSoft_MSU(1), "MSU(1)")
	  Var(45, MSoft_MSU(2), "MSU(2)")
	  Var(46, MSoft_MSU(3), "MSU(3)")
	  Var(47, MSoft_MSD(1), "MSD(1)")
	  Var(48, MSoft_MSD(2), "MSD(2)")
	  Var(49, MSoft_MSD(3), "MSD(3)")
	enddo

	do im = 1, SLHAExist(BlockAe(1), LengthAe)
	  Block(ae), Ae_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, Ae_Af, "Af")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockAu(1), LengthAu)
	  Block(au), Au_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, Au_Af, "Af")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockAd(1), LengthAd)
	  Block(ad), Ad_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, Ad_Af, "Af")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockYe(1), LengthYe)
	  Block(ye), Ye_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, Ye_Yf, "Yf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockYu(1), LengthYu)
	  Block(yu), Yu_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, Yu_Yf, "Yf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockYd(1), LengthYd)
	  Block(yd), Yd_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, Yd_Yf, "Yf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVLamLLEIn(1), LengthRVLamLLEIn)
	  Block(lllei)
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVLamLLEIn_lamLLE, "lamLLE")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVLamLQDIn(1), LengthRVLamLQDIn)
	  Block(llqdi)
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVLamLQDIn_lamLQD, "lamLQD")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVLamUDDIn(1), LengthRVLamUDDIn)
	  Block(luddi)
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVLamUDDIn_lamUDD, "lamUDD")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVLamLLE(1), LengthRVLamLLE)
	  Block(llle), RVLamLLE_Q
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVLamLLE_lamLLE, "lamLLE")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVLamLQD(1), LengthRVLamLQD)
	  Block(llqd), RVLamLQD_Q
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVLamLQD_lamLQD, "lamLQD")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVLamUDD(1), LengthRVLamUDD)
	  Block(ludd), RVLamUDD_Q
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVLamUDD_lamUDD, "lamUDD")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVTLLEIn(1), LengthRVTLLEIn)
	  Block(tllei)
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVTLLEIn_TLLE, "TLLE")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVTLQDIn(1), LengthRVTLQDIn)
	  Block(tlqdi)
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVTLQDIn_TLQD, "TLQD")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVTUDDIn(1), LengthRVTUDDIn)
	  Block(tuddi)
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVTUDDIn_TUDD, "TUDD")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVTLLE(1), LengthRVTLLE)
	  Block(tlle), RVTLLE_Q
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVTLLE_TLLE, "TLLE")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVTLQD(1), LengthRVTLQD)
	  Block(tlqd), RVTLQD_Q
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVTLQD_TLQD, "TLQD")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVTUDD(1), LengthRVTUDD)
	  Block(tudd), RVTUDD_Q
	  do i = 1, 3
	    do j = 1, 3
	      do k = 1, 3
	        Triple(i,j,k, RVTUDD_TUDD, "TUDD")
	      enddo
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVKappaIn(1), LengthRVKappaIn)
	  Block(rki)
	  do i = 1, 3
	    Vector(i, RVKappaIn_kappa, "kappa")
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVKappa(1), LengthRVKappa)
	  Block(rk), RVKappa_Q
	  do i = 1, 3
	    Vector(i, RVKappa_kappa, "kappa")
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVDIn(1), LengthRVDIn)
	  Block(rdi)
	  do i = 1, 3
	    Vector(i, RVDIn_D, "D")
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVD(1), LengthRVD)
	  Block(rd), RVD_Q
	  do i = 1, 3
	    Vector(i, RVD_D, "D")
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVSnVEVIn(1), LengthRVSnVEVIn)
	  Block(rvevi)
	  do i = 1, 3
	    Vector(i, RVSnVEVIn_VEV, "VEV")
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVSnVEV(1), LengthRVSnVEV)
	  Block(rvev), RVSnVEV_Q
	  do i = 1, 3
	    Vector(i, RVSnVEV_VEV, "VEV")
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVM2LH1In(1), LengthRVM2LH1In)
	  Block(rmlhi)
	  do i = 1, 3
	    Vector(i, RVM2LH1In_M2LH1, "M2LH1")
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVM2LH1(1), LengthRVM2LH1)
	  Block(rmlh), RVM2LH1_Q
	  do i = 1, 3
	    Vector(i, RVM2LH1_M2LH1, "MLH12")
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVNMix(1), LengthRVNMix)
	  Block(rnmix)
	  do i = 1, 7
	    do j = 1, 7
	      Matrix(i,j, RVNMix_ZNeu, "ZNeu")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVUMix(1), LengthRVUMix)
	  Block(rumix)
	  do i = 1, 5
	    do j = 1, 5
	      Matrix(i,j, RVUMix_UCha, "UCha")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVVMix(1), LengthRVVMix)
	  Block(rvmix)
	  do i = 1, 5
	    do j = 1, 5
	      Matrix(i,j, RVVMix_VCha, "VCha")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVHMix(1), LengthRVHMix)
	  Block(rhmix)
	  do i = 1, 5
	    do j = 1, 5
	      Matrix(i,j, RVHMix_UH, "UH")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVAMix(1), LengthRVAMix)
	  Block(ramix)
	  do i = 1, 4
	    do j = 1, 5
	      Matrix(i,j, RVAMix_UA, "UA")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockRVLMix(1), LengthRVLMix)
	  Block(rlmix)
	  do i = 1, 7
	    do j = 1, 8
	      Matrix(i,j, RVLMix_CLep, "CLep")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockVCKMIn(1), LengthVCKMIn)
	  Block(vckmi)
	  Var(1, VCKMIn_lambda, "lambda")
	  Var(2, VCKMIn_A, "A")
	  Var(3, VCKMIn_rhobar, "rhobar")
	  Var(4, VCKMIn_etabar, "etabar")
	enddo

	do im = 1, SLHAExist(BlockVCKM(1), LengthVCKM)
	  Block(vckm), VCKM_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, VCKM_VCKM, "VCKM")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockUPMNSIn(1), LengthUPMNSIn)
	  Block(umnsi)
	  Var(1, UPMNSIn_theta12, "theta12")
	  Var(2, UPMNSIn_theta23, "theta23")
	  Var(3, UPMNSIn_theta13, "theta13")
	  Var(4, UPMNSIn_delta13, "delta13")
	  Var(5, UPMNSIn_alpha1, "alpha1")
	  Var(6, UPMNSIn_alpha2, "alpha2")
	enddo

	do im = 1, SLHAExist(BlockUPMNS(1), LengthUPMNS)
	  Block(umns), UPMNS_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, UPMNS_UPMNS, "UPMNS")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockMSL2In(1), LengthMSL2In)
	  Block(msl2i)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, MSL2In_MSL2, "MSL2")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockMSE2In(1), LengthMSE2In)
	  Block(mse2i)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, MSE2In_MSE2, "MSE2")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockMSQ2In(1), LengthMSQ2In)
	  Block(msq2i)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, MSQ2In_MSQ2, "MSQ2")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockMSU2In(1), LengthMSU2In)
	  Block(msu2i)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, MSU2In_MSU2, "MSU2")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockMSD2In(1), LengthMSD2In)
	  Block(msd2i)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, MSD2In_MSD2, "MSD2")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockMSL2(1), LengthMSL2)
	  Block(msl2), MSL2_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, MSL2_MSL2, "MSL2")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockMSE2(1), LengthMSE2)
	  Block(mse2), MSE2_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, MSE2_MSE2, "MSE2")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockMSQ2(1), LengthMSQ2)
	  Block(msq2), MSQ2_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, MSQ2_MSQ2, "MSQ2")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockMSU2(1), LengthMSU2)
	  Block(msu2), MSU2_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, MSU2_MSU2, "MSU2")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockMSD2(1), LengthMSD2)
	  Block(msd2), MSD2_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, MSD2_MSD2, "MSD2")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockTeIn(1), LengthTeIn)
	  Block(tei)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, TeIn_Tf, "Tf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockTuIn(1), LengthTuIn)
	  Block(tui)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, TuIn_Tf, "Tf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockTdIn(1), LengthTdIn)
	  Block(tdi)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, TdIn_Tf, "Tf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockTe(1), LengthTe)
	  Block(te), Te_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, Te_Tf, "Tf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockTu(1), LengthTu)
	  Block(tu), Tu_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, Tu_Tf, "Tf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockTd(1), LengthTd)
	  Block(td), Td_Q
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, Td_Tf, "Tf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockSnuMix(1), LengthSnuMix)
	  Block(snmix)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, SnuMix_UASf, "UASf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockSelMix(1), LengthSelMix)
	  Block(slmix)
	  do i = 1, 6
	    do j = 1, 6
	      Matrix(i,j, SelMix_UASf, "UASf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockUSqMix(1), LengthUSqMix)
	  Block(usmix)
	  do i = 1, 6
	    do j = 1, 6
	      Matrix(i,j, USqMix_UASf, "UASf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockDSqMix(1), LengthDSqMix)
	  Block(tsmix)
	  do i = 1, 6
	    do j = 1, 6
	      Matrix(i,j, DSqMix_UASf, "UASf")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockSnSMix(1), LengthSnSMix)
	  Block(ssmix)
	  do i = 1, 6
	    do j = 1, 6
	      Matrix(i,j, SnSMix_US, "US")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockSnAMix(1), LengthSnAMix)
	  Block(samix)
	  do i = 1, 6
	    do j = 1, 6
	      Matrix(i,j, SnAMix_UA, "UA")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockCVHMix(1), LengthCVHMix)
	  Block(hcmix)
	  do i = 1, 4
	    do j = 1, 4
	      Matrix(i,j, CVHMix_UH, "UH")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockNMNMix(1), LengthNMNMix)
	  Block(nnmix)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, NMNMix_ZNeu, "ZNeu")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockNMHMix(1), LengthNMHMix)
	  Block(nhmix)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, NMHMix_UH, "UH")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockNMAMix(1), LengthNMAMix)
	  Block(namix)
	  do i = 1, 3
	    do j = 1, 3
	      Matrix(i,j, NMAMix_UA, "UA")
	    enddo
	  enddo
	enddo

	do im = 1, SLHAExist(BlockPrecObs(1), LengthPrecObs)
	  Block(prcob)
	  Var(1, PrecObs_DeltaRho, "DeltaRho")
	  Var(2, PrecObs_MWMSSM, "MWMSSM")
	  Var(3, PrecObs_MWSM, "MWSM")
	  Var(4, PrecObs_SW2effMSSM, "SW2effMSSM")
	  Var(5, PrecObs_SW2effSM, "SW2effSM")
	  Var(11, PrecObs_gminus2mu, "gminus2mu")
	  Var(21, PrecObs_EDMeTh, "EDMeTh")
	  Var(22, PrecObs_EDMn, "EDMn")
	  Var(23, PrecObs_EDMHg, "EDMHg")
	  Var(31, PrecObs_bsgammaMSSM, "bsgammaMSSM")
	  Var(32, PrecObs_bsgammaSM, "bsgammaSM")
	  Var(33, PrecObs_DeltaMsMSSM, "DeltaMsMSSM")
	  Var(34, PrecObs_DeltaMsSM, "DeltaMsSM")
	  Var(35, PrecObs_BsmumuMSSM, "BsmumuMSSM")
	  Var(36, PrecObs_BsmumuSM, "BsmumuSM")
	enddo

	decay = 1
	do while( Decay(decay) .ne. invalid )
	  id = Decay_Id(decay)
	  call SLHAPDGName(id, name)
	  namelhs = len_trim(name(1:PDGLen))
	  write(unit,300) id, Decay(decay + 1), name(1:namelhs)
	  next = decay + Decay_Next(decay)
	  decay = decay + 2

	  name(namelhs+1:namelhs+4) = " ->"
	  namelhs = namelhs + 4

	  do while( decay .lt. next )
	    namelen = namelhs
	    idslen = 1
	    n = Decay_Next(decay) - 1
	    do i = 1, n
	      id = Decay_Id(decay)

	      write(ids(idslen:), 301) id
	      idslen = idslen + len_trim(ids(idslen:idslen+PDGDigits))

	      name(namelen:namelen) = ' '
	      call SLHAPDGName(id, name(namelen+1:))
	      namelen = namelen + len_trim(name(namelen:namelen+PDGLen))

	      decay = decay + 1
	    enddo
	    write(unit,302) Decay(decay), n,
     &        ids(1:idslen-1), name(1:namelen-1)
	    decay = decay + 1
	  enddo
	enddo

	close(unit)
	return

997	error = 1
	end

