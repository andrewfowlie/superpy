! This file is part of HiggsBounds
!  -KW
!******************************************************************
module theory_colliderSfunctions
!******************************************************************
!currently, all tevS are all valid in range [50:400]
!if these changes, change comment at the beginning of theory_XS_SM_functions.f90

 contains

!xsection files from OB (his calculations)
#include "cs-ratios_sigma-bg-Hb/Tevatron.R_H_bb.h"
#include "cs-ratios_sigma-bg-Hb/Tevatron.R_H_gg.h"
#include "cs-ratios_sigma-bg-Hb/Tevatron.R_HWp_csb.h"
#include "cs-ratios_sigma-bg-Hb/Tevatron.R_HWm_dub.h"
#include "cs-ratios_sigma-bg-Hb/Tevatron.R_HWm_scb.h"
#include "cs-ratios_sigma-bg-Hb/Tevatron.R_HWp_udb.h"
#include "cs-ratios_sigma-bg-Hb/Tevatron.R_HZ_bbb.h"
#include "cs-ratios_sigma-bg-Hb/Tevatron.R_HZ_ccb.h"
#include "cs-ratios_sigma-bg-Hb/Tevatron.R_HZ_ddb.h"
#include "cs-ratios_sigma-bg-Hb/Tevatron.R_HZ_ssb.h"
#include "cs-ratios_sigma-bg-Hb/Tevatron.R_HZ_uub.h"

#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_H_bb.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_H_gg.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HVBF_WW.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HVBF_ZZ.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HWp_csb.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HWm_dub.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HWm_scb.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HWp_udb.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HZ_gg.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HZ_bbb.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HZ_ccb.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HZ_ddb.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HZ_ssb.h"
#include "cs-ratios_sigma-bg-Hb/LHC7TeV.R_HZ_uub.h"

! VBF ratios for 8 TeV not yet included - 7 TeV ratios used in the code !
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_H_bb.h"
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_H_gg.h"
!#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HVBF_WW.h"
!#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HVBF_ZZ.h"
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HWp_csb.h"
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HWm_dub.h"
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HWm_scb.h"
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HWp_udb.h"
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HZ_gg.h"
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HZ_bbb.h"
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HZ_ccb.h"
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HZ_ddb.h"
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HZ_ssb.h"
#include "cs-ratios_sigma-bg-Hb/LHC8TeV.R_HZ_uub.h"


end module theory_colliderSfunctions
!******************************************************************
