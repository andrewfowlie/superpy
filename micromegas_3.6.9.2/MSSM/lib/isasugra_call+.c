#include"pmodel.h"
#include"pmodel_aux.h"
#include"../../sources/micromegas_aux.h"
#include"../../sources/micromegas.h"

#define xxxSUGRAc  isajetSUGRAc
#define xxxSUGRA   isajetSUGRA
#define xxxsugra_  isajetsugra_
#include "sugrac.inc"

#define xxxEwsbMSSMc  isajetEwsbMSSMc
#define xxxEwsbMSSM   isajetEwsbMSSM
#define xxxewsbmssm_  isajetewsbmssm_
#include "pmssm.inc"

#define xxxAMSBc  isajetAMSBc
#define xxxAMSB   isajetAMSB
#define xxxamsb_  isajetamsb_
#include "amsb.inc"

#define xxxSUGRAnuhc  isajetSUGRAnuhc
#define xxxSUGRAnuh   isajetSUGRAnuh 
#define xxxsugranuhs_  isajetsugranuh_
#include "sugracHiggs.inc"
