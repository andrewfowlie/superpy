#include"../../sources/micromegas_aux.h"
#include"../../sources/micromegas.h"
#include"pmodel.h"
#include"pmodel_aux.h"

#define xxxSUGRAc softSusySUGRAc
#define xxxSUGRA  softSusySUGRA
#define xxxsugra_ softsusysugra_ 
#include "sugrac.inc"

#define xxxAMSBc softSusyAMSBc
#define xxxAMSB  softSusyAMSB
#define xxxamsb_ softsusyamsb_ 
#include "amsb.inc"

#define xxxEwsbMSSMc softSusyEwsbMSSMc
#define xxxEwsbMSSM  softSusyEwsbMSSM
#define xxxewsbmssm_ softsusyewsbmssm_ 
#include "pmssm.inc"

#define xxxSUGRAnuhc  softSusySUGRAnuhc
#define xxxSUGRAnuh   softSusySUGRAnuh
#define xxxsugranuhs_  softSusysugranuh_
#include "sugracHiggs.inc"
