module PDGnumbering
 implicit none

 integer,parameter :: tquark = 6
 integer,parameter :: bquark = 5
 integer,parameter :: cquark = 4
 integer,parameter :: squark = 3
 integer,parameter :: uquark = 2
 integer,parameter :: dquark = 1

 integer,parameter :: tbar = -6
 integer,parameter :: bbar = -5
 integer,parameter :: cbar = -4
 integer,parameter :: sbar = -3
 integer,parameter :: ubar = -2
 integer,parameter :: dbar = -1

 integer,parameter :: em    =  11
 integer,parameter :: ep    = -11
 integer,parameter :: mum   =  13
 integer,parameter :: mup   = -13
 integer,parameter :: taum  =  15
 integer,parameter :: taup  = -15

 integer,parameter :: nue   =  12
 integer,parameter :: numu  =  14
 integer,parameter :: nutau =  16

 integer,parameter :: nuebar   =  -12
 integer,parameter :: numubar  =  -14
 integer,parameter :: nutaubar =  -16

 integer,parameter :: h0     =  25
 integer,parameter :: HH     =  35
 integer,parameter :: A0     =  36

 integer,parameter :: h03     =  45
 integer,parameter :: A02     =  46

 integer,parameter :: Hp     =  37
 integer,parameter :: Hm     = -37
 integer,parameter :: Z0     =  23
 integer,parameter :: Wp     =  24
 integer,parameter :: Wm     = -24
 integer,parameter :: photon =  22
 integer,parameter :: gluon  =  21

 integer,parameter :: neut1 = 1000022
 integer,parameter :: neut2 = 1000023
 integer,parameter :: neut3 = 1000025
 integer,parameter :: neut4 = 1000035
 integer,parameter :: neut5 = 1000045

 integer,parameter :: char1p =  1000024
 integer,parameter :: char1m = -1000024
 integer,parameter :: char2p =  1000037
 integer,parameter :: char2m = -1000037
 
 integer,parameter :: s_t1 = 1000006
 integer,parameter :: s_b1 = 1000005
 integer,parameter :: s_cL = 1000004
 integer,parameter :: s_sL = 1000003
 integer,parameter :: s_uL = 1000002
 integer,parameter :: s_dL = 1000001

 integer,parameter :: s_t2 = 2000006
 integer,parameter :: s_b2 = 2000005
 integer,parameter :: s_cR = 2000004
 integer,parameter :: s_sR = 2000003
 integer,parameter :: s_uR = 2000002
 integer,parameter :: s_dR = 2000001

 integer,parameter :: s_emL    =  1000011
 integer,parameter :: s_mumL   =  1000013
 integer,parameter :: s_taumL  =  1000015

 integer,parameter :: s_emR    =  2000011
 integer,parameter :: s_mumR   =  2000013
 integer,parameter :: s_taumR  =  2000015

 integer,parameter :: s_nueL   =  1000012
 integer,parameter :: s_numuL  =  1000014
 integer,parameter :: s_nutauL =  1000016

 integer,parameter :: s_nutau1A =  1000017
 integer,parameter :: s_nutau2A =  1000018
 integer,parameter :: s_nutau3A =  1000019

 integer,parameter :: gluino  = 1000021

end module PDGnumbering
