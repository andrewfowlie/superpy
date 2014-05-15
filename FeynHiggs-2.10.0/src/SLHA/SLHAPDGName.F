* SLHAPDGName.F
* look up the particle name of the PDG code
* last modified 9 Aug 11 th

#include "PDG.h"


	subroutine SLHAPDGName(code, name)
	implicit none
	integer code
	character*(*) name

	integer npdg
	parameter (npdg = 56)

	integer i, c, pdgcode(npdg)
	character*(PDGLen) pdgname(npdg)

#define Data(i,code,name) data pdgcode(i) /code/, pdgname(i) /name/

	Data( 1, PDG_down,        "down")
	Data( 2, PDG_up,          "up")
	Data( 3, PDG_strange,     "strange")
	Data( 4, PDG_charm,       "charm")
	Data( 5, PDG_bottom,      "bottom")
	Data( 6, PDG_top,         "top")
	Data( 7, PDG_electron,    "electron")
	Data( 8, PDG_nu_e,        "nu_e")
	Data( 9, PDG_muon,        "muon")
	Data(10, PDG_nu_mu,       "nu_mu")
	Data(11, PDG_tau,         "tau")
	Data(12, PDG_nu_tau,      "nu_tau")
	Data(13, PDG_gluon,       "gluon")
	Data(14, PDG_photon,      "photon")
	Data(15, PDG_Z,           "Z")
	Data(16, PDG_W,           "W")
	Data(17, PDG_h0,          "h0")
	Data(18, PDG_HH,          "HH")
	Data(19, PDG_A0,          "A0")
	Data(20, PDG_Hp,          "Hp")
	Data(21, PDG_graviton,    "graviton")
	Data(22, PDG_H3,          "H3")
	Data(23, PDG_A2,          "A2")
	Data(24, PDG_sdown1,      "sdown1")
	Data(25, PDG_sup1,        "sup1")
	Data(26, PDG_sstrange1,   "sstrange1")
	Data(27, PDG_scharm1,     "scharm1")
	Data(28, PDG_sbottom1,    "sbottom1")
	Data(29, PDG_stop1,       "stop1")
	Data(30, PDG_selectron1,  "selectron1")
	Data(31, PDG_snu_e1,      "snu_e1")
	Data(32, PDG_smuon1,      "smuon1")
	Data(33, PDG_snu_mu1,     "snu_mu1")
	Data(34, PDG_stau1,       "stau1")
	Data(35, PDG_snu_tau1,    "snu_tau1")
	Data(36, PDG_gluino,      "gluino")
	Data(37, PDG_neutralino1, "neutralino1")
	Data(38, PDG_neutralino2, "neutralino2")
	Data(39, PDG_chargino1,   "chargino1")
	Data(40, PDG_neutralino3, "neutralino3")
	Data(41, PDG_neutralino4, "neutralino4")
	Data(42, PDG_chargino2,   "chargino2")
	Data(43, PDG_gravitino,   "gravitino")
	Data(44, PDG_neutralino5, "neutralino5")
	Data(45, PDG_sdown2,      "sdown2")
	Data(46, PDG_sup2,        "sup2")
	Data(47, PDG_sstrange2,   "sstrange2")
	Data(48, PDG_scharm2,     "scharm2")
	Data(49, PDG_sbottom2,    "sbottom2")
	Data(50, PDG_stop2,       "stop2")
	Data(51, PDG_selectron2,  "selectron2")
	Data(52, PDG_snu_e2,      "snu_e2")
	Data(53, PDG_smuon2,      "smuon2")
	Data(54, PDG_snu_mu2,     "snu_mu2")
	Data(55, PDG_stau2,       "stau2")
	Data(56, PDG_snu_tau2,    "snu_tau2")

	c = abs(code)
	do i = 1, npdg
	  if( c .eq. pdgcode(i) ) then
	    name = pdgname(i)
	    if( code .lt. 0 )
     &        name(1:1) = char(iand(ichar(name(1:1)), 223))
	    return
	  endif
	enddo

	name = "invalid"
	end
