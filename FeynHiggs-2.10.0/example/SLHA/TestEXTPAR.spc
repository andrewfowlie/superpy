# SUSY Les Houches Accord - MSSM spectrum + Decays
# SPheno2.1
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.unizh.ch
Block MODSEL  # Select model
    1    1
Block MINPAR  # Input parameters
    1   0.100000000E+03  # m0      
    2   0.250000000E+03  # m12     
    3   0.100000000E+02  # tanb    
    4   0.100000000E+01  # Sign(mu)
    5  -0.100000000E+03  # A0
#
Block SMINPUTS  # Mass spectrum
          4   0.911870000E+02  # MZ
          5   0.425000000E+01  # mb(mb)
          6   0.175000000E+03  # t
#
Block EXTPAR  # Mass spectrum
#   
   0    0.484786694E+03  # M_input
   1    0.300000000E+03  # M_1
   2    0.600000000E+03  # M_2
   3    0.100000000E+04  # M_3
  11    0.400000000E+03  # At
  12    0.300000000E+03  # Ab
  13    0.200000000E+03  # Atau
  23    0.500000000E+03  # mu
  24    0.250000000E+03  # MA
  31    0.206630723E+03  # ~e_L-
  32    0.206645846E+03  # ~mu_L-
  33    0.134514453E+03  # ~tau_L-
  34    0.143872558E+03  # ~e_R-
  35    0.143838140E+03  # ~mu_R-
  36    0.210401949E+03  # ~tau_R-
  41    0.564892619E+03  # ~u_L
  42    0.564902784E+03  # ~c_L
  43    0.398749215E+03  # ~t_L
  44    0.547790210E+03  # ~u_R
  45    0.547775859E+03  # ~c_R
  46    0.589079372E+03  # ~t_R
  47    0.547601268E+03  # ~d_R
  48    0.547594947E+03  # ~s_R
  49    0.547471349E+03  # ~b_R
    
