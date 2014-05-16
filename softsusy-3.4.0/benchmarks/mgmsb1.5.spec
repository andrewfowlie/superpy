# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-1.61333832e+03
# dtanbeta/dmu=8.03655618e-02
# mu=4.14538157e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=2.80700461e-01
# dtanbeta/dm3sq=-6.72272047e-04
# m3sq= 2.52168530e+04
# SOFTSUSY3.4.0 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.4.0       # version number
Block MODSEL  # Select model
     1    2   # gmsb
Block SMINPUTS             # Standard Model inputs
     1    1.27934000e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637000e-05   # G_Fermi
     3    1.17200000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.25000000e+00   # mb(mb)
     6    1.73300000e+02   # Mtop(pole)
     7    1.77700000e+00   # Mtau(pole)
Block MINPAR               # SUSY breaking input parameters
     3    1.50000000e+01   # tanb, DRbar, Feynman gauge
     4    1.00000000e+00   # sign(mu)
     1    5.50000000e+04   # lambda
     2    1.10000000e+05   # M_mess
     5    3.00000000e+00   # N5
     6    1.00000000e+00   # cgrav
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.28763193e-05
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.04003439e+01   # MW
        25     1.13487701e+02   # h0
        35     5.21240971e+02   # H0
        36     5.20984036e+02   # A0
        37     5.27357284e+02   # H+
   1000021     1.26173534e+03   # ~g
   1000022     2.28536271e+02   # ~neutralino(1)
   1000023     3.75235022e+02   # ~neutralino(2)
   1000024     3.72564527e+02   # ~chargino(1)
   1000025    -4.13392302e+02   # ~neutralino(3)
   1000035     4.99916346e+02   # ~neutralino(4)
   1000037     4.99685539e+02   # ~chargino(2)
   1000039     1.42360821e-09   # ~gravitino
   1000001     1.20490135e+03   # ~d_L
   1000002     1.20241263e+03   # ~u_L
   1000003     1.20489984e+03   # ~s_L
   1000004     1.20241111e+03   # ~c_L
   1000005     1.14763315e+03   # ~b_1
   1000006     1.07335349e+03   # ~t_1
   1000011     3.57185814e+02   # ~e_L
   1000012     3.47962450e+02   # ~nue_L
   1000013     3.57184368e+02   # ~mu_L
   1000014     3.47960967e+02   # ~numu_L
   1000015     1.70022181e+02   # ~stau_1
   1000016     3.47284835e+02   # ~nu_tau_L
   2000001     1.15609473e+03   # ~d_R
   2000002     1.15885507e+03   # ~u_R
   2000003     1.15609262e+03   # ~s_R
   2000004     1.15885400e+03   # ~c_R
   2000005     1.16431471e+03   # ~b_2
   2000006     1.17679656e+03   # ~t_2
   2000011     1.75959408e+02   # ~e_R
   2000013     1.75953454e+02   # ~mu_R
   2000015     3.58201560e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -7.44795660e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.79039382e-01   # N_{1,1}
  1  2    -3.93646181e-02   # N_{1,2}
  1  3     1.70909351e-01   # N_{1,3}
  1  4    -1.03548580e-01   # N_{1,4}
  2  1    -1.83718063e-01   # N_{2,1}
  2  2    -5.63235989e-01   # N_{2,2}
  2  3     5.90960741e-01   # N_{2,3}
  2  4    -5.47520134e-01   # N_{2,4}
  3  1    -4.48702270e-02   # N_{3,1}
  3  2     6.02092041e-02   # N_{3,2}
  3  3     7.00739147e-01   # N_{3,3}
  3  4     7.09454835e-01   # N_{3,4}
  4  1    -7.56057137e-02   # N_{4,1}
  4  2     8.23158855e-01   # N_{4,2}
  4  3     3.61275580e-01   # N_{4,3}
  4  4    -4.31477962e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1    -5.08300256e-01   # U_{1,1}
  1  2     8.61179917e-01   # U_{1,2}
  2  1    -8.61179917e-01   # U_{2,1}
  2  2    -5.08300256e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1    -6.07680457e-01   # V_{1,1}
  1  2     7.94181631e-01   # V_{1,2}
  2  1    -7.94181631e-01   # V_{2,1}
  2  2    -6.07680457e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     2.84784691e-01   # F_{11}
  1  2     9.58591508e-01   # F_{12}
  2  1     9.58591508e-01   # F_{21}
  2  2    -2.84784691e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     4.68302538e-01   # F_{11}
  1  2     8.83568182e-01   # F_{12}
  2  1     8.83568182e-01   # F_{21}
  2  2    -4.68302538e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.09729572e-01   # F_{11}
  1  2     9.93961479e-01   # F_{12}
  2  1     9.93961479e-01   # F_{21}
  2  2    -1.09729572e-01   # F_{22}
Block gauge Q= 1.09054749e+03  # SM gauge couplings
     1     3.63267416e-01   # g'(Q)MSSM DRbar
     2     6.43394094e-01   # g(Q)MSSM DRbar
     3     1.05266439e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.09054749e+03  
  3  3     8.58368627e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.09054749e+03  
  3  3     2.03132806e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.09054749e+03  
  3  3     1.51776003e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.09054749e+03 # Higgs mixing parameters
     1     4.04626369e+02    # mu(Q)MSSM DRbar
     2     1.44867359e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43591260e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     3.09091104e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.09054749e+03  # MSSM DRbar SUSY breaking parameters
     1     2.38808714e+02      # M_1(Q)
     2     4.47862588e+02      # M_2(Q)
     3     1.20463538e+03      # M_3(Q)
    21     1.05466834e+05      # mH1^2(Q)
    22    -1.42379092e+05      # mH2^2(Q)
    31     3.47483338e+02      # meL(Q)
    32     3.47481858e+02      # mmuL(Q)
    33     3.47025974e+02      # mtauL(Q)
    34     1.64108719e+02      # meR(Q)
    35     1.64102351e+02      # mmuR(Q)
    36     1.62130086e+02      # mtauR(Q)
    41     1.16012225e+03      # mqL1(Q)
    42     1.16012070e+03      # mqL2(Q)
    43     1.12329067e+03      # mqL3(Q)
    44     1.11707459e+03      # muR(Q)
    45     1.11707348e+03      # mcR(Q)
    46     1.04218595e+03      # mtR(Q)
    47     1.11313381e+03      # mdR(Q)
    48     1.11313165e+03      # msR(Q)
    49     1.10893786e+03      # mbR(Q)
Block au Q= 1.09054749e+03  
  1  1    -3.61154646e+02      # Au(Q)MSSM DRbar
  2  2    -3.61154139e+02      # Ac(Q)MSSM DRbar
  3  3    -3.40565913e+02      # At(Q)MSSM DRbar
Block ad Q= 1.09054749e+03  
  1  1    -3.83972346e+02      # Ad(Q)MSSM DRbar
  2  2    -3.83971640e+02      # As(Q)MSSM DRbar
  3  3    -3.76216579e+02      # Ab(Q)MSSM DRbar
Block ae Q= 1.09054749e+03  
  1  1    -3.76468028e+01      # Ae(Q)MSSM DRbar
  2  2    -3.76465286e+01      # Amu(Q)MSSM DRbar
  3  3    -3.75621024e+01      # Atau(Q)MSSM DRbar
