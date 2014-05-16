# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-3.15576124e+03
# dtanbeta/dmu=-5.63349500e-01
# mu=9.24428175e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=9.68132262e-02
# dtanbeta/dm3sq=-1.93595673e-03
# m3sq= -3.30258886e+05
# SOFTSUSY3.4.0 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.4.0       # version number
Block MODSEL  # Select model
     1    1   # sugra
Block SMINPUTS             # Standard Model inputs
     1    1.27934000e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637000e-05   # G_Fermi
     3    1.17200000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.25000000e+00   # mb(mb)
     6    1.73300000e+02   # Mtop(pole)
     7    1.77700000e+00   # Mtau(pole)
Block MINPAR               # SUSY breaking input parameters
     3    4.00000000e+01   # tanb, DRbar, Feynman gauge
     4    1.00000000e+00   # sign(mu)
     1    3.75000000e+02   # m0
     2    6.50000000e+02   # m12
     5   -5.00000000e+02   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.71624629e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=9.24472884e-05
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03886878e+01   # MW
        25     1.18305576e+02   # h0
        35     7.58115772e+02   # H0
        36     7.58120747e+02   # A0
        37     7.62665582e+02   # H+
   1000021     1.46727872e+03   # ~g
   1000022     2.72324218e+02   # ~neutralino(1)
   1000023     5.18125071e+02   # ~neutralino(2)
   1000024     5.18266350e+02   # ~chargino(1)
   1000025    -8.60813291e+02   # ~neutralino(3)
   1000035     8.68428624e+02   # ~neutralino(4)
   1000037     8.69328718e+02   # ~chargino(2)
   1000001     1.38375681e+03   # ~d_L
   1000002     1.38158897e+03   # ~u_L
   1000003     1.38371872e+03   # ~s_L
   1000004     1.38155080e+03   # ~c_L
   1000005     1.16924269e+03   # ~b_1
   1000006     1.01167541e+03   # ~t_1
   1000011     5.75708761e+02   # ~e_L
   1000012     5.69951387e+02   # ~nue_L
   1000013     5.75567636e+02   # ~mu_L
   1000014     5.69808940e+02   # ~numu_L
   1000015     2.72396715e+02   # ~stau_1
   1000016     5.20576039e+02   # ~nu_tau_L
   2000001     1.32750300e+03   # ~d_R
   2000002     1.33228506e+03   # ~u_R
   2000003     1.32742722e+03   # ~s_R
   2000004     1.33228024e+03   # ~c_R
   2000005     1.23480439e+03   # ~b_2
   2000006     1.23564719e+03   # ~t_2
   2000011     4.48955223e+02   # ~e_R
   2000013     4.48587738e+02   # ~mu_R
   2000015     5.46000158e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.60025225e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.98034185e-01   # N_{1,1}
  1  2    -6.88743766e-03   # N_{1,2}
  1  3     5.88848043e-02   # N_{1,3}
  1  4    -2.03201517e-02   # N_{1,4}
  2  1     1.69708832e-02   # N_{2,1}
  2  2     9.85825166e-01   # N_{2,2}
  2  3    -1.42129624e-01   # N_{2,3}
  2  4     8.75208578e-02   # N_{2,4}
  3  1    -2.69116375e-02   # N_{3,1}
  3  2     3.92671486e-02   # N_{3,2}
  3  3     7.04925384e-01   # N_{3,3}
  3  4     7.07682173e-01   # N_{3,4}
  4  1    -5.39955409e-02   # N_{4,1}
  4  2     1.62970536e-01   # N_{4,2}
  4  3     6.92395806e-01   # N_{4,3}
  4  4    -7.00794644e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.79329383e-01   # U_{1,1}
  1  2    -2.02271992e-01   # U_{1,2}
  2  1     2.02271992e-01   # U_{2,1}
  2  2     9.79329383e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.92159334e-01   # V_{1,1}
  1  2    -1.24979425e-01   # V_{1,2}
  2  1     1.24979425e-01   # V_{2,1}
  2  2     9.92159334e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     4.17867892e-01   # F_{11}
  1  2     9.08507801e-01   # F_{12}
  2  1     9.08507801e-01   # F_{21}
  2  2    -4.17867892e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     8.08911623e-01   # F_{11}
  1  2     5.87930257e-01   # F_{12}
  2  1    -5.87930257e-01   # F_{21}
  2  2     8.08911623e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     3.02959119e-01   # F_{11}
  1  2     9.53003553e-01   # F_{12}
  2  1     9.53003553e-01   # F_{21}
  2  2    -3.02959119e-01   # F_{22}
Block gauge Q= 1.08706473e+03  # SM gauge couplings
     1     3.62685064e-01   # g'(Q)MSSM DRbar
     2     6.41137383e-01   # g(Q)MSSM DRbar
     3     1.04857635e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.08706473e+03  
  3  3     8.48084659e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.08706473e+03  
  3  3     4.89904379e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.08706473e+03  
  3  3     4.27865070e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.08706473e+03 # Higgs mixing parameters
     1     8.57816598e+02    # mu(Q)MSSM DRbar
     2     3.91757375e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43731050e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     7.28517712e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.08706473e+03  # MSSM DRbar SUSY breaking parameters
     1     2.76577979e+02      # M_1(Q)
     2     5.10367908e+02      # M_2(Q)
     3     1.42573129e+03      # M_3(Q)
    21    -1.63092332e+05      # mH1^2(Q)
    22    -7.21531147e+05      # mH2^2(Q)
    31     5.68931038e+02      # meL(Q)
    32     5.68787834e+02      # mmuL(Q)
    33     5.22709238e+02      # mtauL(Q)
    34     4.43880501e+02      # meR(Q)
    35     4.43508773e+02      # mmuR(Q)
    36     3.07313024e+02      # mtauR(Q)
    41     1.34026217e+03      # mqL1(Q)
    42     1.34022334e+03      # mqL2(Q)
    43     1.16072027e+03      # mqL3(Q)
    44     1.29200584e+03      # muR(Q)
    45     1.29200095e+03      # mcR(Q)
    46     1.01384746e+03      # mtR(Q)
    47     1.28614727e+03      # mdR(Q)
    48     1.28607026e+03      # msR(Q)
    49     1.17580537e+03      # mbR(Q)
Block au Q= 1.08706473e+03  
  1  1    -1.78771102e+03      # Au(Q)MSSM DRbar
  2  2    -1.78766910e+03      # Ac(Q)MSSM DRbar
  3  3    -1.25732017e+03      # At(Q)MSSM DRbar
Block ad Q= 1.08706473e+03  
  1  1    -2.06510903e+03      # Ad(Q)MSSM DRbar
  2  2    -2.06500229e+03      # As(Q)MSSM DRbar
  3  3    -1.74605091e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.08706473e+03  
  1  1    -6.86158575e+02      # Ae(Q)MSSM DRbar
  2  2    -6.85825543e+02      # Amu(Q)MSSM DRbar
  3  3    -5.74035200e+02      # Atau(Q)MSSM DRbar
