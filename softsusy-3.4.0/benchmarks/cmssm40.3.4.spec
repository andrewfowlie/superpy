# SuperPy: Jacobian for naturalness priors.
# J = 2.26427561e+04
# b = 7.70026336e+04
# Mu = 9.11876000e+01
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
     1    1.15000000e+03   # m0
     2    4.25000000e+02   # m12
     5   -5.00000000e+02   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    2.12821336e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=3.39353755e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03913921e+01   # MW
        25     1.16569409e+02   # h0
        35     9.02126219e+02   # H0
        36     9.02235834e+02   # A0
        37     9.05965595e+02   # H+
   1000021     1.04501969e+03   # ~g
   1000022     1.77577883e+02   # ~neutralino(1)
   1000023     3.39304280e+02   # ~neutralino(2)
   1000024     3.39407075e+02   # ~chargino(1)
   1000025    -5.96236956e+02   # ~neutralino(3)
   1000035     6.06779432e+02   # ~neutralino(4)
   1000037     6.07957039e+02   # ~chargino(2)
   1000001     1.43785735e+03   # ~d_L
   1000002     1.43578685e+03   # ~u_L
   1000003     1.43780732e+03   # ~s_L
   1000004     1.43573676e+03   # ~c_L
   1000005     1.12660414e+03   # ~b_1
   1000006     9.33638416e+02   # ~t_1
   1000011     1.18111797e+03   # ~e_L
   1000012     1.17810028e+03   # ~nue_L
   1000013     1.18082683e+03   # ~mu_L
   1000014     1.17780847e+03   # ~numu_L
   1000015     9.52396766e+02   # ~stau_1
   1000016     1.08259900e+03   # ~nu_tau_L
   2000001     1.41775108e+03   # ~d_R
   2000002     1.41882954e+03   # ~u_R
   2000003     1.41765548e+03   # ~s_R
   2000004     1.41882302e+03   # ~c_R
   2000005     1.25096650e+03   # ~b_2
   2000006     1.15882086e+03   # ~t_2
   2000011     1.16036966e+03   # ~e_R
   2000013     1.15977477e+03   # ~mu_R
   2000015     1.08902162e+03   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.59457841e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.95885729e-01   # N_{1,1}
  1  2    -1.38329293e-02   # N_{1,2}
  1  3     8.51460407e-02   # N_{1,3}
  1  4    -2.77563787e-02   # N_{1,4}
  2  1     3.30694423e-02   # N_{2,1}
  2  2     9.74414736e-01   # N_{2,2}
  2  3    -1.91889643e-01   # N_{2,3}
  2  4     1.12252832e-01   # N_{2,4}
  3  1    -3.95176753e-02   # N_{3,1}
  3  2     5.81001487e-02   # N_{3,2}
  3  3     7.02508333e-01   # N_{3,3}
  3  4     7.08198255e-01   # N_{3,4}
  4  1    -7.45411300e-02   # N_{4,1}
  4  2     2.16677051e-01   # N_{4,2}
  4  3     6.80007763e-01   # N_{4,3}
  4  4    -6.96479804e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.61785096e-01   # U_{1,1}
  1  2    -2.73805460e-01   # U_{1,2}
  2  1     2.73805460e-01   # U_{2,1}
  2  2     9.61785096e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.86996949e-01   # V_{1,1}
  1  2    -1.60738989e-01   # V_{1,2}
  2  1     1.60738989e-01   # V_{2,1}
  2  2     9.86996949e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.06514029e-01   # F_{11}
  1  2     9.51866141e-01   # F_{12}
  2  1     9.51866141e-01   # F_{21}
  2  2    -3.06514029e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.83046197e-01   # F_{11}
  1  2     1.83358052e-01   # F_{12}
  2  1    -1.83358052e-01   # F_{21}
  2  2     9.83046197e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.58207136e-01   # F_{11}
  1  2     9.87405946e-01   # F_{12}
  2  1     9.87405946e-01   # F_{21}
  2  2    -1.58207136e-01   # F_{22}
Block gauge Q= 1.01597243e+03  # SM gauge couplings
     1     3.62147645e-01   # g'(Q)MSSM DRbar
     2     6.41963771e-01   # g(Q)MSSM DRbar
     3     1.05665861e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.01597243e+03  
  3  3     8.55903590e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.01597243e+03  
  3  3     5.16761766e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.01597243e+03  
  3  3     4.12747408e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.01597243e+03 # Higgs mixing parameters
     1     5.88385690e+02    # mu(Q)MSSM DRbar
     2     3.92188422e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43665976e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     9.80407253e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.01597243e+03  # MSSM DRbar SUSY breaking parameters
     1     1.79435908e+02      # M_1(Q)
     2     3.32931664e+02      # M_2(Q)
     3     9.45470536e+02      # M_3(Q)
    21     4.90056325e+05      # mH1^2(Q)
    22    -3.31661565e+05      # mH2^2(Q)
    31     1.17826158e+03      # meL(Q)
    32     1.17797059e+03      # mmuL(Q)
    33     1.08516356e+03      # mtauL(Q)
    34     1.15825157e+03      # meR(Q)
    35     1.15765647e+03      # mmuR(Q)
    36     9.57645291e+02      # mtauR(Q)
    41     1.41174351e+03      # mqL1(Q)
    42     1.41169262e+03      # mqL2(Q)
    43     1.10916229e+03      # mqL3(Q)
    44     1.39583395e+03      # muR(Q)
    45     1.39582732e+03      # mcR(Q)
    46     9.19974302e+02      # mtR(Q)
    47     1.39406890e+03      # mdR(Q)
    48     1.39397149e+03      # msR(Q)
    49     1.22581549e+03      # mbR(Q)
Block au Q= 1.01597243e+03  
  1  1    -1.30190658e+03      # Au(Q)MSSM DRbar
  2  2    -1.30187395e+03      # Ac(Q)MSSM DRbar
  3  3    -8.79626040e+02      # At(Q)MSSM DRbar
Block ad Q= 1.01597243e+03  
  1  1    -1.50810296e+03      # Ad(Q)MSSM DRbar
  2  2    -1.50801989e+03      # As(Q)MSSM DRbar
  3  3    -1.24469241e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.01597243e+03  
  1  1    -5.82671435e+02      # Ae(Q)MSSM DRbar
  2  2    -5.82369937e+02      # Amu(Q)MSSM DRbar
  3  3    -4.88523454e+02      # Atau(Q)MSSM DRbar
