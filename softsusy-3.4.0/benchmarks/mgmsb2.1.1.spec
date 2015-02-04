# SuperPy: Jacobian for naturalness priors.
# J = inf
# b = 1.33941584e+04
# Mu = 9.11876000e+01
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
     1    7.20000000e+04   # lambda
     2    8.00000000e+04   # M_mess
     5    1.00000000e+00   # N5
     6    1.00000000e+00   # cgrav
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.16707875e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.04112039e+01   # MW
        25     1.10609543e+02   # h0
        35     3.73635936e+02   # H0
        36     3.73219731e+02   # A0
        37     3.82168990e+02   # H+
   1000021     7.07369845e+02   # ~g
   1000022     1.14265585e+02   # ~neutralino(1)
   1000023     2.06541896e+02   # ~neutralino(2)
   1000024     2.05380648e+02   # ~chargino(1)
   1000025    -3.04618852e+02   # ~neutralino(3)
   1000035     3.38713696e+02   # ~neutralino(4)
   1000037     3.38982927e+02   # ~chargino(2)
   1000039     1.36512000e-09   # ~gravitino
   1000001     8.40579509e+02   # ~d_L
   1000002     8.36984993e+02   # ~u_L
   1000003     8.40578357e+02   # ~s_L
   1000004     8.36983836e+02   # ~c_L
   1000005     7.95991395e+02   # ~b_1
   1000006     7.47813634e+02   # ~t_1
   1000011     2.56424329e+02   # ~e_L
   1000012     2.43498501e+02   # ~nue_L
   1000013     2.56423248e+02   # ~mu_L
   1000014     2.43497365e+02   # ~numu_L
   1000015     1.22807892e+02   # ~stau_1
   1000016     2.42996153e+02   # ~nu_tau_L
   2000001     8.04159211e+02   # ~d_R
   2000002     8.05699322e+02   # ~u_R
   2000003     8.04157608e+02   # ~s_R
   2000004     8.05698506e+02   # ~c_R
   2000005     8.11379283e+02   # ~b_2
   2000006     8.27438075e+02   # ~t_2
   2000011     1.29436453e+02   # ~e_R
   2000013     1.29432125e+02   # ~mu_R
   2000015     2.58355079e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -7.96317228e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.75217145e-01   # N_{1,1}
  1  2    -6.35461325e-02   # N_{1,2}
  1  3     1.93870154e-01   # N_{1,3}
  1  4    -8.56024106e-02   # N_{1,4}
  2  1     1.61784223e-01   # N_{2,1}
  2  2     8.50470065e-01   # N_{2,2}
  2  3    -4.05482990e-01   # N_{2,3}
  2  4     2.93445186e-01   # N_{2,4}
  3  1    -6.88318965e-02   # N_{3,1}
  3  2     9.77347653e-02   # N_{3,2}
  3  3     6.92438772e-01   # N_{3,3}
  3  4     7.11504486e-01   # N_{3,4}
  4  1    -1.34311411e-01   # N_{4,1}
  4  2     5.12942952e-01   # N_{4,2}
  4  3     5.64381480e-01   # N_{4,3}
  4  4    -6.32711243e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     7.99959131e-01   # U_{1,1}
  1  2    -6.00054489e-01   # U_{1,2}
  2  1     6.00054489e-01   # U_{2,1}
  2  2     7.99959131e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.00657260e-01   # V_{1,1}
  1  2    -4.34530206e-01   # V_{1,2}
  2  1     4.34530206e-01   # V_{2,1}
  2  2     9.00657260e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.27989852e-01   # F_{11}
  1  2     9.44681246e-01   # F_{12}
  2  1     9.44681246e-01   # F_{21}
  2  2    -3.27989852e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     5.43480840e-01   # F_{11}
  1  2     8.39421573e-01   # F_{12}
  2  1     8.39421573e-01   # F_{21}
  2  2    -5.43480840e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.54170368e-01   # F_{11}
  1  2     9.88044279e-01   # F_{12}
  2  1     9.88044279e-01   # F_{21}
  2  2    -1.54170368e-01   # F_{22}
Block gauge Q= 7.68693064e+02  # SM gauge couplings
     1     3.62512220e-01   # g'(Q)MSSM DRbar
     2     6.46107286e-01   # g(Q)MSSM DRbar
     3     1.07638115e+00   # g3(Q)MSSM DRbar
Block yu Q= 7.68693064e+02  
  3  3     8.73707037e-01   # Yt(Q)MSSM DRbar
Block yd Q= 7.68693064e+02  
  3  3     2.06920356e-01   # Yb(Q)MSSM DRbar
Block ye Q= 7.68693064e+02  
  3  3     1.51625828e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 7.68693064e+02 # Higgs mixing parameters
     1     2.95115539e+02    # mu(Q)MSSM DRbar
     2     1.45437300e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43895830e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.54862745e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 7.68693064e+02  # MSSM DRbar SUSY breaking parameters
     1     1.21037427e+02      # M_1(Q)
     2     2.29808727e+02      # M_2(Q)
     3     6.41114056e+02      # M_3(Q)
    21     5.42288402e+04      # mH1^2(Q)
    22    -7.85627109e+04      # mH2^2(Q)
    31     2.49846917e+02      # meL(Q)
    32     2.49845806e+02      # mmuL(Q)
    33     2.49507368e+02      # mtauL(Q)
    34     1.18261004e+02      # meR(Q)
    35     1.18256260e+02      # mmuR(Q)
    36     1.16803613e+02      # mtauR(Q)
    41     8.15912995e+02      # mqL1(Q)
    42     8.15911802e+02      # mqL2(Q)
    43     7.87292174e+02      # mqL3(Q)
    44     7.83708481e+02      # muR(Q)
    45     7.83707636e+02      # mcR(Q)
    46     7.25286539e+02      # mtR(Q)
    47     7.80682986e+02      # mdR(Q)
    48     7.80681320e+02      # msR(Q)
    49     7.77420916e+02      # mbR(Q)
Block au Q= 7.68693064e+02  
  1  1    -1.99093016e+02      # Au(Q)MSSM DRbar
  2  2    -1.99092730e+02      # Ac(Q)MSSM DRbar
  3  3    -1.87348569e+02      # At(Q)MSSM DRbar
Block ad Q= 7.68693064e+02  
  1  1    -2.12230137e+02      # Ad(Q)MSSM DRbar
  2  2    -2.12229737e+02      # As(Q)MSSM DRbar
  3  3    -2.07802528e+02      # Ab(Q)MSSM DRbar
Block ae Q= 7.68693064e+02  
  1  1    -1.94808715e+01      # Ae(Q)MSSM DRbar
  2  2    -1.94807272e+01      # Amu(Q)MSSM DRbar
  3  3    -1.94367974e+01      # Atau(Q)MSSM DRbar
