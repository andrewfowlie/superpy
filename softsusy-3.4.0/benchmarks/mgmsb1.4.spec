# SuperPy: Jacobian for naturalness priors.
# J = inf
# b = 3.35211396e+04
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
     1    6.00000000e+04   # lambda
     2    1.00000000e+05   # M_mess
     5    3.00000000e+00   # N5
     6    1.00000000e+00   # cgrav
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=2.66873410e-05
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03989120e+01   # MW
        25     1.14110037e+02   # h0
        35     5.59348391e+02   # H0
        36     5.59108643e+02   # A0
        37     5.65053015e+02   # H+
   1000021     1.39075811e+03   # ~g
   1000022     2.55824688e+02   # ~neutralino(1)
   1000023     4.07203265e+02   # ~neutralino(2)
   1000024     4.04086839e+02   # ~chargino(1)
   1000025    -4.37206004e+02   # ~neutralino(3)
   1000035     5.44210311e+02   # ~neutralino(4)
   1000037     5.44075455e+02   # ~chargino(2)
   1000039     1.42200000e-09   # ~gravitino
   1000001     1.31156595e+03   # ~d_L
   1000002     1.30928827e+03   # ~u_L
   1000003     1.31156437e+03   # ~s_L
   1000004     1.30928668e+03   # ~c_L
   1000005     1.25021462e+03   # ~b_1
   1000006     1.17118816e+03   # ~t_1
   1000011     3.88660935e+02   # ~e_L
   1000012     3.80180391e+02   # ~nue_L
   1000013     3.88659421e+02   # ~mu_L
   1000014     3.80178846e+02   # ~numu_L
   1000015     1.85140354e+02   # ~stau_1
   1000016     3.79466237e+02   # ~nu_tau_L
   2000001     1.25840200e+03   # ~d_R
   2000002     1.26155741e+03   # ~u_R
   2000003     1.25839979e+03   # ~s_R
   2000004     1.26155628e+03   # ~c_R
   2000005     1.26802361e+03   # ~b_2
   2000006     1.27953766e+03   # ~t_2
   2000011     1.90782068e+02   # ~e_R
   2000013     1.90775811e+02   # ~mu_R
   2000015     3.89420558e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -7.38339629e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.78798766e-01   # N_{1,1}
  1  2    -3.68058036e-02   # N_{1,2}
  1  3     1.70073586e-01   # N_{1,3}
  1  4    -1.08042975e-01   # N_{1,4}
  2  1    -1.90864347e-01   # N_{2,1}
  2  2    -4.82722353e-01   # N_{2,2}
  2  3     6.20688754e-01   # N_{2,3}
  2  4    -5.87618415e-01   # N_{2,4}
  3  1    -4.15506350e-02   # N_{3,1}
  3  2     5.52856883e-02   # N_{3,2}
  3  3     7.01620739e-01   # N_{3,3}
  3  4     7.09186418e-01   # N_{3,4}
  4  1    -6.16224127e-02   # N_{4,1}
  4  2     8.73251370e-01   # N_{4,2}
  4  3     3.05857459e-01   # N_{4,3}
  4  4    -3.74280560e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1    -4.29949350e-01   # U_{1,1}
  1  2     9.02853009e-01   # U_{1,2}
  2  1    -9.02853009e-01   # U_{2,1}
  2  2    -4.29949350e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1    -5.26611510e-01   # V_{1,1}
  1  2     8.50106062e-01   # V_{1,2}
  2  1    -8.50106062e-01   # V_{2,1}
  2  2    -5.26611510e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     2.64092776e-01   # F_{11}
  1  2     9.64497281e-01   # F_{12}
  2  1     9.64497281e-01   # F_{21}
  2  2    -2.64092776e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     4.15408091e-01   # F_{11}
  1  2     9.09635156e-01   # F_{12}
  2  1     9.09635156e-01   # F_{21}
  2  2    -4.15408091e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     9.82154832e-02   # F_{11}
  1  2     9.95165172e-01   # F_{12}
  2  1     9.95165172e-01   # F_{21}
  2  2    -9.82154832e-02   # F_{22}
Block gauge Q= 1.18797304e+03  # SM gauge couplings
     1     3.63457025e-01   # g'(Q)MSSM DRbar
     2     6.42920981e-01   # g(Q)MSSM DRbar
     3     1.04797875e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.18797304e+03  
  3  3     8.55355249e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.18797304e+03  
  3  3     2.02428660e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.18797304e+03  
  3  3     1.51696910e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.18797304e+03 # Higgs mixing parameters
     1     4.28302464e+02    # mu(Q)MSSM DRbar
     2     1.44725627e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43490948e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     3.57068062e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.18797304e+03  # MSSM DRbar SUSY breaking parameters
     1     2.66937366e+02      # M_1(Q)
     2     4.99442100e+02      # M_2(Q)
     3     1.33309543e+03      # M_3(Q)
    21     1.25658111e+05      # mH1^2(Q)
    22    -1.57187013e+05      # mH2^2(Q)
    31     3.78269079e+02      # meL(Q)
    32     3.78267530e+02      # mmuL(Q)
    33     3.77790177e+02      # mtauL(Q)
    34     1.78696594e+02      # meR(Q)
    35     1.78689931e+02      # mmuR(Q)
    36     1.76624759e+02      # mtauR(Q)
    41     1.26284555e+03      # mqL1(Q)
    42     1.26284392e+03      # mqL2(Q)
    43     1.22428826e+03      # mqL3(Q)
    44     1.21596303e+03      # muR(Q)
    45     1.21596187e+03      # mcR(Q)
    46     1.13761953e+03      # mtR(Q)
    47     1.21167688e+03      # mdR(Q)
    48     1.21167461e+03      # msR(Q)
    49     1.20728331e+03      # mbR(Q)
Block au Q= 1.18797304e+03  
  1  1    -3.84515131e+02      # Au(Q)MSSM DRbar
  2  2    -3.84514612e+02      # Ac(Q)MSSM DRbar
  3  3    -3.63468538e+02      # At(Q)MSSM DRbar
Block ad Q= 1.18797304e+03  
  1  1    -4.07920614e+02      # Ad(Q)MSSM DRbar
  2  2    -4.07919892e+02      # As(Q)MSSM DRbar
  3  3    -3.99995804e+02      # Ab(Q)MSSM DRbar
Block ae Q= 1.18797304e+03  
  1  1    -4.02928046e+01      # Ae(Q)MSSM DRbar
  2  2    -4.02925219e+01      # Amu(Q)MSSM DRbar
  3  3    -4.02054065e+01      # Atau(Q)MSSM DRbar
