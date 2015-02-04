# SuperPy: Jacobian for naturalness priors.
# J = 1.20527899e+05
# b = 1.42769266e+05
# Mu = 9.11876000e+01
# SOFTSUSY3.4.0 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.4.0       # version number
     3   # Warning: stau LSP
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
     1    3.90000000e+02   # m0
     2    7.00000000e+02   # m12
     5   -5.00000000e+02   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.67496263e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.20638183e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03883200e+01   # MW
        25     1.18668458e+02   # h0
        35     8.06859619e+02   # H0
        36     8.06862947e+02   # A0
        37     8.11134333e+02   # H+
   1000021     1.57094184e+03   # ~g
   1000022     2.94221603e+02   # ~neutralino(1)
   1000023     5.59329036e+02   # ~neutralino(2)
   1000024     5.59473256e+02   # ~chargino(1)
   1000025    -9.12542128e+02   # ~neutralino(3)
   1000035     9.20022537e+02   # ~neutralino(4)
   1000037     9.20855440e+02   # ~chargino(2)
   1000001     1.47782236e+03   # ~d_L
   1000002     1.47580169e+03   # ~u_L
   1000003     1.47778264e+03   # ~s_L
   1000004     1.47576191e+03   # ~c_L
   1000005     1.25367713e+03   # ~b_1
   1000006     1.08845475e+03   # ~t_1
   1000011     6.10397271e+02   # ~e_L
   1000012     6.04954973e+02   # ~nue_L
   1000013     6.10253714e+02   # ~mu_L
   1000014     6.04810221e+02   # ~numu_L
   1000015     2.93077860e+02   # ~stau_1
   1000016     5.54505766e+02   # ~nu_tau_L
   2000001     1.41690713e+03   # ~d_R
   2000002     1.42222740e+03   # ~u_R
   2000003     1.41682809e+03   # ~s_R
   2000004     1.42222236e+03   # ~c_R
   2000005     1.31847809e+03   # ~b_2
   2000006     1.31752722e+03   # ~t_2
   2000011     4.71619681e+02   # ~e_R
   2000013     4.71242306e+02   # ~mu_R
   2000015     5.78488590e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.59709267e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.98236697e-01   # N_{1,1}
  1  2    -6.14358629e-03   # N_{1,2}
  1  3     5.57096641e-02   # N_{1,3}
  1  4    -1.95495677e-02   # N_{1,4}
  2  1     1.53966010e-02   # N_{2,1}
  2  2     9.86734844e-01   # N_{2,2}
  2  3    -1.36971642e-01   # N_{2,3}
  2  4     8.57674847e-02   # N_{2,4}
  3  1    -2.52706646e-02   # N_{3,1}
  3  2     3.67768258e-02   # N_{3,2}
  3  3     7.05185656e-01   # N_{3,3}
  3  4     7.07617163e-01   # N_{3,4}
  4  1    -5.14571095e-02   # N_{4,1}
  4  2     1.58000218e-01   # N_{4,2}
  4  3     6.93432328e-01   # N_{4,3}
  4  4    -7.01098926e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.80821423e-01   # U_{1,1}
  1  2    -1.94908534e-01   # U_{1,2}
  2  1     1.94908534e-01   # U_{2,1}
  2  2     9.80821423e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.92473798e-01   # V_{1,1}
  1  2    -1.22457175e-01   # V_{1,2}
  2  1     1.22457175e-01   # V_{2,1}
  2  2     9.92473798e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     4.01228718e-01   # F_{11}
  1  2     9.15977901e-01   # F_{12}
  2  1     9.15977901e-01   # F_{21}
  2  2    -4.01228718e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     8.06978492e-01   # F_{11}
  1  2     5.90580827e-01   # F_{12}
  2  1    -5.90580827e-01   # F_{21}
  2  2     8.06978492e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     2.87747645e-01   # F_{11}
  1  2     9.57706267e-01   # F_{12}
  2  1     9.57706267e-01   # F_{21}
  2  2    -2.87747645e-01   # F_{22}
Block gauge Q= 1.16478036e+03  # SM gauge couplings
     1     3.62839382e-01   # g'(Q)MSSM DRbar
     2     6.40796362e-01   # g(Q)MSSM DRbar
     3     1.04507689e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.16478036e+03  
  3  3     8.45800851e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.16478036e+03  
  3  3     4.89469472e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.16478036e+03  
  3  3     4.27982909e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.16478036e+03 # Higgs mixing parameters
     1     9.09602922e+02    # mu(Q)MSSM DRbar
     2     3.91596198e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43656116e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     8.26567305e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.16478036e+03  # MSSM DRbar SUSY breaking parameters
     1     2.98760417e+02      # M_1(Q)
     2     5.50316849e+02      # M_2(Q)
     3     1.52737185e+03      # M_3(Q)
    21    -1.79446242e+05      # mH1^2(Q)
    22    -8.09970358e+05      # mH2^2(Q)
    31     6.03178278e+02      # meL(Q)
    32     6.03032606e+02      # mmuL(Q)
    33     5.56022335e+02      # mtauL(Q)
    34     4.66340535e+02      # meR(Q)
    35     4.65958815e+02      # mmuR(Q)
    36     3.26028829e+02      # mtauR(Q)
    41     1.43151695e+03      # mqL1(Q)
    42     1.43147649e+03      # mqL2(Q)
    43     1.24317383e+03      # mqL3(Q)
    44     1.37920886e+03      # muR(Q)
    45     1.37920375e+03      # mcR(Q)
    46     1.08720239e+03      # mtR(Q)
    47     1.37283118e+03      # mdR(Q)
    48     1.37275088e+03      # msR(Q)
    49     1.25718344e+03      # mbR(Q)
Block au Q= 1.16478036e+03  
  1  1    -1.88756441e+03      # Au(Q)MSSM DRbar
  2  2    -1.88752071e+03      # Ac(Q)MSSM DRbar
  3  3    -1.33467387e+03      # At(Q)MSSM DRbar
Block ad Q= 1.16478036e+03  
  1  1    -2.17630612e+03      # Ad(Q)MSSM DRbar
  2  2    -2.17619485e+03      # As(Q)MSSM DRbar
  3  3    -1.84323054e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.16478036e+03  
  1  1    -7.06310052e+02      # Ae(Q)MSSM DRbar
  2  2    -7.05970662e+02      # Amu(Q)MSSM DRbar
  3  3    -5.91870605e+02      # Atau(Q)MSSM DRbar
