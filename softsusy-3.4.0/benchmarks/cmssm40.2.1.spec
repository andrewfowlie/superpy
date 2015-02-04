# SuperPy: Jacobian for naturalness priors.
# J = 5.39142259e+04
# b = 6.73447033e+04
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
     1    5.50000000e+02   # m0
     2    4.50000000e+02   # m12
     5   -5.00000000e+02   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    2.02197545e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=7.41071404e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03910848e+01   # MW
        25     1.16557282e+02   # h0
        35     6.30939485e+02   # H0
        36     6.30911473e+02   # A0
        37     6.36450164e+02   # H+
   1000021     1.06318928e+03   # ~g
   1000022     1.86425625e+02   # ~neutralino(1)
   1000023     3.55657080e+02   # ~neutralino(2)
   1000024     3.55774157e+02   # ~chargino(1)
   1000025    -6.44992765e+02   # ~neutralino(3)
   1000035     6.53676897e+02   # ~neutralino(4)
   1000037     6.54956723e+02   # ~chargino(2)
   1000001     1.09460464e+03   # ~d_L
   1000002     1.09186764e+03   # ~u_L
   1000003     1.09456964e+03   # ~s_L
   1000004     1.09183254e+03   # ~c_L
   1000005     8.95282839e+02   # ~b_1
   1000006     7.50184654e+02   # ~t_1
   1000011     6.27371898e+02   # ~e_L
   1000012     6.21995122e+02   # ~nue_L
   1000013     6.27202288e+02   # ~mu_L
   1000014     6.21824154e+02   # ~numu_L
   1000015     4.23392991e+02   # ~stau_1
   1000016     5.64287449e+02   # ~nu_tau_L
   2000001     1.06129841e+03   # ~d_R
   2000002     1.06354062e+03   # ~u_R
   2000003     1.06122973e+03   # ~s_R
   2000004     1.06353620e+03   # ~c_R
   2000005     9.69667192e+02   # ~b_2
   2000006     9.60080450e+02   # ~t_2
   2000011     5.76242579e+02   # ~e_R
   2000013     5.75870074e+02   # ~mu_R
   2000015     5.84034103e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -2.60620201e-02       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.96575587e-01   # N_{1,1}
  1  2    -1.21175620e-02   # N_{1,2}
  1  3     7.79436564e-02   # N_{1,3}
  1  4    -2.48002163e-02   # N_{1,4}
  2  1     2.78881840e-02   # N_{2,1}
  2  2     9.79609997e-01   # N_{2,2}
  2  3    -1.73022298e-01   # N_{2,3}
  2  4     9.82333343e-02   # N_{2,4}
  3  1    -3.67057599e-02   # N_{3,1}
  3  2     5.42119517e-02   # N_{3,2}
  3  3     7.03052257e-01   # N_{3,3}
  3  4     7.08118123e-01   # N_{3,4}
  4  1    -6.86442714e-02   # N_{4,1}
  4  2     1.93076367e-01   # N_{4,2}
  4  3     6.85350710e-01   # N_{4,3}
  4  4    -6.98787439e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.69100494e-01   # U_{1,1}
  1  2    -2.46666236e-01   # U_{1,2}
  2  1     2.46666236e-01   # U_{2,1}
  2  2     9.69100494e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.90069020e-01   # V_{1,1}
  1  2    -1.40582130e-01   # V_{1,2}
  2  1     1.40582130e-01   # V_{2,1}
  2  2     9.90069020e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     4.51814689e-01   # F_{11}
  1  2     8.92111813e-01   # F_{12}
  2  1     8.92111813e-01   # F_{21}
  2  2    -4.51814689e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     8.88887293e-01   # F_{11}
  1  2     4.58125944e-01   # F_{12}
  2  1    -4.58125944e-01   # F_{21}
  2  2     8.88887293e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     3.13495695e-01   # F_{11}
  1  2     9.49589622e-01   # F_{12}
  2  1     9.49589622e-01   # F_{21}
  2  2    -3.13495695e-01   # F_{22}
Block gauge Q= 8.23117242e+02  # SM gauge couplings
     1     3.61863494e-01   # g'(Q)MSSM DRbar
     2     6.42453061e-01   # g(Q)MSSM DRbar
     3     1.06387661e+00   # g3(Q)MSSM DRbar
Block yu Q= 8.23117242e+02  
  3  3     8.58488303e-01   # Yt(Q)MSSM DRbar
Block yd Q= 8.23117242e+02  
  3  3     4.97304720e-01   # Yb(Q)MSSM DRbar
Block ye Q= 8.23117242e+02  
  3  3     4.23565555e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 8.23117242e+02 # Higgs mixing parameters
     1     6.40358662e+02    # mu(Q)MSSM DRbar
     2     3.92490005e+01    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.44080520e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     5.09552386e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 8.23117242e+02  # MSSM DRbar SUSY breaking parameters
     1     1.88921430e+02      # M_1(Q)
     2     3.51467000e+02      # M_2(Q)
     3     1.01087160e+03      # M_3(Q)
    21    -3.01384451e+03      # mH1^2(Q)
    22    -4.04645265e+05      # mH2^2(Q)
    31     6.23756061e+02      # meL(Q)
    32     6.23585513e+02      # mmuL(Q)
    33     5.68121813e+02      # mtauL(Q)
    34     5.73354028e+02      # meR(Q)
    35     5.72980047e+02      # mmuR(Q)
    36     4.41258926e+02      # mtauR(Q)
    41     1.06299423e+03      # mqL1(Q)
    42     1.06295797e+03      # mqL2(Q)
    43     8.85958081e+02      # mqL3(Q)
    44     1.03452788e+03      # muR(Q)
    45     1.03452331e+03      # mcR(Q)
    46     7.58311527e+02      # mtR(Q)
    47     1.03117994e+03      # mdR(Q)
    48     1.03110893e+03      # msR(Q)
    49     9.26480711e+02      # mbR(Q)
Block au Q= 8.23117242e+02  
  1  1    -1.37478656e+03      # Au(Q)MSSM DRbar
  2  2    -1.37475216e+03      # Ac(Q)MSSM DRbar
  3  3    -9.37347461e+02      # At(Q)MSSM DRbar
Block ad Q= 8.23117242e+02  
  1  1    -1.60153425e+03      # Ad(Q)MSSM DRbar
  2  2    -1.60144665e+03      # As(Q)MSSM DRbar
  3  3    -1.33779444e+03      # Ab(Q)MSSM DRbar
Block ae Q= 8.23117242e+02  
  1  1    -6.02572760e+02      # Ae(Q)MSSM DRbar
  2  2    -6.02265536e+02      # Amu(Q)MSSM DRbar
  3  3    -5.01609859e+02      # Atau(Q)MSSM DRbar
