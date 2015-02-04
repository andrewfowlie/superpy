# SuperPy: Jacobian for naturalness priors.
# J = 4.01167492e+04
# b = 1.28531915e+05
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
     3    1.00000000e+01   # tanb, DRbar, Feynman gauge
     4    1.00000000e+00   # sign(mu)
     1    2.50000000e+02   # m0
     2    6.00000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.79640086e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.10524604e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03941527e+01   # MW
        25     1.16111784e+02   # h0
        35     8.71235995e+02   # H0
        36     8.70996582e+02   # A0
        37     8.74940543e+02   # H+
   1000021     1.35966331e+03   # ~g
   1000022     2.48142184e+02   # ~neutralino(1)
   1000023     4.69490826e+02   # ~neutralino(2)
   1000024     4.69553690e+02   # ~chargino(1)
   1000025    -7.45756848e+02   # ~neutralino(3)
   1000035     7.58194664e+02   # ~neutralino(4)
   1000037     7.58435563e+02   # ~chargino(2)
   1000001     1.26232745e+03   # ~d_L
   1000002     1.25998342e+03   # ~u_L
   1000003     1.26232433e+03   # ~s_L
   1000004     1.25998030e+03   # ~c_L
   1000005     1.15391482e+03   # ~b_1
   1000006     9.66350523e+02   # ~t_1
   1000011     4.75467735e+02   # ~e_L
   1000012     4.68672709e+02   # ~nue_L
   1000013     4.75463063e+02   # ~mu_L
   1000014     4.68667971e+02   # ~numu_L
   1000015     3.32173805e+02   # ~stau_1
   1000016     4.67065753e+02   # ~nu_tau_L
   2000001     1.20942982e+03   # ~d_R
   2000002     1.21380017e+03   # ~u_R
   2000003     1.20942657e+03   # ~s_R
   2000004     1.21379682e+03   # ~c_R
   2000005     1.20456791e+03   # ~b_2
   2000006     1.19225806e+03   # ~t_2
   2000011     3.38923361e+02   # ~e_R
   2000013     3.38910056e+02   # ~mu_R
   2000015     4.75570939e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05842778e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.96958087e-01   # N_{1,1}
  1  2    -1.26834113e-02   # N_{1,2}
  1  3     7.07619831e-02   # N_{1,3}
  1  4    -3.01072218e-02   # N_{1,4}
  2  1     2.87255093e-02   # N_{2,1}
  2  2     9.76042895e-01   # N_{2,2}
  2  3    -1.77865455e-01   # N_{2,3}
  2  4     1.21979473e-01   # N_{2,4}
  3  1    -2.81046190e-02   # N_{3,1}
  3  2     4.07279267e-02   # N_{3,2}
  3  3     7.04490995e-01   # N_{3,3}
  3  4     7.07985737e-01   # N_{3,4}
  4  1    -6.67798434e-02   # N_{4,1}
  4  2     2.13355651e-01   # N_{4,2}
  4  3     6.83409877e-01   # N_{4,3}
  4  4    -6.94960977e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.67163581e-01   # U_{1,1}
  1  2    -2.54154694e-01   # U_{1,2}
  2  1     2.54154694e-01   # U_{2,1}
  2  2     9.67163581e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.84594265e-01   # V_{1,1}
  1  2    -1.74854608e-01   # V_{1,2}
  2  1     1.74854608e-01   # V_{2,1}
  2  2     9.84594265e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.76440821e-01   # F_{11}
  1  2     9.26440667e-01   # F_{12}
  2  1     9.26440667e-01   # F_{21}
  2  2    -3.76440821e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.85007101e-01   # F_{11}
  1  2     1.72513799e-01   # F_{12}
  2  1    -1.72513799e-01   # F_{21}
  2  2     9.85007101e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.17727051e-01   # F_{11}
  1  2     9.93045992e-01   # F_{12}
  2  1     9.93045992e-01   # F_{21}
  2  2    -1.17727051e-01   # F_{22}
Block gauge Q= 1.04161479e+03  # SM gauge couplings
     1     3.62710172e-01   # g'(Q)MSSM DRbar
     2     6.41986406e-01   # g(Q)MSSM DRbar
     3     1.05179422e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.04161479e+03  
  3  3     8.55317355e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.04161479e+03  
  3  3     1.34167848e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.04161479e+03  
  3  3     1.00383725e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.04161479e+03 # Higgs mixing parameters
     1     7.40241401e+02    # mu(Q)MSSM DRbar
     2     9.65267456e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43896717e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     7.85371268e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.04161479e+03  # MSSM DRbar SUSY breaking parameters
     1     2.53210588e+02      # M_1(Q)
     2     4.67952583e+02      # M_2(Q)
     3     1.31890253e+03      # M_3(Q)
    21     1.95209038e+05      # mH1^2(Q)
    22    -5.27291667e+05      # mH2^2(Q)
    31     4.68054896e+02      # meL(Q)
    32     4.68050140e+02      # mmuL(Q)
    33     4.66612734e+02      # mtauL(Q)
    34     3.32675955e+02      # meR(Q)
    35     3.32662390e+02      # mmuR(Q)
    36     3.28542172e+02      # mtauR(Q)
    41     1.21952883e+03      # mqL1(Q)
    42     1.21952564e+03      # mqL2(Q)
    43     1.12139917e+03      # mqL3(Q)
    44     1.17420464e+03      # muR(Q)
    45     1.17420123e+03      # mcR(Q)
    46     9.60717033e+02      # mtR(Q)
    47     1.16867959e+03      # mdR(Q)
    48     1.16867629e+03      # msR(Q)
    49     1.16264257e+03      # mbR(Q)
Block au Q= 1.04161479e+03  
  1  1    -1.34254409e+03      # Au(Q)MSSM DRbar
  2  2    -1.34253811e+03      # Ac(Q)MSSM DRbar
  3  3    -1.03867288e+03      # At(Q)MSSM DRbar
Block ad Q= 1.04161479e+03  
  1  1    -1.63918489e+03      # Ad(Q)MSSM DRbar
  2  2    -1.63917936e+03      # As(Q)MSSM DRbar
  3  3    -1.53251998e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.04161479e+03  
  1  1    -3.56268858e+02      # Ae(Q)MSSM DRbar
  2  2    -3.56262527e+02      # Amu(Q)MSSM DRbar
  3  3    -3.54349961e+02      # Atau(Q)MSSM DRbar
