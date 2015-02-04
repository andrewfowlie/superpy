# SuperPy: Jacobian for naturalness priors.
# J = 5.25890889e+03
# b = 1.57702828e+04
# Mu = 9.11876000e+01
# SOFTSUSY3.4.0 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.4.0       # version number
Block MODSEL  # Select model
     1    3   # amsb
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
     1    4.50000000e+02   # m0
     2    6.00000000e+04   # m3/2
Block EXTPAR               # scale of SUSY breaking BCs
     0    2.30022864e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.20212975e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03918518e+01   # MW
        25     1.15406705e+02   # h0
        35     1.06453297e+03   # H0
        36     1.06437187e+03   # A0
        37     1.06764246e+03   # H+
   1000021    -1.28226775e+03   # ~g
   1000022     1.97416072e+02   # ~neutralino(1)
   1000023     5.49091484e+02   # ~neutralino(2)
   1000024     1.97585917e+02   # ~chargino(1)
   1000025    -1.02069863e+03   # ~neutralino(3)
   1000035     1.02661882e+03   # ~neutralino(4)
   1000037     1.02559434e+03   # ~chargino(2)
   1000039     6.00000000e+04   # ~gravitino
   1000001     1.27602355e+03   # ~d_L
   1000002     1.27371975e+03   # ~u_L
   1000003     1.27601582e+03   # ~s_L
   1000004     1.27371200e+03   # ~c_L
   1000005     1.10853326e+03   # ~b_1
   1000006     9.26016590e+02   # ~t_1
   1000011     3.87388640e+02   # ~e_L
   1000012     3.79011610e+02   # ~nue_L
   1000013     3.87371796e+02   # ~mu_L
   1000014     3.78994413e+02   # ~numu_L
   1000015     3.48477511e+02   # ~stau_1
   1000016     3.73707966e+02   # ~nu_tau_L
   2000001     1.28948825e+03   # ~d_R
   2000002     1.27993113e+03   # ~u_R
   2000003     1.28947997e+03   # ~s_R
   2000004     1.27992404e+03   # ~c_R
   2000005     1.27225775e+03   # ~b_2
   2000006     1.14253169e+03   # ~t_2
   2000011     3.74881416e+02   # ~e_R
   2000013     3.74846514e+02   # ~mu_R
   2000015     3.96445294e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05614165e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1    -3.76985466e-03   # N_{1,1}
  1  2     9.96454350e-01   # N_{1,2}
  1  3    -8.08618303e-02   # N_{1,3}
  1  4     2.29321027e-02   # N_{1,4}
  2  1     9.97183950e-01   # N_{2,1}
  2  2     9.81478551e-03   # N_{2,2}
  2  3     6.34808239e-02   # N_{2,3}
  2  4    -3.87043248e-02   # N_{2,4}
  3  1    -1.78544014e-02   # N_{3,1}
  3  2     4.08983472e-02   # N_{3,2}
  3  3     7.05428476e-01   # N_{3,3}
  3  4     7.07374873e-01   # N_{3,4}
  4  1    -7.27404861e-02   # N_{4,1}
  4  2     7.28678524e-02   # N_{4,2}
  4  3     7.01286115e-01   # N_{4,3}
  4  4    -7.05405474e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.93500804e-01   # U_{1,1}
  1  2    -1.13825094e-01   # U_{1,2}
  2  1     1.13825094e-01   # U_{2,1}
  2  2     9.93500804e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.99477490e-01   # V_{1,1}
  1  2    -3.23225461e-02   # V_{1,2}
  2  1     3.23225461e-02   # V_{2,1}
  2  2     9.99477490e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1    -3.61119852e-01   # F_{11}
  1  2     9.32519411e-01   # F_{12}
  2  1     9.32519411e-01   # F_{21}
  2  2     3.61119852e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.98813191e-01   # F_{11}
  1  2     4.87053376e-02   # F_{12}
  2  1    -4.87053376e-02   # F_{21}
  2  2     9.98813191e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     5.57692407e-01   # F_{11}
  1  2     8.30047697e-01   # F_{12}
  2  1     8.30047697e-01   # F_{21}
  2  2    -5.57692407e-01   # F_{22}
Block gauge Q= 9.91542975e+02  # SM gauge couplings
     1     3.62517713e-01   # g'(Q)MSSM DRbar
     2     6.44082625e-01   # g(Q)MSSM DRbar
     3     1.05399668e+00   # g3(Q)MSSM DRbar
Block yu Q= 9.91542975e+02  
  3  3     8.57777944e-01   # Yt(Q)MSSM DRbar
Block yd Q= 9.91542975e+02  
  3  3     1.47581614e-01   # Yb(Q)MSSM DRbar
Block ye Q= 9.91542975e+02  
  3  3     9.92462253e-02   # Ytau(Q)MSSM DRbar
Block hmix Q= 9.91542975e+02 # Higgs mixing parameters
     1     1.02020169e+03    # mu(Q)MSSM DRbar
     2     9.65972704e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.44235971e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.13330253e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 9.91542975e+02  # MSSM DRbar SUSY breaking parameters
     1     5.59602231e+02      # M_1(Q)
     2     1.90422220e+02      # M_2(Q)
     3    -1.22154507e+03      # M_3(Q)
    21     7.36616811e+04      # mH1^2(Q)
    22    -1.01740341e+06      # mH2^2(Q)
    31     3.83166359e+02      # meL(Q)
    32     3.83149343e+02      # mmuL(Q)
    33     3.78142931e+02      # mtauL(Q)
    34     3.67230901e+02      # meR(Q)
    35     3.67195378e+02      # mmuR(Q)
    36     3.56660352e+02      # mtauR(Q)
    41     1.23783675e+03      # mqL1(Q)
    42     1.23782877e+03      # mqL2(Q)
    43     1.07285286e+03      # mqL3(Q)
    44     1.24507893e+03      # muR(Q)
    45     1.24507161e+03      # mcR(Q)
    46     9.07645894e+02      # mtR(Q)
    47     1.25424937e+03      # mdR(Q)
    48     1.25424080e+03      # msR(Q)
    49     1.23650688e+03      # mbR(Q)
Block au Q= 9.91542975e+02  
  1  1     1.92812767e+03      # Au(Q)MSSM DRbar
  2  2     1.92811140e+03      # Ac(Q)MSSM DRbar
  3  3     1.09361033e+03      # At(Q)MSSM DRbar
Block ad Q= 9.91542975e+02  
  1  1     2.72566757e+03      # Ad(Q)MSSM DRbar
  2  2     2.72565252e+03      # As(Q)MSSM DRbar
  3  3     2.42743392e+03      # Ab(Q)MSSM DRbar
Block ae Q= 9.91542975e+02  
  1  1     5.88174849e+02      # Ae(Q)MSSM DRbar
  2  2     5.88136919e+02      # Amu(Q)MSSM DRbar
  3  3     5.76926619e+02      # Atau(Q)MSSM DRbar
