# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dmu=-3.50136476e+03
# dtanbeta/dmu=2.83114442e-02
# mu=8.00390928e+02
# SuperPy. Derivatives for naturalness priors.
# dMZ^2/dm3sq=4.43263913e-01
# dtanbeta/dm3sq=-8.97292783e-05
# m3sq= 9.24077652e+04
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
     1    4.50000000e+02   # m0
     2    6.75000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.77127428e+16   # MX scale
# SOFTSUSY-specific non SLHA information:
# MIXING=0 Desired accuracy=1.00000000e-03 Achieved accuracy=1.64596919e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.03930678e+01   # MW
        25     1.16892796e+02   # h0
        35     1.02955355e+03   # H0
        36     1.02929761e+03   # A0
        37     1.03273363e+03   # H+
   1000021     1.52386690e+03   # ~g
   1000022     2.81702438e+02   # ~neutralino(1)
   1000023     5.33076780e+02   # ~neutralino(2)
   1000024     5.33161529e+02   # ~chargino(1)
   1000025    -8.26025010e+02   # ~neutralino(3)
   1000035     8.37848356e+02   # ~neutralino(4)
   1000037     8.38034734e+02   # ~chargino(2)
   1000001     1.44819379e+03   # ~d_L
   1000002     1.44617583e+03   # ~u_L
   1000003     1.44819007e+03   # ~s_L
   1000004     1.44617210e+03   # ~c_L
   1000005     1.31711735e+03   # ~b_1
   1000006     1.10245501e+03   # ~t_1
   1000011     6.37623878e+02   # ~e_L
   1000012     6.32479349e+02   # ~nue_L
   1000013     6.37616851e+02   # ~mu_L
   1000014     6.32472270e+02   # ~numu_L
   1000015     5.10122063e+02   # ~stau_1
   1000016     6.30169455e+02   # ~nu_tau_L
   2000001     1.39075765e+03   # ~d_R
   2000002     1.39574628e+03   # ~u_R
   2000003     1.39075381e+03   # ~s_R
   2000004     1.39574226e+03   # ~c_R
   2000005     1.38423098e+03   # ~b_2
   2000006     1.34936238e+03   # ~t_2
   2000011     5.17264789e+02   # ~e_R
   2000013     5.17247246e+02   # ~mu_R
   2000015     6.36568209e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.05289518e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.97502163e-01   # N_{1,1}
  1  2    -1.03313721e-02   # N_{1,2}
  1  3     6.41238012e-02   # N_{1,3}
  1  4    -2.77639282e-02   # N_{1,4}
  2  1     2.40920824e-02   # N_{2,1}
  2  2     9.78818715e-01   # N_{2,2}
  2  3    -1.66647671e-01   # N_{2,3}
  2  4     1.16456208e-01   # N_{2,4}
  3  1    -2.52423173e-02   # N_{3,1}
  3  2     3.64374034e-02   # N_{3,2}
  3  3     7.05003462e-01   # N_{3,3}
  3  4     7.07817251e-01   # N_{3,4}
  4  1    -6.14152409e-02   # N_{4,1}
  4  2     2.01195183e-01   # N_{4,2}
  4  3     6.86357641e-01   # N_{4,3}
  4  4    -6.96176598e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.71224339e-01   # U_{1,1}
  1  2    -2.38166504e-01   # U_{1,2}
  2  1     2.38166504e-01   # U_{2,1}
  2  2     9.71224339e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.85953607e-01   # V_{1,1}
  1  2    -1.67019414e-01   # V_{1,2}
  2  1     1.67019414e-01   # V_{2,1}
  2  2     9.85953607e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     3.29978107e-01   # F_{11}
  1  2     9.43988585e-01   # F_{12}
  2  1     9.43988585e-01   # F_{21}
  2  2    -3.29978107e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.92189041e-01   # F_{11}
  1  2     1.24743366e-01   # F_{12}
  2  1    -1.24743366e-01   # F_{21}
  2  2     9.92189041e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.03820750e-01   # F_{11}
  1  2     9.94596024e-01   # F_{12}
  2  1     9.94596024e-01   # F_{21}
  2  2    -1.03820750e-01   # F_{22}
Block gauge Q= 1.18405695e+03  # SM gauge couplings
     1     3.62870428e-01   # g'(Q)MSSM DRbar
     2     6.41220984e-01   # g(Q)MSSM DRbar
     3     1.04544721e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.18405695e+03  
  3  3     8.51410163e-01   # Yt(Q)MSSM DRbar
Block yd Q= 1.18405695e+03  
  3  3     1.33444524e-01   # Yb(Q)MSSM DRbar
Block ye Q= 1.18405695e+03  
  3  3     1.00201855e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 1.18405695e+03 # Higgs mixing parameters
     1     8.20159274e+02    # mu(Q)MSSM DRbar
     2     9.63744620e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.43734992e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     1.09377880e+06    # mA^2(Q)MSSM DRbar
Block msoft Q= 1.18405695e+03  # MSSM DRbar SUSY breaking parameters
     1     2.86443518e+02      # M_1(Q)
     2     5.27790444e+02      # M_2(Q)
     3     1.47073568e+03      # M_3(Q)
    21     3.64916420e+05      # mH1^2(Q)
    22    -6.42942643e+05      # mH2^2(Q)
    31     6.31135489e+02      # meL(Q)
    32     6.31128373e+02      # mmuL(Q)
    33     6.28977203e+02      # mtauL(Q)
    34     5.12681160e+02      # meR(Q)
    35     5.12663464e+02      # mmuR(Q)
    36     5.07294221e+02      # mtauR(Q)
    41     1.40090569e+03      # mqL1(Q)
    42     1.40090187e+03      # mqL2(Q)
    43     1.28028791e+03      # mqL3(Q)
    44     1.35157176e+03      # muR(Q)
    45     1.35156766e+03      # mcR(Q)
    46     1.08876339e+03      # mtR(Q)
    47     1.34553934e+03      # mdR(Q)
    48     1.34553543e+03      # msR(Q)
    49     1.33835426e+03      # mbR(Q)
Block au Q= 1.18405695e+03  
  1  1    -1.49018067e+03      # Au(Q)MSSM DRbar
  2  2    -1.49017407e+03      # Ac(Q)MSSM DRbar
  3  3    -1.15487502e+03      # At(Q)MSSM DRbar
Block ad Q= 1.18405695e+03  
  1  1    -1.81701863e+03      # Ad(Q)MSSM DRbar
  2  2    -1.81701253e+03      # As(Q)MSSM DRbar
  3  3    -1.69933298e+03      # Ab(Q)MSSM DRbar
Block ae Q= 1.18405695e+03  
  1  1    -3.98528285e+02      # Ae(Q)MSSM DRbar
  2  2    -3.98521255e+02      # Amu(Q)MSSM DRbar
  3  3    -3.96398295e+02      # Atau(Q)MSSM DRbar
