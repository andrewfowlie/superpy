# SuperPy: Jacobian for naturalness priors.
# J = 4.30175644e+05
# b = 1.03564365e+05
# Mu = 9.11876000e+01
# SOFTSUSY3.4.0 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
# B.C. Allanach and M.A. Bernhardt, Comput. Phys. Commun. 181 (2010) 232,
# arXiv:0903.1805
# B.C. Allanach, M. Hanussek and C.H. Kom, arXiv:1109.3735
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.4.0       # version number
Block MODSEL  # Select model
     1    1   # sugra
     4    1   # R-parity violating
Block SMINPUTS             # Standard Model inputs
     1    1.27934000e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637000e-05   # G_Fermi
     3    1.17200000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.25000000e+00   # mb(mb)
     6    1.73300000e+02   # Mtop(pole)
     7    1.77700000e+00   # Mtau(pole)
    21    4.75000000e-03   # Mdown(2 GeV) MSbar
    22    2.40000000e-03   # Mup(2 GeV) MSbar
    23    1.04000000e-01   # Mstrange(2 GeV) MSbar
    24    1.27000000e+00   # Mcharm(Mcharm) MSbar
    11    5.10998902e-04   # Me(pole)
    13    1.05658357e-01   # Mmu(pole)
Block VCKMIN               # input CKM mixing matrix parameters
     1    2.27200000e-01   # lambda_W
     2    8.18000000e-01   # A
     3    2.21000000e-01   # rhobar
     4    3.40000000e-01   # etabar (no phases used in SOFTSUSY yet though)
Block MINPAR               # SUSY breaking input parameters
     3    1.00000000e+01   # tanb, DRbar, Feynman gauge
     4    1.00000000e+00   # sign(mu)
     1    5.00000000e+01   # m0
     2    5.50000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.50850607e+16   # MX scale
Block RVLAMLQDIN           # input LLE couplings at MSUSY
  1 1 2   1.00000000e-03   # lambda'_{112}
# SOFTSUSY-specific non SLHA information:
# MIXING=1 Desired accuracy=1.00000000e-03 Achieved accuracy=2.73815457e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.05114093e+01   # MW
        25     1.15380156e+02   # CP even neutral scalar
        35     3.64986005e+02   # CP even neutral scalar
        36     3.64986005e+02   # CP odd neutral scalar
        37     2.09337919e+02   # charged scalar
   1000021     1.25106511e+03   # ~g
   1000022     2.26755488e+02   # ~neutralino(1)
   1000023     4.27900258e+02   # ~neutralino(2)
   1000024     4.27945396e+02   # ~chargino(1)
   1000025    -6.97918653e+02   # ~neutralino(3)
   1000035     7.10605534e+02   # ~neutralino(4)
   1000037     7.10892697e+02   # ~chargino(2)
   1000011     2.16984815e+02   # charged scalar
   1000013     2.17008594e+02   # charged scalar
   1000015     3.74564361e+02   # charged scalar
   2000011     3.74569844e+02   # charged scalar
   2000013     3.74677683e+02   # charged scalar
   2000015     7.75147054e+02   # charged scalar
   1000012     3.65957087e+02   # CP even neutral scalar
   1000014     3.65959173e+02   # CP even neutral scalar
   1000016     7.70743718e+02   # CP even neutral scalar
   1000017     3.65957087e+02   # CP odd neutral scalar
   1000018     3.65959173e+02   # CP odd neutral scalar
   1000019     7.70743718e+02   # CP odd neutral scalar
   1000001     1.05045729e+03   # ~d_1
   1000003     1.08771829e+03   # ~d_2
   1000005     1.09555270e+03   # ~d_3
   2000001     1.09555665e+03   # ~d_4
   2000003     1.14471281e+03   # ~d_5
   2000005     1.14473318e+03   # ~d_6
   1000002     8.80397780e+02   # ~u_1
   1000004     1.08263010e+03   # ~u_2
   1000006     1.09942761e+03   # ~u_3
   2000002     1.09943316e+03   # ~u_4
   2000004     1.14212047e+03   # ~u_5
   2000006     1.14212461e+03   # ~u_6
        12     0.00000000e+00   # Mnu1 inverted hierarchy output
        14     5.79701731e-23   # Mnu2 inverted hierarchy output
        16     0.00000000e+00   # Mnu3 inverted hierarchy output
Block RVNMIX  # neutrino-neutralino mixing matrix 
  1 1    0.00000000e+00   # N_{11}
  1 2    1.00000000e+00   # N_{12}
  1 3    0.00000000e+00   # N_{13}
  1 4    0.00000000e+00   # N_{14}
  1 5    0.00000000e+00   # N_{15}
  1 6    0.00000000e+00   # N_{16}
  1 7    0.00000000e+00   # N_{17}
  2 1    0.00000000e+00   # N_{21}
  2 2    0.00000000e+00   # N_{22}
  2 3    1.00000000e+00   # N_{23}
  2 4    0.00000000e+00   # N_{24}
  2 5    0.00000000e+00   # N_{25}
  2 6    0.00000000e+00   # N_{26}
  2 7    0.00000000e+00   # N_{27}
  3 1    1.00000000e+00   # N_{31}
  3 2    0.00000000e+00   # N_{32}
  3 3    0.00000000e+00   # N_{33}
  3 4   -1.52378561e-13   # N_{34}
  3 5    1.76647974e-13   # N_{35}
  3 6   -5.93762083e-13   # N_{36}
  3 7    2.69014505e-15   # N_{37}
  4 1    1.99406298e-13   # N_{41}
  4 2    0.00000000e+00   # N_{42}
  4 3    0.00000000e+00   # N_{43}
  4 4    9.96529381e-01   # N_{44}
  4 5   -1.47304411e-02   # N_{45}
  4 6    7.55678781e-02   # N_{46}
  4 7   -3.16496847e-02   # N_{47}
  5 1   -2.76601044e-13   # N_{51}
  5 2    0.00000000e+00   # N_{52}
  5 3    0.00000000e+00   # N_{53}
  5 4    3.22473212e-02   # N_{54}
  5 5    9.74698234e-01   # N_{55}
  5 6   -1.83582767e-01   # N_{56}
  5 7    1.23372728e-01   # N_{57}
  6 1    4.03769217e-13   # N_{61}
  6 2    0.00000000e+00   # N_{62}
  6 3    0.00000000e+00   # N_{63}
  6 4   -3.02509171e-02   # N_{64}
  6 5    4.39821532e-02   # N_{65}
  6 6    7.04075107e-01   # N_{66}
  6 7    7.08116301e-01   # N_{67}
  7 1   -3.57328687e-13   # N_{71}
  7 2    0.00000000e+00   # N_{72}
  7 3    0.00000000e+00   # N_{73}
  7 4    7.05279018e-02   # N_{74}
  7 5   -2.18659407e-01   # N_{75}
  7 6   -6.81810169e-01   # N_{76}
  7 7    6.94513334e-01   # N_{77}
Block RVUMIX  # lepton-chargino mixing matrix U
  1 1    1.00000000e+00   # U_{11}
  1 2    0.00000000e+00   # U_{12}
  1 3    0.00000000e+00   # U_{13}
  1 4    1.72218139e-13   # U_{14}
  1 5   -1.98293394e-14   # U_{15}
  2 1    0.00000000e+00   # U_{21}
  2 2    1.00000000e+00   # U_{22}
  2 3    0.00000000e+00   # U_{23}
  2 4    0.00000000e+00   # U_{24}
  2 5    0.00000000e+00   # U_{25}
  3 1    0.00000000e+00   # U_{31}
  3 2    0.00000000e+00   # U_{32}
  3 3    1.00000000e+00   # U_{33}
  3 4    0.00000000e+00   # U_{34}
  3 5    0.00000000e+00   # U_{35}
  4 1    1.71364140e-13   # U_{41}
  4 2    0.00000000e+00   # U_{42}
  4 3    0.00000000e+00   # U_{43}
  4 4   -9.64732328e-01   # U_{44}
  4 5    2.63232853e-01   # U_{45}
  5 1    2.62034673e-14   # U_{51}
  5 2    0.00000000e+00   # U_{52}
  5 3    0.00000000e+00   # U_{53}
  5 4   -2.63232853e-01   # U_{54}
  5 5   -9.64732328e-01   # U_{55}
Block RVVMIX  # lepton-chargino mixing matrix V
  1 1    1.00000000e+00   # V_{11}
  1 2    0.00000000e+00   # V_{12}
  1 3    0.00000000e+00   # V_{13}
  1 4    0.00000000e+00   # V_{14}
  1 5    9.69315873e-35   # V_{15}
  2 1    0.00000000e+00   # V_{21}
  2 2    1.00000000e+00   # V_{22}
  2 3    0.00000000e+00   # V_{23}
  2 4    0.00000000e+00   # V_{24}
  2 5    0.00000000e+00   # V_{25}
  3 1    0.00000000e+00   # V_{31}
  3 2    0.00000000e+00   # V_{32}
  3 3    1.00000000e+00   # V_{33}
  3 4    0.00000000e+00   # V_{34}
  3 5    0.00000000e+00   # V_{35}
  4 1   -1.68618867e-35   # V_{41}
  4 2    0.00000000e+00   # V_{42}
  4 3    0.00000000e+00   # V_{43}
  4 4   -9.84753325e-01   # V_{44}
  4 5    1.73956572e-01   # V_{45}
  5 1    9.54537029e-35   # V_{51}
  5 2    0.00000000e+00   # V_{52}
  5 3    0.00000000e+00   # V_{53}
  5 4   -1.73956572e-01   # V_{54}
  5 5   -9.84753325e-01   # V_{55}
Block gauge Q= 9.51752323e+02  # SM gauge couplings
     1     3.62936260e-01   # g'(Q)MSSM DRbar
     2     6.41937799e-01   # g(Q)MSSM DRbar
     3     1.05995676e+00   # g3(Q)MSSM DRbar
Block yu Q= 9.51752323e+02   # diagonal Up Yukawa matrix
  1  1     7.31143324e-06    # YU_{11}(Q)MSSM DRbar
  2  2     3.34917051e-03    # YU_{22}(Q)MSSM DRbar
  3  3     8.54686173e-01    # YU_{33}(Q)MSSM DRbar
Block yd Q= 9.51752323e+02   # diagonal down Yukawa matrix
  1  1     1.40710994e-04    # YD_{11}(Q)MSSM DRbar
  2  2     3.08089076e-03    # YD_{22}(Q)MSSM DRbar
  3  3     1.34260137e-01    # YD_{33}(Q)MSSM DRbar
Block ye Q= 9.51752323e+02   # diagonal lepton Yukawa matrix
  1  1     2.78741930e-05    # YE_{11}(Q)MSSM DRbar
  2  2     5.76350705e-03    # YE_{22}(Q)MSSM DRbar
  3  3     1.00227909e-01    # YE_{33}(Q)MSSM DRbar
Block hmix Q= 9.51752323e+02 # Higgs mixing parameters
     1     6.92578239e+02    # mu(Q)MSSM DRbar
     2     9.66570166e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.44438972e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     6.16418595e+05    # mA^2(Q)MSSM DRbar
Block RVLAMLLE Q= 9.51752323e+02 # non-zero R-Parity violating LLE couplings 
  1 2 2   -2.33515020e-21   # lambda_{122}
  1 3 3   -4.06106976e-20   # lambda_{133}
  2 1 2    2.33515020e-21   # lambda_{212}
  3 1 3    4.06106976e-20   # lambda_{313}
Block RVLAMLQD Q= 9.51752323e+02 # non-zero R-Parity violating LQD couplings 
  1 1 1   -1.17162894e-22   # lambda'_{111}
  1 1 2    9.99999882e-04   # lambda'_{112}
  1 1 3   -9.85512692e-18   # lambda'_{113}
  1 2 1    1.01467611e-28   # lambda'_{121}
  1 2 2   -2.79882545e-15   # lambda'_{122}
  1 2 3    1.67140779e-23   # lambda'_{123}
  1 3 1   -2.39466401e-27   # lambda'_{131}
  1 3 2    6.56345728e-14   # lambda'_{132}
  1 3 3   -5.37998956e-20   # lambda'_{133}
Block RVLAMUDD Q= 9.51752323e+02 # non-zero R-Parity violating UDD couplings 
Block RVTLLE Q= 9.51752323e+02 # non-zero R-Parity violating LLE soft terms 
  1 2 2   -3.16335467e-17   # T_{122}
  1 3 3   -5.50145248e-16   # T_{133}
  2 1 2    3.16335467e-17   # T_{212}
  3 1 3    5.50145248e-16   # T_{313}
Block RVTLQD Q= 9.51752323e+02 # non-zero R-Parity violating LQD soft terms 
  1 1 1   -1.13579886e-18   # T'_{111}
  1 1 2    2.52837758e-07   # T'_{112}
  1 1 3   -6.92597107e-14   # T'_{113}
  1 2 1    3.98376611e-24   # T'_{121}
  1 2 2    5.50971303e-12   # T'_{122}
  1 2 3    6.45434083e-19   # T'_{123}
  1 3 1   -9.32090582e-23   # T'_{131}
  1 3 2   -1.29183109e-10   # T'_{132}
  1 3 3   -6.88049708e-16   # T'_{133}
Block RVTUDD Q= 9.51752323e+02 # non-zero R-Parity violating UDD soft terms 
Block RVKAPPA Q= 9.51752323e+02 # R-Parity violating kappa 
     1   -2.85003144e-16   # kappa_{1}
     2    0.00000000e+00   # kappa_{2}
     3    0.00000000e+00   # kappa_{3}
Block RVD Q= 9.51752323e+02 # R-Parity violating D 
     1   -3.84430335e-12   # D_{1}
     2    0.00000000e+00   # D_{2}
     3    0.00000000e+00   # D_{3}
Block RVSNVEV Q= 9.51752323e+02 # sneutrino VEVs D 
     1   -1.60002565e-10   # SneutrinoVev_{1}
     2    0.00000000e+00   # SneutrinoVev_{2}
     3    0.00000000e+00   # SneutrinoVev_{3}
Block RVM2LH1 Q= 9.51752323e+02 # M2LH1 
     1    1.19353999e-11   # M2LH1_{1}
     2    0.00000000e+00   # M2LH1_{2}
     3    0.00000000e+00   # M2LH1_{3}
Block UPMNS Q= 9.11876000e+01 # neutrino mixing matrix (inverted  hierarchy)
  1  1     0.00000000e+00   # UPMNS_{11} matrix element
  1  2     1.00000000e+00   # UPMNS_{12} matrix element
  1  3     0.00000000e+00   # UPMNS_{13} matrix element
  2  1     0.00000000e+00   # UPMNS_{21} matrix element
  2  2     0.00000000e+00   # UPMNS_{22} matrix element
  2  3     1.00000000e+00   # UPMNS_{23} matrix element
  3  1     1.00000000e+00   # UPMNS_{31} matrix element
  3  2     0.00000000e+00   # UPMNS_{32} matrix element
  3  3     0.00000000e+00   # UPMNS_{33} matrix element
Block msq2 Q= 9.51752323e+02 # super CKM squark mass^2 matrix - DRbar
  1  1     1.22019611e+06    # (m^_Q^2)_{11}
  1  2     4.11585017e+01    # (m^_Q^2)_{12}
  1  3    -9.73233166e+02    # (m^_Q^2)_{13}
  2  1     4.11585017e+01    # (m^_Q^2)_{21}
  2  2     1.21989812e+06    # (m^_Q^2)_{22}
  2  3     7.16513513e+03    # (m^_Q^2)_{23}
  3  1    -9.73233166e+02    # (m^_Q^2)_{31}
  3  2     7.16513513e+03    # (m^_Q^2)_{32}
  3  3     1.04214523e+06    # (m^_Q^2)_{33}
Block msl2 Q= 9.51752323e+02 # super MNS slepton mass^2 matrix - DRbar
  1  1     1.34222647e+05    # (m^_L^2)_{11}
  1  2     0.00000000e+00    # (m^_L^2)_{12}
  1  3     0.00000000e+00    # (m^_L^2)_{13}
  2  1     0.00000000e+00    # (m^_L^2)_{21}
  2  2     1.34221148e+05    # (m^_L^2)_{22}
  2  3     0.00000000e+00    # (m^_L^2)_{23}
  3  1     0.00000000e+00    # (m^_L^2)_{31}
  3  2     0.00000000e+00    # (m^_L^2)_{32}
  3  3     1.33519679e+05    # (m^_L^2)_{33}
Block msd2 Q= 9.51752323e+02 # super CKM squark mass^2 matrix - DRbar
  1  1     1.11827922e+06    # (m^_d^2)_{11}
  1  2    -3.12939351e-06    # (m^_d^2)_{12}
  1  3     3.09950842e-03    # (m^_d^2)_{13}
  2  1    -3.12939351e-06    # (m^_d^2)_{21}
  2  2     1.11827254e+06    # (m^_d^2)_{22}
  2  3    -4.99618889e-01    # (m^_d^2)_{23}
  3  1     3.09950842e-03    # (m^_d^2)_{31}
  3  2    -4.99618889e-01    # (m^_d^2)_{32}
  3  3     1.10715766e+06    # (m^_d^2)_{33}
Block msu2 Q= 9.51752323e+02 # super CKM squark mass^2 matrix - DRbar
  1  1     1.12916299e+06    # (m^_u^2)_{11}
  1  2    -1.48354502e-08    # (m^_u^2)_{12}
  1  3    -1.36941680e-05    # (m^_u^2)_{13}
  2  1    -1.48354502e-08    # (m^_u^2)_{21}
  2  2     1.12915671e+06    # (m^_u^2)_{22}
  2  3    -6.66474553e-02    # (m^_u^2)_{23}
  3  1    -1.36941680e-05    # (m^_u^2)_{31}
  3  2    -6.66474553e-02    # (m^_u^2)_{32}
  3  3     7.75196253e+05    # (m^_u^2)_{33}
Block mse2 Q= 9.51752323e+02 # super MNS slepton mass^2 matrix - DRbar
  1  1     4.31434578e+04    # (m^_e^2)_{11}
  1  2     0.00000000e+00    # (m^_e^2)_{12}
  1  3     0.00000000e+00    # (m^_e^2)_{13}
  2  1     0.00000000e+00    # (m^_e^2)_{21}
  2  2     4.31386917e+04    # (m^_e^2)_{22}
  2  3     0.00000000e+00    # (m^_e^2)_{23}
  3  1     0.00000000e+00    # (m^_e^2)_{31}
  3  2     0.00000000e+00    # (m^_e^2)_{32}
  3  3     4.17042781e+04    # (m^_e^2)_{33}
Block tu Q= 9.51752323e+02   # super CKM trilinear matrix - DRbar
  1  1    -9.10660665e-03    # (T^_u)_{11}
  1  2    -1.92239475e-08    # (T^_u)_{12}
  1  3    -8.62645507e-08    # (T^_u)_{13}
  2  1    -8.80598447e-06    # (T^_u)_{21}
  2  2    -4.17152553e+00    # (T^_u)_{22}
  2  3    -4.19471345e-04    # (T^_u)_{23}
  3  1    -1.19021262e-02    # (T^_u)_{31}
  3  2    -1.26319541e-01    # (T^_u)_{32}
  3  3    -8.26621280e+02    # (T^_u)_{33}
Block td Q= 9.51752323e+02   # super CKM trilinear matrix - DRbar
  1  1    -2.13517915e-01    # (T^_d)_{11}
  1  2    -3.13928124e-06    # (T^_d)_{12}
  1  3     7.46294062e-05    # (T^_d)_{13}
  2  1    -6.87348001e-05    # (T^_d)_{21}
  2  2    -4.67452401e+00    # (T^_d)_{22}
  2  3    -1.20299763e-02    # (T^_d)_{23}
  3  1     7.10718387e-02    # (T^_d)_{31}
  3  2    -5.23245821e-01    # (T^_d)_{32}
  3  3    -1.90629287e+02    # (T^_d)_{33}
Block te Q= 9.51752323e+02   # super CKM trilinear matrix - DRbar
  1  1    -9.07847016e-03    # (T^_e)_{11}
  1  2     0.00000000e+00    # (T^_e)_{12}
  1  3     0.00000000e+00    # (T^_e)_{13}
  2  1     0.00000000e+00    # (T^_e)_{21}
  2  2    -1.87710714e+00    # (T^_e)_{22}
  2  3     0.00000000e+00    # (T^_e)_{23}
  3  1     0.00000000e+00    # (T^_e)_{31}
  3  2     0.00000000e+00    # (T^_e)_{32}
  3  3    -3.24678525e+01    # (T^_e)_{33}
Block VCKM Q= 9.51752323e+02 # DRbar CKM mixing matrix
  1  1     9.73840728e-01    # CKM_{11} matrix element
  1  2     2.27197389e-01    # CKM_{12} matrix element
  1  3     3.94754488e-03    # CKM_{13} matrix element
  2  1    -2.27161577e-01    # CKM_{21} matrix element
  2  2     9.72961872e-01    # CKM_{22} matrix element
  2  3     4.17470254e-02    # CKM_{23} matrix element
  3  1     5.64400451e-03    # CKM_{31} matrix element
  3  2    -4.15516841e-02    # CKM_{32} matrix element
  3  3     9.99120415e-01    # CKM_{33} matrix element
Block msoft Q= 9.51752323e+02 # MSSM DRbar SUSY breaking parameters
     1     2.32454746e+02     # M_1(Q)
     2     4.28959088e+02     # M_2(Q)
     3     1.21684310e+03     # M_3(Q)
    21     1.15648575e+05     # mH1^2(Q)
    22    -4.49251326e+05     # mH2^2(Q)
Block USQMIX  # super CKM squark mass^2 matrix
  1  1     2.12118866e-05   # (USQMIX)_{11}
  1  2     2.24598782e-04   # (USQMIX)_{12}
  1  3     4.12240027e-01   # (USQMIX)_{13}
  1  4     1.23860495e-10   # (USQMIX)_{14}
  1  5     6.01409193e-07   # (USQMIX)_{15}
  1  6     9.11075249e-01   # (USQMIX)_{16}
  2  1     1.49995173e-04   # (USQMIX)_{21}
  2  2     1.58744722e-03   # (USQMIX)_{22}
  2  3     9.11073965e-01   # (USQMIX)_{23}
  2  4     6.98821930e-09   # (USQMIX)_{24}
  2  5     3.38862147e-05   # (USQMIX)_{25}
  2  6    -4.12239841e-01   # (USQMIX)_{26}
  3  1     1.34715137e-07   # (USQMIX)_{31}
  3  2     7.92565601e-03   # (USQMIX)_{32}
  3  3    -4.33164882e-05   # (USQMIX)_{33}
  3  4     2.16313963e-08   # (USQMIX)_{34}
  3  5     9.99968590e-01   # (USQMIX)_{35}
  3  6     1.69857588e-05   # (USQMIX)_{36}
  4  1     1.73022378e-05   # (USQMIX)_{41}
  4  2     1.22676057e-10   # (USQMIX)_{42}
  4  3    -8.93314834e-09   # (USQMIX)_{43}
  4  4     1.00000000e+00   # (USQMIX)_{44}
  4  5    -2.16358254e-08   # (USQMIX)_{45}
  4  6     3.50323771e-09   # (USQMIX)_{46}
  5  1     1.63078723e-01   # (USQMIX)_{51}
  5  2     9.86580764e-01   # (USQMIX)_{52}
  5  3    -1.54169245e-03   # (USQMIX)_{53}
  5  4    -2.82193241e-06   # (USQMIX)_{54}
  5  5    -7.81964177e-03   # (USQMIX)_{55}
  5  6     4.50575289e-04   # (USQMIX)_{56}
  6  1     9.86613048e-01   # (USQMIX)_{61}
  6  2    -1.63073634e-01   # (USQMIX)_{62}
  6  3     1.07454646e-04   # (USQMIX)_{63}
  6  4    -1.70705646e-05   # (USQMIX)_{64}
  6  5     1.29237839e-03   # (USQMIX)_{65}
  6  6    -3.13911150e-05   # (USQMIX)_{66}
Block DSQMIX  # super CKM squark mass^2 matrix
  1  1     4.58910896e-03   # (DSQMIX)_{11}
  1  2    -3.37863748e-02   # (DSQMIX)_{12}
  1  3     9.74684519e-01   # (DSQMIX)_{13}
  1  4     9.54065843e-07   # (DSQMIX)_{14}
  1  5    -1.53805180e-04   # (DSQMIX)_{15}
  1  6     2.20969421e-01   # (DSQMIX)_{16}
  2  1    -1.70102345e-03   # (DSQMIX)_{21}
  2  2     1.25245925e-02   # (DSQMIX)_{22}
  2  3    -2.20659687e-01   # (DSQMIX)_{23}
  2  4    -2.20237648e-06   # (DSQMIX)_{24}
  2  5     3.55191883e-04   # (DSQMIX)_{25}
  2  6     9.75268895e-01   # (DSQMIX)_{26}
  3  1     1.79991756e-06   # (DSQMIX)_{31}
  3  2     4.07109163e-03   # (DSQMIX)_{32}
  3  3     3.74137454e-04   # (DSQMIX)_{33}
  3  4     4.71006867e-06   # (DSQMIX)_{34}
  3  5     9.99991588e-01   # (DSQMIX)_{35}
  3  6    -3.31823905e-04   # (DSQMIX)_{36}
  4  1     1.86544395e-04   # (DSQMIX)_{41}
  4  2     6.29777134e-08   # (DSQMIX)_{42}
  4  3    -2.32250950e-06   # (DSQMIX)_{43}
  4  4     9.99999983e-01   # (DSQMIX)_{44}
  4  5    -4.70914818e-06   # (DSQMIX)_{45}
  4  6     2.05901419e-06   # (DSQMIX)_{46}
  5  1    -1.37739022e-01   # (DSQMIX)_{51}
  5  2     9.89792595e-01   # (DSQMIX)_{52}
  5  3     3.60449183e-02   # (DSQMIX)_{53}
  5  4     2.57066491e-05   # (DSQMIX)_{54}
  5  5    -4.04439927e-03   # (DSQMIX)_{55}
  5  6    -4.79452312e-03   # (DSQMIX)_{56}
  6  1     9.90456447e-01   # (DSQMIX)_{61}
  6  2     1.37824748e-01   # (DSQMIX)_{62}
  6  3     1.17633432e-04   # (DSQMIX)_{63}
  6  4    -1.84775127e-04   # (DSQMIX)_{64}
  6  5    -5.62932986e-04   # (DSQMIX)_{65}
  6  6    -1.56394547e-05   # (DSQMIX)_{66}
Block RVLMIX  # charged higgs-slepton mixing matrix 
  1 1    0.00000000e+00   # C_{11}
  1 2    0.00000000e+00   # C_{12}
  1 3    0.00000000e+00   # C_{13}
  1 4    0.00000000e+00   # C_{14}
  1 5    1.32480073e-01   # C_{15}
  1 6    0.00000000e+00   # C_{16}
  1 7    0.00000000e+00   # C_{17}
  1 8    9.91185669e-01   # C_{18}
  2 1    0.00000000e+00   # C_{21}
  2 2    0.00000000e+00   # C_{22}
  2 3    0.00000000e+00   # C_{23}
  2 4    7.88980026e-03   # C_{24}
  2 5    0.00000000e+00   # C_{25}
  2 6    0.00000000e+00   # C_{26}
  2 7    9.99968875e-01   # C_{27}
  2 8    0.00000000e+00   # C_{28}
  3 1   -4.60678352e-18   # C_{31}
  3 2    5.84664998e-17   # C_{32}
  3 3    3.81626013e-05   # C_{33}
  3 4    0.00000000e+00   # C_{34}
  3 5    0.00000000e+00   # C_{35}
  3 6    9.99999999e-01   # C_{36}
  3 7    0.00000000e+00   # C_{37}
  3 8    0.00000000e+00   # C_{38}
  4 1    1.61974673e-15   # C_{41}
  4 2    8.26965448e-16   # C_{42}
  4 3    9.99999999e-01   # C_{43}
  4 4    0.00000000e+00   # C_{44}
  4 5    0.00000000e+00   # C_{45}
  4 6   -3.81626013e-05   # C_{46}
  4 7    0.00000000e+00   # C_{47}
  4 8    0.00000000e+00   # C_{48}
  5 1    0.00000000e+00   # C_{51}
  5 2    0.00000000e+00   # C_{52}
  5 3    0.00000000e+00   # C_{53}
  5 4    9.99968875e-01   # C_{54}
  5 5    0.00000000e+00   # C_{55}
  5 6    0.00000000e+00   # C_{56}
  5 7   -7.88980026e-03   # C_{57}
  5 8    0.00000000e+00   # C_{58}
  6 1    0.00000000e+00   # C_{61}
  6 2    0.00000000e+00   # C_{62}
  6 3    0.00000000e+00   # C_{63}
  6 4    0.00000000e+00   # C_{64}
  6 5    9.91185669e-01   # C_{65}
  6 6    0.00000000e+00   # C_{66}
  6 7    0.00000000e+00   # C_{67}
  6 8   -1.32480073e-01   # C_{68}
  7 1    9.94690746e-01   # C_{71}
  7 2    1.02909279e-01   # C_{72}
  7 3   -1.69624955e-15   # C_{73}
  7 4    0.00000000e+00   # C_{74}
  7 5    0.00000000e+00   # C_{75}
  7 6   -1.36968709e-18   # C_{76}
  7 7    0.00000000e+00   # C_{77}
  7 8    0.00000000e+00   # C_{78}
Block RVHMIX  # CP-even neutral scalar mixing matrix 
  1 1    4.92740646e-01   # curlyN_{11}
  1 2    8.70176221e-01   # curlyN_{12}
  1 3    2.45324816e-11   # curlyN_{13}
  1 4    0.00000000e+00   # curlyN_{14}
  1 5    0.00000000e+00   # curlyN_{15}
  2 1    0.00000000e+00   # curlyN_{21}
  2 2    0.00000000e+00   # curlyN_{22}
  2 3    0.00000000e+00   # curlyN_{23}
  2 4    0.00000000e+00   # curlyN_{24}
  2 5    1.00000000e+00   # curlyN_{25}
  3 1    0.00000000e+00   # curlyN_{31}
  3 2    0.00000000e+00   # curlyN_{32}
  3 3    0.00000000e+00   # curlyN_{33}
  3 4    1.00000000e+00   # curlyN_{34}
  3 5    0.00000000e+00   # curlyN_{35}
  4 1   -4.32796109e-11   # curlyN_{41}
  4 2   -3.68529736e-12   # curlyN_{42}
  4 3    1.00000000e+00   # curlyN_{43}
  4 4    0.00000000e+00   # curlyN_{44}
  4 5    0.00000000e+00   # curlyN_{45}
  5 1    8.70176221e-01   # curlyN_{51}
  5 2   -4.92740646e-01   # curlyN_{52}
  5 3    3.58449925e-11   # curlyN_{53}
  5 4    0.00000000e+00   # curlyN_{54}
  5 5    0.00000000e+00   # curlyN_{55}
Block RVAMIX  # CP-odd neutral scalar mixing matrix 
  1 1    0.00000000e+00   # curlyN~_{11}
  1 2    0.00000000e+00   # curlyN~_{12}
  1 3    0.00000000e+00   # curlyN~_{13}
  1 4    0.00000000e+00   # curlyN~_{14}
  1 5    1.00000000e+00   # curlyN~_{15}
  2 1    0.00000000e+00   # curlyN~_{21}
  2 2    0.00000000e+00   # curlyN~_{22}
  2 3    0.00000000e+00   # curlyN~_{23}
  2 4    1.00000000e+00   # curlyN~_{24}
  2 5    0.00000000e+00   # curlyN~_{25}
  3 1   -4.39765256e-11   # curlyN~_{31}
  3 2    3.80810553e-12   # curlyN~_{32}
  3 3    1.00000000e+00   # curlyN~_{33}
  3 4    0.00000000e+00   # curlyN~_{34}
  3 5    0.00000000e+00   # curlyN~_{35}
  4 1    8.71992550e-01   # curlyN~_{41}
  4 2    4.89519144e-01   # curlyN~_{42}
  4 3    3.64830621e-11   # curlyN~_{43}
  4 4    0.00000000e+00   # curlyN~_{44}
  4 5    0.00000000e+00   # curlyN~_{45}
