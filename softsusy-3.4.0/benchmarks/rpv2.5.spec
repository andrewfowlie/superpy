# SuperPy: Jacobian for naturalness priors.
# J = 5.07371281e+05
# b = 1.22544186e+05
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
     2    6.00000000e+02   # m12
     5    0.00000000e+00   # A0
Block EXTPAR               # scale of SUSY breaking BCs
     0    1.46222145e+16   # MX scale
Block RVLAMLQDIN           # input LLE couplings at MSUSY
  1 1 2   1.00000000e-03   # lambda'_{112}
# SOFTSUSY-specific non SLHA information:
# MIXING=1 Desired accuracy=1.00000000e-03 Achieved accuracy=3.57798181e-04
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.05220000e+01   # MW
        25     1.15932660e+02   # CP even neutral scalar
        35     3.97929570e+02   # CP even neutral scalar
        36     3.97929570e+02   # CP odd neutral scalar
        37     2.27191106e+02   # charged scalar
   1000021     1.35583764e+03   # ~g
   1000022     2.48582033e+02   # ~neutralino(1)
   1000023     4.69257887e+02   # ~neutralino(2)
   1000024     4.69324413e+02   # ~chargino(1)
   1000025    -7.53671881e+02   # ~neutralino(3)
   1000035     7.65770057e+02   # ~neutralino(4)
   1000037     7.66023368e+02   # ~chargino(2)
   1000011     2.34806675e+02   # charged scalar
   1000013     2.34830214e+02   # charged scalar
   1000015     4.06599165e+02   # charged scalar
   2000011     4.06905719e+02   # charged scalar
   2000013     4.06910336e+02   # charged scalar
   2000015     8.38378793e+02   # charged scalar
   1000012     3.98974613e+02   # CP even neutral scalar
   1000014     3.98976852e+02   # CP even neutral scalar
   1000016     8.34298292e+02   # CP even neutral scalar
   1000017     3.98974613e+02   # CP odd neutral scalar
   1000018     3.98976852e+02   # CP odd neutral scalar
   1000019     8.34298292e+02   # CP odd neutral scalar
   1000001     1.13781765e+03   # ~d_1
   1000003     1.17687186e+03   # ~d_2
   1000005     1.18546075e+03   # ~d_3
   2000001     1.18546487e+03   # ~d_4
   2000003     1.23923629e+03   # ~d_5
   2000005     1.23925858e+03   # ~d_6
   1000002     9.56915942e+02   # ~u_1
   1000004     1.16648297e+03   # ~u_2
   1000006     1.18988770e+03   # ~u_3
   2000002     1.18989322e+03   # ~u_4
   2000004     1.23685612e+03   # ~u_5
   2000006     1.23686105e+03   # ~u_6
        12     0.00000000e+00   # Mnu1 inverted hierarchy output
        14     5.40002589e-23   # Mnu2 inverted hierarchy output
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
  3 4   -1.38250363e-13   # N_{34}
  3 5    1.60518071e-13   # N_{35}
  3 6   -5.86148318e-13   # N_{36}
  3 7    2.25773047e-15   # N_{37}
  4 1    1.80971790e-13   # N_{41}
  4 2    0.00000000e+00   # N_{42}
  4 3    0.00000000e+00   # N_{43}
  4 4    9.97025615e-01   # N_{44}
  4 5   -1.26302783e-02   # N_{45}
  4 6    7.00135608e-02   # N_{46}
  4 7   -2.96395197e-02   # N_{47}
  5 1   -2.55154973e-13   # N_{51}
  5 2    0.00000000e+00   # N_{52}
  5 3    0.00000000e+00   # N_{53}
  5 4    2.81031681e-02   # N_{54}
  5 5    9.77243906e-01   # N_{55}
  5 6   -1.73860466e-01   # N_{56}
  5 7    1.18224778e-01   # N_{57}
  6 1    4.00999408e-13   # N_{61}
  6 2    0.00000000e+00   # N_{62}
  6 3    0.00000000e+00   # N_{63}
  6 4   -2.79140717e-02   # N_{64}
  6 5    4.04956598e-02   # N_{65}
  6 6    7.04526855e-01   # N_{66}
  6 7    7.07970915e-01   # N_{67}
  7 1   -3.60275988e-13   # N_{71}
  7 2    0.00000000e+00   # N_{72}
  7 3    0.00000000e+00   # N_{73}
  7 4    6.61130872e-02   # N_{74}
  7 5   -2.07833889e-01   # N_{75}
  7 6   -6.84479766e-01   # N_{76}
  7 7    6.95644726e-01   # N_{77}
Block RVUMIX  # lepton-chargino mixing matrix U
  1 1    1.00000000e+00   # U_{11}
  1 2    0.00000000e+00   # U_{12}
  1 3    0.00000000e+00   # U_{13}
  1 4    1.56283148e-13   # U_{14}
  1 5   -1.44908680e-14   # U_{15}
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
  4 1    1.54972146e-13   # U_{41}
  4 2    0.00000000e+00   # U_{42}
  4 3    0.00000000e+00   # U_{43}
  4 4   -9.68534978e-01   # U_{44}
  4 5    2.48877471e-01   # U_{45}
  5 1    2.48604422e-14   # U_{51}
  5 2    0.00000000e+00   # U_{52}
  5 3    0.00000000e+00   # U_{53}
  5 4   -2.48877471e-01   # U_{54}
  5 5   -9.68534978e-01   # U_{55}
Block RVVMIX  # lepton-chargino mixing matrix V
  1 1    1.00000000e+00   # V_{11}
  1 2    0.00000000e+00   # V_{12}
  1 3    0.00000000e+00   # V_{13}
  1 4    0.00000000e+00   # V_{14}
  1 5    9.69149783e-35   # V_{15}
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
  4 1   -1.61240830e-35   # V_{41}
  4 2    0.00000000e+00   # V_{42}
  4 3    0.00000000e+00   # V_{43}
  4 4   -9.86062808e-01   # V_{44}
  4 5    1.66373489e-01   # V_{45}
  5 1    9.55642556e-35   # V_{51}
  5 2    0.00000000e+00   # V_{52}
  5 3    0.00000000e+00   # V_{53}
  5 4   -1.66373489e-01   # V_{54}
  5 5   -9.86062808e-01   # V_{55}
Block gauge Q= 1.03033242e+03  # SM gauge couplings
     1     3.63109772e-01   # g'(Q)MSSM DRbar
     2     6.41501676e-01   # g(Q)MSSM DRbar
     3     1.05584083e+00   # g3(Q)MSSM DRbar
Block yu Q= 1.03033242e+03   # diagonal Up Yukawa matrix
  1  1     7.29208536e-06    # YU_{11}(Q)MSSM DRbar
  2  2     3.34030789e-03    # YU_{22}(Q)MSSM DRbar
  3  3     8.51789055e-01    # YU_{33}(Q)MSSM DRbar
Block yd Q= 1.03033242e+03   # diagonal down Yukawa matrix
  1  1     1.40211742e-04    # YD_{11}(Q)MSSM DRbar
  2  2     3.06996124e-03    # YD_{22}(Q)MSSM DRbar
  3  3     1.33701830e-01    # YD_{33}(Q)MSSM DRbar
Block ye Q= 1.03033242e+03   # diagonal lepton Yukawa matrix
  1  1     2.78419993e-05    # YE_{11}(Q)MSSM DRbar
  2  2     5.75685067e-03    # YE_{22}(Q)MSSM DRbar
  3  3     1.00158550e-01    # YE_{33}(Q)MSSM DRbar
Block hmix Q= 1.03033242e+03 # Higgs mixing parameters
     1     7.48344811e+02    # mu(Q)MSSM DRbar
     2     9.65652312e+00    # tan beta(Q)MSSM DRbar Feynman gauge
     3     2.44381404e+02    # higgs vev(Q)MSSM DRbar Feynman gauge
     4     7.21969961e+05    # mA^2(Q)MSSM DRbar
Block RVLAMLLE Q= 1.03033242e+03 # non-zero R-Parity violating LLE couplings 
  1 2 2   -2.68154920e-21   # lambda_{122}
  1 3 3   -4.66566444e-20   # lambda_{133}
  2 1 2    2.68154920e-21   # lambda_{212}
  3 1 3    4.66566444e-20   # lambda_{313}
Block RVLAMLQD Q= 1.03033242e+03 # non-zero R-Parity violating LQD couplings 
  1 1 1   -1.34822869e-22   # lambda'_{111}
  1 1 2    9.99999869e-04   # lambda'_{112}
  1 1 3   -1.13799395e-17   # lambda'_{113}
  1 2 1    1.19414260e-28   # lambda'_{121}
  1 2 2   -3.11895523e-15   # lambda'_{122}
  1 2 3    1.96523139e-23   # lambda'_{123}
  1 3 1   -2.81818094e-27   # lambda'_{131}
  1 3 2    7.31415181e-14   # lambda'_{132}
  1 3 3   -6.15719257e-20   # lambda'_{133}
Block RVLAMUDD Q= 1.03033242e+03 # non-zero R-Parity violating UDD couplings 
Block RVTLLE Q= 1.03033242e+03 # non-zero R-Parity violating LLE soft terms 
  1 2 2   -3.76775385e-17   # T_{122}
  1 3 3   -6.55561178e-16   # T_{133}
  2 1 2    3.76775385e-17   # T_{212}
  3 1 3    6.55561178e-16   # T_{313}
Block RVTLQD Q= 1.03033242e+03 # non-zero R-Parity violating LQD soft terms 
  1 1 1   -1.34300435e-18   # T'_{111}
  1 1 2    3.05481340e-07   # T'_{112}
  1 1 3   -8.15499474e-14   # T'_{113}
  1 2 1    4.91727655e-24   # T'_{121}
  1 2 2    6.64505460e-12   # T'_{122}
  1 2 3    7.95763302e-19   # T'_{123}
  1 3 1   -1.15037222e-22   # T'_{131}
  1 3 2   -1.55797769e-10   # T'_{132}
  1 3 3   -8.15270223e-16   # T'_{133}
Block RVTUDD Q= 1.03033242e+03 # non-zero R-Parity violating UDD soft terms 
Block RVKAPPA Q= 1.03033242e+03 # R-Parity violating kappa 
     1   -3.54210214e-16   # kappa_{1}
     2    0.00000000e+00   # kappa_{2}
     3    0.00000000e+00   # kappa_{3}
Block RVD Q= 1.03033242e+03 # R-Parity violating D 
     1   -4.95290629e-12   # D_{1}
     2    0.00000000e+00   # D_{2}
     3    0.00000000e+00   # D_{3}
Block RVSNVEV Q= 1.03033242e+03 # sneutrino VEVs D 
     1   -1.58854635e-10   # SneutrinoVev_{1}
     2    0.00000000e+00   # SneutrinoVev_{2}
     3    0.00000000e+00   # SneutrinoVev_{3}
Block RVM2LH1 Q= 1.03033242e+03 # M2LH1 
     1    1.56323882e-11   # M2LH1_{1}
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
Block msq2 Q= 1.03033242e+03 # super CKM squark mass^2 matrix - DRbar
  1  1     1.43111176e+06    # (m^_Q^2)_{11}
  1  2     4.80677411e+01    # (m^_Q^2)_{12}
  1  3    -1.13695727e+03    # (m^_Q^2)_{13}
  2  1     4.80677411e+01    # (m^_Q^2)_{21}
  2  2     1.43076375e+06    # (m^_Q^2)_{22}
  2  3     8.37050596e+03    # (m^_Q^2)_{23}
  3  1    -1.13695727e+03    # (m^_Q^2)_{31}
  3  2     8.37050596e+03    # (m^_Q^2)_{32}
  3  3     1.22307006e+06    # (m^_Q^2)_{33}
Block msl2 Q= 1.03033242e+03 # super MNS slepton mass^2 matrix - DRbar
  1  1     1.58687831e+05    # (m^_L^2)_{11}
  1  2     0.00000000e+00    # (m^_L^2)_{12}
  1  3     0.00000000e+00    # (m^_L^2)_{13}
  2  1     0.00000000e+00    # (m^_L^2)_{21}
  2  2     1.58686074e+05    # (m^_L^2)_{22}
  2  3     0.00000000e+00    # (m^_L^2)_{23}
  3  1     0.00000000e+00    # (m^_L^2)_{31}
  3  2     0.00000000e+00    # (m^_L^2)_{32}
  3  3     1.57862468e+05    # (m^_L^2)_{33}
Block msd2 Q= 1.03033242e+03 # super CKM squark mass^2 matrix - DRbar
  1  1     1.31010361e+06    # (m^_d^2)_{11}
  1  2    -3.62429023e-06    # (m^_d^2)_{12}
  1  3     3.58155820e-03    # (m^_d^2)_{13}
  2  1    -3.62429023e-06    # (m^_d^2)_{21}
  2  2     1.31009582e+06    # (m^_d^2)_{22}
  2  3    -5.77321698e-01    # (m^_d^2)_{23}
  3  1     3.58155820e-03    # (m^_d^2)_{31}
  3  2    -5.77321698e-01    # (m^_d^2)_{32}
  3  3     1.29715706e+06    # (m^_d^2)_{33}
Block msu2 Q= 1.03033242e+03 # super CKM squark mass^2 matrix - DRbar
  1  1     1.32307270e+06    # (m^_u^2)_{11}
  1  2    -1.75187438e-08    # (m^_u^2)_{12}
  1  3    -1.57316372e-05    # (m^_u^2)_{13}
  2  1    -1.75187438e-08    # (m^_u^2)_{21}
  2  2     1.32306536e+06    # (m^_u^2)_{22}
  2  3    -7.65648369e-02    # (m^_u^2)_{23}
  3  1    -1.57316372e-05    # (m^_u^2)_{31}
  3  2    -7.65648369e-02    # (m^_u^2)_{32}
  3  3     9.09448048e+05    # (m^_u^2)_{33}
Block mse2 Q= 1.03033242e+03 # super MNS slepton mass^2 matrix - DRbar
  1  1     5.08063264e+04    # (m^_e^2)_{11}
  1  2     0.00000000e+00    # (m^_e^2)_{12}
  1  3     0.00000000e+00    # (m^_e^2)_{13}
  2  1     0.00000000e+00    # (m^_e^2)_{21}
  2  2     5.08007350e+04    # (m^_e^2)_{22}
  2  3     0.00000000e+00    # (m^_e^2)_{23}
  3  1     0.00000000e+00    # (m^_e^2)_{31}
  3  2     0.00000000e+00    # (m^_e^2)_{32}
  3  3     4.91163895e+04    # (m^_e^2)_{33}
Block tu Q= 1.03033242e+03   # super CKM trilinear matrix - DRbar
  1  1    -9.81938960e-03    # (T^_u)_{11}
  1  2    -2.04884561e-08    # (T^_u)_{12}
  1  3    -9.19176654e-08    # (T^_u)_{13}
  2  1    -9.38522663e-06    # (T^_u)_{21}
  2  2    -4.49803448e+00    # (T^_u)_{22}
  2  3    -4.46961970e-04    # (T^_u)_{23}
  3  1    -1.26643621e-02    # (T^_u)_{31}
  3  2    -1.34409895e-01    # (T^_u)_{32}
  3  3    -8.91976823e+02    # (T^_u)_{33}
Block td Q= 1.03033242e+03   # super CKM trilinear matrix - DRbar
  1  1    -2.29766809e-01    # (T^_d)_{11}
  1  2    -3.36193903e-06    # (T^_d)_{12}
  1  3     7.99440683e-05    # (T^_d)_{13}
  2  1    -7.36099528e-05    # (T^_d)_{21}
  2  2    -5.03026443e+00    # (T^_d)_{22}
  2  3    -1.28866893e-02    # (T^_d)_{23}
  3  1     7.60882203e-02    # (T^_d)_{31}
  3  2    -5.60177568e-01    # (T^_d)_{32}
  3  3    -2.05069033e+02    # (T^_d)_{33}
Block te Q= 1.03033242e+03   # super CKM trilinear matrix - DRbar
  1  1    -9.85261530e-03    # (T^_e)_{11}
  1  2     0.00000000e+00    # (T^_e)_{12}
  1  3     0.00000000e+00    # (T^_e)_{13}
  2  1     0.00000000e+00    # (T^_e)_{21}
  2  2    -2.03717333e+00    # (T^_e)_{22}
  2  3     0.00000000e+00    # (T^_e)_{23}
  3  1     0.00000000e+00    # (T^_e)_{31}
  3  2     0.00000000e+00    # (T^_e)_{32}
  3  3    -3.52536855e+01    # (T^_e)_{33}
Block VCKM Q= 1.03033242e+03 # DRbar CKM mixing matrix
  1  1     9.73840738e-01    # CKM_{11} matrix element
  1  2     2.27197366e-01    # CKM_{12} matrix element
  1  3     3.94627287e-03    # CKM_{13} matrix element
  2  1    -2.27161577e-01    # CKM_{21} matrix element
  2  2     9.72962448e-01    # CKM_{22} matrix element
  2  3     4.17335952e-02    # CKM_{23} matrix element
  3  1     5.64218761e-03    # CKM_{31} matrix element
  3  2    -4.15383168e-02    # CKM_{32} matrix element
  3  3     9.99120981e-01    # CKM_{33} matrix element
Block msoft Q= 1.03033242e+03 # MSSM DRbar SUSY breaking parameters
     1     2.54573204e+02     # M_1(Q)
     2     4.68778478e+02     # M_2(Q)
     3     1.31944421e+03     # M_3(Q)
    21     1.37066780e+05     # mH1^2(Q)
    22    -5.22726073e+05     # mH2^2(Q)
Block USQMIX  # super CKM squark mass^2 matrix
  1  1     2.01392561e-05   # (USQMIX)_{11}
  1  2     2.13239016e-04   # (USQMIX)_{12}
  1  3     3.92332894e-01   # (USQMIX)_{13}
  1  4     1.13198004e-10   # (USQMIX)_{14}
  1  5     5.49673684e-07   # (USQMIX)_{15}
  1  6     9.19823273e-01   # (USQMIX)_{16}
  2  1     1.38644664e-04   # (USQMIX)_{21}
  2  2     1.46727706e-03   # (USQMIX)_{22}
  2  3     9.19822172e-01   # (USQMIX)_{23}
  2  4     4.63781313e-09   # (USQMIX)_{24}
  2  5     2.24875893e-05   # (USQMIX)_{25}
  2  6    -3.92332767e-01   # (USQMIX)_{26}
  3  1     1.21017549e-07   # (USQMIX)_{31}
  3  2     7.17361159e-03   # (USQMIX)_{32}
  3  3    -3.11816867e-05   # (USQMIX)_{33}
  3  4     1.96641084e-08   # (USQMIX)_{34}
  3  5     9.99974269e-01   # (USQMIX)_{35}
  3  6     1.10393440e-05   # (USQMIX)_{36}
  4  1     1.56601958e-05   # (USQMIX)_{41}
  4  2     1.23151289e-10   # (USQMIX)_{42}
  4  3    -6.43101075e-09   # (USQMIX)_{43}
  4  4     1.00000000e+00   # (USQMIX)_{44}
  4  5    -1.96676188e-08   # (USQMIX)_{45}
  4  6     2.27706666e-09   # (USQMIX)_{46}
  5  1     1.47032180e-01   # (USQMIX)_{51}
  5  2     9.89105140e-01   # (USQMIX)_{52}
  5  3    -1.43744652e-03   # (USQMIX)_{53}
  5  4    -2.30282420e-06   # (USQMIX)_{54}
  5  5    -7.09570549e-03   # (USQMIX)_{55}
  5  6     3.80599772e-04   # (USQMIX)_{56}
  6  1     9.89131699e-01   # (USQMIX)_{61}
  6  2    -1.47028443e-01   # (USQMIX)_{62}
  6  3     7.67553715e-05   # (USQMIX)_{63}
  6  4    -1.54899567e-05   # (USQMIX)_{64}
  6  5     1.05463500e-03   # (USQMIX)_{65}
  6  6    -2.03108721e-05   # (USQMIX)_{66}
Block DSQMIX  # super CKM squark mass^2 matrix
  1  1     4.61397475e-03   # (DSQMIX)_{11}
  1  2    -3.39694087e-02   # (DSQMIX)_{12}
  1  3     9.77345899e-01   # (DSQMIX)_{13}
  1  4     9.02851694e-07   # (DSQMIX)_{14}
  1  5    -1.45549051e-04   # (DSQMIX)_{15}
  1  6     2.08853447e-01   # (DSQMIX)_{16}
  2  1    -1.58575067e-03   # (DSQMIX)_{21}
  2  2     1.16756516e-02   # (DSQMIX)_{22}
  2  3    -2.08566458e-01   # (DSQMIX)_{23}
  2  4    -1.88624174e-06   # (DSQMIX)_{24}
  2  5     3.04199711e-04   # (DSQMIX)_{25}
  2  6     9.77937168e-01   # (DSQMIX)_{26}
  3  1     1.59641154e-06   # (DSQMIX)_{31}
  3  2     3.70195818e-03   # (DSQMIX)_{32}
  3  3     3.38090123e-04   # (DSQMIX)_{33}
  3  4     4.06313812e-06   # (DSQMIX)_{34}
  3  5     9.99993050e-01   # (DSQMIX)_{35}
  3  6    -2.83150678e-04   # (DSQMIX)_{36}
  4  1     1.69615803e-04   # (DSQMIX)_{41}
  4  2     5.78278437e-08   # (DSQMIX)_{42}
  4  3    -2.09854551e-06   # (DSQMIX)_{43}
  4  4     9.99999986e-01   # (DSQMIX)_{44}
  4  5    -4.06244420e-06   # (DSQMIX)_{45}
  4  6     1.75684548e-06   # (DSQMIX)_{46}
  5  1    -1.36834723e-01   # (DSQMIX)_{51}
  5  2     9.89923621e-01   # (DSQMIX)_{52}
  5  3     3.59847756e-02   # (DSQMIX)_{53}
  5  4     2.32203297e-05   # (DSQMIX)_{54}
  5  5    -3.67786509e-03   # (DSQMIX)_{55}
  5  6    -4.36495634e-03   # (DSQMIX)_{56}
  6  1     9.90581863e-01   # (DSQMIX)_{61}
  6  2     1.36920705e-01   # (DSQMIX)_{62}
  6  3     8.45794935e-05   # (DSQMIX)_{63}
  6  4    -1.68028128e-04   # (DSQMIX)_{64}
  6  5    -5.08490451e-04   # (DSQMIX)_{65}
  6  6    -1.02540236e-05   # (DSQMIX)_{66}
Block RVLMIX  # charged higgs-slepton mixing matrix 
  1 1    0.00000000e+00   # C_{11}
  1 2    0.00000000e+00   # C_{12}
  1 3    0.00000000e+00   # C_{13}
  1 4    0.00000000e+00   # C_{14}
  1 5    1.21350041e-01   # C_{15}
  1 6    0.00000000e+00   # C_{16}
  1 7    0.00000000e+00   # C_{17}
  1 8    9.92609776e-01   # C_{18}
  2 1    0.00000000e+00   # C_{21}
  2 2    0.00000000e+00   # C_{22}
  2 3    0.00000000e+00   # C_{23}
  2 4    7.19151420e-03   # C_{24}
  2 5    0.00000000e+00   # C_{25}
  2 6    0.00000000e+00   # C_{26}
  2 7    9.99974141e-01   # C_{27}
  2 8    0.00000000e+00   # C_{28}
  3 1   -4.09741393e-18   # C_{31}
  3 2    5.23344320e-17   # C_{32}
  3 3    3.47844614e-05   # C_{33}
  3 4    0.00000000e+00   # C_{34}
  3 5    0.00000000e+00   # C_{35}
  3 6    9.99999999e-01   # C_{36}
  3 7    0.00000000e+00   # C_{37}
  3 8    0.00000000e+00   # C_{38}
  4 1    1.37831259e-15   # C_{41}
  4 2    6.93925425e-16   # C_{42}
  4 3    9.99999999e-01   # C_{43}
  4 4    0.00000000e+00   # C_{44}
  4 5    0.00000000e+00   # C_{45}
  4 6   -3.47844614e-05   # C_{46}
  4 7    0.00000000e+00   # C_{47}
  4 8    0.00000000e+00   # C_{48}
  5 1    0.00000000e+00   # C_{51}
  5 2    0.00000000e+00   # C_{52}
  5 3    0.00000000e+00   # C_{53}
  5 4    9.99974141e-01   # C_{54}
  5 5    0.00000000e+00   # C_{55}
  5 6    0.00000000e+00   # C_{56}
  5 7   -7.19151420e-03   # C_{57}
  5 8    0.00000000e+00   # C_{58}
  6 1    0.00000000e+00   # C_{61}
  6 2    0.00000000e+00   # C_{62}
  6 3    0.00000000e+00   # C_{63}
  6 4    0.00000000e+00   # C_{64}
  6 5    9.92609776e-01   # C_{65}
  6 6    0.00000000e+00   # C_{66}
  6 7    0.00000000e+00   # C_{67}
  6 8   -1.21350041e-01   # C_{68}
  7 1    9.94680729e-01   # C_{71}
  7 2    1.03006053e-01   # C_{72}
  7 3   -1.44245954e-15   # C_{73}
  7 4    0.00000000e+00   # C_{74}
  7 5    0.00000000e+00   # C_{75}
  7 6   -1.26496943e-18   # C_{76}
  7 7    0.00000000e+00   # C_{77}
  7 8    0.00000000e+00   # C_{78}
Block RVHMIX  # CP-even neutral scalar mixing matrix 
  1 1    4.93154436e-01   # curlyN_{11}
  1 2    8.69941781e-01   # curlyN_{12}
  1 3    2.50128197e-11   # curlyN_{13}
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
  4 1   -4.52543327e-11   # curlyN_{41}
  4 2   -3.09841975e-12   # curlyN_{42}
  4 3    1.00000000e+00   # curlyN_{43}
  4 4    0.00000000e+00   # curlyN_{44}
  4 5    0.00000000e+00   # curlyN_{45}
  5 1    8.69941781e-01   # curlyN_{51}
  5 2   -4.93154436e-01   # curlyN_{52}
  5 3    3.78406353e-11   # curlyN_{53}
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
  3 1   -4.58838875e-11   # curlyN~_{31}
  3 2    3.18627954e-12   # curlyN~_{32}
  3 3    1.00000000e+00   # curlyN~_{33}
  3 4    0.00000000e+00   # curlyN~_{34}
  3 5    0.00000000e+00   # curlyN~_{35}
  4 1    8.71470395e-01   # curlyN~_{41}
  4 2    4.90448112e-01   # curlyN~_{42}
  4 3    3.84237448e-11   # curlyN~_{43}
  4 4    0.00000000e+00   # curlyN~_{44}
  4 5    0.00000000e+00   # curlyN~_{45}
