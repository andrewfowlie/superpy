Conjugate[m_Mass] ^= m

M$CouplingMatrices = M$CouplingMatrices /.
    { Af[t_, g1_, g2_] -> Kf[g1, g2, t]/Mass[F[t, {g1}]],
      AfC[t_, g1_, g2_] -> KfC[g1, g2, t]/Mass[F[t, {g1}]] } //.
    a_/m_Mass + b_ -> (a + m b)/m


