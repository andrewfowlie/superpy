* Eigen2x2.f
* diagonalization of a Hermitian 2-by-2 matrix
* this file is part of FormCalc
* last modified 15 Jun 04 th


	subroutine Eigen2x2(m, msq, U, M11, M22, M12, *)
	implicit none
	double precision m(*), msq(*)
	double complex U(2,*)
	double precision M11, M22
	double complex M12

	double precision m1, m2, delta, tau
	double complex s

	call Jacobi(delta, s, tau, M11 - M22, M12, abs(M12))

	m1 = M11 + delta
	m2 = M22 - delta

	if( m1 .gt. m2 ) then
	  msq(1) = m2
	  msq(2) = m1
	  U(1,1) = -dconjg(s)
	  U(2,2) = s
	  U(1,2) = 1 - tau
	  U(2,1) = U(1,2)
	else
	  msq(1) = m1
	  msq(2) = m2
	  U(1,1) = 1 - tau
	  U(2,2) = U(1,1)
	  U(1,2) = s
	  U(2,1) = -dconjg(s)
	endif

	if( msq(1) .lt. 0 ) then
	  print *, "Negative mass-squares."
	  return 1
	endif
	m(1) = sqrt(msq(1))
	m(2) = sqrt(msq(2))
	end

