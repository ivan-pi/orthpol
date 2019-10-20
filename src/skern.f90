!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:34:36

SUBROUTINE skern(n,nu0,numax,z,eps,a,b,ker,nu,ierr,rold)

! This routine generates the kernels in the Gauss quadrature remainder
! term, namely

!           K(k)(z)=rho(k)(z)/pi(k)(z), k=0,1,2,...,n,

! where  rho(k)  are the output quantities of the routine  knum, and
! pi(k)  the (monic) orthogonal polynomials. The results are returned
! in the array  ker  as ker(k)=K(k-1)(z), k=1,2,...,n+1. All the other
! input and output parameters have the same meaning as in the routine
! knum.


INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: nu0
INTEGER, INTENT(IN OUT)                  :: numax
COMPLEX, INTENT(IN OUT)                  :: z
REAL, INTENT(IN OUT)                     :: eps
REAL, INTENT(IN)                         :: a(numax)
REAL, INTENT(IN)                         :: b(numax)
COMPLEX, INTENT(OUT)                     :: ker(*)
INTEGER, INTENT(IN OUT)                  :: nu
INTEGER, INTENT(IN OUT)                  :: ierr
COMPLEX, INTENT(IN OUT)                  :: rold(*)
COMPLEX :: p0,p,pm1


! The arrays  ker,rold  are assumed to have dimension  n+1.

CALL knum(n,nu0,numax,z,eps,a,b,ker,nu,ierr,rold)
p0=(0.,0.)
p=(1.,0.)
DO  k=1,n
  pm1=p0
  p0=p
  p=(z-a(k))*p0-b(k)*pm1
  ker(k+1)=ker(k+1)/p
END DO
RETURN
END SUBROUTINE skern

