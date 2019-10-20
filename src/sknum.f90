!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:34:32

SUBROUTINE sknum(n,nu0,numax,z,eps,a,b,rho,nu,ierr,rold)

! This routine generates

!   rho(k)(z)=integral pi(k)(t)dlambda(t)/(z-t), k=0,1,2,...,n,

! where  pi(k)(t)  is the (monic) k-th degree orthogonal polynomial
! with respect to the measure  dlambda(t), and the integral is extended
! over the support of  dlambda. It is assumed that  z  is a complex
! number outside the smallest interval containing the support of
! dlambda. The quantities  rho(k)(z)  are computed as the first  n+1
! members of the minimal solution of the basic three-term recurrence
! relation

!      y(k+1)(z)=(z-a(k))y(k)(z)-b(k)y(k-1)(z), k=0,1,2,...,

! satisfied by the orthogonal polynomials  pi(k)(z).

!   Input:  n  - -  the largest integer  k  for which  rho(k)  is
!                   desired
!           nu0  -  an estimate of the starting backward recurrence
!                   index; if no better estimate is known, set
!                   nu0 = 3*n/2; for Jacobi, Laguerre and Hermite
!                   weight functions, estimates of  nu0  are generated
!                   respectively by the routines  nu0jac,nu0lag  and
!                   nu0her
!           numax - an integer larger than  n  cutting off backward
!                   recursion in case of nonconvergence; if  nu0
!                   exceeds  numax, then the routine aborts with the
!                   error flag  ierr  set equal to  nu0
!           z - - - the variable in  rho(k)(z); type complex
!           eps - - the relative accuracy to which the  rho(k)  are
!                   desired
!           a,b - - arrays of dimension  numax  to be supplied with the
!                   recurrence coefficients  a(k-1), b(k-1), k=1,2,...,
!                   numax.

!   Output: rho - - an array of dimension  n+1  containing the results
!                   rho(k)=rho(k-1)(z), k=1,2,...,n+1; type complex
!           nu  - - the starting backward recurrence index that yields
!                   convergence
!           ierr  - an error flag equal to zero on normal return, equal
!                   to  nu0  if  nu0 > numax, and equal to  numax in
!                   case of nonconvergence.

! The complex array  rold  of dimension  n+1  is used for working space.


INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: nu0
INTEGER, INTENT(IN)                      :: numax
COMPLEX, INTENT(IN OUT)                  :: z
REAL, INTENT(IN)                         :: eps
REAL, INTENT(IN OUT)                     :: a(numax)
REAL, INTENT(IN OUT)                     :: b(numax)
COMPLEX, INTENT(OUT)                     :: rho(*)
INTEGER, INTENT(OUT)                     :: nu
INTEGER, INTENT(OUT)                     :: ierr
COMPLEX, INTENT(OUT)                     :: rold(*)
COMPLEX :: r


! The arrays  rho,rold  are assumed to have dimension  n+1.

ierr=0
np1=n+1
IF(nu0 > numax) THEN
  ierr=nu0
  RETURN
END IF
IF(nu0 < np1) nu0=np1
nu=nu0-5
DO  k=1,np1
  rho(k)=(0.,0.)
END DO
20 nu=nu+5
IF(nu > numax) THEN
  ierr=numax
  GO TO 60
END IF
DO  k=1,np1
  rold(k)=rho(k)
END DO
r=(0.,0.)
DO  j=1,nu
  j1=nu-j+1
  r=CMPLX(b(j1),0.)/(z-CMPLX(a(j1),0.)-r)
  IF(j1 <= np1) rho(j1)=r
END DO
DO  k=1,np1
  IF(ABS(rho(k)-rold(k)) > eps*ABS(rho(k))) GO TO 20
END DO
60 IF(n == 0) RETURN
DO  k=2,np1
  rho(k)=rho(k)*rho(k-1)
END DO
RETURN
END SUBROUTINE sknum

