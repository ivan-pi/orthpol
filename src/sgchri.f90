!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:34:40

SUBROUTINE sgchri(n,iopt,nu0,numax,eps,a,b,x,y,alpha,beta,nu,  &
    ierr,ierrc,fnu,rho,rold,s,s0,s1,s2)

! This routine implements the generalized Christoffel theorem, using
! the method of modified moments (cf. Section 4 of W. Gautschi,
! Minimal solutions of three-term recurrence relations and orthogonal
! polynomials'', Math. Comp. 36, 1981, 547-554). Given the recursion
! coefficients  a(k), b(k), k=0,1,...n, for the (monic) orthogonal
! polynomials with respect to some measure  dlambda(t), it generates
! the recursion coefficients  alpha(k), beta(k), k=0,1,2,...,n-1 for
! the measure

!         dlambda(t)/(t-x)        if iopt=1
!         dlambda(t)/{(t-x)**2+y**2} if iopt=2

!   Input:  n   - - the number of recurrence coefficients desired;
!                   type integer
!           iopt  - an integer selecting the desired weight distribution
!           nu0   - an integer estimating the starting backward
!                   recurrence index; in the absence of any better
!                   choice, take  nu0 = 3*n
!           numax - an integer controlling termination of backward
!                   recursion in case of nonconvergence; a conservative
!                   choice is  numax = 500
!           eps - - a relative error tolerance; type real
!           a,b - - arrays of dimension numax to be supplied with the
!                   recursion coefficients a(k)=alpha(k-1),b(k)=beta(k),
!                   k=1,2,...,numax, for the measure  dlambda
!           x,y - - real parameters defining the linear and quadratic
!                   divisors of  dlambda

!   Output: alpha,beta - arrays of dimension  n  containing the desired
!                   recursion coefficients  alpha(k-1), beta(k-1), k=1,
!                   2,...,n
!           nu  - - the backward recurrence index yielding convergence;
!                   in case of nonconvergence,  nu  will have the value
!                   numax
!           ierr  - an error flag, where
!                   ierr=0     on normal return
!                   ierr=1     if  iopt  is neither 1 nor 2
!                   ierr=nu0   if  nu0 > numax
!                   ierr=numax if the backward recurrence algorithm does
!                              not converge
!                   ierr=-1    if  n  is not in range
!           ierrc - an error flag inherited from the routine  cheb

! The arrays  fnu,s,s0,s1,s2  are working space. The routine calls
! upon the routines  knum  and  cheb.


INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: iopt
INTEGER, INTENT(IN OUT)                  :: nu0
INTEGER, INTENT(IN OUT)                  :: numax
REAL, INTENT(IN OUT)                     :: eps
REAL, INTENT(IN OUT)                     :: a(numax)
REAL, INTENT(IN OUT)                     :: b(numax)
REAL, INTENT(IN OUT)                     :: x
REAL, INTENT(OUT)                        :: y
REAL, INTENT(IN OUT)                     :: alpha(n)
REAL, INTENT(IN OUT)                     :: beta(n)
INTEGER, INTENT(IN OUT)                  :: nu
INTEGER, INTENT(OUT)                     :: ierr
INTEGER, INTENT(IN OUT)                  :: ierrc
REAL, INTENT(OUT)                        :: fnu(*)
COMPLEX, INTENT(IN OUT)                  :: rho(*)
COMPLEX, INTENT(IN OUT)                  :: rold(*)
REAL, INTENT(IN OUT)                     :: s(n)
REAL, INTENT(IN OUT)                     :: s0(*)
REAL, INTENT(IN OUT)                     :: s1(*)
REAL, INTENT(IN OUT)                     :: s2(*)
COMPLEX :: z


! The arrays  fnu,rho,rold,s0,s1,s2  are assumed to have dimension  2*n.

IF(n < 1) THEN
  ierr=-1
  RETURN
END IF
ierr=0
nd=2*n
ndm1=nd-1

! Linear divisor

IF(iopt == 1) THEN
  
! Generate the modified moments of  dlambda.
  
  z=CMPLX(x,0.)
  CALL sknum(ndm1,nu0,numax,z,eps,a,b,rho,nu,ierr,rold)
  DO  k=1,nd
    fnu(k)=-REAL(rho(k))
  END DO
  
! Compute the desired recursion coefficients by means of the modified
! Chebyshev algorithm.
  
  CALL scheb(n,a,b,fnu,alpha,beta,s,ierrc,s0,s1,s2)
  RETURN
  
! Quadratic divisor
  
ELSE IF(iopt == 2) THEN
  
! Generate the modified moments of  dlambda.
  
  y=ABS(y)
  z=CMPLX(x,y)
  CALL sknum(ndm1,nu0,numax,z,eps,a,b,rho,nu,ierr,rold)
  DO  k=1,nd
    fnu(k)=-AIMAG(rho(k))/y
  END DO
  
! Compute the desired recursion coefficients by means of the modified
! Chebyshev algorithm.
  
  CALL scheb(n,a,b,fnu,alpha,beta,s,ierrc,s0,s1,s2)
  RETURN
ELSE
  ierr=1
  RETURN
END IF
END SUBROUTINE sgchri

