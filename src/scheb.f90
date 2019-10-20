!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:34:56

SUBROUTINE scheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)

! Given a set of polynomials  p(0),p(1),...,p(2*n-1)  satisfying

!        p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
!                        k=0,1,...,2*n-2,

!        p(-1)(x)=0,  p(0)(x)=1,

! and associated modified moments

!           fnu(k)=integral of p(k)(x)*dlambda(x),
!                        k=0,1,...,2*n-1,

! this subroutine uses the modified Chebyshev algorithm (see, e.g.,
! Section 2.4 of W. Gautschi,On generating orthogonal polynomials'',
! SIAM J. Sci. Statist. Comput. 3, 1982, 289-317) to generate the
! recursion coefficients  alpha(k),beta(k), k=0,1,...,n-1, for the
! polynomials  pi(k)  orthogonal with respect to the integration
! measure  dlambda(x), i.e.,

!        pi(k+1)(x)=(x-alpha(k))*pi(k)(x)-beta(k)*pi(k-1)(x),
!                          k=0,1,...,n-1,

!        pi(-1)(x)=0,  pi(0)(x)=1.

!     Input:    n - - the number of recursion coefficients desired
!               a,b-- arrays of dimension 2*n-1 to be filled with the
!                     values of  a(k-1),b(k-1), k=1,2,...,2*n-1
!               fnu-- array of dimension  2*n  to be filled with the
!                     values of the modified moments  fnu(k-1), k=1,2,
!                     ...,2*n
!     Output:   alpha,beta-- arrays containing, respectively, the
!                     recursion coefficients  alpha(k-1),beta(k-1),
!                     k=1,2,...,n, where  beta(0)  is the total mass.
!               s - - array containing the normalization factors
!                     s(k)=integral [pi(k)(x)]**2 dlambda(x), k=0,1,
!                     2,...,n-1.
!               ierr- an error flag, equal to  0  on normal return,
!                     equal to  1  if  abs(fnu(0))  is less than the
!                     machine zero, equal to  2  if  n  is out of range,
!                     equal to  -k  if  s(k), k=0,1,2,...,n-1, is about
!                     to underflow, and equal to  +k  if it is about to
!                     overflow.

! The arrays  s0,s1,s2  are needed for working space.

! On machines with limited exponent range, the occurrence of underflow
! [overflow] in the computation of the  alpha's  and  beta's  can often
! be avoided by multiplying all modified moments by a sufficiently large
! [small] scaling factor and dividing the new  beta(0)  by the same
! scaling factor.

! The routine uses the function subroutine  r1mach.



INTEGER, INTENT(IN)                      :: n
REAL, INTENT(IN)                         :: a(*)
REAL, INTENT(IN OUT)                     :: b(*)
REAL, INTENT(IN)                         :: fnu(*)
REAL, INTENT(OUT)                        :: alpha(n)
REAL, INTENT(OUT)                        :: beta(n)
REAL, INTENT(OUT)                        :: s(n)
INTEGER, INTENT(OUT)                     :: ierr
REAL, INTENT(OUT)                        :: s0(*)
REAL, INTENT(OUT)                        :: s1(*)
REAL, INTENT(OUT)                        :: s2(*)


! The arrays  a,b  are assumed to have dimension  2*n-1, the arrays
! fnu,s0,s1,s2  dimension  2*n.

nd=2*n
tiny=10.*r1mach(1)
huge=.1*r1mach(2)
ierr=0
IF(ABS(fnu(1)) < tiny) THEN
  ierr=1
  RETURN
END IF
IF(n < 1) THEN
  ierr=2
  RETURN
END IF

! Initialization

alpha(1)=a(1)+fnu(2)/fnu(1)
beta(1)=fnu(1)
IF(n == 1) RETURN
s(1)=fnu(1)
DO  l=1,nd
  s0(l)=0.
  s1(l)=fnu(l)
END DO

! Continuation

DO  k=2,n
  lk=nd-k+1
  DO  l=k,lk
    
! The quantities  s2(l)  for l > k are auxiliary quantities which may
! be zero or may become so small as to underflow, without however
! causing any harm.
    
    s2(l)=s1(l+1)-(alpha(k-1)-a(l))*s1(l)-beta(k-1)*s0(l) +b(l)*s1(l-1)
    IF(l == k) s(k)=s2(k)
    
! Check impending underflow or overflow
    
    IF(ABS(s(k)) < tiny) THEN
      ierr=-(k-1)
      RETURN
    ELSE IF(ABS(s(k)) > huge) THEN
      ierr=k-1
      RETURN
    END IF
  END DO
  
! Compute the alpha- and beta-coefficient
  
  alpha(k)=a(k)+(s2(k+1)/s2(k))-(s1(k)/s1(k-1))
  beta(k)=s2(k)/s1(k-1)
  DO  l=k,lk
    s0(l)=s1(l)
    s1(l)=s2(l)
  END DO
END DO
RETURN
END SUBROUTINE scheb

