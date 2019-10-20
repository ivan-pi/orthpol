!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:34:23

SUBROUTINE slob(n,alpha,beta,aleft,right,zero,weight,ierr,e,a,b)

! Given  n  and a measure  dlambda, this routine generates the
! (n+2)-point Gauss-Lobatto quadrature formula

!   integral over supp(dlambda) of f(x)dlambda(x)

!      = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k))

!              + w(n+1)f(x(n+1)) + R(n;f).

! The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
! =w(k), k=0,1,...,n,n+1. The user has to supply the recursion
! coefficients  alpha(k), beta(k), k=0,1,...,n,n+1, for the measure
! dlambda. The nodes and weights are computed in terms of the
! eigenvalues and first component of the normalized eigenvectors of
! a slightly modified Jacobi matrix of order  n+2. The routine calls
! upon the subroutine  gauss  and the function subroutine  r1mach.

!   Input:  n - -  the number of interior points in the Gauss-Lobatto
!                  formula; type integer
!           alpha,beta - arrays of dimension  n+2  to be supplied with
!                  the recursion coefficients  alpha(k-1), beta(k-1),
!                  k=1,2,...,n+2, of the underlying measure; the
!                  routine does not use  alpha(n+2), beta(n+2)
!           aleft,right - the prescribed left and right endpoints
!                  x(0)  and  x(n+1)  of the Gauss-Lobatto formula

!   Output: zero - an array of dimension  n+2  containing the nodes (in
!                  increasing order)  zero(k)=x(k), k=0,1,...,n,n+1
!           weight-an array of dimension  n+2  containing the weights
!                  weight(k)=w(k), k=0,1,...,n,n+1
!           ierr - an error flag inherited from the routine  gauss

! The arrays  e,a,b  are needed for working space.



INTEGER, INTENT(IN)                      :: n
REAL, INTENT(IN)                         :: alpha(*)
REAL, INTENT(IN)                         :: beta(*)
REAL, INTENT(IN OUT)                     :: aleft
REAL, INTENT(IN)                         :: right
REAL, INTENT(IN OUT)                     :: zero(*)
REAL, INTENT(IN OUT)                     :: weight(*)
INTEGER, INTENT(IN OUT)                  :: ierr
REAL, INTENT(IN OUT)                     :: e(*)
REAL, INTENT(OUT)                        :: a(*)
REAL, INTENT(OUT)                        :: b(*)


! The arrays  alpha,beta,zero,weight,e,a,b  are assumed to have
! dimension  n+2.

epsma=r1mach(3)

! epsma is the machine single precision.

np1=n+1
np2=n+2
DO  k=1,np2
  a(k)=alpha(k)
  b(k)=beta(k)
END DO
p0l=0.
p0r=0.
p1l=1.
p1r=1.
DO  k=1,np1
  pm1l=p0l
  p0l=p1l
  pm1r=p0r
  p0r=p1r
  p1l=(aleft-a(k))*p0l-b(k)*pm1l
  p1r=(right-a(k))*p0r-b(k)*pm1r
END DO
det=p1l*p0r-p1r*p0l
a(np2)=(aleft*p1l*p0r-right*p1r*p0l)/det
b(np2)=(right-aleft)*p1l*p1r/det
CALL sgauss(np2,a,b,epsma,zero,weight,ierr,e)
RETURN
END SUBROUTINE slob

