!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:34:20

SUBROUTINE smccheb(n,ncapm,mc,mp,xp,yp,squad,eps,iq,idelta,  &
    finl,finr,endl,endr,xfer,wfer,a,b,fnu,alpha,beta,ncap,kount,  &
    ierr,be,x,w,xm,wm,s,s0,s1,s2)

! This is a multiple-component discretized modified Chebyshev
! algorithm, basically a modified Chebyshev algorithm in which the
! modified moments are discretized in the same manner as the inner
! product in the discretization procedure  mcdis. The input and
! output parameters are as in  mcdis. In addition, the arrays  a,b
! must be filled with the recursion coefficients  a(k-1),b(k-1),
! k=1,2,...,2*n-1, defining the modified moments. The arrays
! be,x,w,xm,wm,s,s0,s1,s2  are used for working space. The routine
! calls upon the subroutine  cheb. The routine exits immediately with
! ierr=-1  if  n  is not in range.


INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: ncapm
INTEGER, INTENT(IN)                      :: mc
INTEGER, INTENT(IN)                      :: mp
REAL, INTENT(IN)                         :: xp(*)
REAL, INTENT(IN)                         :: yp(*)
REAL, INTENT(IN OUT)                     :: squad
REAL, INTENT(IN)                         :: eps
INTEGER, INTENT(IN OUT)                  :: iq
INTEGER, INTENT(OUT)                     :: idelta
LOGICAL, INTENT(IN OUT)                  :: finl
LOGICAL, INTENT(IN OUT)                  :: finr
REAL, INTENT(IN OUT)                     :: endl(mc)
REAL, INTENT(IN OUT)                     :: endr(mc)
REAL, INTENT(IN OUT)                     :: xfer(ncapm)
REAL, INTENT(IN OUT)                     :: wfer(ncapm)
REAL, INTENT(IN)                         :: a(*)
REAL, INTENT(IN)                         :: b(*)
REAL, INTENT(OUT)                        :: fnu(*)
REAL, INTENT(IN OUT)                     :: alpha(n)
REAL, INTENT(OUT)                        :: beta(n)
INTEGER, INTENT(OUT)                     :: ncap
INTEGER, INTENT(OUT)                     :: kount
INTEGER, INTENT(OUT)                     :: ierr
REAL, INTENT(OUT)                        :: be(n)
REAL, INTENT(IN)                         :: x(ncapm)
REAL, INTENT(IN)                         :: w(ncapm)
REAL, INTENT(OUT)                        :: xm(*)
REAL, INTENT(OUT)                        :: wm(*)
REAL, INTENT(IN OUT)                     :: s(n)
REAL, INTENT(IN OUT)                     :: s0(*)
REAL, INTENT(IN OUT)                     :: s1(*)
REAL, INTENT(IN OUT)                     :: s2(*)



! The arrays  xp,yp  are assumed to have dimension  mp  if mp > 0,
! the arrays  a,b  dimension 2*n-1, the arrays  fnu,s0,s1,s2  dimension
! 2*n, and the arrays  xm,wm  dimension  mc*ncapm+mp.

nd=2*n
IF(idelta <= 0) idelta=1
IF(n < 1) THEN
  ierr=-1
  RETURN
END IF

! Initialization

incap=1
kount=-1
ierr=0
DO  k=1,n
  beta(k)=0.
END DO
ncap=(nd-1)/idelta
20 DO  k=1,n
  be(k)=beta(k)
END DO
kount=kount+1
IF(kount > 1) incap=2**(kount/5)*n
ncap=ncap+incap
IF(ncap > ncapm) THEN
  ierr=ncapm
  RETURN
END IF

! Discretization of the modified moments

mtncap=mc*ncap
DO  i=1,mc
  im1tn=(i-1)*ncap
  IF(iq == 1) THEN
    CALL squad(ncap,x,w,i,ierr)
  ELSE
    CALL sqgp(ncap,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,wfer)
  END IF
  IF(ierr /= 0) THEN
    ierr=i
    RETURN
  END IF
  DO  k=1,ncap
    xm(im1tn+k)=x(k)
    wm(im1tn+k)=w(k)
  END DO
END DO
IF(mp /= 0) THEN
  DO  k=1,mp
    xm(mtncap+k)=xp(k)
    wm(mtncap+k)=yp(k)
  END DO
END IF
mtnpmp=mtncap+mp
DO  k=1,nd
  km1=k-1
  sum=0.
  DO  i=1,mtnpmp
    p1=0.
    p=1.
    IF(k > 1) THEN
      DO  l=1,km1
        pm1=p1
        p1=p
        p=(xm(i)-a(l))*p1-b(l)*pm1
      END DO
    END IF
    sum=sum+wm(i)*p
  END DO
  fnu(k)=sum
END DO

! Computation of the desired recursion coefficients

CALL scheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)

! In the following statement, the absolute value of the beta's is
! used to guard against failure in cases where the routine is applied
! to variable-sign weight functions and hence the positivity of the
! beta's is not guaranteed.

DO  k=1,n
  IF(ABS(beta(k)-be(k)) > eps*ABS(beta(k))) GO TO 20
END DO
RETURN
END SUBROUTINE smccheb

