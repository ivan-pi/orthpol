!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:34:28

SUBROUTINE slancz(n,ncap,x,w,alpha,beta,ierr,p0,p1)

! This routine carries out the same task as the routine  sti, but
! uses the more stable Lanczos method. The meaning of the input
! and output parameters is the same as in the routine  sti. (This
! routine is adapted from the routine RKPW in W.B. Gragg and
! W.J. Harrod,The numerically stable reconstruction of Jacobi
! matrices from spectral data'', Numer. Math. 44, 1984, 317-335.)



INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                      :: ncap
REAL, INTENT(IN)                         :: x(ncap)
REAL, INTENT(IN)                         :: w(ncap)
REAL, INTENT(OUT)                        :: alpha(n)
REAL, INTENT(OUT)                        :: beta(n)
INTEGER, INTENT(OUT)                     :: ierr
REAL, INTENT(OUT)                        :: p0(ncap)
REAL, INTENT(OUT)                        :: p1(ncap)

IF(n <= 0 .OR. n > ncap) THEN
  ierr=1
  RETURN
ELSE
  ierr=0
END IF
DO  i=1,ncap
  p0(i)=x(i)
  p1(i)=0.
END DO
p1(1)=w(1)
DO  i=1,ncap-1
  pi=w(i+1)
  gam=1.
  sig=0.
  t=0.
  xlam=x(i+1)
  DO  k=1,i+1
    rho=p1(k)+pi
    tmp=gam*rho
    tsig=sig
    IF(rho <= 0.) THEN
      gam=1.
      sig=0.
    ELSE
      gam=p1(k)/rho
      sig=pi/rho
    END IF
    tk=sig*(p0(k)-xlam)-gam*t
    p0(k)=p0(k)-(tk-t)
    t=tk
    IF(sig <= 0.) THEN
      pi=tsig*p1(k)
    ELSE
      pi=(t**2)/sig
    END IF
    tsig=sig
    p1(k)=tmp
  END DO
END DO
DO  k=1,n
  alpha(k)=p0(k)
  beta(k)=p1(k)
END DO
RETURN
END SUBROUTINE slancz

