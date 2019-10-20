!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:35:08

FUNCTION nu0jac(n,z,eps)

! This is an auxiliary function routine providing a starting backward
! recurrence index for the Jacobi measure that can be used in place
! of  nu0  in the routines  knum  and  dknum.



INTEGER, INTENT(IN OUT)                  :: n
COMPLEX, INTENT(IN OUT)                  :: z
REAL, INTENT(IN OUT)                     :: eps

pi=4.*ATAN(1.)
x=REAL(z)
y=ABS(AIMAG(z))
IF(x < 1.) THEN
  IF(x < -1.) angle=.5*(2.*pi+ATAN(y/(x-1.))+ATAN(y/(x+1.)))
  IF(x == -1.) angle=.5*(1.5*pi-ATAN(.5*y))
  IF(x > -1.) angle=.5*(pi+ATAN(y/(x-1.))+ATAN(y/(x+1.)))
ELSE
  IF(x == 1.) angle=.5*(.5*pi+ATAN(.5*y))
  IF(x > 1.) angle=.5*(ATAN(y/(x-1.))+ATAN(y/(x+1.)))
END IF
x2=x*x
y2=y*y
r=((x2-y2-1.)**2+4.*x2*y2)**.25
r=SQRT((x+r*COS(angle))**2+(y+r*SIN(angle))**2)
!       nu0jac=real(n+1)+.5*alog(1./eps)/alog(r)
nu0jac=INT(REAL(n+1)+.5*ALOG(1./eps)/ALOG(r))
RETURN
END FUNCTION nu0jac

