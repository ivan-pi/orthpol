!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:35:05

FUNCTION nu0lag(n,z,al,eps)

! This is an auxiliary function routine providing a starting backward
! recurrence index for the Laguerre measure that can be used in place
! of  nu0  in the routines  knum  and  dknum.



INTEGER, INTENT(IN OUT)                  :: n
COMPLEX, INTENT(IN OUT)                  :: z
REAL, INTENT(IN OUT)                     :: al
REAL, INTENT(IN OUT)                     :: eps

pi=4.*ATAN(1.)
x=REAL(z)
y=AIMAG(z)
phi=.5*pi
IF(y < 0.) phi=1.5*pi
IF(x == 0.) GO TO 10
phi=ATAN(y/x)
IF(y > 0. .AND. x > 0.) GO TO 10
phi=phi+pi
IF(x < 0.) GO TO 10
phi=phi+pi
!    10 nu0lag=(sqrt(real(n+1)+.5*(al+1.))+alog(1./eps)/(4.*(x*x+
!      *  y*y)**.25*cos(.5*(phi-pi))))**2-.5*(al+1.)
10 nu0lag=INT((SQRT(REAL(n+1)+.5*(al+1.))+ALOG(1./eps)/(4.*(x*x+  &
    y*y)**.25*COS(.5*(phi-pi))))**2-.5*(al+1.))
RETURN
END FUNCTION nu0lag

