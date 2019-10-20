!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:35:13

FUNCTION nu0her(n,z,eps)

! This is an auxiliary function routine providing a starting backward
! recurrence index for the Hermite measure that can be used in place
! of  nu0  in the routines  knum  and  dknum.



INTEGER, INTENT(IN OUT)                  :: n
COMPLEX, INTENT(IN OUT)                  :: z
REAL, INTENT(IN OUT)                     :: eps

nu0her=2.*(SQRT(.5*REAL(n+1))+.25*ALOG(1./eps)/ ABS(AIMAG(z)))**2
RETURN
END FUNCTION nu0her

