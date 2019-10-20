!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:30:50

SUBROUTINE srecur(n,ipoly,al,be,a,b,ierr)

! This subroutine generates the coefficients  a(k),b(k), k=0,1,...,n-1,
! in the recurrence relation

!       p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
!                            k=0,1,...,n-1,

!       p(-1)(x)=0,  p(0)(x)=1,

! for some classical (monic) orthogonal polynomials, and sets  b(0)
! equal to the total mass of the weight distribution. The results are
! stored in the arrays  a,b,  which hold, respectively, the coefficients
! a(k-1),b(k-1), k=1,2,...,n.

!       Input:  n - - the number of recursion coefficients desired
!               ipoly-integer identifying the polynomial as follows:
!                     1=Legendre polynomial on (-1,1)
!                     2=Legendre polynomial on (0,1)
!                     3=Chebyshev polynomial of the first kind
!                     4=Chebyshev polynomial of the second kind
!                     5=Jacobi polynomial with parameters  al=-.5,be=.5
!                     6=Jacobi polynomial with parameters  al,be
!                     7=generalized Laguerre polynomial with
!                       parameter  al
!                     8=Hermite polynomial
!               al,be-input parameters for Jacobi and generalized
!                     Laguerre polynomials

!       Output: a,b - arrays containing, respectively, the recursion
!                     coefficients  a(k-1),b(k-1), k=1,2,...,n.
!               ierr -an error flag, equal to  0  on normal return,
!                     equal to  1  if  al  or  be  are out of range
!                     when  ipoly=6  or  ipoly=7, equal to  2  if  b(0)
!                     overflows when  ipoly=6  or  ipoly=7, equal to  3
!                     if  n  is out of range, and equal to  4  if  ipoly
!                     is not an admissible integer. In the case  ierr=2,
!                     the coefficient  b(0)  is set equal to the largest
!                     machine-representable number.

! The subroutine calls for the function subroutines  r1mach,gamma  and
! alga. The routines  gamma  and  alga , which are included in this
! file, evaluate respectively the gamma function and its logarithm for
! positive arguments. They are used only in the cases  ipoly=6  and
! ipoly=7.


INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: ipoly
REAL, INTENT(IN)                         :: al
REAL, INTENT(IN)                         :: be
REAL, INTENT(OUT)                        :: a(n)
REAL, INTENT(OUT)                        :: b(n)
INTEGER, INTENT(OUT)                     :: ierr
EXTERNAL gamma


IF(n < 1) THEN
  ierr=3
  RETURN
END IF
almach=LOG(r1mach(2))
ierr=0
DO  k=1,n
  a(k)=0.
END DO
IF(ipoly == 1) THEN
  b(1)=2.
  IF (n == 1) RETURN
  DO  k=2,n
    fkm1=REAL(k-1)
    b(k)=1./(4.-1./(fkm1*fkm1))
  END DO
  RETURN
ELSE IF (ipoly == 2) THEN
  a(1)=.5
  b(1)=1.
  IF(n == 1) RETURN
  DO  k=2,n
    a(k)=.5
    fkm1=REAL(k-1)
    b(k)=.25/(4.-1./(fkm1*fkm1))
  END DO
  RETURN
ELSE IF(ipoly == 3) THEN
  b(1)=4.*ATAN(1.)
  IF(n == 1) RETURN
  b(2)=.5
  IF(n == 2) RETURN
  DO  k=3,n
    b(k)=.25
  END DO
  RETURN
ELSE IF(ipoly == 4) THEN
  b(1)=2.*ATAN(1.)
  IF(n == 1) RETURN
  DO  k=2,n
    b(k)=.25
  END DO
  RETURN
ELSE IF(ipoly == 5) THEN
  b(1)=4.*ATAN(1.)
  a(1)=.5
  IF(n == 1) RETURN
  DO  k=2,n
    b(k)=.25
  END DO
  RETURN
ELSE IF(ipoly == 6) THEN
  IF(al <= -1. .OR. be <= -1.) THEN
    ierr=1
    RETURN
  ELSE
    alpbe=al+be
    a(1)=(be-al)/(alpbe+2.)
    t=(alpbe+1.)*LOG(2.)+salga(al+1.)+salga(be+1.)- salga(alpbe+2.)
    IF(t > almach) THEN
      ierr=2
      b(1)=r1mach(2)
    ELSE
      b(1)=EXP(t)
    END IF
    IF(n == 1) RETURN
    al2=al*al
    be2=be*be
    a(2)=(be2-al2)/((alpbe+2.)*(alpbe+4.))
    b(2)=4.*(al+1.)*(be+1.)/((alpbe+3.)*(alpbe+2.)**2)
    IF(n == 2) RETURN
    DO  k=3,n
      fkm1=REAL(k-1)
      a(k)=.25*(be2-al2)/(fkm1*fkm1*(1.+.5*alpbe/fkm1)*  &
          (1.+.5*(alpbe+2.)/fkm1))
      b(k)=.25*(1.+al/fkm1)*(1.+be/fkm1)*(1.+alpbe/fkm1)/  &
          ((1.+.5*(alpbe+1.)/fkm1)*(1.+.5*(alpbe-1.)/fkm1) *(1.+.5*alpbe/fkm1)**2)
    END DO
    RETURN
  END IF
ELSE IF(ipoly == 7) THEN
  IF(al <= -1.) THEN
    ierr=1
    RETURN
  ELSE
    a(1)=al+1.
    b(1)=sgamma(al+1.,ierr)
    IF(ierr == 2) b(1)=r1mach(2)
    IF(n == 1) RETURN
    DO  k=2,n
      fkm1=REAL(k-1)
      a(k)=2.*fkm1+al+1.
      b(k)=fkm1*(fkm1+al)
    END DO
    RETURN
  END IF
ELSE IF(ipoly == 8) THEN
  b(1)=SQRT(4.*ATAN(1.))
  IF(n == 1) RETURN
  DO  k=2,n
    b(k)=.5*REAL(k-1)
  END DO
  RETURN
ELSE
  ierr=4
END IF
END SUBROUTINE srecur

FUNCTION salga(x)

! This is an auxiliary function subroutine (not optimized in any
! sense) evaluating the logarithm of the gamma function for positive
! arguments  x. It is called by the subroutine  gamma. The integer  m0
! in the first executable statement is the smallest integer  m  such
! that  1*3*5* ... *(2*m+1)/(2**m)  is greater than or equal to the
! largest machine-representable number. The routine is based on a
! rational approximation valid on [.5,1.5] due to W.J. Cody and
! K.E. Hillstrom; see Math. Comp. 21, 1967, 198-203, in particular the
! case  n=7  in Table II. For the computation of  m0  it calls upon the
! function subroutines  t  and  r1mach. The former, appended below,
! evaluates the inverse function  t = t(y)  of  y = t ln t.


REAL, INTENT(IN)                         :: x
DIMENSION cnum(8),cden(8)
DATA cnum/4.120843185847770,85.68982062831317,243.175243524421,  &
    -261.7218583856145,-922.2613728801522,-517.6383498023218,  &
    -77.41064071332953,-2.208843997216182/,  &
    cden/1.,45.64677187585908,377.8372484823942,951.323597679706,  &
    846.0755362020782,262.3083470269460,24.43519662506312, .4097792921092615/

! The constants in the statement below are  exp(1.)  and  .5*log(8.).

m0=2.71828*t((LOG(r1mach(2))-1.03972)/2.71828)
xi=AINT(x)
IF(x-xi > .5) xi=xi+1.
m=IFIX(xi)-1

! Computation of log gamma on the standard interval (1/2,3/2]

xe=x-REAL(m)
snum=cnum(1)
sden=cden(1)
DO  k=2,8
  snum=xe*snum+cnum(k)
  sden=xe*sden+cden(k)
END DO
salga=(xe-1.)*snum/sden

! Computation of log gamma on (0,1/2]

IF(m == -1) THEN
  salga=salga-LOG(x)
  RETURN
ELSE IF(m == 0) THEN
  RETURN
ELSE
  
! Computation of log gamma on (3/2,5/2]
  
  p=xe
  IF(m == 1) THEN
    salga=salga+LOG(p)
    RETURN
  ELSE
    
! Computation of log gamma for arguments larger than 5/2
    
    mm1=m-1
    
! The else-clause in the next statement is designed to avoid possible
! overflow in the computation of  p  in the if-clause, at the expense
! of computing many logarithms.
    
    IF(m < m0) THEN
      DO  k=1,mm1
        p=(xe+REAL(k))*p
      END DO
      salga=salga+LOG(p)
      RETURN
    ELSE
      salga=salga+LOG(xe)
      DO  k=1,mm1
        salga=salga+LOG(xe+REAL(k))
      END DO
      RETURN
    END IF
  END IF
END IF
END FUNCTION salga

FUNCTION sgamma(x,ierr)

! This evaluates the gamma function for real positive  x, using the
! function subroutines  salga  and  r1mach. In case of overflow, the
! routine returns the largest machine-representable number and the
! error flag  ierr=2.

almach=LOG(r1mach(2))
ierr=0
t=salga(x)
IF(t >= almach) THEN
  ierr=2
  gamma=r1mach(2)
  RETURN
ELSE
  gamma=EXP(t)
  RETURN
END IF
END FUNCTION sgamma

FUNCTION t(y)

! This evaluates the inverse function  t = t(y)  of y = t ln t  for
! nonnegative  y  to an accuracy of about one percent. For the
! approximation used, see pp. 51-52 in W. Gautschi,Computational
! aspects of three-term recurrence relations'', SIAM Rev. 9, 1967,
! 24-82.

IF(y <= 10.) THEN
  p=.000057941*y-.00176148
  p=y*p+.0208645
  p=y*p-.129013
  p=y*p+.85777
  t=y*p+1.0125
ELSE
  z=LOG(y)-.775
  p=(.775-LOG(z))/(1.+z)
  p=1./(1.+p)
  t=y*p/z
END IF
RETURN
END FUNCTION t

