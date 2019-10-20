!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:34:54

SUBROUTINE schri(n,iopt,a,b,x,y,hr,hi,alpha,beta,ierr)

! This subroutine implements the Christoffel or generalized Christoffel
! theorem. In all cases except  iopt=7, it uses nonlinear recurrence
! algorithms described in W. Gautschi,An algorithmic implementation
! of the generalized Christoffel theorem'', Numerical Integration
! (G. Haemmerlin, ed.), Birkhaeuser, Basel, 1982, pp. 89-106. The case
! iopt=7  incorporates a QR step with shift  x  in the manner of
! J. Kautsky and G.H. Golub, On the calculation of Jacobi matrices'',
! Linear Algebra Appl. 52/53, 1983, 439-455, using the algorithm of
! Eq. (67.11) on p. 567 in J.H. Wilkinson,The Algebraic Eigenvalue
! Problem'', Clarendon Press, Oxford, 1965. Given the recursion
! coefficients  a(k),b(k), k=0,1,...,n, for the (monic) orthogonal
! polynomials with respect to some measure  dlambda(t), it generates
! the recursion coefficients  alpha(k),beta(k), k=0,1,...,n-1, for the
! measure

!              (t-x)dlambda(t)               if  iopt=1
!              [(t-x)**2+y**2]dlambda(t)     if  iopt=2
!              (t**2+y**2)dlambda(t) with    if  iopt=3
!                dlambda(t) and supp(dlambda)
!                symmetric  with respect to
!                the origin
!              dlambda(t)/(t-x)              if  iopt=4
!              dlambda(t)/[(t-x)**2+y**2]    if  iopt=5
!              dlambda(t)/(t**2+y**2) with   if  iopt=6
!                dlambda(t) and supp(dlambda)
!                symmetric with respect to
!                the origin
!              [(t-x)**2]dlambda(t)          if  iopt=7


!      Input:   n  - - - the number of recurrence coefficients
!                        desired; type integer
!               iopt - - an integer selecting the desired weight
!                        distribution
!               a,b  - - arrays of dimension  n+1  containing the
!                        recursion coefficients a(k-1),b(k-1),k=1,2,
!                        ...,n+1, of the polynomials orthogonal with
!                        respect to the given measure  dlambda(t)
!               x,y  - - real parameters defining the linear and
!                        quadratic factors, or divisors, of  dlambda(t)
!               hr,hi  - the real and imaginary part, respectively, of
!                        the integral of dlambda(t)/(z-t), where z=x+iy;
!                        the parameter  hr  is used only if  iopt=4 or
!                        5, the parameter  hi  only if  iopt=5 or 6

!      Output:  alpha,beta - - arrays of dimension  n  containing the
!                         desired recursion coefficients  alpha(k-1),
!                         beta(k-1), k=1,2,...,n

! It is assumed that  n  is larger than or equal to 2. Otherwise, the
! routine exits immediately with the error flag  ierr  set equal to 1.
! If  iopt  is not between 1 and 7, the routine exits with  ierr=2.

! The routine uses the function subroutine  r1mach  to evaluate the
! constant  eps , which is used only if  iopt=7.


INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: iopt
REAL, INTENT(IN)                         :: a(*)
REAL, INTENT(IN)                         :: b(*)
REAL, INTENT(IN)                         :: x
REAL, INTENT(IN)                         :: y
REAL, INTENT(IN)                         :: hr
REAL, INTENT(IN)                         :: hi
REAL, INTENT(OUT)                        :: alpha(n)
REAL, INTENT(OUT)                        :: beta(n)
INTEGER, INTENT(OUT)                     :: ierr

REAL :: gamma

! The arrays  a,b  are assumed to have dimension  n+1.

eps=5.*r1mach(3)

! The quantity  eps  is a constant slightly larger than the machine
! precision.

ierr=0
IF(n < 2) THEN
  ierr=1
  RETURN
END IF

! What follows implements Eq. (3.7) of W. Gautschi, op. cit.

IF (iopt == 1) THEN
  e=0.
  DO  k=1,n
    q=a(k)-e-x
    beta(k)=q*e
    e=b(k+1)/q
    alpha(k)=x+q+e
  END DO
  
! Set the first beta-coefficient as discussed in Section 5.1 of the
! companion paper.
  
  beta(1)=b(1)*(a(1)-x)
  RETURN
  
! What follows implements Eq. (4.7) of W. Gautschi, op. cit.
  
ELSE IF(iopt == 2) THEN
  s=x-a(1)
  t=y
  eio=0.
  DO  k=1,n
    d=s*s+t*t
    er=-b(k+1)*s/d
    ei=b(k+1)*t/d
    s=x+er-a(k+1)
    t=y+ei
    alpha(k)=x+t*er/ei-s*ei/t
    beta(k)=t*eio*(1.+(er/ei)**2)
    eio=ei
  END DO
  
! Set the first beta-coefficient.
  
  beta(1)=b(1)*(b(2)+(a(1)-x)**2+y*y)
  RETURN
  
! What follows implements Eq. (4.8) of W. Gautschi, op. cit.
  
ELSE IF(iopt == 3) THEN
  t=y
  eio=0.
  DO  k=1,n
    ei=b(k+1)/t
    t=y+ei
    alpha(k)=0.
    beta(k)=t*eio
    eio=ei
  END DO
  
! Set the first beta-coefficient.
  
  beta(1)=b(1)*(b(2)+y*y)
  RETURN
  
! What follows implements Eqs. (5.1),(5.2) of W. Gautschi, op. cit.
  
ELSE IF(iopt == 4) THEN
  alpha(1)=x-b(1)/hr
  beta(1)=-hr
  q=-b(1)/hr
  DO  k=2,n
    e=a(k-1)-x-q
    beta(k)=q*e
    q=b(k)/e
    alpha(k)=q+e+x
  END DO
  RETURN
  
! What follows implements Eq. (5.8) of W. Gautschi, op. cit.
  
ELSE IF(iopt == 5) THEN
  nm1=n-1
  d=hr*hr+hi*hi
  eroo=a(1)-x+b(1)*hr/d
  eioo=-b(1)*hi/d-y
  alpha(1)=x+hr*y/hi
  beta(1)=-hi/y
  alpha(2)=x-b(1)*hi*eroo/(d*eioo)+hr*eioo/hi
  beta(2)=y*eioo*(1.+(hr/hi)**2)
  IF(n == 2) RETURN
  so=b(2)/(eroo**2+eioo**2)
  ero=a(2)-x-so*eroo
  eio=so*eioo-y
  alpha(3)=x+eroo*eio/eioo+so*eioo*ero/eio
  beta(3)=-b(1)*hi*eio*(1.+(eroo/eioo)**2)/d
  IF(n == 3) RETURN
  DO  k=3,nm1
    s=b(k)/(ero**2+eio**2)
    er=a(k)-x-s*ero
    ei=s*eio-y
    alpha(k+1)=x+ero*ei/eio+s*eio*er/ei
    beta(k+1)=so*eioo*ei*(1.+(ero/eio)**2)
    eroo=ero
    eioo=eio
    ero=er
    eio=ei
    so=s
  END DO
  RETURN
  
! What follows implements Eq. (5.9) of W. Gautschi, op. cit.
  
ELSE IF(iopt == 6) THEN
  nm1=n-1
  eoo=-b(1)/hi-y
  eo=b(2)/eoo-y
  alpha(1)=0.
  beta(1)=-hi/y
  alpha(2)=0.
  beta(2)=y*eoo
  IF(n == 2) RETURN
  alpha(3)=0.
  beta(3)=-b(1)*eo/hi
  IF(n == 3) RETURN
  DO  k=3,nm1
    e=b(k)/eo-y
    beta(k+1)=b(k-1)*e/eoo
    alpha(k+1)=0.
    eoo=eo
    eo=e
  END DO
  RETURN
  
! What follows implements a QR step with shift  x.
  
ELSE IF(iopt == 7) THEN
  u=0.
  c=1.
  c0=0.
  DO  k=1,n
    gamma=a(k)-x-u
    cm1=c0
    c0=c
    IF(ABS(c0) > eps) THEN
      p2=(gamma**2)/c0
    ELSE
      p2=cm1*b(k)
    END IF
    IF(k > 1) beta(k)=s*(p2+b(k+1))
    s=b(k+1)/(p2+b(k+1))
    c=p2/(p2+b(k+1))
    u=s*(gamma+a(k+1)-x)
    alpha(k)=gamma+u+x
  END DO
  beta(1)=b(1)*(b(2)+(x-a(1))**2)
  RETURN
ELSE
  ierr=2
  RETURN
END IF
END SUBROUTINE schri

