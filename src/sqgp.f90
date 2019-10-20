!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:33:45

SUBROUTINE sqgp(n,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,wfer)

! This is a general-purpose discretization routine that can be used
! as an alternative to the routine  quad  in the multiple-component
! discretization procedure  mcdis. It takes no account of the special
! nature of the weight function involved and hence may result in slow
! convergence of the discretization procedure. This routine, therefore,
! should be used only as a last resort, when no better, more natural
! discretization can be found.

! It is assumed that there are  mc.ge.1  disjoint component intervals.
! The discretization is effected by the Fejer quadrature rule,
! suitably transformed to the respective interval. An interval that
! extends to minus infinity has to be indexed by  1; one that extends
! to plus infinity has to be indexed by  mc.

! The output variable  ierr  is given the value  0. Additional input
! parameters and working space used by this routine are as follows:

!          mc      - the number of component intervals; type integer
!          finl    - a logical variable to be set .true. if the
!                    extreme left interval is finite, and .false.
!                    otherwise
!          finr    - a logical variable to be set .true. if the
!                    extreme right interval is finite, and .false.
!                    otherwise
!          endl    - an array of dimension  mc  containing the left
!                    endpoints of the component intervals; if the
!                    first of these extends to minus infinity,  endl(1)
!                    can be set to an arbitrary value
!          endr    - an array of dimension  mc  containing the right
!                    endpoints of the component intervals; if the
!                    last of these extends to plus infinity,  endr(mc)
!                    can be set to an arbitrary value
!          xfer,wfer-working arrays holding the Fejer nodes and
!                    weights, respectively, for the interval [-1,1].

! The user has to supply the routine

!                     function wf(x,i),

! which evaluates the weight function at the point  x  on the i-th
! component interval. The routine also uses the subroutines  fejer,
! symtr  and  tr, which are appended.


INTEGER, INTENT(IN)                      :: n
REAL, INTENT(OUT)                        :: x(n)
REAL, INTENT(OUT)                        :: w(n)
INTEGER, INTENT(IN)                      :: i
INTEGER, INTENT(OUT)                     :: ierr
INTEGER, INTENT(IN OUT)                  :: mc
LOGICAL, INTENT(IN OUT)                  :: finl
LOGICAL, INTENT(IN OUT)                  :: finr
REAL, INTENT(IN)                         :: endl(mc)
REAL, INTENT(IN)                         :: endr(mc)
REAL, INTENT(IN)                         :: xfer(*)
REAL, INTENT(IN)                         :: wfer(*)



! The arrays  xfer,wfer  are dimensioned in the routine  mcdis.

ierr=0
IF(i == 1) CALL sfejer(n,xfer,wfer)
IF(i > 1 .AND. i < mc) GO TO 60
IF(mc == 1) THEN
  IF(finl.AND.finr) GO TO 60
  IF(finl) GO TO 20
  IF(finr) GO TO 40
  DO  k=1,n
    CALL ssymtr(xfer(k),phi,phi1)
    x(k)=phi
    w(k)=wfer(k)*wf(phi,i)*phi1
  END DO
  RETURN
ELSE
  IF((i == 1.AND.finl).OR.(i == mc.AND.finr)) GO TO 60
  IF(i == 1) GO TO 40
END IF
20 DO  k=1,n
  CALL str(xfer(k),phi,phi1)
  x(k)=endl(mc)+phi
  w(k)=wfer(k)*wf(x(k),mc)*phi1
END DO
RETURN
40 DO  k=1,n
  CALL str(-xfer(k),phi,phi1)
  x(k)=endr(1)-phi
  w(k)=wfer(k)*wf(x(k),1)*phi1
END DO
RETURN
60 DO  k=1,n
  x(k)=.5*((endr(i)-endl(i))*xfer(k)+endr(i)+endl(i))
  w(k)=.5*(endr(i)-endl(i))*wfer(k)*wf(x(k),i)
END DO
RETURN
END SUBROUTINE sqgp

SUBROUTINE ssymtr(t,phi,phi1)

! This implements a particular transformation  x=phi(t)  mapping
! the t-interval [-1,1] to the x-interval [-oo,oo].

!        input:   t
!        output:  phi=phi(t)
!                 phi1=derivative of phi(t)

t2=t*t
phi=t/(1.-t2)
phi1=(t2+1.)/(t2-1.)**2
RETURN
END SUBROUTINE ssymtr

SUBROUTINE str(t,phi,phi1)

! This implements a particular transformation  x=phi(t)  mapping
! the t-interval [-1,1] to the x-interval [0,oo].

!         input:   t
!         output:  phi=phi(t)
!                  phi1=derivative of phi(t)

phi=(1.+t)/(1.-t)
phi1=2./(t-1.)**2
RETURN
END SUBROUTINE str

SUBROUTINE sfejer(n,x,w)

! This routine generates the n-point Fejer quadrature rule.

!         input:   n   - the number of quadrature nodes
!         output:  x,w - arrays of dimension  n  holding the quadrature
!                        nodes and weights, respectively; the nodes
!                        are ordered increasingly



INTEGER, INTENT(IN)                      :: n
REAL, INTENT(OUT)                        :: x(n)
REAL, INTENT(OUT)                        :: w(n)

pi=4.*ATAN(1.)
nh=n/2
np1h=(n+1)/2
fn=REAL(n)
DO  k=1,nh
  x(n+1-k)=COS(.5*REAL(2*k-1)*pi/fn)
  x(k)=-x(n+1-k)
END DO
IF(2*nh /= n) x(np1h)=0.
DO  k=1,np1h
  c1=1.
  c0=2.*x(k)*x(k)-1.
  t=2.*c0
  sum=c0/3.
  DO  m=2,nh
    c2=c1
    c1=c0
    c0=t*c1-c2
    sum=sum+c0/REAL(4*m*m-1)
  END DO
  w(k)=2.*(1.-2.*sum)/fn
  w(n+1-k)=w(k)
END DO
RETURN
END SUBROUTINE sfejer

