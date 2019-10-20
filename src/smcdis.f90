!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:34:09

SUBROUTINE smcdis(n,ncapm,mc,mp,xp,yp,squad,eps,iq,idelta,irout,  &
    finl,finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie,  &
    be,x,w,xm,wm,p0,p1,p2)

! This is a multiple-component discretization procedure as described in
! Section 4.3 of the companion paper. It generates to a relative
! accuracy of  eps  the recursion coefficients  alpha(k), beta(k),
! k=0,1,...,n-1, for the polynomials orthogonal with respect to a
! weight distribution consisting of the sum of  mc  continuous
! components and a discrete component with  mp  points. The continuous
! part of the spectrum is made up of  mc  weight functions, each
! supported on its own interval. These intervals may or may not be
! disjoint. The discretization of the inner product on the i-th
! interval is furnished either by a user-supplied subroutine  quad,
! or by the general-purpose subroutine  qgp  provided in this package,
! depending on whether  iq  is equal, or not equal, to  1, respectively.
! The user-supplied routine must have the form  quad(n,x,w,i,ierr)  and
! is assumed to supply the abscissas  x(k)  and weights  w(k), k=1,2,
! ...,n, to be used in approximating the i-th inner product

!               integral of p(x)*q(x)*wf(x,i)dx

! by the

!       sum over k from 1 to n of w(k)*p(x(k))*q(x(k)),

!                                        i=1,2,...,mc.

! The desired recurrence coefficients are then approximated by the
! recursion coefficients of the discrete orthogonal polynomials
! belonging to the discretized inner product, which in turn are
! computed by either the Stieltjes procedure or the Lanczos algorithm
! according as  irout  is equal to, or not equal to  1, respectively.
! Two error flags  ierr,ie  are provided which signal the occurrence
! of an error condition in the quadrature process, or in the routine
! sti  or  lancz  (whichever is used), respectively. The point spectrum
! is given through its abscissas  xp  and jumps  yp.

! If the quadrature routine  quad  has polynomial degree of exactness
! at least  id(n)  for each i, and if  id(n)/n = idelta + O(1/n)  as
! n  goes to infinity, then the procedure is designed to converge after
! one iteration, provided  idelta  is set with the appropriate
! integer. Normally,  idelta=1 (for interpolatory rules) or  idelta=2
! (for Gaussian rules). The default value is  idelta=1.

!    Input:  n    - - the number of recursion coefficients desired;
!                     type integer
!            ncapm  - a discretization parameter indicating an upper
!                     limit of the fineness of the discretization;
!                     ncapm=500  will usually be satisfactory; type
!                     integer
!            mc  - -  the number of disjoint intervals in the
!                     continuous part of the spectrum; type integer
!            mp  - -  the number of points in the discrete part of
!                     the spectrum; type integer. If there is no
!                     point spectrum, set  mp=0.
!            xp  - -  an array of dimension  mp  containing the
!                     abscissas of the point spectrum
!            yp  - -  an array of dimension  mp  containing the jumps
!                     of the point spectrum
!            quad  -  a subroutine determining the discretization of
!                     the inner product on each component interval,
!                     or a dummy routine if  iq  is not equal to  1
!                     (see below)
!            eps  - - the desired relative accuracy of the nonzero
!                     recursion coefficients; type real
!            iq   - - an integer selecting a user-supplied quadrature
!                     routine  quad  if  iq=1  or the ORTHPOL routine
!                     qgp  otherwise
!            idelta - a nonzero integer, typically  1  or  2, inducing
!                     fast convergence in the case of special quadrature
!                     routines
!            irout  - an integer selecting the routine for generating
!                     the recursion coefficients from the discrete
!                     inner product. Specifically,  irout=1  selects the
!                     routine  sti, whereas any other value selects the
!                     routine  lancz

! The logical variables  finl,finr, the arrays  endl,endr  of
! dimension  mc, and the arrays  xfer,wfer  of dimension  ncapm  are
! input variables to the subroutine  qgp  and are used (and hence need
! to be properly dimensioned) only if  iq  is not equal to  1.

!    Output:  alpha,beta - arrays of dimension n, holding as k-th
!                     element  alpha(k-1), beta(k-1), k=1,2,...,n,
!                     respectively
!             ncap  - an integer indicating the fineness of the
!                     discretization that yields convergence within
!                     the eps-tolerance
!             kount - the number of iterations used
!             ierr  - an error flag, equal to  0  on normal return,
!                     equal to  -1  if  n  is not in the proper range,
!                     equal to  i  if there is an error condition in
!                     the discretization of the i-th interval,
!                     and equal to  ncapm  if the discretized
!                     Stieltjes procedure does not converge within the
!                     discretization resolution specified by  ncapm
!             ie - -  an error flag inherited from the routine  sti
!                     or  lancz  (whichever is used)

! The array  be  of dimension  n, the arrays  x,w  of dimension  ncapm,
! and the arrays  xm,wm,p0,p1,p2  of dimension mc*ncapm+mp  are used
! for working space. The routine calls upon the subroutine  sti  or
! lancz, depending on the choice of  irout.


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
INTEGER, INTENT(IN OUT)                  :: irout
LOGICAL, INTENT(IN OUT)                  :: finl
LOGICAL, INTENT(IN OUT)                  :: finr
REAL, INTENT(IN OUT)                     :: endl(mc)
REAL, INTENT(IN OUT)                     :: endr(mc)
REAL, INTENT(IN OUT)                     :: xfer(ncapm)
REAL, INTENT(IN OUT)                     :: wfer(ncapm)
REAL, INTENT(IN OUT)                     :: alpha(n)
REAL, INTENT(OUT)                        :: beta(n)
INTEGER, INTENT(OUT)                     :: ncap
INTEGER, INTENT(OUT)                     :: kount
INTEGER, INTENT(OUT)                     :: ierr
INTEGER, INTENT(IN OUT)                  :: ie
REAL, INTENT(OUT)                        :: be(n)
REAL, INTENT(IN)                         :: x(ncapm)
REAL, INTENT(IN)                         :: w(ncapm)
REAL, INTENT(OUT)                        :: xm(*)
REAL, INTENT(OUT)                        :: wm(*)
REAL, INTENT(IN OUT)                     :: p0(*)
REAL, INTENT(IN OUT)                     :: p1(*)
REAL, INTENT(IN OUT)                     :: p2(*)



! The arrays  xp,yp  are assumed to have dimension  mp  if mp > 0, and
! the arrays  xm,wm,p0,p1,p2  dimension  mc*ncapm+mp.

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
ncap=(2*n-1)/idelta
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

! Discretization of the inner product

mtncap=mc*ncap
DO  i=1,mc
  im1tn=(i-1)*ncap
  IF(iq == 1) THEN
    CALL squad(ncap,x,w,i,ierr)
  ELSE
    CALL sqgp(ncap,x,w,i,ierr,mc,finl,finr,endl,endr,xfer, wfer)
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

! Computation of the desired recursion coefficients

IF(irout == 1) THEN
  CALL ssti(n,mtncap+mp,xm,wm,alpha,beta,ie,p0,p1,p2)
ELSE
  CALL slancz(n,mtncap+mp,xm,wm,alpha,beta,ie,p0,p1)
END IF

! In the following statement, the absolute value of the beta's is
! used to guard against failure in cases where the routine is applied
! to variable-sign weight functions and hence the positivity of the
! beta's is not guaranteed.

DO  k=1,n
  IF(ABS(beta(k)-be(k)) > eps*ABS(beta(k))) GO TO 20
END DO
RETURN
END SUBROUTINE smcdis

