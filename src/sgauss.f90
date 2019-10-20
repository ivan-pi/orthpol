!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:34:46

SUBROUTINE sgauss(n,alpha,beta,eps,zero,weight,ierr,e)

! Given  n  and a measure  dlambda, this routine generates the n-point
! Gaussian quadrature formula

!     integral over supp(dlambda) of f(x)dlambda(x)

!        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).

! The nodes are returned as  zero(k)=x(k) and the weights as
! weight(k)=w(k), k=1,2,...,n. The user has to supply the recursion
! coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure
! dlambda. The routine computes the nodes as eigenvalues, and the
! weights in term of the first component of the respective normalized
! eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
! It uses a translation and adaptation of the algol procedure  imtql2,
! Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
! by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
! Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
! routine  imtql2.

!        Input:  n - - the number of points in the Gaussian quadrature
!                      formula; type integer
!                alpha,beta - - arrays of dimension  n  to be filled
!                      with the values of  alpha(k-1), beta(k-1), k=1,2,
!                      ...,n
!                eps - the relative accuracy desired in the nodes
!                      and weights

!        Output: zero- array of dimension  n  containing the Gaussian
!                      nodes (in increasing order)  zero(k)=x(k), k=1,2,
!                      ...,n
!                weight - array of dimension  n  containing the
!                      Gaussian weights  weight(k)=w(k), k=1,2,...,n
!                ierr- an error flag equal to  0  on normal return,
!                      equal to  i  if the QR algorithm does not
!                      converge within 30 iterations on evaluating the
!                      i-th eigenvalue, equal to  -1  if  n  is not in
!                      range, and equal to  -2  if one of the beta's is
!                      negative.

! The array  e  is needed for working space.



INTEGER, INTENT(IN)                      :: n
REAL, INTENT(IN)                         :: alpha(n)
REAL, INTENT(IN)                         :: beta(n)
REAL, INTENT(IN)                         :: eps
REAL, INTENT(OUT)                        :: zero(n)
REAL, INTENT(OUT)                        :: weight(n)
INTEGER, INTENT(OUT)                     :: ierr
REAL, INTENT(OUT)                        :: e(n)

IF(n < 1) THEN
  ierr=-1
  RETURN
END IF
ierr=0
zero(1)=alpha(1)
IF(beta(1) < 0.) THEN
  ierr=-2
  RETURN
END IF
weight(1)=beta(1)
IF (n == 1) RETURN
weight(1)=1.
e(n)=0.
DO  k=2,n
  zero(k)=alpha(k)
  IF(beta(k) < 0.) THEN
    ierr=-2
    RETURN
  END IF
  e(k-1)=SQRT(beta(k))
  weight(k)=0.
END DO
DO  l=1,n
  j=0
  
! Look for a small subdiagonal element.
  
  105   DO  m=l,n
    IF(m == n) EXIT
    IF(ABS(e(m)) <= eps*(ABS(zero(m))+ABS(zero(m+1)))) EXIT
  END DO
  120   p=zero(l)
  IF(m == l) CYCLE
  IF(j == 30) GO TO 400
  j=j+1
  
! Form shift.
  
  g=(zero(l+1)-p)/(2.*e(l))
  r=SQRT(g*g+1.)
  g=zero(m)-p+e(l)/(g+SIGN(r,g))
  s=1.
  c=1.
  p=0.
  mml=m-l
  
! For i=m-1 step -1 until l do ...
  
  DO  ii=1,mml
    i=m-ii
    f=s*e(i)
    b=c*e(i)
    IF(ABS(f) < ABS(g)) GO TO 150
    c=g/f
    r=SQRT(c*c+1.)
    e(i+1)=f*r
    s=1./r
    c=c*s
    GO TO 160
    150     s=f/g
    r=SQRT(s*s+1.)
    e(i+1)=g*r
    c=1./r
    s=s*c
    160     g=zero(i+1)-p
    r=(zero(i)-g)*s +2.*c*b
    p=s*r
    zero(i+1)=g+p
    g=c*r-b
    
! Form first component of vector.
    
    f=weight(i+1)
    weight(i+1)=s*weight(i)+c*f
    weight(i)=c*weight(i)-s*f
  END DO
  zero(l)=zero(l)-p
  e(l)=g
  e(m)=0.
  GO TO 105
END DO

! Order eigenvalues and eigenvectors.

DO  ii=2,n
  i=ii-1
  k=i
  p=zero(i)
  DO  j=ii,n
    IF(zero(j) >= p) CYCLE
    k=j
    p=zero(j)
  END DO
  IF(k == i) CYCLE
  zero(k)=zero(i)
  zero(i)=p
  p=weight(i)
  weight(i)=weight(k)
  weight(k)=p
END DO
DO  k=1,n
  weight(k)=beta(1)*weight(k)*weight(k)
END DO
RETURN

! Set error - no convergence to an eigenvalue after 30 iterations.

400 ierr=l
RETURN
END SUBROUTINE sgauss

