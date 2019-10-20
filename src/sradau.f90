!

! Code converted using TO_F90 by Alan Miller
! Date: 2019-10-20  Time: 15:33:24

subroutine sradau(n,alpha,beta,end,zero,weight,ierr,e,a,b)

! Given  n  and a measure  dlambda, this routine generates the
! (n+1)-point Gauss-Radau quadrature formula

!   integral over supp(dlambda) of f(t)dlambda(t)

!     = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).

! The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
! =w(k), k=0,1,2,...,n. The user has to supply the recursion
! coefficients  alpha(k), beta(k), k=0,1,2,...,n, for the measure
! dlambda. The nodes and weights are computed as eigenvalues and
! in terms of the first component of the respective normalized
! eigenvectors of a slightly modified Jacobi matrix of order  n+1.
! To do this, the routine calls upon the subroutine  gauss. It also
! uses the function subroutine  r1mach.

!    Input:  n - -  the number of interior points in the Gauss-Radau
!                   formula; type integer
!            alpha,beta - arrays of dimension  n+1  to be supplied with
!                   the recursion coefficients  alpha(k-1), beta(k-1),
!                   k=1,2,...,n+1; the coefficient  alpha(n+1)  is not
!                   used by the routine
!            end -  the prescribed endpoint  x(0)  of the Gauss-Radau
!                   formula; type real

!    Output: zero - array of dimension  n+1  containing the nodes (in
!                   increasing order)  zero(k)=x(k), k=0,1,2,...,n
!            weight-array of dimension  n+1  containing the weights
!                   weight(k)=w(k), k=0,1,2,...,n
!            ierr - an error flag inherited from the routine  gauss

! The arrays  e,a,b  are needed for working space.

    integer, intent(in)                      :: n
    real, intent(in)                         :: alpha(n+1)
    real, intent(in)                         :: beta(n+1)
    real, intent(in)                         :: end
    real, intent(in out)                     :: zero(*)
    real, intent(in out)                     :: weight(*)
    integer, intent(in out)                  :: ierr
    real, intent(in out)                     :: e(*)
    real, intent(out)                        :: a(*)
    real, intent(out)                        :: b(*)


! The arrays  alpha,beta,zero,weight,e,a,b  are assumed to have
! dimension  n+1.
    real :: epsma, p0, p1, pm1
    integer :: np1, k

    ! epsma is the machine single precision.
    epsma = r1mach(3)

    do k = 1, n+1
        a(k) = alpha(k)
        b(k) = beta(k)
    end do
    p0 = 0.
    p1 = 1.
    do k = 1, n
        pm1 = p0
        p0 = p1
        p1 = (end - a(k))*p0 - b(k)*pm1
    end do
    a(n+1) = end - b(n+1)*p0/p1
    call sgauss(n+1,a,b,epsma,zero,weight,ierr,e)
end subroutine sradau

