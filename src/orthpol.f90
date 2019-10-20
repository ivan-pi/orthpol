module orthpol

    implicit none
    private

    public :: sp, dp
    public :: r1mach, d1mach
    public :: recur
    public :: cheb
    public :: sti, lancz, qgp, mcdis, mccheb
    public :: chri, knum, kern, gchri
    public :: nu0jac, nu0lag, nu0her 
    public :: gauss, radau, lob

    public :: schri, dchri, sgchri, dgchri
    public :: srecur, drecur
    public :: scheb, dcheb

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)

    real(sp), parameter :: pi = 4.0_sp*atan(1.0_sp)

    interface
        real(sp) function r1mach(i)
            import sp
            integer, intent(in) :: i
        end function
        real(dp) function d1mach(i)
            import dp
            integer, intent(in) :: i
        end function
    end interface

    !--------------------------------------------------------
    ! 2. Classical weight functions
    !--------------------------------------------------------


    interface recur
        subroutine srecur(n,ipoly,al,be,a,b,ierr)
            import sp
            integer, intent(in) :: n
            integer, intent(in) :: ipoly
            real(sp), intent(in) :: al, be
            real(sp), intent(out) :: a(n), b(n)
            integer, intent(out) :: ierr
        end subroutine
        subroutine drecur(n,ipoly,dal,dbe,da,db,iderr)
            import dp
            integer, intent(in) :: n
            integer, intent(in) :: ipoly
            real(dp), intent(in) :: dal, dbe
            real(dp), intent(out) :: da(n), db(n)
            integer, intent(out) :: iderr
        end subroutine
    end interface

    !--------------------------------------------------------
    ! 3. Moment-related methods
    !--------------------------------------------------------
        
    ! Modified Chebyshev algorithm
    interface cheb
        subroutine scheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
            import sp
            integer, intent(in) :: n
            real(sp), intent(in) :: a(*), b(*) ! 2*n-1
            real(sp), intent(in) :: fnu(*) ! 2*n-1
            real(sp), intent(out) :: alpha(n), beta(n), s(n)
            integer, intent(out) :: ierr
            real(sp), intent(inout) :: s0(*), s1(*), s2(*) ! 2*n
        end subroutine
        subroutine dcheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
            import dp
            integer, intent(in) :: n
            real(dp), intent(in) :: a(*), b(*) ! 2*n-1
            real(dp), intent(in) :: fnu(*) ! 2*n-1
            real(dp), intent(out) :: alpha(n), beta(n), s(n)
            integer, intent(out) :: ierr
            real(dp), intent(inout) :: s0(*), s1(*), s2(*) ! 2*n
        end subroutine
    end interface


    !--------------------------------------------------------
    ! 4. Stieltjes, orthogonal reduction, and discretization procedures
    !--------------------------------------------------------

    ! Stieltjes's procedure
    interface sti
        subroutine ssti(n,ncap,x,w,alpha,beta,ierr,p0,p1,p2)
            import sp
            integer, intent(in) :: n
            integer, intent(in) :: ncap
            real(sp), intent(in) :: x(ncap), w(ncap)
            real(sp), intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: ierr
            real(sp), intent(inout) :: p0(ncap), p1(ncap), p2(ncap)
        end subroutine
        subroutine dsti(n,ncap,x,w,alpha,beta,ierr,p0,p1,p2)
            import dp
            integer, intent(in) :: n
            integer, intent(in) :: ncap
            real(dp), intent(in) :: x(ncap), w(ncap)
            real(dp), intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: ierr
            real(dp), intent(inout) :: p0(ncap), p1(ncap), p2(ncap)
        end subroutine
    end interface

    ! Lanczos algorithm
    interface lancz
        subroutine slancz(n,ncap,x,w,alpha,beta,ierr,p0,p1)
            import sp
            integer, intent(in) :: n
            integer, intent(in) :: ncap
            real(sp), intent(in) :: x(ncap), w(ncap)
            real(sp), intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: ierr
            real(sp), intent(inout) :: p0(ncap), p1(ncap)
        end subroutine
        subroutine dlancz(n,ncap,x,w,alpha,beta,ierr,p0,p1)
            import dp
            integer, intent(in) :: n
            integer, intent(in) :: ncap
            real(dp), intent(in) :: x(ncap), w(ncap)
            real(dp), intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: ierr
            real(dp), intent(inout) :: p0(ncap), p1(ncap)
        end subroutine
    end interface

    interface qgp
        ! General-purpose quadrature routine
        subroutine sqgp(n,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,wfer)
            import sp
            integer, intent(in) :: n
            real(sp), intent(out) :: x(n), w(n)
            integer, intent(in) :: i
            integer, intent(out) :: ierr
            integer, intent(in) :: mc
            logical, intent(in) :: finl, finr
            real(sp), intent(in) :: endl(mc), endr(mc)
            real(sp), intent(inout) :: xfer(*), wfer(*) ! dimensioned by the routine mcdis
        end subroutine
        subroutine dqgp(n,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,wfer)
            import dp
            integer, intent(in) :: n
            real(dp), intent(out) :: x(n), w(n)
            integer, intent(in) :: i
            integer, intent(out) :: ierr
            integer, intent(in) :: mc
            logical, intent(in) :: finl, finr
            real(dp), intent(in) :: endl(mc), endr(mc)
            real(dp), intent(inout) :: xfer(*), wfer(*) ! dimensioned by the routine mcdis
        end subroutine
    end interface

    abstract interface
        subroutine squad_interface(n,x,w,i,ierr)
            import sp
            integer, intent(in) :: n
            real(sp), intent(out) :: x(n), w(n)
            integer, intent(in) :: i
            integer, intent(out) :: ierr
        end subroutine
        subroutine dquad_interface(n,x,w,i,ierr)
            import dp
            integer, intent(in) :: n
            real(dp), intent(out) :: x(n), w(n)
            integer, intent(in) :: i
            integer, intent(out) :: ierr
        end subroutine
    end interface

    ! Multiple-component discretization
    interface mcdis
        subroutine smcdis(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta,irout, &
            finl,finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie, &
            be,x,w,xm,wm,p0,p1,p2)
            integer, intent(in) :: n
            integer, intent(in) :: ncapm
            integer, intent(in) :: mc
            integer, intent(in) :: mp
            real, intent(in) :: xp(mp), yp(mp)
            procedure(squad_interface) :: quad
            real, intent(in) :: eps
            integer, intent(in) :: iq
            integer, intent(in) :: idelta
            integer, intent(in) :: irout
            ! Input variables for qgp if iq /= 1
            logical, intent(in) :: finl, finr
            real, intent(in) :: endl(mc), endr(mc)
            real, intent(inout) :: xfer(ncapm), wfer(ncapm)

            ! Output
            real, intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: ncap
            integer, intent(out) :: kount
            integer, intent(out) :: ierr
            integer, intent(out) :: ie
            ! Work arrays
            real, intent(inout) :: be(n), x(ncapm), w(ncapm)
            real, intent(inout) :: xm(*), wm(*), p0(*), p1(*), p2(*) ! mc*ncapm+mp
        end subroutine
        subroutine dmcdis(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta,irout, &
            finl,finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie, &
            be,x,w,xm,wm,p0,p1,p2)
            integer, intent(in) :: n
            integer, intent(in) :: ncapm
            integer, intent(in) :: mc
            integer, intent(in) :: mp
            double precision, intent(in) :: xp(mp), yp(mp)
            procedure(dquad_interface) :: quad
            double precision, intent(in) :: eps
            integer, intent(in) :: iq
            integer, intent(in) :: idelta
            integer, intent(in) :: irout
            ! Input variables for qgp if iq /= 1
            logical, intent(in) :: finl, finr
            double precision, intent(in) :: endl(mc), endr(mc)
            double precision, intent(inout) :: xfer(ncapm), wfer(ncapm)

            ! Output
            double precision, intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: ncap
            integer, intent(out) :: kount
            integer, intent(out) :: ierr
            integer, intent(out) :: ie
            ! Work arrays
            double precision, intent(inout) :: be(n), x(ncapm), w(ncapm)
            double precision, intent(inout) :: xm(*), wm(*), p0(*), p1(*), p2(*) ! mc*ncapm+mp
        end subroutine
    end interface

    ! Discretized modified Chebyshev algorithm
    interface mccheb
        subroutine smccheb(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta, &
            finl,finr,endl,endr,xfer,wfer,a,b,fnu,alpha,beta,ncap,kount, &
            ierr,be,x,w,xm,wm,s,s0,s1,s2)
            integer, intent(in) :: n
            integer, intent(in) :: ncapm
            integer, intent(in) :: mc
            integer, intent(in) :: mp
            real, intent(in) :: xp(mp), yp(mp)
            procedure(squad_interface) :: quad
            real, intent(in) :: eps
            integer, intent(in) :: iq
            integer, intent(in) :: idelta
     
            ! Input variables for qgp if iq /= 1
            logical, intent(in) :: finl, finr
            real, intent(in) :: endl(mc), endr(mc)
            real, intent(inout) :: xfer(ncapm), wfer(ncapm)
            real, intent(in) :: a(*), b(*) ! 2*n-1
            ! Output
            real, intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: ncap
            integer, intent(out) :: kount
            integer, intent(out) :: ierr
            ! Work arrays
            real, intent(inout) :: fnu(*), s(n),  s0(*), s1(*), s2(*)  ! 2*n
            real, intent(inout) :: be(n), x(ncapm), w(ncapm)
            real, intent(inout) :: xm(*), wm(*) ! mc*ncapm+mp
        end subroutine
        subroutine dmccheb(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta, &
            finl,finr,endl,endr,xfer,wfer,a,b,fnu,alpha,beta,ncap,kount, &
            ierr,be,x,w,xm,wm,s,s0,s1,s2)
            integer, intent(in) :: n
            integer, intent(in) :: ncapm
            integer, intent(in) :: mc
            integer, intent(in) :: mp
            double precision, intent(in) :: xp(mp), yp(mp)
            procedure(dquad_interface) :: quad
            double precision, intent(in) :: eps
            integer, intent(in) :: iq
            integer, intent(in) :: idelta
     
            ! Input variables for qgp if iq /= 1
            logical, intent(in) :: finl, finr
            double precision, intent(in) :: endl(mc), endr(mc)
            double precision, intent(inout) :: xfer(ncapm), wfer(ncapm)
            double precision, intent(in) :: a(*), b(*) ! 2*n-1
            ! Output
            double precision, intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: ncap
            integer, intent(out) :: kount
            integer, intent(out) :: ierr
            ! Work arrays
            double precision, intent(inout) :: fnu(*), s(n),  s0(*), s1(*), s2(*)  ! 2*n
            double precision, intent(inout) :: be(n), x(ncapm), w(ncapm)
            double precision, intent(inout) :: xm(*), wm(*) ! mc*ncapm+mp
        end subroutine
    end interface

    !--------------------------------------------------------
    ! 5. Modification algorithms
    !--------------------------------------------------------

    ! Nonlinear Recurrence Algorithms
    interface chri
        subroutine schri(n,iopt,a,b,x,y,hr,hi,alpha,beta,ierr)
            import sp
            integer, intent(in) :: n
            integer, intent(in) :: iopt
            real(sp), intent(in) :: a(n+1), b(n+1)
            real(sp), intent(in) :: x, y
            real(sp), intent(in) :: hr, hi
            real(sp), intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: ierr
        end subroutine
        subroutine dchri(n,iopt,a,b,x,y,hr,hi,alpha,beta,ierr)
            import dp
            integer, intent(in) :: n
            integer, intent(in) :: iopt
            real(dp), intent(in) :: a(n+1), b(n+1)
            real(dp), intent(in) :: x, y
            real(dp), intent(in) :: hr, hi
            real(dp), intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: ierr
        end subroutine
    end interface


    ! interface
    !     ! Jacobi measure estimate
    !     integer function nu0jac(n,z,eps)
    !         import sp
    !         integer, intent(in) :: n
    !         complex(sp), intent(in) :: z
    !         real(sp), intent(in) :: eps
    !     end function
    !     ! Generalized Laguerre measure estimate
    !     integer function nu0lag(n,z,al,eps)
    !         import sp
    !         integer, intent(in) :: n
    !         complex(sp), intent(in) :: z
    !         real(sp), intent(in) :: al
    !         real(sp), intent(in) :: eps            
    !     end function
    !     ! Hermite measure estimate
    !     integer function nu0her(n,z,eps)
    !         import sp
    !         integer, intent(in) :: n
    !         complex(sp), intent(in) :: z
    !         real(sp), intent(in) :: eps
    !     end function
    ! end interface

    interface knum
        subroutine sknum(n,nu0,numax,z,eps,a,b,rho,nu,ierr,rold)
            import sp
            integer, intent(in) :: n
            integer, intent(in) :: nu0
            integer, intent(in) :: numax
            complex(sp), intent(in) :: z
            real(sp), intent(in) :: eps
            real(sp), intent(in) :: a(numax), b(numax)
            complex(sp), intent(out) :: rho(n+1)
            integer, intent(out) :: nu
            integer, intent(out) :: ierr
            complex(sp), intent(inout) :: rold(n+1)
        end subroutine
        subroutine dknum(n,nu0,numax,z,eps,a,b,rho,nu,ierr,rold)
            import dp
            integer, intent(in) :: n
            integer, intent(in) :: nu0
            integer, intent(in) :: numax
            complex(kind=dp), intent(in) :: z
            real(dp), intent(in) :: eps
            real(dp), intent(in) :: a(numax), b(numax)
            complex(kind=dp), intent(out) :: rho(n+1)
            integer, intent(out) :: nu
            integer, intent(out) :: ierr
            complex(kind=dp), intent(inout) :: rold(n+1)
        end subroutine
    end interface

    interface kern
        subroutine skern(n,nu0,numax,z,eps,a,b,ker,nu,ierr,rold)
            import sp
            integer, intent(in) :: n
            integer, intent(in) :: nu0
            integer, intent(in) :: numax
            complex(sp), intent(in) :: z
            real(sp), intent(in) :: eps
            real(sp), intent(in) :: a(numax), b(numax)
            complex(sp), intent(out) :: ker(n+1)
            integer, intent(out) :: nu
            integer, intent(out) :: ierr
            complex(sp), intent(inout) :: rold(n+1)
        end subroutine
        subroutine dkern(n,nu0,numax,z,eps,a,b,ker,nu,ierr,rold)
            import dp
            integer, intent(in) :: n
            integer, intent(in) :: nu0
            integer, intent(in) :: numax
            complex(kind=dp), intent(in) :: z
            real(dp), intent(in) :: eps
            real(dp), intent(in) :: a(numax), b(numax)
            complex(kind=dp), intent(out) :: ker(n+1)
            integer, intent(out) :: nu
            integer, intent(out) :: ierr
            complex(kind=dp), intent(inout) :: rold(n+1)
        end subroutine
    end interface

    interface gchri
        subroutine sgchri(n,iopt,nu0,numax,eps,a,b,x,y,alpha,beta,nu, &
            ierr,ierrc,fnu,rho,rold,s,s0,s1,s2)
            integer, intent(in) :: n
            integer, intent(in) :: iopt
            integer, intent(in) :: nu0
            integer, intent(in) :: numax
            real, intent(in) :: eps
            real, intent(in) :: a(numax), b(numax)
            real, intent(in) :: x, y
            real, intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: nu
            integer, intent(out) :: ierr
            integer, intent(out) :: ierrc
            real, intent(inout) :: fnu(2*n), s(2*n), s0(2*n), s1(2*n), s2(2*n)
            complex, intent(inout) :: rho(2*n), rold(2*n)
        end subroutine
        subroutine dgchri(n,iopt,nu0,numax,eps,a,b,x,y,alpha,beta,nu, &
            ierr,ierrc,fnu,rho,rold,s,s0,s1,s2)
            import dp
            integer, intent(in) :: n
            integer, intent(in) :: iopt
            integer, intent(in) :: nu0
            integer, intent(in) :: numax
            double precision, intent(in) :: eps
            double precision, intent(in) :: a(numax), b(numax)
            double precision, intent(in) :: x, y
            double precision, intent(out) :: alpha(n), beta(n)
            integer, intent(out) :: nu
            integer, intent(out) :: ierr
            integer, intent(out) :: ierrc
            double precision, intent(inout) :: fnu(2*n), s(2*n), s0(2*n), s1(2*n), s2(2*n)
            complex(kind=dp), intent(inout) :: rho(2*n), rold(2*n)
        end subroutine
    end interface

    !--------------------------------------------------------
    ! 6. Gauss-Type Quadrature Rules
    !--------------------------------------------------------

    ! Gaussian Quadrature
    interface gauss
        ! subroutine sgauss(n,alpha,beta,eps,zero,weight,ierr,e)
        !     integer, intent(in) :: n
        !     real, intent(in) :: alpha(n), beta(n)
        !     real, intent(in) :: eps
        !     real, intent(out) :: zero(n)
        !     real, intent(out) :: weight(n)
        !     integer, intent(out) :: ierr
        !     real, intent(inout) :: e(n)
        ! end subroutine
        module procedure sgauss
        subroutine dgauss(n,alpha,beta,eps,zero,weight,ierr,e)
            integer, intent(in) :: n
            double precision, intent(in) :: alpha(n), beta(n)
            double precision, intent(in) :: eps
            double precision, intent(out) :: zero(n)
            double precision, intent(out) :: weight(n)
            integer, intent(out) :: ierr
            double precision, intent(inout) :: e(n)
        end subroutine
    end interface

    ! Gauss - Radau Quadrature
    interface radau
        ! subroutine sradau(n,alpha,beta,end,zero,weight,ierr,e,a,b)
        !     integer, intent(in) :: n
        !     real, intent(in) :: alpha(n+1), beta(n+1)
        !     real, intent(in) :: end
        !     real, intent(out) :: zero(n+1)
        !     real, intent(out) :: weight(n+1)
        !     integer, intent(out) :: ierr
        !     real, intent(inout) :: e(n+1), a(n+1), b(n+1)
        ! end subroutine
        module procedure sradau
        subroutine dradau(n,alpha,beta,end,zero,weight,ierr,e,a,b)
            integer, intent(in) :: n
            double precision, intent(in) :: alpha(n+1), beta(n+1)
            double precision, intent(in) :: end
            double precision, intent(out) :: zero(n+1)
            double precision, intent(out) :: weight(n+1)
            integer, intent(out) :: ierr
            double precision, intent(inout) :: e(n+1), a(n+1), b(n+1)
        end subroutine
    end interface

    ! Gauss - Lobatto Quadrature
    interface lob
        ! subroutine slob(n,alpha,beta,aleft,right,zero,weight,ierr,e,a,b)
        !     integer, intent(in) :: n
        !     real, intent(in) :: alpha(n+2), beta(n+2)
        !     real, intent(in) :: aleft, right
        !     real, intent(out) :: zero(n+2)
        !     real, intent(out) :: weight(n+2)
        !     integer, intent(out) :: ierr
        !     real, intent(inout) :: e(n+2), a(n+2), b(n+2)
        ! end subroutine
        module procedure slob
        subroutine dlob(n,alpha,beta,aleft,right,zero,weight,ierr,e,a,b)
            integer, intent(in) :: n
            double precision, intent(in) :: alpha(n+2), beta(n+2)
            double precision, intent(in) :: aleft, right
            double precision, intent(out) :: zero(n+2)
            double precision, intent(out) :: weight(n+2)
            integer, intent(out) :: ierr
            double precision, intent(inout) :: e(n+2), a(n+2), b(n+2)
        end subroutine
    end interface

contains

!   Given  n  and a measure  dlambda, this routine generates the
!   (n+1)-point Gauss-Radau quadrature formula
    
!     integral over supp(dlambda) of f(t)dlambda(t)
    
!       = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
    
!   The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
!   =w(k), k=0,1,2,...,n. The user has to supply the recursion
!   coefficients  alpha(k), beta(k), k=0,1,2,...,n, for the measure
!   dlambda. The nodes and weights are computed as eigenvalues and
!   in terms of the first component of the respective normalized
!   eigenvectors of a slightly modified Jacobi matrix of order  n+1.
!   To do this, the routine calls upon the subroutine  gauss. It also
!   uses the function subroutine  r1mach.
    subroutine sradau(n,alpha,beta,end,zero,weight,ierr)
        implicit none

    !    Input:  n - -  the number of interior points in the Gauss-Radau
    !                   formula; type integer
        integer, intent(in)                      :: n
    !            alpha,beta - arrays of dimension  n+1  to be supplied with
    !                   the recursion coefficients  alpha(k-1), beta(k-1),
    !                   k=1,2,...,n+1; the coefficient  alpha(n+1)  is not
    !                   used by the routine
        real(sp), intent(in) :: alpha(n+1), beta(n+1)
    !            end -  the prescribed endpoint  x(0)  of the Gauss-Radau
    !                   formula; type real
        real(sp), intent(in) :: end
        ! Output
    !    Output: zero - array of dimension  n+1  containing the nodes (in
    !                   increasing order)  zero(k)=x(k), k=0,1,2,...,n
        real(sp), intent(out) :: zero(n+1)
    !            weight-array of dimension  n+1  containing the weights
    !                   weight(k)=w(k), k=0,1,2,...,n
        real(sp), intent(out) :: weight(n+1)
    !            ierr - an error flag inherited from the routine  gauss
        integer, intent(out), optional :: ierr
     
        ! The arrays  a,b  are needed for working space.
        real(sp) :: a(n+1), b(n+1)

        real(sp) :: epsma, p0, p1, pm1
        integer :: k, ierr_

        ! epsma is the machine single precision.
        epsma = r1mach(3)

        do k = 1, n+1
            a(k) = alpha(k)
            b(k) = beta(k)
        end do

        p0 = 0._sp
        p1 = 1._sp
        do k = 1, n
            pm1 = p0
            p0 = p1
            p1 = (end - a(k))*p0 - b(k)*pm1
        end do
        a(n+1) = end - b(n+1)*p0/p1
        call sgauss(n+1,a,b,epsma,zero,weight,ierr_)
        if (present(ierr)) ierr = ierr_
    end subroutine sradau

    SUBROUTINE sgauss(n,alpha,beta,eps,zero,weight,ierr)
        implicit real(a-h,o-z)
        implicit integer(i-n)
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

        integer, intent(in) :: n
        real, intent(in) :: alpha(n)
        real, intent(in) :: beta(n)
        real, intent(in) :: eps
        real, intent(out) :: zero(n)
        real, intent(out) :: weight(n)
        integer, intent(out), optional :: ierr

        ! The array  e  is needed for working space.
        real :: e(n)

        if (n < 1) then
            ierr = -1 ! n is not in range
            return
        end if

        ierr = 0 ! normal return

        zero(1) = alpha(1)
        if (beta(1) < 0.) then
          ierr = -2 ! beta(1) is negative
          return
        end if
        
        weight(1) = beta(1)
        if (n == 1) return
        weight(1) =1.

        e(n) = 0.
        do k = 2, n
            zero(k) = alpha(k)
            if(beta(k) < 0.) then
                ierr = -2 ! beta(k) is negative
                return
            end if
            e(k-1) = sqrt(beta(k))
            weight(k) = 0.
        end do

        j = 0
        do l = 1, n
          
            ! Look for a small subdiagonal element.
          
            do m = l, n
                if (m == n) exit
                if (abs(e(m)) <= eps*(abs(zero(m)) + abs(zero(m+1)))) exit
            end do
            p = zero(l)
            
            if (m == l) then
                j = 0
                cycle
            end if

            if (j == 30) then
                ierr = l ! no convergence to an eigenvalue after 30 iterations.
                return
            end if
            j = j + 1
          
            ! Form shift.
          
            g = (zero(l+1) - p)/(2.*e(l))
            r = sqrt(g*g + 1.)
            g = zero(m) - p + e(l)/(g + sign(r,g))
            s = 1.
            c = 1.
            p = 0.
            mml = m-l
          
            ! For i=m-1 step -1 until l do ...
          
            do ii = 1, mml
                i = m-ii
                f = s*e(i)
                b = c*e(i)
                if (abs(f) < abs(g)) then
                    s = f/g
                    r = sqrt(s*s+1.)
                    e(i+1) = g*r
                    c = 1./r
                    s = s*c
                else
                    c = g/f
                    r = sqrt(c*c+1.)
                    e(i+1) = f*r
                    s = 1./r
                    c = c*s
                end if

                g = zero(i+1)-p
                r = (zero(i)-g)*s +2.*c*b
                p = s*r
                zero(i+1) = g+p
                g = c*r-b
            
                ! Form first component of vector.
            
                f = weight(i+1)
                weight(i+1) = s*weight(i) + c*f
                weight(i) = c*weight(i) - s*f
            end do
            zero(l) = zero(l) - p
            e(l) = g
            e(m) = 0.
        end do

        ! Order eigenvalues and eigenvectors.

        do ii = 2, n
            i = ii-1
            k = i
            p = zero(i)
            do  j = ii,n
                if (zero(j) >= p) cycle
                k = j
                p = zero(j)
            end do
            if(k == i) cycle
            zero(k) = zero(i)
            zero(i) = p
            p = weight(i)
            weight(i) = weight(k)
            weight(k) = p
        end do
        do k = 1, n
            weight(k) = beta(1)*weight(k)*weight(k)
        end do

    end subroutine sgauss



    subroutine slob(n,alpha,beta,aleft,right,zero,weight,ierr)

    ! Given  n  and a measure  dlambda, this routine generates the
    ! (n+2)-point Gauss-Lobatto quadrature formula

    !   integral over supp(dlambda) of f(x)dlambda(x)

    !      = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k))

    !              + w(n+1)f(x(n+1)) + R(n;f).

    ! The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
    ! =w(k), k=0,1,...,n,n+1. The user has to supply the recursion
    ! coefficients  alpha(k), beta(k), k=0,1,...,n,n+1, for the measure
    ! dlambda. The nodes and weights are computed in terms of the
    ! eigenvalues and first component of the normalized eigenvectors of
    ! a slightly modified Jacobi matrix of order  n+2. The routine calls
    ! upon the subroutine  gauss  and the function subroutine  r1mach.

    !   Input:  n - -  the number of interior points in the Gauss-Lobatto
    !                  formula; type integer
    !           alpha,beta - arrays of dimension  n+2  to be supplied with
    !                  the recursion coefficients  alpha(k-1), beta(k-1),
    !                  k=1,2,...,n+2, of the underlying measure; the
    !                  routine does not use  alpha(n+2), beta(n+2)
    !           aleft,right - the prescribed left and right endpoints
    !                  x(0)  and  x(n+1)  of the Gauss-Lobatto formula

    !   Output: zero - an array of dimension  n+2  containing the nodes (in
    !                  increasing order)  zero(k)=x(k), k=0,1,...,n,n+1
    !           weight-an array of dimension  n+2  containing the weights
    !                  weight(k)=w(k), k=0,1,...,n,n+1
    !           ierr - an error flag inherited from the routine  gauss

        integer, intent(in) :: n
        real(sp), intent(in) :: alpha(n+2)
        real(sp), intent(in) :: beta(n+2)
        real(sp), intent(in) :: aleft
        real(sp), intent(in) :: right
        real(sp), intent(out) :: zero(n+2)
        real(sp), intent(out) :: weight(n+2)
        integer, intent(out) :: ierr
        
        ! The arrays  e,a,b  are needed for working space.
        real(sp) :: a(n+2), b(n+2)

        ! The arrays  alpha,beta,zero,weight,e,a,b  are assumed to have
        ! dimension  n+2.

        integer :: np1, np2, k
        real(sp) :: epsma,p0l,p0r,p1l,p1r,pm1l,pm1r,det

        ! epsma is the machine single precision.
        epsma=r1mach(3)


        np1 = n+1
        np2 = n+2
        do k = 1, np2
            a(k) = alpha(k)
            b(k) = beta(k)
        end do

        p0l = 0.
        p0r = 0.
        p1l = 1.
        p1r = 1.
        do k = 1, np1
            pm1l = p0l
            p0l = p1l
            pm1r = p0r
            p0r = p1r
            p1l = (aleft-a(k))*p0l-b(k)*pm1l
            p1r = (right-a(k))*p0r-b(k)*pm1r
        end do

        det = p1l*p0r - p1r*p0l
        a(np2) = (aleft*p1l*p0r - right*p1r*p0l)/det
        b(np2) = (right - aleft)*p1l*p1r/det
        
        call sgauss(np2,a,b,epsma,zero,weight,ierr)

    end subroutine slob


! This is an auxiliary function routine providing a starting backward
! recurrence index for the Hermite measure that can be used in place
! of  nu0  in the routines  knum  and  dknum.
    integer function nu0her(n,z,eps)
        integer, intent(in) :: n
        complex(sp), intent(in) :: z
        real(sp), intent(in) :: eps

        nu0her = int(2.0_sp*(sqrt(0.5_sp*real(n+1)) + &
            0.25_sp*log(1.0_sp/eps)/abs(aimag(z)))**2)
    end function nu0her

! This is an auxiliary function routine providing a starting backward
! recurrence index for the Jacobi measure that can be used in place
! of  nu0  in the routines  knum  and  dknum.
    integer function nu0jac(n,z,eps)
        integer, intent(in) :: n
        complex(sp), intent(in) :: z
        real(sp), intent(in) :: eps

        real(sp) :: x,y,angle,r,x2,y2

        x = real(z,sp)
        y = abs(aimag(z))
        if(x < 1.) then
          if(x < -1.) angle=.5*(2.*pi+atan(y/(x-1.))+atan(y/(x+1.)))
          if(x == -1.) angle=.5*(1.5*pi-atan(.5*y))
          if(x > -1.) angle=.5*(pi+atan(y/(x-1.))+atan(y/(x+1.)))
        else
          if(x == 1.) angle=.5*(.5*pi+atan(.5*y))
          if(x > 1.) angle=.5*(atan(y/(x-1.))+atan(y/(x+1.)))
        end if
        x2 = x*x
        y2 = y*y
        r = ((x2-y2-1.)**2+4.*x2*y2)**.25
        r = sqrt((x+r*cos(angle))**2 + (y+r*sin(angle))**2)
        nu0jac = int(real(n+1)+.5*log(1./eps)/log(r))
    end function nu0jac

! This is an auxiliary function routine providing a starting backward
! recurrence index for the Laguerre measure that can be used in place
! of  nu0  in the routines  knum  and  dknum.
    integer function nu0lag(n,z,al,eps)
        integer, intent(in) :: n
        complex(sp), intent(in) :: z
        real(sp), intent(in) :: al
        real(sp), intent(in) :: eps

        real(sp) :: x,y,phi

        x = real(z,sp)
        y = aimag(z)
        phi = 0.5_sp*pi
        if(y < 0.) phi = 1.5_sp*pi
        if(x == 0.) go to 10
        phi = atan(y/x)
        if(y > 0. .and. x > 0.) go to 10
        phi = phi + pi
        if(x < 0.) go to 10
        phi = phi + pi

  10    nu0lag=int((sqrt(real(n+1)+.5*(al+1.))+alog(1./eps)/(4.*(x*x+  &
            y*y)**.25*cos(.5*(phi-pi))))**2-.5*(al+1.))
    end function nu0lag

end module