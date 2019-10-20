c
c
      program test7
c
c
      external quad,dquad
      dimension endl(4),endr(4),xp(1),yp(1),xfer(100),wfer(100),
     *alpha(40),beta(40),be(40),x(100),w(100),xm(400),wm(400),p0(400),
     *p1(400),p2(400)
      double precision d1mach,depsma,dendl(4),dendr(4),dxp(1),dyp(1),
     *dxfer(250),dwfer(250),dalpha(40),dbeta(40),dbe(40),dx(250),
     *dw(250),dxm(1000),dwm(1000),dp0(1000),dp1(1000),dp2(1000),di,deps
      logical finl,finr,finld,finrd
c
c This test generates in single and double precision the first 40
c recurrence coefficients of the orthogonal polynomials for the half-
c range Hermite weight function
c
c                   exp(-x**2)  on  (0,oo).
c
c Printed are the double-precision values of the alpha- and beta-
c coefficients along with the respective relative errors of the single-
c precision values.
c
      finl=.true.
      finr=.false.
      finld=.true.
      finrd=.false.
      epsma=r1mach(3)
      depsma=d1mach(3)
c
c epsma and depsma are the machine single and double precision.
c
      iq=2
      idelta=1
      irout=1
      n=40
      mc=4
      mcd=4
      mp=0
      ncapm=100
      ncpmd=250
c
c Set up the partition for the discretization of the inner product.
c
      do 10 i=1,4
        fi=real(i)
        di=dble(fi)
        endl(i)=3.*(fi-1.)
        endr(i)=3.*fi
        dendl(i)=3.d0*(di-1.d0)
        dendr(i)=3.d0*di
   10 continue
      eps=50.*epsma
      deps=1000.*depsma
c
c Compute the desired recursion coefficients by the multiple-component
c discretization procedure. On the third and fourth subinterval of the
c partition, the quadrature weights produced by  qgp  resp.  dqgp  may
c underflow, on the fourth subinterval even on machines with large
c exponent range.
c
      call smcdis(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta,irout,finl,
     *  finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie,be,x,w,
     *  xm,wm,p0,p1,p2)
      call dmcdis(n,ncpmd,mcd,mp,dxp,dyp,dquad,deps,iq,idelta,irout,
     *  finld,finrd,dendl,dendr,dxfer,dwfer,dalpha,dbeta,ncapd,kountd,
     *  ierrd,ied,dbe,dx,dw,dxm,dwm,dp0,dp1,dp2)
      write(*,1) ncap,kount,ierr,ie
    1 format(/1x,'ncap = ',i3,' kount = ',i2,' ierr = ',i3,' ie = ',i3)
      write(*,2) ncapd,kountd,ierrd,ied
    2 format(1x,'ncapd =',i3,' kountd =',i2,' ierrd =',i3,' ied =',i3/)
c
c Print the results.
c
      write(*,3)
    3 format(/5x,'k',9x,'dalpha(k)',15x,'dbeta(k)')
      write(*,4)
    4 format(10x,'erra',20x,'errb')
      do 20 k=1,n
        km1=k-1
        erra=abs(sngl((dble(alpha(k))-dalpha(k))/dalpha(k)))
        errb=abs(sngl((dble(beta(k))-dbeta(k))/dbeta(k)))
        write(*,5) km1,dalpha(k),dbeta(k)
    5   format(1x,i5,2d24.16)
        write(*,6) erra,errb
    6   format(6x,e12.4,12x,e12.4)
   20 continue
      stop
      end 

      subroutine quad(n,x,w,i,ierr)
      dimension x(n),w(n)
      print *,' User has selected the wrong SP-quadrature routine.'
      stop
      end

      subroutine dquad(n,dx,dw,i,ierr)
      double precision dx(n),dw(n)
      print *,' User has selected the wrong DP-quadrature routine.'
      stop
      end

      function wf(x,i)
      wf=exp(-x*x)
      return
      end

      double precision function dwf(dx,i)
      double precision dx
      dwf=dexp(-dx*dx)
      return
      end

