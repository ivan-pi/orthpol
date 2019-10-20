c
c
      program test6
c
c
      external qlag,dqlag
      dimension a(100),b(100),e(100),xp(1),yp(1),endl(1),endr(1),
     *xfer(1),wfer(1),alpha(40),beta(40),be(40),x(100),w(100),
     *xm(200),wm(200),p0(200),p1(200),p2(200)
      double precision d1mach,da(300),db(300),de(300),depsma,dxp(1),
     *dyp(1),dendl(1),dendr(1),dxfer(1),dwfer(1),dalpha(40),dbeta(40),
     *dbe(40),dx(300),dw(300),dxm(600),dwm(600),dp0(600),dp1(600),
     *dp2(600),deps
      logical finl,finr,finld,finrd
      common/s/a,b,e,epsma
      common/d/da,db,de,depsma
c
c This test generates in single and double precision the first 40
c recursion coefficients of the orthogonal polynomials belonging to
c the logistic density function
c
c              exp(-x)/((1+exp(-x))**2)  on (-oo, oo).
c
c It prints the double-precision beta-coefficients (the alpha's being
c all zero) along with the absolute and relative errors of the alpha-
c coefficients resp. beta-coefficients.
c
      write(*,1)
    1 format(/)
c
c epsma and depsma are the machine single and double precision.
c
      iq=1
      idelta=1
      irout=1
      epsma=r1mach(3)
      depsma=d1mach(3)
      n=40
      mc=2
      mp=0
      ncapm=100
      ncapmm=mc*ncapm+mp
      ncpmd=300
      ncpmmd=mc*ncpmd+mp
      eps=5000.*epsma
      deps=1000.*depsma
c
c Compute the desired coefficients. On machines with limited exponent
c range, some of the weights in the Gauss-Laguerre quadrature rule may
c underflow. 
c
      call smcdis(n,ncapm,mc,mp,xp,yp,qlag,eps,iq,idelta,irout,finl,
     *  finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie,be,x,w,
     *  xm,wm,p0,p1,p2)
      call dmcdis(n,ncpmd,mc,mp,dxp,dyp,dqlag,deps,iq,idelta,irout,
     *  finld,finrd,dendl,dendr,dxfer,dwfer,dalpha,dbeta,ncapd,kountd,
     *  ierrd,ied,dbe,dx,dw,dxm,dwm,dp0,dp1,dp2)
      write(*,2) ncap,kount,ierr,ie
    2 format(/1x,'ncap = ',i3,' kount = ',i2,' ierr = ',i3,' ie = ',i5)
      write(*,3) ncapd,kountd,ierrd,ied
    3 format(1x,'ncapd= ',i3,' kountd= ',i2,' ierrd= ',i3,' ied= ',i5/) 
c
c Print the results.
c
      write(*,4)
    4 format(/5x,'k',8x,'dbeta(k)',13x,'erra',8x,'errb'/)
      do 10 k=1,n
        km1=k-1
        erra=abs(alpha(k))
        errb=abs(sngl((dble(beta(k))-dbeta(k))/dbeta(k)))
        if(ied.eq.0 .or. km1.lt.ied) then
          if(ie.eq.0 .or. km1.lt.ie) then
            write(*,5) km1,dbeta(k),erra,errb
    5       format(1x,i5,d24.16,2e12.4)
          else
            write(*,6) km1,dbeta(k)
    6       format(1x,i5,d24.16)
          end if
        end if
   10 continue
      stop
      end 

      subroutine qlag(n,x,w,i,ierr)
      dimension x(n),w(n),a(100),b(100),e(100)
      common/s/a,b,e,epsma
      call srecur(n,7,0.,0.,a,b,ierr)
      call sgauss(n,a,b,epsma,x,w,ierr,e)
      do 10 k=1,n
        w(k)=w(k)/((1.+exp(-x(k)))**2)
        if(i.eq.1) x(k)=-x(k)
   10 continue
      return
      end  

      subroutine dqlag(n,dx,dw,i,ierr)
      double precision dx(n),dw(n),da(300),db(300),de(300),depsma
      common/d/da,db,de,depsma
      call drecur(n,7,0.d0,0.d0,da,db,ierr)
      call dgauss(n,da,db,depsma,dx,dw,ierr,de)
      do 10 k=1,n
        dw(k)=dw(k)/((1.d0+dexp(-dx(k)))**2)
        if(i.eq.1) dx(k)=-dx(k)
   10 continue
      return
      end

      function wf(x,i)
c
c This is a dummy function. It is never called, since the routine
c qgp  in  mcdis  which requires  wf  is not activated in this test.
c
      wf=0.
      return
      end

      double precision function dwf(dx,i)
      double precision dx
c
c This is a dummy function. It is never called, since the routine
c dqgp  in  dmcdis  which requires  dwf  is not activated in this test.
c
      dwf=0.d0
      return
      end 

