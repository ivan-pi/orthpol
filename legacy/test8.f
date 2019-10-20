c
c
      program test8
c
c
      external qcheb
      dimension oom2(7),xp(1),yp(1),endl(1),endr(1),xfer(1),wfer(1),
     *a(79),b(79),fnu(80),alpha(40),beta(40),be(40),x(500),w(500),
     *xm(500),wm(500),s(40),s0(80),s1(80),s2(80)
      logical finl,finr
      common/s/om2
      data oom2/.1,.3,.5,.7,.9,.99,.999/
c
c This test reproduces the results of  test1  in single precision,
c for n=40, using the routine  mccheb  in place of  cheb  and 
c Gauss-Chebyshev quadrature to discretize the modified moments.
c
      write(*,1)
    1 format(/)
      epsma=r1mach(3)
c
c epsma is the machine single precision.
c
      iq=1
      idelta=1
      n=40
      ndm1=2*n-1
      mc=1
      mp=0
      ncapm=500
      eps=100.*epsma
c
c Generate the recurrence coefficients for the (Chebyshev) polynomials
c defining the modified moments.
c
      call recur(ndm1,3,0.,0.,a,b,ierr)
c
c Compute the desired recursion coefficients by the discretized
c Chebyshev algorithm.
c
      do 20 iom=1,7
        om2=oom2(iom)
        call mccheb(n,ncapm,mc,mp,xp,yp,qcheb,eps,iq,idelta,finl,finr,
     *    endl,endr,xfer,wfer,a,b,fnu,alpha,beta,ncap,kount,ierr,be,x,
     *    w,xm,wm,s,s0,s1,s2)
c
c On machines with limited single-precision exponent range, the routine
c cheb  may have generated an underflow exception, which however is
c harmless and can be ignored.
c
        write(*,2) ncap,kount,ierr
    2   format(/'ncap=',i3,'  kount=',i3,'  ierr=',i3/)
c
c Print the results.
c
        write(*,3)
    3   format(5x,'k',4x,'beta(k)'/)
        do 10 k=1,n
          km1=k-1
          if(k.eq.1) then
            write(*,4) km1,beta(k),om2
    4       format(1x,i5,e14.6,'   om2 =',f6.3)
          else
            write(*,5) km1,beta(k)
    5     format(1x,i5,e14.6)
          end if
   10   continue
        write(*,1)
   20 continue      
      stop
      end

      subroutine qcheb(n,x,w,i,ierr)
      dimension x(n),w(n) 
      common/s/om2
      fn=real(n)
      pi=4.*atan(1.)
      do 10 k=1,n
        fk=real(k)
        x(k)=cos((2.*fk-1.)*pi/(2.*fn))
        w(k)=pi/(fn*sqrt(1.-om2*x(k)**2))
   10 continue
      return
      end

      function wf(x,i)
c
c This is a dummy function. It is never called, since the routine
c qgp  in  mccheb  which requires  wf  is not activated in this test.
c
      wf=0.
      return
      end

