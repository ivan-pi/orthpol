c
c
      program test4
      use orthpol, only: mcdis
c
c
      external qchle
      dimension a(81),b(81),e(81),xp(1),yp(1),endl(1),endr(1),
     *xfer(1),wfer(1),alpha(80),beta(80),be(80),x(81),w(81),
     *xm(162),wm(162),p0(162),p1(162),p2(162),betap(80,3)
      logical finl,finr
      common/s/c,a,b,e,epsma
c
c This is a test of the routine  mcdis, which is applied to generate
c the first 80 recurrence coefficients for the orthogonal polynomials
c belonging to the weight function
c
c       (1-x**2)**(-1/2) + c   on (-1,1),  c = 1, 10, 100.
c
c The corresponding inner product is discretized by applying the
c Gauss-Chebyshev quadrature rule to the first term of the weight
c function, and the Gauss-Legendre rule to the second term. In
c addition to the beta-coefficients (all alpha's are zero), the
c routine prints the variables  ncap  and  kount  to confirm 
c convergence after one iteration.
c 
      write(*,1)
    1 format(/)
c
c epsma is the machine single precision.
c
      epsma=r1mach(3)
      iq=1
      idelta=2
      irout=1
      n=80
      mc=2
      mp=0
      ncapm=81
      eps=5000.*epsma
      iem=0
      c=.1
      do 20 ic=1,3
        c=10.*c
c
c Compute the desired recursion coefficients. On machines with limited
c exponent range, harmless underflow may occur in the routine  gauss
c used in  qchle  to generate the Gauss-Legendre quadrature rule.
c
        call mcdis(n,ncapm,mc,mp,xp,yp,qchle,eps,iq,idelta,irout,finl,
     *    finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie,be,
     *    x,w,xm,wm,p0,p1,p2)
        write(*,2) ncap,kount,ierr,ie,c
    2   format(3x,'ncap = ',i3,' kount = ',i2,' ierr = ',i3,
     *    ' ie = ',i3,' for c = ',f5.1)
        if(abs(ie).gt.iem) iem=abs(ie)
        if(ie.ne.0 .and. abs(ie).le.n) then
          call mcdis(abs(ie)-1,ncapm,mc,mp,xp,yp,qchle,eps,iq,idelta,
     *    irout,finl,finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,
     *    ierr,ie,be,x,w,xm,wm,p0,p1,p2)
          write(*,2) ncap,kount,ierr,ie,c
          write(*,1)
        end if
c
c Assemble the results in an array.
c
        do 10 k=1,n
          if(ie.eq.0 .or. k-1.lt.abs(ie)) then
            betap(k,ic)=beta(k)
          else
            betap(k,ic)=0.
          end if
   10   continue
   20 continue
c
c Print the results.
c
      write(*,3)
    3 format(//3x,'k',2x,'beta(k), c=1',3x,'beta(k), c=10',2x,
     *  'beta(k), c=100'/)
      do 30 k=1,n
        km1=k-1
        if(iem.eq.0 .or. km1.lt.iem) write(*,4) km1,betap(k,1),
     *    betap(k,2),betap(k,3)
    4   format(1x,i3,3e15.7)
   30 continue
      stop
      end 

      subroutine qchle(n,x,w,i,ierr)
      dimension x(n),w(n),a(81),b(81),e(81)
      common/s/c,a,b,e,epsma
      fn=real(n)
      pi=4.*atan(1.)
      if(i.eq.1) then
        do 10 k=1,n
          fk=real(k)
          x(k)=cos((2.*fk-1.)*pi/(2.*fn))
          w(k)=pi/fn
   10   continue
      else
        call srecur(n,1,0.,0.,a,b,ierr)
        call sgauss(n,a,b,epsma,x,w,ierr,e)
        do 20 k=1,n
          w(k)=c*w(k)
   20   continue
      end if
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

