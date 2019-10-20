c
c
      program test1
c
c
      dimension fnu(160),f(80),f0(80),rr(80),a(159),b(159),
     *alpha(80),beta(80),s(80),s0(160),s1(160),s2(160)
      double precision doom2(7),deps,d1mach,dom2,dnu(160),d(80),
     *d0(80),drr(80),da(159),db(159),dalpha(80),dbeta(80),ds(80),
     *ds0(160),ds1(160),ds2(160)
      logical modmom
      data doom2/.1d0,.3d0,.5d0,.7d0,.9d0,.99d0,.999d0/
c
c This test generates the first n beta-coefficients in the recurrence
c relation for the orthogonal polynomials relative to the weight
c function
c
c         ((1-om2*x**2)*(1-x**2))**(-1/2)  on (-1,1)
c
c for om2=.1(.2).9,.99,.999, both in single and double precision,
c using modified moments if  modmom=.true.  and ordinary moments
c otherwise. In the former case, n=80, in the latter, n=20. Printed
c are the double-precision values of the coefficients along with the
c relative errors of the single-precision values.
c
      write(*,1)
    1 format(/)
      modmom=.true.
      eps=r1mach(3)
      deps=d1mach(3)
      if(modmom) then
        n=80
      else
        n=20
      end if
      ndm1=2*n-1
      do 30 iom=1,7
        dom2=doom2(iom)
        om2=sngl(dom2)
c
c Compute the modified resp. ordinary moments using Eqs. (3.7) and (3.9)
c of the companion paper. On machines with limited exponent range, some
c of the high-order modified moments may underflow, without this having
c any deteriorating effect on the accuracy.
c
        call fmm(n,eps,modmom,om2,fnu,ierr,f,f0,rr) 
        call dmm(n,deps,modmom,dom2,dnu,iderr,d,d0,drr)
        if(ierr.ne.0 .or. iderr.ne.0) then
          write(*,2) ierr,iderr,om2
    2     format(/5x,'ierr in fmm = ',i1,'  iderr in dmm = ',i1,
     *      '  for om2 = ',f8.4/)
          goto 30
        end if 
c
c Generate the recursion coefficients for the polynomials defining the
c modified resp. ordinary moments.
c
        if(modmom) then
          call recur(ndm1,3,0.,0.,a,b,ierr)
          call drecur(ndm1,3,0.d0,0.d0,da,db,iderr)
        else
          do 10 k=1,ndm1
            a(k)=0.
            b(k)=0.
            da(k)=0.d0
            db(k)=0.d0    
   10     continue
        end if
c
c Compute the desired recursion coefficients by means of the modified
c Chebyshev algorithm; for the latter, see, e.g., Section 2.4 of
c W. Gautschi, On generating orthogonal polynomials'', SIAM J. Sci.
c Statist. Comput. 3, 1982, 289-317.
c
        call cheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
c
c On machines with limited single-precision exponent range, the routine
c cheb  may generate an underflow exception, which however is harmless
c and can be ignored.
c
        call dcheb(n,da,db,dnu,dalpha,dbeta,ds,iderr,ds0,ds1,ds2)
        write(*,3) ierr,iderr
    3   format(/5x,'ierr in cheb = ',i3,'  iderr in dcheb = ',i3/)
        write(*,4)
    4   format(/5x,'k',9x,'dbeta(k)'/)
        do 20 k=1,n
          km1=k-1
          if(iderr.eq.0 .or. km1.lt.abs(iderr)) then
            if(ierr.eq.0 .or. km1.lt.abs(ierr)) then
              errb=sngl(dabs(dble(beta(k))-dbeta(k))/dbeta(k))
              if(k.eq.1) then
                write(*,5) km1,dbeta(k),errb,om2
    5           format(1x,i5,d24.16,e12.4,'   om2 =',f6.3)
              else
                write(*,6) km1,dbeta(k),errb
    6           format(1x,i5,d24.16,e12.4)
              end if
            else
              write(*,7) km1,dbeta(k)
    7         format(1x,i5,d24.16)
            end if
          end if
   20   continue
        write(*,1)
   30 continue
      stop
      end

      subroutine fmm(n,eps,modmom,om2,fnu,ierr,f,f0,rr)
c
c This routine generates the modified (Chebyshev) resp. ordinary
c moments of the weight function
c
c          ((1-om2*x**2)*(1-x**2))**(-1/2)  on (-1,1)
c
c using Eqs. (3.7) resp. (3.9) of the companion paper.
c
      dimension fnu(*),f(n),f0(n),rr(n)
      logical modmom
c
c The array  fnu  is assumed to have dimension  2*n.
c
      ierr=0
      nd=2*n
      ndm1=nd-1
      pi=4.*atan(1.)
c
c Compute the Fourier coefficients of ((1-om2*sin(theta)**2))**(-1/2)
c as minimal solution of a three-term recurrence relation as described
c on pp.310-311 of W. Gautschi,On generating orthogonal polynomials'',
c SIAM J. Sci. Statist. Comput. 3, 1982, 289-317.
c
      q=om2/(2.-om2+2.*sqrt(1.-om2))
      q1=(1.+q*q)/q
      do 10 k=1,n
        f(k)=0.
   10 continue
      nu=nd
   20 nu=nu+10
      do 30 k=1,n
        f0(k)=f(k) 
   30 continue
      if(nu.gt.500) then
        ierr=1
        return
      end if
      r=0.
      s=0.
      do 40 k=1,nu
        n1=nu-k+1
        fn1=real(n1)
        r=-(fn1-.5)/(fn1*q1+(fn1+.5)*r)
        s=r*(2.+s)
        if(n1.le.n) rr(n1)=r
   40 continue
      c0=1./(1.+s)
      f(1)=rr(1)*c0
      if(n.gt.1) then
        do 50 k=2,n
          f(k)=rr(k)*f(k-1)
   50   continue
      end if
      do 60 k=1,n
        if(abs(f(k)-f0(k)).gt.eps*abs(f(k))) goto 20
   60 continue
c
c Compute the desired modified resp. ordinary moments in term of
c the above Fourier coefficients.
c
      fnu(1)=pi*c0
      if(n.eq.1) return
      fnu(2)=0.
      if(n.eq.2) return
      if(modmom) then
        c=2.*pi
        do 70 k=3,ndm1,2
          k1=(k-1)/2
          c=-.25*c
          fnu(k)=c*f(k1)
          fnu(k+1)=0.
   70   continue
      else
        c=.5*pi
        fnu(3)=c*(c0-f(1))
        fnu(4)=0.
        c=-c
        do 90 k=5,ndm1,2
          k1=(k-1)/2
          k1m1=k1-1
          c=-.25*c
          c1=1.
          sum=f(k1)
          do 80 i=1,k1m1
            c1=-c1*real(2*k1-i+1)/real(i)
            sum=sum+c1*f(k1-i)
   80     continue
          c1=-c1*real(k1+1)/real(2*k1)
          sum=sum+c1*c0
          fnu(k)=c*sum
          fnu(k+1)=0.
   90   continue
      end if
      end

      subroutine dmm(n,deps,modmom,dom2,dnu,ierrd,d,d0,drr)
c
c This is a double-precision version of the routine  fmm.
c
      double precision deps,dom2,dnu(*),d(n),d0(n),drr(n),dpi,dq,
     *dq1,dr,ds,dn1,dc0,dc,dc1,dsum
      logical modmom
c
c The array  dnu  is assumed to have dimension  2*n.
c
      ierrd=0
      nd=2*n
      ndm1=nd-1
      dpi=4.d0*datan(1.d0)
      dq=dom2/(2.d0-dom2+2.d0*dsqrt(1.d0-dom2))
      dq1=(1.d0+dq*dq)/dq
      do 10 k=1,n
        d(k)=0.d0
   10 continue
      nud=nd
   20 nud=nud+10
      do 30 k=1,n
        d0(k)=d(k)
   30 continue
      if(nud.gt.1000) then
        ierrd=1
        return
      end if
      dr=0.d0
      ds=0.d0
      do 40 k=1,nud
        n1=nud-k+1
        dn1=dble(n1)
        dr=-(dn1-.5d0)/(dn1*dq1+(dn1+.5d0)*dr)
        ds=dr*(2.d0+ds)
        if(n1.le.n) drr(n1)=dr
   40 continue
      dc0=1.d0/(1.d0+ds)
      d(1)=drr(1)*dc0
      if(n.gt.1) then
        do 50 k=2,n
          d(k)=drr(k)*d(k-1)
   50   continue
      end if
      do 60 k=1,n
        if(dabs(d(k)-d0(k)).gt.deps*dabs(d(k))) goto 20
   60 continue
      dnu(1)=dpi*dc0
      if(n.eq.1) return
      dnu(2)=0.d0
      if(n.eq.2) return
      if(modmom) then
        dc=2.d0*dpi
        do 70 k=3,ndm1,2
          k1=(k-1)/2
          dc=-.25d0*dc
          dnu(k)=dc*d(k1)
          dnu(k+1)=0.d0
   70   continue
      else
        dc=.5d0*dpi
        dnu(3)=dc*(dc0-d(1))
        dnu(4)=0.d0
        dc=-dc
        do 90 k=5,ndm1,2
          k1=(k-1)/2
          k1m1=k1-1
          dc=-.25d0*dc
          dc1=1.d0
          dsum=d(k1)
          do 80 i=1,k1m1
            dc1=-dc1*dble(2*k1-i+1)/dble(i)
            dsum=dsum+dc1*d(k1-i)
   80     continue
          dc1=-dc1*dble(k1+1)/dble(2*k1)
          dsum=dsum+dc1*dc0
          dnu(k)=dc*dsum
          dnu(k+1)=0.d0
   90   continue
      end if
      end

