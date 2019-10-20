c
c
      program test2
      use orthpol, only: srecur, drecur, scheb, dcheb
c
c
      dimension a(199),b(199),fnu(200),alpha(100),beta(100),s(100),
     *s0(200),s1(200),s2(200)
      double precision dsigma,da(199),db(199),dnu(200),dalpha(100),
     *dbeta(100),ds(100),ds0(200),ds1(200),ds2(200)
      logical modmom,intexp
c
c This test generates the first n recursion coefficients for the
c orthogonal polynomials relative to the weight function
c
c        (x**sigma)*ln(1/x)  on (0,1],  sigma = -.5, 0, .5,
c
c where n=100 when using modified (Legendre) moments, and n=12 when
c using ordinary moments. It prints the double-precision values of the
c coefficients as well as the relative errors of the respective single-
c precision values and the maximum relative errors.
c
      modmom=.true.
c
c Generate the recursion coefficients for the polynomials defining the
c modified resp. ordinary moments.
c
      if(modmom) then
        n=100
        ndm1=2*n-1
        call srecur(ndm1,2,0.,0.,a,b,ierr)
        call drecur(ndm1,2,0.d0,0.d0,da,db,iderr)
      else 
        n=12
        ndm1=2*n-1
        do 10 k=1,ndm1
          a(k)=0.
          b(k)=0.
          da(k)=0.d0
          db(k)=0.d0
   10   continue
      end if
      do 30 is=1,3
        dsigma=-.5d0+.5d0*dble(is-1)
        sigma=sngl(dsigma)
        if(is.eq.2) then
          intexp=.true.
        else
          intexp=.false.
        end if
c
c Compute the modified resp. ordinary moments using Eqs. (3.12) and
c (3.11) of the companion paper. On machines with limited exponent
c range, some of the high-order modified moments may underflow, without
c this having any deteriorating effect on the accuracy.
c
        call fmm(n,modmom,intexp,sigma,fnu)
        call dmm(n,modmom,intexp,dsigma,dnu)
c
c Compute the desired recursion coefficients by means of the modified
c Chebyshev algorithm; for the latter, see, e.g., Section 2.4 of
c W. Gautschi, On generating orthogonal polynomials'', SIAM J. Sci.
c Statist. Comput. 3, 1982, 289-317.
c
        call scheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
c
c On machines with limited single-precision exponent range, the routine
c cheb  may generate an underflow exception, which however is harmless
c and can be ignored.
c
        call dcheb(n,da,db,dnu,dalpha,dbeta,ds,iderr,ds0,ds1,ds2)
        write(*,1) ierr,iderr
    1   format(/6x,'ierr in cheb = ',i4,' iderr in dcheb = ',i4/)
c
c Compute and print the relative errors and their maxima.
c
        eamax=0.
        ebmax=0.
        write(*,2) sigma
    2   format(/21x,'sigma =',f5.1)
        write(*,3)
    3   format(/3x,'k',9x,'dalpha(k)',15x,'dbeta(k)'/)
        do 20 k=1,n
          km1=k-1
          if(iderr.eq.0 .or. km1.lt.abs(iderr)) then
            write(*,4) km1,dalpha(k),dbeta(k)
    4       format(1x,i3,2d24.16)
            if(ierr.eq.0 .or. km1.lt.abs(ierr)) then
              erra=sngl(dabs(dble(alpha(k))-dalpha(k))/dalpha(k))
              errb=sngl(dabs(dble(beta(k))-dbeta(k))/dbeta(k))
              write(*,5) erra,errb
    5         format(4x,e12.4,21x,e12.4)
              if(erra.gt.eamax) then
                eamax=erra
                kamax=km1
              end if
              if(errb.gt.ebmax) then
                ebmax=errb
                kbmax=km1
              end if
            end if
          end if
   20   continue
        write(*,6) eamax,kamax,ebmax,kbmax
    6   format(/2x,'eamax =',e11.4,' at',i3,4x,'ebmax =',e11.4,
     *    ' at',i3//)
   30 continue
      stop
      end

      subroutine fmm(n,modmom,intexp,sigma,fnu)
c
c This generates the first  2*n  modified moments (if modmom=.true.)
c relative to shifted monic Legendre polynomials, using Eq. (3.12) of
c the companion paper, and the first  2*n  ordinary moments (if modmom
c =.false.) by Eq. (3.11), of the weight function
c
c          (x**sigma)*ln(1/x)  on (0,1],   sigma > -1,
c
c for sigma an integer (if intexp=.true.) or a real number (if intexp
c =.false.). In either case, the input variable  sigma  is of type real.
c
      dimension fnu(*)
      logical modmom,intexp
c
c The array  fnu  is assumed to have dimension  2*n.
c
      nd=2*n
      sigp1=sigma+1.
      if(modmom) then
        isigma=int(sigma)
        isigp1=isigma+1
        isigp2=isigma+2
        isigp3=isigma+3
        if(intexp .and. isigp1.lt.nd) then
          kmax=isigp1
        else
          kmax=nd
        end if
        c=1.
        do 20 k=1,kmax
          km1=k-1
          fk=real(k)
          p=1.
          s=1./sigp1
          if(kmax.gt.1) then
            do 10 i=1,km1
              fi=real(i)
              p=(sigp1-fi)*p/(sigp1+fi)
              s=s+1./(sigp1+fi)-1./(sigp1-fi)
   10       continue
          end if
          fnu(k)=c*s*p/sigp1
          c=fk*c/(4.*fk-2.)
   20   continue
        if(.not.intexp .or. isigp1.ge.nd) return
        q=-.5
        if(isigma.gt.0) then
          do 30 iq=1,isigma
            fiq=real(iq)
            q=fiq*fiq*q/((2.*fiq+1.)*(2.*fiq+2.))
   30     continue
        end if
        fnu(isigp2)=c*q
        if(isigp2.eq.nd) return
        do 40 k=isigp3,nd
          km1=k-1
          fkm1=real(km1)
          fnu(k)=-fkm1*(fkm1-sigp1)*fnu(km1)/((4.*fkm1-2.)*
     *      (fkm1+sigp1))
   40   continue
        return
      else
        do 50 k=1,nd
          fkm1=real(k-1)
          fnu(k)=(1./(sigp1+fkm1))**2
   50   continue
      end if
      end

      subroutine dmm(n,modmom,intexp,dsigma,dnu)
c
c This is a double-precision version of the routine  fmm.
c
      double precision dsigma,dnu(*),dsigp1,dc,dk,dp,ds,di,dq,diq,dkm1
      logical modmom,intexp
c
c The array  dnu  is assumed to have dimension  2*n.
c
      nd=2*n
      dsigp1=dsigma+1.d0
      if(modmom) then
        isigma=idint(dsigma)
        isigp1=isigma+1
        isigp2=isigma+2
        isigp3=isigma+3
        if(intexp .and. isigp1.lt.nd) then
          kmax=isigp1
        else
          kmax=nd
        end if
        dc=1.d0
        do 20 k=1,kmax
          km1=k-1
          dk=dble(k)
          dp=1.d0
          ds=1.d0/dsigp1
          if(kmax.gt.1) then
            do 10 i=1,km1
              di=dble(i)
              dp=(dsigp1-di)*dp/(dsigp1+di)
              ds=ds+1.d0/(dsigp1+di)-1.d0/(dsigp1-di)
   10       continue
          end if
          dnu(k)=dc*ds*dp/dsigp1
          dc=dk*dc/(4.d0*dk-2.d0)
   20   continue
        if(.not.intexp .or. isigp1.ge.nd) return
        dq=-.5d0
        if(isigma.gt.0) then
          do 30 iq=1,isigma
            diq=dble(iq)
            dq=diq*diq*dq/((2.d0*diq+1.d0)*(2.d0*diq+2.d0))
   30     continue
        end if
        dnu(isigp2)=dc*dq
        if(isigp2.eq.nd) return
        do 40 k=isigp3,nd
          km1=k-1
          dkm1=dble(km1)
          dnu(k)=-dkm1*(dkm1-dsigp1)*dnu(km1)/((4.d0*dkm1-2.d0)*
     *      (dkm1+dsigp1))
   40   continue
        return
      else
        do 50 k=1,nd
          dkm1=dble(k-1)
          dnu(k)=(1.d0/(dsigp1+dkm1))**2
   50   continue
      end if
      end

