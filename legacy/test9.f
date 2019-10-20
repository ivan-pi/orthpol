c
c
      program test9
c
c
      dimension a(199),b(199),fnu(200),alpha(100),beta(100),s(100),
     *s0(200),s1(200),s2(200),alphc(100),betc(100)
      double precision dsigma,da(199),db(199),dnu(200),dalpha(100),
     *dbeta(100),ds(100),ds0(200),ds1(200),ds2(200),dalphc(100),
     *dbetc(100)
c
c This test recomputes the results of  test2  for  sigma=.5  by applying
c the routine  chri  with  iopt=1, x=0  to the weight function with
c parameter  sigma=-.5. Printed are the relative discrepancies in both
c single and double precision between these results and those obtained
c in  test 2  by the modified Chebyshev algorithm. The test is embedded
c in the routine  test2, from which all print statements have been
c removed.
c
      logical modmom,intexp
      modmom=.true.
c
c Generate the recursion coefficients for the polynomials defining the
c modified resp. ordinary moments.
c
      if(modmom) then
        n=100
        ndm1=2*n-1
        call recur(ndm1,2,0.,0.,a,b,ierr)
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
        call cheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
c
c On machines with limited single-precision exponent range, the routine
c cheb  may generate an underflow exception, which however is harmless
c and can be ignored.
c
        call dcheb(n,da,db,dnu,dalpha,dbeta,ds,iderr,ds0,ds1,ds2)
c
c Up to this point the code is identical with the one of  test2.
c
        if(is.eq.1) then
          write(*,1) ierr,iderr
    1     format(/1x,'ierr in cheb = ',i4,' iderr in dcheb = ',i4/)
          if(ierr.ne.0) then
            nc=abs(ierr)
          else
            nc=n
          end if
          if(iderr.ne.0) then
            ncd=abs(iderr)
          else
            ncd=n
          end if 
c
c Compute the desired recursion coefficients by a modification
c algorithm.
c
          nm1=nc-1
          call chri(nm1,1,alpha,beta,0.,0.,0.,0.,alphc,betc,ierr)
          nm1=ncd-1
          call dchri(nm1,1,dalpha,dbeta,0.d0,0.d0,0.d0,0.d0,dalphc,
     *      dbetc,iderr)
        end if
        if(is.eq.3) then
          write(*,2)
    2     format(/1x,'test of the results for sigma=1/2'/)
          np=nc
          if(ncd.lt.nc) np=ncd
          nm1=np-1
c
c Compute and print the relative discrepancies between the results of
c the modified Chebyshev algorithm and the modification algorithm.
c
          write(*,3)
    3     format(3x,'k',2x,'err alpha',4x,'err beta',12x,'err dalpha',
     *      2x,'err dbeta'/)
          do 15 k=1,nm1
            km1=k-1
            errac=abs(alpha(k)-alphc(k))/alpha(k)
            errbc=abs(beta(k)-betc(k))/beta(k)
            errdac=sngl(dabs(dalpha(k)-dalphc(k))/dalpha(k))
            errdbc=sngl(dabs(dbeta(k)-dbetc(k))/dbetc(k))
            write(*,4) km1,errac,errbc,errdac,errdbc
    4       format(1x,i3,2e12.4,9x,2e12.4)
   15     continue
          write(*,5)
    5     format(/1x,'end of test'/)
        end if
c
c The rest of the code is essentially the same as the corresponding
c piece of code in  test2  with all print statements removed.
c
        eamax=0.
        ebmax=0.
        do 20 k=1,n
          km1=k-1
          erra=sngl(dabs(dble(alpha(k))-dalpha(k))/dalpha(k))
          errb=sngl(dabs(dble(beta(k))-dbeta(k))/dbeta(k))
          if(erra.gt.eamax) then
            eamax=erra
            kamax=km1
          end if
          if(errb.gt.ebmax) then
            ebmax=errb
            kbmax=km1
          end if
   20   continue
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
c =.false.). In either case, the input variable sigma is of type real.
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

