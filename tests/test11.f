c
c
      program test11
      use orthpol, only: sgchri, dgchri, schri, dchri
c
c
      complex rho,rold,z,e
      dimension xx(5),rr(5),a(500),b(500),alpha(40),beta(40),alphc(40),
     *betc(40),fnu(80),rho(80),rold(80),s(40),s0(80),s1(80),s2(80),
     *alphr(40),betr(40),alphcr(40),betcr(40)

      complex(kind(1.0d0)) drho(80), drold(80), dz, de
      double precision d1mach,depsma,deps,dal,dbe,dhi,da(800),db(800),
     *dx,dy,dalpha(40),dbeta(40),dnu(80),ds(40),ds0(80),ds1(80),ds2(80),
     *dhr,dalphc(40),dbetc(40),dalphr(40),dbetr(40),dalcr(40),dbecr(40)
      data xx/1.001,1.01,1.04,1.07,1.1/
      data rr/1.05,1.1625,1.275,1.3875,1.5/
c
c This test is to illustrate the dissimilar performance of the routines
c chri  and  gchri  in the case of division of the Jacobi weight
c function  w(t;alj,bej)  with parameters  alj,bej  by either a linear
c divisor  t-x  or a quadratic divisor  (t-x)**2 + y**2 . In either
c case, the parameters selected are  alj=-.8(.4).8, bej=alj(.4).8.  In
c the former case, x = -1.001, -1.01, -1.04, -1.07 and -1.1, whereas in
c the latter case,  x and  y  are chosen to lie, regularly spaced, on
c the upper half of an ellipse with foci at +1 and -1 and sum of the
c semiaxes equal to  rho = 1.05, 1.1625, 1.275, 1.3875 and 1.5. The
c routines are run in both single and double precision with n=40, the 
c results of the latter being used to calculate, and print, the maximum
c absolute and relative error of the single-precision alpha- and beta-
c coefficients, respectively. Also printed are the starting recurrence 
c indexes required in the backward recurrence schemes of  gchri,dgchri  
c to achieve single- resp. double-precision accuracy. This information 
c is contained in the first line of each 3-line block of the output,
c where in the case of quadratic divisors only average values (averaged
c over the upper half of the respective ellipse) are shown. The second
c and third line of each 3-line block display the maximum 
c reconstruction error'', that is, the maximum errors in the alpha's
c and beta's if the coefficients produced by  gchri,chri  and  dgchri,
c dchri  are fed back to the routines  chri  and  dchri  with  iopt=1
c to recover the original recursion coefficients in single and double
c precision.
c
      write(*,1)
    1 format(/)
      epsma=r1mach(3)
      depsma=d1mach(3)
c
c epsma and depsma are the machine single and double precision.
c
      n=40
      np1=n+1
      nm1=n-1
      nd=2*n
      ndm1=nd-1
      numax=500
      numaxd=800
      eps=10.*epsma
      deps=100.d0*depsma
      epsd=sngl(deps)
      ipoly=6
      do 70 ial=1,5
        al=-.8+.4*real(ial-1)
        dal=dble(al)
        ibemax=6-ial
        do 60 ibe=1,ibemax
          be=al+.4*real(ibe-1)
          dbe=dble(be)
          write(*,2) al,be
    2     format(///1x,'al = ',f6.2,'  be = ',f6.2//)
          hi=0.
          dhi=0.d0
c
c Generate the Jacobi recurrence coefficients to be used in the
c backward recurrence algorithm of the routines  gchri  and  dgchri.
c
          call srecur(numax,ipoly,al,be,a,b,ierr)
          call drecur(numaxd,ipoly,dal,dbe,da,db,ierrd)
          write(*,3)
    3     format(30x,'gchri',20x,'chri')
          write(*,4)
    4     format(5x,'x',5x,'nu0',2x,'nud0',4x,'erra',8x,'errb',10x,
     *      'erra',7x,'errb'/)

          do 20 ix=1,5
            x=-xx(ix)
            dx=dble(x)
            y=0.
            dy=0.d0
            z=cmplx(x,y)
c
c Compute the starting index for backward recurrence.
c
            nu0=nu0jac(ndm1,z,eps)
            nu0d=nu0jac(ndm1,dz,epsd)
c
c Generate the recurrence coefficients for the Jacobi weight function
c divided by a linear divisor, using the routines  gchri,dgchri.
c
            call sgchri(n,1,nu0,numax,eps,a,b,x,y,alpha,beta,nu,ierrg,
     *        ierrc,fnu,rho,rold,s,s0,s1,s2)
c
c On machines with limited single-precision exponent range, the routine
c cheb  used in  gchri  may have generated an underflow exception,
c which however is harmless and can be ignored.
c
            call dgchri(n,1,nu0d,numaxd,deps,da,db,dx,dy,dalpha,dbeta,
     *        nud,ierrgd,ierrcd,dnu,drho,drold,ds,ds0,
     *        ds1,ds2)
            if(ierrg.ne.0 .or. ierrc.ne.0 .or. ierrgd.ne.0 .or.ierrcd
     *        .ne.0) then
              write(*,5) ierrg,ierrgd,al,be,x
    5         format(/1x,'ierrg in gchri = ',i4,' ierrg in dgchri = ',
     *          i4,' for al = ',f6.2,' be = ',f6.2,' x = ',f7.4)
              write(*,6) ierrc,ierrcd,al,be,x
    6         format(1x,'ierrc in gchri = ',i4,' ierrc in dgchri = ',
     *          i4,' for al = ',f6.2,' be = ',f6.2,' x = ',f7.4/) 
              goto 20
            end if
c
c Generate the recurrence coefficients for the Jacobi weight function
c divided by a linear divisor, using the routines  chri,dchri.
c
            hr=real(rho(1))
            dhr=real(drho(1))
            call schri(n,4,a,b,x,y,hr,hi,alphc,betc,ierr)
            call dchri(n,4,da,db,dx,dy,dhr,dhi,dalphc,dbetc,ierr)
c
c Do the reconstruction.
c
            call schri(nm1,1,alpha,beta,x,y,0.,0.,alphr,betr,ierr)
            call dchri(nm1,1,dalpha,dbeta,dx,dy,0.d0,0.d0,dalphr,dbetr,
     *        ierr)
            call schri(nm1,1,alphc,betc,x,y,0.,0.,alphcr,betcr,ierr)
            call dchri(nm1,1,dalphc,dbetc,dx,dy,0.d0,0.d0,dalcr,dbecr,
     *        ierr)
c
c Compute and print the maximum errors.
c
            erragm=0.
            errbgm=0.
            erracm=0.
            errbcm=0.
            erram=0.
            errbm=0.
            errdam=0.
            errdbm=0.
            eracrm=0.
            erbcrm=0.
            edacrm=0.
            edbcrm=0.
            do 10 k=1,n
              km1=k-1
              errag=abs(sngl(dble(alpha(k))-dalpha(k)))
              errbg=abs(sngl((dble(beta(k))-dbeta(k))/dbeta(k)))
              errac=abs(sngl(dble(alphc(k))-dalphc(k)))
              errbc=abs(sngl((dble(betc(k))-dbetc(k))/dbetc(k)))
              if(k.lt.n) then
                erra=abs(alphr(k)-a(k))
                errb=abs((betr(k)-b(k))/b(k))
                errda=abs(sngl(dalphr(k)-da(k)))
                errdb=abs(sngl((dbetr(k)-db(k))/db(k)))
                eracr=abs(alphcr(k)-a(k))
                erbcr=abs((betcr(k)-b(k))/b(k))
                edacr=abs(sngl(dalcr(k)-da(k)))
                edbcr=abs(sngl((dbecr(k)-db(k))/db(k)))
                if(erra.gt.erram) erram=erra
                if(errb.gt.errbm) errbm=errb
                if(errda.gt.errdam) errdam=errda
                if(errdb.gt.errdbm) errdbm=errdb
                if(eracr.gt.eracrm) eracrm=eracr
                if(erbcr.gt.erbcrm) erbcrm=erbcr
                if(edacr.gt.edacrm) edacrm=edacr
                if(edbcr.gt.edbcrm) edbcrm=edbcr
              end if
              if(errag.gt.erragm) erragm=errag
              if(errbg.gt.errbgm) errbgm=errbg
              if(errac.gt.erracm) erracm=errac
              if(errbc.gt.errbcm) errbcm=errbc
   10       continue
            write(*,7) x,nu0,nu0d,erragm,errbgm,erracm,errbcm
    7       format(/1x,f7.4,2i6,2e12.4,2x,2e12.4)
            if(ix.eq.1) then
              write(*,8) erram,errbm,errdam,errdbm
    8         format(11x,'reconstr.',2e12.4,2x,2e12.4)
              write(*,9) eracrm,erbcrm,edacrm,edbcrm
    9         format(12x,'errors',2x,2e12.4,2x,2e12.4)
            else
              write(*,11) erram,errbm,errdam,errdbm
   11         format(20x,2e12.4,2x,2e12.4)
              write(*,11) eracrm,erbcrm,edacrm,edbcrm
            end if
   20     continue
          write(*,1)
          ndiv=20
          ndivm1=ndiv-1
          fndiv=real(ndiv)
          fndm1=real(ndivm1)
          pi=4.*atan(1.)
          write(*,3)
          write(*,12)
   12     format(4x,'rho',4x,'nu0',2x,'nud0',4x,'erra',8x,'errb',10x,
     *      'erra',7x,'errb'/)
          do 50 ir=1,5
            r=rr(ir)
            agmv=0.
            bgmv=0.
            acmv=0.
            bcmv=0.
            amv=0.
            bmv=0.
            damv=0.
            dbmv=0.
            acrmv=0.
            bcrmv=0.
            dacrmv=0.
            dbcrmv=0.
            nu0v=0
            nu0dv=0
            do 40 ith=1,ndivm1
c
c Generate the points on the ellipse.
c
              theta=pi*real(ith)/fndiv
              e=cmplx(cos(theta),sin(theta)) 
              z=.5*(r*e+1./(r*e))
              x=real(z)
              y=aimag(z)
              dx=dble(x)
              dy=dble(y)

c Compute the starting index for backward recurrence.
c
              nu0=nu0jac(ndm1,z,eps)
              nu0d=nu0jac(ndm1,z,epsd)
c
c Generate the recurrence coefficients for the Jacobi weight function
c divided by a quadratic divisor, using the routines  gchri,dgchri.
c
              call sgchri(n,2,nu0,numax,eps,a,b,x,y,alpha,beta,nu,ierrg,
     *          ierrc,fnu,rho,rold,s,s0,s1,s2)
              call dgchri(n,2,nu0d,numaxd,deps,da,db,dx,dy,dalpha,dbeta,
     *          nud,ierrgd,ierrcd,dnu,drho,drold,ds,ds0,
     *          ds1,ds2)
              if(ierrg.ne.0 .or.ierrc.ne.0 .or. ierrgd.ne.0 .or. ierrcd
     *          .ne.0) then
                write(*,5) ierrg,ierrgd,al,be,x
                write(*,6) ierrc,ierrcd,al,be,x
                goto 40
              end if
              nu0v=nu0v+nu0
              nu0dv=nu0dv+nu0d
c
c Generate the recurrence coefficients for the Jacobi weight function
c divided by a quadratic divisor, using the routines  chri,dchri.
c
              hr=real(rho(1))
              hi=aimag(rho(1))
              dhr=real(drho(1))
              dhi=aimag(drho(1))
              call schri(n,5,a,b,x,y,hr,hi,alphc,betc,ierr)
              call dchri(n,5,da,db,dx,dy,dhr,dhi,dalphc,dbetc,ierr)
c
c Do the reconstruction.
c
              call schri(nm1,2,alpha,beta,x,y,0.,0.,alphr,betr,ierr)
              call dchri(nm1,2,dalpha,dbeta,dx,dy,0.d0,0.d0,dalphr,
     *          dbetr,ierr)
              call schri(nm1,2,alphc,betc,x,y,0.,0.,alphcr,betcr,ierr)
              call dchri(nm1,2,dalphc,dbetc,dx,dy,0.d0,0.d0,dalcr,
     *          dbecr,ierr)
c
c Compute and print the maximum average errors.
c
              erragm=0.
              errbgm=0.
              erracm=0.
              errbcm=0.
              erram=0.
              errbm=0.
              errdam=0.
              errdbm=0.
              eracrm=0.
              erbcrm=0.
              edacrm=0.
              edbcrm=0.
              do 30 k=1,n
                km1=k-1
                errag=abs(sngl(dble(alpha(k))-dalpha(k)))
                errbg=abs(sngl((dble(beta(k))-dbeta(k))/dbeta(k)))
                errac=abs(sngl(dble(alphc(k))-dalphc(k)))
                errbc=abs(sngl((dble(betc(k))-dbetc(k))/dbetc(k)))
                if(k.lt.n) then
                  erra=abs(alphr(k)-a(k))
                  errb=abs((betr(k)-b(k))/b(k))
                  errda=abs(sngl(dalphr(k)-da(k)))
                  errdb=abs(sngl((dbetr(k)-db(k))/db(k)))
                  eracr=abs(alphcr(k)-a(k))
                  erbcr=abs((betcr(k)-b(k))/b(k))
                  edacr=abs(sngl(dalcr(k)-da(k)))
                  edbcr=abs(sngl((dbecr(k)-db(k))/db(k)))
                  if(erra.gt.erram) erram=erra
                  if(errb.gt.errbm) errbm=errb
                  if(errda.gt.errdam) errdam=errda
                  if(errdb.gt.errdbm) errdbm=errdb
                  if(eracr.gt.eracrm) eracrm=eracr
                  if(erbcr.gt.erbcrm) erbcrm=erbcr
                  if(edacr.gt.edacrm) edacrm=edacr
                  if(edbcr.gt.edbcrm) edbcrm=edbcr
                end if
                if(errag.gt.erragm) erragm=errag
                if(errbg.gt.errbgm) errbgm=errbg
                if(errac.gt.erracm) erracm=errac
                if(errbc.gt.errbcm) errbcm=errbc
   30         continue
              agmv=agmv+erragm
              bgmv=bgmv+errbgm
              acmv=acmv+erracm
              bcmv=bcmv+errbcm
              amv=amv+erram
              bmv=bmv+errbm
              damv=damv+errdam
              dbmv=dbmv+errdbm
              acrmv=acrmv+eracrm
              bcrmv=bcrmv+erbcrm
              dacrmv=dacrmv+edacrm
              dbcrmv=dbcrmv+edbcrm
   40       continue
            nu0=real(nu0v)/fndm1
            nu0d=dble(nu0dv)/fndm1
            erragm=agmv/fndm1
            errbgm=bgmv/fndm1
            erracm=acmv/fndm1
            errbcm=bcmv/fndm1
            erram=amv/fndm1
            errbm=bmv/fndm1
            errdam=damv/fndm1
            errdbm=dbmv/fndm1
            eracrm=acrmv/fndm1
            erbcrm=bcrmv/fndm1
            edacrm=dacrmv/fndm1
            edbcrm=dbcrmv/fndm1
            if(ir.eq.1) then
              write(*,7) r,nu0,nu0d,erragm,errbgm,erracm,errbcm
              write(*,8) erram,errbm,errdam,errdbm
              write(*,9) eracrm,erbcrm,edacrm,edbcrm
            else
              write(*,7) r,nu0,nu0d,erragm,errbgm,erracm,errbcm
              write(*,11) erram,errbm,errdam,errdbm
              write(*,11) eracrm,erbcrm,edacrm,edbcrm
            end if
   50     continue
   60   continue
   70 continue
      stop
      end

