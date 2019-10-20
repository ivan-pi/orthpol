c
c
      program test5
c
c
      external qjac
      dimension a(41),b(41),e(41),xp(1),yp(1),endl(1),endr(1),xfer(1),
     *wfer(1),alpha(40), beta(40),be(40),x(41),w(41),xm(42),wm(42),
     *p0(42),p1(42),p2(42)
      double precision dy,dalj,dbej,dalpbe,da(41),db(41),daex(40),
     *dbex(40),dnum,dd,dkm1,dden
      logical finl,finr
      common/s/a,b,e,epsma
c
c This test generates the first 40 recursion coefficients for 
c polynomials orthogonal with respect to the Jacobi weight function
c with parameters  alj = -.8(.2)1., bej = -.8(.2)1.  and an added mass
c point of strength  y = .5, 1, 2, 4, 8  at the left end point. It 
c also computes the maximum relative errors (absolute errors for alpha-
c coefficients near zero) of the computed coefficients by comparing
c them against the exact coefficients known analytically.
c
      write(*,1)
    1 format(/)
c
c epsma is the machine single precision
c
      epsma=r1mach(3)
      iq=1
      idelta=2
      irout=1
      n=40
      mc=1
      mp=1
      xp(1)=-1.
      ncapm=41
      eps=5000.*epsma
      y=.25
      do 40 iy=1,5
        y=2.*y
        dy=dble(y)
        yp(1)=y
        write(*,2) y
    2   format(/1x,'y = ',f6.2/)
        write(*,3)
    3   format(2x,'alj',2x,'bej',5x,'erra',7x,'errb',6x,'alpha',
     *    4x,'beta',4x,'ka',2x,'kb',1x,'ierr',1x,'ie',1x,'it'/)
        do 30 ia=1,10
          alj=-1.+.2*real(ia)
          dalj=dble(alj)
          do 20 ib=1,10
            bej=-1.+.2*real(ib)
            alpbe=alj+bej
            dbej=dble(bej)
            dalpbe=dalj+dbej
c
c Generate the Jacobi recurrence coefficients.
c
            call recur(ncapm,6,alj,bej,a,b,ierr)
            call drecur(ncapm,6,dalj,dbej,da,db,iderr)
c
c Compute the desired recursion coefficients.
c
            call mcdis(n,ncapm,mc,mp,xp,yp,qjac,eps,iq,idelta,irout,
     *        finl,finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,
     *        ierr,ie,be,x,w,xm,wm,p0,p1,p2)
c
c Compute the exact coefficients by Eqs. (4.19)-(4.21) of the companion
c paper along with the relative errors (absolute errors for alpha-
c coefficients close to zero).
c
            daex(1)=(da(1)-dy)/(1.d0+dy)
            dbex(1)=1.d0+dy
            erra=abs(sngl(dble(alpha(1))-daex(1)))
            if(abs(sngl(daex(1))).gt.eps) erra=erra/abs(sngl(daex(1)))
            errb=abs(sngl((dble(beta(1))-dbex(1))/dbex(1)))
            erram=erra
            alpham=alpha(1)
            errbm=errb
            betam=beta(1)
            kam=0
            kbm=0
            dnum=1.d0+dy
            dd=1.d0
            do 10 k=2,n
              km1=k-1
              dkm1=dble(km1)
              dden=dnum
              if(k.gt.2) dd=(dbej+dkm1)*(dalpbe+dkm1)*dd/((dalj+dkm1
     *          -1.d0)*(dkm1-1.d0))
              dnum=(1.d0+(dbej+dkm1+1.d0)*(dalpbe+dkm1+1.d0)*dy*dd/
     *          (dkm1*(dalj+dkm1)))/(1.d0+dy*dd)
              daex(k)=da(k)+2.d0*dkm1*(dkm1+dalj)*(dnum-1.d0)/
     *          ((dalpbe+2.d0*dkm1)*(dalpbe+2.d0*dkm1+1.d0))
     *          +2.d0*(dbej+dkm1+1.d0)*(dalpbe+dkm1+1.d0)*((1.d0/dnum)
     *          -1.d0)/((dalpbe+2.d0*dkm1+1.d0)*(dalpbe+2.d0*dkm1+2.d0))
              dbex(k)=dnum*db(k)/dden
              erra=abs(sngl(dble(alpha(k))-daex(k)))
              if(abs(sngl(daex(k))).gt.eps) erra=erra/abs(sngl(daex(k)))
              errb=abs(sngl((dble(beta(k))-dbex(k))/dbex(k)))
              if(erra.gt.erram) then
                erram=erra
                alpham=alpha(k)
                kam=km1
              end if
              if(errb.gt.errbm) then
                errbm=errb
                betam=beta(k)
                kbm=km1
              end if
   10       continue
c
c Print the results.
c
            write(*,4) alj,bej,erram,errbm,alpham,betam,kam,kbm,ierr,
     *        ie,kount
    4       format(1x,2f5.2,2e11.4,2f9.6,4i4,i2)
   20     continue
          write(*,1)
   30   continue
   40 continue
      stop
      end 

      subroutine qjac(n,x,w,i,ierr)
      dimension x(n),w(n),a(41),b(41),e(41)
      common/s/a,b,e,epsma
      call gauss(n,a,b,epsma,x,w,ierr,e)
      do 10 k=1,n
        w(k)=w(k)/b(1)
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

