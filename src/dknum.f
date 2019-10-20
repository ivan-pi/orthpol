c
c
C       subroutine dknum(n,nu0,numax,dx,dy,deps,da,db,drhor,drhoi,nu,
C      *ierr,droldr,droldi)
      subroutine dknum(n,nu0,numax,dz,deps,da,db,drho,nu,
     *ierr,drold)
c
c This is a double-precision version of the routine  knum.
c
C       double precision dx,dy,deps,da(numax),db(numax),drhor(*),
C      *drhoi(*),droldr(*),droldi(*),drr,dri,dden,dt
      double precision deps,da(numax),db(numax)
      complex(kind(1.0d0)) dz,drho,drold,dr
      dimension drho(*), drold(*)
c
c The arrays  drhor,drhoi,droldr,droldi  are assumed to have
c dimension  n+1.
c
      ierr=0
      np1=n+1
      if(nu0.gt.numax) then
        ierr=nu0
        return
      end if
      if(nu0.lt.np1) nu0=np1
      nu=nu0-5
      do 10 k=1,np1
C         drhor(k)=0.d0
C         drhoi(k)=0.d0
        drho(k) = (0.0d0,0.0d0)
   10 continue
   20 nu=nu+5
      if(nu.gt.numax) then
        ierr=numax
        goto 60
      end if
      do 30 k=1,np1
C         droldr(k)=drhor(k)
C         droldi(k)=drhoi(k)
        drold(k)=drho(k)
   30 continue
C       drr=0.d0
C       dri=0.d0
      dr=(0.d0,0.0d0)
      do 40 j=1,nu
        j1=nu-j+1
C         dden=(dx-da(j1)-drr)**2+(dy-dri)**2
C         drr=db(j1)*(dx-da(j1)-drr)/dden
C         dri=-db(j1)*(dy-dri)/dden
C         if(j1.le.np1) then
C           drhor(j1)=drr
C           drhoi(j1)=dri
C         end if
        dr=cmplx(db(j1),0.d0,kind=kind(1.0d0))/(dz-cmplx(da(j1),0.d0,
     *    kind=kind(1.0d0))-dr)
        if(j1.le.np1) drho(j1)=dr
   40 continue
      do 50 k=1,np1
C         if((drhor(k)-droldr(k))**2+(drhoi(k)-droldi(k))**2.gt.
C      *    deps*(drhor(k)**2+drhoi(k)**2)) goto 20
        if(abs(drho(k)-drold(k)).gt.deps*abs(drho(k))) goto 20
   50 continue
   60 if(n.eq.0) return
      do 70 k=2,np1
C         dt=drhor(k)*drhor(k-1)-drhoi(k)*drhoi(k-1)
C         drhoi(k)=drhor(k)*drhoi(k-1)+drhoi(k)*drhor(k-1)
C         drhor(k)=dt
        drho(k)=drho(k)*drho(k-1)
   70 continue
      return
      end 

