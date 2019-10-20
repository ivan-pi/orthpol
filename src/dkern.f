c
c
C       subroutine dkern(n,nu0,numax,dx,dy,deps,da,db,dkerr,dkeri,
C      *  nu,ierr,droldr,droldi)
      subroutine dkern(n,nu0,numax,dz,deps,da,db,dker,nu,ierr,drold)
c
c This is a double-precision version of the routine  kern.
c
C       double precision dx,dy,deps,da(numax),db(numax),dkerr(*),
C      *  dkeri(*),droldr(*),droldi(*),dp0r,dp0i,dpr,dpi,dpm1r,
C      *  dpm1i,dden,dt
      double precision deps, da(numax), db(numax)
      complex(kind(1.0d0)) dz, dker, drold, dp0, dp, dpm1
      dimension dker(*), drold(*)
c
c The arrays  dkerr,dkeri,droldr,droldi  are assumed to have
c dimension  n+1.
c
C       call dknum(n,nu0,numax,dx,dy,deps,da,db,dkerr,dkeri,nu,ierr,
C      *  droldr,droldi)
      call dknum(n,nu0,numax,dz,deps,da,db,dker,nu,ierr,drold)
      if(ierr.ne.0) return
C       dp0r=0.d0
C       dp0i=0.d0
      dp0=(0.0d0,0.0d0)
C       dpr=1.d0
C       dpi=0.d0
      dp=(1.0d0,0.0d0)
      do 10 k=1,n
C         dpm1r=dp0r
C         dpm1i=dp0i
        dpm1=dp0
C         dp0r=dpr
C         dp0i=dpi
        dp0=dp
C         dpr=(dx-da(k))*dp0r-dy*dp0i-db(k)*dpm1r
C         dpi=(dx-da(k))*dp0i+dy*dp0r-db(k)*dpm1i
        dp=(dz-da(k))*dp0-db(k)*dpm1
C         dden=dpr**2+dpi**2
C         dt=(dkerr(k+1)*dpr+dkeri(k+1)*dpi)/dden
C         dkeri(k+1)=(dkeri(k+1)*dpr-dkerr(k+1)*dpi)/dden
C         dkerr(k+1)=dt
        dker(k+1)=dker(k+1)/p
   10 continue
      return
      end

