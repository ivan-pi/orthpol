c
c
      program test10
c
c
      dimension a(31),b(31),alpha(31),beta(31),z(11),w(11),e(11),
     *a1(31),b1(31),betap(20,12),erram(12),errbm(12)
      double precision depsma,d1mach,da(31),db(31),dalpha(31),dbeta(31),
     *dz(11),dw(11),de(11),da1(31),db1(31)
c
c epsma and depsma are the machine single and double precision.
c
      epsma=r1mach(3)
      depsma=d1mach(3)
c
c This test applies the routines  indp  and  dindp  to generate the
c first 20 recursion coefficients of the induced Legendre polynomials
c pind(k,m)(.), m=0,1,2,...,11, that is, of the polynomials orthogonal
c relative to the weight function
c
c                 [p(m)(x)]**2    on [-1,1],
c
c where  p(m)(.)  is the (monic) Legendre polynomial of degree m.
c (When m=0, then  pind(k,0)(.)=p(k)(.).) The routine also prints the
c absolute and relative errors, respectively, of the alpha- and beta-
c coefficients.
c
      n=20
      do 20 im=1,12
        m=im-1
        npm=n+m
c
c Generate the Legendre recurrence coefficients required in the
c routines  indp  and dindp.
c
        call recur(npm,1,0.,0.,a,b,ierr) 
        call drecur(npm,1,0.d0,0.d0,da,db,ierr)
c
c Compute the desired recursion coefficients.
c
        call indp(n,m,a,b,epsma,alpha,beta,ierr,z,w,e,a1,b1)
        call dindp(n,m,da,db,depsma,dalpha,dbeta,ierr,dz,dw,
     *    de,da1,db1)
c
c Compute and print the respective errors.
c
        erram(im)=0.
        errbm(im)=0.
        do 10 k=1,n
          erra=sngl(dabs(dble(alpha(k))-dalpha(k)))
          errb=sngl(dabs((dble(beta(k))-dbeta(k))/dbeta(k)))
          if(erra.gt.erram(im)) erram(im)=erra
          if(errb.gt.errbm(im)) errbm(im)=errb
          betap(k,im)=sngl(dbeta(k))
   10   continue
   20 continue
      do 40 ip=1,3
        ip4=1+4*(ip-1)
        write(*,1) ip4-1,ip4,ip4+1,ip4+2
    1   format(5x,'k',5x,'m=',i1,'  beta(k)',3x,'m=',i1,'  beta(k)',3x,
     *    'm=',i2,' beta(k)',3x,'m=',i2,' beta(k)'/)
        do 30 k=1,n
          km1=k-1
          write(*,2) km1,betap(k,ip4),betap(k,ip4+1),betap(k,ip4+2),
     *      betap(k,ip4+3)
    2     format(1x,i5,4f15.7)
   30   continue
        write(*,3) erram(ip4),erram(ip4+1),erram(ip4+2),erram(ip4+3)
    3   format(/4x,'erra',e14.4,3e15.4)
        write(*,4) errbm(ip4),errbm(ip4+1),errbm(ip4+2),errbm(ip4+3)
    4   format(4x,'errb',e14.4,3e15.4//)
   40 continue
      stop
      end

      subroutine indp(n,m,a,b,eps,alpha,beta,ierr,z,w,e,a1,b1)
c
c If  p(m)(.)  denotes the (monic) orthogonal polynomial of degree  m
c relative to the weight function  w(x), then the corresponding m-th 
c induced orthogonal polynomials  pind(k,m)(.), k=0,1,2,..., are those
c orthogonal with respect to the weight function
c
c                   (p(m)(x)**2)*w(x).  
c
c (For background on induced orthogonal polynomials, including an
c algorithm for generating their recursion coefficients, see W. Gautschi
c and S. Li,A set of orthogonal polynomials induced by a given
c orthogonal polynomial'', Aequationes Math., to appear.) This routine
c obtains the first n recurrence coefficients of the m-th induced
c orthogonal polynomials by an m-fold application of the routine  chri
c with  iopt=7, the shifts taken being, in succession, the zeros of
c p(m)(.).
c
      dimension a(*),b(*),alpha(*),beta(*),z(m),w(m),e(m),
     *a1(*),b1(*)
c
c The arrays  a,b,alpha,beta,a1,b1  are assumed to have dimension  n+m.
c
      npm=n+m
      do 10 k=1,npm
        alpha(k)=a(k)
        beta(k)=b(k)
   10 continue
      if(m.eq.0) return
      call gauss(m,a,b,eps,z,w,ierr,e)
      do 30 imu=1,m
        mi=npm-imu
        do 20 k=1,mi+1
          a1(k)=alpha(k)
          b1(k)=beta(k)
   20   continue
        x=z(imu)
        call chri(mi,7,a1,b1,x,0.,0.,0.,alpha,beta,ierrc)
   30 continue
      return
      end

      subroutine dindp(n,m,da,db,deps,dalpha,dbeta,ierr,dz,dw,de,
     *da1,db1)
c
c This is a double-precision version of the routine  indp.
c
      double precision da(*),db(*),deps,dalpha(*),dbeta(*),
     *dz(m),dw(m),de(m),da1(*),db1(*),dx
c
c The arrays  da,db,dalpha,dbeta,da1,db1  are assumed to have
c dimension  n+m.
c
      npm=n+m
      do 10 k=1,npm
        dalpha(k)=da(k)
        dbeta(k)=db(k)
   10 continue
      if(m.eq.0) return
      call dgauss(m,da,db,deps,dz,dw,ierr,de)
      do 30 imu=1,m
        mi=npm-imu
        do 20 k=1,mi+1
          da1(k)=dalpha(k)
          db1(k)=dbeta(k)
   20   continue
        dx=dz(imu)
        call dchri(mi,7,da1,db1,dx,0.d0,0.d0,0.d0,dalpha,
     *    dbeta,ierrc)
   30 continue
      return
      end

