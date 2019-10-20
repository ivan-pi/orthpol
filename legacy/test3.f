c
c
      program test3
c
c
      dimension x(320),w(320),bexact(320),p0(320),p1(320),p2(320),
     *als(320),bes(320),all(320),bel(320)
c
c This test applies both the Stieltjes procedure (cf. Section 2.1 of
c W. Gautschi, On generating orthogonal polynomials'', SIAM J. Sci.
c Statist. Comput. 3, 1982, 289-317) and the Lanczos algorithm (cf.
c W.B. Gragg and W.J. Harrod, The numerically stable reconstruction of
c Jacobi matrices from spectral data'', Numer. Math. 44, 1984, 317-335)
c to generate the  N  recursion coefficients, N = 40, 80, 160, 320,
c for the (monic) orthogonal polynomials relative to the discrete inner
c product supported on  N  equally spaced points on [-1,1] (including
c the end points) and having equal weights 2/N. The routine prints the
c absolute errors in the alpha-coefficients and the relative errors in
c the beta-coefficients, these being computed using known formulae for
c the coefficients. The maxima of these errors are also printed.
c
      ncap=20
      do 30 incap=1,4
        ncap=2*ncap
        fncap=real(ncap)
        fncm1=real(ncap-1)
c
c Generate the abscissae and weights of the discrete inner product and
c the exact beta-coefficients.
c
        x(1)=-1.
        w(1)=2./fncap
        bexact(1)=2.
        do 10 k=2,ncap
          fkm1=real(k-1)
          x(k)=-1.+2.*fkm1/fncm1
          w(k)=w(1)
          bexact(k)=(1.+1./fncm1)**2*(1.-(fkm1/fncap)**2)/
     *      (4.-(1./fkm1)**2)
   10   continue
c
c Compute the desired coefficients, first by the Stieltjes procedure,
c and then by the Lanczos algorithm. Indicate via the error flag  ierrs
c whether a critical underflow condition has arisen in Stieltjes's
c procedure. (There may, in addition, occur harmless underflow, which
c the routine  sti  does not test for.)
c
        call sti(ncap,ncap,x,w,als,bes,ierrs,p0,p1,p2)
        call lancz(ncap,ncap,x,w,all,bel,ierrl,p0,p1)
        write(*,1) ierrs,ierrl
    1   format(/5x,'ierr in sti = ',i4,11x,'ierr in lancz = ',i3/)
c
c Compute and print the absolute errors of the alpha-coefficients and 
c the relative errors of the beta-coefficients as well as the maximum 
c respective errors.
c
        erralm=0.
        errblm=0.
        write(*,2)
    2   format(5x,'k',4x,'erra',8x,'errb',10x,'erra',8x,'errb'/)
        do 20 in=1,ncap
          inm1=in-1
          erras=abs(als(in))
          errbs=abs((bes(in)-bexact(in))/bexact(in))
          erral=abs(all(in))
          errbl=abs((bel(in)-bexact(in))/bexact(in))
          if(erral.gt.erralm) erralm=erral
          if(errbl.gt.errblm) errblm=errbl
          if(ierrs.eq.0 .or. inm1.lt.abs(ierrs)) then
            if(in.eq.1) then
              write(*,3) inm1,erras,errbs,erral,errbl,ncap
    3         format(1x,i5,2e12.4,2x,2e12.4,'   N =',i4)
            else
              write(*,4) inm1,erras,errbs,erral,errbl
    4         format(1x,i5,2e12.4,2x,2e12.4)
            end if
          else
            if(in.eq.1) then
              write(*,5) inm1,erral,errbl,ncap
    5         format(1x,i5,26x,2e12.4,'   N =',i4)
            else
              write(*,6) inm1,erral,errbl
    6         format(1x,i5,26x,2e12.4)
            end if
          end if
   20   continue
        write(*,7) erralm,errblm
    7   format(/32x,2e12.4//)
   30 continue
      stop
      end 

