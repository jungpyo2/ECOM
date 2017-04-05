!--------------------------------------------------------
! elliptic.f: This subroutine is for a Poisson solver in an unit disk
!--------------------------------------------------------

      subroutine elliptic(nsub,k,nth,rnd,rad,u,ur,uth,urr,urt,utt,
     1                     f,blength,bnodes,ibc,g)
c
c     fast solver for Poisson, modified Helmholtz or Helmholtz equation 
c     in disk with Dirichlet or outgoing boundary conditions at r=RAD.
c
c     A polar grid is used with arbitrary subdivisions in the r direction
c     and equispaced points in the theta direction.
c
c     we solve:    lap(u)             = f(r,theta)
c  
c
c
c     If ibc = 0, then we impose outgoing boundary conditions at r=RAD.
c     If ibc = 1, then we impose Dirichlet boundary conditions at r=RAD,
c
c     with boundary condition u = g at r=R.
c
c
c     INPUT:
c
c     nsub     = number of subintervals in r on [0,R].
c     k        = number of (scaled) Chebyshev nodes used on each subinterval
c     nth      = number of points in theta discretization.
c     rnd      = (interior) grid pts in r - does not include r=RAD.
c     rad      = outer radius of domain (disk).

c     f        = right-hand side of PDE
c     blength  = subinterval length array
c     bnodes   = subinterval boundary array
c     ibc      = boundary condition flag. 0= outgoing, 1= Dirichlet
c     g        = Dirichlet data needed if ibc = 1, otherwise not used.
c
c     OUTPUT:
c
c     u        = solution on same grid as f
c     ur       = r-deriv of solution 
c     uth      = theta-deriv of solution 
c     urr      = r-second deriv of solution 
c     urt      = r-deriv of t-deriv solution 
c     utt      = theta-second deriv of solution c

      use arrays, only: isolver, Halpha, ucoeff, urcoeff, urrcoeff
      use arrays, only: zw,gcoeff, pt, qt, fr, fi, un, unr

      implicit real*8 (a-h,o-z)
      real*8  blength(nsub)
      real*8  bnodes(nsub+1)
      real *8  u(nsub*k,nth),ur(nsub*k,nth),uth(nsub*k,nth)
      real *8  urr(nsub*k,nth),urt(nsub*k,nth), utt(nsub*k,nth)
      real *8  rnd(nsub*k),f(nsub*k,nth),g(nth)
 
   
      real *8 wsave(6*nth+100)
c
c     get arrays for FFTs, ODE solver, etc.
c
c     determine length of work array zw
c     allocate zw, wsave for ffts, wbvp/iwbvp for bvp solver
c
c     isolve   = flag specifying which PDE is to be solved.
c     alpha    = Helmholtz coefficient for isolve=0,1,2

      isolve = isolver !from namelist
      alpha  = halpha  !from namelist

      nr = nsub*k
      ntot = nr*nth
      ifcoeff = 1 
      iu = ifcoeff + ntot
      iur = iu + ntot
      iurr = iur + ntot
      izc = iurr + ntot
      izcp = izc + nth
      itot = izcp + nth
c
      write (*,*) 'zw',itot,nth,nth,nr
c
c     compute FFT of the rhs of PDE at each chebychev node
c
      call fftrhs(nsub,k,nr,nth,rnd,f,zw(ifcoeff),
     *                zw(izc),ibc,g,gcoeff,wsave)
c
c   solve ODE in r for each fourier mode
c
         call ellode(nsub,k,blength,
     *                bnodes,nr,nth,rnd,isolve,alpha,pt,qt,fr,fi,
     *                un,unr,zw(ifcoeff),gcoeff,zw(iu),zw(iur),zw(iurr))

         ucoeff(1:ntot)=zw(iu:iu+ntot-1)  !saved coefficient for interpolation off grid points
         urcoeff(1:ntot)=zw(iur:iur+ntot-1)
         urrcoeff(1:ntot)=zw(iurr:iurr+ntot-1)
c
c   do a backwards FFT on ucoeff and urcoeff to get u, ur, uth and utt
c
         call fftsoln(nsub,k,nr,nth,rnd,zw(iu),zw(iur),u,ur,urt,
     *                 uth,utt,zw(izc),zw(izcp),wsave)
         

c   calculate urr
         do l = 1,nr
            r = rnd(l)
            urr(l,:) = f(l,:)-ur(l,:)/r-utt(l,:)/r**2
c       write(*,*) 'f,u,ur,r,utt,urr'
c     1           ,f(l,1),u(l,1),ur(l,1),r,utt(l,1),urr(l,1)
         end do
      return
      end 
c
c
c***********************************************************************
      subroutine fftrhs(nsub,k,nr,nth,rnd,f,fcoeff,zc,ibc,
     *                   g,gcoeff,wsave)
c***********************************************************************
c     This subroutine computes and stores the Fourier transform of 
c     rhs of PDE.
c
c     For each fixed radius - compute FFT in theta.
c
c     INPUT:
c
c     nsub     = number of subintervals in r on [0,R].
c     k        = number of (scaled) Chebyshev nodes used on each subinterval
c     nr       = number of points in r discretization.
c     nth      = number of points in theta discretization.
c     rnd      = (interior) grid pts in r - does not include r=RAD.
c     f        = right-hand side of PDE
c     zc       = workspace
c     ibc      = boundary condition flag. 0= outgoing, 1= Dirichlet
c     g        = Dirichlet data needed if ibc = 1, otherwise not used.
c     wsave    = FFT workspace
c
c     OUTPUT:
c
c     fcoeff   = Fourier transform in theta of f data
c     gcoeff   = Fourier transform in theta of g data
c
      implicit real*8 (a-h,o-z)
      real *8 rnd(k*nsub),f(nr,nth),g(nth)
      complex *16 zc(nth),fcoeff(nr,nth)
      complex *16 gcoeff(nth)
      complex *16 wsave(*)
c
      pi = 4.d0*datan(1.d0)
c
c     initialize FFT
c
      call dcffti(nth,wsave)
c
c     Fourier transform and scale for each fixed r
c
      do l = 1,nr
         r = rnd(l)
         do j = 1,nth
            zc(j) = f(l,j)
         enddo
         call dcfftf(nth,zc,wsave)
         do j = 1,nth
            fcoeff(l,j) = zc(j)/nth
         enddo
      enddo
c
c     Fourier transform and scale Dirichlet boundary condition
c
      if (ibc.eq.1) then
         do j = 1,nth
            zc(j) = g(j)
         enddo
         call dcfftf(nth,zc,wsave)
         do j = 1,nth
            gcoeff(j) = zc(j)/nth
         enddo
      endif
      return
      end
c
c---------------
      subroutine ellode(nsub,k,blength,bnodes,nr,nth,rnd,isolve,
     *                  alpha,pt,qt,fr,fi,u,ur,
     *                  fcoeff,gcoeff,ucoeff,urcoeff, urrcoeff)
c---------------
c
c     For each fourier mode N, -NFFT <= N <= NFFT, solve the following 
c     BVP:
c
c     If ISOLVE = 0:
c      (u_n)'' + (u_n)'/r - n^2 (u_n)/r^2 = F_n(r)
c           u_n (0) = 0 for n<>0  (u_0)' (0) = 0
c           u_n (1) = 0
c
c     If ISOLVE = 1: 
c      (u_n)'' + (u_n)'/r - (alpha^2 + n^2/r^2) u_n = F_n(r)
c           u_n (0) = 0 for n<>0  (u_0)' (0) = 0
c           u_n (1) = 0
c
c     If ISOLVE = 2: 
c      (u_n)'' + (u_n)'/r + (alpha^2 - n^2/r^2) u_n = F_n(r)
c           u_n (0) = 0 for n<>0  (u_0)' (0) = 0
c           u_n (1) = 0
c
      implicit real*8 (a-h,o-z)
      real *8 rnd(k*nsub), pt(k*nsub), qt(k*nsub), fr(k*nsub),
     *          fi(k*nsub),u(k*nsub), ur(k*nsub) 
      complex*16 ucoeff(nr,nth),urcoeff(nr,nth),fcoeff(nr,nth),eye
      complex*16 urrcoeff(nr,nth), gcoeff(nth)
      real *8 blength(nsub)
      real *8 bnodes(nsub+1)
c
      eye = dcmplx(0.d0,1.d0)
c
c     depending on parity, carry out fftshift procedure:
c
c     if (n.eq.6), for example, the FTT produces six numbers:
c     uhat(1,2,3,4,5,6)  <->  uhat(0,1,2,3,-2,-1) so that the 4th
c     component should be ignored since its conjugate pair is missing
c     for a real-valued function.
c     if (n.eq.5), the FTT produces five numbers:
c     uhat(1,2,3,4,5)  <->  uhat(0,1,2,-2,-1) so that all
c     components can be used.
c
      if (mod(nth,2).eq.0) then
         nfft = (nth-2)/2
         nfftr = nfft+1
      else
         nfft = (nth-1)/2
         nfftr = nfft
      end if
      do n = -nfft,nfftr
         if (n.lt.0) then
            jfcoeff = nth+n+1
         else
            jfcoeff = n+1
         end if
c
c     Set the coefficients of the ODE according to:
c     u'' + pt(x) u'(x) + q(x) u(x) = f(x)
c
      do l = 1 , k*nsub
         r = rnd(l)
         pt(l) = 1d0 / r
         fr(l) = dreal(fcoeff(l,jfcoeff))
         fi(l) = dimag(fcoeff(l,jfcoeff))
         if (isolve.eq.0) then
            qt(l) = -n**2 / r**2
         else if (isolve.eq.1) then
            qt(l) = -(alpha**2 + n**2 / r**2)
         else if (isolve.eq.2) then
            qt(l) = (alpha**2 - n**2 / r**2)
         end if
      end do
c
c     Set coefficients of boundary condition at origin according to:
c     za1 u'(0) + za0 u(0) = 0
c
      a = 0.d0
      if (n.eq.0) then
         za1 = 1.d0
         za0 = 0.d0
      else
         za1 = 0.d0
         za0 = 1.d0
      end if
      ga = 0.d0
      c = 1.d0
      zc1 = 0.d0
      zc0 = 1.d0
c
c     Solve BVP for the real part of the fourier coefficient
c
      gc = dreal(gcoeff(jfcoeff))
      call BVP (nsub,k,a,za1,za0,ga,c,zc1,zc0,gc,rnd,pt,qt,fr,
     *          bnodes,u,ur,ura,urc)
c
c     Store fourier coefficients in 2d array
c
      do l = 1 , k*nsub
         ucoeff(l,jfcoeff) = u(l)
         urcoeff(l,jfcoeff) = ur(l)
         urrcoeff(l,jfcoeff) = fr(l)-pt(l)*ur(l)-qt(l)*u(l)  
       
      end do
c      write(*,*) 'inext55before',  urcoeff(nr,jfcoeff)
c      urcoeff(nr,jfcoeff) = urc  !JPLEE: removed due to dirichlet B.C. 
c       write(*,*) 'inext55',  urcoeff(nr,jfcoeff)
c
c     Solve BVP for the imaginary part of the fourier coefficient
c
      gc = dimag(gcoeff(jfcoeff))
      call BVP (nsub,k,a,za1,za0,ga,c,zc1,zc0,gc,rnd,pt,qt,fi,
     *          bnodes,u,ur,ura,urc)
c
c     Increment fourier coefficients in 2D array
c
      do l = 1 , k*nsub
         ucoeff(l,jfcoeff) = ucoeff(l,jfcoeff) + eye*u(l)
         urcoeff(l,jfcoeff) = urcoeff(l,jfcoeff) + eye*ur(l)
         urrcoeff(l,jfcoeff) = urrcoeff(l,jfcoeff) +
     1        eye*(fi(l)-pt(l)*ur(l)-qt(l)*u(l))  
      end do
c      urcoeff(nr,jfcoeff) = urcoeff(nr,jfcoeff) + eye*urc !JPLEE: removed due to dirichlet B.C. 
      end do
      return
      end
c
c---------------
      subroutine fftsoln(nsub,k,nr,nth,rnd,ucoeff,urcoeff,
     *                   u,ur,urt,uth,utt,zc,zcp,wsave)
c---------------
c
c     Compute inverse FFT of ucoeff and urcoeff to get the solution
c     U and its radial derivative UR.
c
      implicit real*8 (a-h,o-z)
      real *8 rnd(k*nsub),u(nr,nth),ur(nr,nth),urt(nr,nth)
      real *8 uth(nr,nth),utt(nr,nth)
      complex *16 wsave(*)
      complex*16 ucoeff(nr,nth),urcoeff(nr,nth),zc(nth)
      complex*16 zcp(nth),zcpp(nth)
c
      nr = k*nsub
      do l = 1,nr
         r = rnd(l)
         do j = 1,nth
            zc(j) = ucoeff(l,j)
         end do

         if(l.eq.nr) then
c            write(*,*) 'r,ucoeff',ucoeff
         end if
         
c
c     Note: fdiffc takes Fourier transform of function and returns
c     derivative values of function itself (i.e. it incorporates
c     dcfftb).
c
         call fdiffc(zc,zcp,nth,wsave)
         call fdiff2c(zc,zcpp,nth,wsave)
         call dcfftb(nth,zc,wsave)
         do j = 1,nth
            u(l,j) = dreal(zc(j))
            uth(l,j) = dreal(zcp(j))
            utt(l,j) = dreal(zcpp(j))
            zc(j) = urcoeff(l,j)
         end do
  
         call fdiffc(zc,zcp,nth,wsave)
         call dcfftb(nth,zc,wsave)
         do j = 1,nth
         
            urt(l,j) = dreal(zcp(j))
            ur(l,j) = dreal(zc(j))
c            inext = l + (j-1)*nr
c         write(*,*) 'inext5', inext, ur(l,j), uth(l,j)
         end do
    
      end do
c
      return
      end

c---------------
      subroutine ftsolnoffth(thin,nth,ucin,urcin,urrcin,
     *                   uout,urout,urtout,uthout,uttout,urrout)
c---------------
c
c     Compute inverse Fourier transform of ucoeff,urcoeff and urrcoeff at a specific r
c     to get the solution at a specific theta
c
      implicit real*8 (a-h,o-z)
      real *8 thin,uout,urout,urtout,urrout
      real *8 uthout,uttout
      complex *16 zout,zpout,zppout
      complex*16 ucin(nth),urcin(nth),urrcin(nth)
      complex*16 zc(nth),zcp(nth),zcpp(nth)


      zc(1:nth)=ucin(1:nth)

  
      call fcoffth(zc,nth,thin,zout)
c      write(*,*) 'zc',thin,zout,zc(1:nth)
      call fdiffcoffth(zc,zcp,nth,thin,zpout)
      call fdiff2coffth(zc,zcpp,nth,thin,zppout)

      uout = dreal(zout)
      uthout = dreal(zpout)
      uttout = dreal(zppout)

      zc(1:nth)=urcin(1:nth)
    
      call fcoffth(zc,nth,thin,zout)
      call fdiffcoffth(zc,zcp,nth,thin,zpout)

      urout = dreal(zout)
      urtout = dreal(zpout)

      zc(1:nth)=urrcin(1:nth)
    
      call fcoffth(zc,nth,thin,zout)

      urrout = dreal(zout)

c      write(*,*) 'u, ur, ut2',uout, urout, uthout
c
      return
      end

c---------------
      subroutine ftsolnoffth0pi(thin,nth,ucin,urcin,urrcin,
     *                   uout,urout,urrout)
c---------------
c
c     Compute inverse Fourier transform of ucoeff,urcoeff and urrcoeff at a specific r
c     to get the solution at a specific theta
c
      implicit real*8 (a-h,o-z)
      real *8 thin,uout,urout,urtout,urrout
      real *8 uthout,uttout
      complex *16 zout,zpout,zppout
      complex*16 ucin(nth),urcin(nth),urrcin(nth)
      complex*16 zc(nth),zcp(nth),zcpp(nth)


      zc(1:nth)=ucin(1:nth)

  
      call fcoffth(zc,nth,thin,zout)

      uout = dreal(zout)

      zc(1:nth)=urcin(1:nth)
    
      call fcoffth(zc,nth,thin,zout)

      urout = dreal(zout)

      zc(1:nth)=urrcin(1:nth)
    
      call fcoffth(zc,nth,thin,zout)

      urrout = dreal(zout)

c      write(*,*) 'u, ur, ut2',uout, urout, uthout
c
      return
      end
