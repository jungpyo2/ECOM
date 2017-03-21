!--------------------------------------------------------
! bvp.f: This subroutine is to solve an ODE using its Greeen function
!--------------------------------------------------------

c---------------
      subroutine bvp (nsub,nnd,a,za1,za0,ga,c,zc1,zc0,gc,xnd,pt,qt,ft,
     *                bnodes,ucom,upcom,uprimea,uprimec)
ccc     *                bnodes,ucom,upcom,uprimea,uprimec,
ccc     *                leafid,leaflist,parentlist)
c---------------
c     Solves ODE using fast direct solver.
c
c     The governing ODE is assumed to be of the form
c      u'' + pt(x) u'(x) + q(x) u(x) = f(x)
c     with boundary conditions at x=a and x=c:
c
c     za1 u'(a) + za0 u(a) = ga.
c     zc1 u'(c) + zc0 u(c) = gc.
c
c     INPUT:
c
c     nsub       = number of subintervals in r on [0,R].
c     nnd        = number of (scaled) Chebyshev nodes used on each subinterval
c     a          = left endpoint
c     za1,za0    = coeffs of linear b.c. at x=a
c     ga         = b.c. value at x=a
c     c          = left endpoint
c     zc1,zc0    = coeffs of linear b.c. at x=c
c     gc         = b.c. value at x=c
c     xnd        = (interior) grid pts - does not include a,c.
c     pt         = value of p(x) at xnd grid points
c     qt         = value of q(x) at xnd grid points
c     ft         = value of f(x) at xnd grid points
c     bnodes     = interval boundaries
c     leafid     = tree data structure for subintervals in r.
c     leaflist   = tree data structure for subintervals in r.
c     parentlist = tree data structure for subintervals in r.
c
c     OUTPUT:
c
c     ucom       = value of solution at xnd grid points
c     upcom      = value of derivative of solution at xnd grid points
c     uprimea    = value of u' at endpoint x=a
c     uprimec    = value of u' at endpoint x=c
c
      implicit real*8 (a-h,o-z)
ccc      integer leafid(nsub), leaflist(nsub),parentlist(nsub)
      real *8 xnd(nsub*nnd),pt(nsub*nnd),qt(nsub*nnd),ft(nsub*nnd),
     *          ucom(nsub*nnd),upcom(nsub*nnd),mergcond
      real *8 bnodes(nsub+1)
c
      integer, allocatable :: leafid(:)
      integer, allocatable :: leaflist(:)
      integer, allocatable :: parentlist(:)
c
      real *8, allocatable :: blength(:)
      real *8, allocatable :: xndsort(:)
      real *8, allocatable :: ptsort(:)
      real *8, allocatable :: qtsort(:)
      real *8, allocatable :: ftsort(:)
      real *8, allocatable :: ucomsort(:)
      real *8, allocatable :: upcomsort(:)
      real *8, allocatable :: gl(:)
      real *8, allocatable :: gr(:)
      real *8, allocatable :: phil(:)
      real *8, allocatable :: phir(:)
      real *8, allocatable :: eta(:)
      real *8, allocatable :: sigma(:)
c
      allocate(leafid(nsub))
      allocate(leaflist(nsub))
      allocate(parentlist(nsub))
c
      call gentree(nsub, leafid, leaflist, parentlist)
c
      npnts = nsub*nnd
      allocate(gl(npnts))
      allocate(gr(npnts))
      allocate(phil(npnts))
      allocate(phir(npnts))
      allocate(eta(npnts))
      allocate(sigma(npnts))
      allocate(xndsort(npnts))
      allocate(ptsort(npnts))
      allocate(qtsort(npnts))
      allocate(ftsort(npnts))
      allocate(ucomsort(npnts))
      allocate(upcomsort(npnts))
      allocate(blength(nsub))
c
c     various arrays needed by ODE solver.
c     blength    = interval lengths corresponding to leafid ordering
c
c     reorder data to correspond to subinterval ordering dictated by leafid.
c
      call bcond_lin(a, c, za1, za0, zc1, zc0, ga, gc)
      do kp = 1, nsub
         j = leafid(kp)
         blength(j) = bnodes(kp+1) - bnodes(kp)
         istart = nnd*(j-1)
         istartsort = nnd*(kp-1)
         do ii=1,nnd
            xndsort(istart+ii) = xnd(istartsort+ii)
            ptsort(istart+ii) = pt(istartsort+ii)
            qtsort(istart+ii) = qt(istartsort+ii)
            ftsort(istart+ii) = ft(istartsort+ii)
         enddo
      end do
c
c -------------------------------
c     SOLVE ODE
c -------------------------------
         do kp = 1, nsub
            j = leafid(kp)
            call int_eqn(j,nnd,a,c,xndsort,ptsort,qtsort,ftsort,gl,gr,
     1             phil,phir,eta)
         end do
         call solve_sigma(nnd, nsub, leaflist, leafid, nsub, 
     &        leaflist, leafid, parentlist,
     &        blength,gl,gr,phil,phir,eta, mergcond, discond, sigma)
c
c
c   set output parameter ieval (I can't find documentation, so
c   leave this set to 3, which seems to make sure that both
c   u and u' are returned at Chebyshev nodes.
c
         ieval = 3
         call solve_eval(ieval, nnd, nsub, leaflist, leafid,
     &        blength,xndsort,gl,gr, sigma,upcomsort, ucomsort)
c
c     Resort solution into linearly ordered arrays.
c
      do kp = 1, nsub
         j = leafid(kp)
         istart = nnd*(j-1)
         istartsort = nnd*(kp-1)
         do ii=1,nnd
            ucom(istartsort+ii) = ucomsort(istart+ii)
            upcom(istartsort+ii) = upcomsort(istart+ii)
         enddo
      end do
      return
      end
c
***********************************************************************
c
      subroutine chsetup(nnd)

c     

      integer nnd, nmxnd
      parameter (nmxnd=256)
      real*8 pi
      parameter (pi=3.141592653589793238d0)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
      real*8 cftm, citm, spdef, spbx, spxb
      common /chmat/ cftm(nmxnd,nmxnd), citm(nmxnd,nmxnd)
      common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
c
c------------------------------------------------------------
c input  :
c    nnd = number of nodes specified by a user
c output :
c    theta = Chebychev nodes on [0,pi] were stored in decreasing order
c            so that corresponding nodes in [-1,1] are in increasing order
c    chnd01 = Chebyshev nodes between [0,1]
c    cftm : ch(1:n) = cftm(1:n,1:n) * fn(1:n)
c    citm : fn(1:n) = citm(1:n,1:n)' * fn(1:n)
c    SPBX sigma  = Cn^{-1} SPINT_BXn Cn = \INT_{-1}^{x} sigma(t) dt
c    SPXB sigma  = Cn^{-1} SPINT_XBn Cn = \INT_{x}^{1} sigma(t) dt
***********************************************************************
        integer i,k
        real*8 work(nmxnd), defint

        do i = 1, nnd
           theta(i) = (nnd-i+.5d0)/nnd * pi
           ns(i) = dsin(theta(i)) * dsqrt(2d0)
           c(i) = dcos(theta(i))
           chnd01(i) = ( dcos(theta(i)) + 1d0 ) / 2d0
        end do
ccc        write(6,*) ' chnodes are ',(chnd01(i),i=1,nnd)

c------------------------------------------------------------
c computes terms in cosine inverse transform
c and stores them into matrix citm in the transposed format
c to facilitate data access during matrix multiplications.
c------------------------------------------------------------
        do  k = 1,nnd
           do i = 1, nnd
              citm(i,k) = dcos((i-1)*theta(k))
           end do
        end do

c------------------------------------------------------------
c computes terms in cosine forward transform with scaling factors
c and stores them into matrix cftm.
c    cftm(i,k) = 2d0/dble(nnd) * dcos((i-1)*theta(k))
c------------------------------------------------------------
        do k = 1,nnd
           cftm(1,k) = 1d0/dble(nnd)
           do i = 2, nnd
              cftm(i,k) = 2d0/dble(nnd)*citm(i,k)
c              cftmsub((i-1)*nnd+k) = cftm(i,k)
            end do
        end do

c------------------------------------------------------------
c Spectral integration matrix
c   ctfm(:,i) = Cn En where En = ( 0, 0, ... , 1=i-th, ... , 0 )
c   SPBX Sigma  = Cn^{-1} SPINT_BXn Cn = \INT_{-1}^{x} Sigma(t) dt
c   SPXB Sigma  = Cn^{-1} SPINT_XBn Cn = \INT_{x}^{1} Sigma(t) dt
c------------------------------------------------------------
        do i = 1, nnd
           call chindef_bx(work,cftm(1,i),nnd)
           call chbtrans(spbx(1,i),work,nnd)
           call chindef_xb(work,cftm(1,i),nnd)
           call chbtrans(spxb(1,i),work,nnd)
           spdef(i) = defint(cftm(1,i),nnd)
        end do
      end
       
***********************************************************************
      subroutine chnodes(nnd,j,b,blength,xnd)
      integer nmxnd, nmxsub, nnd, j
      parameter (nmxnd=256,nmxsub=1024)
      real*8 xnd(nnd), b(nmxsub),blength(nmxsub)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
c --------------------------------------------------------------------
c return the scaled chebychev nodes on a subinterval [b(j),b(j)+blength(j)]
c input : b(j), blength(j)
c common : chnd01(nnd)
c output : xnd(nnd)
***********************************************************************
      integer i, l
        do i = 1, nnd
           l = i + nnd*(j-1)
           xnd(l) = b(j) + blength(j) * chnd01(i)
        end do
      end

*************************************************************************
      subroutine chftrans(ch, fn, nnd)
      integer nmxnd, nnd, i, j
      parameter (nmxnd=256)
      real*8 ch(nnd), fn(nnd)
      real*8 cftm, citm
      common /chmat/ cftm(nmxnd,nmxnd), citm(nmxnd,nmxnd)
c------------------------------------------------------------------------
c  Given function values at the chebyshev nodes (physical space)
c  calculate the corresponding chebyshev coefficients (fourier space).
*************************************************************************
        do i = 1, nnd
           ch(i) = 0d0
           do j = 1, nnd
              ch(i) = ch(i) + cftm(i,j)*fn(j) 
           end do
        end do
      end

*********************************************************************
      subroutine chbtrans(fn, ch, nnd)
      integer nmxnd, nnd, i, j
      parameter (nmxnd=256)
      real*8 ch(nnd), fn(nnd)
      real*8 cftm, citm
      common /chmat/ cftm(nmxnd,nmxnd), citm(nmxnd,nmxnd)
c------------------------------------------------------------------------
c  from the chebyshev coefficients to function values at the chebyshev nodes
*************************************************************************
        do i = 1, nnd
           fn(i) = 0d0
           do j = 1, nnd
              fn(i) = fn(i) + citm(j,i)*ch(j)
           end do
        end do
      end

*********************************************************************
      function defint(coef,nnd)
      integer  nnd
      real*8 defint,coef(0:nnd-1)
*********************************************************************
      integer i
        defint =  2d0 * coef(0)
        do i = 2, nnd-1, 2
           defint = defint - 2d0 * coef(i) / dble(i*i-1)
        end do
      end
	
*********************************************************************
      subroutine chindef_bx(chintfl,coeffl,nnd)
      integer i,nnd
      real*8 coeffl(0:nnd-1), chintfl(0:nnd-1)
c  ------------------------------------------------------------------
c this routine calculates the chebyshev coefficients of
c the left (normal) indefinite integrals
*********************************************************************

        chintfl(1) = coeffl(0) - coeffl(2)/2d0
        do i = 2,nnd-2
           chintfl(i) = (coeffl(i-1)-coeffl(i+1)) / dble(2*i)
        end do
        chintfl(nnd-1) = coeffl(nnd-2) / dble(2*nnd-2)

        chintfl(0) = (-1)**(nnd+1) * coeffl(nnd-1) / dble(2*nnd)
        do i = nnd-1, 1, -1
           chintfl(0) = chintfl(0) - (-1)**i * chintfl(i)
        end do
      end

*********************************************************************
      subroutine chindef_xb(chintfr,coeffr,nnd)
      integer i,nnd
      real*8 coeffr(0:nnd-1),chintfr(0:nnd-1)
c  ------------------------------------------------------------------
c this routine calculates the chebyshev coefficients of
c the right (backward) indefinite integrals
*********************************************************************

        chintfr(1) = coeffr(2)/2d0 - coeffr(0)
        do i = 2, nnd-2
           chintfr(i) = (coeffr(i+1)-coeffr(i-1)) / dble(2*i)
        end do
        chintfr(nnd-1) = -coeffr(nnd-2) / dble(2*nnd-2)

        chintfr(0) = coeffr(nnd-1) / dble(2*nnd)
        do i = 1, nnd-1
           chintfr(0) = chintfr(0) -  chintfr(i)
        end do
      end
c
c
c---------------
      subroutine chrad (nsub,nnd,a,c,blength,bnodes,rnd,delrmin)
c---------------
c
c     Computes subinterval boundaries and radial grid points
c     on composite Chebyshev mesh.
c
c     INPUT:
c
c     nsub       - number of subintervals
c     nnd        - number of Chebyshev nodes on each subinterval
c     a          - left endpoint
c     c          - right endpoint
c
c     OUTPUT:
c
c     blength    - blength(i) is length of ith subinterval
c     bnodes     - bnodes(1)=a,...bnodes(nsub+1)=c  are the 
c                  subinterval boundary points
c     rnd        - rnd(i) is ith grid point in composite grid
c                - rnd(nsub*nnd+1) is set to c.
c
      implicit real*8 (a-h,o-z)
      real *8 blength(nsub)
      real *8 bnodes(nsub+1)
      real *8 rnd(nsub*nnd+1)

      save
c
c     chsetup sets up important arrays in comon blocks of ODE solver.
c
      write(*,*) 'beforechsetup'
      call chsetup(nnd)
      write(*,*) 'afterchsetup'
c
c     ------------------------------------
c     generate subinterval boundaries
c     ------------------------------------
c
      bnodes(1) = a
      bnodes(nsub+1) = c
      nweight=4
      if(nsub .gt. 800) then
         bnodes(2) = bnodes(1) + (c-a)/dble(10*nweight*nsub)
         bnodes(3) = bnodes(2) + (c-a)/dble(nweight*nsub)
         bnodes(4) = bnodes(3) + (c-a)/dble(nweight*nsub)
         bnodes(nsub) = bnodes(nsub+1)-(c-a)/dble(10*nweight*nsub)
         bnodes(nsub-1) = bnodes(nsub)-(c-a)/dble(nweight*nsub)
         bnodes(nsub-2) = bnodes(nsub-1)-(c-a)/dble(nweight*nsub)
         do kp = 5, nsub-3
            bnodes(kp) = bnodes(4) + 
     1           (kp-4)*(bnodes(nsub-2)-bnodes(4))/dble(nsub-6)
         end do
      else 
         do kp = 1, nsub
            bnodes(kp) = a + (kp-1)*(c-a)/dble(nsub)
         end do
         bnodes(2) = (bnodes(1) + bnodes(2))/2.0d0
      end if

      do kp = 1, nsub
         blength(kp) = bnodes(kp+1) - bnodes(kp)
      end do
c
c     ------------------------------------
c      generate composite chebychev mesh
c     ------------------------------------ 
c
      do kp = 1, nsub
         j = kp
         call chnodes(nnd,j,bnodes,blength,rnd)
      end do
      rnd(nsub*nnd+1) = c
c
c     ------------------------------------
c     calculate the minimum mesh spacing
c     ------------------------------------
c
      delrmin = rnd(2)-rnd(1)
      do i = 1,nsub*nnd
         dr = rnd(i+1)-rnd(i)
         delrmin = min(delrmin,dr)
      end do
      return
      end
 
      subroutine chsetupsub(nnd, psiich, cftm)
      integer nnd, nmxnd
      parameter (nmxnd=256)
      real*8 pi
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      real*8 psiich(*), cftm(*), citm(nmxnd,nmxnd)
      
      integer i,k
      real*8 work(nmxnd), defint


      pi=4.0d0*datan(1.0d0)
      do i = 1, nnd
         theta(i) = (nnd-i+.5d0)/nnd * pi
         ns(i) = dsin(theta(i)) * dsqrt(2d0)
         c(i) = dcos(theta(i))
         chnd01(i) = ( dcos(theta(i)) + 1d0 ) / 2d0
         
         psiich(i)=(1.0d0-chnd01(i))
      end do
cc        write(6,*) ' chnodes are ',(chnd01(i),i=1,nnd)

c------------------------------------------------------------
c computes terms in cosine inverse transform
c and stores them into matrix citm in the transposed format
c to facilitate data access during matrix multiplications.
c------------------------------------------------------------
      do  k = 1,nnd
         do i = 1, nnd
            citm(i,k) = dcos((i-1)*theta(k))
         end do
      end do

c------------------------------------------------------------
c computes terms in cosine forward transform with scaling factors
c and stores them into matrix cftm.
c    cftm(i,k) = 2d0/dble(nnd) * dcos((i-1)*theta(k))
c------------------------------------------------------------
      do k = 1,nnd
         cftm((k-1)*nnd+1) = 1d0/dble(nnd)
         do i = 2, nnd
            cftm((k-1)*nnd+i) = 2d0/dble(nnd)*citm(i,k)
         end do
      end do
 
      return
      end 

      subroutine chsetupq(nnd, psiichq, cftm, spbx, spxb, spdef)
      integer nnd, nmxnd
      parameter (nmxnd=256)
      real*8 pi
c      parameter (pi=3.141592653589793238d0)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      real*8 psiichq(nnd), cftm(nnd), spbx(nnd)
      real*8 spxb(nnd*nnd), spdef(nnd*nnd)
      real*8 citm(nmxnd,nmxnd)
      real*8 work1(nmxnd), work2(nmxnd)
c      common /geofn/ theta, ns, c, chnd01
c      real*8 cftm, citm, spdef, spbx, spxb
c      common /chmat/ cftm(nmxnd,nmxnd), citm(nmxnd,nmxnd)
c      common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
c
c------------------------------------------------------------
c input  :
c    nnd = number of nodes specified by a user
c output :
c    theta = Chebychev nodes on [0,pi] were stored in decreasing order
c            so that corresponding nodes in [-1,1] are in increasing order
c    chnd01 = Chebyshev nodes between [0,1]
c    cftm : ch(1:n) = cftm(1:n,1:n) * fn(1:n)
c    citm : fn(1:n) = citm(1:n,1:n)' * fn(1:n)
c    SPBX sigma  = Cn^{-1} SPINT_BXn Cn = \INT_{-1}^{x} sigma(t) dt
c    SPXB sigma  = Cn^{-1} SPINT_XBn Cn = \INT_{x}^{1} sigma(t) dt
***********************************************************************
        integer i,k
        real*8 work(nmxnd), defint

        save

        pi=4.0d0*datan(1.0d0)
        do i = 1, nnd
           theta(i) = (nnd-i+.5d0)/nnd * pi
           ns(i) = dsin(theta(i)) * dsqrt(2d0)
           c(i) = dcos(theta(i))
           chnd01(i) = ( dcos(theta(i)) + 1d0 ) / 2d0
           
           psiichq(i)=(1.0d0-chnd01(i))
        end do
c        write(*,*) 'psiichq0', psiichq(1:nnd)
cc        write(6,*) ' chnodes are ',(chnd01(i),i=1,nnd)

c------------------------------------------------------------
c computes terms in cosine inverse transform
c and stores them into matrix citm in the transposed format
c to facilitate data access during matrix multiplications.
c------------------------------------------------------------
        do  k = 1,nnd
           do i = 1, nnd
              citm(i,k) = dcos((i-1)*theta(k))
           end do
        end do

c------------------------------------------------------------
c computes terms in cosine forward transform with scaling factors
c and stores them into matrix cftm.
c    cftm(i,k) = 2d0/dble(nnd) * dcos((i-1)*theta(k))
c------------------------------------------------------------
        do k = 1,nnd
           cftm((k-1)*nnd+1) = 1d0/dble(nnd)
           do i = 2, nnd
              cftm((k-1)*nnd+i) = 2d0/dble(nnd)*citm(i,k)
            end do
        end do
c------------------------------------------------------------
c Spectral integration matrix
c   ctfm(:,i) = Cn En where En = ( 0, 0, ... , 1=i-th, ... , 0 )
c   SPBX Sigma  = Cn^{-1} SPINT_BXn Cn = \INT_{-1}^{x} Sigma(t) dt
c   SPXB Sigma  = Cn^{-1} SPINT_XBn Cn = \INT_{x}^{1} Sigma(t) dt
c------------------------------------------------------------
        do k = 1, nnd
           call chindef_bx(work1,cftm((k-1)*nnd+1),nnd)
           call chindef_xb(work2,cftm((k-1)*nnd+1),nnd)
           do i = 1, nnd
              spbx((k-1)*nnd+i) = 0d0
              spxb((k-1)*nnd+i) = 0d0
              do j = 1, nnd
                 spbx((k-1)*nnd+i)=spbx((k-1)*nnd+i)+citm(j,i)*work1(j)
                 spxb((k-1)*nnd+i)=spxb((k-1)*nnd+i)+citm(j,i)*work2(j)
              end do
           end do
           spdef(k) = defint(cftm((k-1)*nnd+1),nnd)
        end do

        return
        end 


      subroutine chftransq(ch, fn, nnd, cftm)

      integer nmxnd, nnd, i, j
      parameter (nmxnd=256)
      real*8 ch(*), fn(*), cftm(*)

      save
c      real*8 cftm, citm
c      common /chmat/ cftm(nmxnd,nmxnd), citm(nmxnd,nmxnd)
c------------------------------------------------------------------------
c  Given function values at the chebyshev nodes (physical space)
c  calculate the corresponding chebyshev coefficients (fourier space).
*************************************************************************
      do i = 1, nnd
         ch(i) = 0d0
         do j = 1, nnd
            ch(i) = ch(i) + cftm((j-1)*nnd+i)*fn(j) 
         end do
      end do
c      write(*,*) 'fn',fn(1:nnd)
c      write(*,*) 'ch',ch(1:nnd)

      return
      end

      subroutine chderiv(x, nnd, ch, dydx)

      implicit real*8 (a-h,o-z)
      integer nmxnd, nnd, i, j
      parameter (nmxnd=256)
      real*8 ch(*)

      save
c      real*8 cftm, citm
c      common /chmat/ cftm(nmxnd,nmxnd), citm(nmxnd,nmxnd)
c------------------------------------------------------------------------
c  Calculate the first derivative of the chebyshev polynomials at x
c  dTn/dx=n/(1-x**2)*(Tn+1-x*Tn)
c  Using the given chebyshev coefficient,
c  df/dx=sum_{n=1}^{n=nnd} ch_n* dTn/dx
*************************************************************************
     
      xx=x*2.0d0-1.0d0
      dydx=0.0d0
    
      if (xx.eq.1.0d0) then
         do i=1,nnd
            dTndx=(i-1)**2
            dydx=dydx+ch(i)*dTndx
         end do
         dydx = -dydx
      else if (xx.eq.-1.0d0) then
         do i=1,nnd
            dTndx=(-1)**(i)*(i-1)**2
            dydx=dydx+ch(i)*dTndx
         end do
         dydx = -dydx
      else
         Tn=1
         do i = 1, nnd
            Tnp1=dcos((i)*dacos(xx))
            dTndx=(i-1)*(Tnp1-xx*Tn)
            dydx=dydx+ch(i)*dTndx
            Tn=Tnp1
         end do
         dydx=dydx/(1-xx**2)
      end if
      dydx=dydx*2 !scaling of x from [-1.1] to [0,1] 
c      write(*,*) x,dydx
      return
      end
   
      subroutine chfit(x, nnd, ch, y)

      implicit real*8 (a-h,o-z)
      integer nmxnd, nnd, i, j
      parameter (nmxnd=256)
      real*8 ch(nnd), x,xx, y, Tn
    
      xx=x*2.0d0-1.0d0
      y=0.0d0
      do i = 1, nnd
         Tn=dcos((i-1)*dacos(xx))
         y=y+ch(i)*Tn
      end do
      return
      end

**************************************************************************
      subroutine discret(nnd,j,blength,v_l,v_r,localcond,eta,phil,phir)
      integer nmxnd, nmxsub, nnd, j
      parameter (nmxnd=256,nmxsub=1024)
      real*8 blength(nmxsub), localcond(nmxsub)
      real*8 v_l(*), v_r(*), phil(*), phir(*), eta(*)
      real*8 spdef, spbx, spxb
      common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
c------------------------------------------------------------
c this subroutine implements the
c Step 1 of the adaptive ode solver by Lee and Greengard, 1995 SISC 
c------------------------------------------------------------
c for each subinterval at the finest level, it sets up a nnd*nnd
c matrix that is derived from discretizing the integral equation 
c solution. then linpack subroutines are used to solve three linear
c systems with the matrix.
c------------------------------------------------------------
c input  : nmxnd = number of maximum nodes per subinterval
c          nnd = number of nodes specified by user
c          v_l, v_r, blength
c in/out : eta = solution with the right hand side f
c          phil = solution with the right hand side psil
c          phir = solution with the right hand side phir
c output : localcond = inverse condition number of local integral operator
**************************************************************************
      real*8 pcc(nmxnd,nmxnd), z(nmxnd)
      integer i, ipvt(nmxnd), k, kp, ip

        do i = 1, nnd
           ip = i + nnd*(j-1)
           do k = 1, nnd
              kp = k + nnd*(j-1)
              pcc(k,i) = blength(j)/2d0 * 
     &         ( phil(kp)*spbx(k,i)*v_l(ip)
     &         + phir(kp)*spxb(k,i)*v_r(ip) )
           end do
           pcc(i,i) = 1d0 + pcc(i,i)
        end do

c ----------------------------------------------------------
c sets up the right hand side vectors for f, psil, psir
c    solves P_i(eta_i) = f_i
c    solves P_i(phil_i) = psil_i
c    solves P_i(phir_i) = phir_i
c and copy it to the right place
c ----------------------------------------------------------
c       call dgefa(pcc,nmxnd,nnd,ipvt,info)
        call dgeco(pcc,nmxnd,nnd,ipvt,localcond(j),z)
        call dgesl(pcc,nmxnd,nnd,ipvt,eta(1+nnd*(j-1)),0)
        call dgesl(pcc,nmxnd,nnd,ipvt,phil(1+nnd*(j-1)),0)
        call dgesl(pcc,nmxnd,nnd,ipvt,phir(1+nnd*(j-1)),0)
      end

****************************************************************************
      subroutine mkcoef(nnd,k,j,blength,v_l,v_r,eta,phil,phir,
     &            alpha11,alpha21,alpha12,alpha22,delta_l,delta_r)
      integer nmxnd, nmxsub, nnd, k, j
      parameter (nmxnd=256,nmxsub=1024)
      real*8 blength(nmxsub)
      real*8 v_r(*), v_l(*), phil(*),phir(*),eta(*)
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      real*8 spdef, spbx, spxb
      common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
c------------------------------------------------------------
c this subroutine implements the
c Step 1 of the adaptive ode solver by Lee and Greengard, 1995 SISC 
c------------------------------------------------------------
c it uses chebychev quadrature method to evaluate inner products on
c each subinterval
c------------------------------------------------------------
c input  : nmxnd = number of maximum nodes per subinterval
c          nnd = number of nodes specified by user
c          nsub = number of subintervals specified by user
c          k = node number in hierarchy tree
c          j = node id to update
c          blength, v_l, v_r
c          eta = stores solution with f as the right hand side 
c          phil = stores solution with psil as the right hand side 
c          phir= stores solution with psir as the right hand side 
c output : alpha11 = inner product of V_l and phi_l
c          alpha12 = inner product of V_l and phi_r
c          alpha21 = inner product of V_r and phi_l
c          alpha22 = inner product of V_r and phi_r
c          delta_l = inner product of V_l and eta
c          delta_r = inner product of V_r and eta
****************************************************************************
      integer i, l
      real*8 sc

      sc = blength(j) / 2d0
      alpha11(k) = 0d0
      alpha12(k) = 0d0
      alpha21(k) = 0d0
      alpha22(k) = 0d0
      delta_l(k) = 0d0
      delta_r(k) = 0d0

      do i = 1, nnd
         l = i + nnd*(j-1)
         alpha11(k) = alpha11(k) + v_l(l)* phil(l)*spdef(i)*sc
         alpha12(k) = alpha12(k) + v_l(l)* phir(l)*spdef(i)*sc
         alpha21(k) = alpha21(k) + v_r(l)* phil(l)*spdef(i)*sc
         alpha22(k) = alpha22(k) + v_r(l)* phir(l)*spdef(i)*sc
         delta_l(k) = delta_l(k) + v_l(l)* eta(l) *spdef(i)*sc
         delta_r(k) = delta_r(k) + v_r(l)* eta(l) *spdef(i)*sc
      end do
      end
**********************************************************************
      subroutine bcond_lin(a,c, za1,za0,zc1,zc0, ga,gc)
      real*8 a, c, za1, za0, zc1, zc0, ga, gc
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c  --------------------------------------------------------------
c  Boundary conditions (all parameter in common /bc/)
c         should be specfied in the subroutine
c  [a,c] : two coordinates of boundary points
c         za1 u'(a) + za0 u(a) = ga
c         zc1 u'(c) + zc0 u(c) = gc
c  G(x,t) = gr(x) gl(t) / wron  or  gr(t) gl(x) / wron,
c      where wron(zl,zr) = zl(x) zrp(x) - zlp(x) zr(x) = constant
c  ---------------------------
c  How to choose Green's functions:
c
c  If za0 * zc0 * (c-a) + ( zc1*za0 - za1*zc0 ) > 1e-4 then use
c    # linear green's functions (q0=0) and linear inhomogeneous term
c    gl = gl1*x + gl0
c    gr = gr1*x + gr0
c    ui = ui1*x + ui0
c
c  otherwise (for Neumann type problems)
c    # cosh,sinh for green's functions (q0=-1) and quadratic u_i(x)
c    gl = gl1 * dcosh(x) - gl0 * dsinh(x)
c    gr = gr1 * dcosh(x) - gr0 * dsinh(x)
c
c    ui = ui2*x**2 + ui1*x + ui0 where
c    ( 2 a za1 + a^2 za0 , za1 + a za0 , za0 ) (ui2)   (ga)
c    (                                       ) (ui1) = (  )
c    ( 2 c za1 + c^2 zc0 , zc1 + c zc0 , zc0 ) (ui0)   (gc)
c
c    - It solves (ui2,ui0) first to get (ui0) and then
c         solves (ui2,ui1) with (ga - za0 ui0, gc - zc0 ui0)
c    - This will be performed only when | det(ui1,ui0) | < 1d-4
c         i.e. (za1 + a za0, zc1 + c zc0) // (za0, zc0)
c      So det(ui2,ui1)=0 implies that
c         (2a za1 + a^2 za0, 2c za1 + c^2 zc0 ) // (za0, zc0)
c         and the rank of the linear system is ONE.
c      Therefore it is solvable only when (ga, gc) is in the range space.
c
c  TO DO -------------------------------
c  1. BETTER WAY TO FIND TO LINEARLY INDEPENDENT GREEN'S FUNCTION
c  2. Using a SVD routine to find a best possible quadratic inhomogeneous term
**********************************************************************
      real*8 det

      wron = za0 * zc0 * (c-a) + ( zc1*za0 - za1*zc0 )
      if ( abs(wron) .ge. 1d-4 ) then
         q0 = 0d0
         gl1 = za0
         gr1 = zc0
         gl0 = -a * za0 - za1
         gr0 = -c * zc0 - zc1

         ui2 = 0d0
         ui1 = ( gc*za0 - ga*zc0 ) / wron
         ui0 = ( (c*zc0+zc1)*ga - (a*za0+za1)*gc ) / wron
         return
      end if

      wron = (za0*zc0-za1*zc1)*dsinh(c-a) + (zc1*za0-za1*zc0)*dcosh(c-a)
      if ( abs(wron) .ge. 1d-4 ) then
         q0 = -1d0
         gl1 = za1*dcosh(a) + za0*dsinh(a)
         gl0 = za1*dsinh(a) + za0*dcosh(a)
         gr1 = zc1*dcosh(c) + zc0*dsinh(c)
         gr0 = zc1*dsinh(c) + zc0*dcosh(c)

         det = za0*zc0*(a*a-c*c) + ( 2d0*a*za1*zc0 - 2d0*c*zc1*za0 )
         if ( abs(det) .lt. 1d-4 ) then
            ui0 = 0d0
         else
            ui0 = ((a*a*za0+2d0*a*za1)*gc-(c*c*zc0+2d0*c*zc1)*ga) / det
         end if
         det = a*c*(a-c)*za0*zc0 + 2d0*(a-c)*za1*zc1
     &       + c*(2d0*a-c)*za1*zc0 + a*(a-2d0*c)*za0*zc1
         ui2 = (  (zc1+c*zc0)*(ga-ui0*za0)
     &          - (za1+a*za0)*(gc-ui0*zc0) ) / det
         ui1 = (  (2d0*a*za1+a*a*za0)*(gc-ui0*zc0)
     &          - (2d0*c*zc1+c*c*zc0)*(ga-ui0*za0) ) / det
         return
      end if

      print*, 'Given boundary conditions are Neumann Type', wron
      print*, 'This mixed boundary condition make some troubles'
      print*, 'to form two linearly independent Green''s functions'
      print*, 'using q0 = 0 or -1. Please make a new choice'
      stop

      end

**********************************************************************
c   gl(x), gr(x), gl'(x), gr'(x) :
c   Green's functions and it's derivatives
c   will be used for function evaluation (int_density.f)
**********************************************************************
      real*8 function glp_x(x)
      real*8 x, q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      if ( q0 .eq. 0 ) then
         glp_x = gl1
      else
         glp_x = gl1 * dsinh(x) - gl0 * dcosh(x)
      end if
      end

      real*8 function grp_x(x)
      real*8 x, q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      if ( q0 .eq. 0 ) then
         grp_x = gr1
      else
         grp_x = gr1 * dsinh(x) - gr0 * dcosh(x)
      end if
      end

      real*8 function gl_x(x)
      real*8 x, q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      if ( q0 .eq. 0 ) then
         gl_x = gl1*x + gl0
      else
         gl_x = gl1 * dcosh(x) - gl0 * dsinh(x)
      end if
      end

      real*8 function gr_x(x)
      real*8 x, q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      if ( q0 .eq. 0 ) then
         gr_x = gr1*x + gr0
      else
         gr_x = gr1 * dcosh(x) - gr0 * dsinh(x)
      end if
      end

**************************************************************************
      subroutine int_eqn(j,nnd,a,c,xnd,pt,qt,ft,gl,gr,phil,phir,eta)
      integer j, nnd, l
      real*8 x, a, c, xnd(*), pt(*), qt(*), ft(*)
      real*8 gl(*), gr(*), phil(*), phir(*), eta(*)
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c  -----------------------------------------------------------------------
c  This routine defines a local integral equation from ODE
c  int_eqn: with inhomogeneous term, int_eig: without inhomogeneous term
c
c  all arguments except j and nnd are defind between 1+nnd*(j-1) and nnd*j
c  xnd, pt, qt, ft is a mandatory input
c  gl, gr, eta, phil, phir are mandatory outputs
c  -----------------------------------------------------------------------
c    ur(x) * wron = vl(x) = gl1 x + gl0 or gl1 * dcosh(x) - gl0 * dsinh(x)
c    ul(x) * wron = vr(x) = gr1 x + gr0 or gr1 * dcosh(x) - gr0 * dsinh(x)
**************************************************************************
      if ( q0 .eq. 0 ) then
         do l = 1 + nnd*(j-1), nnd*j
            x = xnd(l)
            gl(l) = gl1*x + gl0
            gr(l) = gr1*x + gr0
            phil(l) = ( pt(l) * gr1 + qt(l) * gr(l) ) / wron
            phir(l) = ( pt(l) * gl1 + qt(l) * gl(l) ) / wron
            eta(l) = ft(l) - ( pt(l)*ui1 + qt(l)*(ui1*x+ui0) )
         end do
      else
         do l = 1 + nnd*(j-1), nnd*j
            x = xnd(l)
            gl(l) = gl1 * dcosh(x) - gl0 * dsinh(x)
            gr(l) = gr1 * dcosh(x) - gr0 * dsinh(x)
            phil(l) = ( pt(l) * ( gl1*dsinh(x) - gl0*dcosh(x) )
     &                  + (qt(l)-q0) * gr(l) ) / wron
            phir(l) = ( pt(l) * ( gr1*dsinh(x) - gr0*dcosh(x) )
     &                  + (qt(l)-q0) * gl(l) ) / wron
            eta(l) = ft(l) - ( ui2 + pt(l)*(2d0*ui2*x+ui1)
     &                             + qt(l)*(ui2*x**2+ui1*x+ui0) )
         end do
      end if

      end

**************************************************************************
      subroutine int_eig(j,nnd,a,c,xnd,pt,qt,wcopy,gl,gr,psil,psir)
      integer j, nnd, l
      real*8 x, a, c, xnd(*), pt(*), qt(*)
      real*8 wcopy, gl(*), gr(*), psil(*), psir(*)
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c  -----------------------------------------------------------------------
c  This routine defines a local integral equation from ODE
c  int_eqn: with inhomogeneous term, int_eig: without inhomogeneous term
c
c  all arguments except j and nnd are defind between 1+nnd*(j-1) and nnd*j
c  xnd,pt, qt, ft is a mandatory input
c  wcopy = a copy of wron for shifted eigenfunction iteration.
c  gl, gr, eta, phil, phir are mandatory outputs
c  -----------------------------------------------------------------------
c    ur(x) * wron = vl(x) = gl1 x + gl0 or gl1 * dcosh(x) - gl0 * dsinh(x)
c    ul(x) * wron = vr(x) = gr1 x + gr0 or gr1 * dcosh(x) - gr0 * dsinh(x)
**************************************************************************
      wcopy = wron
      if ( q0 .eq. 0 ) then
         do l = 1 + nnd*(j-1), nnd*j
            x = xnd(l)
            gl(l) = gl1*x + gl0
            gr(l) = gr1*x + gr0
            psil(l) = ( pt(l) * gr1 + qt(l) * gr(l) ) / wron
            psir(l) = ( pt(l) * gl1 + qt(l) * gl(l) ) / wron
         end do
      else
         do l = 1 + nnd*(j-1), nnd*j
            x = xnd(l)
            gl(l) = gl1 * dcosh(x) - gl0 * dsinh(x)
            gr(l) = gr1 * dcosh(x) - gr0 * dsinh(x)
            psil(l) = ( pt(l) * ( gl1*dsinh(x) - gl0*dcosh(x) )
     &                  + (qt(l)-q0) * gr(l) ) / wron
            psil(l) = ( pt(l) * ( gr1*dsinh(x) - gr0*dcosh(x) )
     &                  + (qt(l)-q0) * gl(l) ) / wron
         end do
      end if

      end

**********************************************************************
      subroutine initmesh(meshtype, br, nsub, leafid, a, c, b, blength)
      integer meshtype, nsub, nmxsub
      parameter (nmxsub=1024)
      integer leafid(nmxsub)
      real*8 br, a, c, b(nmxsub), blength(nmxsub)
c -----------------------------------------
c  Mesh type:
c    0 - Equal space
c    1,2 - Polynomial order 1 and 2
c    -1(br) - centered, br = b(center+1)/b(center)
c    -2(br) - 0 divides the intervals (nsub=even) : br = b(center+1)/b(center)
**********************************************************************
      integer kp, j
      real*8 bsum

      do kp = 1, nsub
         j = leafid(kp)
         if ( meshtype .eq. 0 ) then
            b(j) = ( dble(kp-1)*c + dble(nsub-kp+1)* a ) / dble(nsub)
            blength(j) = (c-a) / dble(nsub)
         else if (meshtype .eq. 1 ) then
            b(j) = dble(kp*(kp-1)) / (nsub*(nsub+1))
            blength(j) = 2d0 * kp / (nsub*(nsub+1))
         else if (meshtype .eq. 2 ) then
            b(j) = dble((kp-1)*kp*(2*kp-1)) / (nsub*(nsub+1)*(2*nsub+1))
            blength(j) = 6d0 * kp**2 / (nsub*(nsub+1)*(2*nsub+1))
         else if ( meshtype.eq.-1 .and. nsub/2*2.eq.nsub ) then
            bsum = 2d0 * br * (br**nsub-1d0) / (br**2-1d0)
            blength(j) = (c-a)/bsum * br**iabs(2*kp-nsub-1)
            if ( kp .le. nsub/2 ) then
               b(j) = (c+a)/2d0 - (c-a)/bsum
     &               * (br**(nsub+3-2*kp)-br)/(br**2-1d0)
            else
               b(j) = (c+a)/2d0 + (c-a)/bsum
     &               * (br**(2*kp-nsub-1)-br)/(br**2-1d0)
            end if
         else if ( meshtype.eq.-1 .and. nsub/2*2.ne.nsub ) then
            bsum = 2d0 * (br**(nsub+1)-1d0) / (br**2-1d0) - 1d0
            blength(j) = (c-a)/bsum * br**iabs(2*kp-nsub-1)
            if ( kp .le. (nsub+1)/2 ) then
               b(j) = (c+a)/2d0 - (c-a)/bsum
     &               * ( (br**(nsub+3-2*kp)-1d0)/(br**2-1d0) -.5d0)
            else
               b(j) = (c+a)/2d0 + (c-a)/bsum
     &               * ( (br**(2*kp-nsub-1)-1d0)/(br**2-1d0) -.5d0)
            end if
         else if ( meshtype.eq.-2 .and. nsub/2*2.eq.nsub ) then
            bsum = br * (br**nsub-1d0) / (br**2-1d0)
            if ( kp .le. nsub/2 ) then
               blength(j) = -a/bsum * br**iabs(2*kp-nsub-1)
               b(j) = a/bsum * (br**(nsub+3-2*kp)-br)/(br**2-1d0)
            else
               blength(j) = c/bsum * br**iabs(2*kp-nsub-1)
               b(j) = c/bsum * (br**(2*kp-nsub-1)-br)/(br**2-1d0)
            end if
         else
            print*,'ERROR: Unknown Mesh Initialization' 
         end if
      end do

      end
*************************************************************************
      subroutine evalsigma(nnd,sigma,located,leaflist,leafid,
     &           lambda1,lambda2,phil,phir,eta)
      integer nmxsub, nnd
      parameter (nmxsub=1024)
      integer located, leaflist(nmxsub), leafid(nmxsub)
      real*8 lambda1(2*nmxsub), lambda2(2*nmxsub)
      real*8 sigma(*), phil(*), phir(*), eta(*)
c -------------------------------------------------
c   Step 2.C of the algorithm by Lee and Greengard, SISC 1995
c -------------------------------------------------
c     sigma = (def) eta + lambda1*phi_l + lambda2*phi_r
*************************************************************************
      integer i, kp, k, j

      do kp = 1, located
           k = leaflist(kp)
           j = leafid(kp)
           do i = (j-1)*nnd+1, j*nnd
              sigma(i) = eta(i)+lambda1(k)*phil(i)+lambda2(k)*phir(i)
           end do
      end do
      end

*************************************************************************
      subroutine defint_intv(nsub,leaflist, deffl,deffr,
     &  lambda1,lambda2,alpha11,alpha12,alpha21,alpha22,delta_l,delta_r)
      integer nmxsub, nsub
      parameter (nmxsub=1024)
      integer leaflist(nmxsub)
      real*8 deffl(nmxsub+1), deffr(nmxsub+1)
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      real*8 lambda1(2*nmxsub),lambda2(2*nmxsub)
c -------------------------------------------------
c   Step 3.A of the algorithm by Lee and Greengard, SISC 1995
c -------------------------------------------------
c  deffl(j) = \int_A^B(j-1) { vl sigma }
c  deffr(j) = \int_B(j)^C { vr sigma }
c      where sigma = eta + lambda1*phi_l + lambda2*phi_r
*************************************************************************
      integer k, kp

        deffl(1) = 0d0
        do kp = 1, nsub
           k = leaflist(kp)
           deffl(kp+1) = deffl(kp) +  delta_l(k)
     &            + lambda1(k)*alpha11(k) + lambda2(k)*alpha12(k)
        end do

        deffr(nsub+1) = 0d0
        do kp = nsub, 1, -1
           k = leaflist(kp)
           deffr(kp) = deffr(kp+1) + delta_r(k)
     &            + lambda1(k)*alpha21(k) + lambda2(k)*alpha22(k)
        end do
      end

*************************************************************************
      subroutine int_cheby(ieval,nnd,located,leafid,
     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
      integer ieval, nmxnd, nmxsub, nnd
      parameter (nmxnd=256,nmxsub=1024)
      integer located, leafid(nmxsub)
      real*8 blength(nmxsub), xeval(*)
      real*8 glp_x, grp_x
      real*8 gl(*), gr(*), ucom(*), upcom(*), sigma(*)
      real*8 deffl(nmxsub+1), deffr(nmxsub+1)
      real*8 spdef, spbx, spxb
      common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c------------------------------------------------------------------------
c   Step 3.B of the algorithm by Lee and Greengard, SISC 1995
c------------------------------------------------------------------------
c    Given Sigma (a density funtion) and G0 (through ul,ur,vl,vr)
c        Evaluate phi on all of chebyshev nodal points (nnd*located)
c    Phi(x) = ul(x) \int_a^x{vl sigma} + ur(x) \int_x^b{vr sigma}
c    Sol(x) = Ui(x) + Uh(x) = ui2 * x**2 + ui1 * x + ui0 + Phi(x)
c
c  G(x,t) = gr(x) gl(t) / wron  or  gr(t) gl(x) / wron,
c      where wron(gl,gr) = gl(x) grp(x) - glp(x) gr(x) = constant
c------------------------------------------------------------------------
c   input variables: nmxsub, nnd, located, leafid
c                   sigma, gl, gr, xnd, blength, deffl, deffr
c   (optional) output variables : ucom, upcom
*************************************************************************
      integer i, j, k, kp, ip, lp
      real*8 scale, sum_r, sum_l

      do kp = 1, located
           j = leafid(kp)

c ----------------------------------------------------------------
c  ucom(x) = \int G0(x,t) sigma(t) dt
c         = ul(x) * int_A^x {vl*sigma} + ur(x) int_x^C {vr*sigma}
c  upcom(x) = \int G1(x,t) sigma(t) dt
c         = ulp(x) * int_A^x {vl*sigma} + urp(x) int_x^C {vr*sigma}
c ----------------------------------------------------------------
          scale = 0.5d0 * blength(j)
          do i = 1, nnd
             ip = i + nnd * (j-1)
             sum_l = 0d0
             sum_r = 0d0
             do k = 1, nnd
                lp = k + nnd * (j-1)
                sum_l = sum_l + spbx(i,k)*gl(lp)* sigma(lp)
                sum_r = sum_r + spxb(i,k)*gr(lp)* sigma(lp)
             end do
             if ( ieval/2 - ieval/4*2 .eq. 1 ) then
                upcom(ip) = ( grp_x(xeval(ip))*(deffl(kp)+scale*sum_l)
     &                      + glp_x(xeval(ip))*(scale*sum_r+deffr(kp+1))
     &                      ) / wron  + 2d0*ui2*xeval(ip) + ui1
             end if
             if ( ieval - ieval/2*2 .eq. 1 ) then
                ucom(ip) = ( gr(ip)*(deffl(kp)+scale*sum_l)
     &                     + gl(ip)*(scale*sum_r+deffr(kp+1)) ) / wron
     &                   + ui2*xeval(ip)**2 + ui1*xeval(ip) + ui0
             end if
	  end do
      end do
      end

*************************************************************************
      subroutine int_bnodes(ieval,located,leafid,
     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
      integer ieval, nmxsub
      parameter (nmxsub=1024)
      integer located, leafid(nmxsub)
      real*8 blength(nmxsub), xeval(*)
      real*8 glp_x, grp_x, gl_x, gr_x
      real*8 gl(*), gr(*), ucom(*), upcom(*), sigma(*)
      real*8 deffl(nmxsub+1), deffr(nmxsub+1)
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c------------------------------------------------------------------------
c   Step 3.B of the algorithm by Lee and Greengard, SISC 1995
c------------------------------------------------------------------------
c    Given Sigma (a density funtion) and G0 (through ul,ur,vl,vr)
c        Evaluate phi on each nodal point b(1)=a, ... , b(i), ... , c
c    Phi(x) = ul(x) \int_a^x{vl sigma} + ur(x) \int_x^b{vr sigma}
c    Sol(x) = Ui(x) + Uh(x) = ui2 * x**2 + ui1 * x + ui0 + Phi(x)
c
c  G(x,t) = gr(x) gl(t) / wron  or  gr(t) gl(x) / wron,
c      where wron(gl,gr) = gl(x) grp(x) - glp(x) gr(x) = constant
c------------------------------------------------------------------------
c   input variables: nmxsub, nnd, located, leafid
c                   sigma, gl, gr, xnd, blength, deffl, deffr
c   (optional) output variables : ucom, upcom
*************************************************************************
      integer j, kp
      real*8 x

      do kp = 1, located+1
         if ( kp .ne. located+1) then
            j = leafid(kp)
            x = xeval(j)
         else
            j = leafid(located)
            x = xeval(j) + blength(j)
         end if

c ----------------------------------------------------------------
c  ucom(x) = \int G0(x,t) sigma(t) dt
c         = ul(x) * int_A^x {vl*sigma} + ur(x) int_x^C {vr*sigma}
c  upcom(x) = \int G1(x,t) sigma(t) dt
c         = ulp(x) * int_A^x {vl*sigma} + urp(x) int_x^C {vr*sigma}
c ----------------------------------------------------------------
         if ( ieval/2 - ieval/4*2 .eq. 1 ) then
            upcom(kp) = ( grp_x(x)*deffl(kp) + glp_x(x)*deffr(kp) )
     &                  / wron + 2d0*ui2*x + ui1
         end if
         if ( ieval - ieval/2*2 .eq. 1 ) then
            ucom(kp) = ( gr_x(x)*deffl(kp) + gl_x(x)*deffr(kp) )
     &                 / wron + ui2*x**2 + ui1*x + ui0
         end if
      end do
      end
*************************************************************************
      subroutine int_bnodes1(ieval,located,leafid,
     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
c
c  DOES END POINTS ONLY!
c
      integer ieval, nmxsub
      parameter (nmxsub=1024)
      integer located, leafid(nmxsub)
      real*8 blength(nmxsub), xeval(*)
      real*8 glp_x, grp_x, gl_x, gr_x
      real*8 gl(*), gr(*), ucom(*), upcom(*), sigma(*)
      real*8 deffl(nmxsub+1), deffr(nmxsub+1)
      real*8 q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
      common /bc/ q0, wron, gl1, gr1, gl0, gr0, ui2, ui1, ui0
c------------------------------------------------------------------------
c   Step 3.B of the algorithm by Lee and Greengard, SISC 1995
c------------------------------------------------------------------------
c    Given Sigma (a density funtion) and G0 (through ul,ur,vl,vr)
c        Evaluate phi on each nodal point b(1)=a, ... , b(i), ... , c
c    Phi(x) = ul(x) \int_a^x{vl sigma} + ur(x) \int_x^b{vr sigma}
c    Sol(x) = Ui(x) + Uh(x) = ui2 * x**2 + ui1 * x + ui0 + Phi(x)
c
c  G(x,t) = gr(x) gl(t) / wron  or  gr(t) gl(x) / wron,
c      where wron(gl,gr) = gl(x) grp(x) - glp(x) gr(x) = constant
c------------------------------------------------------------------------
c   input variables: nmxsub, nnd, located, leafid
c                   sigma, gl, gr, xnd, blength, deffl, deffr
c   (optional) output variables : ucom, upcom
*************************************************************************
      integer j, kp
      real*8 x

      do kp = 1, located+1,located
         if ( kp .ne. located+1) then
            j = leafid(kp)
            x = xeval(j)
         else
ccc            kp = located+1
            j = leafid(located)
            x = xeval(j) + blength(j)
         end if

c ----------------------------------------------------------------
c  ucom(x) = \int G0(x,t) sigma(t) dt
c         = ul(x) * int_A^x {vl*sigma} + ur(x) int_x^C {vr*sigma}
c  upcom(x) = \int G1(x,t) sigma(t) dt
c         = ulp(x) * int_A^x {vl*sigma} + urp(x) int_x^C {vr*sigma}
c ----------------------------------------------------------------
         if ( ieval/2 - ieval/4*2 .eq. 1 ) then
            upcom(kp) = ( grp_x(x)*deffl(kp) + glp_x(x)*deffr(kp) )
     &                  / wron + 2d0*ui2*x + ui1
         end if
         if ( ieval - ieval/2*2 .eq. 1 ) then
            ucom(kp) = ( gr_x(x)*deffl(kp) + gl_x(x)*deffr(kp) )
     &                 / wron + ui2*x**2 + ui1*x + ui0
         end if
      end do
      end

********************************
* dgeco, dgeco -> dgesl
* all come from linpack
* except dasum.f from blas
********************************
      subroutine dgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(1)
      double precision a(lda,1),z(1)
      double precision rcond
c
c     dgeco factors a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgefa is slightly faster.
c     to solve  a*x = b , follow dgeco by dgesl.
c     to compute  inverse(a)*c , follow dgeco by dgesl.
c     to compute  determinant(a) , follow dgeco by dgedi.
c     to compute  inverse(a) , follow dgeco by dgedi.
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgefa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,dx(1)
      integer i,incx,m,mp1,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end


      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end


      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision dx(1),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end
**********************************************************************
      function monfn(nnd, j, sigma)
      integer nnd, j, nmxnd
      parameter (nmxnd=256)
      real*8 monfn, sigma(*), ch(nmxnd)
c------------------------------------------------------------------------
c   A monitor function (phi) is defined as
c   phi(x) ^ (nnd-2)
c   = (d/dx)^(nnd-1) \int^x sigma(t) dt
c   = (d/dx)&(nnd-1) { ch_(nnd-1) T_(nnd-1) + ... }
c             where ch_(nnd-1) = sigma_(nnd-2) / (2*nnd-2)
c   = ch_(nnd-1) * (nnd-1)! * 2^(nnd-1)
c   \propto sigma_(nnd-2)
c ------------------------------------------------------------
c  TO DO :
c  1. find better monitor function
**********************************************************************
        call chftrans(ch, sigma(1+nnd*(j-1)), nnd)
        monfn = dabs(ch(nnd-1)) + dabs(ch(nnd)-ch(nnd-2))
*       print*,j, dabs(ch(nnd-1)),dabs(ch(nnd)-ch(nnd-2))
      end


**********************************************************************
      subroutine doublemesh(b,blength,
     &       divided, updated,updateid,updatelist,updatedone,
     &       located, leafid,leaflist)
      integer nmxsub
      parameter (nmxsub=1024)
      real*8 b(nmxsub),blength(nmxsub)
      integer divided,updated,located, leafid(nmxsub),leaflist(nmxsub)
      integer updateid(nmxsub), updatelist(nmxsub), updatedone(nmxsub)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
c------------------------------------------------------------------------
c Divide all intervals halves
c ---------------------------
c   in/out : b, blength
c   input :  located, leafid, leaflist
c   output : divided, updated,updateid,updatelist,updatedone
c   common /hierarchy/ updated though make/free children
**********************************************************************
      integer j, k, kp

c ----------------------
c  Divide
c ----------------------
      updated = 0
      divided = 0
      if ( 2*located .gt. nmxsub ) then
         print*, 'EXHAUST all subintervals, can''t double mesh'
         return
      end if
      do kp = 1, located
         k = leaflist(kp)
         j = leafid(kp)
             call makechildren(k)
             updateid(updated+1) = j
             updateid(updated+2) = tree(3,tree(3,k))
             updatelist(updated+1) = tree(2,k)
             updatelist(updated+2) = tree(3,k)
             updatedone(updated+1) = 0
             updatedone(updated+2) = 0
             blength(tree(3,tree(3,k))) = blength(j) / 2d0
             b(tree(3,tree(3,k))) = b(j) + blength(j) / 2d0
             blength(j) = blength(j) / 2d0
             updated = updated + 2
             divided = divided + 1
      end do
             
      end

**********************************************************************
      subroutine topmesh(tolratio, sumtrun, nnd, b,blength,sigma,
     &       divided, updated,updateid,updatelist,updatedone,
     &       located, leafid,leaflist)
      real*8 sigma(*)
      integer nmxsub, nnd
      parameter (nmxsub=1024)
      real*8 tolratio,maxtol, mintol, maxtrun, sumtrun
      real*8 b(nmxsub), blength(nmxsub), monitor(nmxsub), monfn
      integer located, leafid(nmxsub), leaflist(nmxsub)
      integer divided, updated
      integer updateid(nmxsub), updatelist(nmxsub), updatedone(nmxsub)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
c------------------------------------------------------------------------
c Divide if Mi >= max(Mmax*tolratio)
c Merge if Mi+Mj <= min(M_div) / 2**nnd
c ---------------------------
c   input : nnd, tolratio, sigma
c   output : sumtrun
c   in/out : b, blength
c   input : nnd, located, leafid, leaflist
c   output : divided, updated,updateid,updatelist,updatedone
c   common /hierarchy/ updated though make/free children
c ---------------------------
c  TO DO :
c  1. Can you make an automatic mesh refinement routine without
c     user specified constant TOLRATIO ?
**********************************************************************
      integer j, k, kp, prnted

      prnted = 0

c ----------------------
c  Scan Monitor function
c ----------------------
      maxtrun = 0d0
      sumtrun = 0d0
      do kp = 1, located
         j = leafid(kp)
         monitor(j) = monfn(nnd, j, sigma)
         maxtrun = dmax1(maxtrun, monitor(j))
         sumtrun = sumtrun + blength(j)**(nnd-1)*monitor(j)
*        print'(f12.5,g10.3,i6,f12.7)',b(j),blength(j),j,monitor(j)
      end do

      updated = 0
      divided = 0
      mintol = 1d38
      maxtol = maxtrun * tolratio

c ----------------------
c  Divide
c ----------------------
      do kp = 1, located
         k = leaflist(kp)
         j = leafid(kp)
         if ( monitor(j).ge.maxtol ) then
            mintol = dmin1(mintol, monitor(j)/2d0**nnd)
            if ( prnted.eq.0 .and. leaf(0).gt.nmxsub ) then
               print*,'EXHAUST all subintervals, can''t divide it'
               prnted = 1
            else
               call makechildren(k)
               updateid(updated+1) = j
               updateid(updated+2) = tree(3,tree(3,k))
               updatelist(updated+1) = tree(2,k)
               updatelist(updated+2) = tree(3,k)
               updatedone(updated+1) = 0
               updatedone(updated+2) = 0
               blength(tree(3,tree(3,k))) = blength(j) / 2d0
               b(tree(3,tree(3,k))) = b(j) + blength(j) / 2d0
               blength(j) = blength(j) / 2d0
               updated = updated + 2
               divided = divided + 1
            end if
         end if
      end do

c ----------------------
c  Merge
c ----------------------
      do kp = 1, located
         k = leaflist(kp)
         j = leafid(kp)
         if ( tree(3,tree(1,k)).eq.k .and.
     &      tree(2,tree(1,k)).eq.leaflist(kp-1) .and.
     &      monitor(j)+monitor(leafid(kp-1)) .lt. mintol ) then
            call freechildren(tree(1,k))
            updateid(updated+1) = leafid(kp-1)
            updatelist(updated+1) = tree(1,k)
            updatedone(updated+1) = j
            blength(leafid(kp-1)) = blength(leafid(kp-1)) + blength(j)
            updated = updated + 1
         end if
      end do
             
      end

**********************************************************************
      subroutine writefn(nsub,nnd,ieval,xnd,u,ugvn,up,sigma,leafid)
      integer nsub, nnd, ieval, nmxsub
      parameter (nmxsub=1024)
      integer leafid(nmxsub)
      real*8 xnd(*), u(*), ugvn(*), up(*), sigma(*)
**********************************************************************
      integer kp, j, l
         open(file='exout.plt',unit=23,status='unknown')

c        write (23,'(a,i4,a9,i3)') '# nsub =',nsub,'nnd =',nnd
c        write (23,'(a,4a15)') '# x','Sol','Real','Err','Sigma'
         do kp = 1, nsub
            j = leafid(kp)
            do l = 1 + nnd*(j-1), nnd*j
               if ( ieval/2 - ieval/4*2 .eq. 1 ) then
                  write (23,'(f10.5,f15.9,f15.9,e15.6,e15.7)')
     &              xnd(l),u(l), ugvn(l), up(l), sigma(l)
               else
                  write (23,'(f10.5,f15.9,f15.9,e15.6,e15.7)')
     &              xnd(l),u(l), ugvn(l), u(l)-ugvn(l), sigma(l)
               end if
            end do
         end do

         close(unit=23, status='keep')
      end

**********************************************************************
      subroutine writemesh(countup, located,leafid,b,c)
      integer nmxsub
      parameter (nmxsub=1024)
      integer countup, located, leafid(nmxsub), inited, i
      real*8 b(nmxsub), c
      data inited /0/
**********************************************************************
      if ( inited .eq. 0 ) then
         open(file='exout.tree',unit=25,status='unknown')
         inited = 1
      end if

      do i = 1, located
         write (25,'(f8.4,i3)')  b(leafid(i)), countup
      end do
      write (25,'(f8.4,i3)')  c, countup

      end

**********************************************************************
      subroutine writestatus(countup, l2err, errest, trun,
     &                divided, updated, located, discond, mergcond)
      integer countup, divided, updated, located, inited
      real*8 l2err, errest, trun, discond, mergcond
      data inited /0/
**********************************************************************
      if ( inited .eq. 0 ) then
         open(file='exout.conv',unit=27,status='unknown')
c        write(27,'(a2,2a4,a4,3a11,2a15)'),
c    &        '#C','Div','Del','Loc',
c    &        'P_i','Merge','Trun', 'Conv Err','E: Real,Conv'
         inited = 1
      end if

      write(27,'(i2,2i4,i4,3e11.4,2e15.8)')
     &        countup,divided,max0(updated-2*divided,0),located,
     &        discond, mergcond, trun, errest, l2err
      end

**********************************************************************
      subroutine prerr2(nsub,nnd,blength,u,ugvn,leafid)
      integer nsub,nnd, nmxnd, nmxsub
      parameter (nmxnd=256,nmxsub=1024)
      integer leafid(nmxsub)
      real*8 dsqrt, blength(nmxsub), u(*), ugvn(*)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
c -------------------------------------------------
c  print Rel_2_Error = |u-ugvn|_2 / |ugvn|_2
**********************************************************************
      integer kp, j, l, lp
      real*8 err2, l2norm
         err2 = 0d0
         l2norm = 0d0
         do kp = 1, nsub
            j = leafid(kp)
            do lp = 1, nnd
                l = lp + nnd*(j-1)
                err2 = err2 + (ugvn(l)-u(l))**2 * blength(j) * ns(lp)
                l2norm = l2norm + ugvn(l)**2 * blength(j)
            end do
         end do
         print'(a,i4,a5,i2,a5,i5,a,e15.8)',
     &       'm=',nsub,'p=',nnd,'n=',nsub*nnd,
     &       ',   Relative L2 Error =',dsqrt(err2/l2norm)
      end

**********************************************************************
      function err2(nsub,nnd,blength,u,ugvn,leafid)
      integer nmxnd, nmxsub, nsub, nnd
      parameter (nmxnd=256,nmxsub=1024)
      integer leafid(nmxsub)
      real*8 err2, dsqrt, blength(nmxsub), u(*), ugvn(*)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
c -------------------------------------------------
c  print Rel_2_Error = |u-ugvn|_2 / |ugvn|_2
**********************************************************************
      integer kp, j, l, lp
      real*8 err, l2norm
         err = 0d0
         l2norm = 0d0
         do kp = 1, nsub
            j = leafid(kp)
            do lp = 1, nnd
                l = lp + nnd*(j-1)
                err = err + (ugvn(l)-u(l))**2 * blength(j) * ns(lp)
                l2norm = l2norm + ugvn(l)**2 * blength(j)
            end do
         end do
         err2 = dsqrt(err/l2norm)
      end

************************************************************************
      subroutine two2one(nnd,combined,first,second)
      integer nmxnd, nnd
      real*8 combined(nnd), first(nnd), second(nnd)
      parameter (nmxnd=256)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
c  -------------------------------------------------------------------
c   when two functions are defined on adjacent intervals with same length,
c   evaluate functional values on chevshev nodes of merged interval 
************************************************************************
      real*8 work(nmxnd), ti, sum
      integer i,k

c ----------------------
c f(x<0) = fl(2x+1)
c ----------------------
      call chftrans(work, first, nnd)
      do i = 1, nnd/2
         ti = dacos(2d0*c(i)+1d0)
         sum = 0d0
         do k = 1, nnd
            sum = sum + dcos( (k-1) * ti ) * work(k)
         end do
         combined(i) = sum
      end do
          
c -------------------------------------------------------
c f(x_mid = boundary between two intevals) = (fl(1)+fr(-1)) / 2d0
c f(x=0) = (fl(1)+fr(-1)) / 2d0 when nnd = odd
c And chebychev transfromation for second interval
c -------------------------------------------------------
      if ( (nnd+1)/2*2 .eq. nnd+1 ) then
         sum = 0d0
         do k = 1, nnd
            sum = sum + work(k)
         end do
      end if

      call chftrans(work, second, nnd)

      if ( (nnd+1)/2*2 .eq. nnd+1 ) then
         do k = 1, nnd*3/4
            sum = sum + (-1)**(k-1) * work(k)
         end do
         combined((nnd+1)/2) = sum / 2d0
      end if

c ----------------------
c f(x>0) = fr(2x-1)
c ----------------------
      do i = (nnd+1)/2+1, nnd
         ti = dacos(2d0*c(i)-1d0)
         sum = 0d0
         do k = 1, nnd*3/4
            sum = sum + dcos( (k-1) * ti ) * work(k)
         end do
         combined(i) = sum
      end do

      end
         
**********************************************************************
      function updateerr(located,leafid,updated,updateid,updatedone,
     &        blength,nnd,eval,seval)
      integer nmxnd, nmxsub, nnd
      parameter (nmxnd=256,nmxsub=1024)
      integer located, leafid(nmxsub)
      integer updated, updateid(nmxsub), updatedone(nmxsub)
      real*8 updateerr, eval(*), seval(*)
      real*8 combined(nmxnd), blength(nmxsub)
c -------------------------------------------------
c  updateerr = |eval-seval|_2 / |(eval+seval)/2|_2
c  Unfortuately, if there were updates in mesh then
c  eval(l) and seval(l) represent different nodal points
c  So use two2one subroutine for updated intervals
**********************************************************************
      integer lp, up, j, jl, jr
      real*8 bj, num, den

      lp = 1
      up = 1
      num = 0d0
      den = 0d0

      do while (lp .le. located)
         j = leafid(lp)
         jl = (j-1)*nnd + 1
c wasn't updated
         if ( up.gt.updated .or.
     &        j.ne.updateid(up) ) then
            bj = blength(j)
            call l2sum(num, den, seval(jl), eval(jl), nnd, bj)
            lp = lp + 1
c divided
         else if ( leafid(lp).eq.updateid(up) .and.
     &             updatedone(up).le.0 ) then
            bj = blength(j) * 2d0
            jr = (leafid(lp+1)-1)*nnd +1
            call two2one(nnd,combined, eval(jl),eval(jr))
            call l2sum(num, den, seval(jl), combined, nnd, bj)
            up = up + 2
            lp = lp + 2
         else
c merged
            bj = blength(j)
            jr = (updatedone(up)-1)*nnd +1
            call two2one(nnd,combined, seval(jl), seval(jr))
            call l2sum(num, den, combined, eval(jl), nnd, bj)
            up = up + 1
            lp = lp + 1
         end if
      end do

      updateerr = dsqrt(num/den)
      end

**********************************************************************
      subroutine l2sum(num, den, s1, s2, nnd, blength)
      integer nnd, l, nmxnd
      parameter (nmxnd=256)
      real*8 theta(nmxnd), ns(nmxnd), c(nmxnd), chnd01(nmxnd)
      common /geofn/ theta, ns, c, chnd01
      real*8 s1(nnd), s2(nnd), num, den, blength
c -------------------------------------
c  Called by updateerr
**********************************************************************
      do l = 1, nnd
         num = num + (s1(l)-s2(l))**2 * blength * ns(l)
         den = den + (s1(l)+s2(l))**2 / 4d0 * blength 
      end do
      end

**********************************************************************
      subroutine savenodeval(located, leafid, nnd, val, sval)
      integer located, nmxsub, nnd
      parameter (nmxsub=1024)
      integer leafid(nmxsub)
      real*8 val(*), sval(*)
**********************************************************************
      integer kp, j, l

      do kp = 1, located
         j = leafid(kp)
         do l = 1 + nnd*(j-1), nnd*j
            sval(l) = val(l)
         end do
      end do
      end

*****************************************************************
      subroutine upward_sweep(nsub,parentlist,mergcond,
     1        alpha11,alpha12,alpha21,alpha22,delta_l,delta_r)
      integer nmxsub, nsub
      parameter (nmxsub=1024)
      integer parentlist(nmxsub)
      real*8 mergcond
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
      common /hierarchy/ tree, node, leaf
c  -------------------------------------------------------------
c     Merge Condition Number = Max (1-A21b*A12a) in absolute value 
c  -------------------------------------------------------------
c    Step 2.A in the paper by Lee and Greengard, SISC 1995
c       Given alpha, delta in the finest level (nlevel)
c       recusively calculate alpha and delta in all levels
*****************************************************************
      integer k, pc, pa, pb
      real*8 deter

      mergcond = 0d0
c    ------------------------------------------------
c              pc = parentlist(k)
c              alpha_Coarse = alpha(pc)
c              alpha_finer_A = alpha(tree(2),pc)
c              alpha_finer_B = alpha(tree(3),pc)
c              See Eq. 67 & 61
c    ------------------------------------------------

      do k = nsub-1, 1, -1
         pc = parentlist(k)
         pa = tree(2,pc)
         pb = tree(3,pc)

         deter = 1d0 - alpha21(pb) * alpha12(pa)
         mergcond = dmax1(mergcond, dabs(deter))

         alpha11(pc) = alpha11(pb) + (1d0-alpha11(pb)) *
     &         (alpha11(pa)-alpha12(pa)*alpha21(pb)) / deter
         alpha21(pc) = alpha21(pa) +
     &         alpha21(pb)*(1d0-alpha22(pa))*(1d0-alpha11(pa))/deter
         alpha12(pc) = alpha12(pb) +
     &         alpha12(pa)*(1d0-alpha22(pb))*(1d0-alpha11(pb))/deter
         alpha22(pc) = alpha22(pa) + (1d0-alpha22(pa)) *
     &         (alpha22(pb)-alpha12(pa)*alpha21(pb)) / deter

         delta_l(pc) = (1d0-alpha11(pb)) / deter * delta_l(pa) +
     &         delta_l(pb) +
     &         (alpha11(pb)-1d0)*alpha12(pa)/deter * delta_r(pb)
         delta_r(pc) = (1d0-alpha22(pa)) / deter * delta_r(pb) +
     &         delta_r(pa) +
     &         (alpha22(pa)-1d0)*alpha21(pb)/deter * delta_l(pa)

      end do
      end

*****************************************************************
	subroutine down_sweep(nsub,parentlist,alpha11,alpha12,
     1     alpha21,alpha22, delta_l,delta_r, lambda1,lambda2)
      integer nmxsub, nsub
      parameter (nmxsub=1024)
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      real*8 lambda1(2*nmxsub),lambda2(2*nmxsub)
      integer parentlist(nmxsub)
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
      common /hierarchy/ tree, node, leaf
c  -------------------------------------------------------------
c    Step 2.B in the paper by Lee and Greengard, SISC 1995
c       Given alpha, delta in all levels
c       recusively calculate lambda in all levels
*****************************************************************
	integer k, pc, pa, pb
	real*8 deter

c    ---------------------
c      initial data
c      lambda^{0th-level}_{1,2} = 0
c    ---------------------
        lambda1(1+1) = 0d0
        lambda2(1+1) = 0d0

c    ------------------------------------------------
c              pc = parentlist(k)
c              lambda_Coarse = lambda(pc)
c              lambda_finer_A = lambda(tree(2),pc)
c              lambda_finer_B = lambda(tree(3),pc)
c              See Eq. 73 & 75
c    ------------------------------------------------
        do k = 1, nsub-1
           pc = parentlist(k)
           pa = tree(2,pc)
           pb = tree(3,pc)
           deter = 1d0 - alpha21(pb) * alpha12(pa)

           lambda1(pa) = lambda1(pc)
           lambda2(pa) = -(1d0-alpha11(pa))*alpha21(pb)/deter *
     &         lambda1(pc) + (1d0-alpha22(pb))/deter*lambda2(pc) +
     &         alpha21(pb)/deter*(delta_l(pa)-alpha12(pa)*delta_r(pb))
     &             - delta_r(pb)
           lambda1(pb) = +(1d0-alpha11(pa))/deter*lambda1(pc) -
     &         (1d0-alpha22(pb))*alpha12(pa)/deter*lambda2(pc) +
     &         alpha12(pa)/deter*(delta_r(pb)-alpha21(pb)*delta_l(pa))
     &             - delta_l(pa)
           lambda2(pb) = lambda2(pc)
        end do
      end

************************************************************************
      subroutine solve_sigma(nnd, updated, updatelist, updateid,
     &                    located, leaflist, leafid, parentlist,
     &                    blength,gl,gr,phil,phir,eta,
     &                    mergcond, discond, sigma)
      integer nmxsub, nnd, updated, located
      parameter (nmxsub=1024)
      integer updatelist(nmxsub), updateid(nmxsub)
      integer leafid(nmxsub), leaflist(nmxsub), parentlist(nmxsub)
      real*8 blength(nmxsub), discond, mergcond
      real*8 gl(*), gr(*), phil(*), phir(*), eta(*), sigma(*)
c ----------------------------------------------------------------
c
c    Solving an integeral equaion with unknown SIGMA
c    To get solution U(*), use solve_eval
c
c input:
c    nnd : order of solver
c    updated, updatelist, updateid for local integral equation solver
c    located, leaflist, leafid : list of subintervals
c        may be linked with the variables for updated subintervals
c    parentlist for hierarchical solver of global integral equation
c    blength(j), j = leafid(kp), kp = 1..located
c    gl, gr, phil, phir, eta at xnd(l+nnd*(j-1)), l = 1..nnd
c
c output :
c    discond, mergcond : Pivoting Not implement
c    sigma : the solution of the problem at xnd(l+nnd*(j-1)), l = 1..nnd
c
c common variables :
c    Used in discrete.f & initialized in chsetup(nnd) in chtools.f
c        common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
c    Used in recur.f & initialized in make*, free*, gentree in tree.f
c        common /hierarchy/ tree, node, leaf
c ----------------------------------------------------------------
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      real*8 lambda1(2*nmxsub),lambda2(2*nmxsub)
      common /inner/ alpha11, alpha12, alpha21, alpha22,
     &                delta_l, delta_r, lambda1, lambda2
************************************************************************

      integer kp, k, j
      real*8 localcond(nmxsub)

c ---------------------------------------------
c  Solve local integral equations for eta, phi
c ---------------------------------------------
         do kp = 1, updated
            k = updatelist(kp)
            j = updateid(kp)
            call discret(nnd,j,blength, gl,gr,localcond,eta,phil,phir)
            call mkcoef(nnd,k,j,blength,gl,gr,eta,phil,phir,
     &              alpha11,alpha21,alpha12,alpha22,delta_l,delta_r)
         end do

c --------------------------------
c  Condtion numer of local solver
c --------------------------------
         discond = 1d0
         do kp = 1, located
            discond = dmin1( discond, localcond(leafid(kp)) )
         end do

c -------------------------------
c  Solve global integral equation
c -------------------------------
         call upward_sweep(located,parentlist,mergcond,
     &           alpha11,alpha12,alpha21,alpha22,delta_l,delta_r)
         call down_sweep(located,parentlist,alpha11,alpha12,
     &           alpha21,alpha22, delta_l,delta_r, lambda1,lambda2)
         call evalsigma(nnd,sigma,located,leaflist,leafid,
     &           lambda1,lambda2,phil,phir,eta)
      end

************************************************************************
      subroutine solve_eval(ieval, nnd, located, leaflist, leafid,
     &                    blength,xeval,gl,gr, sigma, upcom, ucom)
      integer nmxsub
      parameter (nmxsub=1024)
      integer ieval, nnd, located, leafid(nmxsub), leaflist(nmxsub)
      real*8 blength(nmxsub), xeval(*)
      real*8 gl(*), gr(*), ucom(*), upcom(*), sigma(*)
c ----------------------------------------------------------------
c
c    Using Sigma and G0,G1 (through gl,gr),
c    Evaluate U(*) / U'(*) at specified point(s)
c 
c input:
c -> ieval : solution to evaluate (applying Green's function with sigma)
c        (xeval:bnodes=1,chnodes=0)(upcom:yes=1,no=0)(ucom:yes=1,no=0)_2
c        Examples) 3 = both u and u' at chebyshev nodes
c        Examples) 4 + XY_2 = u/u' at b nodes
c    located, leaflist, leafid : list of subintervals
c        may be linked with the variables for updated subintervals
c    blength(j), j = leafid(kp), kp = 1..located
c -> xeval : nodal points where solution will be evaluated
c        ieval=0 then won't be referenced
c        btest(ieval,3) = 0 then xnd(l+nnd*(j-1)) = chebyshev nodes
c        btest(ieval,3) = 1 then b(j), j = leafid(1..located) and c
c        cf. btest(i,j) = mod(i/2**(j-1),2)
c    gl, gr at xnd(l+nnd*(j-1)), l = 1..nnd
c    sigma : the solution of the problem at xnd(l+nnd*(j-1)), l = 1..nnd
c
c output :
c -> upcom, ucom : solution at xeval()
c        will be referenced only when requested by ieval
c
c common variables :
c    Used in int_density.f & initialized in chsetup(nnd) in chtools.f
c        common /spint/ spdef(nmxnd), spbx(nmxnd,nmxnd), spxb(nmxnd,nmxnd)
c    Used in int_density.f & initialized in set_examples(a,b) in examples.f
c        common /bc/ q0, wron, zl1, zr1, zl0, zr0, ui2, ui1, ui0
c ----------------------------------------------------------------
      real*8 alpha11(2*nmxsub),alpha12(2*nmxsub)
      real*8 alpha21(2*nmxsub),alpha22(2*nmxsub)
      real*8 delta_l(2*nmxsub),delta_r(2*nmxsub)
      real*8 lambda1(2*nmxsub),lambda2(2*nmxsub)
      common /inner/ alpha11, alpha12, alpha21, alpha22,
     &                delta_l, delta_r, lambda1, lambda2
************************************************************************
      real*8 deffl(nmxsub+1), deffr(nmxsub+1)

c -------------------------------
c   Get the solution of ode
c -------------------------------
         if ( ieval .ne. 0 ) then
            call defint_intv(located, leaflist,
     &           deffl,deffr,lambda1,lambda2,
     &           alpha11,alpha12,alpha21,alpha22,delta_l,delta_r)
         end if

         if ( ieval/4 - ieval/8*2 .eq. 0 ) then
            call int_cheby(ieval,nnd,located,leafid,
     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
         else if ( ieval/4 - ieval/8*2 .eq. 1 ) then
ccc            call int_bnodes(ieval,located,leafid,
ccc     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
c
c  EVALUATE AT END POINTS ONLY!!
c
            call int_bnodes1(ieval,located,leafid,
     &           blength,xeval,gl,gr,sigma,deffl,deffr,ucom,upcom)
         endif

	end

**********************************************************************
      block data nodetree
      integer nmxsub
      parameter (nmxsub=1024)
c ---------------------------------------------------------
c  node(1) : next look up pointer, number of used nodes+2 (initial=2)
c  node(2) : the whole interval
c  ordinary node : 2 <= j <= 2*nmxsub
c      node(j) = 0   : never used, allocate j-th node
c      non-zero(ge2) : next available node number
c ---------------------------------------------------------
c  grand-parent: node number = 2
c     tree(1,1) = 0, tree(2,1) = 0, tree(3,1) = 0
c     tree(1,2) = no parent(1), tree(2,2) = 3, tree(3,2) = 4
c  a node with children: node number = j
c     tree(1,j) = parent of j-th node
c     tree(2,j) = left child of j-th node
c     tree(3,j) = right child of j-th node
c  a leaf, node without child: node number = j
c     tree(1,j) = parent of j-th node
c     tree(2,j) = 0
c     tree(3,j) = leaf number
c ---------------------------------------------------------
c  leaf(0) : next look up pointer, number of used leaves
c  leaf( leaf(0) ) : current avaible leaf if not zero
**********************************************************************
      integer tree(3,2*nmxsub),node(2*nmxsub), leaf(0:nmxsub)
      common /hierarchy/ tree, node, leaf
      data node(1),node(2),node(3) /3,0,0/
      data leaf(0),leaf(1),leaf(2) /2,0,0/
      data tree(1,2),tree(2,2),tree(3,2) /1,0,1/
      end

**********************************************************************
      subroutine gentree(n, updateid, updatelist, parentlist)
      integer nmxsub, n, j
      parameter (nmxsub=1024)
      integer updateid(nmxsub), updatelist(nmxsub), parentlist(nmxsub)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
c ----------------------------------------------
c   Example:n=1    Example:n=4     Example:n=7
c   2              2 - 3 - 5       2 - 3 - 5 - 9
c                    |     6         |   |     10
c                    - 4 - 7         |   - 6 - 11
c                          8         |         12
c                                    - 4 - 7 - 13
c                                        |     14
c                                        - 8
**********************************************************************
        do j = 2, n
           tree(1,j) = (j+1) / 2
           tree(2,j) = 2*j-1
           tree(3,j) = 2*j
        end do

        do j = n+1, 2*n
           tree(1,j) = (j+1) / 2 
           tree(2,j) = 0
           tree(3,j) = j - n
        end do

        node(1) = 2*n + 1
        leaf(0) = n + 1
        call listtree(2, n, updateid, updatelist, parentlist)
      end

**********************************************************************
      subroutine makechildren(j)
      integer nmxsub, j
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub),node(2*nmxsub), leaf(0:nmxsub)
c ----------------------------------------------
c  A LEAF becomes a parent from geting two children
c  one leaf comes from leaf stack and the other from parent
**********************************************************************
      integer leftchild, rightchild, newleaf

        leftchild = node(node(1))
           if ( leftchild .le. 0 ) leftchild = node(1) 
        rightchild = node(node(1)+1)
           if ( rightchild .le. 0 ) rightchild = node(1) + 1
        newleaf = leaf(leaf(0))
           if ( newleaf .le. 0 ) newleaf = leaf(0)

        node(1) = node(1) + 2
        leaf(0) = leaf(0) + 1

        tree(2,j) = leftchild
        tree(1,leftchild) = j
        tree(2,leftchild) = 0
        tree(3,leftchild) = tree(3,j)

        tree(3,j) = rightchild
        tree(1,rightchild) = j
        tree(2,rightchild) = 0
        tree(3,rightchild) = newleaf
      end

**********************************************************************
      subroutine makeoffspring(allocate, n)
      integer nmxsub, n, nextnode, allocate
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub),node(2*nmxsub), leaf(0:nmxsub)
c ----------------------------------------------
c  if there is a node with children kill it
c  use depth first search to find a type 3 node
c  deletion order is inverse order of search
**********************************************************************
      integer j, jold, next, kind, allocated

        allocated = 0

        do while ( allocated .lt. allocate)
           jold = tree(1,n)
           j = n
           do while ( allocated .lt. allocate .and.
     &                j.ne.tree(1,n) .and. j.ne.tree(3,tree(1,n)) )
              next = nextnode(kind, j, jold)
              if ( abs(kind) .eq. 4 ) then
                 call makechildren(j)
                 allocated = allocated + 1
              end if
              jold = j
              j = next
           end do
        end do
      end

**********************************************************************
      subroutine freechildren(j)
      integer nmxsub, j
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub),node(2*nmxsub), leaf(0:nmxsub)
c ----------------------------------------------
c  parent release two children which SHOULD BE leaves
c  push 2 children nodes goes back to the node stack
c  left child's leaf goes back to the leaf stack
**********************************************************************
        leaf(0) = leaf(0) - 1
        leaf(leaf(0)) = tree(3,tree(3,j))

        node( node(1)-1 ) = tree(3,j)
        node( node(1)-2 ) = tree(2,j)
        node(1) = node(1) - 2

        tree(3,j) = tree(3,tree(2,j))
        tree(2,j) = 0
      end

**********************************************************************
      subroutine freeoffspring(released, n)
      integer nmxsub, n, nextnode, released
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub),node(2*nmxsub), leaf(0:nmxsub)
c ----------------------------------------------
c  if there is a node with children kill it
c  use depth first search to find a type 3 node
c  deletion order is inverse order of search
**********************************************************************
      integer j, jold, next, kind

        released = 0
        jold = tree(1,n)
        j = n

        do while ( j.ne.tree(1,n) .and. j.ne.tree(3,tree(1,n)) )
           next = nextnode(kind, j, jold)
           if ( abs(kind) .eq. 3 ) then
              call freechildren(j)
              released = released + 1
           end if
           jold = j
           j = next
        end do
      end

**********************************************************************
      function nextnode(kind,j,jold)
      integer nmxsub, nextnode, kind, j, jold
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
c -----------------------------------------------------------------
c  finding next node in "depth first search algorithm"
c  brife description of the algorithm:
c     1. parent with children            (downward pass)
c     3. parent visited both children    (upward pass)
c     4. a leaf                          (finest level pass)
c.   -3. If you returned root node, Stop
c.   -4. If you have only root node (a leaf), Stop
c -----------------------------------------------------------------
c   Example:n=1
c   2               jold  : 1*
c                   j     : 2
c                   next  : 1*
c                   kind  : -4
c   Example:n=4
c   2 - 3 - 5       jold  : 1* 2  3  5  6  3  4  7  8  4
c     |     6       j     : 2  3  5  6  3  4  7  8  4  2
c     - 4 - 7       next  : 3  5  6  3  4  7  8  4  2  1*
c           8       kind  : 1  1  4  4  3  1  4  4  3  -3
**********************************************************************
c ----- A leaf so return to parent, TYPE 4 or -4
      if ( tree(2,j).eq.0 ) then
         if ( j .eq. 2 ) then
            kind = -4
            nextnode = 1
         else if ( j .eq. tree(2,tree(1,j)) ) then
            kind = 4
            nextnode = tree(3,tree(1,j))
         else
            kind = 4
            nextnode = tree(1,j)
         end if
c ----- Move to parent from right child, TYPE 3 or -3
      else if ( jold.eq.tree(3,j) ) then
         if ( j.eq.2 ) then
            kind = -3
            nextnode = 1
         else if ( j .eq. tree(2,tree(1,j)) ) then
            kind = 3
            nextnode = tree(3,tree(1,j))
         else
            kind = 3
            nextnode = tree(1,j)
         end if
c ----- Move to left child, TYPE 1
      else
         kind = 1
         nextnode = tree(2,j)
      end if
      end

**********************************************************************
      subroutine listtree(n, located, leafid, leaflist, parentlist)
      integer nmxsub, n, located, nextnode
      parameter (nmxsub=1024)
      integer leafid(nmxsub), leaflist(nmxsub), parentlist(nmxsub)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
**********************************************************************
      integer j, jold, next, kind, parentlocate

        located = 0
        parentlocate = 0

        jold = tree(1,n)
        j = n
        do while ( j.ne.tree(1,n) .and. j.ne.tree(3,tree(1,n)) )
           next = nextnode(kind, j, jold)

           if ( kind .eq. 1 ) then
              parentlocate = parentlocate + 1
              parentlist( parentlocate ) = j
           else if ( abs(kind) .eq. 3 ) then
           else
              located = located + 1
              leaflist(located) = j
              leafid(located) = tree(3,j)
           end if

           jold = j
           j = next
        end do

        if ( located-1 .ne. parentlocate ) then
           print*,'Internal Error: located-1 .ne. parentlocate'
        end if
      end

**********************************************************************
      subroutine printtree(n)
      integer nmxsub, n, nextnode
      parameter (nmxsub=1024)
      common /hierarchy/ tree, node, leaf
      integer tree(3,2*nmxsub), node(2*nmxsub), leaf(0:nmxsub)
**********************************************************************
      integer k, j, jold, depth_pr, depth, next, kind

        print'(i3,a,a,i3)',node(1)-2,' nodes used,   ',
     &            'next node is',node(node(1))
        print'(i3,a,a,i3)',leaf(0)-1,' leaves used,   ',
     &            'next leaf is',leaf(leaf(0))

        jold = tree(1,n)
        j = n
        depth = 0
        depth_pr = 0

        do while ( j.ne.tree(1,n) .and. j.ne.tree(3,tree(1,n)) )
           next = nextnode(kind, j, jold)

           if ( kind .eq. 1 ) then
              do k = depth_pr, depth-1
                 print'(a4,$)',':'
              end do
              depth = depth + 1
              depth_pr = depth
              print'(i4,$)',j
           else if ( abs(kind) .eq. 3 ) then
              depth = depth - 1
           else
              do k = depth_pr, depth-1
                 print'(a4,$)',':'
              end do
              print'(i4,a,i3,a)',j,'   (',tree(3,j),' )'
              depth_pr = 0
           end if

           jold = j
           j = next
        end do
      end

