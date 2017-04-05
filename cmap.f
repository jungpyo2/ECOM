!--------------------------------------------------------
! cmap.f: This subroutine is for conformal mapping 
! from the given boundary data to the unit disk using Szego kernel
!--------------------------------------------------------
      subroutine cmapbf 

      use arrays, only: allocarray2
      use arrays, only:  Rmap0, Zmap0, iprintmap 
      use arrays, only:  nt1, ksamp,  nt2, eps, epsLag
      use arrays, only:  Rt1, Zt1, T1, dRt1, dZt1, w1, dwdz1, dwdzz1
      use arrays, only:  Rt2, Zt2, T2, dRt2, dZt2, dwdz2, dwdzz2

      use arrays, only:  dfdzbf,dfdzzbf,fbf,zbf,zmbf, dwdz1m, dwdzz1m 
c
c     read the fixed boundary from a file
      implicit double precision(a-h, o-z)
c      real *8 T2temp(1000 000), Rt2temp(1000 000), Zt2temp(1000 000)
c      real *8 dRt2temp(1000 000), dZt2temp(1000 000)
      !complex *16 dfdz(nt1+1), dfdzz(nt1+1)
      real * 8 curv(6,nt1+1), fr(nt1*ksamp+1), fi(nt1*ksamp+1)
      integer m, nnt2
      complex *16 ima, fixed, w(4*nt1*nt1+10000)
!      complex *16 f(nt1*ksamp+1), w(4*nt1*nt1+10000)
!      complex *16 z(nt1*ksamp+1), zm(nt1*ksamp+1), dwdz1m(nt1*ksamp+1)
!      complex *16 dwdzz1m(nt1*ksamp+1)
c      complex *16 dwdz2temp(1000 000), dwdzz2temp(1000 000)
      
      
 
      write(*,*) 'nt',nt1,ksamp
      m = nt1*ksamp 
 
 
      ! nt2 is the number of theta grid in the unit disk 
      ! to be used for Poisson solver. 
      ! The grid is used for interpolation of the inverse map,
      ! so it may need (nt2 < m) 



      done=1
      pi=atan(done)*4
      ima=(0,1.0d0)
      rl=pi*2
      eps7=1.0d-14

c       ...set the point to be mapped to the origin, make sure it's
c       inside the curve!!!
c

      a=Rmap0 !fixed
      b=Zmap0 !fixed
      fixed=Rmap0+ima*Zmap0

      write(*,*) 'mapping center',fixed
      if (iprintmap.eq.1) then
         itype1=1
         itype2=2
         itype3=3
         iw=11
         call quaplot2(iw,Rt1,Zt1,nt1,itype2,a,b,1,itype2,
     1        'equispaced nodes and zero-mapped point*')
      end if
c------------------------------------
c
c       now compute the forward map
c
c      allocate (curv(6,nt1+1)); 
      curv = 0.0
      curv(1,:)=Rt1
      curv(2,:)=dRt1
      curv(4,:)=Zt1
      curv(5,:)=dZt1
      rltot = t1(nt1+1)
      lw=4*nt1*nt1+10000 !need to be larger than 4*nt1**2 !20 000 000
      write(*,*) 'nt1-2',nt1,ksamp,lw,w(lw)
!      allocate (f(m), f2(m), w(lw))
      call szego_fft0(nt1,curv,rltot,ksamp,fixed,crowd)
      write(*,*) 'dfdzbf',dfdzbf(1:nt1)
      do i=1,nt1
         w1(i)=fbf(1+(i-1)*ksamp)
      end do
      do i=1,nt1
         dwdz1(i)=dfdzbf(i)
         dwdzz1(i) = dfdzzbf(i)
      end do
      d=-1
      do 2400 i=1,m
         dd=1-abs(fbf(i))
         dd=abs(dd)
         if (dd .gt. d) d=dd
           
 2400 continue
c
      write (*,*) 'farthest distance off unit disc=',d
      write (*,*) 'odonnell-rokhlin crowding factor=',crowd
c

      hmin=1
      hmax=-1
c
      do 3100 i=1,nt1-1
         d1=abs(fbf(i+1)-fbf(i))
         if (d1 .lt. hmin) hmin=d1
         if (d1 .gt. hmax) hmax=d1
 3100 continue
c
      write (*,*) 'hmin=*',hmin
      write (*,*) 'hmax=*',hmax
      write (*,*) 'crowding factor hmax/hmin=', hmax/hmin
c
      
!      allocate(fr(m),fi(m))
      do 3200 i=1,m
         fr(i)=fbf(i)
         fi(i)=-ima*fbf(i)
 3200 continue
c
      if (iprintmap.eq.1) then
         iw=12     
       write(*,*) 'test1'

         call quaplot(iw,fr,fi,m,itype2,'conformal map via fft*')
      end if

c------------------------------------
c       compute the inverse map by interpolation
c       zero pad and fft the original curve so that it is of length m
c

      do 4800 i=1,nt1+1
         zbf(i)=curv(1,i)+ima*curv(4,i)
 4800 continue

      m=nt1*ksamp

      call interpfft(nt1,zbf,m,zmbf,w)
      call interpfft(nt1,dwdz1,m,dwdz1m,w)
      call interpfft(nt1,dwdzz1,m,dwdzz1m,w)
c
c       interpolate the whole inverse map, z=h(w), w on the unit disc
c
!        eps=1.0d-10
        k=7
        nmax=1000
 
 
        nnt2=nt2+1

        call cnfinv_init(
     1       nt2+1,T2, Rt2, Zt2, dRt2, dZt2,dwdz2, dwdzz2
     2       ,ier,eps,m,k,nmax,errmax,w,lw,lused)

        write(*,*) 'after cnfinc_init, ier=',ier
        write(*,*) 'after cnfinc_init, errmax=*',errmax

      if (iprintmap.eq.1) then
         iw=20
         call quaplot2(iw,Rt1,Zt1,nt1+1,itype2,Rt2,Zt2,nt2+1,itype2,
     1        'Lagrange interpolation for uniform arc length*')
      end if
c      write(*,*) 'test312', RT2(1:nt2+1)
c      write(*,*) 'test313', T2
c
c       produce a bunch of rays and rings on the disc, then plot 
c       their image under the inverse map
c
      if (iprintmap.eq.1) then
        iw=31
        nr1=5
        ntheta1=100
        nr2=50
        ntheta2=30
        npts=750
   
 !       npts=750
        do i=1,nt2+1
           w(i)=Rt2(i)+ima*Zt2(i)
        end do

        call plotgrid(iw,nt1+1,z,nr1,ntheta1,nr2,ntheta2,nt2,w)
c
      end if

      return
      end !end of the subrouine cmap

c---------------------------------------------------------------
      subroutine cmapin

      use arrays, only:  allocarray3, iprintmap
      use arrays, only:  Rt3, Zt3, rnd, tnd, dzdw3, dzdww3
      use arrays, only:  nt1, ksamp,  nt2, eps, epsiter, maxiter
      use arrays, only:  nsub, kcheb, kLag, nr, ntot 
      use arrays, only:  Rt1, Zt1, T1, dRt1, dZt1, Rmap0, Zmap0
      use arrays, only:  Rt2, Zt2, T2, dRt2, dZt2, dwdz2, dwdzz2
      use arrays, only:  zk, dzdw2k, dzdww2k, dwdz3, dwdzz3
      use arrays, only:  blength, bnodes
      use arrays, only:  epslogzk, epslogzwk, epslogzwwk 
c
c     read the fixed boundary from a file
      implicit double precision(a-h, o-z)
c      parameter (nmxsub = 200, kmax = 32)
c      parameter (nr = nsub * kcheb)
c      parameter (ntmax = 300)

      
 
      
c      parameter (ntcon = ntheta*ncon)

c      real *8 g(nt2)
 

!      real *8 wsave(4*nt2+100)
      real *8 zw(3*ntot+6*nt2+100)
  


c      real *8 rconst(nsub*nt2)
c      real *8 tconst(nsub*nt2)
c      real *8 uconst(nsub*nt2)

      real *8 wsave(10 000)

      complex *16 ima, wz, cint, wcon, z(nt2+1), dzdw2(nt2+1)
      complex *16 dzdww2(nt2+1), temp1(nt2), temp2(nt2), temp3(nt2)
      integer sol_unit
      real *8 maxdpsi,  maxpsi, lambda
      character(3) ko, ns
      character(4) nt
      logical checkerr
      
      save

      checkerr = .false.
      write(*,*) 'start poisson, nt2',nt2

      write(*,*) 'start poisson, korder,nr,ntot',kcheb,nr,ntot
      pi = 4.0d0*datan(1.0d0)
      ima = (0,1)

      ntheta=nt2

      korder = kcheb
c
c     Set number of points in theta discretization 
c
c      ntheta = 16
c
      call dcffti(ntheta,wsave)
      do i = 1,ntheta
         tnd(i) = (i-1)*2*pi/ntheta
      enddo


c
c     set grid in r.
c
c      nsub = 4
c      korder = 16
c      nr = korder*nsub
      rad = 1.0d0
      a = 0.0d0
      c = rad
c      write(*,*) 'start poisson, nr',nr
      call chrad(nsub,korder,a,c,blength,bnodes,rnd,delrmin)

c      do inext = 1, ntot
c      if(int(inext/nr).eq.10) then
c         write(*,*) 'inext, nr', inext,nr      
c         write(*,*) 'rnd0(inext) (inext)',rnd(inext)
c      end if
c      end do
      write(6,*)' bnodes are',(bnodes(j),j=1,nsub+1)

c     conformal backward mapping to the original domain
!      write(*,*) 'Rt2', rt2(1:nt2+1)
      do i=1,nt2+1
         z(i)=Rt2(i)+ima*Zt2(i)
         dzdw2(i) = (1,0)/dwdz2(i)
         dzdww2(i) = -dwdzz2(i)/(dwdz2(i)**3)
      end do
!      write(*,*) 'dzdw2', dzdw2(1:nt2+1)

      call fft_wk(nt2, z, zk, wsave)
  
      call fft_wk(nt2, dzdww2, dzdww2k, wsave)
!      write(*,*) 'dzdww2',dzdww2(1:nt2)
!      write(*,*) 'dzdww2k',abs(dzdww2k(1:nt2))

      call fft_wk(nt2, dzdw2, dzdw2k, wsave)
!      write(*,*) 'dzdw2',dzdw2(1:nt2)
!      write(*,*) 'dzdw2k',abs(dzdw2k(1:nt2))
      if (iprintmap.eq.1) then
         iw=27
         itype=3
         call quagraph2(iw,t2,realpart(dzdww2),nt2,
     1        itype,t2,imagpart(dzdww2),nt2,itype,
     2       'real and imag parts of the interpolated inverse map*')
         
         iw=26
         itype=3
         call quagraph2(iw,t2,realpart(dzdw2),nt2,
     1        itype,t2,imagpart(dzdw2),nt2,itype,
     2        'real and imag parts of the interpolated inverse map*')

         iw=25
         itype=3
         call quagraph2(iw,t2,abs(dzdw2k),nt2,
     1      itype,t2,abs(dzdww2k),nt2,itype,
     2       'real and imag parts of the interpolated inverse map*')
c      write(*,*) 'wk',wk(1:nt2+1)
      end if
      call runtime(' started fft cauchy ')

      if (nt2.gt.256) then
         epsmapin=1.0d-30
      else if (nt2.gt.128) then
         epsmapin=1.0d-25
      else
         epsmapin=1.0d-20
      end if
      epslogzk=log(epsmapin/maxval(abs(zk)))
      epslogzwk=log(epsmapin/maxval(abs(dzdw2k)))
      epslogzwwk=log(epsmapin/maxval(abs(dzdww2k)))

c      write(*,*) 'zkabs',abs(zk(1:nt2+1))
c      write(*,*) 'zwwkabs',abs(dzdww2k(1:nt2+1))
c      abszk=maxval(abs(zk(1:nt2+1)))               !abs(wk)
c      absdzdwk=maxval(abs(dzdw2k(1:nt2+1)))
c      absdzdwwk=maxval(abs(dzdww2k(1:nt2+1)))

   
      do i = 1,nr
c         nlimz=0
c         nlimzw=0
c         nlimzww=0
         nlimz=epslogzk/log(rnd(i))+1
         nlimzw=epslogzwk/log(rnd(i))+1
         nlimzww=epslogzwwk/log(rnd(i))+1
c         do k=1,nt2 
c            ntmp=log(epsmapin/abs(zk(k)))/log(rnd(i))+1
c            if (ntmp.gt.nlmz) then
c               nlimz=ntmp
c            end if
c            ntmp=log(epsmapin/abs(dzdw2k(k)))/log(rnd(i))+1
c            if (ntmp.gt.nlmzw) then
c               nlimzw=ntmp
c            end if
c            ntmp=log(epsmapin/abs(dzdww2k(k)))/log(rnd(i))+1
c            if (ntmp.gt.nlmzww) then
c               nlimzww=ntmp
c            end if
c         end do
         if (nlimz.gt.nt2) then
            nlimz=nt2
         end if
         if (nlimzw.gt.nt2) then
            nlimzw=nt2
         end if
         if (nlimzww.gt.nt2) then
            nlimzww=nt2
        end if

        call fft_cauchy1D(nt2,rnd(i),zk,wsave,temp1)
        call fft_cauchy1D(nt2,rnd(i),dzdw2k,wsave,temp2)
        call fft_cauchy1D(nt2,rnd(i),dzdww2k,wsave,temp3)
        do j = 1,ntheta
           inext = i + (j-1)*nr
           Rt3(inext)=temp1(j)
           Zt3(inext)=-ima*temp1(j)
           dzdw3(inext) = temp2(j)
           dwdz3(inext) =  (1.0d0,0.0d0)/dzdw3(inext)
           dzdww3(inext) = temp3(j)
           dwdzz3(inext) =-dzdww3(inext)/(dzdw3(inext)**3)
        end do

c        do j = 1,ntheta
c            inext = i + (j-1)*nr
    
c         write(*,*) 'i,j',i,j,rnd(i),tnd(j)
c            x(inext) = rnd(i)*cos(tnd(j))
c            y(inext) = rnd(i)*sin(tnd(j))
c         write(*,*) 'inext,x,y',inext,x(inext),y(inext)
c            wz =rnd(i)*exp(ima*tnd(j))

c            call fft_cauchy(nlimz,wz,zk,cint)

c            write(*,*) 'diff1',cint-complex(Rt3(inext),Zt3(inext))

c            Rt3(inext)=cint
c            Zt3(inext)=-ima*cint
c            call fft_cauchy(nlimzw,wz,dzdw2k,cint)

c            write(*,*) 'diff2',cint-dzdw3(inext)
c            dzdw3(inext) = cint
c            dwdz3(inext) = (1.0d0,0.0d0)/dzdw3(inext)
c            call fft_cauchy(nlimzww,wz,dzdww2k,cint)
c            dzdww3(inext) = cint
c            dwdzz3(inext) = -dzdww3(inext)/(dzdw3(inext)**3)

c         enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c         This is the end of the debugging code and the beginning
c         of the conformal mapping code
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine szego_fft(n,curv,rltot,ksamp,fixed,
     1      crowd,w,lw)

        use arrays, only:  dfdzbf,dfdzzbf, fbf
        implicit real *8 (a-h,o-z)
        real *8 curv(6,1)
        complex *16 fixed,w(lw)
c
c       this routine conformally maps equispaced points in arclength
c       on the boundary of a curve to the boundary of the unit disc
c       using fancy szego kernel - kerzman-stein machinary, and then
c       integrating the resulting derivative of the map via an fft.
c       this is a memory management routine, the real work is done
c       in szego_fft0
c
c       input:
c
c         n - the number of points in the curve
c         curv - a 6 x n array that contains information about the
c             curve to map as follows:
c             curv(1,i) - the x coordinate of the ith point
c             curv(2,i) - the value of dxdt, where t is the
c                         parameterization variable of the curve
c             curv(3,i) - not used
c             curv(4,i) - the y coordinate of the ith point
c             curv(5,i) - the value of dydt, where t is the
c                         parameterization variable of the curve
c             curv(6,i) - not used
c
c             NOTE: it is assumed that the points x,y in curv are
c             sampled in terms of arclength, it is NOT ASSUMED that
c             the tangent vectors dxdt,dydt are of unit length
c
c             NOTE: curv should not be closed, that is to say that
c             curv(j,1) should NOT be equal to curv(j,n)
c
c         rltot - the length of the curve
c         ksamp - the factor by which to oversample upon integration
c             via fft, start with 100 points end with ksamp*100 points
c         fixed - the point inside curv which to map to the origin
c         lw - the length of the work array w in real *8 words
c
c       output:
c
c         f - the conformally mapped points that lie on the unit disc,
c             i.e, f(i) is the image of curv(1,i),curv(4,i)
c         crowd - the crowding factor, as defined in odonnell-rokhlin,
c             min(f')/max(f')
c         w - work array, it's up to you to make sure it's long enough,
c             it should be something like 4*n**2+20*n real *8 words long
c
c
c       memory manage...
c
        write(*,*) 'nt1-3',n,ksamp,lw
        iz=1
        lz=n+10
c
        idzds=iz+lz
        ldzds=n+10
c
        irhs=idzds+ldzds
        lrhs=n+10
c
        isol=irhs+lrhs
        lsol=n+10
c
        iw7=isol+lsol
        lw7=2*n**2+10*n+10000
c
        iw2=iw7+lw7
        lw2=lw-iw2
c
        call runtime(' before szego_fft0 ')
c        call szego_fft0(n,curv,rltot,ksamp,fixed,
C     1     crowd,w(iz),w(idzds),w(irhs),w(isol),w(iw7),
C     2      lw7,w(iw2),lw2)
c
c        write(*,*) 'dfdz0',dfdz(1:n)
        return
        end
c
c
c
        subroutine szego_fft0(n,curv,rltot,ksamp,fixed,
     1     crowd)

        use arrays, only:  dfdzbf,dfdzzbf, fbf
        use arrays, only: iprintmap, iiter, wsaved, w7saved
        use arrays, only: w7, w8, w9

        implicit real *8 (a-h,o-z)
        real *8 curv(6,1),ttt(n),ys(n),xs(n)
        complex *16 fixed,z(n)
        complex *16 dzds(n),rhs(n), dfdss(n),dsdz(n),dsdzs(n) 
     1       ,sol(n),ima,offset
c
c
         write(*,*) 'nt1-4',n,ksamp
        done=1
        pi=4*atan(done)
        ima=(0,1)
c
c       since already sampled by arclength, just make sure |dfdt|=1
c
        do 1200 i=1,n
        dd=curv(2,i)**2+curv(5,i)**2
        dd=sqrt(dd)
        curv(2,i)=curv(2,i)/dd
        curv(5,i)=curv(5,i)/dd
 1200   continue
c
c       get z and dzds for rhs construction
c
        do 1400 i=1,n
        ttt(i)=i
        z(i)=curv(1,i)+ima*curv(4,i)
        dzds(i)=curv(2,i)+ima*curv(5,i)
 1400   continue
c
        do 2000 i=1,n
        rhs(i)=conjg((dzds(i))/2/pi/ima/(z(i)-fixed))
 2000   continue

c       if 'cmapbf' routine is initially called (iiter.eq.1)
c       construct matrix and solve the thing directly 
c       using old-faithful
c
        if (iiter.eq.1) then

        do 2800 i=1,n
        do 2600 j=1,n
c
        if (i .eq. j) goto 2600
c
        w9(i,j)=conjg(dzds(i)/(z(i)-z(j))/2/pi/ima)
        w9(i,j)=w9(i,j)-dzds(j)/(z(j)-z(i))/2/pi/ima
        w9(i,j)=rltot*w9(i,j)/n
c
 2600   continue
c
        w9(i,i)=1
 2800   continue
c
        call runtime(' before cqrdecom   ')
        call cqrdecom(w9,n,w7,rcond)
        call runtime(' after cqrdecom   ')
        write(*,*) 'rcond=*',rcond
c
        do i=1,n
        do j=1,n
           wsaved((i-1)*n+j)=w9(i,j)
           w7saved((i-1)*n+j)=w7((i-1)*n+j)
        end do
        end do

        else   !if 'cmapbf' routine is recalled (iiter.ne.1)
               ! skip the QR decomposition, instead restore the saved w and w7
        do i=1,n
        do j=1,n
           w9(i,j)=wsaved((i-1)*n+j)
           w7((i-1)*n+j)=w7saved((i-1)*n+j)
        end do
        end do

        end if

        call runtime(' before cqrsolve   ')
        call cqrsolve(n,w7,rhs,sol)
c
        call cmatvec(n,w9,sol,w7)
        call runtime(' after cqrsolve    ')
c
        cd1=0
        cd2=0
        do 3000 i=1,n
        cd1=cd1+abs(w7(i)-rhs(i))**2
        cd2=cd2+abs(rhs(i))**2
 3000   continue
c
        cd1=sqrt(cd1/cd2)
        write(*,*) 'residual from qr=*',cd1
c
c       calculate fprime (w7 below) and then integrate it to get f
c
        saa=0
        do 3400 i=1,n
        saa=saa+sol(i)*conjg(sol(i))
 3400   continue
c
        saa=saa*rltot/n
c
        do 3800 i=1,n
        w7(i)=2*pi*sol(i)**2/saa
        xs(i)=w7(i)
        ys(i)=-ima*w7(i)
        dfdzbf(i)= w7(i)
 3800   continue

        
c
        if (iprintmap.eq.1) then
           iw=15
           itype=3
           call quagraph2(iw,ttt,xs,n,itype,ttt,ys,n,itype,
     1       'real and imag parts of dfdz*')
        end if
c
c       calculate the crowding factor
c
        cmin=1.0d20
        cmax=-1
        do 4000 i=1,n
        if (abs(w7(i)) .lt. cmin) cmin=abs(w7(i))
        if (abs(w7(i)) .gt. cmax) cmax=abs(w7(i))
 4000   continue
c
        crowd=cmin/cmax

c
c       and now scale w7 to get dfds...
c
        do 4400 i=1,n
        w7(i)=w7(i)*dzds(i)*rltot/2/pi
        xs(i)=w7(i)
        ys(i)=-ima*w7(i)
        dsdz(i)=(1.0,0.0)/(dzds(i)*rltot/2/pi)
         
 4400   continue
c        write(*,*) 'dzds=',dzds(1:n)
c        write(*,*) 'dsdz=',dsdz
c
        if (iprintmap.eq.1) then
           iw=16
           itype=3
        call quagraph2(iw,ttt,xs,n,itype,ttt,ys,n,itype,
     1       'real and imag parts of dfds*')
      end if
c
c       ...differentiate dfds to find dfdss via an fft
c       and use it to find dfdzz=dfdss*dsdz**2+dfds*d(dsdz)ds*dsdz

c        write(*,*) 'w7=',w7(1:n)
        call cdifffft(n,w7,dfdss,w8)
c        write(*,*) 'dfdss=',dfdss
        call cdifffft(n,dsdz,dsdzs,w8)
c        write(*,*) 'dsdzs=',dsdzs
 
        do i=1,n
        dfdzzbf(i)=dfdss(i)*(dsdz(i)*dsdz(i))
     1       +w7(i)*dsdzs(i)*dsdz(i)
        end do
c        write(*,*) 'dfdzz=',dfdzz
c
c       ...and now integrate dfds via an fft
c
        m=n*ksamp
        call cintfft(n,w7,m,fbf,w8)

c
c       compensate for off center cicle - 
c       use formula for center of a circle given 3 points on the circle
c
        x1=fbf(1)
        y1=-ima*fbf(1)
        r1=fbf(1)*conjg(fbf(1))
c
        do 5000 i=2,m
        dis2=(fbf(1)-fbf(i))*conjg(fbf(1)-fbf(i))
        if (dis2 .ge. 1.0d0) goto 5200
 5000   continue
 5200   continue
c
        x2=fbf(i)
        y2=-ima*fbf(i)
        r2=fbf(i)*conjg(fbf(i))
c
        do 5600 j=i,m
        dis2=(fbf(j)-fbf(i))*conjg(fbf(j)-fbf(i))
        if (dis2 .ge. 1.0d0) goto 5800
 5600   continue
 5800   continue
c
        x3=fbf(j)
        y3=-ima*fbf(j)
        r3=fbf(j)*conjg(fbf(j))
c
        d=2*(x1*y2-x1*y3+y1*x3-y1*x2+x2*y3-y2*x3)
        tx=(r1*y2-r1*y3+y1*r3-y1*r2+r2*y3-y2*r3)
        ty=(x1*r2-x1*r3+r1*x3-r1*x2+x2*r3-r2*x3)
c
        offset=tx/d+ima*ty/d
c
        do 6200 i=1,m
        fbf(i)=fbf(i)-offset
 6200   continue
c
        return
        end
c
c
c
c
c
        subroutine cmatvec(n,a,x,y)
        implicit real *8 (a-h,o-z)
        complex *16 a(n,n),x(1),y(1),cd
c
        do 1600 i=1,n
        cd=0
        do 1400 j=1,n
        cd=cd+a(i,j)*x(j)
 1400   continue
        y(i)=cd
 1600   continue
c
        return
        end
c
c
c
c
c
       subroutine cintfft(n,dfdt,m,f,work)
        implicit real *8 (a-h,o-z)
        dimension work(1)
        complex *16 dfdt(n),f(m)
c
c       This subroutine takes the FFT of the input, zero pads it with
c       m-n zeros such that the new length is m, and then divides by
c       i*wave_number and takes the inverse FFT. This results in the
c       integration and oversampling of the input.
c       
c              input parameters:
c
c       n - current number of sample points
c       m - desired number of sample points
c       dfdt - input to be integrated
c       work - work array that must be at least length 4*m+2*n+15 real *8
c
c              output parameters:
c
c       f - integrated dfdt (zero padded with m-n zeros)
c
c
        ncir=2
        iwsave=1
        lwsave=4*m+15
c
        idfdt2=iwsave+lwsave
        ldfdt2=ncir*n
c
        call cintfft0(n,dfdt,m,f,work(iwsave),work(idfdt2))
c
        return
        end
c
c
c
        subroutine cintfft0(n,dfdt,m,f,wsave,dfdt2)
        implicit real *8 (a-h,o-z)
        dimension wsave(1)
        complex *16 dfdt(n),f(m),ima,dfdt2(n)
        integer n,m
c
c       initialize f
c
        do 1200 i=1,m
        f(i)=0
 1200   continue
c
        ima=(0,1)
        done=1
        pi=4*atan(done)
c
        do 1400 i=1,n
        dfdt2(i)=dfdt(i)
 1400   continue
c
        call dcffti(n,wsave)
        call dcfftf(n,dfdt2,wsave)
cccc        call prin2('fft of dfds=*',dfdt2,2*n)
c
c       divide by i*wave_number
c
        do 1600 i=1,n/2-1
c
        j=i+1
        dfdt2(j)=dfdt2(j)/ima/i
      
c
        j=n-i+1
        dfdt2(j)=-dfdt2(j)/ima/i
      
c        write(*,*) dfdt3(j)
 1600 continue
c
c       zero pad the result and take the inverse fft now
c
        do 1700 i=1,n/2+1
        f(i)=dfdt2(i)/n
 1700   continue
c
        do 1750 i=1,n/2
        f(m-i+1)=dfdt2(n-i+1)/n
      
 1750 continue
c      write(*,*) dfdtt
c
        call dcffti(m,wsave)
        call dcfftb(m,f,wsave)
c
        return
        end

        subroutine cdifffft(n,f,dfdt,work)
        implicit real *8 (a-h,o-z)
        dimension work(1)
        complex *16 f(n),f2(n), dfdt(n)
c
c       This subroutine takes the FFT of the input,and FFT differentiate it. 
c       Then, take the inverse FFT to return dfdt
c
c              input parameters:
c
c       n - number of sample points
c       f - input to be differentiated
c       work - work array that must be at least length 4*n+15 real *8
c
c              output parameters:
c
c       dfdt - differntiated f 
c
c
        ncir=2
        iwsave=1
        lwsave=4*n+15
c
        idfdt2=iwsave+lwsave
        ldfdt2=ncir*n

        call cdifffft0(n,f,dfdt,work(iwsave),work(idfdt2))
c
        return
        end
c
c
c
        subroutine cdifffft0(n,f,dfdt,wsave,f2)
        implicit real *8 (a-h,o-z)
        dimension wsave(1)
        complex *16 dfdt(n),f(n),ima,f2(n)
        integer n
c
c       initialize f
c
c
        do i=1,n
        dfdt(i) = (0.0d0, 0.0d0)
        f2(i)=f(i)
        end do

c        write(*,*) 'f2=',f2
c
        ima=(0,1)
        done=1
        pi=4*atan(done)
c
 
c
        call dcffti(n,wsave)
        call dcfftf(n,f2,wsave)
cccc        call prin2('fft of dfds=*',dfdt2,2*n)
c
c       divide by i*wave_number
c
        do i=1,n/2-1
c
        j=i+1
        f2(j)=f2(j)*(ima*i)
        j=n-i+1
        f2(j)=-f2(j)*(ima*i)
        
        end do

        f2(1) = (0.0d0,0.0d0)
 
        do i=1,n
        dfdt(i)=f2(i)/n
        end do
c
        call dcffti(n,wsave)
        call dcfftb(n,dfdt,wsave)
c
        
c        write(*,*) 'dfdt=',dfdt
        return
        end
c
c
c
c
c

        subroutine interpfft(n,fin,m,f,work)
        implicit real *8 (a-h,o-z)
        dimension work(1)
        complex *16 fin(n),f(m)
c
c       This subroutine takes the FFT of the input, zero pads it with
c       m-n zeros such that the new length is m, and takes the inverse 
c       FFT. This results in the oversampling of the input.
c       
c              input parameters:
c
c       n - current number of sample points
c       m - desired number of sample points
c       fin - input to be overesampled
c       work - work array that must be at least length 4*m+2*n+15 real *8
c
c              output parameters:
c
c       f - oversampled fin
c
c
        ncir=2
        iwsave=1
        lwsave=4*m+15
c
        ifin2=iwsave+lwsave
        lfin2=ncir*n
c
        call interpfft0(n,fin,m,f,work(iwsave),work(ifin2))
c
        return
        end
c
c
c

        subroutine interpfft0(n,fin,m,f,wsave,fin2)
        implicit real *8 (a-h,o-z)
        real *8 wsave(1)
        complex *16 fin(n),f(m),ima,fin2(n)
c
c       initialize f
c
        do 1200 i=1,m
        f(i)=0
 1200   continue
c
        ima=(0,1)
        done=1
        pi=4*atan(done)
c
        do 1400 i=1,n
        fin2(i)=fin(i)
 1400   continue
c
        call dcffti(n,wsave)
        call dcfftf(n,fin2,wsave)
c
c       zero pad and take the inverse fft
c
        fre=0.0d0
        fim=0.0d0
        do 1800 i=1,n/2+1
        f(i)=fin2(i)/n  !*(1.0d0,0.0d0)
        fre=fre+real(fin2(i))
        fim=fim+real(fin2(i)*(0.0d0,-1.0d0))
 1800   continue
c
c        write(*,*) 'fre,fim',fre,fim
        do 2200 i=1,n/2
        f(m-i+1)=fin2(n-i+1)/n
 2200 continue
c
        call dcffti(m,wsave)
        call dcfftb(m,f,wsave)
c
        return
        end

        subroutine interpffte(n,fin,m,f,work)
        implicit real *8 (a-h,o-z)
        dimension work(1)
        complex *16 fin(n),f(m)
c
c       This subroutine takes the FFT of the input, zero pads it with
c       m-n zeros such that the new length is m, and takes the inverse 
c       FFT. This results in the oversampling of the input.
c       
c              input parameters:
c
c       n - current number of sample points
c       m - desired number of sample points
c       fin - input to be overesampled
c       work - work array that must be at least length 4*m+2*n+15 real *8
c
c              output parameters:
c
c       f - oversampled fin
c
c
        ncir=2
        iwsave=1
        lwsave=4*m+15
c
        ifin2=iwsave+lwsave
        lfin2=ncir*n
c
        call interpffte0(n,fin,m,f,work(iwsave),work(ifin2))
c
        return
        end
c
c
c

        subroutine interpffte0(n,fin,m,f,wsave,fin2)
        implicit real *8 (a-h,o-z)
        real *8 wsave(1)
        complex *16 fin(n),f(m),ima,fin2(n)
c
c       initialize f
c
        do 1200 i=1,m
        f(i)=0
 1200   continue
c
        ima=(0,1)
        done=1
        pi=4*atan(done)
c
        do 1400 i=1,n
        fin2(i)=fin(i)
 1400   continue
c
        call dcffti(n,wsave)
        call dcfftf(n,fin2,wsave)
c
c       zero pad and take the inverse fft
c
        fre=0.0d0
        fim=0.0d0
        do 1800 i=1,n/2+1
        f(i)=real(fin2(i)/n)*(1.0d0,0.0d0)
c        fre=fre+real(f(i))
c        fim=fim+real(f(i)*(0.0d0,-1.0d0))
 1800   continue
c
c        write(*,*) 'fre,fim',fre,fim
        do 2200 i=1,n/2
        f(m-i+1)=fin2(n-i+1)/n
 2200 continue
c
        call dcffti(m,wsave)
        call dcfftb(m,f,wsave)
c
        return
        end

      subroutine interpffto(n,fin,m,f,work)
        implicit real *8 (a-h,o-z)
        dimension work(1)
        complex *16 fin(n),f(m)
c
c       This subroutine takes the FFT of the input, zero pads it with
c       m-n zeros such that the new length is m, and takes the inverse 
c       FFT. This results in the oversampling of the input.
c       
c              input parameters:
c
c       n - current number of sample points
c       m - desired number of sample points
c       fin - input to be overesampled
c       work - work array that must be at least length 4*m+2*n+15 real *8
c
c              output parameters:
c
c       f - oversampled fin
c
c
        ncir=2
        iwsave=1
        lwsave=4*m+15
c
        ifin2=iwsave+lwsave
        lfin2=ncir*n
c
        call interpffto0(n,fin,m,f,work(iwsave),work(ifin2))
c
        return
        end
c
c
c

        subroutine interpffto0(n,fin,m,f,wsave,fin2)
        implicit real *8 (a-h,o-z)
        real *8 wsave(1)
        complex *16 fin(n),f(m),ima,fin2(n)
c
c       initialize f
c
        do 1200 i=1,m
        f(i)=0
 1200   continue
c
        ima=(0,1)
        done=1
        pi=4*atan(done)
c
        do 1400 i=1,n
        fin2(i)=fin(i)
 1400   continue
c
        call dcffti(n,wsave)
        call dcfftf(n,fin2,wsave)
c
c       zero pad and take the inverse fft
c
        fre=0.0d0
        fim=0.0d0
        do 1800 i=1,n/2+1
        f(i)=aimag(fin2(i)/n)*(0.0d0,1.0d0)
c        fre=fre+real(f(i))
c        fim=fim+real(f(i)*(0.0d0,-1.0d0))
 1800   continue
c
c        write(*,*) 'fre,fim',fre,fim
        do 2200 i=1,n/2
        f(m-i+1)=fin2(n-i+1)/n
 2200 continue
c
        call dcffti(m,wsave)
        call dcfftb(m,f,wsave)
c
        return
        end

      subroutine cnfinv_init(nt2temp,T2temp,Rt2temp,Zt2temp,
     1      dRt2temp,dZt2temp, dwdz2temp,dwdzz2temp,
     2            ier,eps,m,k,nmax,errmax,w,lw,lused)
      
        use arrays, only:  fbf, zmbf, dwdz1m, dwdzz1m
c        use arrays, only:  Rt2, Zt2, T2, dRt2, dZt2
        use arrays, only:  kLag, epsLag, iprintmap
c        use arrays, only:  IntLag_periodic2 

        implicit real *8 (a-h,o-z)
        save
        real *8 w(lw), t(m+1)
        integer nt2temp, i2pi
!        real *8 T2o(n2o), R2o(n2o), Z2o(n2o), dR2o(n2o), dZ2o(n2o)
        real *8 T2temp(nt2temp),Rt2temp(nt2temp),Zt2temp(nt2temp)
        real *8 dRt2temp(nt2temp),dZt2temp(nt2temp)
        complex *16 ima,zz
        !complex *16 f(m+1),z(m+1), dwdz1(m+1)
        !complex *16 dwdzz1(m+1),
        complex *16  dwdz2temp(m+1),dwdzz2temp(m+1)
        complex *16 ftemp(m+1),ztemp(m+1),dwdz1temp(m+1),dwdzz1temp(m+1)
!        real *8, allocatable :: x(:),y(:), tin(:)
        real *8 ttemp(m+1),min2pi
        real *8 x(m+1),y(m+1), tin(m+1),dxdz(m+1),dydz(m+1)
        real *8 dxdzz(m+1),dydzz(m+1)
        real *8 dxdz2(m+1),dydz2(m+1)
        real *8 dxdzz2(m+1),dydzz2(m+1)
        real *8 ddxdzdt2(m+1),ddydzdt2(m+1)
        real *8 ddxdzzdt2(m+1),ddydzzdt2(m+1)
        integer *4 isort(m+1)
c
c       this subroutine takes a more or less random, but well
c       oversampled, discritization of a forward conformal map
c       and interpolates it properly so as to compute the
c       inverse map
c
c       input:
c
c         eps - the precision to which the inverse map should
c             be interpolated
c         m - length of f and z
c         f - the image of the forward map on the unit circle
c         z - the curve that was mapped to the unit circle
c         k - the degree polynomial that should be used on each
c             subinterval to interpolate the inverse map, somethingc
c             like 5 or 7 or 8 should be used, k should not be large
c         nmax - used in subdivision of f, the maximum number of points
c             that should be interpolated using a k-th degree
c             polynomial
c         w - large work array
c         lw - length of w in real *8 words
c         kLag - the order of Lagrange interpolation
c         nt2 - the number of grid points to be interpolated
c
c       output:
c
c         ier - ier=0 is good, anything else is bad
c         errmax - the maximum interpolation error over the whole
c             unit circle
c         w - contains info to be used by cnfinv_eval, etc
c         lused - the first lused real *8 elements of w should NOT
c             be touched afterwards
c
c       NOTE: f and z should NOT be closed, i.e. f(1) should not
c           equal f(m) and z(1) should not equal z(m)
c       NOTE: f and z may have an additional item added to the end
c
c 
        write(*,*) 'Start mapping'
        ier=0
        ima=(0,1)
        pi=4*atan(1.0d0)

 

        if (m .gt. 1000 000) then
            ier=1024
            write (*,*) 'need to increase static memory',a
            return
        endif
c
c       turn f into thetas on the unit disc and ensure that the
c       inverse map is defined on exactly [0+a,2*pi+a], for some a
c
        zz=log(fbf(1))
        theta=-ima*zz
        ifpos=0        
        if (theta .ge. 0) ifpos=1
c

        do 1000 i=1,m
c
        yy=-ima*(fbf(i))
        zz=log(fbf(i))
        theta=-ima*zz
c
        if ((ifpos .eq. 1) .and. (yy .lt. 0)) then
            theta=theta+2*pi
        endif
c
        if ((ifpos .eq. 0) .and. (yy .ge. 0)) then
            ifpos=1
        endif
c
        t(i)=theta
 1000   continue

c
        d=t(m)-t(1)
        if (d .ge. 2*pi) goto 1200
c
c        write(*,*) "tcmap0",t(1:m)
        min2pi=10.0d0
        do i=1, m
           isort(i)=i
           diff2pi=dabs(t(i)-2*pi)
c            write(*,*) 'isort',isort(i),tmaxp(isort(i)),diff2pi
           if (diff2pi .lt. min2pi) then
               min2pi=diff2pi
               i2pi=i
c               write(*,*) 'i2pi',i2pi,diff2pi
           endif
        end do
        if(i2pi.ne.m) then
           do i=1, m-i2pi
              isort(i)=i2pi+i
           end do
           do i=1, i2pi
              isort(i+m-i2pi)=i
           end do
        end if
c        call qsortinc(m,t(1:m),isort(1:m))  
  
        do i=1, m
           ttemp(i)=t(isort(i))
           ftemp(i)=fbf(isort(i))
           ztemp(i)=zmbf(isort(i))
           dwdz1temp(i)=dwdz1m(isort(i))
           dwdzz1temp(i)=dwdzz1m(isort(i))
        end do
        t(1:m)=ttemp(1:m)
        fbf(1:m)=ftemp(1:m)
        zmbf(1:m)=ztemp(1:m)
        dwdz1m(1:m)=dwdz1temp(1:m)
        dwdzz1m(1:m)=dwdzz1temp(1:m)
c        write(*,*) "tcmap",t(1:m)
c        write(*,*) "isortcmap",isort(1:m)
        mm=m+1
        t(mm)=t(1)+2*pi
        fbf(mm)=fbf(1)
        zmbf(mm)=zmbf(1)
        dwdz1m(mm)=dwdz1m(1)
        dwdzz1m(mm)=dwdzz1m(1)
 1200   continue

c
c       ...interpolate the real part of the inverse map
c
!        allocate (x(mm),y(mm),tin(mm))
        do 1400 i=1,mm
        x(i)=zmbf(i)
        y(i)=-ima*(zmbf(i))
        dxdz(i)=dwdz1m(i)
        dydz(i)=-ima*(dwdz1m(i))
        dxdzz(i)=dwdzz1m(i)
        dydzz(i)=-ima*(dwdzz1m(i))
        tin(i) = t(i)
 1400   continue
 

      !calculate the derivatives for the first point
      ! for Tt1(1)=Tin(1)
      ! using formula given in Berrut and Trefethen, SIAM Review 2004

        do i=1,nt2temp
           T2temp(i)=(i-1)*(2*pi/(nt2temp-1))+tin(1) 
        end do
 
        ! Lagrange interpolation for 
        !epsLag=1e-3
       write(*,*), '123',mm,nt2temp
c        nt2temp = 513
       call IntLag(mm,tin, x, y, kLag, nt2temp, T2temp, Rt2temp, 
     1      Zt2temp,  dRt2temp, dZt2temp, epsLag)
       write(*,*), '456, nnt',nt2temp
       write(*,*), '456, nnt',dxdz2(1:nt2temp)
       call IntLag(mm,tin, dxdz, dydz, kLag, nt2temp, T2temp, dxdz2, 
     1      dydz2,  ddxdzdt2, ddydzdt2, epsLag)
c       write(*,*), '789, nnt',nt2temp
c       write(*,*), '456, nnt',dxdzz2(1:nt2)
       call IntLag(mm,tin, dxdzz, dydzz, kLag, nt2temp, T2temp, dxdzz2, 
     1      dydzz2,  ddxdzzdt2, ddydzzdt2, epsLag)
c       write(*,*), '789, nnt',nt2temp

       if (iprintmap.eq.1) then
          iw=16
          itype=3
c          nt2temp = 513
          call quagraph2(iw,t,x,mm,itype,t,y,mm,itype,
     1       'real and imag parts of the inverse map*')
 
c          iw=17
c          itype=3
c          call quagraph2(iw,t2temp(1:nt2temp),Rt2temp(1:nt2temp),nt2temp,
c     1         itype,t2temp(1:nt2temp),Zt2temp(1:nt2temp),nt2temp,itype,
c     2         'real and imag parts of the interpolated inverse map*')
c          iw=18
c          itype=3
          
          call quagraph2(iw,t,dxdz,mm,itype,t,dydz,mm,itype,
     1         'real and imag parts of the inverse map*')
          
          iw=19
          itype=3
          call quagraph2(iw,t2temp,dxdz2,nt2temp,
     1         itype,t2temp,dydz2,nt2temp,itype,
     2         'real and imag parts of the interpolated inverse map*')
          iw=28
          itype=3
          
          call quagraph2(iw,t,dxdzz,mm,itype,t,dydzz,mm,itype,
     1         'real and imag parts of the inverse map*')
          
          iw=29
          itype=3
          call quagraph2(iw,t2temp,dxdzz2,nt2temp,
     1         itype,t2temp,dydzz2,nt2temp,itype,
     2         'real and imag parts of the interpolated inverse map*')
       end if

c        write(*,*), '10, nnt',nt2temp
        do i=1,nt2temp
           dwdz2temp(i) = dxdz2(i)+ ima*dydz2(i)
           dwdzz2temp(i) = dxdzz2(i)+ ima*dydzz2(i)
        end do
c        write(*,*) 'dwdz2'!, dwdz2temp(1:nt2temp)

        write(*,*) 'Start mapping2'
        return
        end
c
c
c
c
c
        subroutine plotgrid(iw,n,z,nr1,ntheta1,nr2,ntheta2,npts,w)
        implicit real *8 (a-h,o-z)
        real *8 x(10 000),y(10 000),fr(10 000),fi(10 000)
        real *8 wk(10 000), wsave(10 000)!, h(10 000)
        complex *16 w(1),z(1),ima,wz,cint
c
c       plots the image of a regular polar grid under the
c       inverse conformal map stored in w
c
c       input:
c            iw - file number
c            n - the length of z
c            z - complex array of the boundary points of the original
c                domain
c            nr1,ntheta1 - for rays, number of radial points 
c                angle points
c            nr2,ntheta2 - for rings, number of radial points
c                angle points
c            npts - number of points to use when evaluating the
c                cauchy integral
c            w - the array created by cnfinv_init
c
c
        ima=(0,1)
        pi=4*atan(1.0d0)
c
c       plot the rays
c
c        nr1 =2
c        h(1)=0.1
c        h(2)= 0.99
        h=1.0d0/(nr1+1)
        h2=2*pi/(ntheta1)

        ntot=0
        do 2000 j=1,ntheta1

        theta=h2*(j-1)

        call fft_wk(npts, w, wk, wsave)
c        write(*,*) 'wk',wk(1:npts)

        do 1800 i=1,nr1
c           wz=h(i)*exp(ima*theta)
           wz=h*i*exp(ima*theta)
           call fft_cauchy(npts,wz,wk,cint)
c          call cnfinv_cauchy(npts,wz,w,cint)
           x(ntot+i)=cint
           y(ntot+i)=-ima*cint
 1800   continue

        ntot=ntot+nr1
 2000   continue
c
c       plot the rings
c
        hr=1.0d0/(nr2+1)
        ht=2*pi/ntheta2

c
        do 3000 j=1,nr2
c
        r=hr*j
c
        do 2800 i=1,ntheta2
           wz=r*exp(ima*ht*i)
           call fft_cauchy(npts,wz,wk,cint)
c           call cnfinv_cauchy(npts,wz,w,cint)
           x(ntot+i)=cint
           y(ntot+i)=-ima*cint
 2800   continue
c
        ntot=ntot+ntheta2
 3000   continue

c
        do 4000 i=1,n
        fr(i)=z(i)
        fi(i)=-ima*z(i)
 4000   continue

        itype2=2
        itype3=3
c        write(*,*) 'iw',iw
c        write(*,*) 'x',x
c        write(*,*) 'y',y
        call quaplot2(iw,fr,fi,n,itype3,x,y,ntot,itype2,
     1       'image of uniform polar grid under inverse map*')

c
        return
        end
c
c
c
c
        subroutine fft_wk(n,w,wk,wsave)
       
        implicit real *8 (a-h,o-z)
        complex *16 w(1),wk(1)
        real *8 wsave(1)
        complex *16 win(n),ima

        ima=(0,1)
c
        do i=1,n
           win(i)=w(i)
        end do
c
        call dcffti(n,wsave)
        call dcfftf(n,win,wsave)
c
c       zero pad and take the inverse fft
c
        do i=1,n
           wk(i)=win(i)/n
        end do

        return
        end
 
        subroutine fft_cauchy1D(n,r,wk,wsave, val)
        implicit real *8 (a-h,o-z)
        dimension wsave(1)
        real *8 r,phi
        complex *16 wk(1),rr
        complex *16 val(1),ima,z
c
c       this subroutine evaluates the integral
c
c           val = \frac{1}{2 \pi \i} \int \frac{w(z)}{z-\zeta} dz
c               = sum_n=0^n=npts (r^n) Z_n e^(ima*n*phi)
c                  where zeta=r*e^(ima*phi)
c
c
        ima=(0,1)
        done=1
        pi=4*atan(done)

        do i=1,n
           val(i)=r**(i-1)*wk(i)
        end do
        call dcffti(n,wsave)
        call dcfftb(n,val,wsave)

        return
        end

        subroutine fft_cauchy(n,zeta,wk, val)
        implicit real *8 (a-h,o-z)

        real *8 r,phi
        complex *16 wk(1),rr
        complex *16 zeta,val,ima,z
c
c       this subroutine evaluates the integral
c
c           val = \frac{1}{2 \pi \i} \int \frac{w(z)}{z-\zeta} dz
c               = sum_n=0^n=npts (r^n) Z_n e^(ima*n*phi)
c                  where zeta=r*e^(ima*phi)
c
c

        ima=(0,1)
c
c       get r and phi
c
        r=abs(zeta)
        rr=log(zeta)
        phi=-ima*rr

        z=0

        do i=1,n
           z=z+r**(i-1)*wk(i)*exp(ima*(i-1)*phi)
        end do
c
        val=z
        return
        end

        subroutine cnfinv_cauchy(npts,w0,z,val)
        implicit real *8 (a-h,o-z)
        complex *16 z(1)
        complex *16 w0,val,ima,w,wwv,w2
c
c       this subroutine evaluates the integral
c
c           val = \frac{1}{2 \pi \i} \int \frac{z(w)}{w-w0} dz
c
c       where the integral is taken over the boundary of the unit disc,
c       and z(w) is the conformal map taking the unit disc to
c       another region, all of which has been precomputed and stored
c       in z from calls to other subroutines
c
c       input:
c
c         npts - the number of points to use in the trapezoidal rule
c         zeta - as in the integral above
c         z - contains info on the conformal map on the boundary
c           as created by a previous call to cnfinv_init
c
c       output:
c
c         val - the value of the cauchy integral above
c
c
        ima=(0,1)
        pi=4*atan(1.0d0)
c
        thresh=1.0d-8
        d=1-abs(zeta)
        if (d .lt. thresh) goto 3000

c
c       use normal trapezoidal quadrature here, zeta far enough
c       from the boundary
c
        w=0
        n=npts
        h=2*pi/n
c
        do 1400 i=1,n
        ttt=h*(i-1)
c        call cnfinv_eval(ier,ttt,w,zzv)
        wwv=z(i) ! assume w is the conformap map for z (nt2=npts)
c        write(*,*) 'zzv(i),i',i,zzv
        w=w+wwv*h*ima*exp(ima*ttt)/(exp(ima*ttt)-w0)
c        z2=zzv*ima*exp(ima*ttt)/(exp(ima*ttt)-zeta)
c        write(*,*) 'z(i),i',i,z
 1400   continue

        w=w/2/pi/ima
        val=w
        return
c

 3000   continue
c
c       the point zeta is very close to the unit disc, use other method
c
        write (*,*) 'zeta too close to boundary, bomb!',a
        stop

c
        return
        end

        subroutine cnf_cauchy(npts,dt,z0,zin,dzdtin,win,val)
        implicit real *8 (a-h,o-z)
        complex *16 win(1),zin(1), dzdtin(1)
        complex *16 z0,val,ima,z,zzv,z2
c
c       this subroutine evaluates the integral
c
c           val = \frac{1}{2 \pi \i} \int \frac{w(z)}{z-z0} dz
c
c       where the integral is taken over the boundary of the unit disc,
c       and w(z) is the conformal map taking the original doamin to the unit disc
c       , all of which has been precomputed and stored
c       in w from calls to other subroutines
c
c       input:
c
c         npts - the number of points to use in the trapezoidal rule
c         zeta - as in the integral above
c         w - contains info on the conformal map on the boundary
c           as created by a previous call to cnfinv_init
c
c       output:
c
c         val - the value of the cauchy integral above
c
c
        ima=(0,1)
        pi=4.0d0*atan(1.0d0)
c
c
c       use normal trapezoidal quadrature here, zeta far enough
c       from the boundary
c
c       val  = \frac{1}{2 \pi \i} sum_{i=1}^{nt} \frac{w_i}{z_i-z0} dzdt_i dt

        z=0
        n=npts
        h=dt
c
c        write(*,*) 'n,h',n,h
        do 1400 i=1,n
c           write(*,*) 'z',z
           z=z+h*win(i)*dzdtin(i)/(zin(i)-z0)
c           write(*,*) 'zin',zin(i),z0
c           write(*,*) 'dzdtin',dzdtin(i),win(i)
 1400   continue

        z=z/2.0d0/pi/ima
        val=z
        return
        end
