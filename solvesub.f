
      subroutine Solovev(x,y,c,d1,d2,d3,
     1        uex,uxex,uyex,uxxex,uxyex, uyyex,f)
      implicit real*8 (a-h,o-z)
c
c      c=1.0d0
c      epsilon=0.5d0!.32d0
c      rkappa=1.0d0 !7d0
c      delta=0.0!.33d0
c      alpha=asin(delta)


      f=c*x**2

     
c      write(*,*) 'd',d1,d2,d3
c      write(*,*) 'rhs1',rhs(1), d1*w(1,1)+d2*w(1,2)+d3*w(1,3)
c      write(*,*) 'rhs2',rhs(2), d1*w(2,1)+d2*w(2,2)+d3*w(2,3)
c      write(*,*) 'rhs3',rhs(3), d1*w(3,1)+d2*w(3,2)+d3*w(3,3)
      uex = c/8.0d0*x**4+d1+d2*x**2+d3*(x**4-4.0d0*x**2*y**2)
      uxex = c/2.0d0*x**3+2.0d0*d2*x+d3*(4.0d0*x**3-8.0d0*x*y**2)
      uyex = -8.0d0*d3*x**2*y
      uxxex = c*3.0d0/2.0d0*x**2+2.0d0*d2 + d3*(12.0d0*x**2-8.0d0*y**2)
      uxyex = -16.0d0*d3*x*y 
      uyyex = -8.0d0*d3*x**2

      return
      end

      subroutine uexact_org(xi,yi,x0,y0,
     1        uex,uxex,uyex,uxxex,uxyex,uyyex,f)
      implicit real*8 (a-h,o-z)
c
      pi = 4.0d0*datan(1.0d0)
      sigma = 0.2
      x0s=x0 !+0.03d0
      x=(xi-x0s)/sigma ; y=(yi-y0)/sigma
      uex = exp(-x**2-y**2)
      uxex = (-2*x)/sigma*exp(-x**2-y**2)
      uxxex = (4*x**2-2)/sigma**2*exp(-x**2-y**2)
      uyex = (-2*y)/sigma*exp(-x**2-y**2)
      uyyex = (4*y**2-2)/sigma**2*exp(-x**2-y**2)
      uxyex = (4*y*x)/sigma**2*exp(-x**2-y**2)
c      if (isolve.eq.0) then
         f = uxxex-1.0d0/xi*uxex+uyyex 

      return
      end

      subroutine uexact(r,theta,
     1        uex,urex,utex,f)
      implicit real*8 (a-h,o-z)
c
      pi = 4.0d0*datan(1.0d0)
      uex = sin(pi*r)*sin(theta)
      urex = pi*cos(pi*r)*sin(theta)
      urrex = -pi*pi*sin(pi*r)*sin(theta)
      utex = sin(pi*r)*cos(theta)
      uttex = -sin(pi*r)*sin(theta)
c      if (isolve.eq.0) then
         f = urrex + urex/r + uttex/(r*r)
c      else if (isolve.eq.1) then
c         f = (urrex + urex/r + uttex/(r*r) - alpha*alpha*uex)
c      else if (isolve.eq.2) then
c         f = (urrex + urex/r + uttex/(r*r) + alpha*alpha*uex)
c      endif
c
      return
      end

      subroutine findinitmax(inext_maxpsi, rndi, tndi
     1     , dpsidri, dpsidti, dpsidrri ,dpsidrti, dpsidtti
     2     , Rt3i, Zt3i, psimi)

      use arrays, only: nt2, nr, rnd, tnd, Rt3, Zt3
      use arrays, only: u, ur, uth, urr, urt, utt, dzdw3, dzdww3

      implicit real*8 (a-h,o-z)
     
      complex *16 ima

      ima = (0.0d0,1.0d0)
      pi = 4.0d0*datan(1.0d0)
      twopi=2.0d0*pi
      teps=0.1   ! criteria of theta to determine up-down symmetry

      if (isymud.eq.1) then     ! for up-down symmetric, tndm=0.0d0 or tndm=pi
         ith=inext_maxpsi/nr+1
         ir=inext_maxpsi-(ith-1)*nr
         rndi=rnd(ir)
         tndi=tnd(ith)
          
         write (*,*) 'original theta for maximum psi=',tndi
         if ((dabs(mod(tndi,twopi)).lt.teps).or.
     1        (dabs(dabs(mod(tndi,twopi))-twopi).lt.teps)) then
            tndi=0.0d0
         else if (dabs(dabs(mod(tndi,twopi))-pi).lt.teps) then
            tndi=pi
         else
            write(*,*) 'maximum psi is not on the midplane'
            write(*,*) 'isymud is enforced to be 0'
            isymud=0
         end if
         write (*,*) 'theta for maximum psi for isym=1 is',tndi
         dzdwr=dzdw3(inext_maxpsi)
         dzdwi=(-ima)*dzdw3(inext_maxpsi)
         dzdwwr=dzdww3(inext_maxpsi)
         dzdwwi=(-ima)*dzdww3(inext_maxpsi)
         
         dRdri=dzdwr*dcos(tndi)-dzdwi*dsin(tndi)
         dRdrri=dzdwwr*dcos(2.0d0*tndi)
     1        -dzdwwi*dsin(2.0d0*tndi)
      
         Rt3i=Rt3(inext_maxpsi)
         Zt3i=Zt3(inext_maxpsi)
         ui=u(inext_maxpsi)
         uri=ur(inext_maxpsi)
         urri=urr(inext_maxpsi)
  
         psimi=ui*dsqrt(Rt3i)
         dpsidri = uri*dsqrt(Rt3i)
     1        +ui/(2.0d0*dsqrt(Rt3i))*dRdri
      
         dpsidrri = urri*dsqrt(Rt3i)
     1        +uri*dRdri/dsqrt(Rt3i)
     1        +ui*dRdrri/(2.0d0*sqrt(Rt3i))

      end if

      if (isymud.ne.1) then     ! for up-down asymmetric
         ith=inext_maxpsi/nr+1
         ir=inext_maxpsi-(ith-1)*nr
         rndi=rnd(ir)
         tndi=tnd(ith)

         dzdwr=dzdw3(inext_maxpsi)
         dzdwi=(-ima)*dzdw3(inext_maxpsi)
         dzdwwr=dzdww3(inext_maxpsi)
         dzdwwi=(-ima)*dzdww3(inext_maxpsi)
         
         dRdri=dzdwr*dcos(tndi)-dzdwi*dsin(tndi)
         dRdrri=dzdwwr*dcos(2.0d0*tndi)
     1        -dzdwwi*dsin(2.0d0*tndi)
         dRdti=dzdwr*rndi*dsin(tndi)+dzdwi*rndi*dcos(tndi)
         dRdrti=dzdwwr*rndi*dsin(2.0d0*tndi)
     1        +dzdwwi*rndi*dcos(2.0d0*tndi)
     2        +dzdwr*(-dsin(tndi))+dzdwi*dcos(tndi)
         dRdtti=dzdwwr*rndi**2*(-dcos(2.0d0*tndi))
     1        +dzdwwi*rndi**2*dsin(2.0d0*tndi)
     2        +dzdwr*rndi*dcos(tndi)-dzdwi*rndi*dsin(tndi)
         
         Rt3i=Rt3(inext_maxpsi)
         Zt3i=Zt3(inext_maxpsi)
         ui=u(inext_maxpsi)
         uri=ur(inext_maxpsi)
         uti=uth(inext_maxpsi)
         urri=urr(inext_maxpsi)
         urti=urt(inext_maxpsi)
         urti=utt(inext_maxpsi)
         
         psimi=ui*dsqrt(Rt3i)
         dpsidri = uri*dsqrt(Rt3i)
     1        +ui/(2.0d0*dsqrt(Rt3i))*dRdri
         dpsidti = uti*dsqrt(Rt3i)
     1        +ui/(2.0d0*dsqrt(Rt3i))*dRdti
         
         dpsidrri = urri*dsqrt(Rt3i)
     1        +uri*dRdri/dsqrt(Rt3i)
     1        +ui*dRdrri/(2.0d0*sqrt(Rt3i))
         dpsidrti=  urti*sqrt(Rt3i)
     1        +uti*(dRdri+dRdti)/(2.0d0*sqrt(Rt3i))
     1        +ui*dRdrti/(2.0d0*sqrt(Rt3i))
         dpsidtti = utti*sqrt(Rt3i)
     1     +uti*dRdti/sqrt(Rt3i)
     1        +ui*dRdtti/(2.0d0*sqrt(Rt3i))

      end if
c      write(*,*) 'u, ur, ut',ui, uri, uti
c      write(*,*) 'rndi, tndi, Rt3i, Zt3i', rndi, tndi, Rt3i, Zt3i

      return
      end

      subroutine findpsim(rndm, tndm, psim)

      use arrays, only: nt2, nr, rnd, kcheb, ucoeff, urcoeff, urrcoeff
      use arrays, only: cftmsub, bnodes, zk
      use arrays, only: isymud, nsub
      
      implicit real*8 (a-h,o-z)

      

      real *8  ucinr(nt2), urcinr(nt2), urrcinr(nt2)
      real *8  ucini(nt2), urcini(nt2), urrcini(nt2)  
      real *8  chcoeff(kcheb)
      real *8  uch(kcheb), urch(kcheb), uthch(kcheb)
      real *8  urrch(kcheb)
      complex *16 ima, wzm, cint
      complex *16 ucin(nt2)
      complex *16 urcin(nt2), urrcin(nt2)

      ima = (0.0d0,1.0d0)
      pi = 4.0d0*datan(1.0d0)
      twopi=2.0d0*pi

         if (dabs(rndm).le.1.0d-8) then
            rndm=1.0d-8
         end if
      
        
         bnodes(1)=0.0d0

         do iisub=2,nsub+1
            if (rndm.lt.bnodes(iisub)) then
               isub=iisub-1
               xch=(rndm-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
               exit
            end if
         end do

         
         do i=1,kcheb
            do j=1,nt2
               inext=((isub-1)*kcheb+i)+(j-1)*nr
               ucin(j)=ucoeff(inext)
               urcin(j)=urcoeff(inext)
               urrcin(j)=urrcoeff(inext)
            end do
            call ftsolnoffth0pi(tndm,nt2,ucin,urcin,urrcin,
     1           uch(i),urch(i),urrch(i))
         end do
         call chftransq(chcoeff,uch,kcheb,cftmsub)
         call chfit(xch, kcheb, chcoeff, um)



         wzm =rndm*exp(ima*tndm)
         call fft_cauchy(nt2,wzm,zk,cint)
         Rt3m=cint
         Zt3m=-ima*cint


      
         psim = um*dsqrt(Rt3m)

         
      return
      end
   



      subroutine findmagaxis(inext_maxpsi, rndi, tndi, dpsidri, dpsidti 
     1     ,dpsidrri ,dpsidrti, dpsidtti, Rt3i, Zt3i, psii
     2     , dpsi2)

      use arrays, only: nt2, nr, rnd, kcheb, ucoeff, urcoeff, urrcoeff
      use arrays, only: cftmsub, bnodes, zk, dzdw2k, dzdwW2k, epsmag
      use arrays, only: isymud, nsub
      
      implicit real*8 (a-h,o-z)

      real *8  ucinr(nt2), urcinr(nt2), urrcinr(nt2)
      real *8  ucini(nt2), urcini(nt2), urrcini(nt2)  
      real *8  chcoeff(kcheb)
      real *8  uch(kcheb), urch(kcheb), uthch(kcheb)
      real *8  urtch(kcheb), uttch(kcheb), urrch(kcheb)
      complex *16 ima, wzm, cint
      complex *16 ucin(nt2)
      complex *16 urcin(nt2), urrcin(nt2)

      ima = (0.0d0,1.0d0)
      pi = 4.0d0*datan(1.0d0)
      twopi=2.0d0*pi

c      eps7=1.0d-18 !14
      dpsi2= dpsidri**2+dpsidti**2 !/rndi**2
      rndm=rndi
      tndm=tndi
      dpsidrm=dpsidri
      dpsidtm=dpsidti
      dpsidrrm=dpsidrri
      dpsidrtm=dpsidrti
      dpsidttm=dpsidtti

      Rt3m=Rt3i
      Zt3m=Zt3i
      psim=psii


      ith=inext_maxpsi/nr+1
      ir=inext_maxpsi-(ith-1)*nr
      a = rnd(ir-1)
      b = rndm
      c = rnd(ir+1)

      hnt=nt2/2
      if (abs(ith-hnd).lt.10) then
         d = a
         a = -c
         b = -b
         c = -d
      end if


      e = 0

      iimag=0
      maxiimag=10
      do while ((iimag.eq.0).or.
     1     (dabs(dpsi2).gt.epsmag).and.(iimag.lt.maxiimag))

         iimag=iimag+1

         if (isymud.eq.1) then  ! for up-down symmetric, tndm=0.0d0 or tndm=pi
                   
            call findpsim(abs(a), tndm, psia) !rnd to psi
            call findpsim(abs(b), tndm, psib)
            call findpsim(abs(c), tndm, psic)

            pbc = psib-psic
            pba = psib-psia
            ba = b-a
            bc = b-c
            d=b-0.5d0*(pbc*ba**2.0d0-pba*bc**2.0d0)/(ba*pbc-bc*pba)

            if (((ba*pbc-bc*pba).ne.0).and.(e.eq.0)) then !2nd interpolating

               drm = abs(d) - rndm
               rndm = abs(d)
               if (d.ge.c) then
                  a = b
                  b = c
                  c = d
               else if (d.lt.a) then
                  c = b
                  b = a
                  a = d
               else if ((d.ge.a).and.(d.lt.b)) then
                  c = b
                  b = d
               else if ((d.ge.b).and.(d.lt.c)) then
                  a = b
                  b = d
               end if
            else if (((ba*pbc-bc*pba).eq.0).or.(e.ne.0)) then !Nelder Mead method
               e=e+1
               if (psia.le.psic) then
                  b = a
                  psib = psia
                  a = c
                  psia = psic
                  c = b
                  psic = psib
               end if

               r = 2.0d0*a - c
               call findpsim(abs(r), tndm, psir)
               if (psir.ge.psia) then !r good > choose good bt s r
                  s = a + 2.0d0*(a-c)
                  call findpsim(abs(s), tndm, psis)
                  if (psis.ge.psir) then
                     c = a
                     psic = psia
                     a = s
                     psia = psis
                  else if (psis.lt.psir) then
                     c = a
                     psic = psia
                     a = r
                     psia = psir
                  end if
                  drm = a-c
               else if (psir.lt.psia) then ! a good > compare r c > choose c close to good r c
                  if (psir.ge.psic) then
                     cc = a + (r-a)/2.0d0
                  else if (psir.lt.psic) then
                     cc = a + (c-a)/2.0d0
                  end if
                  call findpsim(abs(cc), tndm, psicc)
                  if (psicc.ge.psia) then ! c good than a > new a = cc
                     c = a
                     psic = psia
                     a = cc
                     psia = psicc
                     drm = a-c
                  else if (psicc.lt.psia) then ! a good than cc > c = cc
                     c = cc
                     psic = psicc
                     drm = 0.0d0
                  end if
               end if
               rndm = abs(a)
            end if

c            tndm=0.0d0
            if (dabs(rndm).le.1.0d-8) then
               rndm=1.0d-8
               write(*,*) 'rndm',rndm
               exit !return
            end if
         else                   ! for up-down asymmetric
            
            if (rndm.gt.rnd(2)) then !if the magnetic axis is not so close to the unitdisk center 
               drdtjac=dpsidrrm*dpsidttm-dpsidrtm**2
               drm=(-dpsidttm*dpsidrm+dpsidrtm*dpsidtm)/drdtjac
               dtm=(dpsidrtm*dpsidrm-dpsidrrm*dpsidtm)/drdtjac
               
               rndm=rndm+drm
               tndm=tndm+dtm
            else if (rndm.le.1.0d-8) then
               write(*,*) 'rndm',rndm
               exit !return
            else
               drm=-dpsidrm/dpsidrrm 
               rndm=rndm+drm

               if (rndm.le.0.0d0) then
                  rndm=1.0d-8
                  write(*,*) 'rndm',rndm
                  exit !rndm=1.0d-8
               end if
            end if
         end if
         
         write (*,*) 'findmagaxis', iimag, rndm, tndm, drm, dtm

   
c         do i=1, nr-1
c            if (rndm.le.rnd(1)) then
c               isub=1
c               xch=(rndm-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
c               write(*,*) 'very close to origin: smaller than rnd(1)'
c               exit
c            end if
c            if ((rnd(i).le.rndm).and.(rnd(i+1).gt.rndm)) then
c               isub=(i-1)/kcheb+1
c               xch=(rndm-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
c               if (xch.gt.1.0d0) then
c                  isub=isub+1
c                  xch=(rndm-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
c               end if
c               exit
c            end if
c         end do

         bnodes(1)=0.0d0

         do iisub=2,nsub+1
            if (rndm.lt.bnodes(iisub)) then
               isub=iisub-1
               xch=(rndm-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
               exit
            end if
         end do
         write(*,*) 'bnodes',bnodes
         write(*,*) 'isub',isub,xch,bnodes(isub),bnodes(isub+1)

         if (isymud.eq.1) then  ! for up-down symmetric, tndm=0.0d0 or tndm=pi
    
            do i=1,kcheb
               do j=1,nt2
                  inext=((isub-1)*kcheb+i)+(j-1)*nr
                  ucin(j)=ucoeff(inext)
                  urcin(j)=urcoeff(inext)
                  urrcin(j)=urrcoeff(inext)
               end do
               call ftsolnoffth0pi(tndm,nt2,ucin,urcin,urrcin,
     1              uch(i),urch(i),urrch(i))
            end do
            call chftransq(chcoeff,uch,kcheb,cftmsub)
            call chfit(xch, kcheb, chcoeff, um)
            call chftransq(chcoeff,urch,kcheb,cftmsub)
            call chfit(xch, kcheb, chcoeff, urm)
            call chftransq(chcoeff,urrch,kcheb,cftmsub)
            call chfit(xch, kcheb, chcoeff, urrm)

            wzm =rndm*exp(ima*tndm)
            call fft_cauchy(nt2,wzm,zk,cint)
            Rt3m=cint
            Zt3m=-ima*cint
            call fft_cauchy(nt2,wzm,dzdw2k,cint)
            dzdwr=cint
            dzdwi=-ima*cint
            call fft_cauchy(nt2,wzm,dzdww2k,cint)
            dzdwwr=cint
            dzdwwi=-ima*cint
            
            
            dRdrm=dzdwr*dcos(tndm)-dzdwi*dsin(tndm)
            dRdrrm=dzdwwr*dcos(2.0d0*tndm)
     1           -dzdwwi*dsin(2.0d0*tndm)
            
            psim = um*dsqrt(Rt3m)
            
            dpsidrm = urm*dsqrt(Rt3m)
     1           +um/(2.0d0*dsqrt(Rt3m))*dRdrm
            
            dpsidrrm = urrm*dsqrt(Rt3m)
     1           +urm*dRdrm/dsqrt(Rt3m)
     1           +um*dRdrrm/(2.0d0*sqrt(Rt3m))
            
         else  ! for up-down asymmetric

            do i=1,kcheb
               do j=1,nt2
                  inext=((isub-1)*kcheb+i)+(j-1)*nr
                  ucin(j)=ucoeff(inext)
                  urcin(j)=urcoeff(inext)
                  urrcin(j)=urrcoeff(inext)
               end do
               call ftsolnoffth(tndm,nt2,ucin,urcin,urrcin,
     1              uch(i),urch(i),urtch(i),uthch(i),uttch(i),urrch(i))
               
            end do
            call chftransq(chcoeff,uch,kcheb,cftmsub)
            call chfit(xch, kcheb, chcoeff, um)
            call chftransq(chcoeff,urch,kcheb,cftmsub)
            call chfit(xch, kcheb, chcoeff, urm)
            call chftransq(chcoeff,urtch,kcheb,cftmsub)
            call chfit(xch, kcheb, chcoeff, urtm)
            call chftransq(chcoeff,uthch,kcheb,cftmsub)
            call chfit(xch, kcheb, chcoeff, uthm)
            call chftransq(chcoeff,uttch,kcheb,cftmsub)
            call chfit(xch, kcheb, chcoeff, uttm)
            call chftransq(chcoeff,urrch,kcheb,cftmsub)
            call chfit(xch, kcheb, chcoeff, urrm)
            
            
          write(*,*) 'ucin',ucin(1:nt2)
           write(*,*) 'urcin',urcin(1:nt2)
           write(*,*) 'urrcin',urrcin(1:nt2)
           write(*,*) 'xch',xch,um,urm,urtm,uthm,urrm
             write(*,*) 'rndm,tndm',rndm,tndm 
            wzm =rndm*exp(ima*tndm)
            call fft_cauchy(nt2,wzm,zk,cint)
            Rt3m=cint
            Zt3m=-ima*cint
            call fft_cauchy(nt2,wzm,dzdw2k,cint)
            dzdwr=cint
            dzdwi=-ima*cint
            call fft_cauchy(nt2,wzm,dzdww2k,cint)
            dzdwwr=cint
            dzdwwi=-ima*cint
            
            write(*,*) 'Rt3m',Rt3m,zt3m,dzdwr,dzdwi

            dRdrm=dzdwr*dcos(tndm)-dzdwi*dsin(tndm)
            dRdrrm=dzdwwr*dcos(2.0d0*tndm)
     1           -dzdwwi*dsin(2.0d0*tndm)
            dRdtm=dzdwr*rndm*dsin(tndm)+dzdwi*rndm*dcos(tndm)
            dRdrtm=dzdwwr*rndm*dsin(2.0d0*tndm)
     1           +dzdwwi*rndm*dcos(2.0d0*tndm)
     2           +dzdwr*(-dsin(tndm))+dzdwi*dcos(tndm)
            dRdttm=dzdwwr*rndm**2*(-dcos(2.0d0*tndm))
     1           +dzdwwi*rndm**2*dsin(2.0d0*tndm)
     2           +dzdwr*rndm*dcos(tndm)-dzdwi*rndm*dsin(tndm)
            write(*,*) 'dRdrm',drdrm,drdrrm,drdtm,drdrtm
            psim = um*dsqrt(Rt3m)
            
            dpsidrm = urm*dsqrt(Rt3m)
     1           +um/(2.0d0*dsqrt(Rt3m))*dRdrm
            dpsidtm = uthm*dsqrt(Rt3m)
     1           +um/(2.0d0*dsqrt(Rt3m))*dRdtm
            
            dpsidrrm = urrm*dsqrt(Rt3m)
     1           +urm*dRdrm/dsqrt(Rt3m)
     1           +um*dRdrrm/(2.0d0*sqrt(Rt3m))
            dpsidrtm=  urtm*sqrt(Rt3m)
     1           +uthm*(dRdrm+dRdtm)/(2.0d0*sqrt(Rt3m))
     1           +um*dRdrtm/(2.0d0*sqrt(Rt3m))
            dpsidttm = uttm*sqrt(Rt3m)
     1           +uthm*dRdtm/sqrt(Rt3m)
     1           +um*dRdttm/(2.0d0*sqrt(Rt3m))

            write(*,*) 'dpsdrm',dpsitdrm,dpsidtm,dpsidrrm
            dpsi2= dpsidrm**2+dpsidtm**2 !/rndm**2
         end if
      write (*,*) 'findmagaxis2',psim, dpsidrm, dpsidtm

      end do

      rndi=rndm
      tndi=tndm
      dpsidri=dpsidrm
      dpsidti=dpsidtm
      dpsidrri=dpsidrrm
      dpsidrti=dpsidrtm
      dpsidtti=dpsidttm

      Rt3i=Rt3m
      Zt3i=Zt3m
      psii=psim

      return
      end
