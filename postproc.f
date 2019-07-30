!--------------------------------------------------------
! postproc.f: This subroutine is to find contour integral values 
! (e.g. saftey factor, magnetic shear), check the errors of soultion,
! and print the results in some files
!--------------------------------------------------------


!--------------------------------------------------------
!     This routine (postproc2) is called after the iteration converges 
!     to do postprocessing of the solutions
!---------------------------------------------------------     
      subroutine postproc2

      use arrays, only:  Rt1, Zt1, nt1, file_bc, lambdaex, itftype
      use arrays, only:  nsub, kcheb, kLag, nr, ntot, nt2, ksamp2
      use arrays, only:  Rt3, Zt3, rnd, tnd, dzdw3, dzdww3
      use arrays, only:  dwdz3, dwdzz3,  zk, dzdw2k, dzdww2k
      use arrays, only:  g, psii, f, lambda, maxpsi, maxdpsi, psiex      
      use arrays, only:  dpsidr, dpsidrr, u, ur, uth, urr, urt, utt
      use arrays, only:  npsi, psi0, psiB, psiEF, psiiEF, qpsiEF, qpsich
      use arrays, only:  Rmid, epsmaxdist, Rmaxpsi, Zmaxpsi, Rmax, Rmin
      use arrays, only:  csol,reps,rkappa,delta,d1,d2,d3,Rmap0 
      use arrays, only:  nchq, psiichq, fpolsgn, fpolch, shat, cftmq
      use arrays, only:  pprimch, ffprimch, fpolEF, presch, rhoch
      use arrays, only:  spbx, spxb, spdef, qpsich1, jpar0, phiichq
      use arrays, only: q0, F0, p0, pin, pout, FF0, FFin, FFout, phitot 
      use arrays, only: ici, icr, ic_max, iptype, fiter, fout, iprintmap
      use arrays, only: iprintsol,iprintsoldisk, iprintcon, iprintQS
      use arrays, only: iprintEFIT, ifpol, Fedge, fiter22, iogyropsi
      use arrays, only: ifitMil, torcur, iscale, iecom, ijtype
      use arrays, only:  Rt3, Zt3, istability,R0,ibscur, B0vac
      use arrays, only: nchy, psiichy, cftmy, spybx, spyxb,spydef
      use arrays, only: fpass, etai, Zi, nue_star,nui_star, jBSch
      use arrays, only: ne, Te, Vloop, jOhmch, Rmaxis, Zmaxis
      use arrays, only: itrinity, ntpsi, tpsi, nttin, ttin, trho
      use arrays, only: gradpar,Rtrin,Ztrin,btot,bpol,gbdrift,gbdrift0
      use arrays, only: shatlocal1,shatlocal2,shatlocal3
      use arrays, only: cvdrift,cvdrift0, dpsidrho, irho, ptflowch
      use arrays, only: ptfprimch, dpsiidrhoch

      implicit double precision(a-h, o-z)

      parameter(ncon0=200)
      real *8 psi(ntot)
      real *8 psiR(ntot)
      real *8 psiZ(ntot)
      real *8 psiRR(ntot)
      real *8 psiRZ(ntot)
      real *8 psiZZ(ntot)
 
      real *8 uex(ntot)
      real *8 urex(ntot)
      real *8 utex(ntot)

      real *8 shattemp1(ntot)
      real *8 shattemp2(ntot)

c      real *8 psiex(ntot)
      real *8 psiRex(ntot)
      real *8 psiZex(ntot)
      real *8 psiRRex(ntot)
      real *8 psiRZex(ntot)
      real *8 psiZZex(ntot)
      real *8 fsol(ntot) 

      real *8 psicon(ncon0)
      real *8 psicontmp(ncon0)
      real *8 fpolcon(ncon0),FFp(ncon0)
      real *8 qcon(ncon0)
      real *8 qpsiconex(ncon0)
      real *8 rndcon(3*nt2*ncon0)
      real *8 tcon(3*nt2*ncon0)
      real *8 urcon(3*nt2*ncon0)
      real *8 utcon(3*nt2*ncon0)

      real *8 Rcon(3*ksamp2*ncon0)
      real *8 Zcon(3*ksamp2*ncon0)
      real *8 psiRcon(3*ksamp2*ncon0)
      real *8 psiZcon(3*ksamp2*ncon0)
      real *8 tmaxp(3*ksamp2*ncon0)
      real *8 psiconex(3*ksamp2*ncon0)
      real *8 psiRconex(3*ksamp2*ncon0)
      real *8 psiZconex(3*ksamp2*ncon0)
      real *8 Ravg(ncon0), Zavg(ncon0)

      real *8 trnd(3*nt2*ntpsi)  !for trinity option
      real *8 ttnd(3*nt2*ntpsi)
      real *8 tur(3*nt2*ntpsi),tut(3*nt2*ntpsi)
      real *8 turr(3*nt2*ntpsi),turt(3*nt2*ntpsi),tutt(3*nt2*ntpsi)
      real *8 btotch(nttin*ncon0),alphapsi0ch(nttin*ncon0)
      real *8 dbtotdpsi(nttin*ntpsi),dalphadpsi0(nttin*ntpsi)
      real *8 ddalphadpsidt0(nttin*ntpsi),alphatheta0ch(nttin*ncon0)
      real *8 dalphatt(nttin*ntpsi),bpR(nttin*ntpsi)
      real *8 RtrinT(nttin*ntpsi),ZtrinT(nttin*ntpsi)
      real *8 phitrin(nttin*ntpsi) !JPL

      real *8 Rgeom(ntpsi),ageom(ntpsi),gdpsi_av(ntpsi)
      real *8 tzeros1(ntpsi), tzeros2(nttin*ntpsi)

      real *8 tRcon(3*ksamp2*ntpsi)
      real *8 tZcon(3*ksamp2*ntpsi)
      real *8 psiRtcon(3*ksamp2*ntpsi),psiZtcon(3*ksamp2*ntpsi)
      real *8 psiRRtcon(3*ksamp2*ntpsi),psiRZtcon(3*ksamp2*ntpsi)
      real *8 psiZZtcon(3*ksamp2*ntpsi),ttmaxp(3*ksamp2*ntpsi)
      real *8 psitconex(3*ksamp2*ntpsi),psiRtconex(3*ksamp2*ntpsi)
      real *8 psiZtconex(3*ksamp2*ntpsi),psiRRtconex(3*ksamp2*ntpsi)
      real *8 psiRZtconex(3*ksamp2*ntpsi),psiZZtconex(3*ksamp2*ntpsi)
      real *8 shatlocal0(3*ksamp2*ntpsi)

      real *8 tRavg(ntpsi), tZavg(ntpsi), tpsireal(ntpsi)
      real *8 R0mil(ncon0),iaspmil(ncon0)
      real *8 kappamil(ncon0),deltamil(ncon0)
      real *8 Ia(ncon0),Ib(ncon0),Ic(ncon0)
      real *8 Lp(ncon0),Id(ncon0),Ie(ncon0)
      real *8 Ig(ncon0),Ih(ncon0)

 
      real *8 nDi(ncon0),dqdpsi(ncon0)

      real *8 plmvol(ncon0),plmArea(ncon0),acmcur(ncon0),parcur(ncon0)
      real *8 betaP(ncon0),intind(ncon0)
      real *8 bfpol(npsi) ,cfpol(npsi),dfpol(npsi)
      real *8 qoutEF(npsi) ,soutEf(npsi)

      real *8 Rmido(ncon0), Rmidi(ncon0)
      real *8 Btheta2(1:ncon0),betatf(1:ncon0)
      real *8 betar1(ncon0),betar2(ncon0), tempIb(ncon0)
      real *8 tRmido(ntpsi), tRmidi(ntpsi)
      real *8 chcoeff0(ncon0),chcoeff1(ncon0),chcoeff2(ncon0)
      real *8 chcoeff3(ncon0),chcoeff4(ncon0),chcoeff5(ncon0)
      real *8 chcoeff6(ncon0),chcoeff7(ncon0),chcoeff8(ncon0)

      real *8 tfpolcon(ntpsi),tpcon(ntpsi), tdqdpsi(ntpsi)
      real *8 tdfdpsi(ntpsi), tdpdpsi(ntpsi), tdVdpsi(ntpsi)
      real *8 tR0mil(ntpsi),tiaspmil(ntpsi),tkappamil(ntpsi)
      real *8 tdeltamil(ntpsi), tVol(ntpsi), tshatmil(ntpsi)
      real *8 talphamil(ntpsi), tqmil(ntpsi)
      real *8 tdR0mil(ntpsi),tdkappamil(ntpsi)
      real *8 tddeltamil(ntpsi)

      integer *4 icon(ncon0),itcon(ntpsi)
      complex *16 ima
      real *8 pi,mu, intind0
      integer sol_unit, icheckerr, fogyro
      character(15) file_name, sformat

      ima = (0,1)
      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi

      select case (iptype)

      case (0)  !Solovev solutions
         call scalesol(lambda, psiB, psi)
         call calDer(psiR, psiZ, psiRR, psiRZ, psiZZ)
         ncon=nchq
         psicon(1:ncon)=psiichq(1:nchq)/lambda +psiB
   
         fpolcon(1:ncon)= fpolsgn*dabs(F0)
         fpolch(1:ncon)= fpolcon(1:ncon)
         call findcontour(ncon,psicon,nr,nt2,rnd,tnd,psi,dpsidr,dpsidrr
     1        ,rndcon,tcon,urcon,utcon)

         call sortcontour(ncon,psicon,zk,dzdw2k,rndcon
     1        ,tcon,urcon,utcon,Rcon,Zcon,icon,Ravg,Zavg
     2        ,tmaxp,psiRcon,psiZcon,psiconex,psiRconex,psiZconex
     3        ,Rmido, Rmidi)      

    
         write(*,*) 'contour finished'

         call calfint(ncon,psicon,Rcon,Zcon,tmaxp,icon
     1        ,psiRcon,psiZcon,Ia,Ib,Ic,Lp,Id)

         Ia=Ia*dabs(lambda) !rescaling to dimensionless values
         Ib=Ib/dabs(lambda)
         Ic=Ic/dabs(lambda)

         do i=1,ntot
          
            call Solovev(Rt3(i),Zt3(i),csol,d1,d2,d3,
     1           psiex(i),psiRex(i),psiZex(i),psiRRex(i),psiRZex(i), 
     1           psiZZex(i),fsolex)      
         end do

         call findq_solovev(ncon,psicon,fpolcon,qpsiconex)
   
         icheckerr = 1

      case (1,11,2,4)
         ncon=nchq
         psicon(1:ncon)=psiichq(1:nchq)

         call findcontour(ncon,psicon,nr,nt2,rnd,tnd,psii,dpsidr,dpsidrr
     1        ,rndcon,tcon,urcon,utcon)

         call sortcontour(ncon,psicon,zk,dzdw2k,rndcon
     1        ,tcon,urcon,utcon,Rcon,Zcon,icon,Ravg,Zavg
     2        ,tmaxp,psiRcon,psiZcon,psiconex,psiRconex,psiZconex
     3        ,Rmido, Rmidi)
   
         write(*,*) 'contour finished'
         call calfint(ncon,psicon,Rcon,Zcon,tmaxp,icon
     1        ,psiRcon,psiZcon,Ia,Ib,Ic,Lp,Id)

         if (istability.eq.1) then
            call calfint2(ncon,psicon,Rcon,Zcon,tmaxp,icon
     1        ,psiRcon,psiZcon,Ie,Ig,Ih)
c            write(*,*) 'Ie',Ie(1:ncon)
c            write(*,*) 'Ig',Ig(1:ncon)
c            write(*,*) 'Ih',Ih(1:ncon)
         end if
         

         select case (iscale)
            
         case (0) ! p0, FF0 and F0 are used
            call findq0(ncon,lambda,Ic,cftmq,F0,q0)  ! calculate q0 to match F0 
         case (1) ! p0, FF0 and q0 are used
            call findF0(ncon,lambda,Ic,cftmq,q0,F0) ! calculate F0 to match q0 
         case (2) ! p0/FF0, torcur and F0 are used (default)
            call matchtorcur(ncon,lambda,Ia,cftmq,torcur,factorl)
            lambda=lambda*factorl
            FF0=FF0/factorl
            p0=p0/factorl
            pprimch=pprimch/factorl
            ffprimch=ffprimch/factorl
            call findq0(ncon,lambda,Ic,cftmq,F0,q0)
            write(fiter,*) 'After normalization by toroidal current:'
            write(fiter,*) 'lambda,factor,psi0', lambda,factorl,
     1        (1.0d0/lambda)+psiB
         case (3) ! p0/FF0, torcur and q0 are used 
            call matchtorcur(ncon,lambda,Ia,cftmq,torcur,factorl)
            lambda=lambda*factorl
            FF0=FF0/factorl
            p0=p0/factorl
            pprimch=pprimch/factorl
            ffprimch=ffprimch/factorl
            call findF0(ncon,lambda,Ic,cftmq,q0,F0)
            write(fiter,*) 'After normalization by toroidal current:'
            write(fiter,*) 'lambda,factor,psi0', lambda,factorl,
     1        (1.0d0/lambda)+psiB
         case (4) ! p0/FF0, F0 and q0 are used 
            call normbyF0q0(ncon,lambda,Ic,cftmq,q0,F0,factorl)
            lambda=lambda*factorl
            FF0=FF0/factorl
            p0=p0/factorl
            pprimch=pprimch/factorl
            ffprimch=ffprimch/factorl
            write(fiter,*) 'After normalization to match F0 and q0:'
            write(fiter,*) 'lambda,factor,psi0', lambda,factorl,
     1        (1.0d0/lambda)+psiB
         case(5,7) ! dp/dpsii, FF0, torcur and F0 are used where psii=0 at the core and psii=1 at the edge
            call matchtorcur(ncon,lambda,Ia,cftmq,torcur,factorl)
            lambda=lambda*factorl
            FF0=FF0/factorl
            jpar0=jpar0/factorl
c            pprimch=pprimch/factorl
            ffprimch=ffprimch/factorl
            if (ifpol.eq.1) then
               call findF0byFedge(ncon,lambda,psiichq,ffprimch,
     1           spdef,Fedge,F0)
            end if
            call findq0(ncon,lambda,Ic,cftmq,F0,q0)
            write(fiter,*) 'After normalization by toroidal current:'
            write(fiter,*) 'lambda,factor,psi0', lambda,factorl,
     1        (1.0d0/lambda)+psiB
         case(6,8) ! dp/dpsii, FF0, torcur and q0 are used where psii=0 at the core and psii=1 at the edge
            call matchtorcur(ncon,lambda,Ia,cftmq,torcur,factorl)
            lambda=lambda*factorl
            FF0=FF0/factorl
            jpar0=jpar0/factorl
c            pprimch=pprimch/factorl
            ffprimch=ffprimch/factorl
            call findF0(ncon,lambda,Ic,cftmq,q0,F0)
            write(fiter,*) 'After normalization by toroidal current:'
            write(fiter,*) 'lambda,factor,psi0', lambda,factorl,
     1        (1.0d0/lambda)+psiB
         end select

         write(*,*) 'p0,FF0',p0,FF0
         write(*,*) 'q0,F0',q0,F0
        

         !call scalesol(lambda, psiB, psi)
         !call calDer(psiR, psiZ, psiRR, psiRZ, psiZZ)

         
c         FFp(1:ncon)=FF0*(1-(1-psiichq(1:ncon))**FFin)**FFout
         if (iecom.ne.2) then !F in ECOMQ is determined
            if (ifpol.eq.0) then
               call findFbyFFp0(ncon,lambda,psiichq,ffprimch,
     1              spbx,F0,fpolch)
            else if (ifpol.eq.1) then
               call findFbyFFp1(ncon,lambda,psiichq,ffprimch,
     1              spxb,Fedge,fpolch)
            end if
         end if
         psicon(1:ncon)=psiichq(1:nchq)/lambda +psiB
         fpolcon(1:ncon) = fpolch(1:nchq) 
         
         icheckerr = 0

      case (3) !EFIT
         
         ncon=nchq
         psicon(1:ncon)=psiichq(1:nchq)

         call spline(-psiiEF, fpolEF     
     1           , bfpol, cfpol, dfpol, npsi)
         do i= 1, ncon
            call ispline(-psicon(i),-psiiEF
     1              , fpolEF, bfpol ,cfpol,dfpol,npsi,fpolch(i))
         end do
         call ispline(-1.0d0,-psiiEF
     1              , fpolEF, bfpol ,cfpol,dfpol,npsi,F0)
c         write(*,*) 'fpolcon',fpolcon(1:ncon)
c         ncon=npsi-2
c         psicon(1:ncon)=psiEF(3:npsi)
c         qpsiconex(1:ncon)=qpsiEF(3:npsi)
      
         call findcontour(ncon,psicon,nr,nt2,rnd,tnd,psii,dpsidr,dpsidrr
     1        ,rndcon,tcon,urcon,utcon)
    
         call sortcontour(ncon,psicon,zk,dzdw2k,rndcon
     1        ,tcon,urcon,utcon,Rcon,Zcon,icon,Ravg,Zavg
     2        ,tmaxp,psiRcon,psiZcon,psiconex,psiRconex,psiZconex
     3        ,Rmido, Rmidi)

         write(*,*) 'contour finished'
         call calfint(ncon,psicon,Rcon,Zcon,tmaxp,icon
     1        ,psiRcon,psiZcon,Ia,Ib,Ic,Lp,Id)
c         write(*,*) 'Ic0',lambda,Ic(1:ncon)
!         call scalesol(lambda, psiB, psi)
!         call calDer(psiR, psiZ, psiRR, psiRZ, psiZZ)
 
         if (ifpol.eq.0) then
            call findFbyFFp0(ncon,lambda,psiichq,ffprimch,
     1           spbx,F0,fpolch)
         else if (ifpol.eq.1) then
            call findFbyFFp1(ncon,lambda,psiichq,ffprimch,
     1           spxb,Fedge,fpolch)
         end if
         fpolcon(1:ncon) = fpolch(1:nchq) 

c         if (ibtype.eq.3) then
c            icheckerr = 1
c         else
            icheckerr = 0
c         end if
         
      end select

      qcon(1:ncon)=Ic(1:ncon)*dabs(lambda/(2.0d0*pi)
     1        *fpolcon(1:ncon)) !*fpolsgn
c      write(*,*) 'qcon',qcon(1:ncon)
  
      call findtorflux(ncon,lambda,psiichq,qcon,
     1     spbx,spdef, phitot,phiichq)

      call findrhoch(ncon,Rmido, Rmidi)

      if (itrinity.eq.1) then 
         call findtpsifromtrho(ncon)
         call findcontour2der(ntpsi,tpsi,nr,nt2,rnd,tnd,psii
     1        ,dpsidr,dpsidrr,trnd,ttnd,tur,tut
     2        ,turr,turt,tutt)
         
         if (iptype.eq.0) then  !Solovev solutions
            tpsireal(1:ntpsi)=tpsi(1:ntpsi)/lambda +psiB
         else
            tpsireal(1:ntpsi)=tpsi(1:ntpsi)
         end if
         
         call sortcontour2der (ntpsi,tpsireal,zk,dzdw2k,dzdww2k
     1        ,trnd,ttnd,tur,tut,turr,turt,tutt ,tRcon,tZcon,itcon
     2        ,tRavg,tZavg,ttmaxp,psiRtcon,psiZtcon,psiRRtcon
     3        ,psiRZtcon,psiZZtcon,psitconex,psiRtconex,psiZtconex
     4        ,psiRRtconex,psiRZtconex,psiZZtconex,tRmido,tRmidi)
         
         if (iptype.eq.0) then 
            itotc = 0
            errmax =0.0d0;errpsiRmax=0.0d0;errpsiZmax=0.0d0
            errpsiRRmax=0.0d0;errpsiRZmax=0.0d0;errpsiZZmax=0.0d0
            
            do k=1, ntpsi
               istart=itcon(k)
               iend=itcon(k+1)-1
               nin=iend-istart+1
               do j=1, nin
                  itotc=itotc+1
                  errtmp = abs(psitconex(itotc)-tpsireal(k))
                  errpsiRtmp = abs(psiRtconex(itotc)-psiRtcon(itotc))
                  errpsiZtmp = abs(psiZtconex(itotc)-psiZtcon(itotc))
                  errpsiRRtmp = abs(psiRRtconex(itotc)-psiRRtcon(itotc))
                  errpsiRZtmp = abs(psiRZtconex(itotc)-psiRZtcon(itotc))
                  errpsiZZtmp = abs(psiZZtconex(itotc)-psiZZtcon(itotc))
                  
                  if (errmax.lt.errtmp) then
                     errmax= errtmp
                  end if
                  if (errpsiRmax.lt.errpsiRtmp) then
                     errpsiRmax= errpsiRtmp
                  end if
                  if (errpsiZmax.lt.errpsiZtmp) then
                     errpsiZmax= errpsiZtmp
                  end if
                  if (errpsiRRmax.lt.errpsiRRtmp) then
                     errpsiRRmax= errpsiRRtmp
                  end if
                  
                  if (errpsiRZmax .lt. errpsiRZtmp) then
                     errpsiRZmax = errpsiRZtmp
                  end if
                  
                  if (errpsiZZmax.lt.errpsiZZtmp) then
                     errpsiZZmax= errpsiZZtmp
                  end if
               end do
            end do
            
            write(*,*) 'error in contours',errmax,errpsiRmax
     1              ,errpsiZmax,errpsiRRmax,errpsiRZmax,errpsiZZmax
         end if                 !iptype=0
         
      end if                    !itrinity=1


      call scalesol(lambda, psiB, psi)
      call calDer(psiR, psiZ, psiRR, psiRZ, psiZZ)

      if (icheckerr.eq.1) then
         if ((iptype.eq.3) .and.(ibtype.eq.3)) then
            call findqprofile(ncon,lambda,psiichq,qcon,cftmq,
     1        npsi,psiiEF,qoutEF,soutEF)
            call checkerr(psiRex, psiZex, psiRRex, psiRZex, psiZZex,
     1     psi, psiR, psiZ, psiRR, psiRZ, psiZZ, npsi,psicon, psiRcon, 
     2      psiZcon, psiconex, psiRconex, psiZconex, qoutEF, 
     3     qpsiEF, icon)
         else 
            call checkerr(psiRex, psiZex, psiRRex, psiRZex, psiZZex,
     1     psi, psiR, psiZ, psiRR, psiRZ, psiZZ, ncon,psicon, psiRcon, 
     2      psiZcon, psiconex, psiRconex, psiZconex, qcon, 
     3     qpsiconex, icon)
         end if
         write(*,*) 'finished checking errors'

         !after error checking, rescaling to normalized psi for iptype=0 
         psiRcon=psiRcon*lambda  
         psiZcon=psiZcon*lambda
         psiRtcon=psiRtcon*lambda
         psiZtcon=psiZtcon*lambda
         psiRRtcon=psiRRtcon*lambda
         psiRZtcon=psiRZtcon*lambda
         psiZZtcon=psiZZtcon*lambda

      end if

      do i=1,nchq
           presch(i)=0.0d0
           do k=1,nchq
              presch(i)=presch(i)+spxb((k-1)*nchq+i)*pprimch(k)
           end do
      end do
      presch(1:nchq)=presch(1:nchq)/(2.0d0*(lambda)) !final scaling
      if (itftype.ne.0) then
          ptflowch(1:nchq)=ptflowch(1:nchq)/(-lambda) !final scaling
      end if 

      if (iecom.eq.2) then   !ECOMQ
c         if(iqtype.ne.1) then
c            if (ifpol.eq.0) then
c               call findFbyFFp0(ncon,lambda,psiichq,ffprimch,
c     1              spbx,F0,fpolch)
c            else if (ifpol.eq.1) then
c               call findFbyFFp1(ncon,lambda,psiichq,ffprimch,
c     1              spxb,Fedge,fpolch)
c            end if
c         end if
c         write(*,*) 'qpsich input',qpsich(1:ncon)
         call findq(ncon,psicon,Rcon,Zcon,tmaxp,icon,fpolch
     1        ,psiRcon,psiZcon, qpsich1)
         qpsich1=qpsich1*dabs(lambda)
c     qpsich=-qpsich
      write(*,*) 'qpsich calculated',qpsich1(1:ncon)
         open(2,file='ECOMQ_profiles.out')
         write(2,*) 'normalized psi'
         write(2,*) psiichq(1:nchq)
         write(2,*) 'r/a'
         write(2,*) rhoch(1:nchq)
         write(2,*) 'FdF/dpsi'
         write(2,*) ffprimch(1:nchq)  !exact ffprimch
         write(2,*) 'F'
         write(2,*) fpolch(1:nchq)
         write(2,*) 'dp/dpsi'
         write(2,*) pprimch(1:nchq)
         write(2,*) 'p [Pa]'
         write(2,*) presch(1:nchq)
         write(2,*)  'pressure due to flow :nmR^2omega^2/2 [Pa]'
         write(2,*) ptflowch(1:nchq)
         write(2,*) 'flow Mach number : sqrt(2*ptf/p)'
         write(2,*) dsqrt(2*ptflowch(1:nchq)/presch(1:nchq))
         write(2,*) 'q profile'
         write(2,*) qpsich1(1:nchq)
         close(2)
      else if (iecom.eq.1) then !ECOMJ
         if(ijtype.ne.1) then
           
            if (ifpol.eq.0) then
               call findFbyFFp0(ncon,lambda,psiichq,ffprimch,
     1              spbx,F0,fpolch)
            else if (ifpol.eq.1) then
               call findFbyFFp1(ncon,lambda,psiichq,ffprimch,
     1              spxb,Fedge,fpolch)
            end if
         end if
         call findq(ncon,psicon,Rcon,Zcon,tmaxp,icon,fpolch
     1        ,psiRcon,psiZcon,qpsich1)
         qpsich1=qpsich1*dabs(lambda)
c     qpsich=-qpsich
      write(*,*) 'qpsich calculated',qpsich1(1:ncon)
         open(2,file='ECOMJ_profiles.out')
         write(2,*) 'normalized psi'
         write(2,*) psiichq(1:nchq)
         write(2,*) 'r/a'
         write(2,*) rhoch(1:nchq)
         write(2,*) 'FdF/dpsi'
         write(2,*) ffprimch(1:nchq)  !exact ffprimch
         write(2,*) 'F'
         write(2,*) fpolch(1:nchq)
         write(2,*) 'dp/dpsi'
         write(2,*) pprimch(1:nchq)
         write(2,*) 'p [Pa]'
         write(2,*) presch(1:nchq)
         write(2,*)  'pressure due to flow :nmR^2omega^2/2 [Pa]'
         write(2,*) ptflowch(1:nchq)
         write(2,*) 'flow Mach number : sqrt(2*ptf/p)'
         write(2,*) dsqrt(2*ptflowch(1:nchq)/presch(1:nchq))
         write(2,*) 'q profile'
         write(2,*) qpsich1(1:nchq)
    
         close(2)
      else    !ECOM

        open(2,file='ECOM_profiles.out')
         write(2,*) 'normalized psi'
         write(2,*) 1.0d0-psiichq(1:nchq)
         write(2,*) 'r/a'
         write(2,*) rhoch(1:nchq)
         write(2,*) 'FdF/dPsi'
         write(2,*) ffprimch(1:nchq) 
         write(2,*) 'F'
         write(2,*) fpolch(1:nchq)
         write(2,*) 'dp/dPsi'
         write(2,*) pprimch(1:nchq)
         write(2,*) 'p [Pa]'
         write(2,*) presch(1:nchq)
        
         if (itftype.ne.0) then
            rmin=R0*reps
            write(*,*) 'rmin',rmin,lambda
            tempIb(1:ncon)=0.0d0
            do i=1,ncon
               do k=1,ncon
               tempIb(i)=tempIb(i)+spbx((k-1)*ncon+i)*Ib(k)
               end do
            end do
         
            do i=1,ncon
               temp1=0.0d0
               temp2=0.0d0
               temp1d=0.0d0
               do k=1,ncon
               temp1=temp1+spbx((k-1)*ncon+i)*ptfprimch(k)*tempIb(k)
               temp2=temp2+spbx((k-1)*ncon+i)*pprimch(k)*tempIb(k)
               end do
c               write(*,*) 'temp1,temp2',temp1,temp2
               temp1=temp1/lambda/2.0d0
               temp2=temp2/lambda/2.0d0*mu
               Btheta2(i)=dpsiidrhoch(i)**2/rmin**2/(lambda**2)/R0**2
               temp3=tempIb(i)*Btheta2(i)/2.0d0
               betatf(i)=temp1/temp3


               betaP(i)=temp2/temp3
            end do
            
            write(2,*) 'Poloidal beta'
            write(2,*) betaP(1:ncon)
            write(2,*) 'beta toroidal flow 1'
            write(2,*) betatf(1:ncon)

            do i=1,ncon
               temp1=0.0d0
               temp2=0.0d0
               do k=1,ncon
               temp1=temp1+spbx((k-1)*ncon+i)*ptfprimch(k)*rhoch(k)**2
               temp2=temp2+spbx((k-1)*ncon+i)*pprimch(k)*rhoch(k)**2
               end do
c               write(*,*) 'temp1,temp2',temp1,temp2
               temp1=temp1*rmin**2/lambda/2.0d0
               temp2=temp2*rmin**2/lambda/2.0d0*mu
               Btheta2(i)=dpsiidrhoch(i)**2/rmin**2/(lambda**2)/R0**2
               temp3=rhoch(i)**2*rmin**2*Btheta2(i)/2.0d0
               betatf(i)=temp1/temp3


               betaP(i)=temp2/temp3
            end do
       
            do i=1,ncon
               temp1=0.0d0
               temp2=0.0d0
               do k=1,ncon
               drdpsii=1.0d0/(dpsiidrhoch(k))*rmin
               temp1=temp1+spxb((k-1)*ncon+i)
     1             *betatf(k)*rhoch(k)*rmin*drdpsii
               temp2=temp2+spbx((k-1)*ncon+i)*Btheta2(k)
     1             *betatf(k)*rhoch(k)**3*rmin**3*drdpsii
               end do
               betar1(i)=temp1/2.0d0
               betar2(i)=temp2/(Btheta2(i)*rhoch(i)**2*rmin**2)/2.0d0
            end do
            write(2,*) 'dp_omega/dPsi'
            write(2,*) ptfprimch(1:nchq)
            write(2,*)  'pressure due to flow :nmR^2omega^2/2 [Pa]'
            write(2,*) ptflowch(1:nchq)
            write(2,*) 'flow Mach number : sqrt(2*ptf/p)'
            write(2,*) dsqrt(2*ptflowch(1:nchq)/presch(1:nchq))
            write(2,*) 'Btheta^2 [T^2]'
            write(2,*) Btheta2(1:ncon) 
            write(2,*) 'Poloidal beta'
            write(2,*) betaP(1:ncon)
            write(2,*) 'beta toroidal flow'
            write(2,*) betatf(1:ncon)
            write(2,*) 'betar1'
            write(2,*) betar1(1:ncon)
            write(2,*) 'betar2'
            write(2,*) betar2(1:ncon)
         end if

         close(2)
      end if

      call findshat(ncon,lambda,psiichq
     1    ,qcon,cftmq,shat)
   
c      write(*,*) 'Ia',Ia(1:ncon)
c      write(*,*) 'Ib',Ib(1:ncon)
c      write(*,*) 'Ic',Ic(1:ncon)
c      write(*,*) 'Id',Id(1:ncon)
      call findPhyVal(ncon,lambda,psiichq,Ia,Ib,Ic,Id,Lp
     1    ,pprimch,fpolch,cftmq,spdef,spbx,spxb
     2     ,plmvol,plmArea,acmcur,presch,betaP,intind,B0vac
     3     ,presavg,Wpol,betaTot,plmvol0,plmArea0,acmcur0,betaP0
     4     ,intind0,edgecur,parcur)
 
   

      if (istability.eq.1) then 
         call Mercier(ncon,lambda,psiichq,Ia,Ib,Ic,Ie,Ig,Ih
     1    ,pprimch,fpolch,qcon,cftmq,dqdpsi,nDi)
         call Troyon(acmcur0,B0vac,reps,R0,rkappa,plmArea0
     1    ,betaTroy1,betaTroy2)
      end if

      if (ibscur.eq.1) then 
         call chsetupq(nchy, psiichy, cftmy, spybx, spyxb, spydef) 
         call findfpass(ncon,lambda,psiichq,Rcon,Zcon,tmaxp,icon
     1    ,psiRcon,psiZcon
     2        ,fpolch,nchy,psiichy,spydef,fpass)
       
        call findBScur(ncon,lambda,Ib,Ic,pprimch
     1    ,fpolch,fpass,etai,Zi,nue_star,nui_star,jBSch)
        call findOhmcur(ncon,lambda,ne,Te,fpass,Zi
     1           ,nue_star,Vloop,jOhmch)
      end if

      if (ifitMil.eq.1) then
         epsFIT = 1.0d-13
         call fitMiller(ncon,icon,Rcon,Zcon,Ravg,Zavg
     1        ,epsFIT,R0mil,iaspmil,kappamil,deltamil)

         call printMiller(ncon,psicon,R0mil,iaspmil
     1        ,kappamil,deltamil, Ib)
      end if
      
      if (itrinity.eq.1) then 
        
c         call sortcontour(ntpsi,tpsi,zk,dzdw2k,trnd
c     1        ,ttnd,tur,tut,tRcon,tZcon,itcon,tRavg,tZavg
c     2        ,ttmaxp,psiRtcon,psiZtcon,psitconex,psiRtconex,psiZtconex
c     3        ,tRmido, tRmidi)
         write(*,*) 'contour finished for trinity'
c         write(*,*) 'psiRZtcon',psiRztcon(1:ntpsi*ksamp2)
c         write(*,*) 'psiRZtcon',ntpsi,ncon

         call calmetric0(ncon,psicon,Rcon,Zcon,tmaxp,icon
     1        ,psiRcon,psiZcon, fpolch
     2     ,nttin,ttin, btotch,alphapsi0ch,alphatheta0ch)  ! chebyshev psii grid

c         write(*,*) 'btotch',btotch(1:ncon*nttin)
c         write(*,*) 'alphapsi0ch',alphapsi0ch(1:ncon*nttin)
c         write(*,*) 'alphatheta0ch',alphatheta0ch(1:ncon*nttin)

         call chftransq(chcoeff0, rhoch, ncon, cftmq)
         call chftransq(chcoeff1, fpolch, ncon, cftmq)
         call chftransq(chcoeff2, presch, ncon, cftmq) 
         call chftransq(chcoeff3, plmvol, ncon, cftmq) 
         call chftransq(chcoeff4, qcon, ncon, cftmq)
         call chftransq(chcoeff5, R0mil, ncon, cftmq)
         call chftransq(chcoeff6, iaspmil, ncon, cftmq)
         call chftransq(chcoeff7, kappamil, ncon, cftmq)
         call chftransq(chcoeff8, deltamil, ncon, cftmq)
 
c         write(*,*) 'fpolch',fpolch(1:ncon),presch(1:ncon)
c         write(*,*) 'fpolch2',chcoeff1(1:ncon),chcoeff2(1:ncon)
     
         do i=1,ntpsi
c            call chfit(1.0d0-tpsi(i),ncon,chcoeff0,trho(i))
            call chfit(1.0d0-tpsi(i),ncon,chcoeff3,tvol(i))
            if (tvol(i).lt.0) tvol(i)=0.0d0
            call chfit(1.0d0-tpsi(i),ncon,chcoeff4,tqmil(i))
            call chfit(1.0d0-tpsi(i),ncon,chcoeff5,tR0mil(i))
            call chfit(1.0d0-tpsi(i),ncon,chcoeff6,tiaspmil(i))
            call chfit(1.0d0-tpsi(i),ncon,chcoeff7,tkappamil(i))
            call chfit(1.0d0-tpsi(i),ncon,chcoeff8,tdeltamil(i))

            call chderiv(1.0d0-tpsi(i),ncon,chcoeff0,drhodpsii)
            dpsidrho(i)=1.0d0/drhodpsii/lambda !used non-normalized psi

            call chfit(1.0d0-tpsi(i),ncon,chcoeff1,tfpolcon(i))
            call chderiv(1.0d0-tpsi(i),ncon,chcoeff1,tdfdpsi(i))
            tdfdpsi(i)=tdfdpsi(i)*lambda

            call chfit(1.0d0-tpsi(i),ncon,chcoeff2,tpcon(i))
            call chderiv(1.0d0-tpsi(i),ncon,chcoeff2,tdpdpsi(i))
            tdpdpsi(i)=tdpdpsi(i)*lambda

            call chderiv(1.0d0-tpsi(i),ncon,chcoeff3,tdvdpsi(i))
            tdvdpsi(i)=tdvdpsi(i)*lambda
c            write(*,*) 'tdvdps',tdvdpsi(i),1.0d0-tpsi(i)
            call chderiv(1.0d0-tpsi(i),ncon,chcoeff4,tdqdpsi(i))
            tdqdpsi(i)=tdqdpsi(i)*lambda

            tshatmil(i)=tdqdpsi(i)*2.0d0*tvol(i)/tdvdpsi(i)


            talphamil(i)=-tdpdpsi(i)*2.0d0*tdvdpsi(i)/(2.0d0*pi)**2
     1          *mu*dsqrt(tvol(i)/(2.0d0*pi**2*tR0mil(i)))
            
            call chderiv(1.0d0-tpsi(i),ncon,chcoeff5,tdR0mil(i))
            tdR0mil(i)=tdR0mil(i)/(drhodpsii*reps*R0)

            call chderiv(1.0d0-tpsi(i),ncon,chcoeff6,diaspdpsii)

            call chderiv(1.0d0-tpsi(i),ncon,chcoeff7,tdkappamil(i))
            tdkappamil(i)=tdkappamil(i)/diaspdpsii*tiaspmil(i)
     1           /tkappamil(i)


            call chderiv(1.0d0-tpsi(i),ncon,chcoeff8,tddeltamil(i))
            tddeltamil(i)=tddeltamil(i)/diaspdpsii*tiaspmil(i)
     1           /dsqrt(1-tdeltamil(i)**2)
         end do
c         write(*,*) 'dpsidrho1',dpsiidrhoch(1:ncon)
c         write(*,*) 'dpsidrho2',dpsidrho(1:ntpsi)*lambda
c         write(*,*) 'tdfdpsi',tdfdpsi(1:ntpsi),tfpolcon(1:ntpsi)
c         write(*,*) 'tdpdpsi',tdpdpsi(1:ntpsi),tpcon(1:ntpsi)
c         write(*,*) 'plmvol',plmvol(1:ncon),chcoeff3(1:ncon)
         

         do j=1,nttin
            call chftransq(chcoeff3, btotch((j-1)*ncon+1), ncon, cftmq)
            call chftransq(chcoeff4, alphapsi0ch((j-1)*ncon+1)
     1           , ncon, cftmq)
            call chftransq(chcoeff5, alphatheta0ch((j-1)*ncon+1)
     1           , ncon, cftmq)
            do i=1,ntpsi
               call chderiv(1.0d0-tpsi(i),ncon,chcoeff3,tmp3)
               dbtotdpsi((i-1)*nttin+j)=tmp3*lambda
               call chderiv(1.0d0-tpsi(i),ncon,chcoeff4,tmp4)
               dalphadpsi0((i-1)*nttin+j)=tmp4*lambda
               call chderiv(1.0d0-tpsi(i),ncon,chcoeff5,tmp4)
               ddalphadpsidt0((i-1)*nttin+j)=tmp4*lambda
            end do
         end do
c         write(*,*) 'tfpolcon',tfpolcon(1:ntpsi),tpcon(1:ntpsi)

c         write(*,*) 'dbtotdpsi',dbtotdpsi(1:ntpsi*nttin)
c         write(*,*) 'dalphapsi0ch',dalphadpsi0(1:ntpsi*nttin)
c         write(*,*) 'dalphapsidt0',ddalphadpsidt0(1:ntpsi*nttin)

         write(*,*) 'calculate local mertics for trinity'
         call calmetric(ntpsi,tpsi,tRcon,tZcon,ttmaxp,itcon,psiRtcon
     1        ,psiZtcon,tfpolcon,tdfdpsi,tdpdpsi,dbtotdpsi,dalphadpsi0
     2        ,nttin,ttin,gradpar,Rtrin,Ztrin,bpR,btot,bpol,gbdrift
     3        ,gbdrift0,cvdrift,cvdrift0,shatlocal1,dalphatt,phitrin)



         gdpsi_av=0.0d0
         do i=1,ntpsi
            do j=1,nttin
               gdpsi_av(i)=gdpsi_av(i)+bpR((i-1)*nttin+j)
         !phitrin((i-1)*nttin+j)=dalphadpsi0((i-1)*nttin+j)*tfpolcon(i) !JPL 
            end do
            gdpsi_av(i)=gdpsi_av(i)/nttin
         end do

c         shatlocal2(1:ntpsi*nttin)=ddalphadpsidt0(1:ntpsi*nttin)
c     1       *(gradpar(1:ntpsi*nttin)*btot(1:ntpsi*nttin)) 
         shatlocal2(1:ntpsi*nttin)=
     1        (ddalphadpsidt0(1:ntpsi*nttin)+dalphatt(1:ntpsi*nttin))
     1       *(gradpar(1:ntpsi*nttin)*btot(1:ntpsi*nttin)) 
c         write(*,*) 'shatlocal201',ddalphadpsidt0(1:nttin)
c         write(*,*) 'shatlocal202',    dalphatt(1:nttin)
c         write(*,*) 's201',ddalphadpsidt0(ntpsi*nttin-nttin:ntpsi*nttin)
c         write(*,*) 'shat202',dalphatt(ntpsi*nttin-nttin:ntpsi*nttin)
c         write(*,*) 'shatlocal2',shatlocal2(1:ntpsi*ntpsi)

         call findlocalshear(ntpsi,tpsi,itcon,ttmaxp,tRcon,tZcon
     1        ,psiRtcon,psiZtcon,psiRRtcon,psiRZtcon,psiZZtcon
     2        ,tfpolcon,tdfdpsi,shatlocal0,nttin,ttin,shatlocal3)

c         write(*,*) 'shatlocal0',shatlocal0(1:ksamp2*ntpsi)
         ! calculate local shear except fpol and df/dpsi part in 2d domain
         call chftransq(chcoeff0, ffprimch, nchq, cftmq)
         call chftransq(chcoeff1, fpolch, nchq, cftmq)
         do i=1,ntot
            call chfit(1.0d0-psii(i), nchq, chcoeff0, ffprime3)
            call chfit(1.0d0-psii(i), nchq, chcoeff1, fpol3)
            dpsi2=psiR(i)**2+psiZ(i)**2
            dpsi2part1=(psiZ(i)**2-psiR(i)**2)*(psiRR(i)-psiZZ(i))
     1           /dpsi2**2
            Reffect1=-psiR(i)/Rt3(i)/dpsi2
            dpsidRdZpart1=-4.0d0*psiR(i)*psiZ(i)*psiRZ(i)/dpsi2**2
            shattemp1(i)=(dpsi2part1+Reffect1+dpsidRdZpart1)
     1           *fpol3/Rt3(i)**2+ffprime3/fpol3/Rt3(i)**2

            dR=Rt3(i)-Rmaxis
            dZ=Zt3(i)-Zmaxis
            dnum1=(dR*psiR(i)+dZ*psiZ(i))
            dpsi2part2=(-dR**2*psiRR(i)-dZ**2*psiZZ(i))/dnum1**2
            Reffect2=-dR/Rt3(i)/dnum1
            dpsidRdZpart2=-2.0d0*dR*dz*psiRZ(i)/dnum1**2
            shattemp2(i)=dpsi2part2+Reffect2+dpsidRdZpart2

c            if (i>0.99*ntot) then
c               write(*,*) 'shattemp1',dpsi2part1,Reffect1,dpsidRdZpart1
c     1              ,fpol3,ffprime3, shattemp1(i)
c               write(*,*) 'shattemp2',dpsi2part2,Reffect2,dpsidRdZpart2
c     1              ,shattemp2(i)
c            end if
         end do
         sol_unit = 10
         open(sol_unit, file='shattemp.out') !, status='new')
       
         do j=1, nt2
            do i=1, nr
               inext = i + (j-1)*nr
               write (sol_unit, *) Rt3(inext), 
     $              Zt3(inext), shattemp1(inext)
            end do
         end do
         close(sol_unit)

         open(sol_unit, file='fitmiller_trin.out')
         write(sol_unit, *) 'rho, R0, dR0/dr, r/R0, kappa
     1, s_kappa, delta, s_delta, q , s, alpha'
         do k=1, ntpsi
            write (sol_unit, *) trho(k), tR0mil(k), tdR0mil(k)
     1           ,tiaspmil(k),tkappamil(k),tdkappamil(k)
     2           ,tdeltamil(k),tddeltamil(k),tqmil(k),tshatmil(k)
     3           ,talphamil(k)
         end do  
         close(sol_unit)

         write(fiter22,*) 'ntpsi,nttin'
         write(fiter22,*) ntpsi,nttin
         write(fiter22,*) 'tpsi'
         write(fiter22,*) tpsi
         write(fiter22,*) 'ttin'
         write(fiter22,*) ttin

         write(fiter22,*) 'trho'
         write(fiter22,*) trho(1:ntpsi)
         write(fiter22,*) 'tvol(1:ntpsi)'
         write(fiter22,*) tvol(1:ntpsi)
         write(fiter22,*) 'dpsidrho(1:ntpsi)'
         write(fiter22,*) dpsidrho(1:ntpsi)
         write(fiter22,*) 'dVdpsi(1:ntpsi)'
         write(fiter22,*) tdVdpsi(1:ntpsi)
         write(fiter22,*) 'shat(1:ntpsi)'
         write(fiter22,*)  tshatmil(1:ntpsi)
         write(fiter22,*) 'alpha(1:ntpsi)'
         write(fiter22,*)  talphamil(1:ntpsi)

         write(fiter22,*) 'Phitrin(1:ntpsi*nttin)'
         write(fiter22,*) phitrin
         write(fiter22,*) 'Rtrin(1:ntpsi*nttin)'
         write(fiter22,*) Rtrin
         write(fiter22,*) 'Ztrin(1:ntpsi*nttin)'
         write(fiter22,*) Ztrin
         write(fiter22,*) 'Btot(1:ntpsi*nttin)'
         write(fiter22,*) btot
         write(fiter22,*) 'Bpol(1:ntpsi*nttin)'
         write(fiter22,*) bpol
         write(fiter22,*) 'gradpar(1:ntpsi*nttin)'
         write(fiter22,*) gradpar
         write(fiter22,*) 'gbdrift(1:ntpsi*nttin)'
         write(fiter22,*) gbdrift
         write(fiter22,*) 'gbdrift0(1:ntpsi*nttin)'
         write(fiter22,*) gbdrift0
         write(fiter22,*) 'cvdrift(1:ntpsi*nttin)'
         write(fiter22,*) cvdrift
         write(fiter22,*) 'cvdrift0(1:ntpsi*nttin)'
         write(fiter22,*) cvdrift0
         write(fiter22,*) 'local shear by d2alpha/dtdp:S(1:ntpsi*nttin)'
         write(fiter22,*) shatlocal1
         write(fiter22,*) 'local shear by d2alpha/dpdt:S(1:ntpsi*nttin)'
         write(fiter22,*) shatlocal2
         write(fiter22,*) 'local shear from direct cal:S(1:ntpsi*nttin)'
         write(fiter22,*) shatlocal3

c         epsFIT = 1.0d-13
c         file_name='fitMillerTr.out'
c         call fitMiller(ntpsi,itcon,tRcon,tZcon,tRavg,tZavg
c     1        ,epsFIT,tR0mil,tiaspmil,tkappamil,tdeltamil)

c         call printMiller(ntpsi,tpsi,tR0mil,tiaspmil
c     1        ,tkappamil,tdeltamil,file_name)

         if (iogyropsi.eq.1) then
            fogyro=23

            tzeros1=0.0d0
            tzeros2=0.0d0
            sformat='(5E20.10)'
            open(fogyro, file='ogyropsi.dat') !, status='new')
            write(fogyro,*) 'NPSI'
            write(fogyro,*) ntpsi
            write(fogyro,*) 'NCHI'
            write(fogyro,*) nttin
            write(fogyro,*) 'R0EXP'
            write(fogyro,*) R0
            write(fogyro,*) 'B0EXP'
            write(fogyro,*) dabs(F0/R0)
            write(fogyro,*) 'PSI'
            write(fogyro,sformat) (tpsi-1.0d0)/lambda
            write(fogyro,*) 'CHI'
            write(fogyro,sformat) ttin 
         
            do k=1,ntpsi
               Rgeom(k)=(Rtrin((k-1)*nttin+1)
     1              +Rtrin((k-1)*nttin+nttin/2+1))/2.0d0 ! assume ttin is uniform and ntpsi is even
               ageom(k)=Rtrin((k-1)*nttin+1)-Rgeom(k)
            end do
            write(fogyro,*) 'Rgeom'
            write(fogyro,sformat) Rgeom
            write(fogyro,*) 'ageom'
            write(fogyro,sformat) ageom
            write(fogyro,*) 'q'
            write(fogyro,sformat) tqmil
            write(fogyro,*) 'dqdpsi'
            write(fogyro,sformat) tdqdpsi
            write(fogyro,*) 'd2qdpsi2'
            write(fogyro,sformat) tzeros1    !not needed

            write(fogyro,*) 'p'
            write(fogyro,sformat) tpcon
            write(fogyro,*) 'dpdpsi'
            write(fogyro,sformat) tdpdpsi
            write(fogyro,*) 'f'
            write(fogyro,sformat) tfpolcon
            write(fogyro,*) 'fdfdpsi'
            write(fogyro,sformat) tdfdpsi*tfpolcon
        
            write(fogyro,*) 'V'
            write(fogyro,sformat) tvol
            write(fogyro,*) 'rho_t' !For normalized toroidal flux, make irho=3, and then rho=sqrt(tor_flux)
           frho_t=dsqrt(dabs(phitot*2.0/(tfpolcon(ntpsi)/Rgeom(ntpsi))))
            write(*,*) 'frho_t',frho_t

            write(fogyro,sformat) trho
            
            write(fogyro,*) 'Shear' ! magnetic shear defined in Miller PoP1998
            write(fogyro,sformat) tshatmil
            write(fogyro,*) 'dSheardpsi' ! magnetic shear defined in Miller PoP1998
            write(fogyro,sformat)  tzeros1    !not needed
            write(fogyro,*)  'kappa'
            write(fogyro,sformat)  tkappamil
            write(fogyro,*)  'delta_lower'  !assumed up-down symmetry in Miller fitting
            write(fogyro,sformat)  tdeltamil
            write(fogyro,*)  'delta_upper'
            write(fogyro,sformat)   tdeltamil
            write(fogyro,*)  'dVdpsi'
            write(fogyro,sformat)   tdvdpsi
            write(fogyro,*)  'dpsidrhotor'
            write(fogyro,sformat)   dpsidrho

            
            write(fogyro,*)  'GDPSI_av'
            write(fogyro,sformat)   gdpsi_av
            write(fogyro,*)  'radius_av'
            write(fogyro,sformat)  tzeros1   !not needed
            write(fogyro,*)  'R_av'
            write(fogyro,sformat)  tzeros1   !not needed
            write(fogyro,*)  'TE'
            write(fogyro,sformat)  tzeros1   !not needed
            write(fogyro,*)  'DTEDPSI'
            write(fogyro,sformat)  tzeros1   !not needed

            write(fogyro,*)  'NE'
            write(fogyro,sformat)  tzeros1   !not needed
            write(fogyro,*)  'DNEDPSI'
            write(fogyro,sformat)  tzeros1   !not needed

            write(fogyro,*)  'TI'
            write(fogyro,sformat)  tzeros1   !not needed
            write(fogyro,*)  'DTIDPSI'
            write(fogyro,sformat)  tzeros1   !not needed
            write(fogyro,*)  'NI'
            write(fogyro,sformat)  tzeros1   !not needed
            write(fogyro,*)  'DNIDPSI'
            write(fogyro,sformat)  tzeros1   !not needed

            write(fogyro,*)  'ZEFF'
            write(fogyro,sformat)  tzeros1   !not needed
            write(fogyro,*)  'SIGNEO'
            write(fogyro,sformat)  tzeros1   !not needed
            write(fogyro,*)  'JASBAV'
            write(fogyro,sformat)  tzeros1   !not needed

            !2D array
            write(fogyro,*) 'g11'
            write(fogyro,sformat) tzeros2  !not needed
            write(fogyro,*) 'g12'
            write(fogyro,sformat) tzeros2  !not needed
            write(fogyro,*) 'g22'
            write(fogyro,sformat) tzeros2  !not needed
            write(fogyro,*) 'g33'
            write(fogyro,sformat) tzeros2  !not needed
            write(fogyro,*) 'B'
            write(fogyro,sformat) tzeros2  !not needed
            write(fogyro,*) 'dBdpsi'
            write(fogyro,sformat) tzeros2  !not needed
            write(fogyro,*) 'dBdchi'
            write(fogyro,sformat) tzeros2  !not needed
            write(fogyro,*) 'dPsidR'
            write(fogyro,sformat) tzeros2  !not needed
            write(fogyro,*) 'dPsidZ'
            write(fogyro,sformat) tzeros2  !not needed
            write(fogyro,*) 'dChidR'
            write(fogyro,sformat) tzeros2  !not needed
            write(fogyro,*) 'dChidZ'
            write(fogyro,sformat) tzeros2  !not needed
       
            write(fogyro,*) 'Jacobian'
            write(fogyro,sformat) tzeros2  !not needed

            do i=1,ntpsi
               do j=1,nttin
                  inext1=(i-1)*nttin+j
                  inext2=(j-1)*ntpsi+i
                  RtrinT(inext2)=Rtrin(inext1)
                  ZtrinT(inext2)=Ztrin(inext1)
               end do
            end do

            write(fogyro,*) 'R'
            write(fogyro,sformat) RtrinT  !not needed
            write(fogyro,*) 'Z'
            write(fogyro,sformat) ZtrinT  !not needed

            close(fogyro)
         end if

      end if

      write(*,*) 'finished finding physical quantities'
     
      if (iprintmap.eq.1) then
         iw=41
         itype2=2
         itype3=3
         call quaplot2(iw,Rt1,Zt1,nt1,itype3,Rt3,Zt3,nr*nt2,itype2,
     1       'image of uniform polar grid under inverse map*')
      end if

      
     
      if (iprintsoldisk.eq.1) then
         call printSolDisk(psi)
      end if
  

      if (iprintsol.eq.1) then
         select case (iptype)
         case (0,3)               !If there are exact solutions or EFIT solutions
            call printSolex(psi, psiR, psiZ, psiRex,psiZex)
         case(1,2,4)
            call printSol(psi, psiR, psiZ)
         end select
      end if

      if (iprintcon.eq.1) then
         select case (iptype)
         case (0)             !If there are exact solutions
            call printContourEx(ncon, icon, rndcon, tcon, psicon, Rcon,
     1     Zcon, psiRcon, psiZcon, psiconex, psiRconex, psiZconex)
         case (1,11,2,3,4)
            call printContour(ncon, icon, rndcon, tcon, psicon, Rcon,
     1     Zcon, psiRcon, psiZcon)
         end select
         write(*,*) 'finished printing contours'
      end if

      if (iprintQS.eq.1) then
         select case (iptype)
         case (0,3)               !If there are exact solutions or EFIT solutions
            call printQSex(ncon, psicon, qpsiconex, qcon, shat)
         case (1,2,4)
            call printQS(ncon, psicon, qcon, shat)
         end select
         write(*,*) 'finished printing Q and S profiles'
      end if

 

      if (iprintEFIT.ge.1) then

         write(*,*) 'before printEFIT'
         call printEfitG (npsi, ncon, qcon, acmcur0)
c printEfitG (npsi,ncon,cftmq,pprimch,ffprimch,presch,fpolch
c     1     ,qcon)
         Write(*,*) 'finished printing gfile for EFIT output format'
      end if

      return
      end
!---------------------------------------------------------
!--------------------------------------------------------
!     This routine postproc1 is called during nonlinar iterations 
!     to find some contour intergral values, 
!     especially for q solver ECOMQ and j solver ECOMJ
!---------------------------------------------------------
      subroutine postproc1


      use arrays, only:  nsub, kcheb, kLag, nr, nt2, ntot, ksamp2 
      use arrays, only:  zk, dzdw2k, dzdww2k, rnd, tnd, lambda, itftype

      use arrays, only:  psii, dpsidr, dpsidrr, psiex, iiter, iLag
      use arrays, only:  u, ur, uth, urr, urt, utt, torcur, iecom
      use arrays, only:  wffprim, fpolcon2, ffprimpre, ffprimpre1, q0
      use arrays, only:  ffprimgoal, fpolch, ffprimpre2, ffprimch, F0
      use arrays, only:  nchq, psiichq, cftmq, qpsich, pprimch, qpsich1
 
      use arrays, only: epsmaxdist, Rmaxpsi, Zmaxpsi, fpolsgn, jpsich 
      use arrays, only: spdef, spbx, spxb, iscale, lambdaex,maxdpsi
      use arrays, only: maxdpsiq, maxdpsij, FF0, maxiterF, Rmax, Rmaxis
      use arrays, only: fiter1, fiter2, fiter3,fiter4,fiter5, fiter6
      use arrays, only: nchy, psiichy, cftmy, spybx, spyxb,spydef, fpass
      use arrays, only: fpass,etai,Zi,nue_star,nui_star, jBSch, ijtype
      use arrays, only: ne, Te, Vloop, jOhmch, ifpol, Fedge, iptable
      use arrays, only: rhoch, dpsiidrhoch, iptable,  iqconstraint
      use arrays, only: ptflowch,ptfprimch,presch,phitot,phiichq

      implicit double precision(a-h, o-z)
      
      parameter (ncon0 = 100)
      real *8 psicon(ncon0)
      real *8 psicontmp(ncon0)
      real *8 fpolcon(ncon0),alphafpol(ncon0)
      real *8 fpoltmp1(ncon0),fpoltmp2(ncon0)
      real *8 ffprimnum(ncon0)
      real *8 Ia(ncon0),Ib(ncon0),Ic(ncon0)
      real *8 Lp(ncon0),Id(ncon0),Ie(ncon0)
      real *8 Itf1(ncon0),Itf2(ncon0),Ibtf(ncon0)
      real *8 rndcon(3*nt2*ncon0)
      real *8 tcon(3*nt2*ncon0)
      real *8 urcon(3*nt2*ncon0)
      real *8 utcon(3*nt2*ncon0)
      real *8 Rcon(3*ksamp2*ncon0)
      real *8 Zcon(3*ksamp2*ncon0)
      real *8 Ravg(ncon0), Zavg(ncon0)
      real *8 psiRcon(3*ksamp2*ncon0)
      real *8 psiZcon(3*ksamp2*ncon0)
      real *8 tmaxp(3*ksamp2*ncon0)
      real *8 psiconex(3*ksamp2*ncon0)
      real *8 psiRconex(3*ksamp2*ncon0)
      real *8 psiZconex(3*ksamp2*ncon0)
      real *8 maxF, errF, epsF

      real *8 Rmido(ncon0), Rmidi(ncon0)
      integer *4 icon(ncon0), sol_unit, iiterF
 
      complex *16 ima
      real *8 ffprimnext(ncon0), chcoeff(ncon0), mu0 
      real *8 dZdtmido(ncon0), dZdtmidi(ncon0)
      save

 
      pi = 4.0d0*datan(1.0d0)
      mu0=(4*pi*1.0d-7)
      ima = (0,1)
      ncon=nchq
      
      psicon(1:ncon)=psiichq(1:nchq)
c      write(*,*) 'psiichq',psiichq(1:nchq)
c      write(*,*) 'psicon',psicon(1:ncon)
      call findcontour(ncon,psicon,nr,nt2,rnd,tnd,psii,dpsidr,dpsidrr
     1        ,rndcon,tcon,urcon,utcon)
c      write(*,*) 'psicon',psicon(1:ncon)
      call sortcontour(ncon,psicon,zk,dzdw2k,rndcon
     1        ,tcon,urcon,utcon,Rcon,Zcon,icon,Ravg,Zavg
     2        ,tmaxp,psiRcon,psiZcon,psiconex,psiRconex,psiZconex
     3        ,Rmido, Rmidi)
c      write(*,*) 'psicon',psicon(1:ncon)
      write(*,*) 'contour finished'
    
            
      call calfint(ncon,psicon,Rcon,Zcon,tmaxp,icon
     1     ,psiRcon,psiZcon,Ia,Ib,Ic,Lp,Id)

      if ((itftype.ne.0).and.((iLag.eq.1).or.(iecom.ne.0))) then
         call calfinttf(ncon,psicon,Rcon,Zcon,tmaxp,icon
     1        ,psiRcon,psiZcon,Itf1, Itf2)

         do k=1,ncon
            ptfbyp=ptflowch(k)/presch(k)
            dptfbydp=ptfprimch(k)/pprimch(k)
            Ibtf(k)=Itf1(k)+(dptfbydp-ptfbyp)*Itf2(k)
         end do
c         write(*,*) 'Itf1',Itf1(1:ncon)
c         write(*,*) 'Itf2',Itf2(1:ncon)
c         write(*,*) 'Ibtf',Ibtf(1:ncon)
      end if

c      write(*,*) 'Ia',Ia(1:ncon)
c      write(*,*) 'Ib',Ib(1:ncon) 
 
      select case (iecom)
      case (0)  !F constraint
c         write(*,*) 'maxdpsi',maxdpsi,maxdpsij
         if ((iscale.eq.5).or.(iscale.eq.6)) then
            if (maxdpsi.lt.maxdpsij) then
               call matchtorcur(ncon,lambda,Ia,cftmq,torcur,factorl)
               write(*,*) 'torcur, factorl',torcur, factorl
               FF0=FF0/factorl
               write(*,*) 'FF0', FF0
c               lambdaex=lambda*factorl
c               write(*,*) 'lambdaex', lambdaex
            end if
         end if
         
c         if (iLag.eq.1) then
            if (ifpol.eq.0) then
               F0tmp=dabs(F0)*fpolsgn
               call findFbyFFp0(ncon,lambda,psiichq,ffprimch,
     1              spbx,F0tmp,fpolch)
            else if (ifpol.eq.1) then
               Fedgetmp=dabs(Fedge)*fpolsgn
               call findFbyFFp1(ncon,lambda,psiichq,ffprimch,
     1              spxb,Fedgetmp,fpolch)
            end if
c         end if

         qpsich(1:ncon)=Ic(1:ncon)*dabs(lambda/(2.0d0*pi)
     1        *fpolch(1:ncon)) !*fpolsgn

c         write(*,*) 'qpsich',qpsich(1:ncon)
      case (1)  !j_par constraint
c         write(*,*) 'maxdpsi',maxdpsi,maxdpsij
         if (maxdpsi.lt.maxdpsij) then
            if (iscale.eq.5) then
               call matchtorcur(ncon,lambda,Ia,cftmq,torcur,factorl)
               write(*,*) 'torcur, factorl',torcur, factorl
               lambdaex=lambda*factorl
               write(*,*) 'lambdaex', lambdaex
            else if (iscale.eq.6) then
               call matchtorcur(ncon,lambda,Ia,cftmq,torcur,factorl)
               write(*,*) 'torcur, factorl',torcur, factorl
               lambdaex=lambda*factorl
               write(*,*) 'lambdaex', lambdaex
               call findF0(ncon,lambdaex,Ic,cftmq,q0,F0)
               write(*,*) 'q0, F0',q0, F0
            end if
            
         end if

         if (ijtype.eq.3) then
            call findfpass(ncon,lambda,psiichq,Rcon,Zcon,tmaxp,icon
     1           ,psiRcon,psiZcon,fpolch,nchy,psiichy,spydef,fpass)
            call findBScur(ncon,lambda,Ib,Ic,pprimch !Evaluate Boostrap current
     1           ,fpolch,fpass,etai,Zi,nue_star,nui_star,jBSch)
            call findOhmcur(ncon,lambda,ne,Te,fpass,Zi
     1           ,nue_star,Vloop,jOhmch)
            jpsich(1:ncon)=(jBSch(1:ncon)+jOhmch(1:ncon))*mu0  !update Jpar by boostrap and Ohmic current
c            write(*,*) 'JBsch',jBSch(1:ncon)
c            write(*,*) 'JOhmch',jOhmch(1:ncon)
c            write(*,*) 'Jpsich',jpsich(1:ncon)

         else if (ijtype.eq.5) then
            call findfpass(ncon,lambda,psiichq,Rcon,Zcon,tmaxp,icon
     1           ,psiRcon,psiZcon,fpolch,nchy,psiichy,spydef,fpass)
            call findOhmcur(ncon,lambda,ne,Te,fpass,Zi
     1           ,nue_star,Vloop,jOhmch)
            jpsich(1:ncon)=(jOhmch(1:ncon))*mu0  !update Jpar by  Ohmic current
c            write(*,*) 'JOhmch',jOhmch(1:ncon)
c            write(*,*) 'Jpsich',jpsich(1:ncon)
         end if

         iiterF=0
         maxiterF=1 !inner iteration for self-consistent F and dF/dpsi. If maxiterF=1, no inner iteration 
         epsF=1.0d-5
         maxF = 1.0d5
         errF = 1.0d5
         do while ((dabs(errF).gt.dabs(epsF*maxF)).and.
     1        (iiterF.lt.maxiterF))
            iiterF=iiterF+1
            if (ifpol.eq.0) then
               call findFbyFFp0(ncon,lambda,psiichq,ffprimch,
     1              spbx,F0,fpolch)
            else if (ifpol.eq.1) then
               call findFbyFFp1(ncon,lambda,psiichq,ffprimch,
     1              spxb,Fedge,fpolch)
            end if
c            write(*,*) 'fpolch',fpolch(1:ncon)
            if (itftype.eq.0) then
               call findJparfac(ncon,lambda,psicon,Ia,Ib,Ic
     1              ,fpolch,jpsich,cftmq,ffprimpre1,ffprimpre2)
            else if (itftype.ne.0) then

               call findJparfac(ncon,lambda,psicon,Ia,Ibtf,Ic
     1              ,fpolch,jpsich,cftmq,ffprimpre1,ffprimpre2)
            end if

            ffprimnext(1:ncon)=ffprimpre1(1:ncon)
     1           -ffprimpre2(1:ncon)*mu0*pprimch(1:ncon)
            maxF=0.0d0
            errF=0.0d0
            do i=1,ncon
               if (dabs(ffprimnext(i)).gt.maxF) then
                  maxF=ffprimnext(i)
               end if
               errF1=ffprimnext(i)-ffprimch(i)
               if (dabs(errF1).gt.errF) then
                  errF=errF1
               end if
            end do
            ffprimch(1:ncon)=ffprimnext(1:ncon)
            write(*,*) 'iiterF,fpprimch',iiterF,ffprimch(1:ncon)
            write(*,*) 'maxF,errF',maxF,errF
         end do

         write(fiter3,*) ffprimpre1(1:ncon)
         write(fiter5,*) ffprimpre2(1:ncon) 
         write(fiter6,*) ffprimch(1:ncon) 
         qpsich(1:ncon)=Ic(1:ncon)*dabs(lambda/(2.0d0*pi)
     1        *fpolch(1:ncon)) !*fpolsgn

         write(*,*) 'qpsich',qpsich(1:ncon)
      case (2)  ! q constraint
c         write(*,*) 'maxdpsi',maxdpsi,maxdpsiq
         if ((iscale.eq.2) .and. (maxdpsi.lt.maxdpsiq)) then
            call matchtorcur(ncon,lambda,Ia,cftmq,torcur,factorl)
            write(*,*) 'torcur, factorl',torcur, factorl
            write(*,*) 'factorl'
            lambdaex=lambdaex*factorl
            write(*,*) 'lambdaex=',lambdaex
         end if

c         if 
         call findfpol(ncon,lambda,psicon,qpsich,Ic
     1        ,fpoltmp1)
c         factorl=fpoltmp1(ncon)/Fedge
c         lambdaex=lambdaex*factorl
c         write(*,*) 'factorl', factorl
c         write(*,*) 'lambdaex=',lambdaex
         if (iqconstraint.eq.0) then
c            lambda=lambdaex
            call findFbyFFp1(ncon,lambda,psicon,ffprimnum,
     1              spxb,fpoltmp1(ncon),fpoltmp2)
c            wfpoltmp1=0.2d0  !1.0d0 is better?
            wfpoltmp1=1.0d0  !1.0d0 is better for the convergence
         else if (iqconstraint.eq.1) then
            call findFbyFFp1(ncon,lambda,psicon,ffprimnum,
     1           spxb,fedge,fpoltmp2)
            !wfpoltmp1=0.0d0
            wfpoltmp1=0.2d0 ! to make fpoltmp1 and fpoltmp2 matched
         end if
        
         fpolcon(1:ncon)=wfpoltmp1*fpoltmp1(1:ncon)
     1        +(1.0d0-wfpoltmp1)*fpoltmp2(1:ncon)
         fpolch(1:ncon)=fpolcon(1:ncon)
         fpolcon2(1:ncon)=fpolcon(1:ncon)*fpolcon(1:ncon)/2.0d0
c         write(*,*) 'fpoltmp1',fpoltmp1(1:ncon)
c         write(*,*) 'fpoltmp2',fpoltmp2(1:ncon)
      write(*,*) 'fpolcon',fpolcon(1:ncon)

c      write(*,*) 'Ic',Ic(1:ncon)
      write(*,*) '<R-2>',Ic(1:ncon)/Ib(1:ncon)
      call findffprimnext(ncon,lambda,psicon,cftmq,fpolcon2
     1        ,ffprimnext)
  
         if (itftype.eq.0) then
            call findfpol0(ncon,lambda,psicon,Ia,Ib,Ic
     1           ,pprimch,cftmq,ffprimpre)

            call findfpol1(ncon,lambda,psicon,Ia,Ib,Ic
     1           ,pprimch,fpolcon,qpsich,cftmq,ffprimpre1
     1           ,ffprimpre2)
            call findfpol11(ncon,lambda,psicon,Ia,Ib,Ic
     1           ,pprimch,fpolcon,qpsich,cftmq,ffprimnum)

         else if (itftype.ne.0) then
 

            call findfpol0(ncon,lambda,psicon,Ia,Ibtf,Ic
     1           ,pprimch,cftmq,ffprimpre)

            call findfpol1(ncon,lambda,psicon,Ia,Ibtf,Ic
     1           ,pprimch,fpolcon,qpsich,cftmq,ffprimpre1
     1           ,ffprimpre2)
            call findfpol11(ncon,lambda,psicon,Ia,Ibtf,Ic
     1           ,pprimch,fpolcon,qpsich,cftmq,ffprimnum)
         end if

       
         
c     call findfpolalpha(ncon,lambda,fpolcon,Ia,Ic,alphafpol)
c     ffprimch(1:ncon)=ffprimnum(1:ncon) 
         ffprimgoal(1:ncon)= ffprimnext(1:ncon)
         write(fiter1,*) ffprimpre(1:ncon)
         write(fiter2,*) ffprimnext(1:ncon) !/lambda**2
         write(fiter3,*) ffprimpre1(1:ncon)
c     write(fiter4,*) alphafpol(1:ncon)
         write(fiter5,*) ffprimpre2(1:ncon) 
c         write(fiter6,*) 'iiter=',iiter,ffprimnum(1:ncon) !-ffprimpre(1:ncon)! ffprimnum(1:ncon)
         ffprimch(1:ncon)=ffprimnum(1:ncon)

 
      end select

      if (iLag.eq.1) then
         call findLagrange(ncon,lambda,psiichq,Ia,Ib,Ic
     1        ,Itf1,pprimch,fpolch,cftmq,spdef,spbx,spxb
     2        ,WLag)
      end if

      if (iptable.eq.2) then
c         rhoch(1:ncon)=(Rmido(1:ncon)-Rmaxis)/(Rmax-Rmaxis)
         call findtorflux(ncon,lambda,psiichq,qpsich,
     1     spbx,spdef, phitot,phiichq)
         call findrhoch(ncon,Rmido, Rmidi)
         write(*,*) 'rhoch',rhoch(1:ncon)
         call chftransq(chcoeff, rhoch, ncon, cftmq)
         do i=1,nchq
            psin=1.0d0-psiichq(i)
            call chderiv(psin,nchq,chcoeff,drhodpsi)
            dpsiidrhoch(i)=-1.0d0/drhodpsi !used normalized psi
         end do
c         write(*,*) 'dpsiidrhoch',dpsiidrhoch(1:ncon)
      end if

      return
      end
!-----------------
      subroutine scalesol(lambda, psiB, psi)

      use arrays, only:  ntot, psii
      use arrays, only: dpsidr, dpsidrr
      use arrays, only: u, ur, uth, urr, urt, utt

      implicit double precision(a-h, o-z)

      real *8 psi(*), lambda, psiB

      psi(1:ntot) = psii(1:ntot)/lambda +psiB

      u(1:ntot) = u(1:ntot)/lambda
      ur(1:ntot) = ur(1:ntot)/lambda
      uth(1:ntot) = uth(1:ntot)/lambda
      urr(1:ntot) = urr(1:ntot)/lambda
      urt(1:ntot) = urt(1:ntot)/lambda
      utt(1:ntot) = utt(1:ntot)/lambda
      dpsidr(1:ntot) = dpsidr(1:ntot)/lambda
      dpsidrr(1:ntot) = dpsidrr(1:ntot)/lambda

      return
      end

      subroutine calDer(psiR, psiZ, psiRR, psiRZ, psiZZ)

      use arrays, only:  nr, ntot, nt2, kcheb, nsub  
      use arrays, only:  dpsidr, dpsidrr, u, ur, uth, urr, urt, utt
      use arrays, only:  Rt3, Zt3, rnd, tnd, dzdw3, dzdww3
      use arrays, only:  dwdz3, dwdzz3

      implicit double precision(a-h, o-z)

      real *8 psiR(*), psiZ(*)
      real *8 psiRR(*),psiRZ(*),psiZZ(*)

      real *8 ux(ntot)
      real *8 uy(ntot)
      real *8 uxx(ntot)
      real *8 uxy(ntot)
      real *8 uyy(ntot)
      real *8 uxdisk(ntot)
      real *8 uydisk(ntot)
      real *8 uxxdisk(ntot)
      real *8 uxydisk(ntot)
      real *8 uyydisk(ntot)
      real *8 dRdr(ntot)
      complex *16 ima
      integer *4 i,j, inext, ntheta

      ima = (0,1)
      ntheta=nt2

      do j=1, ntheta
         do i=1, nr
            inext = i + (j-1)*nr
            uxdisk(inext)=ur(inext)*dcos(tnd(j))
     1           +1/rnd(i)*uth(inext)*(-dsin(tnd(j)))
            uydisk(inext)=ur(inext)*dsin(tnd(j))
     1           +1/rnd(i)*uth(inext)*(dcos(tnd(j)))
            uxxdisk(inext)=urr(inext)*dcos(tnd(j))**2
     1           +1/rnd(i)*ur(inext)*dsin(tnd(j))**2
     1           +1/rnd(i)**2*utt(inext)*dsin(tnd(j))**2
     1           +2/rnd(i)**2*uth(inext)*dsin(tnd(j))*dcos(tnd(j))
     1           -2/rnd(i)*urt(inext)*dsin(tnd(j))*dcos(tnd(j))
           
            uyydisk(inext)=urr(inext)*dsin(tnd(j))**2
     1           +1/rnd(i)*ur(inext)*dcos(tnd(j))**2
     1           +1/rnd(i)**2*utt(inext)*dcos(tnd(j))**2
     1           -2/rnd(i)**2*uth(inext)*dsin(tnd(j))*dcos(tnd(j))
     1           +2/rnd(i)*urt(inext)*dsin(tnd(j))*dcos(tnd(j))

            uxydisk(inext)=urr(inext)*dsin(tnd(j))*dcos(tnd(j))
     1           -1/rnd(i)*ur(inext)*dsin(tnd(j))*dcos(tnd(j))
     1           -1/rnd(i)**2*utt(inext)*dsin(tnd(j))*dcos(tnd(j))
     1           -1/rnd(i)**2*uth(inext)*dcos(2*tnd(j))
     1           +1/rnd(i)*urt(inext)*dcos(2*tnd(j))

            dzdwr=dzdw3(inext)
            dzdwi=(-ima)*dzdw3(inext)
            dRdr(inext)=dzdwr*dcos(tnd(j))-dzdwi*dsin(tnd(j))

            dwdzr=dwdz3(inext)
            dwdzi=(-ima)*dwdz3(inext)
            ux(inext)=uxdisk(inext)*dwdzr+uydisk(inext)*dwdzi
            uy(inext)=uxdisk(inext)*(-dwdzi)+uydisk(inext)*dwdzr


            dwdzzr=dwdzz3(inext)
            dwdzzi=(-ima)*dwdzz3(inext)

            uxx(inext)=uxxdisk(inext)*dwdzr**2+uxdisk(inext)*dwdzzr
     1           +2*uxydisk(inext)*dwdzr*dwdzi
     1           +uyydisk(inext)*dwdzi**2+uydisk(inext)*dwdzzi

            uyy(inext)=uxxdisk(inext)*dwdzi**2-uxdisk(inext)*dwdzzr
     1           -2*uxydisk(inext)*dwdzr*dwdzi
     1           +uyydisk(inext)*dwdzr**2-uydisk(inext)*dwdzzi

            uxy(inext)=-uxxdisk(inext)*dwdzr*dwdzi-uxdisk(inext)*dwdzzi
     1           +uxydisk(inext)*(dwdzr**2-dwdzi**2)
     1           +uyydisk(inext)*dwdzr*dwdzi+uydisk(inext)*dwdzzr

            
            psiR(inext) = u(inext)/(2.0d0*sqrt(Rt3(inext)))
     1           + ux(inext)*sqrt(Rt3(inext))
            psiZ(inext) = uy(inext)*sqrt(Rt3(inext))
            psiRR(inext) = -u(inext)/(4.0d0*Rt3(inext)**1.5)
     1           +ux(inext)/(sqrt(Rt3(inext)))
     1           + uxx(inext)*sqrt(Rt3(inext))
            psiRZ(inext) = uy(inext)/(2.0d0*sqrt(Rt3(inext)))
     1           + uxy(inext)*sqrt(Rt3(inext))
            psiZZ(inext) = uyy(inext)*sqrt(Rt3(inext))

         end do
      end do

      return
      end

      subroutine printSolDisk(psi)

      use arrays, only:  nr, ntot, nt2, dpsidr  
      use arrays, only:  rnd, tnd, u, ur, uth, urr, urt, utt
      implicit double precision(a-h, o-z)

      real *8 x(ntot)
      real *8 y(ntot)
      real *8 psi(*)

      integer *4 sol_unit
      integer *4 i,j, inext, ntheta

      ntheta = nt2
      sol_unit = 10
      open(sol_unit, file='sol_disk.out')!, status='new')
      do j=1, ntheta
         do i=1, nr
            inext = i + (j-1)*nr

            write (sol_unit, *) tnd(j), 
     $           rnd(i), u(inext), ur(inext), uth(inext),
     $           urr(inext), urt(inext),utt(inext)

         end do
      end do
      close(sol_unit)


      open(sol_unit, file='sol_diskxy.out') !, status='new')
      do j=1, ntheta
         do i=1, nr

            inext = i + (j-1)*nr
            x(inext) = rnd(i)*cos(tnd(j))
            y(inext) = rnd(i)*sin(tnd(j))
            write (sol_unit, *) x(inext), 
     $           y(inext),  u(inext), psi(inext), dpsidr(inext)

         end do
      end do
      close(sol_unit)

      iw=40
      itype2=3
      call quaplot(iw,x,y,nr*ntheta,itype2,
     1       'image of uniform polar grid under inverse map*')
      
      return
      end
    
      subroutine printSolex(psi, psiR, psiZ, psiRex,psiZex)

      use arrays, only: Rt1,Zt1,nt1
      use arrays, only: nr, ntot, nt2, Rt3, Zt3, psiex
      use arrays, only:  rnd, tnd, u, ur, uth, urr, urt, utt
      implicit double precision(a-h, o-z)

      real *8 x(ntot)
      real *8 y(ntot)
      real *8 psi(*), psiR(*), psiZ(*), psiRex(*),psiZex(*)

      integer *4 sol_unit
      integer *4 i,j, inext, ntheta

      ntheta = nt2
      sol_unit = 10
c      open(sol_unit, file='solu_orig.out') !, status='new')
c      do j=1, ntheta
c         do i=1, nr
c            inext = i + (j-1)*nr
c            write (sol_unit, *) Rt3(inext), 
c     $           Zt3(inext), u(inext), ux(inext),uy(inext),
c     $           uxx(inext), uxy(inext),uyy(inext)
c         end do
c      end do
c      close(sol_unit)

      open(sol_unit, file='sol_orig.out') !, status='new')
      do j=1, ntheta
         do i=1, nr
            inext = i + (j-1)*nr
            write (sol_unit, *) Rt3(inext), 
     $           Zt3(inext), psi(inext), psiex(inext)
         end do
      end do
      close(sol_unit)

      open(sol_unit, file='sol_origder.out') !, status='new')
      do j=1, ntheta
         do i=1, nr
            inext = i + (j-1)*nr
            write (sol_unit, *) Rt3(inext), 
     $           Zt3(inext), psiR(inext), psiRex(inext),
     $           psiZ(inext), psiZex(inext)
         end do
      end do
      close(sol_unit)

      return
      end

      subroutine printContourEx(ncon, icon, rndcon, tcon, psicon, Rcon,
     1     Zcon, psiRcon, psiZcon, psiconex, psiRconex, psiZconex)

      use arrays, only:  nr, ntot, nt2, kcheb, nsub
      implicit double precision(a-h, o-z)

      real *8 rndcon(*), tcon(*), psicon(*), Rcon(*), Zcon(*),
     1   psiRcon(*), psiZcon(*), psiconex(*), psiRconex(*), psiZconex(*)

      integer *4 sol_unit, ncon, icon(*)
      integer *4 i,j, inext, ntheta, istart, iend, nin

      ntheta = nt2
      sol_unit = 10

      open(sol_unit, file='contour_disk.out') !, status='new')
      itotc = 0
      do k=1, ncon
         do j=1, ntheta*ksamp2
               itotc=itotc+1
               write (sol_unit,*) rndcon(itotc), 
     $           tcon(itotc), psicon(k)
          
         end do
      end do
      close(sol_unit)

      open(sol_unit, file='contour_orig.out') !, status='new')
      itotc = 0
      do k=1, ncon
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         do j=1, nin !ntheta*ksamp2
               itotc=itotc+1
               write (sol_unit, *) Rcon(itotc), 
     $           Zcon(itotc), psicon(k), psiconex(itotc),
     $           psiRcon(itotc), psiRconex(itotc)
         end do
      end do
      close(sol_unit)

      return
      end

      subroutine printMiller(ncon,psicon,R0,iasp
     1     ,kappa,delta, Ib)

      use arrays, only:  nchq, psiichq,  cftmq
      use arrays, only: reps ,spbx, spxb,dpsiidrhoch
      use arrays, only: ptfprimch,psiichq,lambda,rhoch
      implicit double precision(a-h, o-z)

      real *8 psicon(*),R0(*),iasp(*),kappa(*),delta(*)
      real *8 chcoeffR(ncon),chcoeffk(ncon),chcoeffd(ncon)
      real *8 dR0(ncon),dkappa(ncon),ddelta(ncon)

      real *8 chcoeff0(ncon),dptfprimch(ncon)
      real *8 Btheta2(ncon),betatf2(ncon),tempIb(ncon),Ib(*)
      integer *4 sol_unit, ncon
 
      call chftransq(chcoeffR, R0, ncon, cftmq)
      call chftransq(chcoeffk, kappa, ncon, cftmq)
      call chftransq(chcoeffd, delta, ncon, cftmq)
      do i=1,ncon
         xch=1.0d0-psiichq(i)
         call chderiv(xch,ncon,chcoeffR,dR0(i))
         call chderiv(xch,ncon,chcoeffk,dkappa(i))
         call chderiv(xch,ncon,chcoeffd,ddelta(i))
      end do
      dR0(1:ncon)=-dR0(1:ncon)
      dkappa(1:ncon)=-dkappa(1:ncon)
      ddelta(1:ncon)=-ddelta(1:ncon)

      R00=R0(ncon)
      rmin=R00*reps
      tempIb(1:ncon)=0.0d0
      do i=1,ncon
c         do k=1,ncon
c            tempIb(i)=tempIb(i)+spbx((k-1)*ncon+i)*Ib(k)
c         end do
         tempIb(i)=ptfprimch(i)*dpsiidrhoch(i)/lambda
      end do
      call chftransq(chcoeff0, tempIb, ncon, cftmq)
      do i=1,ncon
         psin=1.0d0-psiichq(i)
         call chderiv(psin,ncon,chcoeff0,dptfdr2)
         dptfprimch(i)=-dptfdr2/dpsiidrhoch(i)
      end do

      do i=1,ncon
         temp1=0.0d0
         do k=1,ncon
            temp1=temp1+spbx((k-1)*ncon+i)*dptfprimch(k)*rhoch(k)**2
     1           *(R0(k)-R0(ncon))
         end do
         temp1=temp1/2.0d0
         Btheta2(i)=dpsiidrhoch(i)**2/rmin**2/(lambda**2)/R00**2
         temp3=rhoch(i)*R00*Btheta2(i)
         betatf2(i)=temp1/temp3/2.0d0
c         write(*,*) 'temp13',temp1,temp3
      end do

c      write(*,*) 'betatf2',betatf2

      open(sol_unit, file='fitMiller.out')

      write(sol_unit, *) 'psi, R0, dR0/dpsi, r/R0, kappa
     1 , dkappa/dpsi, delta, ddelta/dpsi'
      do k=1, ncon
         write (sol_unit, *) psicon(k), R0(k),dR0(k),iasp(k),
     $        kappa(k), dkappa(k),delta(k),ddelta(k)
      end do  
      close(sol_unit)

      return
      end
      
      subroutine printQSEx(ncon, psicon, qpsiconex, qcon, shat)

      use arrays, only:  nr, ntot, nt2, kcheb, nsub, f0, iptype, fpolsgn
      use arrays, only: nflx, psiflx, lambda,psiB,psiichq,cftmq
      use arrays, only:  npsi,psiEF, psiiEF, qpsiEF
      implicit double precision(a-h, o-z)

      real *8 psicon(*), qcon(*), shat(*),
     1   qpsiconex(*)

      real *8 psiiflx(nflx),qout(nflx),psiout(nflx)
      real *8 sout(nflx),qoutex(nflx), fpolflx(nflx)
      real *8 qoutEF(npsi) ,soutEf(npsi)

      integer *4 sol_unit, ncon
      integer *4 k, ntheta
      character(3) ko, ns
      character(4) nt

      ntheta = nt2
      sol_unit = 10
      write (nt,'(i4.4)') ntheta
      write (ko,'(i3.3)') kcheb
      write (ns,'(i3.3)') nsub

      open(sol_unit, file='qcontour'//nt//'_ko'//ko
     $     //'_nsub'//ns//'.out')

      write(sol_unit, *) 'psi, psii, q(exact), q(out), shat(out)'
      do k=1, ncon
         write (sol_unit, *) psicon(k), psiichq(k), qpsiconex(k),
     $        qcon(k),shat(k)
      end do  
      close(sol_unit)

c      write (*,*) 'nflx=',nflx
      if (nflx.gt.0) then
       
         select case (iptype)    
         case (0)
            do i=1,nflx
               psiiflx(i)=1-psiflx(i)  
               psiout(i)=psiiflx(i)/lambda +psiB
            end do
            call findqprofile(ncon,lambda,psiichq,qcon,cftmq,
     1           nflx,psiiflx,qout,sout)
            fpolflx(1:nflx)=fpolsgn*dabs(F0)
c            write(*, *) 'fpol', fpolflx(1:nflx)
            call findq_solovev(nflx,psiout,fpolflx,qoutex)
            open(sol_unit, file='qcontour2_'//nt//'_ko'//ko
     $           //'_nsub'//ns//'.out')
            write(*, *) 'psi,qout(ex),qout(out),sout(out)'
            do k=1, nflx
               write (sol_unit, *) psiout(k), 
     $              qoutex(k),qout(k),sout(k)
            end do
            close(sol_unit)
         case (3)
            call findqprofile(ncon,lambda,psiichq,qcon,cftmq,
     1        npsi,psiiEF,qoutEF,soutEF)
            open(sol_unit, file='qcontour2_'//nt//'_ko'//ko
     $           //'_nsub'//ns//'.out')
            write(sol_unit, *) 'psi,qout(EFIT), qout(out), sout(out)'
            do k=1, npsi
               write (sol_unit, *) psiEF(k),qpsiEF(k),
     $              qoutEF(k),soutEF(k)
            end do
            close(sol_unit)
         end select
      end if
      return
      end

      subroutine printSol(psi, psiR, psiZ)

      use arrays, only: Rt1,Zt1,nt1
      use arrays, only:  nr, ntot, nt2, Rt3, Zt3
      use arrays, only:  rnd, tnd, u, ur, uth, urr, urt, utt
      implicit double precision(a-h, o-z)

      real *8 x(ntot)
      real *8 y(ntot)
      real *8 psi(*), psiR(*), psiZ(*)

      integer *4 sol_unit
      integer *4 i,j, inext, ntheta

      ntheta = nt2
      sol_unit = 10
c      open(sol_unit, file='solu_orig.out') !, status='new')
c      do j=1, ntheta
c         do i=1, nr
c            inext = i + (j-1)*nr
c            write (sol_unit, *) Rt3(inext), 
c     $           Zt3(inext), u(inext), ux(inext),uy(inext),
c     $           uxx(inext), uxy(inext),uyy(inext)
c         end do
c      end do
c      close(sol_unit)

      open(sol_unit, file='sol_orig.out') !, status='new')
      do j=1, ntheta
         do i=1, nr
            inext = i + (j-1)*nr
            write (sol_unit, *) Rt3(inext), 
     $           Zt3(inext), psi(inext)
         end do
      end do
      close(sol_unit)

      open(sol_unit, file='sol_origder.out') !, status='new')
      do j=1, ntheta
         do i=1, nr
            inext = i + (j-1)*nr
            write (sol_unit, *) Rt3(inext), 
     $           Zt3(inext), psiR(inext), 
     $           psiZ(inext)
         end do
      end do
      close(sol_unit)

      iw=41
      itype2=2
      itype3=3
      call quaplot2(iw,Rt1,Zt1,nt1,itype3,Rt3,Zt3,nr*ntheta,itype2,
     1       'image of uniform polar grid under inverse map*')

      return
      end

      subroutine printContour(ncon, icon, rndcon, tcon, psicon, Rcon,
     1     Zcon, psiRcon, psiZcon)

      use arrays, only:  nr, ntot, nt2, kcheb, nsub
      implicit double precision(a-h, o-z)

      real *8 rndcon(*), tcon(*), psicon(*), Rcon(*), Zcon(*),
     1   psiRcon(*), psiZcon(*)

      integer *4 sol_unit, ncon, icon(*)
      integer *4 i,j, inext, ntheta, istart, iend, nin

      ntheta = nt2
      sol_unit = 10

      open(sol_unit, file='contour_disk.out') !, status='new')
      itotc = 0
      do k=1, ncon
         do j=1, ntheta*ksamp2
               itotc=itotc+1
               write (sol_unit,*) rndcon(itotc), 
     $           tcon(itotc), psicon(k)
         end do
      end do
      close(sol_unit)

      open(sol_unit, file='contour_orig.out') !, status='new')
      itotc = 0
      do k=1, ncon
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         do j=1, nin !ntheta*ksamp2
               itotc=itotc+1
               write (sol_unit, *) Rcon(itotc), 
     $           Zcon(itotc), psicon(k),
     $           psiRcon(itotc)
         end do
      end do
      close(sol_unit)

      return
      end

      subroutine printQS(ncon, psicon,  qcon, shat)

      use arrays, only:  nr, ntot, nt2, kcheb, nsub
      use arrays, only: nflx, psiflx, lambda,psiB,psiichq,cftmq

      implicit double precision(a-h, o-z)

      real *8 psicon(*), qcon(*), shat(*)
      real *8 psiiflx(nflx),psiout(nflx), qout(nflx), sout(nflx)

      integer *4 sol_unit, ncon
      integer *4 k, ntheta
      character(3) ko, ns
      character(4) nt

      ntheta = nt2
      sol_unit = 10
      write (nt,'(i4.4)') ntheta
      write (ko,'(i3.3)') kcheb
      write (ns,'(i3.3)') nsub

      open(sol_unit, file='qcontour'//nt//'_ko'//ko
     $     //'_nsub'//ns//'.out')

      write(sol_unit, *) 'psi,psii,  q(out), shat(out)'
      do k=1, ncon
         write (sol_unit, *) psicon(k),psiichq(k), 
     $        qcon(k),shat(k)
      end do
      close(sol_unit)
 
      if (nflx.gt.0) then
         do i=1,nflx
            psiiflx(i)=1-psiflx(i)  
            psiout(i)=psiiflx(i)/lambda +psiB
         end do
         call findqprofile(ncon,lambda,psiichq,qcon,cftmq,
     1        nflx,psiiflx,qout,sout)

         open(sol_unit, file='qcontour2_'//nt//'_ko'//ko
     $     //'_nsub'//ns//'.out')

         write(sol_unit, *) 'psi, qout(out), sout(out)'
         do k=1, nflx
            write (sol_unit, *) psiout(k), 
     $        qout(k),sout(k)
         end do
         close(sol_unit)
      end if
  
      return
      end

    
      subroutine printEfitG (npsi, nchq, qcon, acmcur0)

      use arrays, only:  cftmq, pprimch,ffprimch, presch,fpolch
      use arrays, only:  psi0, psiB, lambda, iscale
      use arrays, only:  nin, Rin, Zin, R0, Z0
      use arrays, only:  nt1, Rt1, Zt1, T1, dRt1, dZt1, dT1, w1
      use arrays, only:  Rt3, Zt3, ntot
      use arrays, only:  B0vac, Rmid, Rmaxis, Zmaxis
      use arrays, only:  iptype, csol,reps,rkappa,delta,d1,d2,d3
      use arrays, only:  iprintEFIT, T1_mil, T1_mil2, kLag, epsLag
      
      implicit double precision(a-h, o-z)
      character*20 case_out(6), file_out

      integer *4 kvtor,nmaxx 

      real *8 pprimchtmp(nchq), chcoeff(nchq), qcon(*)
      real *8 REFo(npsi), ZEFo(npsi), PSIEFo(npsi), psiiEFo(npsi)
      real *8 pprimEFo(npsi), ffprimEFo(npsi), fpolEFo(npsi)
      real *8 presEFo(npsi), qpsiEFo(npsi), rhovnEFo(npsi)
      real *8 rndEFo(npsi*npsi), tndEFo(npsi*npsi), psiex(npsi*npsi)
      real *8 psiRZ(npsi*npsi)
      real *8 rbbbs(nin+1), zbbbs(nin+1)
      real *8 rlim(nin+1), zlim(nin+1), dsdt(nt1)

      complex *16 ima, zc1(nt1), dzdtc1(nt1)

      ima = (0.0d0,1.0d0)
      pi=4.0d0*datan(1.0d0)
      
      neqdsk=21
      write(*,*) 'EFITin0', npsi,nin,nchq
c      write(*,*) 'EFITin',pprimch(1:nchq),ffprimch(1:nchq)
c      write(*,*) 'EFITin',fpolch(1:nchq), presch(1:nchq)
      nw=npsi
      nh=npsi

      Rleft = 1000.0
      Rright = -1000.0
      Zbott = 1000.0
      Ztop = -1000.0
      do i=1, nt1
         if (Rt1(i).lt.Rleft) then
            Rleft=Rt1(i)
         end if
         if (Rt1(i).gt.Rright) then
            Rright=Rt1(i)
         end if
         if (Zt1(i).lt.Zbott) then
            Zbott=Zt1(i)
         end if
         if (Zt1(i).gt.Ztop) then
            Ztop=Zt1(i)
         end if
      end do
      rdim=(Rright-Rleft)*1.4d0
      delR=(Rright-Rleft)*0.2d0
      Rleft=Rleft-delR
      Rright=Rright-delR

      Zmid=Z0
      zdim=(Ztop-Zbott)*1.4d0
      rzdimmin=min((Rright-Rleft)/2.0,(Ztop-Zbott)/2.0)
      do j=1, nh
         ZEFo(j) = Zmid-Zdim/2.0d0+Zdim*(j-1)/(nh-1)
      end do

      do i=1, npsi
         PSIEFo(i)= PSI0+(PSIB-PSI0)*(i-1)/(npsi-1)
         psiiEFo(i)=(psiEFo(i)-psiB)/(psi0-psiB)
         REFo(i) = rleft+rdim*(i-1)/(nw-1)
      end do

      if (iscale.ge.5) then
         pprimchtmp(1:nchq)=pprimch(1:nchq)*(-lambda) !from dp/dpsii to dp/dpsi where psii=0 at the center and psii=1 at the edge
      else
         pprimchtmp(1:nchq)=pprimch(1:nchq)
      end if

      call chftransq(chcoeff, pprimchtmp, nchq, cftmq)
      do i=1,npsi
         call chfit(1.0d0-psiiEFo(i), nchq, chcoeff, pprimEFo(i))
      end do
      call chftransq(chcoeff, ffprimch, nchq, cftmq)
      do i=1,npsi
         call chfit(1.0d0-psiiEFo(i), nchq, chcoeff, ffprimEFo(i))
      end do
      call chftransq(chcoeff, fpolch, nchq, cftmq)
      do i=1,npsi
         call chfit(1.0d0-psiiEFo(i), nchq, chcoeff, fpolEFo(i))
      end do
      call chftransq(chcoeff, presch, nchq, cftmq)
      do i=1,npsi
         call chfit(1.0d0-psiiEFo(i), nchq, chcoeff, presEFo(i))
      end do
      call chftransq(chcoeff, qcon, nchq, cftmq)
      do i=1,npsi
         call chfit(1.0d0-psiiEFo(i), nchq, chcoeff, qpsiEFo(i))
      end do
      write(*,*) 'EFITdone1'

      idum = 1
      simag = psi0
      sibry = psib
      bcentr = B0vac*(-1.0d0) !Cmod sign
      Rcentr=Rmid
      current=acmcur0*(-1.0d0) ! assumed the opposite current direction to B field
      xdum = 0.0d0
      
      zc1(1:nt1)=Rt1(1:nt1)+Zt1(1:nt1)*ima
      dsdt(1:nt1)=dsqrt(dRt1(1:nt1)**2+dZt1(1:nt1)**2)
     
      dzdtc1(1:nt1)=(dRt1(1:nt1)+dZt1(1:nt1)*ima)/dsdt(1:nt1)
c      write(*,*) 'error',nt1, dt1, rzdimmin,iexterior
c      write(*,*) 'error',zc1(1:nt1)
c      write(*,*) 'error3', dzdtc1(1:nt1), w1(1:nt1)
      errSolomax=0.0d0
      do i=1, nw
         do j=1, nh        
 1          ij=(i-1)*nh+j

            call findpsioffaxis(nt1, dt1, rzdimmin, zc1, dzdtc1, w1, 
     1           REFo(i), ZEFo(j), iexterior, rndEfo(ij), 
     2           tndEfo(ij), psirz(ij))

            write(*,*) 'iexterior',iprintEFIT,iexterior,psirz(ij)
            if (iexterior.eq.0) then
               psirz(ij)=psirz(ij)/lambda+psiB
            else
               if (iprintEFIT.eq.1) then
                  psirz(ij)=(-1.0d-3)/lambda+psiB !sign(1.0d0,lambda)
               else if (iprintEFIT.eq.2) then !use the artificial contours outside last closed flux surface
                  
                  tRZ =atan2(ZEFo(j)-Z0,REFo(i)-R0)
                  rRz =sqrt((ZEFo(j)-Z0)**2+(REFo(i)-R0)**2)
                  if (tRZ.lt.0) then 
                     tRZ=2.0d0*pi+tRZ
                  end if

                  !write(*,*) 't1_mil',t1_mil(1:nt1+1)
                  !write(*,*) 't1_mil2',t1_mil2(1:nt1+1)
                  nout=1
                  call IntLag1np(nt1+1,t1_mil2, t1_mil, kLag
     1                 ,nout,tRZ, tRZ2,  dtRZ2, epsLag)
                  
                  write(*,*) 'tRZ,tRZ2',tRZ,tRZ2
                  
                  alpha=dasin(delta)
                  rlast=R0*reps*sqrt(dcos(tRZ2+alpha*dsin(tRZ2))**2+
     1                 (rkappa*dsin(tRZ2))**2)

                  write(*,*) 'R,Z',REFo(i), ZEFo(j),tRZ,rRZ, rlast
                  if (rRz.lt.rlast) then
                      psirz(ij)=psiB !(1.0d0-(rRz/rlast))/lambda+psiB
           write(*,*) 'wrong in calculating radius',rRz,rlast,psirz(ij) 
                  else
                     psirz(ij)=(1.0d0-(rRz/rlast))/lambda+psiB
             write(*,*) 'good in calculating radius',rRz,rlast,psirz(ij) 
                     
                  end if

                  
               end if
            end if
            if ((iptype.eq.0).and.(iexterior.eq.0)) then
               call Solovev(REFo(i),ZEFo(j),csol,d1,d2,d3,
     1           psiex(ij),psiRex,psiZex,psiRRex,psiRZex, 
     1           psiZZex,fsolex)
               if (dabs(psirz(ij)-psiex(ij)).gt.errSolomax) then
                  errSolomax=dabs(psirz(ij)-psiex(ij))
                  write(*,*) 'error EFIT1 output',errSolomax
c               write(*,*) 'error psi1 output',psirz(i,j),psiex(i,j)
               end if
            end if
         end do
      end do

      
      write(*,*) 'EFITdone2'
      nbbbs = nin+1
      rbbbs(1:nin)=Rin(1:nin)
      zbbbs(1:nin)=Zin(1:nin)
      rbbbs(nin+1)=Rin(1)
      zbbbs(nin+1)=Zin(1)

      limitr = nin+1
      rlim(1:nin)=Rin(1:nin)
      zlim(1:nin)=Zin(1:nin)
      rlim(nin+1)=Rin(1)
      zlim(nin+1)=Zin(1)
      kvtor = 0
      rvtor = 0.0d0
      nmass = 0
      rhovnEFo =0.0d0

      write(*,*) 'writing Efit output g file'
      write(*,*) 'EFITdone3'

      case_out(1)='ECOM'
      case_out(2)='results'
      case_out(3)='results'
      case_out(4)='results'
      case_out(5)='results'
      case_out(6)='results'

      open(neqdsk, file='gfromECOM', status='unknown')
      write (neqdsk,2000) (case_out(i),i=1,6),idum,nw,nh
      write (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid
      write (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr
      write (neqdsk,2020) current,simag,xdum,rmaxis,xdum
      write (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
      write (neqdsk,2020) (fpolEFo(i),i=1,nw)
      write (neqdsk,2020) (presEFo(i),i=1,nw)
      write (neqdsk,2020) (ffprimEFo(i),i=1,nw)
      write (neqdsk,2020) (pprimEFo(i),i=1,nw)
      write (neqdsk,2020) ((psirz((i-1)*nh+j),i=1,nw),j=1,nh)
      write (neqdsk,2020) (qpsiEFo(i),i=1,nw)
      write (neqdsk,2022) nbbbs,limitr
      write (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
      write (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)
      write (neqdsk,2024) kvtor,rvtor,nmass
c      if (kvtor.gt.0) then
c         write (neqdsk,2020) (pressw(i),i=1,nw)
c         write (neqdsk,2020) (pwprim(i),i=1,nw)
c      endif
c      if (nmass.gt.0) then
c         write (neqdsk,2020) (dmion(i),i=1,nw)
c      endif
      write (neqdsk,2020) (rhovnEFo(i),i=1,nw)
c     
 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)
 2024 format (i5,e16.9,i5)
      close(neqdsk)
 
  
      return
      end
      

      subroutine checkerr(psiRex, psiZex, psiRRex, psiRZex, psiZZex,
     1     psi, psiR, psiZ, psiRR, psiRZ, psiZZ, ncon,psicon, psiRcon, 
     2     psiZcon, psiconex, psiRconex, psiZconex, qcon, 
     3     qpsiconex, icon)

      use arrays, only:  nr, ntot, nt2, psiex, kcheb, nsub  

      implicit double precision(a-h, o-z)

      real *8 psiRex(*), psiZex(*), psiRRex(*)
      real *8 psiRZex(*),psiZZex(*)

      real *8 psi(*), psiR(*), psiZ(*)
      real *8 psiRR(*),psiRZ(*),psiZZ(*)

      real *8 psicon(*), psiRcon(*), psiZcon(*)
      real *8 psiconex(*),psiRconex(*), psiZconex(*)
      
      real *8 qcon(*),qpsiconex(*)
     
      integer *4 icon(*)
      integer *4 i,j, inext, ntheta
      character(3) ko, ns
      character(4) nt

      errmax = 0.0d0
      errRmax = 0.0d0
      errZmax = 0.0d0
      errRRmax = 0.0d0
      errRZmax = 0.0d0
      errZZmax = 0.0d0
      errconmax = 0.0d0
      errpsiRconmax = 0.0d0
      errpsiZconmax = 0.0d0
      errqconmax = 0.0d0

      err = 0.0d0
      errR = 0.0d0
      errZ = 0.0d0
      errRR = 0.0d0
      errZZ = 0.0d0
      errRZ = 0.0d0
      errcon = 0.0d0
      errpsiRcon = 0.0d0
      errpsiZcon = 0.0d0
      errqcon = 0.0d0

      sol = 0.0d0
      solR = 0.0d0
      solZ = 0.0d0
      solRR = 0.0d0
      solZZ = 0.0d0
      solRZ = 0.0d0
      solcon = 0.0d0
      solpsiRcon = 0.0d0
      solpsiZcon = 0.0d0
      solqcon = 0.0d0
 
      ntheta = nt2
      do j = 1,ntheta
      do i = 1,nr
         inext = i + (j-1)*nr
         errtmp = abs(psiex(inext)-psi(inext))
         errRtmp = abs(psiRex(inext)-psiR(inext))
         errZtmp = abs(psiZex(inext)-psiZ(inext))
         errRRtmp = abs(psiRRex(inext)-psiRR(inext))
         errRZtmp = abs(psiRZex(inext)-psiRZ(inext))
         errZZtmp = abs(psiZZex(inext)-psiZZ(inext))

         if (errmax.lt.errtmp) then
            errmax= errtmp
         end if
         if (errRmax.lt.errRtmp) then
            errRmax= errRtmp
         end if
         if (errZmax.lt.errZtmp) then
            errZmax= errZtmp
         end if
         if (errRRmax.lt.errRRtmp) then
            errRRmax= errRRtmp
         end if
         if (errRZmax.lt.errRZtmp) then
            errRZmax= errRZtmp
         end if
         if (errZmax.lt.errZtmp) then
            errZZmax= errZZtmp
         end if
         err = err + errtmp**2
         errR = errR + errRtmp**2
         errZ = errZ + errZtmp**2
         errRR = errRR + errRRtmp**2
         errRZ = errRZ + errRZtmp**2
         errZZ = errZZ + errZZtmp**2
         
         sol = sol + abs(psiex(inext))**2
         solR = solR + abs(psiRex(inext))**2
         solZ = solZ + abs(psiZex(inext))**2
         solRR = solRR + abs(psiRRex(inext))**2
         solRZ = solRZ + abs(psiRZex(inext))**2
         solZZ = solZZ + abs(psiZZex(inext))**2
      enddo
      enddo
      relerr = sqrt(err/sol)
      relerrR = sqrt(errR/solR)
      relerrZ = sqrt(errZ/solZ)
      err = sqrt(err)
      errR = sqrt(errR)
      errZ = sqrt(errZ)
      
      relerrRR = sqrt(errRR/solRR)
      relerrRZ = sqrt(errRZ/solRZ)
      relerrZZ = sqrt(errZZ/solZZ)
      errRR = sqrt(errRR)
      errRZ = sqrt(errRZ)
      errZZ = sqrt(errZZ)

      itotc = 0

      do k=1, ncon-1
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         do j=1, nin
               itotc=itotc+1
               errcontmp = abs(psiconex(itotc)-psicon(k))
               errpsiRcontmp = abs(psiRconex(itotc)-psiRcon(itotc))
               errpsiZcontmp = abs(psiZconex(itotc)-psiZcon(itotc))
           

               if (errconmax.lt.errcontmp) then
                  errconmax= errcontmp
                  errpsiRconmax= errpsiRcontmp
                  errpsiZconmax= errpsiZcontmp  
               end if
               errcon = errcon + errcontmp**2
               solcon = solcon + abs(psiconex(itotc))**2
               errpsiRcon = errpsiRcon + errpsiRcontmp**2
               solpsiRcon = solpsiRcon + abs(psiRconex(itotc))**2
               errpsiZcon = errpsiZcon + errpsiZcontmp**2
               solpsiZcon = solpsiZcon + abs(psiZconex(itotc))**2
c            end do
         end do
         errqtmp = abs(qpsiconex(k)-qcon(k))
         if (errqconmax.lt.errqtmp) then
            errqconmax= errqtmp
         end if
         errqcon = errqcon + errqtmp**2
         solqcon = solqcon + abs(qpsiconex(k))**2
      end do
      relerrcon = sqrt(errcon/solcon)
      errcon = sqrt(errcon)
      relerrpsiRcon = sqrt(errpsiRcon/solpsiRcon)
      errpsiRcon = sqrt(errpsiRcon)
      relerrpsiZcon = sqrt(errpsiZcon/solpsiZcon)
      errpsiZcon = sqrt(errpsiZcon)
      errqcon = sqrt(errqcon)
      relerrqcon = sqrt(errqcon/solqcon)

      write (nt,'(i4.4)') ntheta
      write (ko,'(i3.3)') kcheb
      write (ns,'(i3.3)') nsub

      open(7, file='err_nt'//nt//'_ko'//ko
     $     //'_nsub'//ns//'.out') !, status='new')
      write(7,*) ' nt1, klag, ntheta, nr,  kcheb, nsub ', 
     $     nt1, klag, ntheta, nr, kcheb, nsub
      write(7,*) ' Solution, err, rel err, max(err),',err,','
     $     ,relerr,',',errmax,','
      write(7,*) 'dpsi/dR,err, rel err, max(errR),',errR,',',relerrR,','
     $     ,errRmax,','
      write(7,*) 'dpsi/dZ,err, rel err, max(errZ),',errZ,',',relerrZ,','
     $     ,errZmax,','
      write(7,*) 'd2psi/dRR,  err, rel err, max(errRR),'
     1     ,errRR,',',relerrRR,',',errRRmax,','
      write(7,*) 'd2psi/dRdZ,  err, rel err, max(errRZ),'
     1     ,errRZ,',',relerrRZ,',',errRZmax,','
      write(7,*) 'd2psi/dZdZ,  err, rel err, max(errZZ),'
     1     ,errZZ,',',relerrZZ,',',errZZmax,','
      write(7,*) 'contours,  err, rel err, max(err),'
     1     ,errcon,',',relerrcon,',',errconmax,','
      write(7,*) 'contours dpsi/dR,  err, rel err, max(errR),'
     1     ,errpsiRcon,',',relerrpsiRcon,',',errpsiRconmax,','
      write(7,*) 'contours dpsi/dZ,  err, rel err, max(errR),'
     1     ,errpsiZcon,',',relerrpsiZcon,',',errpsiZconmax,','

      write(7,*) 'contours q, err, rel err, max(errR),'
     1     ,errqcon,',',relerrqcon,',',errqconmax,','
      close(7)
      write(*,*) 'finished errorcheck'
     
      return
      end

      subroutine sortcontour(ncon,psicon,zk,dzdw2k,rndcon
     1     ,tcon,urcon,utcon,Rcon,Zcon,icon,Ravg,Zavg
     2     ,tmaxp,psiRcon,psiZcon,psiconex,psiRconex,psiZconex
     3     ,Rmido, Rmidi)

      use arrays, only:  csol,d1,d2,d3,iptype, Rmaxis, Zmaxis
      use arrays, only:  Rmaxpsi,Zmaxpsi, ksamp2, nt2
      use arrays, only:  epsLag2, psiB
      use arrays, only: ici, icr, ic_max, kLag, isymud
      use arrays, only:  epslogzk, epslogzwk, epslogzwwk 

      implicit real*8 (a-h,o-z)
 
      real *8 psiconex(*),psiRconex(*),psiZconex(*)
      real *8 psicon(*),psiRcon(*),psiZcon(*)
      real *8 rndcon(*),tcon(*),urcon(*),utcon(*),tmaxp(*)
      real *8 Rcon(*), Zcon(*), Ravg(*), Zavg(*)
      real *8 Rmido(*), Rmidi(*)
      complex *16 zk(*),dzdw2k(*)
      integer *4 ncross,icon(*), im(ksamp2), jlag(2*ksamp2+2*klag)
      integer *4 mlag(2*ksamp2+2*klag),jarr(ksamp2)

      real *8 Rconm(2*nt2*ksamp2), Zconm(2*nt2*ksamp2)
      real *8 psiRconm(2*nt2*ksamp2),psiZconm(2*nt2*ksamp2)
      real *8 tmaxpm(2*nt2*ksamp2), tlag(2*ksamp2+2*klag)
      real *8 wbary(ksamp2),mlagtmp(ksamp2+klag)
      real *8 Rcontemp(2*nt2*ksamp2),Zcontemp(2*nt2*ksamp2)
      real *8 tmaxptemp(2*nt2*ksamp2),wlag(klag)
      real *8 tlagtmp(ksamp2+klag),Ztmp(ksamp2),Rtmp(ksamp2)
      real *8 psiRcontemp(2*nt2*ksamp2),psiZcontemp(2*nt2*ksamp2)
      real *8 min2pi, mindth 

      real *8 rndconm(2*nt2*ksamp2),tconm(2*nt2*ksamp2)
      real *8 urconm(2*nt2*ksamp2), utconm(2*nt2*ksamp2)


      real *8 tmaxpmid(2),Rmid,dRdt

      complex *16 rndconc(2*nt2)
      complex *16 urconc(2*nt2),utconc(2*nt2)
      complex *16 rndconmc(2*nt2*ksamp2)
      complex *16 wm(10*nt2*ksamp2)
      complex *16 urconmc(2*nt2*ksamp2), utconmc(2*nt2*ksamp2)
   
      complex *16 wcon,ima,c1,cint,dzdwcon,dzdwwcon,dwdzzcon

      integer *4 isort(2*nt2*ksamp2), imindth(2*ksamp2)
      integer *4 isortstart,itotc,isortend,nsort
      integer i,j,k,imaxu,jmaxt, irinsub, nr, i2pi, isame
      

      save 

      kLag=8
      ima = (0.0d0,1.0d0)
      c1 =  (1.0d0,0.0d0)
      pi = 4.0d0*datan(1.0d0)

      itotc = 0
      itotc2 = 0
      do k=1, ncon
        
         Ravg(k) = 0.0d0
         Zavg(k) = 0.0d0
         m=nt2*ksamp2
c         write(*,*) 'nt2,ksamp2,m',nt2,ksamp2,m
         iorgstart=(k-1)*nt2+1
         rndconc(1:nt2)=rndcon(iorgstart:iorgstart+nt2-1)*c1
         urconc(1:nt2)=urcon(iorgstart:iorgstart+nt2-1)*c1
         utconc(1:nt2)=utcon(iorgstart:iorgstart+nt2-1)*c1
         if (isymud.eq.1) then     !oversampling using FFT pading for backward conformal mapping
            call interpffte(nt2,rndconc,m,rndconmc,wm)
         else
            call interpfft(nt2,rndconc,m,rndconmc,wm)
         end if
         rndconm(1:m)=(rndconmc(1:m))
         do j=1,m
            if(mod((j-1),ksamp2).eq.0) then
               tconm(j)=tcon(iorgstart+int((j-1)/ksamp2))
            else
               tconm(j)=tcon(iorgstart+int((j-1)/ksamp2))+
     1              (2.0d0*pi/m)*mod((j-1),ksamp2)
            end if
c            tcon(1:m)=tcon(iorgstart)+(2.0d0*pi/m)*(j-1)
         end do
         
        
         if (isymud.eq.1) then
            call interpffte(nt2,urconc,m,urconmc,wm)
         else
            call interpfft(nt2,urconc,m,urconmc,wm)
         end if
         call system_clock(ict, icr, ic_max)
         if (ict.gt.ici) then
            tsec=real(ict-ici)/real(icr)
         else
            tsec=real(ict-ici+ic_max)/real(icr)
         end if

         if (isymud.eq.1) then
            call interpffto(nt2,utconc,m,utconmc,wm)
         else
            call interpfft(nt2,utconc,m,utconmc,wm)
         end if
         urconm(1:m)=(urconmc(1:m))
         utconm(1:m)=(utconmc(1:m))

    
         if (psicon(k).lt.psicon(1)*0.2) then
            ksamp3= ksamp2
         else
            ksamp3= ksamp2
         end if

 

         idsamp3=m/ksamp3
         do j=1, ksamp3     !initial guess of the uniform theta grid in plasma domain
            im(j)=(j-1)*idsamp3+1
            wcon=rndconm(im(j))*exp(ima*tconm(im(j)))      

            nlimz=epslogzk/log(rndconm(im(j)))+1
            if (nlimz.gt.nt2) then
               nlimz=nt2
            end if
            call fft_cauchy(nlimz,wcon,zk,cint)
            Rconm(j)=cint
            Zconm(j)=-ima*cint
            Ravg(k)=Ravg(k)+Rconm(j)
            Zavg(k)=Zavg(k)+Zconm(j)
         end do
         Ravg(k)=Ravg(k)/ksamp3
         Zavg(k)=Zavg(k)/ksamp3

c         write(*,*) "(Ravg,Zavg)=",Ravg(k),Zavg(k)
c         write(*,*) "(Rmax,Zmax)=",Rmaxpsi,Zmaxpsi
         tmin=100.0
         do j=1, ksamp3
c            tmaxpm(j)=datan2((Zconm(j)-Zavg(k)),
c     1           (Rconm(j)-Ravg(k)))
            tmaxpm(j)=datan2((Zconm(j)-Zmaxis),
     1           (Rconm(j)-Rmaxis))
            if (tmaxpm(j).lt.0.0d0) then 
               tmaxpm(j)=2.0d0*pi+tmaxpm(j)
            end if 
            if (tmaxpm(j).lt.tmin) then
               tmin=tmaxpm(j)
               jtmin=j
            end if
         end do
         do j=1,ksamp3
            jarr(j)=j
         end do
         jlag(1:klag)=jarr(ksamp3-klag+1:ksamp3)
         mlag(1:klag)=im(ksamp3-klag+1:ksamp3)-m      
         tlag(1:klag)=tmaxpm(ksamp3-klag+1:ksamp3)-2.0d0*pi
         
         jlag(klag+1:klag+ksamp3)=jarr(1:ksamp3)
         mlag(klag+1:klag+ksamp3)=im(1:ksamp3)      
         tlag(klag+1:klag+jtmin-1)=tmaxpm(1:jtmin-1)-2.0d0*pi
         tlag(klag+jtmin:klag+ksamp3)=tmaxpm(jtmin:ksamp3)

         jlag(klag+ksamp3+1:klag+ksamp3+jtmin-1)=jarr(1:jtmin-1)
         mlag(klag+ksamp3+1:klag+ksamp3+jtmin-1)=im(1:jtmin-1)+m      
         tlag(klag+ksamp3+1:klag+ksamp3+jtmin-1)=tmaxpm(1:jtmin-1)

         jlag(klag+ksamp3+jtmin:2*klag+ksamp3+jtmin-1)=
     1        jarr(jtmin:jtmin+klag-1)
         mlag(klag+ksamp3+jtmin:2*klag+ksamp3+jtmin-1)=
     1        im(jtmin:jtmin+klag-1)+m      
         tlag(klag+ksamp3+jtmin:2*klag+ksamp3+jtmin-1)=
     1        tmaxpm(jtmin:jtmin+klag)+2.0d0*pi
         isortstart = itotc +1
         icon(k)=isortstart
         dth=0.0d0 !2*pi/ksamp3*0.001d0
         !asuume dth=0 ---> thksamp3(j=1)=0.0d0
         istart=jtmin-1+klag

         klagl=floor(klag/2.0)
         klagr=klag-klagl

         eps13=1.0d-13
         do j=1, ksamp3
c            isort(i)=isortstart+i-1
            itotc=itotc+1
            thksamp3= 2.0d0*pi/ksamp3*(j-1)+dth
            ifind = 0
            jj=0
                        
            do while ((ifind.eq.0) .and. (jj.lt.1)) !find the theta grid which is equispaced as much as possible (not perfectly) 
               jj=jj+1
               if (jj.eq.1) then
                  do while (thksamp3.gt.tlag(istart)) 
                     istart=istart+1
                  end do
                  tlagtmp(1:klag)=tlag(istart-klagl:istart+klagr)
                  mlagtmp(1:klag)=mlag(istart-klagl:istart+klagr)  !index of m

                  call IntLagAll(klag, tlagtmp, mlagtmp,
     1              thksamp3,y1,wbary,eps13,isame)
                  if (isame.ne.0) then
                     ifind=1
                     itry2=modulo(nint(y1-1),m)+1
                     ttry2=tlagtmp(isame)
                     Rtry=Rconm(jlag(istart-klagl-1+isame))
                     Ztry=Zconm(jlag(istart-klagl-1+isame))
                     exit
                  else                                 
                     call CheckSameInt(klag,mlag(istart-klagl)
     1                   ,nint(y1),isame2)
                     if (isame2.ne.0) then
                        ifind=1
                        itry2=modulo(nint(y1-1),m)+1
                        ttry2=tlagtmp(isame2)
                        Rtry=Rconm(jlag(istart-klagl-1+isame2))
                        Ztry=Zconm(jlag(istart-klagl-1+isame2))
                        exit 
                     else
                        itry1=modulo(nint(y1-1),m)+1
                        wcon=rndconm(itry1)*exp(ima*tconm(itry1))
                        nlimz=epslogzk/log(rndconm(itry1))+1
                        if (nlimz.gt.nt2) then
                           nlimz=nt2
                        end if
                        call fft_cauchy(nlimz,wcon,zk,cint)
                        Rtry=cint
                        Ztry=-ima*cint
c                        ttry1=datan2((Ztry-Zavg(k)),(Rtry-Ravg(k)))
                        ttry1=datan2((Ztry-Zmaxis),(Rtry-Rmaxis))
                        if (ttry1.lt.0.0d0) then 
                           ttry1=2.0d0*pi+ttry1
                        end if 
                     end if
                  end if
               else
                  call CheckSameInt(klag,mlag(istart-klagl)
     1                   ,nint(y2),isame2)
                  if (isame2.ne.0) then
                     ifind=1
                     itry2=modulo(nint(y2-1),m)+1
                     ttry2=tlagtmp(isame2)
                     Rtry=Rconm(jlag(istart-klagl-1+isame2))
                     Ztry=Zconm(jlag(istart-klagl-1+isame2))
                     exit 
                  else
                     itry1=itry2
                     ttry1=ttry2
                     y1=y2
                  end if
               end if
               tlagtmp(klag+jj)=ttry1
               mlagtmp(klag+jj)=nint(y1)
               call IntLagAdd(klag+jj, tlagtmp, mlagtmp,  
     1              thksamp3,y2,wbary,eps13,isame) 
               
               itry2=modulo(nint(y2-1),m)+1
               if (itry1.eq.itry2) then
                  ifind=1
                  ttry2=ttry1
                  exit
               else
                  wcon=rndconm(itry2)*exp(ima*tconm(itry2))
                  nlimz=epslogzk/log(rndconm(itry2))+1
                  if (nlimz.gt.nt2) then
                     nlimz=nt2
                  end if
                  call fft_cauchy(nlimz,wcon,zk,cint)
                  Rtry=cint
                  Ztry=-ima*cint
c                  ttry2=datan2((Ztry-Zavg(k)),(Rtry-Ravg(k)))
                  ttry2=datan2((Ztry-Zmaxis),(Rtry-Rmaxis))
                  if ((ttry2.lt.0.0d0).and.(j.ne.1)) then 
                     ttry2=2.0d0*pi+ttry2
                  end if
               end if
            end do
  
           
            tmaxp(itotc)=ttry2
            if ((j.gt.1).and.(tmaxp(itotc-1).gt.tmaxp(itotc))) then 
               tmaxp(itotc-1)=tmaxp(itotc-1)-2.0d0*pi
            end if
            Rcon(itotc)=Rtry
            Zcon(itotc)=Ztry

            wcon=rndconm(itry2)*exp(ima*tconm(itry2))
            nlimzw=epslogzwk/log(rndconm(itry2))+1
            if (nlimzw.gt.nt2) then
               nlimzw=nt2
            end if
            call fft_cauchy(nlimzw,wcon,dzdw2k,cint)
            dzdwcon = cint
            dwdzr = (1.0d0,0.0d0)/cint
            dwdzi = -ima*(1.0d0,0.0d0)/cint
            
            dct=dcos(tconm(itry2))
            dst=dsin(tconm(itry2))
            uxdisk=urconm(itry2)*dct
     1           +1/rndconm(itry2)*utconm(itry2)*(-dst)
            uydisk=urconm(itry2)*dst
     1           +1/rndconm(itry2)*utconm(itry2)*(dct) 
            psiRcon(itotc)=(psicon(k))/(2.0d0*Rtry)   !scaled psicon (not normalized psii)
     1              +dsqrt(Rtry)*(uxdisk*dwdzr+uydisk*dwdzi)
            psiZcon(itotc)=dsqrt(Rtry)
     1              *(uxdisk*(-dwdzi)+uydisk*dwdzr)

       
            select case (iptype)
            case(0)
               call Solovev(Rcon(itotc),Zcon(itotc),csol,d1,d2,d3,
     1        psiconex(itotc),psiRconex(itotc),psiZconex(itotc)
     2              ,psiconRR,psiconRZ,psiconZZ,fsolcon)

            end select
         end do
   
         call system_clock(ict, icr, ic_max)
         if (ict.gt.ici) then
            tsec=real(ict-ici)/real(icr)
         else
            tsec=real(ict-ici+ic_max)/real(icr)
         end if
c        write(*,*) "CPU time sortcontour4: t(sec)= ",tsec 

         !find outer mid-plane point and inner mid-plane point
         tmaxpmid(1)=tmaxp(1) !(isortstart)
         tmaxpmid(2)=tmaxp(ksamp3/2+1) !isortstart)
c         write(*,*) 'tmaxpmid',tmaxpmid
c         write(*,*) 'tmaxp',tmaxp(1:ksamp3)
         do j=1,klag
            jj=modulo(jtmin-klagl-2+j,m)+1
            wcon=rndconm(jj)*exp(ima*tconm(jj))

            nlimz=epslogzk/log(rndconm(jj))+1
            if (nlimz.gt.nt2) then
               nlimz=nt2
            end if

            call fft_cauchy(nlimz,wcon,zk,cint)
            Rtmp(j)=cint
            Ztmp(j)=-ima*cint
c            tlagtmp(j)=datan2((Ztmp(j)-Zavg(k)),(Rtmp-Ravg(k)))
            tlagtmp(j)=datan2((Ztmp(j)-Zmaxis),(Rtmp(j)-Rmaxis))
         end do
         epsLag22=epsLag2
         call IntLagAll(klag,tlagtmp,Rtmp,
     1           tmaxpmid(1),Rmid,wlag,epsLag22, isame)
c         write(*,*) 'tlagtmp',tlagtmp(1:klag),Rtmp(1:klag),Rmid
         Rmido(k)=Rmid
         
 
c         dZdt= 0.0d0
c         epsLag22=epsLag2
c         do while (dZdt.eq.0.0d0)
c            call IntLag_sameptsAll(klag,tlagtmp,Ztmp,
c     1           tmaxpmid(1),Zmid(1),dZdt,epsLag22)
c            write(*,*) 'epsLag2 is changed to ',epsLag22
c            if (epsLag22.gt.1d-8) then
c               write(*,*) 'no z grid in outer-midplane'
c               if (isymud.eq.1) then
c                  write(*,*) 'contour integrals may not be accurate'
c                  isymud=0
c               end if
c               exit
c            else if (dZdt.eq.0.0d0) then
c               epsLag22=epsLag22*10
c               write(*,*) 'epsLag2 is changed to ',epsLag22
c            end if
c         end do
c         epsLag2=epsLag22
c         dZdtmido(k)=dZdt

         do j=1,klag
            jj=modulo(jtmin-klagl-2+m/2+j,m)+1
c            write(*,*) 'tconma',tconm(jj)-pi
            wcon=rndconm(jj)*exp(ima*tconm(jj))

            nlimz=epslogzk/log(rndconm(jj))+1
            if (nlimz.gt.nt2) then
               nlimz=nt2
            end if

            call fft_cauchy(nlimz,wcon,zk,cint)
            Rtmp(j)=cint
            Ztmp(j)=-ima*cint
            tlagtmp(j)=datan2((Ztmp(j)-Zmaxis),(Rtmp(j)-Rmaxis))
            if (tlagtmp(j).lt.0.0d0) then 
               tlagtmp(j)=2.0d0*pi+tlagtmp(j)
            end if
         end do

         call IntLagAll(klag,tlagtmp,Rtmp,
     1           tmaxpmid(2),Rmid,wlag,epsLag22, isame)

         Rmidi(k)=Rmid
c         dZdt= 0.0d0
c         epsLag22=epsLag2
c         do while (dZdt.eq.0.0d0)
c            call IntLag_sameptsAll(klag,tlagtmp,Ztmp,
c     1           tmaxpmid(2),Zmid(2),dZdt,epsLag22)         
c            if (epsLag22.gt.1d-8) then
c               write(*,*) 'no z grid in inner-midplane'
c               if (isymud.eq.1) then
c                  write(*,*) 'contour integrals may not be accurate'
c                  isymud=0
c               end if
c               exit
c            else if (dZdt.eq.0.0d0) then
c               epsLag22=epsLag22*10
c               write(*,*) 'epsLag2 is changed to ',epsLag22
c            end if
c         end do
c         dZdtmidi(k)=dZdt
c         write(*,*) 'klag',klag,tlagtmp(1:klag),Ztmp(1:klag)
c         write(*,*) 'klag2',tmaxpmid(2), Zmid(2),dZdt

c         write(*,*) "finished sortcontour0 for icon=",k
      end do
      icon(ncon+1)=itotc+1
     
      write(*,*) "finished sortcontour1"
      return
      end

      subroutine sortcontour2der(ncon,psicon,zk,dzdw2k,dzdww2k
     1     ,rndcon,tcon,urcon,utcon,urrcon,urtcon,uttcon
     2     ,Rcon,Zcon,icon,Ravg,Zavg,tmaxp,psiRcon,psiZcon
     3     ,psiRRcon,psiRZcon,psiZZcon,psiconex,psiRconex,psiZconex
     4     ,psiRRconex,psiRZconex,psiZZconex,Rmido, Rmidi)

      use arrays, only:  csol,d1,d2,d3,iptype, Rmaxis, Zmaxis
      use arrays, only:  Rmaxpsi,Zmaxpsi, ksamp2, nt2
      use arrays, only:  epsLag2,psiB,psii
      use arrays, only: ici, icr, ic_max, kLag, isymud
      use arrays, only:  epslogzk, epslogzwk, epslogzwwk 

      implicit real*8 (a-h,o-z)
 
      real *8 psiconex(*),psiRconex(*),psiZconex(*)
      real *8 psicon(*),psiRcon(*),psiZcon(*)
      real *8 psiRRconex(*),psiRZconex(*),psiZZconex(*)
      real *8 psiRRcon(*),psiRZcon(*),psiZZcon(*)
      real *8 urrcon(*),urtcon(*),uttcon(*)

      real *8 rndcon(*),tcon(*),urcon(*),utcon(*),tmaxp(*)
      real *8 Rcon(*), Zcon(*), Ravg(*), Zavg(*)
      real *8 Rmido(*), Rmidi(*)
      complex *16 zk(*),dzdw2k(*),dzdww2k(*)
      integer *4 ncross,icon(*), im(ksamp2), jlag(2*ksamp2+2*klag)
      integer *4 mlag(2*ksamp2+2*klag),jarr(ksamp2)

      real *8 Rconm(2*nt2*ksamp2), Zconm(2*nt2*ksamp2)
      real *8 psiRconm(2*nt2*ksamp2),psiZconm(2*nt2*ksamp2)
      real *8 tmaxpm(2*nt2*ksamp2), tlag(2*ksamp2+2*klag)
      real *8 wbary(ksamp2),mlagtmp(ksamp2+klag)
      real *8 Rcontemp(2*nt2*ksamp2),Zcontemp(2*nt2*ksamp2)
      real *8 tmaxptemp(2*nt2*ksamp2),wlag(klag)
      real *8 tlagtmp(ksamp2+klag),Ztmp(ksamp2),Rtmp(ksamp2)
      real *8 psiRcontemp(2*nt2*ksamp2),psiZcontemp(2*nt2*ksamp2)
      real *8 min2pi, mindth 

      real *8 rndconm(2*nt2*ksamp2),tconm(2*nt2*ksamp2)
      real *8 urconm(2*nt2*ksamp2), utconm(2*nt2*ksamp2)
      real *8 urrconm(2*nt2*ksamp2), urtconm(2*nt2*ksamp2)
      real *8 uttconm(2*nt2*ksamp2)

      real *8 tmaxpmid(2),Rmid,dRdt

      complex *16 rndconc(2*nt2)
      complex *16 urconc(2*nt2),utconc(2*nt2)
      complex *16 urrconc(2*nt2),urtconc(2*nt2),uttconc(2*nt2)
      complex *16 rndconmc(2*nt2*ksamp2)
      complex *16 wm(10*nt2*ksamp2)
      complex *16 urconmc(2*nt2*ksamp2), utconmc(2*nt2*ksamp2)
      complex *16 urrconmc(2*nt2*ksamp2), urtconmc(2*nt2*ksamp2)
      complex *16 uttconmc(2*nt2*ksamp2)

      complex *16 wcon,ima,c1,cint,dzdwcon,dzdwwcon,dwdzzcon

      integer *4 isort(2*nt2*ksamp2), imindth(2*ksamp2)
      integer *4 isortstart,itotc,isortend,nsort
      integer i,j,k,imaxu,jmaxt, irinsub, nr, i2pi, isame
      

      save 

      kLag=8
      ima = (0.0d0,1.0d0)
      c1 =  (1.0d0,0.0d0)
      pi = 4.0d0*datan(1.0d0)

      itotc = 0
      itotc2 = 0
      do k=1, ncon
         !if (psicon(k).lt.psicon(1)*0.2) then
            ksamp3= ksamp2
         !else
         !   ksamp3= ksamp2
         !end if
c            write(*,*) 'ksamp3',ksamp3,ksamp2
c            write(*,*) 'psicon11',psicon(k),psii(1)
         if (psicon(k).gt.psii(1)) then
            Rcon((k-1)*ksamp3+1:k*ksamp3)=Rmaxis
            Zcon((k-1)*ksamp3+1:k*ksamp3)=Zmaxis
            icon(k)=itotc+1
            itotc=itotc+ksamp3
            Ravg(k)=Rmaxis;Zavg(k)=Rmaxis;
            do i=1,ksamp3
               tmaxp((k-1)*ksamp3+i)=(i-1)*2.0*pi/ksamp3
            end do
            psiRcon((k-1)*ksamp3+1:k*ksamp3)=0.0d0
            psiZcon((k-1)*ksamp3+1:k*ksamp3)=0.0d0
            psiRRcon((k-1)*ksamp3+1:k*ksamp3)=0.0d0
            psiRZcon((k-1)*ksamp3+1:k*ksamp3)=0.0d0
            psiZZcon((k-1)*ksamp3+1:k*ksamp3)=0.0d0
            Rmido(k)=Rmaxis
            Rmidi(k)=Rmaxis
        
         else

         
         Ravg(k) = 0.0d0
         Zavg(k) = 0.0d0
         m=nt2*ksamp2
c         write(*,*) 'nt2,ksamp2,m',nt2,ksamp2,m
         iorgstart=(k-1)*nt2+1
         rndconc(1:nt2)=rndcon(iorgstart:iorgstart+nt2-1)*c1
         urconc(1:nt2)=urcon(iorgstart:iorgstart+nt2-1)*c1
         utconc(1:nt2)=utcon(iorgstart:iorgstart+nt2-1)*c1
         urrconc(1:nt2)=urrcon(iorgstart:iorgstart+nt2-1)*c1
         urtconc(1:nt2)=urtcon(iorgstart:iorgstart+nt2-1)*c1
         uttconc(1:nt2)=uttcon(iorgstart:iorgstart+nt2-1)*c1
         if (isymud.eq.1) then     !oversampling using FFT pading for backward conformal mapping
            call interpffte(nt2,rndconc,m,rndconmc,wm)
         else
            call interpfft(nt2,rndconc,m,rndconmc,wm)
         end if
         rndconm(1:m)=(rndconmc(1:m))
         do j=1,m
            if(mod((j-1),ksamp2).eq.0) then
               tconm(j)=tcon(iorgstart+int((j-1)/ksamp2))
            else
               tconm(j)=tcon(iorgstart+int((j-1)/ksamp2))+
     1              (2.0d0*pi/m)*mod((j-1),ksamp2)
            end if
c            tcon(1:m)=tcon(iorgstart)+(2.0d0*pi/m)*(j-1)
         end do
         
        
         if (isymud.eq.1) then
            call interpffte(nt2,urconc,m,urconmc,wm)
         else
            call interpfft(nt2,urconc,m,urconmc,wm)
         end if
         call system_clock(ict, icr, ic_max)
         if (ict.gt.ici) then
            tsec=real(ict-ici)/real(icr)
         else
            tsec=real(ict-ici+ic_max)/real(icr)
         end if

         if (isymud.eq.1) then
            call interpffto(nt2,utconc,m,utconmc,wm)
         else
            call interpfft(nt2,utconc,m,utconmc,wm)
         end if
         urconm(1:m)=(urconmc(1:m))
         utconm(1:m)=(utconmc(1:m))
c         write(*,*) 'psiRin, temp1-0',urconmc(1:m)
         if (isymud.eq.1) then
            call interpffte(nt2,urrconc,m,urrconmc,wm)
            call interpffto(nt2,urtconc,m,urtconmc,wm)
            call interpffte(nt2,uttconc,m,uttconmc,wm)
         else
            call interpfft(nt2,urrconc,m,urrconmc,wm)
            call interpfft(nt2,urtconc,m,urtconmc,wm)
            call interpfft(nt2,uttconc,m,uttconmc,wm)
         end if
         urrconm(1:m)=urrconmc(1:m)
         urtconm(1:m)=urtconmc(1:m)
         uttconm(1:m)=uttconmc(1:m)
c
c         write(*,*) 'psiRin, temp1-1',urtconmc(1:m)
    
        

         idsamp3=m/ksamp3
         do j=1, ksamp3     !initial guess of the uniform theta grid in plasma domain
            im(j)=(j-1)*idsamp3+1
            wcon=rndconm(im(j))*exp(ima*tconm(im(j)))      
            nlimz=epslogzk/log(rndconm(im(j)))+1
            if (nlimz.gt.nt2) then
               nlimz=nt2
            end if
            call fft_cauchy(nlimz,wcon,zk,cint)
            Rconm(j)=cint
            Zconm(j)=-ima*cint
            Ravg(k)=Ravg(k)+Rconm(j)
            Zavg(k)=Zavg(k)+Zconm(j)
         end do
         Ravg(k)=Ravg(k)/ksamp3
         Zavg(k)=Zavg(k)/ksamp3

c         write(*,*) "(Ravg,Zavg)=",Ravg(k),Zavg(k)
c         write(*,*) "(Rmax,Zmax)=",Rmaxpsi,Zmaxpsi
         tmin=100.0
         do j=1, ksamp3
c            tmaxpm(j)=datan2((Zconm(j)-Zavg(k)),
c     1           (Rconm(j)-Ravg(k)))
            tmaxpm(j)=datan2((Zconm(j)-Zmaxis),
     1           (Rconm(j)-Rmaxis))
            if (tmaxpm(j).lt.0.0d0) then 
               tmaxpm(j)=2.0d0*pi+tmaxpm(j)
            end if 
            if (tmaxpm(j).lt.tmin) then
               tmin=tmaxpm(j)
               jtmin=j
            end if
         end do
c        write(*,*)'tt1',tmaxpm(1:ksamp3),Rconm(1:ksamp3),Zconm(1:ksamp3)
         do j=1,ksamp3
            jarr(j)=j
         end do
         jlag(1:klag)=jarr(ksamp3-klag+1:ksamp3)
         mlag(1:klag)=im(ksamp3-klag+1:ksamp3)-m      
         tlag(1:klag)=tmaxpm(ksamp3-klag+1:ksamp3)-2.0d0*pi
         
         jlag(klag+1:klag+ksamp3)=jarr(1:ksamp3)
         mlag(klag+1:klag+ksamp3)=im(1:ksamp3)      
         tlag(klag+1:klag+jtmin-1)=tmaxpm(1:jtmin-1)-2.0d0*pi
         tlag(klag+jtmin:klag+ksamp3)=tmaxpm(jtmin:ksamp3)

         jlag(klag+ksamp3+1:klag+ksamp3+jtmin-1)=jarr(1:jtmin-1)
         mlag(klag+ksamp3+1:klag+ksamp3+jtmin-1)=im(1:jtmin-1)+m      
         tlag(klag+ksamp3+1:klag+ksamp3+jtmin-1)=tmaxpm(1:jtmin-1)

         jlag(klag+ksamp3+jtmin:2*klag+ksamp3+jtmin-1)=
     1        jarr(jtmin:jtmin+klag-1)
         mlag(klag+ksamp3+jtmin:2*klag+ksamp3+jtmin-1)=
     1        im(jtmin:jtmin+klag-1)+m      
         tlag(klag+ksamp3+jtmin:2*klag+ksamp3+jtmin-1)=
     1        tmaxpm(jtmin:jtmin+klag)+2.0d0*pi
         isortstart = itotc +1
         icon(k)=isortstart
         dth=0.0d0 !2*pi/ksamp3*0.001d0
         !asuume dth=0 ---> thksamp3(j=1)=0.0d0
         istart=jtmin-1+klag

         klagl=floor(klag/2.0)
         klagr=klag-klagl

         eps13=1.0d-13
         do j=1, ksamp3
c            isort(i)=isortstart+i-1
            itotc=itotc+1
            thksamp3= 2.0d0*pi/ksamp3*(j-1)+dth
            ifind = 0
            jj=0
                        
            do while ((ifind.eq.0) .and. (jj.lt.10)) !find the theta grid which is equispaced as much as possible (not perfectly) 
               jj=jj+1
               if (jj.eq.1) then
                  do while (thksamp3.gt.tlag(istart)) 
                     istart=istart+1
                  end do
                  tlagtmp(1:klag)=tlag(istart-klagl:istart+klagr)
                  mlagtmp(1:klag)=mlag(istart-klagl:istart+klagr)  !index of m

                  call IntLagAll(klag, tlagtmp, mlagtmp,
     1              thksamp3,y1,wbary,eps13,isame)
                  if (isame.ne.0) then
                     ifind=1
                     itry2=modulo(nint(y1-1),m)+1
                     ttry2=tlagtmp(isame)
                     Rtry=Rconm(jlag(istart-klagl-1+isame))
                     Ztry=Zconm(jlag(istart-klagl-1+isame))
                     exit
                  else                                 
                     call CheckSameInt(klag,mlag(istart-klagl)
     1                   ,nint(y1),isame2)
                     if (isame2.ne.0) then
                        ifind=1
                        itry2=modulo(nint(y1-1),m)+1
                        ttry2=tlagtmp(isame2)
                        Rtry=Rconm(jlag(istart-klagl-1+isame2))
                        Ztry=Zconm(jlag(istart-klagl-1+isame2))
                        exit 
                     else
                        itry1=modulo(nint(y1-1),m)+1
                        wcon=rndconm(itry1)*exp(ima*tconm(itry1))

                        nlimz=epslogzk/log(rndconm(itry1))+1
                        if (nlimz.gt.nt2) then
                           nlimz=nt2
                        end if
                        call fft_cauchy(nlimz,wcon,zk,cint)
                        Rtry=cint
                        Ztry=-ima*cint
c                        ttry1=datan2((Ztry-Zavg(k)),(Rtry-Ravg(k)))
                        ttry1=datan2((Ztry-Zmaxis),(Rtry-Rmaxis))
                        if (ttry1.lt.0.0d0) then 
                           ttry1=2.0d0*pi+ttry1
                        end if 
                     end if
                  end if
               else
                  call CheckSameInt(klag,mlag(istart-klagl)
     1                   ,nint(y2),isame2)
                  if (isame2.ne.0) then
                     ifind=1
                     itry2=modulo(nint(y2-1),m)+1
                     ttry2=tlagtmp(isame2)
                     Rtry=Rconm(jlag(istart-klagl-1+isame2))
                     Ztry=Zconm(jlag(istart-klagl-1+isame2))
                     exit 
                  else
                     itry1=itry2
                     ttry1=ttry2
                     y1=y2
                  end if
               end if
               tlagtmp(klag+jj)=ttry1
               mlagtmp(klag+jj)=nint(y1)
               call IntLagAdd(klag+jj, tlagtmp, mlagtmp,  
     1              thksamp3,y2,wbary,eps13,isame) 
               
               itry2=modulo(nint(y2-1),m)+1
               if (itry1.eq.itry2) then
                  ifind=1
                  ttry2=ttry1
                  exit
               else
                  wcon=rndconm(itry2)*exp(ima*tconm(itry2))

                  nlimz=epslogzk/log(rndconm(itry2))+1
                  if (nlimz.gt.nt2) then
                     nlimz=nt2
                  end if

                  call fft_cauchy(nlimz,wcon,zk,cint)
                  Rtry=cint
                  Ztry=-ima*cint
c                  ttry2=datan2((Ztry-Zavg(k)),(Rtry-Ravg(k)))
                  ttry2=datan2((Ztry-Zmaxis),(Rtry-Rmaxis))
                  if ((ttry2.lt.0.0d0).and.(j.ne.1)) then 
                     ttry2=2.0d0*pi+ttry2
                  end if
               end if
            end do
  
           
            tmaxp(itotc)=ttry2
            if ((j.gt.1).and.(tmaxp(itotc-1).gt.tmaxp(itotc))) then 
               tmaxp(itotc-1)=tmaxp(itotc-1)-2.0d0*pi
            end if
            Rcon(itotc)=Rtry
            Zcon(itotc)=Ztry

            wcon=rndconm(itry2)*exp(ima*tconm(itry2))

            nlimzw=epslogzwk/log(rndconm(itry2))+1
            if (nlimzw.gt.nt2) then
               nlimzw=nt2
            end if
            call fft_cauchy(nlimzw,wcon,dzdw2k,cint)
            dzdwcon = cint
            dwdzr = (1.0d0,0.0d0)/cint
            dwdzi = -ima*(1.0d0,0.0d0)/cint

            nlimzww=epslogzwwk/log(rndconm(itry2))+1
            if (nlimzww.gt.nt2) then
               nlimzww=nt2
            end if
            call fft_cauchy(nlimzww,wcon,dzdww2k,cint)
            dzdwwcon = cint
            dwdzzcon = -dzdwwcon/(dzdwcon**3)
            dwdzzr=dwdzzcon
            dwdzzi=(-ima)*dwdzzcon

            dct=dcos(tconm(itry2))
            dst=dsin(tconm(itry2))
            dc2t=dcos(2.0d0*tconm(itry2))

            uxdisk=urconm(itry2)*dct
     1           +1/rndconm(itry2)*utconm(itry2)*(-dst)
            uydisk=urconm(itry2)*dst
     1           +1/rndconm(itry2)*utconm(itry2)*(dct) 
            uxxdisk=urrconm(itry2)*dct**2
     1           +1/rndconm(itry2)*urconm(itry2)*dst**2
     1           +1/rndconm(itry2)**2*uttconm(itry2)*dst**2
     1           +2/rndconm(itry2)**2*utconm(itry2)*dst*dct
     1           -2/rndconm(itry2)*urtconm(itry2)*dst*dct
            
            uyydisk=urrconm(itry2)*dst**2
     1           +1/rndconm(itry2)*urconm(itry2)*dct**2
     1           +1/rndconm(itry2)**2*uttconm(itry2)*dct**2
     1           -2/rndconm(itry2)**2*utconm(itry2)*dst*dct
     1           +2/rndconm(itry2)*urtconm(itry2)*dst*dct
            
            uxydisk=urrconm(itry2)*dst*dct
     1           -1/rndconm(itry2)*urconm(itry2)*dst*dct
     1           -1/rndconm(itry2)**2*uttconm(itry2)*dst*dct
     1           -1/rndconm(itry2)**2*utconm(itry2)*dc2t
     1           +1/rndconm(itry2)*urtconm(itry2)*dc2t

            ux=uxdisk*dwdzr+uydisk*dwdzi
            uy=uxdisk*(-dwdzi)+uydisk*dwdzr
            uxx=uxxdisk*dwdzr**2+uxdisk*dwdzzr
     1           +2*uxydisk*dwdzr*dwdzi
     1           +uyydisk*dwdzi**2+uydisk*dwdzzi

            uyy=uxxdisk*dwdzi**2-uxdisk*dwdzzr
     1           -2*uxydisk*dwdzr*dwdzi
     1           +uyydisk*dwdzr**2-uydisk*dwdzzi

            uxy=-uxxdisk*dwdzr*dwdzi-uxdisk*dwdzzi
     1           +uxydisk*(dwdzr**2-dwdzi**2)
     1           +uyydisk*dwdzr*dwdzi+uydisk*dwdzzr

            
            psiRcon(itotc)=(psicon(k))/(2.0d0*Rtry)
     1              +ux*dsqrt(Rtry)
            psiZcon(itotc)=uy*dsqrt(Rtry)

            psiRRcon(itotc) = -(psicon(k))/(4.0d0*Rtry**2)
     1           +ux/(sqrt(Rtry))
     1           + uxx*sqrt(Rtry)
            psiRZcon(itotc) = uy/(2.0d0*sqrt(Rtry))
     1           + uxy*sqrt(Rtry)
            psiZZcon(itotc) = uyy*sqrt(Rtry)

c            write(*,*) 'dwdzz',dwdzzr,dwdzzi
c            write(*,*) 'urrconm',urrconm(itry2),urconm(itry2)
c            write(*,*) 'urtconm',urtconm(itry2),uttconm(itry2)
          
c            write(*,*) 'uxxdisk',uxxdisk,uyydisk,uxydisk
c            write(*,*) 'uxx',uxx,uyy,uxy
 
c           write(*,*) 'psiRRcon(itotc)',psiRRcon(itotc),psiZZcon(itotc)
c     1           ,psiRZcon(itotc)
            select case (iptype)
            case(0)
               call Solovev(Rcon(itotc),Zcon(itotc),csol,d1,d2,d3
     1       ,psiconex(itotc),psiRconex(itotc),psiZconex(itotc)
     2       ,psiRRconex(itotc),psiRZconex(itotc),psiZZconex(itotc)
     3              ,fsolcon)
c           write(*,*) 'psiRRconex',psiRRconex(itotc),psiZZconex(itotc)
c     1           ,psiRZconex(itotc),uy/(2.0d0*sqrt(Rtry)),uxy*sqrt(Rtry)
            end select
         end do
   
         call system_clock(ict, icr, ic_max)
         if (ict.gt.ici) then
            tsec=real(ict-ici)/real(icr)
         else
            tsec=real(ict-ici+ic_max)/real(icr)
         end if
c        write(*,*) "CPU time sortcontour4: t(sec)= ",tsec 

         !find outer mid-plane point and inner mid-plane point
         tmaxpmid(1)=tmaxp(1) !(isortstart)
         tmaxpmid(2)=tmaxp(ksamp3/2+1) !isortstart)
c         write(*,*) 'tmaxpmid',tmaxpmid
c         write(*,*) 'tmaxp',tmaxp(1:ksamp3)
         do j=1,klag
            jj=modulo(jtmin-klagl-2+j,m)+1
            wcon=rndconm(jj)*exp(ima*tconm(jj))
            nlimz=epslogzk/log(rndconm(jj))+1
            if (nlimz.gt.nt2) then
               nlimz=nt2
            end if
            call fft_cauchy(nlimz,wcon,zk,cint)
            Rtmp(j)=cint
            Ztmp(j)=-ima*cint
c            tlagtmp(j)=datan2((Ztmp(j)-Zavg(k)),(Rtmp-Ravg(k)))
            tlagtmp(j)=datan2((Ztmp(j)-Zmaxis),(Rtmp(j)-Rmaxis))
         end do
         epsLag22=epsLag2
         call IntLagAll(klag,tlagtmp,Rtmp,
     1           tmaxpmid(1),Rmid,wlag,epsLag22, isame)
c         write(*,*) 'tlagtmp',tlagtmp(1:klag),Rtmp(1:klag),Rmid
         Rmido(k)=Rmid
         
 
c         dZdt= 0.0d0
c         epsLag22=epsLag2
c         do while (dZdt.eq.0.0d0)
c            call IntLag_sameptsAll(klag,tlagtmp,Ztmp,
c     1           tmaxpmid(1),Zmid(1),dZdt,epsLag22)
c            write(*,*) 'epsLag2 is changed to ',epsLag22
c            if (epsLag22.gt.1d-8) then
c               write(*,*) 'no z grid in outer-midplane'
c               if (isymud.eq.1) then
c                  write(*,*) 'contour integrals may not be accurate'
c                  isymud=0
c               end if
c               exit
c            else if (dZdt.eq.0.0d0) then
c               epsLag22=epsLag22*10
c               write(*,*) 'epsLag2 is changed to ',epsLag22
c            end if
c         end do
c         epsLag2=epsLag22
c         dZdtmido(k)=dZdt

         do j=1,klag
            jj=modulo(jtmin-klagl-2+m/2+j,m)+1
c            write(*,*) 'tconma',tconm(jj)-pi
            wcon=rndconm(jj)*exp(ima*tconm(jj))
            nlimz=epslogzk/log(rndconm(jj))+1
            if (nlimz.gt.nt2) then
               nlimz=nt2
            end if
            call fft_cauchy(nlimz,wcon,zk,cint)
            Rtmp(j)=cint
            Ztmp(j)=-ima*cint
            tlagtmp(j)=datan2((Ztmp(j)-Zmaxis),(Rtmp(j)-Rmaxis))
            if (tlagtmp(j).lt.0.0d0) then 
               tlagtmp(j)=2.0d0*pi+tlagtmp(j)
            end if
         end do

         call IntLagAll(klag,tlagtmp,Rtmp,
     1           tmaxpmid(2),Rmid,wlag,epsLag22, isame)

         Rmidi(k)=Rmid
c         dZdt= 0.0d0
c         epsLag22=epsLag2
c         do while (dZdt.eq.0.0d0)
c            call IntLag_sameptsAll(klag,tlagtmp,Ztmp,
c     1           tmaxpmid(2),Zmid(2),dZdt,epsLag22)         
c            if (epsLag22.gt.1d-8) then
c               write(*,*) 'no z grid in inner-midplane'
c               if (isymud.eq.1) then
c                  write(*,*) 'contour integrals may not be accurate'
c                  isymud=0
c               end if
c               exit
c            else if (dZdt.eq.0.0d0) then
c               epsLag22=epsLag22*10
c               write(*,*) 'epsLag2 is changed to ',epsLag22
c            end if
c         end do
c         dZdtmidi(k)=dZdt
c         write(*,*) 'klag',klag,tlagtmp(1:klag),Ztmp(1:klag)
c         write(*,*) 'klag2',tmaxpmid(2), Zmid(2),dZdt

c         write(*,*) "finished sortcontour0 for icon=",k

         end if
c         write(*,*) 'Rcon3',Rcon((k-1)*ksamp3+1:k*ksamp3)
      end do
      icon(ncon+1)=itotc+1
     
      write(*,*) "finished sortcontour1"
      return
      end

      subroutine checkcontour(ncon,psicon,nr,ntheta,rnd,tnd,psi,icheck)

      use arrays, only: kcheb,cftmsub, bnodes

      implicit real*8 (a-h,o-z)
      save

      real *8 maxu
c      real *8 rndin(nsub),
      real *8 psicon(*),rnd(*),tnd(*),psi(*)

      integer i,j,k,imaxu,jmaxt, isub, nr
      integer inextsub

c      nr= korder*nsub
c     pick up the position in the chebycheb grid of each subinterval for maximum u in the contour
c      irinsub=1 


      itotc = 0   
      itotcp=1
      icheck =1

    

c     check if the conformal mapping center is fine 
c     to construct the coutour close to the magnetic axis.
! number of contour plots is ncon
c         imaxu=(i-1)*korder+irinsub  ! the r index for maxsimum u at each controur
c         maxu = -1000.0d0
c        Find the crossing points of (psi-psicon(k)) at theta=0
      k=1  !first contour
         psiiB = 0.0d0
         do j=1,ntheta
            inext=1+(j-1)*nr
            temp1 = psi(inext)-psicon(k)
            ncross = 0
            do i=2,nr+1
               inext=i+(j-1)*nr
               if (i.eq.nr+1) then
                  temp2=psiiB-psicon(k)
               else
                  temp2=psi(inext)-psicon(k)
               end if
               if (temp1*temp2 .lt. 0.0d0) then                  
                  ncross = ncross+1
               end if
               temp1=temp2
            end do
            if (ncross.ne.1)  then
               icheck = 0
               exit
            end if
         
         end do
 
            
      return
      end

      subroutine findlocalshear(ncon,psicon,icon,tmaxp,Rcon,Zcon
     1     ,psiR,psiZ,psiRR,psiRZ,psiZZ,fpolcon,dfdpsi,shatlocal0
     2     ,nttin,ttin,shatlocal2)

      use arrays, only:  nt2, ksamp2,epstri,isymud,lambda

      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*),tmaxp(*),ttin(*),psiR(*),psiZ(*)
      real *8 Rcon(*), Zcon(*),psiRR(*),psiRZ(*),psiZZ(*)
      real *8 fpolcon(*),dfdpsi(*),shatlocal0(*),shatlocal2(*)
      real *8 thin(nt2*ksamp2+1),temp(nttin)
      real *8 thintmp(nt2*ksamp2+1), ttintmp(nttin)

      integer ncon,nttin,i,j,k,icon(*)

      !psi in this routine is normalized one to the range from 0 to 1
      pi=4.0d0*datan(1.0d0)
      do k=1, ncon
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         thin(1:nin)=tmaxp(istart:iend)
         thin(nin+1)=tmaxp(istart)+2.0d0*pi

         if (isymud.eq.0) then
            thintmp(1:nint+1)=thin(1:nint+1)/2.0d0
            ttintmp(1:nttin)=ttin(1:nttin)/2.0d0
         end if

         do j=1,nin
            kj=istart+j-1
            dpsi2=(psiR(kj)**2+psiZ(kj)**2)/lambda**2
            cd2psidR2=fpolcon(k)/Rcon(kj)/dpsi2
     1           *(1.0d0-2.0d0*psiR(kj)**2/lambda**2/dpsi2)
            cd2psidZ2=fpolcon(k)/Rcon(kj)/dpsi2
     1           *(1.0d0-2.0d0*psiZ(kj)**2/lambda**2/dpsi2)
            cdpsidR=(dfdpsi(k)*psiR(kj)/lambda/Rcon(kj)
     1           -fpolcon(k)/Rcon(kj)**2)/dpsi2
c            write(*,*) 'shatlocal00k',dpsi2,psiR(kj),psiZ(kj)
c     1           ,psiRR(kj),psiRZ(kj),psiZZ(kj)
            cdpsidZ=dfdpsi(k)*psiZ(kj)/lambda/Rcon(kj)/dpsi2
            cd2psidRdZ=-4.0d0*psiR(kj)*psiZ(kj)/lambda**2
     1           *fpolcon(k)/Rcon(kj)/dpsi2**2

            shatlocal0(kj)=1/Rcon(kj)/lambda*(cd2psidR2*psiRR(kj)
     1           +cd2psidZ2*psiZZ(kj)+cdpsidR*psiR(kj)
     2           +cdpsidZ*psiZ(kj)+cd2psidRdZ*psiRZ(kj))

c            write(*,*) 'shatlocal00',cd2psidR2,psiRR(kj) 
c     1            ,cd2psidZ2,psiZZ(kj),cdpsidR,psiR(kj)
c     2           ,cdpsidZ,psiZ(kj),cd2psidRdZ,psiRZ(kj)


            dpsi2part1=(psiZ(kj)**2-psiR(kj)**2)*(psiRR(kj)-psiZZ(kj))
     1           /dpsi2**2/lambda**3
            Reffect1=-psiR(kj)/Rcon(kj)/dpsi2/lambda
            dpsidRdZpart1=-4.0d0*psiR(kj)*psiZ(kj)*psiRZ(kj)
     1           /dpsi2**2/lambda**3

            shattemp1=(dpsi2part1+Reffect1+dpsidRdZpart1)
     1           *fpolcon(k)/Rcon(kj)**2+dfdpsi(k)/Rcon(kj)**2
c            write(*,*) 'shatlocalcomp1', dpsi2part1,Reffect1
c     1           ,dpsidRdZpart1,fpolcon(k)/Rcon(kj)**2,
c     2           dfdpsi(k)/Rcon(kj)**2
c            write(*,*) 'shatlocalcomp', shatlocal0(kj),
c     1           shattemp1, shatlocal0(kj)-shattemp1
         end do
         
         if (isymud.eq.1) then
  
            call IntTriCos(nin/2+1,thin,shatlocal0(istart),nttin/2+1,
     1           ttin,temp,epstri)
            ist=(k-1)*nttin+1
            shatlocal2(ist:ist+nttin/2)=temp(1:nttin/2+1)
            do i=2,nttin/2
               shatlocal2(ist+nttin+1-i)=temp(i)
            end do
c            write(*,*) 'shatlocal20',ist,istart
c            write(*,*) 'shatlocal20',shatlocal0(istart:istart+nin-1)
c            write(*,*) 'shatlocal2', shatlocal2(ist:ist+nttin-1)
         else
            
            call IntTriCos(nin+1,thintmp,shatlocal0(istart),nttin,
     1           ttintmp,temp,epstri)
            ist=(k-1)*nttin+1
            shatlocal2(ist:ist+nttin-1)=temp(1:nttin)
         end if
      end do

      return
      end

      subroutine findcontour2der(ncon,psicon,nr,ntheta,rnd,tnd,psi
     1     ,dpsidr,dpsidrr,rcon,tcon,urcon,utcon
     2     ,urrcon,urtcon,uttcon)
      use arrays, only: kcheb,cftmsub, bnodes, ur,nsub
      use arrays, only: ur,uth,urr,urt,utt
      implicit real*8 (a-h,o-z)
      save

      real *8 maxu
c      real *8 rndin(nsub),
      real *8 psicon(*),rnd(*),tnd(*),psi(*),dpsidr(*)
      real *8 urcon(*), utcon(*),urrcon(*), urtcon(*),uttcon(*)
      real *8 dpsidrr(*),rcon(*),tcon(*)
      real *8 chcoeff(kcheb), psisub(kcheb)
      integer *4 ncross

      integer i,j,k,imaxu,jmaxt, isub, nr
      integer inextsub

c      nr= korder*nsub
c     pick up the position in the chebycheb grid of each subinterval for maximum u in the contour
c      irinsub=1 

      write(*,*)' bnodes are',(bnodes(j),j=1,nsub+1)

      eps10=1.0d-5
      itotc = 0   
      itotcp=1
      ncross = 0
      do k=1,ncon              ! number of contour plots is ncon

         write(*,*) 'psicon(k)',psicon(k)
         psiiB = 0.0d0
         do j=1,ntheta
            temp1 = 1.0d0-psicon(k)
            
c            inext=1+(j-1)*nr
c            temp1 = psi(inext)-psicon(k)
c             write(*,*) 'temp1',temp1,psi(inext),psicon(k)
c            if (temp1.eq.0.0d0) then
c               isub=1
c               inextsub=((isub-1)*kcheb+1)+(j-1)*nr
c               itotc = itotc + 1
c               rl=0.0d0
c               rr=rnd(1)
c            else
            do i=1,nr+1
               inext=i+(j-1)*nr
               if (i.eq.nr+1) then
                  temp2=psiiB-psicon(k)
               else
                  temp2=psi(inext)-psicon(k)
               end if
c               write(*,*) 'temp1,temp2',temp1,temp2
               if (temp1*temp2 .le. 0.0d0) then
                  
                  if (i.eq.1) then
                     rl=0.0d0+eps10
                     rr=rnd(1)
                     isub=1
                  else 
                     rl=rnd(i-1)
                     if (i.eq.nr+1) then
                        rr=1.0d0-eps10
                     else 
                        rr=rnd(i)
                     end if         
                     isub=(i-2)/kcheb+1
                  end if
                  
                  inextsub=((isub-1)*kcheb+1)+(j-1)*nr
                  itotc = itotc + 1
                  exit
               end if
               temp1=temp2
            end do
c            write(*,*) 'isub(k)',isub,j
c                  if(psicon(k).lt.0.05) then
c                     write(*,*) 'psi(inext,psicon(k),itotc',itotc,j
c                     write(*,*) psi(inext-1), psi(inext),psicon(k)
c                     write(*,*) rnd(i-1),rnd(i)
c                     write(*,*) i,isub,bnodes(isub),bnodes(isub+1)
c                     write(*,*) psi(inextsub:inextsub+kcheb)
c                  end if
                  psisub(1:kcheb)=psi(inextsub:inextsub+kcheb-1)
                  call chftransq(chcoeff,psisub,kcheb,cftmsub)

                  eps7=1.0d-14
                  dpsic=temp1 
                  rc=rl
                  xch=(rc-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
                  idomax=20
                  ido=0
                  do while ((dabs(dpsic).gt.eps7).and.(ido.lt.idomax))
                     ido=ido+1
c                     write(*,*) 'rc',rc,xch
                     call chderiv(xch, kcheb, chcoeff, dpsidrch)
                     dpsidrch= -dpsidrch/(bnodes(isub+1)-bnodes(isub))

                     rc=rc-dpsic/dpsidrch
                     if (rc.gt.rr) then
                        rc=rr               
c                        exit
                     else if (rc.lt.rl) then
                        rc=rl
                        exit
                     end if

                   
                     xch=(rc-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
                     if (xch.gt.1.0d0) exit
                    
                     
                     call chfit(xch, kcheb, chcoeff,psic)
                     dpsic=psic-psicon(k)
                     if(psicon(k).lt.0.05) then
                      end if
                  end do

c     If the interpolated point is outside the subnode,  
                  if (xch.gt.1.0d0) then 
                     write (*,*) 'xch>1'
                     isub=isub+1
                     inextsub=((isub-1)*kcheb+1)+(j-1)*nr
                     psisub(1:kcheb)=psi(inextsub:inextsub+kcheb-1)
                     call chftransq(chcoeff,psisub,kcheb,cftmsub)
                     dpsic=temp2
                     rc=rr
                     xch=(rc-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
                     do while (dabs(dpsic).gt.eps7)
                 
                        call chderiv(xch, kcheb, chcoeff, dpsidrch)
                        dpsidrch=-dpsidrch/(bnodes(isub+1)-bnodes(isub))
                        rc=rc-dpsic/dpsidrch
                        xch=(rc-bnodes(isub))
     1                      /(bnodes(isub+1)-bnodes(isub))
                           
                        if (xch.lt.0.0d0) exit
                        call chfit(xch, kcheb, chcoeff,psic)
                        dpsic=psic-psicon(k)
                         
                     end do
                  end if
                     
                  if (xch.lt.0.0d0) then  !unlikely happend
                     dr1=temp1/(-dpsidr(inext-1))
                     dr2=-dr1**2.0d0*dpsidrr(inext-1)/2.0d0
     1                 /(dpsidr(inext-1)+dpsidrr(inext-1)*dr1)
                     rcon(itotc) = rnd(i-1)+dr1+dr2
                     urcon(itotc) = ur(inext-1)
     1                 +(rcon(itotc)-rnd(i-1))*urr(inext-1)
                     utcon(itotc) = uth(inext-1)
     1                 +(rcon(itotc)-rnd(i-1))*urt(inext-1)
                     urrcon(itotc) = urr(inext-1)
     1                 +(rcon(itotc)-rnd(i-1))/(rnd(i)-rnd(i-1))
     2                    *(urr(inext)-urr(inext-1))
                     urtcon(itotc) = urt(inext-1)
     1                 +(rcon(itotc)-rnd(i-1))/(rnd(i)-rnd(i-1))
     2                    *(urt(inext)-urt(inext-1))
                     uttcon(itotc) = utt(inext-1)
     1                 +(rcon(itotc)-rnd(i-1))/(rnd(i)-rnd(i-1))
     2                    *(utt(inext)-utt(inext-1))
                  else 
                     rcon(itotc) = rc
                     call chftransq(chcoeff,ur(inextsub),kcheb,cftmsub)
                     call chfit(xch, kcheb, chcoeff,urcon(itotc))
                     call chftransq(chcoeff,uth(inextsub),kcheb,cftmsub)
                     call chfit(xch, kcheb, chcoeff,utcon(itotc))
                     call chftransq(chcoeff,urr(inextsub),kcheb,cftmsub)
                     call chfit(xch, kcheb, chcoeff,urrcon(itotc))
                     call chftransq(chcoeff,urt(inextsub),kcheb,cftmsub)
                     call chfit(xch, kcheb, chcoeff,urtcon(itotc))
                     call chftransq(chcoeff,utt(inextsub),kcheb,cftmsub)
                     call chfit(xch, kcheb, chcoeff,uttcon(itotc))
                  end if

                  tcon(itotc) = tnd(j)
                  ncross= ncross+1

c               end if
c               end if
c               temp1=temp2
c            end do
               
         end do
         itotcp=itotc+1

c         write(*,*) 'rcon1D',rcon((k-1)*ntheta+1:k*ntheta) 
      end do
 
            
      return
      end
c
c
c
c
      subroutine findcontour(ncon,psicon,nr,ntheta,rnd,tnd,psi,dpsidr,
     1     dpsidrr,rcon,tcon,urcon,utcon)

      use arrays, only: kcheb,cftmsub, bnodes, nsub
      use arrays, only: ur,uth,urr,urt
      implicit real*8 (a-h,o-z)
      save

      real *8 maxu
c      real *8 rndin(nsub),
      real *8 psicon(*),rnd(*),tnd(*),psi(*),dpsidr(*)
      real *8 urcon(*), utcon(*)
      real *8 dpsidrr(*),rcon(*),tcon(*)
      real *8 chcoeff(kcheb), psisub(kcheb)
      integer *4 ncross

      integer i,j,k,imaxu,jmaxt, isub, nr
      integer inextsub

c      nr= korder*nsub
c     pick up the position in the chebycheb grid of each subinterval for maximum u in the contour
c      irinsub=1 

c      bnodes(1)=0.0d0
      write(*,*)' bnodes are',(bnodes(j),j=1,nsub+1)


      itotc = 0   
      itotcp=1
      ncross = 0
      do k=1,ncon              ! number of contour plots is ncon

         write(*,*) 'psicon(k)',psicon(k)
         psiiB = 0.0d0
         do j=1,ntheta
            inext=1+(j-1)*nr
            temp1 = psi(inext)-psicon(k)
            do i=2,nr+1
               inext=i+(j-1)*nr
               if (i.eq.nr+1) then
                  temp2=psiiB-psicon(k)
               else
                  temp2=psi(inext)-psicon(k)
               end if

               if (temp1*temp2 .lt. 0.0d0) then
                  isub=(i-2)/kcheb+1
                  inextsub=((isub-1)*kcheb+1)+(j-1)*nr
                  itotc = itotc + 1
c                  if(psicon(k).lt.0.05) then
c                     write(*,*) 'psi(inext,psicon(k),itotc',itotc,j
c                     write(*,*) psi(inext-1), psi(inext),psicon(k)
c                     write(*,*) rnd(i-1),rnd(i)
c                     write(*,*) i,isub,bnodes(isub),bnodes(isub+1)
c                     write(*,*) psi(inextsub:inextsub+kcheb)
c                  end if
                  psisub(1:kcheb)=psi(inextsub:inextsub+kcheb-1)
                  call chftransq(chcoeff,psisub,kcheb,cftmsub)

                  eps7=1.0d-14
                  dpsic=temp1 
                  rc=rnd(i-1)
                  xch=(rc-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))

                  do while (dabs(dpsic).gt.eps7)
                     call chderiv(xch, kcheb, chcoeff, dpsidrch)
                     dpsidrch= -dpsidrch/(bnodes(isub+1)-bnodes(isub))

                     rc=rc-dpsic/dpsidrch
                     if (rc.gt.rnd(i)) then
                        rc=rnd(i)               
c                        exit
                     else if (rc.lt.rnd(i-1)) then
                        rc=rnd(i-1)
                        exit
                     end if

                   
                     xch=(rc-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
                     if (xch.gt.1.0d0) exit
                    
                     
                     call chfit(xch, kcheb, chcoeff,psic)
                     dpsic=psic-psicon(k)
                     if(psicon(k).lt.0.05) then
                      end if
                  end do

c     If the interpolated point is outside the subnode,  
                  if (xch.gt.1.0d0) then 
                     write (*,*) 'xch>1'
                     isub=isub+1
                     inextsub=((isub-1)*kcheb+1)+(j-1)*nr
                     psisub(1:kcheb)=psi(inextsub:inextsub+kcheb-1)
                     call chftransq(chcoeff,psisub,kcheb,cftmsub)
                     dpsic=temp2
                     rc=rnd(i)
                     xch=(rc-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
                     do while (dabs(dpsic).gt.eps7)
                 
                        call chderiv(xch, kcheb, chcoeff, dpsidrch)
                        dpsidrch=-dpsidrch/(bnodes(isub+1)-bnodes(isub))
                        rc=rc-dpsic/dpsidrch
                        xch=(rc-bnodes(isub))
     1                      /(bnodes(isub+1)-bnodes(isub))
                           
                        if (xch.lt.0.0d0) exit
                        call chfit(xch, kcheb, chcoeff,psic)
                        dpsic=psic-psicon(k)
                         
                     end do
                  end if
                     
                  if (xch.lt.0.0d0) then
                     dr1=temp1/(-dpsidr(inext-1))
                     dr2=-dr1**2.0d0*dpsidrr(inext-1)/2.0d0
     1                 /(dpsidr(inext-1)+dpsidrr(inext-1)*dr1)
                     rcon(itotc) = rnd(i-1)+dr1+dr2
                     urcon(itotc) = ur(inext-1)
     1                 +(rcon(itotc)-rnd(i-1))*urr(inext-1)
                     utcon(itotc) = uth(inext-1)
     1                 +(rcon(itotc)-rnd(i-1))*urt(inext-1)
                  else 
                     rcon(itotc) = rc
                     call chftransq(chcoeff,ur(inextsub),kcheb,cftmsub)
                     call chfit(xch, kcheb, chcoeff,urcon(itotc))
                     call chftransq(chcoeff,uth(inextsub),kcheb,cftmsub)
                     call chfit(xch, kcheb, chcoeff,utcon(itotc))
                  end if

                  tcon(itotc) = tnd(j)
                  ncross= ncross+1

c               end if
               end if
               temp1=temp2
            end do
               
         
         end do
         itotcp=itotc+1
      end do
 
            
      return
      end
c
c
c
c
      subroutine findq0(ncon,lambda,Ic,cftmq,F0, q0)

      use arrays, only: fpolsgn
      implicit real*8 (a-h,o-z)
      save

      real *8 Ic(*), cftmq(*), Ic0
      real *8 chcoeff1(ncon)
      real *8 lambda,pi,q0,F0

      pi=4.0d0*datan(1.0d0)
      call chftransq(chcoeff1, Ic, ncon, cftmq)
      call chfit(0.0d0,ncon,chcoeff1,Ic0)
      q0=Ic0/(2*pi)*dabs(F0*lambda)
      return
      end

      subroutine findF0(ncon,lambda,Ic,cftmq,q0,F0)

      use arrays, only: fpolsgn
      implicit real*8 (a-h,o-z)
      save

      real *8 Ic(*), cftmq(*), Ic0
      real *8 chcoeff1(ncon)
      real *8 lambda,pi,q0,F0

      pi=4.0d0*datan(1.0d0)
      call chftransq(chcoeff1, Ic, ncon, cftmq)
      call chfit(0.0d0,ncon,chcoeff1,Ic0)
      F0=q0/Ic0*2*pi*fpolsgn/dabs(lambda)

      return
      end

      subroutine normbyF0q0(ncon,lambda,Ic,cftmq,q0,F0,factorl)

      use arrays, only: fpolsgn
      implicit real*8 (a-h,o-z)
      save

      real *8 Ic(*), cftmq(*), Ic0
      real *8 chcoeff1(ncon)
      real *8 lambda,lambda0,pi,q0,F0

      pi=4.0d0*datan(1.0d0)
      call chftransq(chcoeff1, Ic, ncon, cftmq)
      call chfit(0.0d0,ncon,chcoeff1,Ic0)
      lambda0=2*pi/Ic0*dabs(q0/F0)
      factorl=lambda0/dabs(lambda)
      return
      end

      subroutine findtorflux(ncon,lambda,psiichq,qcon,
     1     spbx,spdef, phitot,phiich)

      use arrays, only: fpolsgn
      implicit real*8 (a-h,o-z)
      save

      real *8 psiichq(*),spbx(*), qcon(*)
      real *8 phiich(*),spdef(*)
      real *8 lambda,pi,FF0,FFin,ffout,F0

      pi=4.0d0*datan(1.0d0)
      phiich(1:ncon)=0.0d0
      phitot=0.0d0
      do i=1,ncon
         do k=1,ncon
            ! int_{-1}^{x} FFp((x+1)/2) dx
            phiich(i)=phiich(i)+spbx((k-1)*ncon+i)*qcon(k)
         end do
         phitot=phitot+spdef(i)*qcon(i)
      end do
      phiich(1:ncon)=phiich(1:ncon)/dabs(lambda*2.0d0)
      phitot=phitot/dabs(lambda*2.0d0)

      phiich(1:ncon)=phiich(1:ncon)/phitot !normalized toroidal flux
      write(*,*) 'phitot',phitot
      write(*,*) 'phiich',phiich(1:ncon)
      return
      end
      subroutine findFbyFFp0(ncon,lambda,psiichq,FFp,
     1     spbx,F0,fpolch)

      use arrays, only: fpolsgn
      implicit real*8 (a-h,o-z)
      save

      real *8 psiichq(*),spbx(*), fpolch(*)
      real *8 FFp(*), FFint(ncon)
      real *8 lambda,pi,FF0,FFin,ffout,F0

      pi=4.0d0*datan(1.0d0)
      FFint(1:ncon)=0.0d0
      do i=1,ncon
         do k=1,ncon
            ! int_{-1}^{x} FFp((x+1)/2) dx
            FFint(i)=FFint(i)+spbx((k-1)*ncon+i)*FFp(k)
         end do
      end do
      FFint(1:ncon)=FFint(1:ncon)/dabs(lambda*2.0d0)
      fpolch(1:ncon)=dsqrt(F0**2+2.0d0*FFint(1:ncon))*fpolsgn
c      write(*,*) 'fpolch',fpolch(1:ncon)
      return
      end

      subroutine findFbyFFp1(ncon,lambda,psiichq,FFp,
     1     spxb,Fedge,fpolch)

      use arrays, only: fpolsgn
      implicit real*8 (a-h,o-z)
      save

      real *8 psiichq(*),spxb(*), fpolch(*)
      real *8 FFp(*), FFint(ncon)
      real *8 lambda,pi,FF0,FFin,ffout,Fedge

      pi=4.0d0*datan(1.0d0)
      FFint(1:ncon)=0.0d0
      do i=1,ncon
         do k=1,ncon
            ! int_{x}^{1} FFp((x+1)/2) dx
            FFint(i)=FFint(i)+spxb((k-1)*ncon+i)*FFp(k)
         end do
      end do
      FFint(1:ncon)=FFint(1:ncon)/dabs(lambda*2.0d0)
      fpolch(1:ncon)=dsqrt(Fedge**2-2.0d0*FFint(1:ncon))*fpolsgn
c      write(*,*) 'fpolch',fpolch(1:ncon)
      return
      end

      subroutine findF0byFedge(ncon,lambda,psiichq,FFp,
     1     spdef,Fedge,F0)

      use arrays, only: fpolsgn
      implicit real*8 (a-h,o-z)
      save

      real *8 psiichq(*),spdef(*)
      real *8 FFp(*), FFint
      real *8 lambda,pi,FF0,FFin,ffout,Fedge,F0

      pi=4.0d0*datan(1.0d0)
      FFint=0.0d0
      do i=1,ncon
            ! int_{-1}^{1} FFp((x+1)/2) dx
            FFint=FFint+spdef(i)*FFp(i)
      end do
      FFint=FFint/dabs(lambda*2.0d0)
      F0=dsqrt(Fedge**2-2.0d0*FFint)*fpolsgn
c      write(*,*) 'fpolch',fpolch(1:ncon)
      return
      end

      subroutine findfpol(ncon,lambda,psicon,qpsich,Ic
     1    ,fpolcon)

      use arrays, only: fpolsgn
      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*),fpolcon(*), Ic(*)
      real *8 qpsich(*), lambda
      real *8 bfpol(2000), cfpol(2000), dfpol(2000)
 
      integer i,j,k

      pi=4.0d0*datan(1.0d0)
      do k=1,ncon
         fpolcon(k)=qpsich(k)/Ic(k)
      end do
      fpolcon(1:ncon)=fpolcon(1:ncon)*2*pi*fpolsgn/dabs(lambda)

c      write(*,*) "fpolcon2 before interp", fpolcon(1:ncon)
c      call spline(-psicon(2:ncon-1), fpolcon(2:ncon-1) 
c     1           , bfpol, cfpol, dfpol, ncon-2)
c      call ispline(-psicon(1),-psicon(2:ncon-1),fpolcon(2:ncon-1) 
c     1           , bfpol, cfpol, dfpol, ncon-2, fpolcon(1))
c      call ispline(-psicon(ncon),-psicon(2:ncon-1),fpolcon(2:ncon-1) 
c     1     , bfpol, cfpol, dfpol, ncon-2, fpolcon(ncon))
c      write(*,*) "fpolcon2 after interp", fpolcon(1:ncon)
      return
      end

      subroutine findffprimnext(ncon,lambda,psicon,cftmq,fpolcon2
     1    ,ffprimnext)
      implicit real*8 (a-h,o-z)
      save

      real *8 fpolcon2(*),ffprimnext(*),cftmq(*),psicon(*)
      real *8 lambda
      real *8 chcoeff(2000)
      integer i,j,k

      call chftransq(chcoeff, fpolcon2, ncon, cftmq)
     
      do i=1,ncon
         xch=1.0d0-psicon(i)
         call chderiv(xch, ncon, chcoeff, ffprimnext(i)) 
         ffprimnext(i)=ffprimnext(i)*lambda
c         call chfit(xch, ncon, chcoeff,fpolcon2fit)
c         write(*,*) 'fpolcon2fit',fpolcon2fit        
      end do

      return
      end

      subroutine fitMiller(ncon,icon,Rcon,Zcon,Ravg,Zavg
     1    ,epsFIT,R0mil,iaspmil,kappamil,deltamil)

      use arrays, only:  nt2, ksamp2
      

      implicit real*8 (a-h,o-z)
      save

      real *8 Rcon(*), Zcon(*), Ravg(*),Zavg(*)
      real *8 R0mil(*),iaspmil(*),kappamil(*),deltamil(*)
      real *8 Rin(nt2*ksamp2), Zin(nt2*ksamp2),theta(nt2*ksamp2)
      real *8 Rfit(nt2*ksamp2),t0(nt2*ksamp2)
      real *8 dely(nt2*ksamp2),Omega(3),Hessian(3,3)
      complex *16 A(3,3),y(3),delb(3),w7(1000)
      integer *4  icon(*)
      integer i,j,k, istart,iend

c      ksamp3=ksamp2/4 
      pi=4.0d0*datan(1.0d0)

c      write(*,*) 'psiRintot',psiRcon(1:ncon*ksamp3)
      do k=1, ncon
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         Rin(1:nin)=Rcon(istart:iend)
         Zin(1:nin)=Zcon(istart:iend)
         iiterf=0
         imaxiterf=30
         do while (iiterf .lt.imaxiterf) !maximum iterations for fiiting =30
            iiterf=iiterf+1
            if (iiterf.eq.1) then ! initial guess
               R0mil(k)=Ravg(k)
               Z0=Zavg(k)
          
               RatRmax = 0.0d0
               ZatZmax = 0.0d0
               ZatZmin = 0.0d0
               do i=1,nin
                  if (Rin(i) .gt. RatRmax) then
                     RatRmax = Rin(i)
                  end if
                  if (Zin(i) .gt. ZatZmax) then
                     ZatZmax = Zin(i)
                     RatZmax = Rin(i)
                     iZmax=i
                  end if
                  if (Zin(i) .lt. ZatZmin) then
                     ZatZmin = Zin(i)
                     iZmin=i
                  end if
               end do
               iaspmil(k)=(RatRmax-Ravg(k))/Ravg(k)
               kappamil(k)=ZatZmax/(Ravg(k)*iaspmil(k))
               deltamil(k)=(Ravg(k)-RatZmax)/(Ravg(k)*iaspmil(k))
c               write(*,*) 'deltamil',deltamil(k),Ravg(k),RatZmax
               b1=R0mil(k)
               b2=R0mil(k)*iaspmil(k)
               b3=dasin(deltamil(k))
c               b4=R0mil(k)*iaspmil(k)*kappamil(k)
c               b1=1.0d0;b2=0.32d0;b3=dasin(0.33d0);b4=1.7d0*b2
               sumy2=0.0d0

               do i=1,nin
               Z0=0.0d0 !symmetric
               b4=(ZatZmax-Z0)
               if (i.le.iZmax) then
                  theta(i)=dasin((Zin(i)-Z0)/(ZatZmax-Z0))
               else if (i.le.nin/2) then
                  theta(i)=pi-dasin((Zin(i)-Z0)/(ZatZmax-Z0))
               else if (i.lt.iZmin) then
                  theta(i)=pi-dasin((Zin(i)-Z0)/(ZatZmax-Z0))
               else if (i.eq.iZmin) then
                  theta(i)=pi*1.5d0
               else
                  theta(i)=2.0d0*pi+dasin((Zin(i)-Z0)/(ZatZmax-Z0))
               end if
c               datan2((Zin(i)-Z0),
c     1              (Rin(i)-b1))
               end do

            else
               fac=1.0d0
               b1=b1-fac*real(delb(1))
               b2=b2-fac*real(delb(2))
               b3=b3-fac*real(delb(3))
c               b4=b4-fac*real(delb(4))
               
            end if

         
c            write(*,*) 'theta',theta(1:nin)
c            write(*,*) 'b1,b2,b3,b4',b1,b2,b3,b4
            sumy1=sumy2
            sumy2=0.0d0
            do i=1,nin
               t0(i)=theta(i)+b3*dsin(theta(i))
               Rfit(i)=b1+b2*dcos(t0(i))
c               write(*,*) 'theta,t0',theta(i),t0(i)
c               Zfit(i)=b4*dsin(theta(i))
               dely(i)=((Rfit(i)-Rin(i))**2)/2.0d0 !(Zfit(i)-Zin(i))**2
c                write(*,*) 'R,z',Rfit(i),Rin(i),Zfit(i),Zin(i)
               sumy2=sumy2+dely(i)        
            end do
c            write(*,*) 'dely',dely(1:nin)
c            write(*,*) 'sum2,sum1',sumy2,sumy1
            dsumy=sumy2-sumy1
            if (dabs(dsumy).lt.epsFIT) then
               write(*,*) 'finished iteration for fitting 
     1             to miller parameters when iiterf=', iiterf
               exit  !exit iteration for iiterf
            end if
            Omega(1:3)=0.0d0
            Hessian=0.0d0
            do i=1,nin

               dRdbeta1=1.0d0
               dRdbeta2=dcos(t0(i))
               dRdbeta3=-b2*dsin(t0(i))*dsin(theta(i))

               Omega(1)=Omega(1)+dRdbeta1*(Rfit(i)-Rin(i))
               Omega(2)=Omega(2)+dRdbeta2*(Rfit(i)-Rin(i))
               Omega(3)=Omega(3)+dRdbeta3*(Rfit(i)-Rin(i))
               d2Rdb2b3=-dsin(t0(i))*dsin(theta(i))
               d2Rdb3=-b2*dcos(t0(i))*dsin(theta(i))**2
               Hessian(1,1)=Hessian(1,1)+dRdbeta1*dRdbeta1
               Hessian(1,2)=Hessian(1,2)+dRdbeta1*dRdbeta2
               Hessian(1,3)=Hessian(1,3)+dRdbeta1*dRdbeta3
               Hessian(2,3)=Hessian(2,3)+dRdbeta2*dRdbeta3
     1              +d2Rdb2b3*(Rfit(i)-Rin(i))
               Hessian(2,2)=Hessian(2,2)+dRdbeta2*dRdbeta2
               Hessian(3,3)=Hessian(3,3)+dRdbeta3*dRdbeta3
     1              +d2Rdb3*(Rfit(i)-Rin(i))
            end do
            Hessian(2,1)=Hessian(1,2)
            Hessian(3,1)=Hessian(1,3)
            Hessian(3,2)=Hessian(2,3)

            regulp=1.0d-2 ! regularization parameter to reduce condition number
            do j1=1,3
               do j2=1,3
                  A(j1,j2)=(1.0d0,0.0d0)*Hessian(j1,j2)
               end do
               y(j1)=(1.0d0,0.0d0)*Omega(j1)
               ! Tikhonov's regularization
c               JTJ(j1,j1)=JTJ(j1,j1)+regulp*(1.0d0,0.0d0)
c               write(*,*) 'jTj',   JTJ(j1,1:4)
            end do
            call cqrdecom(A,3,w7,rcond)
c            write(*,*) 'condition number', rcond
            call cqrsolve(3,w7,y,delb)  
c            write(*,*) 'iiterf,delb', iiterf,real(delb(1:3))
         end do

         R0mil(k)=b1
         iaspmil(k)=b2/b1
         kappamil(k)=b4/b2
         deltamil(k)=dsin(b3)
      end do

      return
      end 

      subroutine calfint(ncon,psicon,Rcon,Zcon,tmaxp,icon
     1    ,psiRcon,psiZcon,Ia,Ib,Ic,Lp,Id)

      use arrays, only:  nt2, ksamp2

      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*),psiRcon(*),psiZcon(*)
      real *8 Rcon(*), Zcon(*),tmaxp(*)
      real *8 psiRin(nt2*ksamp2+1),psiZin(nt2*ksamp2+1)
      real *8 Rin(nt2*ksamp2+1), Zin(nt2*ksamp2+1),thin(nt2*ksamp2+1)
      real *8 Ia(*),Ib(*),Ic(*),Lp(*),Id(*)
c      real *8 fint1(ncon),fint2(ncon), fint3(ncon), chcoeff(ncon)
c      real *8 dfint1(ncon), lambda

      integer *4  icon(*)
      integer i,j,k

c      ksamp3=ksamp2/4 
      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi
c      write(*,*) 'psicon1',psicon(1:ncon)
c      write(*,*) 'psiRintot',psiRcon(1:ncon*ksamp3)
      do k=1, ncon
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         Rin(1:nin)=Rcon(istart:iend)
         Zin(1:nin)=Zcon(istart:iend)
         thin(1:nin)=tmaxp(istart:iend)
         psiRin(1:nin)=psiRcon(istart:iend)
         psiZin(1:nin)=psiZcon(istart:iend)
         Rin(nin+1)=Rcon(istart)
         Zin(nin+1)=Zcon(istart)
         thin(nin+1)=tmaxp(istart)+2.0d0*pi
         psiRin(nin+1)=psiRcon(istart)
         psiZin(nin+1)=psiZcon(istart)

         call findsinglefint(nin,Rin, Zin
     1     ,thin, psiRin, psiZin
     2     ,Ia(k),Ib(k),Ic(k),Lp(k),Id(k))
      end do
c      write(*,*) 'psicon2',psicon(1:ncon)
      return
      end 

      subroutine findsinglefint(nin,Rin,Zin,thin,psiRin,psiZin
     1    ,fout1,fout2,fout3,fout4,fout5)

      use arrays, only:  nt2, ksamp2, epstri
      use arrays, only: isymud, Rmaxis, Zmaxis

      implicit real*8 (a-h,o-z)
c      save
c      parameter (nint=nt2)
      real *8 Rin(*), Zin(*),thin(*)
    
      real *8 Rint(nin+1),Zint(nin+1), dldt(nin+1)
 
      real *8 Rintmp(nin+1), Zintmp(nin+1),thintmp(nin+1)
      real *8 psiRin(*),psiZin(*), BpR2(nin+1), tttmp(nin)

      real *8 temp0(nin+1), temp1(nin+1), temp2(nin+1)
      real *8 temp3(nin+1), temp4(nin+1), temp5(nin+1)
      real *8 temp6(nin+1), temp7(nin+1), temp8(nin+1)
  
      real *8 thint(nin+1),integ11(nin+1)
      real *8 integ12(nin+1),integ13(nin+1)
      real *8 integ14(nin+1),integ15(nin+1)

      integer *4 i,j,k, kLag2, iw, itype, nint

      save
      
      nint=nin !nt2*ksamp2/2
      pi = 4.0d0*datan(1.0d0)
c      epstri= 1.0d-10

      delth = 0.0d0
      do i=1,nint
         thint(i)=2.0d0*pi*(i-1)/nint+thin(1) !delth
      end do


      kLag2=4
      if (isymud.eq.1) then
         do i=1,nin+1
            thintmp(i)=thin(i)
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
      else
         do i=1,nin+1
            thintmp(i)=thin(i)/2.0d0
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
         tttmp(1:nint)=thint(1:nint)/2.0d0
      end if
 

      BpR2(1:nin+1)=(psiRin(1:nin+1)**2.0d0+psiZin(1:nin+1)**2.0d0)
!      temp0(1:nin+1)=1.0d0/
!     1     (((Zin(1:nin+1)-Zmaxis)/(Rin(1:nin+1)-Rmaxis))**2.0d0+1.0d0)
!      temp1(1:nin+1)=-temp0(1:nin+1)*(Zin(1:nin+1)-Zmaxis)
!     1     /(Rin(1:nin+1)-Rmaxis)**2.0d0   !dtheta/dR
!      temp2(1:nin+1)=temp0(1:nin+1)
!     1     /(Rin(1:nin+1)-Rmaxis)   !dtheta/dZ


      temp0(1:nin+1)=(Rin(1:nin+1)-Rmaxis)**2.0d0/
     1     ((Zin(1:nin+1)-Zmaxis)**2.0d0+(Rin(1:nin+1)-Rmaxis)**2.0d0)
      temp1(1:nin+1)=-(Zin(1:nin+1)-Zmaxis)/
     1     ((Zin(1:nin+1)-Zmaxis)**2.0d0+(Rin(1:nin+1)-Rmaxis)**2.0d0)   !dtheta/dR
      temp2(1:nin+1)=(Rin(1:nin+1)-Rmaxis)/
     1     ((Zin(1:nin+1)-Zmaxis)**2.0d0+(Rin(1:nin+1)-Rmaxis)**2.0d0)  !dtheta/dZ

      temp3(1:nin+1)=dabs((psiRin(1:nin+1)*temp2(1:nin+1)
     1     -psiZin(1:nin+1)*temp1(1:nin+1))/Rin(1:nin+1)) !(grad phi times grad psi) cdot grad theta)

      temp4(1:nin+1)=BpR2(1:nin+1)/temp3(1:nin+1)/Rin(1:nin+1)**2 

c      write(*,*) 'temp4 ',temp4(1:nin+1),thintmp(1:nin/2+1)
      temp5(1:nin+1)=1.0d0/temp3(1:nin+1)
      temp6(1:nin+1)=1.0d0/temp3(1:nin+1)/Rin(1:nin+1)**2 
      temp7(1:nin+1)=dsqrt(BpR2(1:nin+1))/temp3(1:nin+1)/Rin(1:nin+1)
      temp8(1:nin+1)=1.0d0/temp3(1:nin+1)/Rin(1:nin+1)

      if (isymud.eq.1) then

         call IntTriCos(nin/2+1,thintmp,temp4,nint/2,
     1        thint,integ11,epstri)
c     write(*,*) 'intlag2'
         call IntTriCos(nin/2+1,thintmp,temp5,nint/2,
     1        thint,integ12,epstri)

c     write(*,*) 'intlag3'
         call IntTriCos(nin/2+1,thintmp,temp6,nint/2,
     1        thint,integ13,epstri) 
c     write(*,*) 'intlag4'
         call IntTriCos(nin/2+1,thintmp,temp7,nint/2,
     1        thint,integ14,epstri)

c     write(*,*) 'intlag5'
         call IntTriCos(nin/2+1,thintmp,temp8,nint/2,
     1        thint,integ15,epstri)

         integ11(nint/2+1)=temp4(nin/2+1)
         integ12(nint/2+1)=temp5(nin/2+1)
         integ13(nint/2+1)=temp6(nin/2+1)
         integ14(nint/2+1)=temp7(nin/2+1)
         integ15(nint/2+1)=temp8(nin/2+1)

         do i=1,nint/2
            integ11(nint+2-i)=integ11(i)
            integ12(nint+2-i)=integ12(i)
            integ13(nint+2-i)=integ13(i)
            integ14(nint+2-i)=integ14(i)
            integ15(nint+2-i)=integ15(i)
         end do
c         write(*,*) 'integ11t2',integ11(1:nint+1)
      else                      ! if isymud is not 1 (up-down assymmetric)
        
         call IntTriCos(nin+1,thintmp,temp4,nint+1,
     1        tttmp,integ11,epstri)
  
c     write(*,*) 'intlag2'
         call IntTriCos(nin+1,thintmp,temp5,nint+1,
     1        tttmp,integ12,epstri)

c      write(*,*) 'intlag3'
         call IntTriCos(nin+1,thintmp,temp6,nint+1,
     1        tttmp,integ13,epstri)
c      write(*,*) 'intlag4'
         call IntTriCos(nin+1,thintmp,temp7,nint+1,
     1        tttmp,integ14,epstri)

c      write(*,*) 'intlag5'
         call IntTriCos(nin+1,thintmp,temp8,nint+1,
     1        tttmp,integ15,epstri)

      end if                    ! end if isymud.eq.1

      fout1=0.0d0
      fout2=0.0d0
      fout3=0.0d0
      fout4=0.0d0
      fout5=0.0d0
      do i=1,nint
         fout1=fout1+integ11(i)
         fout2=fout2+integ12(i)  
         fout3=fout3+integ13(i)  
         fout4=fout4+integ14(i) 
         fout5=fout5+integ15(i) 
      end do
      fout1=fout1*2.0d0*pi/nint
      fout2=fout2*2.0d0*pi/nint
      fout3=fout3*2.0d0*pi/nint
      fout4=fout4*2.0d0*pi/nint
      fout5=fout5*2.0d0*pi/nint
c      write(*,*) 'integ11',integ11(1:nin+1)
c      write(*,*) 'integ12',integ12(1:nin+1)
c      write(*,*) 'integ13',integ13(1:nin+1)
c      write(*,*) 'fout2',fout2
      return 

      end

      subroutine calfint2(ncon,psicon,Rcon,Zcon,tmaxp,icon
     1    ,psiRcon,psiZcon,Ie,Ig,Ih)

      use arrays, only:  nt2, ksamp2

      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*),psiRcon(*),psiZcon(*)
      real *8 Rcon(*), Zcon(*),tmaxp(*)
      real *8 psiRin(nt2*ksamp2+1),psiZin(nt2*ksamp2+1)
      real *8 Rin(nt2*ksamp2+1), Zin(nt2*ksamp2+1),thin(nt2*ksamp2+1)
      real *8 Ie(*),Ig(*),Ih(*)
c      real *8 fint1(ncon),fint2(ncon), fint3(ncon), chcoeff(ncon)
c      real *8 dfint1(ncon), lambda

      integer *4  icon(*)
      integer i,j,k

c      ksamp3=ksamp2/4 
      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi
c      write(*,*) 'psicon1',psicon(1:ncon)
c      write(*,*) 'psiRintot',psiRcon(1:ncon*ksamp3)
      do k=1, ncon
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         Rin(1:nin)=Rcon(istart:iend)
         Zin(1:nin)=Zcon(istart:iend)
         thin(1:nin)=tmaxp(istart:iend)
         psiRin(1:nin)=psiRcon(istart:iend)
         psiZin(1:nin)=psiZcon(istart:iend)
         Rin(nin+1)=Rcon(istart)
         Zin(nin+1)=Zcon(istart)
         thin(nin+1)=tmaxp(istart)+2.0d0*pi
         psiRin(nin+1)=psiRcon(istart)
         psiZin(nin+1)=psiZcon(istart)

         call findsinglefint2(nin,Rin, Zin
     1     ,thin, psiRin, psiZin
     2     ,Ie(k),Ig(k),Ih(k))
      end do
c      write(*,*) 'psicon2',psicon(1:ncon)
      return
      end 

      subroutine findsinglefint2(nin,Rin,Zin,thin,psiRin,psiZin
     1    ,fout6,fout7,fout8)

      use arrays, only:  nt2, ksamp2, epstri
      use arrays, only: isymud, Rmaxis, Zmaxis


      implicit real*8 (a-h,o-z)
c      save
c      parameter (nint=nt2)
      real *8 Rin(*), Zin(*),thin(*)
    
      real *8 Rint(nin+1),Zint(nin+1), dldt(nin+1)
 
      real *8 Rintmp(nin+1), Zintmp(nin+1),thintmp(nin+1)
      real *8 psiRin(*),psiZin(*), BpR2(nin+1), tttmp(nin)

      real *8 temp0(nin+1), temp1(nin+1), temp2(nin+1)
      real *8 temp3(nin+1), temp4(nin+1), temp5(nin+1)
      real *8 temp6(nin+1), temp7(nin+1), temp8(nin+1)
  
      real *8 thint(nin+1),integ11(nin+1)
      real *8 integ12(nin+1),integ13(nin+1)
      real *8 integ14(nin+1),integ15(nin+1)
      real *8 integ16(nin+1),integ17(nin+1)
      real *8 integ18(nin+1),integ19(nin+1)
      integer *4 i,j,k, kLag2, iw, itype, nint

      save
      
      nint=nin !nt2*ksamp2/2
      pi = 4.0d0*datan(1.0d0)
c      epstri= 1.0d-10

      delth = 0.0d0
      do i=1,nint
         thint(i)=2.0d0*pi*(i-1)/nint+thin(1) !delth
      end do


      kLag2=4
      if (isymud.eq.1) then
         do i=1,nin+1
            thintmp(i)=thin(i)
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
      else
         do i=1,nin+1
            thintmp(i)=thin(i)/2.0d0
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
         tttmp(1:nint)=thint(1:nint)/2.0d0
      end if
 

      BpR2(1:nin+1)=(psiRin(1:nin+1)**2.0d0+psiZin(1:nin+1)**2.0d0)
      temp0(1:nin+1)=1.0d0/
     1     (((Zin(1:nin+1)-Zmaxis)/(Rin(1:nin+1)-Rmaxis))**2.0d0+1.0d0)
      temp1(1:nin+1)=-temp0(1:nin+1)*(Zin(1:nin+1)-Zmaxis)
     1     /(Rin(1:nin+1)-Rmaxis)**2.0d0   !dtheta/dR
      temp2(1:nin+1)=temp0(1:nin+1)
     1     /(Rin(1:nin+1)-Rmaxis)   !dtheta/dZ

      temp3(1:nin+1)=dabs((psiRin(1:nin+1)*temp2(1:nin+1)
     1     -psiZin(1:nin+1)*temp1(1:nin+1))/Rin(1:nin+1)) !(grad phi times grad psi) cdot grad theta)

      temp4(1:nin+1)=1.0d0/BpR2(1:nin+1)/temp3(1:nin+1)/Rin(1:nin+1)**2 
      temp5(1:nin+1)=1.0d0/BpR2(1:nin+1)/temp3(1:nin+1)
      temp6(1:nin+1)=1.0d0/BpR2(1:nin+1)/temp3(1:nin+1)*Rin(1:nin+1)**2 

      if (isymud.eq.1) then

         call IntTriCos(nin/2+1,thintmp,temp4,nint/2,
     1        thint,integ16,epstri)
  
         call IntTriCos(nin/2+1,thintmp,temp5,nint/2,
     1        thint,integ17,epstri)

         call IntTriCos(nin/2+1,thintmp,temp6,nint/2,
     1        thint,integ18,epstri) 

         integ16(nint/2+1)=temp4(nin/2+1)
         integ17(nint/2+1)=temp5(nin/2+1)
         integ18(nint/2+1)=temp6(nin/2+1)

         do i=1,nint/2
            integ16(nint+2-i)=integ16(i)
            integ17(nint+2-i)=integ17(i)
            integ18(nint+2-i)=integ18(i)
         end do
      else                      ! if isymud is not 1 (up-down assymmetric)
        
         call IntTriCos(nin+1,thintmp,temp4,nint+1,
     1        tttmp,integ16,epstri)
  
         call IntTriCos(nin+1,thintmp,temp5,nint+1,
     1        tttmp,integ17,epstri)

         call IntTriCos(nin+1,thintmp,temp6,nint+1,
     1        tttmp,integ18,epstri)

      end if                    ! end if isymud.eq.1

      fout6=0.0d0
      fout7=0.0d0
      fout8=0.0d0
      do i=1,nint
         fout6=fout6+integ16(i)
         fout7=fout7+integ17(i)  
         fout8=fout8+integ18(i)  
      end do
      fout6=fout6*2.0d0*pi/nint
      fout7=fout7*2.0d0*pi/nint
      fout8=fout8*2.0d0*pi/nint
    
      return 

      end

      subroutine calfinttf(ncon,psicon,Rcon,Zcon,tmaxp,icon
     1    ,psiRcon,psiZcon,Itf1, Itf2)

      use arrays, only:  nt2, ksamp2
      use arrays, only:  itftype, presch, ptflowch, mach

      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*),psiRcon(*),psiZcon(*)
      real *8 Rcon(*), Zcon(*),tmaxp(*)
      real *8 psiRin(nt2*ksamp2+1),psiZin(nt2*ksamp2+1)
      real *8 Rin(nt2*ksamp2+1), Zin(nt2*ksamp2+1),thin(nt2*ksamp2+1)
      real *8 Itf1(*),Itf2(*)
c      real *8 fint1(ncon),fint2(ncon), fint3(ncon), chcoeff(ncon)
c      real *8 dfint1(ncon), lambda

      integer *4  icon(*)
      integer i,j,k

c      ksamp3=ksamp2/4 
      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi
c      write(*,*) 'psicon1',psicon(1:ncon)
c      write(*,*) 'psiRintot',psiRcon(1:ncon*ksamp3)
      do k=1, ncon
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         Rin(1:nin)=Rcon(istart:iend)
         Zin(1:nin)=Zcon(istart:iend)
         thin(1:nin)=tmaxp(istart:iend)
         psiRin(1:nin)=psiRcon(istart:iend)
         psiZin(1:nin)=psiZcon(istart:iend)
         Rin(nin+1)=Rcon(istart)
         Zin(nin+1)=Zcon(istart)
         thin(nin+1)=tmaxp(istart)+2.0d0*pi
         psiRin(nin+1)=psiRcon(istart)
         psiZin(nin+1)=psiZcon(istart)

         if (itftype.eq.1) then
            ptfbyp=mach**2/(2.0d0)
         else 
            ptfbyp=ptflowch(k)/presch(k)
         end if
c         write(*,*) 'ptfbyp', ptfbyp
         call findsinglefinttf(nin,Rin, Zin
     1     ,thin, psiRin, psiZin
     2     ,ptfbyp ,Itf1(k), Itf2(k))
      end do
c      write(*,*) 'psicon2',psicon(1:ncon)
      return
      end 

      subroutine findsinglefinttf(nin,Rin,Zin,thin,psiRin,psiZin
     1    , ptfbyp, fouttf1, fouttf2)

      use arrays, only:  nt2, ksamp2, epstri
      use arrays, only: isymud, Rmid, Rmaxis, Zmaxis

      implicit real*8 (a-h,o-z)
c      save
c      parameter (nint=nt2)
      real *8 Rin(*), Zin(*),thin(*)
    
      real *8 Rint(nin+1),Zint(nin+1), dldt(nin+1)
 
      real *8 Rintmp(nin+1), Zintmp(nin+1),thintmp(nin+1)
      real *8 psiRin(*),psiZin(*), pexp(nin+1), tttmp(nin)

      real *8 temp0(nin+1), temp1(nin+1), temp2(nin+1)
      real *8 temp3(nin+1), temp4(nin+1), temp5(nin+1)
      real *8 temp6(nin+1), temp7(nin+1), temp8(nin+1)
  
      real *8 thint(nin+1),integ1tf(nin+1),integ2tf(nin+1)

      integer *4 i,j,k, kLag2, iw, itype, nint

      save
      
      nint=nin !nt2*ksamp2/2
      pi = 4.0d0*datan(1.0d0)
c      epstri= 1.0d-10

      delth = 0.0d0
      do i=1,nint
         thint(i)=2.0d0*pi*(i-1)/nint+thin(1) !delth
      end do


      kLag2=4
      if (isymud.eq.1) then
         do i=1,nin+1
            thintmp(i)=thin(i)
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
      else
         do i=1,nin+1
            thintmp(i)=thin(i)/2.0d0
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
         tttmp(1:nint)=thint(1:nint)/2.0d0
      end if
 
      pexp(1:nin+1)=dexp(ptfbyp*((Rin(1:nin+1)/Rmid)**2-1)) !exponential factor
      temp0(1:nin+1)=1.0d0/
     1     (((Zin(1:nin+1)-Zmaxis)/(Rin(1:nin+1)-Rmaxis))**2.0d0+1.0d0)
      temp1(1:nin+1)=-temp0(1:nin+1)*(Zin(1:nin+1)-Zmaxis)
     1     /(Rin(1:nin+1)-Rmaxis)**2.0d0   !dtheta/dR
      temp2(1:nin+1)=temp0(1:nin+1)
     1     /(Rin(1:nin+1)-Rmaxis)   !dtheta/dZ

      temp3(1:nin+1)=dabs((psiRin(1:nin+1)*temp2(1:nin+1)
     1     -psiZin(1:nin+1)*temp1(1:nin+1))/Rin(1:nin+1)) !(grad phi times grad psi) cdot grad theta)

      temp4(1:nin+1)=pexp(1:nin+1)/temp3(1:nin+1)
      temp5(1:nin+1)=pexp(1:nin+1)*((Rin(1:nin+1)/Rmid)**2-1)
     1     /temp3(1:nin+1)
     
      if (isymud.eq.1) then

         call IntTriCos(nin/2+1,thintmp,temp4,nint/2+1,
     1     thint,integ1tf,epstri)
         call IntTriCos(nin/2+1,thintmp,temp5,nint/2+1,
     1     thint,integ2tf,epstri)

         do i=1,nint/2
            integ1tf(nint+2-i)=integ1tf(i)
            integ2tf(nint+2-i)=integ2tf(i)
         end do
      else                      ! if isymud is not 1 (up-down assymmetric)
        
         call IntTriCos(nin+1,thintmp,temp4,nint+1,
     1        tttmp,integ1tf,epstri)
         call IntTriCos(nin+1,thintmp,temp5,nint+1,
     1        tttmp,integ2tf,epstri)

      end if                    ! end if isymud.eq.1

      fouttf1=0.0d0
      fouttf2=0.0d0
      do i=1,nint
         fouttf1=fouttf1+integ1tf(i)  
         fouttf2=fouttf2+integ2tf(i)  
      end do
      fouttf1=fouttf1*2.0d0*pi/nint
      fouttf2=fouttf2*2.0d0*pi/nint

c      write(*,*) 'integtf',integ1tf(1:nin+1)
c      write(*,*) 'fouttf1',fouttf1
      return 

      end


      subroutine calmetric(ncon,psicon,Rcon,Zcon,tmaxp,icon,psiRcon
     1        ,psiZcon,fpolcon,dfdpsi,dpdpsi,dBdpsi,dalphadpsi0
     2     ,nttin,ttin,gradpar,Rtrin,Ztrin,bpR,btot,bpol,gbdrift
     3     ,gbdrift0,cvdrift,cvdrift0,tshatlocal,dalphatt,phitrin)

      use arrays, only:  nt2, ksamp2, dpsidrho,epsLag, kLag

      implicit real*8 (a-h,o-z)
      save

      parameter (ndt=8)

      real *8 psicon(*),psiRcon(*),psiZcon(*), BpR(*)
      real *8 Rcon(*), Zcon(*),tmaxp(*),ttin(*), Rtrin(*), Ztrin(*)
      real *8 psiRin(nt2*ksamp2+1),psiZin(nt2*ksamp2+1)
      real *8 Rin(nt2*ksamp2+1), Zin(nt2*ksamp2+1),thin(nt2*ksamp2+1)
      real *8 fpolcon(*),dfdpsi(*),dpdpsi(*),dBdpsi(*),dalphadpsi0(*)
      real *8 gradpar(*), btot(*),bpol(*), gbdrift(*),gbdrift0(*)
      real *8 cvdrift(*),cvdrift0(*),tshatlocal(*),dalphatt(*)
      real *8 phitrin(*)
      real *8 sgradpar(nttin),sbtot(nttin),sbpR(nttin)
      real *8 sRtrin(nttin),sZtrin(nttin),sbpol(nttin),salphapsi0(nttin)
      real *8 spsidott(nttin),sdBdtheta(nttin),sdalphat(nttin)
      real *8 dalphadpsi(nttin), ddalphadpsidt(nttin),temp(nttin)

      real *8 alphapsi0(nttin*ncon)
      real *8 psidott(nttin*ncon),dBdtheta(nttin*ncon)

      real *8 tn(ndt),cftm(ndt*ndt)
      real *8 spbx(ndt*ndt), spxb(ndt*ndt), spdef(ndt)
c      real *8 fint1(ncon),fint2(ncon), fint3(ncon), chcoeff(ncon)
c      real *8 dfint1(ncon), lambda
      real *8 pi, mu
      integer *4  icon(*)
      integer i,j,k

c      ksamp3=ksamp2/4 
      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi

      !to calculate integral from 0 to theta
      call chsetupq(ndt, tn, cftm, spbx, spxb, spdef) 
      tn(1:ndt)=1.0d0-tn(1:ndt) 
c      write(*,*) 'tn',tn(1:ncon)
c      write(*,*) 'psiRintot',psiRcon(1:ncon*ksamp3)
      do k=1, ncon
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         Rin(1:nin)=Rcon(istart:iend)
         Zin(1:nin)=Zcon(istart:iend)
         thin(1:nin)=tmaxp(istart:iend)
         psiRin(1:nin)=psiRcon(istart:iend)
         psiZin(1:nin)=psiZcon(istart:iend)
         Rin(nin+1)=Rcon(istart)
         Zin(nin+1)=Zcon(istart)
         thin(nin+1)=tmaxp(istart)+2.0d0*pi
         psiRin(nin+1)=psiRcon(istart)
         psiZin(nin+1)=psiZcon(istart)

         call findsinglemetric(nin,Rin, Zin
     1     ,thin, psiRin, psiZin, fpolcon(k)
     2     ,nttin, ttin, ndt, tn, spdef, cftm
     3     ,sgradpar, sRtrin, sZtrin,sbpR,sbtot,sbpol
     4     ,salphapsi0,sdBdtheta,spsidott,sdalphat)

         Rtrin((k-1)*nttin+1:k*nttin)=sRtrin(1:nttin)
         Ztrin((k-1)*nttin+1:k*nttin)=sZtrin(1:nttin)
         BpR((k-1)*nttin+1:k*nttin)=sbpR(1:nttin)
         bpol((k-1)*nttin+1:k*nttin)=sbpol(1:nttin)
         gradpar((k-1)*nttin+1:k*nttin)=sgradpar(1:nttin)
         btot((k-1)*nttin+1:k*nttin)=sbtot(1:nttin)
         alphapsi0((k-1)*nttin+1:k*nttin)=salphapsi0(1:nttin)
         dBdtheta((k-1)*nttin+1:k*nttin)=sdBdtheta(1:nttin)
         psidott((k-1)*nttin+1:k*nttin)=spsidott(1:nttin)
         dalphatt((k-1)*nttin+1:k*nttin)=sdalphat(1:nttin)

      end do
      do k=1, ncon
         do j=1,nttin
            kj=(k-1)*nttin+j
            if (isNaN(gradpar(kj))) gradpar(kj)=0.0 !JPL
            if (isNaN(alphapsi0(kj))) alphapsi0(kj)=0.0 !JPL
            if (isNaN(dalphadpsi0(kj))) dalphadpsi0(kj)=0.0 !JPL           
            gbdriftphi=-bpol(kj)*dBdpsi(kj)
     1           -(psidott(kj)*dBdtheta(kj)/Rtrin(kj)**2)
            gbdrifttheta=fpolcon(k)*dBdpsi(kj)*(gradpar(kj)*btot(kj))
            gbdriftpsi=-fpolcon(k)*dBdtheta(kj)*(gradpar(kj)*btot(kj))

            alphaphi=1.0d0
            alphatheta=-fpolcon(k)/(gradpar(kj)*btot(kj)*Rtrin(kj)**2)
            alphapsi=-dfdpsi(k)*alphapsi0(kj)-fpolcon(k)*dalphadpsi0(kj)
   
            phitrin(kj)=alphapsi0(kj)*abs(fpolcon(k))/(2.0d0*pi) !JPL
            if (isNaN(phitrin(kj))) phitrin(kj)=0.0 
            dalphadpsi(j)=alphapsi
c           write(*,*) 'alpahpsi',dfdpsi(k),alphapsi0(kj),fpolcon(k)
c     1           ,dalphadpsi0(kj)
            gbdrift(kj)=(gbdriftphi*alphaphi+gbdrifttheta*alphatheta
     1           +gbdriftpsi*alphapsi)/dpsidrho(k)
            gbdrift0(kj)=gbdriftpsi/dpsidrho(k)

c            write(*,*) 'gbdrift',gbdriftphi,alphaphi,gbdrifttheta
c     1          ,alphatheta ,gbdriftpsi,alphapsi

            cvdriftphi=-bpol(kj)*(dBdpsi(kj)+mu*dpdpsi(k)/btot(kj))
     1           -(psidott(kj)*dBdtheta(kj)/Rtrin(kj)**2)
            cvdrifttheta=fpolcon(k)*(dBdpsi(kj)+mu*dpdpsi(k)/btot(kj))
     1           *(gradpar(kj)*btot(kj))
            cvdriftpsi=-fpolcon(k)*dBdtheta(kj)*(gradpar(kj)*btot(kj))

            cvdrift(kj)=(cvdriftphi*alphaphi+cvdrifttheta*alphatheta
     1           +cvdriftpsi*alphapsi)/dpsidrho(k)
            cvdrift0(kj)=cvdriftpsi/dpsidrho(k)
c            write(*,*) 'cvdrift',cvdriftphi,alphaphi,cvdrifttheta
c     1          ,alphatheta ,cvdriftpsi,alphapsi

         end do

         !write(*,*) 'phitrin',k,phitrin((k-1)*nttin+1:(k)*nttin)
c         write(*,*) 'temp5lag', dalphadpsi(1:nttin)

         call IntLag1np(nttin,ttin,dalphadpsi,kLag,nttin,ttin,
     1        temp,ddalphadpsidt,epsLag)   !alpha is not periodic 
c         write(*,*) 'temp5lagder',ddalphadpsidt(1:nttin)
c         call findddalphadpsidt(nin,thin,nttin, ttin
c     1     , ndt, tn, spdef, cftm
c     2    ,dalphadpsi ,ddalphadpsidt)

         tshatlocal((k-1)*nttin+1:k*nttin)=ddalphadpsidt(1:nttin)
     1     *(gradpar((k-1)*nttin+1:k*nttin)*btot((k-1)*nttin+1:k*nttin))
      end do


      return
      end 


      subroutine findsinglemetric(nin,Rin,Zin
     1     ,thin,psiRin,psiZin, fpolin
     2     ,nttin, ttin, ndt, tn, spdef, cftm
     3      ,sgradpar, sR,sZ,sBpR,sbtot,sbpol
     4     ,alphapsi0,dBdtheta,spsidott, sdalphat)

      use arrays, only: nt2, ksamp2, epstri
      use arrays, only: isymud, Rmid, lambda, Rmaxis, Zmaxis

      implicit real*8 (a-h,o-z)
c      save
      real *8 Rin(*), Zin(*),thin(*)
    
      real *8 Rint(nin+1),dRdt(nin+1) ,Zint(nin+1)
      real *8 dZdt(nin+1),dldt(nin+1), ttintmp(nttin+1)
 
      real *8 Rintmp(nin+1), Zintmp(nin+1),thintmp(nin+1)
      real *8 psiRin(*),psiZin(*), BpR2(nin+1),BR2(nin+1),BpR(nin+1)
      real *8 temp0(nin+1), temp1(nin+1), temp2(nin+1)
      real *8 temp3(nin+1), temp4(nin+1), temp5(nin+1)
      real *8 temp6(nin+1), temp7(nin+1), temp8(nin+1), temp9(nin+1)
      real *8 gradparin(nin+1), Jac(nin+1)

      real *8 ttin(*),sgradpar(*), sbtot(*)
      real *8 alphapsi0(*),dBdtheta(*),spsidott(*),sdalphat(*)

      real *8 tn(*), cftm(*), spdef(*),dttin(nttin),tfine1(nttin*ndt)
      real *8 tfine2(nttin*ndt)
      real *8 tempfine5(nttin*ndt), tempfine8(nttin*ndt), chcoeff5(ndt)
      real *8 tempfine9(nttin*ndt),chcoeff9(ndt)

      real *8 sR(*),sZ(*),sBpR(*),sbpol(*)

      integer *4 i,j,k, kLag2, iw, itype, nint

      save

      !psi in this routine is normalized one to the range from 0 to 1
      nint=nin !nt2*ksamp2/2
      pi = 4.0d0*datan(1.0d0)

      do i=1,nttin ! assumed ttin(1)=0
         if (i.ne.nttin) then
            dttin(i)=ttin(i+1)-ttin(i)
         else
            dttin(i)=ttin(1)+2.0d0*pi-ttin(nttin)
         end if

         tfine1((i-1)*ndt+1:i*ndt)=ttin(i)+tn(1:ndt)*dttin(i)
c         tfine2((i-1)*ndt+1:i*ndt)=ttin(i)+tn(1:ndt)*dttin(i)
c     1        -dttin(i)/2.0d0
      end do
      
c      write(*,*) 'tfine1',tfine1(1:nttin*ndt)
c      do j=1,ndt
c         if (tfine2(j).lt.0.0d0) then
c            tfine2(j)=tfine2(j)+2.0d0*pi
c         end if
c      end do
     
      kLag2=4
      if (isymud.eq.1) then
         do i=1,nin+1
            thintmp(i)=thin(i)
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
      else
         do i=1,nin+1
            thintmp(i)=thin(i)/2.0d0
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
         tfine1(1:nttin*ndt)=tfine1(1:nttin*ndt)/2.0d0
         ttintmp(1:nttin)=ttin(1:nttin)/2.0d0
c         write(*,*) 'ttintmp',ttintmp(1:nttin),ttin(1:nttin)
      end if
 
      BpR2(1:nin+1)=(psiRin(1:nin+1)**2.0d0+psiZin(1:nin+1)**2.0d0)
     1     /lambda**2           

      BR2(1:nin+1)= BpR2(1:nin+1)+fpolin**2.0d0
      temp0(1:nin+1)=1.0d0/
     1     (((Zin(1:nin+1)-Zmaxis)/(Rin(1:nin+1)-Rmaxis))**2.0d0+1.0d0)
      temp1(1:nin+1)=-temp0(1:nin+1)*(Zin(1:nin+1)-Zmaxis)
     1     /(Rin(1:nin+1)-Rmaxis)**2.0d0   !dtheta/dR
      temp2(1:nin+1)=temp0(1:nin+1)
     1     /(Rin(1:nin+1)-Rmaxis)   !dtheta/dZ
      temp3(1:nin+1)=(psiRin(1:nin+1)*temp2(1:nin+1)
     1     -psiZin(1:nin+1)*temp1(1:nin+1))/Rin(1:nin+1)/lambda !(grad phi times grad Psi) cdot grad theta)
      BpR(1:nin+1)=dsqrt(BpR2(1:nin+1)) !magnitude of BpR
      temp4(1:nin+1)=dsqrt(BpR2(1:nin+1))/Rin(1:nin+1) !magnitude of Bp
      temp5(1:nin+1)=dsqrt(BR2(1:nin+1))/Rin(1:nin+1) !magnitude of Btot

      temp6(1:nin+1)=temp3(1:nin+1)/temp5  ! b dot grad theta

      temp7(1:nin+1)=(psiRin(1:nin+1)*temp1(1:nin+1)
     1     +psiZin(1:nin+1)*temp2(1:nin+1))/lambda     !grad Psi cdot grad theta

c      write(*,*) 'psiRin, temp1-1',psiRin(1:nin+1)
c      write(*,*) 'psiRin, temp1-2',temp2(1:nin+1)
c      write(*,*) 'psiRin, temp1-31',psiZin(1:nin+1)
c      write(*,*) 'psiRin, temp1-32',temp1(1:nin+1)
      temp8(1:nin+1)=1.0d0/temp3(1:nin+1)
     1       /Rin(1:nin+1)**2                    ! 1/(B dot grad theta)/R^2

      do i=1,nin+1
         if (BpR2(i).eq.0) then
            temp9(i)=0.0d0
         else
            temp9(i)=fpolin*temp8(i)*temp7(i)/BpR2(i)
         end if
      end do

      
c      write(*,*) 'BpR23',(BpR2(1:nin+1))
c      write(*,*) 'Rin3',Rin(1:nin+1)
c      write(*,*) 'Bp3',dsqrt(BpR2(1:nin+1))/Rin(1:nin+1)
c      write(*,*) 'Btot3',temp1(1:nin+1),lambda

      if (isymud.eq.1) then
  
         call IntTriCos(nin/2+1,thintmp,Rin,nttin/2+1,
     1        ttin,sR,epstri)
c         write(*,*) 'Rin-sR1',thintmp(1:nin/2+1)
c         write(*,*) 'Rin-sR2',Rin(1:nin/2+1)
c         write(*,*) 'Rin-sR3',ttin(1:nttin/2+1)
c         write(*,*) 'Rin-sR4',sR(1:nttin/2+1)
         call IntTriCos(nin/2+1,thintmp,Zin,nttin/2+1,
     1        ttin,sZ,epstri)
c         write(*,*) 'Zin-sR1',thintmp(1:nin/2+1)
c         write(*,*) 'Zin-sR2',Zin(1:nin/2+1)
c         write(*,*) 'Zin-sR3',ttin(1:nttin/2+1)
c         write(*,*) 'Zin-sR4',sZ(1:nttin/2+1)
         call IntTriCos(nin/2+1,thintmp,BpR,nttin/2+1,
     1        ttin,sBpR,epstri) 
         call IntTriCos(nin/2+1,thintmp,temp4,nttin/2+1,
     1        ttin,sbpol,epstri) 
         call IntTriCos(nin/2+1,thintmp,temp5,nttin/2+1,
     1        ttin,sbtot,epstri) !assumed symmetric ttin
         
         call IntTriCos(nin/2+1,thintmp,temp6,nttin/2+1,
     1        ttin,sgradpar,epstri) 
         call IntTriCos(nin/2+1,thintmp,temp7,nttin/2+1,
     1        ttin,spsidott,epstri) 
         

         call IntTriCos(nin/2+1,thintmp,temp8,nttin*ndt/2+1,
     1        tfine1,tempfine8,epstri) 
         call IntTriCos(nin/2+1,thintmp,temp5,nttin*ndt/2+1,
     1        tfine1,tempfine5,epstri) 
         call IntTriCos(nin/2+1,thintmp,temp9,nttin*ndt/2+1,
     1        tfine1,tempfine9,epstri) 

         do i=2,nttin/2
            sR(nttin+2-i)=sR(i)
            sZ(nttin+2-i)=-sZ(i)
            sbpR(nttin+2-i)=sbpR(i)
            sbpol(nttin+2-i)=sbpol(i)
            sbtot(nttin+2-i)=sbtot(i)
            sgradpar(nttin+2-i)=sgradpar(i)
            spsidott(nttin+2-i)=spsidott(i)
         end do

         do i=2,nttin*ndt/2+1
            tempfine5(nttin*ndt+2-i)=tempfine5(i)
            tempfine8(nttin*ndt+2-i)=tempfine8(i)
            tempfine9(nttin*ndt+2-i)=tempfine9(i)
         end do

         do i=1,nttin/2
            call chftransq(chcoeff5,tempfine5((i-1)*ndt+1), ndt, cftm)
            call chderiv(0.0d0, ndt, chcoeff5,dtempfine5)
            dBdtheta(i)=dtempfine5/dttin(i)
         end do
         call chftransq(chcoeff5,tempfine5((nttin/2-1)*ndt+1),ndt, cftm)
         call chderiv(1.0d0, ndt, chcoeff5,dtempfine5)
         dBdtheta(nttin/2+1)=dtempfine5/dttin(i)
         do i=2,nttin/2
            dBdtheta(nttin+2-i)=-dBdtheta(i)
         end do

         do i=1,nttin/2
            call chftransq(chcoeff9,tempfine9((i-1)*ndt+1), ndt, cftm)
            call chderiv(0.0d0, ndt, chcoeff9,dtempfine9)
            sdalphat(i)=dtempfine9/dttin(i)
         end do
         call chftransq(chcoeff9,tempfine9((nttin/2-1)*ndt+1),ndt, cftm)
         call chderiv(1.0d0, ndt, chcoeff9,dtempfine9)
         sdalphat(nttin/2+1)=dtempfine9/dttin(i)
         do i=2,nttin/2
            sdalphat(nttin+2-i)=sdalphat(i)
         end do
         sdalphat(1:nttin)=-sdalphat(1:nttin)
c         write(*,*) 'temp9',(temp9(1:nin+1))
c         write(*,*) 'tempfine9',(tempfine9(1:nttin))
c         write(*,*) 'sdalphat',sdalphat(1:nttin)
      else  ! if isymud is not 1 (up-down assymmetric)

         call IntTriCos(nin+1,thintmp,Rin,nttin,
     1        ttintmp,sR,epstri)
         call IntTriCos(nin+1,thintmp,Zin,nttin,
     1        ttintmp,sZ,epstri)
c         write(*,*) 'Rin-sR1',thintmp(1:nin+1)
c         write(*,*) 'Rin-sR2',Rin(1:nin+1)
c         write(*,*) 'Rin-sR3',ttintmp(1:nttin)
c         write(*,*) 'Rin-sR4',sR(1:nttin)

c         write(*,*) 'Zin-sR1',thintmp(1:nin+1)
c         write(*,*) 'Zin-sR2',Zin(1:nin+1)
c         write(*,*) 'Zin-sR3',ttintmp(1:nttin)
c         write(*,*) 'Zin-sR4',sZ(1:nttin)

         call IntTriCos(nin+1,thintmp,BpR,nttin,
     1        ttintmp,sBpR,epstri) 
         call IntTriCos(nin+1,thintmp,temp4,nttin,
     1        ttintmp,sbpol,epstri) 
         call IntTriCos(nin+1,thintmp,temp5,nttin,
     1        ttintmp,sbtot,epstri) !assumed symmetric ttin
         
         call IntTriCos(nin+1,thintmp,temp6,nttin,
     1        ttintmp,sgradpar,epstri) 
         call IntTriCos(nin+1,thintmp,temp7,nttin,
     1        ttintmp,spsidott,epstri) 
         
         call IntTriCos(nin+1,thintmp,temp8,nttin*ndt,
     1        tfine1,tempfine8,epstri) 
         call IntTriCos(nin+1,thintmp,temp5,nttin*ndt,
     1        tfine1,tempfine5,epstri) 
         call IntTriCos(nin+1,thintmp,temp9,nttin*ndt,
     1        tfine1,tempfine9,epstri) 
         do i=1,nttin
            call chftransq(chcoeff5,tempfine5((i-1)*ndt+1), ndt, cftm)
            call chderiv(0.0d0, ndt, chcoeff5,dtempfine5)
            dBdtheta(i)=dtempfine5/dttin(i)
         end do
         do i=1,nttin 
            call chftransq(chcoeff9,tempfine9((i-1)*ndt+1), ndt, cftm)
            call chderiv(0.0d0, ndt, chcoeff9,dtempfine9)
            sdalphat(i)=dtempfine9/dttin(i)
         end do
      end if ! end if isymud.eq.1

c      Jac(1:nttin)=1.0d0/(sbtot(1:nttin)*sgradpar(1:nttin))
c      alphat(1:nttin)=fpolin*Jac(1:nttin)/sR(1:nttin)**2.0d0
      alphapsi0(1)=0.0d0
      do i=1,nttin
         dalphapsi0=0.0d0
         do j=1,ndt
            dalphapsi0=dalphapsi0+spdef(j)*tempfine8((i-1)*ndt+j)
         end do
         dalphapsi0=dalphapsi0*dttin(i)/2.0d0*(2.0d0*pi)
         alphapsi0(i+1)= alphapsi0(i)+ dalphapsi0
      end do


c      write(*,*) 'sbtot',sbtot(1:nttin)
c      write(*,*) 'sgradpar',sgradpar(1:nttin)
c      write(*,*) 'fout2',fout2
      return 

      end

      subroutine calmetric0(ncon,psicon,Rcon,Zcon,tmaxp,icon
     1        ,psiRcon,psiZcon,fpolcon
     2     ,nttin,ttin, btot,alphapsi0,alphatheta0)

      use arrays, only:  nt2, ksamp2

      implicit real*8 (a-h,o-z)
      save

      parameter (ndt=8)

      real *8 psicon(*),psiRcon(*),psiZcon(*)
      real *8 Rcon(*), Zcon(*),tmaxp(*),fpolcon(*)
      real *8 psiRin(nt2*ksamp2+1),psiZin(nt2*ksamp2+1)
      real *8 Rin(nt2*ksamp2+1), Zin(nt2*ksamp2+1),thin(nt2*ksamp2+1)
      real *8 ttin(*), btot(*),alphapsi0(*),alphatheta0(*)
      real *8 sbtot(nttin), salphapsi0(nttin),salphatheta0(nttin)

      real *8 tn(ndt),cftm(ndt*ndt)
      real *8 spbx(ndt*ndt), spxb(ndt*ndt), spdef(ndt)
c      real *8 fint1(ncon),fint2(ncon), fint3(ncon), chcoeff(ncon)
c      real *8 dfint1(ncon), lambda

      integer *4  icon(*)
      integer i,j,k

c      ksamp3=ksamp2/4 
      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi

      !to calculate integral from 0 to theta
      call chsetupq(ndt, tn, cftm, spbx, spxb, spdef) 
      tn(1:ndt)=1.0d0-tn(1:ndt) 
c      write(*,*) 'psicon1',psicon(1:ncon)
c      write(*,*) 'psiRintot',psiRcon(1:ncon*ksamp3)
      do k=1, ncon
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         Rin(1:nin)=Rcon(istart:iend)
         Zin(1:nin)=Zcon(istart:iend)
         thin(1:nin)=tmaxp(istart:iend)
         psiRin(1:nin)=psiRcon(istart:iend)
         psiZin(1:nin)=psiZcon(istart:iend)
         Rin(nin+1)=Rcon(istart)
         Zin(nin+1)=Zcon(istart)
         thin(nin+1)=tmaxp(istart)+2.0d0*pi
         psiRin(nin+1)=psiRcon(istart)
         psiZin(nin+1)=psiZcon(istart)


         call findsinglemetric0(nin,Rin, Zin
     1     ,thin, psiRin, psiZin, fpolcon(k)
     2     ,nttin, ttin, ndt, tn, spdef, cftm
     3     ,sbtot,salphapsi0,salphatheta0)
         
         do j=1,nttin
            btot((j-1)*ncon+k)=sbtot(j)
            alphapsi0((j-1)*ncon+k)=salphapsi0(j)
            alphatheta0((j-1)*ncon+k)=salphatheta0(j)
         end do
      end do
     
      return
      end 

      subroutine findsinglemetric0(nin,Rin,Zin
     1     ,thin,psiRin,psiZin, fpolin
     2     ,nttin, ttin, ndt, tn, spdef, cftm
     3     ,sbtot,salphapsi0,salphatheta0)

      use arrays, only: nt2, ksamp2, epstri
      use arrays, only: isymud, Rmid, lambda, Rmaxis, Zmaxis

      implicit real*8 (a-h,o-z)
c      save
      real *8 Rin(*), Zin(*),thin(*)
    
      real *8 Rint(nin+1),dRdt(nin+1) ,Zint(nin+1)
      real *8 dZdt(nin+1),dldt(nin+1), ttintmp(nttin+1)
 
      real *8 Rintmp(nin+1), Zintmp(nin+1),thintmp(nin+1)
      real *8 psiRin(*),psiZin(*), BpR2(nin+1),BR2(nin+1)
      real *8 temp0(nin+1), temp1(nin+1), temp2(nin+1)
      real *8 temp3(nin+1), temp4(nin+1), temp5(nin+1)
      real *8 temp6(nin+1), temp7(nin+1), temp8(nin+1)
      real *8 gradparin(nin+1), Jac(nin+1), temp88(nin+1)

      real *8 ttin(*), sbtot(*),salphapsi0(*),salphatheta0(*)

      real *8 tn(*), cftm(*), spdef(*),dttin(nttin),tfine(nttin*ndt)
      real *8 tempfine5(nttin*ndt), tempfine8(nttin*ndt), chcoeff5(ndt)


      integer *4 i,j,k, kLag2, iw, itype, nint

      save
       !psi in this routine is normalized one to the range from 0 to 1
      nint=nin !nt2*ksamp2/2
      pi = 4.0d0*datan(1.0d0)

      do i=1,nttin ! assumed ttin(1)=0
         if (i.ne.nttin) then
            dttin(i)=ttin(i+1)-ttin(i)
         else
            dttin(i)=ttin(1)+2.0d0*pi-ttin(nttin)
         end if

         tfine((i-1)*ndt+1:i*ndt)=ttin(i)+tn(1:ndt)*dttin(i)
      end do

      kLag2=4
      if (isymud.eq.1) then
         do i=1,nin+1
            thintmp(i)=thin(i)
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
      else
         do i=1,nin+1
            thintmp(i)=thin(i)/2.0d0
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
         tfine(1:nttin*ndt)=tfine(1:nttin*ndt)/2.0d0
         ttintmp(1:nttin)=ttin(1:nttin)/2.0d0
      end if
 
      BpR2(1:nin+1)=(psiRin(1:nin+1)**2.0d0+psiZin(1:nin+1)**2.0d0) 
     1     /lambda**2 
c      write(*,*) 'BpR2(1:nin+1)',BpR2(1:nin+1)
      BR2(1:nin+1)= BpR2(1:nin+1)+fpolin**2.0d0  
c      BR2(1:nin+1)= BpR2(1:nin+1)+fpolin**2.0d0  
      temp0(1:nin+1)=1.0d0/
     1     (((Zin(1:nin+1)-Zmaxis)/(Rin(1:nin+1)-Rmaxis))**2.0d0+1.0d0)
      temp1(1:nin+1)=-temp0(1:nin+1)*(Zin(1:nin+1)-Zmaxis)
     1     /(Rin(1:nin+1)-Rmaxis)**2.0d0   !dtheta/dR
      temp2(1:nin+1)=temp0(1:nin+1)
     1     /(Rin(1:nin+1)-Rmaxis)   !dtheta/dZ
      temp3(1:nin+1)=(psiRin(1:nin+1)*temp2(1:nin+1)
     1     -psiZin(1:nin+1)*temp1(1:nin+1))/Rin(1:nin+1)/lambda !(grad phi times grad Psi) cdot grad theta)
      temp5(1:nin+1)=dsqrt(BR2(1:nin+1))/Rin(1:nin+1) !magnitude of Btot

c      temp8(1:nin+1)=1.0d0/dabs(temp3(1:nin+1))
c     1       /Rin(1:nin+1)**2 
      temp8(1:nin+1)=1.0d0/temp3(1:nin+1)
     1       /Rin(1:nin+1)**2                    ! 1/(B dot grad theta)/R^2

c      write(*,*) 'Rin3',Rin(1:nin+1)
c      write(*,*) 'Bp3',dsqrt(BpR2(1:nin+1))/Rin(1:nin+1)
c      write(*,*) 'Btot3',temp1(1:nin+1),lambda

      if (isymud.eq.1) then
  
     
         call IntTriCos(nin/2+1,thintmp,temp5,nttin/2+1,
     1        ttin,sbtot,epstri) !assumed symmetric ttin
         call IntTriCos(nin/2+1,thintmp,temp8,nttin/2+1,
     1        ttin,temp88,epstri) 
         call IntTriCos(nin/2+1,thintmp,temp8,nttin*ndt/2+1,
     1        tfine,tempfine8,epstri) 
 

         do i=2,nttin/2
            sbtot(nttin+2-i)=sbtot(i)
            temp88(nttin+2-i)=temp88(i)
         end do

         do i=2,nttin*ndt/2+1
            tempfine8(nttin*ndt+2-i)=tempfine8(i)
         end do

      else  ! if isymud is not 1 (up-down assymmetric)

         call IntTriCos(nin+1,thintmp,temp5,nttin,
     1        ttintmp,sbtot,epstri) !assumed symmetric ttin
         call IntTriCos(nin+1,thintmp,temp8,nttin,
     1        ttintmp,temp88,epstri) 

         call IntTriCos(nin+1,thintmp,temp8,nttin*ndt,
     1        tfine,tempfine8,epstri) 

      end if ! end if isymud.eq.1

      salphapsi0(1)=0.0d0
      do i=1,nttin
         dalphapsi0=0.0d0
         do j=1,ndt
            dalphapsi0=dalphapsi0+spdef(j)*tempfine8((i-1)*ndt+j)
         end do
         dalphapsi0=dalphapsi0*dttin(i)/2.0d0*(2.0d0*pi)
         salphapsi0(i+1)= salphapsi0(i)+ dalphapsi0
      end do
      salphatheta0(1:nttin)=temp88(1:nttin)*fpolin
      
      return 

      end

      subroutine findfpass(ncon,lambda,psicon,Rcon,Zcon,tmaxp,icon
     1    ,psiRcon,psiZcon
     2     ,fpolch,nchy,psiichy,spydef,fpass)

      use arrays, only:  nt2, ksamp2

      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*),psiRcon(*),psiZcon(*),fpolch(*)
      real *8 Rcon(*), Zcon(*),tmaxp(*)
      real *8 psiRin(nt2*ksamp2+1),psiZin(nt2*ksamp2+1)
      real *8 Rin(nt2*ksamp2+1), Zin(nt2*ksamp2+1),thin(nt2*ksamp2+1)
      real *8 fpass(*),psiichy(*),spydef(*),lambda


      integer *4  icon(*)
      integer i,j,k

c      ksamp3=ksamp2/4 
      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi
c      write(*,*) 'psicon1',psicon(1:ncon)
c      write(*,*) 'psiRintot',psiRcon(1:ncon*ksamp3)

c     call chsetupq(nchy, psiichy, cftmy, spybx, spyxb, spydef)

      do k=1, ncon
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         Rin(1:nin)=Rcon(istart:iend)
         Zin(1:nin)=Zcon(istart:iend)
         thin(1:nin)=tmaxp(istart:iend)
         psiRin(1:nin)=psiRcon(istart:iend)
         psiZin(1:nin)=psiZcon(istart:iend)
         Rin(nin+1)=Rcon(istart)
         Zin(nin+1)=Zcon(istart)
         thin(nin+1)=tmaxp(istart)+2.0d0*pi
         psiRin(nin+1)=psiRcon(istart)
         psiZin(nin+1)=psiZcon(istart)

         call findsinglefpass(nin,lambda,Rin, Zin
     1     ,thin, psiRin, psiZin
     2     ,fpolch(k),nchy,psiichy,spydef,fpass(k))
      end do
c      write(*,*) 'psicon2',psicon(1:ncon)
      return
      end 

      subroutine findsinglefpass(nin,lambda,Rin,Zin,thin,psiRin,psiZin
     1    , fpol,nchy,psiichy,spydef, fpassout)

      use arrays, only:  nt2, ksamp2, epstri
      use arrays, only: isymud, Rmaxis, Zmaxis
 
      implicit real*8 (a-h,o-z)
c      save
c      parameter (nint=nt2)
      real *8 Rin(*), Zin(*),thin(*),psiichy(*),spydef(*)
    
      real *8 Rint(nin+1),Zint(nin+1), dldt(nin+1)
 
      real *8 Rintmp(nin+1), Zintmp(nin+1),thintmp(nin+1)
      real *8 psiRin(*),psiZin(*), BpR2(nin+1), tttmp(nin)
      real *8 B2(nin+1),B(nin+1)

      real *8 temp0(nin+1), temp1(nin+1), temp2(nin+1)
      real *8 temp3(nin+1), temp4(nin+1), temp5(nin+1)
      real *8 temp6(nin+1), temp7(nin+1), temp8(nin+1)
  
      real *8 thint(nin+1),integ11(nin+1)
      real *8 integ12(nin+1),integ13(nin+1)

      real *8 y(nchy), yden(nchy),lambda
      integer *4 i,j,k, kLag2, iw, itype, nint

      save
      
      nint=nin !nt2*ksamp2/2
      pi = 4.0d0*datan(1.0d0)
c      epstri= 1.0d-10

      delth = 0.0d0
      do i=1,nint
         thint(i)=2.0d0*pi*(i-1)/nint+thin(1) !delth
      end do


      kLag2=4
      if (isymud.eq.1) then
         do i=1,nin+1
            thintmp(i)=thin(i)
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
      else
         do i=1,nin+1
            thintmp(i)=thin(i)/2.0d0
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
         tttmp(1:nint)=thint(1:nint)/2.0d0
      end if
 

      BpR2(1:nin+1)=(psiRin(1:nin+1)**2.0d0+psiZin(1:nin+1)**2.0d0)

      B2(1:nin+1)=(BpR2(1:nin+1)/lambda**2+fpol**2)/Rin(1:nin+1)**2
c      write(*,*) 'B2',fpol**2, B2(1:nin+1), BpR2(1:nin+1)/lambda**2
      B(1:nin+1)=dsqrt(B2(1:nin+1))
  
      Bmax=0.0d0
      do i=1,nin+1
         if (B(i).gt.Bmax) then
            Bmax=B(i)
         end if
      end do

      temp0(1:nin+1)=1.0d0/
     1     (((Zin(1:nin+1)-Zmaxis)/(Rin(1:nin+1)-Rmaxis))**2.0d0+1.0d0)
      temp1(1:nin+1)=-temp0(1:nin+1)*(Zin(1:nin+1)-Zmaxis)
     1     /(Rin(1:nin+1)-Rmaxis)**2.0d0   !dtheta/dR
      temp2(1:nin+1)=temp0(1:nin+1)
     1     /(Rin(1:nin+1)-Rmaxis)   !dtheta/dZ

      temp3(1:nin+1)=dabs((psiRin(1:nin+1)*temp2(1:nin+1)
     1     -psiZin(1:nin+1)*temp1(1:nin+1))/Rin(1:nin+1)) !(grad phi times grad psi) cdot grad theta)

      temp4(1:nin+1)=1.0d0/temp3(1:nin+1)*B2(1:nin+1)
      temp5(1:nin+1)=1.0d0/temp3(1:nin+1)

      if (isymud.eq.1) then

         call IntTriCos(nin/2+1,thintmp,temp4,nint/2+1,
     1        thint,integ11,epstri)
  
         call IntTriCos(nin/2+1,thintmp,temp5,nint/2+1,
     1        thint,integ12,epstri)

         do i=1,nint/2
            integ11(nint+2-i)=integ11(i)
            integ12(nint+2-i)=integ12(i)
         end do
      else                      ! if isymud is not 1 (up-down assymmetric)
        
         call IntTriCos(nin+1,thintmp,temp4,nint+1,
     1        tttmp,integ11,epstri)
  
         call IntTriCos(nin+1,thintmp,temp5,nint+1,
     1        tttmp,integ12,epstri)

      end if                    ! end if isymud.eq.1

      fout1=0.0d0
      fout2=0.0d0
      do i=1,nint
         fout1=fout1+integ11(i)
         fout2=fout2+integ12(i)  
      end do
      B2avg=fout1/fout2
c      write(*,*) 'integ12',integ12
c      write(*,*) 'B2avg',fout2,fout1,B2avg
      y(1:nchy)=psiichy(1:nchy)/Bmax  
c      write(*,*) 'y',y(1:nchy)

      do j=1,nchy   !find denominator of the integrand
         temp6(1:nin+1)=1.0d0/temp3(1:nin+1)
     1        *dsqrt(1.0d0-y(j)*B(1:nin+1))

         if (isymud.eq.1) then

            call IntTriCos(nin/2+1,thintmp,temp6,nint/2+1,
     1        thint,integ13,epstri)
            do i=1,nint/2
               integ13(nint+2-i)=integ13(i)
            end do
         else                   ! if isymud is not 1 (up-down assymmetric)
        
            call IntTriCos(nin+1,thintmp,temp6,nint+1,
     1        tttmp,integ13,epstri)
         end if                 ! end if isymud.eq.1
         fout3=0.0d0
         do i=1,nint
            fout3=fout3+integ13(i)
         end do
         yden(j)=fout3/fout2
      end do

c      write(*,*) 'integ03',integ03
c      write(*,*) 'yden', yden
      fpint=0
      do j=1,nchy
         fpint=fpint+spydef(j)*y(j)/yden(j)
      end do
      fpint=fpint/(2.0d0*Bmax)
      fpassout=0.75d0*B2avg*fpint

c      write(*,*) 'fpass',fpassout,fpint, B2avg

      return 

      end

      subroutine findLagrange(ncon,lambda,psiichq,Ia,Ib,Ic
     1     ,Itf,pprimch,fpolcon,cftmq,spdef,spbx,spxb
     2     ,WLag)

      use arrays, only: iiter,fout, itftype
      implicit real*8 (a-h,o-z)
      save

      real *8 psiichq(*),Ia(*),Ib(*),Ic(*),Itf(*),fpolcon(*)
      real *8 cftmq(*), spdef(*), spbx(*), spxb(*),pprimch(*)
      real *8 presch(ncon), pi, mu, lambda
 
      integer i,j,k

      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi

!     <Integration using Chebyshev grid>

    
      presch(1:ncon)=0.0d0
      do i=1,ncon     
         do k=1,ncon           
            presch(i)=presch(i)+spxb((k-1)*ncon+i)*pprimch(k)
         end do
      end do
      if ((iscale.eq.5) .or.(iscale.eq.6)) then !pprimch=dp/dpsi (psi is normalized flux)
         presch(1:ncon)=presch(1:ncon)/(-2.0d0)
      else                      !pprimch=dp/dPsi=lambda*dp/dpsi (Psi is actual scale flux)
         presch(1:ncon)=presch(1:ncon)/(lambda*2.0d0)
      end if

c      write(*,*) 'presch',presch(1:ncon)
c      write(*,*) 'fpolcon',fpolcon(1:ncon)
c      write(*,*) 'Ia',Ia(1:ncon)
c      write(*,*) 'Ib',Ib(1:ncon)
c      write(*,*) 'Ic',Ic(1:ncon)
c      write(*,*) 'Itf',Itf(1:ncon)
      Wpol=0.0d0
      Wtor=0.0d0
      Wpre=0.0d0
      do i=1,ncon     
         Wpol=Wpol+spdef(i)*Ia(i)
         if (itftype.eq.0) then !toroidal flow is not considered
            Wpre=Wpre+spdef(i)*Ib(i)*presch(i)
         else !toroidal flow is considered
            Wpre=Wpre+spdef(i)*Itf(i)*presch(i)
         end if
         Wtor=Wtor+spdef(i)*Ic(i)*fpolcon(i)**2
 
      end do
      Wpol=Wpol/(lambda**2*4.0d0)*2*pi
      Wpre=Wpre/(2.0d0)*mu*2*pi
      Wtor=Wtor/(4.0d0)*2*pi
      
      WLag=Wpol-Wtor-Wpre
      Wtot=Wpol+Wtor+Wpre
c      write(*,*) 'Wpol,Wpre,Wtor',Wpol,Wpre,Wtor
 
c      write(fout,*) 'iiter, WLag, Wtot', iiter, WLag, Wtot
      return
      end 

      subroutine findPhyVal(ncon,lambda,psiichq,Ia,Ib,Ic,Id,Lp
     1    ,pprimch,fpolcon,cftmq,spdef,spbx,spxb
     2     ,plmvol,plmArea,acmcur,presch,betaP,intind,B0
     3     ,presavg,Wpol,betaTot,plmvol0,plmArea0,acmcur0,betaP0
     4     ,intind0,edgecur,parcur)

      use arrays, only: Rmid,iptype, Rmaxis, Zmaxis, fout
      use arrays, only: psi0, psiB, iscale, ijtype,fpoledge
      use arrays, only: jpar0, jpin, jpout, Jomax, Joloc, Jow
      implicit real*8 (a-h,o-z)
      save

      real *8 psiichq(*),Ia(*),Ib(*),Ic(*),Id(*),Lp(*),fpolcon(*)
      real *8 cftmq(*), spdef(*), spbx(*), spxb(*),pprimch(*)
      real *8 plmvol(*),plmArea(*),acmcur(*),presch(*),betaP(*)
      real *8 intind(*),parcur(*)

      real *8 chcoeff1(ncon),chcoeff2(ncon)
      real *8 chcoeff3(ncon),chcoeff4(ncon),chcoeff5(ncon)
      real *8 chcoefftmp(ncon),temp(ncon)
      real *8 lambda,Lpedge,mu,intind0
  
      integer i,j,k

      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi

!     Accumulated toroidal current inside the flux surface       
      acmcur(1:ncon)=Ia(1:ncon)/(dabs(lambda)*mu)
   

      call chftransq(chcoeff1, Ia, ncon, cftmq)
      call chftransq(chcoeff2, Id, ncon, cftmq)
c      call chftransq(chcoeff3, Ic, ncon, cftmq)
      call chftransq(chcoeff4, Lp, ncon, cftmq)
      call chftransq(chcoeff5, fpolcon, ncon, cftmq)
c      write(*,*) 'fpolcon', fpolcon(1:ncon)

!     Total toroidal current in [A]  
      call chfit(1.0d0,ncon,chcoeff1,acmcur0)
      acmcur0=acmcur0/(dabs(lambda)*mu)


!     Parallel current density at many flux surfaces   

      do i=1,ncon
         psin=1.0d0-psiichq(i)
         call chderiv(psin, ncon, chcoeff1,dIa1)
         call chderiv(psin, ncon, chcoeff5,dF1)
         parcur(i)=(dIa1-dF1/fpolcon(i)*Ia(i))/(Ic(i)*lambda) 
c         write(*,*) 'parcur1',dIa1,dF1/fpolcon(i)*Ia(i)
         ! parallel current profile includes mu0
      end do

      if (ijtype.eq.4) then
         acmparcur1=0.0d0
         acmparcur2=0.0d0
         do i=1,ncon
            psin=1.0d0-psiichq(i)
            acmparcur1=acmparcur1+spdef(i)*Id(i)*
     1           jpar0*(1.0d0-psin**jpin)**jpout
            acmparcur2=acmparcur2+spdef(i)*Id(i)*
     1           jpar0*Jomax*exp(-((psin-Joloc)/Jow)**2)
         end do
         acmparcur1=acmparcur1/2.0d0/mu
         acmparcur2=acmparcur2/2.0d0/mu
      end if

!     Toroidal surface current at the last flux surface  
      call chderiv(1.0d0,ncon, chcoeff1,edgecur)
      edgecur=edgecur*lambda
      edgecur=edgecur/(lambda*mu)

!     Perimeter of the last flux surface [m]  
      call chfit(1.0d0,ncon,chcoeff4,Lpedge)
     

!     <Integration using Chebyshev grid>
!     Plasma volume inside a flux surface [m^3]
!     : plmvol=(2*pi)*(int_{-1}^{x} Ib dx)/2.0
!
!     Total plasma volume [m^3]:
!     : plmvol0=(2*pi)*(int_{-1}^{1} Ib dx)/2.0
!
!     Poloidal cross section inside a flux surface [m^2]
!     : plmArea=(int_{-1}^{x} Id dx)/2.0
!
!     Total plasma volume [m^3]:
!     : plmArea0=(int_{-1}^{1} Id dx)/2.0
!
!     Pressure [Pa]
!     : p(x)=int_{x}^{1} dp/dpsi dx/(lambda*2.0)
!
!     internal inductance []
!     : intind(x)=mu*int_{-1}^{1} Ia((x+1)/2) dx*(lambda/2.0)
!               /(Ia(x=1)**2*Rmid)
   
      plmvol(1:ncon)=0.0d0
      plmvol0=0.0d0
      plmArea(1:ncon)=0.0d0
      plmArea0=0.0d0
      presch(1:ncon)=0.0d0
     
      intind(1:ncon)=0.0d0
      intind0=0.0d0
      Wpol=0.0d0
      do i=1,ncon     
         do k=1,ncon           
            plmvol(i)=plmvol(i)+spbx((k-1)*ncon+i)*Ib(k)        
            plmArea(i)=plmArea(i)+spbx((k-1)*ncon+i)*Id(k) 
            presch(i)=presch(i)+spxb((k-1)*ncon+i)*pprimch(k)        
            intind(i)=intind(i)+spbx((k-1)*ncon+i)*Ia(k)
         end do
         plmvol0=plmvol0+spdef(i)*Ib(i)
         plmArea0=plmArea0+spdef(i)*Id(i)
         intind0=intind0+spdef(i)*Ia(i)
         Wpol=Wpol+spdef(i)*Ia(i)
      end do

      plmvol(1:ncon)=plmvol(1:ncon)/2.0d0*(2*pi)
      plmvol0=plmvol0/2.0d0*(2*pi)

      plmArea(1:ncon)=plmArea(1:ncon)/2.0d0
      plmArea0=plmArea0/2.0d0

      if ((iscale.eq.5) .or.(iscale.eq.6)) then !pprimch=dp/dpsi (psi is normalized flux)
         presch(1:ncon)=presch(1:ncon)/(-2.0d0)
      else  !pprimch=dp/dPsi=lambda*dp/dpsi (Psi is actual scale flux)
         presch(1:ncon)=presch(1:ncon)/(lambda*2.0d0)
      end if
      intind(1:ncon)=intind(1:ncon)*(2*pi)
     1     /(mu**2*acmcur0**2*Rmid*(lambda**2))
      intind0=intind0*(2*pi)
     1     /(mu**2*acmcur0**2*Rmid*(lambda**2))
      Wpol=Wpol/(lambda**2*4.0d0) !*(2*pi)


!     Volume averaged Pressure [Pa]
!     : presavg=int_{-1}^{1} (p*Ib) dx/(lambda*2.0)
!              =- int {-1)^{1} (dp/dpsi)*vol dx
!     internal inductance []
!     : intind(x)=mu*int_{-1}^{1} Ia((x+1)/2) dx*(lambda/2.0)
!               /(Ia(x=1)**2*Rmid)
 

      betaP(1:ncon)=0.0d0
      Presavg=0.0d0
      Presavg2=0.0d0
      Pres2avg=0.0d0
      do i=1,ncon     
         do k=1,ncon
            ! <p>=(int_{-1}^{x} p((x+1)/2)*Ib((x+1)/2) dx)*(lambda/2.0)/Vol
            betaP(i)=betaP(i)+spbx((k-1)*ncon+i)*Ib(k)*presch(k)
         end do
         Presavg=Presavg+spdef(i)*pprimch(i)*plmvol(i)
         Presavg2=Presavg2+spdef(i)*Ib(i)*presch(i)
         Pres2avg=Pres2avg+spdef(i)*Ib(i)*presch(i)**2
      end do

      if ((iscale.eq.5) .or.(iscale.eq.6)) then !pprimch=dp/dpsi (psi is normalized flux)
         Presavg=(-1.0d0)*Presavg/(2.0d0)/plmvol0
      else
         Presavg=(-1.0d0)*Presavg/dabs(lambda*2.0d0)/plmvol0
      end if

      Presavg2=Presavg2/2.0d0*(2*pi)/plmvol0
      Pres2avg=Pres2avg/2.0d0*(2*pi)/plmvol0
         
c      betaP(1:ncon)=2.0d0*Lpedge**2/(acmcur0**2*mu)
c     1     *betaP(1:ncon)*dabs(lambda)/2.0d0/plmvol(1:ncon)
c      betaP0 =2.0d0*Lpedge**2/(acmcur0**2*mu)*Pressavg
      betaP(1:ncon)=2.0d0/(mu*acmcur0**2*Rmid)
     1     *betaP(1:ncon)/2.0d0*(2*pi)*2.0d0
      betaP0 =2.0d0/(mu*acmcur0**2*Rmid)*Presavg*plmvol0*2.0d0


      call chfit(1.0d0,ncon,chcoeff5,fpoledge)
      call chfit(0.0d0,ncon,chcoeff5,fpol0)
c      write(*,*) 'chcoeff5',chcoeff5(1:ncon)
      write(*,*) 'fpoledge',fpoledge,fpol0

      B0=dabs(fpoledge/Rmid) !*2*pi
      write(*,*) 'Rmid',Rmid,B0
c      B0=dabs(fpol0/Rmid)
      betaTot= 2.0d0*mu*Presavg/B0**2
      psi0 = (1.0d0/lambda)+psiB
      write(fout,*) '*************************************'
      write(fout,*) '*****< Output: Physical Values >*****'
      write(fout,*) 'Magnetic axis :(R,Z)=', Rmaxis, Zmaxis
      write(fout,*) 'Psi at the magnetic axis (Psi0) : '
     1     , psi0
      write(fout,*) 'Psi at the boundary (PsiB): ', psiB
      write(fout,*) 'Normalized Psi (Psi-PsiB)/(Psi0-PsiB) :  '
     1     ,psiichq(1:ncon)
      write(fout,*) 'Total toroidal current [A] : ',acmcur0

      if (ijtype.eq.4) then
         write(fout,*) 'Total parallel current1 [A](%) : ',acmparcur1
     1        ,'( ',acmparcur1/(acmparcur1+acmparcur2)*100,' %)'
         write(fout,*) 'Total parallel current2 [A] : ',acmparcur2
     1        ,'( ',acmparcur2/(acmparcur1+acmparcur2)*100,' %)'
      end if

      write(fout,*) 'Tor. current inside psi :'
      write(fout,*) acmcur(1:ncon)
      write(fout,*) 'Par. current density at psi [A/m^2*mu0]:'
      write(fout,*) parcur(1:ncon)
      write(fout,*) 'Edge surface. current density dI/dpsi : ',lpedge
      write(fout,*) 'Total plasma volume [m^3] : ',plmvol0
      write(fout,*) 'Voulme inside psi :'
      write(fout,*) plmvol(1:ncon)

      write(fout,*) 'Total poloidal cross section [m^2] : ',plmArea0
      write(fout,*) 'cross section inside psi :'
      write(fout,*) plmArea(1:ncon)

      write(fout,*) 'Volume averaged Pressure[Pa]',Presavg
      write(fout,*) 'Volume averaged Pressure[Pa](integralB)',Presavg2
      write(fout,*) 'Pressure at psi [Pa]:'
      write(fout,*) presch(1:ncon)
      write(fout,*) 'Volume averaged Pressure**2 [Pa**2]',Pres2avg
      write(fout,*) 'Total internal inductance : ',intind0
      write(fout,*) 'Internal inductance inside psi :'
      write(fout,*)  intind(1:ncon)
      write(fout,*) 'Poloidal field energy : ', Wpol

      write(fout,*) 'Vacuum magnetic field  at R=Rmid (B0): ',B0
      write(fout,*) 'Beta (2.0d0*mu*Presavg/B0**2)',betaTot
     
      write(fout,*) 'Poloidal beta',betaP0
      write(fout,*) 'Poloidal beta inside psi :'
      write(fout,*) betaP(1:ncon)
      return
      end

      subroutine Troyon(torcur,B0,reps,R0,rkappa,plmArea0
     1    ,betaTroy1,betaTroy2)

      use arrays, only: fout
      implicit real*8 (a-h,o-z)
      save

      real *8 mu,pi

      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi
      qstar = 2*B0*plmArea0/(mu*R0*torcur)

      torcur=torcur*1.0d-6 ! [MA]
      betaTroy1=0.028d0*torcur/(B0*R0*reps)
      betaTroy2=0.14d0*(reps*rkappa/qstar)

      write(fout,*) '*************************************'
      write(fout,*) '*****< Troyon beta limit >*****'
      write(fout,*) 'beta_t=0.028*I0/(a*B0)=',betaTroy1
      write(fout,*) 'beta_t=0.14*eps*kappa/q=',betaTroy2
      return
      end

      subroutine Mercier(ncon,lambda,psiichq,Ia,Ib,Ic,Ie,Ig,Ih
     1    ,pprimch,fpolcon,qpsich,cftmq,dqdpsi,nDi)

      use arrays, only: fout
      implicit real*8 (a-h,o-z)
      save

      real *8 psiichq(*),Ia(*),Ib(*),Ic(*),Ie(*),Ig(*),Ih(*)
      real *8 pprimch(*),fpolcon(*),qpsich(*),cftmq(*)

      real *8 chcoeff1(ncon),chcoeff2(ncon)
      real *8 J1(ncon),J2(ncon),J3(ncon),J4(ncon),J5(ncon),J6(ncon)
      real *8 dJ5(ncon),dqdpsi(*),nDi(*)
      real *8 lambda,mu,pi
  
      integer i,j,k

      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi

!     scaling     
      J1(1:ncon)=Ie(1:ncon)*(dabs(lambda)**3)/(2.d0*pi)
      J2(1:ncon)=Ig(1:ncon)*(dabs(lambda)**3)/(2.d0*pi)
      J3(1:ncon)=Ih(1:ncon)*(dabs(lambda)**3)/(2.d0*pi)
      J4(1:ncon)=Ic(1:ncon)*(dabs(lambda))/(2.d0*pi)
      J5(1:ncon)=Ib(1:ncon)*(dabs(lambda))/(2.d0*pi)
      J6(1:ncon)=Ia(1:ncon)/(dabs(lambda))/(2.d0*pi)

c      write(*,*) 'J1',J1(1:ncon)
c      write(*,*) 'J2',J2(1:ncon)
c      write(*,*) 'J3',J3(1:ncon)
c      write(*,*) 'J4',J4(1:ncon)
c      write(*,*) 'J5',J5(1:ncon)
c      write(*,*) 'J6',J6(1:ncon)

      call chftransq(chcoeff1, J5, ncon, cftmq)
      call chftransq(chcoeff2, qpsich, ncon, cftmq)

      do i=1,ncon
         psin=1.0d0-psiichq(i)
         call chderiv(psin, ncon, chcoeff1,dJ5(i))
         call chderiv(psin, ncon, chcoeff2,dqdpsi(i))
      end do
      dJ5(1:ncon)=dJ5(1:ncon)*lambda
      dqdpsi(1:ncon)=dqdpsi(1:ncon)*lambda

      do i=1,ncon
         D1=(mu*pprimch(i)*fpolcon(i)*J2(i)/dqdpsi(i)-0.5d0)**2
         D2=mu*pprimch(i)/dqdpsi(i)**2*(dJ5(i)-mu*pprimch(i)*J3(i))
     1        *(fpolcon(i)**2*J1(i)+J4(i))
         nDi(i)=D1+D2
c         write(*,*) 'D1,D2',i,D1,D2,dqdpsi(i),mu*pprimch(i),dJ5(i)
c         write(*,*) 'test',mu*pprimch(i)*fpolcon(i)*J2(i)/dqdpsi(i)
      end do

      write(fout,*) '*************************************'
      write(fout,*) '*****< Mercier stability condition >*****'
      write(fout,*) 'necessary condition for stability: -Di>0'

      write(fout,*) 'Normalized Psi (Psi-PsiB)/(Psi0-PsiB) :  '
     1     ,psiichq(1:ncon)
      write(fout,*) '-Di : ',nDi(1:ncon)

   
      return
      end

      subroutine matchtorcur(ncon,lambda,Ia,cftmq,torcur,factorl)

      implicit real*8 (a-h,o-z)
      save

      real *8 Ia(*), cftmq(*)
      real *8 chcoeff1(ncon)
      real *8 lambda,mu, torcur0, factorl
  
      integer i,j,k

      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi

      call chftransq(chcoeff1, Ia, ncon, cftmq)

c      write(*,*) 'chcoe',chcoeff1(1:ncon), Ia(1:ncon)
      call chfit(1.0d0,ncon,chcoeff1,torcur0)

      torcur0=dabs(torcur0/(lambda*mu))
    
      factorl=torcur0/torcur
      write(*,*) 'factorl',factorl,torcur0,torcur
      return
      end

      subroutine findIpar(ncon,lambda,psicon,Ia,Ib,Ic
     1    ,fpol,ffprimch,pprimch,cftmq,spdef,curpar)

      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*),Ia(*),Ib(*),Ic(*)
      real *8 fpol(*), curpar(*)
      real *8 ffprimch(*),pprimch(*), cftmq(*), spdef(*)
      real *8 chcoeff1(ncon),chcoeff2(ncon)
      real *8 chcoeff3(ncon),chcoeff4(ncon)
      real *8 chcoefftmp(ncon),temp(ncon)
      real *8 mu
  
      integer i,j,k

      pi=4.0d0*datan(1.0d0)
      mu=4.0d-7*pi

      do i=1, ncon
         curpar(i)=-Ib(i)/Ic(i)*pprimch(i)*mu
     1        -ffprimch(i)*(1+Ia(i)/Ic(i)/fpol(i)**2)
      end do
c      curpar=curpar/(lambda*mu)

c      write(*,*) 'curpar', curpar(1:ncon)
c      write(*,*) 'ffrimch',ffprimch(1:ncon)

      return
      end

      subroutine findfpol0(ncon,lambda,psicon,Ia,Ib,Ic
     1    ,pprimch,cftmq,ffprimpre)

      
      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*),Ia(*),Ib(*),Ic(*),ffprimpre(*)
      real *8 pprimch(*), cftmq(*)
      real *8 chcoeff(ncon),ffprimch(ncon),temp(ncon)
      real *8 dIa(ncon), lambda, Ibedge, Icedge

      integer i,j,k

      pi=4.0d0*datan(1.0d0)

c      write(*,*) 'cftmq',cftmq(1:ncon*ncon)
c      write(*,*) 'Ia',Ia(1:ncon)
      call chftransq(chcoeff, Ia, ncon, cftmq)
c      write(*,*) 'chcoeff',chcoeff(1:ncon)
      do i=1,ncon
         call chderiv(1.0d0-psicon(i), ncon, chcoeff, dIa(i)) 
         dIa(i)=dIa(i)*lambda
      end do
c      write(*,*) 'dIa',dIa(1:ncon)/lambda**2
c      ffprimch(1:ncon)=(-1.0d1)*(psicon(1:ncon))
c      write(*,*) 'dIa2',-Ic(1:ncon)*ffprimch(1:ncon)
c     1     -(4*pi*1.0d-7)*pprimch(1:ncon)*Ib(1:ncon)
c      write(*,*) 'pprimch',(4*pi*1.0d-7)*pprimch(1:ncon)
c      call chderiv(1.0d0, ncon, chcoeff, dIaedge)
c      write(*,*) 'dIaedge',dIaedge/lambda
c      call chftransq(chcoeff, pprimch, ncon, cftmq)
c      call chfit(1.0d0, ncon, chcoeff, dpedge)
c      write(*,*) 'dpprimedge',(4*pi*1.0d-7)*dpedge
c      call chftransq(chcoeff, ffprimch, ncon, cftmq)
c      call chfit(1.0d0, ncon,  chcoeff, dfedge)
c      write(*,*) 'dpprimedge',dfedge
c      temp(1:ncon)=Ib(1:ncon)
c      call chftransq(chcoeff, temp, ncon, cftmq)
c      call chfit(1.0d0, ncon, chcoeff, Ibedge)
c      write(*,*) 'dpprimedge',(4*pi*1.0d-7)*dpedge*Ibedge
c      temp(1:ncon)=Ic(1:ncon)
c      call chftransq(chcoeff, temp, ncon, cftmq)
c      call chfit(1.0d0, ncon,  chcoeff, Icedge)
c      write(*,*) 'dpprimedge',dfedge*Icedge
c      write(*,*) 'sumedge',(-(4*pi*1.0d-7)*dpedge*Ibedge-dfedge*Icedge)

c      temp(1:ncon)=-Ic(1:ncon)*ffprimch(1:ncon)
c     1     -(4*pi*1.0d-7)*pprimch(1:ncon)*Ib(1:ncon)
c      call chftransq(chcoeff, temp, ncon, cftmq)
c      call chfit(1.0d0, ncon,  chcoeff, dfedge3)
c      write(*,*) 'dfedge3',dfedge3

    
      ffprimpre(1:ncon)=-(dIa(1:ncon)/dabs(lambda)**2+       
     1     (4*pi*1.0d-7)*pprimch(1:ncon)*Ib(1:ncon))/Ic(1:ncon)
 
      return
      end

      subroutine findshat(ncon,lambda,psicon
     1    ,qpsi,cftmq,shat)
      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*), shat(*)
      real *8 qpsi(*), cftmq(*)
      real *8 chcoeff(ncon)
      real *8 dIa(ncon), lambda

      integer i,j,k

      pi=4.0d0*datan(1.0d0)

  
      call chftransq(chcoeff, qpsi, ncon, cftmq)
      do i=1,ncon
         call chderiv(1.0d0-psicon(i), ncon, chcoeff, shat(i)) 
         shat(i)=shat(i)*lambda
      end do
      write(*,*) 'shat',shat(1:ncon)
  
      return
      end

      subroutine findqprofile(nin,lambda,psiin
     1    ,qin,cftmq,nout,psiout,qout,sout)
      implicit real*8 (a-h,o-z)
      save

      real *8 psiin(*), psiout(*)
      real *8 qin(*), qout(*), sout(*), cftmq(*)
      real *8 chcoeff(nin)
      real *8 lambda

      integer i, nin, nout
  
      call chftransq(chcoeff, qin, nin, cftmq)
      do i=1,nout
         call chfit(1.0d0-psiout(i), nin, chcoeff, qout(i))
         call chderiv(1.0d0-psiout(i), nin, chcoeff, sout(i)) 
         sout(i)=sout(i)*lambda
      end do
c      write(*,*) 'psiout',psiout(1:nout)
c      write(*,*) 'qout',qout(1:nout)
c      write(*,*) 'sout',sout(1:nout)
  
      return
      end

      subroutine findfpolalpha(ncon,lambda,fpolcon,Ia,Ic
     1    ,alphaout)

      implicit real*8 (a-h,o-z)
      save

      real *8 fpolcon(*),Ia(*),Ic(*),alphaout(*)
      real *8 lambda

  
      alphaout(1:ncon)=Ia(1:ncon)
     1     /(fpolcon(1:ncon)**2*lambda*Ic(1:ncon))
 
      return
      end
    


      subroutine findfpol1(ncon,lambda,psicon,Ia,Ib,Ic
     1    ,pprimch, fpolcon,qpsich,cftmq,ffprimpre1
     1    ,ffprimpre2)
  
      use arrays, only: fpolsgn, lambdaf, lambdaex, iiter

      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*), fpolcon(*),qpsich(*)
      real *8 ffprimpre1(*),ffprimpre2(*)
      real *8 pprimch(*), cftmq(*),mu0
      real *8 Ia(*), Ib(*),Ic(*), chcoeff2(ncon), IaoF(ncon)
      real *8 dIaoF(ncon), den(ncon),  num1(ncon),lambda, num2(ncon)
      real *8 alpha(ncon),IaIcbyq(ncon), dIaIcbyq(ncon)
      real *8 ffprimpre1tmp(ncon),ffprimpre2tmp(ncon)

      integer i,j,k

      pi=4.0d0*datan(1.0d0)
      mu0=pi*4.0d-7
c      alpha(1:ncon)=Ia(1:ncon)/
c     1     (qpsich(1:ncon)*2.0d0*pi*fpolcon(1:ncon)*lambda)
      alpha(1:ncon)=Ia(1:ncon)*Ic(1:ncon)/(qpsich(1:ncon)*2.0d0*pi)**2
      !write(*,*) 'alpha',alpha(1:ncon)
      den(1:ncon)=-(1.0d0+alpha(1:ncon))
      !IaIcbyq(1:ncon)=Ia(1:ncon)*(2.0d0*pi)/(lambda*fpolcon(1:ncon))
      IaIcbyq(1:ncon)=Ia(1:ncon)*Ic(1:ncon)/qpsich(1:ncon)
      call chftransq(chcoeff2, IaIcbyq, ncon, cftmq)
      do i=1,ncon
         call chderiv(1.0d0-psicon(i), ncon, chcoeff2, dIaIcbyq(i))         
      end do
      !write(*,*) 'dIaIcbyq',dIaIcbyq(1:ncon)
c      num1(1:ncon)=fpolcon(1:ncon)/(2.0d0*pi*Ic(1:ncon))
c     1     *dIaIcbyq(1:ncon)
c      num1(1:ncon)=fpolcon(1:ncon)**2*lambda/(2.0d0*pi)**2
c     1    /qpsich(1:ncon)*dIaIcbyq(1:ncon)
      num1(1:ncon)=qpsich(1:ncon)/(lambda*Ic(1:ncon)**2)
     1     *dIaIcbyq(1:ncon)
c      num2(1:ncon)=mu0*fpolcon(1:ncon)*lambda/(qpsich(1:ncon)*2.0d0*pi)
c     1     *pprimch(1:ncon)*Ib(1:ncon)
      num2(1:ncon)=mu0/Ic(1:ncon)*pprimch(1:ncon)*Ib(1:ncon)
      !write(*,*) 'num1',num1(1:ncon)
      !write(*,*) 'lambdaic',lambda*Ic(1:ncon)**2
      !write(*,*) 'num2',num2(1:ncon)
      ffprimpre1(1:ncon)=num1(1:ncon)/den(1:ncon)
      ffprimpre2(1:ncon)=num2(1:ncon)/den(1:ncon)

      IaoF(1:ncon)=Ia(1:ncon)/fpolcon(1:ncon)
      call chftransq(chcoeff2, IaoF, ncon, cftmq)
      do i=1,ncon
         call chderiv(1.0d0-psicon(i), ncon, chcoeff2, dIaoF(i))         
      end do
      dIaoF(1:ncon)=dIaoF(1:ncon)*lambda

      den(1:ncon)=2*pi*qpsich(1:ncon)*fpolcon(1:ncon)
     1     +Ia(1:ncon)/dabs(lambda)*fpolsgn
c
      num1(1:ncon)= dIaoF(1:ncon)*fpolcon(1:ncon)/dabs(lambda)!*lambda
     1     *fpolsgn 

      num2(1:ncon)= (4*pi*1.0d-7)*pprimch(1:ncon)*Ib(1:ncon)
     1     *dabs(lambda)*fpolsgn 

      ffprimpre1tmp(1:ncon)=-fpolcon(1:ncon)**2*num1(1:ncon) 
     1     /den(1:ncon)
      ffprimpre2tmp(1:ncon)=-fpolcon(1:ncon)**2*num2(1:ncon) 
     1     /den(1:ncon)

c      write(*,*)'ffprim_pre1',ffprimpre1(1:ncon)
c      write(*,*)'ffprim_pre1tmp',ffprimpre1tmp(1:ncon)

c      write(*,*)'ffprim_pre2',ffprimpre2(1:ncon)
c      write(*,*)'ffprim_pre2tmp',ffprimpre2tmp(1:ncon)
      return
      end

      subroutine findfpol11(ncon,lambda,psicon,Ia,Ib,Ic
     1    ,pprimch, fpolcon,qpsich,cftmq,ffprimpre1)
  
      use arrays, only: fpolsgn, lambdaf, lambdaex, iiter

      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*), fpolcon(*),qpsich(*),mu0
      real *8 ffprimpre1(*)
      real *8 pprimch(*), cftmq(*)
      real *8 Ia(*), Ib(*),Ic(*), chcoeff2(ncon), IaoF(ncon)
      real *8 dIaoF(ncon), den(ncon),  num1(ncon),lambda, num2(ncon)
      real *8 alpha(ncon),IaIcbyq(ncon), dIaIcbyq(ncon)
      real *8 ffprimpre1tmp(ncon)

      integer i,j,k

      pi=4.0d0*datan(1.0d0)
      mu0=pi*4.0d-7
      alpha(1:ncon)=Ia(1:ncon)*Ic(1:ncon)/(qpsich(1:ncon)*2.0d0*pi)**2
c      write(*,*) 'alpha',alpha(1:ncon)
      den(1:ncon)=-(1.0d0+alpha(1:ncon))

      IaIcbyq(1:ncon)=Ia(1:ncon)*Ic(1:ncon)/qpsich(1:ncon)
      call chftransq(chcoeff2, IaIcbyq, ncon, cftmq)
      do i=1,ncon
         call chderiv(1.0d0-psicon(i), ncon, chcoeff2, dIaIcbyq(i))         
      end do
c      write(*,*) 'dIaIcbyq',dIaIcbyq(1:ncon)
      num1(1:ncon)=qpsich(1:ncon)/(lambda*Ic(1:ncon)**2)
     1     *dIaIcbyq(1:ncon)
      num2(1:ncon)=mu0/Ic(1:ncon)*pprimch(1:ncon)*Ib(1:ncon)
c      write(*,*) 'num1',num1(1:ncon)
c      write(*,*) 'lambdaic',lambda*Ic(1:ncon)**2
c     write(*,*) 'num2',num2(1:ncon)
      ffprimpre1(1:ncon)=(num1(1:ncon)+num2(1:ncon))/den(1:ncon)


      IaoF(1:ncon)=Ia(1:ncon)/fpolcon(1:ncon)
      call chftransq(chcoeff2, IaoF, ncon, cftmq)
      do i=1,ncon
         call chderiv(1.0d0-psicon(i), ncon, chcoeff2, dIaoF(i))         
      end do
      dIaoF(1:ncon)=dIaoF(1:ncon)*lambda

      den(1:ncon)=2*pi*qpsich(1:ncon)*fpolcon(1:ncon)
     1     +Ia(1:ncon)/dabs(lambda)*fpolsgn

      num1(1:ncon)=
     1     (dIaoF(1:ncon)*fpolcon(1:ncon)/dabs(lambda) +
     1     mu0*pprimch(1:ncon)*Ib(1:ncon)*dabs(lambda)
     1     )*fpolsgn
      ffprimpre1tmp(1:ncon)=-fpolcon(1:ncon)**2*num1(1:ncon) 
     1     /den(1:ncon)

c     write(*,*)'ffprim_num1',ffprimpre1(1:ncon)
c     write(*,*)'ffprim_num1tmp',ffprimpre1tmp(1:ncon)

      return
      end

      subroutine findBScur(ncon,lambda,Ib,Ic,pprimcon
     1     ,fpolcon,fpass,etai,Zi,nue_star,nui_star,jBScon)
  
      use arrays, only: fout, ijBSmodel

      implicit real*8 (a-h,o-z)
      save

      real *8 Ib(*),Ic(*),pprimcon(*),fpolcon(*),fpass(*)
      real *8 etai(*), Zi(*), nue_star(*), nui_star(*), jBScon(*)
      real *8 x(ncon),D(ncon),L31(ncon),L32(ncon),alpha(ncon)
      real *8 ftrap(ncon), fteff31(ncon),fteff32ee(ncon),fteff32ei(ncon)
      real *8 alpha0(ncon), F32ee(ncon), F32ei(ncon),L34(ncon)
      real *8 A2e(ncon),A2i(ncon),jBSdotB(ncon), Bdotgradphi(ncon)
      real *8 lambda, mu0, pi,fteff34(ncon)

      pi=4.0d0*datan(1.0d0)
      mu0=(4*pi*1.0d-7)
      
c     write (*,*) "Ic=",Ic(1:ncon)

      select case (ijBSmodel)
        
      case (1) ! model by Hirshman [PF 31(1998), 3150]

      x(1:ncon)=(1.0d0-fpass(1:ncon))/fpass(1:ncon)
      D(1:ncon)=1.414d0*Zi(1:ncon)+Zi(1:ncon)**2
     1     +x(1:ncon)*(0.754d0+2.657d0*Zi(1:ncon)+2.0d0*Zi(1:ncon)**2)
     2     +x(1:ncon)**2*(0.348d0+1.243d0*Zi(1:ncon)+Zi(1:ncon)**2)
      L31(1:ncon)=fpolcon(1:ncon)*x(1:ncon)
     1     *(0.754d0+2.21d0*Zi(1:ncon)+Zi(1:ncon)**2
     2     +x(1:ncon)*(0.348d0+1.243d0*Zi(1:ncon)+Zi(1:ncon)**2))
     3     /D(1:ncon)
      L32(1:ncon)=-fpolcon(1:ncon)*x(1:ncon)
     1   *(0.884d0+2.074d0*Zi(1:ncon))/D(1:ncon)
      alpha(1:ncon)=-1.172d0/(1.0d0+0.462d0*x(1:ncon))

      !Assumed Ti=Te
      A2e(1:ncon)=etai(1:ncon)/(1+etai(1:ncon))
     1     *Zi(1:ncon)/(Zi(1:ncon)+1)
      A2i(1:ncon)=A2e(1:ncon)
      jBSdotB(1:ncon)=pprimcon(1:ncon)
     1     *(-L31(1:ncon)*(1.0d0+alpha(1:ncon)/Zi(1:ncon)*A2i(1:ncon))
     2    -L32(1:ncon)*A2e(1:ncon))
      Bdotgradphi(1:ncon)=fpolcon(1:ncon)
     1     *Ic(1:ncon)/Ib(1:ncon)
      jBScon(1:ncon)=jBSdotB(1:ncon)/Bdotgradphi(ncon)


      case(2) ! model by Sauter et. al.[PoP 6(1999), 2834]


      ftrap(1:ncon)= 1.0d0-fpass(1:ncon)
      fteff31(1:ncon)=ftrap(1:ncon)
     1     /(1.0d0+(1.0d0-0.1d0*ftrap(1:ncon))*dsqrt(nue_star(1:ncon))
     2     +0.5d0*(1.0d0-ftrap(1:ncon))*nue_star(1:ncon)/Zi(1:ncon))
      fteff32ee(1:ncon)=ftrap(1:ncon)
     1     /(1.0d0+0.26d0*(1.0d0-ftrap(1:ncon))*dsqrt(nue_star(1:ncon))
     2     +0.18d0*(1.0d0-0.37d0*ftrap(1:ncon))*nue_star(1:ncon)
     3     /dsqrt(Zi(1:ncon)))
      fteff32ei(1:ncon)=ftrap(1:ncon)
     1     /(1+(1+0.6d0*ftrap(1:ncon))*dsqrt(nue_star(1:ncon))
     2     +0.85d0*(1.0d0-0.37d0*ftrap(1:ncon))*nue_star(1:ncon)
     3     *(1.0d0+Zi(1:ncon)))

      alpha0(1:ncon)=-1.17d0*(1.0d0-ftrap(1:ncon))
     1     /(1.0d0-0.22d0*ftrap(1:ncon)-0.19d0*ftrap(1:ncon)**2)

      alpha(1:ncon)=((alpha0(1:ncon)+0.25d0*(1.0d0-ftrap(1:ncon)**2)
     1   *dsqrt(nui_star(1:ncon)))/(1.0d0+0.5d0*dsqrt(nui_star(1:ncon)))
     2   +0.315d0*nui_star(1:ncon)**2*ftrap(1:ncon)**6)
     3   /(1.0d0+0.15d0*nui_star(1:ncon)**2*ftrap(1:ncon)**6)

      F32ee(1:ncon)=(0.05d0+0.62d0*Zi(1:ncon))/Zi(1:ncon)
     1     /(1.0d0+0.44d0*Zi(1:ncon))
     2     *(fteff32ee(1:ncon)-fteff32ee(1:ncon)**4)
     3    +1.0d0/(1.0d0+0.22d0*Zi(1:ncon))
     4     *(fteff32ee(1:ncon)**2-fteff32ee(1:ncon)**4
     5     -1.2d0*(fteff32ee(1:ncon)**3-fteff32ee(1:ncon)**4))
     6    +1.2d0/(1.0d0+0.5d0*Zi(1:ncon))*fteff32ee(1:ncon)**4

      F32ei(1:ncon)=-(0.56d0+1.93d0*Zi(1:ncon))/Zi(1:ncon)
     1     /(1.0d0+0.44d0*Zi(1:ncon))
     2     *(fteff32ei(1:ncon)-fteff32ei(1:ncon)**4)
     3    +4.95d0/(1.0d0+2.48d0*Zi(1:ncon))
     4     *(fteff32ei(1:ncon)**2-fteff32ei(1:ncon)**4
     5     -0.55d0*(fteff32ei(1:ncon)**3-fteff32ei(1:ncon)**4))
     6    -1.2d0/(1.0d0+0.5d0*Zi(1:ncon))*fteff32ei(1:ncon)**4

      L31(1:ncon)=(1.0d0+1.4d0/(1.0d0+Zi(1:ncon)))*fteff31(1:ncon)
     1     -1.9d0/(1.0d0+Zi(1:ncon))*fteff31(1:ncon)**2
     2     +0.3d0/(1.0d0+Zi(1:ncon))*fteff31(1:ncon)**3
     3     +0.2d0/(1.0d0+Zi(1:ncon))*fteff31(1:ncon)**4
      
      L32(1:ncon)=F32ee(1:ncon)+F32ei(1:ncon)

      fteff34(1:ncon)=ftrap(1:ncon)
     1     /(1.0d0+(1.0d0-0.1d0*ftrap(1:ncon))*dsqrt(nue_star(1:ncon))
     2   +0.5d0*(1.0d0-0.5d0*ftrap(1:ncon))*nue_star(1:ncon)/Zi(1:ncon))

      L34(1:ncon)=(1.0d0+1.4d0/(1.0d0+Zi(1:ncon)))*fteff34(1:ncon)
     1     -1.9d0/(1.0d0+Zi(1:ncon))*fteff34(1:ncon)**2
     2     +0.3d0/(1.0d0+Zi(1:ncon))*fteff34(1:ncon)**3
     3     +0.2d0/(1.0d0+Zi(1:ncon))*fteff34(1:ncon)**4
c     write(*,*) 'fteff',fteff31(1:ncon),fteff32ee(1:ncon)
c    1     ,fteff32ei(1:ncon)
c     write(*,*) 'alpha',alpha(1:ncon),F32ee(1:ncon)
c    1     ,F32ei(1:ncon)
c     write(*,*) 'L31',L31(1:ncon),L32(1:ncon)
c    1     ,alpha0(1:ncon)
      !Assumed ne=ni and Ti=Te making Rpe=pe/p=0.5
      Rpe=0.5d0
      jBSdotB(1:ncon)=-fpolcon(1:ncon)*pprimcon(1:ncon)
     1    *etai(1:ncon)/(etai(1:ncon)+1.0d0)
     2    *(L31(1:ncon)/etai(1:ncon)+Rpe*(L31(1:ncon)+L32(1:ncon))
     3     +(1.0d0-Rpe)*(L31(1:ncon)+alpha(1:ncon)*L34(1:ncon)))
      Bdotgradphi(1:ncon)=fpolcon(1:ncon)
     1     *Ic(1:ncon)/Ib(1:ncon)
      jBScon(1:ncon)=jBSdotB(1:ncon)/Bdotgradphi(ncon)

      end select

      write(fout,*) '*************************************'
      write(fout,*) '*****< Boostrap current >*****'
      write(fout,*) 'fraction of passing particles'
      write(fout,*) fpass(1:ncon)
      write(fout,*) '<j_BS dot B>='
      write(fout,*) jBSdotB(1:ncon)
      write(fout,*) '<j_BS dot B>/(B dot grad phi>='
      write(fout,*) jBScon(1:ncon)

      return
      end

      subroutine findOhmcur(ncon,lambda,ne,Te,fpass,Zi
     1     ,nue_star,Vloop,johmcon)
  
      use arrays, only: fout, ijBSmodel

      implicit real*8 (a-h,o-z)
      save

      real *8 ne(*),Te(*),fpass(*)
      real *8 Zi(*),nue_star(*),Vloop(*),johmcon(*)
      real *8 ftrap(ncon), loglambdae(ncon),fteff33(ncon)
      real *8 alpha0(ncon), F32ee(ncon), F32ei(ncon)
      real *8 sigma_sptz(ncon),sigma_neo(ncon)
      real *8 lambda, mu0, pi

      pi=4.0d0*datan(1.0d0)
      mu0=(4*pi*1.0d-7)
      

      select case (ijBSmodel)
        
      case (1) ! model by Hirshman [PF 31(1998), 3150]
         ftrap(1:ncon)= 1.0d0-fpass(1:ncon)

      case(2) ! model by Sauter et. al.[PoP 6(1999), 2834]


      ftrap(1:ncon)= 1.0d0-fpass(1:ncon)
      fteff33(1:ncon)=ftrap(1:ncon)
     1     /(1.0d0+(0.55d0-0.1d0*ftrap(1:ncon))*dsqrt(nue_star(1:ncon))
     2     +0.45d0*(1.0d0-ftrap(1:ncon))*nue_star(1:ncon)
     3     /Zi(1:ncon)**1.5d0)

      ! ne is in [m^-3] and Te is in [eV]
      loglambdae(1:ncon)=31.3-dlog(dsqrt(ne(1:ncon))/Te(1:ncon))
      sigma_sptz(1:ncon)=1.9012d4*Te(1:ncon)**1.5d0/loglambdae(1:ncon)
     1     /(Zi(1:ncon)*(0.58d0+0.74d0/(0.76d0+Zi(1:ncon))))
      sigma_neo(1:ncon)=(1.0d0-(1.0d0+0.36d0/Zi(1:ncon))*fteff33(1:ncon)
     1     +0.59d0/Zi(1:ncon)*fteff33(1:ncon)**2
     2     -0.23d0/Zi(1:ncon)*fteff33(1:ncon)**3)*sigma_sptz(1:ncon)

      johmcon(1:ncon)=sigma_neo(1:ncon)*Vloop(1:ncon)/(2.0d0*pi)
   
      end select

      write(fout,*) '*************************************'
      write(fout,*) '*****< Ohmic current >*****'
      write(fout,*) 'Spitzer conductivity [m^-1 Ohm^-1]='
      write(fout,*) sigma_sptz(1:ncon)
      write(fout,*) 'Neoclassical conductivity [m^-1 Ohm^-1]='
      write(fout,*) sigma_neo(1:ncon)
      write(fout,*) '<j_ohm dot B>/(B dot grad phi>='
      write(fout,*) johmcon(1:ncon)

      return
      end

      subroutine findrhoch(ncon,Rmido, Rmidi)
      use arrays, only: phiichq,psiichq
      use arrays, only: irho, rhoch, dpsiidrhoch, cftmq
      use arrays, only: Rmaxis,Rmax,Rmin
      implicit real*8 (a-h,o-z)
      real *8 chcoeff0(ncon), Rmido(*), Rmidi(*)

     
      if (irho.eq.0) then
         rhoch(1:ncon)=(Rmido(1:ncon)-Rmidi(1:ncon))/(Rmax-Rmin)
      else if (irho.eq.1) then
         rhoch(1:ncon)=(Rmido(1:ncon)-Rmaxis)/(Rmax-Rmaxis)
      else if (irho.eq.2) then
         rhoch(1:ncon)=dsqrt(1.0d0-psiichq(1:ncon))
      else if (irho.eq.3) then
         rhoch(1:ncon)= dsqrt(phiichq(1:ncon))
      end if

      call chftransq(chcoeff0, rhoch, ncon, cftmq)
      do i=1,ncon
         psin=1.0d0-psiichq(i)
         call chderiv(psin,ncon,chcoeff0,drhodpsi)
         dpsiidrhoch(i)=-1.0d0/drhodpsi !used normalized psi
      end do

      write(*,*) 'dpsiidrhoch',dpsiidrhoch(1:ncon)
c     write(*,*) 'Rmido,rhomid',Rmido(1:ncon),rhoch(1:ncon)

      return
      end

      subroutine findtpsifromtrho(ncon)
 
      use arrays, only: psiichq, rhoch, cftmq
      use arrays, only: trho, tpsi, ntpsi
      implicit real*8 (a-h,o-z)
      real *8 chcoeff0(ncon)

      call chftransq(chcoeff0, rhoch, ncon, cftmq)
      eps7=1.0d-14
c     write (*,*) 'trho',trho(1:ntpsi),rhoch(1:ncon)
      do i=1,ntpsi
         if (dabs(trho(i)).lt.eps7) then
            tpsi(i)=1.0d0
         else if (dabs(trho(i)-1.0d0).lt.eps7) then
            tpsi(i)=0.0d0
            exit
         else

         do j=1,ncon
            drho = rhoch(j)-trho(i)
            if (drho.gt.0.0d0)  exit
c            write(*,*) 'j',j,psiichq(1:ncon)
         end do
         tpsitmp=1.0d0-psiichq(j) !initial guess
         iloop=0
         do while (dabs(drho).gt.eps7)
c           write(*,*) 'tpsitmp',iloop,tpsitmp
            call chderiv(tpsitmp,ncon,chcoeff0,drhodpsitmp)
c           write(*,*) 'drhodpsitmp',drho,drhodpsitmp
            tpsitmp=tpsitmp+drho/drhodpsitmp

            if (tpsitmp.lt.0.0d0) then
               tpsitmp=eps7
               trhotmp=eps7
            else if (tpsitmp.gt.1.0d0) then
               tpsitmp=1.0d0-eps7
               trhotmp=1.0d0-eps7
            end if

            call chfit(tpsitmp,ncon,chcoeff0,trhotmp)
c           write(*,*) 'tpsitmp',trhotmp
            drho=trhotmp-trho(i)
            iloop=iloop+1
            if (iloop.gt.10) then
               write(*,*) 'error in finding trho'
               exit
            end if
         end do
         tpsi(i)=1.0d0-tpsitmp
         write(*,*) 'exit of the loop of findtpsi in',iloop

         end if
      end do
      write(*,*) 'tpsi',tpsi(1:ntpsi)
         
      return
      end 

      subroutine findJparfac(ncon,lambda,psicon,Ia,Ib,Ic
     1     ,fpolcon,jpsicon,cftmq,ffprimpre1,ffprimpre2)
  
      use arrays, only: fpolsgn, lambdaf, lambdaex, iiter

      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*),fpolcon(*), jpsicon(*)
      real *8 ffprimpre1(*),ffprimpre2(*)
      real *8 cftmq(*)
      real *8 Ia(*), Ib(*),Ic(*),lambda, mu0, pi

      pi=4.0d0*datan(1.0d0)
      mu0=(4*pi*1.0d-7)
      
c      write (*,*) "Ia=",Ia(1:ncon)
c      write (*,*) "Ib=",Ib(1:ncon)
c      write (*,*) "Ic=",Ic(1:ncon)
      ffprimpre1(1:ncon)=1.0d0/(1.0d0+Ia(1:ncon)/(fpolcon(1:ncon)**2
     1     *Ic(1:ncon)*lambda**2))
c
      ffprimpre2(1:ncon)=ffprimpre1(1:ncon)*Ib(1:ncon)/Ic(1:ncon)

      ffprimpre1(1:ncon)=-ffprimpre1(1:ncon)*jpsicon(1:ncon)
      return
      end

      subroutine findq(ncon,psicon,Rcon,Zcon,tmaxp,icon,fpolcon
     1    ,psiRcon,psiZcon,qcon)

      use arrays, only:  nt2, ksamp2

      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*),psiRcon(*),psiZcon(*)
      real *8 Rcon(*), Zcon(*),tmaxp(*),fpolcon(*)
      real *8 psiRin(nt2*ksamp2+1),psiZin(nt2*ksamp2+1)
      real *8 Rin(nt2*ksamp2+1), Zin(nt2*ksamp2+1),thin(nt2*ksamp2+1)
      real *8 qcon(*)
      integer *4  icon(*)
      integer i,j,k

      do k=1, ncon
         istart=icon(k)
         iend=icon(k+1)-1
         nin=iend-istart+1
         Rin(1:nin)=Rcon(istart:iend)
         Zin(1:nin)=Zcon(istart:iend)
         thin(1:nin)=tmaxp(istart:iend)
         psiRin(1:nin)=psiRcon(istart:iend)
         psiZin(1:nin)=psiZcon(istart:iend)
         
         Rin(nin+1)=Rcon(istart)
         Zin(nin+1)=Zcon(istart)
         thin(nin+1)=tmaxp(istart)+8.0d0*datan(1.0d0)
         psiRin(nin+1)=psiRcon(istart)
         psiZin(nin+1)=psiZcon(istart)

c      write(*,*) 'Rin,Zin',Rin(1:nin+1),Zin(1:nin+1)
c      write(*,*) 'psiRin',psiRin(1:nin+1)
c      write(*,*) 'psiZin',psiZin(1:nin+1)
c      write(*,*) 'dzdt',dZdtmido(k), dZdtmidi(k)
         call findsingleq(nin,Rin, Zin
     1     ,thin,fpolcon(k), psiRin
     2     ,psiZin, qcon(k))
      end do

      return
      end

      subroutine findsingleq(nin,Rin,Zin,thin,fpolin
     1    ,psiRin,psiZin,  qout)

      use arrays, only:  nt2, fpolsgn, ksamp2
      use arrays, only: csol,d1,d2,d3, epstri
      use arrays, only:  psiichq,lambdaex,psiB
      use arrays, only: isymud, Rmaxis, Zmaxis
  
      implicit real*8 (a-h,o-z)
c      save
c      parameter (nint=nt2)
      real *8 Rin(*), Zin(*),thin(*), qout
    
      real *8 Rint(nin+1),Zint(nin+1), dldt(nin+1)
 
      real *8 Rintmp(nin+1), Zintmp(nin+1),thintmp(nin+1)
      real *8 psiRin(*),psiZin(*), BpR2(nin+1), tttmp(nin)

      real *8 temp0(nin+1), temp1(nin+1), temp2(nin+1)
      real *8 temp3(nin+1), temp4(nin+1), temp5(nin+1)
      real *8 temp6(nin+1), temp7(nin+1), temp8(nin+1)
  
      real *8 thint(nin+1),integ11(nin+1)
      real *8 integ12(nin+1),integ13(nin+1)
      real *8 integ14(nin+1),integ15(nin+1)

      integer *4 i,j,k, kLag2, iw, itype, nint

      save
      
      nint=nin !nt2*ksamp2/2
      pi = 4.0d0*datan(1.0d0)
c      epstri= 1.0d-10

      delth = 0.0d0
      do i=1,nint
         thint(i)=2.0d0*pi*(i-1)/nint+thin(1) !delth
      end do


      kLag2=4
      if (isymud.eq.1) then
         do i=1,nin+1
            thintmp(i)=thin(i)
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
      else
         do i=1,nin+1
            thintmp(i)=thin(i)/2.0d0
            Rintmp(i)=Rin(i)
            zintmp(i)=zin(i)
         end do
         tttmp(1:nint)=thint(1:nint)/2.0d0
      end if
 

      BpR2(1:nin+1)=(psiRin(1:nin+1)**2.0d0+psiZin(1:nin+1)**2.0d0)
      temp0(1:nin+1)=1.0d0/
     1     (((Zin(1:nin+1)-Zmaxis)/(Rin(1:nin+1)-Rmaxis))**2.0d0+1.0d0)
      temp1(1:nin+1)=-temp0(1:nin+1)*(Zin(1:nin+1)-Zmaxis)
     1     /(Rin(1:nin+1)-Rmaxis)**2.0d0   !dtheta/dR
      temp2(1:nin+1)=temp0(1:nin+1)
     1     /(Rin(1:nin+1)-Rmaxis)   !dtheta/dZ

      temp3(1:nin+1)=dabs((psiRin(1:nin+1)*temp2(1:nin+1)
     1     -psiZin(1:nin+1)*temp1(1:nin+1))/Rin(1:nin+1)) !(grad phi times grad psi) cdot grad theta)

      temp6(1:nin+1)=1.0d0/temp3(1:nin+1)/Rin(1:nin+1)**2 

      if (isymud.eq.1) then

         call IntTriCos(nin/2+1,thintmp,temp6,nint/2+1,
     1        thint,integ13,epstri) 
     
         do i=1,nint/2
            integ13(nint+2-i)=integ13(i)
         end do

      else                      ! if isymud is not 1 (up-down assymmetric)
        
         call IntTriCos(nin+1,thintmp,temp6,nint+1,
     1        tttmp,integ13,epstri)

      end if                    ! end if isymud.eq.1

      qout=0.0d0
      do i=1,nint
         qout=qout+integ13(i)/nint            
      end do
c      write(*,*) 'nint,qout',nint,qout,fpolin
      qout=qout*fpolin*fpolsgn

c      write(*,*) 'qout',qout

      return
      end

      subroutine findpsioffaxis(nt1, dt, rzdimmin, zc1, dzdtc1, 
     1     wc1,  Roff, Zoff, iexterior, rndm, tndm, psi)

      use arrays, only: nt2, nr, rnd, kcheb, ucoeff, nsub
      use arrays, only: cftmsub, bnodes, Rmap0, Zmap0
      use arrays, only: R0,Z0,reps,rkappa, delta, ibtype

      implicit real*8 (a-h,o-z)
      
      integer *4 iexterior, nt1
      real *8  dt, Roff, Zoff, rzdimmin, rndm, tndm, um, psi
      real *8  ucinr(nt2), ucini(nt2)  
      real *8  ucoeffsubj(kcheb),chcoeff(kcheb), uch(kcheb)
      complex *16 zc1(*), dzdtc1(*), wc1(*), wc2(nt1)
      complex *16 ima, zc0, w0, w2pi,zout
      complex *16 ucin(nt2)
     
      ima = (0.0d0,1.0d0)
      
      pi = 4.0d0*datan(1.0d0)
      
      zc0=Roff+ Zoff*ima 
      dist=dsqrt((Roff-Rmap0)**2+(Zoff-Zmap0)**2)
c      write(*,*) 'dist,rzdimmin',dist,rzdimmin
c      write(*,*) 'zc0',zc0
c      write(*,*) 'zc1',zc1(1:nt1)
c      write(*,*) 'dzdtc1',dzdtc1
c      write(*,*) 'wc1',wc1
      wc2 = 1.0d0
      call cnf_cauchy(nt1,dt,zc0,zc1,dzdtc1,wc1,w0)
      w0r=w0
      w0i=-ima*w0
      call cnf_cauchy(nt1,dt,zc0,zc1,dzdtc1,wc2,w2pi)
      w2pir=w2pi
      w2pii=-ima*w2pi
c      write(*,*) 'w0',w0
c      write(*,*) 'w2pi',w2pi
      !if ((isNaN(w0r).eq.1).or.(isNan(w0i).eq.1)) then
      !   psi=-1.0d0
      !   rndm=-1.0d0
      !   tndm=0.0d0
      !   iexterior= 1
c         write(*,*) 'exterior point1'
      !   return
      !end if

      rndm=dsqrt(w0r**2+w0i**2)
      tndm=datan2(w0i,w0r) 

          write(*,*) 'zc0',zc0
      write(*,*) 'rndm,tndm',w0r,w0i,rndm,tndm
     
      if (ibtype.eq.1) then
         alpha=dasin(delta)
c        write(*,*) 'alpha',R0,Z0,reps,rkappa, delta, alpha 
         if ((zoff.lt.(Z0+rkappa*reps*R0)).and.
     1      (zoff.gt.(Z0-rkappa*reps*R0))) then
         t=dasin((Zoff-Z0)/R0)/reps/rkappa
         Rmax=R0*(1+reps*dcos(t+alpha*dsin(t)))
         Rmin=R0*(1+reps*dcos(pi-t+alpha*dsin(pi-t)))

         write(*,*) 't',t, Rmax, Rmin
         if ((Roff.gt.Rmin).and.(Roff.lt.Rmax)) then
              iexterior= 0 !interior

              write(*,*) 'possible interior point: rnd,tnd',rndm,tndm
              w0r2=w0/w2pi
              w0i2=-ima*w0/w2pi
              rndm=dsqrt(w0r2**2+w0i2**2)
              tndm=datan2(w0i2,w0r2)
              if (rndm.ge.1.0d0) then
                 psi=-1.0d0
                 rndm=-1.0d0
                 tndm=0.0d0
                 iexterior= 4
                 write(*,*) 'exterior point4' !make cauchy integral 0
                 return
              end if

         write(*,*) 'interior point(corrected): rnd,tnd',rndm,tndm

              do i=1, nr-1
                 if (rndm.le.rnd(1)) then
               isub=1
               xch=(rndm-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
c               write(*,*) 'very close to origin: smaller than rnd(1)'
               exit
            else if((rnd(i).le.rndm).and.(rnd(i+1).gt.rndm)) then
               isub=(i-1)/kcheb+1
               xch=(rndm-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
               if (xch.gt.1.0d0) then
                  isub=isub+1
                  xch=(rndm-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
               end if
               exit
            else if ((rnd(nr).lt.rndm)) then
               isub=nsub
               xch=(rndm-bnodes(nsub))/(bnodes(nsub+1)-bnodes(nsub))
               exit
            end if
              end do
    
         else ! R is outside [Rmin,Rmax]
            psi=0.0d0
         rndm=-1.0d0
         tndm=0.0d0
         iexterior= 5
         write(*,*) 'exterior point5'
         return
         end if
      
        else ! z is outside [zmin,zmax]

         psi=0.0d0
         rndm=-1.0d0
         tndm=0.0d0
         iexterior= 6
         write(*,*) 'exterior point6'
         return
           end if


        else ! if not ibtype.eq.1
           
      if (rndm.ge.1.0d0) then
         psi=-1.0d0
         rndm=-1.0d0
         tndm=0.0d0
         iexterior= 2
c         write(*,*) 'check: w2pir=0',w2pi
c         write(*,*) 'exterior point2'
         return
      else if ((rndm.lt.3.0d-1).and.(dist.gt. rzdimmin*0.6)) then 
         psi=-1.0d0
         rndm=-1.0d0
         tndm=0.0d0
         iexterior= 3
c         write(*,*) 'check: w2pir=0',w2pi
c         write(*,*) 'exterior point3' !make cauchy integral 0
         return
      else
c         write(*,*) 'check: (w2pir-1.0d0)=0',w2pi
         write(*,*) 'possible interior point: rnd,tnd',rndm,tndm
         w0r2=w0/w2pi
         w0i2=-ima*w0/w2pi
         rndm=dsqrt(w0r2**2+w0i2**2)
         tndm=datan2(w0i2,w0r2)
     
         if (rndm.ge.1.0d0) then
            psi=-1.0d0
            rndm=-1.0d0
            tndm=0.0d0
            iexterior= 4
c            write(*,*) 'exterior point4' !make cauchy integral 0
            return
         end if

         write(*,*) 'interior point(corrected): rnd,tnd',rndm,tndm
         iexterior= 0 !interior
         do i=1, nr-1
            if (rndm.le.rnd(1)) then
               isub=1
               xch=(rndm-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
c               write(*,*) 'very close to origin: smaller than rnd(1)'
               exit
            else if((rnd(i).le.rndm).and.(rnd(i+1).gt.rndm)) then
               isub=(i-1)/kcheb+1
               xch=(rndm-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
               if (xch.gt.1.0d0) then
                  isub=isub+1
                  xch=(rndm-bnodes(isub))/(bnodes(isub+1)-bnodes(isub))
               end if
               exit
            else if ((rnd(nr).lt.rndm)) then
               isub=nsub
               xch=(rndm-bnodes(nsub))/(bnodes(nsub+1)-bnodes(nsub))
               exit
            end if
         end do
      end if


      end if                    ! ibtype.eq.1
      
      write(*,*) 'fine1',isub,xch
      do i=1,kcheb
         do j=1,nt2
            inext=((isub-1)*kcheb+i)+(j-1)*nr
            ucin(j)=ucoeff(inext)
         end do
         write(*,*) 'fine2',ucin(1:nt2)
         call fcoffth(ucin,nt2,tndm,zout)
         uch(i) = dreal(zout)
      end do
      write(*,*) 'psi00',uch(1:kcheb)
      call chftransq(chcoeff,uch,kcheb,cftmsub)
      call chfit(xch, kcheb, chcoeff, um)


      write(*,*) 'psi00',uch(1:kcheb),um, Roff
      psi = um*dsqrt(Roff)
c      write(*,*) 'psi0',psi,rndm,xch !,ucin(1:nt2)
     
      return
      end

      subroutine findq_solovev(ncon,psicon,fpolcon,qcon)
   
      use arrays, only:  csol,reps,rkappa,delta,d1,d2,d3
      use arrays, only: fpolsgn
      implicit real*8 (a-h,o-z)
      save

      real *8 psicon(*),fpolcon(*),qcon(*)
      real *8 c0(ncon),rmin2(ncon),rmax2(ncon)
      real *8 rmin(ncon),rmax(ncon)

      pi = 4.0d0*datan(1.0d0)
      do i=1,ncon
         c0(i) = d1 - psicon(i);
      end do
      c4 = csol/8.0d0 + d3
    
      do i=1,ncon
         if ((d2**2 - 4*c0(i)*c4).ge.0.0d0) then
            rmin2(i)= (-d2 - dsqrt(d2**2 - 4*c0(i)*c4))
     1           /(2.0d0*c4)
            rmax2(i)= (-d2 + dsqrt(d2**2 - 4*c0(i)*c4))
     1           /(2.0d0*c4)
         else
            rmin2(i)=-d2/(2.0d0*c4)
            rmax2(i)=-d2/(2.0d0*c4)
         end if
      end do
      rmin(1:ncon) = dsqrt(rmin2(1:ncon))
      rmax(1:ncon) = dsqrt(rmax2(1:ncon))
      errtol=1d-15
      do i=1,ncon
         ekorder=dsqrt(1-rmin2(i)/rmax2(i))
         factor = dabs(fpolcon(i))/
     1        (4*pi*dsqrt(-d3*c4)*rmax(i)*rmin(i)**2)
c         write(*,*) 'factor',factor,korder
         ellipticE=RF(0.0d0,1-ekorder**2,1.0d0,errtol,ierr)
     1        -ekorder**2*RD(0.0d0,1-ekorder**2,1.0d0,errtol,ierr)/3.0d0
         qcon(i)=factor*ellipticE !*fpolsgn
      end do
c      write(*,*) 'ellipt',Rf(0.0d0,0.5d0,0.0d0,1d-3,ierr)
c      write(*,*) 'qcon',qcon(1:ncon)
c      write(*,*) 'fpol',fpolcon(1:ncon)

      return
      end

C
C
C     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C ACM ALGORITHM 577
C
C ALGORITHMS FOR INCOMPLETE ELLIPTIC INTEGRALS
C
C BY B.C. CARLSON AND E.M. NOTIS
C
C ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE, SEPTEMBER, 1981.
C          ********************************************************
C
      DOUBLE PRECISION FUNCTION RF(X,Y,Z,ERRTOL,IERR)
C
C          THIS FUNCTION SUBROUTINE COMPUTES THE INCOMPLETE ELLIPTIC
C          INTEGRAL OF THE FIRST KIND,
C          RF(X,Y,Z) = INTEGRAL FROM ZERO TO INFINITY OF
C
C                                -1/2     -1/2     -1/2
C                      (1/2)(T+X)    (T+Y)    (T+Z)    DT,
C
C          WHERE X, Y, AND Z ARE NONNEGATIVE AND AT MOST ONE OF THEM
C          IS ZERO.  IF ONE OF THEM IS ZERO, THE INTEGRAL IS COMPLETE.
C          THE DUPLICATION THEOREM IS ITERATED UNTIL THE VARIABLES ARE
C          NEARLY EQUAL, AND THE FUNCTION IS THEN EXPANDED IN TAYLOR
C          SERIES TO FIFTH ORDER.  REFERENCE: B. C. CARLSON, COMPUTING
C          ELLIPTIC INTEGRALS BY DUPLICATION, NUMER. MATH. 33 (1979),
C          1-16.  CODED BY B. C. CARLSON AND ELAINE M. NOTIS, AMES
C          LABORATORY-DOE, IOWA STATE UNIVERSITY, AMES, IOWA 50011.
C          MARCH 1, 1980.
C
C          CHECK BY ADDITION THEOREM: RF(X,X+Z,X+W) + RF(Y,Y+Z,Y+W)
C          = RF(0,Z,W), WHERE X,Y,Z,W ARE POSITIVE AND X * Y = Z * W.
C
      INTEGER IERR,PRINTR
      DOUBLE PRECISION C1,C2,C3,E2,E3,EPSLON,ERRTOL,LAMDA
      DOUBLE PRECISION LOLIM,MU,S,UPLIM,X,XN,XNDEV,XNROOT
      DOUBLE PRECISION Y,YN,YNDEV,YNROOT,Z,ZN,ZNDEV,ZNROOT
C          INTRINSIC FUNCTIONS USED: DABS,DMAX1,DMIN1,DSQRT
C
C          PRINTR IS THE UNIT NUMBER OF THE PRINTER.
      DATA PRINTR/6/
C
C          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
C          LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
C          UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.
C
C          ACCEPTABLE VALUES FOR:   LOLIM      UPLIM
C          IBM 360/370 SERIES   :   3.D-78     1.D+75
C          CDC 6000/7000 SERIES :   1.D-292    1.D+321
C          UNIVAC 1100 SERIES   :   1.D-307    1.D+307
C
C          WARNING: IF THIS PROGRAM IS CONVERTED TO SINGLE PRECISION,
C          THE VALUES FOR THE UNIVAC 1100 SERIES SHOULD BE CHANGED TO
C          LOLIM = 1.E-37 AND UPLIM = 1.E+37 BECAUSE THE MACHINE
C          EXTREMA CHANGE WITH THE PRECISION.
C
      DATA LOLIM/3.D-78/, UPLIM/1.D+75/
C
C          ON INPUT:
C
C          X, Y, AND Z ARE THE VARIABLES IN THE INTEGRAL RF(X,Y,Z).
C
C          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
C          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
C          ERRTOL ** 6 / (4 * (1 - ERRTOL)).
C
C          SAMPLE CHOICES:  ERRTOL   RELATIVE TRUNCATION
C                                    ERROR LESS THAN
C                           1.D-3    3.D-19
C                           3.D-3    2.D-16
C                           1.D-2    3.D-13
C                           3.D-2    2.D-10
C                           1.D-1    3.D-7
C
C          ON OUTPUT:
C
C          X, Y, Z, AND ERRTOL ARE UNALTERED.
C
C          IERR IS THE RETURN ERROR CODE:
C               IERR = 0 FOR NORMAL COMPLETION OF THE SUBROUTINE,
C               IERR = 1 FOR ABNORMAL TERMINATION.
C
C          ********************************************************
C          WARNING: CHANGES IN THE PROGRAM MAY IMPROVE SPEED AT THE
C          EXPENSE OF ROBUSTNESS.
C
      IF (DMIN1(X,Y,Z) .LT. 0.D0) GO TO 100
      IF (DMIN1(X+Y,X+Z,Y+Z) .LT. LOLIM) GO TO 100
      IF (DMAX1(X,Y,Z) .LE. UPLIM) GO TO 112
  100 WRITE(PRINTR,104)
  104 FORMAT(1H0,42H*** ERROR - INVALID ARGUMENTS PASSED TO RF)
      WRITE(PRINTR,108) X,Y,Z
  108 FORMAT(1H ,4HX = ,D23.16,4X,4HY = ,D23.16,4X,4HZ = ,D23.16)
      IERR = 1
      GO TO 124
C
  112 IERR = 0
      XN = X
      YN = Y
      ZN = Z
C
  116 MU = (XN + YN + ZN) / 3.D0
      XNDEV = 2.D0 - (MU + XN) / MU
      YNDEV = 2.D0 - (MU + YN) / MU
      ZNDEV = 2.D0 - (MU + ZN) / MU
      EPSLON = DMAX1(DABS(XNDEV),DABS(YNDEV),DABS(ZNDEV))
      IF (EPSLON .LT. ERRTOL) GO TO 120
      XNROOT = DSQRT(XN)
      YNROOT = DSQRT(YN)
      ZNROOT = DSQRT(ZN)
      LAMDA = XNROOT * (YNROOT + ZNROOT) + YNROOT * ZNROOT
      XN = (XN + LAMDA) * 0.25D0
      YN = (YN + LAMDA) * 0.25D0
      ZN = (ZN + LAMDA) * 0.25D0
      GO TO 116
C
  120 C1 = 1.D0 / 24.D0
      C2 = 3.D0 / 44.D0
      C3 = 1.D0 / 14.D0
      E2 = XNDEV * YNDEV - ZNDEV * ZNDEV
      E3 = XNDEV * YNDEV * ZNDEV
      S = 1.D0 + (C1 * E2 - 0.1D0 - C2 * E3) * E2 + C3 * E3
      RF = S / DSQRT(MU)
C
  124 RETURN
      END
C
C
C     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
      DOUBLE PRECISION FUNCTION RD(X,Y,Z,ERRTOL,IERR)
C
C          THIS FUNCTION SUBROUTINE COMPUTES AN INCOMPLETE ELLIPTIC
C          INTEGRAL OF THE SECOND KIND,
C          ellipticE(X,Y,Z) = INTEGRAL FROM ZERO TO INFINITY OF
C
C                                -1/2     -1/2     -3/2
C                      (3/2)(T+X)    (T+Y)    (T+Z)    DT,
C
C          WHERE X AND Y ARE NONNEGATIVE, X + Y IS POSITIVE, AND Z IS
C          POSITIVE.  IF X OR Y IS ZERO, THE INTEGRAL IS COMPLETE.
C          THE DUPLICATION THEOREM IS ITERATED UNTIL THE VARIABLES ARE
C          NEARLY EQUAL, AND THE FUNCTION IS THEN EXPANDED IN TAYLOR
C          SERIES TO FIFTH ORDER.  REFERENCE: B. C. CARLSON, COMPUTING
C          ELLIPTIC INTEGRALS BY DUPLICATION, NUMER. MATH. 33 (1979),
C          1-16.  CODED BY B. C. CARLSON AND ELAINE M. NOTIS, AMES
C          LABORATORY-DOE, IOWA STATE UNIVERSITY, AMES, IOWA 50011.
C          MARCH 1, 1980..
C
C          CHECK: RD(X,Y,Z) + RD(Y,Z,X) + RD(Z,X,Y)
C          = 3 / DSQRT(X * Y * Z), WHERE X, Y, AND Z ARE POSITIVE.
C
      INTEGER IERR,PRINTR
      DOUBLE PRECISION C1,C2,C3,C4,EA,EB,EC,ED,EF,EPSLON,ERRTOL,LAMDA
      DOUBLE PRECISION LOLIM,MU,POWER4,SIGMA,S1,S2,UPLIM,X,XN,XNDEV
      DOUBLE PRECISION XNROOT,Y,YN,YNDEV,YNROOT,Z,ZN,ZNDEV,ZNROOT
C          INTRINSIC FUNCTIONS USED: DABS,DMAX1,DMIN1,DSQRT
C
C          PRINTR IS THE UNIT NUMBER OF THE PRINTER.
      DATA PRINTR/6/
C
C          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
C          LOLIM IS NOT LESS THAN 2 / (MACHINE MAXIMUM) ** (2/3).
C          UPLIM IS NOT GREATER THAN (0.1 * ERRTOL / MACHINE
C          MINIMUM) ** (2/3), WHERE ERRTOL IS DESCRIBED BELOW.
C          IN THE FOLLOWING TABLE IT IS ASSUMED THAT ERRTOL WILL
C          NEVER BE CHOSEN SMALLER THAN 1.D-5.
C
C          ACCEPTABLE VALUES FOR:   LOLIM      UPLIM
C          IBM 360/370 SERIES   :   6.D-51     1.D+48
C          CDC 6000/7000 SERIES :   5.D-215    2.D+191
C          UNIVAC 1100 SERIES   :   1.D-205    2.D+201
C
C          WARNING: IF THIS PROGRAM IS CONVERTED TO SINGLE PRECISION,
C          THE VALUES FOR THE UNIVAC 1100 SERIES SHOULD BE CHANGED TO
C          LOLIM = 1.E-25 AND UPLIM = 2.E+21 BECAUSE THE MACHINE
C          EXTREMA CHANGE WITH THE PRECISION.
C
      DATA LOLIM/6.D-51/, UPLIM/1.D+48/
C
C          ON INPUT:
C
C          X, Y, AND Z ARE THE VARIABLES IN THE INTEGRAL RD(X,Y,Z).
C
C          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
C          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
C          3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2.
C
C          SAMPLE CHOICES:  ERRTOL   RELATIVE TRUNCATION
C                                    ERROR LESS THAN
C                           1.D-3    4.D-18
C                           3.D-3    3.D-15
C                           1.D-2    4.D-12
C                           3.D-2    3.D-9
C                           1.D-1    4.D-6
C
C          ON OUTPUT:
C
C          X, Y, Z, AND ERRTOL ARE UNALTERED.
C
C          IERR IS THE RETURN ERROR CODE:
C               IERR = 0 FOR NORMAL COMPLETION OF THE SUBROUTINE,
C               IERR = 1 FOR ABNORMAL TERMINATION.
C
C          ********************************************************
C          WARNING: CHANGES IN THE PROGRAM MAY IMPROVE SPEED AT THE
C          EXPENSE OF ROBUSTNESS.
C
      IF (DMIN1(X,Y) .LT. 0.D0) GO TO 100
      IF (DMIN1(X+Y,Z) .LT. LOLIM) GO TO 100
      IF (DMAX1(X,Y,Z) .LE. UPLIM) GO TO 112
  100 WRITE(PRINTR,104)
  104 FORMAT(1H0,42H*** ERROR - INVALID ARGUMENTS PASSED TO RD)
      WRITE(PRINTR,108) X,Y,Z
  108 FORMAT(1H ,4HX = ,D23.16,4X,4HY = ,D23.16,4X,4HZ = ,D23.16)
      IERR = 1
      GO TO 124
C
  112 IERR = 0
      XN = X
      YN = Y
      ZN = Z
      SIGMA = 0.D0
      POWER4 = 1.D0
C
  116 MU = (XN + YN + 3.D0 * ZN) * 0.2D0
      XNDEV = (MU - XN) / MU
      YNDEV = (MU - YN) / MU
      ZNDEV = (MU - ZN) / MU
      EPSLON = DMAX1(DABS(XNDEV),DABS(YNDEV),DABS(ZNDEV))
      IF (EPSLON .LT. ERRTOL) GO TO 120
      XNROOT = DSQRT(XN)
      YNROOT = DSQRT(YN)
      ZNROOT = DSQRT(ZN)
      LAMDA = XNROOT * (YNROOT + ZNROOT) + YNROOT * ZNROOT
      SIGMA = SIGMA + POWER4 / (ZNROOT * (ZN + LAMDA))
      POWER4 = POWER4 * 0.25D0
      XN = (XN + LAMDA) * 0.25D0
      YN = (YN + LAMDA) * 0.25D0
      ZN = (ZN + LAMDA) * 0.25D0
      GO TO 116
C
  120 C1 = 3.D0 / 14.D0
      C2 = 1.D0 / 6.D0
      C3 = 9.D0 / 22.D0
      C4 = 3.D0 / 26.D0
      EA = XNDEV * YNDEV
      EB = ZNDEV * ZNDEV
      EC = EA - EB
      ED = EA - 6.D0 * EB
      EF = ED + EC + EC
      S1 = ED * (- C1 + 0.25D0 * C3 * ED - 1.5D0 * C4 * ZNDEV * EF)
      S2 = ZNDEV * (C2 * EF + ZNDEV * (- C3 * EC + ZNDEV * C4 * EA))
      RD = 3.D0 * SIGMA + POWER4 * (1.D0 + S1 + S2) / 
     1     (MU * DSQRT(MU))
C
  124 RETURN
      END
C
C
C     $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

