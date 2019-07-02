!--------------------------------------------------------
! solvej.f: This subroutine is to solve non-linear Grad-Shafranov equation 
! by iterating the Poisson solver in an unit disk using jpar inputs
!--------------------------------------------------------

      subroutine solvej

      use arrays, only:  allocarray3, iptype, file_bc, csol
      use arrays, only:  Rt3, Zt3, rnd, tnd, dzdw3, dzdww3
      use arrays, only:  nt1, ksamp,  nt2, eps, epsiter, maxiter
      use arrays, only:  nsub, kcheb, kLag, nr, ntot, blength, bnodes
      use arrays, only:  Rt1, Zt1, T1, dRt1, dZt1, Rmid, Rmaxis, Zmaxis
      use arrays, only:  Rt2, Zt2, T2, dRt2, dZt2, dwdz2, dwdzz2
      use arrays, only:  dwdz3, dwdzz3, cftmsub, Rmap0, Zmap0

      use arrays, only:  csol,reps,rkappa,delta,d1,d2,d3
      use arrays, only:  npsi, psi0, psiB, psiEF, psiiEF 
      use arrays, only:  rhoEF, presEF, chcoeffpp, chcoeffpp0, chcoeffff 
      use arrays, only:  pprimEF,ffprimEF,chcoeffprho,chcoefff2rho
      use arrays, only:  fpolEF, qpsiEF, bpprim, cpprim, dpprim
      use arrays, only:  bffprim, cffprim, dffprim, fpolch, rhoch
      use arrays, only:  dpsiidrhoch, iptable, isplch
      use arrays, only:  nchq, psiichq, cftmq, spbx, spxb, spdef
      use arrays, only:  nchq0, qpsich0, psiichq0, cftmq0
      use arrays, only:  spbx0, spxb0, spdef0, pprimch0, fpolch0
      use arrays, only:  qpsich, pprimch,ffprimch, qpsich1
      use arrays, only:  jpsich, jpsich0
      use arrays, only:  wffprim, fpolcon2, ffprimpre, ffprimpre1
      use arrays, only:  ffprimpre2, ffprimgoal

      use arrays, only:  iiter, iremap, iiterl, lambdaex
      use arrays, only:  g, psii, f, lambda, maxpsi, maxdpsi, lambdaf
      use arrays, only:  epsmaxdist, Rmaxpsi, Zmaxpsi, fpolsgn, psiex 
      use arrays, only:  dpsidr, dpsidrr, u, ur, uth, urr, urt, utt
      use arrays, only: q0, F0, p0, pin, pout, FF0, FFin, FFout
      use arrays, only: qfac,qpin, qpout,Ffac, Fpin, Fpout
      use arrays, only: jpar0, jpin, jpout, Jomax, Joloc, Jow
      use arrays, only: itftype, mach, ptf0, ptfin, ptfout
      use arrays, only: fiter, ftime, fout, isymud
      use arrays, only: chcoeffp, presch, chcoefftf, ptfprimch, ptflowch
      use arrays, only: ptflowEF, btflow, ctflow, dtflow, iscale 
      use arrays, only: ptfmax, ptfloc, ptfw
      use arrays, only: rndm, tndm, dpsidrm, dpsidtm, dpsidrrm ,dpsidrtm
      use arrays, only: dpsidttm, Rt3m, Zt3m, psim

      use arrays, only: iqtype, nqpsi, psiiq, qpsiin, bqpsi,cqpsi,dqpsi
      use arrays, only: ijtype, njpsi, psiij, jpsiin, bjpsi,cjpsi,djpsi
      use arrays, only: iconvj, maxdpsiq, relaxiter
      use arrays, only: nchy, psiichy, cftmy, spybx, spyxb,spydef

      use arrays, only: isolver, halpha, ifindmagaxis, dpsistartmag
      use arrays, only: dense0, denseb, rn1de, rn2de, ne, etai
      use arrays, only: ti0, tib, rn1ti, rn2ti
      use arrays, only: te0, teb, rn1te, rn2te, Te, nue_star, nui_star 
c
c     read the fixed boundary from a file
      implicit double precision(a-h, o-z)
c      parameter (nmxsub = 200, kmax = 32)
c      parameter (nr = nsub * kcheb)
c      parameter (ntmax = 300)

       
 
c      parameter (ntcon = ntheta*ncon)
      real *8 fsol(ntot)
      real *8 psif(ntot)
   
 
      real *8 dRdr(ntot), dRdrr(ntot)
c      real *8 presch(nchq)
      real *8 chcoeff(nchq),chcoeffpre1(nchq),chcoeffpre2(nchq)
      real *8 chcoeffgoal(nchq)
      real *8 bfpol(npsi),cfpol(npsi),dfpol(npsi)
      real *8 bpres(npsi),cpres(npsi),dpres(npsi) 
      real *8 preschrho(nchq),fpol2chrho(nchq)
      real *8 psiiEF2(npsi)

      real *8 psicontmp(nchq), psicontmp0(nchq0), psiichsub(kcheb)
      real *8 fpolcon(nchq)
      real *8 chcoeff0(nchq0)
      real *8 mu0
      complex *16 ima
      
c      save


      write(*,*) 'start poisson, nt2',nt2

      write(*,*) 'start poisson, korder,nr,ntot',kcheb,nr,ntot
      pi = 4.0d0*datan(1.0d0)
      mu0=4.0d0*pi*1.0d-7
      ima = (0,1)

      ntheta=nt2

      korder = kcheb
    

      if (iiter.eq.1) then !initial condition for the first iteration
         lambda = -10.0          ! Initial guess for lambda. Sign should be correct.
         do i=1,ntheta
            g(i)=0.0d0          !psii at the fixed Boundary  
         end do
         do j = 1,ntheta
            do i = 1,nr
               inext = i + (j-1)*nr
               psii(inext)=1.0d0-rnd(i)**2 !intial guess for psii
            end do
         end do

         call chsetupsub(kcheb, psiichsub, cftmsub)  
         call chsetupq(nchq0, psiichq0, cftmq0, spbx0, spxb0, spdef0)
         write(*,*) 'psiichq0',psiichq0(1:nchq0)
c     For several intial iterations, try to find contour integrals 
c     only for 6 flux surfaces because the solutions are not accurate 
c     The number of flux surfaces will be reconvered after iterations iiter>iiterq
         if (nchq.gt.6) nchq=6            
         call chsetupq(nchq, psiichq, cftmq, spbx, spxb, spdef)
         write(*,*) 'psiichq',psiichq(1:nchq)
         select case (iptype)   ! Option for pressure and current profiles
         case(0)                ! solovev solution 
            pprimch0(1:nchq0)=-csol/mu0
            fpolch0(1:nchq0)= dabs(F0)*fpolsgn
            ffprimch(1:nchq0)= 0.0d0

            pprimch(1:nchq)=-csol/mu0
            fpolch(1:nchq)= dabs(F0)*fpolsgn
            lambdaex=1.0d0/(d1-d2**2.0d0/(4.0d0*(csol/8.0d0+d3)))
            write(*,*) 'lambdaex=',lambdaex
            iqtype = 0
           
         case(1)
            do i=1,nchq0
               psinchq=1.0d0-psiichq0(i)
               pprimch0(i)=p0*(1.0d0-psinchq**pin)**pout/(mu0*Rmid**2)
               ffprimch(i)=FF0*(1-psinchq**pin)**pout
c               fpolch0(i)=F0*(1.0d0+Ffac*psinchq**Fpin)**Fpout
            end do
            do i=1,nchq
               psinchq=1.0d0-psiichq(i)
               pprimch(i)=p0*(1.0d0-psinchq**pin)**pout/(mu0*Rmid**2)
c               fpolch(i)=F0*(1.0d0+Ffac*psinchq**Fpin)**Fpout
            end do
            fpolch(1:nchq)= dabs(F0)*fpolsgn !initial guess
         case(2,3) ! profile by Table or EFIT
            if (iptable.eq.0) then
               call spline(-psiiEF, pprimEF     
     1              , bpprim, cpprim, dpprim, npsi) 
            
               call spline(-psiiEF, ffprimEF
     1              , bffprim, cffprim, dffprim,npsi)   
           
               do i= 1, nchq0
                  call ispline(-psiichq0(i),-psiiEF
     1              , pprimEF, bpprim ,cpprim,dpprim,npsi,pprimch0(i))
                  call ispline(-psiichq0(i),-psiiEF
     1           , ffprimEF, bffprim ,cffprim,dffprim,npsi,ffprimch(i))
               end do
               do i= 1, nchq
                  call ispline(-psiichq(i),-psiiEF
     1                 , pprimEF, bpprim ,cpprim,dpprim,npsi,pprimch(i))
               end do
            
               call ispline(-1.0d0,-psiiEF
     1              , pprimEF, bpprim ,cpprim,dpprim,npsi,p0)
               call ispline(-1.0d0,-psiiEF
     1              , ffprimEF, bffprim ,cffprim,dffprim,npsi,FF0)
c               call chftransq(chcoeff0,qpsich0,nchq0,cftmq0)
c               do i=1,nchq
c                  call chfit(1-psiichq(i), nchq0, chcoeff0,qpsich(i))
c               end do
            else if (iptable.eq.1) then
               isplch=1  !interpolate presure proflie to chebyshev points 
                         ! and then use it to find dpdpsi
               call spline(-psiiEF, presEF     
     1              , bpres, cpres, dpres, npsi) 
            
               call spline(-psiiEF, ffprimEF
     1              , bffprim, cffprim, dffprim,npsi)

               do i= 1, nchq
                  call ispline(-psiichq(i),-psiiEF
     1                 , presEF, bpres ,cpres,dpres,npsi,presch(i))
                  call ispline(-psiichq(i),-psiiEF
     1            , ffprimEF, bffprim ,cffprim,dffprim,npsi,ffprimch(i))
               end do
               call chftransq(chcoeffp, presch, nchq, cftmq)
               do i= 1, nchq
                  psin=1.0d0-psiichq(i)
                  call chderiv(psin,nchq, chcoeffp, pprimch0tmp) !dpdpsi where psi is normalized flux
                  pprimch0(i)=pprimch0tmp
               end do

               call chftransq(chcoeffpp0, pprimch0, nchq, cftmq) !pprimch, checoeffpp will be multiplied by lamdba
               call chftransq(chcoeffff, ffprimch, nchq, cftmq)

            else if (iptable.eq.2) then
               isplch=1  !interpolate presure proflie to chebyshev points 
                         ! and then use it to find dpdpsi
               call spline(rhoEF, presEF     
     1              , bpres, cpres, dpres, npsi) 
            
               call spline(rhoEF, fpolEF
     1              , bfpol, cfpol, dfpol,npsi)

               do i= 1, nchq
                  call ispline(1.0d0-psiichq(i),rhoEF
     1                 , presEF, bpres ,cpres,dpres,npsi,preschrho(i)) !pressure at the chebyshev grids in rho
                  call ispline(1.0d0-psiichq(i),rhoEF
     1            , fpolEF, bfpol ,cfpol,dfpol,npsi,fpolchrho)  !F at the chebyshev grids in rho
                  fpol2chrho(i)=fpolchrho**2
               end do
               call chftransq(chcoeffprho, preschrho, nchq, cftmq)
               call chftransq(chcoefff2rho, fpol2chrho, nchq, cftmq)
               
               rhoch(1:nchq)=dsqrt(1.0d0-psiichq(1:nchq))  !initial guess for the relataion of rho and psi
               dpsiidrhoch(1:nchq)=2*rhoch(1:nchq)
            end if
         case (4)               !Genray, TORIC analytic form
!     dense(i)=(dense0(i)-denseb(i))*(1-rho**rn1de(i))**rn2de(i)+denseb(i)
!     tempe(i)=(te0(i)-teb(i))*(1-rho**rn1te(i))**rn2te(i)+teb(i)
!     rho=sqrt(1-psii)
            write(*,*) 'psiichq',psiichq(1:nchq)
            do i=1,nchq
               psinchq=1.0d0-psiichq(i)
               rho=dsqrt(psinchq)
               dens=(dense0-denseb)*(1.0d0-rho**rn1de)**rn2de
     1              +denseb ! in [10^19 m^-3]
               tempe=(te0-teb)*(1.0d0-rho**rn1te)**rn2te+teb  ! in [keV]
               tempi=(ti0-tib)*(1.0d0-rho**rn1ti)**rn2ti+tib 
               dndpsi=(dense0-denseb)*(1.0d0-rho**rn1de)**(rn2de-1)
     1              *(-rn1de*rn2de)*rho**(rn1de-2)/2.0d0
               dTdpsie=(te0-teb)*(1.0d0-rho**rn1te)**(rn2te-1)
     1              *(-rn1te*rn2te)*rho**(rn1te-2)/2.0d0
               dTdpsii=(ti0-tib)*(1.0d0-rho**rn1ti)**(rn2ti-1)
     1              *(-rn1ti*rn2ti)*rho**(rn1ti-2)/2.0d0

               const=1.60217657d+3 ! [10^19 keV/m^3 -> J/m^3]
               presch(i)=const*dens*(tempe+tempi)
               ! dp/dpsii where psii=0 at the core and psii=1 at the edge
               pprimch(i)=const*(dndpsi*(tempe+tempi)
     1              +dens*(dTdpsie+dTdpsii))      !assumed Te=Ti
               ne(i)=dens*0.1d0       ! [10^20 m^-3]
               Te(i)=tempe*1.0d+3 ! [eV]
               etai(i)=(dTdpsii/dndpsi)*(dens/tempi)
            end do
            fpolch(1:nchq)= dabs(F0)*fpolsgn !initial guess
            write(*,*) ' presch',presch(1:nchq)
            write(*,*) ' pprimch',pprimch(1:nchq)
         end select

         select case (ijtype)   ! Option for J par profiles
         case(0)                ! solovev solution 
          
            psicontmp0(1:nchq0)=psiichq0(1:nchq0)/lambdaex 
            call findq_solovev(nchq0,psicontmp0,fpolch0,qpsich0)
            psicontmp(1:nchq)=psiichq(1:nchq)/lambdaex 
            call findq_solovev(nchq,psicontmp,fpolch,qpsich)
           
         case(1)
            do i=1,nchq0
               psinchq=1.0d0-psiichq0(i)
               jpsich0(i)=jpar0*(1.0d0-psinchq**jpin)**jpout
            end do
            do i=1,nchq
               psinchq=1.0d0-psiichq(i)
               jpsich(i)=jpar0*(1.0d0-psinchq**jpin)**jpout
            end do
            write(*,*) 'jpsich',jpsich(1:nchq)
            write(*,*) 'jpsich0',jpsich0(1:nchq0)
         case(2) ! profile by Table
            call spline(-psiij, jpsiin     
     1           , bjpsi, cjpsi, djpsi, njpsi)        
            write(*,*) 'bjpsi',bjpsi(1:njpsi)
            do i= 1, nchq0
               call ispline(-psiichq0(i),-psiij
     1              , jpsiin , bjpsi, cjpsi, djpsi,njpsi, jpsich0(i))
            end do
            write(*,*) 'jpsich0',jpsich0(1:nchq0)
            call chftransq(chcoeff0,jpsich0,nchq0,cftmq0)
            do i=1,nchq
               call chfit(1-psiichq(i), nchq0, chcoeff0,jpsich(i))
            end do
            write(*,*) 'jpsich',jpsich(1:nchq)
            
         case(3) ! Ohmic + off-axis external current
            do i=1,nchq0
               psinchq=1.0d0-psiichq0(i)
               jpsich0(i)=jpar0*((1.0d0-psinchq**jpin)**jpout
     1              +Jomax*exp(-((psinchq-Joloc)/Jow)**2))
            end do
            do i=1,nchq
               psinchq=1.0d0-psiichq(i)
               jpsich(i)=jpar0*((1.0d0-psinchq**jpin)**jpout
     1             +Jomax*exp(-((psinchq-Joloc)/Jow)**2))
            end do
         case(4) ! Only Boostrap current, which will be evaluated during iterations   
            jpsich(1:nchq)=-mu0*Rmid**2*pprimch(1:nchq) !initial guess from pressure profile    
            write(*,*) 'Jpsich',jpsich(1:nchq)
         !  jpsich0 will be determined by interpolating jpsich during iteration

            call chsetupq(nchy, psiichy, cftmy, spybx, spyxb, spydef) !setup for pitchangle grid to evalute Jbs

         case(5) ! Only Ohmic current, which will be evaluated during iterations   
            if (iconvj.eq.1) then
               jpsich(1:nchq)=mu0*Rmid**2*lambda*pprimch(1:nchq) !initial guess from normalized pressure profile
            else
               jpsich(1:nchq)=-mu0*Rmid**2*pprimch(1:nchq) !initial guess from pressure profile    
            end if
            write(*,*) 'Jpsich',jpsich(1:nchq)
         !  jpsich0 will be determined by interpolating jpsich during iteration

            call chsetupq(nchy, psiichy, cftmy, spybx, spyxb, spydef) !setup for pitchangle grid to evalute Jbs
         end select

         ffprimpre1(1:nchq)=-jpsich(1:nchq) !initialization
         ffprimpre2(1:nchq)=Rmid**2.0d0
         write(*,*) 'pprim',pprimch(1:nchq)
         write(*,*) 'qpsich',qpsich(1:nchq)
         write(*,*) 'fpolch',fpolch(1:nchq)
         write(*,*) 'lambdaex',lambdaex
c         lambda=lambdaex

         if (itftype.ne.0) then  !If torflow is considered
            do i=1,nchq
               presch(i)=0.0d0
               do k=1,nchq
                  presch(i)=presch(i)+spxb((k-1)*nchq+i)*pprimch(k)
               end do
            end do
            presch(1:nchq)=-presch(1:nchq)/(2.0d0) !lambda is assumed to be -1
            write(*,*) 'presch',presch(1:nchq)
            call chftransq(chcoeffp, presch, nchq, cftmq)

            select case (itftype) ! option for the pressure by toroidal flow
            case(1)             !mach**2*pressure profile
          
               ptflowch(1:nchq)=presch(1:nchq)*mach**2/(2.0d0) 
               ptfprimch(1:nchq)=pprimch(1:nchq)*mach**2/(2.0d0) 
               ! presch is due to total pressure by electrons and ions
               ! but ptflow is only due to ions
c               write(*,*) 'ptflowch',ptflowch(1:nchq)
            case(2)             !ptf0*(1-psin**ptfin)**ptfout
               do i=1,nchq
                  psinchq=1.0d0-psiichq(i)
                  ptflowch(i)=ptf0*(1-psinchq**ptfin)**ptfout
     1                 /(mu0*Rmid**2)    !ptfflowch includes  (mu0*Rmid**2)
                  ptfprimch(i)=ptf0*ptfout*(1-psinchq**ptfin)
     1                 **(ptfout-1)*(-ptfin)*psinchq**(ptfin-1)
               end do
            case(22)             !ptf0*(1-psin**ptfin)**ptfout
               do i=1,nchq
                  psinchq=1.0d0-psiichq(i)
                  ptflowch(i)=ptfmax*exp(-((psinchq-ptfloc)/ptfw)**2)
     1                 /(mu0*Rmid**2)    !ptfflowch includes  (mu0*Rmid**2)
                  ptfprimch(i)=ptfmax*exp(-((psinchq-ptfloc)/ptfw)**2)
     1                 *(-2.0d0)*(psinchq-ptfloc)/ptfw**2
               end do
            case (3)
               call spline(-psiiEF, ptflowEF   
     1              , btflow, ctflow, dtflow, npsi) 
            
               do i= 1, nchq
                  call ispline(-psiichq(i),-psiiEF, ptflowEF 
     1                , btflow, ctflow, dtflow, npsi , ptflowch(i))
               end do

               call chftransq(chcoefftf, ptflowch, nchq, cftmq)
               do i= 1, nchq
                  psinchq=1.0d0-psiichq(i)
                  call chderiv(psinchq,nchq, chcoefftf,ptfprimch(i))
               end do

            end select
         end if
      end if
 

!-------------------------------------
! The following is for every iteration

      call chftransq(chcoeffpre1, ffprimpre1, nchq, cftmq)
      call chftransq(chcoeffpre2, ffprimpre2, nchq, cftmq)
c      call chftransq(chcoeffgoal, ffprimgoal, nchq, cftmq)
      write(*,*) 'ffprimpre1',ffprimpre1(1:nchq)
      write(*,*) 'ffprimpre2',ffprimpre2(1:nchq)

      if ((iptype.eq.2).and.(iptable.eq.1)) then ! make pprimch from pressure for iptable=1 or 2
         pprimch(1:nchq)=pprimch0(1:nchq)*lambda
         chcoeffpp(1:nchq)=chcoeffpp0(1:nchq)*lambda
         write(*,*) 'pprimch',pprimch(1:nchq)
      else if  ((iptype.eq.2).and.(iptable.eq.2)) then
         do i= 1, nchq
            psichrho=rhoch(i)
            call chderiv(psichrho,nchq,chcoeffprho,pprimchrho) !dpdpsi where psi is normalized flux
            pprimch(i)=pprimchrho/dpsiidrhoch(i)*lambda
         end do
         write(*,*) 'pprimch',pprimch(1:nchq)
         call chftransq(chcoeffpp, pprimch, nchq, cftmq)
      end if
         
c      write(*,*) 'ffprimgoal',ffprimgoal(1:nchq)
      if (itftype.eq.0) then     ! If toroflow is not considered
         do inext = 1, ntot
            psin=1.0d0-psii(inext)
            if (psin.lt.0.0d0) then
               psin=0.0d0
               write(*,*) 'negative psin !!'
            end if

            call chfit(psin, nchq, chcoeffpre1,ffprim_pre1)
            call chfit(psin, nchq, chcoeffpre2,ffprim_pre2)
c            call chfit(psin, nchq, chcoeffgoal,ffprim_next)
            ffprim_inext=ffprim_pre1
c            ffprim_inext=ffprim_next*wffprim+
c     1           (1.0d0-wffprim)*ffprim_pre1

            select case (iptype)
            case(0)             ! solovev solution      
               pprim_inext=-csol/mu0
            case(1)             ! profiles by parameters 
               pprim_inext=p0*(1-psin**pin)**pout !p0 includes (mu0*Rmid**2)
               pprim_inext=pprim_inext/(mu0*Rmid**2)
            
            case (2,3)          ! profile by table or EFIT
               if (isplch.eq.0) then
                  call ispline(-psii(inext),-psiiEF
     1                , pprimEF, bpprim ,cpprim,dpprim,npsi,pprim_inext)
               else if (isplch.eq.1) then
                  psin=1-psii(inext)
                  call chfit(psin,nchq,chcoeffpp,pprim_inext)
               end if
            case (4)            !Genray, TORIC analytic form
          
               rho=dsqrt(psin)
               if (rho.eq.0.0d0) then
                  pprim_inext=0.0d0
               else
                  dens=(dense0-denseb)*(1.0d0-rho**rn1de)**rn2de
     1                 +denseb  ! in [10^19 m^-3]
                  tempi=(ti0-tib)*(1.0d0-rho**rn1ti)**rn2ti+tib 
                  dndpsi=(dense0-denseb)*(1.0d0-rho**rn1de)**(rn2de-1)
     1                 *(-rn1de*rn2de)*rho**(rn1de-2)/2.0d0
                  dTdpsie=(te0-teb)*(1.0d0-rho**rn1te)**(rn2te-1)
     1                 *(-rn1te*rn2te)*rho**(rn1te-2)/2.0d0
                  dTdpsii=(ti0-tib)*(1.0d0-rho**rn1ti)**(rn2ti-1)
     1                 *(-rn1ti*rn2ti)*rho**(rn1ti-2)/2.0d0

                  const=1.60217657d+3 ! [10^19 keV/m^3 -> J/m^3]

               ! dp/dpsii where psii=0 at the core and psii=1 at the edge
                  pprim_inext=const*(dndpsi*(temp+tempi)
     1              +dens*(dTdpsie+dTdpsii))      !assumed Te=Ti
               end if
            end select
      
            select case (iconvj)
            case (0)   ! fixed dP/dPsi (actual flux, Psi)
               fsol(inext)=lambda*((-Rt3(inext)**2.0d0+ffprim_pre2)*mu0
     1              *pprim_inext 
     2              -ffprim_inext)
            case (1)  ! fixed dP/dpsi (normalized flux, psi)
               fsol(inext)=-lambdaex**2*(-Rt3(inext)**2.0d0+ffprim_pre2)
     1              *mu0*pprim_inext
     2          -ffprim_inext*lambda
c               write(*,*) 'pprim_inext',fsol(inext),pprim_inext
c               write(*,*) 'frim_inext',ffprim_pre2,ffprim_pre2
            case (2)
               fsol(inext)=lambda*(-Rt3(inext)**2.0d0*mu0
     1           *pprim_inext -ffprim_pre2)!*lambdaex/lambda)
     2          -ffprim_inext*lambda*(lambda/lambdaex)
            case (3)
               fsol(inext)=lambda*(-Rt3(inext)**2.0d0*mu0
     1           *pprim_inext -ffprim_pre2*lambdaex/lambda)
     2          -ffprim_inext*lambda
            end select
         end do

      else                      ! If toroflow is considered
        
         do inext = 1, ntot
            psin=1.0d0-psii(inext)
            call chfit(psin, nchq, chcoeffpre1,ffprim_pre1)
            call chfit(psin, nchq, chcoeffpre2,ffprim_pre2)
c            call chfit(psin, nchq, chcoeffgoal,ffprim_next)
            ffprim_inext=ffprim_pre1
c            ffprim_inext=ffprim_next*wffprim+
c     1           (1.0d0-wffprim)*ffprim_pre1

            select case (iptype)
            case(0)             ! solovev solution      
                  !will be added
            case(1)             ! profiles by parameters 
               select case (itftype) 
               case(1)          !mach**2*pressure profile
                  ptfbyp=mach**2/(2.0d0) 
                  pprim_inext=p0*(1-psin**pin)**pout !p0 includes (mu0*Rmid**2)

                  pexp=dexp(ptfbyp*((Rt3(inext)/Rmid)**2-1))
                  pprim_inext1=pexp*pprim_inext
               case(2)          !ptf0*(1-psin**ptfin)**ptfout
                  pprim_inext=p0*(1-psin**pin)**pout !p0 includes (mu0*Rmid**2)
                  ptflow =ptf0*(1-psin**ptfin)**ptfout !ptf0 includes (mu0*Rmid**2)
                  call chfit(psin,nchq,chcoeffp,p_inext)
                  p_inext=p_inext*mu0*Rmid**2  !chcoeffp does not include (mu0*Rmid**2)

                  ptfbyp = ptflow/p_inext
                  if (ptfbyp.gt.10) then
                     ptfbyp = 10.0
                  end if
                  pexp=dexp(ptfbyp*((Rt3(inext)/Rmid)**2-1))
c                  write(*,*) 'ptflow',ptflow,p_inext,ptfbyp,pexp

                  ptfprim =ptf0*ptfout*(1-psin**ptfin)
     1                 **(ptfout-1)*(-ptfin)*psin**(ptfin-1)
                  pprim_inext1=pexp*(pprim_inext+
     1                 ((Rt3(inext)/Rmid)**2-1)
     2                 *(ptfprim-ptfbyp*pprim_inext))

               case(22)         !ptfmax*exp(-((psin-ptfloc)/ptfw)**2))

                  pprim_inext=p0*(1-psin**pin)**pout !p0 includes (mu0*Rmid**2)
                  ptflow =ptfmax*exp(-((psin-ptfloc)/ptfw)**2) !ptfmax includes (mu0*Rmid**2)
                  call chfit(psin,nchq,chcoeffp,p_inext)
                  p_inext=p_inext*mu0*Rmid**2  !chcoeffp does not include (mu0*Rmid**2)

                  ptfbyp = ptflow/p_inext
                  if (ptfbyp.gt.10) then
                     ptfbyp = 10.0
                  end if
                  pexp=dexp(ptfbyp*((Rt3(inext)/Rmid)**2-1))
c                  write(*,*) 'ptflow',ptflow,p_inext,ptfbyp,pexp

                  ptfprim =ptfmax*exp(-((psin-ptfloc)/ptfw)**2)
     1                 *(-2.0d0)*(psin-ptfloc)/ptfw**2

                  pprim_inext1=pexp*(pprim_inext+
     1                 ((Rt3(inext)/Rmid)**2-1)
     2                 *(ptfprim-ptfbyp*pprim_inext))
               case (3)         !from table
                  pprim_inext=p0*(1-psin**pin)**pout !p0 includes (mu0*Rmid**2)
                  call ispline(-psii(inext),-psiiEF, ptflowEF 
     1                 , btflow, ctflow, dtflow, npsi , ptflow)
                  call chfit(psin,nchq,chcoeffp,p_inext)
                  p_inext=p_inext*mu0*Rmid**2 !chcoeffp does not include (mu0*Rmid**2)
                  ptfbyp = ptflow/p_inext
                  pexp=dexp(ptfbyp*((Rt3(inext)/Rmid)**2-1))

                  call chderiv(psin,nchq, chcoefftf,ptfprim)
                  pprim_inext1=pexp*(pprim_inext+
     1                 ((Rt3(inext)/Rmid)**2-1)
     2                 *(ptfprim-ptfbyp*pprim_inext))

               end select
               pprim_inext=pprim_inext/(mu0*Rmid**2)
               pprim_inext1=pprim_inext1/(mu0*Rmid**2)
            case (2,3)          ! profile by table or EFIT
                !will be added
            end select



            select case (iconvj)
            case (0)   ! fixed dP/dPsi (actual flux, Psi)
               fsol(inext)=lambda*((-Rt3(inext)**2.0d0
     1             *pprim_inext1*mu0 +ffprim_pre2*pprim_inext*mu0)
     2              -ffprim_inext)
            case (1)  ! fixed dP/dpsi (normalized flux, psi)
               fsol(inext)=-lambdaex**2*(-Rt3(inext)**2.0d0
     1             *pprim_inext1*mu0 +ffprim_pre2*pprim_inext*mu0)
     2              -ffprim_inext*lambda
            end select
c           write(*,*) 'pprim_inext',fsol(inext),pprim_inext,pprim_inext1
c               write(*,*) 'frim_inext',ffprim_pre2,ffprim_inext
c            case (2)
c               fsol(inext)=lambda*(-Rt3(inext)**2.0d0*mu0
c     1           *pprim_inext -ffprim_pre2)!*lambdaex/lambda)
c     2          -ffprim_inext*lambda*(lambda/lambdaex)
c            case (3)
c               fsol(inext)=lambda*(-Rt3(inext)**2.0d0*mu0
c     1           *pprim_inext -ffprim_pre2*lambdaex/lambda)
c     2          -ffprim_inext*lambda
        
         end do
      end if

      write(*,*) 'iterataion index =', iiter
      do inext = 1, ntot
    
c         if (fsol(inext).gt.0.0d0) then
c            fsol(inext)=0.0d0
c         end if
 
         f(inext)=(fsol(inext)/(Rt3(inext)**(0.5d0))
     1        +3.0d0/(4.0d0*(Rt3(inext)**(2.5d0)))*psii(inext))
     2        *(dzdw3(inext)*conjg(dzdw3(inext)))
     3       -halpha**2.0d0*psii(inext)/dsqrt(Rt3(inext))
      end do


      ibc = 1
      call elliptic(nsub,korder,ntheta,rnd,rnd(nr+1),
     1     u,ur,uth,urr,urt,utt,f,blength,bnodes,ibc,g)

      psif(1:ntot) =u(1:ntot)*sqrt(Rt3(1:ntot))
     
      maxpsi = 0.0d0
      do inext = 1, ntot
         if (abs(psif(inext)) > maxpsi) then
            maxpsi = abs(psif(inext))
            inext_maxpsi=inext
         end if
      end do
 
      maxdpsi = 0.0d0
      do inext = 1, ntot
         if (abs(psif(inext)-psii(inext)) > maxdpsi) then
            maxdpsi = abs(psif(inext)-psii(inext))
         end if
      end do

      ! find magnetic axis off the grid if needed
      if ((maxdpsi.lt.dpsistartmag).and.(ifindmagaxis.eq.1)) then
       
         if ((rndm.lt.0.0d0).or.(iremap.eq.1)) then !For initial value, find values on grid 
            call findinitmax(inext_maxpsi, rndm, tndm
     1           , dpsidrm, dpsidtm, dpsidrrm ,dpsidrtm, dpsidttm
     2           , Rt3m, Zt3m, psim)
         end if
         call findmagaxis(inext_maxpsi, rndm, tndm, dpsidrm, dpsidtm 
     1        ,dpsidrrm ,dpsidrtm, dpsidttm, Rt3m, Zt3m, psim
     2        ,dpsim2)

         maxpsiongird=maxpsi/psim
         maxpsi=psim
         Rmaxis=Rt3m
         Zmaxis=Zt3m

         dist2=(Rmaxis-Rmap0)**2+(Zmaxis-Zmap0)**2
       else  ! magnetic axis on the grid
         Rmaxpsi=Rt3(inext_maxpsi)
         Zmaxpsi=Zt3(inext_maxpsi)
         Rmaxis=Rmaxpsi
         Zmaxis=Zmaxpsi

         dist2=(Rmaxpsi-Rmap0)**2+(Zmaxpsi-Zmap0)**2
       
      end if
c      Rmid=Rmaxis
      psif(1:ntot) = psif(1:ntot)/maxpsi
c      if (iconvq.ne.0) then
c         if ((iconvq.ne.1).or.(iiter.gt.relaxiter)) then 
            lambda = lambda/maxpsi
c         end if
c      end if

      psii(1:ntot) = psif(1:ntot)

      write(*,*) "(Rmap0,Zmap0): disk center=",Rmap0,Zmap0
      write(*,*) "(Rmaxpsi,Zmaxpsi): maxpsi on grid =",Rmaxpsi,Zmaxpsi
      write(*,*) "(Rmaxis,Zmaxis): maxpsi off grid =",Rmaxis,Zmaxis
      write(*,*) "dist2,epsmaxdist=",dist2,epsmaxdist
      diffpsi=1.0d0-psii(1)
      diffpsich=(1.0d0-dcos(pi/2.0d0/(nchq+1)))/2.0d0
      write(*,*) "distpsi,diffpsich=",diffpsi,diffpsich
c      if ((dist2.lt.epsmaxdist).and.
c     1     (dabs(diffpsi).lt.0.9*dabs(diffpsich))) then
         if (dabs(diffpsi).lt.0.8*dabs(diffpsich)) then
c         if (iremap.eq.1) then
            call checkcontour(nchq,psiichq,nr,ntheta,rnd,tnd
     1          ,psii,icheck)
         write(*,*) 'icheck',icheck
         if (icheck.ne.1) then
            iremap = 1
            Rmap0 = Rmaxis      !-0.8*diffpsich
            Zmap0 = Zmaxis 
         else 
c            epsmaxdist=epsmaxdist*0.1d0
c            iremap = 1
c            Rmap0 = Rmaxis      !-0.8*diffpsich
c            Zmap0 = Zmaxis      !maxpsi
            
c         end if

         if ((iremap.eq.0).and.(nchq.lt.nchq0)) then
            nchq=nchq+1
            call chsetupq(nchq, psiichq, cftmq, spbx, spxb, spdef)
            write(*,*) 'nchq=',nchq
            select case (iptype) ! Option for pressure and current profiles
            case(0)             ! solovev solution 
               pprimch(1:nchq)=-csol/mu0
               fpolch(1:nchq)= dabs(F0)*fpolsgn
               lambdaex=1.0d0/(d1-d2**2.0d0/(4.0d0*(csol/8.0d0+d3)))
               write(*,*) 'lambdaex=',lambdaex
               iqtype = 0
               
            case(1)
               do i=1,nchq
                  psinchq=1.0d0-psiichq(i)
c          pprimch(i)=p0*psinchq*(1.0d0-psinchq**pin)**pout/(mu0*Rmid**2) !test
                  pprimch(i)=p0*(1.0d0-psinchq**pin)**pout/(mu0*Rmid**2) !test
c     fpolch(i)=F0*(1.0d0+Ffac*psinchq**Fpin)**Fpout
               end do
               fpolch(1:nchq)= dabs(F0)*fpolsgn !initial guess
            case(2,3)           ! profile by Table or EFIT
               call spline(-psiiEF, pprimEF     
     1              , bpprim, cpprim, dpprim, npsi) 
            
               call spline(-psiiEF, fpolEF
     1              , bffprim, cffprim, dffprim,npsi)   
            
               do i= 1, nchq
                  call ispline(-psiichq(i),-psiiEF
     1                 , pprimEF, bpprim ,cpprim,dpprim,npsi,pprimch(i))
               end do

               call ispline(-1.0d0,-psiiEF
     1              , pprimEF, bpprim ,cpprim,dpprim,npsi,p0)
               call ispline(-1.0d0,-psiiEF
     1              , ffprimEF, bffprim ,cffprim,dffprim,npsi,FF0)

            case (4)            !Genray, TORIC analytic form
               do i=1,nchq
                  psinchq=1.0d0-psiichq(i)
                  rho=dsqrt(psinchq)
                  dens=(dense0-denseb)*(1.0d0-rho**rn1de)**rn2de
     1                 +denseb  ! in [10^19 m^-3]
                  tempe=(te0-teb)*(1.0d0-rho**rn1te)**rn2te+teb ! in [keV]
                  tempi=(ti0-tib)*(1.0d0-rho**rn1ti)**rn2ti+tib 
                  dndpsi=(dense0-denseb)*(1.0d0-rho**rn1de)**(rn2de-1)
     1                 *(-rn1de*rn2de)*rho**(rn1de-2)/2.0d0
                  dTdpsie=(te0-teb)*(1.0d0-rho**rn1te)**(rn2te-1)
     1                 *(-rn1te*rn2te)*rho**(rn1te-2)/2.0d0
                  dTdpsii=(ti0-tib)*(1.0d0-rho**rn1ti)**(rn2ti-1)
     1                 *(-rn1ti*rn2ti)*rho**(rn1ti-2)/2.0d0
                  
                  const=1.60217657d+3 ! [10^19 keV/m^3 -> J/m^3]

                  ! dp/dpsii where psii=0 at the core and psii=1 at the edge
                  pprimch(i)=const*(dndpsi*(temp+tempi)
     1                 +dens*(dTdpsie+dTdpsii)) !assumed Te=Ti
                  ne(i)=dens*0.1d0 ! [10^20 m^-3]
                  Te(i)=tempe*1.0d+3 ! [eV]
                  etai(i)=(dTdpsii/dndpsi)*(dens/tempi)

                  !nue_star(i)
                  !nui_star(i)
c     fpolch(i)=F0*(1.0d0+Ffac*psinchq**Fpin)**Fpout
               end do
               write(*,*) 'pprimch',pprimch(1:nchq)
               fpolch(1:nchq)= dabs(F0)*fpolsgn !initial guess
   
            end select

            select case (ijtype) ! Option for J par profiles
            case(0)             ! solovev solution 
               
               psicontmp(1:nchq)=psiichq(1:nchq)/lambdaex 
               call findq_solovev(nchq,psicontmp,fpolch,qpsich)
           
            case(1)
              
               do i=1,nchq
                  psinchq=1.0d0-psiichq(i)
                  jpsich(i)=jpar0*(1.0d0-psinchq**jpin)**jpout
               end do
               write(*,*) 'jpsich',jpsich(1:nchq)
            case(2)             ! profile by Table
               
               call chftransq(chcoeff0,jpsich0,nchq0,cftmq0)
               do i=1,nchq
                  call chfit(1-psiichq(i), nchq0, chcoeff0,jpsich(i))
               end do
               write(*,*) 'jpsich',jpsich(1:nchq)
            case(3)             ! Ohmic + off-axis external current
            
               do i=1,nchq
                  psinchq=1.0d0-psiichq(i)
                  jpsich(i)=jpar0*((1.0d0-psinchq**jpin)**jpout
     1                 +Jomax*exp(-((psinchq-Joloc)/Jow)**2))
               end do
               write(*,*) 'jpsich',jpsich(1:nchq)
            case(4)             ! Only Boostrap current, which will be evaluated during iterations   
c               jpsich(1:nchq)=-mu0*Rmid**2*pprimch(1:nchq) !initial guess from pressure profile    
c               write(*,*) 'Jpsich',jpsich(1:nchq)
!  jpsich0 will be determined by interpolating jpsich during iteration

            case (5)

            end select
      
            write(*,*) 'psiichq', psiichq(1:nchq)

c         epsmaxdist = 0.0001d0
            write(*,*) 'fpolch', fpolch(1:nchq)
            write(*,*) 'pprimch', pprimch(1:nchq)
     
            if (itftype.ne.0) then
               do i=1,nchq
                  presch(i)=0.0d0
                  do k=1,nchq
                     presch(i)=presch(i)+spxb((k-1)*nchq+i)*pprimch(k)
                  end do
               end do
               presch(1:nchq)=-presch(1:nchq)/(2.0d0) !lambda is assumed to be -1
               write(*,*) 'presch',presch(1:nchq)
               call chftransq(chcoeffp, presch, nchq, cftmq)
               
               select case (itftype) ! option for the pressure by toroidal flow
               case(1)          !mach**2*pressure profile
                  
                  ptflowch(1:nchq)=presch(1:nchq)*mach**2/(2.0d0) 
                  ptfprimch(1:nchq)=pprimch(1:nchq)*mach**2/(2.0d0) 
! presch is due to total pressure by electrons and ions
! but ptflow is only due to ions
c     write(*,*) 'ptflowch',ptflowch(1:nchq)
               case(2)          !ptf0*(1-psin**ptfin)**ptfout
                  do i=1,nchq
                     psinchq=1.0d0-psiichq(i)
                     ptflowch(i)=ptf0*(1-psinchq**ptfin)**ptfout
     1                    /(mu0*Rmid**2) !ptfflowch includes  (mu0*Rmid**2)
                     ptfprimch(i)=ptf0*ptfout*(1-psinchq**ptfin)
     1                    **(ptfout-1)*(-ptfin)*psinchq**(ptfin-1)
                  end do
               case(22)         !ptf0*(1-psin**ptfin)**ptfout
                  do i=1,nchq
                     psinchq=1.0d0-psiichq(i)
                     ptflowch(i)=ptfmax*exp(-((psinchq-ptfloc)/ptfw)**2)
     1                    /(mu0*Rmid**2) !ptfflowch includes  (mu0*Rmid**2)
                    ptfprimch(i)=ptfmax*exp(-((psinchq-ptfloc)/ptfw)**2)
     1                    *(-2.0d0)*(psinchq-ptfloc)/ptfw**2
                  end do
               case (3)
                  call spline(-psiiEF, ptflowEF   
     1                 , btflow, ctflow, dtflow, npsi) 
                  
                  do i= 1, nchq
                     call ispline(-psiichq(i),-psiiEF, ptflowEF 
     1                    , btflow, ctflow, dtflow, npsi , ptflowch(i))
                  end do
                  
                  call chftransq(chcoefftf, ptflowch, nchq, cftmq)
                  do i= 1, nchq
                     psinchq=1.0d0-psiichq(i)
                     call chderiv(psinchq,nchq, chcoefftf,ptfprimch(i))
                  end do
                  
               end select
            end if

         end if 
         iremap = 0
         end if
         
      else
         iremap = 1
         Rmap0 = Rmaxis !+0.5d0*(Rmap0-Rmaxis) !0.8*diffpsich !Rmaxpsi
         if (isymud.eq.1) then
            Zmap0=0.0d0
         else
            Zmap0 = Zmaxis      !Zmaxpsi
         end if
      end if


      do j=1, ntheta
         do i=1, nr
            inext = i + (j-1)*nr
      
            dzdwr=dzdw3(inext)
            dzdwi=(-ima)*dzdw3(inext)
            dzdwwr=dzdww3(inext)
            dzdwwi=(-ima)*dzdww3(inext)
            dRdr(inext)=dzdwr*dcos(tnd(j))-dzdwi*dsin(tnd(j))
            dRdrr(inext)=dzdwwr*(dcos(tnd(j)))**2
     1           -dzdwwr*(dsin(tnd(j)))**2
     1           -dzdwwi*dsin(tnd(j))*dcos(tnd(j))
         end do
      end do
       
      u(1:ntot)=u(1:ntot)/maxpsi
      ur(1:ntot)=ur(1:ntot)/maxpsi
      urr(1:ntot)=urr(1:ntot)/maxpsi
      uth(1:ntot)=uth(1:ntot)/maxpsi
      urt(1:ntot)=urt(1:ntot)/maxpsi
      utt(1:ntot)=utt(1:ntot)/maxpsi
      
      dpsidr(1:ntot)=ur(1:ntot)*sqrt(Rt3(1:ntot))
     1     +u(1:ntot)/(2*sqrt(Rt3(1:ntot)))*dRdr(1:ntot)

      dpsidrr(1:ntot) = urr(1:ntot)*sqrt(Rt3(1:ntot))
     1     +ur(1:ntot)*dRdr(1:ntot)/sqrt(Rt3(1:ntot))
     1     +u(1:ntot)*dRdrr(1:ntot)/(2*sqrt(Rt3(1:ntot)))

      write(*,*) '||psif(inext)|| ',maxpsi
      write(*,*) '||psif(inext)-psii(inext)|| ',maxdpsi
      write(*,*) 'lambda', lambda
      write(*,*) 'iremap ', iremap
      write(fiter,*) 'iter: lambda,maxdpsi,psi0',iiter,lambda,
     1        maxdpsi,(1.0d0/lambda)+psiB, iremap

      return
      end 
 
