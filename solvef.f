!--------------------------------------------------------
! solvef.f: This subroutine is to solve non-linear Grad-Shafranov equation 
! by iterating the Poisson solver in an unit disk using FF' input
!--------------------------------------------------------

      subroutine solvef

      use arrays, only:  iptype, file_bc, csol
      use arrays, only:  Rt3, Zt3, rnd, tnd, dzdw3, dzdww3
      use arrays, only:  nt2, eps, epsiter, maxiter
      use arrays, only:  nsub, kcheb, kLag, nr, ntot, blength, bnodes
      use arrays, only:  Rmid, Rmaxis, Zmaxis, cftmsub, Rmap0, Zmap0

      use arrays, only:  csol,reps,rkappa,delta,d1,d2,d3
      use arrays, only:  npsi, psi0, psiB, psiEF, psiiEF 
      use arrays, only:  rhoEF, presEF, chcoeffpp, chcoeffpp0, chcoeffff 
      use arrays, only:  pprimEF,ffprimEF,chcoeffprho,chcoefff2rho
      use arrays, only:  fpolEF, qpsiEF, bpprim, cpprim, dpprim
      use arrays, only:  bffprim, cffprim, dffprim, fpolch, rhoch
      use arrays, only:  bpres,cpres,dpres,bfpol,cfpol,dfpol

      use arrays, only:  nchq, psiichq, cftmq, spbx, spxb, spdef
      use arrays, only:  qpsich, pprimch, pprimch0,ffprimch, qpsich1
      use arrays, only:  wffprim, fpolcon2, ffprimpre, ffprimpre1
      use arrays, only:  ffprimpre2, ffprimgoal, dpsiidrhoch

      use arrays, only:  iiter, iremap, iremapon, iiterl, lambdaex
      use arrays, only:  g, psii, f, lambda, maxpsi, maxdpsi, lambdaf
      use arrays, only:  epsmaxdist, Rmaxpsi, Zmaxpsi, fpolsgn, psiex 
      use arrays, only:  dpsidr, dpsidrr, u, ur, uth, urr, urt, utt
      use arrays, only: q0, F0, p0, pin, pout, FF0, FFin, FFout
      use arrays, only: iftype,itftype, mach, ptf0, ptfin, ptfout
      use arrays, only: fiter, ftime, fout, isplch, iptable
      use arrays, only: chcoeffp, presch, chcoefftf, ptfprimch, ptflowch
      use arrays, only: ptflowEF, btflow, ctflow, dtflow, iscale
      use arrays, only: ptfmax, ptfloc, ptfw
      use arrays, only: rndm, tndm, dpsidrm, dpsidtm, dpsidrrm ,dpsidrtm
      use arrays, only: dpsidttm, Rt3m, Zt3m, psim
      use arrays, only: isolver, halpha, ifindmagaxis, dpsistartmag

      implicit double precision(a-h, o-z)
      real *8 fsol(ntot)
      real *8 psif(ntot) 
      real *8 dRdr(ntot), dRdrr(ntot)
      real *8 chcoeff(nchq),chcoeffpre1(nchq),chcoeffpre2(nchq)
      real *8 chcoeffgoal(nchq)
      real *8 preschrho(nchq),fpol2chrho(nchq)
      real *8 psicontmp(nchq), psiichsub(kcheb)
      real *8 fpolcon(nchq)
      real *8 mu0
      complex *16 ima
      

      write(*,*) 'start poisson, nt2',nt2

      write(*,*) 'start poisson, korder,nr,ntot',kcheb,nr,ntot
      pi = 4.0d0*datan(1.0d0)
      mu0=4.0d0*pi*1.0d-7
      ima = (0.0d0,1.0d0)

      ntheta=nt2

      korder = kcheb
      !--------------------------------------------------
      if (iiter.eq.1) then !initial condition for the first iteration
         
         do i=1,ntheta
            g(i)=0.0d0     !psii at the fixed Boundary for u 
         end do
         do j = 1,ntheta
            do i = 1,nr
               inext = i + (j-1)*nr
               psii(inext)=1.0d0-rnd(i)**2 !intial guess for psii
            end do
         end do

         lambda = lambdaex !-10.0d0   ! Initial guess for lambda. Sign should be correct.
         call chsetupsub(kcheb, psiichsub, cftmsub)  
         call chsetupq(nchq, psiichq, cftmq, spbx, spxb, spdef)

         select case (iptype)   ! Option for pressure profiles
         case(0)                ! solovev solution 
            pprimch(1:nchq)=-csol/mu0
         case(1)
            do i=1,nchq
               psinchq=1.0d0-psiichq(i)
               pprimch(i)=p0*(1-psinchq**pin)**pout/(mu0*Rmid**2)
            end do
         case(11)
            do i=1,nchq
               psinchq=1.0d0-psiichq(i)
               pprimch(i)=-p0*pin*pout*(1-psinchq**pin)**(pout-1.0d0)
     1              *psinchq**(pin-1.0d0)/(mu0*Rmid**2)
            end do
         case(2,3) ! profile by Table or EFIT
            if (iptable.eq.0) then
               call spline(-psiiEF, pprimEF     
     1              , bpprim, cpprim, dpprim, npsi) 
               do i= 1, nchq
                  call ispline(-psiichq(i),-psiiEF
     1                 , pprimEF, bpprim ,cpprim,dpprim,npsi,pprimch(i))
               end do

               call chftransq(chcoeffpp, pprimch, nchq, cftmq)
               
               call ispline(-1.0d0,-psiiEF
     1              , pprimEF, bpprim ,cpprim,dpprim,npsi,p0)
            else if (iptable.eq.1) then
               isplch=1  !interpolate presure proflie to chebyshev points 
                         ! and then use it to find dpdpsi
               call spline(-psiiEF, presEF     
     1              , bpres, cpres, dpres, npsi) 

               do i= 1, nchq
                  call ispline(-psiichq(i),-psiiEF
     1                 , presEF, bpres ,cpres,dpres,npsi,presch(i))
               end do
               call chftransq(chcoeffp, presch, nchq, cftmq)
               do i= 1, nchq
                  psin=1.0d0-psiichq(i)
                  call chderiv(psin,nchq, chcoeffp, pprimch0tmp) !dpdpsi where psi is normalized flux
                  pprimch0(i)=pprimch0tmp
               end do

               call chftransq(chcoeffpp0, pprimch0, nchq, cftmq) !pprimch, checoeffpp will be multiplied by lamdba

            else if (iptable.eq.2) then
               isplch=1  !interpolate presure proflie to chebyshev points 
                         ! and then use it to find dpdpsi
               call spline(rhoEF, presEF     
     1              , bpres, cpres, dpres, npsi) 

               rhoch(1:nchq)=dsqrt(1.0d0-psiichq(1:nchq))  !initial guess for the relataion of rho and psi
               dpsiidrhoch(1:nchq)=2*rhoch(1:nchq)
            end if
         end select

         select case (iftype)   ! Option for current profiles
         case(0)                ! solovev solution 
            ffprimch(1:nchq) = 0.0d0
         case(1)
            do i=1,nchq
               psinchq=1.0d0-psiichq(i)
               ffprimch(i)=FF0*(1-psinchq**ffin)**ffout
            end do
         case(11)
            do i=1,nchq
               ffprimch(i)=-FF0*ffin*ffout*(1-psinchq**ffin)
     1              **(ffout-1.0d0)*psinchq**(ffin-1.0d0)
            end do
         case(2,3) ! profile by Table or EFIT
            if (iptable.eq.0) then
               call spline(-psiiEF, ffprimEF
     1              , bffprim, cffprim, dffprim,npsi)   
            
               do i= 1, nchq
                  call ispline(-psiichq(i),-psiiEF
     1            , ffprimEF, bffprim ,cffprim,dffprim,npsi,ffprimch(i))
               end do

               call chftransq(chcoeffff, ffprimch, nchq, cftmq)
               call ispline(-1.0d0,-psiiEF
     1              , ffprimEF, bffprim ,cffprim,dffprim,npsi,FF0)
            else if (iptable.eq.1) then
               isplch=1  !interpolate presure proflie to chebyshev points 
                         ! and then use it to find dpdpsi
               call spline(-psiiEF, ffprimEF
     1              , bffprim, cffprim, dffprim,npsi)

               do i= 1, nchq
                  call ispline(-psiichq(i),-psiiEF
     1            , ffprimEF, bffprim ,cffprim,dffprim,npsi,ffprimch(i))
               end do
               call chftransq(chcoeffff, ffprimch, nchq, cftmq)

            else if (iptable.eq.2) then
               isplch=1  !interpolate presure proflie to chebyshev points 
                         ! and then use it to find dpdpsi
               call spline(rhoEF, fpolEF
     1              , bfpol, cfpol, dfpol,npsi)
            end if
         end select
         

         if ((itftype.ne.0).and.(iptable.eq.0)) then !If torflow is considered
            do i=1,nchq
               presch(i)=0.0d0
               do k=1,nchq
                  presch(i)=presch(i)+spxb((k-1)*nchq+i)*pprimch(k)
               end do
            end do
            presch(1:nchq)=-presch(1:nchq)/(2.0d0) !lambda is assumed to be -1 just for ratio between presch and ptflowch
            write(*,*) 'presch',presch(1:nchq)
            call chftransq(chcoeffp, presch, nchq, cftmq)

            select case (itftype) ! option for the pressure by toroidal flow
            case(1)             !mach**2*pressure profile
               ptflowch(1:nchq)=presch(1:nchq)*mach**2/(2.0d0) 
               ptfprimch(1:nchq)=pprimch(1:nchq)*mach**2/(2.0d0) 
               ! presch is due to total pressure by electrons and ions
               ! but ptflow is only due to ions
            case(2)             !ptf0*(1-psin**ptfin)**ptfout
               do i=1,nchq
                  psinchq=1.0d0-psiichq(i)
                  ptflowch(i)=ptf0*(1-psinchq**ptfin)**ptfout
     1                 /(mu0*Rmid**2)    !ptfflowch includes  (mu0*Rmid**2)
                  ptfprimch(i)=ptf0*ptfout*(1-psinchq**ptfin)
     1                 **(ptfout-1)*(-ptfin)*psinchq**(ptfin-1) 
                  !lambda is assumed to be -1 just for ratio between ptflowch and ptfprimch
               end do
            case(22)             !ptf0*(1-psin**ptfin)**ptfout
               do i=1,nchq
                  psinchq=1.0d0-psiichq(i)
                  ptflowch(i)=ptfmax*exp(-((psinchq-ptfloc)/ptfw)**2)
     1                 /(mu0*Rmid**2)    !ptfflowch includes  (mu0*Rmid**2)
                  ptfprimch(i)=ptfmax*exp(-((psinchq-ptfloc)/ptfw)**2)
     1                 *(-2.0d0)*(psinchq-ptfloc)/ptfw**2  !lambda is assumed to be -1
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
      endif
      !---------------------------------------------
      ! The following is for every iteration
      if (iptable.eq.1) then    ! make pprimch from pressure for iptable=1 
         pprimch(1:nchq)=pprimch0(1:nchq)*lambda
         chcoeffpp(1:nchq)=chcoeffpp0(1:nchq)*lambda
      else if (iptable.eq.2) then ! make pprimch from pressure for iptable=2

         call spline(rhoEF
     1           , fpolEF, bfpol ,cfpol,dfpol,npsi) ! JHSONG
         do i= 1, nchq
            rhochpsi=rhoch(i)
            call ispline(rhochpsi,rhoEF
     1           , presEF, bpres ,cpres,dpres,npsi,preschrho(i)) !pressure at rho corresponding to the chebyshev grids of psi
            call ispline(rhochpsi,rhoEF
     1           , fpolEF, bfpol ,cfpol,dfpol,npsi,fpolchrho)
            fpol2chrho(i)=fpolchrho**2
         end do
         call chftransq(chcoeffprho, preschrho, nchq, cftmq)
         call chftransq(chcoefff2rho, fpol2chrho, nchq, cftmq)
         do i =1, nchq
            psinchq=1.0d0-psiichq(i)
            call chderiv(psinchq,nchq,chcoeffprho,pprimchrho) !dpdpsi where psi is normalized flux
            call chderiv(psinchq,nchq,chcoefff2rho,f2primchrho)
            pprimch(i)=pprimchrho*lambda
            ffprimch(i)=f2primchrho/2.0d0*lambda
         end do
         call chftransq(chcoeffpp, pprimch, nchq, cftmq)
         call chftransq(chcoeffff, ffprimch, nchq, cftmq)
      end if

      if ((itftype.ne.0).and.(iptable.ne.0)) then !If pprimch changes, update presch
         do i=1,nchq
            presch(i)=0.0d0
            do k=1,nchq
               presch(i)=presch(i)+spxb((k-1)*nchq+i)*pprimch(k)
            end do
         end do
         presch(1:nchq)=-presch(1:nchq)/(2.0d0) !lambda is assumed to be -1 just for ratio between presch and ptflowch
         write(*,*) 'presch',presch(1:nchq)
         call chftransq(chcoeffp, presch, nchq, cftmq)
         
         select case (itftype)  ! option for the pressure by toroidal flow
         case(1)                !mach**2*pressure profile
            ptflowch(1:nchq)=presch(1:nchq)*mach**2/(2.0d0) 
            ptfprimch(1:nchq)=pprimch(1:nchq)*mach**2/(2.0d0) 
            ! presch is due to total pressure by electrons and ions
            ! but ptflow is only due to ions
         case(2)                !ptf0*(1-psin**ptfin)**ptfout
            do i=1,nchq
               psinchq=1.0d0-psiichq(i)
               ptflowch(i)=ptf0*(1-psinchq**ptfin)**ptfout
     1              /(mu0*Rmid**2) !ptfflowch includes  (mu0*Rmid**2)
               ptfprimch(i)=ptf0*ptfout*(1-psinchq**ptfin)
     1              **(ptfout-1)*(-ptfin)*psinchq**(ptfin-1) 
               !lambda is assumed to be -1 just for ratio between ptflowch and ptfprimch
            end do
         case(22)               !ptf0*(1-psin**ptfin)**ptfout
            do i=1,nchq
               psinchq=1.0d0-psiichq(i)
               ptflowch(i)=ptfmax*exp(-((psinchq-ptfloc)/ptfw)**2)
     1              /(mu0*Rmid**2) !ptfflowch includes  (mu0*Rmid**2)
               ptfprimch(i)=ptfmax*exp(-((psinchq-ptfloc)/ptfw)**2)
     1              *(-2.0d0)*(psinchq-ptfloc)/ptfw**2 !lambda is assumed to be -1
            end do
         case (3)
            call spline(-psiiEF, ptflowEF   
     1           , btflow, ctflow, dtflow, npsi) 
            
            do i= 1, nchq
               call ispline(-psiichq(i),-psiiEF, ptflowEF 
     1              , btflow, ctflow, dtflow, npsi , ptflowch(i))
            end do
            
            call chftransq(chcoefftf, ptflowch, nchq, cftmq)
            do i= 1, nchq
               psinchq=1.0d0-psiichq(i)
               call chderiv(psinchq,nchq, chcoefftf,ptfprimch(i))
            end do

         end select
      end if

      ! evaluate pprim, ffprim
      if (itftype.eq.0) then    ! If toroflow is not considered

         do inext = 1, ntot
            psin=1-psii(inext)
	    if (psin.lt.0.0d0) then
		psin=0.0d0
	    else if (psin.gt.1.0d0) then
		psin=1.0d0
	    end if

            select case (iptype)
            case(0)             ! solovev solution      
               pprim_inext=-csol/mu0
            case(1)             ! profiles by parameters 
               pprim_inext=p0*(1-psin**pin)**pout/(mu0*Rmid**2)
               !p0 includes (mu0*Rmid**2)
            case(11)            ! profiles by parameters   
               pprim_inext=-p0*pin*pout*(1-psin**pin)**(pout-1.0d0)
     1              *psin**(pin-1.0d0)/(mu0*Rmid**2)
            case (2,3)          ! profile by table or EFIT
               if (isplch.eq.0) then
                  call ispline(-psii(inext),-psiiEF
     1                , pprimEF, bpprim ,cpprim,dpprim,npsi,pprim_inext)
                        !pprimEF does not includes 4*pi*1.0d-7 
               else if (isplch.eq.1) then
                  call chfit(psin,nchq,chcoeffpp,pprim_inext)
               end if
            end select

            select case (iftype)
            case(0)             ! solovev solution      
               ffprim_inext=0.0d0
            case(1)             ! profiles by parameters 
               ffprim_inext=FF0*(1-psin**ffin)**ffout 
            case(11)            ! profiles by parameters   
               ffprim_inext=-FF0*ffin*ffout*(1-psin**ffin)
     1              **(ffout-1.0d0)*psin**(ffin-1.0d0)
            case (2,3)          ! profile by table or EFIT
               if (isplch.eq.0) then
                  call ispline(-psii(inext),-psiiEF, ffprimEF, bffprim 
     1                 ,cffprim, dffprim,npsi,ffprim_inext)
               else if (isplch.eq.1) then
                  call chfit(psin,nchq,chcoeffff,ffprim_inext)
               end if
            end select

            fsol(inext)=lambda*(-Rt3(inext)**2.0d0*mu0
     1           *pprim_inext-ffprim_inext)
         end do
    
      else                      ! If toroflow is considered
         
         do inext = 1, ntot
            psin=1-psii(inext)
            if (psin.lt.0.0d0) then
                psin=0.0d0
            else if (psin.gt.1.0d0) then
                psin=1.0d0
            end if

            select case (iptype)
            case(0)             ! solovev solution      
               pprim_inext=-csol/mu0
               p_inext=-csol/mu0*(psin-1)
           
            case(1)             ! profiles by parameters 
               pprim_inext=p0*(1-psin**pin)**pout/(mu0*Rmid**2) !p0 includes (mu0*Rmid**2)
               call chfit(psin,nchq,chcoeffp,p_inext)
            case (2,3)          ! profile by table or EFIT          
               if (isplch.eq.0) then
                  call ispline(-psii(inext),-psiiEF
     1                , pprimEF, bpprim ,cpprim,dpprim,npsi,pprim_inext)
                        !pprimEF does not includes 4*pi*1.0d-7 
               else if (isplch.eq.1) then
                  call chfit(psin,nchq,chcoeffpp,pprim_inext)
               end if
               call chfit(psin,nchq,chcoeffp,p_inext)
            end select

            select case (itftype) 
            case(1)             !mach**2*pressure profile
               ptfbyp=mach**2/(2.0d0)          
               pexp=dexp(ptfbyp*((Rt3(inext)/Rmid)**2-1))
               pprim_inext1=pexp*pprim_inext
            case(2)             !ptf0*(1-psin**ptfin)**ptfout
               ptflow =ptf0*(1-psin**ptfin)**ptfout/(mu0*Rmid**2) !ptf0 includes (mu0*Rmid**2)                
               ptfbyp = ptflow/p_inext
               if (ptfbyp.gt.10) then ! put a physical lmit of ptfbyp
                  ptfbyp = 10.0
               end if
               pexp=dexp(ptfbyp*((Rt3(inext)/Rmid)**2-1))
               ptfprim =ptf0*ptfout*(1-psin**ptfin)
     1              **(ptfout-1)*(-ptfin)*psin**(ptfin-1)
     2              /(mu0*Rmid**2) !ptf0 includes (mu0*Rmid**2)
               pprim_inext1=pexp*(pprim_inext+
     1              ((Rt3(inext)/Rmid)**2-1)
     2              *(ptfprim-ptfbyp*pprim_inext))
            case(22)            !ptfmax*exp(-((psin-ptfloc)/ptfw)**2))
               ptflow =ptfmax*exp(-((psin-ptfloc)/ptfw)**2)
     1              /(mu0*Rmid**2) !ptfmax includes (mu0*Rmid**2)
               ptfbyp = ptflow/p_inext
               if (ptfbyp.gt.10) then
                  ptfbyp = 10.0
               end if
               pexp=dexp(ptfbyp*((Rt3(inext)/Rmid)**2-1))
               ptfprim =ptfmax*exp(-((psin-ptfloc)/ptfw)**2)
     1              *(-2.0d0)*(psin-ptfloc)/ptfw**2
     1              /(mu0*Rmid**2) !ptfmax includes (mu0*Rmid**2)
               pprim_inext1=pexp*(pprim_inext+
     1              ((Rt3(inext)/Rmid)**2-1)
     2              *(ptfprim-ptfbyp*pprim_inext))
            case (3)            !from table            
               call ispline(-psii(inext),-psiiEF, ptflowEF 
     1              , btflow, ctflow, dtflow, npsi , ptflow)
               ptfbyp = ptflow/p_inext
               pexp=dexp(ptfbyp*((Rt3(inext)/Rmid)**2-1))
               call chderiv(psin,nchq, chcoefftf,ptfprim)
               ptfprim=ptfprim/(mu0*Rmid**2) !chcoefftf includes (mu0*Rmid**2)
               pprim_inext1=pexp*(pprim_inext+
     1                 ((Rt3(inext)/Rmid)**2-1)
     2              *(ptfprim-ptfbyp*pprim_inext))
            end select

            select case (iftype)
            case(0)             ! solovev solution      
               ffprim_inext=0.0d0
            case(1)             ! profiles by parameters 
               ffprim_inext=FF0*(1-psin**ffin)**ffout 
            case(11)            ! profiles by parameters   
               ffprim_inext=-FF0*ffin*ffout*(1-psin**ffin)
     1              **(ffout-1.0d0)*psin**(ffin-1.0d0)
            case (2,3)          ! profile by table or EFIT
               if (isplch.eq.0) then
                  call ispline(-psii(inext),-psiiEF, ffprimEF, bffprim 
     1                 ,cffprim, dffprim,npsi,ffprim_inext)
               else if (isplch.eq.1) then
                  call chfit(psin,nchq,chcoeffff,ffprim_inext)
               end if
            end select
            
            if ((iscale.eq.5).or.(iscale.eq.6)) then
               fsol(inext)=lambda**2*(-mu0*Rt3(inext)**2.0d0
     1              *pprim_inext1)-lambda*ffprim_inext
            else
               fsol(inext)=lambda*(-mu0*Rt3(inext)**2.0d0
     1              *pprim_inext1-ffprim_inext)
            end if
         end do
      end if

      write(*,*) 'iterataion index =', iiter

      ! right hand side of poisson solver
      do inext = 1, ntot
 
         f(inext)=(fsol(inext)/(Rt3(inext)**(0.5d0))
     1           +3.0d0/4.0d0*(psii(inext)/Rt3(inext)**2.5d0)) 
     2           *(dzdw3(inext)*conjg(dzdw3(inext)))
     3           -halpha**2.0d0*psii(inext)/dsqrt(Rt3(inext))
      end do

      ! call the poisson solver
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
     2        ,dpsim2)   !JunHSong

         maxpsiongird=maxpsi/psim
         maxpsi=psim
         Rmaxis=Rt3m
         Zmaxis=Zt3m

         dist2=(Rmaxis-Rmap0)**2+(Zmaxis-Zmap0)**2
       else ! magnetic axis on the grid
         Rmaxpsi=Rt3(inext_maxpsi)
         Zmaxpsi=Zt3(inext_maxpsi)
         Rmaxis=Rmaxpsi
         Zmaxis=Zmaxpsi

         dist2=(Rmaxpsi-Rmap0)**2+(Zmaxpsi-Zmap0)**2
       
      end if
      psif(1:ntot) = psif(1:ntot)/maxpsi
      lambda = lambda/maxpsi

      psii(1:ntot) = psif(1:ntot)

      write(*,*) "(Rmap0,Zmap0): disk center=",Rmap0,Zmap0
      write(*,*) "(Rmaxpsi,Zmaxpsi): maxpsi on grid =",Rmaxpsi,Zmaxpsi
      write(*,*) "(Rmaxis,Zmaxis): maxpsi off grid =",Rmaxis,Zmaxis
      write(*,*) "dist2,epsmaxdist=",dist2,epsmaxdist
      diffpsi=1.0d0-psii(1)
      diffpsich=(1.0d0-dcos(pi/2.0d0/(nchq+1)))/2.0d0
      write(*,*) "diffpsi,diffpsich=",diffpsi,diffpsich

      if (iremapon.eq.1) then
         if ((dist2.gt.epsmaxdist).or.
     1        (dabs(diffpsi).gt.0.9*dabs(diffpsich))) then
             iremap = 1
         else if (iremap.eq.1) then
            call checkcontour(nchq,psiichq,nr,ntheta,rnd,tnd
     1           ,psii,icheck)
            write(*,*) 'icheck',icheck
            if (icheck.eq.1) then
               iremap = 0
            else 
               epsmaxdist=epsmaxdist*0.1d0
               iremap = 1
            end if
         end if
      else
         iremap=0
      end if

      if (iremap.eq.1) then
         Rmap0 = Rmaxis   !
!       JunHsong: remove the shift by findmagaxis algorithm corrections
!        -0.1*diffpsich !Rmaxpsi: put some distance to find magaxis routine
         Zmap0 = Zmaxis         !Zmaxpsi
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
     1     -u(1:ntot)*(dRdr(1:ntot)**2)/(4*sqrt(Rt3(1:ntot)**3))   
     1     +u(1:ntot)*dRdrr(1:ntot)/(2*sqrt(Rt3(1:ntot)))

      write(*,*) '||psif(inext)|| ',maxpsi
      write(*,*) '||psif(inext)-psii(inext)|| ',maxdpsi
      write(*,*) 'lambda', lambda
      write(*,*) 'iremap ', iremap
      write(fiter,*) 'iter: lambda,maxdpsi,psi0',iiter,lambda,
     1        maxdpsi,(1.0d0/lambda)+psiB, iremap

      return
      end

