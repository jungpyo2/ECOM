!--------------------------------------------------------
! readbf.f: This subroutine is to read the given boundary in terms of R,z
! and give the equispaced arclength data for conformal mapping
!--------------------------------------------------------
      subroutine readbf    
      
      use arrays, only:  allocarray1, write_namelists,ibtype
      use arrays, only:  nt2, iiter, iremap, file_bc, iptype
      use arrays, only:  ibtype, nt1, kLag, eps, epsLag
      use arrays, only:  R0, Z0, Rmid, Rmax, Rmin
      use arrays, only:  nin, Rin, Zin, Tin
      use arrays, only:  Rt1, Zt1, T1, dRt1, dZt1, dT1, T1_mil, T1_mil2
      use arrays, only:  csol,reps,rkappa,delta,q0,F0,d1,d2,d3
      use arrays, only:  npsi, psiEF, psiiEF, pprimEF, ffprimEF
      use arrays, only:  allocarrayProf, file_prof, iprintmap
      use arrays, only:  iecom, iqtype, allocarrayQprof
      use arrays, only:  file_qprof,nqpsi,psiiq,qpsiin
      use arrays, only:  ijtype,file_jprof,allocarrayJprof
      use arrays, only:  njpsi,psiij,jpsiin, isymud, Rmap0, Zmap0
      use arrays, only:  iptable, rhoEF, presEF, fpolEF
c
c     read the fixed boundary from a file
      implicit double precision(a-h, o-z)
      integer i, j, k,  unitf, ij, ik
      real *8 tfunc(100 000), thin(10 000), coeff(10 000)
      real *8 Rint(10 000), Zint(10 000),Tint(10 000)
      real *8 dRint(10 000), dZint(10 000),Thint(10 000)
c      real *8 temp1(10 000), temp2(10 000),temp3(10 000)
      real *8 Rintemp(10 000), Zintemp(10 000), qpsitmp(1000)
      real *8 jpsitmp(1000),psinj(1000), rho(1000), pres(1000)
      real *8 psin(1000), psinq(1000), pprim(1000), ffprim(1000)
      real *8 psiichr(100), cftmr(10000), spbxr(10000), spxbr(10000)
      real *8 spdefr(100),thintr(100), Rintr(100), Zintr(100)
      real *8 dRintr(100), dZintr(100), dl(100), fpol(1000)

      integer *4 isort(10 000)
   
      

      done=1
      pi=atan(done)*4.0d0
      ima=(0,1)
c
      call write_namelists 
      write(*,*) 'iptype, ibtype',iptype,ibtype
      if ((iptype.eq.3) .or. (ibtype.eq.3)) then
         call readEfitG  !read from EFIT output
      end if

      if (iptype.eq.0) then !Solovev solution
         rfac = rkappa*F0/(2.0d0*R0**3*q0)
         d3 = rfac/rkappa**2/(-4.0d0)
         d2 = rfac*(-R0**2/2.0d0)
         d1 = rfac*(-reps**2*R0**4+R0**4/4.0d0)
         csol = (rfac/4.0d0-d3)*8.0d0 
         write(*,*) 'd1,d2,d3',d1,d2,d3
         write(*,*) 'csol',csol,reps,rkappa,delta,d1,d2,d3
      end if
      
      select case (ibtype)
         case (0) !Solovev solution
c           call Solovev_param(csol,reps,rkappa,delta,d1,d2,d3)
        
        
            rl=pi*2
            eps7=1.0d-14 

            nin=nt1
            call allocarray1

            call anaresa(ier,funcurve_solovev,rl,nt1,eps7,tfunc,h,rltot)
      
            write(*,*) 'h and rltot', h, rltot
        
            if(ier .ne. 0) stop
c 
c       load up the array curv
c
            dt1=h
            do i=1,nt1+1 
c               tfunc(i)=(i-1)*2*pi/nt1
               call funcurve_solovev(tfunc(i), Rt1(i), Zt1(i)
     1             , dRt1(i), dZt1(i))
                t1(i)=(i-1)*h 
c               t1(i)=tfunc(i)
            end do
        ! note that dRdt1 and dZt1 is derivative in terms of the partmeter t 
        ! instead of the arclength t1
        ! The derivative in terms of the arclengh is given 
        ! by dR/ds=dR/dt/(ds/dt) in szego_fft0

            tin(1:nt1)=t1(1:nt1)
            Rin(1:nt1)=Rt1(1:nt1) 
            Zin(1:nt1)=Zt1(1:nt1) 

            Rmax=R0*dsqrt(1.0d0+2.0d0*reps)
            Rmin=R0*dsqrt(1.0d0-2.0d0*reps)
            Rmid= (Rmax+Rmin)/2.0d0

         case (1) !Miller equilibrium           
c
c       resample the user-specified curve funcurve and get total arc length
c 
            rl=pi*2
            eps7=1.0d-14
            nin=nt1
            call allocarray1
            call anaresa(ier,funcurve,rl,nt1,eps7,tfunc,h,rltot)
            write(*,*) 'conformal mapping: h and rltot', h, rltot
            if(ier .ne. 0) stop
c 
c       load up the array curv
c
            dt1=h
            do i=1,nt1+1 
               call funcurve(tfunc(i), Rt1(i), Zt1(i), dRt1(i), dZt1(i))
               t1(i)=(i-1)*h
            end do
        ! note that dRdt1 and dZt1 is derivative in terms of the partmeter t 
        ! instead of the arclength t1
        ! The derivative in terms of the arclengh is given 
        ! by dR/ds=dR/dt/(ds/dt) in szego_fft0

            tin(1:nt1)=t1(1:nt1)
            Rin(1:nt1)=Rt1(1:nt1) 
            Zin(1:nt1)=Zt1(1:nt1) 
c
            t1_mil(1)=0.0
            t1_mil2(1)=0.0
            do i=2, nt1
               t1_mil(i)=2.0*pi/nt1*(i-1)
               x=dcos(t1_mil(i)+alpha*dsin(t1_mil(i)))
               y=rkappa*dsin(t1_mil(i))

               t1_mil2(i) =atan2(y,x)
               if (t1_mil2(i).lt.0) then 
                  t1_mil2(i)=2.0d0*pi+t1_mil2(i)
               end if
            end do
            t1_mil(nt1+1)=2.0*pi
            t1_mil2(nt1+1)=2.0*pi

            write(*,*) 'theta for miller',t1_mil2(1:nt1+1)
            
            Rmax=R0*(1.0d0+reps)  
            Rmin=R0*(1.0d0-reps) 
            Rmid=R0  
 
         case (2)  ! read from the table in file_bc
            
            unitf=7
            open(unitf, file=file_bc, status='unknown')
            read(unitf,*) 
            read(unitf,*) nin
            read(unitf,*)
            read(unitf,*) R0, Z0
            read(unitf,*)
         
            call allocarray1
       
            Rmap0=R0
            Zmap0=Z0
            Rmax=-1.0d0         !Change Z0 to have maximum R point in midplane 
            do i=1,nin+1
               if (Rin(i).gt.Rmax) then
                  Rmax=Rin(i)
                  Z0=Zin(i)
               end if
            end do

            Tin(1)=0.0 
            do i=1,nin+1
              read(unitf,*) Rin(i), Zin(i)
!             calculate arc length
               thin(i)=atan2((Zin(i)-Z0),(Rin(i)-R0))
               if (thin(i).lt.0) then 
                  thin(i)=2.0d0*pi+thin(i)
               end if
               isort(i)=i
c               if (i.eq.nin+1) then
c                  thin(i)=thin(1)+2.0d0*pi
c               end if
               
c               if (i>1) then
c                  Tin(i)=Tin(i-1)+dsqrt((Rin(i)-Rin(i-1))**2 
c     1              +(Zin(i)-Zin(i-1))**2)
c               end if
            end do
            close(unitf)

            call qsortinc(nin,thin(1:nin),isort(1:nin))  
            do i=1, nin
               Rintemp(i)=Rin(isort(i))
               Zintemp(i)=Zin(isort(i))
            end do
            Rin(1:nin)=Rintemp(1:nin)
            Zin(1:nin)=Zintemp(1:nin)
            Rin(nin+1)=Rin(1)
            Zin(nin+1)=Zin(1)
            thin(nin+1)=thin(1)+2*pi

c            nint=400
            nint=2*nt2          !nt2*20
            do i=1,nint+1
               thint(i)=2.0d0*pi*(i-1)/nint+thin(1)
            end do
            kLag2=8
            call IntLag(nin+1,thin,Rin,Zin,kLag,nint+1,thint,
     1           Rint,Zint,dRint, dZint, epsLag)
        
            Rmax=Rint(1)
            Rmin=Rint(nint/2+1) !Assumed even nint
            Rmid=(Rmax+Rmin)/2.0d0 

            Tint(1)=0.0d0
            kchebr=8
            call chsetupq(kchebr, psiichr, cftmr, spbxr, spxbr
     1           , spdefr)  
            psiichr(1:kchebr)=1.0d0-psiichr(1:kchebr)
            do i=2,nint+1
               thintr(1:kchebr)=thint(i-1)
     1              +(thint(i)-thint(i-1))*psiichr(1:kchebr)
               call IntLag(nin+1,thin,Rin,Zin,kLag,kchebr,thintr,
     1              Rintr,Zintr,dRintr, dZintr, epsLag)
               dl(1:kchebr)=dsqrt(dRintr(1:kchebr)**2
     1              +dZintr(1:kchebr)**2)
c     call chftransq(chcoeff,dl,kchebr,cftmr)
               dtint=0.0d0
               do j=1,kchebr
                  dtint=dtint+spdefr(j)*dl(j)
               end do
               dtint= dtint/2.0d0*(thint(i)-thint(i-1))
               Tint(i)=Tint(i-1)+dtint
            end do

            dt1=Tint(nint+1)/nt1
            do i=1,nt1+1
               T1(i)=(i-1)*dt1
            end do

            if (isymud.eq.1) then !assymed up-down symmetry
               nt1half=nt1/2
               call IntLag(nint+1,Tint,Rint,Zint,kLag,nt1half+1, T1,
     1              Rt1,Zt1,dRt1, dZt1, epsLag) 
               dRt1(1)=0.0d0 
               dRt1(nt1half+1)=0.0d0
               Zt1(1)=0.0d0
               Zt1(nt1half+1)=Zt1(1)
               dZt1(1)=1.0d0
               dZt1(nt1half+1)=-dZt1(1)
               
               
               do i=1,nt1half
                  Zt1(nt1half+1+i)=-Zt1(nt1half+1-i)
                  dZt1(nt1half+1+i)=dZt1(nt1half+1-i)
                  Rt1(nt1half+1+i)=Rt1(nt1half+1-i)
                  dRt1(nt1half+1+i)=-dRt1(nt1half+1-i)
               end do
               
            else !assymed up-down asymmetry
               call IntLag(nint+1,Tint,Rint,Zint,kLag,nt1+1,T1,
     1              Rt1,Zt1,dRt1, dZt1, epsLag)
            end if
         
            open(7, file='RZ_input.dat', status='unknown')
            write(7,*) ' R,  Z, T '
            do i=1,nin+1
               write(7,*) Rin(i), Zin(i), Tin(i)
            end do
            close(7)

         case (3)  !EFIT
        
            write(*,*) 'R0,Z0',R0,Z0
            write(*,*) 'Rin',Rin

            Rmax=-1.0d0         !Change Z0 to have maximum R point in midplane 
            do i=1,nin+1
               if (Rin(i).gt.Rmax) then
                  Rmax=Rin(i)
                  Z0=Zin(i)
               end if
            end do

            Tin(1)=0.0 
            do i=1,nin+1
               thin(i)=atan2((Zin(i)-Z0),(Rin(i)-R0))
               isort(i)=i
               if (thin(i).lt.0) then 
                  thin(i)=2.0d0*pi+thin(i)
               end if
            end do
c            write(*,*) 'thin',thin(1:nin+1)

            call qsortinc(nin,thin(1:nin),isort(1:nin))  
            do i=1, nin
               Rintemp(i)=Rin(isort(i))
               Zintemp(i)=Zin(isort(i))
            end do
            Rin(1:nin)=Rintemp(1:nin)
            Zin(1:nin)=Zintemp(1:nin)
 
            Rin(nin+1)=Rin(1)
            Zin(nin+1)=Zin(1)
            thin(nin+1)=thin(1)+2*pi

            nint=2*nt2 !nt2*20
            do i=1,nint+1
               thint(i)=2*pi*(i-1)/nint+thin(1)
            end do

            kLag2=8
            call IntLag(nin+1,thin,Rin,Zin,kLag,nint+1,thint,
     1           Rint,Zint,dRint, dZint, epsLag)

            Rmax=Rint(1)
            Rmin=Rint(nint/2+1) !Assumed even nint
            Rmid=(Rmax+Rmin)/2.0d0 

            Tint(1)=0.0d0
            kchebr=8
            call chsetupq(kchebr, psiichr, cftmr, spbxr, spxbr
     1           , spdefr)  
            psiichr(1:kchebr)=1.0d0-psiichr(1:kchebr)
            do i=2,nint+1
               thintr(1:kchebr)=thint(i-1)
     1              +(thint(i)-thint(i-1))*psiichr(1:kchebr)
               call IntLag(nin+1,thin,Rin,Zin,kLag,kchebr,thintr,
     1              Rintr,Zintr,dRintr, dZintr, epsLag)
               dl(1:kchebr)=dsqrt(dRintr(1:kchebr)**2
     1              +dZintr(1:kchebr)**2)
               dtint=0.0d0
               do j=1,kchebr
                  dtint=dtint+spdefr(j)*dl(j)
               end do
               dtint= dtint/2.0d0*(thint(i)-thint(i-1))
               Tint(i)=Tint(i-1)+dtint
            end do

            dt1=Tint(nint+1)/nt1
            do i=1,nt1+1
               T1(i)=(i-1)*dt1
            end do

    
            call IntLag(nint+1,Tint,Rint,Zint,kLag,nt1+1,T1,
     1           Rt1,Zt1,dRt1, dZt1, epsLag)

        
         
      end select             !end select ibtype

c-----------------------------
      open(8, file='RZ_intp.gen', status='unknown')
      write(8,*) 'nt1'
      write(8,*)  nt1
      write(8,*) ' R0, Z0 '
      write(8,*) R0, Z0
      write(8,*) ' Rmid '
      write(8,*) Rmid
      write(8,*) ' R1, Z1 '
      do i=1,nt1+1
         write(8,*) Rt1(i), Zt1(i)
      end do
      close(8)

      open(8, file='RZ_intp1.dat', status='unknown')
      write(8,*) ' T, R, Z, dR/dT, dZ/dT '
      do i=1,nt1+1
         write(8,*) T1(i), Rt1(i), Zt1(i), dRt1(i), dZt1(i)
      end do
      close(8)

      if (iprintmap.eq.1) then
         itype1=1
         itype2=2
         itype3=3
         iw=10
         call quaplot2(iw,Rt1,Zt1,nt1+1,itype2,Rin,Zin,nin+1,itype2,
     1        'Lagrange interpolation for uniform arc length*')
      end if

      if (iptype.eq.2) then  !psin, p', FF' are given by a table in 'file_prof'
         if (iptable.eq.0) then
            unitf=7
            open(unitf, file=file_prof, status='unknown')
            read(unitf,*) 
            read(unitf,*) npsi
            read(unitf,*)
            read(unitf,*) (psin(i),i=1,npsi)
            read(unitf,*)
            read(unitf,*) (pprim(i),i=1,npsi)
            read(unitf,*)
            read(unitf,*) (ffprim(i),i=1,npsi)
            close (unitf)

            call allocarrayProf
            psiiEF(1:npsi)=1-psin(1:npsi)
            pprimEF(1:npsi)=pprim(1:npsi)
            ffprimEF(1:npsi)=ffprim(1:npsi)

            write(*,*) 'npsi',npsi
            write(*,*) 'pprimEF',pprimEF
            write(*,*) 'ffprimEF',ffprimEF

         else if (iptable.eq.1) then
            unitf=7
            open(unitf, file=file_prof, status='unknown')
            read(unitf,*) 
            read(unitf,*) npsi
            read(unitf,*)
            read(unitf,*) (psin(i),i=1,npsi)
            read(unitf,*)
            read(unitf,*) (pres(i),i=1,npsi)
            read(unitf,*)
            read(unitf,*) (ffprim(i),i=1,npsi)
            close (unitf)

            call allocarrayProf
            psiiEF(1:npsi)=1-psin(1:npsi)
            presEF(1:npsi)=pres(1:npsi)
            ffprimEF(1:npsi)=ffprim(1:npsi)

            write(*,*) 'npsi',npsi
            write(*,*) 'presEF',presEF
            write(*,*) 'ffprimEF',ffprimEF

         else if (iptable.eq.2) then
            unitf=7
            open(unitf, file=file_prof, status='unknown')
            read(unitf,*) 
            read(unitf,*) npsi
            read(unitf,*)
            read(unitf,*) (rho(i),i=1,npsi)
            read(unitf,*)
            read(unitf,*) (pres(i),i=1,npsi)
            read(unitf,*)
            read(unitf,*) (fpol(i),i=1,npsi)
            close (unitf)

            call allocarrayProf
            rhoEF(1:npsi)=rho(1:npsi)
            presEF(1:npsi)=pres(1:npsi)
            fpolEF(1:npsi)=fpol(1:npsi)

            write(*,*) 'rho',rhoEF
            write(*,*) 'presEF',presEF
            write(*,*) 'fpolEF',fpolEF
         end if
      end if

      if ((iecom.eq.2) .and. (iqtype.eq.2)) then  !psinq, q are given by a table in 'file_qprof'
         unitf=7
         open(unitf, file=file_qprof, status='unknown')
         read(unitf,*) 
         read(unitf,*) nqpsi
         read(unitf,*)
         read(unitf,*) (psinq(i),i=1,nqpsi)
         read(unitf,*)
         read(unitf,*) (qpsitmp(i),i=1,nqpsi)
         close (unitf)

         call allocarrayQprof
         psiiq(1:nqpsi)=1-psinq(1:nqpsi)
         qpsiin(1:nqpsi)=qpsitmp(1:nqpsi)
         write(*,*) 'nqpsi',nqpsi
         write(*,*) 'psiiq',psiiq
         write(*,*) 'qpsiin',qpsiin

      else if ((iecom.eq.1) .and. (ijtype.eq.2)) then  !psinq, q are given by a table in 'file_jprof'
         unitf=7
         open(unitf, file=file_jprof, status='unknown')
         read(unitf,*) 
         read(unitf,*) njpsi
         read(unitf,*)
         read(unitf,*) (psinj(i),i=1,njpsi)
         read(unitf,*)
         read(unitf,*) (jpsitmp(i),i=1,njpsi)
         close (unitf)

         call allocarrayJprof
         psiij(1:njpsi)=1-psinj(1:njpsi)
         jpsiin(1:njpsi)=jpsitmp(1:njpsi)
         write(*,*) 'njpsi',njpsi
         write(*,*) 'psiij',psiij
         write(*,*) 'jpsiin',jpsiin
      end if

      write(*,*) 'Finish reading'
      return
      end   !end of the subrouine readbf




       subroutine funcurve(t,x,y,dxdt,dydt)

       use arrays, only: R0,Z0,reps,rkappa, delta 
       implicit real *8 (a-h,o-z)
c        complex *16 rhs(3),w(3,3),w7(1000), sol(3)
 
        save
c
        done=1
        pi=4*atan(done)
c
c       barbell
c
c        a=6
c
c        x=a*cos(t)
c        y=sin(t)*(5*cos(2*t)+6)
c        dxdt=-a*sin(t)
c        dydt=cos(t)*(5*cos(2*t)+6)-sin(t)*sin(2*t)*5*2
c
cccc        y=y/8
cccc        return
c
c       cloud
c
c        x=5*cos(t)+.3*cos(10*t)
c        y=5*sin(t)+.3*sin(10*t)
c        dxdt=-5*sin(t)-3*sin(10*t)
c        dydt=5*cos(t)+3*cos(10*t)
c
cccc        return
c
c       ellipse
c
c        a=1
c        b=1.5d0
c
c        x=a*cos(t)
c        y=b*sin(t)
c        dxdt=-a*sin(t)
c        dydt=b*cos(t)

cccc        return
c
c       the plasma cross section
c

c        reps=0.32d0
c        rkappa=1.7d0
c        delta=0.33d0
c        goto 2000
c
c        reps=.78d0
c        rkappa=2.0d0
c        delta=.35d0
c        alpha=asin(delta)
c        goto 2000
c
c        reps=.95d0
c        rkappa=1.0d0
c        delta=.2d0
c        alpha=asin(delta)
c        goto 2000
c
c 2000   continue 

        
        alpha=dasin(delta)
c        write(*,*) 'alpha',R0,Z0,reps,rkappa, delta, alpha 
        x=R0*(1+reps*dcos(t+alpha*dsin(t)))
        y=Z0+R0*reps*rkappa*dsin(t)

        dxdt=-R0*reps*dsin(t+alpha*dsin(t))*(1+alpha*dcos(t))
        dydt=R0*reps*rkappa*dcos(t)
c
c        write(*,*) x,y,dxdt,dydt
        return
        end

      subroutine funcurve_Solovev(t,x,y,dxdt,dydt)

      use arrays, only:  csol,reps,rkappa,delta,d1,d2,d3
      use arrays, only: reps, rkappa, delta, R0, Z0
      implicit real *8 (a-h,o-z)

c        integer, parameter :: QR_K = selected_real_kind (32)
      real (kind=10) :: temp1, temp2 !  MyReal
c        real * 16 temp1, temp2
         
c        complex *16 rhs(3),w(3,3),w7(1000), sol(3)    
      save
c
c
c        eps=1d-7
c        t=t+eps 
c        done=1
c        pi=4*atan(done)
c       
     
c         c=csol
c        p=d2/(2*(c/8.0d0+d3))
c        q=dsqrt((-4.0d0*d3)/(c/8.0d0+d3))
c        r=dsqrt((d2**2/(4*(c/8.0d0+d3))-d1)/(c/8.0d0+d3))
c        write(*,*) 'p,q,r',p,q,r
c        x=dsqrt(r*dcos(t)-p)
c        y=r*dsin(t)/(q*x)

c        dxdt=-r*dsin(t)/(2*x)
c        dydt=(r*dcos(t)*x+r**2*dsin(t)**2/(2*x))/(q*x**2)

        x=R0*dsqrt(1.0d0+2.0d0*reps*dcos(t))
        y=Z0+rkappa*reps*R0**2/x*dsin(t)

        dxdt=-reps*R0**2*dsin(t)/(x)
        dydt=-y*dxdt/x+rkappa*reps*R0**2/x*dcos(t)

c        write(*,*) 'd1,d2,d3csd',d1,d2,d3

c        x=1.0d0+reps*dcos(t)
c        dxdt=-reps*dsin(t)
c        dxdtt=-reps*dcos(t)
       

c        xo=(1+reps)
          
c        yo2d1=(2.0d0*(c/8.0+d3)*xo-2.0d0*d1/xo**3)/(4.0d0*d3)
c        yo2d2=(2.0d0*(c/8.0+d3)+6.0d0*d1/xo**4)/(4.0d0*d3)
c        yo2d3=(-24.0d0*d1/xo**5)/(4.0d0*d3)
c        yo2d4=(120.0d0*d1/xo**6)/(4.0d0*d3)

c        xi=(1+reps)

c        yi2d1=(2.0d0*(c/8.0+d3)*xi-2.0d0*d1/xi**3)/(4.0d0*d3)
c        yi2d2=(2.0d0*(c/8.0+d3)+6.0d0*d1/xi**4)/(4.0d0*d3)
c        yi2d3=(-24.0d0*d1/xi**5)/(4.0d0*d3)
c        yi2d4=(120.0d0*d1/xi**6)/(4.0d0*d3)
      

c        if(t.lt.pi) then
c           ytemp2=((c/8.0d0+d3)*x**2+d2+d1/x**2)/(4.0d0*d3)
c           if(ytemp2.ge.1e-8) then
c              y=dsqrt(ytemp2) !((c/8.0d0+d3)*x**2+d2+d1/x**2)/(4.0d0*d3))
c           else if (ytemp2.ge.1e-26) then
c              dx=x-xo   
c              ytay2=y2d1*dx+y2d2*dx**2/2.0d0+y2d3*dx**3/6.0d0
c     1             +y2d4*dx**4/24.0d0
c              if(ytay2.gt.0) then
c                 y=dsqrt(ytay2)
c              else 
c                 y=0.0d0
c              end if
c           else
c              y=0.0d0
c           end if
c           dy2dx=((c/8.0d0+d3)*x-d1/x**3)/(2.0d0*d3)
c           if (abs(y).gt.1.0e-7) then
c              dydt=dy2dx/(2.0d0)*dxdt/y!dsqrt(dxdt**2/

c           else
c              dydt=sign(1.0d0,dcos(t))*dsqrt(dy2dx*dxdtt/2.0d0)
c              write(*,*) 'y=0,t,dy2dx,dxdtt,dydt',t,dy2dx,dxdtt,dydt
c           end if
c        else
c           if(((c/8.0d0+d3)*x**2+d2+d1/x**2)/(4.0d0*d3).ge.0) then
c              y=-dsqrt(((c/8.0d0+d3)*x**2+d2+d1/x**2)/(4.0d0*d3))
c           else
c              y=0.0d0
c           end if
c          dy2dx=((c/8.0d0+d3)*x-d1/x**3)/(2.0d0*d3)!(c/8.0+d3)/(2.0*d3)*x-d1/(2.0*d3)*x**(-3)
c           if (abs(y).gt.1.0e-7) then
c              dydt=dy2dx/(2.0d0)*dxdt/y!
c           else 
c              dydt=sign(1.0d0,dcos(t))*dsqrt(dy2dx*dxdtt/2.0d0)
c              write(*,*) 'y=0,t,dy2dx,dxdtt,dydt',t,dy2dx,dxdtt,dydt
c           end if
c        end if

        return
        end


        subroutine Solovev_param(c,reps,rkappa,delta,d1,d2,d3)
        implicit real *8 (a-h,o-z)
        complex *16 rhs(3),w(3,3),w7(1000), sol(3)
        save
c
 
        rhs(1)=(1,0)*(-c/8)*(1+reps)**4
        rhs(2)=(1,0)*(-c/8)*(1-reps)**4
        rhs(3)=(1,0)*(-c/8)*(1-delta*reps)**4

        w(1,1)=(1,0); w(1,2)=(1,0)*(1+reps)**2
        w(1,3)=(1,0)*(1+reps)**4

        w(2,1)=(1,0); w(2,2)=(1,0)*(1-reps)**2
        w(2,3)=(1,0)*(1-reps)**4

        w(3,1)=(1,0);w(3,2)=(1,0)*(1-delta*reps)**2
        w(3,3)=(1,0)*((1-delta*reps)**4-
     1     4*rkappa**2*reps**2*(1-delta*reps)**2)

        call cqrdecom(w,3,w7,rcond)
        write (*,*) 'rcond for conformal mapping=',rcond
c     
        call cqrsolve(3,w7,rhs,sol)
        d1=sol(1)
        d2=sol(2)
        d3=sol(3)

        write(*,*) 'd1,d2,d3',d1,d2,d3
        return
        end

      subroutine readEfitG

      use arrays, only:  allocarray1, allocarrayEFIT
      use arrays, only:  ibtype,iptype
      use arrays, only:  npsi, psi0, psiB, psiEF, pprimEF, ffprimEF
      use arrays, only:  fpolEF, qpsiEF, psiiEF, file_efit
      use arrays, only:  nt1, kLag, eps, epsLag
      use arrays, only:  nin, Rin, Zin, Tin, R0, Z0, Rmap0, Zmap0
      use arrays, only:  Rt1, Zt1, T1, dRt1, dZt1, F0, q0
      use arrays, only:  nrEF, nzEF, REf, ZEF
      implicit double precision(a-h, o-z)
      character*20 case(6), file_out

      real *8 psirz(1000,1000),fpol(1000),pres(1000),ffprim(1000)
      real *8 pprime(1000),qpsi(1000),rbbbs(1000),zbbbs(1000)
      real *8 pressw(1000),pwprim(1000),dmion(1000),rhovn(1000)
      real *8 rlim(1000),zlim(1000)

      neqdsk=10
c      file_efit='g133221.01000'
      write(*,*) file_efit
      open(neqdsk, file=file_efit, status='unknown')
      read (neqdsk,2000) (case(i),i=1,6),idum,nw,nh
c      read (neqdsk,2010) idum,nw,nh
      write(*,*) 'idum',idum,nw,nh
      read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid
      read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr
      read (neqdsk,2020) current,simag,xdum,rmaxis,xdum
      read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
      read (neqdsk,2020) (fpol(i),i=1,nw)
      read (neqdsk,2020) (pres(i),i=1,nw)
      read (neqdsk,2020) (ffprim(i),i=1,nw)
      read (neqdsk,2020) (pprime(i),i=1,nw)
      read (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
      read (neqdsk,2020) (qpsi(i),i=1,nw)
      read (neqdsk,2022) nbbbs,limitr
      write(*,*) 'nbbs', nbbbs,limitr
      read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
      read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)
c      read (neqdsk,2024) kvtor,rvtor,nmass
c      read (neqdsk,2023) kvtor,nmass
c      write(*,*) 'kvtor', kvtor,rvtor,nmass
c      if (kvtor.gt.0) then
c         read (neqdsk,2020) (pressw(i),i=1,nw)
c         read (neqdsk,2020) (pwprim(i),i=1,nw)
c      endif
c      if (nmass.gt.0) then
c         read (neqdsk,2020) (dmion(i),i=1,nw)
c      endif
c      read (neqdsk,2020) (rhovn(i),i=1,nw)
c     
 2000 format (6a8,3i4)
 2010 format (3i4)
 2020 format (5e16.9)
 2022 format (2i5)
 2023 format (i5,i5)
 2024 format (i5,e16.9,i5)
      close(neqdsk)
 
      nefito=11
      file_out='readEfit.out'
      open(nefito, file=file_out, status='unknown')
      write(nefito,*) 'idum, nw, nh',idum, nw, nh
      write(nefito,*) 'nbbbs, limitr', nbbs, limitr
      write(nefito,*) rdim,zdim,rcentr,rleft,zmid
      write(nefito,*) rmaxis,zmaxis,simag,sibry,bcentr
      write(nefito,*) current,simag,xdum,rmaxis,xdum
      write(nefito,*) zmaxis,xdum,sibry,xdum,xdum
      write(nefito,*) 'ffrime'
      write(nefito,*) (ffprim(i),i=1,nw)
      write(nefito,*) 'pprime'
      write(nefito,*) (pprime(i),i=1,nw)
      write(nefito,*) 'qsi'
      write(nefito,*) (qpsi(i),i=1,nw)

      write(nefito,*) 'rbbs,zbbs'
      write(nefito,*) (rbbbs(i),zbbbs(i),i=1,nbbbs)
      write(nefito,*) 'rlim,zlim'
      write(nefito,*) (rlim(i),zlim(i),i=1,limitr)
      write(nefito,*) 'fpol'
      write(nefito,*) (fpol(i),i=1,nw)
      write(nefito,*) 'pressure'
      write(nefito,*) (pres(i),i=1,nw)
   
      PSI0 = SIMAG
      PSIB = SIBRY
      write(nefito,*) 'SIMAG, SIBRY'
      write(nefito,*) SIMAG,SIBRY
      write(nefito,*) 'current [A]'
      write(nefito,*)  current
      if (iptype.eq.3) then
         npsi = nw
         nREF = nw
         nZEF = nh
         call allocarrayEFIT
         do i=1, npsi
            PSIEF(i)= PSI0+(PSIB-PSI0)*(i-1)/(npsi-1)
            psiiEF(i)=(psiEF(i)-psiB)/(psi0-psiB)
            pprimEF(i)=pprime(i)
            ffprimEF(i)=ffprim(i)
            fpolEF(i)=fpol(i)
            qpsiEF(i)=qpsi(i)
            REF(i) = rleft+rdim*(i-1)/(nw-1)
         end do
         do j=1, nh
            ZEF(j) = zmid-zdim/2.0d0+zdim*(j-1)/(nh-1)
         end do
         f0=fpol(1)
         q0=qpsi(1)
      end if
      write(nefito,*) 'psi'
      write(nefito,*) (PSIEF(i),i=1,nw)
      write(nefito,*) 'R'
      write(nefito,*) (REF(i),i=1,nw)
      write(nefito,*) 'Z'
      write(nefito,*) (ZEF(j),j=1,nh)
      write(nefito,*) 'psiRZ'
      write(nefito,*) ((PSIrz(i,j),i=1,nw),j=1,nh)

      close(nefito)

      if (ibtype.eq.3) then
         R0=RMAXIS
         Rmid=RMAXIS
         Z0=ZMAXIS
         Rmap0=RMAXIS
         Zmap0=ZMAXIS
         nin = nbbbs-1
         write (*,*) 'nin',nin
         call allocarray1         
         Tin(1)=0.0 
         do i=1,nin+1
            Rin(i)=rbbbs(i)
            Zin(i)=zbbbs(i) 
         end do
      end if

      return
      end
      
