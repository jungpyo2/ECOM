! For public variables of ECOM

module arrays 
     
  public :: sort_samevalue 
  public :: deallocarray
! Lagder_periodic, Lagder_periodic2  
  public :: IntLag_periodic, IntLag_periodic2, barycentric_w, fpoledge, B0vac
  public :: file_prof, file_bc, file_efit, nt1, ksamp,ksamp2, kLag, nt2, eps, nsub, kcheb
  public :: epsLag, epsLag2, epsiter, epsmag, maxiter, nr, ntot, blength, bnodes
  public :: nin, Rin, Zin, Tin, R0, Z0, Rmap0, Zmap0, Rmid,Rmax,Rmin, Rmaxis, Zmaxis
  public :: iremap, iiter, isymud, iecom, iLag, iprintEFIT
  public :: Rt1, Zt1, T1, dT1, dRt1, dZt1, w1, dwdz1, dwdzz1, wsaved, w7saved
  public :: Rt2, Zt2, T2, dRt2, dZt2, dwdz2, dwdzz2, zk, dzdw2k, dzdww2k
  public :: Rt3, Zt3, rnd, tnd, dzdw3, dzdww3, fp3, gp3, dwdz3, dwdzz3
  public :: csol, reps, rkappa, delta, d1, d2, d3
  public :: npsi, psi0, psiB, psiEF, psiiEF, pprimEF, ffprimEF, fpolEF, qpsiEF, rhoEF, presEF
  public :: nREF, nZEF, REF, ZEF, isplch, irho
  public :: bpprim, cpprim, dpprim, bffprim, cffprim, dffprim, bpres,cpres,dpres,bfpol,cfpol,dfpol

  public :: nchq, psiichq, phiichq, cftmq, spbx, spxb, spdef, qpsich, pprimch, ffprimch,qpsich1
  public :: cftmsub, pprimch0, fpolch0, jpsich, phitot
  public :: wffprim, fpolcon2, ffprimpre, ffprimpre1, ffprimpre2, ffprimgoal, fpolch
  public :: nchq0, qpsich0, psiichq0, cftmq0, spbx0, spxb0, spdef0, jpsich0
!  public :: nchqs, psiichqs, cftmqs, spbxs, spxbs, spdefs, qpsichs, pprimchs, qpsich1s
!  public :: fpolcon2s, ffprimpres, ffprimpre1s, ffprimpre2s, ffprimgoals, fpolchs
!  public :: dZdtmidos, dZdtmidis

  public :: g, psii, f, lambda, maxpsi, maxdpsi, lambdaf, iiterl, psiex, lambdaex
  public :: epsmaxdist, Rmaxpsi, Zmaxpsi, fpolsgn , iiterq, relaxiter, iqconstraint
  public :: dpsidr, dpsidrr, u, ur, uth, urr, urt, utt, ucoeff, urcoeff, urrcoeff
  public :: q0, F0, Fedge, p0, pin, pout, FF0, FFin, FFout,qfac,qpin,qpout,ffac,Fpin,Fpout
  public :: Jpar0, Jpin, Jpout, Jomax, Joloc, Jow, pgmax, pgloc, pgw
  public :: itftype, file_tflow, ptf0, ptfin, ptfout, mach, ptfmax, ptfloc, ptfw
  public :: nflx, psiflx, maxiterF, dpsistartmag, dpsistartnchq
  public :: torcur, betaP,  curden, dlnqdpsi, dlnIcdpsi, shat, curpar 
  public :: ici, icr, ic_max, iptype, iftype, iscale, iconvq, iconvj, ifpol, maxdpsiq, maxdpsij
  public :: ftime, fiter, fout,fiter1, fiter2, fiter3,fiter4,fiter5,fiter6,fiter7,fiter22
  public :: iprintsol, iprintsoldisk, iprintcon, iprintQS, ifitMil, iprintmap
  public :: istability,ibscur, ijBSmodel, isolver, Halpha, ifindmagaxis, iptable,iguessRmap0
  public :: chcoeffp, presch, chcoefftf, ptfprimch, ptflowch, rhoch, chcoeffpp, chcoeffff
  public :: ptflowEF, btflow, ctflow, dtflow, chcoeffpp0, dpsiidrhoch,chcoeffprho,chcoefff2rho
 
  public :: iqtype, file_qprof, nqpsi, psiiq, qpsiin, bqpsi, cqpsi, dqpsi
  public :: ijtype, file_jprof, njpsi, psiij, jpsiin, bjpsi, cjpsi, djpsi
  public :: nchy, psiichy, cftmy, spybx, spyxb,spydef,fpass, etai, Zi, jBSch
  public :: nue_star, nui_star, ne, Te, Vloop, Vloop0, jOhmch
  public :: rndm, tndm, dpsidrm, dpsidtm, dpsidrrm ,dpsidrtm, dpsidttm, Rt3m, Zt3m, psim

  public :: dense0, denseb, rn1de, rn2de, te0, teb, rn1te, rn2te, ti0, tib, rn1ti, rn2ti
  public :: itrinity, iogyropsi, ntpsi, tpsi, trho, nttin, ttin, dpsidrho
  public :: gradpar,Rtrin,Ztrin,btot,bpol,gbdrift,gbdrift0, cvdrift,cvdrift0
  public :: shatlocal1,shatlocal2,shatlocal3
  public :: zw,gcoeff, pt, gt, fr, fi, un, unr
  public :: dfdzbf,dfdzzbf,fbf,zbf,zmbf, dwdz1m, dwdzz1m, w7, w8, w9 

  character (20) :: file_prof, file_qprof,file_jprof,file_bc, file_efit, file_tflow
  integer :: nin, nt1, ksamp,ksamp2, kLag, nt2, nsub, kcheb, maxiter 
  integer :: nr, ntot, npsi, iremap, iremapon, iiter, nchq, iiterl, isymud, nchq0
  integer :: iiterq, ici, icr, ic_max, iptype, iftype, ibtype,itftype, nflx
  integer :: ftime, fiter, fout, iscale, iecom, iqtype, nqpsi, njpsi, isplch
  integer :: iconvq,iconvj, ijtype, maxiterF, iLag, ifindmagaxis, ifpol, irho
  integer :: istability,ibscur, nchy, ijBSmodel, isolver, iprintEFIT, itrinity,iogyropsi
  integer :: fiter1, fiter2, fiter3,fiter4,fiter5,fiter6,fiter7,fiter22
  integer :: iprintsol,iprintsoldisk, iprintcon, iprintQS, ifitMil, iprintmap
  integer :: relaxiter, nrEF, nzEF, ntpsi, nttin, iptable, iqconstraint,iguessRmap0

  real * 8 :: R0, Z0, Rmaxis, Zmaxis, Rmid, Rmax, Rmin
  real * 8 :: eps, epsLag, epsLag2, epstri, epsiter, fpolsig, ptfmax, ptfloc, ptfw, pgmax, pgloc, pgw
  real * 8 :: q0, F0, Fedge, p0, pin, pout, FF0, FFin, FFout,qfac, qpin, qpout,Ffac,Fpin, Fpout
  real * 8 :: Jpar0, Jpin, Jpout, Jomax, Joloc, Jow, Halpha, Vloop0, dRmap0, phitot
  real * 8 :: csol, reps, rkappa, delta, d1, d2, d3, torcur,betaP, fpoledge, B0vac
  real * 8 :: curden, dlnqdpsi, dlnIcdpsi, Rmaxpsi, Zmaxpsi, Rmap0, Zmap0, dT1
  real * 8 :: mach,ptf0, ptfin, ptfout, maxdpsiq, maxdpsij, epsmag, dpsistartmag, dpsistartnchq
  real * 8 :: psi0, psiB, wffprim, epsmaxdist, lambda, maxpsi, maxdpsi, lambdaf, lambdaex
  real * 8 :: rndm, tndm, dpsidrm, dpsidtm, dpsidrrm ,dpsidrtm, dpsidttm, Rt3m, Zt3m, psim
  real * 8 :: epslogzk, epslogzwk, epslogzwwk
  real * 8 ::  dense0, denseb, rn1de, rn2de, te0, teb, rn1te, rn2te, ti0, tib, rn1ti, rn2ti
  
  real * 8, dimension(1000) :: psiflx, Vloop

  real * 8, dimension(:), allocatable :: Rin, Zin, Tin
  real * 8, dimension(:), allocatable :: Rt1, Zt1, T1, dRt1, dZt1, T1_mil, T1_mil2
  real * 8, dimension(:), allocatable :: psiEF,psiiEF, pprimEF, ffprimEF, fpolEF,qpsiEF, REF, ZEF
  real * 8, dimension(:), allocatable :: Rt2, Zt2, T2, dRt2, dZt2
  real * 8, dimension(:), allocatable :: Rt3, Zt3, rnd, tnd, blength, bnodes
  real * 8, dimension(:), allocatable ::  bpprim, cpprim, dpprim, bffprim, cffprim, dffprim, rhoEF, presEF
  real * 8, dimension(:), allocatable :: psiichq, phiichq,cftmq, spbx, spxb, spdef, qpsich, pprimch, ffprimch
  real * 8, dimension(:), allocatable :: cftmsub, shat, fpolcon2,chcoeffpp,chcoeffpp0,chcoeffff
  real * 8, dimension(:), allocatable :: ffprimpre, ffprimpre1, ffprimgoal, qpsich1, fpolch, ffprimpre2
  real * 8, dimension(:), allocatable :: qpsich0, curpar, jpsich, jpsich0,chcoeffprho,chcoefff2rho
  real * 8, dimension(:), allocatable :: psiichq0, cftmq0, spbx0, spxb0, spdef0, fpolch0, pprimch0
!  real * 8, dimension(:), allocatable :: fpolcon2s, dZdtmidos, dZdtmidis
!  real * 8, dimension(:), allocatable :: ffprimpres, ffprimpre1s, ffprimgoals, qpsich1s, fpolchs, ffprimpre2s
  real * 8, dimension(:), allocatable ::  g, psii, f,dpsidr, dpsidrr, u, ur, uth, urr, urt, utt, psiex
  real * 8, dimension(:), allocatable :: chcoeffp, presch, rhoch, dpsiidrhoch,chcoefftf, ptfprimch, ptflowch
  real * 8, dimension(:), allocatable :: ptflowEF, btflow, ctflow, dtflow
  real * 8, dimension(:), allocatable :: psiiq, qpsiin, bqpsi, cqpsi, dqpsi
  real * 8, dimension(:), allocatable :: psiij, jpsiin, bjpsi, cjpsi, djpsi
  real * 8, dimension(:), allocatable :: psiichy, cftmy, spybx, spyxb,spydef
  real * 8, dimension(:), allocatable :: fpass, etai, Zi, jBSch, nue_star, nui_star, ne, Te, jOhmch
  real * 8, dimension(:), allocatable :: tpsi, trho, dpsidrho,ttin, gradpar,Rtrin,Ztrin,btot,bpol !for trinity
  real * 8, dimension(:), allocatable :: gbdrift,gbdrift0,cvdrift,cvdrift0,shatlocal1,shatlocal2,shatlocal3 !for trinity
  real * 8, dimension(:), allocatable :: bpres,cpres,dpres, bfpol,cfpol,dfpol
  complex * 16, dimension(:), allocatable :: w1, dwdz1,dwdzz1, dwdz2, dwdzz2, dzdw3, dzdww3, dwdz3, dwdzz3
  complex * 16, dimension(:), allocatable :: wsaved, w7saved, zk, dzdw2k, dzdww2k
  complex * 16, dimension(:), allocatable :: ucoeff, urcoeff, urrcoeff

  complex * 16, dimension(:), allocatable :: zw,gcoeff  !arrays for elliptic.f
  real * 8, dimension(:), allocatable :: pt, qt, fr, fi, un, unr  !
  complex * 16, dimension(:), allocatable :: dfdzbf,dfdzzbf,fbf,zbf,zmbf, dwdz1m, dwdzz1m, w7, w8  ! arrays for cmap.f
  complex * 16, dimension(:,:), allocatable :: w9
  !real * 8, dimension(:), allocatable :: pt, qt, fr, fi, un, unr  
contains

  subroutine initarray

    implicit none
    
    logical :: file_exists
    character (20) :: filein

    call system_clock(ici, icr, ic_max)   

!   Default values of the global variables
    iecom = 0 ! 0-ECOM (FF' input), 1-ECOMQ (q input), 2-ECOMJ (J par input) 
    iremap = 1
    iremapon = 1
    iiter = 0
    file_bc='ITER.dat' 
    file_efit='g133221.01000'
    file_prof='profile.dat'
    file_tflow='torflow.dat'
    file_qprof='qprofile.dat'
    file_jprof='jprofile.dat'
    nt1 = 32
    kLag=8  ! 8 th order Lagrange interpolation
    eps = 1.0d-14
    epsLag = 1.0d-13
    epsLag2 = 1.0d-13
    epstri =1.0d-10
    epsiter = 1.0d-14
    epsmag = 1.0d-21 !gradient square for magnetic axis
    maxiter = 50
    nt2 = 64
    ksamp = 8
    ksamp2 = 256
    nsub = 8
    kcheb = 16

    !If file_bf.eq. 'Solovev', use parameter
    csol=1.0
    reps=0.32d0
    rkappa=1.0 !1.7
    delta=0.0 !0.33
    R0=1.0d0
    Z0=0.0d0
    rndm = -1.0d0 !initially negative value

    nREF= 65  ! number of R grid for EFIT output format
    nZEF= 65  ! number of Z grid for EFIT output format

    q0=1.0d0
    F0=-1.0d0
    Fedge=-1.0d0
    p0=-1.0d0
    pin=2.0d0
    pout=1.0d0
    FF0=-1.0d0
    FFin=2.0d0
    FFout=1.0d0
    qfac=1.0d0
    qpin=2.0d0
    qpout=1.0d0
    Ffac =-0.1d0
    Fpin=2.0d0
    Fpout=1.0d0
    Jpar0 = 2.0d0 ! included mu0 (i.e. mu0*jpar)
    Jpin=2.0d0
    Jpout=1.0d0
    Jomax=0.0d0
    Joloc=0.7d0
    Jow=0.1d0
    nchq=10
    nflx=nchq
    psiflx=0.0d0
    wffprim=0.0d0
    epsmaxdist=1e-3
    maxdpsi = 1.0d0
    psiB = 0.0d0
    fpolsgn = -1.0d0
    iiterl = 1000
    iiterq = 1000
    lambdaex = -10.0d0
    torcur = 1.0d0
    Vloop0 = 1.0d0 !loop voltage [V]
    Vloop = -100.0d0 !will be overlapped

    ptf0=0.0d0
    ptfin=1.0d0 
    ptfout=1.0d0 
    ptfmax=0.0d0
    ptfloc=0.5d0 
    ptfw=0.1d0
    pgmax=0.0d0
    pgloc=0.5d0 
    pgw=0.1d0

    mach=0.0d0
    isymud = 1 !up-down symmetric tokamak :1
    iptype = 1
    iftype = 1
    ibtype = 1
    itftype= 0 
    iscale = 2
    isplch = 0 ! 0-use spline interpolation on psi of plasma domain from 1d-table for iptype=2
               ! 1-use spline interpolation on psi of nchq chebyshev point
               !   ,then use chebyshev interpolation on psi of plasma domain
    iptable = 0! 0- the table in "file_prof" gives dpdPsi in terms of psi (normalized pol.flx) in MKS unit
               ! 1- the table in "file_prof" gives pressure in terms of psi in MKS unit
               ! 2- the table in "file_prof" gives pressure in terms of rho in MKS unit
               ! 3- the table in "file_prof" gives pressure in terms of R in MKS unit

    irho  = 1  ! normalized radial coordintae rho=r/a 
               ! 0- r is the half-radius between the outer-midplane and the inner-midplane
               ! 1- r is the minor radius on outer-midplane from the magnetic axis
               ! 2- r is the sqrt of the poloidal flux (sqrt(Psi))
               ! 3- r is the sqrt of the toroidal flux (sqrt(phi))

    ifpol = 0  ! 0-use F0 at the center, 1-use Fedge at the edge for F profile
    iqtype = 1 ! 1-proifle, 2-table
    ijtype = 1 ! 1-proifle, 2-table
    iconvq = 2 ! 0-fixed lambda, 1-relax lambda, 2-method 2, 2-method 3
    iconvj = 0 
    maxiterF = 1
    relaxiter = 30 
    maxdpsiq = 1.0d-4
    maxdpsij = 1.0d-4

    iprintmap = 0
    iprintsol = 0
    iprintsoldisk = 0
    iprintcon = 1
    iprintQS = 1
    ifitMil = 1
    itftype=0

    itrinity =0   !Evaluation of local metrics for trinity inputs
    iogyropsi=1

    ntpsi = 8     !number of flux surface for trinity
    nttin = 12    !number of theta grid in a flux suface of trinity

    istability=0  !Anaysis for stability
    iLag = 0      !Evaluation of Lagrange (Grad-Hirshman variational form)
    ibscur =0     !Evaluation of boostrap current
    ijBSmodel =2  !Boostrap current model 1: Hirshman, 2: Sauter
    ifindmagaxis = 1 !0: use the maximum on grid as magnetic axis, 1: interpolate to find magnetic axis
    iprintEFIT = 0
    dpsistartmag =1.0d-3  ! finding magnetic axis off grid starts when dpsi is below dpsistartmag
    dpsistartnchq =1.0d-4  ! increase nchq over iterations when dpsi is below dpsistartnchq
    isolver=0   !unit disk PDE solver type: 
                !0-Poisson (laplacian)u=f
                !1-Helmholtz (laplacian-alpha*2)u=f, 
                !2-Helmholtz (laplacian-alpha*2)u=f
 
    iqconstraint =1  ! The additional constraint of q solver
                     ! 0: coneserved magnetic flux (lambda=lambdaex). iconvq=2 is recommended
                     ! 1: conserved magnetic field at the edge (Fedge is conserved). iconvq=1 is required


    iguessRmap0 = 1  ! guess the mapping center to be close the magnetic axis and reduce the number of mappings
                     ! 0: fix mapping center by R0
                     ! 1: shift mapping center by pressure term (shafranov shift)

    Halpha = 1.0d0 !alpha coefficeints in Helmholtz equation 
!   choose isolve =0 and alpha=0 for Poisson solver (default)
!   choose isolve =1 and alpha>0 for convergence in case of low current at core
 
   
    nchy = 16 ! Picthangle velocity grid to evaluate boostrap current
    npsi = 31 ! Number of flux surface for EFIT ouput format
    filein ="ecom.in"
    INQUIRE(FILE=filein, EXIST=file_exists)
    if(file_exists) then
       call get_namelists(filein)
    end if
    
    nr=nsub*kcheb
    ntot = nr*nt2


    call allocarray2
    call allocarray3
    call allocarray4

!    if ((ibscur.eq.1) .or. (ijtype.ge.4)) then
       call allocarrayBScur
!    end if
       if (itrinity.eq.1)  call allocarrayTRINITY
  end subroutine initarray
  subroutine get_namelists(filein)
 
    implicit none
!    logical :: list, accelx=.false., accelv=.false.
!    logical:: debug=.false.
    integer :: i,j
    real * 8 :: sflx, dRShaf
    character (20) :: filein
     
    namelist /GSparameter/ iecom, iptype, iftype, ibtype,file_prof, file_bc, file_efit, file_tflow,&
          nt1, nt2, kLag, ksamp, ksamp2, nsub, kcheb,&
         eps, epsiter, epsLag, epsLag2, maxiter, epsmaxdist,nchq, nflx, psiflx, &
        reps, rkappa, delta, R0, Rmid, Z0, q0, F0, Fedge, p0, pin, pout, FF0, FFin, FFout, &
        ptf0, ptfin, ptfout, ptfmax, ptfloc, ptfw, mach, isymud, ifpol, torcur, lambdaex,&
        iprintsol,iprintsoldisk, iprintcon, iprintQS, ifitMil, iscale, iprintmap, &
         itftype, qfac, qpin, qpout, Ffac, Fpin, Fpout, iqtype, file_qprof, &
         jpar0, jpin, jpout, ijtype, file_jprof, iconvq, iconvj, maxdpsiq, relaxiter, &
         maxdpsij, Jomax, Joloc, Jow, maxiterF, istability, ibscur, nchy, ijBSmodel, &
         Vloop0, Vloop, ifindmagaxis, epsmag, dpsistartmag, iprintEFIT, &
         dense0, denseb, rn1de, rn2de, te0, teb, rn1te, rn2te, ti0, tib, rn1ti, rn2ti, npsi, &
         itrinity, ntpsi, nttin, dpsistartnchq, isplch, iptable, irho, iqconstraint, iguessRmap0, &
         pgmax, pgloc, pgw, iremapon, iogyropsi

    ! NAMELIST characteristics
	! need the delim, else some implementations will not surround
	! character strings with delimiters
	! recl limits the I/O to 80 character lines
!	open(6, recl=80, delim='APOSTROPHE')
        open(8,file=filein, status='OLD', recl=80, delim='APOSTROPHE')
	read(8,nml=GSparameter)
! using the KEYWORD feature of routine
! calls adds clarity
!	call diffstuff(fullout=.TRUE.)
!        close(6)
        close(8)

        select case (iptype)
        case (0)
           !              ibtype = 0
           fpolsgn= -1.0d0 !F0/dabs(F0) 
        case (1,2)
           torcur = torcur*1.0d6  !from [MA] to [A] 
           fpolsgn= -1.0d0
           
           if ((iptable.eq.2).and.(isolver.eq.0)) then  !If r/a instead of psi is used for input profile,
                                      !helmoltz solver is recommended to make sure the negative right hand side 
              isolver=1
              halpha=5.0d0
           end if
        end select
           
        if (psiflx(2).eq.0.0d0) then
           do i=1,nflx
              sflx=real(i-1)/(nflx-1)
              psiflx(i)=sflx**2
              
           end do
        end if
        if (Vloop(1).eq.-100.0d0) then
           Vloop(1:nchq) =Vloop0*1.0d0
        end if

        select case (ibtype)
        case (0,2)
           Rmap0=R0
           Zmap0=Z0
           Rmid = R0 
           Rmaxpsi=R0 
           Zmaxpsi=Z0
        case (1)
           dRShaf=reps**2*dabs(p0)*R0*0.25  !estimate shafranov shit
           if (iguessRmap0.eq.0) then
              Rmap0=R0        ! fix the mapping center by R0
           else if (iguessRmap0.eq.1) then
              Rmap0=R0+dRShaf ! guess it using the magnetic axis shift
           end if
           Zmap0=Z0
           Rmid = R0 
           Rmaxpsi=R0 
           Zmaxpsi=Z0
        case (3) 
           isymud = 0 !assumed up-down asymmetric boundary from EFIT 
           Rmap0=R0
           Zmap0=Z0
           Rmid = R0 
           Rmaxpsi=R0 
           Zmaxpsi=Z0
        end select

        dRmap0=3.0d-3*R0*reps ! To find magnetic axis, put small distance from the center of unit disk to magnetic axis 
        if ((ifindmagaxis.eq.1).and.(iguessRmap0.ne.0)) then
           Rmap0=Rmap0-dRmap0
        end if

        write(*,*) 'psiflx',psiflx(1:nflx)
        if ((iscale.eq.5) .or.(iscale.eq.6)) then
           iconvj=1
        end if
        if (isolver.eq.0) Halpha=0.0d0 
      end subroutine get_namelists

    subroutine write_namelists
 
    implicit none

    
    write(fout,*) '**********< namelist >********' 
    write(fout,*) 'iecom=',iecom
    write(fout,*) 'iptype=',iptype
    write(fout,*) 'iftype=',iftype
    write(fout,*) 'iqtype=',iqtype 
    write(fout,*) 'ijtype=',ijtype 
    write(fout,*) 'ibtype=',ibtype 
    write(fout,*) 'itftype=',itftype
    write(fout,*) 'iptable=',iptable
    write(fout,*) 'irho=',irho
    write(fout,*) 'iscale=',iscale 
    write(fout,*) 'iconvq',iconvq
    write(fout,*) 'iconvj',iconvj
!    write(fout,*) 'isolver',isolver
    write(fout,*) 'file_prof=',file_prof
    write(fout,*) 'file_qprof=',file_qprof
    write(fout,*) 'file_jprof=',file_jprof
    write(fout,*) 'file_bc=',file_bc
    write(fout,*) 'file_efit=',file_efit
    write(fout,*) 'file_tflow=',file_tflow
    write(fout,*) 'nt1=',nt1
    write(fout,*) 'nt2=',nt2
    write(fout,*) 'nsub=',nsub
    write(fout,*) 'kcheb=',kcheb
    write(fout,*) 'nchq=',nchq
    write(fout,*) 'nflx=',nflx
    write(fout,*) 'kLag=',kLag
    write(fout,*) 'ksamp=',ksamp
    write(fout,*) 'ksamp2=',ksamp2
    write(fout,*) 'eps=',eps
    write(fout,*) 'epsiter=',epsiter
    write(fout,*) 'maxiter=',maxiter
    write(fout,*) 'epsmaxdist=',epsmaxdist
    write(fout,*) 'maxdpsiq', maxdpsiq
    write(fout,*) 'maxdpsij', maxdpsij
    write(fout,*) 'relaxiter', relaxiter
    write(fout,*) 'reps=',reps
    write(fout,*) 'rkappa=',rkappa
    write(fout,*) 'delta=',delta
    write(fout,*) 'R0=',R0
    write(fout,*) 'Z0=',Z0
    write(fout,*) 'q0=',q0
    write(fout,*) 'F0=',F0
    write(fout,*) 'p0=',p0
    write(fout,*) 'pin=',pin
    write(fout,*) 'pout=',pout
    write(fout,*) 'FF0=',FF0
    write(fout,*) 'FFin=',FFin
    write(fout,*) 'FFout=',FFout
    write(fout,*) 'jpar0=',jpar0
    write(fout,*) 'jpin=',jpin
    write(fout,*) 'jpout=',jpout
    write(fout,*) 'jomax=',Jomax
    write(fout,*) 'joloc=', Joloc
    write(fout,*) 'jow=', Jow 
    write(fout,*) 'qfac=',qfac
    write(fout,*) 'qpin=',qpin
    write(fout,*) 'qpout=',qpout
    write(fout,*) 'ptf0=',ptf0
    write(fout,*) 'ptfin=',ptfin
    write(fout,*) 'ptfout=',ptfout
    write(fout,*) 'mach=',mach
    write(fout,*) 'isymud=',isymud
    write(fout,*) 'iprintsol',iprintsol
    write(fout,*) 'iprintsoldisk',iprintsoldisk
    write(fout,*) 'iprintcon', iprintcon
    write(fout,*) 'iprintQS', iprintQS
    write(fout,*) 'ifitMil', ifitMil
    write(fout,*) 'iprintmap',iprintmap
    write(fout,*) 'iprintEFIT',iprintEFIT
    write(fout,*) 'istability',istability
    write(fout,*) 'ibscur', ibscur
    write(fout,*) 'nchy', nchy
    write(fout,*) 'ijBSmodel', ijBSmodel
    write(fout,*) 'iogyropsi', iogyropsi

  end subroutine write_namelists

  subroutine allocarray1

    implicit none

    allocate (Rin(nin+1), Zin(nin+1), Tin(nin+1))
    allocate (Rt1(nt1+1), Zt1(nt1+1), T1(nt1+1), T1_mil(nt1+1))  
    allocate (dRt1(nt1+1), dZt1(nt1+1),T1_mil2(nt1+1))  
    allocate (dWdZ1(nt1+1), dWdZZ1(nt1+1)) 
    allocate (w1(nt1), wsaved(nt1*nt1), w7saved(nt1*nt1))

    allocate (dfdzbf(nt1+1),dfdzzbf(nt1+1))
    allocate (fbf(nt1*ksamp+1), zbf(nt1+1), zmbf(nt1*ksamp+1))
    allocate (dwdz1m(nt1*ksamp+1),dwdzz1m(nt1*ksamp+1))
    allocate (w7(2*nt1*nt1+10*nt1+10000), w8(10*nt1*ksamp), w9(nt1,nt1)) 
    Rt1 = 0.0; Zt1 = 0.0; T1 = 0.0; T1_mil=0.0; T1_mil2=0.0
    dRt1 = 0.0; dZt1 = 0.0; w1= 0.0

  end subroutine allocarray1
  
  subroutine allocarrayProf 

     implicit none

     allocate (psiEF(npsi),psiiEF(npsi), pprimEF(npsi), ffprimEF(npsi))
     allocate (bpprim(npsi),cpprim(npsi),dpprim(npsi))
     allocate (bpres(npsi),cpres(npsi),dpres(npsi))
     allocate (bfpol(npsi),cfpol(npsi),dfpol(npsi))
     allocate (rhoEF(npsi), presEF(npsi), fpolEF(npsi))
     allocate (bffprim(npsi),cffprim(npsi),dffprim(npsi))
     allocate (btflow(npsi),ctflow(npsi),dtflow(npsi))
   end subroutine allocarrayProf

  subroutine allocarrayQProf

     implicit none

     allocate (psiiq(nqpsi),qpsiin(nqpsi))
     allocate (bqpsi(nqpsi),cqpsi(nqpsi),dqpsi(nqpsi))
   end subroutine allocarrayQProf

   subroutine allocarrayJProf

     implicit none

     allocate (psiij(njpsi),jpsiin(njpsi))
     allocate (bjpsi(njpsi),cjpsi(njpsi),djpsi(njpsi))
   end subroutine allocarrayJProf

   subroutine allocarrayBScur

     implicit none

     allocate (psiichy(nchy),cftmy(nchy*nchy),spybx(nchy*nchy))
     allocate (spyxb(nchy*nchy),spydef(nchy))
     allocate (fpass(nchq),etai(nchq),Zi(nchq),jBSch(nchq))
     allocate (nue_star(nchq), nui_star(nchq))
     allocate (ne(nchq), Te(nchq), jOhmch(nchq)) !Vloop(nchq)
     etai(1:nchq)=1.0d0
     Zi(1:nchq)=1.0d0
     ne(1:nchq)=0.1d0  ! [10^20 m^-3]
     Te(1:nchq)=5.0d2  ! [eV]
     !Vloop(1:nchq) =Vloop0*1.0d0 ![Volt]
     nue_star(1:nchq)=0.1d0
     nui_star(1:nchq)=nue_star(1:nchq)/60.0d0

   end subroutine allocarrayBScur

   subroutine allocarrayEFIT

     implicit none

     allocate (psiEF(npsi),psiiEF(npsi), pprimEF(npsi), ffprimEF(npsi))
     allocate (fpolEF(npsi), qpsiEF(npsi))
     allocate (bpprim(npsi),cpprim(npsi),dpprim(npsi))
     allocate (bffprim(npsi),cffprim(npsi),dffprim(npsi))
     allocate (REF(nREF), ZEF(nZEF))

   end subroutine allocarrayEFIT

   subroutine allocarrayTrinity
     implicit none
     integer :: i,j

     allocate (tpsi(ntpsi),ttin(nttin), trho(ntpsi),dpsidrho(ntpsi))
     
     do i=1,ntpsi
!        tpsi(i)=1.0d0-real(i)/(ntpsi+1)
!        tpsi(i)=1.0d0-real(2*i-1)/(2*ntpsi)
        trho(i)=1.0d0/(ntpsi-1)*(i-1)
!        trho(i)=1.0d0/(ntpsi)*real(i) !JPL:remove psi=0
     end do
     do j=1,nttin
        ttin(j)=8.0d0*datan(1.0d0)/nttin*(j-1)             
     end do

     allocate (gradpar(ntpsi*nttin),btot(ntpsi*nttin),bpol(ntpsi*nttin))
     allocate (Rtrin(ntpsi*nttin),Ztrin(ntpsi*nttin),gbdrift(ntpsi*nttin), gbdrift0(ntpsi*nttin))
     allocate (cvdrift(ntpsi*nttin), cvdrift0(ntpsi*nttin),shatlocal1(ntpsi*nttin))
     allocate (shatlocal2(ntpsi*nttin),shatlocal3(ntpsi*nttin))

   end subroutine allocarrayTrinity

  subroutine allocarray2

    implicit none

    allocate (Rt2(nt2+1), Zt2(nt2+1), T2(nt2+1))  
    allocate (dRt2(nt2+1), dZt2(nt2+1))
    allocate (zk(nt2+1), dzdw2k(nt2+1), dzdww2k(nt2+1))
    allocate (dWdZ2(nt2+1))
    allocate (dWdZZ2(nt2+1))
  
end subroutine allocarray2

  subroutine allocarray3

    implicit none

    save

    allocate (rnd(nr), tnd(nt2))  
    allocate (Rt3(ntot),Zt3(ntot))  
    allocate (dzdw3(ntot)) 
    allocate (dzdww3(ntot)) 
    allocate (dwdz3(ntot)) 
    allocate (dwdzz3(ntot))
    allocate (blength(nsub),bnodes(nsub+1))
 
  end subroutine allocarray3

  subroutine allocarray4 

    implicit none

    integer ::  ifcoeff,  iu ,iur, iurr, izc, izcp, itot
    allocate (cftmsub(kcheb*kcheb))
    allocate (psiichq(nchq), phiichq(nchq), cftmq(nchq*nchq),spbx(nchq*nchq),spxb(nchq*nchq))
    allocate (spdef(nchq), qpsich(nchq), jpsich(nchq), pprimch(nchq),ffprimch(nchq))  
    allocate (qpsich1(nchq),fpolcon2(nchq), ffprimpre(nchq), ffprimpre1(nchq))
    allocate (ffprimpre2(nchq) , ffprimgoal(nchq))
    allocate (fpolch(nchq),shat(nchq),chcoeffpp(nchq),chcoeffpp0(nchq),chcoeffff(nchq)) 
    allocate (curpar(nchq), ptfprimch(nchq), ptflowch(nchq),chcoeffprho(nchq),chcoefff2rho(nchq))
    allocate (chcoeffp(nchq), presch(nchq), rhoch(nchq),dpsiidrhoch(nchq), chcoefftf(nchq))

    nchq0=nchq
    allocate (qpsich0(nchq0), jpsich0(nchq0),fpolch0(nchq0), pprimch0(nchq0))
    allocate (psiichq0(nchq0), cftmq0(nchq0*nchq0),spbx0(nchq0*nchq0))
    allocate (spxb0(nchq0*nchq0),spdef0(nchq0))

!    allocate (spdefs(nchqs), qpsichs(nchqs), pprimchs(nchqs), qpsich1s(nchqs))  
!    allocate (fpolcon2s(nchqs), ffprimpres(nchqs), ffprimpre1s(nchqs))
!    allocate (ffprimpre2s(nchqs) , ffprimgoals(nchqs))
!    allocate (fpolchs(nchqs), dZdtmidos(nchqs), dZdtmidis(nchqs)) 

    allocate (g(nt2), f(ntot), psii(ntot), psiex(ntot))
    allocate (dpsidr(ntot), dpsidrr(ntot), u(ntot), ur(ntot), uth(ntot), utt(ntot))
    allocate (urr(ntot), urt(ntot), ucoeff(ntot), urcoeff(ntot), urrcoeff(ntot))


      !nr = nsub*k
      !ntot = nr*nth
      ifcoeff = 1 
      iu = ifcoeff + ntot
      iur = iu + ntot
      iurr = iur + ntot
      izc = iurr + ntot
      izcp = izc + nt2
      itot = izcp + nt2

      allocate(zw(itot))
      allocate(gcoeff(nt2))
      allocate(pt(nr+10))
      allocate(qt(nr+10))
      allocate(fr(nr+10))
      allocate(fi(nr+10))
      allocate(un(nr+10))
      allocate(unr(nr+10))
 
    fpolcon2 = 0.0d0
    ffprimpre = 0.0d0 
    ffprimpre1 = 0.0d0
    ffprimpre2 = 0.0d0 
    ffprimgoal = 0.0d0 
  end subroutine allocarray4

  subroutine deallocarray

    implicit none

    deallocate (Rin, Zin, Tin,  Rt1, Zt1, T1, dRt1, dZt1, dwdz1)
    deallocate (Rt2, Zt2, T2, dRt2, dZt2, dwdz2)
    deallocate (rnd, tnd, Rt3, Zt3, dzdw3,T1_mil, T1_mil2)
    deallocate (dpsidr, dpsidrr, u, ur, uth, utt)
    deallocate (urr, urt, ucoeff, urcoeff, urrcoeff)
    deallocate (g, f, psii, psiex) 
    deallocate (zw, gcoeff, pt, qt, fr, fi, un, unr) 
    deallocate (dfdzbf,dfdzzbf, fbf,zbf,zmbf, dwdz1m, dwdzz1m )
    deallocate (w7, w8, w9)

  end subroutine deallocarray
end module arrays
