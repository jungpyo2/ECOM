
      program ecom
c----------------------------------------------------------
c
c             EEEEE    CCC     OOO     MM   MM
c             E       C       O   O    M M M M
c             EEEEE  C        O   O    M  M  M
c             E       C       O   O    M  M  M
c             EEEEE    CCC     OOO     M  M  M 
c
c-----------------------------------------------------------
c       ECOM (Equilibria via COnformal Mapping) :
c       A program to solve Grad-Shafnov equation 
c       for a static magnetic equilibrium of a toroidal axisymmetric system
c       with a fixed boundary using a conformal mapping to unit disk.
c
c       date: Mar.1. 2014
c       version : 1.0
c       contact: Jungpyo Lee (jungpyo@cims.nyu.edu)
c                or Antoine Cerfon (cerfon@cims.nyu.edu)  
c------------------------------------------------------------
c
c     Files:
c     ecom.f--arrays.f90 (Initialization, Global variables,...)
c            -readbf.f (Reading fixed boundary)
c            -cmap.f (conformal mapping)
c            -solve.f (Grad-Sharanov solver in an unit disk)
c            -postproc.f (Contour integrals, Error analysis,...)
c
c     readbf.f - anaresa.f (Calculation of arc-length)
c              - interpol.f (Lagrange or Spline interpolalation)
c              - cqrsolve.f (QR decomposition)  
c
c     cmap.f - interpol.f (Lagrange interpolalation)
c            - durdecplus.f (Discrete Fast Fourier Transform)
c            - quaplot.f (Plotting for Conformal mapping)
c            - cqrsolve.f (QR decomposition) 
c            < These conformal mapping routines are written
c              by Micheal Oneil>
c
c     solvef.f (iecom=0) predfined FF' (poloidal current) profile 
c     solvej.f (iecom=1) predefined j_par (parallel current) profile
c     solveq.f (iecom=2) predfined q  (safety factor) profile
c             - elliptic.f (Poisson solver in an unit disk )
c             - bvpadapfix.f (ODE solver using Green function)
c             < These Poisson solver routines are written 
c               by Jun-Yup Lee and Leslie Greengard>
c
c     postproc.f - interpol.f (Trigonometric or Lagrange interpolalation)
c                - durdecplus.f (Discrete Fast Fourier Transform)
c
c---------------------------------------------------------------
c
c
c    Namelist: declared in arrays.f90 and read from "ecom.in"
c    iecom  - Predefined profile along with the pressure profile in G-S eq.
c           - 0 (default): Poloidal current profile FdF/dpsi
c           - 1 : Parallel current profile J_par
c           - 2 : Safety factor profile q  
c    iptype - type of pressure (p') profile 
c           - 0: p'=constant for Solovev's analytic solution
c           - 1 (default): mu0*R0**2*p'=p0*(1-psin**pin)**pout 
c                where psin is normalized flux to [0,1] (0:axis, 1:edge)
c           - 2: discrete values of p' or p in terms of psin or rho
c                are given by a table in 'file_prof'
c           - 3: discrete values of p' in terms of psin
c                are given by EFIT output 'file_efit'
c    iftype - type of poloidal current (FF') profile (activated for iecom=0) 
c           - 0: FF'=0 for Solovev's analytic solution
c           - 1 (default): FF'=ff0*(1-psin**ffin)**ffout 
c                where psin is normalized flux to [0,1] (0:axis, 1:edge)
c           - 2: discrete values of FF' or F in terms of psin or rho
c                are given by a table in 'file_prof'
c           - 3: discrete values of FF' in terms of psin
c                are given by EFIT output 'file_efit'
c    ijtype - type of parallel current (jpar) profile (activated for iecom=1) 
c           - 1 (default): jpar=jpar0*(1-psin**jpin)**jpout 
c                where psin is normalized flux to [0,1] (0:axis, 1:edge)
c           - 2: discrete values of jpar in terms of psin 
c                are given by a table in 'file_jprof'
c           - 3: jpar as evaluated from ohmic and boostrap current model
c                 at eachc iteration
c           - 4: jpar=jpar0*((1-psin**jpin)**jpout
c                            +jomax*exp(-((psin-joloc)/jow)**2))  
c    iqtype - type of safety factor (q) profile (activated for iecom=2) 
c           - 1 (default): q=q0*(1+qfac*psin**qin)**qout 
c                where psin is normalized flux to [0,1] (0:axis, 1:edge)
c           - 2: discrete values of jpar in terms of psin 
c                are given by a table in 'file_qprof' 
c    ibtype - type of the fixed boundary contour (last closed flux surface) 
c           - 0: For Solovev's analytic solution (only for iptype=0),
c                the contour is given by the equations with the parameter, t 
c                R=R0*sqrt(1+2*reps*cos(t))
c                Z=Z0+reps*rkappa*R0**2/R*sin(t)
c            -1 (default): For Miller local equilibrium,
c                the contour is given by the equations with the parameter, t 
c                R=R0*(1+reps*cos(t+asin(delta)*sin(t)))
c                Z=Z0+R0*reps*rkappa*sin(t)
c            -2: (R,Z) are given by a table in 'file_bc' 
c            -3: (R,Z) are given by EFIT output file 'file_efit'
c    itftype - type of the toroidal flow pressure (ptorflow=1/2*rho*R0**2*omega**2) 
c              Modified pressure in G-S equation is 
c              pmod(R,psi)=p(psi)*exp(ptorflow(psi)/p(psi)*(R**2-R0**2)/R0**2),  
c              which is valid for negligible poloidal flow and strong toroidal flow 
c              (comparable to thermal flow).
c            - 0 (default): no toroidal flow
c            - 1: ptorflow=mach**2*pressure profile
c            - 2: ptorflow=ptf0*(1-psin**ptfin)**ptfout  
c            - 3: ptorflow is given by a table in 'file_tflow'
c    iptable - type of pressure profile table in "file_prof" 
c            - 0 (default): dp/dPsi in terms of psin in MKS unit
c            - 1: p in terms of psin in MKS unit
c            - 2: p in terms of rho in MKS unit
c    irho    - type of normalized radial coordintae rho 
c            - 0 (default) : half-radius between the outer-midplane and the inner-midplane
c                 rho=(Rout(psin)-Rin(psin))/(Rout(psin=1)-Rin(psin=1))
c            - 1: rho=(Rout(psin)-R0)/(Rout(psin=1)-R0)
c            - 2: rho=sqrt of the normalized polodail flux (sqrt(Psi))
c            - 3: rho=sqrt of the normalized toroidal flux (sqrt(phi))
c    iscale - control parameter to scale the solutions by normalizing
c             pressure and toroidal magnetic field profiles (only for iptype=1 and iptype=2)
c           - Among p0, FF0, torcur, F0 and q0 in the namelist, some are selected and the others are ignored
c           - 0: Absolute values of p0, FF0 are used to solve psi, and F0 is used for contour integral.
c                torcur and q0 are ignored  
c           - 1: Absolute values of p0, FF0 are used to solve psi, and q0 is used for contour integral.
c                torcur and F0 are ignored  
c           - 2: torcur is used to normalize p0 and FF0, so the ratio p0/FF0 and torcur is considered to solve psi. 
c                F0 is used for contour integrals. q0 is ignored. (default)
c           - 3: torcur is used to normalize p0 and FF0, so the ratio p0/FF0 and torcur is considered to solve psi. 
c                q0 is used for contour integrals. F0 is ignored. 
c           - 4: Both F0 and q0 are used to normalize p0 and FF0, so the ratio p0/FF0, F0, q0 are considered 
c                to solve psi and find contour integrals. torcur is ignored. 
c    ifpol  - 0: use F0 at the center (psin=0) to reconstruct F profile from FF' 
c           - 1: use Fedge at the edge (psin=1) to reconstruct F profile from FF' 
c
c    file_prof - File to give 1-D profiles of (psin, p', FF')
c    file_bc   - File to give (R0, Z0) and table for (R,Z) of the boundary
c    file_efit - File to give EFIT input and output information
c
c    nt1 - number of grid points in the boundary for conformal mapping
c    nt2 - number of theta points in the unit disk (theta grid of solutions)
c    nsub - number of radial subgroups in the unit disk
c    kcheb - number of chebyshev points in a radial subgroup
c            (radial grid of solutions is nr=nsub*kcheb)
c    nchq - number of flux surface in chebyshev grid of psin
c    nflx - number of flux surface of specific psi. Interpolated by nchq data 
c    psiflx - 1-d array of specific normalized psi of the flux surface
c                  
c
c    kLag - order of Lagrange interpolation 
c    ksamp - factor of oversampling for FFT pading for conformal mapping
c    ksamp2 - factor of oversampling for FFT pading for contour integrals
c   
c    R0  - R (Major radius) at the center point of the boundary.  
c    Z0  - Z at the center point of the boundary.
c          (R0,Z0) is used for the center point of conformal maping at the first iteration
c    q0  - safety factor at the magnetic axis
c    F0  - constant toroidal flux at the magnetic axis for Solovev solution, iptype=0  
c    p0  - pressure gradient at the magnetic axis multipled by permeability and R0**2 for iptype=1
c    pin - inner exponent of the pressure gradient profile for iptype=1
c    pout - outer exponent of the pressure gradient profile for iptype=1
c    FF0  - FF' at the magnetic axis for iptype=1
c    FFin - inner exponent of the FF' profile for iptype=1
c    Ffout - outer exponent of the FF' profile for iptype=1
c    torcur - total toroidal current in [MA] for iptype=1 and iptype=2
c    reps - ratio of minor radius to major radius (R0) for ibtype=0 and ibtype=1
c    rkappa - elongation of the boundary for ibtype=0 and ibtype=1
c    delta  - triangularity of the boundary for ibtype=1
c   
c    epsiter - small real constant to determine the convergence of iteration
c    maxiter - maximum number of iterations
c    epsmaxdist - allowed maximum distance between (R0,Z0) and (R,Z) of magnetic axis (minimum psi).
c                 If the distance is larger than epsmaxdist,(R,Z) of magnetic axis is assigned 
c                 to be the center of new conformal mapping during iterations.
c 
c    isymud - index for up-down symmetry of the solution and boundary in z=0 
c           - 0: up-down assymetric
c           - 1: up-down symmetric 
c             isymud=1 is used to improve accuracy and computation speed of solver
c

c
c    mach  - ratio of toroidal flow to thermal flow for itftype=0
c    ptf0 - torodial flow pressure at the magnetic axis for itftype=1  
c    ptfin - inner exponent of the toroidal flow pressure profile for itftype=1
c    ptfout - outer exponent of the toroidal flow pressure profile for itftype=1
c    file_tflow - File to give 1-D profiles of toroidal flow pressure (ptorflow) for itftype=2
c
c    iprintmap     - Print Ascii file for conformal mapping information (0: no print, 1: print)
c    iprintsol     - Print Ascii file for psi and its derivatives (0: no print, 1: print)
c    iprintsoldisk - Print Ascii file for psi and its derivatives in unit disk (0: no print, 1: print)
c    iprintcon     - Print Ascii file for Contours of constant psi (0: no print, 1: print)
c    iprintQS      - Print Ascii file for q (safety factor) and s (magnetic shear) profiles (0: no print, 1: print)
c    iprintEFIT    - Print Ascii file in EFIT g file format (0: no print, 1: print)
c    npsi - Number of flux surface for EFIT ouput format
c
c    ifitMil   - fitting the psi contour of nchq chebyshev grid 
c                using local Miller equilibrium parameter (0: no fitting, 1: fitting)  
c    istability - Control parameter for stability analysis (0: off, 1: on)
c    ibscur     - Control parameter for evaluating boostrap current (0: off, 1: on)
c    ijBSmodel  - Boostrap current model 1: Hirshman, 2 (default): Sauter 
c    ifindmagaxis - Control paramter for finding maximum psi (magnetic axis) (0: on grod, 1 (default) off grid)
c    itrinity   - Control parameter for evaluating local metrics for Trinity code inputs
c    ntpsi      - number of flux surface for trinity
c    nttin      - number of theta grid in a flux suface of trinity

c
c---------------------------------------------------------------

      use arrays, only: iecom, initarray, deallocarray
      use arrays, only: iiter, iremap, epsiter, maxiter, maxdpsi, fout
      use arrays, only: iscale, iLag, iptable


      call initarray
      call initfiles

      call readbf

      call runtime('   reading boundary   ')

      iiter=0
     
      do while ((maxdpsi > epsiter) .and. (iiter < maxiter)) 

         iiter = iiter +1
         write(fout,*) 'iteration: ',iiter,', new mapping ',iremap

c        conformal mapping during iterations if iremap=1
         if (iremap.eq.1) then
            call cmapbf
            call runtime(' started backward mapping ')
            call cmapin
            call runtime(' finished backward mapping ')
         end if
        

c        Poisson solver for the unit disk
         select case (iecom)

         case (0) ! Predefined poloidal current profile (dF/dpsi)
            call solvef
            if ((iremap.ne.1).and.((iscale.eq.5).or.(iscale.eq.6)
     1           .or.(iLag.eq.1).or.(iptable.ge.2))) then
               call postproc1
            end if  
 
         case (1) ! Predefined parallel current profile (j_par)
            call solvej
            if (iremap.ne.1) then !Find FF' by contour integrals
               call postproc1
            end if
        
         case (2) ! Predefined safety factor profile (q) and lambda
            call solveq
            if (iremap.ne.1) then !Find FF' by contour integrals
               call postproc1
            end if

         end select
      end do  
     
      write(*,*) 'Iteration finished '
      call runtime('   poisson solver    ')

c     Post-processing
      call postproc2

      write(*,*) 'Postprocess finished '
      call runtime('   postprocessing    ')

      call deallocarray
      call closefiles

      stop
      end   
