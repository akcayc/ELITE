
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ELITE (Edge Localized Instabilities in Tokamak Experiments)
!
! Authors:  P.B. Snyder, H.R. Wilson, R.L. Miller
!
! Copyright(c) 1999-2005, All rights reserved
!
! Code is under development, and the source code documentation and 
!  user's guide are presently in a preliminary form. The authors 
!  therefore request that users consult with the code authors, 
!  particularly before publishing ELITE results.   We also request 
!  that users report any errors in the documentation or code to the 
!  authors, and share any modifications made to the code or accompanying 
!  idl routines with the authors.  Please request permission from
!  P.B. Snyder or H.R. Wilson before distributing the code, as we wish 
!  to keep track of all users and changes to the code, and maintain a 
!  master version of the code to prevent unnecessary branching.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module elite_data
!  Philip B. Snyder 1/31/00

  implicit none

      real ppmult   ! multiplier for pprime, to allow variation of
                       ! p' on all surfaces
      real alphamult  ! multiplier for alpha on all surfaces, should not
            ! be used with ppmult
      real shearmult  ! multiplier for shear on all surfaces, should not
            ! be used with ppmult
      logical sa_fix ! set true to fix shear and alpha on outer surface
         ! (and same multiplier on inner surfaces), requires
         !   ppmult=alphamult=shearmult=1.
      real shear_fix ! fixed value for outer surface shear
      real alpha_fix ! fixed value for outer surface alpha
      logical analytic_peel ! if true, code outputs the analytic peeling
        ! criterion lhspeel-rhspeel (>0 for stability) rather than doing
        ! the actual stability calculation
      integer surfplot ! tag for outputting data for surface plots in fun2d
        ! if surfplot=0 not surface plotting data is written.  If
        ! surfplot>0 then surface plotting data is written with
        ! nxplt=surfplot*nxinterp+1 radial points 
      integer wsmodel ! tag for using various omega* models to assign
        !  the value of gamsq for the stability check
        !  default is zero.  wsmodel=1 uses 1/2 of the maximum
        !  omega_star_p model to determine gamsq
      integer ncrit_ws  ! critical n above which diamagnetic stabilization
        !  rolls over for wsmodel=10 and wsmodel=11.  default=13
        !  these rollover models are based on BOUT++ results

      real ion_mass ! ion mass relative to proton mass. used for rotation
           ! normalization and omega*

      real qref     !reference q surface
      real dmercier ! D_M=alpha*dmercier
                                ! >0 is stabilizing ?
      real aspect   !aspect ratio
!      real xwid     !half-width of fourier modelet
      real dx
      real dw
      real lamdist
      real q0
      integer verbose           !integer 0-10 indicating amount of text and
                                ! file output.  0=no text, essential files
                                ! defaul=4
      integer nmwinhalf         !total number of modes in mode window
                                ! is 2*nmwinhalf+1
      integer nn                !toroidal mode number
      integer mmax              !maximum poloidal mode number
      integer nm                !number of poloidal mode numbers
      integer nmodes            ! useful for length of vectors
                                !=nm for updownsym=.t., =2*nm for updownsym=.f.
      integer nmres             ! number of m's resonant in the plasma
      integer nxinterp          !number of x points at which matrices
                                !are evalutated; used for interpolation
                                !to other values
      logical dens              ! user supplied n,T profiles if .true.
      logical rotation          ! rotation not supported in this version
                                !  but keep full namelist for compatibility
      integer nmvac             ! number of vacuum modes, should be one for now
      integer nmlow             ! number of modes with m<resonant values,
                                !   added on so nm=nmvac+nmres+nmlow
      real psimin               ! minimum psi value requested in qref namelist
      integer nd
      integer ndist
      integer ns                !number of points on reference surface
      integer meshtype          !1 use dist variables, 2 is wilson's orig
                                ! for determining xx(i)
      integer mindex            !indexing parameter =1 for updownsym
                                ! =2 for .not.updownsym
      integer n1term            !if non-zero, additional 1/n correction
                                !terms are included
      logical autorun           !if true then last line of output is 
                                !dW for least stable or most unstable mode
                                !useful for scans using scripts
      logical usegamr2          !if true, then iterate to convergence using gamr2
                                !(the eigenvalue without forced hermitian symmetry)
                                ! in addition to gamr.  This option is useful for testing
                                ! the impact of equilibrium and other errors on the
                                ! growth rate
      logical updownsym         !true for updown symmetric equilibria
      logical vacuum            !true if dW_vacuum to be included
      logical funcal            !true if want to calculate eigenfunction as well
                                !as growth rate
      logical kinkterms         !true by default, set false to eliminate kink terms
      logical nowindow          !true if every value of x calculates every
                                !poloidal mode number
      logical splineeq          !1/02 true to spline au* & iu* arrays over
              ! the nxinterp grid rather than using linear interpolation
              ! use only when nxinterp is small!!!(<~20) or it will be slow
      real newtstep  ! parameter to set initial step in Newton iteration to
              ! calculate growth rate.  First step will multiply or divide
              ! gamsq by 1.+newtstep depending on sign of gamr
      real btnorm
      real rnorm
      real bug(10)            !spare parameters during program development
      integer mmin              !minimum poloidal mode number=mmax-nm+1
      integer mres              !m value for primary resonant surface
                                !so x=mres-n*q(psi)
      integer nedge             !# of modes non-zero at edge--effect of 
                                ! m window
      integer mwindow           ! if windowning is on, max size of window
                                !   if not, mwindow=nmodes
      real, dimension(:), allocatable :: rl,zl,bpl,uprimel,cosul,sinul
! 10/20/00 poloidal differentials of psi1
      real, dimension(:), allocatable :: dpsi1l,d2psi1l
!      real rl(ks)   !R on reference surface
!      real zl(ks)   !Z on reference surface
!      real bpl(ks)  !poloidal mag. field on ref. surf.
!      real uprimel(ks) !du/dl on ref. surf.
!      real cosul(ks) !cos(u) on ref. surf
!      real sinul(ks) !sin(u) on ref. surf.
      real, dimension(:), allocatable :: xinterp
!      real xinterp(kxinterp) !x values that au,aup,aupp defined at
!      real, dimension(:,:,:), allocatable :: au,aup,aupp,audd,aupdd,auppdd, &
!            aud,aupd,auppd,aund,aupnd,auppnd
! 1/02 replace au.. arrays with one large array, 
!   auarr(19,lwindow,lwindow,nxinterp) where elements 1-19 correspond
!   to the old array set in the order below
      real, dimension(:,:,:,:), allocatable :: auarr
! 1/02 auspline to hold spline coefficients for auarr
      real, dimension(:,:,:,:), allocatable :: auspline
!      real, dimension(:,:,:), allocatable :: aumk,au2mk,aum2k,au2m, &
!             au2k,aum,auk,aund
!      real, dimension(:,:,:), allocatable :: aupmk,aup2mk,aupm2k,aup2m, &
!             aup2k,aupm,aupk,aupnd
!      real, dimension(:,:,:), allocatable :: auppmk,aupp2mk,auppm2k
! 1/02 also replace iu arrays with one large array 1-9 correspond to below
      real, dimension(:,:,:,:), allocatable :: iuarr,iuspline
!      real, dimension(:,:,:), allocatable :: iuk,ium,iund
!      real, dimension(:,:,:), allocatable :: iupk,iupm,iupnd
!      real, dimension(:,:,:), allocatable :: iuppk,iuppm,iuppnd
!      real au(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um
!      real aup(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um'
!      real aupp(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um''
!      real audd(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um
!      real aupdd(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um'
!      real auppdd(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um''
!      real aud(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um
!      real aupd(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um'
!      real auppd(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um''
!      real aund(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um
!      real aupnd(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um'
!      real auppnd(kmodes,kmodes,kxinterp) 
!                                !ballon eq. matrix--coeff of um''
      real, dimension(:,:), allocatable :: asurfu,asurfup, &
                                           isurfu,isurfup
!      real agamma(kmodes,kmodes) 
!                                !surface eq. matrix--ffprime piece
!      real asurfu(kmodes,kmodes) 
!                                !surface eq. matrix--coeff of um
!      real asurfup(kmodes,kmodes) 
!                                !surface eq. matrix--coeff of um'
      real circumf !total circumf of ref. surf.
      real f        !f=B_tor*R
!      real ffprime  !ff' piece of toroidal current, use ffprime_eq&_adj
!      real pprime   !p'  PBS: actually 4 Pi p', use pprime_eq&pprime_adj   
      real pi       !3.14...
      real psisep   !psi value where b_pol -> 0 based on
                                ! linear expansion from ref. surf.
      real gamsq    !the current guess at the square of the growth rate
                    !the code will iterate on this if igam >/=1
      real gamscl   !scales gamsq to give it normalised to Alfven frequency,
                    !ie, to omega_A=B/sqrt(mu_0*rho)*R, with B and R an
                    !average on outer surface.
      real gamscl2  ! new factor using R averaged from in and out and
                    !  f=R*B_phi on outer surface to represent B
      integer igam  !if >/= 1, code will iterate to find the true growth
                    !rate-squared, gamsq...otherwise it solves for the 
                    !eigenvalues for the given gamsq; igam=no. iterations
      real gamim ! imaginary part of eigenvalue for rotational case
      real del      !del=mres-n*qedge
      real xmin     !minimum x=mres-n*q
      real xmax     !maximum x
      real qrefp    !dq/dpsi at the reference surface
      real vprime
      real volume
      real vpp
      real abishop
      real sigbishop
      real shbishop
      real rmajor
      real b2inv_edge !<B^-2> at edge defined in nustuff
                                !strictly speaking it is not dl/B_p average
                                !but rather dl/(R^2*B_p) average
                                !may want to change this some time
      character*24 cdatetime
      character*14 codver
!-----------------------------------------------------------------------
!c variables from original edgsa code
!-----------------------------------------------------------------------
      integer kx
      real gam
      real dm
      real qa
      integer m0,m1up,m1lo
      real, dimension(:,:), allocatable :: ffun
      real, dimension(:), allocatable :: aval
      real, dimension(:), allocatable :: xx
      real, dimension(:,:), allocatable :: pmat,qmat,smat
      real, dimension(:,:,:), allocatable :: aa,ainv
      real, dimension(:,:), allocatable :: bb
      real, dimension(:,:,:), allocatable :: fmat
      real, dimension(:,:), allocatable :: dw_vac
!      real dw_vac             !vacuum matrix from vac_edge
      integer offset_v
      integer nx
      real tau         !Te/Ti
      real te          !Te (eV)
      real etai        !ratio of n to T lenth scale (=etae)
      real ne          !electron density (=ni)
      real vloop       !loop voltage
      real bstrap      !bootstrap current <J.B>
      real indj        !inductive current <j.B>
      real zeff        !Zeff

      character runname*50  ! Name of this run.  Associated i/o files
                        ! will be named RUNNAME.in, etc.
      integer lrunname      ! length of runname
      integer outfile    ! file number for .out file
! flux surface quantities from EFIT read in from .eqdat file
      real, dimension(:), allocatable :: psiv_eq,pprime_eq,f_eq, &
           ffprime_eq,vprime_eq,q_eq,ppp_eq,ffpp_eq,qp_eq,qpp_eq
! ppp_eq is the equilibrium p'', ffpp_eq=(ff')'
! qp_eq is the equilibrium q', qpp_eq=q''    
      real, dimension(:), allocatable :: ne_eq,neprime_eq,te_eq,teprime_eq, &
           ti_eq,tiprime_eq
      real, dimension(:), allocatable :: q_calc ! calculated q(ixinterp)       
      real alpha_edge  ! edge value of alpha, used in matgen
      real, dimension(:), allocatable :: alpha ! calculated alpha(ixinterp)
      real, dimension(:), allocatable :: shear ! calculated shear(ixinterp)
      real, dimension(:), allocatable :: pprime_adj ! adjusted value
         ! of 4 Pi p', including ppmult or alphamult
      real, dimension(:), allocatable :: ffprime_adj ! adjusted ff'
      real f_in,q_in ! edge f and q values, read from .eq file
      real qmin ! qmin value for calculating nmres, read from .eq
      integer maxdxinterp  ! maximum delta_x between interpolation
         ! surfaces.  Used to pad size of nustuff arrays properly.
      logical qmonotonic ! flag for whether q profile is monotonic on
             ! the nxinterp surfaces
      real, dimension(:), allocatable :: xinterp_equil,psigrid,rgrid
        ! x and psi and R grid read from .xtopsi file
      real, dimension(:), allocatable :: psixx  ! psi values for xx grid
      real, dimension(:), allocatable :: omegas  ! omega_*p
      real, dimension(:), allocatable :: omegasn  ! omega_*_i_n
      real, dimension(:), allocatable :: omegaspi  ! omega_*_pi
      real omegas_max, omegasn_max, omegaspi_max  ! max values of omega*
      real jedgen1store   ! store jedge, normalization 1, at outer surf
      real, dimension(:), allocatable :: jedgen1_arr ! store jn1 array
      real, dimension(:), allocatable :: jedgen2_arr ! store jn2 array
      real lhspeelstore   ! store lhs of peeling crit for outer surf
      real rhspeelstore   !  "    rhs " "
      real dmstore        ! store Mercier D_m for outer surface
      real tmatasym       ! value of T matrix asymmetry (2nd norm) after
                           !  final iteration
      integer nsurfdat    ! file unit number for *.surfdat file used
        ! by surfcalc and mercy
! 8/02 pbs, new vars for compression and rotation
!  in this version these are not used, just here for namelist compatibility
      real g_comp  ! gamma value for compression (1=isothermal,5/3=adiabat)
      logical compression  ! flag for using compression (default=.false.)
      complex, dimension(:,:), allocatable :: comp_ai,comp_api,comp_appi
!   matrices to hold compressional & rotation terms (v_m and vprime_m) 
!      which enter smat,pmat,qmat
      complex ii  ! to hold sqrt(-1)=cmplx(0.,1.)

      real rot_const  ! multiplier for rotation in simple rotation models
      integer rot_model  !   -1 : tanh rotation file (fixed width for now)
                     !   1-20 : polynomial ~ (1-psinorm**rot_model)
! 6/2020 variables for reading in info from .eq file - needed to ensure
!     nmvac is set consistently here and in the vac code
      real, dimension(:), allocatable :: tsurf,rsurf


contains   

  subroutine init
!-----------------------------------------------------------------------
! sets default values for parameters
! reads in input parameters from inbalpeel
!  Same version works for updown symmetric or antisymmetric
!-----------------------------------------------------------------------
      use command_line
      implicit none
      integer ninput,iosa,m,nmvac_inp
!      namelist/qref_modes/nn,nm,nmvac,nmwinhalf,xwid,nxinterp,dens
      namelist/qref_modes/nn,psimin,nmlow,nmvac,nmwinhalf,nxinterp, &
           dens,rotation,verbose
      namelist/plas/tau,te,etai,ne,vloop,zeff, &
          dx,nd,dw,ndist,lamdist, &
          dmercier,ns,meshtype,autorun,updownsym,n1term,gamsq,igam,funcal, &
          vacuum,nowindow,bug,ppmult,alphamult,shearmult, &
          sa_fix,shear_fix,alpha_fix,analytic_peel,surfplot,splineeq, &
          wsmodel,g_comp,compression,ion_mass,gamim,rot_const,rot_model, &
          kinkterms,newtstep,ncrit_ws,usegamr2
      character*24 ctime
!      character*8 date
!      character*9 ctime
!      character*6 zone
!      integer values(8)
      character*1 ystring
      integer time
      integer i,nunit
!      real f_in,q_in
      integer icount,ierr,ios,nsurfpts
      real rmin, rmaj  ! minor and major radius for use in nmvac 
      real fn1, fn2, fn3  ! functions for setting nmvac
!     code version
      codver='v5.9sym Jun20'
      pi=acos(-1.d0)
!-----------------------------------------------------------------------
! default qref_modes values
!-----------------------------------------------------------------------
      nn=10                     !toroidal mode number
      nmvac=-20  ! will default to formula after nn, q, eps are read
      nmlow=1
      verbose=4  ! default will produce .out file and moderate level of text output
      nmwinhalf=-1 ! will default it to nmvac+3 once nmvac set
      nxinterp=100
      psimin=0.7
      dens=.false.
      rotation=.false.
      rot_const=0.
      rot_model=1
!-----------------------------------------------------------------------
! default plas namelist values
!-----------------------------------------------------------------------
      tau=1.                    !Te/Ti
      te=100.                   !Te (eV)
      etai=1.                   !ratio of n to T lenth scale (=etae)
      ne=0.1                    !electron density (=ni)
      vloop=1.                  !loop voltage
      zeff=1.                   !Zeff
      dx=0.02                   !at the moment, only used in solveit
                                ! to define dxmin=dx/nd
      nd=50                     !at the moment, only used in solveit
                                ! to define dxmin=dx/nd
      dw=0.1                    !not used--originally dx=dx except 
                                ! dw width around rat. surf. 
                                ! where dx=dx/nd
! in qref_modes   xwid=3.       !estimated half-width of fourier modes
                                ! added to innermost x to determine xmax
! in qref_modes   nmwinhalf=3       !total number of modes in mode window
                                ! is 2*nmwinhalf+1
!-----------------------------------------------------------------------
! next two params determine x mesh-see routine getxx
!-----------------------------------------------------------------------
      ndist=100                 ! # of x points betwee rational surf.
      lamdist=0.3               ! packing of points--smaller lamdist is
                                !  more packing
! note: should use either ppmult or (alphamult & shearmult), not both.
!   ppmult will change both alpha and shear, because it leaves ffprime
!   unchanged.  alphamult & shmult adjust both pprime & ffprime to
!   achieve desired multipliers of the equilibrium shear and alpha
      ppmult=1.0                ! multiplier for pprime_eq on all surfaces
      alphamult=1.0             ! multiplier for alpha on all surfaces
      shearmult=1.0             ! multiplier for shear on all surfaces
      sa_fix=.false.            ! whether to fix outer s&alpha
      shear_fix=1.              !  fixed value for outer surface shear
      alpha_fix=1.              ! fixed value for outer surface alpha
      analytic_peel=.false.         ! if true, code quits after outputting
         ! the analytic peeling stability criterion lhspeel-rhspeel
         ! for the outer flux surface
      surfplot=2                ! default to output surface plot data
         ! with nxplt=2*nxinterp+1 radial points
! note: using sa_fix requires that the multipliers (alphamult etc.)
!   be set to one.
!      dmercier=0.1              ! D_M=alpha*dmercier
                                ! >0 is stabilizing ?
      dmercier=0.0   ! should be zero (is this for comp to s-alpha code?)
      wsmodel=0                 ! don't use omegastar to set gamsq
      ncrit_ws=13               ! default rollover from BOUT++ cbm18 results
      ion_mass=2.0              ! default=deuterium
      g_comp=1.                 ! isothermal compressionality
      compression=.false.
      ns=501                    !number of poloidal points on ref. surf.
!      nxinterp=10               !# of x points used for interpolation
                                ! of non-resonant quantities across
                                ! full range of x: xmax>=x>=xmin
      meshtype=1                !1 uses dist variables, 2 is wilson's orig
                                !   xx(i)
      autorun=.false.           !if true then last line of output is 
                                !dW for least stable or most unstable mode
                                !useful for scans using scripts
      usegamr2=.false.          !if true, then iterate to convergence using gamr2
                                !(the eigenvalue without forced hermitian symmetry)
                                ! in addition to gamr.  This option is useful for testing
                                ! the impact of equilibrium and other errors on the
                                ! growth rate
      updownsym=.false.         !set true for updown symmetric equilibria
                                ! to save storage and cpu time
      n1term=2                  !set to zero to switch off additional 1/n terms
      gamsq=0.                  !initial set of square of growth rate
      igam=0                    !set to max no. iterations to find 
                                !self-consistent gamsq...if </=1 code solves 
                                !for given gamsq
      gamim=0.                  ! not used in this version
      vacuum=.true.             !set false to make dW_vac=0
      funcal=.true.             !set to false if eigenfunction not required
      kinkterms=.true.          !set to false to eliminate kink terms
      nowindow=.false.          !if true then full range of m used at 
                                !every x location
      splineeq=.false.          ! true to spline au* & iu* arrays over
                                ! xinterp grid rather than using linear
                                ! interpolation.  use only when nxinterp
                                ! is small (<~20) or it will be _very_ slow
      newtstep=0.2              ! default initial Newton step of 20%
      do i=1,10
         bug(i)=0.              !spare parameters for program development
      end do
!-----------------------------------------------------------------------
! read new namelist values from file inbalpeel
!-----------------------------------------------------------------------
      ninput=20
! determine the name of the input file:
      icount=iargc()
      if(icount  >=  1) then
          call cl_getarg(1,runname,lrunname,ierr)
      else
          write(6,*) 'What is the name of the *.in input file?'
          write(6,*) '(without the .in suffix):'
          read(5,*) runname
          lrunname=index(runname,' ')-1
      endif
     



      open(unit=ninput,file=runname(1:lrunname)//'.in',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem opening file',runname(1:lrunname)//'.in'
         stop
      endif

      read(ninput,qref_modes)
      read(ninput,plas)
      close(ninput)

      
      if (verbose .ge. 1) write(6,*) 'Running elite for ', runname(1:lrunname)
      ! open runname.out file for standard text output
      outfile=39  
      if (verbose .ge. 2) open(unit=outfile,file=runname(1:lrunname)//'.out',iostat=ios)

!  c*24
!      call date_and_time(date,ctime,zone,values)
!      cdatetime=ctime(1:2)//':'//ctime(3:4)//' '//date(7:8)//'/'//date(5:6)// &
!            '/'//date(1:4)
      cdatetime =  ctime (time ())
      if (verbose .ge. 2) then
         write(6,*) cdatetime
         write(6,*) "code version: ",codver
         write(6,qref_modes)
         write(6,plas)

         write(outfile,*) cdatetime
         write(outfile,*) "code version: ",codver
         write(outfile,qref_modes)
         write(outfile,plas)
      endif

! read f_in and q_in only from the name.eq file (before used only by vac)
!   add qmin, 9/01
!  6/20 also read rsurf the name.eq file so that nmvac options can be setup
!    identically as in the vac code
      nunit=71
      open(unit=nunit,file=runname(1:lrunname)//'.eq',status='old',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'does file', runname(1:lrunname)//'.eq',' exist?'
         stop
      endif
      read(nunit,*) nsurfpts
      allocate ( tsurf(nsurfpts), rsurf(nsurfpts) )
      tsurf=0.
      rsurf=0.
      read(nunit,*) f_in,q_in,qmin
      read(nunit,*) (tsurf(i),i=1,nsurfpts)
      read(nunit,*) (rsurf(i),i=1,nsurfpts)
      close(nunit)

      rmaj=(rsurf(nsurfpts/2+1)+rsurf(1))/2.
      rmin=(rsurf(1)-rsurf(nsurfpts/2+1))/2.
      if (verbose .ge. 2) then
         write(*,*) 'read surfdat, q_in=',q_in,'rmaj=',rmaj, &
              'rmin=',rmin,'eps=',rmin/rmaj,'nn*q_in=',real(nn)*q_in
      endif

! set default values of nmvac and nmwinhalf based on nn:
    if (nmvac.lt.0) then
       nmvac_inp=nmvac
       if (nmvac.eq.-1) then
          nmvac=min(nn+3,12)
       else if (nmvac.eq.-2) then
          nmvac=min(nn+2,11)
       else if (nmvac.eq.-3) then 
          nmvac=min(nn+4,15)
       else if (nmvac.eq.-4) then 
          nmvac=min(int(sqrt(nn+1.)),12)
       else if (nmvac.eq.-5) then
          if (nn .lt. 10.) then
             nmvac=nn+2
          else
             nmvac=min(int((nn*2.)**0.6667+6.),50)
          endif
       else if (nmvac.eq.-6) then
          if (nn .lt. 7) then
             nmvac=nn+1
          else
             nmvac=min(nn+2,15)
          endif
       else if (nmvac.eq.-15) then 
          nmvac=max(int(nn*2.5),10)
       else if (nmvac.eq.-16) then 
          nmvac=nn*1.5
       else if (nmvac.eq.-17) then 
          nmvac=nn*3.5
       else if (nmvac.eq.-18) then 
          nmvac=min(int(nn*1.5),30)
       else if (nmvac.eq.-19) then 
          nmvac=min(nn*2,30)
       else if (nmvac .eq. -20) then
! 6/20 new default based on assessment of nmvac convergence in 11 cases with a wide range of q_in and aspect ratio
          fn1=2.+0.4333*real(nn)*q_in
          fn2=15.+0.0343*(rmin/rmaj)*q_in*(real(nn)*q_in-30.)
          fn3=15.+6.*(rmin/rmaj)*q_in
          nmvac=min(int(fn1),int(fn3))
          if (nn*q_in .ge. 30) nmvac=min(nmvac,int(fn2))
          if (verbose .ge. 1) then
             write (*,*) 'nmvac=-20 set,', &
                  'fn1=',int(fn1),' fn2=',int(fn2),' fn3=',int(fn3)
             write (outfile,*) 'nmvac=-20 set,', &
                  'fn1=',int(fn1),' fn2=',int(fn2),' fn3=',int(fn3)
          endif
       else if (nmvac .eq. -21) then
! 6/20 more relaxed, faster estimate based on assessment of nmvac convergence in 11 cases with a wide range of q_in and aspect ratio
          fn1=1.+0.4*real(nn)*q_in
          fn2=13.+0.03*(rmin/rmaj)*q_in*(real(nn)*q_in-30.)
          fn3=13.+4.*(rmin/rmaj)*q_in
          nmvac=min(int(fn1),int(fn3))
          if (nn*q_in .ge. 30) nmvac=min(nmvac,int(fn2))
          if (verbose .ge. 1) then
             write (*,*) 'nmvac=-21 set,', &
                  'fn1=',int(fn1),' fn2=',int(fn2),' fn3=',int(fn3)
             write (outfile,*) 'nmvac=-21 set,', &
                  'fn1=',int(fn1),' fn2=',int(fn2),' fn3=',int(fn3)
          endif
       else if (nmvac .eq. -22) then
! 6/20 more agressive, slow estimate based on assessment of nmvac convergence in 11 cases with a wide range of q_in and aspect ratio
          fn1=3.+0.43*real(nn)*q_in
          fn2=20.+0.06*(rmin/rmaj)*q_in*(real(nn)*q_in-30.)
          fn3=20.+10.*(rmin/rmaj)*q_in
          nmvac=min(int(fn1),int(fn3))
          if (nn*q_in .ge. 40) nmvac=min(nmvac,int(fn2))
          if (verbose .ge. 1) then 
             write (*,*) 'nmvac=-22 set,', &
                  'fn1=',int(fn1),' fn2=',int(fn2),' fn3=',int(fn3)
             write (outfile,*) 'nmvac=-22 set,', &
                  'fn1=',int(fn1),' fn2=',int(fn2),' fn3=',int(fn3)
           endif
       else
          write(*,*) ' Entered unsupported nmvac value=',nmvac
          if (verbose .ge. 2) write(outfile,*) &
               ' Entered unsupported nmvac value=',nmvac
          stop
       endif
       if (verbose .ge. 2) then
          write(*,*) ' nmvac=',nmvac_inp,' input, set to nmvac=',nmvac
          write(outfile,*) ' nmvac=',nmvac_inp,' input, set to nmvac=',nmvac
       endif
    endif

    if (nmwinhalf.eq.-1) then
       nmwinhalf=nmvac+3
       if (verbose .ge. 2) then
          write(*,*) 'setting nmwinhalf to default value of ',nmwinhalf
          write(outfile,*) 'setting nmwinhalf to default value of', &
            nmwinhalf
       endif
    endif

! check that sa_fix is true only when multipliers are one
     if ((sa_fix).and.((ppmult.ne.1.).or.(alphamult.ne.1.).or. &
          (shearmult.ne.1.)) ) then
         write(*,*) 'fixing outer s&alpha with sa_fix requires'
         write(*,*) '  ppmult=alphamult=shearmult=1.'
         stop
      endif

! check to make sure that ppmult and (alphamult and shearmult) aren't
!  mistakenly being used together
      if ((ppmult.ne.1.0).and.((alphamult.ne.1.).or.(shearmult.ne.1))) then
         write(*,*) 'should adjust params either by using ppmult or'
         write(*,*) ' alphamult and shearmult, not both'
         write(*,*) 'ppmult=',ppmult,' alphamult=',alphamult, &
              'shearmult=',shearmult
         stop
      endif

      if (.not. kinkterms) then
         write(*,*) 'WARNING!! Running with kink terms turned off.  This is', &
              ' unphysical and should be used for physics studies only!!'
         if (verbose .ge. 2) write(outfile,*) &
              'WARNING!! Running with kink terms turned off.  This is', &
              ' unphysical and should be used for physics studies only!!'
      endif



! must read or calculate qsurf from .rzbp before doing this        
!!  edge q....
!      qa=q0-del/nn
!      qref=qa
!      m0=nn*q0                  !find first rational surface outside plasma

      qa=q_in
      qref=q_in
      m0=int(nn*q_in+1.)
      q0=real(m0)/real(nn)
      del=(q0-q_in)*nn
      if (verbose .ge. 2) then
         write(*,*) 'q_in= ',q_in,' m0=',m0,' q0=',q0,' del=',del,'qmin=',qmin
         write(outfile,*) 'q_in= ',q_in,' m0=',m0,' q0=',q0,' del=',del, &
              'qmin=',qmin
      endif
      mmax=m0+nmvac-1
! 9/5/01 qmin now read in in .eq file

      nmres=(qa-qmin)*nn+1   ! 9/01 could make this more precise,
         ! but as long as it's consistent with vac, should be fine
      nm=nmres+nmvac+nmlow ! total m's in calculation
      if (verbose .ge. 2) write(*,*) &
           'number of resonant ms in plasma, nmres=',nmres, &
           'total nm=nmres+nmvac+nmlow=',nm
      
      mmin=mmax-nm+1
      if (verbose .ge. 2) write(*,*) 'mmin=',mmin,'mmax=',mmax

      if (nmvac < 1) then
         write(*,*) 'must have at least one vacuum mode (nmvac>=1)'
         stop
      endif

      if(mmin < 1) then
!         write(*,*)'mmin=m0+nmvac-nm must be >=1'
!         stop
         write(*,*) 'WARNING, mmin<1, mmin=',mmin
         if (verbose .ge. 2) write(outfile,*) 'WARNING, mmin<1, mmin=',mmin
      endif

!      xmax=m0-mmin+xwid
!      xmax=m0-mmin
      xmax=nn*(q0-qmin)
      xmin=del
      if(updownsym) then
         nmodes=nm
         mindex=1
      else
         nmodes=2*nm            ! useful because re and im has doubled 
         mindex=2               ! length of vectors
      endif
      if (verbose .ge. 2) then
         write(*,*)'xmax=',xmax,' xmin=',xmin,' nmodes=',nmodes
         write(outfile,*)'xmax=',xmax,' xmin=',xmin,' nmodes=',nmodes
      endif
! more compatibility tests
      if((.not.nowindow) .and. (nmwinhalf.lt.(nmvac+1))) then
         write(6,*) 'WARNING--resetting nmwinhalf=nmvac+1'
         if (verbose .ge. 2) write(outfile,*) 'WARNING--resetting nmwinhalf=nmvac+1'
         nmwinhalf=nmvac+1
      endif
! finally need number of modes that can be non-zero at edge
      mwindow=min0(mindex*(2*nmwinhalf+1),nmodes)
      if (nowindow) then
         mwindow=nmodes
      endif
! 3/02 turn windowing off if it's not going to be used 
      if ((.not.nowindow).and.(mwindow.eq.nmodes)) then
         if (verbose .ge. 2) then
            write(*,*) 'window size larger than m domain'
            write(*,*) 'windowing turned off'
            write(outfile,*) 'window size larger than m domain'
            write(outfile,*) 'windowing turned off'
         endif
         nowindow=.true.
      endif

      if (verbose .ge. 2) then
         write(6,*)' mwindow=',mwindow
         write(outfile,*)' mwindow=',mwindow
      endif

! 6/01 with new windowning algorithm nedge should always be mwindow
      nedge=mwindow

! check for rotation or compression, not supported by this version
      if (rotation) then
         write(*,*) '***rotation is not supported in this version'
         write(*,*) 'set rotation to false or use version not marked sym'
         stop
      endif
      if (compression) then
         write(*,*) '***compression is not supported in this version'
         write(*,*) 'set compression to false or use version not marked sym'
         stop
      endif
      if (.not.updownsym) then
         write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write(*,*) '***this version is intended for efficient treatment of'
         write(*,*) ' updown symmetric equilibria.  For updown asymmetric'
         write(*,*) ' cases please use version not labelled sym.'
         write(*,*) 'Using the script elite.x will choose correct version'
         write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write(*,*) ' If you wish to continue for testing purposes, enter y'
         read(*,*) ystring
         if (ystring .ne. 'y') stop
!         stop
      endif

!  read the flux surface arrays from name.eqdat
      call rdeqdat
!  read .xtopsi file to get psigrid and xinterp_equil
      call rdxtopsi
      
!  calculate maxdxinterp=max change in x coord between interpolation
!    surfaces.  needed to properly size the nustuff arrays for interpolation
!    in matsetup  7/01
!    9/01 use more accurate xinterp_equil grid for this now
      maxdxinterp=1
      do i=2,nxinterp
!         maxdxinterp=max0(maxdxinterp,int(nn*(q_eq(i-1)-q_eq(i))+1.05))
         maxdxinterp=max0(maxdxinterp, &
              int(abs(xinterp_equil(i-1)-xinterp_equil(i))+1.01) )
         ! use 1.05 to allow for small differences between q_eq and
         !  qsurf calculated by nustuff for the x grid
!         write(*,*) 'nn dq=',nn*(q_eq(i-1)-q_eq(i)), &
!              int(nn*(q_eq(i-1)-q_eq(i))+1.05)
      enddo
!      write(*,*) 'maxdxinterp=',maxdxinterp
      if (verbose .ge. 3) write(outfile,*) 'maxdxinterp=',maxdxinterp

      call alloc_nustuff

      return
 
  end subroutine init


  subroutine alloc_arrays

!     integer mwindow  ! now a global
    implicit none

     kx=(xmax+1)*ndist ! max number of xx points (perhaps could be precise?)
     if (meshtype == 2 .or. meshtype == 4) kx=ndist  
     ! new equal space option uses nx=ndist 6/01
     if (meshtype == 3) kx=ndist*(nxinterp-1)  ! grid using ndist points
               ! between each nxinterp surface

     allocate( rl(ns),zl(ns),bpl(ns),uprimel(ns),cosul(ns),sinul(ns) )
     rl=0.; zl=0.; bpl=0.; uprimel=0.; cosul=0.; sinul=0.
     allocate( dpsi1l(ns),d2psi1l(ns) )
     dpsi1l=0.; d2psi1l=0.

     allocate( xinterp(nxinterp),dw_vac(nmodes,nmodes),xx(kx),psixx(kx))
     xinterp=0.; dw_vac=0.; xx=0.; psixx=0.

     allocate( psiv_eq(nxinterp),pprime_eq(nxinterp),f_eq(nxinterp), &
          ffprime_eq(nxinterp),vprime_eq(nxinterp),q_eq(nxinterp), &
          ppp_eq(nxinterp),ffpp_eq(nxinterp),qp_eq(nxinterp),qpp_eq(nxinterp) )
     psiv_eq=0.; pprime_eq=0.; f_eq=0.; ffprime_eq=0.
     vprime_eq=0.; q_eq=0.; ppp_eq=0.; ffpp_eq=0.; qp_eq=0.; qpp_eq=0.

     allocate(ne_eq(nxinterp),neprime_eq(nxinterp),te_eq(nxinterp), &
          teprime_eq(nxinterp),ti_eq(nxinterp),tiprime_eq(nxinterp) )
     ne_eq=0.; neprime_eq=0.; te_eq=0.; teprime_eq=0.; ti_eq=0.; tiprime_eq=0.

     allocate( q_calc(nxinterp),alpha(nxinterp),shear(nxinterp), &
          pprime_adj(nxinterp),ffprime_adj(nxinterp), &
          jedgen1_arr(nxinterp),jedgen2_arr(nxinterp) )
     q_calc=0.; alpha=0.; shear=0.; pprime_adj=0.; ffprime_adj=0.
     jedgen1_arr=0.; jedgen2_arr=0.

     allocate( omegas(nxinterp),omegasn(nxinterp),omegaspi(nxinterp) )
     omegas=0.; omegasn=0.; omegaspi=0.

     allocate( xinterp_equil(nxinterp), psigrid(nxinterp), rgrid(nxinterp) )
     xinterp_equil=0.; psigrid=0.; rgrid=0.

  end subroutine alloc_arrays


  subroutine alloc_nustuff
  !  allocate nustuff arrays.  separate routine is needed because must
  !   use q_eq to calculate maxdxinterp to get needed size for these 7/01.

     implicit none
     integer lwindow ! mwindow+2*mindex*maxdxinterp, larger window with 
                     !   room on both sides
                     !   used to allow proper offset of nustuff arrays in 
                     !   matgen  7/01
     lwindow=min0(mwindow+2*mindex*maxdxinterp,nmodes)

! 1/02 replace au.. arrays with one large array, 
!   auarr(19,lwindow,lwindow,nxinterp) where elements 1-19 correspond
!   to the old array set in the order below
     allocate( auarr(19,lwindow,lwindow,nxinterp) )
     auarr=0.

!     allocate( aumk(lwindow,lwindow,nxinterp), &
!          au2mk(lwindow,lwindow,nxinterp), &
!          aum2k(lwindow,lwindow,nxinterp),au2m(lwindow,lwindow,nxinterp), &
!          au2k(lwindow,lwindow,nxinterp),aum(lwindow,lwindow,nxinterp), &
!          auk(lwindow,lwindow,nxinterp),aund(lwindow,lwindow,nxinterp), &
!          aupmk(lwindow,lwindow,nxinterp),aup2mk(lwindow,lwindow,nxinterp), &
!          aupm2k(lwindow,lwindow,nxinterp),aup2m(lwindow,lwindow,nxinterp), &
!          aup2k(lwindow,lwindow,nxinterp),aupm(lwindow,lwindow,nxinterp), &
!          aupk(lwindow,lwindow,nxinterp),aupnd(lwindow,lwindow,nxinterp), &
!          auppmk(lwindow,lwindow,nxinterp),aupp2mk(lwindow,lwindow,nxinterp),&
!          auppm2k(lwindow,lwindow,nxinterp) )
!     aumk=0.; au2mk=0.; aum2k=0.; au2m=0.; au2k=0.; aum=0.; auk=0.
!     aund=0.; aupmk=0.; aup2mk=0.; aupm2k=0.; aup2m=0.; aup2k=0.
!     aupm=0.; aupk=0.; aupnd=0.; auppmk=0.; aupp2mk=0.; auppm2k=0.

! 1/02 replace iu.. arrays with one large array, 
!   iuarr(9,lwindow,lwindow,nxinterp) where elements 1-9 correspond
!   to the old array set in the order below
     allocate( iuarr(9,lwindow,lwindow,nxinterp) )
     iuarr=0.

! 1/02 allocate arrays for spline coefficients only in splineeq is true
     if (splineeq) then
        allocate( auspline(19,nmodes,nmodes,nxinterp), &
             iuspline(7:9,nmodes,nmodes,nxinterp) )
        auspline=0.; iuspline=0.
     endif

!     allocate( iuk(lwindow,lwindow,nxinterp),ium(lwindow,lwindow,nxinterp), &
!          iund(lwindow,lwindow,nxinterp),iupk(lwindow,lwindow,nxinterp), &
!          iupm(lwindow,lwindow,nxinterp),iupnd(lwindow,lwindow,nxinterp), &
!          iuppk(lwindow,lwindow,nxinterp),iuppm(lwindow,lwindow,nxinterp), &
!          iuppnd(lwindow,lwindow,nxinterp) )
!     iuk=0.; ium=0.; iund=0.; iupk=0.; iupm=0.; iupnd=0.
!     iuppk=0.; iuppm=0.; iuppnd=0.

     allocate( asurfu(lwindow,lwindow), &
             asurfup(lwindow,lwindow),isurfu(lwindow,lwindow), &
             isurfup(lwindow,lwindow) )
     asurfu=0.; asurfup=0.; isurfu=0.; isurfup=0.

  end subroutine alloc_nustuff


  subroutine rdxtopsi

! 9/01 read psigrid from the .xtopsi file, also reads x
!   grid, but don't use this - recalculated in nustuff

! 6/02 also read r grid (R on outer midplane) for use in
!   estimating ped width for omega_star models

    implicit none
    integer nxtopsi,nxinterp_equil,ios

    nxtopsi=96
    if (verbose .ge. 3) write(outfile,*) 'reading ',runname(1:lrunname)//'.xtopsi'
    open(unit=nxtopsi,file=runname(1:lrunname)//'.xtopsi', &
         iostat=ios)
    if(ios.ne.0) then
       write(6,*) 'could not find file ',runname(1:lrunname)//'.xtopsi'
       stop
    endif

    read(nxtopsi,*) nxinterp_equil
    if (nxinterp_equil.ne.nxinterp) then
       write(*,*) 'nxinterp value read from .xtopsi file', &
            ' does not equal nxinterp',nxinterp_equil,nxinterp
       stop
    endif

    read(nxtopsi,*) xinterp_equil
!    write(*,*) 'xinterp grid read from .xtopsi',xinterp_equil
    if (abs(xinterp_equil(1)-del) > 1.d-3) then
       write(*,*) 'x(1) read from .xtopsi file does not match del', &
            xinterp_equil(1),del
       stop
    endif

    read(nxtopsi,*) psigrid
!    write(*,*) 'psigrid read from .xtopsi',psigrid
    if (verbose .ge. 3) write(outfile,*) 'psigrid read from .xtopsi',psigrid
    
    read(nxtopsi,*) rgrid
 !   write(*,*) 'rgrid read from .xtopsi',rgrid
    if (verbose .ge. 3) write(outfile,*) 'rgrid read from .xtopsi',rgrid

    close(nxtopsi)


  end subroutine rdxtopsi

  subroutine rdeqdat

! reads flux surface data for the nxinterp reference flux surfaces
!  from the file runname.eqdat created by the equilibrium code (eliteeq)
!   PBS 2/8/00

      implicit none
      integer i,neqdat,ios
      character*30 dumstring
      integer npts_loc,nxinterp_loc

      neqdat=97

      if (verbose .ge. 3) write(outfile,*) 'reading ',runname(1:lrunname)//'.eqdat'
      open(unit=neqdat,file=runname(1:lrunname)//'.eqdat', &
           iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'could not find file ',runname(1:lrunname)//'.eqdat'
         stop
      endif
      
      read(neqdat,*) dumstring
      read(neqdat,*) nxinterp_loc,npts_loc
      if (verbose .ge. 4) write(outfile,*) dumstring,nxinterp_loc,npts_loc
      if (nxinterp_loc.ne.nxinterp) then
!         write(*,*) 'nxinterp in ',runname(1:lrunname)//'.eqdat', &
!              ' does not match nxinterp',nxinterp_loc,nxinterp
         if (verbose .ge. 2) write(outfile,*) &
              'nxinterp in ',runname(1:lrunname)//'.eqdat', &
              ' does not match nxinterp',nxinterp_loc,nxinterp
         ! use nxinterp from .eqdat and give warning
            if (verbose .ge. 1) write(*,*) ' Using nxinterp value from ', &
                 runname(1:lrunname)//'.eqdat,', 'nxinterp= ',nxinterp_loc
            if (verbose .ge. 2) write(outfile,*) ' Using nxinterp value from ', &
                 runname(1:lrunname)//'.eqdat,', 'nxinterp= ',nxinterp_loc
            nxinterp=nxinterp_loc
!            stop
      endif
! now nxinterp has its final value, allocate arrays
      if (verbose .ge. 5) write(outfile,*)' Allocating arrays'
      call alloc_arrays
      if (verbose .ge. 5) write(outfile,*)' Done allocating'

      read(neqdat,*) dumstring
      read(neqdat,*) (psiv_eq(i),i=1,nxinterp)
      if (verbose .ge. 6) write(*,*) dumstring,psiv_eq
      read(neqdat,*) dumstring
      read(neqdat,*) (pprime_eq(i),i=1,nxinterp)
!  2/29/00 Miller's pprime is actually 4*pi*p'!!!!!!!!
      pprime_eq=4.*pi*pprime_eq
      read(neqdat,*) dumstring
      read(neqdat,*) (f_eq(i),i=1,nxinterp)
      read(neqdat,*) dumstring
      read(neqdat,*) (ffprime_eq(i),i=1,nxinterp)
      read(neqdat,*) dumstring
      read(neqdat,*) (vprime_eq(i),i=1,nxinterp)
      read(neqdat,*) dumstring
      read(neqdat,*) (q_eq(i),i=1,nxinterp)
      read(neqdat,*) dumstring
      read(neqdat,*) (ppp_eq(i),i=1,nxinterp)
! Also multiply p'' by 4*pi !!!!
      ppp_eq=4.*pi*ppp_eq
      read(neqdat,*) dumstring
      read(neqdat,*) (ffpp_eq(i),i=1,nxinterp)
      read(neqdat,*) dumstring
      read(neqdat,*) (qp_eq(i),i=1,nxinterp)
      read(neqdat,*) dumstring
      read(neqdat,*) (qpp_eq(i),i=1,nxinterp)

      if (verbose .ge. 4) then
         write(outfile,*) 'ppp_eq=',ppp_eq
         write(outfile,*) 'ffpp_eq=',ffpp_eq
         write(outfile,*) 'qp_eq=',qp_eq
         write(outfile,*) 'qpp_eq=',qpp_eq
      endif
      if (dens) then
        read(neqdat,*) dumstring
        read(neqdat,*) (ne_eq(i),i=1,nxinterp)
        read(neqdat,*) dumstring
        read(neqdat,*) (neprime_eq(i),i=1,nxinterp)
        read(neqdat,*) dumstring
        read(neqdat,*) (te_eq(i),i=1,nxinterp)
        read(neqdat,*) dumstring
        read(neqdat,*) (teprime_eq(i),i=1,nxinterp)
        read(neqdat,*) dumstring
        read(neqdat,*) (ti_eq(i),i=1,nxinterp)
        read(neqdat,*) dumstring
        read(neqdat,*) (tiprime_eq(i),i=1,nxinterp)
      end if
      close(neqdat)

      if (verbose .ge. 2) write(outfile,*) 'q_eq=',q_eq

  end subroutine rdeqdat


end module elite_data








