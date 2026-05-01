
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ELITE (Edge Localized Instabilities in Tokamak Experiments)
!
! Authors:  P.B. Snyder, H.R. Wilson, R.L. Miller
!
! Copyright(c) 1999-2007, All rights reserved
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


module eliteeq_data
!  Philip B. Snyder 1/26/00

  implicit none

!   variables needed in equilibrium code

      real q_surf
      real, dimension(:), allocatable :: theta,rpts,zpts,drdt,dzdt,bppts,omega
      real, dimension(:), allocatable :: om_edge

! Need ns variable from namelist plas to spline functions into equal arc length
      integer ns
! variables for namelist equil
!      real k,gamma,rmajor,rminor,f_surf,kappa,triang,kappap,triangp,deltap
      real rminor,rmajor,f_surf
      integer npts
      character*4 shape
      character runname*50  ! Name of this run.  Associated i/o files
                        ! will be named RUNNAME.in, etc.
      integer lrunname      ! length of runname

! new input variables for reading eqdsk files
      integer npsi              !size of mapped psi mesh
      real percenflux           !outermost flux surface is 
                                !psiv(npsi)=psiaxis+(psilim-psiaxis)*percenflux
      integer rotate            !=1 for rotation to read eqdsk correctly
      real alpsi
      integer npfit             !threshhold pts for using furpl
      logical delmin            ! whether to search for nn which minimizes del
      logical dens
      logical rotation
      integer nnmin,nnmax       ! range to search for nn which minimizes del
      logical setdel            ! whether to calculte q_new/q_old needed
         ! to get del=del_fix.  designed for use with scripts which adjust
         ! q in the efit
      real del_fix              ! desired value of del
      logical qafix    ! if true, and setdel is true, elite will modify
         !  bps and bppts so that calculated q will lead to the specified
         !  del value  (meaningless if setdel is false)
      logical wrtgato  ! if true write out dskgato file filename_out.dskgato
         !  for use in GATO (used to convert eqin or dskbal files to dskgato)

! variables for eqdsk info
      real, dimension(:), allocatable :: xgrid,zgrid
      real xaxis                !x of mag axis
      real zaxis                !z of mag axis
      real psiaxis              !psi at mag axis
      real psilim               !psi at limiter
      real dx                   !x(2)-x(1)
      real dz                   !z(2)-z(1)
      real, dimension(:,:), allocatable :: psixz,bpsq
!      real psixz(nxd,nzd)         !complete psi array
!      real bpsq(nxd,nzd)        !calculated in readeqdsk
      real, dimension(:), allocatable :: sp,spp,sf,sffp,qpsi
!      real sp(nxd)              !pressure on equally spaced psi
!      real spp(nxd)             !pprime on equally spaced psi
!      real sf(nxd)              !f on equally spaced psi
!      real sffp(nxd)            !ffprime on equally spaced psi
!      real qpsi(nxd)            !qpsi on equally spaced psi
      real, dimension(:), allocatable :: xbndry,zbndry
!      real xbndry(kbnd)         !x position of boundary points
!      real zbndry(kbnd)         !z position of boundary points
      real, dimension(:), allocatable :: xlim,zlim
!      real xlim(klim)           !x position of limiter points
!      real zlim(klim)           !z position of limiter points
      real, dimension(:), allocatable :: pressw,pwprim,rho0,rho0p
!      real pressw(nxd)          !see Lang's defn for press with rotation
!      real pwprim(nxd)          !see Lang's defn for press with rotation
!      real rho0(nxd)            !see Lang's defn for press with rotation
!      real rho0p(nxd)           !d(rho0)/d(psi)
      real rvtor                !reference r for rotation
      integer nbndry            !number of boundary points
      integer nlim              !number of limiter points
      integer kvtor             !>0 means rotation
      integer nmass             !>0 means rho0 specified
      integer nx                !size of x mesh
      integer nz                !size of z mesh

! ===== variables read in from dskbal or obtained from mapper=====
      real, dimension(:), allocatable :: psic,psiv,pprime,fval, &
           ffprime,chipsi,qsfin,ppp,ffpp,qsfinp,qsfinpp,pressval,qcalc
      real, dimension(:), allocatable :: nel,nprime,tel,teprime,tion,tiprime
! ppp is p'', ffpp=(ff')'
!      real psic(kpsi)           !the psi coordinates
!      real psiv(kpsi)           !the real poloidal flux array (a.k.a. chi)
!      real pprime(kpsi)         !pprime
!      real fval(kpsi)           !this is f as in Btor=f/R
!      real ffprime(kpsi)        !this is ff'
!      real chipsi(kpsi)         !d(psiv)/d(psic)
!      real qsfin(kpsi)          !the safety factor from input
!      real nel(kpsi)            !electron density
!      real nprime(kpsi)         !and its psi-derivative
!      real tel(kpsi)            !electron temperature
!      real teprime(kpsi)        !and its psi derivative
!      real tion(kpsi)           !ion temperature
!      real tiprime(kpsi)        !and its psi derivative
      real, dimension(:,:), allocatable :: xs,zs,bps
!      real xs(kpsi,kthet)       !R-coordinate of pts. on psiv contour
!      real zs(kpsi,kthet)       !Z-coordinate
!      real bps(kpsi,kthet)      !bp in flux coordinates
      real, dimension(:), allocatable :: press,pw,pwp,rho,rhop
      real, dimension(:,:,:), allocatable :: csplpsi 

! vars needed for in equalarc
      real, dimension(:,:), allocatable :: arcsur 

      real pi   
      integer nn,mmax,nmwinhalf,nm,nxinterp,nmvac
      integer nmlow  ! number of below resonant m's added
      real psimin  ! minimum psi value for cutoff of equilibrium
      real q0,del
      integer verbose  ! sets how much output 0=minimal 10=everything,default=4
!,xwid    

! new toq file reading vars 4/3/00
      real dthe,dpsi
! allow getting rid of outer toq surfaces
      integer garbage
      integer outfile ! file number for .eqout file to replace text output
! 9/01 need to write qmin for vac and plas to use to get nm
      real qmin


contains

    subroutine init
!      use command_line
      integer nunit,ios,nmvac_inp
      namelist/qref_modes/nn,psimin,nmlow,nmvac,nmwinhalf,nxinterp, &
           dens,rotation,verbose
      integer icount,ierr,iargc
      character*5 test_end
    real tau,te,etai,ne,vloop,zeff,dx,dw,lamdist,dmercier,gamsq, &
         bug(10),ppmult,alphamult,shearmult,shear_fix,alpha_fix, &
         g_comp,ion_mass,gamim,rot_const,newtstep,rot_pedmid, &
         rot_pedwid,rot_sep,rot_axis,rot_expin,rot_expout,shootpsi,psimax
    integer  nd,ndist,meshtype,n1term,igam,surfplot, &
         wsmodel,rot_model,ncrit_ws,ific
    logical autorun,updownsym,funcal,vacuum,nowindow,analytic_peel,splineeq, &
         compression,sa_fix,kinkterms,usegamr2,idealwall
      external iargc
!      namelist/equil/k,gamma,rminor,rmajor,npts,f_surf,shape,&
!          kappa,triang,kappap,triangp,deltap,percenflux,alpsi,npsi, &
!          delmin,nnmin,nnmax,setdel,del_fix,qafix
      namelist/equil/npts,shape,&
          percenflux,alpsi,npsi, &
          delmin,nnmin,nnmax,setdel,del_fix,qafix,wrtgato
    namelist/plas/tau,te,etai,ne,vloop,zeff, &
         dx,nd,dw,ndist,lamdist, &
         dmercier,ns,meshtype,autorun,updownsym,n1term,gamsq,igam,funcal, &
         vacuum,nowindow,bug,ppmult,alphamult,shearmult, &
         sa_fix,shear_fix,alpha_fix,analytic_peel,surfplot,splineeq, &
         wsmodel,g_comp,compression,ion_mass,gamim,rot_const,rot_model, &
         kinkterms,newtstep,rot_pedmid,rot_pedwid,rot_sep,rot_axis, &
         rot_expin,rot_expout,ncrit_ws,usegamr2,idealwall,ific, &
         shootpsi,psimax
!    namelist/plas/tau,te,etai,ne,vloop,zeff, &
!         dx,nd,dw,ndist,lamdist, &
!         dmercier,ns,meshtype,autorun,updownsym,n1term,gamsq,igam,funcal, &
!         ific,vacuum,nowindow,bug,ppmult,alphamult,shearmult, &
!         sa_fix,shear_fix,alpha_fix,analytic_peel,surfplot,splineeq, &
!         wsmodel,g_comp,compression,ion_mass,gamim,rot_const,rot_model, &
!         kinkterms,newtstep,ncrit_ws,usegamr2,idealwall

      pi=acos(-1.)
      garbage=0
!-----------------------------------------------------------------------
! default values for namelist equil
!-----------------------------------------------------------------------
!      k=0.96
!      gamma= 3.1416
!      rminor=1.
!      rmajor=1000.
!      kappa=1.8                 !elongation--shape='dee'
!      kappap=0.                 !d(elongation)/dr--shape='dee'
!      triang=0.5                !triangularity--shape='dee'
!      triangp=0.0               !d(triang)/dr--shape='dee'
!      deltap=0.                 !d(rmajor)/dr--shape='dee'
      npts=501 ! # of points in equil boundary
!      f_surf=1.
      shape='toq'
!-----------------------------------------------------------------------
! default values for namelist qref_modes
!-----------------------------------------------------------------------
      nn=10
!      nm= 12
!      nmvac=5
      nmvac=-1  ! will default it to min(nn+3,12) after nn is read
      nmlow=1
!      nmwinhalf=8
      nmwinhalf=-1 ! will default it to nmvac+3 once nmvac set
!      xwid=3.
      nxinterp=100
      psimin=0.7
      dens=.false.     ! default constant density profile
      rotation=.false.
      verbose=4   ! some text output, not too much

! defaults for eqdsk input vars
      rotate=0
      percenflux=0.99
      npsi=170
      alpsi=-0.97  ! extreme concentration of surfaces toward edge
      npfit=25
      delmin=.false.   ! default not to search for nn which minimizes del
      nnmin=nn
      nnmax=nn
      setdel=.false.
      del_fix=0.1
      qafix=.false.
      wrtgato=.false.

!-----------------------------------------------------------------------
! read namelists
!-----------------------------------------------------------------------


! determine the name of the input file:
      icount=iargc()
      if(icount  >=  1) then
          call getarg(1,runname)
       else
          write(6,*) 'What is the name of the *.in input file?'
          write(6,*) '(without the .in suffix):'
          read(5,*) runname
       endif
       lrunname=index(runname,' ')-1
!       if (verbose .ge. 1) write(6,*) 'Running eliteeq for ', runname(1:lrunname)


      nunit=17
      open(unit=nunit,file=runname(1:lrunname)//'.in',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem opening file',runname(1:lrunname)//'.in'
         stop
      endif
      read(nunit,equil)
      read(nunit,qref_modes)
      if (verbose .ge. 1) write(6,*) 'Running eliteeq for ', runname(1:lrunname)
      if (verbose .ge. 2) write(6,equil)
      if (verbose .ge. 2) write(6,qref_modes)
      read(nunit,plas)
      close(nunit)
!-----------------------------------------------------------------------
! remaining initializations
!-----------------------------------------------------------------------
!      q_surf=q0-del/nn ! here is the reason we need to read qref_modes
                                !to determine q of the equil surface

      if (qafix .and. (.not. setdel)) then
         write(*,*) 'qafix option not defined without setdel'
         stop
      endif
      
!  Open the .eqout file for text output
      outfile=39
      if (verbose .ge. 2) then
         open(unit=outfile,file=runname(1:lrunname)//'.eqout',iostat=ios)
         write(outfile,equil)
         write(outfile,qref_modes)
      endif

! set default values of nmvac and nmwinhalf based on nn:
!  6/2020 don't need this, and new options require q_surf and eps
!      if (nmvac.lt.0) then
!         nmvac_inp=nmvac
!         if (nmvac.eq.-1) then
!            nmvac=min(nn+3,12)
!         else if (nmvac.eq.-2) then
!            nmvac=min(nn+2,11)
!         else if (nmvac.eq.-3) then 
!            nmvac=min(nn+4,15)
!         else if (nmvac.eq.-4) then 
!            nmvac=min(int(sqrt(nn+1.)),12)
!         else if (nmvac.eq.-5) then
!            if (nn .lt. 10.) then
!               nmvac=nn+2
!            else
!               nmvac=min(int((nn*2.)**0.6667+6.),50)
!            endif
!         else if (nmvac.eq.-6) then
!            if (nn .lt. 7) then
!               nmvac=nn+1
!            else
!               nmvac=min(nn+2,15)
!            endif
!         else if (nmvac.eq.-15) then 
!            nmvac=max(int(nn*2.5),10)
!         else if (nmvac.eq.-16) then 
!            nmvac=nn*1.5
!         else if (nmvac.eq.-17) then 
!            nmvac=nn*3.5
!         else if (nmvac.eq.-18) then 
!            nmvac=min(int(nn*1.5),30)
!         else if (nmvac.eq.-19) then 
!            nmvac=min(nn*2,30)
!         else
!            write(*,*) ' Entered unsupported nmvac value=',nmvac
!            if (verbose .gt. 3) then
!               write(outfile,*) ' Entered unsupported nmvac value=',nmvac
!            endif
!            stop
!         endif
!         if (verbose .ge. 1) then
!            write(*,*) ' nmvac=',nmvac_inp,' input, set to nmvac=',nmvac
!         endif
!         if (verbose .gt. 3) then
!            write(outfile,*) ' nmvac=',nmvac_inp,' input, set to nmvac=',nmvac
!         endif
!      endif

!      if (nmwinhalf.eq.-1) then
!         nmwinhalf=nmvac+3
!         if (verbose .ge. 1) then
!            write(*,*) 'setting nmwinhalf to default value of ',nmwinhalf
!         endif
!         if (verbose .gt. 3) then
!            write(outfile,*) 'setting nmwinhalf to default value of', &
!                 nmwinhalf
!         endif
!      endif

! make sure npsi is larger than nxinterp for efit files
      if (shape.eq.'eqds'.and.(npsi.le.nxinterp)) then
         npsi=nxinterp+20
         if (verbose .gt. 1) then
            write(*,*) 'increased npsi to be larger than nxinterp, npsi=',npsi
         endif
         if (verbose .gt. 3) then
            write(outfile,*) 'increased npsi to be larger than nxinterp, npsi=',npsi
         endif
      endif


!  Allocate the main set of common arrays
      call alloc_arrays

      return
   end subroutine init

   subroutine alloc_arrays

     allocate( theta(npts),rpts(npts),zpts(npts),drdt(npts), &
       dzdt(npts),bppts(npts),omega(npts),om_edge(npts) )

     theta=0.; rpts=0.; zpts=0.; drdt=0.; dzdt=0.; bppts=0.; omega=0.
     om_edge=0.

   end subroutine alloc_arrays

   subroutine dealloc_arrays

     deallocate(theta,rpts,zpts,drdt,dzdt,bppts,omega,om_edge)

   end subroutine dealloc_arrays

   subroutine alloc_eqdsk

     allocate( xgrid(nx),zgrid(nz) )
     xgrid=0.; zgrid=0.

     allocate( psixz(nx,nz),bpsq(nx,nz) )
     psixz=0.; bpsq=0.

     allocate( sp(nx),spp(nx),sf(nx),sffp(nx),qpsi(nx) )
     sp=0.; spp=0.; sf=0.; sffp=0.; qpsi=0.

     allocate( pressw(nx),pwprim(nx),rho0(nx),rho0p(nx) )
     pressw=0.; pwprim=0.; rho0=0.; rho0p=0.

   end subroutine alloc_eqdsk

   subroutine alloc_bndry

     allocate( xbndry(nbndry),zbndry(nbndry) )
     xbndry=0.; zbndry=0.

     allocate( xlim(nlim),zlim(nlim) )
     xlim=0.; zlim=0.

   end subroutine alloc_bndry

   subroutine alloc_mapper

     allocate( psic(npsi),psiv(npsi),pprime(npsi),fval(npsi), &
           ffprime(npsi),chipsi(npsi),qsfin(npsi),ppp(npsi), &
           ffpp(npsi),qsfinp(npsi),qsfinpp(npsi),pressval(npsi),qcalc(npsi) )
     psic=0.; psiv=0.; pprime=0.; fval=0.; qsfinp=0.; qsfinpp=0.
     ffprime=0.; chipsi=0.; qsfin=0.; ppp=0.; ffpp=0.; pressval=0.; qcalc=0.
     
     allocate( xs(npsi,npts),zs(npsi,npts),bps(npsi,npts) )
     xs=0.; zs=0.; bps=0.

     allocate( press(npsi),pw(npsi),pwp(npsi),rho(npsi),rhop(npsi) )
     press=0.; pw=0.; pwp=0.; rho=0.; rhop=0.

     allocate( csplpsi(4,nx,nz) )
     csplpsi=0.

     allocate( arcsur(npsi,npts) )  ! used in equalarc
     arcsur=0.

     if (dens) then
        allocate( nel(npsi),nprime(npsi),tel(npsi),teprime(npsi),  &
             tion(npsi),tiprime(npsi) )
        nel=0.; nprime=0.; tel=0.; teprime=0.; tion=0.; tiprime=0.
     endif

   end subroutine alloc_mapper
   

   subroutine dealloc_mapper

      deallocate(psic,psiv,pprime,fval,ffprime,chipsi,qsfin,ppp, &
           ffpp,qsfinp,qsfinpp,pressval,qcalc)
      deallocate(xs,zs,bps)
      deallocate(press,pw,pwp,rho,rhop,csplpsi,arcsur)

      if (dens) deallocate( nel,nprime,tel,teprime,tion,tiprime )

   end subroutine dealloc_mapper

   subroutine alloc_toq

     allocate( psiv(npsi),pprime(npsi),fval(npsi), &
           ffprime(npsi),chipsi(npsi),qsfin(npsi),ppp(npsi), &
           press(npsi),ffpp(npsi),qsfinp(npsi),qsfinpp(npsi), &
           pressval(npsi),qcalc(npsi) )
     psiv=0.; pprime=0.; fval=0.; ffpp=0.; qsfinp=0.; qsfinpp=0.
     ffprime=0.; chipsi=0.; qsfin=0.; ppp=0.; press=0.; pressval=0.; qcalc=0.
!
     allocate( nel(npsi),nprime(npsi),tel(npsi),teprime(npsi),  &
               tion(npsi),tiprime(npsi) )
     nel=0.; nprime=0.; tel=0.; teprime=0.; tion=0.; tiprime=0.
     
     allocate( xs(npsi,npts),zs(npsi,npts),bps(npsi,npts),arcsur(npsi,npts) )
     xs=0.; zs=0.; bps=0.; arcsur=0.

   end subroutine alloc_toq

   subroutine dealloc_toq

     deallocate(psiv,pprime,fval,ffprime,chipsi,qsfin,ppp,press, &
          ffpp,qsfinp,qsfinpp,pressval,qcalc)
     deallocate( nel,nprime,tel,teprime,tion,tiprime )
     deallocate(xs,zs,bps,arcsur)

   end subroutine dealloc_toq

end module eliteeq_data









