      
subroutine wrteigfunc

!-----------------------------------------------------------------------
! store eigenfunction in runname.plas
!  Same version works for updown symmetric or antisymmetric
! rlm--probably need to patch up by passing nmodes as well as nm
!-----------------------------------------------------------------------
      use elite_data, only: runname,lrunname,cdatetime,codver, &
           nowindow,nmwinhalf,meshtype,mmax,nd,ndist,nmodes,nn,ns, &
           nxinterp,m0,nx,aspect,btnorm,del,dmercier,dx,dw,etai, &
           lamdist,ne,ppmult,q0,qref,rnorm,tau,te,vloop,qmin,zeff, &
           gam,qa,xx,ffun,psixx,psimin,gamsq,gamscl2,jedgen1store, &
           lhspeelstore,rhspeelstore,alpha,shear,pprime_eq,tmatasym, &
           nmlow,dmstore,rotation,psigrid,rotnorm,gamim,mind_new,gamre, &
           verbose
      implicit none
      integer i,m
      integer ios,nunit
! variables for namelist equil
!      real k
!      real gamma
!!      real rmajor
!      real rminor
!      real f_surf
!      real kappa
!      real triang
!      real kappap
!      real triangp
!      real deltap
      integer npts
      character*4 shape
      integer npsi              !size of mapped psi mesh
      real percenflux           !outermost flux surface is 
                                !psiv(npsi)=psiaxis+(psilim-psiaxis)*percenflux
      real alpsi,del_fix
      logical delmin,setdel,qafix
      integer nnmin,nnmax
      namelist/equil/npts,shape,&
          percenflux,alpsi,npsi, &
          delmin,nnmin,nnmax,setdel,del_fix,qafix
!      namelist/equil/k,gamma,rminor,rmajor,npts,f_surf,shape,&
!          kappa,triang,kappap,triangp,deltap,percenflux,alpsi,npsi, &
!          delmin,nnmin,nnmax,setdel,del_fix,qafix
!      namelist/equil/k,gamma,rminor,rmajor,npts,f_surf,shape,&
!          kappa,triang,kappap,triangp,deltap,percenflux,alpsi,npsi, &
!          delmin,nnmin,nnmax
!      namelist/equil/k,gamma,rminor,rmajor,npts,f_surf,shape, &
!          kappa,triang,kappap,triangp,deltap
      integer ng
      namelist/vac/npts,ng !beware npts is in two different namelists
      integer nmwinhlf
      integer npts_equil
      real rotnormxx(nx),rotspline(nxinterp),Gamr,Gamimg,Gamimgnorm
      complex, dimension(:,:), allocatable :: gfun
      complex Gam_c
!      real gfun(nmodes,nx)
      open(unit=12,file=runname(1:lrunname)//'.in',status='old',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'could not open ',runname(1:lrunname)//'.in in wrteigfunc'
         stop
      endif
      read(12,equil)
      npts_equil=npts
      read(12,vac)
      close(12)
      nunit=30
      open(unit=nunit,file=runname(1:lrunname)//'.plas', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.plas'
         stop
      endif
      write(nunit,*) cdatetime
      write(nunit,*) "code version: ",codver
! write the equil namelist quantities first
      write(nunit,'(a4)') shape
      write(nunit,'(i6)') npts_equil
! 12 reals
      write(nunit,'(1p,6e13.5)') psimin,gamre*sqrt(gamscl2),gamsq, &
           jedgen1store,lhspeelstore,rhspeelstore,alpha(1),shear(1), &
           pprime_eq(1),tmatasym,percenflux,alpsi
!      write(nunit,'(1p,6e13.5)') k,gamma,rminor,rmajor,f_surf, &
!          kappa,triang,kappap,triangp,deltap,percenflux,alpsi
!     namelist & calculated quantities
      if(nowindow) then
         nmwinhlf=-1
      else
         nmwinhlf=nmwinhalf
      endif
      write(nunit,'(13i6)') meshtype,mmax, &
          nmlow,ndist,nmodes*mind_new,nn,ns,nxinterp
      write(nunit,'(5i6)') m0,nx,nmwinhlf,npts,ng 
      write(nunit,'(1p,6e13.5)') dmstore,btnorm,del, &
          dmercier,dx,dw,etai,lamdist,ne,ppmult, &
          q0,qref,rnorm,tau,te,vloop,qmin,zeff
      write(nunit,'(1p,6e13.5)') gam,qa ! calculated
      write(nunit,'(1p,6e13.5)') (xx(i),i=1,nx)
      if (rotation) then
   ! need to renormalize, to output X=Gamma*F rather than F eigfct
   !  spline rotation onto xx grid
         allocate( gfun(nmodes,nx) )
         call spline(-psigrid,rotnorm,nxinterp,-1.d30,-1.d30,rotspline)
         do i=1,nx
            call splint(-psigrid,rotnorm,rotspline, &
                  nxinterp,-psixx(i),rotnormxx(i))
            if (verbose .ge. 1) write(85,*) 'i=',i,' rotnormxx=',rotnormxx(i)
         enddo
!         Gamr=sqrt(gamsq)
         Gamr=gamre
!         do m=1,nmodes
!            if (mod(m,2)==1) then  ! real parts 
               ! g_r=f_r-gamimgnorm*f_i
         do i=1,nx
            Gamimgnorm=(gamim+nn*rotnormxx(i))/Gamr
            Gam_c=cmplx(1.,Gamimgnorm)
            gfun(:,i)=Gam_c*ffun(:,i)
!                  gfun(m,i)=ffun(m,i)-Gamimgnorm*ffun(m+1,i)
         enddo
!         enddo
         ffun=gfun  ! renormalizing ffun to get proper X eigenfunction
              ! in fun2d as well as here
      endif
      if (mind_new.eq.1) then
         do m=1,nmodes
            write(nunit,'(1p,6e13.5)') (real(ffun(m,i)),i=1,nx)
         enddo
      else
         do  m=1,nmodes
            write(nunit,'(1p,6e13.5)') (real(ffun(m,i)),i=1,nx)
            write(nunit,'(1p,6e13.5)') (aimag(ffun(m,i)),i=1,nx)
         end do
      endif
! 9/01 add psixx because xx may not be monotonic
      write(nunit,'(1p,6e13.5)') (psixx(i),i=1,nx)
      close(nunit)

    return
      
end subroutine wrteigfunc
