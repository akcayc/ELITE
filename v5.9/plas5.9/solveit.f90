
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ELITE (Edge Localized Instabilities in Tokamak Experiments)
!
! Authors:  P.B. Snyder, H.R. Wilson, R.L. Miller
!
! Copyright(c) 1999-2014, All rights reserved
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



subroutine solveit

  ! translated to F90 1/31/00
  !-----------------------------------------------------------------------
  ! this routine contains main logic of the algorithm
  !  Same version works for updown symmetric or antisymmetric
  !-----------------------------------------------------------------------
  use elite_data, only: mwindow,dx,nd,meshtype,lamdist,ndist, &
       outfile,xmin,xmax,xx,nx,kx,del,dw,xinterp,igam,nedge, &
       gamsq,aspect,qref,nn,gamscl,alpha,shear,ppmult,alphamult, &
       shearmult,autorun,runname,lrunname,gam,aval,auarr,iuarr, &
       qmonotonic,nxinterp,psixx,psigrid,splineeq,tmatasym,pi, &
       ne_eq,wsmodel,omegas,pprime_eq,rgrid,rmajor,circumf,q_eq, &
       rotation,gamim,dens,ion_mass,rot_eq,mind_new,gamre, &
       rav,cssq1,cssq2,rva1,rva2,psiv_eq,funcal,gamscl2, &
       rot_const,newtstep,omegasn,omegaspi,omegas_max,omegasn_max, &
       omegaspi_max,ncrit_ws,verbose
  implicit none
!  real vr(mwindow,mwindow),vi(mwindow,mwindow)
  complex vcomp(mwindow,mwindow)  ! complex array for eigenvectors
  complex vcdum(mwindow,mwindow)
  real gamr(mind_new*mwindow),gami(mind_new*mwindow)
  real gamr2(mind_new*mwindow),gami2(mind_new*mwindow)
  real dxmin,gammin,gammin2
  real sigr,sigi,dwid,st,stp,alph,dedge
  real yold,xold,errgam,errnow,gamnew,gamax,gamrtol
  integer i,m,ni,ilab,nunit,ios,j,igmax,interp,nxpersurf
  real niws(nxinterp),niwsspline(nxinterp),niwsmax,niwsloc
  real ppspline(nxinterp),ppmax,pploc,rppmax,rpphalf
  real rspline(nxinterp),delta,qppmax,qspline(nxinterp)
  real circum_approx,ktheta,rminor 
  integer ippmax,ipphalf
  real gamimnew,gnmax,gamnorm(mind_new*mwindow),x2old,y2old,errimnow
  integer ignmax,ii
  real rot_norm
  integer lilabarg,iargc,icount
  character*6 ilabarg
  real wsspline(nxinterp),wsloc,psiwsmax,omegas_eff(nxinterp)
  real rpsispline(nxinterp),psirot,matchrot,maxfrac
  integer unst_count    ! when igam=1, =0 stable, =2 unstable, =1 ambiguous
  !-----------------------------------------------------------------------      
  !     1.0 generate x-mesh
  !-----------------------------------------------------------------------
  ! 9/01 first check for monotonicity of q profile - different gridding options
  xold=xinterp(1)
  qmonotonic=.true.
  do i=2,nxinterp
     if (xinterp(i) .le. xold) then
        qmonotonic=.false.
     endif
     xold=xinterp(i)
  enddo
  if (.not.qmonotonic) then
     write(*,*) 'q profile is not monotonic'
     if (verbose .ge. 2) write(outfile,*) 'q profile is not monotonic'
  endif

  xold=0.
  dxmin=dx/nd
!  5/14 pbs: New meshtype=6 option which will revert to type 1 if q is
!   monotonic and type 5 if it is not
  if ((meshtype .eq. 6) .and. qmonotonic) then
     meshtype=1
     if (verbose .ge. 2) write(*,*) 'meshtype 6 changed to 1 as q is monotonic'
     if (verbose .ge. 2) write(outfile,*) 'meshtype 6 changed to 1 as q is monotonic'
  endif
  if ((meshtype .eq. 6) .and. (.not. qmonotonic)) then
     meshtype=5
     if (verbose .ge. 2) write(*,*) 'meshtype 6 changed to 5 as q is non-monotonic'
     if (verbose .ge. 2) write(outfile,*) 'meshtype 6 changed to 5 as q is non-monotonic'
     if ((gamsq .lt. 5.) .and. (.not. rotation)) then
        gamsq=2500.
        if (verbose .ge. 2) write(*,*) 'gamsq set to large initial guess value for q reversed case, gamsq=',gamsq
        if (verbose .ge. 2) write(outfile,*) 'gamsq set to large initial guess value for q reversed case, gamsq=',gamsq
     endif
  endif

  if (meshtype.eq.1) then
     if (.not.qmonotonic) then
        write(*,*) 'meshtype=1 not supported for non-monotonic q'
        stop
     endif
     dxmin=lamdist*(200./ndist)*0.0004
     if (verbose .ge. 3) write(outfile,*) 'dxmin=',dxmin
     call getxx(ndist,dxmin,xmin,xmax,xx,nx)
!     write(*,*) 'return from getxx'
     if (nx.ne.kx) then
        if (verbose .ge. 3) then
           write(outfile,*) 'nx .ne. kx',nx,kx
           write(outfile,*) 'reallocating xx and psixx'
        endif
        deallocate(xx,psixx)
        kx=nx
        allocate( xx(kx), psixx(kx) )
        xx=0.; psixx=0.
        call getxx(ndist,dxmin,xmin,xmax,xx,nx)

        !            stop
     endif

  else if (meshtype.eq.0) then
     if (.not.qmonotonic) then
        write(*,*) 'meshtype=0 not supported for non-monotonic q'
        stop
     endif
     i=1
     xx(1)=del
30   continue
     st=dx
     m=int(xx(i)+dx/nd)
!     dwid=abs(float(m)-xx(i))
     dwid=abs(dble(m)-xx(i))
     if (dwid.gt.0.5) dwid=abs(1.-dwid)
     if (dwid.lt.dw/2.) st=dx/nd
     if ((xx(i)-del).lt.dw/2.) st=dx/nd
     !     write(6,*)' m=',m,' x=',xx(i),' dwid=',dwid,' dw=',dw/2.
     i=i+1
     if (i.gt.kx) then
        write(6,*)' Input error***dx too small for array size'
        write(6,*)' Max achievable x for this mesh=',xx(i-1),' i=',i
        write(6,*)' Max x asked for through inputs=',xmax
        stop
     end if
     xx(i)=xx(i-1)+st
     !     write(6,*)' i=',i,' x=',xx(i)
     if (xx(i).lt.xmax) goto 30
     nx=i-1
  else if (meshtype.eq.2) then
     if (.not.qmonotonic) then
        write(*,*) 'meshtype=2 not supported for non-monotonic q'
        stop
     endif
     !  equal spaced mesh, using ndist as total number of x-mesh points...
     st=(xmax-del)/(ndist-1.)
     nx=ndist
     do i=1,nx
        xx(i)=del+(i-1)*st
     end do
     ! can't have last point exactly on last xinterp surface
     if (xx(nx) .ge. xinterp(nxinterp) ) then
        xx(nx)=0.99*xx(nx)+0.01*xx(nx-1)
     endif
  else if (meshtype == 3) then
     ! nieve mesh, meshing evenly in x between each ixinterp surface
     if (verbose .ge. 3) write(*,*) &
          'meshtype=3, nieve grid with ndist evenly spaced', &
          ' linear points between each of the nxinterp equil surfaces'
     nx=(nxinterp-1)*ndist
     do i=2,nxinterp-1
        st=(xinterp(i)-xinterp(i-1))/(ndist)
        stp=(psigrid(i)-psigrid(i-1))/(ndist)
        do j=1,ndist
           xx((i-2)*ndist+j)=xinterp(i-1)+(j-1)*st
           psixx((i-2)*ndist+j)=psigrid(i-1)+(j-1)*stp
        enddo
     enddo
     !  for last surface, grid evenly to xinterp(nxinterp), include end point
     !  9/28 - can't include end point as this messes up interpolation in matgen
     st=(xinterp(nxinterp)-xinterp(nxinterp-1))/(ndist-1.)
     stp=(psigrid(nxinterp)-psigrid(nxinterp-1))/(ndist-1.)
     do j=1,ndist-1
        xx((nxinterp-2)*ndist+j)=xinterp(nxinterp-1)+(j-1)*st
        psixx((nxinterp-2)*ndist+j)=psigrid(nxinterp-1)+(j-1)*stp
     enddo
     xx((nxinterp-2)*ndist+ndist)=xinterp(nxinterp-1)+(ndist-1.01)*st
     psixx((nxinterp-2)*ndist+ndist)=psigrid(nxinterp-1)+(ndist-1.01)*stp
     !         write(*,*) 'xx=',xx
     !         write(*,*) 'psixx=',psixx

!!$  else if (meshtype == 4) then
!!$     ! attempt to evenly space in x even for non-monotonic q profiles
!!$     !  should be similar to meshtype=2, use ndist total points
!!$     !  start with psii at outermost point, then grid inward until
!!$     !   psigrid(nxinterp) is reached, adjust dxguess until the total
!!$     !   number of points needed = ndist
!!$     call getxx4
!!$     write(24,*) 'xx=',xx
!!$     write(24,*) 'psixx=',psixx
!!$
!!$     !        stop
  else if (meshtype == 5) then
  ! nieve mesh evenly spaced between ixinterp surfaces, but with
  !  a total number of points = ndist*number of rational surfaces
     nxpersurf=ndist*xmax/(nxinterp-1)
     nx=nxpersurf*(nxinterp-1)
     if (verbose .ge. 3) write(*,*) &
          'meshtype=5, nieve grid with ',nxpersurf,' evenly spaced', &
          ' linear points between each of the nxinterp equil surfaces'
     if (verbose .ge. 3) write(*,*) 'nx=',nx
     do i=2,nxinterp-1
        st=(xinterp(i)-xinterp(i-1))/(nxpersurf)
        stp=(psigrid(i)-psigrid(i-1))/(nxpersurf)
        do j=1,nxpersurf
           xx((i-2)*nxpersurf+j)=xinterp(i-1)+(j-1)*st
           psixx((i-2)*nxpersurf+j)=psigrid(i-1)+(j-1)*stp
        enddo
     enddo
     !  for last surface, grid evenly to xinterp(nxinterp), include end point
     !  9/28 - can't include end point as this messes up interpolation in matgen
     st=(xinterp(nxinterp)-xinterp(nxinterp-1))/(nxpersurf-1.)
     stp=(psigrid(nxinterp)-psigrid(nxinterp-1))/(nxpersurf-1.)
     do j=1,nxpersurf-1
        xx((nxinterp-2)*nxpersurf+j)=xinterp(nxinterp-1)+(j-1)*st
        psixx((nxinterp-2)*nxpersurf+j)=psigrid(nxinterp-1)+(j-1)*stp
     enddo
     xx((nxinterp-2)*nxpersurf+nxpersurf)=xinterp(nxinterp-1)+(nxpersurf-1.01)*st
     psixx((nxinterp-2)*nxpersurf+nxpersurf)=psigrid(nxinterp-1)+(nxpersurf-1.01)*stp
  else
     write(6,*)' meshtype=',meshtype,' is not a permissible choice'
     stop
  endif

  ! 9/01 now need to map x grid onto psi for the monotonic grid choices, 
  !   assume linear relation between x and psi between interp surfaces
  if (meshtype .le. 2) then
     interp=1
     do i=1,nx
        do while (xx(i) > xinterp(interp+1))
           interp=interp+1
           if (interp .ge. nxinterp) then
              write(*,*) 'problem mapping xx onto psi', &
                   '  interp out of range',interp, &
                   ' xx=',xx(i),'xinterp(nxinterp=',xinterp(nxinterp)
              stop
           endif
        enddo
        alph=(xx(i)-xinterp(interp))/(xinterp(interp+1)-xinterp(interp))
        if ( (alph>1.) .or. (alph<0.)) then
           write(*,*) 'problem mapping onto xx onto psi'
           write(*,*) 'i=',i,'alph=',alph,'interp=',interp
           stop
        endif
        psixx(i)=(1.-alph)*psigrid(interp)+(alph)*psigrid(interp+1)
     enddo
 !    write(24,*) 'xx=',xx
 !    write(24,*) 'psixx=',psixx
  endif


 ! do i=1,nx
 !    write(996,*)' i=',i,' x=',xx(i)
 ! end do
  if (verbose .ge. 3) then 
     write(outfile,*)' nx=',nx,' kx=',kx
     write(6,*)' nx=',nx,' meshtype=',meshtype
  endif
  !      do  ni=1,nalpha
  !         if (nalpha.gt.1) then
  !            alpha=almin+(ni-1)*dalf
  !            call refsurf        ! calculate reference surface
  !         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 5/29/02  add options for setting gamsq according to omegas
  !       for stability check
  !  make this a separate subroutine?

!  6/06 calculate omegastar related quantities for all cases with
!       density information
  if (dens) then
     call spline(-psigrid,omegas,nxinterp,-1.d30,-1.d30,wsspline)
     omegas_max=0.
     do i=1,nx
        call splint(-psigrid,omegas,wsspline,nxinterp,-psixx(i), &
             wsloc)
        if (abs(wsloc) > omegas_max) then
           omegas_max=abs(wsloc)
           psiwsmax=psixx(i)
        endif
     enddo
     if (verbose .ge. 2) then
        write(*,*) 'omegas_max=',omegas_max,' at psi=',psiwsmax
        write(outfile,*) 'omegas_max=',omegas_max,' at psi=',psiwsmax
     endif
     call spline(-psigrid,omegasn,nxinterp,-1.d30,-1.d30,wsspline)
     omegasn_max=0.
     psiwsmax=0.
     do i=1,nx
        call splint(-psigrid,omegasn,wsspline,nxinterp,-psixx(i), &
             wsloc)
        if (abs(wsloc) > omegasn_max) then
           omegasn_max=abs(wsloc)
           psiwsmax=psixx(i)
        endif
     enddo
     if (verbose .ge. 3) then
        write(*,*) 'omegasn_max=',omegasn_max,' at psi=',psiwsmax
        write(outfile,*) 'omegasn_max=',omegasn_max,' at psi=',psiwsmax
     endif
     call spline(-psigrid,omegaspi,nxinterp,-1.d30,-1.d30,wsspline)
     omegaspi_max=0.
     psiwsmax=0.
     do i=1,nx
        call splint(-psigrid,omegaspi,wsspline,nxinterp,-psixx(i), &
             wsloc)
        if (abs(wsloc) > omegaspi_max) then
           omegaspi_max=abs(wsloc)
           psiwsmax=psixx(i)
        endif
     enddo
     if (verbose .ge. 2) then
        write(*,*) 'omegaspi_max=',omegaspi_max,' at psi=',psiwsmax
        write(outfile,*) 'omegaspi_max=',omegaspi_max,' at psi=',psiwsmax
     endif
  endif

! wsmodels don't make sense with rotation/complex gamma
  if (wsmodel .ne. 0 .and. rotation) then
     write(*,*) 'wsmodel=',wsmodel,' not compatible with', &
          'rotation=.true.'
     stop
  endif
     
!! 4/06 PBS found bug for wsmodel cases with varying density
!!   gamsq is defined in terms of edge density, not local density
!!   so should be maximizing over omegas^2 not (ni * omegas^2)
!!   gamsq=4 pi ion_mass 1.6726e-24 ne_eq(1) gamma^2
!!   gamma > omegas/2 ->  gamma^2 > omegas^2/4
!!   gamsq > pi ion_mass 1.6726e-24 ne_eq(1) omegas^2
  if (wsmodel.eq.1) then    
     !       use half of max value of omega_*_p_i**2 to determine gamsq
     !         spline onto the fine grid to find maximum omegas**2*ne
     niws=omegas**2*ne_eq(1) 
     call spline(-psigrid,niws,nxinterp,-1.d30,-1.d30,niwsspline)
     niwsmax=0.
     do i=1,nx
        call splint(-psigrid,niws,niwsspline,nxinterp,-psixx(i), &
             niwsloc)
        if (niwsloc > niwsmax) niwsmax=niwsloc
     enddo
     gamsq=pi*ion_mass*1.6726e-24*(0.5*niwsmax) ! 4 pi m_i 1/4 n_i oms^2
     if (verbose .ge. 3) then
        write(*,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(outfile,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(*,*) 'for wsmodel=1 using half this value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        write(outfile,*) 'for wsmodel=1 using half this value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
     
        write(*,*) 'assigned gamsq=',gamsq
        write(outfile,*) 'assigned gamsq=',gamsq
     endif
  endif

  if (wsmodel.eq.6) then    
     !       use 1/4 of max value of omega_*_p_i**2 to determine gamsq
     !         spline onto the fine grid to find maximum omegas**2*ne
     niws=omegas**2*ne_eq(1) 
     call spline(-psigrid,niws,nxinterp,-1.d30,-1.d30,niwsspline)
     niwsmax=0.
     do i=1,nx
        call splint(-psigrid,niws,niwsspline,nxinterp,-psixx(i), &
             niwsloc)
        if (niwsloc > niwsmax) niwsmax=niwsloc
     enddo
     gamsq=pi*ion_mass*1.6726e-24*(0.25*niwsmax) ! 4 pi m_i 1/4 n_i oms^2
     if (verbose .ge. 3) then
        write(*,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(outfile,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(*,*) 'for wsmodel=6 using one fourth this value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        write(outfile,*) 'for wsmodel=6 using one fourth this value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
     
        write(*,*) 'assigned gamsq=',gamsq
        write(outfile,*) 'assigned gamsq=',gamsq
     endif

  endif
  if (wsmodel.eq.7) then    
     !       use 1/8 of max value of omega_*_p_i**2 to determine gamsq
     !         spline onto the fine grid to find maximum omegas**2*ne
     niws=omegas**2*ne_eq(1) 
     call spline(-psigrid,niws,nxinterp,-1.d30,-1.d30,niwsspline)
     niwsmax=0.
     do i=1,nx
        call splint(-psigrid,niws,niwsspline,nxinterp,-psixx(i), &
             niwsloc)
        if (niwsloc > niwsmax) niwsmax=niwsloc
     enddo
     gamsq=pi*ion_mass*1.6726e-24*(0.125*niwsmax) ! 4 pi m_i 1/4 n_i oms^2
     if (verbose .ge. 3) then
        write(*,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(outfile,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(*,*) 'for wsmodel=7 using one eighth this value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        write(outfile,*) 'for wsmodel=7 using one eighth this value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        gamsq=pi*ion_mass*1.6726e-24*(0.125*niwsmax) ! 4 pi m_i 1/4 n_i oms^2
        write(*,*) 'assigned gamsq=',gamsq
        write(outfile,*) 'assigned gamsq=',gamsq
     endif
  endif

  if (wsmodel.eq.8) then    
     !       use 1/4 of max value of omega_*_n_i**2 to determine gamsq
     !         spline onto the fine grid to find maximum omegas**2*ne
     !      (like wsmodel 6 except use omegesn rather than omegas
     niws=omegasn**2*ne_eq(1) 
     call spline(-psigrid,niws,nxinterp,-1.d30,-1.d30,niwsspline)
     niwsmax=0.
     do i=1,nx
        call splint(-psigrid,niws,niwsspline,nxinterp,-psixx(i), &
             niwsloc)
        if (niwsloc > niwsmax) niwsmax=niwsloc
     enddo
     gamsq=pi*ion_mass*1.6726e-24*(0.25*niwsmax) ! 4 pi m_i 1/4 n_i oms^2
     if (verbose .ge. 3) then
        write(*,*) 'maximum ni*omegasn^2 value is ',niwsmax
        write(outfile,*) 'maximum ni*omegasn^2 value is ',niwsmax
        write(*,*) 'for wsmodel=8 using one fourth this value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        write(outfile,*) 'for wsmodel=8 using one fourth this value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
     
        write(*,*) 'assigned gamsq=',gamsq
        write(outfile,*) 'assigned gamsq=',gamsq
     endif
  endif

  if (wsmodel.eq.9) then    
     !       use 1/32 of max value of omega_*_p_i**2 to determine gamsq
     !         spline onto the fine grid to find maximum omegas**2*ne
     !  8/17/09 pbs  This is meant to be a minimal omegastar model motivated
     !    by a mode width twice the pedestal width (hence use omegasmax/4)
     !    and by the Hastie/Catto/Ramos result of around sqrt(2) weakening
     !    of omegastar stabilization due to profiles of omegastar
     !   Hence use omegastar_max/(4*sqrt(2))  or omegastar_max^2/32
     niws=omegas**2*ne_eq(1) 
     call spline(-psigrid,niws,nxinterp,-1.d30,-1.d30,niwsspline)
     niwsmax=0.
     do i=1,nx
        call splint(-psigrid,niws,niwsspline,nxinterp,-psixx(i), &
             niwsloc)
        if (niwsloc > niwsmax) niwsmax=niwsloc
     enddo
     gamsq=pi*ion_mass*1.6726e-24*(0.03125*niwsmax) ! 4 pi m_i 1/4 n_i oms^2
     if (verbose .ge. 3) then
        write(*,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(outfile,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(*,*) 'for wsmodel=9 using 1/32nd of this value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        write(outfile,*) 'for wsmodel=9 using 1/32nd of this value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
     
        write(*,*) 'assigned gamsq=',gamsq
        write(outfile,*) 'assigned gamsq=',gamsq
     endif
  endif

  if (wsmodel.eq.10) then    
     !  pbs 4/2010:  bi-linear model based on BOUT++ results
     !       use 1/2 of max value of omega_*_p_i**2 to determine gamsq
     !        up to nn=ncrit_ws, then for nn>ncrit_ws increases by only 6%
     !           that is astar=omega_star/nn
     !         nn<=ncrit_ws:  omega_star_eff = omegas
     !         nn>ncrit_ws:   omega_star_eff = astar*ncrit_ws + (nn-ncrit_ws)*.06*astar
     !         spline onto the fine grid to find maximum omegas**2*ne
     if (nn .le. ncrit_ws) then
        omegas_eff=omegas
     else
        omegas_eff=omegas*ncrit_ws/nn + 0.168*(nn-ncrit_ws)*omegas/nn
     endif
     niws=omegas_eff**2*ne_eq(1) 
     call spline(-psigrid,niws,nxinterp,-1.d30,-1.d30,niwsspline)
     niwsmax=0.
     do i=1,nx
        call splint(-psigrid,niws,niwsspline,nxinterp,-psixx(i), &
             niwsloc)
        if (niwsloc > niwsmax) niwsmax=niwsloc
     enddo
     gamsq=pi*ion_mass*1.6726e-24*(0.5*niwsmax) ! 4 pi m_i 1/2 n_i oms^2
     if (verbose .ge. 3) then
        write(*,*) 'maximum ni*omegas_eff^2 value is ',niwsmax
        write(outfile,*) 'maximum ni*omegas_eff^2 value is ',niwsmax
        write(*,*) 'for wsmodel=10 using one half this value (rollover model) ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        write(outfile,*) 'for wsmodel=10 using one half this value (rollover model) ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
     
        write(*,*) 'assigned gamsq=',gamsq
        write(outfile,*) 'assigned gamsq=',gamsq
     endif
  endif

  if (wsmodel.eq.11) then    
     !  pbs 4/2010:  bi-linear model based on BOUT++ results
     !       use 1/4 of max value of omega_*_p_i**2 to determine gamsq
     !        up to nn=ncrit_ws, then for nn>ncrit_ws increases by only 6%
     !           that is astar=omega_star/nn
     !         nn<=ncrit_ws:  omega_star_eff = omegas
     !         nn>ncrit_ws:   omega_star_eff = astar*ncrit_ws + (nn-ncrit_ws)*.06*astar
     !         spline onto the fine grid to find maximum omegas**2*ne
     if (nn .le. ncrit_ws) then
        omegas_eff=omegas
     else
        omegas_eff=omegas*ncrit_ws/nn + 0.06*(nn-ncrit_ws)*omegas/nn
     endif
     niws=omegas_eff**2*ne_eq(1) 
     call spline(-psigrid,niws,nxinterp,-1.d30,-1.d30,niwsspline)
     niwsmax=0.
     do i=1,nx
        call splint(-psigrid,niws,niwsspline,nxinterp,-psixx(i), &
             niwsloc)
        if (niwsloc > niwsmax) niwsmax=niwsloc
     enddo
     gamsq=pi*ion_mass*1.6726e-24*(0.25*niwsmax) ! 4 pi m_i 1/4 n_i oms^2
     if (verbose .ge. 3) then
        write(*,*) 'maximum ni*omegas_eff^2 value is ',niwsmax
        write(outfile,*) 'maximum ni*omegas_eff^2 value is ',niwsmax
        write(*,*) 'for wsmodel=11 using one fourth this value (rollover model) ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        write(outfile,*) 'for wsmodel=10 using one fourth this value (rollover model) ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        
        write(*,*) 'assigned gamsq=',gamsq
        write(outfile,*) 'assigned gamsq=',gamsq
     endif
  endif

  if (wsmodel.eq.12) then    
     !  pbs 4/2010:  bi-linear model based on BOUT++ results
     !       use .84^2~0.7 of max value of omega_*_p_i**2 to determine gamsq
     !        up to nn=ncrit_ws, then for nn>ncrit_ws increases by only 6%
     !           that is astar=omega_star/nn
     !         nn<=ncrit_ws:  omega_star_eff = omegas
     !         nn>ncrit_ws:   omega_star_eff = astar*ncrit_ws + (nn-ncrit_ws)*.06*astar
     !         spline onto the fine grid to find maximum omegas**2*ne
     if (nn .le. ncrit_ws) then
        omegas_eff=omegas
     else
        omegas_eff=omegas*ncrit_ws/nn + 0.06*(nn-ncrit_ws)*omegas/nn
     endif
     niws=omegas_eff**2*ne_eq(1) 
     call spline(-psigrid,niws,nxinterp,-1.d30,-1.d30,niwsspline)
     niwsmax=0.
     do i=1,nx
        call splint(-psigrid,niws,niwsspline,nxinterp,-psixx(i), &
             niwsloc)
        if (niwsloc > niwsmax) niwsmax=niwsloc
     enddo
     gamsq=pi*ion_mass*1.6726e-24*(0.7*niwsmax) ! 4 pi m_i 0.7 n_i oms^2
     if (verbose .ge. 3) then
        write(*,*) 'maximum ni*omegas_eff^2 value is ',niwsmax
        write(outfile,*) 'maximum ni*omegas_eff^2 value is ',niwsmax
        write(*,*) 'for wsmodel=12 using 0.7 this value (rollover model) ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        write(outfile,*) 'for wsmodel=12 using 0.7 this value (rollover model) ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
     
        write(*,*) 'assigned gamsq=',gamsq
        write(outfile,*) 'assigned gamsq=',gamsq
     endif
  endif

  if (wsmodel.eq.13) then
     !  PS  7/27/2011 add combination model with 0.02w_A + weak w_s threshold    
     !       use 1/16 of max value of omega_*_p_i**2 +.02w_A to determine gamsq
     !         spline onto the fine grid to find maximum omegas**2*ne
     niws=omegas**2*ne_eq(1) 
     call spline(-psigrid,niws,nxinterp,-1.d30,-1.d30,niwsspline)
     niwsmax=0.
     do i=1,nx
        call splint(-psigrid,niws,niwsspline,nxinterp,-psixx(i), &
             niwsloc)
        if (niwsloc > niwsmax) niwsmax=niwsloc
     enddo
     gamsq=pi*ion_mass*1.6726e-24*(0.0625*niwsmax) ! 4 pi m_i 1/4 n_i oms^2
     gamsq=( sqrt(gamsq)+ sqrt(((0.02)**2)/gamscl2) )**2   ! add .02omega_A to this
     if (verbose .ge. 3) then
        write(*,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(outfile,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(*,*) 'for wsmodel=13 using 1/16 this value +.02w_A ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        write(outfile,*) 'for wsmodel=6 using one fourth this value +.02w_A ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
     
        write(*,*) 'assigned gamsq=',gamsq
        write(outfile,*) 'assigned gamsq=',gamsq
     endif
  endif

  if (wsmodel.eq.2) then    
     !       use the max value of omega_*_p_i**2 to determine gamsq
     !         spline onto the fine grid to find maximum omegas**2*ne
     !          this provides something of an upper limit to w* stabilization
     niws=omegas**2*ne_eq(1)
     call spline(-psigrid,niws,nxinterp,-1.d30,-1.d30,niwsspline)
     niwsmax=0.
     do i=1,nx
        call splint(-psigrid,niws,niwsspline,nxinterp,-psixx(i), &
             niwsloc)
        if (niwsloc > niwsmax) niwsmax=niwsloc
     enddo
     gamsq=pi*ion_mass*1.6726e-24*(niwsmax)
     if (verbose .ge. 3) then
        write(*,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(outfile,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(*,*) 'for wsmodel=2 using full maximum value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        write(outfile,*) 'for wsmodel=2 using full maximum value ', &
             'to determine gamsq with ion mass ion_mass=',ion_mass
        write(*,*) 'assigned gamsq=',gamsq
        write(outfile,*) 'assigned gamsq=',gamsq
     endif
  endif

  if (wsmodel.eq.3 .or. wsmodel.eq.4 .or. wsmodel.eq.5) then
     !       use model 1 with modification for finite wavelength similar
     !          to Rogers-Drake 
     !          reduce omega* by 1/( 1+ 1/(k_theta L_p) )    
     !       use half of max value of omega_*_p_i**2 to determine gamsq
     !         spline onto the fine grid to find maximum omegas**2*ne
     niws=omegas**2*ne_eq(1) 
     call spline(-psigrid,niws,nxinterp,-1.d30,-1.d30,niwsspline)
     niwsmax=0.
     do i=1,nx
        call splint(-psigrid,niws,niwsspline,nxinterp,-psixx(i), &
             niwsloc)
        if (niwsloc > niwsmax) niwsmax=niwsloc
     enddo
     if (verbose .ge. 5) then
        write(*,*) 'maximum ni*omegas^2 value is ',niwsmax
        write(outfile,*) 'maximum ni*omegas^2 value is ',niwsmax
     endif
     ! find maximum in |p'|, use this to define approx ped center
     call spline(-psigrid,pprime_eq,nxinterp,-1.d30,-1.d30,ppspline)
     ppmax=0.
     do i=1,nx
        call splint(-psigrid,pprime_eq,ppspline,nxinterp,-psixx(i), &
             pploc)
        if (abs(pploc) > abs(ppmax)) then
           ppmax=abs(pploc)
           ippmax=i
        endif
     enddo
     if (ippmax>(nx-2)) then
        write(*,*) 'found max pprime at inner bound', &
             ' failed to locate pedestal center'
        stop
     endif
     ipphalf=0
     i=ippmax
     do while (ipphalf == 0)
        i=i+1
        call splint(-psigrid,pprime_eq,ppspline,nxinterp,-psixx(i), &
             pploc)
        if (abs(pploc) .le. abs(0.5*ppmax)) then
           ipphalf=i
           if (verbose .ge. 3) write(*,*) 'ipphalf=',ipphalf
        endif
     enddo

     ! spline to get radial locations off ppmax and pphalf
     call spline(-psigrid,rgrid,nxinterp,-1.d30,-1.d30,rspline)
     call splint(-psigrid,rgrid,rspline,nxinterp,-psixx(ippmax), &
          rppmax)
     call splint(-psigrid,rgrid,rspline,nxinterp,-psixx(ipphalf), &
          rpphalf)
     if (verbose .ge. 2) then
        write(*,*) 'radius of ped mid=',rppmax,psixx(ippmax)
        write(*,*) 'radius of ped half=',rpphalf,psixx(ipphalf)
        write(*,*) 'rgrid(1)=',rgrid(1)
        write(outfile,*) 'radius of ped mid=',rppmax,psixx(ippmax)
        write(outfile,*) 'radius of ped half=',rpphalf,psixx(ipphalf)
        write(outfile,*) 'rgrid(1)=',rgrid(1)
     endif
     !         delta=2.*(rppmax-rpphalf)
     delta=min(2.*(rppmax-rpphalf), &
          0.5*(rgrid(1)-(rppmax-2.*rpphalf)))
     if (verbose .ge. 3) then
        write(*,*) 'delta=',delta
        write(outfile,*) 'delta=',delta
     endif
     rminor=rmajor/aspect
     circum_approx=circumf*(rminor-(rgrid(1)-rppmax))/rminor
     if (verbose .ge. 3) then
        write(*,*) 'circum_approx=',circum_approx,circumf
        write(outfile,*) 'circum_approx=',circum_approx,circumf
     endif
     ! spline to get q at ppmax
     call spline(-psigrid,q_eq,nxinterp,-1.d30,-1.d30,qspline)
     call splint(-psigrid,q_eq,qspline,nxinterp,-psixx(ippmax), &
          qppmax)
     if (verbose .ge. 3) then
        write(*,*) 'qppmax=',qppmax
        write(outfile,*) 'qppmax=',qppmax
     endif
     ktheta=2.*pi*nn*qppmax/circum_approx
     if (verbose .ge. 3) then
        write(*,*) 'ktheta=',ktheta,'k d=',ktheta*delta, &
             '1/(1+1/k d)=',1./(1.+1./(ktheta*delta) ) 
        write(outfile,*) 'ktheta=',ktheta,'k d=',ktheta*delta, &
             '1/(1+1/k d)=',1./(1.+1./(ktheta*delta) )
     endif
     niwsmax=niwsmax* (1./(1.+1./(ktheta*delta)))**2.
     if (verbose .ge. 2) then
        write(*,*) 'modified niwsmax=',niwsmax
        write(outfile,*) 'modified niwsmax=',niwsmax
     endif
     if (wsmodel .eq. 3) then
        if (verbose .ge. 3) then
           write(*,*) 'for wsmodel=3 using half this modified value ', &
                'to determine gamsq with ion mass ion_mass=',ion_mass
           write(outfile,*) 'for wsmodel=3 using half this modified', &
                'value to determine gamsq with ion mass ion_mass=',ion_mass
        endif
        gamsq=pi*ion_mass*1.6726e-24*(0.5*niwsmax)
     endif
     if (wsmodel .eq. 4) then
        gamsq=pi*ion_mass*1.6726e-24*(0.125*niwsmax)
        if (verbose .ge. 3) then
           write(*,*) 'for wsmodel=4 using 1/8 this modified value ', &
                'to determine gamsq with ion mass ion_mass=',ion_mass
           write(outfile,*) 'for wsmodel=4 using 1/8 this modified', &
                'value to determine gamsq with ion mass ion_mass=',ion_mass
           write(*,*) 'assigned gamsq=',gamsq
           write(outfile,*) 'assigned gamsq=',gamsq
        endif
     endif
     if (wsmodel .eq. 5) then
        gamsq=pi*ion_mass*1.6726e-24*(0.0625*niwsmax)
        if (verbose .ge. 3) then
           write(*,*) 'for wsmodel=5 using 1/16 this modified value ', &
                'to determine gamsq with ion mass ion_mass=',ion_mass
           write(outfile,*) 'for wsmodel=5 using 1/16 this modified', &
                'value to determine gamsq with ion mass ion_mass=',ion_mass
           write(*,*) 'assigned gamsq=',gamsq
           write(outfile,*) 'assigned gamsq=',gamsq
        endif
     endif
     
  endif
!  if wsmodel is < 0 then set gamsq so that gamma/omega_A = 0.01 * abs(wsmodel)
  if (wsmodel .le. -1) then 
     gamsq=((0.01*abs(wsmodel))**2)/gamscl2
     gamre=sqrt(gamsq)
     if (verbose .ge. 3) then
        write(*,*) 'for wsmodel=',wsmodel,' setting gamma/omega_A=',0.01*abs(wsmodel)
        write(outfile,*) 'for wsmodel=',wsmodel,' setting gamma/omega_A=',0.01*abs(wsmodel)
        write(*,*) 'assigned gamsq=',gamsq
        write(outfile,*) 'assigned gamsq=',gamsq
     endif
  endif
  if (wsmodel .ne. 0) gamre=sqrt(gamsq)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! now do setup in earlier loop
  !         call matsetup          ! calculate the matrices needed to
  ! construct P, Q, and S at the
  ! interpolation x points

  ! allocate arrays needed for shoot,matgen,fcal
  call alloc_solveit

  if (verbose .ge. 2) write(outfile,*) 'xx(1)=',xx(1),'xx(nx)=',xx(nx)
  if (verbose .ge. 2) write(outfile,*) 'xinterp=',xinterp
  if (igam.le.1) then
     errnow=1.
     errgam=1.0e-6
     gamrtol=2.0e-5
     ! add most inertia terms here rather than 
     !  in matgen for efficiency 7/01
     if (rotation) gamsq=gamre**2
     auarr(8,:,:,:)=auarr(8,:,:,:)-gamsq*iuarr(3,:,:,:)
     auarr(6,:,:,:)=auarr(6,:,:,:)-gamsq*iuarr(2,:,:,:)
     auarr(7,:,:,:)=auarr(7,:,:,:)-gamsq*iuarr(1,:,:,:)
     auarr(16,:,:,:)=auarr(16,:,:,:)-gamsq*iuarr(6,:,:,:)
     auarr(14,:,:,:)=auarr(14,:,:,:)-gamsq*iuarr(5,:,:,:)
     auarr(15,:,:,:)=auarr(15,:,:,:)-gamsq*iuarr(4,:,:,:)

     if (splineeq .and. rotation) then
        write(*,*) 'splineeq not compatible with rotation'
        stop
     endif
     if (splineeq) then
        if (verbose .ge. 5) write(*,*) 'call splinearr'
        call splinearr
        if (verbose .ge. 5) write(*,*) 'return from splinearr'
     endif

     call shoot(gamr,gami,vcomp,gamr2,gami2)
  else if (rotation) then
!!     call findzero(gamr,gami,vr,vi,gamr2,gami2)
!!     call newton2d
     call compnewton(gamr,gami,vcomp,gamr2,gami2)
     gamsq=gamre**2
  else
     errnow=1.
     errgam=1.0e-6
     gamrtol=2.0e-5
     do 10 i=1,igam
        if ((errnow.lt.errgam).and. &
             (abs(yold*aspect*qref).lt.gamrtol)) goto 10
        ! add most inertia terms here rather than 
        !  in matgen for efficiency 7/01
        auarr(8,:,:,:)=auarr(8,:,:,:)-(gamsq-xold)*iuarr(3,:,:,:)
        auarr(6,:,:,:)=auarr(6,:,:,:)-(gamsq-xold)*iuarr(2,:,:,:)
        auarr(7,:,:,:)=auarr(7,:,:,:)-(gamsq-xold)*iuarr(1,:,:,:)
        auarr(16,:,:,:)=auarr(16,:,:,:)-(gamsq-xold)*iuarr(6,:,:,:)
        auarr(14,:,:,:)=auarr(14,:,:,:)-(gamsq-xold)*iuarr(5,:,:,:)
        auarr(15,:,:,:)=auarr(15,:,:,:)-(gamsq-xold)*iuarr(4,:,:,:)

        if (splineeq .and. rotation) then
           write(*,*) 'splineeq not compatible with rotation'
           stop
        endif
        if (splineeq) then
           if (verbose .ge. 5) write(*,*) 'call splinearr'
           call splinearr
           if (verbose .ge. 5) write(*,*) 'return from splinearr'
        endif

        call shoot(gamr,gami,vcomp,gamr2,gami2)
! use norm for rotation case with complex eigenvalue
        if (rotation) then
           gamnorm=gamr2
!           gamnorm=gamr
!           do j=1,nedge
!              gamnorm(j)=sqrt(gamr2(j)**2+gami2(j)**2)
!           enddo
        else
           gamnorm=gamr
        endif
!        write(*,*) 'gamr=',gamr,'gamr2=',gamr2
        do j=1,nedge
           if (j.eq.1) then
              gamax=gamnorm(j)
              igmax=j
           else
              if (gamax.gt.gamnorm(j)) then
                 gamax=gamnorm(j)
                 igmax=j
              end if
           end if
        end do
        if (i.eq.1) then
           errnow=1.
           yold=gamnorm(igmax)
           xold=gamsq
           
           if (.not.rotation) then
              if (verbose .ge. 1) write(6,12) &
                   i,gamsq,errnow,aspect*qref*gamnorm(igmax)
!              write(6,*)' '
              if (verbose .ge. 2) write(outfile,12) &
                   i,gamsq,errnow,aspect*qref*gamnorm(igmax)
           else
              if (verbose .ge. 1) write(6,13) &
                   i,gamsq,gamim,errnow,aspect*qref*gamr2(igmax), &
                   aspect*qref*gami2(igmax)
              if (verbose .ge. 2) write(outfile,13) i,gamsq,gamim,errnow, &
                   aspect*qref*gamr2(igmax),aspect*qref*gami2(igmax)
           endif

! 8.05, update to allow newtstep as input value
           if (yold.le.0.) gamsq=gamsq*(1.+newtstep)
           if (yold.gt.0.) gamsq=gamsq/(1.+newtstep)

!           if (gamsq.gt.0.5) then
!              gamsq=gamsq+1.
!           else
!              if (yold.le.0.) then
!                 gamsq=gamsq*2.
!              else
!                 gamsq=gamsq*0.5
!              endif
!           endif
!
!           ! change by more than one if gamsq starts large
!           if (gamsq.gt.5..and.yold.le.0.) gamsq=(gamsq-1.)*1.01
!           if (gamsq.gt.5..and.yold.gt.0.) gamsq=(gamsq-1.)*0.99
        else
           !               gamnew=(yold*gamsq-gamr(igmax)*xold)/(yold-gamr(igmax))
           gamnew=(yold*gamsq-gamnorm(igmax)*xold)/(yold-gamnorm(igmax)) 
           errnow=abs((gamnew-gamsq)/(gamnew+gamsq))
           !               errnow=abs((gamnew-gamsq)/(abs(gamnew)+abs(gamsq)))
           if (.not.rotation) then
              if (verbose .ge. 1) write(6,12)i,gamsq,errnow,aspect*qref*gamr(igmax)
              if (verbose .ge. 2) write(outfile,12)i,gamsq,errnow,aspect*qref*gamr(igmax)

! debug
!              write(*,*) 'yold=',yold,' igmax=',igmax,' gamr=',gamr(igmax), &
!                   ' xold=',xold,' gamnew=',gamnew,' errnow=',errnow

           else
              if (verbose .ge. 1) write(6,13) &
                   i,gamsq,gamim,errnow,aspect*qref*gamr2(igmax), &
                   aspect*qref*gami2(igmax)
              if (verbose .ge. 2) write(outfile,13)i,gamsq,gamim,errnow, &
                   aspect*qref*gamr2(igmax),aspect*qref*gami2(igmax)
           endif
           !don't let it go negative
           if (gamnew.lt.0..and.gamsq.gt.0.) gamnew=0.1*min(gamsq,xold)
           !don't let it increase too fast
           if ((gamnew/gamsq).gt.20. .and. (gamnew/xold).gt.20.) &
                gamnew=20.*max(gamsq,xold)

           ! pbs 12/10 if newtstep set by hand to small value, don't let it change
           !   by more than newtstep
           if (newtstep .lt. 0.2) then
             ! max fractional change - increase with step to avoid loops
              maxfrac=1.+newtstep*(i**0.3)
              if ( (gamnew/gamsq) .gt. maxfrac ) &
                   gamnew=gamsq*maxfrac
              if ( (gamnew/gamsq) .lt. (1./maxfrac) ) &
                   gamnew=gamsq/maxfrac
           endif

           xold=gamsq
           gamsq=gamnew
!           !don't let it go negative 
!           if (gamnew.lt.0..and.xold.gt.0.) gamsq=xold/10.
           yold=gamnorm(igmax)

        end if
      ! set gamre and gamim values based on iterated gamsq
        if (gamsq .ge. 0.) then
           gamre=sqrt(gamsq)
           gamim=0.
        else
           gamre=0.
           gamim=sqrt(-gamsq)
           write(*,*) '!! found negative gamsq for no rotation case', &
                'gamsq=', gamsq,' gamre=',gamre,' gamim=',gamim
           if (verbose .ge. 2) write(outfile,*) &
                '!! found negative gamsq for no rotation case', &
                'gamsq=', gamsq,' gamre=',gamre,' gamim=',gamim          
        endif
10      continue

12      format('Iteration number =',i3,' gamsq=',e14.6,   &
             ' Error=',e13.5,' gamr=',e14.6)
13      format('Iteration #=',i3,' gamsq=',e14.6, &
             ' gamim=',e14.6,' Error=',e13.5, &
             'gamr=',e14.6,'gami2=',e14.6)
     end if
     do i=1,nedge ! renormalize eigenvalues
        gamr(i)=aspect*qref*gamr(i)
        gami(i)=aspect*qref*gami(i)
        gamr2(i)=aspect*qref*gamr2(i)
        gami2(i)=aspect*qref*gami2(i)
     enddo
     !         call sortout(gamr,gami)

! TEMP!!!!
     vcdum=vcomp
     call sortout(gamr2,gami2,vcdum)
! only using the eigvalues associated with gamr2 for rotation case
     if (rotation) then
        vcomp=vcdum ! use sorted "2" eigenvectors
     endif
     !         if (nalpha.gt.1) then
     !            write(6,34)alpha,gamr(1),gami(1)
     ! 34         format('alpha=',e13.5,' gamma=',e13.5,', ',e13.5)
     !         else
     !     Write out eigenvalues
     if (verbose .ge. 2) then
        write(6,*)' '
        write(outfile,*) ' '
        write(6,*)' nn=',nn,' del=',del
        write(outfile,*)' nn=',nn,' del=',del
        write(6,*)'*************Calculated Growth Rates****************'
        write(outfile,*)'*************Calculated Growth Rates****************'
     endif
     if (((errnow.lt.errgam).and.(abs(gamr(1)).lt.gamrtol)).or.(rotation)) then
!  growth rate converged
!   check convergence of fictitious eigenvalues
        if (abs(gamr(1)).lt.1.e-4 .and. abs(gamr2(1)).lt.1.e-4) then
           if (verbose .ge. 1) write(6,*) 'fictitious eigenvalues well converged'
           if (verbose .ge. 2) write(outfile,*) &
                'fictitious eigenvalues well converged'
        else if (abs(gamr(1)).lt.1.e-3 .and. abs(gamr2(1)).lt.1.e-3) then
           if (verbose .ge. 1) write(6,*) &
                'fictitious eigenvalues marginally well converged'
           if (verbose .ge. 2) write(outfile,*) &
                'fictitious eigenvalues marginally well converged'
        else if (abs(gamr(1)).lt.1.e-2 .and. abs(gamr2(1)).lt.1.e-2) then
           if (verbose .ge. 1) write(6,*) &
                '!!fictitious eigenvalue convergence is marginal, growth rates may not be accurate!!'
           if (verbose .ge. 1) write(6,*) &
                ' !!Higher resolution and/or more tightly converged equilibrium should help!!'
           if (verbose .ge. 2) write(outfile,*) &
                '!!fictitious eigenvalue convergence is marginal, growth rates may not be accurate!!'
           if (verbose .ge. 2) write(outfile,*) &
                ' !!Higher resolution and/or more tightly converged equilibrium should help!!'
        else
           if (verbose .ge. 1) then
              write(6,*) '!!!!!WARNING!!!! At least one of the fictitious eigenvalues is not converged to 1e-2!!'
              write(6,*) ' !!!!Growth rates are likely not accurate!!!!'
              write(6,*) ' !!!!Higher resolution and/or more tightly converged equilibrium are needed!!!!!!'
           endif
           if (verbose .ge. 2) then
              write(outfile,*) '!!!!!WARNING!!!! At least one of the fictitious eigenvalues is not converged to 1e-2!!'
              write(outfile,*) ' !!!!Growth rates are likely not accurate!!!!'
              write(outfile,*) ' !!!!Higher resolution and/or more tightly converged equilibrium are needed!!!!!!'
           endif
        endif

!        write(6,*)'WARNING***Check fictitious eigenvalue (below) is close to zero'
        if (verbose .ge. 1) write(6,*) '***Growth rates with geometric center norm***'
        if (verbose .ge. 2) write(outfile,*) '***Growth rates with geometric center norm***'

        if (verbose .ge. 1) write(6,*) &
             'center norm, (gam/om_A)**2=',(gamre**2*gamscl2),' (R/B)**2=',gamscl2
        if (verbose .ge. 2) write(outfile,*) &
             'center norm, (gam/om_A)**2=',(gamre**2*gamscl2),' (R/B)**2=',gamscl2
        
        if (verbose .ge. 1) write(6,*) &
             'center norm, (gam/om_A)=',gamre*sqrt(gamscl2),' (R/B)=',sqrt(gamscl2)
        if (verbose .ge. 2) write(outfile,*) &
             'center norm, (gam/om_A)=',gamre*sqrt(gamscl2),' (R/B)=',sqrt(gamscl2)
        
!        write(6,*) '***Growth rates with old R/B in out average***'
        if (verbose .ge. 3) then
           write(outfile,*) '***Growth rates with old R/B in out average***'
!        write(6,*)'old norm, (gam/om_A)**2=',gamre**2*gamscl,' (R/B)**2=', &
!             gamscl
           write(outfile,*)'old norm, (gam/om_A)**2=',gamre**2*gamscl,' (R/B)**2=',gamscl
!        write(6,*)'old norm, (gam/om_A)=',gamre*sqrt(gamscl),' (R/B)=', &
!             sqrt(gamscl)
           write(outfile,*)'old norm (gam/om_A)=',gamre*sqrt(gamscl),' (R/B)=',sqrt(gamscl)
        endif
        dedge=1.
        if (dens) dedge=ne_eq(1) 
    !! gamsq is in CGS units for all file types
        if (verbose .ge. 1) write(6,*) &
             ' Growth rate for rho=1kgm**-3=',sqrt(gamsq/(4.*pi*1.0d-3)),' s^{-1}'
        if (verbose .ge. 2) write(outfile,*) &
             ' Growth rate for rho=1kgm**-3=',sqrt(gamsq/(4.*pi*1.0d-3)),' s^{-1}'
!!!!PBS 7.04 I'm confused by this- ne_eq is density, not 4pi*rho
!!!!     write(6,*)' Growth rate=',sqrt(gamsq/dedge),' s^{-1}'
!!!! growth rate in s^-1 output below when dens=true
        if (rotation) then
           if (verbose .ge. 1) then
              write(6,*) ' gamim=',gamim,' psi_inner=',psiv_eq(nxinterp)/psiv_eq(1)
              write(6,*)' Omega_sep=',rot_eq(1),' Omega_inner=',rot_eq(nxinterp)
              write(6,*)' Cs (outer,inner)=',sqrt(cssq1),sqrt(cssq2),' cm s^-1'
              write(6,*)' Va (outer,inner)=',rva1/rav,rva2/rav,' cm s^-1'
              write(6,*)' R*Omega_sep/Cs=',rav*rot_eq(1)/sqrt(cssq1),' R*Omega_inner/Cs=',rav*rot_eq(nxinterp)/sqrt(cssq2)
              write(6,*)' rav=',rav
           !        write(6,*) ' gamim=',gamim, ' Omega_sep=',rot_eq(1)

! pbs 6/06 omega_sep no longer subtracted out
              write(6,*) ' in lab frame omega/omega_A =',-gamim*sqrt(gamscl)
           endif
           if (verbose .ge. 2) then
              write(outfile,*) ' gamim=',gamim,' psi_inner=',psiv_eq(nxinterp)/psiv_eq(1)
              write(outfile,*)' Omega_sep=',rot_eq(1),' Omega_inner=',rot_eq(nxinterp)
              write(outfile,*)' Cs (outer,inner)=',sqrt(cssq1),sqrt(cssq2),' cm s^-1'
              write(outfile,*)' Va (outer,inner)=',rva1/rav,rva2/rav,' cm s^-1'
              write(outfile,*)' R*Omega_sep/Cs=',rav*rot_eq(1)/sqrt(cssq1),' R*Omega_inner/Cs=',rav*rot_eq(nxinterp)/sqrt(cssq2)
              write(outfile,*)' rav=',rav
              write(outfile,*) ' gamim=',gamim,' in lab frame omega/omega_A =',-gamim*sqrt(gamscl)
           endif
        endif
        if (dens) then
           rot_norm=sqrt(4.*pi*ion_mass*1.6726e-24*ne_eq(1))
           if (rotation) then
              if (verbose .ge. 1) write(6,*) ' in omega_sep frame omega/omega_A =',-gamim*sqrt(gamscl)-nn*rot_eq(1)*sqrt(gamscl)*rot_norm
              if (verbose .ge. 2) write(outfile,*) ' in omega_sep frame omega/omega_A =',-gamim*sqrt(gamscl)-nn*rot_eq(1)*sqrt(gamscl)*rot_norm
           endif
           if (verbose .ge. 1) write(6,*) 'ne_eq(1)=',ne_eq(1),'sqrt(4pi rho0)=',rot_norm
           if (verbose .ge. 2) write(outfile,*) 'ne_eq(1)=',ne_eq(1),'sqrt(4pi rho0)=',rot_norm
           if (verbose .ge. 1) write(6,*) 'growth rate= ',gamre/rot_norm, ' s^-1'
           if (verbose .ge. 2) write(outfile,*) 'growth rate= ',gamre/rot_norm, ' s^-1'
           if (verbose .ge. 1) write(6,'(a19,E13.6)') ' gam/omegaspi_max=',gamre/rot_norm/omegaspi_max
           if (verbose .ge. 2) write(outfile,'(a19,E13.6)') ' gam/omegaspi_max=',gamre/rot_norm/omegaspi_max
           if (rotation) then
              if (verbose .ge. 1) write(6,*) 'lab frame omega= ',-gamim/rot_norm,' rad/s'
              if (verbose .ge. 2) write(outfile,*) 'lab frame omega= ',-gamim/rot_norm,' rad/s'
              matchrot=-gamim/rot_norm/nn
              call spline(rot_eq,-psigrid,nxinterp,-1.d30,-1.d30,rpsispline)
              call splint(rot_eq,-psigrid,rpsispline,nxinterp, &
                   matchrot,psirot)
              if (verbose .ge. 1) then
                 write(6,*) ' Mode co-rotates with plasma where rot_eq=',matchrot/(2.*pi*1000.),' kHz'
                 write(6,*) '  This occurs near psi=',-psirot
              endif
              if (verbose .ge. 2) then
                 write(outfile,*) ' Mode co-rotates with plasma where rot_eq=', matchrot/(2.*pi*1000.),' kHz'
                 write(outfile,*) '  This occurs near psi=',-psirot
              endif
           endif
        endif
     else
        if (igam.eq.1) then
           if (verbose .ge. 1) then
              write(6,*) '(igam=1, stability check only, no growth rate calculation)'
              write(6,*) '(set igam to higher value to iterate to find growth rate)'
           endif
           if (verbose .ge. 2) then
              write(outfile,*) '(igam=1, stability check only, no growth rate calculation)'
              write(outfile,*) '(set igam to higher value to iterate to find growth rate)'
           endif
        else if ((igam>5).and.(gamsq<1.e-5).and.(gamr(1)>0.) &
             .and.(gamr2(1)>0.).and.(yold>0.)) then
           if (verbose .ge. 1) then
              write(6,*) ' Equilibrium appears to be stable'
              write(6,*) '  Did not converge to a finite growth rate'
           endif
           if (verbose .ge. 2) then
              write(outfile,*) ' Equilibrium appears to be stable'
              write(outfile,*) '  Did not converge to a finite growth rate'
           endif
        else
           if (verbose .ge. 1) then
              write(6,*)' Growth rate not sufficiently well converged'
              write(6,*)' Set igam to higher value (no. iterations) or improve guess of gamsq'
              write(6,*) '  If this fails, may need better converged or higher resolution equilibrium'
           endif
           if (verbose .ge. 2) then
              write(outfile,*)' Growth rate not sufficiently well converged'
              write(outfile,*)' Set igam to higher value (no. iterations) or improve guess of gamsq'
              write(outfile,*) '  If this fails, may need better converged or higher resolution equilibrium'
           endif
        endif
     end if
     
     if (verbose .ge. 2) then
        write(outfile,*)' gamsq=',gamsq
        write(outfile,*)'alpha=',alpha
        write(outfile,*)'shear=',shear
        write(outfile,*)'ppmult=',ppmult,' alphamult=',alphamult,' shearmult=',shearmult
     endif
     if (verbose .ge. 1) then
        write(6,*)' '
        write(6,*)'**************** fictitious eigenvalue****************'
        write(6,*)'************************* Sym eig, Non-sym eig'
        if (verbose .ge. 2) write(outfile,*)'************************* Sym eig, Non-sym eig'
        do i=1,1
           write(6,'(E18.10,E13.5,5X,E18.10,E13.5)') gamr(i),gami(i),gamr2(i),gami2(i)
        end do
        write(6,*)' gamsq=',gamsq
     endif
! 8.04 issue warning if there is large discrepancy between sym and asym eigvals
     !  for now require relative difference less than 20% or absolute < 0.1
     if (abs((gamr(1)-gamr2(1))/gamr(1)).gt.0.2) then
        if (abs(gamr(1)-gamr2(1)).gt.0.01) then
           write(6,*) '!!!!!Large discrepancy between symmetric and asymmetic eigenvalues!!!!'
           write(6,*) ' !!!Suggests that the equilibrium is not sufficiently resolved and/or converged!!!'   
           write(6,*) ' !!!Stability results are likely to be unreliable !!!!!!'
           if (verbose .ge. 2) then 
              write(outfile,*) '!!!!!Large discrepancy between symmetric and asymmetic eigenvalues!!!!'
              write(outfile,*) ' !!!Suggests that the equilibrium is not sufficiently resolved and/or converged!!!'   
              write(outfile,*) ' !!!Stability results are likely to be unreliable !!!!!!'
           endif
        endif
     endif

! 8.04 add feedback on results of stability check
! 5/14 compute unst_count
     unst_count=0
     if (igam.eq.1) then
        if (gamr(1).lt.0.) unst_count=1
        if (gamr2(1).lt.0.) unst_count=unst_count+1
        if (gamsq.le.0.1) then
           if (gamr(1).lt.0.) then
              if (gamr2(1).lt.0.) then
                 if (verbose .ge. 1) write(6,*) 'Stability check finds mode is **UNSTABLE**'
                 if (verbose .ge. 2) write(outfile,*) 'Stability check finds mode is **UNSTABLE**'
              else
                 if (verbose .ge. 1) write(6,*) 'Stability check is **INCONCLUSIVE**'
                 if (verbose .ge. 2) write(outfile,*) 'Stability check is **INCONCLUSIVE**'
                 if (verbose .ge. 1) write(6,*) ' Higher resolution and/or more tightly converged equilibrium is needed'
                 if (verbose .ge. 2) write(outfile,*) ' Higher resolution and/or more tightly converged equilibrium is needed'
                 if (gamsq.lt.1e-6) then
                    if (verbose .ge. 1) write(6,*) &
                         'Increasing gamsq to a small but finite value, eg 0.0001, may also be helpful'
                    if (verbose .ge. 2) write(outfile,*) &
                         'Increasing gamsq to a small but finite value, eg 0.0001, may also be helpful'
                 endif
              endif
           else   ! gamr(1) is positive
              if (gamr2(1).gt.0) then
                 if (verbose .ge. 1) write(6,*) 'Stability check finds mode is **STABLE**'
                 if (verbose .ge. 2) write(outfile,*) 'Stability check finds mode is **STABLE**'
              else
                 if (verbose .ge. 1) write(6,*) 'Stability check is **INCONCLUSIVE**'
                 if (verbose .ge. 2) write(outfile,*) 'Stability check is **INCONCLUSIVE**'
                 if (verbose .ge. 1) write(6,*) ' Higher resolution and/or more tightly converged equilibrium is needed'
                 if (verbose .ge. 2) write(outfile,*) ' Higher resolution and/or more tightly converged equilibrium is needed'
                 if (gamsq.lt.1e-6) then
                    if (verbose .ge. 1) write(6,*) 'Increasing gamsq to a small but finite value, eg 0.0001, may also be helpful'
                    if (verbose .ge. 2) write(outfile,*) 'Increasing gamsq to a small but finite value ,eg 0.0001, may also be helpful'
                 endif
              endif
           endif
        else
           if (verbose .ge. 2) then
              write(6,*)' Need to set gamsq to a smaller value, typically'
              write(6,*)'  around 0.0001 in order to do a stability check'
              write(outfile,*)' Need to set gamsq to a smaller value, typically'
              write(outfile,*)'  around 0.0001 in order to do a stability check'
           endif
        endif
     endif

     if (verbose .ge. 1) write(6,*) &
          '*******************************************************'
     if (verbose .ge. 1) write(6,*)' '
     do  i=1,nedge
        if (verbose .ge. 2) write(outfile,'(I2,1p,E18.10,E13.5,5X,E18.10,E13.5)') &
             i,gamr(i),gami(i),gamr2(i),gami2(i)
     end do

! open and write to .gamma file for convenient one line output
     nunit=49
     open(unit=nunit,file=runname(1:lrunname)//'.gamma', &
          status='unknown',iostat=ios)
     if(ios.ne.0) then
        write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.gamma'
        stop
     endif

     if (igam .gt. 1) then
        if (dens) then
           rot_norm=sqrt(4.*pi*ion_mass*1.6726e-24*ne_eq(1))
           write(49,*) 'nn cent_gam/w_A  gam/(omegspi_max/4)   gamsq     del      gamr        gamr2      tmatasym'
           write(49,2999) nn,sqrt(gamsq*gamscl2), &
                gamre/rot_norm/(omegaspi_max/4.),gamsq,del, &
                gamr(1),gamr2(1),tmatasym
        else
           write(49,*) 'nn cent_gam/w_A  cent_gam**2/w_A**2    gamsq     del      gamr        gamr2      tmatasym'
           write(49,3000) nn,sqrt(gamsq*gamscl2),gamsq*gamscl2,gamsq,del, &
                gamr(1),gamr2(1),tmatasym
        endif
     else
        if (dens) then
           rot_norm=sqrt(4.*pi*ion_mass*1.6726e-24*ne_eq(1))
           write(49,*) 'nn cent_gam/w_A  gam/(omegspi_max/4)   gamsq     del      gamr        gamr2      tmatasym unst_ct'
           write(49,2999) nn,sqrt(gamsq*gamscl2), &
                gamre/rot_norm/(omegaspi_max/4.),gamsq,del, &
                gamr(1),gamr2(1),tmatasym,unst_count
        else
           write(49,*) 'nn cent_gam/w_A  cent_gam**2/w_A**2    gamsq     del      gamr        gamr2      tmatasym unst_ct'
           write(49,3000) nn,sqrt(gamsq*gamscl2),gamsq*gamscl2,gamsq,del, &
                gamr(1),gamr2(1),tmatasym,unst_count
        endif
     end if
!3000     format(i4,2f13.10,e14.6,f9.6,2e12.4,1f10.6)
!2999     format(i4,f13.10,f13.8,e14.6,f9.6,2e12.4,1f10.6)

3000     format(i4,1x,f13.9,1x,f13.8,1x,e14.6,1x,f8.5,2e12.4,1x,1f9.5,1x,i3)
2999     format(i4,1x,f13.9,1x,f13.8,1x,e14.6,1x,f8.5,2e12.4,1x,1f9.5,1x,i3)

! if rotation write .gamim file for compact output for use in scripts etc
     if (rotation) then
        nunit=48
        open(unit=nunit,file=runname(1:lrunname)//'.gamim', &
             status='unknown',iostat=ios)
        if(ios.ne.0) then
           write(6,*) 'problem creating/opening ', &
                runname(1:lrunname)//'.gamim'
           stop
        endif

        write(nunit,*) rot_const
        write(nunit,3002) gamsq
        write(nunit,3002) gamim
        write(nunit,*) 'nn  rot_const  gam/w_A   omega/w_A     gamsq         gamim      del        gamr2'
        write(nunit,3001) nn,rot_const,gamre*sqrt(gamscl2),-gamim*sqrt(gamscl),gamsq,gamim,del,gamr2(1)
3001    format(i4,1x,1f9.3,1x,f11.8,1x,f11.8,1x,f12.6,1x,f12.6,1x,f9.6,e12.4)
3002    format(1f16.8)
        close(nunit)
     endif

     if(autorun) then
        gammin=gamr(1)
        gammin2=gamr2(1)
        do i=2,nedge
           if (gamr(i) < gammin) then
              gammin=gamr(i)
              gammin2=gamr2(i)
              !                 gammin=amin1(gammin,gamr(i))
              !                  gammin2=amin1(gammin2,gamr2(i))
           endif
        end do
        !               write(6,'(f20.14)') gammin
        write(6,*) 'below are int(1000*tmatasym) and then gamr_sym,gamr_asym,tmatasym,gamsq for script use'
        write(6,*) int(1000*tmatasym)
        write(6,'(4f17.10)') gammin,gammin2,tmatasym,gamsq
        if (verbose .ge. 2) then
           write(outfile,*) 'below are int(1000*tmatasym) and then gamr_sym,gamr_asym,tmatasym,gamsq for script use'
           write(outfile,*) int(1000*tmatasym)
           write(outfile,'(4f17.10)') gammin,gammin2,tmatasym,gamsq
        endif
        ! 3/6/00 create file runname.sadel containing s,alpha,del for use
        !   by scripts
        if (verbose .ge. 5) then
           nunit=29
           open(unit=nunit,file=runname(1:lrunname)//'.sadel', &
                status='unknown',iostat=ios)
           if(ios.ne.0) then
              write(6,*) 'problem creating/opening ', &
                   runname(1:lrunname)//'.sadel'
              stop
           endif

           write(nunit,'(3f11.5)') shear(1),alpha(1),del
           
           close(nunit)
        endif
!        stop
        return
     end if
     icount=iargc()
     if(icount  >=  2) then
        call getarg(2,ilabarg)
        lilabarg=index(ilabarg,' ')-1
        write(6,*) 'ilabarg : ',ilabarg,', lilabarg: ',lilabarg
        ilab=0
        do m=1,lilabarg
           ilab=ilab+(iachar(ilabarg(m:m))-48)*10**(lilabarg-m)           
        enddo
        write(6,*) 'Chosen eigenvalue : ',ilab
     else
!        write(6,*) ' Select an eigenvalue label I (a real one please!)'
!        write(6,*) 'Type 0 to end'
!        read(5,*)ilab
        if (funcal) then
           ilab=1
        else
           ilab=0
        end if
     endif
     if(ilab.lt.1.or.ilab.gt.nedge) return
     !     Load up eigenvector for generation of full eigenfunction
     gam=gamr(ilab)
     do  m=1,nedge
        aval(m)=vcomp(m,ilab)
        if (verbose .ge. 6) write(6,*) 'm=',m,'aval=',aval(m)
     end do
     call fcal
     call wrteigfunc
     !            call wrteigfunc2
     call fun2d
     !         end if
     !      end do
     return

   end subroutine solveit
   !
   !----------------------------------------------------------------------
   !
   !

   subroutine compnewton(gamr,gami,vcomp,gamr2,gami2)
! solve newton's method using complex numbers
     use elite_data, only: gamsq,gamim,iuarr,mwindow,igam, &
          auarr,gamre,newtstep,outfile,verbose
     implicit none
     complex z(igam+1),a,b,f(igam+1)
!     real vr(mwindow,mwindow),vi(mwindow,mwindow)
     complex vcomp(mwindow,mwindow)
     real gamr(mwindow),gami(mwindow),gamnorm(mwindow)
     real gamr2(mwindow),gami2(mwindow),multfact
     real gamsqold,gamsqdiff
     integer i,j,igmax
     real gamax,errtol

!     multfact=0.98
!     multfact=0.998
     multfact=sqrt(1.+newtstep) ! 1+newtstep is for gamsq, so take root
     errtol=5.d-10
!     z(1)=cmplx (sqrt(gamsq),gamim)
     z(1)=cmplx(gamre,gamim)
     gamsq=gamre**2
     if (verbose .ge. 1) write(*,*) 'iteration 1   z(1)=',z(1)

     auarr(8,:,:,:)=auarr(8,:,:,:)-gamsq*iuarr(3,:,:,:)
     auarr(6,:,:,:)=auarr(6,:,:,:)-gamsq*iuarr(2,:,:,:)
     auarr(7,:,:,:)=auarr(7,:,:,:)-gamsq*iuarr(1,:,:,:)
     auarr(16,:,:,:)=auarr(16,:,:,:)-gamsq*iuarr(6,:,:,:)
     auarr(14,:,:,:)=auarr(14,:,:,:)-gamsq*iuarr(5,:,:,:)
     auarr(15,:,:,:)=auarr(15,:,:,:)-gamsq*iuarr(4,:,:,:)

!     write(*,*) 'call shoot in compnewton, gamsq=',gamsq,'gamim=',gamim
     if (verbose .ge. 1) write(*,*) '  gamsq=',gamsq,'gamim=',gamim
     call shoot(gamr,gami,vcomp,gamr2,gami2)
!     write(*,*) 'return from shoot'
! use norm for rotation case with complex eigenvalue
     gamnorm=gamr2
!     gamnorm=gamr
     do j=1,mwindow
        if (j.eq.1) then
           gamax=gamnorm(j)
           igmax=j
        else
           if (gamax.gt.gamnorm(j)) then
              gamax=gamnorm(j)
              igmax=j
           end if
        end if
     end do
     f(1)=cmplx(gamr2(igmax),gami2(igmax))
     if (verbose .ge. 1) write(*,*) '  f(1)=',f(1)

     gamsqold=gamsq
!     gamsq=gamsq*(multfact)**2
     gamim=gamim*multfact
     gamre=gamre*multfact
     gamsq=gamre**2
!     z(2)=cmplx(sqrt(gamsq),gamim)
     z(2)=cmplx(gamre,gamim)
     if (verbose .ge. 1) write(*,*) 'iteration 2   z(2)=',z(2)

     do i=2,igam
        gamsqdiff=gamsq-gamsqold
        auarr(8,:,:,:)=auarr(8,:,:,:)-gamsqdiff*iuarr(3,:,:,:)
        auarr(6,:,:,:)=auarr(6,:,:,:)-gamsqdiff*iuarr(2,:,:,:)
        auarr(7,:,:,:)=auarr(7,:,:,:)-gamsqdiff*iuarr(1,:,:,:)
        auarr(16,:,:,:)=auarr(16,:,:,:)-gamsqdiff*iuarr(6,:,:,:)
        auarr(14,:,:,:)=auarr(14,:,:,:)-gamsqdiff*iuarr(5,:,:,:)
        auarr(15,:,:,:)=auarr(15,:,:,:)-gamsqdiff*iuarr(4,:,:,:)

!        write(*,*) 'call shoot in compnewton, gamsq=',gamsq,'gamim=',gamim
        if (verbose .ge. 1) write(*,*) '  gamsq=',gamsq,'gamim=',gamim
        call shoot(gamr,gami,vcomp,gamr2,gami2)
!        write(*,*) 'return from shoot'
! use norm for rotation case with complex eigenvalue
        gamnorm=gamr2
!        gamnorm=gamr
        do j=1,mwindow
           if (j.eq.1) then
              gamax=gamnorm(j)
              igmax=j
           else
              if (gamax.gt.gamnorm(j)) then
                 gamax=gamnorm(j)
                 igmax=j
              end if
           end if
        end do
        f(i)=cmplx(gamr2(igmax),gami2(igmax))
        if (verbose .ge. 1) write(*,*) '  f(',i,')=',f(i)
        if (abs(f(i)) < errtol) then
           if (verbose .ge. 1) write(*,*) 'compnewton converged to ',errtol
           return
        endif
        a=(f(i-1)-f(i))/(z(i-1)-z(i))
        if (verbose .ge. 2) write(outfile,*) 'a=',a
        b=(f(i-1)*z(i)-f(i)*z(i-1))/(z(i)-z(i-1))
        if (verbose .ge. 2) write(outfile,*) 'b=',b
        z(i+1)=-b/a 
        if (i<igam) then
           if (verbose .ge. 1) write(*,*) 'iteration',i+1,'  z(',i+1,')=',z(i+1) 
           gamsqold=gamsq
!           gamsq=real(z(i+1))**2
           gamim=imag(z(i+1))
           gamre=real(z(i+1))
           gamsq=gamre**2
        endif
       
     enddo

   end subroutine compnewton



   subroutine sortout(arr,arrim,arr2)
     !
     !
     !-----------------------------------------------------------------------
     ! this routine orders the elements in an array giving lowest first
     !-----------------------------------------------------------------------
     use elite_data, only: mwindow,nedge
     implicit none
     integer i,j,ij
     real armin
     real arr(mwindow),arr_temp(mwindow),arrim(mwindow),arr_tempim(mwindow)
     complex arr2(mwindow,mwindow),arr2_temp(mwindow,mwindow)
     !
     armin=arr(1)
     arr_temp(1)=armin
     ij=1
     do i=2,nedge
        if (arr(i).lt.armin) then
           arr_temp(1)=arr(i)
           armin=arr(i)
           arr_tempim(1)=arrim(i)
           arr2_temp(:,1)=arr2(:,i)
           ij=i
        end if
     end do
     arr(ij)=arr_temp(1)-100.
     do j=2,nedge
        ij=0
        do i=1,nedge
           if (arr(i).ge.arr_temp(1)) then
              if (ij.eq.0) then
                 ij=i
                 armin=arr(i)
              else
                 if (armin.gt.arr(i)) then
                    armin=arr(i)
                    ij=i
                 end if
              end if
           end if
        end do
        arr_temp(j)=armin
        arr_tempim(j)=arrim(ij)
        arr2_temp(:,j)=arr2(:,ij)
        arr(ij)=arr_temp(1)-100.
     end do
     do i=1,nedge
        arr(i)=arr_temp(i)
        arrim(i)=arr_tempim(i)
        arr2(:,i)=arr2_temp(:,i)
     end do
     return
   end subroutine sortout


   subroutine alloc_solveit
     ! allocate arrays used in solveit,shoot,matgen,fcal
     !   all those not needed until solveit is called
     use elite_data, only: pmat,qmat,smat,aa,ainv,bb,fmat,ffun, &
          aval,mwindow,mindex,kx,nmodes,comp_ai,comp_api, &
          comp_appi,compression
!     integer hmwindow ! half size to conserve memory for mindex=2 case

!     hmwindow=mwindow/mindex

     allocate( pmat(mwindow,mwindow), qmat(mwindow,mwindow), &
          smat(mwindow,mwindow) )
     pmat=0.; qmat=0.; smat=0.

     allocate( aa(mwindow,mwindow,kx), ainv(mwindow,mwindow,kx) )
     aa=0.; ainv=0.

     allocate( bb(mwindow+mindex,kx), fmat(3,mwindow,mwindow) )
     bb=0.; fmat=0.

     allocate( ffun(nmodes,kx),aval(nmodes) )
     ffun=0.; aval=0.

     if (compression) then
        allocate( comp_ai(mwindow,mwindow), comp_api(mwindow,mwindow), &
             comp_appi(mwindow,mwindow) )
        comp_ai=0.; comp_api=0.; comp_appi=0.
     endif

   end subroutine alloc_solveit
