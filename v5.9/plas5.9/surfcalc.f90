
subroutine surfcalc(isurf)

!-----------------------------------------------------------------------
!     calculates surface quantities on surface isurf
!      these are passed to other
!     routines through shared vars in elite_data
!     rl,zl,bpl,uprimel,cosul,sinul
!     circumf,psisep
!-----------------------------------------------------------------------
! 2/8/00 for now, f,ffprime,pprime are read from eqdat, not calculated here
      use elite_data, only: ns,outfile,rl,zl,bpl,uprimel,cosul,sinul, &
           runname,lrunname,nxinterp,f,f_eq,q_eq,rmajor,aspect,rnorm, &
           btnorm,circumf,vprime,pi,qref,vprime_eq,volume,dpsi1l, &
           d2psi1l,ppmult,pprime_eq,q_calc,ffprime_eq,sa_fix,shearmult, &
           shear_fix,alphamult,alpha_fix,alpha_edge,alpha,shear, &
           pprime_adj,ffprime_adj,qrefp,vpp,psisep,abishop,sigbishop, &
           shbishop,jedgen1store,nsurfdat,psigrid,jedgen1_arr, &
           jedgen2_arr,verbose
      implicit none
      integer nplot
      real wav,dav,lambishop
      real dtheta,t1,t2,t3
      real denom
      real x1,x2,x3
      real drdl,dr2dl2,dzdl,dz2dl2
      real, dimension(:), allocatable :: rcinvps,rpts,zpts,bppts,work, &
           uprimepts,cosupts,sinupts,r,theta,dummy,leng
!-----------------------------------------------------------------
!  10/20/00 to evaluate nupp in nustuff, we need additional quantities on
!  a flux surface: poloidal differentials of psi1 (see Miller PoP 5 (1998) 973)
      real, dimension(:), allocatable :: dpsi1dl,d2psi1dl2
      real qsurf,rmax,rmin,rminor
      real rnewton
      real leq(ns)
      real bp0bishop,r0bishop
      real qv
      real dr,drt
      integer i,ip,im,j,ir
      integer npts,nunit,ios
      real i1,i2,i3
      real dl
      real rhomax
      real area,intbp,vp,bsqinv,bpoldsq,bpnewsq
      real jedge_n1,jedge_n1pp,jedge_n1ffp
      real jedge_n2,jedge_n2pp,jedge_n2ffp
      real jedge_n3,jedge_n3pp,jedge_n3ffp
      real bsq,bpsqinv,bsqbpsqinv,bpave
      real vpp_c,vpp_pp,vpp_ffp
      real ffprime_c,ffprime_pp,ffprime_s,qp_s
      real pprime_loc,ffprime_loc,alpha_loc,sh_loc,ffprime_hirsh
      integer isurf,nxinterp_loc
      real bnd_set,dum2(ns),check,dpsi1dl_s(ns),d2psi1dl2_s(ns),dum
      real alpha_max,pprime_max,jedgen1_max
      real psi_alpha_max,psi_pprime_max,psi_jn1_max
      real jedgen2_max,sh_min,sh_al_max,sh_j1_max,sh_j2_max
      real psi_sh_min,psi_jn2_max
      real psifine,varfine,j2spline(nxinterp),ppspline(nxinterp)
      real aspline(nxinterp),sspline(nxinterp),j1spline(nxinterp)

      save npts
      if (verbose .ge. 2) write(outfile,*) &
           "==================begin output from surfcalc surface=",isurf
!-----------------------------------------------------------------------
!     1.0 read in surface info from *.surf file, read up to desired surface
!-----------------------------------------------------------------------

! zero out arrays to be careful
      rl=0.; zl=0.; bpl=0.; uprimel=0.; cosul=0.; sinul=0.

! 11/30/00 modify so file is only opened once and read continuously
      nunit=19
      if (isurf==1) then
         open(unit=nunit,file=runname(1:lrunname)//'.surf',iostat=ios)
         if(ios.ne.0) then
            write(6,*) 'could not find file ',runname(1:lrunname)//'.surf'
            stop
         endif
         read(nunit,*) npts,nxinterp_loc
         if(nxinterp_loc.ne.nxinterp) then
            write(*,*) 'in surfcalc nxinterp_loc must equal nxinterp'
            stop
         endif

!  can now allocate the surface arrays
!!  first time through only
!      if (isurf==1) then
      endif

      if (verbose .ge. 5) write(outfile,*) 'allocating surface arrays, npts= ',npts
      allocate( rcinvps(npts),rpts(npts),zpts(npts),bppts(npts), &
           work(npts),uprimepts(npts),cosupts(npts),sinupts(npts), &
           r(npts),theta(npts),dummy(npts),leng(npts), &
           dpsi1dl(npts),d2psi1dl2(npts) )
!      endif



      rcinvps=0.; rpts=0.; zpts=0.; bppts=0.; work=0.; uprimepts=0. 
      cosupts=0.; sinupts=0.; r=0.; theta=0.; dummy=0.; leng=0.

!      do i=1,isurf   ! read up to desired surface
         do j=1,npts  ! read in backwards as before
            ir=npts+1-j
            read(nunit,*) rpts(ir),zpts(ir),bppts(ir)
         enddo
!      enddo
!      close(nunit)      

      f=f_eq(isurf)
      qsurf=q_eq(isurf)

      leng(1)=0.
      rmax=rpts(1)
      rmin=rpts(1)
      do i=1,npts
         if (i.gt.1) then
            dr=rpts(i)-rpts(i-1)
            drt=zpts(i)-zpts(i-1)
            leng(i)=leng(i-1)+sqrt(dr**2+drt**2)
         end if
         rmax=amax1(rmax,rpts(i))
         rmin=amin1(rmin,rpts(i))
      end do
      rmajor=(rmax+rmin)*0.5
      rminor=(rmax-rmin)*0.5
      aspect=rmajor/rminor
      if (verbose .ge. 3) write(outfile,*)&
           ' Major radius=',rmajor,' Minor radius=',rminor, &
           ' surface=',isurf
      do i=1,npts
         theta(i)=atan2(zpts(i),(rpts(i)-rmajor))
      end do
      rnorm=rminor
      btnorm=f/rmajor
      circumf=leng(npts)
      if (verbose .ge. 2) write(outfile,*) 'btnorm=',btnorm,' circumf=',circumf

!  1/8/02 create name.surfdat file to store equilibrium info for
!     plotting.  File written to by both surfcalc and mercy
      if (isurf == 1 .and. (verbose .ge. 1)) then
         nsurfdat=62
         open(unit=nsurfdat,file=runname(1:lrunname)//'.surfdat', &
              status="unknown",iostat=ios)
         if(ios.ne.0) then
            write(6,*) 'problem creating/opening ', &
                 runname(1:lrunname)//'.surfdat'
            stop
         endif
         write(nsurfdat,*) nxinterp,rmajor,rminor,btnorm
      endif

!-----------------------------------------------------------------------
!     2.0 can now get correct normalization for bppts above
!-----------------------------------------------------------------------
!     we take bt_norm to be 1--the toroidal field at R=aspect so
!-----------------------------------------------------------------------
      qv=0.
      vprime=0.
      do i=2,npts
         dl=(leng(i)-leng(i-1))
         qv=qv+0.5*dl*(1./(rpts(i)**2*bppts(i))+ &
             1./(rpts(i-1)**2*bppts(i-1)))
         vprime=vprime+0.5*dl*(1./bppts(i)+ &
             1./bppts(i-1))
      end do
      qv=qv*f/(2.*pi)
      if (verbose .ge. 2) write(outfile,*) &
           'qv= ',qv,' q_eq=',q_eq(isurf),' qref=',qref
      if( (isurf==1) .and. ( abs(qv-qref)>.01 ) ) then
         write(*,*) 'for ref surface q_calc should = qa=qref',qv,qref
         stop
      endif
      q_calc(isurf)=qv
!      bp0bishop=qv/qref
      bp0bishop=1.
!-----------------------------------------------------------------------
!     note bishop has normalization bp0 in formula for bp(theta) and
!     it is just the above quantity--his r0 we take to be 1 so:
!-----------------------------------------------------------------------
      r0bishop=rnorm
      vprime=2.*pi*vprime/bp0bishop
      if (verbose .ge. 4) then
         write(outfile,*)'vprime=',vprime,' vprime_eq=',vprime_eq(isurf)
         write(outfile,'(1p,3(a,e12.4))') "bp0bishop=",bp0bishop
      endif

!-----------------------------------------------------------------------
!     3.0 calculate volume enclosed by a flux surface
!-----------------------------------------------------------------------
      volume=0.
      dtheta=2.*pi/real(npts-1.)
      do i=2,npts
         volume=volume+0.25*(zpts(i)+zpts(i-1))* &
             (rpts(i-1)**2-rpts(i)**2)
      enddo
      volume=-2.*pi*volume
!-----------------------------------------------------------------------
!     4.0 need combination rcinv+sinu/r for ffprime calculation
!-----------------------------------------------------------------------
      do i=1,npts
         if(i.eq.1) then
            x1=leng(npts-1)-leng(npts)
            x2=leng(i)
            x3=leng(i+1)
            ip=2
            im=npts-1
         else if(i.eq.npts) then
            x1=leng(i-1)
            x2=leng(i)
            x3=leng(i)+leng(2)
            ip=2
            im=npts-1
         else
            x1=leng(i-1)
            x2=leng(i)
            x3=leng(i+1)
            ip=i+1
            im=i-1
         endif
         denom=(x1-x2)*(x1-x3)*(x2-x3)
         drdl=(-rpts(ip)*(x1-x2)**2+rpts(im)*(x2-x3)**2+ &
             rpts(i)*(x1-x3)*(x1-2.*x2+x3))/denom
         dr2dl2=2.*(rpts(ip)*(x1-x2)+rpts(im)*(x2-x3)+rpts(i)*(x3-x1))/ &
             denom
         dzdl=(-zpts(ip)*(x1-x2)**2+zpts(im)*(x2-x3)**2+ &
             zpts(i)*(x1-x3)*(x1-2.*x2+x3))/denom
         dz2dl2=2.*(zpts(ip)*(x1-x2)+zpts(im)*(x2-x3)+zpts(i)*(x3-x1))/ &
             denom
         dpsi1dl(i)=(-rpts(ip)*bppts(ip)*(x1-x2)**2+rpts(im)*bppts(im) &
           *(x2-x3)**2+rpts(i)*bppts(i)*(x1-x3)*(x1-2.*x2+x3))/denom
         d2psi1dl2(i)=2.*(rpts(ip)*bppts(ip)*(x1-x2)+ &
           rpts(im)*bppts(im)*(x2-x3)+rpts(i)*bppts(i)*(x3-x1))/denom
!-----------------------------------------------------------------------
!     the following quantity is -u'+sinu/rs
!-----------------------------------------------------------------------
         rcinvps(i)=(drdl*dz2dl2-dzdl*dr2dl2)-dzdl/rpts(i)
         uprimepts(i)=-drdl*dz2dl2+dzdl*dr2dl2
         cosupts(i)=drdl
         sinupts(i)=-dzdl
      enddo
!-----------------------------------------------------------------------
!     6.0 spline r,z, and rbp onto an equal arc length mesh
!-----------------------------------------------------------------------
! perhaps stop doing this?  2/00  why is it done twice?
      dl=circumf/(ns-1.)
      do i=1,ns
         leq(i)=(i-1)*dl
      end do
      leq(ns)=circumf
      call spline1d(rl     ,leq,ns,rpts     ,leng,npts,work)
      call spline1d(zl     ,leq,ns,zpts     ,leng,npts,work)
      call spline1d(bpl    ,leq,ns,bppts    ,leng,npts,work)
      call spline1d(uprimel,leq,ns,uprimepts,leng,npts,work)
      call spline1d(cosul  ,leq,ns,cosupts  ,leng,npts,work)
      call spline1d(sinul  ,leq,ns,sinupts  ,leng,npts,work)
      call spline1d(dpsi1l  ,leq,ns,dpsi1dl  ,leng,npts,work)
      call spline1d(d2psi1l  ,leq,ns,d2psi1dl2  ,leng,npts,work)


!! 11/7/00 check above derivatives by calculating spline derivatives
!      bnd_set=-2.d30  ! calculate boundary condition from nearest 4 pts
!      call spline(leq,rl*bpl,ns,bnd_set,bnd_set,dum2)
!        ! dum2 contains the interpolating function
!      do i=1,ns
!         call zsplint(leq,rl*bpl,dum2,dum2,ns,leq(i),check,dpsi1dl_s(i),dum)
!!         write(*,*) 'i=',i,' qsfin=',qsfin(i),' check=',check, &
!!              'qsfinp=',qsfinp(i)
!!         write(74,*) 'i=',i,' dpsi1dl_s(i)=',dpsi1dl_s(i), &
!!              'dpsi1l=',dpsi1l(i),' uprimel=',uprimel(i), &
!!              'sinul=',sinul(i),'cosul=',cosul(i), &
!!              ' u=',acos(cosul(i)),' u=',asin(sinul(i))
!!         if (i.gt.1) write(74,*) 'dudl=', &
!!              (acos(cosul(i))-acos(cosul(i-1)))/dl, &
!!              (asin(sinul(i))-asin(sinul(i-1)))/dl
!!         write(74,*) 'uprimel=',uprimel(i)
!      enddo
!      call spline(leq,dpsi1dl_s,ns,bnd_set,bnd_set,dum2)
!        ! dum2 contains the interpolating function
!      do i=1,ns
!         call zsplint(leq,dpsi1dl_s,dum2,dum2,ns,leq(i),check, &
!              d2psi1dl2_s(i),dum)
!!         write(*,*) 'i=',i,' qsfin=',qsfin(i),' check=',check, &
!!              'qsfinp=',qsfinp(i)
!!         write(74,*) 'i=',i,' d2psi1dl2_s(i)=',d2psi1dl2_s(i), &
!!              'd2psi1l=',d2psi1l(i)
!! TEMPORARY CHANGE!!!!! use spline derivatives instead of 3 pt
!!         d2psi1l(i)=d2psi1dl2_s(i)
!!         dpsi1l(i)=dpsi1dl_s(i)
!      enddo

!-----------------------------------------------------------------------
! print quantities out for debug
!-----------------------------------------------------------------------
!      open(44,file='debug.bspline',status='unknown')
!      write(44,*) (leng(i),rpts(i),zpts(i),bppts(i), &
!          uprimepts(i),cosupts(i),sinupts(i),i=1,npts)
!      close(44)
!      open(45,file='debug.aspline',status='unknown')
!      write(45,*) (leq(i),rl(i),zl(i),bpl(i), &
!          uprimel(i),cosul(i),sinul(i),i=1,ns)
!      close(45)

! the rest should go 2/00 - maybe not all - grab ffprime & pprime 
!     from eqdat (ffprime_eq(isurf) & pprime_eq(isurf)
!-----------------------------------------------------------------------
!     5.0 calculate  pprime, and ffprime
!-----------------------------------------------------------------------
!      pprime=-alpha*(2.*pi)**2/(2.*vprime)* &
!          sqrt(2.*pi**2*rmajor/volume)
!      write(*,*) 'alpha= ',alpha,' pi=',pi,' vprime=',vprime, &
!           'volume=',volume,' rmajor=',rmajor
!  include adjustments due to ppmult
      if (ppmult.ne.1.) then
         if (verbose .ge. 2) write(outfile,*) &
              'Adjusting pprime by factor of ',ppmult
      endif
      pprime_loc=ppmult*pprime_eq(isurf)
!--------------------------------------------------------------------
! 3/2/00  first go through and calculate alpha and sh which correspond
!     to the equilibtium (*ppmult) pprime and ffprime.  Then incorporate
!      adjustment factors alphamult & shearmult, and calculate adjusted
!      pprime_adj & ffprime_adj that go along with these.  Or if sa_fix
!      is true, use fixed values shear_fix and alpha_fix
!---------------------------------------------------------------------
!      pprime_adj(isurf)=pprime_loc  ! assign adjusted p' value
!      pprime=pprime_loc  ! used in mercy
      if (verbose .ge. 2) write(outfile,*) 'pprime_eq*ppmult=',pprime_loc
      alpha_loc=-pprime_loc/((2.*pi)**2/(2.*vprime)* &
          sqrt(2.*pi**2*rmajor/volume) )
!      alpha(isurf)=alpha_loc 
      if (verbose .ge. 2) write(outfile,*) &
           'alpha corresponding to pprime_eq*ppmult=',alpha_loc
!      if (isurf==1) then
!         alpha_edge=alpha_loc
!         write(*,*) 'alpha_edge= ',alpha_loc
!      endif
      i1=0.
      i2=0.
      i3=0.
      wav=0.
      dav=0.
!-----------------------------------------------------------------------
!     in terms of the following integrals
!     q'=f'*q/f+f/(2*pi)*(-2*i2+p'*i3+ff'*i1)
!     since we know s,alpha,and f--we can solve for ff'
!-----------------------------------------------------------------------
      do i=2,npts
         dl=0.5*(leng(i)-leng(i-1))
         i1=i1+dl*(1./(rpts(i)**4*bppts(i)**3)+ &
             1./(rpts(i-1)**4*bppts(i-1)**3))
         i2=i2+dl*(1./(rpts(i)**3*bppts(i)**2)*rcinvps(i) &
             +1./(rpts(i-1)**3*bppts(i-1)**2)*rcinvps(i-1))
         i3=i3+dl*(1./(rpts(i)**2*bppts(i)**3)+ &
             1./(rpts(i-1)**2*bppts(i-1)**3))
         wav=wav+dl*bp0bishop*(r(i)*cos(theta(i))/bppts(i)+ &
             r(i-1)*cos(theta(i-1))/bppts(i-1))
         dav=dav+dl*bp0bishop*(1./bppts(i)+1./bppts(i-1))
      enddo
      wav=wav/(2.*pi)
      dav=dav/(2.*pi)
      denom=qv/f**2+f*i1/(2.*pi)
      if (verbose .ge. 4 .or. (verbose .ge. 2 .and. isurf==1)) then
         write(outfile,*)' denom=',denom,' f=',f,' pprime_loc=',pprime_loc
         write(outfile,*)' i2=',i2,' i3=',i3,' volume=',volume,' qv=',qv
         write(outfile,*)' vprime=',vprime,' rmajor=',rmajor
      endif
! no hirsh option for now
!      if (sh.lt.0.) then
!!     calculate ffprime and sh from b/strap+Ohmic+P-S currents
!         call hirsh
!!     write(6,*)' ffprime=',ffprime
!         t1=ffprime*denom
!         t2=f*i2/pi
!         t3=f/(2.*pi)*pprime_loc*i3
!         sh=2*volume*(ffprime*denom-f*i2/pi+f/(2.*pi)*pprime_loc*i3) &
!             /(vprime*qv)
!         write(6,*)' sh after hirsh=',sh
!      else
!     calculate ffprime from sh

! 2/00 just use the EFIT ffprime for this surface
     ffprime_loc=ffprime_eq(isurf)
!     ffprime_adj(isurf)=ffprime_loc ! store adjusted value 
!     ffprime=ffprime_loc  ! for use in mercy
     if (verbose .ge. 2) write(outfile,*) 'equilibrium ffprime=',ffprime_loc
!         ffprime=(f*i2/pi+vprime*qv/(2*volume)*sh- &
!             f/(2.*pi)*pprime_loc*i3)/denom
!      end if

     sh_loc=2.*volume*(ffprime_loc*denom-f*i2/pi+f/(2.*pi)*pprime_loc*i3) &
          /(vprime*qv)
!     shear(isurf)=sh_loc
     if (verbose .ge. 2) write(outfile,*)'unadjusted local shear=',sh_loc
     if (sa_fix) then
        write(*,*) 'fixing outer s&alpha using shear_fix & alpha_fix'
        if (isurf==1) then
           shearmult=shear_fix/sh_loc
           alphamult=alpha_fix/alpha_loc
        endif
     endif
     if ((alphamult.ne.1.).or.(shearmult.ne.1.)) then
        alpha_loc=alpha_loc*alphamult
        pprime_loc=pprime_loc*alphamult
        write(outfile,*) 'adjusted alpha=',alpha_loc,' pprime=',pprime_loc
        sh_loc=sh_loc*shearmult
        ffprime_loc=(f*i2/pi+vprime*qv/(2*volume)*sh_loc- &
             f/(2.*pi)*pprime_loc*i3)/denom
! 3/2/00 if shearmult=-1. determine shear and ffprime from hirsh
        if (shearmult==-1.) then
!!     calculate ffprime and sh from b/strap+Ohmic+P-S currents
          call hirsh(isurf,pprime_loc,ffprime_hirsh)
          ffprime_loc=ffprime_hirsh
          sh_loc=2.*volume*(ffprime_hirsh*denom-f*i2/pi+ &
               f/(2.*pi)*pprime_loc*i3)/(vprime*qv)
          if (verbose .ge. 2) write(outfile,*)' shear after hirsh=',sh_loc
       endif
       if (verbose .ge. 2) write(outfile,*) &
            'adjusted shear=',sh_loc,' ffprime=',ffprime_loc
     endif

      if (isurf==1) then
         alpha_edge=alpha_loc
         if (verbose .ge. 2) write(outfile,*) 'alpha_edge= ',alpha_loc
      endif
! can now assign adjusted values for use in mercy & nustuff
     alpha(isurf)=alpha_loc
     shear(isurf)=sh_loc
     pprime_adj(isurf)=pprime_loc
     ffprime_adj(isurf)=ffprime_loc
      ffprime_c=f*i2/pi/denom
      ffprime_s=vprime*qv/(2*volume)/denom
      ffprime_pp=-f/(2.*pi)*i3/denom
      qrefp=vprime*qv/(2*volume)*sh_loc
      qp_s=vprime*qv/(2*volume)
      if (verbose .ge. 4) write(outfile,'(1p,4(a,e12.4))') &
          'pprime_loc=',pprime_loc, &
          '  ffprime_loc=',ffprime_loc, &
          ' check=',ffprime_c+ffprime_pp*pprime_loc+ffprime_s*sh_loc, &
          '  f=',f
      if (verbose .ge. 2) write(outfile,'(1p,3(a,e12.4))')'qv=',qv,' qrefp=',qrefp, &
          "  circumf=",circumf
!-----------------------------------------------------------------------
!     5.5 can now calculate V'' for Mercier calculation
!     I am calculating individual coefficients below to generate
!     full peeling mode criterion--need to know how vpp changes
!     with s and alpha
!-----------------------------------------------------------------------
      vpp=0.
      vpp_c=0.
      vpp_pp=0.
      vpp_ffp=0.
      do i=2,npts
         dl=0.5*(leng(i)-leng(i-1))
         vpp_c=vpp_c+dl*( &
             1./(rpts(i)**2*bppts(i)**3)* &
             (2.*uprimepts(i)*rpts(i)*bppts(i))+ &
             1./(rpts(i-1)**2*bppts(i-1)**3)* &
             (2.*uprimepts(i-1)*rpts(i-1)*bppts(i-1)))
         vpp_pp=vpp_pp+dl*( &
             1./(rpts(i)**2*bppts(i)**3)* &
             rpts(i)**2+ &
             1./(rpts(i-1)**2*bppts(i-1)**3)* &
             rpts(i-1)**2)
         vpp_ffp=vpp_ffp+dl*( &
             1./(rpts(i)**2*bppts(i)**3)+ &
             1./(rpts(i-1)**2*bppts(i-1)**3))
         vpp=vpp+dl*( &
             1./(rpts(i)**2*bppts(i)**3)* &
             (2.*uprimepts(i)*rpts(i)*bppts(i)+ &
             rpts(i)**2*pprime_loc+ffprime_loc)+ &
             1./(rpts(i-1)**2*bppts(i-1)**3)* &
             (2.*uprimepts(i-1)*rpts(i-1)*bppts(i-1)+ &
             rpts(i-1)**2*pprime_loc+ffprime_loc))
      enddo
      vpp=2.*pi*vpp
      vpp_c=2.*pi*vpp_c
      vpp_pp=2.*pi*vpp_pp
      vpp_ffp=2.*pi*vpp_ffp
      
      if (verbose .ge. 3) write(outfile,'(1p,3(a,e12.4))')&
           'volume=',volume,'  vprime=', &
          vprime,'  vpp=',vpp,' check=', &
          vpp_c+pprime_loc*vpp_pp+ffprime_loc*vpp_ffp
! write info to dskpeel to reconstruct peeling mode in s-alpha space
!      write(65,'(1p,2e16.8)') pprime_loc,ffprime_loc
!      write(65,'(1p,2e16.8)') alpha_loc,sh_loc
!      write(65,'(1p,3e16.8)') ffprime_c,ffprime_s,ffprime_pp
!      write(65,'(1p,3e16.8)') vpp_c,vpp_pp,vpp_ffp
!      write(65,'(1p,3e16.8)') qp_s
!-----------------------------------------------------------------------
! 5.7 calculate j_edge/j_ave
!     for this we need area, <B^-2> and Integrate[Bp,dl]
!-----------------------------------------------------------------------
      area=0.
      intbp=0.
      vp=0.
      bsqinv=0.
      bsq=0.
      bpsqinv=0.
      bsqbpsqinv=0.
      do i=2,npts
         dl=(leng(i)-leng(i-1))
         bpoldsq=(f/rpts(i-1))**2+bppts(i-1)**2
         bpnewsq=(f/rpts(i  ))**2+bppts(i  )**2
         bsqinv=bsqinv+0.5*dl*(1./(bppts(i)*bpnewsq)+ &
             1./(bppts(i-1)*bpoldsq))
         bpsqinv=bpsqinv+0.5*dl*(1./(rpts(i)**2*bppts(i)**3)+ &
             1./(rpts(i-1)**2*bppts(i-1)**3))
         bsqbpsqinv=bsqbpsqinv+0.5*dl*(1./(rpts(i)**2*bppts(i)**3) &
             *bpnewsq+1./(rpts(i-1)**2*bppts(i-1)**3)*bpoldsq)
         bsq=bsq+0.5*dl*(1./bppts(i)*bpnewsq+ &
             1./bppts(i-1)*bpoldsq)
         vp=vp+0.5*dl*(1./bppts(i)+ &
             1./bppts(i-1))
         intbp=intbp+0.5*(bppts(i-1)+bppts(i))*dl
         area=area-0.5*(rpts(i-1)-rpts(i))*(zpts(i-1)+zpts(i))
      end do
      bpave=intbp/leng(npts)
      bsqinv=bsqinv/vp
      bsq=bsq/vp
      bpsqinv=bpsqinv/vp
      bsqbpsqinv=bsqbpsqinv/vp
      jedge_n1=-btnorm*area*(pprime_loc*f*bsqinv+ffprime_loc/f)/intbp
      jedge_n1pp=-btnorm*area*(pprime_loc*f*bsqinv)/intbp
      jedge_n1ffp=-btnorm*area*(ffprime_loc/f)/intbp
      jedge_n2=-btnorm*area*(pprime_loc*f+ffprime_loc/f*bsq)/intbp &
          /btnorm**2
      jedge_n2pp=-btnorm*area*(pprime_loc*f)/intbp &
          /btnorm**2
      jedge_n2ffp=-btnorm*area*(ffprime_loc/f*bsq)/intbp &
          /btnorm**2
      jedge_n3=-btnorm*area*rmajor**2* &
          (pprime_loc*f*bpsqinv+ffprime_loc/f*bsqbpsqinv)/intbp &
          *bpave*bpave/btnorm**2
      jedge_n3pp=-btnorm*area*rmajor**2* &
          (pprime_loc*f*bpsqinv)/intbp &
          *bpave*bpave/btnorm**2
      jedge_n3ffp=-btnorm*area*rmajor**2* &
          (ffprime_loc/f*bsqbpsqinv)/intbp &
          *bpave*bpave/btnorm**2
      if (verbose .ge. 4) then
         write(outfile,*) 'area=',area
         write(outfile,*) 'bsqinv=',bsqinv
         write(outfile,*) 'Integrate[Bp,dl]=',intbp
      endif
      
      jedgen1_arr(isurf)=jedge_n1
      jedgen2_arr(isurf)=jedge_n2
      if (verbose .ge. 2) then
         write(outfile,*) 'jedge_n1=',jedge_n1
         write(outfile,*) 'jedge_n2=',jedge_n2
         write(outfile,*) 'jedge_n3=',jedge_n3
      endif
!      open(unit=88,file='dsk_jedge',status='unknown',iostat=ios)
! 3/7/00  .jedge file should contain edge current only for the
!   outermost surface
      if (isurf==1) then
         jedgen1store=jedge_n1
         if (verbose .ge. 1)  then
            open(unit=88,file=runname(1:lrunname)//'.jedge', &
                 status='unknown',iostat=ios)
            write(88,'(1p,e15.6)') jedge_n1
            write(88,'(1p,2e15.6)') jedge_n1pp,jedge_n1ffp
            write(88,'(1p,e15.6)') jedge_n2
            write(88,'(1p,2e15.6)') jedge_n2pp,jedge_n2ffp
            write(88,'(1p,e15.6)') jedge_n3
            write(88,'(1p,2e15.6)') jedge_n3pp,jedge_n3ffp
         endif
!         close(88)
      endif
! 3/24/05  when we get to last surface, write max pedestal
!   region values of alpha,pprime and jedge_n1 to the
!   end of the .jedge file, along with edge jedge_n1
!   value and psi locations
      if (isurf==nxinterp) then
         alpha_max=alpha(1)
         pprime_max=pprime_adj(1)
         jedgen1_max=jedgen1_arr(1)
         jedgen2_max=jedgen2_arr(1)
         sh_min=shear(1)
         psi_alpha_max=psigrid(1)
         psi_pprime_max=psigrid(1)
         psi_jn1_max=psigrid(1)
         psi_jn2_max=psigrid(1)
         psi_sh_min=psigrid(1)
         sh_al_max=shear(1)
         sh_j1_max=shear(1)
         sh_j2_max=shear(1)
! ps 3/13 spline onto finer grid to get more accurate max values
!    for .jedge file
         call spline(-psigrid,alpha,nxinterp,-1.d30,-1.d30,aspline)
         call spline(-psigrid,pprime_adj,nxinterp,-1.d30,-1.d30,ppspline)
         call spline(-psigrid,shear,nxinterp,-1.d30,-1.d30,sspline)
         call spline(-psigrid,jedgen1_arr,nxinterp,-1.d30,-1.d30,j1spline)
         call spline(-psigrid,jedgen2_arr,nxinterp,-1.d30,-1.d30,j2spline)
         do i=0,399
            psifine=0.92+dble(i)*(psigrid(1)-.92)/400. 
!            if (verbose .ge. 6) write(*,*) 'psifine=',psifine
            call splint(-psigrid,alpha,aspline,nxinterp,-psifine, &
             varfine)
            if (varfine > alpha_max) then
               alpha_max=varfine
               psi_alpha_max=psifine               
            endif
            call splint(-psigrid,pprime_adj,ppspline,nxinterp,-psifine, &
                 varfine)
            if (abs(varfine) > abs(pprime_max)) then
               pprime_max=varfine
               psi_pprime_max=psifine
            endif
            call splint(-psigrid,shear,sspline,nxinterp,-psifine, &
                 varfine)
            if (varfine < sh_min) then
               sh_min=varfine
               psi_sh_min=psifine
            endif
            call splint(-psigrid,jedgen1_arr,j1spline,nxinterp,-psifine, &
                 varfine)
            if (varfine > jedgen1_max) then
               jedgen1_max=varfine
               psi_jn1_max=psifine
            endif
            call splint(-psigrid,jedgen2_arr,j2spline,nxinterp,-psifine, &
                 varfine)
            if (varfine > jedgen2_max) then
               jedgen2_max=varfine
               psi_jn2_max=psifine
            endif
         enddo
! determine shear at max of alpha,jedge_n1 and jedge_n2
         call splint(-psigrid,shear,sspline,nxinterp,-psi_alpha_max, &
                 sh_al_max)
         call splint(-psigrid,shear,sspline,nxinterp,-psi_jn1_max, &
                 sh_j1_max)
         call splint(-psigrid,shear,sspline,nxinterp,-psi_jn2_max, &
                 sh_j2_max)

!         do i=2,nxinterp
!! assigning max values for pedestal region, defined to be psi>0.92
!            if (psigrid(i) > 0.92) then
!               if (alpha(i) > alpha_max) then
!                  alpha_max=alpha(i)
!                  psi_alpha_max=psigrid(i)
!                  sh_al_max=shear(i)
!               endif
!               if (abs(pprime_adj(i)) > abs(pprime_max) ) then
!                  pprime_max=pprime_adj(i)
!                  psi_pprime_max=psigrid(i)
!               endif
!               if (jedgen1_arr(i) > jedgen1_max) then
!                  jedgen1_max=jedgen1_arr(i)
!                  psi_jn1_max=psigrid(i)
!                  sh_j1_max=shear(i)
!               endif
!               if (jedgen2_arr(i) > jedgen2_max) then
!                  jedgen2_max=jedgen2_arr(i)
!                  psi_jn2_max=psigrid(i)
!                  sh_j2_max=shear(i)
!               endif
!               if (shear(i) < sh_min) then
!                  sh_min=shear(i)
!                  psi_sh_min=psigrid(i)
!               endif
!            endif
!         enddo
         
         if (verbose .ge. 1) then
            write(88,*) 'pedestal alpha_max, pprime_max, jedge_n1_max, jedge_n1_sep'
            write(88,'(1p,4e15.6)') alpha_max,pprime_max,jedgen1_max,jedgen1_arr(1)
            write(88,*) 'psi locations of alpha_max,pprime_max,jedgen1_max,edge'
            write(88,'(1p,4e15.6)') psi_alpha_max,psi_pprime_max,psi_jn1_max, &
                 psigrid(1)
            write(88,*) 'pedestal jedge_n2_max, jedge_n2_sep, shear_min'
            write(88,'(1p,3e15.6)') jedgen2_max,jedgen2_arr(1),sh_min
            write(88,*) 'psi locations of jedge_n2_max, jedge_n2_sep, shear_min'
            write(88,'(1p,3e15.6)') psi_jn2_max,psigrid(1),psi_sh_min
            write(88,*) 'shear values at alpha_max, jedge_n1_max, jedge_n2 max'
            write(88,'(1p,3e15.6)') sh_al_max,sh_j1_max,sh_j2_max
            close(88)
         endif
      endif

!-----------------------------------------------------------------------
!     6.0 spline r,z, and rbp onto an equal arc length mesh
!-----------------------------------------------------------------------
      dl=circumf/(ns-1.)
      do i=1,ns
         leq(i)=(i-1)*dl
      end do
      leq(ns)=circumf
      call spline1d(rl     ,leq,ns,rpts     ,leng,npts,work)
      call spline1d(zl     ,leq,ns,zpts     ,leng,npts,work)
      call spline1d(bpl    ,leq,ns,bppts    ,leng,npts,work)
      call spline1d(uprimel,leq,ns,uprimepts,leng,npts,work)
      call spline1d(cosul  ,leq,ns,cosupts  ,leng,npts,work)
      call spline1d(sinul  ,leq,ns,sinupts  ,leng,npts,work)
!-----------------------------------------------------------------------
!     calculate distance to separatrix in rho and in psi
!-----------------------------------------------------------------------
      rhomax=1./(uprimel(1)+(rl(1)*pprime_loc+ffprime_loc/rl(1))/bpl(1))
      psisep=rhomax*bpl(1)*rl(1)
      if (verbose .ge. 3) then 
         write(outfile,'(4(a,1p,e12.4))') "uprimel(1)=",uprimel(1), &
              "  rl(1)=",rl(1),"  bpl(1)=",bpl(1)
         write(outfile,'(a,1p,e12.4,a,e12.4)') "rhomax=",rhomax,"  psisep=", &
              psisep
      endif
!-----------------------------------------------------------------------
!     define bishop's alpha and sigma
!-----------------------------------------------------------------------
      abishop=-2.*pprime_loc*r0bishop**2/bp0bishop
      sigbishop=(pprime_loc+ffprime_loc/rmajor**2)*rmajor*r0bishop/ &
          bp0bishop
      lambishop=abishop*wav/dav-sigbishop
      shbishop=rmajor*r0bishop*bp0bishop*qrefp/qv
      if (verbose .ge. 4) then
         write(outfile,'(1p,3(a,e12.4))') "abishop=",abishop, &
              "  sigbishop=",sigbishop,"  shbishop=",shbishop
         write(outfile,'(1p,3(a,e12.4))')'lambishop=',lambishop, &
              '  wav=',wav,'  dav=',dav
         write(outfile,*) "================end output from surfcalc, isurf=",isurf
      endif

! write line in *.surfdat file
      if (verbose .ge. 1) write(nsurfdat,*) &
           isurf,psigrid(isurf),qv,q_eq(isurf),f_eq(isurf), &
           ffprime_eq(isurf),alpha_loc,pprime_loc/(4.*pi),sh_loc,jedge_n1, &
           jedge_n2,jedge_n3


      deallocate( rcinvps,rpts,zpts,bppts, &
              work,uprimepts,cosupts,sinupts, &
              r,theta,dummy,leng )

    return

end subroutine surfcalc








