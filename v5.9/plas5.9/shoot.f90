      
subroutine shoot(gamr,gami,vcomp,gamr2,gami2)

!-----------------------------------------------------------------------
!  Shoots in x-space to generate eigenfunction
!  ainv_i is the inverse of (Q_i+S_i*aa_(i-1))S_i
!  Same version works for updown symmetric or antisymmetric
!-----------------------------------------------------------------------
      use elite_data, only: aa,ainv,mwindow,nmodes,kx,mindex,bb,fmat, &
           nowindow,nmwinhalf,xx,nedge,mmax,nx,m0,pmat,smat,qmat,outfile, &
           mmin,asurfu,gamsq,isurfu,asurfup,isurfup,bug,dw_vac,updownsym, &
           vacuum,nn,pi,tmatasym,splineeq,compression,vsurf1u,ddsurf1u, &
           vuarr,press_eq,rotation,gamim,rotnorm,tmat_det,ii,mind_new, &
           g_comp,gamre,verbose
      implicit none
      real gamr(mind_new*mwindow),gami(mind_new*mwindow), &
           gamr2(mind_new*mwindow),gami2(mind_new*mwindow)
!      real, dimension(:), allocatable :: gamr,gami,gamr2,gami2
!      real gamr(mwindow),gami(mwindow)!,vr(mwindow,mwindow),vi(mwindow,mwindow)
      real gamrsmall(mwindow),gamismall(mwindow)
      complex vcomp(mwindow,mwindow)
!      real gamr2(mwindow),gami2(mwindow)
      complex fnow(mwindow),flast(mwindow)
      complex bmat(mwindow,mwindow),cmat(mwindow,2*mwindow)
! temp test to solve cgesv problem
!      complex cmat1(mwindow,mwindow),cmat2(mwindow,mwindow), &
!           bmat_tmp(mwindow,mwindow)
      complex tmat(mwindow,mwindow),tmat_herm(mwindow,mwindow)
      complex asu(nmodes,nmodes)
      complex asup(nmodes,nmodes)
      complex ans(mwindow)
      complex rhs(mwindow,mwindow)
      complex tmat2(mwindow,mwindow)
      complex wmat(mwindow,mwindow)
!      real bmat_real(mwindow,mwindow),cmat_real(mwindow,2*mwindow)
!      real tmat3(mwindow,mwindow)
      real, dimension(:,:), allocatable :: tmat_real,tmat_rsym,vac_tmp, &
           vr,vi,dum
      complex vmat(mwindow,mwindow),umat(mwindow,mwindow)
      real x1,x2,x3
      real denom1,denom2,denom3
      real d1,d2,d3
      integer m,mp,i,klo,kup,j
      real t1,t2,t3,f1,f2,f3
      real rz2,rz3,rlamda
      real fac
      real ave,avereal,aveim
      integer mmm
      integer neig,iterp,kk
      integer ixv(kx),mwintopv(kx)
      integer ix
      real x
      integer mup,mlo,nmsp,mwintop,index
      integer msp
      real asym0,asym1,sumsym0,sumsym1,maxdiag,rndoff
      real sanorm0,sanorm1
      complex  & !tmat_c(mwindow/mindex,mwindow/mindex), &
!           tmat_c2(mwindow/mindex,mwindow/mindex), &
           det(2),work(mwindow/mindex)
      integer ipiv(mwindow),ipvt(mwindow/mindex)
      integer info
      complex, dimension(mwindow,mwindow) :: dd,ddinv,ddinvav,ddinvbv, &
           av,bv,tempmat
!      real Gamsqr,Gamsqim,asuri,asupri,asuii,asupii
      complex Gamma_c,Gamma2_c
!      real ddri,ddii,avri,avii,bvri,bvii
      real vrsmall(mwindow/mindex,mwindow/mindex), &
           vismall(mwindow/mindex,mwindow/mindex)
      complex vrtest(mwindow,mwindow)

      fnow=0.; flast=0.; bmat=0.; cmat=0.; tmat=0.
      asu=0.; asup=0.; ans=0.; tmat2=0.
      vmat=0.; umat=0.; tempmat=0. 
!      tmat_real=0.; tmat_rsym=0.

!      ii=cmplx(0.,1.)

! 6/01 PBS  This routine is now being called multiple times to iterate
!   gamsq and calculate growth rate.  Best to zero relevant arrays again.
      aa=0.; ainv=0.; bb=0.; fmat=0.


      if (verbose .ge. 4) write(outfile,*) &
           "==================begin output from shoot==========="
!-----------------------------------------------------------------------
!  Run check on x-mesh compatibility
!  int(x+0.5) increases by unit at integer values of x
!-----------------------------------------------------------------------
      if(.not.nowindow) then
         msp=nmwinhalf
         ix=int(xx(3)+0.5)-int(xx(1)+0.5)
         if (ix.ne.0) then
            write(6,*)'ERROR***first three x-mesh', &
                ' points must be<1, ix=',ix,'xx=',xx
            stop
         end if
      end if
!-----------------------------------------------------------------------
!  No. modes which interact with the edge:
!     msp is maximum number of modes on one side of center mode
!     mm(1)=mmax, the maximum mode number, which is in vacuum
!     m0=mres, the minimum mode number in vacuum
!     so nvac=mm(1)-m0+1
!     and nedge=nvac+msp; have already guarenteed in main that
!        msp>=nvac-1
!-----------------------------------------------------------------------
      if (verbose .ge. 4) write(outfile,*)' nedge=',nedge,' nmodes=',nmodes
!-----------------------------------------------------------------------
!  Initialise independent eigenvectors at x=Delta 
!           (the plasma/vac surface)
!  Only consider restricted "window" in m at one time
!     mindex*(2*msp+1) is maximum number of modes in window
!     recall mindex=1 for up/down symmetric and 2 for non up/down sym
!-----------------------------------------------------------------------
!      if(nowindow) then
!         nmodesmax=nmodes
!      else
!         nmodesmax=min0(mindex*(2*msp+1),nmodes)
!      endif
!      nmodesmax=mwindow
!      do  m=1,nmodesmax        
!         do  mp=1,nmodesmax    
!          aa(m,mp,1)=0.
!       end do
!      end do

!fmat already zeroed out      
      do  m=1,nedge             !m labels components of solution
       fmat(1,m,m)=1.           !each solution has different non-zero b.c.
      end do
!-----------------------------------------------------------------------
! we need various indexing arrays for the windowing feature 
! The ixv array:
!   ix=1 -> xx(i) is integer or first value above integer
!   ix=-1-> xx(i) is integer or first value above integer &&
!        close enough to the edge that last interval had
!        nmodes < 2*msp+1 
!   ix=0 -> xx(i) not integer or first value above integer
!-----------------------------------------------------------------------
!   rewrite 6/2001 PBS
!    now window is always the same size = mwindow
!    use ixv=1 only for those cases where the domain of the window
!      must actually shift (reset it to zero when window is up against
!       upper left or lower right and no shift is needed)
!  9/2001 with reversed shear use ix=-1 for cases where window must
!    shift up rather than down

      ixv=0
      mwintopv=mmax
      nmsp=mwindow

      if (mwindow.ne.nedge) then
         write(*,*) 'problem in shoot, mwindow.ne.nedge'
      endif
      if(.not.nowindow) then
         do i=2,nx
! caculate nmsp
            x=xx(i)
            mup=min0(m0+msp-int(x+0.5),mmax)
!            mlo=max0(m0-msp-int(x+0.5),mmin)
            ! don't let mup be less than mmin+mwindow/mindex
            if (mup < (mmin+nmsp/mindex-1) ) mup=mmin+nmsp/mindex-1
!            if (mlo > (mmax-nmsp/mindex+1) ) mlo=mmax-nmsp/mindex+1
           ! mup-mlo should always be 2*msp
            mlo=mup-2*msp
            if (mlo < mmin) then
               write(*,*) 'problem in shoot, mlo<mmin',mlo,mmin,i,mup
               stop
            endif
            if (i.eq.2) then
               if (verbose .ge. 4) then
                  write(outfile,*)' mup=',mup,' mlo=',mlo,' mindex=',mindex
                  write(outfile,*)' m0=',m0,' mmax=',mmax,' mmin=',mmin
                  write(outfile,*)' x=',x,' nmsp=',nmsp,' int(x+0.5)=',int(x+0.5)
               endif
            end if

            mwintopv(i)=mup     ! top m value in current m window
            ix=int(x+0.5)-int(xx(i-1)+0.5)
            if (ix.eq.1) then
     ! if window is up against upper left or previous window was already
     !   against lower right, then reset ix=0, don't move window
               if (mup.eq.mmax) ix=0
               if ((m0-msp-int(xx(i-1)+0.5)).le.mmin) ix=0
            end if
            if (ix.eq.-1) then
     ! 9/01 non-monotonic q now possible
     ! don't move window up if it's up against lower right, or
     !  if previous window was already up against upper left
               if (mlo.eq.mmin) ix=0
               if ((m0+msp-int(xx(i-1)+0.5)).ge.mmax) ix=0
            endif
            if ( (ix > 1) .or. (ix < -1) ) then
               write(*,*) 'problem in shoot, ix=',ix
               write(*,*) 'i=',i,'mup=',mup,'mlow=',mlo,'x=',x
               stop
            endif
            ixv(i)=ix
            if(ix.ne.0) then
               if (verbose .ge. 4) write(outfile,*) &
                    'ix=',ixv(i),'xx(i)=',xx(i),'nmsp=',nmsp, &
                   'mwintop=',mwintopv(i)
            endif
         end do
      endif
!-----------------------------------------------------------------------
!  calculate matrix manipulations required for calculating aa and bb:
!-----------------------------------------------------------------------
      iterp=0 ! iterp is index of xinterp for determining s,q, & p below
      do 100 i=2,nx
         ix=ixv(i)
         mwintop=mwintopv(i)
!         write(*,*) 'call splinematgen, i=',i
         if (.not.splineeq) then
            if (compression) then
               call comp_matgen(i,iterp,nmsp,mwintop) ! get compressional
                 ! and rotational contributions to ai,api,appi
            endif
            call matgen(i,iterp,nmsp,mwintop) ! calculates smat, qmat, pmat
         else
            write(*,*) 'splineeq not yet setup in complex version'
            stop
            call splinematgen(i,iterp,nmsp,mwintop) ! calculates smat, qmat, pmat
         endif

         cmat(:,1:nmsp)=-pmat
         cmat(:,nmsp+1:2*nmsp)=-smat
!  TEMP!!!
!         cmat1=-pmat
!         cmat2=-smat

         if(ix.eq.0) then
            bmat=qmat+matmul(smat,aa(:,:,i-1))
         elseif(ix.eq.1) then
            bmat=qmat   ! bmat will be q+s*aa
! PBS 6/25/01 fix for case where window "overruns" matrix, that is
!   where nmspv(i-1)>nmspv(i)
            do  mp=1,nmsp-mindex
               do  m=1,nmsp
                     do  kk=1,nmsp-mindex
                        bmat(m,mp)=bmat(m,mp)+ &
                             smat(m,kk)*aa(kk+mindex,mp+mindex,i-1)
                     end do
               end do
            end do
! 9/01 need to add back ix=-1 case for reverse shear
         elseif(ix.eq.-1) then
            bmat=qmat
            do mp=1+mindex,nmsp
               do m=1,nmsp
                  do kk=1+mindex,nmsp
                     bmat(m,mp)=bmat(m,mp)+ &
                          smat(m,kk)*aa(kk-mindex,mp-mindex,i-1)
                  enddo
               enddo
            enddo
         endif     
!         write(*,*) 'calling zgesv in shoot'
!         if (i==2) then
!            write(88,*) 'bmat=',bmat,' cmat=',cmat
!         endif
!         if (updownsym) then
!            bmat_real=real(bmat)
!            cmat_real=real(cmat)
!            call dgesv( nmsp,2*nmsp, bmat_real, mwindow, ipiv, &
!                 cmat_real,mwindow, info )
!            bmat=bmat_real
!            cmat=cmat_real
!         else
!! TEMP fix for cgesv

         call zgesv( nmsp,2*nmsp, bmat, mwindow, ipiv, cmat, &
              mwindow, info )
!            bmat_tmp=bmat
!            call cgesv( nmsp,nmsp, bmat, mwindow, ipiv, cmat1, &
!                 mwindow, info )
!            call cgesv( nmsp,nmsp, bmat_tmp, mwindow, ipiv, cmat2, &
!                 mwindow, info )
         if(info.ne.0) then
            write(6,*) "zgesv in shoot failed; info=",info
            stop
         endif

!         if (i==2) then
!            write(88,*) 'after cgesv, bmat=',bmat,' cmat=',cmat
!!         write(*,*) 'return from cgesv in shoot'
!         endif

         aa(:,:,i)=cmat(:,1:nmsp)            ! -(q+s*aa)^-1*p
         ainv(:,:,i)=cmat(:,nmsp+1:2*nmsp)   ! -(q+s*aa)^-1*s
 100  continue

!-----------------------------------------------------------------------
! Now, for each of the nedge independent solutions
!-----------------------------------------------------------------------
      do 200 neig=1,nedge
!         bb(:,1)=fmat(1,:,neig)
         bb(1:nedge,1)=fmat(1,:,neig)   
!-----------------------------------------------------------------------
!  shoot into the plasma core, out to xmax, calculating bb as you go
!  components of eigenvector labeled by m
!-----------------------------------------------------------------------
         do  i=2,nx
            ix=ixv(i)
            bb(nmsp+1:nmsp+mindex,i)=0.
            do  m=1,nmsp
               bb(m,i)=0.
               
              if (ix.eq.0) then
!                 bb(1:nmsp,i)=matmul(ainv(:,:,i),bb(1:nmsp,i-1))
                 do mp=1,nmsp
                    bb(m,i)=bb(m,i)+ainv(m,mp,i)*bb(mp,i-1)
                 end do
              elseif(ix.eq.1) then
                 do mp=1,nmsp
                    bb(m,i)=bb(m,i)+ainv(m,mp,i)*bb(mp+mindex,i-1)
                 end do
              elseif(ix.eq.-1) then
       ! 9/01  double check this
                 do mp=1+mindex,nmsp
                    bb(m,i)=bb(m,i)+ainv(m,mp,i)*bb(mp-mindex,i-1)
                 end do
              endif
           end do
        end do
!-----------------------------------------------------------------------
!  now return to the plasma surface, constructing the eigenfunction
!  as you go
!-----------------------------------------------------------------------
         flast=bb(1:nmsp,nx)
         do 180 i=nx-1,2,-1
            ix=ixv(i+1)
            do  m=1,nmsp
               fnow(m)=bb(m,i)
               if (ix.eq.0) then
                  do  mp=1,nmsp
                     fnow(m)=fnow(m)+aa(m,mp,i)*flast(mp)
                  end do
               elseif(ix.eq.1) then !ix=1
                  do mp=1+mindex,nmsp
                     fnow(m)=fnow(m)+aa(m,mp,i)*flast(mp-mindex)
                  end do
               elseif(ix.eq.-1) then
                  do mp=1,nmsp-mindex
                     fnow(m)=fnow(m)+aa(m,mp,i)*flast(mp+mindex)
                  end do
               end if
            end do
            if (i.ne.2) flast=fnow
 180     continue
!-----------------------------------------------------------------------
!  Store solution at 2 mesh points next to wall
!-----------------------------------------------------------------------
        if (verbose .ge. 4) write(outfile,*)'CHECK nedge=',nedge,' =nmsp=',nmsp
        fmat(2,:,neig)=fnow
        fmat(3,:,neig)=flast
!        do  m=1,nedge
!           fmat(2,m,neig)=fnow(m)
!           fmat(3,m,neig)=flast(m)
!        end do
 200  continue

!-----------------------------------------------------------------------
!     Form T-matrix from set of eigenvectors
!-----------------------------------------------------------------------
      x1=xx(3)
      x2=xx(1)
      x3=xx(2)
      denom1=1./((x1-x2)*(x1-x3))
      denom2=1./((x1-x2)*(x2-x3))
      denom3=1./((x1-x3)*(x2-x3))
      d1=(x2-x3)*denom1
      d2=(x1-2.*x2+x3)*denom2
      d3=(x2-x1)*denom3
! add inertia terms into boundary terms...
! calculate modified local Gamma values for rotational case
      if (rotation) then
!         rot_loc=alph*rotnorm(iterp+1)+(1.-alph)*rotnorm(iterp)
!         Gamsqr=gamsq-(gamim+nn*rotnorm(1))**2
!         Gamsqim=2*sqrt(gamsq)*(gamim+nn*rotnorm(1))
!         Gamma_c=cmplx(sqrt(gamsq),gamim+nn*rotnorm(1))
         Gamma_c=cmplx(gamre,gamim+nn*rotnorm(1))
         Gamma2_c=Gamma_c*Gamma_c
!         write(*,*) 'boundary, rotnorm=',rotnorm(1),'Gamma2_c=', &
!              Gamma2_c
      else
!         Gamsqr=gamsq
!         Gamsqim=0.
         Gamma_c=sqrt(gamsq)
         Gamma2_c=gamsq
      endif
      do m=1,mwindow
        do mp=1,mwindow
          asu(m,mp)=asurfu(m,mp)-Gamma2_c*isurfu(m,mp)
          asup(m,mp)=asurfup(m,mp)-Gamma2_c*isurfup(m,mp)
        end do
      end do

! add compressional surface terms if necessary
      if (compression.and. press_eq(1).ne.0. .and. g_comp.ne.0.) then
         if (verbose .ge. 2) write(outfile,*) &
              'including compressional surface terms'
         do m=1,nedge
            do mp=1,nedge
               dd(m,mp)=ddsurf1u(m,mp)+ &
                    Gamma2_c*vuarr(2,m,mp,1)
            enddo
         enddo
         av=Gamma2_c*vuarr(4,1:nedge,1:nedge,1)
         bv=Gamma2_c*vuarr(3,1:nedge,1:nedge,1)
! invert dd
         ddinv=0.
         do j=1,nedge
            ddinv(j,j)=1.
         enddo
!       ddinv=dd**-1
!       ddtemp=dd
         call zgesv(nedge,nedge,dd,nedge,ipiv,ddinv, &
              nedge,info)
         if(info.ne.0) then
            write(6,*) "cgesv in comp_matgen failed; info=",info
            write(*,*) 'i=',i,'iterp=',iterp,'xx(i)=',xx(i)
            stop
         endif

!         ddinvav=matmul(ddinv,(Gamsqr*vuarr(4,1:nedge,1:nedge,1)))
         ddinvav=matmul(ddinv,av)
!         ddinvbv=matmul(ddinv,(Gamsqr*vuarr(3,1:nedge,1:nedge,1)))
         ddinvbv=matmul(ddinv,bv)
         tempmat=vsurf1u(1:nedge,1:nedge)
         asu=asu+matmul(tempmat,ddinvav)
         asup=asup+matmul(tempmat,ddinvbv)
!         asu=asu+matmul(vsurf1u(1:nedge,1:nedge),ddinvav)
!         asup=asup+matmul(vsurf1u(1:nedge,1:nedge),ddinvbv)
      endif

      do  m=1,nedge
         do  mp=1,nedge
            tmat(m,mp)=0.
            rhs(m,mp)=0.
            do  kk=1,nedge
               tmat(m,mp)=tmat(m,mp)+asu(m,kk)*fmat(1,kk,mp)+ &
                    asup(m,kk)*(d1*fmat(3,kk,mp)+ &
                    d2*fmat(1,kk,mp)+d3*fmat(2,kk,mp))
            end do
         end do
      end do

!-----------------------------------------------------------------------
! symmetrize tmat--8/26/98
!-----------------------------------------------------------------------

      if (updownsym) then
         allocate( tmat_real(mwindow,mwindow), &
              tmat_rsym(mwindow,mwindow))
         tmat_real=0.; tmat_rsym=0.
         tmat_real=real(tmat)
      else
         allocate( tmat_real(2*mwindow,2*mwindow), &
              tmat_rsym(2*mwindow,2*mwindow) )
         tmat_real=0.; tmat_rsym=0.
         do m=1,mwindow
            do mp=1,mwindow
               tmat_real(2*m-1,2*mp-1)=real(tmat(m,mp))
               tmat_real(2*m,2*mp)=real(tmat(m,mp))
               tmat_real(2*m,2*mp-1)=aimag(tmat(m,mp))
               tmat_real(2*m-1,2*mp)=-aimag(tmat(m,mp))
            enddo
         enddo
      endif

      if (verbose .ge. 5) write(outfile,*) 'Checking symmetry of Tmat_real'
      call checksym(tmat_real,mwindow*mind_new,nedge*mind_new)

      do m=1,mind_new*nedge
        do mp=1,mind_new*nedge
         tmat_rsym(m,mp)=tmat_real(m,mp)
        end do
      end do
      if (bug(3)<10. .and. (updownsym) .and. (.not.rotation)) then
!  3/03 ps looks like this routine as written should not be used
!    for undown asymmetric cases - need to enforce Hermitian instead
!         write(*,*) 'SYMMETRIZING TMAT!!!'
!         write(outfile,*) 'SYMMETRIZING TMAT!!!'
         do m=2,nedge
!      do m=2,nmodes
            do mp=1,m-1
               ave=0.5*(tmat_real(m,mp)+tmat_real(mp,m))
!            write(*,*) 'matmmp=',tmat(m,mp),'tmatmpm=',tmat(mp,m),ave,mp,m
               tmat_rsym(m,mp)=ave
               tmat_rsym(mp,m)=ave
            end do
         end do
      else 
!         write(*,*) 'NOT symetrizing TMAT!!!'
!         write(outfile,*) 'NOT symetrizing TMAT!!!'
      endif

      tmat_herm=real(tmat)
      do m=2,nedge
         do mp=1,m-1
            avereal=0.5*(real(tmat(m,mp))+real(tmat(mp,m)))
            aveim=0.5*(aimag(tmat(m,mp))-aimag(tmat(mp,m)))
            tmat_herm(m,mp)=cmplx(avereal,aveim)
            tmat_herm(mp,m)=cmplx(avereal,-aveim)
         enddo
      enddo
!      write(89,*) 'before adding vac tmat_herm=',tmat_herm
!-----------------------------------------------------------------------
! For the moment, we write dW to give answer of original edge code
!  Since tmat actually needs to multiplied by pi/nn to be dW_s
!    we divide dw_vac by this factor to make them consistent
!  Dividing both by 2pi is a particular kinetic normalization
!-----------------------------------------------------------------------
! vac_tmp for old real matrices only
      allocate(vac_tmp(mind_new*nmodes,mind_new*nmodes))
      vac_tmp=0.
      do m=1,nmodes
         do mp=1,nmodes
!            vac_tmp(mindex*(m-1)+1,mindex*(mp-1)+1)=dw_vac(m,mp)
            vac_tmp(mind_new*(m-1)+1,mind_new*(mp-1)+1)=dw_vac(m,mp)
            if (.not.updownsym) &
                vac_tmp(mind_new*(m-1)+2,mind_new*(mp-1)+2)= &
                dw_vac(m,mp) 
         enddo
      enddo

      do m=1,nedge
         do mp=1,nedge
            if(vacuum) then
!               tmat(m,mp)=(nn*dw_vac(m+offset_v,mp+offset_v)/pi+ &
!                    tmat(m,mp))/(2.*pi)
               tmat(m,mp)=(nn*dw_vac(m,mp)/pi+ &
                    tmat(m,mp))/(2.*pi)
               tmat_real(m,mp)=(nn*vac_tmp(m,mp)/pi+ &
                    tmat_real(m,mp))/(2.*pi)
               tmat_rsym(m,mp)=(nn*vac_tmp(m,mp)/pi+ &
                    tmat_rsym(m,mp))/(2.*pi)
               tmat_herm(m,mp)=(nn*dw_vac(m,mp)/pi+ &
                    tmat_herm(m,mp))/(2.*pi)
            else
               tmat_real(m,mp)=tmat_real(m,mp)/(2.*pi)
               tmat_rsym(m,mp)=tmat_rsym(m,mp)/(2.*pi)
               tmat(m,mp)=tmat(m,mp)/(2.*pi)
               tmat_herm(m,mp)=tmat_herm(m,mp)/(2.*pi)
            endif
         end do
      end do

!-----------------------------------------------------------------------
! Finally solve eigenvalue equation for complex eigenvalue (gamr,gami)
!    and complex eigenvectors (vr,vi)
!-----------------------------------------------------------------------
      if (verbose .ge. 5) write(outfile,*) 'call eigen'
!      do m=1,nedge
!         gami(m)=0.
!         do mp=1,nedge
!            tmat2(m,mp)=tmat(m,mp)
!            vi(m,mp)=0.
!         end do
!      end do
!      tmat2=tmat
      allocate( vr(mind_new*mwindow,mind_new*mwindow), &
           vi(mind_new*mwindow,mind_new*mwindow), &
           dum(mind_new*mwindow,mind_new*mwindow))
      gamr=0. 
      vr=0. 
      vi=0. 
      gami=0.
      dum=0.
      tmat2=tmat

!      write(*,*) 'call eigen, mwindow=',mwindow,nedge,'tmat2=',tmat2
!      write(25,*) tmat2
      call eigen(tmat_rsym,mind_new*mwindow,mind_new*nedge,gamr,vr,dum)
      if (verbose .ge. 5) write(outfile,*) 'return from eigen'
      if (updownsym) then
         vcomp=vr
      else
         do j=1,mwindow
            do kk=1,mwindow
               vcomp(j,kk)=cmplx(vr(2*j-1,2*kk-1),vr(2*j,2*kk-1))
            enddo
         enddo
      endif
! 4/20/01 temporarily remove nag call until library is found
!      call eigen_nag(tmat3,mwindow,nedge,gamr2,gami2)
!      call eigen_lapack(tmat_real,mind_new*mwindow,mind_new*nedge,gamr2,gami2)
!      write(outfile,*) 'return from eigen_lapack'
! use complex routine all the time now
!      if (rotation) then
! use complex eigenvalue solver
!         gamr2=100.
!         gami2=0.
!         tmat_c2=tmat_c

! 4/03 first call eigen_comp with hermitian matrix, then with non
!         call eigen_comp(tmat_herm,mwindow,mwindow/mindex, &
!              gamrsmall,gamismall,vrsmall,vismall)
      call eigen_hcomp(tmat_herm,mwindow,mwindow/mindex, &
              gamrsmall,gamismall,vrsmall,vismall,info)
      if (info.ne.0) then
         write(*,*) 'eigen_hcomp failed in shoot, info=',info
         stop
      endif
         if (updownsym) then
            gamr=gamrsmall
            gami=gamismall
         else
            gamr(1:mwindow)=gamrsmall
            gamr(mwindow+1:2*mwindow)=gamrsmall
            gami(1:mwindow)=gamismall
            gami(mwindow+1:2*mwindow)=gamismall
         endif
         vcomp=cmplx(vrsmall,vismall)

! now call eigen_comp with un-hermitianized matrix
         vrsmall=0.; vismall=0.
!         call eigen_comp(tmat,mwindow,mwindow/mindex, &
!              gamrsmall,gamismall,vrsmall,vismall)
         call eigen_comp2(tmat,mwindow,mwindow/mindex, &
              gamrsmall,gamismall,vrsmall,vismall,info)
         if (info.ne.0) then
            write(*,*) 'eigen_comp2 failed in shoot, info=',info
            stop
         endif
         if (updownsym) then
            gamr2=gamrsmall
            gami2=gamismall
         else
            gamr2(1:mwindow)=gamrsmall
            gamr2(mwindow+1:2*mwindow)=gamrsmall
            gami2(1:mwindow)=gamismall
            gami2(mwindow+1:2*mwindow)=gamismall
         endif
         if (rotation) vcomp=cmplx(vrsmall,vismall)


!      endif


!         write(42,*) 'vr(:,1)=',vr(:,1),'vr(:,2)=',vr(:,2)
!         write(42,*) 'vrsmall(:,12)=',vrsmall(:,12),'vrsmall(:,2)=', &
!              vrsmall(:,2),'vrsmall(12,:)=',vrsmall(12,:)
!         write(42,*) 'vismall(:,12)=',vismall(:,12),'vismall(12,:)=', &
!              vismall(12,:)

!      write(42,*) 'vr=',vr,'vrsmall=',vrsmall,'vismall=',vismall, &
!           'vrtest=',vrtest
!  TEMP!!!!!!
!         vr=vrtest
!         vr=0.
!         do i=1,mwindow/mindex
!            vr(2*i-1,1:mwindow/mindex)=vrsmall(i,:)
!            vr(2*i,1:mwindow/mindex)=vismall(i,:)
!         enddo
        

! evaluate determinant of tmat_c for possible use in iterative
!    schemes
          !tmat_c2 is copy that will be factored
         call zgefa(tmat2,mwindow/mindex,mwindow/mindex, &
              ipvt,info)
         if (info.ne.0) then
            write(*,*) 'zgefa in shoot found zeros in u(k,k),', &
                 'cant calculate determinant'
            tmat_det=(0.,0.)
!            stop
         endif
         if (info.eq.0) then
            call zgedi(tmat2,mwindow/mindex,mwindow/mindex, &
                 ipvt,det,work,10)    ! calculate determinant
            tmat_det=det(1)*10.**det(2)
         endif


!         vr(1:mwindow/mindex,1:mwindow/mindex)=vrsmall
!      endif

!-----------------------------------------------------------------------
! check eigenvalues and eigenvectors
!-----------------------------------------------------------------------
      do m=1,nedge
         ans(m)=0.
         do mp=1,nedge
            ans(m)=ans(m)+tmat(m,mp)*vcomp(mp,nedge-1)
         end do
      end do
      if (verbose .ge. 4) then
         write(outfile,'(5e12.4)') (ans(m),m=1,nedge)
         write(outfile,'(5e12.4)') (vcomp(m,nedge-1),m=1,nedge)
      endif
!      write(74,*) nedge,x1,x2,x3
!      write(74,'(5e20.12)') ((asurfu(i,j),i=1,nmodes),j=1,nmodes)
!      write(74,'(5e20.12)') ((asurfup(i,j),i=1,nmodes),j=1,nmodes)
!      write(74,'(5e20.12)') ((fmat(1,i,j),i=1,nmodes),j=1,nmodes)
!      write(74,'(5e20.12)') ((fmat(2,i,j),i=1,nmodes),j=1,nmodes)
!      write(74,'(5e20.12)') ((fmat(3,i,j),i=1,nmodes),j=1,nmodes)
!      write(74,'(5e20.12)') ((tmat(i,j),i=1,nmodes),j=1,nmodes)
!      write(74,'(5e20.12)') ((dw_vac(i,j),i=1,nmodes),j=1,nmodes)
      if (verbose .ge. 5) write(outfile,*) "return from  eigen"

!      write(*,*) 'in shoot gamr=',gamr,' gamr2=',gamr2

    return
      
end subroutine shoot


subroutine checksym(checkmat,ksize,size)

      use elite_data, only: outfile,tmatasym,verbose
      integer size,ksize
      real checkmat(ksize,ksize)
      real asym0,asym1,sumsym0,sumsym1,maxdiag,rndoff
      real sanorm0,sanorm1

! 1/4/01 evaluate asymmetry in matrix using gato's algorithm (for dW_vac)
      maxdiag=0.
      do i=1,size
         if (abs(checkmat(i,i)).gt.maxdiag) maxdiag=abs(checkmat(i,i))
      enddo
      rndoff   = maxdiag*maxdiag/1.d6
!      write(*,*) 'using rndoff of maxdiag^2/1.d6=',rndoff
      sumsym0  = 0.0
      sanorm0  = 0.0
      sumsym1  = 0.0
      sanorm1  = 0.0
      do j=1,size  ! loop to calculate gato-like asymmetry 
         do i=j,size
            if(i .eq. j) sanorm0 = 2.0* checkmat(i,j)*checkmat(i,j)
            if(i .ne. j) sanorm0 = (checkmat(i,j)*checkmat(i,j) + &
                 checkmat(j,i)*checkmat(j,i))
            if(abs(sanorm0) .le. rndoff) sanorm0 = rndoff
            sumsym0 = sumsym0 + 2.0*(checkmat(j,i)-checkmat(i,j))* &
                 (checkmat(j,i)-checkmat(i,j))/sanorm0

            if(i .eq. j) sanorm1 = sanorm1 +  checkmat(i,j)*checkmat(i,j)
            if(i .ne. j) sanorm1 = sanorm1+(checkmat(i,j)*checkmat(i,j) + &
                 checkmat(j,i)*checkmat(j,i))
            sumsym1  = sumsym1  +  2.0*(checkmat(j,i)-checkmat(i,j))* &
                 (checkmat(j,i)-checkmat(i,j))
         enddo
      enddo

      asym0   = sqrt(sumsym0) / size
      asym1   = sqrt(sumsym1/sanorm1)


!      write(*,7020) asym0,asym1
      if (verbose .ge. 2) write(outfile,7020) asym0,asym1
 7020 format(1x,'Matrix asymmetry (first  norm) is: ',e16.6 &
           ,/,1x,'Matrix asymmetry (second norm) is: ',e16.6,/)
      if (asym0.gt.0.5 .and. asym1.gt.0.5) then
         write (*,*) '!!!WARNING, both asym0 and asym1 > 50%'
         if (verbose .ge. 2) write (outfile,*) &
              '!!!WARNING, both asym0 and asym1 > 50%'
      else if (asym0.gt.0.1 .and. asym1.gt.0.1) then
         if (verbose .ge. 2) write (*,*) '!!!WARNING, both asym0 and asym1 exceed 10%'
         if (verbose .ge. 2) write (outfile,*) 'WARNING, both asym0 and asym1 exceed 10%'
      else if (asym0.gt.0.01 .and. asym1.gt.0.01) then 
!         write (*,*) '!WARNING, one of the asymmetry norms exceeds 1%'
!         write (outfile,*) '!WARNING, one of the asymmetry norms exceeds 1%'
      endif
      tmatasym=asym1

end subroutine checksym



