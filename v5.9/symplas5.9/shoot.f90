      
subroutine shoot(gamr,gami,vr,vi,gamr2,gami2)

!-----------------------------------------------------------------------
!  Shoots in x-space to generate eigenfunction
!  ainv_i is the inverse of (Q_i+S_i*aa_(i-1))S_i
!  Same version works for updown symmetric or antisymmetric
!-----------------------------------------------------------------------
      use elite_data, only: aa,ainv,mwindow,nmodes,kx,mindex,bb,fmat, &
           nowindow,nmwinhalf,xx,nedge,mmax,nx,m0,pmat,smat,qmat,outfile, &
           mmin,asurfu,gamsq,isurfu,asurfup,isurfup,bug,dw_vac,updownsym, &
           vacuum,nn,pi,tmatasym,splineeq,verbose
      implicit none
      real gamr(mwindow),gami(mwindow),vr(mwindow,mwindow),vi(mwindow,mwindow)
      real gamr2(mwindow),gami2(mwindow)
      real fnow(mwindow),flast(mwindow)
      real bmat(mwindow,mwindow),cmat(mwindow,2*mwindow)
      real tmat(mwindow,mwindow)
      real asu(nmodes,nmodes)
      real asup(nmodes,nmodes)
      real ans(mwindow)
      real rhs(mwindow,mwindow)
      real tmat2(mwindow,mwindow),wmat(mwindow,mwindow)
      real tmat3(mwindow,mwindow)
      real vmat(mwindow,mwindow),umat(mwindow,mwindow)
      real x1,x2,x3
      real denom1,denom2,denom3
      real d1,d2,d3
      integer m,mp,i,klo,kup,j
      real t1,t2,t3,f1,f2,f3
      real rz2,rz3,rlamda
      real fac
      real ave
      integer mmm
      integer neig,iterp,kk
      integer ixv(kx),mwintopv(kx)
      integer ix
      real x
      integer mup,mlo,nmsp,mwintop,index
      integer msp
      real vac_tmp(nmodes,nmodes)
      real asym0,asym1,sumsym0,sumsym1,maxdiag,rndoff
      real sanorm0,sanorm1
      complex tmat_c(nmodes/mindex,nmodes/mindex),ii
      integer ipiv(mwindow)
      integer info

      fnow=0.; flast=0.; bmat=0.; cmat=0.; tmat=0.
      asu=0.; asup=0.; ans=0.; tmat2=0.; tmat3=0.
      vmat=0.; umat=0.; vac_tmp=0.; tmat_c=0.

      ii=(0.,1.)

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

      if (verbose .ge. 3) write(outfile,*) 'mwindow=',mwindow,'nedge=',nedge
      if (mwindow.ne.nedge) then
         if (verbose .ge. 1) write(*,*) 'problem in shoot, mwindow.ne.nedge'
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
      if (verbose .ge. 5) write(outfile,*) 'start matgen loop'
      do 100 i=2,nx
         ix=ixv(i)
         mwintop=mwintopv(i)
!         write(*,*) 'call splinematgen, i=',i
         if (.not.splineeq) then
            call matgen(i,iterp,nmsp,mwintop) ! calculates smat, qmat, pmat
         else
            call splinematgen(i,iterp,nmsp,mwintop) ! calculates smat, qmat, pmat
         endif

         cmat(:,1:nmsp)=-pmat
         cmat(:,nmsp+1:2*nmsp)=-smat
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
         call dgesv( nmsp,2*nmsp, bmat, mwindow, ipiv, cmat, mwindow, info )
         if(info.ne.0) then
            write(6,*) "dgesv in shoot failed; info=",info
            stop
         endif


         aa(:,:,i)=cmat(:,1:nmsp)            ! -(q+s*aa)^-1*p
         ainv(:,:,i)=cmat(:,nmsp+1:2*nmsp)   ! -(q+s*aa)^-1*s
 100  continue
         if (verbose .ge. 5) write(outfile,*) 'end matgen call loop'

!-----------------------------------------------------------------------
! Now, for each of the nedge independent solutions
!-----------------------------------------------------------------------
!      write(*,*) 'starting shooting in and out loop'
      do 200 neig=1,nedge
         bb(:,1)=fmat(1,:,neig)
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
!$$$            if(ix.ne.0) then
!$$$               write(6,*) 'xx=',xx(i),'ix=',ix
!$$$               write(6,*) 'i=',i
!$$$               write(6,*) (bb(mmm,i),mmm=1,2*msp+2)
!$$$               write(6,*) 'i=',i-1
!$$$               write(6,*) (bb(mmm,i-1),mmm=1,2*msp+2)
!$$$               write(6,*) 'i=',i+1
!$$$               write(6,*) (bb(mmm,i+1),mmm=1,2*msp+2)
!$$$            end if
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
        if (verbose .ge. 5) write(outfile,*) 'end of shooting in and out loop'

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
      do m=1,mwindow
        do mp=1,mwindow
          asu(m,mp)=asurfu(m,mp)-gamsq*isurfu(m,mp)
          asup(m,mp)=asurfup(m,mp)-gamsq*isurfup(m,mp)
        end do
      end do

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
      do m=1,nedge/mindex
         do mp=1,nedge/mindex
            tmat_c(m,mp)=tmat(m*mindex-1,mp*mindex-1)+ &
                 ii*tmat(m*mindex-1,mp*mindex) 
!            write(*,*) 'tmat_c(m,mp)=',tmat_c(m,mp),m,mp
         enddo
      enddo
!      do m=1,nedge/mindex
!         do mp=1,m
!            write(73,*) 'tmat_c(m,mp)=',tmat_c(m,mp), &
!                 'tmat_c(mp,m)=',tmat_c(mp,m),m,mp
!         enddo
!      enddo

!      write(*,*) 'Checking symmetry of Tmat'
      if (verbose .ge. 5) write(outfile,*) 'Checking symmetry of Tmat'
      call checksym(tmat,mwindow,nedge)

      do m=1,nedge
        do mp=1,nedge
         tmat3(m,mp)=tmat(m,mp)
        end do
      end do
      if (bug(3)<10.) then
!         write(*,*) 'SYMMETRIZING TMAT!!!'
         if (verbose .ge. 4) write(outfile,*) 'SYMMETRIZING TMAT!!!'
         do m=2,nedge
!      do m=2,nmodes
            do mp=1,m-1
               ave=0.5*(tmat(m,mp)+tmat(mp,m))
!            write(*,*) 'matmmp=',tmat(m,mp),'tmatmpm=',tmat(mp,m),ave,mp,m
               tmat(m,mp)=ave
               tmat(mp,m)=ave
            end do
         end do
      else 
         if (verbose .ge. 2) then
            write(*,*) 'NOT symetrizing TMAT!!!'
            write(outfile,*) 'NOT symetrizing TMAT!!!'
         endif
      endif


!-----------------------------------------------------------------------
! For the moment, we write dW to give answer of original edge code
!  Since tmat actually needs to multiplied by pi/nn to be dW_s
!    we divide dw_vac by this factor to make them consistent
!  Dividing both by 2pi is a particular kinetic normalization
!-----------------------------------------------------------------------
      vac_tmp=0.
      do m=1,nmodes/mindex
         do mp=1,nmodes/mindex
!            vac_tmp(mindex*(m-1)+1,mindex*(mp-1)+1)=dw_vac(m,mp)
            vac_tmp(mindex*(m-1)+1,mindex*(mp-1)+1)=dw_vac(m,mp)
            if (.not.updownsym) &
                vac_tmp(mindex*(m-1)+2,mindex*(mp-1)+2)=dw_vac(m,mp) 
         enddo
      enddo

      do m=1,nedge
         do mp=1,nedge
            if(vacuum) then
!               tmat(m,mp)=(nn*dw_vac(m+offset_v,mp+offset_v)/pi+ &
!                    tmat(m,mp))/(2.*pi)
               tmat(m,mp)=(nn*vac_tmp(m,mp)/pi+ &
                    tmat(m,mp))/(2.*pi)
               tmat3(m,mp)=(nn*vac_tmp(m,mp)/pi+ &
                    tmat3(m,mp))/(2.*pi)
            else
               tmat(m,mp)=tmat(m,mp)/(2.*pi)
               tmat3(m,mp)=tmat3(m,mp)/(2.*pi)
            endif
         end do
      end do

!-----------------------------------------------------------------------
! Finally solve eigenvalue equation for complex eigenvalue (gamr,gami)
!    and complex eigenvectors (vr,vi)
!-----------------------------------------------------------------------
      if (verbose .ge. 5) write(outfile,*) 'call eigen'
      do m=1,nedge
         gami(m)=0.
         do mp=1,nedge
            tmat2(m,mp)=tmat(m,mp)
            vi(m,mp)=0.
         end do
      end do
      gamr=0.; vr=0.
      bmat=0.

!      write(*,*) 'call eigen, mwindow=',mwindow,nedge,'tmat2=',tmat2
!      write(25,*) tmat2
      call eigen(tmat2,mwindow,nedge,gamr,vr,bmat)
      if (verbose .ge. 5) write(outfile,*) 'return from eigen'
!   temporarily remove call when not needed to check performance
      call eigen_lapack(tmat3,mwindow,nedge,gamr2,gami2)
      if (verbose .ge. 5) write(outfile,*) 'return from eigen_lapack'
!-----------------------------------------------------------------------
! check eigenvalues and eigenvectors
!-----------------------------------------------------------------------
      do m=1,nedge
         ans(m)=0.
         do mp=1,nedge
            ans(m)=ans(m)+tmat(m,mp)*vr(mp,nedge-1)
         end do
      end do
      if (verbose .ge. 4) then
         write(outfile,'(5e12.4)') (ans(m),m=1,nedge)
         write(outfile,'(5e12.4)') (vr(m,nedge-1),m=1,nedge)
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
    return
      
end subroutine shoot


subroutine checksym(checkmat,ksize,size)

      use elite_data, only: outfile,tmatasym
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
         if (verbose .ge. 2) write (outfile,*) '!!!WARNING, both asym0 and asym1 > 50%'
      else if (asym0.gt.0.1 .and. asym1.gt.0.1) then
         if (verbose .ge. 2) write (*,*) '!!!WARNING, both asym0 and asym1 exceed 10%'
         if (verbose .ge. 2) write (outfile,*) 'WARNING, both asym0 and asym1 exceed 10%'
      else if (asym0.gt.0.01 .and. asym1.gt.0.01) then 
!         write (*,*) '!WARNING, one of the asymmetry norms exceeds 1%'
!         write (outfile,*) '!WARNING, one of the asymmetry norms exceeds 1%'
      endif
      tmatasym=asym1

end subroutine checksym



