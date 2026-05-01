      
subroutine fcal

!-----------------------------------------------------------------------
!  Shoots in x-space to generate eigenfunction ffun(m,i)
!  (m labels the harmonic, i the x mesh point)
!  from the m independent eigenfunctions, weighted with the
!  eigenvector aval(m).
!  ainv_i is the inverse of (Q_i+S_i*aa_(i-1))
!  times S_i
!  Same version works for updown symmetric or antisymmetric
!-----------------------------------------------------------------------
      use elite_data, only: mwindow,kx,nmwinhalf,ffun,fmat,mmax, &
           nowindow,nx,xx,m0,mmin,mindex,outfile,nedge,ainv,aa,aval, &
           nmodes,bb,verbose
      implicit none
      complex fnow(mwindow),flast(mwindow)
      complex tmat(mwindow,mwindow)
      complex vmat(mwindow,mwindow),umat(mwindow,mwindow)
      integer i,j,neig,m,mp,ix
      integer index,mup,mlo,nmsp,msp,mup1,mv,offset
      integer ixv(kx),mwintopv(kx)
      real x
      msp=nmwinhalf
!-----------------------------------------------------------------------
!  Initialise eigenfunction
!-----------------------------------------------------------------------

      ffun=0.; fnow=0.; flast=0.; tmat=0.; vmat=0.; umat=0.

      fmat=0.
      ixv=0
      mwintopv=mmax
      nmsp=mwindow

      if(.not.nowindow) then
         do i=2,nx
! caculate nmsp
            x=xx(i)
            mup=min0(m0+msp-int(x+0.5),mmax)
          ! don't let mup be less than mmin+mwindow/mindex
            if (mup < (mmin+nmsp/mindex-1) ) mup=mmin+nmsp/mindex-1
            mlo=mup-2*msp
            if (mlo < mmin) then
               write(*,*) 'problem in fcal, mlo<mmin',mlo,mmin,i,mup
               stop
            endif
            if (i.eq.2) then
               if (verbose .ge. 3) then
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
               write(*,*) 'problem in fcal, ix=',ix
               write(*,*) 'i=',i,'mup=',mup,'mlow=',mlo,'x=',x
               stop
            endif
            ixv(i)=ix
            if(ix.ne.0) then
               if (verbose .ge. 3) write(outfile,*) &
                    'ix=',ixv(i),'xx(i)=',xx(i),'nmsp=',nmsp, &
                   'mwintop=',mwintopv(i)

            endif
         end do
      endif

!-----------------------------------------------------------------------
!  Initialise independent eigenvectors at x=Delta (the plasma/vac
!  surface):
!  mp labels the eigenvector, m are the components of that vector:
!-----------------------------------------------------------------------
      do  m=1,nedge
         fmat(1,m,m)=1.
      end do
!-----------------------------------------------------------------------
! Now, for each independent eigenvector
!-----------------------------------------------------------------------
      do 200 neig=1,nedge
         do  m=1,nedge
            bb(m,1)=fmat(1,m,neig) !Load initial bb-matrix
         end do
         do index=1,mindex
            bb(nedge+index,1)=0
         end do
!-----------------------------------------------------------------------
!  shoot into the plasma core, out to xmax, calculating bb as you go
!-----------------------------------------------------------------------
        do 120 i=2,nx
           do index=1,mindex
             bb(nmsp+index,i)=0.
          end do
          ix=ixv(i)
          do  m=1,nmsp       !  components of eigenvector labeled by m
             bb(m,i)=0.
             if (ix.eq.0) then
                do mp=1,nmsp
                   bb(m,i)=bb(m,i)+ainv(m,mp,i)*bb(mp,i-1)
                end do
             elseif(ix.eq.1) then
                do mp=1,nmsp
                   bb(m,i)=bb(m,i)+ainv(m,mp,i)*bb(mp+mindex,i-1)
                end do
             else if (ix.eq.-1) then
                do mp=1+mindex,nmsp
                   bb(m,i)=bb(m,i)+ainv(m,mp,i)*bb(mp-mindex,i-1)
                end do
             endif
          end do
120       continue
!-----------------------------------------------------------------------
!  now return to the plasma surface, constructing the eigenfunction
!  as you go
!-----------------------------------------------------------------------
        do 180 i=nx,2,-1
           x=xx(i)
           mup=mwintopv(i)
           mlo=mup-nmsp/mindex+1
!           if (mlo == mmin) write(*,*) 'mlo=mmin=',mlo,mup,mup-mlo,nmsp
           if (mlo < mmin) then
              write(*,*) 'problem in fcal mlo<mmin',mlo,mmin,mup,nmsp
              stop
           endif
           if (i.eq.nx) then
!     store solution vector in flast array:
              do  m=1,nmsp
                 flast(m)=bb(m,i)
              end do
           else
             ix=ixv(i+1)
!  calculate new fm fnow, from last f, flast:
             do m=1,nmsp
                fnow(m)=bb(m,i)
                if (ix.eq.0) then
                   do  mp=1,nmsp
                      fnow(m)=fnow(m)+aa(m,mp,i)*flast(mp)
                   end do
                elseif (ix.eq.1) then
                   do mp=1+mindex,nmsp
                      fnow(m)=fnow(m)+aa(m,mp,i)*flast(mp-mindex)
                   end do
                elseif (ix.eq.-1) then
                   do mp=1,nmsp-mindex
                      fnow(m)=fnow(m)+aa(m,mp,i)*flast(mp+mindex)
                   end do
                end if
             end do
             do  m=1,nmsp
                flast(m)=fnow(m)
             end do
          end if
!  Calculate full eigenfunction
          do m=1,mwindow
             offset=(mmax-mup)*mindex
             ffun(m+offset,i)=ffun(m+offset,i)+aval(neig)*flast(m)
          end do
 180   continue
 200  continue
!  Calculate eigenfunction at edge mesh point
      do  m=1,nedge
         ffun(m,1)=aval(m)
!         write(6,*) 'm=',m,'aval(m)=',aval(m)
         if (verbose .ge. 3) write(outfile,*) 'm=',m,'aval(m)=',aval(m)
      end do
      if(nedge.lt.nmodes) then
         do  m=nedge+1,nmodes
            ffun(m,1)=0.
         end do
      endif
    return
      
end subroutine fcal




