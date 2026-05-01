
subroutine splinematgen(i,iterp,nmodeswin,mwintop)

!-----------------------------------------------------------------------
!     Calculates P, Q and S matrices at each x-point
!     smat*U[x-dx]+qmat*U[x]+pmat*U[x+dx]
!     Same version works for updown symmetric or antisymmetric
!     modified 5/19/98 to allow windowing
!     nmodeswin is number of modes carried in current m window
!     mwintop is maximum m visible in current m window
!-----------------------------------------------------------------------
! 1/02 new version uses splines of the au* and iu* arrays along
!    psi (or x) rather than linear interpolation

      use elite_data, only: psixx,psigrid,nxinterp,xx,m0,auarr, &
           iuarr,auspline,iuspline,nx,nowindow,mmax,nmwinhalf, &
           maxdxinterp,mmin,mwindow,mindex,gamsq,smat,qmat,pmat
      implicit none
      integer i,iterp
      integer nmodeswin,mwintop
      integer j,kk,ii,m,mp
      real x1,x2,x3
      real d1,d2,d3
      real dd1,dd2,dd3
      real denom1,denom2,denom3
      integer mtop,offset,jfull,kkfull,jbig,kkbig,bigoffset
      real auloc(19),iuloc(7:9)
      real dmnq,dmnq2,dmpnq,dmpnq2,dmpnq2m,dmpnqm,dmpnqm2
      real ai,api,appi

!-----------------------------------------------------------------------
!  i labels x-point used
!  iterp labels the sparser interpolation points
!   note iterp is not sent =0

!-----------------------------------------------------------------------
!     2.0 use 3 point differencing to compute d/dx and d2/dx2
!     f'(x2)=d1*f(x1)+d2*f(x2)+d3*f(x3)
!     f''(x2)=dd1*f(x1)+dd2*f(x2)+dd3*f(x3)
!     formulas have been tested
!-----------------------------------------------------------------------
! 9/01 there is now potential for these denominators to blow up where
!   the shear goes through zero, unlikely to cause problems
       x2=xx(i)
       if (i.eq.nx) then
          x1=xx(i-1)
          x3=xx(i-2)
       else if (i.eq.1) then
          x3=xx(i+1)
          x1=xx(i+2)
       else
          x1=xx(i-1)
          x3=xx(i+1)
       end if
       denom1=1./((x1-x2)*(x1-x3))
       denom2=1./((x1-x2)*(x2-x3))
       denom3=1./((x1-x3)*(x2-x3))
       d1=(x2-x3)*denom1
       d2=(x1-2.*x2+x3)*denom2
       d3=(x2-x1)*denom3
       dd1=2.*denom1
       dd2=-2.*denom2
       dd3=2.*denom3

      mtop=mmax
      if(.not.nowindow) then
         ! note that mtop is one above usual mup.  need to leave room
         !  for offset in matgen
         mtop=min0(m0+nmwinhalf-int(xx(i)+0.5)+maxdxinterp,mmax)
         if (mtop < (mmin+mwindow/mindex-1+maxdxinterp)) &
              mtop=mmin+mwindow/mindex-1+maxdxinterp
         mtop=min0(mtop,mmax)
      endif

      bigoffset=mindex*(mmax-mwintop)
      offset=mindex*(mtop-mwintop)
      if ( (offset < 0) .or. &
           (offset > 2*maxdxinterp*mindex) ) then
         write(*,*) 'problem in matgen, offset=',offset,&
              'mtop=',mtop,'mmax=',mmax,'mmin=',mmin, &
              'mwintop=',mwintop,'i=',i
         stop
      endif

       do  j=1,nmodeswin          ! j labels rows of matrices
          jfull=j+offset
          jbig=j+bigoffset
!          mp=mmax-(jfull-1)/mindex  !if updownsym mindex=1
          mp=mtop-(jfull-1)/mindex
          dmpnq=xx(i)+mp-m0
          dmpnq2=dmpnq*dmpnq
!          write(6,*) "mp=",mp,"dmpnq=",dmpnq
          do  kk=1,nmodeswin          ! kk labels columns of matrices
             kkfull=kk+offset
             kkbig=kk+bigoffset
!             m=mmax-(kkfull-1)/mindex !if updownsym mindex=1
             m=mtop-(kkfull-1)/mindex !if updownsym mindex=1
             dmnq=xx(i)+m-m0
             dmnq2=dmnq*dmnq
             dmpnqm=dmpnq*dmnq
             dmpnqm2=dmpnq*dmnq2
             dmpnq2m=dmpnq2*dmnq

             do ii=1,19
                !!!!!assuming no windowing for now!!!!!!!!!!!!!!!!!
                !!!need to fix auspline indices for windowing!!!!!
                call splint(-psigrid,auarr(ii,jfull,kkfull,:), &
                  auspline(ii,jbig,kkbig,:), &
                  nxinterp,-psixx(i),auloc(ii))
             enddo
             do ii=7,9
                !!!!!assuming no windowing for now!!!!!!!!!!!!!!!!!
                !!!need to fix iuspline indices for windowing!!!!!
                call splint(-psigrid,iuarr(ii,jfull,kkfull,:), &
                  iuspline(ii,jbig,kkbig,:), &
                  nxinterp,-psixx(i),iuloc(ii))
             enddo

             
             ai   =auloc(1)*dmpnqm &
                 +auloc(2)*dmpnqm2 &
                 +auloc(3)*dmpnq2m &
                 +auloc(4)*dmnq2 &
                 +auloc(5)*dmpnq2 &
                 +auloc(6)*dmnq &
                 +auloc(7)*dmpnq &
                 +auloc(8)

             api=auloc(9)*dmpnqm &
                 +auloc(10)*dmpnqm2 &
                 +auloc(11)*dmpnq2m &
                 +auloc(12)*dmnq2 &
                 +auloc(13)*dmpnq2 &
                 +auloc(14)*dmnq &
                 +auloc(15)*dmpnq &
                 +auloc(16)
              appi =auloc(17)*dmpnqm &
                 +auloc(18)*dmpnqm2 &
                 +auloc(19)*dmpnq2m - &
! add in inertia terms...
                 gamsq*(iuloc(9) &
                 +iuloc(8)*dmnq &
                 +iuloc(7)*dmpnq) 
             smat(j,kk)= d1*api+dd1*appi
             qmat(j,kk)= ai+d2*api+dd2*appi
             pmat(j,kk)= d3*api+dd3*appi

          end do
       enddo

       return


end subroutine splinematgen




subroutine splinearr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2/11/02 splines au* and iu* arrays and stores spline coeffs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       use elite_data, only: nxinterp,mmax,xinterp,m0,nmwinhalf,mmin, &
            mindex,nowindow,auarr,maxdxinterp,mwindow,psigrid, &
            auspline,nmodes,iuarr,iuspline
       implicit none

       integer i,j,k,mi,mj
       integer iterp,mtop(nxinterp),offset,mlo(nxinterp)
       real autemp(19,nxinterp),iutemp(7:9,nxinterp),bnd_set
      
      mtop=mmax
      mlo=mmin

! determine range mtop-mlo of au* arrays if windowing is used
      if(.not.nowindow) then
         do iterp=1,nxinterp
         ! note that mtop is maxdxinterp above usual mup.  need to leave room
         !  for offset in matgen
            mtop(iterp)=min0(m0+nmwinhalf-int(xinterp(iterp)+0.5) &
                 +maxdxinterp,mmax)
            if (mtop(iterp) < (mmin+mwindow/mindex-1+maxdxinterp) ) &
                 mtop(iterp)=mmin+mwindow/mindex-1+maxdxinterp
            mtop(iterp)=min0(mtop(iterp),mmax)
            mlo(iterp)=mtop(iterp)-mwindow/mindex+1-2*maxdxinterp  ! reduce by two 
                        !  to allow room
            if (mlo(iterp) < mmin) mlo(iterp)=mmin
         enddo
      endif

! spline each i,j element in au* matrices vs psigrid, store in auspline
      bnd_set=-2.d30
      do i=1,nmodes
         mi=mmax-(i-1)/mindex
         do j=1,nmodes
            mj=mmax-(j-1)/mindex
            autemp=0.
            iutemp=0.
            do iterp=1,nxinterp
               if (.not. nowindow) then
                  if ( (mi.le.mtop(iterp)) .and. (mi.ge.mlo(iterp)) &
                       .and. (mj.le.mtop(iterp)) .and. (mj.ge.mlo(iterp))) then
                     autemp(:,iterp)=auarr(:,mindex*(mtop(iterp)-mi)+1, &
                          mindex*(mtop(iterp)-mj)+1,iterp)
                     iutemp(:,iterp)=iuarr(7:9,mindex*(mtop(iterp)-mi)+1, &
                          mindex*(mtop(iterp)-mj)+1,iterp)
                  else  !zeroing shouldn't be necessary
                     autemp(:,iterp)=0.
                     iutemp(:,iterp)=0.
                  endif
               else 
                  autemp(:,iterp)=auarr(:,mindex*(mtop(iterp)-mi)+1, &
                       mindex*(mtop(iterp)-mj)+1,iterp)
                  iutemp(:,iterp)=iuarr(7:9,mindex*(mtop(iterp)-mi)+1, &
                          mindex*(mtop(iterp)-mj)+1,iterp)
               endif
            enddo
            do k=1,19 ! index of auarr elements
               call spline(-psigrid,autemp(k,:),nxinterp, &
                    bnd_set,bnd_set,auspline(k,i,j,:))
            enddo
            do k=7,9 ! index of iuarr arguments
               call spline(-psigrid,iutemp(k,:),nxinterp, &
                    bnd_set,bnd_set,iuspline(k,i,j,:))
            enddo

         enddo
      enddo

      return

end subroutine splinearr

