      
subroutine matgen(i,iterp,nmodeswin,mwintop)

!-----------------------------------------------------------------------
!     Calculates P, Q and S matrices at each x-point
!     smat*U[x-dx]+qmat*U[x]+pmat*U[x+dx]
!     Same version works for updown symmetric or antisymmetric
!     modified 5/19/98 to allow windowing
!     nmodeswin is number of modes carried in current m window
!     mwintop is maximum m visible in current m window
!-----------------------------------------------------------------------
      use elite_data, only: xx,psixx,psigrid,pmat,qmat,smat,nxinterp,nx, &
           mindex,mmax,m0,auarr,gamsq, &
           iuarr,nowindow,nmwinhalf,mmin,mwindow,maxdxinterp,xinterp
      implicit none
      integer i,iterp
      integer nmodeswin,mwintop
!      real xval
      real psival
!      real ain,apin,appin
      real alph,afact
      real x1,x2,x3
      real d1,d2,d3
      real dd1,dd2,dd3
      real denom1,denom2,denom3
      real ai,api,appi
      real aic,apic,appic
      integer j,kk,itinit,m,mp
      integer jfull1,jfull2,kkfull1,kkfull2,offset1,offset2
      real dmnq,dmpnq
      real weight
      real dmpnq2,dmpnqm,dmpnqm2,dmpnq2m,dmnq2  !pre-multiply for efficiency
      integer mtop1,mtop2 ! top index for arrays at iterp and iterp+1
      real xi
!-----------------------------------------------------------------------
!  i labels x-point used
!  iterp labels the sparser interpolation points
!  check Mercier stability of each radius:
! how do we do this now?
!        if (root.lt.error) then
!          write(6,*)'ERROR***surface i=',i,' is Mercier unstable'
!          stop
!        end if
!   x-mesh spacing...
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     1.0 find the two xinterp values which xx(i) falls between
!-----------------------------------------------------------------------
! 9/01 new version must now do all interpolation using the monotonic
!   psixx grid rather than the potentially non-monotonix xx grid
!   also change signs because psi decreases moving inward

       if(iterp.le.0) then
          iterp=1
       endif
!       xval=xx(i)
       psival=psixx(i)
       itinit=iterp
!       if(xval.gt.xinterp(iterp+1)) then
       if(psival < psigrid(iterp+1)) then
!          do while(xval-xinterp(iterp).gt.0.)
          do while(psival-psigrid(iterp) < 0.)
             iterp=iterp+1
             if(iterp.gt.nxinterp) go to 100
          end do
          iterp=iterp-1
!       else if(xval.lt.xinterp(iterp)) then
          else if(psival > psigrid(iterp)) then
!          do while(xval-xinterp(iterp).lt.0.)
             do while(psival-psigrid(iterp) > 0.)
             iterp=iterp-1
             if(iterp.le.0) go to 100
          end do
       endif
 100   if(iterp.lt.1.or.iterp.gt.(nxinterp-1)) then
          write(6,*) 'problem in matgen'
          write(6,*) 'iterp out of range; iterp=',iterp
          write(6,*) 'i,xx(i)=,psixx(i)=',i,xx(i),psixx(i)
          write(6,*) 'psigrid(1),psigrid(nxinterp)', &
              psigrid(1),psigrid(nxinterp)
          stop
       endif
       alph=(psival-psigrid(iterp))/(psigrid(iterp+1)-psigrid(iterp))
       if(alph.lt.0..or.alph.gt.1.0) then
          write(6,*) 'problem in matgen'
          write(6,*) 'alph out of range; alph=',alph
          write(6,*) 'i,iterp,itinit,xx(i)=',i,iterp,itinit,xx(i)
          write(6,*) 'psigrid(iterp),psigrid(iterp+1)', &
              psigrid(iterp),psigrid(iterp+1)
          write(6,*) 'psigrid(iterp-1),psigrid(iterp+2)', &
              psigrid(iterp-1),psigrid(iterp+2)
          stop
       endif
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
!-----------------------------------------------------------------------
!     3.0 do the interpolation to complete construction of s,q, and p
!-----------------------------------------------------------------------

      mtop1=mmax
      mtop2=mmax
      if(.not.nowindow) then
         ! note that mtop is one above usual mup.  need to leave room
         !  for offset in matgen
         mtop1=min0(m0+nmwinhalf-int(xinterp(iterp)+0.5)+maxdxinterp,mmax)
         if (mtop1 < (mmin+mwindow/mindex-1+maxdxinterp)) &
              mtop1=mmin+mwindow/mindex-1+maxdxinterp
         mtop1=min0(mtop1,mmax)
         mtop2=min0(m0+nmwinhalf-int(xinterp(iterp+1)+0.5)+maxdxinterp,mmax)
         if (mtop2 < (mmin+mwindow/mindex-1+maxdxinterp) ) &
              mtop2=mmin+mwindow/mindex-1+maxdxinterp
         mtop2=min0(mtop2,mmax)
      endif

!       offset=mindex*(mmax-mwintop)
      offset1=mindex*(mtop1-mwintop)
      offset2=mindex*(mtop2-mwintop)
      if ( (offset1 < 0) .or. (offset2 < 0) .or. &
           (offset1 > 2*maxdxinterp*mindex) .or. &
           (offset2 > 2*maxdxinterp*mindex) ) then
         write(*,*) 'problem in matgen, offset1=',offset1,'offset2=',offset2,&
              'mtop1=',mtop1,'mtop2=',mtop2,'mmax=',mmax,'mmin=',mmin, &
              'mwintop=',mwintop,'iterp=',iterp
         stop
      endif
!  6/01 for efficiency use afact=(1-alph)/alph for 1-alph terms
!     and multiply sums by alph at end
       afact=(1.-alph)/alph

       do  j=1,nmodeswin          ! j labels rows of matrices
!          jfull=j+offset
          jfull1=j+offset1
          jfull2=j+offset2
!          mp=mmax-(jfull-1)/mindex  !if updownsym mindex=1
          mp=mtop1-(jfull1-1)/mindex
          dmpnq=xx(i)+mp-m0
          dmpnq2=dmpnq*dmpnq
!          write(6,*) "mp=",mp,"dmpnq=",dmpnq
          do  kk=1,nmodeswin          ! kk labels columns of matrices
!             kkfull=kk+offset
             kkfull1=kk+offset1
             kkfull2=kk+offset2
!             m=mmax-(kkfull-1)/mindex !if updownsym mindex=1
             m=mtop1-(kkfull1-1)/mindex !if updownsym mindex=1
             dmnq=xx(i)+m-m0
             dmnq2=dmnq*dmnq
             dmpnqm=dmpnq*dmnq
             dmpnqm2=dmpnq*dmnq2
             dmpnq2m=dmpnq2*dmnq
!             write(6,*) "m=",m,"dmnq=",dmnq
             ai   =(auarr(1,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(1,jfull1,kkfull1,iterp))*dmpnqm &
                 +(auarr(2,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(2,jfull1,kkfull1,iterp))*dmpnqm2 &
                 +(auarr(3,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(3,jfull1,kkfull1,iterp))*dmpnq2m &
                 +(auarr(4,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(4,jfull1,kkfull1,iterp))*dmnq2 &
                 +(auarr(5,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(5,jfull1,kkfull1,iterp))*dmpnq2 &
                 +(auarr(6,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(6,jfull1,kkfull1,iterp))*dmnq &
                 +(auarr(7,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(7,jfull1,kkfull1,iterp))*dmpnq &
                 +auarr(8,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(8,jfull1,kkfull1,iterp)
! add in inertia terms...
!   7/01 these terms now added to au.. in solveit for efficiency
!                 gamsq*(iund(jfull2,kkfull2,iterp+1)+ &
!                 afact*iund(jfull1,kkfull1,iterp) &
!                 +(ium(jfull2,kkfull2,iterp+1)+ &
!                 afact*ium(jfull1,kkfull1,iterp))*dmnq &
!                 +(iuk(jfull2,kkfull2,iterp+1)+ &
!                 afact*iuk(jfull1,kkfull1,iterp))*dmpnq) 
             api   =(auarr(9,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(9,jfull1,kkfull1,iterp))*dmpnqm &
                 +(auarr(10,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(10,jfull1,kkfull1,iterp))*dmpnqm2 &
                 +(auarr(11,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(11,jfull1,kkfull1,iterp))*dmpnq2m &
                 +(auarr(12,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(12,jfull1,kkfull1,iterp))*dmnq2 &
                 +(auarr(13,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(13,jfull1,kkfull1,iterp))*dmpnq2 &
                 +(auarr(14,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(14,jfull1,kkfull1,iterp))*dmnq &
                 +(auarr(15,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(15,jfull1,kkfull1,iterp))*dmpnq &
                 +auarr(16,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(16,jfull1,kkfull1,iterp)
! add in inertia terms...
!   7/01 these terms now added to au.. in solveit for efficiency
!                 gamsq*(iupnd(jfull2,kkfull2,iterp+1)+ &
!                 afact*iupnd(jfull1,kkfull1,iterp) &
!                 +(iupm(jfull2,kkfull2,iterp+1)+ &
!                 afact*iupm(jfull1,kkfull1,iterp))*dmnq &
!                 +(iupk(jfull2,kkfull2,iterp+1)+ &
!                 afact*iupk(jfull1,kkfull1,iterp))*dmpnq) 
              appi =(auarr(17,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(17,jfull1,kkfull1,iterp))*dmpnqm &
                 +(auarr(18,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(18,jfull1,kkfull1,iterp))*dmpnqm2 &
                 +(auarr(19,jfull2,kkfull2,iterp+1)+ &
                 afact*auarr(19,jfull1,kkfull1,iterp))*dmpnq2m - &
! add in inertia terms...
                 gamsq*(iuarr(9,jfull2,kkfull2,iterp+1)+ &
                 afact*iuarr(9,jfull1,kkfull1,iterp) &
                 +(iuarr(8,jfull2,kkfull2,iterp+1)+ &
                 afact*iuarr(8,jfull1,kkfull1,iterp))*dmnq &
                 +(iuarr(7,jfull2,kkfull2,iterp+1)+ &
                 afact*iuarr(7,jfull1,kkfull1,iterp))*dmpnq) 
             smat(j,kk)=   alph*(d1*api+dd1*appi)
             qmat(j,kk)=alph*(ai+d2*api+dd2*appi)
             pmat(j,kk)=   alph*(d3*api+dd3*appi)


!!!!!!!!!!!!TEST SPLINES 1/02
!!!!!! assume no windowing for this test nmodeswin=nmodes
!             call splint(-psigrid,auarr(1,jfull1,kkfull1,:), &
!                  auspline(1,j,kk,:), &
!                  nxinterp,-psixx(i),xi)

!             write(*,*) 'i=',i,'j=',j,'kk=',kk,'linear=', &
!                  alph*auarr(1,jfull2,kkfull2,iterp+1)+ &
!                  (1.-alph)*auarr(1,jfull1,kkfull1,iterp), &
!                  'spline=',xi

!             if(j.eq.kk) then
!                qmat(j,j)=qmat(j,j)+alpha*dmercier* &
!                    (-2.*pi*qref/aspect)
!!!!!!!!!!!Need to look at this again!!!!!!!!!!2/10/00
!             qmat(j,j)=qmat(j,j)+alpha_edge*dmercier* &
!                    (-2.*pi*qref/aspect)   

!-----------------------------------------------------------------------
! following coding is to put inertia at rational surfaces
!  using bug(1)
!-----------------------------------------------------------------------
!                if(m.lt.mres) then
!             if(bug(1) > 0.) then
!                   weight=(alph*auarr(17,jfull,kkfull,iterp+1)+ &
!                       (1.-alph)*auarr(17,jfull,kkfull,iterp))
!                   smat(j,j)=smat(j,j)+weight*bug(1)*dd1
!                   qmat(j,j)=qmat(j,j)+weight*bug(1)*dd2
!                   pmat(j,j)=pmat(j,j)+weight*bug(1)*dd3
!                endif
!             endif
          end do
       end do
    return

end subroutine matgen






