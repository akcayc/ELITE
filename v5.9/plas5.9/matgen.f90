      
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
           mindex,mmax,m0,auarr,gamsq,gamre, &
           iuarr,nowindow,nmwinhalf,mmin,mwindow,maxdxinterp,xinterp, &
           comp_ai,comp_api,comp_appi,compression,rotnorm,rotation, &
           nn,gamim,ruarr
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
      complex ai,api,appi
!      real aic,apic,appic
      integer j,kk,itinit,m,mp
      integer jfull1,jfull2,kkfull1,kkfull2,offset1,offset2
      real dmnq,dmpnq
      real weight
      real dmpnq2,dmpnqm,dmpnqm2,dmpnq2m,dmnq2  !pre-multiply for efficiency
      integer mtop1,mtop2 ! top index for arrays at iterp and iterp+1
      real xi
!      real Gamsqr,Gamr,Gamsqim,Gamimg,Gamsqdiff,
      real rot_loc 
      complex Gamma_c,Gamma2_c,Gam2diff_c
!      real invGamr,invGamim
!      real ari(nmodeswin,nmodeswin),apri(nmodeswin),appri(nmodeswin)
!      real aii(nmodeswin,nmodeswin),apii(nmodeswin),appii(nmodeswin)
      !local values of Gamma for rotation case
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

! calculate modified local Gamma values for rotational case
       if (rotation) then
          rot_loc=alph*rotnorm(iterp+1)+(1.-alph)*rotnorm(iterp)
!          Gamsqr=gamsq-(gamim+nn*rot_loc)**2
!          Gamr=sqrt(gamsq)
!          Gamsqim=2*Gamr*(gamim+nn*rot_loc)
!          Gamimg=gamim+nn*rot_loc
!          Gamsqdiff=Gamsqr-gamsq
          Gamma_c=cmplx(gamre,gamim+nn*rot_loc)
          Gamma2_c=Gamma_c*Gamma_c
          Gam2diff_c=Gamma2_c-gamre**2
!          invGamr=Gamr/(gamsq+Gamimg**2)
!          invGamim=-Gamimg/(gamsq+Gamimg**2)
       else
!          Gamsqr=gamsq
!          Gamr=sqrt(gamsq)
          Gamma_c=sqrt(gamsq)
          Gamma2_c=gamsq
       endif


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
                 Gamma2_c*(iuarr(9,jfull2,kkfull2,iterp+1)+ &
                 afact*iuarr(9,jfull1,kkfull1,iterp) &
                 +(iuarr(8,jfull2,kkfull2,iterp+1)+ &
                 afact*iuarr(8,jfull1,kkfull1,iterp))*dmnq &
                 +(iuarr(7,jfull2,kkfull2,iterp+1)+ &
                 afact*iuarr(7,jfull1,kkfull1,iterp))*dmpnq) 


! add imaginary gamsq and gam terms for rotation 
              if (rotation) then
                 ! first add additional real Gamsq part to terms
                 !   which had gamsq added in shoot
!                 ai=ai+ &
!  bug fix attempt 3/03, should be minus sign
                 ai=ai- &
                      Gam2diff_c*(iuarr(3,jfull2,kkfull2,iterp+1)+ &
                      afact*iuarr(3,jfull1,kkfull1,iterp) &
                      +(iuarr(2,jfull2,kkfull2,iterp+1)+ &
                      afact*iuarr(2,jfull1,kkfull1,iterp))*dmnq &
                      +(iuarr(1,jfull2,kkfull2,iterp+1)+ &
                      afact*iuarr(1,jfull1,kkfull1,iterp))*dmpnq)
!                 api=api+ &
                 api=api- &
                      Gam2diff_c*(iuarr(6,jfull2,kkfull2,iterp+1)+ &
                      afact*iuarr(6,jfull1,kkfull1,iterp) &
                      +(iuarr(5,jfull2,kkfull2,iterp+1)+ &
                      afact*iuarr(5,jfull1,kkfull1,iterp))*dmnq &
                      +(iuarr(4,jfull2,kkfull2,iterp+1)+ &
                      afact*iuarr(4,jfull1,kkfull1,iterp))*dmpnq) 

! add real parts of ai and api due to Omega' terms
                 ai=ai+ &
                      Gamma_c*(ruarr(1,jfull2,kkfull2,iterp+1)+ &
                      afact*ruarr(1,jfull1,kkfull1,iterp))
                 api=api+ &
                      Gamma_c*(ruarr(2,jfull2,kkfull2,iterp+1)+ &
                      afact*ruarr(2,jfull1,kkfull1,iterp))

! 2/06 add parts of ai and api due to Omega*Omega' terms
                 ai=ai+ &
                      (ruarr(3,jfull2,kkfull2,iterp+1)+ &
                      afact*ruarr(3,jfull1,kkfull1,iterp) )
                 api=api+ &
                      (ruarr(4,jfull2,kkfull2,iterp+1)+ &
                      afact*ruarr(4,jfull1,kkfull1,iterp) )

! pbs 3/03 add real parts of ai and api due to Omega'/Gamma terms
!   remove 3/18/03
!                 ai=ai+ &
!                      invGamr*( (ruarr(4,jfull2,kkfull2,iterp+1)+ &
!                      afact*ruarr(4,jfull1,kkfull1,iterp))*dmpnqm &
!                      +(ruarr(5,jfull2,kkfull2,iterp+1)+ &
!                      afact*ruarr(5,jfull1,kkfull1,iterp))*dmpnq )
!                 api=api+ &
!                      invGamr*(ruarr(3,jfull2,kkfull2,iterp+1)+ &
!                      afact*ruarr(3,jfull2,kkfull1,iterp))*dmpnqm

              endif

              if (compression) then
                 ai=ai+comp_ai(j,kk)
                 api=api+comp_api(j,kk)
                 appi=appi+comp_appi(j,kk)
              endif

!  TEMPORARY TEST, see if gamma_c matters (can we divide it out?)
!              ai=Gamma_c*ai
!              api=Gamma_c*api
!              appi=Gamma_c*appi

             smat(j,kk)=   alph*(d1*api+dd1*appi)
             qmat(j,kk)=alph*(ai+d2*api+dd2*appi)
             pmat(j,kk)=   alph*(d3*api+dd3*appi)

          end do
       end do
    return

end subroutine matgen






