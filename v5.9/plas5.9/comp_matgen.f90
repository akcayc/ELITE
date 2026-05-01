      
subroutine comp_matgen(i,iterp,nmodeswin,mwintop)

! 8/02 snyder
!   calculates comp_ai,comp_api and comp_appi matrices which
!    contain the compressional and rotational contributions
!    to the P Q and S matrices calculated by matgen

      use elite_data, only: xx,psixx,psigrid,nxinterp,nx, &
           mindex,mmax,m0,gamsq,vuarr, &
           nowindow,nmwinhalf,mmin,mwindow,maxdxinterp,xinterp, &
           comp_ai,comp_api,comp_appi,rotnorm,rotation,nn,gamim,gamre
      implicit none
      integer i,iterp
      integer nmodeswin,mwintop
      real psival
      real alph,afact
      real x1,x2,x3
      integer j,kk,itinit,m,mp
      integer jfull1,jfull2,kkfull1,kkfull2,offset1,offset2
      real dmnq,dmpnq
      real weight
      real dmpnq2,dmpnqm,dmpnqm2,dmpnq2m,dmnq2  !pre-multiply for efficiency
      integer mtop1,mtop2 ! top index for arrays at iterp and iterp+1
      real xi
      complex, dimension(nmodeswin,nmodeswin) :: dd,av,bv,evp,fvp,gvp,hvp, &
           umv,umvp,ddinv,ddinvav,ddinvbv,ddtemp,ddtest,ddinvtemp,avtemp
      integer ipiv(nmodeswin),info
      real rndoff
!      real Gamsqr,Gamr,Gamsqim,Gamimg,Gamsqdiff,
      real rot_loc
!      real invGamr,invGamim
      complex Gamma_c,Gamma2_c
!      real ddri,ddii,avri,avii,bvri,bvii
!      real evpri,evpii,fvpri,fvpii,gvpri,gvpii,hvpri,hvpii
!      real umvri,umvii
      rndoff=1.e-10
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
       psival=psixx(i)
       itinit=iterp
       if(psival < psigrid(iterp+1)) then
          do while(psival-psigrid(iterp) < 0.)
             iterp=iterp+1
             if(iterp.gt.nxinterp) go to 100
          end do
          iterp=iterp-1
          else if(psival > psigrid(iterp)) then
             do while(psival-psigrid(iterp) > 0.)
             iterp=iterp-1
             if(iterp.le.0) go to 100
          end do
       endif
 100   if(iterp.lt.1.or.iterp.gt.(nxinterp-1)) then
          write(6,*) 'problem in comp_matgen'
          write(6,*) 'iterp out of range; iterp=',iterp
          write(6,*) 'i,xx(i)=,psixx(i)=',i,xx(i),psixx(i)
          write(6,*) 'psigrid(1),psigrid(nxinterp)', &
              psigrid(1),psigrid(nxinterp)
          stop
       endif
       alph=(psival-psigrid(iterp))/(psigrid(iterp+1)-psigrid(iterp))
       if(alph.lt.0..or.alph.gt.1.0) then
          write(6,*) 'problem in comp_matgen'
          write(6,*) 'alph out of range; alph=',alph
          write(6,*) 'i,iterp,itinit,xx(i)=',i,iterp,itinit,xx(i)
          write(6,*) 'psigrid(iterp),psigrid(iterp+1)', &
              psigrid(iterp),psigrid(iterp+1)
          write(6,*) 'psigrid(iterp-1),psigrid(iterp+2)', &
              psigrid(iterp-1),psigrid(iterp+2)
          stop
       endif

!-----------------------------------------------------------------------
!     3.0 do the interpolation to complete construction of comp_ai etc.
!-----------------------------------------------------------------------

      mtop1=mmax
      mtop2=mmax
      if(.not.nowindow) then
         ! note that mtop is one above usual mup.  need to leave room
         !  for offset in comp_matgen
         mtop1=min0(m0+nmwinhalf-int(xinterp(iterp)+0.5)+maxdxinterp,mmax)
         if (mtop1 < (mmin+mwindow/mindex-1+maxdxinterp)) &
              mtop1=mmin+mwindow/mindex-1+maxdxinterp
         mtop1=min0(mtop1,mmax)
         mtop2=min0(m0+nmwinhalf-int(xinterp(iterp+1)+0.5)+maxdxinterp,mmax)
         if (mtop2 < (mmin+mwindow/mindex-1+maxdxinterp) ) &
              mtop2=mmin+mwindow/mindex-1+maxdxinterp
         mtop2=min0(mtop2,mmax)
      endif

      offset1=mindex*(mtop1-mwintop)
      offset2=mindex*(mtop2-mwintop)
      if ( (offset1 < 0) .or. (offset2 < 0) .or. &
           (offset1 > 2*maxdxinterp*mindex) .or. &
           (offset2 > 2*maxdxinterp*mindex) ) then
         write(*,*) 'problem in comp_matgen, offset1=', &
              offset1,'offset2=',offset2,&
              'mtop1=',mtop1,'mtop2=',mtop2,'mmax=',mmax,'mmin=',mmin, &
              'mwintop=',mwintop,'iterp=',iterp
         stop
      endif
!  6/01 for efficiency use afact=(1-alph)/alph for 1-alph terms
!     and multiply sums by alph at end
       afact=(1.-alph)/alph

       if (rotation) then
          rot_loc=alph*rotnorm(iterp+1)+(1.-alph)*rotnorm(iterp)
!          Gamsqr=gamsq-(gamim+nn*rot_loc)**2
!          Gamr=sqrt(gamsq)
!          Gamsqim=2*Gamr*(gamim+nn*rot_loc)
!          Gamimg=gamim+nn*rot_loc
!          Gamsqdiff=Gamsqr-gamsq
!          Gamma_c=cmplx(sqrt(gamsq),gamim+nn*rot_loc)
          Gamma_c=cmplx(gamre,gamim+nn*rot_loc)
          Gamma2_c=Gamma_c*Gamma_c
!          invGamr=Gamr/(gamsq+Gamimg**2)
!          invGamim=-Gamimg/(gamsq+Gamimg**2)
       else
!          Gamsqr=gamsq
!          Gamr=sqrt(gamsq)
          Gamma_c=sqrt(gamsq)
          Gamma2_c=gamsq
       endif


       do  j=1,nmodeswin          ! j labels rows of matrices
          jfull1=j+offset1
          jfull2=j+offset2
          mp=mtop1-(jfull1-1)/mindex
          dmpnq=xx(i)+mp-m0
          dmpnq2=dmpnq*dmpnq
          do  kk=1,nmodeswin          ! kk labels columns of matrices
             kkfull1=kk+offset1
             kkfull2=kk+offset2
             m=mtop1-(kkfull1-1)/mindex !if updownsym mindex=1
             dmnq=xx(i)+m-m0
             dmnq2=dmnq*dmnq
             dmpnqm=dmpnq*dmnq
             dmpnqm2=dmpnq*dmnq2
             dmpnq2m=dmpnq2*dmnq

!  dd is D, the denominator in the vm and vmp equations
             dd(j,kk)=(vuarr(1,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(1,jfull1,kkfull1,iterp))*dmpnqm &
                  +Gamma2_c*(vuarr(2,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(2,jfull1,kkfull1,iterp))

! av is the matrix multiplying u_m in the v_m equation
             av(j,kk)=Gamma2_c*(vuarr(4,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(4,jfull1,kkfull1,iterp))
! bv is the matrix multiplying ump in the v_m equation
             bv(j,kk)=Gamma2_c*(vuarr(3,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(3,jfull1,kkfull1,iterp))
             

! evp is matrix multiplying vm on the rhs of the vmp equation
!             evp(j,kk)=(vuarr(5,jfull2,kkfull2,iterp+1)+ &
!                  afact*vuarr(5,jfull1,kkfull1,iterp))* &
!                  (dmpnq+dmnq)- &
!                  sqrt(gamsq)*(vuarr(6,jfull2,kkfull2,iterp+1)+ &
!                  afact*vuarr(6,jfull1,kkfull1,iterp))
!!!!TEMP DEBUG!!!!
             evp(j,kk)=(vuarr(5,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(5,jfull1,kkfull1,iterp))* &
!                  (2.*m*mp/(m0-xx(i))-mp-m)- &
                  (dmpnq+dmnq)- &
                  Gamma_c*(vuarr(6,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(6,jfull1,kkfull1,iterp))


! fvp multiplies um in the vmp equation
             fvp(j,kk)=Gamma_c* &
                  (vuarr(8,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(8,jfull1,kkfull1,iterp))
! gvp multiplies ump in the vmp equation
             gvp(j,kk)=Gamma_c* &
                  (vuarr(7,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(7,jfull1,kkfull1,iterp)) &
                  +Gamma2_c*(vuarr(12,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(12,jfull1,kkfull1,iterp))+ &
                  Gamma2_c*(vuarr(13,jfull2,kkfull2,iterp+1)+ &
                  vuarr(13,jfull1,kkfull1,iterp))

!             gvp(j,kk)=Gamr* &
!                  (vuarr(7,jfull2,kkfull2,iterp+1)+ &
!                  afact*vuarr(7,jfull1,kkfull1,iterp)) &
!                  +(vuarr(12,jfull2,kkfull2,iterp+1)+ &
!                  afact*vuarr(12,jfull1,kkfull1,iterp))+ &
!                  gamsq*(vuarr(13,jfull2,kkfull2,iterp+1)+ &
!                  vuarr(13,jfull1,kkfull1,iterp))
!             gvp(j,kk)=Gamr* &
!                  (vuarr(7,jfull2,kkfull2,iterp+1)+ &
!                  afact*vuarr(7,jfull1,kkfull1,iterp)) &
!                  +av(j,kk)
!                  (vuarr(4,jfull2,kkfull2,iterp+1)+ &
!                  afact*vuarr(4,jfull1,kkfull1,iterp))
! hvp multiplies umpp in the vmp equation
!             hvp(j,kk)=gamsq*(vuarr(3,jfull2,kkfull2,iterp+1)+ &
!                  afact*vuarr(3,jfull1,kkfull1,iterp))
             hvp(j,kk)=Gamma2_c*(vuarr(11,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(11,jfull1,kkfull1,iterp))

! note: this equals bv, change later for efficiency

! umv multiplies v_m in the um equation
             umv(j,kk)=vuarr(9,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(9,jfull1,kkfull1,iterp) 
!                  +invGamr*(vuarr(14,jfull2,kkfull2,iterp+1)+ &
!                  afact*vuarr(14,jfull1,kkfull1,iterp))
! umvp multiplies vmp in the um equation
             umvp(j,kk)=vuarr(10,jfull2,kkfull2,iterp+1)+ &
                  afact*vuarr(10,jfull1,kkfull1,iterp)

          end do
       end do

! invert dd
       ddinv=0.
       do j=1,nmodeswin
          ddinv(j,j)=1.
       enddo
!       ddinv=dd**-1
       ddtemp=dd
       call zgesv(nmodeswin,nmodeswin,ddtemp,nmodeswin,ipiv,ddinv, &
            nmodeswin,info)
!       call cgesv2(nmodeswin,nmodeswin,ddtemp,nmodeswin,ipiv,ddinv, &
!            nmodeswin,info)
       if(info.ne.0) then
          write(6,*) "cgesv in comp_matgen failed; info=",info
          write(*,*) 'i=',i,'iterp=',iterp,'xx(i)=',xx(i)
          stop
       endif


       ddinvav=matmul(ddinv,av)
       ddinvbv=matmul(ddinv,bv)

       comp_ai=matmul(umv,ddinvav) + &
            matmul(umvp,matmul(ddinv,(matmul(evp,ddinvav)+fvp)))
       comp_api=matmul(umv,ddinvbv) + &
            matmul(umvp,matmul(ddinv,(matmul(evp,ddinvbv)+gvp)))
       comp_appi=matmul(umvp,matmul(ddinv,hvp))

    return

end subroutine comp_matgen






