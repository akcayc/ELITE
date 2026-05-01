
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ELITE (Edge Localized Instabilities in Tokamak Experiments)
!
! Authors:  P.B. Snyder, H.R. Wilson, R.L. Miller
!
! Copyright(c) 1999-2005, All rights reserved
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

program elitevac

!-----------------------------------------------------------------------
!  Updated for new version of code and converted to f90
!    Snyder 1/26/00
!     make former main routine into subroutine dwcalc
 
    use elitevac_data, only: init,surfdat,verbose

!    if (verbose .ge. 5) write(6,*) 'call init'
!    write(*,*) 'verbose= ',verbose
    call init   ! set defaults and read in the runname.in namelist
    if (verbose .ge. 5) write(*,*) 'call dwcalc'
!    write(*,*) 'verbose= ',verbose
    call dwcalc   ! main routine to calculate vacuum dW

!  stop


contains

  subroutine dwcalc

!-----------------------------------------------------------------------
! vacuum code for high mode number perturbations
! based upon A. Pletzer's C++ code
! 5/8/98
!-----------------------------------------------------------------------
      use elitevac_data, only: npts,mmax,mmin,ntor,mref,ng,theta,pi,r,  &
           fsurf,qsurf,dt,runname,lrunname,cosi,sini,zti,verbose
 !     use elitevac_data
      implicit none
      real, dimension(npts,npts) :: re_rgijr,im_rgijr,re_rkijr,im_rkijr, &
           gana,frhs
      real, dimension(npts*2,npts*2) :: g,kmat
      real, dimension(npts*2,mmax-mmin+1) :: bn,rhs,an
      real, dimension(mmax-mmin+1,mmax-mmin+1) :: dw,dwim
      real, dimension(mmax-mmin+1) :: coeff
      real rmajor,rminor
      real frgn,frgsing,frgn1
      real arg1,arg2,arg3,arg4
      real hyp
!      real cosi,sini
      real cosmw,sinmw
      real fac
      real rv,bpv,rt,zt,zv
      real tv,ave,diff,maxdiff,avediff,norm
      real asym0,asym1,asymlg0,asymlg1,lgasym0,lgasym1,sumsym0,sumsym1
      real sanorm0,sanorm1,rndoff,maxdiag
      integer i,j,jm,mj,ki
      integer ib,jb
      integer nrows,ncolumns
      integer, dimension(npts*2) :: ipiv
      integer info
      integer nmodes,k
      integer mval
      integer nunit,ios
      integer npts2
   
      re_rgijr=0.; im_rgijr=0.; re_rkijr=0.; im_rkijr=0.
      gana=0.; frhs=0.; g=0.; kmat=0.

!     write(*,*) 'in dwcalc, verbose=',verbose,' ng=',ng

!-----------------------------------------------------------------------
! debug stuff: check all of the functions
!-----------------------------------------------------------------------
!$$$      do i=1,npts
!$$$         tv=theta(i)
!$$$         write(16,*) rv(tv),zv(tv),rt(tv),zt(tv),cosmw(1,tv),
!$$$     $        sinmw(1,tv),bpv(tv)
!$$$      end do
!$$$      stop
!-----------------------------------------------------------------------
! integrate R*G_n and R*K_n with analytic singularity subtracted out
! to smooth integrands
!-----------------------------------------------------------------------
      if (verbose .ge. 5) write(*,*) 'call integrate'
      call integrate(re_rgijr,im_rgijr,re_rkijr,im_rkijr)

!      write(*,*) 're_rgijr',re_rgijr
!      open(unit=72,file='phil.debug',status='unknown')
!      write(72,*) 're_rkijr', re_rkijr
!-----------------------------------------------------------------------
! integrate analytic singularity piece
!-----------------------------------------------------------------------
      if (verbose .ge. 5) write(*,*) 'call int_ana'
      call int_ana(gana)
!      write(72,*) 'gana',gana
!-----------------------------------------------------------------------
! add back in the analytic singularity piece that was subtracted out
!-----------------------------------------------------------------------
      do j=1,npts
         do i=1,npts
!            cosi=cosmw(mref,theta(i))
!            sini=sinmw(mref,theta(i))
            re_rgijr(i,j)=re_rgijr(i,j)-gana(i,j)*cosi(i)/(2.*pi)
            im_rgijr(i,j)=im_rgijr(i,j)-gana(i,j)*sini(i)/(2.*pi)
            re_rkijr(i,j)=-re_rkijr(i,j)- &
                cosi(i)*zti(i)*gana(i,j)/(4.*pi*r(i))
            im_rkijr(i,j)=-im_rkijr(i,j)- &
                sini(i)*zti(i)*gana(i,j)/(4.*pi*r(i))
         end do
!         re_rkijr(j,j)=cosmw(mref,theta(j))+re_rkijr(j,j)
         re_rkijr(j,j)=cosi(j)+re_rkijr(j,j)
!         im_rkijr(j,j)=sinmw(mref,theta(j))+im_rkijr(j,j)
         im_rkijr(j,j)=sini(j)+im_rkijr(j,j)
      end do
!-----------------------------------------------------------------------
! now construct the big matrices
!-----------------------------------------------------------------------
      do jb=1,npts
         do ib=1,npts
            j=(jb-1)*2+1
            i=(ib-1)*2+1
            g(i  ,j  )= re_rgijr(ib,jb)
            g(i+1,j  )= im_rgijr(ib,jb)
            g(i  ,j+1)=-im_rgijr(ib,jb)
            g(i+1,j+1)= re_rgijr(ib,jb)
            kmat(i  ,j  )= re_rkijr(ib,jb)
            kmat(i+1,j  )= im_rkijr(ib,jb)
            kmat(i  ,j+1)=-im_rkijr(ib,jb)
            kmat(i+1,j+1)= re_rkijr(ib,jb)
         end do
      end do
!-----------------------------------------------------------------------
! the "reduced" bfield
!-----------------------------------------------------------------------
      do mval=mmin,mmax
         j=mval+1-mmin
         do ib=1,npts
            i=2*(ib-1)+1
!-----------------------------------------------------------------------
! correcting error found by Wilson, adding
! sqrt(rt(theta(ib))**2+zt(theta(ib))**2) to fac 8/24/98
!-----------------------------------------------------------------------
            fac=sqrt(rt(theta(ib))**2+zti(ib)**2) &
                /(rv(theta(ib))**3*bpv(theta(ib)))
            bn(i,j)=fac*cosmw(mref-mval,theta(ib))
            bn(i+1,j)=-fac*sinmw(mref-mval,theta(ib))
         end do
      end do
!-----------------------------------------------------------------------
! next calculate rhs=-g*bn
!-----------------------------------------------------------------------
      nmodes=mmax+1-mmin
      do jm=1,nmodes
         do i=1,2*npts
            rhs(i,jm)=0.
            do j=1,2*npts
               rhs(i,jm)=rhs(i,jm)-g(i,j)*bn(j,jm)
            end do
         end do
      end do
!-----------------------------------------------------------------------
! solve k*an==rhs; rhs is replaced by an in routine albet
!-----------------------------------------------------------------------
      nrows=2*npts
      ncolumns=nmodes
!      call dgesv( nrows,ncolumns, kmat, kpts2, ipiv, rhs, kpts2, info )
      npts2=2*npts
!      write(72,*) 'beforerhs ',rhs,'kmat ',kmat
      call dgesv ( nrows,ncolumns, kmat, npts2, ipiv, rhs, npts2, info )
!      write(72,*) 'rhs', rhs, 'kmat', kmat
      if(info.ne.0) then
         write(6,*) "dgesv failed; info=",info
         stop
      endif
!-----------------------------------------------------------------------
! now rhs is the potential
! to get dW matrix we need to integrate bn*rhs
!-----------------------------------------------------------------------
      do j=1,nmodes
         mj=mmin-1+j
         coeff(j)=fsurf*(mj-ntor*qsurf)/qsurf
         do i=1,nmodes
            dw(i,j)=0.
            dwim(i,j)=0. ! this is for debug
         end do
      end do
      do k=1,npts
         rmajor=rv(theta(k))
!-----------------------------------------------------------------------
!         rminor=sqrt(rt(theta(k))**2+zt(theta(k))**2) <--wrong
!     correcting error found by Wilson 8/24/98
!-----------------------------------------------------------------------
         rminor=1.
         ki=2*(k-1)+1
         do i=1,nmodes
            do j=1,nmodes
               dw(i,j)=dw(i,j)+rmajor*rminor*(rhs(ki,i)*bn(ki,j) &
                   +rhs(ki+1,i)*bn(ki+1,j))
               dwim(i,j)=dwim(i,j)+rmajor*rminor*(rhs(ki,i)*bn(ki+1,j) &
                   -rhs(ki+1,i)*bn(ki,j))
            end do
         end do
      end do
!      open(unit=71,file='vac.debug',status='unknown')
!      write(71,*) ((dw(i,j),dwim(i,j),i=1,nmodes),j=1,nmodes)
!      close(71)
      do i=1,nmodes
         do j=1,nmodes
            dw(i,j)=-pi*coeff(i)*coeff(j)*dw(i,j)*dt
         end do
      end do
!-----------------------------------------------------------------------
! symmetrize dw
!-----------------------------------------------------------------------
! 12/20/00 reproduce Gato's method of evaluating the asymmetry of dW_vac
      if (verbose .ge. 5) write(*,*) 'symmetrize dW_vac'
      norm=0.
      maxdiag=0.
      do i=1,nmodes
         norm=norm+abs(dw(i,i))
         if (abs(dw(i,i)).gt.maxdiag) maxdiag=abs(dw(i,i))
      enddo
      norm=norm/nmodes 
      rndoff   = maxdiag*maxdiag/1.d6
      if (verbose .ge. 1) write(*,*) 'using rndoff of maxdiag^2/1.d6=',rndoff
      sumsym0  = 0.0
      sanorm0  = 0.0
      sumsym1  = 0.0
      sanorm1  = 0.0
      do j=1,nmodes  ! loop to calculate gato-like asymmetry 
         do i=j,nmodes
            if(i .eq. j) sanorm0 = 2.0* dw(i,j)*dw(i,j)
            if(i .ne. j) sanorm0 = (dw(i,j)*dw(i,j) + dw(j,i)*dw(j,i))
            if(abs(sanorm0) .le. rndoff) sanorm0 = rndoff
            sumsym0 = sumsym0 + 2.0*(dw(j,i)-dw(i,j))*(dw(j,i)-dw(i,j))/sanorm0

            if(i .eq. j) sanorm1 = sanorm1 +  dw(i,j)*dw(i,j)
            if(i .ne. j) sanorm1 = sanorm1+(dw(i,j)*dw(i,j) + dw(j,i)*dw(j,i))
            sumsym1  = sumsym1  +  2.0*(dw(j,i)-dw(i,j))*(dw(j,i)-dw(i,j))
         enddo
      enddo
  
      maxdiff=0.
      avediff=0.
      if (verbose .ge. 1) write(*,*) 'symmetrizing dW_vac'
      do i=1,nmodes
         do j=1,i
            ave=0.5*(dw(i,j)+dw(j,i))
            diff=abs(dw(i,j)-dw(j,i))/norm

            dw(i,j)=ave
            dw(j,i)=ave
!            write(*,*) 'relative diff,ave', diff,ave
            if (diff.gt.maxdiff) maxdiff=diff
            avediff=avediff+diff
         end do
      end do

      asym0   = sqrt(sumsym0) / nmodes
      asym1   = sqrt(sumsym1/sanorm1)


      avediff=avediff/((nmodes-1.)*nmodes/2.)
!      write(*,*) 'vacuum asymmetry maxdiff=',maxdiff,' avediff=',avediff, &
!           ' norm=',norm
      if (verbose .ge. 1) write(*,7020) asym0,asym1
 7020 format(/,1x,'Vacuum matrix asymmetry (first  norm) is: ',e16.6 &
           ,/,1x,'Vacuum matrix asymmetry (second norm) is: ',e16.6,/)
      if (asym0.gt.0.1 .and. asym1.gt.0.1) then
         write (*,*) '!!!!!WARNING, both asym0 and aym1 > 10%'
      else if (asym0.gt.0.01 .and. asym1.gt.0.01) then
         if (verbose .ge. 1) write (*,*) '!WARNING, both asym0 and asym1 exceed 1%'
!      else if (asym0.gt.0.01 .or. asym1.gt.0.01) then 
!         write (*,*) '!WARNING, one of the asymmetry norms exceeds 1%'
      endif
      write(*,*) 'outputting int(1000*asym1) for use in scripts'
      write(*,*) int(1000*asym1)
!-----------------------------------------------------------------------
! write dw matrix to a file so edge code can read it
!-----------------------------------------------------------------------
      nunit=31
      open(unit=nunit,file=runname(1:lrunname)//'.vac',status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem opening/creating ',runname(1:lrunname)//'.vac'
         stop
      endif
      write(nunit,*) nmodes,mmax,mmin 
! current version of plas stores maximum m value first so store
!   in reverse order below
      write(nunit,*) ((dw(nmodes+1-i,nmodes+1-j),i=1,nmodes),j=1,nmodes)
!-----------------------------------------------------------------------
! debug stuff
!-----------------------------------------------------------------------
!$$$      jb=1
!$$$      do ib=1,npts
!$$$         i=2*(ib-1)+1
!$$$         write(32,'(2e15.5)') rhs(i,jb),rhs(i+1,jb)
!$$$      end do
     
  end subroutine dwcalc


!-----------------------------------------------------------------------
  subroutine integrate(re_rgijr,im_rgijr,re_rkijr,im_rkijr)
!-----------------------------------------------------------------------
! integrate R*G_n(theta_i,theta_j) and R*K(theta_i,theta_j) from
!    theta_j-dt/2 to theta_j+dt/2
! re and im refer to real and imaginary parts which contain
!  cos(m*omega) and sin(m*omega) weighting factors respectively
!  trying a revised version with tent function finite elements
!-----------------------------------------------------------------------
      use elitevac_data, only: npts,theta,mref,z,r,pi,ng,xg,dt,wg,n05, &
           cosi,sini,zti,verbose
!      use functions
      implicit none
!      include 'i_vac.f'
      real re_rgijr(npts,npts),im_rgijr(npts,npts)
      real re_rkijr(npts,npts),im_rkijr(npts,npts)
      real frgn,frgsing,frgn1
      integer i,j,k
      real ti,z1,r1,r2,z2
      real tjk,vrgn,vrgsing,s,lam,temp
      real rdgdz,rdgdr,vrkn
      real rv,zv,rt,zt
      real cosjk,sinjk,rgijrp,rkijrp
      real cosmw,sinmw
      real xleft,xright,xgauss,xplus,xminus
      real costjk(npts,ng),sintjk(npts,ng),rvt(npts,ng),zvt(npts,ng)
      real ztt(npts,ng),rtt(npts,ng)
      integer jp,jm
! 1.0 zero out the arrays
      re_rgijr=0.; im_rgijr=0.; re_rkijr=0.; im_rkijr=0.  

! 1.1 cosine and sine arrays--move to init
!      do j=1,npts
         do i=1,npts
            ti=theta(i)
            cosi(i)=cosmw(mref,ti)
            sini(i)=sinmw(mref,ti)
!            rti(i)=rt(ti)
            zti(i)=zt(ti)
! 7/01 calculate and store additional arrays for use in big loop to improve
!     efficiency
            do k=1,ng/2
               tjk=ti+0.5*dt*(1.+xg(k))
               costjk(i,k)=cosmw(mref,tjk)
               sintjk(i,k)=sinmw(mref,tjk)
               rvt(i,k)=rv(tjk)
               zvt(i,k)=zv(tjk)
               rtt(i,k)=rt(tjk)
               ztt(i,k)=zt(tjk)
            enddo
            do k=ng/2+1,ng
               tjk=ti-0.5*dt*(1.-xg(k))
               costjk(i,k)=cosmw(mref,tjk)
               sintjk(i,k)=sinmw(mref,tjk)
               rvt(i,k)=rv(tjk)
               zvt(i,k)=zv(tjk)
               rtt(i,k)=rt(tjk)
               ztt(i,k)=zt(tjk)
            enddo
         end do
!      end do
! 2.0 integrate over tent functions
! integrals are done from -dt/2 to dt/2 for each theta(j) below
! 
      if (verbose .ge. 5) write(*,*) '2.0 integrate over tent functions'   
      do j=1,npts
         jp=j+1
         jm=j-1
         if(jp.gt.npts) jp=1
         if(jm.lt.1) jm=npts
         do i=1,npts
            ti=theta(i)
            z1=z(i)
            r1=r(i)
            if (r1.eq.0.) write(*,*) '0 radius',i,pi
            do k=1,ng/2
! 2.1 left hand side of integral--from -dt/2 to 0
!     y=(b+a)/2+(b-a)/2*xg where xg ranges from -1 to 1
!     take a=-dt and b=0
               xgauss=xg(ng/2+k)
               xleft=0.5+0.5*xgauss ! tent function value of tent_j
               xminus=1.-xleft   !tent function value of tent_j-1
!               xleft=1.         !debug
!               xminus=0.        !debug
               tjk=theta(j)-0.5d0*dt*(1.-xgauss)
               r2=rvt(j,ng/2+k)
               z2=zvt(j,ng/2+k)
               lam=1.+((r1-r2)**2+(z1-z2)**2)/(2.*r1*r2)
               s=lam/(sqrt(lam*lam-1.0))
               vrgn=frgn(ti,tjk,r1,z1,r2,z2,lam,s)
               vrgsing=frgsing(ti,tjk)
!               cosjk=cosmw(mref,tjk)
               cosjk=costjk(j,ng/2+k)
!               sinjk=sinmw(mref,tjk)
               sinjk=sintjk(j,ng/2+k)
               rgijrp=0.5d0*dt*wg(ng/2+k)
               re_rgijr(i,j)=re_rgijr(i,j)+rgijrp* &
                   (xleft*vrgn*cosjk-vrgsing*cosi(i))
               im_rgijr(i,j)=im_rgijr(i,j)+rgijrp* &
                   (xleft*vrgn*sinjk-vrgsing*sini(i))
               re_rgijr(i,jm)=re_rgijr(i,jm)+xminus*rgijrp* &
                   vrgn*cosjk
               im_rgijr(i,jm)=im_rgijr(i,jm)+xminus*rgijrp* &
                   vrgn*sinjk
!               r2=rv(tjk)
!               z2=zv(tjk)
!               lam=1.+((r1-r2)**2+(z1-z2)**2)/(2.*r1*r2)
!               s=lam/(sqrt(lam*lam-1.0))
               temp = (s*s-1.0) * n05 * ( s*vrgn/sqrt(s*s-1.0) &
                   - frgn1(ti,tjk,r1,z1,r2,z2,lam,s) )/( r1 * r2 )
               rdgdz = ( z1 - z2 ) * temp
               rdgdr = - vrgn/(2.0*r2) + ( lam*r1 - r2 ) * temp
               vrkn = ztt(j,ng/2+k) * rdgdr - rtt(j,ng/2+k) * rdgdz
               rkijrp= (dt/2) * wg(ng/2+k)
               re_rkijr(i,j) = re_rkijr(i,j) + rkijrp* &
                   (xleft*vrkn*cosjk+zti(i)* vrgsing/(2.*r1)*cosi(i))
               im_rkijr(i,j) = im_rkijr(i,j) + rkijrp* &
                   (xleft*vrkn*sinjk+zti(i)* vrgsing/(2.*r1)*sini(i))
               re_rkijr(i,jm) = re_rkijr(i,jm) + rkijrp* &
                   xminus*vrkn*cosjk
               im_rkijr(i,jm) = im_rkijr(i,jm) + rkijrp* &
                   xminus*vrkn*sinjk
!               if (re_rkijr(i,j).gt.1e10) write (*,*) 'point 1, i,j',i,j
               if (re_rkijr(i,j).gt.1e14) write (*,*) 'point 1, i,j',i,j
! 2.2 right hand side of integral--from 0 to dt/2
               xgauss=xg(k)
               xright=0.5-0.5*xgauss ! tent function value of tent_j
               xplus=1-xright ! tent function value of tent_j+1
!               xright=1.        ! debug
!               xplus=0.         ! debug
               tjk = theta(j) +0.5d0* dt*(1.+xgauss)
               r2=rvt(j,k)
               z2=zvt(j,k)
               lam=1.+((r1-r2)**2+(z1-z2)**2)/(2.*r1*r2)
               s=lam/(sqrt(lam*lam-1.0))
               vrgn=frgn(ti,tjk,r1,z1,r2,z2,lam,s)
               vrgsing=frgsing(ti,tjk)
!               cosjk=cosmw(mref,tjk)
               cosjk=costjk(j,k)
!               sinjk=sinmw(mref,tjk)
               sinjk=sintjk(j,k)
               rgijrp=0.5d0*dt*wg(k)
               re_rgijr(i,j)=re_rgijr(i,j)+rgijrp* &
                   (xright*vrgn*cosjk-vrgsing*cosi(i))
               im_rgijr(i,j)=im_rgijr(i,j)+rgijrp* &
                   (xright*vrgn*sinjk-vrgsing*sini(i))
               re_rgijr(i,jp)=re_rgijr(i,jp)+rgijrp* &
                   xplus*vrgn*cosjk
               im_rgijr(i,jp)=im_rgijr(i,jp)+rgijrp* &
                   xplus*vrgn*sinjk
!               r2=rv(tjk)
!               z2=zv(tjk)
!               lam=1.+((r1-r2)**2+(z1-z2)**2)/(2.*r1*r2)
!               s=lam/(sqrt(lam*lam-1.0))
               temp = (s*s-1.0) * n05 * ( s*vrgn/sqrt(s*s-1.0) &
                   - frgn1(ti,tjk,r1,z1,r2,z2,lam,s) )/( r1 * r2 )
               rdgdz = ( z1 - z2 ) * temp
               rdgdr = - vrgn/(2.0*r2) + ( lam*r1 - r2 ) * temp
               vrkn = ztt(j,k) * rdgdr - rtt(j,k) * rdgdz
               rkijrp=(dt/2) * wg(k)
               re_rkijr(i,j) = re_rkijr(i,j) + rkijrp* &
                   (xright*vrkn*cosjk+zti(i)* vrgsing/(2.*r1)*cosi(i))
               im_rkijr(i,j) = im_rkijr(i,j) + rkijrp* &
                   (xright*vrkn*sinjk+zti(i)* vrgsing/(2.*r1)*sini(i))
               re_rkijr(i,jp) = re_rkijr(i,jp) + rkijrp* &
                   xplus*vrkn*cosjk
               im_rkijr(i,jp) = im_rkijr(i,jp) + rkijrp* &
                   xplus*vrkn*sinjk
!               if (re_rkijr(i,j).gt.1e10) write (*,*) 'point 2, i,j',i,j
               if (re_rkijr(i,j).gt.1e14) write (*,*) 'point 2, i,j',i,j
            end do
         end do
      end do
      if (verbose .ge. 5) write(*,*) 'returning from integrate'
      return
   end subroutine integrate

!-----------------------------------------------------------------------
   subroutine int_ana(gana)
! analytically integrate the analytic singularity frsing
      use elitevac_data, only: npts,theta,dt,pi
      implicit none
!      include 'i_vac.f'
      real gana(npts,*)
      integer i,j
      real ti,tj,tij_right,tij_left,rij
      do i=1,npts
         ti=theta(i)
         do j=1,npts
            tj=theta(j)
            tij_right = tj+dt*0.5 - ti
            tij_left = tj-dt*0.5 - ti
            rij = 2.* tij_right * log( abs(tij_right) ) - &
                2.* tij_left * log( abs(tij_left) ) - 2.*dt
            tij_right = tij_right - 2.*pi
            tij_left = tij_left - 2.*pi
            rij = rij + &
                2.* tij_right * log( abs(tij_right) ) - &
                2.* tij_left * log( abs(tij_left) ) - 2.*dt
            tij_right = tij_right + 4.*pi
            tij_left = tij_left + 4.*pi
            rij = rij + &
                2.* tij_right * log( abs(tij_right) ) - &
                2.* tij_left * log( abs(tij_left) ) - 2.*dt
            gana(i,j)=rij
         end do
      end do
      return
   end subroutine int_ana



end program elitevac
