      
subroutine mapperb

!-----------------------------------------------------------------------
!     changed dpsi to dpsic  : rlm 7/3/96
!     changed dpsi2 to dpsiv : same
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
      real dpsic,xval(npsi),chipsipsi(npsi)
      integer itype,fail,failcount
      integer i,j
      real psival2(nx),dum2(nx),dum3(nx),dum,check
      real dpsiv
      real dchi,psitest,bnd_set
      real qpsip(nx) ! derivatives of equil
      real rhosq(nx),rho_tor(nx),torfO2pi
      real eps_percenflux

      eps_percenflux=1e-5
      failcount=0
      if (verbose .gt. 4) write(outfile,*) &
           'allocating mapper variables, npsi=',npsi, &
           ' npts=',npts
      call alloc_mapper

! efit psi grid is psival2
      dpsiv=(psilim-psiaxis)/(nx-1)
      do  i=1,nx
         psival2(i)=psiaxis+(i-1.)*dpsiv
      end do
      psival2(nx)=psilim

!  1/20/10 pbs:  calculate rho_tor grid corresponding to efit psi grid
      rhosq(1)=0.
      do i=2,nx
         rhosq(i)=rhosq(i-1)+0.5*(psival2(i)-psival2(i-1))*  &
              (qpsi(i)+qpsi(i-1))
      end do
      torfO2pi=rhosq(nx)
      do i=1,nx
         rhosq(i)=rhosq(i)/torfO2pi
      enddo
      rho_tor(1)=0.
      do i=2,nx
         rho_tor(i)=sqrt(rhosq(i))
      end do
      if (alpsi.eq.-3) then
         write(*,*) 'rho, psi fractions for efit surfaces'
         do i=nx,1,-1
            write(*,*) rho_tor(i),(psival2(i)-psiaxis)/(psilim-psiaxis)
         enddo
         stop
      endif


!     have trouble finding psilim contour on many eqdsks
!     percenflux allows edge scrapeoff--rlm 10/2/96
15    dchi=(psilim-psiaxis)*percenflux
      itype=0                   ! makes equally spaced psic
      if (verbose .gt. 4) write(outfile,*) 'call initpsi'
      call initpsi(psic,psiv,xval,chipsi,chipsipsi,dpsic,npsi, &
           alpsi,dchi,itype)

      do j=2,npsi
         psiv(j)=psiv(j)+psiaxis
      enddo
      psiv(1)=psiaxis
      psiv(npsi)=psiaxis+dchi
      if (verbose .gt. 3) write(outfile,*) 'psilim= ',psilim,'psiaxis= ',psiaxis, &
           'psiv(npsi)=',psiv(npsi)


! 3/8/00 pbs write out psi percentage for relevant flux surfaces
      if (verbose .gt. 3) write(outfile,*) &
           'psi fractions in outer ixinterp surfaces are:'
      if (verbose .gt. 3) write(outfile,*) ((psiv(i)-psiaxis)/(psilim-psiaxis), &
           i=npsi,npsi-nxinterp+1,-1)
      if (alpsi.eq.-2) then
         write(*,*) 'psi fractions for efit surfaces'
         do i=nx,1,-1
            write(*,*) (psival2(i)-psiaxis)/(psilim-psiaxis), &
                 (psiv(npsi-nx+i)-psiaxis)/(psilim-psiaxis)
         enddo
      endif


!
! set up pprime, ffprime, f, and qsfin arrays
      call spline1d(pprime,psiv,npsi,spp,psival2,nx,dum2)
      call spline1d(ffprime,psiv,npsi,sffp,psival2,nx,dum2)
      call spline1d(fval,psiv,npsi,sf,psival2,nx,dum2)
      call spline1d(qsfin,psiv,npsi,qpsi,psival2,nx,dum2)
      call spline1d(press,psiv,npsi,sp,psival2,nx,dum2)
      call spline1d(pw,psiv,npsi,pressw,psival2,nx,dum2)
      call spline1d(pwp,psiv,npsi,pwprim,psival2,nx,dum2)
      call spline1d(rho,psiv,npsi,rho0,psival2,nx,dum2)
      call spline1d(rhop,psiv,npsi,rho0p,psival2,nx,dum2)


!
! 10/00 pbs spline on original grid first to calculate q' q'' p'' (ff')'
      bnd_set=-2.d30  ! calculate boundary condition from nearest 4 pts
      call spline(psival2,qpsi,nx,bnd_set,bnd_set,dum2)
        ! dum2 contains the interpolating function
      do i=1,npsi
         call zsplint(psival2,qpsi,dum2,dum2,nx,psiv(i),check,qsfinp(i),dum)
!         write(*,*) 'i=',i,' qsfin=',qsfin(i),' check=',check, &
!              'qsfinp=',qsfinp(i)
      enddo

! calculate q'' by first splining q on old grid to keep all points
      call spline(psival2,qpsi,nx,bnd_set,bnd_set,dum2)
        ! dum2 contains the interpolating function
      do i=1,nx
         call zsplint(psival2,qpsi,dum2,dum2,nx,psival2(i),check,qpsip(i),dum)
!         write(*,*) 'i=',i,' qsfin=',qsfin(i),' check=',check, &
!              'qsfinp=',qsfinp(i)
      enddo
      call spline(psival2,qpsip,nx,bnd_set,bnd_set,dum2)
        ! dum2 contains the interpolating function
      do i=1,npsi
         call zsplint(psival2,qpsip,dum2,dum2,nx,psiv(i),check,qsfinpp(i),dum)
!         write(*,*) 'i=',i,' qsfinp=',qsfinp(i),' check=',check, &
!              'qsfinpp=',qsfinpp(i)
      enddo



      call spline(psival2,spp,nx,bnd_set,bnd_set,dum2)
        ! dum2 contains the interpolating function
      do i=1,npsi
         call zsplint(psival2,spp,dum2,dum2,nx,psiv(i),check,ppp(i),dum)
!         write(*,*) 'i=',i,' pprime=',pprime(i),' check=',check, &
!              'ppp=',ppp(i)
      enddo

      call spline(psival2,sffp,nx,bnd_set,bnd_set,dum2)
        ! dum2 contains the interpolating function
      do i=1,npsi
         call zsplint(psival2,sffp,dum2,dum2,nx,psiv(i),check,ffpp(i),dum)
!         write(*,*) 'i=',i,' ffprime=',ffprime(i),' check=',check, &
!              'ffpp=',ffpp(i)
      enddo

!     call spline(theta,rpts,npts,bnd_set,bnd_set,y2)
!     do i=1,npts ! second y2 arg is dummy as are rpt_v and int_v
!        call zsplint(theta,rpts,y2,y2,npts,theta(i), &
!            rpt_v,drdt(i),int_v)
!      end do


!     bicubic spline on psi used by cntour which is called by equalarc
      call bcspline(xgrid,zgrid,psixz,nx,nz,nx,csplpsi)
!     first surface is easy
      i=1
!      do j=1,nthet
      do j=1,npts
         xs(i,j)=xaxis
         zs(i,j)=zaxis
         bps(i,j)=0.
      end do
      if (verbose .gt. 4) write (outfile,*)'call equalarc'
      do i=2,npsi
         psitest=psiv(i)
!         write(6,'("i=",i5," psitest=",e12.4)') i,psitest
         psitest=psiv(i)
         call equalarc(psitest,i,fail)
         if (fail .eq. 1) then
            if (failcount .le. 10) then
               percenflux=percenflux-eps_percenflux
               failcount=failcount+1
               write(*,*) 'retrying with slightly reduced percenflux=',percenflux
               if (verbose .ge. 2) write(outfile,*) &
                    'retrying with slightly reduced percenflux=',percenflux
               goto 15
            else
               write(*,*) 'giving up, failure in equalarc, failcount=',failcount
               stop
            endif
         endif
      end do

!! if dens is true call routine to read mtanh coefficients 
!!   and calculate density and temperature profiles
!      if (dens) then
!         call readmtanh
!      endif

!
!     convert from MKS to cgs units
      call convert

    return
      
end subroutine mapperb


subroutine initpsi(psic,psiv,xval,chipsi,chipsipsi,dpsi &
          ,npsi,alpsi,dchi,itype)
!------------------------------------------------------------------------
!     subroutine which determines psi coordinate variables
!     called by init and by wrtcamino
!     itype=0 gives equally spaced psic
!     itype=1 gives equally spaced psiv
!     note that not all quantities are calculated for itype=1
!------------------------------------------------------------------------
      implicit none
!      real psic(1),psiv(1),xval(1),chipsi(1),chipsipsi(1)
      real psic(*),psiv(*),xval(*),chipsi(*),chipsipsi(*)
      real alpsi,dchi,dpsi
      integer npsi,itype
      integer i
      real fpi,xend,pi
      pi=acos(-1.)
      if(itype.eq.0) then
         dpsi=1./(npsi-1.)
         do  i=1,npsi
            psic(i)=(i-1)*dpsi
            if(alpsi.ge.0.) then
               xval(i)=(1.+alpsi)*psic(i)**2/(1.+alpsi*psic(i))
               chipsi(i)=dchi*(1.+alpsi)*(2.+alpsi*psic(i))* &
                   psic(i)/(1.+alpsi*psic(i))**2
               chipsipsi(i)=dchi*2.*(1.+alpsi) &
                   /(1.+alpsi*psic(i))**3
            else if (alpsi.eq.-2.) then
! for alpsi=-2. use equal spaced grid identical to that in efit file
               xval(i)=psic(i)
               chipsi(i)=dchi
               chipsipsi(i)=0.
            else
               if(alpsi.lt.(-1.)) then
                  print *,"error in alpsi"
                  stop
               endif
               fpi=-alpsi*pi
               xend=sin(fpi*0.5)**2
               xval(i)=(sin(psic(i)*fpi*0.5))**2/xend
               chipsi(i)=dchi*fpi*sin(psic(i)*fpi)*0.5/xend
               chipsipsi(i)=dchi*fpi**2*cos(psic(i)*fpi)*0.5/xend
            end if
            psiv(i)=xval(i)*dchi
         end do
         xval(npsi)=1.0
         psic(npsi)=1.
         psiv(npsi)=dchi
      else
         dpsi=1./(npsi-1.)
         do  i=2,npsi
            xval(i)=(i-1)*dpsi
            if(alpsi.ge.0.) then
               psic(i)=(alpsi*xval(i)+sqrt((alpsi*xval(i))**2+ &
                   4.*(1+alpsi)*xval(i)))/(2.*(1.+alpsi))
            else
               if(alpsi.lt.(-1.)) then
                  print *,"error in alpsi"
                  stop
               endif
               fpi=-alpsi*pi
               xend=sin(fpi*0.5)**2
               psic(i)=asin(sqrt(xval(i)*xend))*2./fpi
            end if
            psiv(i)=xval(i)*dchi
         end do
         xval(1)=0.
         psic(1)=0.
         psiv(1)=0.
         xval(npsi)=1.0
         psic(npsi)=1.
         psiv(npsi)=dchi
      end if
    return
      
end subroutine initpsi



subroutine equalarc(psivin,ipsi,fail)

!-----------------------------------------------------------------------
!     stripped down version of chali from mbc
!     rlm 6/22/94
!     this routine finds psi contour psivin=psiv(ipsi) and return through
!     common blocks the equal arc length values of xs(ipsi,j) and zs(ipsi,j)
!     arcsur(ipsi,j), the arc length along the contour, is also returned
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
!      include 'cliche.f'
      real psivin
      integer ipsi,fail

      integer is,js,npc,ntw,npco,npcm1
      real delxz
      real xp(6*(nx+nz)),zp(6*(nx+nz)),bpsqp(6*(nx+nz)),tp(6*(nx+nz))
!      real csx(3,nxzd),csz(3,nxzd),csbp(3,nxzd)
      real csx(6*(nx+nz)),csz(6*(nx+nz)),csbp(6*(nx+nz))
      real arc(6*(nx+nz))
      real dl,ds,bpsqval
      integer j,nlow
!      real sterpl,f2s,bb(4)
!      integer ier
!-----------------------------------------------------------------------
!     the following common is shared with routine cntour
!-----------------------------------------------------------------------
!  shared in new module contour 2/3/00
      real xemin,xemax          !output form cntour
      real yemin,yemax          !output form cntour
      real yxmin,yxmax          !output form cntour
      real xymin,xymax          !output form cntour
      real xaxd,yaxd            !input to cntour
      real dang,arcl,bperr      !input to cntour
      real xmin,xmax            !input to cntour
      real ymin,ymax            !input to cntour
      common/cntd/xaxd,yaxd, &
          xemin,xemax,yemin,yemax,yxmin,yxmax,xymin,xymax, &
          dang,arcl,bperr,xmin,xmax,ymin,ymax
!-----------------------------------------------------------------------
!     attempt to find contour using furpl
!     is and js are starting values for search for flux surfaces
!-----------------------------------------------------------------------

      integer nxzd

      fail=0    ! set to one as flag for failure of equalarc
      nxzd=6*(nx+nz)

      npc=1
      ntw=nxzd
      is=nx/2
      js=nz/2
      if (verbose .ge. 5) write(*,*) 'about to call furplm in equalarc 1',npc
      call furplm(psixz,bpsq,xgrid,zgrid,psivin,xp,zp,bpsqp,npc, &
          nx,nz,ntw,nx,nz,is,js)
!      npc=0
!      npfit=150
      npfit=40
      if (verbose .ge. 5) write(*,*) 'return from furplm npc=',npc,' npfit=',npfit
      if(npc.le.npfit) then
!         write(6,'("used contour",2i5)') ipsi,npc
         xaxd=xaxis
         yaxd=zaxis
         npco=npc
         xmin=xgrid(1)
         xmax=xgrid(nx)
         ymin=zgrid(1)
         ymax=zgrid(nz)
         bperr=0.1
!         bperr=0.001
         if(ipsi.eq.2) then
            arcl=amax1(dx,dz)*0.1
!            arcl=amax1(dx,dz)*0.01
         else
            arcl=arcsur(ipsi-1,npts)/npfit
         endif
         npfit=400
         dang=1./npfit
         call cntourp(xgrid,nx,zgrid,nz,csplpsi,xp,zp,bpsqp, &
             npc,dx,dz,nxzd,psivin,nx)
!         call cntourp(xgrid,nx,zgrid,nz,csplpsi,xp,zp,bpsqp, &
!             npc,dx/5.,dz/5.,nxzd,psivin,nx)
         call sortr(xp,zp,bpsqp,npc,xaxis,zaxis,xaxis,zaxis,0,1)
      else
         if(xp(1).ne.xp(npc).or.zp(1).ne.zp(npc)) then
            write(6,'("error in equalarc finding contour ", &
                &i3," of ",i3)') ipsi,npsi
            write(6,'("if this is final contour," &
                &," try reducing percenflux")')
            fail=1
!            stop
!   3/14 retry failures with modified percenflux value
            return
         endif
         call sortr(xp,zp,bpsqp,npc,xaxis,zaxis,xaxis,zaxis,0,1)
         delxz=0.20*amin1(xgrid(2)-xgrid(1),zgrid(2)-zgrid(1))
!         delxz=0.020*amin1(xgrid(2)-xgrid(1),zgrid(2)-zgrid(1))
         call pack(bpsqp,xp,zp,npc,delxz)
!         write(*,*) 'finished with pack, npc=',npc
      endif
      if(npc.le.5) then
         write(6,'("error in equalarc--npc too small=")') npc
         stop
      endif
!-----------------------------------------------------------------------
!     set up arclengths on surface
!-----------------------------------------------------------------------
      tp(1)=0.
      npcm1=npc-1
      do  j=2,npc
         tp(j)=tp(j-1)+sqrt((xp(j)-xp(j-1))**2+(zp(j)-zp(j-1))**2)
      end do
!-----------------------------------------------------------------------
!     spline xp wrt t 
!-----------------------------------------------------------------------
      call spline(tp,xp,npc,-1.e31,-1.e31,csx)

!-----------------------------------------------------------------------
!     spline zp wrt t 
!-----------------------------------------------------------------------
      call spline(tp,zp,npc,-1.e31,-1.e31,csz)

!-----------------------------------------------------------------------
!     calculate arclength to each point
!-----------------------------------------------------------------------
      call garc(tp,xp,zp,csx,csz,nxzd,arc,npc)
!-----------------------------------------------------------------------
!     spline xp wrt arc
!-----------------------------------------------------------------------
      call spline(arc,xp,npc,-1.e31,-1.e31,csx)

!-----------------------------------------------------------------------
!     spline zp wrt arc
!-----------------------------------------------------------------------
      call spline(arc,zp,npc,-1.e31,-1.e31,csz)

!-----------------------------------------------------------------------
!     spline bpsqp wrt arc
!-----------------------------------------------------------------------
      call spline(arc,bpsqp,npc,-1.e31,-1.e31,csbp)

!-----------------------------------------------------------------------
!     convert to equal arc length along contour
!-----------------------------------------------------------------------
      ds=arc(npc)/(npts-1.)
      do  j=1,npts-1
         dl=(j-1)*ds
         nlow=0
         call splint(arc,xp,csx,npc,dl,xs(ipsi,j),nlow)
         call splint(arc,zp,csz,npc,dl,zs(ipsi,j),nlow)
         call splint(arc,bpsqp,csbp,npc,dl,bpsqval,nlow)
         if(bpsqval.gt.0.) then
            bps(ipsi,j)=sqrt(bpsqval)
         else
            bps(ipsi,j)=0.
         endif
         arcsur(ipsi,j)=dl
      end do
      xs(ipsi,npts)=xs(ipsi,1)
      zs(ipsi,npts)=zs(ipsi,1)
      bps(ipsi,npts)=bps(ipsi,1)
      arcsur(ipsi,npts)=arc(npc)
    return
      
end subroutine equalarc


subroutine convert
      use eliteeq_data
      implicit none
!     to convert eqdsk quantities from MKS to cgs units
!     02/10/97 add conversion of arcsur and bps

      real convpsi,convpp,convf,convffp,convchi,convm,convb
      parameter (convpsi=1.e8,convpp=1.e-7,convf=1.e6,convffp=1.e4, &
          convchi=convpsi,convm=1.e2,convb=1.e4)
      integer i,j

      do i=1,npsi
         psiv(i)=convpsi*psiv(i)
         press(i)=convpp*press(i)*convpsi
         pprime(i)=convpp*pprime(i)
         fval(i)=convf*fval(i)
         ffprime(i)=convffp*ffprime(i)
         chipsi(i)=convchi*chipsi(i)
         qsfinp(i)=qsfinp(i)/convpsi
         qsfinpp(i)=qsfinpp(i)/(convpsi*convpsi)
         ppp(i)=ppp(i)*convpp/convpsi
         ffpp(i)=ffpp(i)*convffp/convpsi
      enddo

      pressval=press

      do j=1,npts
         do i=1,npsi
            xs(i,j)=convm*xs(i,j)
            zs(i,j)=convm*zs(i,j)
            arcsur(i,j)=convm*arcsur(i,j)
            bps(i,j)=convb*bps(i,j)
         enddo
      enddo
!      xaxis=convm*xaxis
!      zaxis=convm*zaxis
    return
      
end subroutine convert


subroutine readmtanh
! pbs: routine for readint mtanh coefficients for ne, Te, Ti
   use eliteeq_data
   implicit none
   real coeffs(5,3),x,f,fp,z,a0,a1,a2,a3,a4
   integer nmtanh,i,j,k,ios

   nmtanh=43

   open(unit=nmtanh,file=runname(1:lrunname)//'.mtanh', &
        status='old',iostat=ios)
   if(ios==0 .and. (verbose .ge. 1)) write(*,*) 'reading mtanh file ', &
        runname(1:lrunname)//'.mtanh'
   if(ios.ne.0) then
      if (verbose .ge. 1) then
         write(6,*) 'problem opening file ',runname(1:lrunname)//'.mtanh'
         write(6,*) 'will try to read file named mtanh instead'
      endif
      open(unit=nmtanh,file='mtanh',status='old',iostat=ios)
      if(ios.ne.0) then
         if (verbose .ge. 1) then
            write(*,*) 'could not open mtanh file'
            write(*,*) 'will attempt to read density and temp from peqdsk file'
         endif
         call readpeqdsk
         return
!         stop
      endif
      if (verbose .ge. 1) write(*,*) 'reading mtanh coefficients from file mtanh'
   endif

   do i=1,3
      read(nmtanh,*) (coeffs(j,i),j=1,5)
   enddo

   if (verbose .ge. 2) then
      write(*,*) 'ane=',coeffs(:,1)
      write(*,*) 'ate=',coeffs(:,2)
      write(*,*) 'ati=',coeffs(:,3)
   endif
   if (verbose .gt. 3) then
      write(outfile,*) 'ane=',coeffs(:,1)
      write(outfile,*) 'ate=',coeffs(:,2)
      write(outfile,*) 'ati=',coeffs(:,3)
   endif

   do k=1,npsi
      x=(psiv(k)/1.d8-psiaxis)/(psilim-psiaxis)
!      x=psiv(k)
      do i=1,3
         a0=coeffs(1,i)
         a1=coeffs(2,i)
         a2=coeffs(3,i)
         a3=coeffs(4,i)
         a4=coeffs(5,i)
         z=(x-a2)/a3
!         write(*,*) 'k=',k,x,z
         f=a0-a1*tanh(z)-a1*a4*z*exp(-z)/(exp(-z)+exp(z))
         fp=-a1/(a3*cosh(z)**2)-a1*a4*exp(-z)* &
              (1-z-z*tanh(z))/(a3*(exp(z)+exp(-z)))
         fp=fp/(1.d8*(psilim-psiaxis))  ! unnormalize d psi
         if (i==1) then
            nel(k)=f*1.d14  ! density in 10^14 cm^-3
            nprime(k)=fp*1.d14  
         else if (i==2) then
            tel(k)=f  ! electron temperature in eV
            teprime(k)=fp
         else if (i==3) then
            tion(k)=f
            tiprime(k)=fp
         endif
      enddo
      if (verbose .ge. 4) write(*,*) k,x,nel(k),tel(k),tion(k)
   enddo

   return

end subroutine readmtanh




subroutine readpeqdsk
! pbs 7/06: routine for reading ne, Te, Ti from peqdsk file
   use eliteeq_data
   implicit none
   integer nmp,i,j,k,ios
   character*8 line(10)
   integer nxmax
!   parameter (nxmax=256)
!   real psiread(nxmax),fread(nxmax),freadp(nxmax),psival2(nxmax)
   real, dimension(:), allocatable :: psiread,fread,freadp,dum2
   real psival2(nx),psinorm2(nx)
   real dpsiv,psierr,psitol  ! ,dum2(nxmax)
   real nil(npsi),niprime(npsi),psivnorm(npsi),pprime_est,pratio(npsi)
   integer pbug,n_ne,n_te,n_ni,n_ti

!   pbug=0
! 6/07  pbug=1 is for Osborne's p file format with normalized psi values
!    and derivatives wrt normalized psi, now standard
   pbug=1

   n_ne=nx; n_te=nx; n_ni=nx; n_ti=nx

   psitol=1.e-5
   nmp=46

   open(unit=nmp,file=runname(1:lrunname)//'.peqdsk', &
        status='old',iostat=ios)
   if(ios==0 .and. (verbose .ge. 1)) write(*,*) 'reading peqdsk file ', &
        runname(1:lrunname)//'.peqdsk'
   if(ios.ne.0) then
      if (verbose .ge. 1) then
         write(6,*) 'problem opening file ',runname(1:lrunname)//'.peqdsk'
         write(6,*) 'will try to read file named peqdsk instead'
      endif
      open(unit=nmp,file='peqdsk',status='old',iostat=ios)
      if(ios.ne.0) then
         write(*,*) 'could not open peqdsk file'
         write(*,*) 'error: unable to find density and temp data and dens=.true.'
         write(*,*) 'WARNING: USING DEFAULT FLAT DENSITY PROFILE!!'
         if (verbose .gt. 2) then
            write(outfile,*) 'could not open peqdsk file'
            write(outfile,*) 'error: unable to find density and temp data and dens=.true.'
            write(outfile,*) 'WARNING: USING DEFAULT FLAT DENSITY PROFILE!!'
         endif
         nel=1.d14
         nprime=0.
         tel=1.d3
         teprime=0.
         tion=1.d3
         tiprime=0.
         nil=1.d14
         niprime=0.
         return
!         stop
      endif
      if (verbose .ge. 1) write(*,*) 'reading density and temp data from file peqdsk'
   endif

! set up g-file psi values as a check
   dpsiv=(psilim-psiaxis)/(nx-1)
   do  i=1,nx
      psival2(i)=psiaxis+(i-1.)*dpsiv
   end do
   psival2(nx)=psilim
! convert psival2 into cgs units (like psiv)
   psival2=psival2*1.d8
   do i=1,nx
      psinorm2(i)=abs(psival2(i)-psival2(1))/abs(psival2(nx)-psival2(1))
   enddo
   if (verbose .gt. 3) write(outfile,*) 'psinorm2=',psinorm2
   do i=1,npsi
      psivnorm(i)=abs(psiv(i)/1d8-psiaxis)/abs(psilim-psiaxis)
   enddo
   if (verbose .gt. 3) write(outfile,*) 'psivnorm=',psivnorm

!!!!! Read n_e and n_e_prime data
   if (pbug==0) then
      read(nmp,'(10a8)') (line(i),i=1,10)
      if (verbose .ge. 4) write(*,*) 'line=',line
   else
      read(nmp,*) n_ne
      if (verbose .ge. 4) write(*,*) 'read n_ne=',n_ne
   endif
   allocate( psiread(n_ne),fread(n_ne),freadp(n_ne),dum2(n_ne))
   do i=1,n_ne
      read(nmp,*) psiread(i),fread(i),freadp(i)
!      write(*,*) i,(psiread(i)-psiaxis)/(psilim-psiaxis),fread(i),freadp(i)
!  convert psiread into cgs units
      if (pbug==0) psiread(i)=psiread(i)*1.d8
      psierr=abs((psiread(i)-psival2(i))/psiread(i))
!      write(*,*) 'psicheck ',psiread(i),psival2(i),psiv(i)
      if (psierr.gt.psitol .and. pbug==0) then
         if (verbose .ge. 1) write(*,*) &
              '!!!psi read from p file differs from g file by', &
              psierr,' at surface=',i,' psi=',psiread(i)
      endif      
   enddo
! spline values onto mapping surfaces
!   first put psival2 into cgs units as psiv is
!   call spline1d(nel,psiv,npsi,fread,psival2,nx,dum2)
!   call spline1d(nprime,psiv,npsi,freadp,psival2,nx,dum2)
   if (pbug==0) then
      call spline1d(nel,psiv,npsi,fread,psiread,n_ne,dum2)
      call spline1d(nprime,psiv,npsi,freadp,psiread,n_ne,dum2)
   else
      call spline1d(nel,psivnorm,npsi,fread,psiread,n_ne,dum2)
      call spline1d(nprime,psivnorm,npsi,freadp,psiread,n_ne,dum2)
   endif
! convert into CGS units
   nel=nel*1.d14
   nprime=nprime*1.d14/1.d8
! derivatives are with respect to normalzed psi for pbug=1,unnormalize
   if (pbug==1) nprime=nprime/(psilim-psiaxis)
!   write(*,*) 'nel=',nel
!   write(*,*) 'nprime=',nprime
   if (verbose .gt. 3) write(outfile,*) 'psival2, psiread, psivnorm, nel, fread'
   do i=1,n_ne
      if (verbose .gt. 3) write(outfile,*) &
           psival2(i),psiread(i),psivnorm(i),nel(i)/1.d14,fread(i)
   enddo

!!!!! Read T_e and T_e_prime data
   if (pbug==0) then
      read(nmp,'(10a8)') (line(i),i=1,10)
      if (verbose .ge. 4) write(*,*) 'line=',line
   else
      read(nmp,*) n_te
      if (verbose .gt. 3) write(outfile,*) 'read n_te=',n_te
      deallocate(psiread,fread,freadp,dum2)
      allocate( psiread(n_te),fread(n_te),freadp(n_te),dum2(n_te))
   endif
   do i=1,n_te
      read(nmp,*) psiread(i),fread(i),freadp(i)
!      write(*,*) i,(psiread(i)-psiaxis)/(psilim-psiaxis),fread(i),freadp(i)
!  convert psiread into cgs units
      if (pbug==0) psiread(i)=psiread(i)*1.d8
      psierr=abs((psiread(i)-psival2(i))/psiread(i))
!      write(*,*) 'psicheck ',psiread(i),psival2(i),psiv(i)
      if (psierr.gt.psitol .and. pbug==0) then
         if (verbose .ge. 1) write(*,*) &
              '!!!psi read from p file differs from g file by', &
              psierr,' at surface=',i,' psi=',psiread(i)
      endif      
   enddo
! spline values onto mapping surfaces
!   first put psival2 into cgs units as psiv is
   if (pbug==0) then
      call spline1d(tel,psiv,npsi,fread,psiread,n_te,dum2)
      call spline1d(teprime,psiv,npsi,freadp,psiread,n_te,dum2)
   else
      call spline1d(tel,psivnorm,npsi,fread,psiread,n_te,dum2)
      call spline1d(teprime,psivnorm,npsi,freadp,psiread,n_te,dum2)
   endif
! convert into CGS units
   tel=tel*1.d3
   teprime=teprime*1.d3/1.d8
! derivatives are with respect to normalzed psi for pbug=1,unnormalize
   if (pbug==1) teprime=teprime/(psilim-psiaxis)
   if (verbose .gt. 3) then
      write(outfile,*) 'tel=',tel
      write(outfile,*) 'teprime=',teprime
   endif

!!!!! Read n_i and n_i_prime data
!  at present these are not used, n_i=n_e still assumed
   if (pbug==0) then
      read(nmp,'(10a8)') (line(i),i=1,10)
      if (verbose .ge. 4) write(*,*) 'line=',line
   else
      read(nmp,*) n_ni
      if (verbose .gt. 3) write(outfile,*) 'read n_ni=',n_ni
      deallocate(psiread,fread,freadp,dum2)
      allocate( psiread(n_ni),fread(n_ni),freadp(n_ni),dum2(n_ni))
   endif
   do i=1,n_ni
      read(nmp,*) psiread(i),fread(i),freadp(i)
!      write(*,*) i,(psiread(i)-psiaxis)/(psilim-psiaxis),fread(i),freadp(i)
!  convert psiread into cgs units
      if (pbug==0) psiread(i)=psiread(i)*1.d8
      psierr=abs((psiread(i)-psival2(i))/psiread(i))
!      write(*,*) 'psicheck ',psiread(i),psival2(i),psiv(i)
      if (psierr.gt.psitol .and. pbug==0) then
         if (verbose .ge. 1) write(*,*) &
              '!!!psi read from p file differs from g file by', &
              psierr,' at surface=',i,' psi=',psiread(i)
      endif      
   enddo
! spline values onto mapping surfaces
!   first put psival2 into cgs units as psiv is
   if (pbug==0) then
      call spline1d(nil,psiv,npsi,fread,psiread,n_ni,dum2)
      call spline1d(niprime,psiv,npsi,freadp,psiread,n_ni,dum2)
   else
      call spline1d(nil,psivnorm,npsi,fread,psiread,n_ni,dum2)
      call spline1d(niprime,psivnorm,npsi,freadp,psiread,n_ni,dum2)
   endif
! convert into CGS units
   nil=nil*1.d14
   niprime=niprime*1.d14/1.d8
! derivatives are with respect to normalzed psi for pbug=1,unnormalize
   if (pbug==1) niprime=niprime/(psilim-psiaxis)
   if (verbose .gt. 3) write(outfile,*) 'nil=',nil
   if (verbose .gt. 3) write(outfile,*) 'niprime=',niprime
   if (verbose .ge. 1) write(*,*) &
        '!!!! ion density data read from p eqdsk file, but is not', &
        ' presently used in elite, which assumes n_i=n_e in mass density', &
        ' (ie use appropriate ion_mass) and omega* calculations'
   if (verbose .gt. 3) write(outfile,*) &
        '!!!! ion density data read from p eqdsk file, but is not', &
        ' presently used in elite, which assumes n_i=n_e in mass density', &
        ' (ie use appropriate ion_mass) and omega* calculations'


!!!!! Read T_i and T_i_prime data
   if (pbug==0) then
      read(nmp,'(10a8)') (line(i),i=1,10)
      if (verbose .ge. 4) write(*,*) 'line=',line
   else
      read(nmp,*) n_ti
      if (verbose .gt. 3) write(outfile,*) 'read n_ti=',n_ti
      deallocate(psiread,fread,freadp,dum2)
      allocate( psiread(n_ti),fread(n_ti),freadp(n_ti),dum2(n_ti))
   endif
   do i=1,n_ti
      read(nmp,*) psiread(i),fread(i),freadp(i)
!      write(*,*) i,(psiread(i)-psiaxis)/(psilim-psiaxis),fread(i),freadp(i)
!  convert psiread into cgs units
      if (pbug==0) psiread(i)=psiread(i)*1.d8
      psierr=abs((psiread(i)-psival2(i))/psiread(i))
!      write(*,*) 'psicheck ',psiread(i),psival2(i),psiv(i)
      if (psierr.gt.psitol .and. pbug==0) then
         if (verbose .ge. 1) write(*,*) '!!!psi read from p file differs from g file by', &
              psierr,' at surface=',i,' psi=',psiread(i)
      endif      
   enddo
! spline values onto mapping surfaces
!   first put psival2 into cgs units as psiv is
   if (pbug==0) then
      call spline1d(tion,psiv,npsi,fread,psiread,n_ti,dum2)
      call spline1d(tiprime,psiv,npsi,freadp,psiread,n_ti,dum2)
   else
      call spline1d(tion,psivnorm,npsi,fread,psiread,n_ti,dum2)
      call spline1d(tiprime,psivnorm,npsi,freadp,psiread,n_ti,dum2)
   endif
! convert into CGS units
   tion=tion*1.d3
   tiprime=tiprime*1.d3/1.d8
! derivatives are with respect to normalzed psi for pbug=1,unnormalize
   if (pbug==1) tiprime=tiprime/(psilim-psiaxis)
   if (verbose .gt. 3) write(outfile,*) 'tion=',tion
   if (verbose .gt. 3) write(outfile,*) 'tiprime=',tiprime
   
   if (verbose .gt. 3) write(outfile,*) 'psilim-psiaxis=',(psilim-psiaxis)
   if (verbose .gt. 3) write(outfile,*) 'pprime  ne*Tep Te*nep ni*Tip Ti*nip  tot'
   do i=1,npsi
      pprime_est=1.6022e-12*(nel(i)*teprime(i)+tel(i)*nprime(i)+ &
           nil(i)*tiprime(i)+tion(i)*niprime(i))
      if (verbose .gt. 3) write(outfile,*) pprime(i),1.6022e-12*nel(i)*teprime(i), &
           1.6022e-12*tel(i)*nprime(i),1.6022e-12*nil(i)*tiprime(i), &
           1.6022e-12*tion(i)*niprime(i),pprime_est 
      pratio(i)=pprime(i)/pprime_est
!      write(*,*) i,'ratio of pprime/pprime_est=',pprime(i)/pprime_est
!      write(outfile,*) i,'ratio of pprime/pprime_est=',pprime(i)/pprime_est
   enddo
   if (verbose .ge. 4) write(*,*) 'ratios of pprime/(ne*Tep+Te*nep+ni*Tip+Ti*nip)',pratio
   if (verbose .gt. 3) write(outfile,*) &
        'ratios of pprime/(ne*Tep+Te*nep+ni*Tip+Ti*nip)',pratio


end subroutine readpeqdsk

