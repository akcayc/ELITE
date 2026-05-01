
subroutine readeq

!-----------------------------------------------------------------------
! New equilibrium  format: HRW 11/11/02
! Load in SI units, except for pressure->mu0*pressure, and T in keV
! Psi mesh is arbitrary, but poloidal mesh should be equal arc length
! Poloidal angle chi=2*pi*l/L (L=full circumference)
! First and last points on the poloidal mesh are identical (ie 1 and npts
! array points)
!
! PBS 07/08/04 convert to CGS units internally for consistency
!   with other formats in unnormalized quantities such as gamr
!   and R/B.  Also add 'eqin' as default equilibrium file name
!   in case name.eqin does not exist.
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
      integer i,j,k
      integer noutbal,neqin,ios
      integer old_npts
      character*30 xtitle
      character*8 parlab
      double precision var1,var2
      real bndset,dum
      real, dimension(:), allocatable :: za,pspline


      noutbal=12
      neqin=22
      old_npts=npts

!      write(*,*) 'lrunname=',lrunname,'runname=',runname
      open(unit=neqin,file=runname(1:lrunname)//'.eqin', &
           status='old',iostat=ios)
      if(ios==0) write(*,*) 'reading eqin file ', &
           runname(1:lrunname)//'.eqin'
      if(ios.ne.0) then
        write(6,*) 'problem opening file ',runname(1:lrunname)//'.eqin'
        write(6,*) 'will try to read file named eqin instead'
         open(unit=neqin,file='eqin',status='old',iostat=ios)
         if(ios.ne.0) then
           write(*,*) 'could not open eqin file'
           stop
         endif
         write(*,*) 'reading equilibrium from file eqin'
      endif

      read(neqin,'(a30)')xtitle
! No. flux surface, no. poloidal mesh points
      read(neqin,'(2i5)') npsi,npts
      if (npts.ne.old_npts) then
         write(*,*) 'npts in input file must =npts from eqin'
         write(*,*) 're-set npts to',npts,' in your .in file'
         stop
      endif

      if (verbose .gt. 4) write(outfile,*) 'allocating equilibrium input variables'

      call alloc_toq
      allocate( pspline(npsi),za(npsi) )
!      allocate( rl(npts),zl(npts),bpl(npts),work(npts),leng(npts),leq(npts) )

!  Psi on mesh: 0 on axis, max at edge, flux per radian
      read(neqin,'(a8)') parlab
      do i=1,npsi/5
         read(neqin,200) psiv((i*5)-4:i*5)
      enddo
      read(neqin,200) psiv((npsi/5)*5+1:npsi)
!      write(noutbal,*) psiv

!  mu0*dp/dpsi
      read(neqin,'(a8)') parlab
      do i=1,npsi/5
         read(neqin,200) pprime((i*5)-4:i*5)
      enddo
      read(neqin,200) pprime((npsi/5)*5+1:npsi)
      pprime=pprime/(4.*pi)
!      write(noutbal,*) pprime

!  mu0*d^2p(psi)/dpsi^2
      read(neqin,'(a8)') parlab
      do i=1,npsi/5
         read(neqin,200) ppp((i*5)-4:i*5)
      enddo
      read(neqin,200) ppp((npsi/5)*5+1:npsi)
      ppp=ppp/(4.*pi)
!      ppp=-ppp
!      write(noutbal,*) ppp

! f(psi)=R*Bphi
      read(neqin,'(a8)') parlab
      do i=1,npsi/5
         read(neqin,200) fval((i*5)-4:i*5)
      enddo
      read(neqin,200) fval((npsi/5)*5+1:npsi)
!      fval=-fval
!      write(noutbal,*) fval

! f*df/dpsi
      read(neqin,'(a8)') parlab
      do i=1,npsi/5
         read(neqin,200) ffprime((i*5)-4:i*5)
      enddo
      read(neqin,200) ffprime((npsi/5)*5+1:npsi)
!      ffprime=-ffprime
!      write(noutbal,*) ffprime

! d/dpsi(f*df/dpsi)
      read(neqin,'(a8)') parlab
      do i=1,npsi/5
         read(neqin,200) ffpp((i*5)-4:i*5)
      enddo
      read(neqin,200) ffpp((npsi/5)*5+1:npsi)
!      ffpp=-ffpp
!      write(noutbal,*) ffpp

!  q(psi): not needed, just used for consistency check on data
      read(neqin,'(a8)') parlab
      do i=1,npsi/5
         read(neqin,200) qsfin((i*5)-4:i*5)
      enddo
      read(neqin,200) qsfin((npsi/5)*5+1:npsi)
!      write(401,*)' q'
!      write(401,*) qsfin
!      write(noutbal,*) qsfin

! R(psi,chi)
      read(neqin,'(a8)') parlab
!      k=1
!      read(neqin,*)var1
      do j=1,npts
!         k=k+1
         do i=1,npsi/5
            read(neqin,200) xs((i*5)-4:i*5,j)
         enddo
         read(neqin,200) xs((npsi/5)*5+1:npsi,j)
!         read(neqin,*)var2
!         if (var1.eq.var2) write(6,*)' k=',k
      enddo
!      write(noutbal,*) xs

! Z(psi,chi)
      read(neqin,'(a8)') parlab
      do j=1,npts
         do i=1,npsi/5
            read(neqin,200) zs((i*5)-4:i*5,j)
         enddo
         read(neqin,200) zs((npsi/5)*5+1:npsi,j)
      enddo
!      write(noutbal,*) zs

! Bp(psi,chi)
      read(neqin,'(a8)') parlab
      do j=1,npts
         do i=1,npsi/5
            read(neqin,200) bps((i*5)-4:i*5,j)
         enddo
         read(neqin,200) bps((npsi/5)*5+1:npsi,j)
      enddo
!      write(noutbal,*) bps
!
!  If dens is set to .true. then read density, Te and Ti profiles
      if (dens) then

! electron density in m**-3
        read(neqin,'(a8)') parlab
        do i=1,npsi/5
           read(neqin,200) nel((i*5)-4:i*5)
        enddo
        read(neqin,200) nel((npsi/5)*5+1:npsi)
!        write(noutbal,*) nel
!
! dn/dpsi
        read(neqin,'(a8)') parlab
        do i=1,npsi/5
           read(neqin,200) nprime((i*5)-4:i*5)
        enddo
        read(neqin,200) nprime((npsi/5)*5+1:npsi)
!        write(noutbal,*) nprime
!
! Te in keV
        read(neqin,'(a8)') parlab
        do i=1,npsi/5
           read(neqin,200) tel((i*5)-4:i*5)
        enddo
        read(neqin,200) tel((npsi/5)*5+1:npsi)
!        write(noutbal,*) tel
!
!  dTe/dpsi
        read(neqin,'(a8)') parlab
        do i=1,npsi/5
           read(neqin,200) teprime((i*5)-4:i*5)
        enddo
        read(neqin,200) teprime((npsi/5)*5+1:npsi)
!        write(noutbal,*) teprime
!
! Ti in keV 
        read(neqin,'(a8)') parlab
        do i=1,npsi/5
           read(neqin,200) tion((i*5)-4:i*5)
        enddo
        read(neqin,200) tion((npsi/5)*5+1:npsi)
!        write(noutbal,*) tion
!
!   dTi/psi
        read(neqin,'(a8)') parlab
        do i=1,npsi/5
           read(neqin,200) tiprime((i*5)-4:i*5)
!           write(*,*) i,tiprime((i*5)-4:i*5)
        enddo
        read(neqin,200) tiprime((npsi/5)*5+1:npsi)
!        write(noutbal,*) tiprime
      end if

      write(6,*) 'finished reading equilibrium file'
      if (verbose .gt. 4) write(outfile,*) 'finished reading equilibrium file'


      close(neqin)
200   format(5g20.12)

!  convert to Gaussian (CGS) units for consistency with other
!    formats.  This is not strictly necessary, but it allows
!    unnormalized quantities such as gamr to be same size

     call convert_eq


!-----------------------------------------------------------------------
! calculate pressure, added from readbal,  PBS 7.8.04
!-----------------------------------------------------------------------
     if (dens) then
        do i=1,npsi
           pressval(i)=nel(i)* &
                (tel(i)+tion(i))*1.6022e-12
        enddo
     else
        bndset=-2.d30
        call spline(psiv,pprime,npsi,bndset,bndset,pspline)
        call zspline(psiv,pprime,pspline,npsi,0,npsi, &
             za )
        do i=1,npsi
           call zsplint(psiv,pprime,pspline,za,npsi,psiv(i), &
                dum,dum,pressval(i))
           if (pressval(i) < 0.) then
              write(*,*) 'warning, splined pressure<0 at i=',i, &
                   ' resetting to zero'
              pressval(i)=0.
           endif
        enddo
     endif

! set surface 1 to axis, check tolerance
     do i=1,npts
        if (abs(xs(1,i)-xs(1,1)).gt.1e-12) then
           write(*,*) 'in readeq, inner surface not on axis, i=,',i, &
                ' xs(1,i)=',xs(1,i),' xs(1,1)=',xs(1,1)
           if (verbose .gt. 3) &
                write(outfile,*) 'in readeq, inner surface not on axis, i=,',i, &
                ' xs(1,i)=',xs(1,i),' xs(1,1)=',xs(1,1)
        endif
        if (abs(zs(1,i)-zs(1,1)).gt.1e-12) then
           write(*,*) 'in readeq, inner surface not on axis, i=,',i, &
                ' zs(1,i)=',zs(1,i),' zs(1,1)=',zs(1,1)
           if (verbose .gt. 3) write(outfile,*) &
                'in readeq, inner surface not on axis, i=,',i, &
                ' zs(1,i)=',zs(1,i),' zs(1,1)=',zs(1,1)
        endif
        if (abs(bps(1,i)-bps(1,1)).gt.1e-12) then
           write(*,*) 'in readeq, bps varies on axis, i=,',i, &
                ' bps(1,i)=',bps(1,i),' bps(1,1)=',bps(1,1)
           if (verbose .gt. 3) write(outfile,*) &
                'in readeq, bps varies on axis, i=,',i, &
                ' bps(1,i)=',bps(1,i),' bps(1,1)=',bps(1,1)
        endif
        xs(1,i)=xs(1,1)
        zs(1,i)=zs(1,1)
        bps(1,i)=bps(1,1)
     enddo


! possibly temporary, (5/15/05) remap onto equal arc grid
     call equal_arc
   
  return
      
end subroutine readeq



subroutine convert_eq
  use eliteeq_data
  implicit none
! converts from MKS to CGS units for eqin files 7/8/04

  real convpsi,convpp,convf,convffp,convchi,convm,convb
  parameter (convpsi=1.e8,convpp=1.e-7,convf=1.e6,convffp=1.e4, &
       convchi=convpsi,convm=1.e2,convb=1.e4)

  psiv=psiv*convpsi
! pprime is mu0*pprime, already has 10^-7 in it
  ppp=ppp/convpsi
  fval=fval*convf
  ffprime=ffprime*convffp
  ffpp=ffpp*convffp/convpsi
  xs=xs*convm
  zs=zs*convm
  bps=bps*convb
  nel=nel/convm**3
  nprime=nprime/(convm**3*convpsi)
  teprime=teprime/convpsi
  tiprime=tiprime/convpsi

  return

end subroutine convert_eq


subroutine equal_arc

  use eliteeq_data
  implicit none
! Convert R,Z,Bp to equal arc length grid via spline

  integer i,j
  real, dimension(npts) :: rl,zl,bpl,work,leng,leq
  real dr,drt,dl

! Remap eqin grid onto an equal arc grid
!   Gato needs this (wrtgato option) 


  write(*,*) '!!re-mapping eqin grid onto equal arc lengths'
  if (verbose .gt. 3) write(outfile,*) '!!re-mapping eqin grid onto equal arc lengths'
  do j=2,npsi 
     leng(1)=0.
     do i=2,npts
        dr=xs(j,i)-xs(j,i-1)
        drt=zs(j,i)-zs(j,i-1)
        leng(i)=leng(i-1)+sqrt(dr**2+drt**2)
     end do
     dl=leng(npts)/(npts-1.)
     do i=1,npts
        leq(i)=(i-1)*dl
     end do
     leq(npts)=leng(npts)
     call spline1d(rl     ,leq,npts,xs(j,:)     ,leng,npts,work)
     call spline1d(zl     ,leq,npts,zs(j,:)     ,leng,npts,work)
     call spline1d(bpl    ,leq,npts,bps(j,:)    ,leng,npts,work)

     xs(j,:)=rl
     zs(j,:)=zl
     bps(j,:)=bpl

  enddo

   
  return

end subroutine equal_arc


