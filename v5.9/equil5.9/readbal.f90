
subroutine readbal

!-----------------------------------------------------------------------
! read dskbal created by toq equilibrium code
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
      character*8 line(10)
      integer j,i
      integer noutbal,ndskbal,ios
      integer old_npts
      real dum,bndset
      real, dimension(:), allocatable :: za,pspline

      noutbal=12
      ndskbal=22
      old_npts=npts

!      open(unit=noutbal,file=runname(1:lrunname)//'.outbal',status='unknown')

      if (verbose .ge. 2) write(*,*) 'lrunname=',lrunname,'runname=',runname
      if (verbose .gt. 3) write(outfile,*) 'lrunname=',lrunname,'runname=',runname
      open(unit=ndskbal,file=runname(1:lrunname)//'.dskbal', &
           status='old',iostat=ios)
      if(ios==0) write(*,*) 'reading dskbal file ', &
           runname(1:lrunname)//'.dskbal'
      if(ios.ne.0) then
         if (verbose .ge. 1) write(6,*) &
              'problem opening file ',runname(1:lrunname)//'.dskbal'
         if (verbose .ge. 1) write(6,*) &
              'will try to read file named dskbal instead'
         open(unit=ndskbal,file='dskbal',status='old',iostat=ios)
         if(ios.ne.0) then
           write(*,*) 'could not open dskbal file'
           stop
         endif
         if (verbose .ge. 1) write(*,*) 'reading equilibrium from file dskbal'
         if (verbose .gt. 3) write(outfile,*) 'reading equilibrium from file dskbal'
      endif

!      open (unit=ndskbal,file=filename,status='old')
! first read namelist quantities at the head of dskbal
!      write(noutbal,50)
 50   format(//,'*** following is the namelist from the equilibrium' &
              ,' code ***')
!      call header(ndskbal,noutbal)
!      write(6,*)' in header'
      call header(ndskbal,0)  ! 0 produces no outbal file
      read(ndskbal,'(10a8)') (line(i),i=1,10)
      read(ndskbal,100) npsi,npts

      if (npts.ne.old_npts) then
         write(*,*) 'npts in input file must =npts from toq',npts,old_npts
         stop
      endif

! 4/3/00 for now the values of npsi & nthet from the dskbal file will 
!  override the input values, and nxinterp will be determined from
!  q (ie the number of surfaces used by elite will be the number from
!  the dskbal file which enclose the specified range of m numbers)

!      write(*,*) '!!!Elite will use the values of npts given in the', &
!      'dskbal file, overriding the value from the input file'

      if (verbose .ge. 1) write(*,*) 'npsi=',npsi,' npts=',npts
      if (verbose .gt. 4) write(outfile,*) 'allocating toq variables'

      call alloc_toq
      if (verbose .gt. 4) write(outfile,*) &
           'deallocate and reallocate surf arrays with new npts'
      call dealloc_arrays
      call alloc_arrays
      allocate( pspline(npsi),za(npsi) )

 100  format(2i5)
      read(ndskbal,'(10a8)') (line(i),i=1,10)
!      write(6,*) line
      do i=1,npsi/5
         read(ndskbal,200) psiv((i*5)-4:i*5)
!         write(*,*)' i, psi=', i,psiv((i*5)-4:i*5)
      enddo
!      if (mod(npsi,5).ne.0) read(ndskbal,200) psiv((npsi/5)*5+1:npsi)
      read(ndskbal,200) psiv((npsi/5)*5+1:npsi)
!      write(noutbal,*) psiv
!      write(6,*)' got psiv'
      read(ndskbal,'(10a8)') (line(i),i=1,10)
!      write(6,*) line
      do i=1,npsi/5
         read(ndskbal,200) pprime((i*5)-4:i*5)
!         write(*,*) i,pprime((i*5)-4:i*5)
      enddo
!      if (mod(npsi,5).ne.0) read(ndskbal,200) pprime((npsi/5)*5+1:npsi)
      read(ndskbal,200) pprime((npsi/5)*5+1:npsi)
!      write(noutbal,*) pprime

      read(ndskbal,'(10a8)') (line(i),i=1,10)
!      write(noutbal,*) line
      do i=1,npsi/5
         read(ndskbal,200) fval((i*5)-4:i*5)
!         write(*,*) i,fval((i*5)-4:i*5)
      enddo
!      if (mod(npsi,5).ne.0) read(ndskbal,200) fval((npsi/5)*5+1:npsi)
      read(ndskbal,200) fval((npsi/5)*5+1:npsi)
!      write(noutbal,*) fval

      read(ndskbal,'(10a8)') (line(i),i=1,10)
!      write(noutbal,*) line
      do i=1,npsi/5
         read(ndskbal,200) ffprime((i*5)-4:i*5)
!         write(*,*) i,ffprime((i*5)-4:i*5)
      enddo
!      if (mod(npsi,5).ne.0) read(ndskbal,200) ffprime((npsi/5)*5+1:npsi)
      read(ndskbal,200) ffprime((npsi/5)*5+1:npsi)
!      write(noutbal,*) ffprime

      read(ndskbal,'(10a8)') (line(i),i=1,10)
!      write(noutbal,*) line
      do i=1,npsi/5
         read(ndskbal,200) chipsi((i*5)-4:i*5)
!         write(*,*) i,chipsi((i*5)-4:i*5)
      enddo
!      if (mod(npsi,5).ne.0) read(ndskbal,200) chipsi((npsi/5)*5+1:npsi)
      read(ndskbal,200) chipsi((npsi/5)*5+1:npsi)
!      write(noutbal,*) chipsi

!      write(6,*)' reading q'
      read(ndskbal,'(10a8)') (line(i),i=1,10)
!      write(noutbal,*) line
      do i=1,npsi/5
         read(ndskbal,200) qsfin((i*5)-4:i*5)
!         write(*,*) i,qsfin((i*5)-4:i*5)
      enddo
!      if (mod(npsi,5).ne.0) read(ndskbal,200) qsfin((npsi/5)*5+1:npsi)
      read(ndskbal,200) qsfin((npsi/5)*5+1:npsi)
!      write(noutbal,*) qsfin
!      write(6,*)' done q'

      read(ndskbal,'(10a8)') (line(i),i=1,10)
!      write(noutbal,*) line
      do j=1,npts
!      write(*,*) j,npts
         do i=1,npsi/5
            read(ndskbal,200) xs((i*5)-4:i*5,j)
!            write(*,*) i,xs((i*5)-4:i*5,j)
         enddo
         read(ndskbal,200) xs((npsi/5)*5+1:npsi,j)
!         if (mod(npsi,5).ne.0) read(ndskbal,200) xs((npsi/5)*5+1:npsi,j)
      enddo
!      write(noutbal,*) xs

      read(ndskbal,'(10a8)') (line(i),i=1,10)
!      write(noutbal,*) line
      do j=1,npts
!      write(*,*) j,npts
         do i=1,npsi/5
            read(ndskbal,200) zs((i*5)-4:i*5,j)
!            write(*,*) i,zs((i*5)-4:i*5,j)
         enddo
!         if (mod(npsi,5).ne.0) read(ndskbal,200) zs((npsi/5)*5+1:npsi,j)
         read(ndskbal,200) zs((npsi/5)*5+1:npsi,j)
      enddo
!      write(noutbal,*) zs

200   format(5g20.12)
!
!      write(6,*)' dens=',dens
!  If dens is set to .true. then read density, Te and Ti profiles
      if (dens) then
        if (verbose .ge. 2) write(6,*)' reading density'
        read(ndskbal,'(10a8)') (line(i),i=1,10) 
        do i=1,npsi/5
           read(ndskbal,200) nel((i*5)-4:i*5)
           if (verbose .gt. 4) write(outfile,*) i,nel((i*5)-4:i*5)
        enddo
!        if (mod(npsi,5).ne.0) read(ndskbal,200) nel((npsi/5)*5+1:npsi)
        read(ndskbal,200) nel((npsi/5)*5+1:npsi)
!        write(noutbal,*) nel
!        write(6,*)' done density'
!
        read(ndskbal,'(10a8)') (line(i),i=1,10)
        if (verbose .gt. 4) write(outfile,*)' reading n-prime'
        do i=1,npsi/5
           read(ndskbal,200) nprime((i*5)-4:i*5)
           if (verbose .gt. 4) write(outfile,*) i,nprime((i*5)-4:i*5)
        enddo
!        if (mod(npsi,5).ne.0) read(ndskbal,200) nprime((npsi/5)*5+1:npsi)
        read(ndskbal,200) nprime((npsi/5)*5+1:npsi)
!        write(noutbal,*) nprime
!        write(6,*)' done n-prime'
!
        read(ndskbal,'(10a8)') (line(i),i=1,10)
        do i=1,npsi/5
           read(ndskbal,200) tel((i*5)-4:i*5)
!           write(*,*) i,tel((i*5)-4:i*5)
        enddo
!        if (mod(npsi,5).ne.0) read(ndskbal,200) tel((npsi/5)*5+1:npsi)
        read(ndskbal,200) tel((npsi/5)*5+1:npsi)
!        write(noutbal,*) tel
!
        read(ndskbal,'(10a8)') (line(i),i=1,10)
        do i=1,npsi/5
           read(ndskbal,200) teprime((i*5)-4:i*5)
!           write(*,*) i,teprime((i*5)-4:i*5)
        enddo
!        if (mod(npsi,5).ne.0) read(ndskbal,200) teprime((npsi/5)*5+1:npsi)
        read(ndskbal,200) teprime((npsi/5)*5+1:npsi)
!        write(noutbal,*) teprime
!
        read(ndskbal,'(10a8)') (line(i),i=1,10)
        do i=1,npsi/5
           read(ndskbal,200) tion((i*5)-4:i*5)
           if (verbose .gt. 4) write(outfile,*) i,tion((i*5)-4:i*5)
        enddo
!        if (mod(npsi,5).ne.0) read(ndskbal,200) tion((npsi/5)*5+1:npsi)
        read(ndskbal,200) tion((npsi/5)*5+1:npsi)
!        write(noutbal,*) tion
!
        read(ndskbal,'(10a8)') (line(i),i=1,10)
        do i=1,npsi/5
           read(ndskbal,200) tiprime((i*5)-4:i*5)
           if (verbose .gt. 4) write(outfile,*) i,tiprime((i*5)-4:i*5)
        enddo
!        if (mod(npsi,5).ne.0) read(ndskbal,200) tiprime((npsi/5)*5+1:npsi)
        read(ndskbal,200) tiprime((npsi/5)*5+1:npsi)
!        write(noutbal,*) tiprime
      end if

      if (verbose .gt. 4) write(outfile,*) 'finished reading dskbal file'

! TEMPORARY CHANGE!!!
      garbage=0
      if (garbage.gt.0) then
         npsi=npsi-garbage
         write(*,*) 'reducing npsi by ',garbage,'!!!!!!'
         if (verbose .gt. 3) write(outfile,*) 'reducing npsi by ',garbage,'!!!!!!'
      endif

      close(ndskbal)
!      close(noutbal)

!-----------------------------------------------------------------------
! calculate pressure
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
               if (verbose .ge. 1) write(*,*) &
                    'warning, splined pressure<0 at i=',i, &
                    ' resetting to zero'
               pressval(i)=0.
            endif
         enddo
      endif

!      call spline(psiv,pprime,npsi,-1.e30,-1.e30,ppp)

! remove as it causes problems on DECs for some reason 5/00
!      do i=1,npsi
!         call zspline(psiv,pprime,ppp,npsi,0.,npsi,press(i)) 
!      end do

   
  return
      
end subroutine readbal


      
subroutine header(nunin,nunout)
!-----------------------------------------------------------------------
!     read and write namelist header on dskxxx
!-----------------------------------------------------------------------
      implicit none
      integer nunin,nunout
      integer i,j
      character*8 ilab(10)
      character*5 iend
      iend=" &end"
      rewind(nunin)
      do  i=1,3
!         write(6,*)' i=',i
         read(nunin,150) (ilab(j),j=1,10)
!         write(6,*) (ilab(j),j=1,10)
 150     format(10a8)
         if(nunout.ne.0) write(nunout,150) (ilab(j),j=1,10)
         if(ilab(1)(1:5).eq.iend) return
      end do
!      write(6,*)ilab(1)(1:5)
!      write(6,*)' not returned!'
     
 300  read(nunin,150) (ilab(j),j=1,10)
      if(nunout.ne.0) write(nunout,150) (ilab(j),j=1,10)
      if(ilab(1)(1:5).eq.iend) return
      go to 300
      
end subroutine header



subroutine getbps
!-----------------------------------------------------------------------
!     calculate bpoloidal if input is dskbal
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
      integer i,j,jp,jm
      real xthe,zthe,xpsi,zpsi
      real xtheta(npts),ztheta(npts)
      real denom,psix,psiz,gpsi2
      real bnd_set,check,dum,dum2(npsi)
      real xmin,zmin,xmax,zmax
      dthe=2.*pi/(npts-1.)
      dpsi=1./(npsi-1.)
!      dpsi=1./(npsi+garbage-1.)

! calculate bps for equilibrium file types which do not
!  already contain it
      if ( (shape.ne.'eqbm') .and. (shape.ne.'ngat') ) then

         do  i=2,npsi
            do  j=1,npts
               if(i.ne.npsi) then
                  xpsi=(xs(i+1,j)-xs(i-1,j))/(2.*dpsi)
                  zpsi=(zs(i+1,j)-zs(i-1,j))/(2.*dpsi)
               else
                  xpsi=2.*(0.75*xs(npsi,j)-xs(npsi-1,j)+ &
                       0.25*xs(npsi-2,j))/dpsi
                  zpsi=2.*(0.75*zs(npsi,j)-zs(npsi-1,j)+ &
                       0.25*zs(npsi-2,j))/dpsi
               endif
               jp=j+1
               jm=j-1
               if(jm.eq.0) jm=npts-1
               if(jp.eq.npts+1) jp=2
               xthe=(xs(i,jp)-xs(i,jm))/(2.*dthe)
               zthe=(zs(i,jp)-zs(i,jm))/(2.*dthe)
               denom=1./(xpsi*zthe-zpsi*xthe)
               psix=zthe*denom
               psiz=-xthe*denom
               gpsi2=psix**2+psiz**2
               bps(i,j)=sqrt(psix**2+psiz**2)*chipsi(i)/xs(i,j)
! correct for garbage points if there
!            bps(i,j)=bps(i,j)*psiv(npsi)/psiv(npsi+garbage)
               xtheta(j)=xthe
               ztheta(j)=zthe
            end do
            arcsur(i,1)=0.
            do j=2,npts
               arcsur(i,j)=arcsur(i,j-1)+0.5*( &
                    sqrt(xtheta(j-1)**2+ztheta(j-1)**2)+ &
                    sqrt(xtheta(j  )**2+ztheta(j  )**2))*dthe
            end do
         end do
      endif

! calculate approximate values for xaxis and zaxis since these
!  aren't read in
      xmax=xs(1,1)
      xmin=xmax
      zmax=zs(1,1)
      zmin=zmax
      do j=2,npts
         if (xs(1,j) > xmax) xmax=xs(1,j)
         if (xs(1,j) < xmin) xmin=xs(1,j)
         if (zs(1,j) > zmax) zmax=zs(1,j)
         if (zs(1,j) < zmin) zmin=zs(1,j)
      enddo
      zaxis=0.5*(zmax+zmin)
      xaxis=0.5*(xmax+xmin)
      if (verbose .ge. 2) write(*,*) 'inner surface, xmax=',xmax,'xmin=',xmin, &
           'zmax=',zmax,'zmin=',zmin
      if (verbose .ge. 2) write(*,*) &
           'approximate axis location, xaxis=',xaxis,'zaxis=',zaxis
      if (verbose .gt. 3) write(outfile,*) &
           'inner surface, xmax=',xmax,'xmin=',xmin, &
           'zmax=',zmax,'zmin=',zmin
      if (verbose .gt. 3) write(outfile,*) &
           'approximate axis location, xaxis=',xaxis,'zaxis=',zaxis

! 10/00 calculate additional equilibrium quantities q',q'',p'',(ff')'
      bnd_set=-2.d30  ! calculate boundary condition from nearest 4 pts
!      write(*,*) 'psiv=',psiv,'qsfin=',qsfin,'npsi=',npsi
      call spline(psiv,qsfin,npsi,bnd_set,bnd_set,dum2)
        ! dum2 contains the interpolating function
!      write(*,*) 'called spline dum2=',dum2
      do i=1,npsi
         call zsplint(psiv,qsfin,dum2,dum2,npsi,psiv(i),check,qsfinp(i),dum)
!         write(*,*) 'i=',i,' qsfin=',qsfin(i),' check=',check, &
!              'qsfinp=',qsfinp(i)
      enddo

      call spline(psiv,qsfinp,npsi,bnd_set,bnd_set,dum2)
        ! dum2 contains the interpolating function
      do i=1,npsi
         call zsplint(psiv,qsfinp,dum2,dum2,npsi,psiv(i),check,qsfinpp(i),dum)
!         write(*,*) 'i=',i,' qsfinp=',qsfinp(i),' check=',check, &
!              'qsfinpp=',qsfinpp(i)
      enddo

!! eqbm format already contains ppp and ffpp 7.04
      if (shape .ne. 'eqbm') then
         call spline(psiv,pprime,npsi,bnd_set,bnd_set,dum2)
         ! dum2 contains the interpolating function
         do i=1,npsi
            call zsplint(psiv,pprime,dum2,dum2,npsi,psiv(i),check,ppp(i),dum)
!         write(*,*) 'i=',i,' pprime=',pprime(i),' check=',check, &
!              'ppp=',ppp(i)
         enddo

         call spline(psiv,ffprime,npsi,bnd_set,bnd_set,dum2)
        ! dum2 contains the interpolating function
         do i=1,npsi
            call zsplint(psiv,ffprime,dum2,dum2,npsi,psiv(i),check,ffpp(i),dum)
!         write(*,*) 'i=',i,' ffprime=',ffprime(i),' check=',check, &
!              'ffpp=',ffpp(i)
         enddo
      endif
 

      return

end subroutine getbps




