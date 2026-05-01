subroutine readgato

!-----------------------------------------------------------------------
! read dskgato created by toq or caltrans/teq equilibrium code
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
      character*8 line(10)
      integer j,i
      integer noutgato,ndskgato,ios
      integer old_npts
      integer neqtyp,nft,neqsym,ntht,it,ii,jj,nmap,isign
      character*8  etitl(8),date
      real, dimension(:), allocatable :: seqdpdr,seqdpdz, &
           psic_loc,dum2
      real, dimension(:,:), allocatable :: b_r,b_z
      real rcnt,xma,zma,btor,totcur,axddxz,dnnorm,rndoff
      integer ikr0,ikr1,ikz0,ikz1,jkm0,jkm1,jkp0,jkp1,iii,jjj,i0,i1
      real dpsi_loc,dchi
      real bnd_set,check,dum
      real convpsi,convpp,convf,convffp,convchi,convm,convb
      parameter (convpsi=1.e8,convpp=1.e-7,convf=1.e6,convffp=1.e4, &
          convchi=convpsi,convm=1.e2,convb=1.e4)

! assume file is a "new style" gato file for shape="gato" or shape="ngat"
!  new style has header and can be sym or asym.  old is sym with no header.
      neqtyp=0 ! "new style" equilibria
      if (shape.eq.'ogat') neqtyp=1  ! use shape="ogat" for old style
      nft=5
      nmap=1  ! reading inverse equilibrium
      rndoff=1.e-8

      noutgato=12
      ndskgato=22
      old_npts=npts

!      write(*,*) 'lrunname=',lrunname,'runname=',runname
      if (verbose .gt. 3) write(outfile,*) 'lrunname=',lrunname,'runname=',runname
      open(unit=ndskgato,file=runname(1:lrunname)//'.dskgato', &
           status='old',iostat=ios)
      if(ios==0) write(*,*) 'reading dskgato file ', &
           runname(1:lrunname)//'.dskgato'
      if(ios==0 .and. (verbose .gt. 3)) write(outfile,*) 'reading dskgato file ', &
           runname(1:lrunname)//'.dskgato'
      if(ios.ne.0) then
         write(6,*) 'problem opening file ',runname(1:lrunname)//'.dskgato'
         write(6,*) 'will try to read file named dskgato instead'
         open(unit=ndskgato,file='dskgato',status='old',iostat=ios)
         if(ios.ne.0) then
           write(*,*) 'could not open dskgato file'
           stop
         endif
         write(*,*) 'reading equilibrium from file dskgato'
         if (verbose .gt. 3) write(outfile,*) 'reading equilibrium from file dskgato'
      endif

!      read(ndskgato,'(10a8)') (line(i),i=1,10)

      if    (neqtyp .eq. 0) then
         read (ndskgato,1000) date
         read (ndskgato,1000) (etitl(i),i=1,nft)
         read (ndskgato,1010) npsi,ntht,neqsym
! 2.2 For old-style up-down symmetric equilibrium
!
      elseif(neqtyp .ne. 0) then
!
! 2.2.1 Set the default heading and date  (old style)
!
         date      = ' '
         do it  = 1,nft
            etitl(it) = ' '
         enddo
!10          continue
!
! 2.2.2 Read the dimensions                (old style)
!
! 2.2.2.1 Read npsi and ntht
        read (ndskgato,1020) npsi,ntht
!
! 2.2.2.2 Set default for up-down symmetric equilibrium
        neqsym  = 1
      endif

      if((neqsym .lt. 0) .or. (neqsym .gt.  1 )) then
         write(*,*) 'neqsym read from dskgato must be 0 or 1',neqsym
         stop
      endif

      if (npsi.le.0 .or. ntht.le.0) then
         write(*,*) 'npsi and ntht must be >0',npsi,ntht
         stop
      endif

      if(neqsym .eq. 0) npts= ntht
      if(neqsym .eq. 1) npts= 2*(ntht-1) + 1

      if (npts.ne.old_npts) then
         write(*,*) 'npts in input file must =npts from dskgato',npts,old_npts
         stop
      endif

! 4/3/00 for now the values of npsi & nthet from the dskbal file will 
!  override the input values, and nxinterp will be determined from
!  q (ie the number of surfaces used by elite will be the number from
!  the dskbal file which enclose the specified range of m numbers)

!      write(*,*) '!!!Elite will use the values of npts given in the', &
!      'dskbal file, overriding the value from the input file'

      write(*,*) 'npsi=',npsi,' npts=',npts,' neqsym=',neqsym
      if (verbose .gt. 3) write(outfile,*) &
           'npsi=',npsi,' npts=',npts,' neqsym=',neqsym
      if (verbose .gt. 4) write(outfile,*) 'allocating toq variables'

      call alloc_toq
      if (verbose .gt. 4) write(outfile,*) &
           'deallocate and reallocate surf arrays with new npts'
      call dealloc_arrays
      call alloc_arrays
      allocate( sp(npsi), seqdpdr(ntht), seqdpdz(ntht), &
           psic_loc(npsi), dum2(npsi) )

           
!  Read the equilibrium scalar data

      read (ndskgato,2000) rcnt,xma,zma,btor
      if(neqtyp .eq. 0) read (ndskgato,2000) totcur,axddxz
      if(neqtyp .eq. 1) read (ndskgato,2000) totcur,axddxz,dnnorm

      write(*,*) 'btor=',btor,' totcur=',totcur
      if (verbose .gt. 3) write(outfile,*) 'btor=',btor,' totcur=',totcur

! 3.1 Read the source functions
!      read (ndskgato,2000) (psimsh (jj), jj = 1,npsi)
      read (ndskgato,2000) (psiv (jj), jj = 1,npsi)
      read (ndskgato,2000) (fval     (jj), jj = 1,npsi)
      read (ndskgato,2000) (ffprime   (jj), jj = 1,npsi)
      read (ndskgato,2000) (sp     (jj), jj = 1,npsi)
      read (ndskgato,2000) (pprime    (jj), jj = 1,npsi)
      read (ndskgato,2000) (qsfin   (jj), jj = 1,npsi)
      read (ndskgato,2000) (nel   (jj), jj = 1,npsi)

! 3.2.1 Read the derivatives of the poloidal flux around the boundary
!         Note: the input may have this switched!

      read (ndskgato,2000) (seqdpdr(ii), ii = 1,ntht)
      read (ndskgato,2000) (seqdpdz(ii), ii = 1,ntht)

! 3.2.2 Read the inverse equilibrium (r,z)

      read (ndskgato,2000) ((xs(jj,ii), jj = 1,npsi), ii = 1,ntht)
      read (ndskgato,2000) ((zs(jj,ii), jj = 1,npsi), ii = 1,ntht)

! read Br and Bz arrays for shape eq ngat - better way to get bpol
      if (shape .eq. 'ngat') then
         write(*,*) 'Reading Br and Bz which should be present ', &
              'at end of t file for shape=ngat'
         allocate ( b_r(npsi,ntht), b_z(npsi,ntht) )
         read (ndskgato,2000) ((b_r(jj,ii), jj = 1,npsi), ii = 1,ntht)
         read (ndskgato,2000) ((b_z(jj,ii), jj = 1,npsi), ii = 1,ntht)
         bps=1.e4*sqrt(b_r**2+b_z**2)
         write(79,2000) ((bps(jj,ii),jj=1,npsi), ii=1,ntht)

      endif




! 4.0 Reverse signs if input has alternate convention

      if(psiv(1) .gt. psiv(npsi)) then

! 4.1 For nmap = 1 reverse signs of psiv and psi derivatives

! 4.1.1 Decide which signs to reverse: btor or totcur

        if(btor .lt. 0.0) isign  = -1
        if(btor .eq. 0.0) isign  =  0
        if(btor .gt. 0.0) isign  = +1

! 4.1.2 Print a warning if btor is zero

        if(isign .eq. 0) then
           write(*,*) 'WARNING!! btor=0 in readgato!!'
           if (verbose .gt. 3) write(outfile,*) 'WARNING!! btor=0 in readgato!!'
        endif

! 4.1.3 Reverse the signs

        if(isign .ne. 0) then
          if    (nmap .eq. 1) then
             write(*,*) 'reversing sign of btor or totcur',isign
             if (verbose .gt. 3) write(outfile,*) &
                  'reversing sign of btor or totcur',isign
            if(isign .gt. 0) totcur     = -totcur
            if(isign .lt. 0) btor       = -btor

            do jj  = 1,npsi
               psiv(jj) = -psiv(jj)
               ffprime(jj) = -ffprime(jj)
               pprime(jj) = -pprime(jj)
               if(isign .lt. 0 .and. fval(1).lt.0.) fval(jj)= -fval(jj)
            enddo


         endif
        endif
     endif




! 7.2 Fill the r and z arrays for up-down symmetry

      if(neqsym .eq. 1) then
!         write(*,*) 'filling in xs and zs for up-down symmetric case'
         if (verbose .gt. 4) write(outfile,*) &
              'filling in xs and zs for up-down symmetric case'
! 7.2.1 Ensure the magnetic axis is on the midplane

! 7.2.1.1 Check for discrepancy

        if(abs(zma) .gt. 0.0) then
          write(*,3000) xma,zma
          if(abs(zma).gt.rndoff) then
             write(*,*) 'WARNING: zma is not zero for neqdsym=1',zma
             if (verbose .gt. 3) write(outfile,*) &
                  'WARNING: zma is not zero for neqdsym=1',zma
          endif

! 7.2.1.2 Reset the axis on the midplane
        zma           = 0.0

! 7.3 Ensure the first flux surface is at the magnetic axis

! 7.3.1 Check for discrepancy
        ikr0          = 0
        ikr1          = 0
        ikz0          = 0
        ikz1          = 0

        do ii= 1,ntht
! 7.3.1.1 Check the r value is the axis value
           if(abs(xs(1,ii)-xma) .ne. 0.0) then
              if(abs(xs(1,ii)-xma) .le. rndoff) ikr0   = ikr0 + 1
              if(abs(xs(1,ii)-xma) .gt. rndoff) ikr1   = ikr1 + 1
           endif

! 7.3.1.2 Check the z value is the axis value
           if(abs(zs(1,ii)-zma) .ne. 0.0) then
              if(abs(zs(1,ii)-zma) .le. rndoff) ikz0   = ikz0 + 1
              if(abs(zs(1,ii)-zma) .gt. rndoff) ikz1   = ikz1 + 1
           endif
        enddo

! 7.3.2 Print error message if any discrepancy

        if(ikr0 .gt. 0  .or.  ikr1 .gt. 0  .or.  ikz0 .gt. 0  .or. &
             ikz1 .gt. 0) then
          write(*,3100) npsi,ntht
          write(*,3110) (iii, xs(1,iii), zs(1,iii), &
               iii=1,ntht)
          write(*,*) 'Discrepancy in r or z read in',ikr0,ikz0, &
                  ikr1,ikz1
          if (verbose .gt. 3) write(outfile,*) &
               'Discrepancy in r or z read in',ikr0,ikz0, &
                  ikr1,ikz1
          endif
       endif

! 7.3.3 Reset the first flux surface identically to the axis values
       
       do ii= 1,ntht
          xs(1,ii) = xma
          zs(1,ii) = zma
       enddo

! 7.4 Ensure that the first and last ray lay along the midplane

! 7.4.1 Check for discrepancy

        jkm0          = 0
        jkm1          = 0
        jkp0          = 0
        jkp1          = 0

        do jj= 1,npsi
! 7.4.1.1 Check the first ray is on the midplane
           if(abs(zs(jj,  1 )) .ne.   0.0 ) then
              if(abs(zs(jj,  1 )) .le. rndoff) jkp0   = jkp0 + 1
              if(abs(zs(jj,  1 )) .gt. rndoff) jkp1   = jkp1 + 1
           endif

! 7.4.1.2 Check the last  ray is on the midplane
           if(abs(zs(jj,ntht)) .ne.   0.0 ) then
              if(abs(zs(jj,ntht)) .le. rndoff) jkm0   = jkm0 + 1
              if(abs(zs(jj,ntht)) .gt. rndoff) jkm1   = jkm1 + 1
           endif
        enddo

! 7.4.2 Print error message if any discrepancy

        if(jkm0 .gt. 0  .or.  jkp0 .gt. 0  .or.  jkm1 .gt. 0  .or. &
             jkp1 .gt. 0) then
           write(*,3200) npsi,ntht
           write(*,3210) (jjj, xs(jjj,  1 ), zs(jjj,  1 ), &
                xs(jjj,ntht), zs(jjj,ntht), jjj=1,npsi)
           if (verbose .gt. 3) write(outfile,3200) npsi,ntht
           if (verbose .gt. 3) write(outfile,3210) (jjj, xs(jjj,  1 ), zs(jjj,  1 ), &
                xs(jjj,ntht), zs(jjj,ntht), jjj=1,npsi)
           write(*,*) 'first or last ray slightly off midplane', &
                jkm0,jkp0,jkm1,jkp1
           if (verbose .gt. 3) write(outfile,*) &
                'first or last ray slightly off midplane', &
                jkm0,jkp0,jkm1,jkp1
        endif

! 7.4.3 Reset the rays along the midplane

        do jj= 1,npsi
           zs(jj,  1 ) = 0.0
           zs(jj,ntht) = 0.0
        enddo

! 7.5 Fill the complete r and z arrays

        do ii= ntht+1,npts
           i1= ii
           i0= npts + 1 - i1
           do jj= 1,npsi
              xs(jj,i1) = +xs(jj,i0)
              zs(jj,i1) = -zs(jj,i0)
           enddo
        enddo
     endif


! Need to calculate chipsi= d chi/d psi in Bob's toq notation
!      d psiv/d psic_loc
     dpsi_loc=1./(npsi-1.)
     dchi=psiv(npsi)-psiv(1)
     if (verbose .gt. 3) write(outfile,*) 'dpsi_loc=',dpsi_loc,' dchi=',dchi
     do ii=1,npsi
        psic_loc(ii)=(ii-1.)*dpsi_loc
     enddo

      bnd_set=-2.d30  ! calculate boundary condition from nearest 4 pts
      call spline(psic_loc,psiv,npsi,bnd_set,bnd_set,dum2)
        ! dum2 contains the interpolating function
      do i=1,npsi
         call zsplint(psic_loc,psiv,dum2,dum2,npsi,psic_loc(i),check, &
              chipsi(i),dum)
      enddo
! enforce positive chipsi
      do i=1,npsi
         if (chipsi(i)<0.) then
!  print warning only if chipsi is non-negligibly positive and
!   if it is used (for shape.eq.ngat br and bz are read directly to get bp)
            if ((abs(chipsi(i))>1.e-12).and.(shape.ne.'ngat')) then
               write(*,*) 'WARNING!! found negative chipsi=',chipsi(i), &
                    ' for surface=' ,i,' reset to 0'
               if (verbose .gt. 3) write(outfile,*) &
                    'WARNING!! found negative chipsi=',chipsi(i), &
                    ' for surface=',i,' reset to 0'
            else
               if (verbose .gt. 3) write(outfile,*) &
                    'found negative chipsi=',chipsi(i), &
                    ' for surface=',i,' reset to 0'
            endif
            chipsi(i)=0.
         endif
      enddo

! Convert everything to CGS units for use by ELITE
      do i=1,npsi
         psiv(i)=convpsi*psiv(i)
         sp(i)=convpp*sp(i)*convpsi
         pprime(i)=convpp*pprime(i)
         fval(i)=convf*fval(i)
         ffprime(i)=convffp*ffprime(i)
         chipsi(i)=convchi*chipsi(i)
      enddo

      pressval=sp

      do j=1,npts
         do i=1,npsi
            xs(i,j)=convm*xs(i,j)
            zs(i,j)=convm*zs(i,j)
         enddo
      enddo

     write(*,*) 'npsi=',npsi,'npts=',npts
     if (verbose .gt. 3) then
        write(outfile,*) 'npsi=',npsi,'npsts=',npts
        write(outfile,*) 'q=',qsfin
        write(outfile,*) 'pprime=',pprime
        write(outfile,*) 'chipsi=',chipsi
     endif


     close(ndskgato)
!     stop


      return

 1000 format(6a8)
 1010 format(3i5)
 1020 format(2i5)
 2000 format(1p4e19.12)
 3000 format(/,1x,'Warning: the axis is not on the midplane',/ &
           ,1x,'xma    = ',e16.9,4x,'zma    = ',e16.9)
 3100 format(/,1x,'Warning: the first flux surface does not correspond' &
           ,1x,'to the magnetic axis',/ &
           ,10x,'npsi = ',i5,2x,'ntht = ',i5,// &
           ,4x,'ii',6x,'seqrps(1,ii)',5x,'seqzps(1,ii)')
 3110 format((1x,i5,2x,2(1x,e16.9)))
 3200 format(/,1x,'Warning: the first ray is not on the midplane' &
             ,1x,'npsi = ',i5,2x,'ntht = ',i5,// &
             ,4x,'jj',6x,'seqrps(jj,  1 )',2x,'seqzps(jj,  1 )' &
             ,2x,'seqrps(jj,ntht)',2x,'seqzps(jj,ntht)')
 3210 format((1x,i5,2x,4(1x,e16.9)))

end subroutine readgato
