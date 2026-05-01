subroutine writegato

!-----------------------------------------------------------------------
! write dskgato file out, for use in GATO for benchmarking cases
!  where only eqin or dskbal file is available
! 5/05 pbs: to start only use old dskgato format (neqtyp=1) which
!  supports only updown symmetric and writes only upper half of
!  the equilibrium
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none

      integer noutgato,ndskgato,ngato,ios,isurf,j,ntht,jj,ii
      real r4,r2,rz,z4,z2,rsq,zsq,determ,psirr,psizz
      real totcur,axddxz,dnnorm,rcnt
      real convpsi,convpp,convf,convffp,convchi,convm,convb,convp
      character*24 cdatetime,ctime
      integer time,isym
      character*64 alab
      parameter (convpsi=1.e-8,convpp=1.e7,convf=1.e-6,convffp=1.e-4, &
          convchi=convpsi,convm=1.e-2,convb=1.e-4,convp=1.e-1)

      noutgato=12
      ndskgato=22
      ngato=28

      if (verbose .gt. 3) write(outfile,*) &
           'writing ',runname(1:lrunname)//'_out.dskgato'
      open(unit=ndskgato,file=runname(1:lrunname)//'_out.dskgato', &
           status='unknown')
      if (verbose .gt. 3) write(outfile,*) &
           'writing ',runname(1:lrunname)//'_out.ndskgato'
      open(unit=ngato,file=runname(1:lrunname)//'_out.ndskgato', &
           status='unknown')


      ntht=(npts-1)/2+1
      totcur=1.5e6  ! need to calculate this
      dnnorm=1.
      rcnt=0.5*(xs(npsi,1)+xs(npsi,ntht))
      cdatetime =  ctime (time ())
      alab = 'dskgato file from elite, '//runname(1:lrunname)//' '//shape

! calculate elongation on axis
      isurf=2
      r4=0.
      r2=0.
      rz=0.
      z4=0.
      z2=0.
      do j=1,ntht
         rsq=(xs(isurf,j)-xs(1,1))**2
         zsq=zs(isurf,j)**2
         r4=r4+rsq*rsq
         r2=r2+rsq
         rz=rz+rsq*zsq
         z2=z2+zsq
         z4=z4+zsq*zsq
      end do
      determ=r4*z4-rz*rz
      psirr=(r2*z4-z2*rz)/determ
      psizz=(r4*z2-rz*r2)/determ
      axddxz=sqrt(psirr/psizz)



      write(ndskgato,1020) npsi,ntht

      write(ndskgato,2000) rcnt*convm,xs(1,1)*convm, &
           zs(1,1)*convm,(fval(npsi)/rcnt)*convb
      write(ndskgato,2000) totcur,axddxz,dnnorm
      write(ndskgato,2000) (psiv(jj)*convpsi, jj = 1,npsi)
      write(ndskgato,2000) (fval(jj)*convf, jj = 1,npsi)
      write(ndskgato,2000) (ffprime(jj)*convffp, jj = 1,npsi)
      write(ndskgato,2000) (pressval(jj)*convp, jj = 1,npsi)
      write(ndskgato,2000) (pprime(jj)*convpp, jj = 1,npsi)
      write(ndskgato,2000) (qsfin(jj), jj = 1,npsi)
      if (nel(npsi/2) > 0.0) then
         write(ndskgato,2000) (nel(jj), jj = 1,npsi)
      else
         write(ndskgato,2000) (dnnorm, jj = 1,npsi)
      endif

! write zeros for seqdpdr and swqdpdz for now, not used in gato
!      read (ndskgato,2000) (seqdpdr(ii), ii = 1,ntht)
!      read (ndskgato,2000) (seqdpdz(ii), ii = 1,ntht)
      write(ndskgato,2000) (0., ii = 1,ntht)
      write(ndskgato,2000) (0., ii = 1,ntht)

      write(ndskgato,2000) ((xs(jj,ii)*convm, jj = 1,npsi), ii = 1,ntht)
      write(ndskgato,2000) ((zs(jj,ii)*convm, jj = 1,npsi), ii = 1,ntht)

      close(ndskgato)


! now write ngato format (whole equilibrium)

! recalculate elongation on axis
      isurf=2
      r4=0.
      r2=0.
      rz=0.
      z4=0.
      z2=0.
      do j=1,npts
         rsq=(xs(isurf,j)-xs(1,1))**2
         zsq=zs(isurf,j)**2
         r4=r4+rsq*rsq
         r2=r2+rsq
         rz=rz+rsq*zsq
         z2=z2+zsq
         z4=z4+zsq*zsq
      end do
      determ=r4*z4-rz*rz
      psirr=(r2*z4-z2*rz)/determ
      psizz=(r4*z2-rz*r2)/determ
      axddxz=sqrt(psirr/psizz)

      isym=0

      write(ngato,*) cdatetime
      write(ngato,*) alab
      write(ngato,1010) npsi,npts,isym
      write(ngato,2000) rcnt*convm,xs(1,1)*convm, &
           zs(1,1)*convm,(fval(npsi)/rcnt)*convb
      write(ngato,2000) totcur,axddxz,dnnorm
      write(ngato,2000) (psiv(jj)*convpsi, jj = 1,npsi)
      write(ngato,2000) (fval(jj)*convf, jj = 1,npsi)
      write(ngato,2000) (ffprime(jj)*convffp, jj = 1,npsi)
      write(ngato,2000) (pressval(jj)*convp, jj = 1,npsi)
      write(ngato,2000) (pprime(jj)*convpp, jj = 1,npsi)
      write(ngato,2000) (qsfin(jj), jj = 1,npsi)
      if (nel(npsi/2) > 0.0) then
         write(ngato,2000) (nel(jj), jj = 1,npsi)
      else
         write(ngato,2000) (dnnorm, jj = 1,npsi)
      endif

! write zeros for seqdpdr and swqdpdz for now, not used in gato
!      read (ngato,2000) (seqdpdr(ii), ii = 1,npts)
!      read (ngato,2000) (seqdpdz(ii), ii = 1,npts)
      write(ngato,2000) (0., ii = 1,npts)
      write(ngato,2000) (0., ii = 1,npts)

      write(ngato,2000) ((xs(jj,ii)*convm, jj = 1,npsi), ii = 1,npts)
      write(ngato,2000) ((zs(jj,ii)*convm, jj = 1,npsi), ii = 1,npts)

      close(ngato)

      return


 1000 format(6a8)
 1010 format(3i5)
 1020 format(2i5)
 2000 format(1p4e19.12)


end subroutine writegato
