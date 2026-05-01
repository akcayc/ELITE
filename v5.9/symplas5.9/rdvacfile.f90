
subroutine rdvacfile

!  translated to F90  1/31/00
!-----------------------------------------------------------------------
! read dskvac_edge file to get deltaW_vac matrix
!  modified 7/22 to allow vac to include more modes than needed by plas
!-----------------------------------------------------------------------
      use elite_data, only: nmodes,runname,lrunname,outfile,mmax,mmin, &
           nedge,dw_vac,offset_v,mindex,verbose
      implicit none
      integer nunit,nmodes_v,mmax_v,mmin_v
      real gamr(nmodes),vr(nmodes,nmodes)
      real dw_temp(nmodes,nmodes),bmat(nmodes,nmodes)
      integer ios,i,j,m,mp
      nunit=31
      open(unit=nunit,file=runname(1:lrunname)//'.vac', &
           status='old',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem opening ',runname(1:lrunname)//'.vac'
         stop
      endif
      if (verbose .ge. 4) write(outfile,*) 'reading ',runname(1:lrunname)//'.vac'
      read(nunit,*) nmodes_v,mmax_v,mmin_v
      if (verbose .ge. 3) write(outfile,*) &
           'nmodes_v=',nmodes_v,' mmax_v=',mmax_v,' mmin_v=',mmin_v
      if(mmax_v.ne.mmax) then
         write(6,*) 'mmax_vac.ne.mmax_plas',mmax_v,mmax
         stop
      endif
      if(mmin_v.ne.mmin) then
         write(6,*) 'mmin_vac.ne.mmin_plas',mmin_v,mmin
         stop
      endif
      if (verbose .ge. 3) write(outfile,*) &
           'nmodes= ',nmodes,' nmodes_v= ',nmodes_v,' nedge=',nedge
      read(nunit,*) ((dw_vac(i,j),i=1,nmodes_v),j=1,nmodes_v)
!-----------------------------------------------------------------------
! need offset if mmax_v ne. mmax
!-----------------------------------------------------------------------
      offset_v=mmax_v-mmax
!-----------------------------------------------------------------------
! dw_vac=Integral[phi* n.Grad[phi] ds]
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! verify that dw_vac is positive definite matrix
!-----------------------------------------------------------------------
!     do m=1,nedge
!        do mp=1,nedge
      do m=1,nedge/mindex
         do mp=1,nedge/mindex
            dw_temp(m,mp)=dw_vac(m,mp)
            bmat(m,mp)=0.
            vr(m,mp)=0.
         end do
      end do
      if (verbose .ge. 5) write(outfile,*) 'call eigen in rdvacfile'
      gamr=0.
      call eigen(dw_temp,nmodes,nedge/mindex,gamr,vr,bmat)
!      call eigen(dw_temp,nmodes_v,nmodes_v,gamr,vr,bmat)
      do i=1,nmodes_v
         if(gamr(i).lt.0.0) then
            write(6,*) 'error dw_vac is not positive definite'
            write(6,*) gamr(i)
            stop
         endif
      end do
    return
      
end subroutine rdvacfile
