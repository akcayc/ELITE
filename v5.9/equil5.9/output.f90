      
subroutine wrteqdat

!-----------------------------------------------------------------------
!     to test mapper's output
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
      integer i,j,ndskeqdat,ndskpol

      ndskeqdat=97
      if (verbose .gt. 4) write(outfile,*) 'writing ',runname(1:lrunname)//'.eqdat'
      open(unit=ndskeqdat,file=runname(1:lrunname)//'.eqdat', &
           status='unknown')

! For now, rewrite to write only those surfaces used by elite,
!   in reverse order (out to in) the way elite uses them  2/7/00

      write(ndskeqdat,'(" npsi used by elite, npts")')
      write(ndskeqdat,'(2i5)') nxinterp,npts
      write(ndskeqdat,'("poloidal flux--psiv")')
!      write(ndskeqdat,200) (psiv(i),i=npsi,npsi-nxinterp+1,-1)
!  write flux with axis defined to be zero  PS 2/11/00
      write(ndskeqdat,200) (psiv(i)-psiv(1),i=npsi,npsi-nxinterp+1,-1)
      write(ndskeqdat,'("pprime")')
      write(ndskeqdat,200) (pprime(i),i=npsi,npsi-nxinterp+1,-1)
      write(ndskeqdat,'("f")')
      write(ndskeqdat,200) (fval(i),i=npsi,npsi-nxinterp+1,-1)
      write(ndskeqdat,'("ffprime")')
      write(ndskeqdat,200) (ffprime(i),i=npsi,npsi-nxinterp+1,-1)
      write(ndskeqdat,'("vprime")')
      write(ndskeqdat,200) (chipsi(i),i=npsi,npsi-nxinterp+1,-1)
      write(ndskeqdat,'("q")')
      write(ndskeqdat,200) (qsfin(i),i=npsi,npsi-nxinterp+1,-1)
      write(ndskeqdat,'("ppp=pprimeprime")')
      write(ndskeqdat,200) (ppp(i),i=npsi,npsi-nxinterp+1,-1)
      write(ndskeqdat,'("ffpp=(ffprime)prime")')
      write(ndskeqdat,200) (ffpp(i),i=npsi,npsi-nxinterp+1,-1)
      write(ndskeqdat,'("qprime")')
      write(ndskeqdat,200) (qsfinp(i),i=npsi,npsi-nxinterp+1,-1)
      write(ndskeqdat,'("qprimeprime")')
      write(ndskeqdat,200) (qsfinpp(i),i=npsi,npsi-nxinterp+1,-1)
      if (dens) then
        write(ndskeqdat,'("electron density")')
        write(ndskeqdat,200) (nel(i),i=npsi,npsi-nxinterp+1,-1)
        write(ndskeqdat,'("neprime")')
        write(ndskeqdat,200) (nprime(i),i=npsi,npsi-nxinterp+1,-1)
        write(ndskeqdat,'("Te")')
        write(ndskeqdat,200) (tel(i),i=npsi,npsi-nxinterp+1,-1)
        write(ndskeqdat,'("Teprime")')
        write(ndskeqdat,200) (teprime(i),i=npsi,npsi-nxinterp+1,-1)
        write(ndskeqdat,'("tion")')
        write(ndskeqdat,200) (tion(i),i=npsi,npsi-nxinterp+1,-1)
        write(ndskeqdat,'("Tiprime")')
        write(ndskeqdat,200) (tiprime(i),i=npsi,npsi-nxinterp+1,-1)
      end if
      write(ndskeqdat,'("pressure")')
      write(ndskeqdat,200) (pressval(i),i=npsi,npsi-nxinterp+1,-1)


!      write(ndskeqdat,'("X")')

! 12/19/00 stop wasting space by writeing x&z to .eqdat, these
      !    are not being read by the plas code
!      do j=1,npts
!         write(ndskeqdat,200) (xs(i,j),i=npsi,npsi-nxinterp+1,-1)
!      enddo
!      write(ndskeqdat,'("Z")')
!      do j=1,npts
!         write(ndskeqdat,200) (zs(i,j),i=npsi,npsi-nxinterp+1,-1)
!      enddo
 200 format(5d24.16)
      close(ndskeqdat)

    return
      
end subroutine wrteqdat


subroutine wrtedge

!-----------------------------------------------------------------------
!     to test mapper's output
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
      integer i,j,ndskrz

      ndskrz=92
      if (verbose .gt. 4) write(outfile,*) 'writing ',runname(1:lrunname)//'.surf'
      open(unit=ndskrz,file=runname(1:lrunname)//'.surf', &
        status='unknown')

      write(ndskrz,*)npts,nxinterp
      do i=npsi,npsi-nxinterp+1,-1
         do j=1,npts
            write(ndskrz,'(1p,3d26.16)') xs(i,j),zs(i,j),bps(i,j)
         enddo
      end do
      close(ndskrz)
      if (verbose .gt. 4) write(outfile,*) &
           'finished writing ',runname(1:lrunname)//'.surf'

    return
      
end subroutine wrtedge


subroutine wrtdsk
!-----------------------------------------------------------------------
! create and write the file dskeq_edge
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
      integer i,nunit,ios
      nunit=27
      open(unit=nunit,file=runname(1:lrunname)//'.eq', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem with opening file ',runname(1:lrunname)//'.eq'
         stop
      endif
      write(nunit,*) npts
      write(nunit,*) f_surf,q_surf,qmin
      write(nunit,*) (theta(i),i=1,npts)
      write(nunit,*) (rpts(i),i=1,npts)
      write(nunit,*) (zpts(i),i=1,npts)
      write(nunit,*) (drdt(i),i=1,npts)
      write(nunit,*) (dzdt(i),i=1,npts)
      write(nunit,*) (bppts(i),i=1,npts)
      write(nunit,*) (om_edge(i),i=1,npts)
      close(nunit)
      return
    
end subroutine wrtdsk


subroutine wrtpitch

!-----------------------------------------------------------------------
!     write pitch angle info on outer midplane for libeam comparison
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
      integer i,j,nunit,ios
      real bpovbt
      nunit=37
      open(unit=nunit,file=runname(1:lrunname)//'.pitch', &
           status='unknown',iostat=ios)
      if(ios.ne.0) then
         write(6,*) 'problem with opening file ', &
              runname(1:lrunname)//'.pitch'
         stop
      endif
      write(nunit,*) 'x,z,Bt,Bp,Bp/Bt,gamma(rad),gamma(deg)'
      do i=npsi,1,-1
         j=1
         bpovbt=bps(i,j)/(fval(i)/xs(i,j))
         write(nunit,250) xs(i,j),zs(i,j),fval(i)/xs(i,j), &
              bps(i,j),bpovbt,atan(bpovbt),180./pi*atan(bpovbt)
      enddo


250   format(7d14.5)

end subroutine wrtpitch
