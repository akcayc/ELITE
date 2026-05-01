      
subroutine rdeqdsk

!-----------------------------------------------------------------------
!     changed dpsi to dpsiv: rlm 7/3/96
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
      character*6 ntitle(5)
      character*6 dat
      integer ipestg
      real xdim,zdim,rc,redge,zmid
      real btor
      real totcur,psimx(2),xax(2)
      real zax(2),psisep,xsep,zsep
      real dpsiv
      integer i,j,isave
      integer neqdsk,ios,ndskeqwrt

      neqdsk=22
      open(unit=neqdsk,file=runname(1:lrunname)//'.eqdsk', &
           status='old',iostat=ios)
      if(ios==0 .and. (verbose .ge. 1)) write(*,*) 'reading eqdsk file ', &
           runname(1:lrunname)//'.eqdsk'
      if(ios.ne.0) then
         if (verbose .ge. 1) then
            write(6,*) 'problem opening file ',runname(1:lrunname)//'.eqdsk'
            write(6,*) 'will try to read file named eqdsk instead'
         endif
         open(unit=neqdsk,file='eqdsk',status='old',iostat=ios)
         if(ios.ne.0) then
           write(*,*) 'could not open eqdsk file'
           stop
         endif
         if (verbose .ge. 1) write(*,*) 'reading equilibrium from file eqdsk'
      endif
!     read eqdsk file
      if (verbose .gt. 4) write(outfile,'("begin read eqdsk")')      
      read(neqdsk,200)(ntitle(i),i=1,5),dat,ipestg,nx,nz
 200  format(6a8,3i4)
      read(neqdsk,300)xdim,zdim,rc,redge,zmid
      read(neqdsk,300)xaxis,zaxis,psiaxis,psilim,btor
      read(neqdsk,300)totcur,psimx(1),psimx(2),xax(1),xax(2)
      read(neqdsk,300)zax(1),zax(2),psisep,xsep,zsep

      if (verbose .gt. 4) write(outfile,*) &
           'xaxis=',xaxis,'zaxis=',zaxis,'xsep=',xsep, &
           'zsep=',zsep,'redge=',redge,'zmid=',zmid,'xax=',xax, &
           'btor=',btor,'totcur=',totcur,'rc=',rc,'shift=',xaxis-rc
!      stop

! 7/07 PBS, write basic efit params for use in post-processing
      ndskeqwrt=87
      if (verbose .gt. 3) write(outfile,*) 'writing ',runname(1:lrunname)//'.efitdat'
      if (verbose .ge. 1) then
         open(unit=ndskeqwrt,file=runname(1:lrunname)//'.efitdat', &
              status='unknown')
         write(ndskeqwrt,*) 'parameters read from g-eqdsk file'
         write(ndskeqwrt,*) '  xaxis   zaxis   btor   totcur'
         write(ndskeqwrt,250) xaxis,zaxis,btor,totcur
         write(ndskeqwrt,*) '  rc   psiaxis   psisep'
         write(ndskeqwrt,250) rc,psiaxis,psisep
         close(ndskeqwrt)
      endif
250   format(7e16.8)
      

      if (verbose .gt. 4) write(outfile,*) &
           'allocating eqdsk arrays, nx= ',nx,' nz= ',nz
      call alloc_eqdsk

! PBS 8/04 issue warning for low resolution efit files
      if ((nx<128) .or. (nz<128)) then
         write(6,*) '********************************************************'
         write(6,*) '**WARNING, the resolution of this eqdsk file is below  *'
         write(6,*) '** the minumum recommended resolution (129x129) for    *'
         write(6,*) '** edge stability studies.  Results are likely to be   *'
         write(6,*) '** unreliable, particularly if there is an X-point     *'
         write(6,*) '********************************************************'
         if (verbose .gt. 3) then
            write(outfile,*) '********************************************************'
            write(outfile,*) '**WARNING, the resolution of this eqdsk file is below  *'
            write(outfile,*) '** the minumum recommended resolution (129x129) for    *'
            write(outfile,*) '** edge stability studies.  Results are likely to be   *'
            write(outfile,*) '** unreliable, particularly if there is an X-point     *'
            write(outfile,*) '********************************************************'
         endif
      endif

      read(neqdsk,300)(sf(i),i=1,nx)
      read(neqdsk,300)(sp(i),i=1,nx)
      read(neqdsk,300)(sffp(i),i=1,nx)
      read(neqdsk,300)(spp(i),i=1,nx)
      read(neqdsk,300)((psixz(i,j),i=1,nx),j=1,nz)
 300  format(5e16.9)
      kvtor=0
      nmass=0
      read(neqdsk,300,end=500) (qpsi(i),i=1,nx)
      read(neqdsk,'(2i5)') nbndry,nlim

      if (verbose .gt. 4) write(outfile,*) &
           'allocating eqdsk boundary arrays, nbndry= ',nbndry, &
          'nlim= ',nlim
      call alloc_bndry

      read(neqdsk,300) (xbndry(i),zbndry(i),i=1,nbndry)
      read(neqdsk,300) (xlim(i),zlim(i),i=1,nlim)
      if (verbose .gt. 4) write(outfile,'("finished eqdsk")')
      if(rotate.ne.0) then
         read(neqdsk,'(i5,e16.9,i5)',end=500,err=500) kvtor,rvtor,nmass
         if(kvtor.gt.0) then
            read(neqdsk,300) (pressw(i),i=1,nx)
            read(neqdsk,300) (pwprim(i),i=1,nx)
!     guard against negative or zero pressw
            isave=1
            do i=1,nx
               if(pressw(i).le.0.) then
                  if (verbose .ge. 1) write(*,*) 'found negative pressw i=',i
                  isave=i
                  go to 75
               endif
            end do
 75         continue
            if(isave.ne.1) then
               do i=isave,nx
                  pressw(i)=pressw(isave-1)+(i-isave)/(nx-isave)* &
                      (1.-pressw(isave-1))
                  pwprim(i)=(pressw(isave-1)-1.)/ &
                      (psilim-psiaxis)*(nx-isave)/(nx-1)
               end do
            endif
         endif
         if(nmass.gt.0) then
            read(neqdsk,300) (rho0(i),i=1,nx)
            dpsiv=(psilim-psiaxis)/(nx-1.)
            do i=2,nx-1
               rho0p(i)=(rho0(i+1)-rho0(i-1))/(2.*dpsiv)
            end do
            rho0p(1)=(-3.*rho0(1)+4.*rho0(2)-rho0(3))/(2.*dpsiv)
            rho0p(nx)=(3.*rho0(nx)-4.*rho0(nx-1)+rho0(nx-2))/(2.*dpsiv)
         endif
      endif
!     generate x,z gridd
 500  continue
      dx=xdim/(nx-1.)
      dz=zdim/(nz-1.)
      do i=1,nx
         xgrid(i)=redge+(i-1.)*dx
      end do
      do i=1,nz
         zgrid(i)=-0.5*zdim+(i-1.)*dz
      end do
!     if sf is negative change it's sign
!     I don't know what the purpose of a negative B field is.
      if(sf(nx).lt.0.) then
         if (verbose .ge. 1) write(*,*) &
              'read negative value for sf, changing sign of sf,not sffp'
         do i=1,nx
            sf(i)=-sf(i)
!            sffp(i)=-sffp(i)
         end do
      endif

      close(neqdsk)
      call getbpsq(psixz,nx,nz,xgrid,dx,dz,nx,nz,bpsq)
!      write(42,'(4e20.12)')((bpsq(i,j),i=1,nx,2),j=1,nz,2)
    return
      
end subroutine rdeqdsk


subroutine getbpsq(psixz,nxd,nzd,xgrid,dx,dz,nx,nz,bpsq)

      integer nxd,nzd,nx,nz
      real psixz(nxd,nzd)
      real bpsq(nxd,nzd)
      real xgrid(nxd)
      real dx,dz
      real dpsixsq
      integer i,j
      do i=2,nx-1
         do j=2,nz-1
            bpsq(i,j)=(((psixz(i+1,j)-psixz(i-1,j))/(2.*dx))**2+ &
               ((psixz(i,j+1)-psixz(i,j-1))/(2.*dz))**2)/ &
                xgrid(i)**2
         end do
      end do
      i=1
      do j=2,nz-1
         bpsq(i,j)=(((-3.*psixz(i,j)+4.*psixz(i+1,j) &
             -psixz(i+2,j))/(2.*dx))**2+ &
             ((psixz(i,j+1)-psixz(i,j-1))/(2.*dz))**2)/ &
             xgrid(i)**2
      end do
      i=nx
      do j=2,nz-1
         bpsq(i,j)=(((3.*psixz(i,j)-4.*psixz(i-1,j) &
             +psixz(i-2,j))/(2.*dx))**2+ &
             ((psixz(i,j+1)-psixz(i,j-1))/(2.*dz))**2)/ &
             xgrid(i)**2
      end do
      j=1
      do i=1,nx
         if(i.eq.1) then
            dpsixsq=((-3.*psixz(i,j)+4.*psixz(i+1,j) &
               -psixz(i+2,j))/(2.*dx))**2
         else if(i.eq.nx) then
            dpsixsq=((3.*psixz(i,j)-4.*psixz(i-1,j) &
                +psixz(i-2,j))/(2.*dx))**2
         else
            dpsixsq=((psixz(i+1,j)-psixz(i-1,j))/(2.*dx))**2
         endif
         bpsq(i,j)=(dpsixsq+ &
             ((-3.*psixz(i,j)+4.*psixz(i,j+1) &
             -psixz(i,j+2))/(2.*dz))**2)/ &
             xgrid(i)**2
      end do
      j=nz
      do i=1,nx
         if(i.eq.1) then
            dpsixsq=((-3.*psixz(i,j)+4.*psixz(i+1,j) &
              -psixz(i+2,j))/(2.*dx))**2
         else if(i.eq.nx) then
            dpsixsq=((3.*psixz(i,j)-4.*psixz(i-1,j) &
               +psixz(i-2,j))/(2.*dx))**2
         else
            dpsixsq=((psixz(i+1,j)-psixz(i-1,j))/(2.*dx))**2
         endif
         bpsq(i,j)=(dpsixsq+ &
             ((3.*psixz(i,j)-4.*psixz(i,j-1) &
             +psixz(i,j-2))/(2.*dz))**2)/ &
             xgrid(i)**2
      end do
 


   return
      
end subroutine getbpsq

