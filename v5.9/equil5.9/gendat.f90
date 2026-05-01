
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

subroutine gendat

  ! Snyder 10/19/00: preparing for addition of 1/n terms, no longer
  !   adjust B_p so that calculated q values equal those from the
  !   eqdsk or dskbal file, add calculation of ppp=p'' and ffpp=(ff')'
  !
  ! Snyder 5/10/00: redo code to use q values from the eqdsk or
  !   dskbal files.  adjust bp on all surfaces so that calculated
  !   q will equal the original equilibrium value
  ! (the fact that they agree so poorly is still worrisome)
  !
  ! Snyder 1/26/00: change code over to f90 and begin preparations
  !   for handling multiple equilibrium flux surfaces.

  !-----------------------------------------------------------------------
  ! generate the following  equilibrium data
  !      (theta(i),i=1,npts) 0-2pi otherwise theta can be anything
  !      (rpts(i),i=1,npts)
  !      (zpts(i),i=1,npts)
  !      (drdt(i),i=1,npts)
  !      (dzdt(i),i=1,npts)
  !      (bppts(i),i=1,npts)
  !      (omega(i),i=1,npts)
  ! note that in the original refsurf, theta was a local variable
  !   and increased in counter-clockwise direction
  !   but the Bishop convention of length increasing in clockwise
  !   direction provided the real poloidal varible
  !-----------------------------------------------------------------------
  use eliteeq_data
  implicit none
  integer nunit,ios,ipass
  real, dimension(npts) :: leng,r,y2,y3
  real bnd_set,rpt_v,zpt_v,int_v
  real dr,drt,dtheta
  real rnewton
  real bpnorm,q_nonorm,dl,fconst
  real aspect,btnorm,rvnorm,s,alpha
  real, dimension(npts) ::  ccrpts,cczpts,ccrbppts
  real rnorm
  real xtemp,ytemp
  integer mks,nsurf
  integer i,ir,ipsi,m0,mmin,j
  integer count,old_npsi1,old_npsi2
  real qold
  real qtemp,del_min,q_calc,psifrac
  real, dimension(:), allocatable :: q_array
  integer nn_min
  real, dimension(npts) :: rl,zl,bpl,work,leq,nu
  real dchi
  !      real, dimension(nxinterp) :: q_plas,psigrid,workq,rmaxgrid
  !   nxinterp can now change for toq files
  real, dimension(:), allocatable :: q_plas,psigrid,workq,rmaxgrid
  integer scalefact,imin
  !      parameter (scalefact=10)
  !      real q_fine(scalefact*nxinterp),psifine(scalefact*nxinterp)
  real q_fine,psiwork,psiqmin,psitotal,rmax
  logical reverse   ! set to true if q reversal found
  real edgeqerr

  pi=acos(-1.d0)
  reverse=.false.

  if (verbose .gt. 4) then 
     write(outfile,*) 'here in gendat, npts= ',npts
  !-----------------------------------------------------------------------
  ! 1.0 read or generate r,z,bp
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Multi-flux surface version 2/3/00
  ! first option will be to start with an EFIT (eqdsk) file, and map it
  !   into flux surfaces using routines based on those in Miller's bal1.6
  !-----------------------------------------------------------------------
     write(outfile,*) 'shape= ',shape
  endif

!------------------------------------------------------------
!  Routine to handle any flux surface gridded equilibrium file
!-----------------------------------------------------------
!   updated 7.8.04 to include eqbm format
!------------------------------------------------------------TOQ
  if (shape.eq.'toq' .or. shape.eq.'gato' .or. shape.eq.'ngat' &
       .or. shape.eq.'ogat' .or. shape.eq.'pest' .or. &
       shape.eq.'eqbm') then   
     ! read a dskbal file produced by toq or dskgato produced by toq/teq
     !   or eqbm file or pest
     if (shape.eq.'toq') then
        if (verbose .gt. 4) write(outfile,*) 'call readbal'
        call readbal  ! routine to read dskbalnew file named runname.dskbal
     else if (shape.eq.'eqbm') then
        if (verbose .gt. 4) write(outfile,*) 'call readeq'
        call readeq
     else if (shape.eq.'pest') then
        if (verbose .gt. 4) write(outfile,*) 'call readpest'
        call readpest
     else
        if (verbose .gt. 4) write(outfile,*) 'call readgato'
        call readgato
     endif
     if (verbose .gt. 4) write(outfile,*) 'call getbps'
     call getbps
     if (verbose .gt. 4) write(outfile,*) 'return from getbps'

! TEMP
!     write(*,*) 'first and last point for each surface'
!     do i=1,npsi
!        write(*,*) i,xs(i,1),xs(i,npts),xs(i,npts-1),zs(i,1),zs(i,npts)
!     enddo


! 5/05 write out dskgato file for use by gato if requested
!  (useful for convering eqin or dskbal files)
     if (wrtgato) then
        if (verbose .gt. 4) write(outfile,*) 'call writegato'
        call writegato
     endif

     allocate( q_array(npsi) )

! go through and calculate q and x for each needed surface

     qold=1e20   ! used to check monotonicity of q
     ipsi=npsi
     ipass=0
     fconst=0.  ! used to shift f**2 if setdel=true
22   continue
     fval(ipsi)=sqrt(fval(ipsi)**2+fconst) ! adjust if setdel true
     ipass=ipass+1
     psifrac=(psiv(ipsi)-psiv(1))/(psiv(npsi)-psiv(1))
!     write(*,*) 'ipsi=',ipsi,'psifrac=',psifrac

     leng(1)=0.
     do i=2,npts
        dr=xs(ipsi,i)-xs(ipsi,i-1)
        drt=zs(ipsi,i)-zs(ipsi,i-1)
        leng(i)=leng(i-1)+sqrt(dr**2+drt**2)
     end do

     omega(1)=0.
     do i=2,npts
        dl=(leng(i)-leng(i-1))
        if ((xs(ipsi,i)==0.) .or. (bps(ipsi,i)==0.)) then
           write(*,*) 'xs or bps is zero, ipsi=',ipsi,'i=',i, &
                'xs=',xs(ipsi,i),'bps=',bps(ipsi,i)
           stop
        endif
        omega(i)=omega(i-1)+0.5*dl*fval(ipsi)* &
             (1./(xs(ipsi,i)**2*bps(ipsi,i))+ &
             1./(xs(ipsi,i-1)**2*bps(ipsi,i-1)))
     end do
     q_calc=omega(npts)/(2.*pi)
     if (q_calc.ge.qold) then
        write(6,*) 'q profile reversed at psi=',psifrac
        write(6,*) 'q=',q_calc,' qold=',qold
        if (verbose .ge. 2) write(outfile,*) 'q profile reversed at psi=',psifrac
        if (verbose .ge. 2) write(outfile,*) 'q=',q_calc,' qold=',qold
        reverse=.true.
        !              stop
     endif
     call qspline(ipsi,q_calc) 
     qold=q_calc

     if (ipass == 1) then
! First pass through is edge flux surface
        if (setdel) then
! Shift f**2 by a constant to get required del (qa)....does not affect ffprime
! and therefore still satisfies GS equation
           if (verbose .ge. 1) write(6,*)' q_calc no spline=',q_calc
           if (verbose .gt. 3) write(outfile,*)' q_calc no spline=',q_calc
           call qspline(npsi,q_calc)
           if (verbose .ge. 1) write(6,*)' q_calc spline=',q_calc
           if (verbose .gt. 3) write(outfile,*)' q_calc spline=',q_calc
           m0=int(nn*q_calc+1.)
! 5/05 for del_fix=-1., adjust f*2 to force q_edge=qsfin(npsi)=edge q read
!      in from equilibrium file
           if (del_fix.ne.-1.) then
              fconst=fval(ipsi)**2*(((m0-del_fix)/(nn*q_calc))**2-1.)
           else
              write(*,*) 'del_fix=-1., adjusting f**2 so that calculated ', &
                   'edge q exactly matches edge q from equilibrium file=',qsfin(npsi)
              if (verbose .gt. 3) write(outfile,*) &
                   'del_fix=-1., adjusting f**2 so that calculated ', &
                   'edge q exactly matches edge q from equilibrium file',qsfin(npsi)
              fconst=fval(ipsi)**2*((qsfin(npsi)/q_calc)**2-1.)
           endif
           qold=1.e20  ! reset to avoid false finding of reversal
           goto 22
        end if
     end if

     if (ipsi == npsi) then
        if (delmin) then  ! find nn which minimizes del
           del_min=1e10
           do i=nnmin,nnmax
              m0=int(i*q_calc+1.)   
              q0=real(m0)/real(i)
              del=(q0-q_calc)*i
              if (del<del_min) then
                 del_min=del
                 nn=i
              endif
           enddo
           write(*,*) 'minimized del over range ',nnmin, &
                '<=nn<=',nnmax
           write(*,*) 'assigned nn=',nn,' del_min=',del_min
           if (verbose .gt. 3) write(outfile,*) 'minimized del over range ',nnmin, &
                '<=nn<=',nnmax
           if (verbose .gt. 3) write(outfile,*) 'assigned nn=',nn,' del_min=',del_min
        endif
        m0=int(nn*q_calc+1.)
        q0=real(m0)/real(nn)
        del=(q0-q_calc)*nn
        if (verbose .gt. 3) write(outfile,*) &
             'qa= ',q_calc,' m0=',m0,' q0=',q0,' del=',del
        if (verbose .ge. 1) write(*,*) 'qa= ',q_calc,' m0=',m0,' q0=',q0,' del=',del
!   5/20  no longer using or defining nmvac and mmax in equil
!        mmax=m0+nmvac-1
        q_surf=q_calc
     endif

     if (verbose .gt. 4) write(outfile,*) &
          'surface= ',ipsi,' q_spline=', qsfin(ipsi), &
          ' q_calc=',q_calc,' x=',nn*(q0-q_calc)

     if (psifrac > psimin ) then
        ipsi=ipsi-1
        if (ipsi <= 1) then
           write(*,*) 'reached inmost surface with psifrac>psimin'
           stop
        endif
        goto 22
     endif

     nsurf=npsi-ipsi+1
!     write(*,*) 'psifrac <= psimin, total equil surfaces=',nsurf
     if (verbose .gt. 3) write(outfile,*) &
          'psifrac <= psimin, total equil surfaces=',nsurf
!     write(*,*) 'For now, ELITE will be forced to use this number'
!     write(*,*) ' of surfaces.  Make sure to reset this in input file'
!     write(*,*) 'nxinterp=nsurf=',nsurf
     nxinterp=nsurf

! now need to allocate nxinterp size arrays
     allocate( q_plas(nxinterp),psigrid(nxinterp), &
          workq(nxinterp),rmaxgrid(nxinterp) )

! Evaluate q as in ELITE, a test of code
     do ipsi=npsi,npsi-nxinterp+1,-1
        omega(1)=0.
        leng(1)=0.
        do i=2,npts
           dr=xs(ipsi,i)-xs(ipsi,i-1)
           drt=zs(ipsi,i)-zs(ipsi,i-1)
           leng(i)=leng(i-1)+sqrt(dr**2+drt**2)
           dl=(leng(i)-leng(i-1))
           omega(i)=omega(i-1)+0.5*dl*fval(ipsi)* &
                (1./(xs(ipsi,i)**2*bps(ipsi,i))+ &
                1./(xs(ipsi,i-1)**2*bps(ipsi,i-1)))
        end do
        q_calc=omega(npts)/(2.*pi)
        bpnorm=q_calc/qsfin(ipsi)
!        write(*,*) 'surface=',ipsi,' bpnorm=',bpnorm
        if (verbose .gt. 4) write(outfile,*) 'surface=',ipsi,' bpnorm=',bpnorm
        !          do i=1,npts
        !             bps(ipsi,i)=bpnorm*bps(ipsi,i)
        !         omega(i)=omega(i)/(bpnorm*q_surf)
        !          end do
        !          write(*,*) 'Not normalizing B_p'
     end do
!     write(*,*) 'Not normalizing B_p'
!     write(outfile,*) 'Not normalizing B_p'

     ! check it
     do ipsi=npsi,npsi-nxinterp+1,-1
        omega(1)=0.
        leng(1)=0.
        do i=2,npts
           dr=xs(ipsi,i)-xs(ipsi,i-1)
           drt=zs(ipsi,i)-zs(ipsi,i-1)
           leng(i)=leng(i-1)+sqrt(dr**2+drt**2)
           dl=(leng(i)-leng(i-1))
           omega(i)=omega(i-1)+0.5*dl*fval(ipsi)* &
                (1./(xs(ipsi,i)**2*bps(ipsi,i))+ &
                1./(xs(ipsi,i-1)**2*bps(ipsi,i-1)))
        end do
        q_calc=omega(npts)/(2.*pi)
        q_array(ipsi)=q_calc
        qcalc(ipsi)=q_calc
!       write(outfile,*) 'after normalization, q_calc=',q_calc, &
 !            ' qsfin(ipsi)=', qsfin(ipsi)
     end do


     !   write output data
     call wrteqdat
     call wrtedge !write dskrzbp_in for edge code

! 6/03 write pitch angle data for libeam comparisons
!     call wrtpitch

     do i=npsi,npsi-nxinterp+1,-1
        psigrid(npsi-i+1)=(psiv(i)-psiv(1))/(psiv(npsi)-psiv(1))
     enddo

         psitotal=psiv(npsi)-psiv(1)
    
  endif

!-----------------------------------------------------------
!  Procedure to read Efit g files
!------------------------------------------------------------EFIT
  if (shape.eq.'eqds') then   ! shape is a 4 chr variable for now
!     allocate( q_plas(nxinterp),psigrid(nxinterp), &
!          workq(nxinterp),rmaxgrid(nxinterp) )
     if (verbose .gt. 4) write(outfile,*) 'call rdeqdsk'
     call rdeqdsk   ! routine to read eqdsk file named runname.eqdsk

     count=0  ! count for number of iterations to get nsurf=nxinterp
     old_npsi1=0 ! store old values to check if search is falling into loop
     old_npsi2=0
5    continue   ! start of big loop to do mapping and iterate until
     !   nsurf=nxinterp (desired # of equil surfaces in range)
     if (verbose .gt. 4) write(outfile,*) 'call mapperb'
     call mapperb

     ! go through and calculate q and x for each needed surface

     !        qold=qsfin(npsi)+1.e-6   ! used to check monotonicity of q
     qold=1.e30
     ipsi=npsi
     fconst=0.
     ipass=0
10   continue
     ipass=ipass+1
     fval(ipsi)=sqrt(fval(ipsi)**2+fconst) ! adjust if setdel true
     !        do ipsi=npsi,1,-1

     psifrac=(psiv(ipsi)/1.e8-psiaxis)/(psilim-psiaxis)
     if (verbose .gt. 4) write(outfile,*) 'ipsi=',ipsi,'psifrac=',psifrac

     leng(1)=0.
     do i=2,npts
        dr=xs(ipsi,i)-xs(ipsi,i-1)
        drt=zs(ipsi,i)-zs(ipsi,i-1)
        leng(i)=leng(i-1)+sqrt(dr**2+drt**2)
     end do

     omega(1)=0.
     do i=2,npts
        dl=(leng(i)-leng(i-1))
        omega(i)=omega(i-1)+0.5*dl*fval(ipsi)* &
             (1./(xs(ipsi,i)**2*bps(ipsi,i))+ &
             1./(xs(ipsi,i-1)**2*bps(ipsi,i-1)))
     end do
     q_calc=omega(npts)/(2.*pi)
     if (q_calc.ge.qold) then
        reverse=.true.
        write(6,*) 'q profile reversed at psi=',psifrac
        write(6,*) 'q=',q_calc,' qold=',qold
        if (verbose .ge. 2) write(outfile,*) 'q profile reversed at psi=',psifrac
        if (verbose .ge. 2) write(outfile,*) 'q=',q_calc,' qold=',qold
        !              stop
     endif
     qold=q_calc

     if (ipass == 1) then
! First pass through is edge flux surface
        if (setdel) then
! Shift f**2 by a constant to get required del (qa)....does not affect ffprime
! and therefore still satisfies GS equation
           if (verbose .ge. 1) write(6,*)' q_calc no spline=',q_calc
           if (verbose .gt. 3) write(outfile,*)' q_calc no spline=',q_calc
           call qspline(npsi,q_calc)
           if (verbose .ge. 1) write(6,*)' q_calc spline=',q_calc
           if (verbose .gt. 3) write(outfile,*)' q_calc spline=',q_calc
           m0=int(nn*q_calc+1.)
! 5/05 for del_fix=-1., adjust f*2 to force q_edge=qsfin(npsi)=edge q read
!      in from equilibrium file
           if (del_fix.ne.-1.) then
              fconst=fval(ipsi)**2*(((m0-del_fix)/(nn*q_calc))**2-1.)
           else
              write(*,*) 'del_fix=-1., adjusting f**2 so that calculated ', &
                   'edge q exactly matches edge q from equilibrium file=',qsfin(npsi)
              if (verbose .gt. 3) write(outfile,*) &
                   'del_fix=-1., adjusting f**2 so that calculated ', &
                   'edge q exactly matches edge q from equilibrium file',qsfin(npsi)
              fconst=fval(ipsi)**2*((qsfin(npsi)/q_calc)**2-1.)
           endif
           qold=1.e20  ! reset to avoid false finding of reversal
           goto 10
        end if
     end if


     if (ipsi == npsi) then
        if (delmin) then  ! find nn which minimizes del
           del_min=1e10
           do i=nnmin,nnmax
              m0=int(i*q_calc+1.)   
              q0=real(m0)/real(i)
              del=(q0-q_calc)*i
              if (del<del_min) then
                 del_min=del
                 nn=i
              endif
           enddo
           write(*,*) 'minimized del over range ',nnmin, &
                '<=nn<=',nnmax
           write(*,*) 'assigned nn=',nn,' del_min=',del_min
           if (verbose .gt. 3) write(outfile,*) 'minimized del over range ',nnmin, &
                '<=nn<=',nnmax
           if (verbose .gt. 3) write(outfile,*) 'assigned nn=',nn,' del_min=',del_min
        endif
        m0=int(nn*q_calc+1.)
        q0=real(m0)/real(nn)
        del=(q0-q_calc)*nn
        if (verbose .gt. 3) write(outfile,*) 'qa=',q_calc,' m0=',m0
        if (verbose .gt. 3) write(outfile,*) 'q0=',q0,' del=',del
        !              mmax=m0+nmvac-1
        !              mmin=mmax-nm+1
        !              xmax=m0-mmin+xwid
        !              qmin=q0-xmax/real(nn)
        q_surf=q_calc
        !              write(outfile,*) 'mmax=',mmax,' mmin=',mmin,'xmax= ',xmax, &
        !                   ' qmin=',qmin
        !              if (mmin < 1) then
        !                 write(*,*) 'mmin=m0+nmvac-nm must be positive'
        !                 stop
        !              endif

     endif

     if (verbose .gt. 4) write(outfile,*) &
          'surface= ',ipsi,' q_spline=', qsfin(ipsi), &
          ' q_calc=',q_calc,' x=',nn*(q0-q_calc)

     !          write(*,*) 'surface=',ipsi,' q_spline=', qsfin(ipsi), &
     !             ' q_calc=',q_calc,' x=',nn*(q0-q_calc),'qprime=',qsfinp(ipsi)

     !          if ((nn*(q0-q_calc)) < xmax) then
     if (psifrac > psimin ) then          
        ipsi=ipsi-1
        if (ipsi <= 1) then
           write(*,*) 'reached inmost surface with psifrac > psimin'
           stop
        endif
        goto 10
     endif

     nsurf=npsi-ipsi+1
     if (verbose .gt. 3) write(outfile,*) &
          'psifrac<psimin , total equil surfaces=',nsurf
!     write(*,*) 'psifrac<psimin , total equil surfaces=',nsurf

     if (nsurf .ne. nxinterp) then
      ! for alpsi=-2 use equal spaced efit-like grid, and do
      !   not iterate npsi
        if (alpsi.eq.-2) then
           if (verbose .ge. 1) write(*,*) 'alpsi=-2, setting nxinterp=nsurf=',nsurf
           if (verbose .gt. 3) write(outfile,*) &
                'alpsi=-2, setting nxinterp=nsurf=',nsurf
           nxinterp=nsurf
        else if (alpsi.eq.-3) then
           if (verbose .ge. 1) write(*,*) 'alpsi=-3, setting nxinterp=nsurf=',nsurf
           if (verbose .gt. 3) write(outfile,*) &
                'alpsi=-3, setting nxinterp=nsurf=',nsurf
           nxinterp=nsurf
        else
           if (verbose .gt. 3) write(outfile,*) 'nxinterp= ',nxinterp
           count=count+1
           if (count.gt.10) then
              write(*,*) 'too many iterations, unable to get nsurf=nxinterp'
              stop
           endif
           old_npsi2=old_npsi1
           old_npsi1=npsi
           npsi=real(npsi)*real(real(nxinterp-1.)/real(nsurf-1.))
           if (npsi == old_npsi1) then
              if (nxinterp > nsurf) npsi=npsi+1
              if (nxinterp < nsurf) npsi=npsi-1
           endif
           if (npsi == old_npsi2) npsi=real(old_npsi1+npsi)/2.
           if (verbose .gt. 3) write(outfile,*) 'rerunning with npsi= ',npsi
           call dealloc_mapper  ! deallocate mapper vars
           goto 5   ! go back and remap with new npsi
        endif
     endif

     allocate( q_plas(nxinterp),psigrid(nxinterp), &
          workq(nxinterp),rmaxgrid(nxinterp) )

     allocate( q_array(npsi) )
     ! must now adjust bp so that calculated q will match equilibrium values
     ! no longer adjusting 10/00, but still calculate bpnorm
     if (verbose .gt. 4) write(outfile,*) &
          'calculating normalizing factors bpnorm=q_calc/qsfin(ipsi)'
     do ipsi=npsi,npsi-nxinterp+1,-1
        omega(1)=0.
        leng(1)=0.
        do i=2,npts
           dr=xs(ipsi,i)-xs(ipsi,i-1)
           drt=zs(ipsi,i)-zs(ipsi,i-1)
           leng(i)=leng(i-1)+sqrt(dr**2+drt**2)
           dl=(leng(i)-leng(i-1))
           omega(i)=omega(i-1)+0.5*dl*fval(ipsi)* &
                (1./(xs(ipsi,i)**2*bps(ipsi,i))+ &
                1./(xs(ipsi,i-1)**2*bps(ipsi,i-1)))
        end do
        q_calc=omega(npts)/(2.*pi)
        q_array(ipsi)=q_calc
        bpnorm=q_calc/qsfin(ipsi)
!        write(*,*) 'surface=',ipsi,' bpnorm=',bpnorm
        if (verbose .gt. 4) write(outfile,*) 'surface=',ipsi,' bpnorm=',bpnorm
        !          do i=1,npts
        !             bps(ipsi,i)=bpnorm*bps(ipsi,i)
        !         omega(i)=omega(i)/(bpnorm*q_surf)
        !          end do
     end do
!     write(*,*) 'not normalizing!'
!     write(*,*) 'using calculated q rather than q from the eqdsk'
!     write(outfile,*) 'not normalizing!'
!     write(outfile,*) 'using calculated q rather than q from the eqdsk'

!     if (ipass == 1) then
! First pass through is edge flux surface


     ! check it
     do ipsi=npsi,npsi-nxinterp+1,-1
        omega(1)=0.
        leng(1)=0.
        do i=2,npts
           dr=xs(ipsi,i)-xs(ipsi,i-1)
           drt=zs(ipsi,i)-zs(ipsi,i-1)
           leng(i)=leng(i-1)+sqrt(dr**2+drt**2)
           dl=(leng(i)-leng(i-1))
           omega(i)=omega(i-1)+0.5*dl*fval(ipsi)* &
                (1./(xs(ipsi,i)**2*bps(ipsi,i))+ &
                1./(xs(ipsi,i-1)**2*bps(ipsi,i-1)))
        end do
        q_calc=omega(npts)/(2.*pi)
        qcalc(ipsi)=q_calc
        if (verbose .gt. 3) write(outfile,*) &
             'after normalization, q_calc=',q_calc, &
             ' qsfin(ipsi)=', qsfin(ipsi)
     end do

     !  read density and temperature mtanh coefficients and
     !    calculate profiles
     if (dens) then
        call readmtanh
     endif

     !   write output data
     call wrteqdat
     call wrtedge !write dskrzbp_in for edge code

     !     normalize xaxis and zaxis
     zaxis=1.d2*zaxis
     xaxis=1.d2*xaxis

     !   write data for converting x coordinate to psi coordinate
     !      write(outfile,*) 'writing ',runname(1:lrunname)//'.xtopsi'
     !      open(unit=25,file=runname(1:lrunname)//'.xtopsi', &
     !           status='unknown')
     !      write(25,*) nxinterp
     !      write(25,*) (nn*(q0-q_array(i)),i=npsi,npsi-nxinterp+1,-1)
     do i=npsi,npsi-nxinterp+1,-1
        psigrid(npsi-i+1)=(psiv(i)/1.d8-psiaxis)/(psilim-psiaxis)
     enddo
     !      write(25,*) ((psiv(i)/1.e8-psiaxis)/(psilim-psiaxis), &
     !           i=npsi,npsi-nxinterp+1,-1)
     !      write(25,*) 'total flux=',(psilim-psiaxis)*1.e8
     psitotal=(psilim-psiaxis)*1.d8
     !      close(25)

     if (verbose .gt. 3) write(outfile,*) &
          'psi fractions of the equilibrium surfaces:'
     if (verbose .gt.34) write(outfile,*) &
          ((psiv(i)/1.e8-psiaxis)/(psilim-psiaxis), &
          i=npsi,npsi-nxinterp+1,-1)

     !        enddo

     ! 2/3/00 .rzbp & .surf files have now been produced from eqdsk
     !   for now, just read single surface as before below
     !   this single surface will be used in vac (name.eq)

  endif
  !-----------------------------------------------------------------------
  ! 12/19/00 no longer allow options other than 'eqds' or 'toq' or 'gato'
  ! 1.1    shape=='file'    ! or 'eqds' or 'toq'
  !-----------------------------------------------------------------------
  if ((shape.eq.'eqds') .or. (shape.eq.'toq') &
       .or. (shape.eq.'gato') .or. (shape.eq.'ngat') &
       .or. (shape.eq.'ogat') .or. (shape.eq.'pest') .or. (shape.eq.'eqbm')) then
     !         nunit=17
     !         write(*,*) 'reading ',runname(1:lrunname)//'.rzbp'
     !         open(unit=nunit,file=runname(1:lrunname)//'.rzbp',status="old")
     !        write(*,*) npts,fval(npsi)
     !        read(nunit,*)npts,f_surf
     f_surf=fval(npsi)
     if (verbose .gt. 4) write(outfile,*) npts,f_surf
     leng(1)=0.
     do i=1,npts
        ir=npts+1-i
        !           read(nunit,*)rpts(ir),zpts(ir),bppts(ir)
        rpts(ir)=xs(npsi,i)
        zpts(ir)=zs(npsi,i)
        bppts(ir)=bps(npsi,i)
     end do
     do i=2,npts
        dr=rpts(i)-rpts(i-1)
        drt=zpts(i)-zpts(i-1)
        leng(i)=leng(i-1)+sqrt(dr**2+drt**2)
     end do
     !        close(nunit)
     rmajor=(rpts(npts/2+1)+rpts(1))/2.
     rminor=(rpts(1)-rpts(npts/2+1))/2.
     !        write(6,*)' Major radius=',rnorm
     do i=1,npts
        theta(i)=atan2(zpts(i),(rpts(i)-rminor))
     end do
     !        write(*,*) 'here is theta, should be equispaced, ',theta
  else
     write(*,*) 'shape=',shape,' is not supported'
     stop
  end if

  !-----------------------------------------------------------------------
  ! 2.0 still need to:
  !     renormalize bppts to yield correct q
  !     calculate omega
  !     make length/2pi the theta variable
  !     calculate theta derivative of r and z
  !-----------------------------------------------------------------------
  omega(1)=0.
  do i=2,npts
     dl=(leng(i)-leng(i-1))
     omega(i)=omega(i-1)+0.5*dl*f_surf* &
          (1./(rpts(i)**2*bppts(i))+ &
          1./(rpts(i-1)**2*bppts(i-1)))
  end do
  q_nonorm=omega(npts)/(2.*pi)
! HRW and S Saarelma (24/3/04) need to write this omega mesh to vac
! code, but it is over-written below. We introduce om_edge for this
  do i=1,npts
    om_edge(i)=omega(i)/q_nonorm
  end do
  if (verbose .gt. 3) write(outfile,*) &
       'q_nonorm=',q_nonorm,' q_surf=',q_surf, &
       ' qsfin(npsi)=',qsfin(npsi)
  if (verbose .ge. 1) write(*,*) 'q_nonorm=',q_nonorm,' qsfin(npsi)=',qsfin(npsi)

  edgeqerr=abs((q_nonorm-qsfin(npsi))/q_nonorm)
!  write(*,*) 'Edge q calculated from flux surface integral differs from the value'
!  write(*,*) ' read from the equilibrium file 

  if (edgeqerr>0.05) then
     if (verbose .ge. 1) write(6,*) ' '
     if (verbose .ge. 1) write(6,*) &
          '***WARNING!!! Difference between edge q calculated from flux surface integral'
     if (verbose .ge. 1) write(6,*) &
          '  and interpolated from equilibrium file is large (',(100*edgeqerr),'%).'
     if (verbose .gt. 3) write(outfile,*) ' '
     if (verbose .gt. 3) write(outfile,*) &
          '***WARNING!!! Difference between edge q calculated from flux surface integral'
     if (verbose .gt. 3) write(outfile,*) &
          '  and interpolated from equilibrium file is large (',(100*edgeqerr),'%).'
     if (shape.eq.'eqds') then
        if (verbose .ge. 1) write(6,*) &
             '  Reducing percenflux and/or using a higher resolution eqdsk file should help'
        if (verbose .ge. 1) write(6,*) &
             '  Also make sure npts in the equil namelist is sufficiently large'
        if (verbose .gt. 3) write(outfile,*) &
             '  Reducing percenflux and/or using a higher resolution eqdsk file should help'
        if (verbose .gt. 3) write(outfile,*) &
             '  Also make sure npts in the equil namelist is sufficiently large'
     else
        if (verbose .gt. 3) write(outfile,*) &
             '  Higher poloidal resolution may be needed in the original equilibrium file'
     endif
     if (verbose .ge. 1) write(6,*) ' '
     if (verbose .gt. 3) write(outfile,*) ' '
  endif

  if (verbose .ge. 1) write(*,*) 'del=',del
  bpnorm=q_nonorm/q_surf
  if (verbose .gt. 3) write(outfile,*) 'old bpnorm=',bpnorm
  ! change 5/9/00, really shouldn't need to do this, need more
  !    precision in files
  if (verbose .gt. 3) write(outfile,*) 'setting bpnorm=1., setting q_surf=q_nonorm'
  q_surf=q_nonorm
  bpnorm=1.

  if (abs(bpnorm-1.0) > 1e-5) then
     write(*,*) 'old bpnorm=',bpnorm,' should be 1'
     stop
  endif

  do i=1,npts
     !         bppts(i)=bpnorm*bppts(i)
     !         omega(i)=omega(i)/(bpnorm*q_surf)
     omega(i)=omega(i)/q_surf
  end do
  if (abs(omega(npts)-2.*pi)>.01) then
     write(*,*) 'something is wrong, omega(npts)=',omega(npts)
     stop
  endif
  omega(npts)=2.*pi
  ! redefine theta
  do i=1,npts
     theta(i)=leng(i)*2.*pi/leng(npts)
  end do
  ! vacuum calculation requires dr/dtheta and dz/dtheta
  bnd_set=-2.d30

  if (verbose .gt. 4) write(outfile,*) 'about to call spline in gendat'
  call spline(theta,rpts,npts,bnd_set,bnd_set,y2)
  do i=1,npts ! second y2 arg is dummy as are rpt_v and int_v
     call zsplint(theta,rpts,y2,y2,npts,theta(i), &
          rpt_v,drdt(i),int_v)
  end do
  call spline(theta,zpts,npts,bnd_set,bnd_set,y2)
  do i=1,npts ! second y2 arg is dummy as are rpt_v and int_v
     !         write(*,*) i,theta(i),zpt_v,dzdt(i),int_v
     call zsplint(theta,zpts,y2,y2,npts,theta(i), &
          zpt_v,dzdt(i),int_v)

  end do


  ! 9/01 go through all steps from surfdat&nustuff to try to get q
  !   values exactly matching those from the plas code.  Use these
  !   for q_surf and to calculate qmin to be written to output files
  !     of course match is exact only if npts = ns

  !      do j=npsi,npsi-nxinterp+1,-1
  do j=npsi-nxinterp+1,npsi  ! go in to out, so omega will be for
     ! for outer surface at end of loop
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
     do i=1,npts
        nu(i)= (fval(j)*leng(npts))/(2.*pi*bpl(i)*rl(i)**2)
     enddo
     omega(1)=0.
     dchi=2.*pi/(npts-1.)
     do i=2,npts
        omega(i)  =omega(i-1)  +dchi*(nu(i)  +nu(i-1)  )*0.5
     enddo
     !         q_array(j)=omega(npts)/(2.*pi)
     q_plas(npsi-j+1)=omega(npts)/(2.*pi)
  enddo
  del=(q0-q_plas(1))*nn
  if ((del<0.) .or. (del>1.)) then
     write(*,*) 'problem with del=',del
     write(*,*) 'q_plas(1)=',q_plas(1),q_plas(2)
     write(*,*) 'q0=',q0
     stop
  endif
  !      write(*,*) 'newly calculated q_plas(1)=',q_plas(1)
 ! write(*,*) 'using q_plas as q_surf to be written to .eq file', &
!       q_plas(1),'del=',del
  if (verbose .gt. 3) write(outfile,*) 'newly calculated q_plas=',q_plas
  if (verbose .gt. 3) write(outfile,*) &
       'using q_plas as q_surf to be written to .eq file', &
       q_plas(1),'del=',del
  !      write(*,*) 'psigrid=',psigrid
  q_surf=q_plas(1)
  qmin=q_plas(1)
  imin=1
  do i=2,nxinterp
     if (q_plas(i) < qmin) then
        qmin=q_plas(i)
        imin=i
     endif
  enddo
!  write(*,*) 'on nxinterp grid, qmin=',qmin,'at ixinterp=',imin, &
!       'psi=',psigrid(imin)

  do i=1,npts
     omega(i)=omega(i)/q_surf
  end do
  if (abs(omega(npts)-2.*pi)>.01) then
     write(*,*) 'something is wrong, omega(npts)=',omega(npts)
     stop
  endif
  omega(npts)=2.*pi


  ! calculate rmaxgrid to go with psigrid for writing to .xtopsi
  do i=npsi,npsi-nxinterp+1,-1
     rmax=xs(i,1)
     do j=2,npts
        if (xs(i,j) > rmax) rmax=xs(i,j)
     enddo
     rmaxgrid(npsi-i+1)=rmax
  enddo

  ! write out .xtopsi file using recalculated q_plas and psigrid
  if (verbose .gt. 5) write(outfile,*) 'writing ',runname(1:lrunname)//'.xtopsi'
  open(unit=25,file=runname(1:lrunname)//'.xtopsi', &
       status='unknown')
  write(25,*) nxinterp
  write(25,*) (nn*(q0-q_plas(i)),i=1,nxinterp)
  write(25,*) psigrid
  write(25,*) rmaxgrid  ! write R_max
  write(25,*) rmaxgrid-xaxis  ! write r_max on each surface
  write(25,*) 'total flux=',psitotal
  close(25)


  ! spline onto finer grid (modeling xx grid in plas) to get more
  !   accurate qmin value 9/01

  !      dpsi=(psigrid(1)-psigrid(nxinterp))/(nxinterp*scalefact-1.d0)
  !      psifine(1)=psigrid(1)
  !      do i=2,nxinterp*scalefact
  !         psifine(i)=psifine(i-1)-dpsi
  !      enddo
  !      psifine(nxinterp*scalefact)=psigrid(nxinterp)
  !      write(*,*) 'psifine(nxinterp*scalefact)=',psifine(nxinterp*scalefact)
!!! spline routine apparently requires increasing x grid (use -psigrid)
  !      call spline1dx(q_fine,-psifine,nxinterp*scalefact, &
  !           q_plas,-psigrid,nxinterp,workq)



  call spline(-psigrid,q_plas,nxinterp,-1.d30,-1.d30,workq)
  !   qmin already set to min on xinterp grid   qmin=q_plas(1)
  psiqmin=psigrid(imin)
  scalefact=1000
  dpsi=(psigrid(1)-psigrid(nxinterp))/(nxinterp*scalefact-1.d0)
  psiwork=psigrid(1)
  do i=1,nxinterp*scalefact
     call splint(-psigrid,q_plas,workq,nxinterp,-psiwork,q_fine)
     if (q_fine .le. qmin) then
        qmin=q_fine
        psiqmin=psiwork
     endif
     if (q_fine > q_plas(1)) then
        write(*,*) 'q_fine > q_plas(1)',q_fine,psiwork
        write(*,*) 'qmax must be at the edge for now'
        write(*,*) '!!!!WARNING found q_fine > q_sep',q_fine,q_plas(1)
        if (verbose .gt. 3) write(outfile,*) &
             '!!!!WARNING found q_fine > q_sep',q_fine,q_plas(1)
!        stop
     endif
     psiwork=psiwork-dpsi
     if (psiwork < psigrid(nxinterp)) psiwork=psigrid(nxinterp)
     !         write(*,*) 'psiwork=',psiwork,'i=',i
     !         if (i == nxinterp*scalefact) write(*,*) 'psiworklow=',psiwork,dpsi
  enddo
  !      write(*,*) 'q_fine=',q_fine
  !      write(*,*) 'psifine=',psifine

  !9/01 in the future, we should just find the actual minimum of the
  !    spline fit, rather than messing with the above iteration


!  write(*,*) 'on splined grid with ',nxinterp*scalefact, &
!       'points, qmin=',qmin,' at psi=',psiqmin
  if (verbose .gt. 3) write(outfile,*) &
       'on splined grid with ',nxinterp*scalefact, &
       'points, qmin=',qmin,' at psi=',psiqmin


  if (reverse) then
     if (verbose .ge. 1) then
        write(6,*) '********************************************************'
        write(6,*) '**WARNING!! Elite works best with monotonic q          *'
        write(6,*) '** profiles.  Best to increase psimin to omit the      *'
        write(6,*) '** reversed q region if in inner core.  Not all        *'
        write(6,*) '** meshytpe values allow reversed q profiles           *'
        write(6,*) '********************************************************'
     endif
     if (verbose .gt. 3) then
        write(outfile,*) '********************************************************'
        write(outfile,*) '**WARNING!! Elite works best with monotonic q          *'
        write(outfile,*) '** profiles.  Best to increase psimin to omit the      *'
        write(outfile,*) '** reversed q region if in inner core.  Not all        *'
        write(outfile,*) '** meshytpe values allow reversed q profiles           *'
        write(outfile,*) '********************************************************'
     endif
  endif

 ! write(outfile,*) 'finished gendat'
 ! write(*,*) 'nxinterp from eliteeq is'
 ! write(*,*) nxinterp
!!$  if (setdel .and. (.not. qafix)) then
!!$     write(*,*) 'to get del=',del_fix, 'need qnew/qold='
!!$     write(*,*) 1.+(del-del_fix)/(q_nonorm*nn)
!!$  end if

!  write(*,*) 'int(1000*del) (for use in scripts) is:'
!  write(*,*) int(1000*del)
  return   

end subroutine gendat

!
!-------------------------------------------------------------------------
!
!      
subroutine testeq

! HRW 12 Nov 2002: code to calculate qprime and qprime-prime from expansion
! of G-S equation compare with numerical differential of q-profile and 
! thereby serve as a test of the equilibrium "quality"
!-----------------------------------------------------------------------
      use eliteeq_data
      implicit none
!
      real, dimension(npts) :: rl,zl,bpl,leq,leng,work
      real, dimension(npts) :: uprimel,cosul,sinul,dpsi1l,d2psi1l
      real nu(npts),nup(npts),nupp(npts),omeg(npts),omegp(npts),omegpp(npts)
      real, dimension(npts) :: dpsi1dl,d2psi1dl2,rcinvps,uprimepts,cosupts,      &
                               sinupts
      real dr,drt,dl,qsurf,qsurfp,qsurfpp
      real x1,x2,x3,denom,dchi,nu1,nu2,rc
      real drdl,dr2dl2,dzdl,dz2dl2
      real psi_1,psi_2,psi_3
      real f_loc,ffprime_loc,pprime_loc,ppp_loc
      integer ip,im,i,j,ii,ir
      real y1,y2,y3,dqdpsi,d2qdpsi,aa,bb,cc,errqp,errqpp,qpmax,qppmax,errq

      if (nxinterp.le.2) then
        write(6,*)' WARNING****nxinterp</=2, so no equilibrium error check done'
        return
      end if
      errqp=0.
      errqpp=0.
      qpmax=0.
      qppmax=0.
    do j=npsi-1,npsi-nxinterp+2,-1
        pprime_loc=4.*pi*pprime(j)
        ppp_loc=4.*pi*ppp(j)
        ffprime_loc=ffprime(j)  ! 3/1/00 use adjusted values
        f_loc=fval(j)
!  spline r,z,Bp onto equal arc length:
         leng(1)=0.
         do i=1,npts
            ir=npts+1-i
            rpts(ir)=xs(j,i)
            zpts(ir)=zs(j,i)
            bppts(ir)=bps(j,i)
         enddo
         do i=2,npts
            dr=rpts(i)-rpts(i-1)
            drt=zpts(i)-zpts(i-1)
            leng(i)=leng(i-1)+sqrt(dr**2+drt**2)
         end do
!  useful eqbm pieces...sinu, etc
!
         do i=1,npts
            if(i.eq.1) then
               x1=leng(npts-1)-leng(npts)
               x2=leng(i)
               x3=leng(i+1)
               ip=2
               im=npts-1
            else if(i.eq.npts) then
               x1=leng(i-1)
               x2=leng(i)
               x3=leng(i)+leng(2)
               ip=2
               im=npts-1
            else
               x1=leng(i-1)
               x2=leng(i)
               x3=leng(i+1)
               ip=i+1
               im=i-1
            endif
            denom=(x1-x2)*(x1-x3)*(x2-x3)
            drdl=(-rpts(ip)*(x1-x2)**2+rpts(im)*(x2-x3)**2+ &
               rpts(i)*(x1-x3)*(x1-2.*x2+x3))/denom
            dr2dl2=2.*(rpts(ip)*(x1-x2)+rpts(im)*(x2-x3)+rpts(i)*(x3-x1))/ &
               denom
            dzdl=(-zpts(ip)*(x1-x2)**2+zpts(im)*(x2-x3)**2+ &
                zpts(i)*(x1-x3)*(x1-2.*x2+x3))/denom
            dz2dl2=2.*(zpts(ip)*(x1-x2)+zpts(im)*(x2-x3)+zpts(i)*(x3-x1))/ &
                denom
            dpsi1dl(i)=(-rpts(ip)*bppts(ip)*(x1-x2)**2+rpts(im)*bppts(im) &
               *(x2-x3)**2+rpts(i)*bppts(i)*(x1-x3)*(x1-2.*x2+x3))/denom
            d2psi1dl2(i)=2.*(rpts(ip)*bppts(ip)*(x1-x2)+ &
                rpts(im)*bppts(im)*(x2-x3)+rpts(i)*bppts(i)*(x3-x1))/denom
!-----------------------------------------------------------------------
!     the following quantity is -u'+sinu/rs
!-----------------------------------------------------------------------
            rcinvps(i)=(drdl*dz2dl2-dzdl*dr2dl2)-dzdl/rpts(i)
            uprimepts(i)=-drdl*dz2dl2+dzdl*dr2dl2
            cosupts(i)=drdl
            sinupts(i)=-dzdl
         enddo
!-----------------------------------------------------------------------
!     6.0 spline r,z, and rbp onto an equal arc length mesh
!-----------------------------------------------------------------------
         dl=leng(npts)/(npts-1.)
         do i=1,npts
            leq(i)=(i-1)*dl
         end do
         leq(npts)=leng(npts)
         call spline1d(uprimel,leq,npts,uprimepts,leng,npts,work)
         call spline1d(cosul  ,leq,npts,cosupts  ,leng,npts,work)
         call spline1d(sinul  ,leq,npts,sinupts  ,leng,npts,work)
         call spline1d(dpsi1l  ,leq,npts,dpsi1dl  ,leng,npts,work)
         call spline1d(d2psi1l  ,leq,npts,d2psi1dl2  ,leng,npts,work)

!
!----------------------------------------------------------------------------------
!
         call spline1d(rl     ,leq,npts,rpts     ,leng,npts,work)
         call spline1d(zl     ,leq,npts,zpts     ,leng,npts,work)
         call spline1d(bpl    ,leq,npts,bppts    ,leng,npts,work)

!-----------------------------------------------------------------------
!     1.0 define nu and first two psi derivatives
!     explicitly extract the psi dependence of nu to facilitate
!     taking psi derivatives
!     nu=f*L/(2*pi*R^2)*(1+nuhat*dpsi)/(bp+dbpp*dpsi)
!     these values have been checked 27/10/97
!     (10/20/00 nupp not being calculated suffieciently accurately...
!     Miller expansion improved to include O(rho^2) corrections for nupp
!     (see Miller et al PoP 5 (1998) 973)
!-----------------------------------------------------------------------
      do i=1,npts
        psi_1=rl(i)*bpl(i)
        psi_2=0.5*((-rl(i)*uprimel(i)+sinul(i))*bpl(i)- &
          rl(i)**2*pprime_loc-ffprime_loc)
        psi_3=(1./6.)*((-2.*bpl(i)*sinul(i)+4.*psi_2+ffprime_loc)* &
          (sinul(i)/rl(i)-uprimel(i))-rl(i)**2*pprime_loc* &
          (sinul(i)/rl(i)+uprimel(i))-d2psi1l(i)+cosul(i)* &
          dpsi1l(i)/rl(i)-rl(i)*bpl(i)*(rl(i)**2*ppp_loc &
          +ffpp(j)))
        Rc=-1./uprimel(i)
        nu(i)= (f_loc*leng(npts))/(2.*pi*bpl(i)*rl(i)**2)
        nu1=-sinul(i)/rl(i)+uprimel(i)-2.*psi_2/psi_1
! 11/28/00 eliminate the dpsi1l term which cancels when the higher
!   order correction to dl_p is properly considered in the nu'' calculation
        nu2=(sinul(i)/rl(i)+2.*psi_2/psi_1)*(sinul(i)/rl(i)-uprimel(i)) &
            +4.*psi_2**2/psi_1**2-3.*psi_3/psi_1
        nup(i)=nu(i)*(ffprime_loc/f_loc**2+nu1/psi_1)
        nupp(i)=nu(i)*((ffpp(j)-(ffprime_loc/f_loc)**2)/ &
          f_loc**2 +2.*ffprime_loc*nu1/(psi_1*f_loc**2)+ &
          2.*(nu2-nu1*psi_2/psi_1)/psi_1**2)
      end do
!-----------------------------------------------------------------------
!     1.1  can use nu to calculate q's and omega's
!     checked 27/10/97
!-----------------------------------------------------------------------
      omeg(1)=0.
      omegp(1)=0.
      omegpp(1)=0.
      dchi=2.*pi/(npts-1.)
      do i=2,npts
         omeg(i)  =omeg(i-1)  +dchi*(nu(i)  +nu(i-1)  )*0.5
         omegp(i) =omegp(i-1) +dchi*(nup(i) +nup(i-1) )*0.5
         omegpp(i)=omegpp(i-1)+dchi*(nupp(i)+nupp(i-1))*0.5
      end do
      qsurf=omeg(npts)/(2.*pi)
      qsurfp=omegp(npts)/(2.*pi)
      qsurfpp=omegpp(npts)/(2.*pi)
      if (verbose .ge. 5) then
         write(6,*)' surface=',j,' psi=',psiv(j),'psinorm=', &
             percenflux*(psiv(j)-psiv(1))/(psiv(npsi)-psiv(1))
         write(6,*)' q_num=',qcalc(j),' q_test=',qsurf, &
              'q_eq=',qsfin(j),'errq=',abs(qsurf-qsfin(j))/qsfin(j)
      endif
      ip=j+1
      if (ip.gt.npsi) ip=npsi
      ii=ip-1
      im=ip-2
      x1=psiv(im)
      x2=psiv(ii)
      x3=psiv(ip)
      y1=qcalc(im)
      y2=qcalc(ii)
      y3=qcalc(ip)
      aa=((y1-y2)/(x1-x2)-(y2-y3)/(x2-x3))/(x1-x3)
      bb=(y1-y2)/(x1-x2)-aa*(x1+x2)
      cc=y1-aa*x1**2-bb*x1
      if(i.ne.npsi) then
         dqdpsi=2.*aa*x2+bb
      else
         dqdpsi=2.*aa*x3+bb
      endif
      d2qdpsi=2.*aa
      if (verbose .ge. 5) then
         write(6,*)' GS qp=',qsurfp,' num qp=',dqdpsi
         write(6,*)' GSqpp=',qsurfpp,' num qpp=',d2qdpsi
         write(6,*) 'errqp=',abs((qsurfp-dqdpsi)/(qsurfp+dqdpsi))
         write(6,*) 'errqpp=',abs((qsurfpp-d2qdpsi)/(qsurfpp+d2qdpsi))
      endif
      
      errqp=errqp+abs((qsurfp-dqdpsi)/(qsurfp+dqdpsi))
      errqpp=errqpp+abs((qsurfpp-d2qdpsi)/(qsurfpp+d2qdpsi))
      errq=abs((qsurfp-dqdpsi)/(qsurfp+dqdpsi))
      if (errq.gt.qpmax) qpmax=errq
      errq=abs((qsurfpp-d2qdpsi)/(qsurfpp+d2qdpsi))
      if (errq.gt.qppmax) qppmax=errq
  end do
  errqp=errqp/(nxinterp-1)
  errqpp=errqpp/(nxinterp-1)
  if (verbose .ge. 1) then
     write(6,*)' '
     write(6,*)'**********************************************************'
     write(6,*)' Following errors indicate quality of equilibrium:'
     write(6,*)' Average error in q-prime=',errqp
     write(6,*)' Maximum error in q-prime=',qpmax
     write(6,*)' If these errors are greater than few percent, consider'
     write(6,*)' converging equilibrium to higher accuracy' 
     write(6,*)' Average error in q-prime-prime=',errqpp
     write(6,*)'*********************************************************'
  endif
  if (verbose .gt. 3) then
     write(outfile,*)' '
     write(outfile,*)'**********************************************************'
     write(outfile,*)' Following errors indicate quality of equilibrium:'
     write(outfile,*)' Average error in q-prime=',errqp
     write(outfile,*)' Maximum error in q-prime=',qpmax
     write(outfile,*)' If these errors are greater than few percent, consider'
     write(outfile,*)' converging equilibrium to higher accuracy' 
     write(outfile,*)' Average error in q-prime-prime=',errqpp
     write(outfile,*)'*********************************************************'
  endif

!  write(6,*)' Maximum error in q-prime-prime=',qppmax
!  write(*,*) 'int(1000*del) (for use in scripts) is:'
!  write(*,*) int(1000*del)

end subroutine testeq
!
!************************************************************************
!
      subroutine qspline(ipsi,qedge)
!     ******************************
!
!  Evaluates q at edge by splining R and Bp onto equal arc length as
!  done in plas. This leads to a consistent del between equil and plas.

! pbs 5/05 edit to avoid redefinition of the global rpts, zpts, bppts here

      use eliteeq_data
      implicit none
      real leng(npts),work(npts)
      real rl(ns),bpl(ns),leq(ns),nu(ns)
      integer i,ipsi
      real dr,drt,dl,qedge,dchi
      real rpts_loc(npts),zpts_loc(npts),bppts_loc(npts)

!
      do i=1,npts
        rpts_loc(i)=xs(ipsi,i)
        zpts_loc(i)=zs(ipsi,i)
        bppts_loc(i)=bps(ipsi,i)
      end do
      leng(1)=0.
      do i=2,npts
        dr=rpts_loc(i)-rpts_loc(i-1)
        drt=zpts_loc(i)-zpts_loc(i-1)
        leng(i)=leng(i-1)+sqrt(dr**2+drt**2)
      end do
      dl=leng(npts)/(ns-1.)
      do i=1,ns
         leq(i)=(i-1)*dl
      end do
      leq(ns)=leng(npts)
      call spline1d(rl     ,leq,ns,rpts_loc     ,leng,npts,work)
      call spline1d(bpl    ,leq,ns,bppts_loc    ,leng,npts,work)
      do i=1,ns
         nu(i)= (fval(ipsi)*leng(npts))/(2.*pi*bpl(i)*rl(i)**2)
      end do
!  Evaluate qedge
      qedge=0.
      dchi=2.*pi/(ns-1.)
      do i=2,ns
         qedge  =qedge+dchi*(nu(i)  +nu(i-1)  )*0.5
      end do
      qedge=qedge/(2.*pi)
end subroutine qspline
