subroutine readpest

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!read pest style equilibrium files
! 7/02 pbs based on routine by turnbull
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use eliteeq_data
  implicit none

  integer noutpest,ndskpest,ios
  integer old_npts

  character*8  etitl(8),date
  real bfieldf,psifact,prsfact,amu
  real, dimension(:), allocatable :: sg1,sg2,sg3,sg4,sqvlp,dum2, &
       psic_loc
  real, dimension(:,:), allocatable :: seqaj0,seqaj3
  integer nft,ntadjst,npadjst,it
  real rndoff,rndofbf,rndofps,rndofpr,rndofr
  integer nthd1,npsd1,ntht
  integer ntt,npp,ntinpt,npinpt,nremt,nremp
  real dlr,dlt,sdummy
  integer ntbegin,ntlast,ntlst2,ntroom,npbegin,nplast,nproom
  integer neqsym,jj,ii,kk,ll,iii,i1,i0,i,j
  integer kmap,nprofl,nthet,nthp1
  real radfact,btffact,psimaxa,twopi,psilima,pslimdf
  real rcnt,prnorm,psilimp,psiminp,psimax,psmaxdf
  real fbnorm,spj1,spj2,sfj1,sfj2,sgj1,sgj2,rrps1,rrps2,rrps3
  real zzps1,zzps2,zzps3,a3ps1,a3ps2,a3ps3,a0ps1,a0ps2,a0ps3
  real xma,zma,seqraxs,seqzaxs,totcur,sfedge,btor,psisep
  real xsep,zsep
  real xax(2),zax(2),psimx(2),seqrdif,seqzdif
  integer jjp,iip,ibsign,itsign,npst
  integer ircnt,izcnt,ikr0,ikr1,ikz0,ikz1
  real suminr,suminz,dscinr,dscinz
  real dpsi_loc,dchi,bnd_set,check,dum
  real convpsi,convpp,convf,convffp,convchi,convm,convb
  parameter (convpsi=1.e8,convpp=1.e-7,convf=1.e6,convffp=1.e4, &
       convchi=convpsi,convm=1.e2,convb=1.e4)
  parameter (nft=5)

  noutpest=12
  ndskpest=22
  old_npts=npts

  write(*,*) 'lrunname=',lrunname,'runname=',runname
  open(unit=ndskpest,file=runname(1:lrunname)//'.dskpest', &
       status='old',iostat=ios)
  if(ios==0) write(*,*) 'reading dskpest file ', &
       runname(1:lrunname)//'.dskpest'
  if(ios.ne.0) then
     write(6,*) 'problem opening file ',runname(1:lrunname)//'.dskpest'
     write(6,*) 'will try to read file named dskpest instead'
     open(unit=ndskpest,file='dskpest',status='old',iostat=ios)
     if(ios.ne.0) then
        write(*,*) 'could not open dskpest file'
        stop
     endif
     write(*,*) 'reading equilibrium from file dskpest'
  endif


! 1.0 Initialization
! 1.1 Set the conversion factors
  amu      = 4.0*pi*1.0e-07
  twopi=2.*pi
  bfieldf   = 1.0
  psifact   = bfieldf
  prsfact   = (bfieldf*bfieldf)/amu
  kmap=0

! 1.2 Set roundoff parameters

  rndoff    = 1.e-12
  rndofbf   = rndoff*bfieldf
  rndofps   = rndoff*psifact
  rndofpr   = rndoff*prsfact

! 1.3 Set the adjustment in labeling used in JSOLVER
  ntadjst   = 2
  npadjst   = 0

! 1.4 Set the default heading and date
  date      = 'unknown'
  do it  = 1,nft
     etitl(it) = ' '
  enddo

! 2.0 Read in the inverse equilibrium scalar data
! 2.1 Read and check the dimensions
! 2.1.1 Read the input values
  read (ndskpest,1000) nthd1,npsd1,ntht,npsi,neqsym,dlr,dlt

! 2.1.2 Set the theta dimension according to up-down symmetry
! 2.1.2.1 Check for invalid symmetry switch
  if((neqsym .lt. 0) .or. (neqsym .gt.  1 )) then
     write(*,*) 'neqsym is not 0 or 1!!!!!!!'
     if(neqsym .lt. 0) neqsym  = 0
     if(neqsym .gt. 1) neqsym  = 1
  endif

! 2.1.2.2 Set the full theta array size
  if(neqsym .eq. 0) npts=ntht
  if(neqsym .ne. 0) npts=2*(ntht-1) + 1

  if (npts.ne.old_npts) then
     write(*,*) 'npts in input file must =npts from dskpest',npts,old_npts
     stop
  endif

  write(*,*) 'npsi=',npsi,' npts=',npts,' neqsym=',neqsym
  if (verbose .gt. 4) write(outfile,*) 'allocating toq variables'

  call alloc_toq
  if (verbose .gt. 4) write(outfile,*) &
       'deallocate and reallocate surf arrays with new npts'
  call dealloc_arrays
  call alloc_arrays
  allocate( sp(npsi), sg1(npsi), sg2(npsi), sg3(npsi), sg4(npsi) )
  allocate( sqvlp(npsi), dum2(npsi), psic_loc(npsi) )
  allocate( seqaj0(npsi,npts), seqaj3(npsi,npts) )
  write(*,*) 'finished allocating'
! 2.1.3 Set sizes for excess zero values included in the input file
! 2.1.3.1 Set the size differences

! ntt and npp are array sizes, here they're dynamic so...
  npp=npsi
  ntt=npts

  ntinpt  = nthd1 - ntht
  npinpt  = npsd1 - npsi
  nremt   = nthd1 - ntt
  nremp   = npsd1 - npp

! 2.1.3.2 Set the last valid theta index from the JSOLVER indexing
  ntbegin = ntadjst + 1
  ntlast  = ntht + ntbegin
  ntlst2  = npts + ntbegin
  ntroom  = ntt  - ntlst2
  npbegin = npadjst + 1
  nplast  = npsi + npbegin - 1
  nproom  = npp  - nplast

  write(*,*) 'nremp=',nremp,' nproom=',nproom

! 2.1.3.3 Check for invalid input dimensions
  if( ntinpt .lt. ntadjst) then
     write(*,*) 'ntinp < ntadjst'
     stop
  endif
  if( npinpt .lt. npadjst) then
     write(*,*) 'npinpt < npadjst'
  endif

! 2.2 Read the equilibrium scalar data
  write(*,*) 'read rcnt etc'
  read (ndskpest,1100) rcnt,prnorm,psilimp,psiminp
  write(*,*) 'rcnt etc read'

! 3.0 Read the equilibrium quantities
! 3.1 Read the profiles
! 3.1.1 Read profiles up to the full input dimension written if the
!       available dimension is sufficient

  if    (nremp .le. 0) then
     read (ndskpest,2000)  (sp     (jj), jj = 1,npsd1)
     read (ndskpest,2000)  (pprime    (jj), jj = 1,npsd1)
     read (ndskpest,2000)  (qsfin   (jj), jj = 1,npsd1)
     read (ndskpest,2000)  (sqvlp  (jj), jj = 1,npsd1)
     read (ndskpest,2000)  (fval     (jj), jj = 1,npsd1)
     read (ndskpest,2000)  (ffprime   (jj), jj = 1,npsd1)
     read (ndskpest,2000)  (sg1    (jj), jj = 1,npsd1)
     read (ndskpest,2000)  (sg2    (jj), jj = 1,npsd1)
     read (ndskpest,2000)  (sg3    (jj), jj = 1,npsd1)
     read (ndskpest,2000)  (sg4    (jj), jj = 1,npsd1)
     read (ndskpest,2000)  (psiv (jj), jj = 1,npsd1)

! 3.1.2 Read the profiles up to the available dimension and discard
!       the rest if the available dimension is insufficient
  elseif(nremp .gt. 0) then
     read (ndskpest,2000)  (sp     (jj), jj = 1,npp  ) &
          ,(sdummy     , kk = 1,nremp)
     write(*,*) 'read sp=',sp,'sdummy=',sdummy
     read (ndskpest,2000)  (pprime    (jj), jj = 1,npp  ) &
          ,(sdummy     , kk = 1,nremp)
     write(*,*) 'read pprime=',pprime,'sdummy=',sdummy
     read (ndskpest,2000)  (qsfin   (jj), jj = 1,npp  ) &
          ,(sdummy     , kk = 1,nremp)
     read (ndskpest,2000)  (sqvlp  (jj), jj = 1,npp  ) &
          ,(sdummy     , kk = 1,nremp)
     read (ndskpest,2000)  (fval     (jj), jj = 1,npp  ) &
          ,(sdummy     , kk = 1,nremp)
     read (ndskpest,2000)  (ffprime   (jj), jj = 1,npp  ) &
          ,(sdummy     , kk = 1,nremp)
     read (ndskpest,2000)  (sg1    (jj), jj = 1,npp  ) &
          ,(sdummy     , kk = 1,nremp)
     read (ndskpest,2000)  (sg2    (jj), jj = 1,npp  ) &
          ,(sdummy     , kk = 1,nremp)
     read (ndskpest,2000)  (sg3    (jj), jj = 1,npp  ) &
          ,(sdummy     , kk = 1,nremp)
     read (ndskpest,2000)  (sg4    (jj), jj = 1,npp  ) &
          ,(sdummy     , kk = 1,nremp)
     read (ndskpest,2000)  (psiv (jj), jj = 1,npp  ) &
          ,(sdummy     , kk = 1,nremp)
  endif
  write(*,*) 'eq arrays read'


! 3.2 Read the inverse mapping r(psii,theta), and z(psi,theta)

! 3.2.1 Read the inverse equilibrium (r,z)
! 3.2.1.1 Read the mapping up to the full input dimensions written if
!         the available dimensions are sufficient
  if    (nremp .le. 0  .and.  nremt .le. 0) then
     read (ndskpest,2000) ((xs(jj,ii), ii = 1,nthd1) &
          , jj = 1,npsd1)
     read (ndskpest,2000) ((zs(jj,ii), ii = 1,nthd1) &
          , jj = 1,npsd1)

! 3.2.1.2 Read the mapping up to the available dimension and discard
!         the rest if the available dimension in psi is insufficient
  elseif(nremp .gt. 0  .and.  nremt .le. 0) then
     read (ndskpest,2000) ((xs(jj,ii), ii = 1,nthd1) &
          , jj = 1,npp  ) &
          ,((sdummy       , ii = 1,nthd1) &
          , kk = 1,nremp)
     read (ndskpest,2000) ((zs(jj,ii), ii = 1,nthd1) &
          , jj = 1,npp  ) &
          ,((sdummy       , ii = 1,nthd1) &
          , kk = 1,nremp)

! 3.2.1.3 Read the mapping up to the available dimension and discard
!         the rest if the available dimension in theta is insufficient
  elseif(nremp .le. 0  .and.  nremt .gt. 0) then
     read (ndskpest,2000) ((xs(jj,ii), ii = 1,ntt  ) &
          , (sdummy       , ll = 1,nremt) &
          , jj = 1,npsd1)
     read (ndskpest,2000) ((zs(jj,ii), ii = 1,ntt  ) &
          , (sdummy       , ll = 1,nremt) &
          , jj = 1,npsd1)

! 3.2.1.4 Read the mapping up to the available dimension and discard
!         the rest if both available dimensions are insufficient
  elseif(nremp .gt. 0  .and.  nremt .gt. 0) then
     read (ndskpest,2000) ((xs(jj,ii), ii = 1,ntt  ) &
          , (sdummy       , ll = 1,nremt) &
          , jj = 1,npp  ) &
          ,((sdummy       , ii = 1,ntt  ) &
          , (sdummy       , ll = 1,nremt) &
          , kk = 1,nremp)
     read (ndskpest,2000) ((zs(jj,ii), ii = 1,ntt  ) &
          , (sdummy       , ll = 1,nremt) &
          , jj = 1,npp  ) &
          ,((sdummy       , ii = 1,ntt  ) &
          , (sdummy       , ll = 1,nremt) &
          , kk = 1,nremp) 
  endif


  write(*,*) 'xs and zs read'

! 3.2.2 Read all the Jacobian for kmap = -2
  if(kmap .eq. -2) then

! 3.2.2.1 Read the Jacobian up to the full input dimensions written if
!         the available dimensions are sufficient
     if    (nremp .le. 0  .and.  nremt .le. 0) then
        read (ndskpest,2000) ((seqaj3(jj,ii), ii = 1,ntt  ) &
             , jj = 1,npp  )
        read (ndskpest,2000) ((seqaj0(jj,ii), ii = 1,ntt  ) &
             , jj = 1,npp  )

! 3.2.2.2 Read the Jacobian up to the available dimension and discard
!         the rest if the available dimension in psi is insufficient
     elseif(nremp .gt. 0  .and.  nremt .le. 0) then
        read (ndskpest,2000) ((seqaj3(jj,ii), ii = 1,ntt  ) &
             , jj = 1,npp  ) &
             ,((sdummy       , ii = 1,ntt  ) &
             , kk = 1,nremp)
        read (ndskpest,2000) ((seqaj0(jj,ii), ii = 1,ntt  ) &
             , jj = 1,npp  ) &
             ,((sdummy       , ii = 1,ntt  ) &
             , kk = 1,nremp)

! 3.2.2.3 Read the Jacobian up to the available dimension and discard
!         the rest if the available dimension in theta is insufficient
     elseif(nremp .le. 0  .and.  nremt .gt. 0) then
        read (ndskpest,2000) ((seqaj3(jj,ii), ii = 1,ntt  ) &
             , (sdummy       , ll = 1,nremt) &
             , jj = 1,npp  )
        read (ndskpest,2000) ((seqaj0(jj,ii), ii = 1,ntt  ) &
             , (sdummy       , ll = 1,nremt) &
             , jj = 1,npp  )
     
! 3.2.2.4 Read the Jacobian up to the available dimension and discard
!         the rest if both available dimensions are insufficient
     elseif(nremp .gt. 0  .and.  nremt .gt. 0) then
        read (ndskpest,2000) ((seqaj3(jj,ii), ii = 1,ntt  ) &
             , (sdummy       , ll = 1,nremt) &
             , jj = 1,npp  ) &
             ,((sdummy       , ii = 1,ntt  ) &
             , (sdummy       , ll = 1,nremt) &
             , kk = 1,nremp)
        read (ndskpest,2000) ((seqaj0(jj,ii), ii = 1,ntt  ) &
             , (sdummy       , ll = 1,nremt) &
             , jj = 1,npp  ) &
             ,((sdummy       , ii = 1,ntt  ) &
             , (sdummy       , ll = 1,nremt) &
             , kk = 1,nremp)
     endif
  endif


! 3.3 Set the quantities that were not explicitly read in

! 3.3.1 Set the basic integer parameters

  nprofl  = npsi
  nthet   = npts
  nthp1   = ntht + 1

! 4.0 Modify the units
! 4.1 Conversion factors

  radfact   = rcnt
  btffact   = radfact*bfieldf

! 4.2 Reset the poloidal flux to Wb/m**2
! 4.2.1 Reset the end values
  psimaxa = bfieldf*psiminp/twopi
  psilima = bfieldf*psilimp/twopi

! 4.2.2 Reset the flux array

  do jj   = 1,nprofl
     psiv(jj)  = psifact*psiv(jj)
  enddo

! 4.2.3 Check the psi values are consistent
  psimax  = psiv(1)
  psilim  = psiv(npsi)
  psmaxdf = abs(psimax - psimaxa)
  pslimdf = abs(psilim - psilima)

  if(psmaxdf .ge. rndofps) then
     write(*,*) 'psmaxdf .ge. rndofps',psmaxdf,rndofps
  endif

  if(pslimdf .ge. rndofps) then
     write(*,*) 'pslimdf .ge. rndofps!!!',pslimdf,rndofps
  endif

! 4.3 Reset the pressure and toroidal field units and shift the pressure
!     and toroidal field to the grid points

! 4.3.1 Reset the pressure and pprime units
  prnorm      = prsfact*prnorm

  do jj   = 1,nprofl
     sp   (jj)   =  prsfact         *sp   (jj)
     pprime  (jj)   = (prsfact/psifact)*pprime  (jj)
  enddo

! 4.3.2 Reset the toroidal field function and units
  fbnorm      = 1.0*btffact

  do jj   = 1,nprofl
     fval   (jj)   =  btffact         *fval   (jj)
     ffprime (jj)   = (btffact/psifact)*ffprime (jj)
  enddo

! 4.3.3 Reset the remaining input profile units

  do jj   = 1,nprofl
     sqvlp(jj)   = (1.0/psifact)              *sqvlp(jj)
     sg1  (jj)   = (btffact/radfact)          *sg1  (jj)
     sg2  (jj)   = (btffact/(psifact*radfact))*sg2  (jj)
     sg3  (jj)   =  psifact                   *sg3  (jj)
     sg4  (jj)   =  1.0                       *sg4  (jj)
   enddo

! 4.3.4 Shift the pressure, toroidal field function, and sg3 
!       half a grid call

   do jj   = 2,nprofl
      jjp         = nprofl - jj + 2
      spj1        = sp (jjp-1)
      spj2        = sp ( jjp )
      sfj1        = fval (jjp-1)
      sfj2        = fval ( jjp )
      sgj1        = sg3(jjp-1)
      sgj2        = sg3( jjp )
      sp (jjp)    = 0.5*(spj1 + spj2)
      fval (jjp)    = 0.5*(sfj1 + sfj2)
      sg3(jjp)    = 0.5*(sgj1 + sgj2)
   enddo

   sp(1)       = prnorm
   fval(1)       = sg1(1)*qsfin(1)

! 4.3.5 Reset the ffprime to f*fprime

   do jj   = 1,nprofl
      ffprime (jj)   = fval   (jj)*ffprime (jj)
   enddo

! 4.4 Readjust the storage of the 2D mapping quantities

! 4.4.1 Readjust the storage of the mapping quantities xs and zs

   do jj     = 1,nprofl

! 4.4.1.1 Store the excess points
      rrps1         = xs(jj,1)
      rrps2         = xs(jj,2)
      rrps3         = xs(jj,ntlast)
      zzps1         = zs(jj,1)
      zzps2         = zs(jj,2)
      zzps3         = zs(jj,ntlast)

! 4.4.1.2 Shift the index
      do ii     = 1,ntht
         iip           = ii + ntadjst
         xs(jj,ii) = xs(jj,iip)
         zs(jj,ii) = zs(jj,iip)
      enddo

! 4.4.1.3 Add the excess points to the end of the array if space exists
      if(ntroom .ge. 0) then
         xs(jj,ntht+1) = rrps3
         xs(jj,ntht+2) = rrps2
         xs(jj,ntlast) = rrps1
         zs(jj,ntht+1) = zzps3
         zs(jj,ntht+2) = zzps2
         zs(jj,ntlast) = zzps1
      endif
   enddo

! 4.4.2 Readjust the storage of the Jacobian quantities

   if(kmap .eq. -2) then
      do jj     = 1,nprofl

! 4.4.2.1 Store the excess points
         a3ps1         = seqaj3(jj,1)
         a3ps2         = seqaj3(jj,2)
         a3ps3         = seqaj3(jj,ntlast)
         a0ps1         = seqaj0(jj,1)
         a0ps2         = seqaj0(jj,2)
         a0ps3         = seqaj0(jj,ntlast)

! 4.4.2.2 Shift the index
         do ii     = 1,ntht
            iip           = ii + ntadjst
            seqaj3(jj,ii) = seqaj3(jj,iip)
            seqaj0(jj,ii) = seqaj0(jj,iip)
         enddo

! 4.4.2.3 Add the excess points to the end of the array if space exists
         if(ntroom .ge. 0) then
            seqaj3(jj,ntht+1) = a3ps3
            seqaj3(jj,ntht+2) = a3ps2
            seqaj3(jj,ntlast) = a3ps1
            seqaj0(jj,ntht+1) = a0ps3
            seqaj0(jj,ntht+2) = a0ps2
            seqaj0(jj,ntlast) = a0ps1
         endif
      enddo
   endif



! 4.4.4.1 Compute the average of xs and zs on axis
   seqraxs   = 0.0
   seqzaxs   = 0.0
   rndofr    = rndoff*rcnt
   do ii = 1,ntht
      seqraxs   = seqraxs + xs(1,ii)
      seqzaxs   = seqzaxs + zs(1,ii)
   enddo

   seqraxs   = seqraxs/float(ntht)
   seqzaxs   = seqzaxs/float(ntht)

! 3.3.2.2 Check for discrepancies
   ircnt     = 0
   izcnt     = 0
   do ii = 1,ntht
      seqrdif   = abs(xs(1,ii) - seqraxs)
      seqzdif   = abs(zs(1,ii) - seqzaxs)
      if(seqrdif .ge. rndofr) ircnt  = ircnt + 1
      if(seqzdif .ge. rndofr) izcnt  = izcnt + 1
   enddo

   if(ircnt .gt. 0) then
      write(*,*) 'ircnt > 0!!'
      write(*,3100) ntht
      write(*,3110) (xs(1,ii),ii=1,ntht)
   endif

   if(izcnt .gt. 0) then
      write(*,*) 'izcnt > 0!!'
      write(*,3200) ntht
      write(*,3210) (zs(1,ii),ii=1,ntht)
   endif


! 4.5.1 Set the axis and limiter positions
   xma     = seqraxs
   zma     = seqzaxs
   write(*,*) 'xma=',xma,'zma=',zma
!   xlim    = xs(npsi,1)
!   write(*,*) 'xlim=',xlim
!   zlim    = zs(npsi,1)

! 4.5.2 Set the toroidal field and current

! 4.5.2.1 Set the total current
!         This needs to be set properly
   totcur  = bfieldf

! 4.5.2.2 Set sfedge and btor
   sfedge  = fval(nprofl)
   btor    = sfedge/rcnt

! 5.0 Reset the conventions if nonstandard
! 5.1 Reverse signs of psiv and psi derivatives if input has
!     the alternate psi convention

   if(psiv(1) .gt. psiv(npsi)) then

! 5.1.1 Print a warning
      write(*,*) 'reversing convention for psi',kmap,npsi

! 5.1.2 Reverse the signs
      do jj  = 1,npsi
         psiv(jj) = -psiv(jj)
         ffprime  (jj) = -ffprime  (jj)
         pprime   (jj) = -pprime   (jj)
      enddo
   endif


! 5.2 Check the sign convention of btor

! 5.2.1 Check the signs of btor and totcur

   if(btor   .lt. 0.0) ibsign  = -1
   if(btor   .eq. 0.0) ibsign  =  0
   if(btor   .gt. 0.0) ibsign  = +1

! 5.2.2 Print a warning if btor is zero


   if(ibsign .eq. 0) write(*,*) 'sign of btor is zero'
   if(totcur .lt. 0.0) itsign  = -1
   if(totcur .eq. 0.0) itsign  =  0
   if(totcur .gt. 0.0) itsign  = +1

! 5.2.3 Switch the sign of totcur if the convention is reversed

   if(ibsign .lt. 0  .and.  itsign .lt. 0) then
      write(*,*) 'rdjsvmp, btor conventions switched'
      totcur     = -totcur
   endif

! 6.0 Set up the auxiliary quantities

! 6.1 Set scalar quantities

! 6.1.1 Set the profile switch npst

      npst     = 0

! 6.1.2 Set the magnetic axis values

      psisep   = psimax
      xsep     = 0.0
      zsep     = 0.0

      xax(1)   = xma
      xax(2)   = 0.0
      zax(1)   = zma
      zax(2)   = 0.0
      psimx(1) = psimax
      psimx(2) = 0.0


! 7.0 Set the number of poloidal angles and complete arrays
!     for up-down symmetry

! 7.1 Fill the r and z arrays for up-down symmetry
      if    (neqsym .eq. 1) then

! 7.1.1 Ensure the magnetic axis is on the midplane

! 7.1.1.1 Check for discrepancy

        if(abs(zma) .gt. 0.0) then
          if(abs(zma) .le. rndofr) then
             write(*,*) 'zma not zero for neqsym=1'
          endif
          if(abs(zma) .gt. rndofr) then
             write(*,*) 'zma not zero for neqsym=1'
             write(*,4000) xma,zma
          endif
       endif
! 7.1.1.2 Reset the axis on the midplane

        zma           = 0.0


! 7.2 Ensure the first flux surface is at the magnetic axis

! 7.2.1 Check for discrepancies

        ikr0          = 0
        ikr1          = 0
        ikz0          = 0
        ikz1          = 0
        suminr        = 0.0
        suminz        = 0.0
        do ii     = 1,ntht

! 7.2.1.1 Check the r value is the axis value
           dscinr        = xs(1,ii) - xma
           suminr        = suminr + abs(dscinr)
           if(abs(dscinr) .ne. 0.0) then
              if(abs(dscinr) .le. rndofr) ikr0   = ikr0 + 1
              if(abs(dscinr) .gt. rndofr) ikr1   = ikr1 + 1
           endif

! 7.2.1.2 Check the z value is the axis value
           dscinz        = zs(1,ii) - zma
           suminz        = suminz + abs(dscinz)
           if(abs(dscinz) .ne. 0.0) then
              if(abs(dscinz) .le. rndofr) ikz0   = ikz0 + 1
              if(abs(dscinz) .gt. rndofr) ikz1   = ikz1 + 1
           endif

        enddo

! 7.2.1.3 Find the mean discrepancy
        suminr    = suminr / float(ntht)
        suminz    = suminz / float(ntht)

! 7.2.2 Print error message if any discrepancy

        if(ikr0 .gt. 0  .or.  ikr1 .gt. 0  .or.  ikz0 .gt. 0  .or. &
             ikz1 .gt. 0) then
           if(ikr0 .gt. 0) write(*,*) &
                'Small Discrepancy in xs(1)  '
           if(ikz0 .gt. 0) write(*,*) &
                'Small Discrepancy in zs(1)  '
           if(ikr1 .gt. 0) write(*,*) &
                'Large Discrepancy in xs(1)  '
           if(ikz1 .gt. 0) write(*,*) &
               'Large Discrepancy in zs(1)  '
   
           write(*,4100) npsi,ntht,suminr,suminz
           write(*,4110) (iii, xs(1,iii), zs(1,iii), &
                iii=1,ntht)
        endif

! 7.2.3 Reset the first flux surface identically to the axis values

        do ii    = 1,ntht
           xs(1,ii) = xma
           zs(1,ii) = zma
        enddo

! 7.3.3 Reset the rays along the midplane

        do jj       = 1,npsi
           zs(jj,  1 ) = 0.0
           zs(jj,ntht) = 0.0
        enddo

! 7.4 Fill the complete r and z arrays

        do ii     = nthp1,npts
           i1            = ii
           i0            = npts + 1 - i1
           
           do jj     = 1,npsi
              xs(jj,i1) = +xs(jj,i0)
              zs(jj,i1) = -zs(jj,i0)
           enddo
        enddo


! 7.5 Fill the complete Jacobian arrays

        if(kmap .eq. -2) then
           do ii     = nthp1,npts
              i1            = ii
              i0            = npts + 1 - i1

              do jj     = 1,npsi
                 seqaj0(jj,i1) = +seqaj0(jj,i0)
                 seqaj3(jj,i1) = +seqaj3(jj,i0)
              enddo
           enddo
        endif
     endif

! Need to calculate chipsi= d chi/d psi in Bob's toq notation
!      d psiv/d psic_loc
     dpsi_loc=1./(npsi-1.)
     dchi=psiv(npsi)-psiv(1)
     write(*,*) 'dpsi_loc=',dpsi_loc,' dchi=',dchi
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
            write(*,*) 'WARNING!! found negative chipsi for surface=' &
                 ,1,' reset to 0'
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

      do j=1,npts
         do i=1,npsi
            xs(i,j)=convm*xs(i,j)
            zs(i,j)=convm*zs(i,j)
         enddo
      enddo

     write(*,*) 'npsi=',npsi,'npts=',npts
     write(*,*) 'q=',qsfin
     write(*,*) 'pprime=',pprime
     write(*,*) 'chipsi=',chipsi


     close(ndskpest)


! 8.0 Return and end

      return

1000  format(5(1x,i4),2(1x,e11.4))
1100  format(4e16.8)
2000  format(5e16.8)
3000  format(/,10x,'psimax  = ',e20.13,4x,'psimaxa = ',e20.13 &
           ,/,10x,'psmaxdf = ',e16.9,8x, 'rndofps = ',e16.9)
3010  format(/,10x,'psilim  = ',e20.13,4x,'psilima = ',e20.13 &
           ,/,10x,'pslimdf = ',e16.9,8x, 'rndofps = ',e16.9)
3100  format(/,5x,'xs(1,ii),ii=1,',i5)
3110  format(5(1x,e16.9))
3200  format(/,5x,'zs(1,ii),ii=1,',i5)
3210  format(5(1x,e16.9))
3500  format(/,10x,'dpsi01  = ',e16.9,4x, 'dpsi02  = ',e16.9 &
           , 4x,'dpsi03  = ',e16.9,/ &
           ,10x,'dpsi12  = ',e16.9,4x, 'dpsi13  = ',e16.9 &
           , 4x,'dpsi23  = ',e16.9)
4000  format(/,10x,'xma     = ',e16.9, 4x,'zma     = ',e16.9)
4100  format(/,10x,'npsi    = ',i5,   15x,'ntht    = ',i5,/ &
           ,10x,'suminr  = ',e16.9, 4x,'suminz  = ',e16.9,// &
           ,4x,'ii',6x,'xs(1,ii)',5x,'zs(1,ii)')
4110  format((1x,i5,2x,2(1x,e16.9)))
4200  format(/,10x,'npsi    = ',i5,   15x,'ntht    = ',i5,/ &
           ,10x,'sumray1 = ',e16.9, 4x,'sumrayn = ',e16.9,// &
           ,4x,'jj',6x,'xs(jj,  1 )',2x,'zs(jj,  1 )' &
           ,2x,'xs(jj,ntht)',2x,'zs(jj,ntht)')
4210  format((1x,i5,2x,4(1x,e16.9)))


end subroutine readpest
