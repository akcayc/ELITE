subroutine readpest

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!read pest style equilibrium files
! 7/02 pbs based on routine by turnbull
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use eliteeq_data
  implicit none

  integer noutpest,ndskpest,ios
  integer old_nptssubroutine readpest

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


  real dpsi_loc,dchi,bnd_set,check,dum,bfieldf_in
  real convpsi,convpp,convf,convffp,convchi,convm,convb,dum5(100)
  logical there
  character*12 dumst(6)

  parameter (convpsi=1.e8,convpp=1.e-7,convf=1.e6,convffp=1.e4, &
       convchi=convpsi,convm=1.e2,convb=1.e4)
  parameter (nft=5)

  noutpest=12
  ndskpest=22
  old_npts=npts

  write(*,*) 'lrunname=',lrunname,'runname=',runname
  open(unit=ndskpest,file=runname(1:lrunname)//'.dskpest', &
       access='direct',status='old',iostat=ios)
  if(ios==0) write(*,*) 'reading dskpest file ', &
       runname(1:lrunname)//'.dskpest'
  if(ios.ne.0) then
     write(6,*) 'problem opening file ',runname(1:lrunname)//'.dskpest'
     write(6,*) 'will try to read file named dskpest instead'
     inquire (file='dskpest',exist=there,access=dumst(1), &
          direct=dumst(2),form=dumst(3),formatted=dumst(4), &
          sequential=dumst(5) )
     write(*,*) 'exist=',there,' access=',dumst(1), &
          'direct=',dumst(2),'form=',dumst(3), &
          'formatted=',dumst(4),'seq=',dumst(5)
     open(unit=ndskpest,file='dskpest', &
          status='old',form='unformatted',iostat=ios)
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
!  bfieldf   = 1.0
  bfieldf=.39712
  psifact   = bfieldf
  prsfact   = (bfieldf*bfieldf)/amu
  kmap=0

  read (ndskpest) nthd1,npsd1,neqsym
  write(*,*) 'read values nthd1=',nthd1,' npsd1=',npsd1,'neqdsy=',neqsym
  npsi=npsd1-1
  ntht=nthd1+1
  nthd1=ntht+2

  read (ndskpest) psilimp,psiminp,rcnt,zcnt,bfieldf_in
  write(*,*) psilimp,psiminp,rcnt,zcnt,bfieldf_in

  write(*,*) 'ntht=',ntht,'npsi=',npsi,'nthd1=',nthd1,'npsd1=',npsd1
  
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

  read (ndskpest)  (sp     (jj), jj = 1,npsi)
  read (ndskpest)  (pprime    (jj), jj = 1,npsi)
  read (ndskpest)  (qsfin   (jj), jj = 1,npsi)
  read (ndskpest)  (sqvlp  (jj), jj = 1,npsi)
  read (ndskpest)  (fval     (jj), jj = 1,npsi)
  read (ndskpest)  (ffprime   (jj), jj = 1,npsi)
  read (ndskpest)  (sg1    (jj), jj = 1,npsi)
  read (ndskpest)  (sg2    (jj), jj = 1,npsi)
  read (ndskpest)  (psiv (jj), jj = 1,npsi)
!  try using zpsim instead
!     read (ndskpest)  (psiv    (jj), jj = 1,npsi)
  read (ndskpest)  (sg4    (jj), jj = 1,npsi)

  write(*,*) 'eq arrays read'
  write(*,*) 'sp=',sp
  write(*,*) 'psiv=',psiv
  write(*,*) 'ffprime=',ffprime
  write(*,*) 'sg1=',sg1
  write(*,*) 'sg2=',sg2
  write(*,*) 'sg3=',sg3

  nremp=nthd1-ntht
  read(ndskpest) ((sdummy,ii=1,nthd1) ,kk=1,nremp) &
       , ((xs(jj,ii), ii=1,ntht),jj=1,npp)
  read(ndskpest) ((sdummy,ii=1,nthd1),kk=1,nremp) &
       , ((zs(jj,ii), ii=1,ntht),jj=1,npp)

