      
subroutine fun2d

  !-----------------------------------------------------------------------
  ! store 2-d eigenfunction in runname.fun2d for updown symmetric eqbria
  !-----------------------------------------------------------------------
  ! 9/01 need to replace x with psi as radial variable, as reversed shear
  !         is now allowed
  ! 2/03 modify to also output imaginary part, f2d_im, needed to make
  !        3d plots

  use elite_data, only: nxinterp,ns,runname,lrunname,cdatetime, &
       codver,nowindow,nmwinhalf,mmax,nm,m0,nn,nx,nmodes, &
       mindex,ffun,meshtype,nd,ndist,aspect,btnorm,del,dmercier, &
       dx,dw,etai,lamdist,ne,ppmult,q0,qref,rmajor,rnorm,tau,te, &
       vloop,qmin,zeff,gam,qa,psigrid,psixx,psimin,gamsq,gamre, &
       gamscl2,jedgen1store,lhspeelstore,rhspeelstore,alpha, &
       shear,pprime_eq,tmatasym,nmlow,dmstore,surfplot,outfile, &
       rotation,verbose
  implicit none
  integer i,m,ibench
  integer ios,nunit
  integer npts
  character*4 shape
  integer npsi              !size of mapped psi mesh
  real percenflux           !outermost flux surface is 
  !psiv(npsi)=psiaxis+(psilim-psiaxis)*percenflux
  real alpsi,del_fix,cc,ss
  complex*16 cexp
  logical delmin,setdel,qafix
  integer nnmin,nnmax
  namelist/equil/npts,shape,&
       percenflux,alpsi,npsi, &
       delmin,nnmin,nnmax,setdel,del_fix,qafix
  integer ng
  namelist/vac/npts,ng !beware npts is in two different namelists
  integer nmwinhlf
  integer npts_equil
  real omegf(nxinterp,ns),rlf(nxinterp,ns),zlf(nxinterp,ns), &
       bplf(nxinterp,ns)
  integer nxplt,ipsigot,ipsig,j,mv,ipsif
  real dpsif,psif,ratf,ome
  real, dimension(:,:),allocatable ::  rlplt,zlplt,bplplt
  complex, dimension(:,:),allocatable :: f2d
  ibench=0
  open(73,file=runname(1:lrunname)//'.omegadat')
  do i=1,nxinterp
     do j=1,ns
        read(73,*)rlf(i,j),zlf(i,j),bplf(i,j),omegf(i,j)
     end do
  end do
  close(73,status='delete')

  if (surfplot .le. 0) then   ! don't generate *.fun2d
     if (verbose .ge. 1) write(6,*) &
          'no surface plot data written, surfplot=',surfplot
     if (verbose .ge. 2) write(outfile,*) &
          'no surface plot data written, surfplot=',surfplot
     return
  endif

  open(unit=12,file=runname(1:lrunname)//'.in',status='old',iostat=ios)
  if(ios.ne.0) then
     write(6,*) 'could not open ',runname(1:lrunname)//'.in in wrteigfunc'
     stop
  endif
  read(12,equil)
  npts_equil=npts
  read(12,vac)
  close(12)
  nunit=30
  open(unit=nunit,file=runname(1:lrunname)//'.fun2d', &
       status='unknown',iostat=ios)
  if(ios.ne.0) then
     write(6,*) 'problem creating/opening ',runname(1:lrunname)//'.fun2d    '
     stop
  endif
  write(nunit,*) cdatetime
  write(nunit,*) "code version: ",codver
  ! write the equil namelist quantities first
  write(nunit,'(a4)') shape
  write(nunit,'(i6)') npts_equil
  ! 12 reals
  if (.not. rotation) gamre=sqrt(gamsq)
!  write(nunit,'(1p,6e13.5)') psimin,sqrt(gamsq*gamscl2),gamsq, &
  write(nunit,'(1p,6e13.5)') psimin,gamre*sqrt(gamscl2),gamre**2, &
       jedgen1store,lhspeelstore,rhspeelstore,alpha(1),shear(1), &
       pprime_eq(1),tmatasym,percenflux,alpsi

  !     namelist & calculated quantities
  if(nowindow) then
     nmwinhlf=-1
  else
     nmwinhlf=nmwinhalf
  endif
  !      nxplt=4*nxinterp+1
  !      nxplt=2*nxinterp+1
  nxplt=surfplot*nxinterp+1
  allocate (f2d(nxplt,ns),rlplt(nxplt,ns),zlplt(nxplt,ns), &
       bplplt(nxplt,ns))
  dpsif=(psigrid(nxinterp)-psigrid(1))/(nxplt-1)
  ipsigot=1
  if (ibench==1) then
     write(76,*)mmax,nm,m0,nn,nx,nmodes,mindex
     do i=1,nx
        write(76,*)psixx(i)
        do j=1,nmodes
           write(76,*)ffun(j,i)
        end do
     end do
     !  close benchmarking data
     close(76)
  endif
  do i=1,nxplt
     psif=psigrid(1)+(i-1)*dpsif
     ipsif=1
     do 5 j=1,nxinterp
        !          if (psigrid(j).gt.xf) goto 5
        if (psigrid(j).lt.psif) goto 5
        ipsif=j
5       continue
        if (i.eq.nxplt) then
           ratf=0.
        else
           ratf=(psif-psigrid(ipsif))/(psigrid(ipsif+1)-psigrid(ipsif))
        end if
        !  find closest x-mesh point
        do 10 j=ipsigot,nx
           !          if (psixx(j).gt.xf) goto 10
           if (psixx(j).lt.psif) goto 10
           ipsig=j
10         continue
           ipsigot=ipsig
           do j=1,ns
              f2d(i,j)=0.
!              f2d_im(i,j)=0.
              if (ipsif.lt.nxinterp) then
                 ome=omegf(ipsif,j)+ratf*(omegf(ipsif+1,j)-omegf(ipsif,j))
                 rlplt(i,j)=rlf(ipsif,j)+ratf*(rlf(ipsif+1,j)-rlf(ipsif,j))
                 zlplt(i,j)=zlf(ipsif,j)+ratf*(zlf(ipsif+1,j)-zlf(ipsif,j))
                 bplplt(i,j)=bplf(ipsif,j)+ratf*(bplf(ipsif+1,j)-bplf(ipsif,j))
              else
                 ome=omegf(ipsif,j)
                 rlplt(i,j)=rlf(ipsif,j)
                 zlplt(i,j)=zlf(ipsif,j)
                 bplplt(i,j)=bplf(ipsif,j)
              end if
              mv=mmax
              do m=1,nmodes/mindex
! ~e^-i mv ome -> cos - i sin
                 cc=cos(mv*ome)
                 ss=-sin(mv*ome)
                 cexp=cmplx(cc,ss)
                 f2d(i,j)=f2d(i,j)+ffun(m,ipsig)*cexp
!                 if (mindex.eq.2) then
!                    f2d(i,j)=f2d(i,j)+ffun(m*mindex-1,ipsig)*cos(mv*ome)  &
!                         +ffun(m*mindex,ipsig)*sin(mv*ome)
!                    f2d_im(i,j)=f2d_im(i,j)-ffun(m*mindex-1,ipsig)*sin(mv*ome)  &
!                         +ffun(m*mindex,ipsig)*cos(mv*ome)
!                 else
!                    f2d(i,j)=f2d(i,j)+ffun(m,ipsig)*cos(mv*ome)
!                    f2d_im(i,j)=f2d_im(i,j)-ffun(m,ipsig)*sin(mv*ome)
!                 end if
                 mv=mv-1
              end do
           end do
        end do
        write(nunit,'(13i6)') meshtype,mmax, &
             nmlow,ndist,nmodes,nn,ns,nxinterp
        write(nunit,'(5i6)') m0,nx,nmwinhlf,npts,ng
        write(nunit,'(1p,6e13.5)') dmstore,btnorm,del, &
             dmercier,dx,dw,etai,lamdist,ne,ppmult, &
             q0,qref,rnorm,tau,te,vloop,qmin,zeff
        write(nunit,'(1p,6e13.5)') gam,qa ! calculated
        write(nunit,*)nxplt,ns
        do  j=1,ns
           write(nunit,'(1p,6e13.5)') (rlplt(i,j),i=1,nxplt)
        end do
        do  j=1,ns
           write(nunit,'(1p,6e13.5)') (zlplt(i,j),i=1,nxplt)
        end do
        do  j=1,ns
           write(nunit,'(1p,6e13.5)') (bplplt(i,j),i=1,nxplt)
        end do
        do  j=1,ns
           write(nunit,'(1p,6e13.5)') (real(f2d(i,j)),i=1,nxplt)
        end do
        do  j=1,ns
           write(nunit,'(1p,6e13.5)') (aimag(f2d(i,j)),i=1,nxplt)
        end do
        close(nunit)

        return

      end subroutine fun2d
!










