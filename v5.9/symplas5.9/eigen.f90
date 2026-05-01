!      program testsytrd
!c-----------------------------------------------------------------------
!c program to test lapack routines for eigenvalues and eigenvectors
!c   of symmetric matrix
!c-----------------------------------------------------------------------
!      implicit none
!      integer kn,knsq,kwork
!      parameter(kn=5)
!      real*8 matrix(kn,kn)
!      real*8 eigvec(kn,kn)
!      real*8 workm(kn,kn)
!      real*8 eigval(kn)
!      integer n
!c
!      integer i,j,ih,k
!c-----------------------------------------------------------------------
!c create matrix
!c-----------------------------------------------------------------------
!      n=kn
!      matrix=0.
!      do i=1,n
!         do j=i,n
!            matrix(i,j)=0.31*i+0.1*j
!         end do
!         matrix(i,i)=i*i
!      end do
!      call eigen(matrix,kn,n,eigval,eigvec,workm)
!      write(6,*) 'eigenvalues'
!      write(6,'(5e14.5)') eigval
!      write(6,*) 'eigenvector test z'
!      write(6,'(5e14.5)') eigvec
!      stop
!      end

!=======================================================================
      
subroutine hmult(hnew,h,kn,n)

!-----------------------------------------------------------------------
! multiply two matrices and store result in h
!-----------------------------------------------------------------------
      implicit none
      integer kn,n
!      integer kmax
!      parameter(kmax=100)
      real h(kn,*),hnew(kn,*)
      real hright(n,n)
!      real hright(kmax,kmax)
      integer i,j,k
!      if(n.gt.kmax) then
!         write(6,*)'n > kmax in qmat'
!         stop
!      endif
      do i=1,n
         do j=1,n
            hright(i,j)=h(i,j)
            h(i,j)=0.
         end do
      end do
      do i=1,n
         do j=1,n
            do k=1,n
               h(i,j)=h(i,j)+hnew(i,k)*hright(k,j)
            end do
         end do
      end do
    return
      
end subroutine hmult

!=======================================================================
      
subroutine qmat(h,n,kn,tau,matrix,hnew)

!-----------------------------------------------------------------------
! reconstruct the orhogonal matrix Q: Q**T A Q = T
!  q will be returned in h
!-----------------------------------------------------------------------
      implicit none
!      integer kmax
!      parameter(kmax=100)
      integer n,kn
      real h(kn,kn)
      real hnew(kn,kn)
      real matrix(kn,*)
      real tau(*)
!      real v(kmax)
      real v(n)
      integer i,j,ih
!      if(n.gt.kmax) then
!         write(6,*)'n > kmax in qmat'
!         stop
!      endif
      h=0.
      hnew=0.
      do i=1,n
         h(i,i)=1.
      end do
      h(1,1)=1.-tau(1)
      do ih=2,n-1
         do i=1,ih-1
            v(i)=matrix(i,ih+1)
         end do
         v(ih)=1.
         do i=1,n
            if(i.le.ih) then
               do j=1,ih
                  hnew(i,j)=-tau(ih)*v(i)*v(j)
               end do
               hnew(i,i)=hnew(i,i)+1.
            else
               hnew(i,i)=1.
            endif
         end do
         call hmult(hnew,h,kn,n)
      end do
    return
      
end subroutine qmat

!=======================================================================
      
subroutine eigen(matrix,kn,n,eigval,eigvec,work_matrix)

      use elite_data, only:outfile,verbose
      implicit none
      integer kn,n
      real matrix(kn,kn)
      real work_matrix(kn,kn)
      real eigval(kn)
      real eigvec(kn,kn)
!      integer kmax,knsq,kwork
!      parameter(kmax=50,knsq=kmax*kmax,kwork=10*knsq)
      character*1 uplo
!      integer iwork(kwork)
      integer iwork(10*n*n)
      character*1 compz
      integer liwork
      integer info,lda,lwork
!      real d(kmax),e(kmax),tau(kmax),work(kwork)
      real d(n),e(n),tau(n),work(10*n*n)
!-----------------------------------------------------------------------
! get tridiagonal form
!-----------------------------------------------------------------------
!      if(n.gt.kmax) then
!         write(6,*) 'n=',n,'  kmax=',kmax
!         write(6,*)'n > kmax in eigen'
!         stop
!      endif
      uplo='u'
!      lwork=knsq
      lwork=n*n
      if (verbose .ge. 5) write(outfile,*) 'in eigen, calling dsytrd'
      call dsytrd(uplo,n,matrix,kn,eigval,e,tau,work,lwork,info)
      if (verbose .ge. 5) write(outfile,*) 'return from dsytrd'
      call qmat(eigvec,n,kn,tau,matrix,work_matrix)
!      write(6,*) 'orthogonal matrix q'
!      write(6,'(5e14.5)') eigvec
!-----------------------------------------------------------------------
! calculate eigenvalues and eigenvectors
!-----------------------------------------------------------------------
      compz='v'
      lwork=10*n*n
      liwork=10*n*n
      call dstedc(compz,n,eigval,e,eigvec,kn,work, &
          lwork,iwork,liwork,info)
      if (verbose .ge. 5) write(outfile,*) 'dstedc info=',info
    return
      
end subroutine eigen


!-------------------------------------------------------------------------

subroutine eigen_lapack(matrix,kn,n,eigvalr,eigvali)
! use open source slatec library

  implicit none
  integer kn,n,lwork,info
  real matrix(kn,kn),work(6*n)
  real eigvalr(kn),eigvali(kn)
!  real eigvecr(kn,kn),eigveci(kn,kn)
  real eigvecl(kn,kn),eigvecr(kn,kn)
  character*1 jobvl,jobvr

!  call sgeev(matrix,kn,n,eigval,eigvec,

!  matz=0 ! set to one to get eigenvectors
!  call rg(kn,n,matrix,eigvalr,eigvali,matz,eigvec,iwork,fwork

  lwork=6*n ! must match dimension of work
  jobvl='N'  ! no left eigenvectors
  jobvr='N'  ! no right eigenvectors
  call dgeev(jobvl,jobvr,n,matrix,kn,eigvalr,eigvali,eigvecl,n,eigvecr,n, &
       work,lwork,info)

  if (info .ne. 0) then
     write(*,*) 'problem in eigen_lapack, info=',info
     stop
  endif

end subroutine eigen_lapack










