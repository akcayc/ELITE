
subroutine getxx(npts,dxmin,xmin,xmax,xx,ixmax)

!-----------------------------------------------------------------------
! this is a completely rewritten version of getxx
! input is:
!   npts: number of points per unit range of x including 1 endpoint
!         i.e., if x(i)=1 then x(i+npts)=2
!   dxmin: the minimum dx value which applies at rational surfaces,i.e, 
!          integer values of x
!   xmin: the plasma boundary (elsewhere referred to as Delta)
!   xmax: the maximum value of x, assumed to be an integer
! output is:
!   xx(ix): an array of x values ranging from xx(1)=xmin to 
!            xx(ixmax)=xmax
!   ixmax: index for maximum x in xx array
!-----------------------------------------------------------------------
      use elite_data, only: outfile,verbose
      implicit none
      integer ixmax
      real dxmin,xmin,xmax,xx(*)
      integer nint,npts,i,ix,ic,m,im
      integer nhalf,npoints
      real df,xnew,eps,bbar,dxbig
      real xdiff,xdiffnew,xmatch
      integer imatch
!      integer kx
!      parameter(kx=505)
!      real x(kx)
      real x(2*(npts+1)/2)
!-----------------------------------------------------------------------
!  1.0 check input values
!-----------------------------------------------------------------------
      if(xmin.lt.0.) then
         write(6,*) "sorry but getxx is too dumb to deal with ", &
             "your del"
         stop
      endif
      npoints=2*((npts+1)/2)
!      if(npoints+1.gt.kx) then
!         write(6,*) 'npts in getxx is too large for kx parameter'
!         write(6,*) 'npts=',npts, '; kx=',kx
!         stop
!      endif
!-----------------------------------------------------------------------
!  2.0 create array of points from 0-1 to be used on each interval between
!      rational surfaces
!-----------------------------------------------------------------------
      x(1)=0.
      x(npoints+1)=1.
      nhalf=npoints/2+1
      x(nhalf)=0.5d0
      if (verbose .ge. 5) write(outfile,*) 'calling findeps'
      call findeps(dxmin,nhalf-1,0.5d0,eps)
      if (verbose .ge. 5) write(outfile,*) 'back from findeps'
      do i=2,nhalf-1
         x(i)=(x(i-1)+dxmin)/(1.-eps)
         x(npoints+2-i)=1.-x(i)
      end do
!      write(6,'(1p,5e13.5)') (x(i),i=1,npoints+1)
!-----------------------------------------------------------------------
!  3.0 require special array for xx between x=xmin and x=1
!-----------------------------------------------------------------------
      xx(1)=xmin
      xmatch=(xmin+1.)*0.5
      imatch=1
      xdiff=abs(x(imatch)-xmatch)
      do i=1,npoints+1
         xdiffnew=abs(x(i)-xmatch)
         if(xdiffnew.lt.xdiff) then
            imatch=i
            xdiff=xdiffnew
         endif
      end do
      dxbig=(x(imatch+1)-x(imatch))
      eps=(dxbig-dxmin)/(x(imatch)-xmin)
      bbar=dxmin-eps*(xmin+dxmin)
      ix=1
      do while (xx(ix).lt.x(imatch))
         ix=ix+1
         xx(ix)=(xx(ix-1)+bbar)/(1.-eps)
      end do
      ix=ix-1
      if(x(imatch)-xx(ix).lt.0.05*dxbig) ix=ix-1
      do i=imatch,npoints+1
         ix=ix+1
         xx(ix)=x(i)
      end do
!      write(6,*) 'imatch=',imatch
!      write(6,*) 'xmatch=',xmatch
!      write(6,*) 'dxbig=',dxbig
!      write(6,*) 'bbar=',bbar
!      write(6,*) 'eps=',eps
!      write(6,*) 'ix=',ix
!      write(6,'(1p,5e13.5)') (xx(i),i=1,ix)
!-----------------------------------------------------------------------
!  4.0 can now construct the rest of the mesh using the x array
!-----------------------------------------------------------------------
      m=xmax
      do im=1,m-1
         do ic=2,npoints+1
            ix=ix+1
            xx(ix)=im+x(ic)
         end do
      end do

! 8/04 add extra grid points to go in to <= psimin
      ic=1
      do while (xx(ix) < xmax)
         ix=ix+1
         ic=ic+1
         xx(ix)=m+x(ic)
!         write(*,*) 'ix=',ix,' xx(ix)=',xx(ix),' xmax=',xmax
      enddo
      ixmax=ix
! put final point exactly at xmax, subtract small amount for roundoff
!   problem which causes psi interpolation issue
      xx(ixmax)=xmax-1.e-8
!      write(6,*) "ixmax=",ixmax
!      write(31,*) (xx(i),i=1,ixmax)

    return

end subroutine getxx


!-----------------------------------------------------------------------
subroutine findeps(b,nwant,xn,eps)
!-----------------------------------------------------------------------
! this routine calculates the required eps for
!    x[i+1]-x[i]= eps*x[i+1]+b
!    when there are nwant+1 points in the closed interval from 0 to xn
!-----------------------------------------------------------------------
      use elite_data, only: outfile,verbose
      implicit none
      real b,xn,eps,de,tol
      integer nwant,i,count
      real nguess
!      write(6,*) 'b=',b
!      write(6,*) 'nwant=',nwant
!      write(6,*) 'xn=',xn
! make initial guess for eps, which will be too large
      eps=(-2.d0*(b*nwant - xn))/(b*(-nwant + nwant**2))
!      write(6,*) 'eps=',eps
      if (eps.ge.1.) then
         if (verbose .ge. 2) write(outfile,*) 'resetting eps=0.95'
         eps=0.95
      endif
! calculate n for this epsguess
      nguess=1. - alog(1. - eps + (eps*xn)/b)/alog(1. - eps)
!      write(6,*) 'nguess=',nguess
      count=0
      do while (nguess.lt.nwant)
!         write(*,*) 'count=',count,' eps=',eps,' nguess=',nguess
         eps=0.5*eps
         nguess=1. - alog(1. - eps + (eps*xn)/b)/alog(1. - eps)
!         write(6,*) 'eps=',eps
!         write(6,*) 'nguess=',nguess
         count=count+1
         if(count.gt.30) then
            write(6,*) 'runaway job in findeps-1'
            stop
         endif
      end do
! now use Newton's Method to get an accurate value of eps
!      tol=1.d-8
      tol=1.d-13
      de=1.0
      count=0
      do while (abs(de).gt.tol)
         de= (alog(1. - eps)**2* &
              (-1. + nwant + alog(1. + eps*(-1. + xn/b))/ &
              alog(1. - eps)))/ &
              (-(((b - xn)*alog(1. - eps))/ &
              (b*(-1. + eps) - eps*xn)) + & 
              alog(1. + eps*(-1. + xn/b))/(-1. + eps))
         eps=eps+de
         count=count+1
         if(count.gt.30) then
            write(6,*) 'runaway job in findeps-2'
            stop
         endif
      end do
      if (verbose .ge. 3) write(outfile,*) 'eps=',eps
    return

end subroutine findeps


subroutine getxx4

  use elite_data, only: xx,psixx,xinterp,psigrid,ndist,nx,nxinterp,del, &
       verbose
  implicit none
  real dxguess,dpsimin,dpsimax,psii,xi,dxdpsii,check,dum,dum2(nxinterp)
  real bnd_set,dpsii,olddxguess,xinew
  integer count,i

! attempt to evenly space in x even for non-monotonic q profiles
!  should be similar to meshtype=2, use ndist total points
!  start with psii at outermost point, then grid inward until
!   psigrid(nxinterp) is reached, adjust dxguess until the total
!   number of points needed = ndist
        nx=ndist
!        dxguess=( abs(xinterp(1)-minval(xinterp)) + abs(xinterp(nxinterp) &
!             -minval(xinterp)))/ndist

        dxguess=0.
        do i=2,nxinterp
           dxguess=dxguess+abs( xinterp(i)-xinterp(i-1) )
        enddo
        dxguess=dxguess/ndist  ! dxguess is positive

        bnd_set=-2.d30
        ! temporarily set a minimum on dpsi equal to .001 of avg
        dpsimin=.001*(psigrid(1)-psigrid(nxinterp))/ndist
        dpsimax=5.*(psigrid(1)-psigrid(nxinterp))/ndist
        if (verbose .ge. 2) write(*,*) 'dpsimin=',dpsimin,'dpsimax=',dpsimax
        ! spline apparantly requires increasing 'x' grid, use -psigrid
        call spline(-psigrid,xinterp,nxinterp,bnd_set,bnd_set,dum2)

50      count=2
        psii=psigrid(1)
        if (verbose .ge. 2) write(*,*) 'min,dxguess=',minval(xinterp),dxguess,dxguess*ndist

        ! spline to get x and dpsi/dx
        ! first call spline to get interpolating fct, contained in dum2

60      continue

        ! call splint to get x value at psii
        call splint(-psigrid,xinterp,dum2,nxinterp,-psii,xi)
        ! call zsplint to get dxdpsi value at psii
        call zsplint(-psigrid,xinterp,dum2,dum2,nxinterp,-psii,check, &
             dxdpsii,dum)  ! note dxdpsi has wrong sign here

!        write(*,*) 'psii=',psii,'xi=',xi,'dxdpsii=',dxdpsii

        dpsii=-abs(dxguess/dxdpsii)
        if (abs(dpsii) < abs(dpsimin) ) then
           if (verbose .ge. 2) write(*,*) &
                'dpsii<dpsimin, count=',count,dpsii,dpsimin,psii
           dpsii=-abs(dpsimin)
        endif
        if (abs(dpsii) > abs(dpsimax) ) then
!           write(*,*) 'dpsii>dpsimax, count=',count,dpsii,dpsimin,psii
           dpsii=-abs(dpsimax)
        endif

! 1/17/02 don't let dx exceed dxguess by more than factor of two
!    this can happen when shear is changing quickly near zero
        call splint(-psigrid,xinterp,dum2,nxinterp,-(psii+dpsii),xinew)
        if ( abs(xinew-xi) > 2.*dxguess ) then
           do while ( abs(xinew-xi) > 2.*dxguess )
              if (verbose .ge. 2) write(*,*) &
                   'dx > 2*dxguess, reducing psii',xinew,xi,dxguess
              dpsii=0.5*dpsii
              call splint(-psigrid,xinterp,dum2,nxinterp,-(psii+dpsii),xinew)
           enddo
        endif

        psii=psii+dpsii
        if (psii < psigrid(nxinterp)) then
           psii=psigrid(nxinterp)
           if ((count).eq.ndist) goto 70
           olddxguess=dxguess
           dxguess=dxguess*(real(count)/real(ndist))
           if (verbose .ge. 3) write(*,*) 'count=',count,'dxguess=',dxguess
           goto 50
        else
           count=count+1
           goto 60
        endif
           

70      continue
        if (verbose .ge. 3) write(*,*) 'count=ndist',count,ndist,dxguess
        ! now that dxguess value gives correct number of points (ndist)
        ! lay out actual psixx and xx grids
        psii=psigrid(1)
        psixx(1)=psigrid(1)
        xx(1)=del
        do i=2,ndist
!           ! call splint to get x value at psii
           call splint(-psigrid,xinterp,dum2,nxinterp,-psii,xi)
           ! call zsplint to get dxdpsi value at psii
           call zsplint(-psigrid,xinterp,dum2,dum2,nxinterp,-psii,check, &
                dxdpsii,dum)  ! note dxdpsi has wrong sign here

!           write(*,*) 'psii=',psii,'xi=',xi,'dxdpsii=',dxdpsii

           dpsii=-abs(dxguess/dxdpsii)
           if (abs(dpsii) < abs(dpsimin) ) then
              if (verbose .ge. 2) write(*,*) &
                   'dpsii<dpsimin, count=',count,dpsii,dpsimin,psii
              dpsii=-abs(dpsimin)
           endif
           if (abs(dpsii) > abs(dpsimax) ) then
              if (verbose .ge. 2) write(*,*) &
                   'dpsii>dpsimax, count=',count,dpsii,dpsimin,psii
              dpsii=-abs(dpsimax)
           endif

! 1/17/02 don't let dx exceed dxguess by more than factor of two
!    this can happen when shear is changing quickly near zero
           call splint(-psigrid,xinterp,dum2,nxinterp,-(psii+dpsii),xinew)
           if ( abs(xinew-xi) > 2.*dxguess ) then
              do while ( abs(xinew-xi) > 2.*dxguess )
                 if (verbose .ge. 2) write(*,*) &
                      'dx > 2*dxguess, reducing psii',xinew,xi,dxguess
                 dpsii=0.5*dpsii
                 call splint(-psigrid,xinterp,dum2,nxinterp, &
                      -(psii+dpsii),xinew)
              enddo
           endif

           psii=psii+dpsii
           call splint(-psigrid,xinterp,dum2,nxinterp,-psii,xi)
           psixx(i)=psii
           xx(i)=xi
        enddo
        if (psixx(ndist)>psigrid(nxinterp) ) then
           write(*,*) 'problem creating grid, psixx>psigrid(nxinterp)', &
                'psixx=',psixx(ndist),'psigrid=',psigrid(nxinterp)
           stop
        endif
        psixx(ndist)=psigrid(nxinterp)+abs(.1*dpsimin)
        if (psixx(ndist)>psixx(ndist-1)) then
           write(*,*) 'problem with grid, psixx(ndist)>psixx(ndist-1)'
           stop
        endif
        call splint(-psigrid,xinterp,dum2,nxinterp,-psixx(ndist),xx(ndist))

! debug!!!
        if (verbose .ge. 4) write(*,*) 'dxavg=',(xx(ndist)-xx(1))/(ndist-1)
        do i=2,ndist
           if (verbose .ge. 5) write(*,*) 'i=',i,'psixx=',psixx(i),'xx=', &
                xx(i),'dx=',xx(i)-xx(i-1), &
                'err=', ((xx(i)-xx(i-1))-(xx(ndist)-xx(1))/(ndist-1))/ &
                ((xx(ndist)-xx(1))/(ndist-1))
        enddo


end subroutine getxx4

