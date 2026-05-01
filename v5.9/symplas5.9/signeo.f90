
subroutine signeo(rnust,ft,zeff,sigfac)

!     ***************************************
!
! calculates the neoclassical enhancement of the parallel
! conductivity
!      implicit double precision(a-h,o-z)
      implicit none
      real rnust,ft,zeff,sigfac
      real pi,xmax,dx,fac,errf,rle,rlt
      integer nb,i
      real xx,xsq,expo,erfp,gox,taunu,rl1,rl2,fs,ftstar
      real arg,errarg

      pi=4.*atan(1.)
      xmax=10.
      nb=1001
      dx=xmax/(nb-1)
      fac=2./sqrt(pi)
! initialise error funtion
      errf=fac*(exp(-(dx/4.)**2)*0.5*dx)
      rle=3.4*((1.13+zeff)/(2.67+zeff))/zeff
      rlt=2.06*((1.38+zeff)/(3.23+4.68*zeff+zeff*zeff))
      sigfac=0.
      do 10 i=1,nb
	xx=(i-0.5)*dx
	xsq=xx*xx
	expo=exp(-xsq)
!  derivative of error function
	erfp=fac*expo
!  chandrasekhar function
	gox=(errf-xx*erfp)/(2.*xsq)
	taunu=(3./(2.*fac))*(errf-gox+zeff)/(xsq*xx)
	rl1=2.5-xsq
	rl2=35./8.-7.*xsq/2.+xsq*xsq/2.
!  spitzer function
	fs=rle-rlt*rl1+8.*(1./zeff-rle+3.*rlt/2.)*rl2/15.
!  collisionality modified trapped particle fraction
	ftstar=ft/(1.+1.75*rnust*taunu/xx)
	arg=(1.-ftstar)*(1.-ftstar*(taunu*fs-1.))*fs
	arg=arg*xsq*xsq*exp(-xsq)
	sigfac=sigfac+arg*dx
	xx=i*dx
	errarg=fac*exp(-xx*xx)
!  up-date error function
	errf=errf+errarg*dx
 10   continue
      sigfac=sigfac*8./(3.*sqrt(pi))
    return
      
end subroutine signeo
