!
!   P.B. Snyder 1/27/00
!     collect functions used by various elitevac routines


!----------------------------------------------------
   function cosmw(mval,theta_v)
! real part of high m dependence which is factored out of b
      use elitevac_data, only: pi,tsurf,omegasurf,y2omegasurf,nsurfpts
      implicit none
      integer mval
      real cosmw,theta_v,omega_v
      integer i
!      include 'i_vac.f'
      real t_v
      t_v=theta_v
      if(t_v.lt.0.) t_v=t_v+2.*pi
      if(t_v.gt.2.*pi) t_v=t_v-2.*pi
!      write(6,*) 'cosmw:t_v=',t_v
!      write(6,*) (tsurf(i),i=1,nsurfpts)
      call splint(tsurf,omegasurf,y2omegasurf,nsurfpts,t_v,omega_v)
      cosmw=cos(mval*omega_v)
      return
   end function cosmw
!-----------------------------------------------------------------------
   function sinmw(mval,theta_v)
! imaginary part of high m dependence which is factored out of b
      use elitevac_data, only: pi,tsurf,omegasurf,y2omegasurf,nsurfpts
      implicit none
      integer mval
      real sinmw,theta_v,omega_v
!      include 'i_vac.f'
      real t_v
      t_v=theta_v
      if(t_v.lt.0.) t_v=t_v+2.*pi
      if(t_v.gt.2.*pi) t_v=t_v-2.*pi
!      write(6,*) 't_v=',t_v
      call splint(tsurf,omegasurf,y2omegasurf,nsurfpts,t_v,omega_v)
      sinmw=sin(mval*omega_v)
      return
   end function sinmw
!-----------------------------------------------------------------------
   function bpv(theta_v)
! poloidal field on the boundary surface
      use elitevac_data, only: pi,tsurf,bpsurf,y2bpsurf,nsurfpts
      implicit none
      real bpv,theta_v
!      include 'i_vac.f'
      real t_v
      t_v=theta_v
      if(t_v.lt.0.) t_v=t_v+2.*pi
      if(t_v.gt.2.*pi) t_v=t_v-2.*pi
!      write(6,*) 't_v=',t_v
      call splint(tsurf,bpsurf,y2bpsurf,nsurfpts,t_v,bpv)
      return
   end function bpv
!-----------------------------------------------------------------------
   function rv(theta_v)
! coordinate R on the boundary surface
      use elitevac_data, only: pi,tsurf,rsurf,y2rsurf,nsurfpts
      implicit none
      real rv,theta_v
!      include 'i_vac.f'
      real t_v
      t_v=theta_v
      if(t_v.lt.0.) t_v=t_v+2.*pi
      if(t_v.gt.2.*pi) t_v=t_v-2.*pi
!      write(6,*) 't_v=',t_v
      call splint(tsurf,rsurf,y2rsurf,nsurfpts,t_v,rv)
      return
   end function rv
!-----------------------------------------------------------------------
   function rt(theta_v)
! dR/dtheta on the boundary surface
      use elitevac_data, only: pi,tsurf,drdtsurf,y2drdtsurf,nsurfpts
      implicit none
      real rt,theta_v
!      include 'i_vac.f'
      real t_v
      t_v=theta_v
      if(t_v.lt.0.) t_v=t_v+2.*pi
      if(t_v.gt.2.*pi) t_v=t_v-2.*pi
!      write(6,*) 't_v=',t_v
      call splint(tsurf,drdtsurf,y2drdtsurf,nsurfpts,t_v,rt)
      return
   end function rt
!-----------------------------------------------------------------------
   function zv(theta_v)
! this is coordinate z on the boundary surface
      use elitevac_data, only: pi,tsurf,zsurf,y2zsurf,nsurfpts
      implicit none
      real zv,theta_v
!      include 'i_vac.f'
      real t_v
      t_v=theta_v
      if(t_v.lt.0.) t_v=t_v+2.*pi
      if(t_v.gt.2.*pi) t_v=t_v-2.*pi
!      write(6,*) 't_v=',t_v
      call splint(tsurf,zsurf,y2zsurf,nsurfpts,t_v,zv)
      return
   end function zv
!-----------------------------------------------------------------------
   function zt(theta_v)
! dz/dtheta on the boundary surface
      use elitevac_data, only: pi,tsurf,dzdtsurf,y2dzdtsurf,nsurfpts
      implicit none
      real zt,theta_v
!      include 'i_vac.f'
      real t_v
      t_v=theta_v
      if(t_v.lt.0.) t_v=t_v+2.*pi
      if(t_v.gt.2.*pi) t_v=t_v-2.*pi
!      write(6,*) 't_v=',t_v
      call splint(tsurf,dzdtsurf,y2dzdtsurf,nsurfpts,t_v,zt)
      return
   end function zt
!-----------------------------------------------------------------------
   function frgn(t1,t2,r1,z1,r2,z2,lam,s)
! this is the green's function R*G_n(theta_1,theta_2)
!      use elitevac_data, only: gamma_n05,ntor,n05,n1,pi,gamma_n1
      use elitevac_data, only: gn05ovgn1,ntor,n05,n1,pi
      implicit none
!      include 'i_vac.f'
      real frgn,t1,t2
      real lam,r1,r2,z1,z2,s,p
!      real rv,zv
      real Hyper2F1
      real hyp2f1f_s
      real onehalf
      onehalf=0.5d0
!      r1=rv(t1)
!      r2=rv(t2)
!      z1=zv(t1)
!      z2=zv(t2)
!      lam=1.+((r1-r2)**2+(z1-z2)**2)/(2.*r1*r2)
!      s=lam/(sqrt(lam*lam-1.0))
      p=(s-1.)/(s+1.)
! use either Hyper2F1 or faster series version (Pletzer) hyp2f1f_s
! at the moment, I am using the slower full version because p>0.999
! gets called when I use ng=8 or 16. Could test for p before calling
! but total execution time is not that bad anyway.
!      frgn = sqrt(r2)*gamma_n05*p**(ntor*0.5d0+0.25d0)* &
!          Hyper2F1(n05,onehalf,n1,p) / &
!          ( sqrt( pi * r1 ) * gamma_n1 )   
      frgn = sqrt(r2)*gn05ovgn1*p**(ntor*0.5d0+0.25d0)* &
          Hyper2F1(n05,onehalf,n1,p) / &
          ( sqrt( pi * r1 ) )   


      return
   end function frgn
!-----------------------------------------------------------------------
   function frgsing(t1,t2)
! this is the analytic singularity which is subtracted from frgn before
! numerical integration
      use elitevac_data, only: pi
      implicit none
!      include 'i_vac.f'
      real frgsing,t1,t2
      frgsing =  log( abs(t2-t1      ) ) + &
          log( abs(t2-t1-2.*pi) ) + &
          log( abs(t2-t1+2.*pi) )
      frgsing=-frgsing/pi
      return
   end function frgsing
!-----------------------------------------------------------------------
   function frgn1(t1,t2,r1,z1,r2,z2,lam,s)
! R*G_(n+1)
!      use elitevac_data, only: pi,gamma_n15,n1,n15,n2,gamma_n2
      use elitevac_data, only: pi,gn15ovgn2,n1,n15,n2
      implicit none
!      include 'i_vac.f'
      real frgn1,t1,t2
      real r1,r2,z1,z2,lam,s,p
!      real rv,zv
      real Hyper2F1
      real hyp2f1f_s
      real onehalf
      onehalf=0.5d0
!      r1=rv(t1)
!      r2=rv(t2)
!      z1=zv(t1)
!      z2=zv(t2)
!      lam=1.+((r1-r2)**2+(z1-z2)**2)/(2.*r1*r2)
!      s=lam/(sqrt(lam*lam-1.0))
      p=(s-1.)/(s+1.)
! use either Hyper2F1 or faster series version (Pletzer) hyp2f1f_s
!      frgn1 =sqrt(r2) *gamma_n15 *p**(n1*0.5+0.25) * &
!          Hyper2F1(n15,onehalf,n2,p) / &
!          ( sqrt( pi*r1 ) * gamma_n2 )
      frgn1 =sqrt(r2) *gn15ovgn2 *p**(n1*0.5+0.25) * &
          Hyper2F1(n15,onehalf,n2,p) / &
          ( sqrt( pi*r1 ) )
      return
   end function frgn1
!-----------------------------------------------------------------------
   function Hyper2F1(a,b,c,x)
! the full hypergeometric coding which calls the c routine hypf1f
      implicit none
      real Hyper2F1,a,b,c,x
      real hyp2f1f
      if ( x.gt.1.0 ) then      !F1 is complex for arg > 1
         write(6,*) "Hyper2F1 error"
         stop
      else if ( x.le. -1.0 ) then
         Hyper2F1 = (1-x)**(-a) *hyp2f1f(a,c-b,c, x/(x-1)) 
      else                                !-1 < x < 1
         Hyper2F1 = hyp2f1f(a,b,c, x )
      endif
      return
   end function Hyper2F1

!-----------------------------------------------------------------------
   function gammln(xx)

!  log of gamma function, from numerical recipes
      implicit none
      integer j
      real cof(6),ser,stp,tmp,xx,x,y,gammln
      data cof /76.18009172947146d0,-86.50532032941677d0, &
                  24.01409824083091d0,-1.231739572450155d0, &
                  .1208650973866179d-2,-.5395239384953d-5/

      x=xx; y=x
      tmp=x+5.5
      tmp=(x+.5)*log(tmp)-tmp
      ser=1.000000000190015d0
      stp=2.5066282746310005d0
      do j=1,6
         y=y+1
         ser=ser+cof(j)/y;
      enddo
      gammln=tmp+log(stp*ser/x);
      return
    end function gammln

!----------------------------------------------------------
    function hyp2f1f(a,b,c,x)
! 4/26/00 replace old hyp2f1f.c with a
! wrapper around numerical recipes hypgeo function
      real a,b,c,x,hyp2f1f
      complex aa,bb,cc,xx,hypgeo

      aa=cmplx(a,0.)
      bb=cmplx(b,0.)
      cc=cmplx(c,0.)
      xx=cmplx(x,0.)
!      write(*,*) 'aa=',aa,bb,cc,xx
! evaluate with numerical recipes
!      hyptemp=hypgeo(aa,bb,cc,xx)
!      write(*,*) 'hyptemp=',hyptemp
      hyp2f1f=real(hypgeo(aa,bb,cc,xx))
!      write(*,*) 'hyp2f1f=',hyp2f1f
      return
    end function hyp2f1f











