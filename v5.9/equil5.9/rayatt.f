      subroutine rayatt(x,nw,y,nh,cspln,nxd,dx,dy,xaxd,yaxd,psiaxd,
     $     xmin,xmax,ymin,ymax,psivl,thet,xc,yc,pds,ifail)
      implicit none
      integer nw,nh,nxd
      real x(nw),y(nh),cspln(2,nxd,1),dx,dy,xaxd,yaxd,psiaxd
!      real x(nw),y(nh),cspln(2,nxd,*),dx,dy,xaxd,yaxd,psiaxd
      real xmin,xmax,ymin,ymax
      real psivl,thet,xc,yc,pds(6)
      integer ifail
c...........................................................................
      integer kmax,nmax
      parameter (kmax=129,nmax=20)
      real pi,rndoff0,serrt,serrs
!      parameter (pi=3.14159265358979,rndoff0=1.e-12)
      parameter (pi=3.14159265358979,rndoff0=1.e-9)
!      parameter (serrt=1.e-12,serrs=10.*serrt)
      parameter (serrt=1.e-9,serrs=10.*serrt)
      integer izone,iflag,kountr,newti,ier
      real cost,sint,sgn
      real a,bincp,x1,y1,x2,y2,psi1,psi2
      real xn,yn,delx,dely,serr,xnorm,pnorm,dpsi,dpsids,rndoff
      real rangex,rangey,dum
      data dum/0./
c
      cost=cos(thet)
      sint=sin(thet)
      izone=int(2.*thet/pi+sign(0.5,thet)) 
      iflag=mod(izone,2)
c
      if (iflag .eq. 0)then
         a=sint/cost
         bincp=yaxd-a*xaxd
         sgn=sign(1.,cost)
      else
         a=cost/sint
         bincp=xaxd-a*yaxd
         sgn=sign(1.,sint)
      endif
c
c     do the sliding window search
      x1=xaxd
      y1=yaxd
      psi1=psiaxd
!      write(*,*) 'x1= ',x1,' y1= ',y1,' psi1=',psi1,' sgn=',sgn, 
!     $    ' dx= ',dx,' bincp= ',bincp,' a=',a,' dy=',dy,' kmax= ',kmax,
!     $    ' xmin= ',xmin,' xmax= ',xmax,' ymin=',ymin,' ymax=',ymax,
!     $    ' iflag= ',iflag    
      do kountr=1,kmax
         if (iflag .eq. 0)then
            x2=x1+sgn*dx
            y2=a*x2+bincp
         else
            y2=y1+sgn*dy
            x2=a*y2+bincp
         endif
         rangex=(x2-xmin)*(x2-xmax)
         rangey=(y2-ymin)*(y2-ymax)
         if (rangex .gt. 0. .or. rangey .gt. 0.)then
            ifail =1
            write(*,*) 'xmin= ',xmin,' xmax= ',xmax,
     $           ' ymin=',ymin,'ymax= ',ymax
            call errray1(1,psivl,thet,x2,y2,dum,dum,dum) ! terminates execution
         endif
c         call dbcevl2(x,nw,y,nh,cspln,nxd,x2,y2,pds,ier,1)
c         call bcsplint(x,nw,y,nh,nxd,cspln,x2,y2,pds,1,ier)
         call bcsplint(x,y,cspln,nw,nh,nxd,x2,y2,pds,1,ier)
         psi2=pds(1)
         if ((psivl-psi1)*(psivl-psi2) .le. 0.)go to 10
         x1=x2
         y1=y2
         psi1=psi2
      enddo
      ifail=2
      call errray1(2,psivl,thet,x2,y2,dum,dum,dum) ! terminates execution
 10   continue                  ! done with sliding search
c     
      xnorm=sqrt(xmax*xmax+xmin*xmin+ymax*ymax+ymin*ymin)
      pnorm=abs(psi1)+abs(psi2)
      rndoff=rndoff0*(pnorm/xnorm)
      if (iflag .eq. 0)then
         xn=x1+0.5*sgn*dx
         yn=a*xn+bincp
      else
         yn=y1+0.5*sgn*dy
         xn=a*yn+bincp
      endif
      do newti=1,nmax
c         call dbcevl2(x,nw,y,nh,cspln,nxd,xn,yn,pds,ier,3)
         call bcsplint(x,y,cspln,nw,nh,nxd,xn,yn,pds,3,ier)
         dpsids=pds(2)*cost+pds(3)*sint
         dpsi=pds(1)-psivl
         if (abs(dpsids) .gt. rndoff)then
            serr=-dpsi/dpsids
         else
            serr=-dpsi/rndoff
         endif
         if (abs(serr) .lt. serrt)go to 20
         delx=serr*cost
         dely=serr*sint
         xn=xn+delx
         yn=yn+dely
      enddo
      ifail=3
      call errray1(3,psivl,thet,xn,yn,serr,serrt,serrs)
 20   continue
      xc=xn
      yc=yn
      ifail=0
c      write (3,'(2e15.6,2i5)')psivl,thet,kountr,newti
c
      return
      end


      subroutine errray1(errcode,psi,thet,x,y,err,tol,tola)
      implicit none
      integer errcode
      real psi,thet,x,y,err,tol,tola
c.......................................................................
      if (errcode .eq. 1)then
         write (6,'(a)')'rayatt:  sliding window search out of range'
         write (6,100)'  psi = ',psi,'  thet = ',thet
         write (6,100)'  x   = ',x,  '  y    = ',y
         return
      endif
c
      if (errcode .eq. 2)then
         write (6,'(a)')'rayatt:  sliding window search fails (10)'
         write (6,100)'  psi = ',psi,'  thet = ',thet
         write (6,100)'  x   = ',x,  '  y    = ',y
         return
      endif
c
      if (errcode .eq. 3)then
         write (6,'(a)')'rayatt:  Newton-Raphson interations fails (20)'
         write (6,100)'  psi = ',psi,'  thet = ',thet
         write (6,100)'  err = ',err,'  tol  = ',tol
         if (abs(err) .gt. abs(tola))then
            write (6,'(a)')'rayatt: err > max.tol'
            return
         else
            write (6,'(a)')'rayatt: tol < err < max. tol (warning)'
            return
         endif
      endif
c
      write(6,*) 'error in errray1--errray999'
      stop
 100  format(a,1pe12.4,a,1pe12.4)
      end
