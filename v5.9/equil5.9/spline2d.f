      subroutine splnyr(x,y,n,yp,wk)
      implicit none
      integer n
      real x(n),y(n),yp(n),wk(n,5)
c...........................................................................
      real t,v,btmp
      real a2,b2,c2,an2,bn2,cn2
      integer j
c
      t=(x(2)-x(1))/(x(3)-x(2))
      v=(x(n)-x(n-2))/(x(n-1)-x(n-2))
c
      do j=1,n-1
         yp(j)=1./(x(j+1)-x(j))
      enddo
c
c     set up the maxtrix elements for y'
      wk(2,2)=(2.+t+1./t)*yp(2)
      wk(2,3)=(1.+t)*yp(2)
      wk(2,4)=(y(2)-y(1))*yp(1)**2+(y(3)-y(2))*(2.+3.*t)*yp(2)**2
      wk(n-1,1)=(v+1./v-2.)*yp(n-2)
      wk(n-1,2)=(v-1.)*yp(n-2)
      wk(n-1,4)=(1./v**2+2.*v-3.)*yp(n-2)**2*(y(n-1)-y(n-2))+
     $     (y(n)-y(n-1))*yp(n-2)**2/v**2
      if (n .gt. 4)then
         do j=3,n-2
            wk(j,1)=yp(j-1)
            wk(j,2)=2.*(yp(j-1)+yp(j))
            wk(j,3)=yp(j)
            wk(j,4)=3.*((y(j)-y(j-1))*yp(j-1)**2+(y(j+1)-y(j))*yp(j)**2)
         enddo
      endif
c
      btmp=wk(2,2)
      yp(2)=wk(2,4)/btmp
      do j=3,n-1
         wk(j,5)=wk(j-1,3)/btmp
         btmp=wk(j,2)-wk(j,1)*wk(j,5)
         yp(j)=(wk(j,4)-wk(j,1)*yp(j-1))/btmp
      enddo
      do j=n-2,2,-1
         yp(j)=yp(j)-yp(j+1)*wk(j+1,5)
      enddo  
c
      a2=yp(2)*(x(3)-x(2))
      b2=3.*(y(3)-y(2))-(x(3)-x(2))*(2.*yp(2)+yp(3))
      c2=-2.*(y(3)-y(2))+(x(3)-x(2))*(yp(2)+yp(3))
      yp(1)=(a2+t*(-2.*b2+t*3.*c2))/(x(3)-x(2))
c
      an2=yp(n-2)*(x(n-1)-x(n-2))
      bn2=3.*(y(n-1)-y(n-2))-(x(n-1)-x(n-2))*(2.*yp(n-2)+yp(n-1))
      cn2=-2.*(y(n-1)-y(n-2))+(x(n-1)-x(n-2))*(yp(n-2)+yp(n-1))
      yp(n)=(an2+v*(2.*bn2+v*3.*cn2))/(x(n-1)-x(n-2))
c
      return
      end
      subroutine bcspline(x,y,f,nx,ny,nxp,c)
      implicit none
      integer nx,ny,nxp
      real f(nxp,ny),x(nx),y(ny),c(4,nxp,ny)
c..........................................................................
      integer kx
      parameter (kx=2051)
      real g(kx),gp(kx),wk(5*kx)
      integer i,j
c
      if (nx .lt. 4 .or. nx .gt. kx)then
         write (6,'(a)')'bcspline:   nx < 4 .or. nx > kx = 2051'
         stop
      endif
      if (ny .lt. 4 .or. ny .gt. kx)then
         write (6,'(a)')'bcspline:   ny < 4 .or. ny > kx = 2051'
         stop
      endif
      if (nx .gt. nxp)then
         write (6,'(a)')'bcspline:   nx > nxp (storage error)'
         stop
      endif
c
      if (x(nx)-x(1) .lt. 0.)then
         write (6,'(a)')'bcspline:  x not monotonically increasing'
         stop
      endif
      if (y(ny)-y(1) .lt. 0.)then
         write (6,'(a)')'bcspline:  y not monotonically increasing'
         stop
      endif
c
      do j=1,ny
         do i=1,nx
            c(1,i,j)=f(i,j)
         enddo
      enddo
c
      do j=1,ny
         do i=1,nx
            g(i)=f(i,j)
         enddo
         call splnyr(x,g,nx,gp,wk)
         do i=1,nx
            c(2,i,j)=gp(i)
         enddo
      enddo
c
      do i=1,nx
         do j=1,ny
            g(j)=f(i,j)
         enddo
         call splnyr(y,g,ny,gp,wk)
         do j=1,ny
            c(3,i,j)=gp(j)
         enddo
      enddo
c
      do i=1,nx
         do j=1,ny
            g(j)=c(2,i,j)
         enddo
         call splnyr(y,g,ny,gp,wk)
         do j=1,ny
            c(4,i,j)=gp(j)
         enddo
      enddo
c
      return
      end
      subroutine bcsplint(x,y,c,nx,ny,nxp,xl,yl,pds,icalc,ier)
      implicit none
      integer nx,ny,nxp,icalc,ier
      real x(nx),y(ny),c(4,nxp,ny),xl,yl,pds(6)
c...........................................................................
      integer i,j
      real hx, hy, sux(2), suy(2), su(2), svx(2), sv(2), sxy(2),
     $     u, v, spln0, spln1, spln2, s0, sh, sp0, sph, h, d
c
      spln0(s0,sh,sp0,sph,h,d) = s0+d*(h*sp0+d*(3.*(sh-s0)-
     $  (sph+2.*sp0)*h+d*(2.*(s0-sh)+(sph+sp0)*h)))
      spln1(s0,sh,sp0,sph,h,d) = sp0+d*(6.*(sh-s0)/h-2.0*
     $  (sph+2.*sp0)+3.*d*(2.*(s0-sh)/h+(sph+sp0)))
      spln2(s0,sh,sp0,sph,h,d) = 6.*(sh-s0)/h**2-2.0*
     $  (sph+2.*sp0)/h+d*(2.*(s0-sh)/h**2+(sph+sp0)/h)*6.0
c
      save i,j
      data i,j/0,0/
c
      ier = 0
c --- correlated table search for xl
      call tableintrp (x, nx, xl, i)
      if ( i .eq. 0 )then
         if (abs(xl-x(1)) .le. abs(x(2)-x(1)))then
            i=1
            ier = 3
            write (6,'(a,i2)')'bcsplint: x < x(1), ier = ',ier
            go to 10
         else
            ier=33
            write (6,'(a,i2)')'bscplint: x < x(1), ier = ',ier
            return
         endif
      endif
      if ( i .eq. nx)then
         if (abs(xl-x(nx)) .le. abs(x(nx-1)-x(nx)))then
            i=nx-1
            ier=5
            write (6,'(a,i2)')'bcsplint: x > x(nx), ier = ',ier
         else
            ier=35
            write (6,'(a,i2)')'bcsplint: x > x(nx), ier = ',ier
            return
         endif
      endif
 10   continue
c --- correlated table search for yl
      call tableintrp (y, ny, yl, j)
      if ( j .eq. 0 )then
         if (abs(yl-y(1)) .le. abs(y(2)-y(1)))then
            j=1
            ier=4
            write (6,'(a,i2)')'bcsplint: y < y(1), ier = ',ier
            go to 20
         else
            ier=34
            write (6,'(a,i2)')'bcsplint: y < y(1), ier = ',ier
            return
         endif
      endif
      if ( j .eq. ny)then
         if (abs(yl-y(ny)) .le. abs(y(ny-1)-y(ny)))then
            j=ny-1
            ier=6
            write (6,'(a,i2)')'bcsplint: y > y(ny), ier = ',ier
         else
            ier=36
            write (6,'(a,i2)')'bcsplint: y > y(ny), ier = ',ier
            return
         endif
      endif
 20   continue
c
      hx   = x(i+1)-x(i)
      hy   = y(j+1)-y(j)
      u    = (xl-x(i))/hx
      v    = (yl-y(j))/hy
c
      sv(1)= spln0(c(1,i,j),c(1,i,j+1),c(3,i,j),c(3,i,j+1),hy,v)
      svx(1)=spln0(c(2,i,j),c(2,i,j+1),c(4,i,j),c(4,i,j+1),hy,v)
      sv(2)= spln0(c(1,i+1,j),c(1,i+1,j+1),c(3,i+1,j),c(3,i+1,j+1),hy,v)
      svx(2)=spln0(c(2,i+1,j),c(2,i+1,j+1),c(4,i+1,j),c(4,i+1,j+1),hy,v)      
      pds(1)=spln0(sv(1),sv(2),svx(1),svx(2),hx,u)
      if (icalc .le. 1)return
c
      su(1)= spln0(c(1,i,j),c(1,i+1,j),c(2,i,j),c(2,i+1,j),hx,u)
      suy(1)=spln0(c(3,i,j),c(3,i+1,j),c(4,i,j),c(4,i+1,j),hx,u)
      su(2)= spln0(c(1,i,j+1),c(1,i+1,j+1),c(2,i,j+1),c(2,i+1,j+1),hx,u)
      suy(2)=spln0(c(3,i,j+1),c(3,i+1,j+1),c(4,i,j+1),c(4,i+1,j+1),hx,u)
      pds(2)=spln1(sv(1),sv(2),svx(1),svx(2),hx,u)
      if (icalc .eq. 2)return
      pds(3)=spln1(su(1),su(2),suy(1),suy(2),hy,v)
      if (icalc .eq. 3)return
c
      sux(1)=spln1(c(1,i,j),c(1,i+1,j),c(2,i,j),c(2,i+1,j),hx,u)
      sxy(1)=spln1(c(3,i,j),c(3,i+1,j),c(4,i,j),c(4,i+1,j),hx,u)
      sux(2)=spln1(c(1,i,j+1),c(1,i+1,j+1),c(2,i,j+1),c(2,i+1,j+1),hx,u)
      sxy(2)=spln1(c(3,i,j+1),c(3,i+1,j+1),c(4,i,j+1),c(4,i+1,j+1),hx,u)
      pds(4) = spln1(sux(1),sux(2),sxy(1),sxy(2),hy,v)
      if (icalc .eq. 4)  return
      pds(5) = spln2(sv(1),sv(2),svx(1),svx(2),hx,u)
      if (icalc .eq. 5)  return
      pds(6) = spln2(su(1),su(2),suy(1),suy(2),hy,v)
c
      return
      end


      subroutine tableintrp (xx, n, x, jlo)
c
c ----------------------------------------------------------------------
c --- correlated table search routine. use jlo from previous call to
c --- get jlow for current value of x. if jlo from previous call is
c --- no good(jlo=0 or jlo=n) then do a binary search.
c --- this routine returns jlo so that
c   xx(jlo) .le. x  .le. xx(jlo+1) if jlo=1
c   xx(jlo) .lt. x  .le. xx(jlo+1) if jlo =2,3...n-1
c   it is assumed that xx(j),j=1,2..n is monotonically
c   increasing,this is NOT checked for.
c  this is a modified version of subroutine HUNT (Numerical Recipes)
c ------------------------------------------------------------------ HSJ
c
      implicit none
c
      integer n,jlo,jhi,jmid,inc
      real    x,xx(n)
c
      if (x .lt. xx(1)) then
        jlo=0                       ! indicates out of range below xx(1)
      else if (x .le. xx(2)) then
        jlo=1                       ! xx(1) .le. x .le. xx(2)
      else if (x .le. xx(n)) then   ! xx(2) .lt. x .le. xx(n)
c
c         check if jlo from previous call is usable
c
          if (jlo .eq. 0 .or. jlo .eq. n) then
              jlo=2
              jhi=n
              go to 15  ! no correlation go directly to bisection
          else          ! 1 .le. jlo .lt. n
c
c             bracket x, then use bisection
c             start with jlo from previous call
c
              inc=1
              if (x .gt. xx(jlo)) then ! search up
    4             jhi=jlo+inc
                  if (jhi .ge. n) then
                      jhi=n
                  else if (x .gt. xx(jhi)) then
                      inc=inc+inc
                      jlo=jhi
                      go to 4
                  end if
              else
    5             jhi=jlo
                  jlo=jlo-inc
                  if (jlo .le. 1) then
                      jlo=1
                  else if (x .le. xx(jlo)) then
                      inc=inc+inc
                      go to 5
                  end if
              end if
c
          end if
c
c         bisection
c
   10     if (jhi-jlo .eq. 1)  return
   15     jmid=(jhi+jlo)/2
          if (x .gt. xx(jmid)) then
              jlo=jmid
          else  ! x .le. xx(jmid)
              jhi=jmid
          end if
          go to 10
c
      else  ! x .gt. xx(n)
          jlo=n
      end if
c
      return
c
      end




ccccccccccccccccccccccc
c=======================================================================
      subroutine spline1dx(ynew,xnew,nnew,yold,xold,nold,y2old)
c-----------------------------------------------------------------------
c     use 1d cubic spline on yold[xold] to produce ynew[xnew]
c     y2old(1:nold) is a work array
c     ynew(1:nnew) is the output
c-----------------------------------------------------------------------
      implicit none
      double precision ynew(*),yold(*),xnew(*),xold(*)
      double precision y2old(*)
      double precision yp1,ypn
      integer nnew,nold,i
      yp1=-2.e30
      ypn=-2.e30
      call splinex(xold,yold,nold,yp1,ypn,y2old)
      do i=1,nnew
         call splintx(xold,yold,y2old,nold,xnew(i),ynew(i))
      end do
      return
      end
c=======================================================================
      subroutine splinex(x,y,n,yp1,ypn,y2)
c-----------------------------------------------------------------------
c     spline routine based upon numerical recipes
c     this is the setup routine which needs to be called only once
c     splines y as a function of x--both arrays have n elements
c     yp1 and ypn are boundary conditions on the spline
c     yp1=y'[x] at x=x[1]
c     if yp1>=1.e30 then y''[x]=0 at x=x[1] is used
c     if yp1<=-1.e30 then y'[x[1]] is calculated from first four points
c     ypn=y'[x] at x=x[n]
c     if ypn>=1.e30 then y''[x]=0 at x=x[n] is used
c     if ypn<=-1.e30 then y'[x[n]] is calculated from last four points
c     y2[1:n] is calculated array of the second derivatives of the
c     interpolating function at the x[i]
c-----------------------------------------------------------------------
      implicit none
      integer nmax,n
      parameter (nmax=5000)
      double precision x(*),y(*),y2(*),yp1,ypn
      double precision u(nmax)
      double precision yp1t,ypnt,sig,p,qn,un
      integer k,i
      if (n .gt. nmax)then
         write(*,*) 'spline: n=',n,' nmax=',nmax
         stop 'spline;  dimensional error '
      endif
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else if(yp1.lt.-.99e30) then
         yp1t=(3*x(1)**2+x(2)*x(3)+
     $        x(2)*x(4)+x(3)*x(4)-2*x(1)*(x(2)+x(3)+x(4)))*
     $ y(1)/((x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4)))+
     $ (-x(1)**2+x(1)*x(3)+x(1)*x(4)-x(3)*x(4))*y(2)/
     $ ((x(1)-x(2))*(x(2)-x(3))*(x(2)-x(4)))+
     $ (x(1)**2-x(1)*x(2)-x(1)*x(4)+x(2)*x(4))*y(3)/
     $ ((x(1)-x(3))*(x(2)-x(3))*(x(3)-x(4)))+
     $ (-x(1)**2+x(1)*x(2)+x(1)*x(3)-x(2)*x(3))*y(4)/
     $ ((x(1)-x(4))*(x(2)-x(4))*(x(3)-x(4)))
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1t)
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do  i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else if(ypn.lt.-.99e30) then
         ypnt=(-(x(-2+n)*x(-1+n))+x(-2+n)*x(n)+x(-1+n)*x(n)-x(n)**2)*
     $  y(-3+n)/
     $  ((-x(-3+n)+x(-2+n))*(-x(-3+n)+x(-1+n))*
     $  (-x(-3+n)+x(n)))+
     $  (x(-3+n)*x(-1+n)-x(-3+n)*x(n)-x(-1+n)*x(n)+x(n)**2)*
     $  y(-2+n)/
     $  ((-x(-3+n)+x(-2+n))*(-x(-2+n)+x(-1+n))*
     $  (-x(-2+n)+x(n)))+
     $  (-(x(-3+n)*x(-2+n))+x(-3+n)*x(n)+x(-2+n)*x(n)-x(n)**2)*
     $  y(-1+n)/
     $  ((-x(-3+n)+x(-1+n))*(-x(-2+n)+x(-1+n))*
     $  (-x(-1+n)+x(n)))+
     $        (x(-3+n)*x(-2+n)+x(-3+n)*x(-1+n)+x(-2+n)*x(-1+n)-
     $        2*(x(-3+n)+x(-2+n)+x(-1+n))*x(n)+3*x(n)**2)*y(n)/
     $        ((-x(-3+n)+x(n))*(-x(-2+n)+x(n))*(-x(-1+n)+x(n)))
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypnt-(y(n)-y(n-1))/(x(n)-x(n-1)))
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do  k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      return
      end
c=======================================================================
      subroutine splintx(xa,ya,y2a,n,x,y)
c-----------------------------------------------------------------------
c     cubic spline evaluator--spline must be called first to evaluate
c     y2a
c     ya is a function of xa--both are arrays of length n
c     ya2[1:n] contains spline coefficients calculated in spline
c     x is the argument of y[x] where y is to be evaluated
c     y=y[x] is the returned value
c-----------------------------------------------------------------------
      implicit none
      double precision xa(*),ya(*),y2a(*)
      double precision x,y
      double precision a,b,h
      integer n,klo,khi,k
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) then
         write(6,*) 'bad xa input.'
         stop
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end




