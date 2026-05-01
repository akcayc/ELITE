      subroutine spline2d(fun,x,y,nx,ny,kx,coef)
c-----------------------------------------------------------------------
c     setup routine for bicubic spline of fun[x,y]
c     output of this routine is nx*ny spline coefficients 
c     stored in coef(kx,ny)
c-----------------------------------------------------------------------
      implicit none
      double precision fun(*),x(*),y(*),coef(*)
      integer kx,nx,ny
      integer ind,j
      do j=1,ny
         ind=(j-1)*kx+1
         call spline(x,fun(ind),nx,-1.e30,-1.e30,coef(ind))
      end do
      return
      end
c=======================================================================
      subroutine spline2dt(fun_new,x_new,y_new,kx_new,nx_new,ny_new,
     $     fun_old,x_old,y_old,kx_old,nx_old,ny_old,coef)
c-----------------------------------------------------------------------
c     evaluates bicubic spline: determines fun_new[x_new,y_new]
c     after spline2d has been called to determine coef
c     in the calling program fun_new,run_old, and coef are 2d arrays:
c     fun_old(kx_old,ky_old),fun_new(kx_new,ky_new),coef(kx_old,ky_old)
c-----------------------------------------------------------------------
      implicit none
      double precision fun_new(*),x_new(*),y_new(*)
      double precision fun_old(*),x_old(*),y_old(*)
      double precision coef(*)
      integer kx_new,nx_new,ny_new
      integer kx_old,nx_old,ny_old
      integer i,j,ky_old,ind
      parameter(ky_old=4000)
      double precision ftemp(ky_old),ctemp(ky_old)
      if(ny_old.gt.ky_old) then
         write(6,'("dimensioning problem in spline2dt")')
         stop
      end if
      do i=1,nx_new
         do j=1,ny_old
            ind=(j-1)*kx_old+1
            call splint(x_old,fun_old(ind),coef(ind),nx_old,x_new(i),
     $           ftemp(j))
         end do
         call  spline(y_old,ftemp,ny_old,-1.e30,-1.e30,ctemp)
         do j=1,ny_new
            call splint(y_old,ftemp,ctemp,ny_old,y_new(j),
     $           fun_new((j-1)*kx_new+i))
         end do
      end do
      return
      end
c=======================================================================
      subroutine spline1d(ynew,xnew,nnew,yold,xold,nold,y2old)
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
      yp1=-1.e30
      ypn=-1.e30
      call spline(xold,yold,nold,yp1,ypn,y2old)
      do i=1,nnew
         call splint(xold,yold,y2old,nold,xnew(i),ynew(i))
      end do
      return
      end
c=======================================================================
      subroutine spline(x,y,n,yp1,ypn,y2)
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
      parameter (nmax=4000)
      double precision x(*),y(*),y2(*),yp1,ypn
      double precision u(nmax)
      double precision yp1t,ypnt,sig,p,qn,un
      integer k,i
      if (n .gt. nmax)then
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
      subroutine splint(xa,ya,y2a,n,x,y)
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
c=======================================================================
      subroutine zspline(xa,ya,y2a,n,zg,ng,za)
c-----------------------------------------------------------------------
c     zspline integrates the cubic spline of ya[xa]
c     it assumes that spline has already been called to evaluate y2a
c     xa[1:n], ya[1:n], y2a[1:n]
c     In Mathematica notation z[i]=zg+Integrate[y[x],{x,x[ng],x[i]}]
c     where both zg and ng are input quantities,
c     za is calculated here. 
c     za can then be used in zsplint to determine z at a 
c     specific x
c-----------------------------------------------------------------------
      implicit none
c
      double precision xa(*),ya(*),y2a(*),zg,za(*)
      integer n,ng
c
      integer j
      double precision const
c
      if (ng .lt. 0 .or. ng .gt. n)stop 'zspline: wrong ng'
c
      za(1)=0.
      do j=2,n
         za(j)=za(j-1)+0.5*(xa(j)-xa(j-1))*(ya(j)+ya(j-1))-
     >        (xa(j)-xa(j-1))**3*(y2a(j)+y2a(j-1))/24.0
      enddo
c
      const=zg-za(ng)
      do j=1,n
         za(j)=za(j)+const
      enddo
c
      return
      end
c=======================================================================
      subroutine zsplint(xa,ya,y2a,za,n,x,y,yp,z)
c-----------------------------------------------------------------------
c     evaluate cubic spline to determine function (y),
c     derivative (yp) and integral (z) at location x.
c     first spline must be called to obtain y2a and
c     zspline must be called to obtain za
c-----------------------------------------------------------------------
      implicit none
      double precision xa(*),ya(*),y2a(*),za(*)
      integer n
      double precision x,y,yp,z
      integer klo,khi,k
      double precision h,a,b
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
     *     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      yp=(ya(khi)-ya(klo))/h-(3.*a**2-1.)/6.*h*y2a(klo)
     >     +(3.*b**2-1.)/6.*h*y2a(khi)
      z=za(klo)+0.5*h*ya(klo)*(1.-a**2)+0.5*ya(khi)*h*b**2+
     >     y2a(klo)*h**3*(2.*a**2-a**4-1.)/24.-
     >     y2a(khi)*h**3*(2.*b**2-b**4)/24.
      return
      end

      SUBROUTINE LZHES(N, A, NA, B, NB, X, NX, WANTX)                   LZH   10
C THIS SUBROUTINE REDUCES THE COMPLEX MATRIX A TO UPPER
C HESSENBERG FORM AND REDUCES THE COMPLEX MATRIX B TO
C TRIANGULAR FORM
C INPUT PARAMETERS..
C N   THE ORDER OF THE A AND B MATRICES
C A   A COMPLEX MATRIX
C NA  THE ROW DIMENSION OF THE A MATRIX
C B   A COMPLEX MATRIX
C NB  THE ROW DIMENSION OF THE B MATRIX
C NX  THE ROW DIMENSION OF THE X MATRIX
C WANTX A LOGICAL VARIABLE WHICH IS SET TO .TRUE. IF
C       THE EIGENVECTORS ARE WANTED. OTHERWISE IT SHOULD
C     BE SET TO .FALSE.
C OUTPUT PARAMETERS..
C A  ON OUTPUT A IS AN UPPER HESSENBERG MATRIX, THE
C    ORIGINAL MATRIX HAS BEEN DESTROYED
C B  AN UPPER TRIANGULAR MATRIX, THE ORIGINAL MATRIX
C    HAS BEEN DESTROYED
C X  CONTAINS THE TRANSFORMATIONS NEEDED TO COMPUTE
C    THE EIGENVECTORS OF THE ORIGINAL SYSTEM
      COMPLEX Y, W, Z, A(NA,N), B(NB,N), X(NX,N)
      double precision C, D
      LOGICAL WANTX
      NM1 = N - 1
C REDUCE B TO TRIANGULAR FORM USING ELEMENTARY
C TRANSFORMATIONS
      DO 80 I=1,NM1
        D = 0.d0
        IP1 = I + 1
        DO 10 K=IP1,N
          Y = B(K,I)
          C = DABS(REAL(Y)) + DABS(IMAG(Y))
          IF (C.LE.D) GO TO 10
          D = C
          II = K
   10   CONTINUE
        IF (D.EQ.0.d0) GO TO 80
        Y = B(I,I)
        IF (D.LE.ABS(REAL(Y))+ABS(IMAG(Y))) GO TO 40
C MUST INTERCHANGE
        DO 20 J=1,N
          Y = A(I,J)
          A(I,J) = A(II,J)
          A(II,J) = Y
   20   CONTINUE
        DO 30 J=I,N
          Y = B(I,J)
          B(I,J) = B(II,J)
          B(II,J) = Y
   30   CONTINUE
   40   DO 70 J=IP1,N
          Y = B(J,I)/B(I,I)
          IF (REAL(Y).EQ.0.d0 .AND. IMAG(Y).EQ.0.d0) GO TO 70
          DO 50 K=1,N
            A(J,K) = A(J,K) - Y*A(I,K)
   50     CONTINUE
          DO 60 K=IP1,N
            B(J,K) = B(J,K) - Y*B(I,K)
   60     CONTINUE
   70   CONTINUE
        B(IP1,I) = (0.d0,0.d0)
   80 CONTINUE
C INITIALIZE X
      IF (.NOT.WANTX) GO TO 110
      DO 100 I=1,N
        DO 90 J=1,N
          X(I,J) = (0.d0,0.d0)
   90   CONTINUE
        X(I,I) = (1.d0,0.d0)
  100 CONTINUE
C REDUCE A TO UPPER HESSENBERG FORM
  110 NM2 = N - 2
      IF (NM2.LT.1) GO TO 270
      DO 260 J=1,NM2
        JM2 = NM1 - J
        JP1 = J + 1
        DO 250 II=1,JM2
          I = N + 1 - II
          IM1 = I - 1
          IMJ = I - J
          W = A(I,J)
          Z = A(IM1,J)
          IF (ABS(REAL(W))+ABS(IMAG(W)).LE.ABS(REAL(Z))
     *     +ABS(IMAG(Z))) GO TO 140
C MUST INTERCHANGE ROWS
          DO 120 K=J,N
            Y = A(I,K)
            A(I,K) = A(IM1,K)
            A(IM1,K) = Y
  120     CONTINUE
          DO 130 K=IM1,N
            Y = B(I,K)
            B(I,K) = B(IM1,K)
            B(IM1,K) = Y
  130     CONTINUE
  140     Z = A(I,J)
          IF (REAL(Z).EQ.0.d0 .AND. IMAG(Z).EQ.0.d0) GO TO 170
          Y = Z/A(IM1,J)
          DO 150 K=JP1,N
            A(I,K) = A(I,K) - Y*A(IM1,K)
  150     CONTINUE
          DO 160 K=IM1,N
            B(I,K) = B(I,K) - Y*B(IM1,K)
  160     CONTINUE
C TRANSFORMATION FROM THE RIGHT
  170     W = B(I,IM1)
          Z = B(I,I)
          IF (ABS(REAL(W))+ABS(IMAG(W)).LE.ABS(REAL(Z))
     *     +ABS(IMAG(Z))) GO TO 210
C MUST INTERCHANGE COLUMNS
          DO 180 K=1,I
            Y = B(K,I)
            B(K,I) = B(K,IM1)
            B(K,IM1) = Y
  180     CONTINUE
          DO 190 K=1,N
            Y = A(K,I)
            A(K,I) = A(K,IM1)
            A(K,IM1) = Y
  190     CONTINUE
          IF (.NOT.WANTX) GO TO 210
          DO 200 K=IMJ,N
            Y = X(K,I)
            X(K,I) = X(K,IM1)
            X(K,IM1) = Y
  200     CONTINUE
  210     Z = B(I,IM1)
          IF (REAL(Z).EQ.0.d0 .AND. IMAG(Z).EQ.0.d0) GO TO 250
          Y = Z/B(I,I)
          DO 220 K=1,IM1
            B(K,IM1) = B(K,IM1) - Y*B(K,I)
  220     CONTINUE
          B(I,IM1) = (0.d0,0.d0)
          DO 230 K=1,N
            A(K,IM1) = A(K,IM1) - Y*A(K,I)
  230     CONTINUE
          IF (.NOT.WANTX) GO TO 250
          DO 240 K=IMJ,N
            X(K,IM1) = X(K,IM1) - Y*X(K,I)
  240     CONTINUE
  250   CONTINUE
        A(JP1+1,J) = (0.d0,0.d0)
  260 CONTINUE
  270 RETURN
      END
      SUBROUTINE LZIT(N, A, NA, B, NB, X, NX, WANTX, ITER, EIGA,        LZI   10
     * EIGB)
C THIS SUBROUTINE SOLVES THE GENERALIZED EIGENVALUE PROBLEM
C              A X  = LAMBDA B X
C WHERE A IS A COMPLEX UPPER HESSENBERG MATRIX OF
C ORDER N AND B IS A COMPLEX UPPER TRIANGULAR MATRIX OF ORDER N
C INPUT PARAMETERS
C N      ORDER OF A AND B
C A      AN N X N UPPER HESSENBERG COMPLEX MATRIX
C NA     THE ROW DIMENSION OF THE A MATRIX
C B      AN N X N UPPER TRIANGULAR COMPLEX MATRIX
C NB     THE ROW DIMENSION OF THE B MATRIX
C X      CONTAINS TRANSFORMATIONS TO OBTAIN EIGENVECTORS OF
C        ORIGINAL SYSTEM. IF EIGENVECTORS ARE REQUESTED AND QZHES
C        IS NOT CALLED, X SHOULD BE SET TO THE IDENTITY MATRIX
C NX     THE ROW DIMENSION OF THE X MATRIX
C WANTX  LOGICAL VARIABLE WHICH SHOULD BE SET TO .TRUE.
C        IF EIGENVECTORS ARE WANTED. OTHERWISE IT
C        SHOULD BE SET TO .FALSE.
C OUTPUT PARAMETERS
C X      THE ITH COLUMN CONTAINS THE ITH EIGENVECTOR
C        IF EIGENVECTORS ARE REQUESTED
C ITER   AN INTEGER ARRAY OF LENGTH N WHOSE ITH ENTRY
C        CONTAINS THE NUMBER OF ITERATIONS NEEDED TO FIND
C        THE ITH EIGENVALUE. FOR ANY I IF ITER(I) =-1 THEN
C        AFTER 30 ITERATIONS THERE HAS NOT BEEN A SUFFICIENT
C        DECREASE IN THE LAST SUBDIAGONAL ELEMENT OF A
C        TO CONTINUE ITERATING.
C EIGA   A COMPLEX ARRAY OF LENGTH N CONTAINING THE DIAGONAL OF A
C EIGB   A COMPLEX ARRAY OF LENGTH N CONTAINING THE DIAGONAL OF B
C THE ITH EIGENVALUE CAN BE FOUND BY DIVIDING EIGA(I) BY
C EIGB(I). WATCH OUT FOR EIGB(I) BEING ZERO
      COMPLEX A(NA,N), B(NB,N), EIGA(N), EIGB(N)
      COMPLEX S, W, Y, Z
!, CSQRT
      COMPLEX X(NX,N)
      INTEGER ITER(N)
      COMPLEX ANNM1, ALFM, BETM, D, SL, DEN, NUM, ANM1M1
      REAL EPSA, EPSB, SS, R, ANORM, BNORM, ANI, BNI, C
      REAL D0, D1, D2, E0, E1
      LOGICAL WANTX
      NN = N
C COMPUTE THE MACHINE PRECISION TIMES THE NORM OF A AND B
      ANORM = 0.d0
      BNORM = 0.d0
      DO 30 I=1,N
        ANI = 0.d0
        IF (I.EQ.1) GO TO 10
        Y = A(I,I-1)
        ANI = ANI + ABS(REAL(Y)) + ABS(IMAG(Y))
   10   BNI = 0.d0
        DO 20 J=I,N
          ANI = ANI + ABS(REAL(A(I,J))) + ABS(IMAG(A(I,J)))
          BNI = BNI + ABS(REAL(B(I,J))) + ABS(IMAG(B(I,J)))
   20   CONTINUE
        IF (ANI.GT.ANORM) ANORM = ANI
        IF (BNI.GT.BNORM) BNORM = BNI
   30 CONTINUE
      IF (ANORM.EQ.0.) ANORM = 1.d0
      IF (BNORM.EQ.0.) BNORM = 1.d0
      EPSB = BNORM
      EPSA = ANORM
   40 EPSA = EPSA/2.d0
      EPSB = EPSB/2.d0
      C = ANORM + EPSA
      IF (C.GT.ANORM) GO TO 40
      IF (N.LE.1) GO TO 320
   50 ITS = 0
      NM1 = NN - 1
C CHECK FOR NEGLIGIBLE SUBDIAGONAL ELEMENTS
   60 D2 = ABS(REAL(A(NN,NN))) + ABS(IMAG(A(NN,NN)))
      DO 70 LB=2,NN
        L = NN + 2 - LB
        SS = D2
        Y = A(L-1,L-1)
        D2 = ABS(REAL(Y)) + ABS(IMAG(Y))
        SS = SS + D2
        Y = A(L,L-1)
        R = SS + ABS(REAL(Y)) + ABS(IMAG(Y))
        IF (R.EQ.SS) GO TO 80
   70 CONTINUE
      L = 1
   80 IF (L.EQ.NN) GO TO 320
      IF (ITS.LT.30) GO TO 90
      ITER(NN) = -1
      IF (ABS(REAL(A(NN,NM1)))+ABS(IMAG(A(NN,NM1))).GT.0.8d0*
     * ABS(REAL(ANNM1))+ABS(IMAG(ANNM1))) RETURN
   90 IF (ITS.EQ.10 .OR. ITS.EQ.20) GO TO 110
C COMPUTE SHIFT AS EIGENVALUE OF LOWER 2 BY 2
      ANNM1 = A(NN,NM1)
      ANM1M1 = A(NM1,NM1)
      S = A(NN,NN)*B(NM1,NM1) - ANNM1*B(NM1,NN)
      W = ANNM1*B(NN,NN)*(A(NM1,NN)*B(NM1,NM1)-B(NM1,NN)*ANM1M1)
      Y = (ANM1M1*B(NN,NN)-S)/2.d0
!      Z = CSQRT(Y*Y+W)
      z=sqrt(y*y+w)
      IF (REAL(Z).EQ.0.d0 .AND. IMAG(Z).EQ.0.d0) GO TO 100
      D0 = REAL(Y/Z)
      IF (D0.LT.0.0) Z = -Z
  100 DEN = (Y+Z)*B(NM1,NM1)*B(NN,NN)
      IF (REAL(DEN).EQ.0.d0 .AND. IMAG(DEN).EQ.0.d0) DEN =
     * CMPLX(EPSA,0.d0)
      NUM = (Y+Z)*S - W
      GO TO 120
C AD-HOC SHIFT
  110 Y = A(NM1,NN-2)
      NUM = CMPLX(ABS(REAL(ANNM1))+ABS(IMAG(ANNM1)),ABS(REAL(Y))
     * +ABS(IMAG(Y)))
      DEN = (1.d0,0.d0)
C CHECK FOR 2 CONSECUTIVE SMALL SUBDIAGONAL ELEMENTS
  120 IF (NN.EQ.L+1) GO TO 140
      D2 = ABS(REAL(A(NM1,NM1))) + ABS(IMAG(A(NM1,NM1)))
      E1 = ABS(REAL(ANNM1)) + ABS(IMAG(ANNM1))
      D1 = ABS(REAL(A(NN,NN))) + ABS(IMAG(A(NN,NN)))
      NL = NN - (L+1)
      DO 130 MB=1,NL
        M = NN - MB
        E0 = E1
        Y = A(M,M-1)
        E1 = ABS(REAL(Y)) + ABS(IMAG(Y))
        D0 = D1
        D1 = D2
        Y = A(M-1,M-1)
        D2 = ABS(REAL(Y)) + ABS(IMAG(Y))
        Y = A(M,M)*DEN - B(M,M)*NUM
        D0 = (D0+D1+D2)*(ABS(REAL(Y))+ABS(IMAG(Y)))
        E0 = E0*E1*(ABS(REAL(DEN))+ABS(IMAG(DEN))) + D0
        IF (E0.EQ.D0) GO TO 150
  130 CONTINUE
  140 M = L
  150 CONTINUE
      ITS = ITS + 1
      W = A(M,M)*DEN - B(M,M)*NUM
      Z = A(M+1,M)*DEN
      D1 = ABS(REAL(Z)) + ABS(IMAG(Z))
      D2 = ABS(REAL(W)) + ABS(IMAG(W))
C FIND L AND M AND SET A=LAM AND B=LBM
C     NP1 = N + 1
      LOR1 = L
      NNORN = NN
      IF (.NOT.WANTX) GO TO 160
      LOR1 = 1
      NNORN = N
  160 DO 310 I=M,NM1
        J = I + 1
C FIND ROW TRANSFORMATIONS TO RESTORE A TO
C UPPER HESSENBERG FORM. APPLY TRANSFORMATIONS
C TO A AND B
        IF (I.EQ.M) GO TO 170
        W = A(I,I-1)
        Z = A(J,I-1)
        D1 = ABS(REAL(Z)) + ABS(IMAG(Z))
        D2 = ABS(REAL(W)) + ABS(IMAG(W))
        IF (D1.EQ.0.d0) GO TO 60
  170   IF (D2.GT.D1) GO TO 190
C MUST INTERCHANGE ROWS
        DO 180 K=I,NNORN
          Y = A(I,K)
          A(I,K) = A(J,K)
          A(J,K) = Y
          Y = B(I,K)
          B(I,K) = B(J,K)
          B(J,K) = Y
  180   CONTINUE
        IF (I.GT.M) A(I,I-1) = A(J,I-1)
        IF (D2.EQ.0.d0) GO TO 220
C THE SCALING OF W AND Z IS DESIGNED TO AVOID A DIVISION BY ZERO
C WHEN THE DENOMINATOR IS SMALL
        Y = CMPLX(REAL(W)/D1,IMAG(W)/D1)/CMPLX(REAL(Z)/D1,IMAG(Z)/
     *   D1)
        GO TO 200
  190   Y = CMPLX(REAL(Z)/D2,IMAG(Z)/D2)/CMPLX(REAL(W)/D2,IMAG(W)/
     *   D2)
  200   DO 210 K=I,NNORN
          A(J,K) = A(J,K) - Y*A(I,K)
          B(J,K) = B(J,K) - Y*B(I,K)
  210   CONTINUE
  220   IF (I.GT.M) A(J,I-1) = (0.d0,0.d0)
C PERFORM TRANSFORMATIONS FROM RIGHT TO RESTORE B TO
C   TRIANGLULAR FORM
C APPLY TRANSFORMATIONS TO A
        Z = B(J,I)
        W = B(J,J)
        D2 = ABS(REAL(W)) + ABS(IMAG(W))
        D1 = ABS(REAL(Z)) + ABS(IMAG(Z))
        IF (D1.EQ.0.d0) GO TO 60
        IF (D2.GT.D1) GO TO 270
C MUST INTERCHANGE COLUMNS
        DO 230 K=LOR1,J
          Y = A(K,J)
          A(K,J) = A(K,I)
          A(K,I) = Y
          Y = B(K,J)
          B(K,J) = B(K,I)
          B(K,I) = Y
  230   CONTINUE
        IF (I.EQ.NM1) GO TO 240
        Y = A(J+1,J)
        A(J+1,J) = A(J+1,I)
        A(J+1,I) = Y
  240   IF (.NOT.WANTX) GO TO 260
        DO 250 K=1,N
          Y = X(K,J)
          X(K,J) = X(K,I)
          X(K,I) = Y
  250   CONTINUE
  260   IF (D2.EQ.0.d0) GO TO 310
        Z = CMPLX(REAL(W)/D1,IMAG(W)/D1)/CMPLX(REAL(Z)/D1,IMAG(Z)/
     *   D1)
        GO TO 280
  270   Z = CMPLX(REAL(Z)/D2,IMAG(Z)/D2)/CMPLX(REAL(W)/D2,IMAG(W)/
     *   D2)
  280   DO 290 K=LOR1,J
          A(K,I) = A(K,I) - Z*A(K,J)
          B(K,I) = B(K,I) - Z*B(K,J)
  290   CONTINUE
        B(J,I) = (0.d0,0.d0)
        IF (I.LT.NM1) A(I+2,I) = A(I+2,I) - Z*A(I+2,J)
        IF (.NOT.WANTX) GO TO 310
        DO 300 K=1,N
          X(K,I) = X(K,I) - Z*X(K,J)
  300   CONTINUE
  310 CONTINUE
      GO TO 60
  320 CONTINUE
      EIGA(NN) = A(NN,NN)
      EIGB(NN) = B(NN,NN)
      IF (NN.EQ.1) GO TO 330
      ITER(NN) = ITS
      NN = NM1
      IF (NN.GT.1) GO TO 50
      ITER(1) = 0
      GO TO 320
C FIND EIGENVECTORS USING B FOR INTERMEDIATE STORAGE
  330 IF (.NOT.WANTX) RETURN
      M = N
  340 CONTINUE
      ALFM = A(M,M)
      BETM = B(M,M)
      B(M,M) = (1.d0,0.d0)
      L = M - 1
      IF (L.EQ.0) GO TO 370
  350 CONTINUE
      L1 = L + 1
      SL = (0.d0,0.d0)
      DO 360 J=L1,M
        SL = SL + (BETM*A(L,J)-ALFM*B(L,J))*B(J,M)
  360 CONTINUE
      Y = BETM*A(L,L) - ALFM*B(L,L)
      IF (REAL(Y).EQ.0.d0 .AND. IMAG(Y).EQ.0.d0) Y =
     * CMPLX((EPSA+EPSB)/2.d0,0.d0)
      B(L,M) = -SL/Y
      L = L - 1
  370 IF (L.GT.0) GO TO 350
      M = M - 1
      IF (M.GT.0) GO TO 340
C TRANSFORM TO ORIGINAL COORDINATE SYSTEM
      M = N
  380 CONTINUE
      DO 400 I=1,N
        S = (0.d0,0.d0)
        DO 390 J=1,M
          S = S + X(I,J)*B(J,M)
  390   CONTINUE
        X(I,M) = S
  400 CONTINUE
      M = M - 1
      IF (M.GT.0) GO TO 380
C NORMALIZE SO THAT LARGEST COMPONENT = 1.
      M = N
  410 CONTINUE
      SS = 0.
      DO 420 I=1,N
        R = ABS(REAL(X(I,M))) + ABS(IMAG(X(I,M)))
        IF (R.LT.SS) GO TO 420
        SS = R
        D = X(I,M)
  420 CONTINUE
      IF (SS.EQ.0.d0) GO TO 440
      DO 430 I=1,N
        X(I,M) = X(I,M)/D
  430 CONTINUE
  440 M = M - 1
      IF (M.GT.0) GO TO 410
      RETURN
      END

      subroutine zgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(1),job
      complex*16 a(lda,1),det(2),work(1)
c
c     zgedi computes the determinant and inverse of a matrix
c     using the factors computed by zgeco or zgefa.
c
c     on entry
c
c        a       complex*16(lda, n)
c                the output from zgeco or zgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from zgeco or zgefa.
c
c        work    complex*16(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     complex*16(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. cabs1(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if zgeco has set rcond .gt. 0.0 or zgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,zswap
c     fortran dabs,dcmplx,mod
c
c     internal variables
c
      complex*16 t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = (1.0d0,0.0d0)
         det(2) = (0.0d0,0.0d0)
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (cabs1(det(1)) .eq. 0.0d0) go to 60
 10                if (cabs1(det(1)) .ge. 1.0d0) go to 20
               det(1) = dcmplx(ten,0.0d0)*det(1)
               det(2) = det(2) - (1.0d0,0.0d0)
            go to 10
 20                continue
 30                       if (cabs1(det(1)) .lt. ten) go to 40
               det(1) = det(1)/dcmplx(ten,0.0d0)
               det(2) = det(2) + (1.0d0,0.0d0)
            go to 30
 40                continue
 50                 continue
 60                     continue
 70                      continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = (1.0d0,0.0d0)/a(k,k)
            t = -a(k,k)
            call zscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = (0.0d0,0.0d0)
               call zaxpy(k,t,a(1,k),1,a(1,j),1)
 80                continue
 90                       continue
 100                       continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = (0.0d0,0.0d0)
 110               continue
            do 120 j = kp1, n
               t = work(j)
               call zaxpy(n,t,a(1,j),1,a(1,k),1)
 120               continue
            l = ipvt(k)
            if (l .ne. k) call zswap(n,a(1,k),1,a(1,l),1)
 130         continue
 140             continue
 150              continue
      return
      end

      subroutine zgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      complex*16 a(lda,1)
c
c     zgefa factors a complex*16 matrix by gaussian elimination.
c
c     zgefa is usually called by zgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
c
c     on entry
c
c        a       complex*16(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that zgesl or zgedi will divide by zero
c                     if called.  use  rcond  in zgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,izamax
c     fortran dabs
c
c     internal variables
c
      complex*16 t
      integer izamax,j,k,kp1,l,nm1
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = izamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
 10                   continue
c
c           compute multipliers
c
            t = -(1.0d0,0.0d0)/a(k,k)
            call zscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
 20                         continue
               call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
 30                continue
         go to 50
 40          continue
            info = k
 50             continue
 60           continue
 70            continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0d0) info = n
      return
      end


