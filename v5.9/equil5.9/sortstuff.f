      subroutine sortr(xp,zp,gsq,npc,xmx,zmx,xset,zset,j,jsep)
      dimension xp(1),zp(1),gsq(1)
c ---------------------------------------
c order points counterclockwise if needed
c ---------------------------------------
      i=0
   10 i=i+1
      dr1=xp(i)-xmx
      dr2=xp(i+1)-xmx
      dz1=zp(i)-zmx
      dz2=zp(i+1)-zmx
      dot=dr1*dr2+dz1*dz2
      cross=dr1*dz2-dr2*dz1
      if(cross.eq.0.) go to 10
      if(dot.eq.0.) go to 10
      if(cross/dot.gt.0.) go to 20
c ------------------------------------------------
c reorder the points to make them counterclockwise
c ------------------------------------------------
      call sorter(xp,zp,gsq,npc)
   20 continue
c ---------------------
c delete the last point
c ---------------------
      npc=npc-1
c --------------------------------------------------------------------
c sort so first point is one with smallest positive angle about center
c xmx, zmx for j .lt. jsep
c pick point closest to xset, zset for j .ge. jsep
c --------------------------------------------------------------------
      call sortog(xp,zp,gsq,npc,xmx,zmx,xset,zset,j,jsep)
c -----------------------------------------
c reset new last point to equal first point
c -----------------------------------------
      npc=npc+1
      xp(npc)=xp(1)
      zp(npc)=zp(1)
      gsq(npc)=gsq(1)
      return
      end
      subroutine sorter(x,y,d,n)
cray  lcm (x),(y),(d)
        dimension x(1),y(1),d(1)
        nsort=n/2
        nl=n+1
        do 100 i=1,nsort
        nl=nl-1
        t=x(i)
        x(i)=x(nl)
        x(nl)=t
        t=y(i)
        y(i)=y(nl)
        y(nl)=t
        t=d(i)
        d(i)=d(nl)
        d(nl)=t
  100   continue
        return
      end
      subroutine sortog(x,y,d,n,xs,ys,xset,zset,j,jsep)
cray  lcm (x),(y),(d)
        dimension x(1),y(1),d(1)
      data pi/3.1415926535898/
      if(j.ge.jsep) go to 110
      ang=1.e10
        do 50 i=1,n
      angp=atan2(y(i)-ys,x(i)-xs)
      angp=abs(angp)
      if(angp.gt.ang) go to 50
   30   isp=i
      ang=angp
   50   continue
   60 if(isp.eq.1) return
        im1=isp-1
        nm1=n-1
        do 100 i=1,im1
        tx=x(1)
        ty=y(1)
        td=d(1)
        do 90 jj=1,nm1
        x(jj)=x(jj+1)
        y(jj)=y(jj+1)
        d(jj)=d(jj+1)
   90   continue
        x(n)=tx
        y(n)=ty
        d(n)=td
  100   continue
        return
c ----------------------------------------------
c sorting procedure for points inside separatrix
c find points closest to xset, zset
c ----------------------------------------------
  110 dist=1.e10
      do 130 i=1,n
      dt=(y(i)-zset)**2+(x(i)-xset)**2
      if(dt.gt.dist) go to 130
      isp=i
      dist=dt
  130 continue
      go to 60
      end

