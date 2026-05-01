
subroutine mercy(isurf)

!  1/31/00 translated to F90
!-----------------------------------------------------------------------
!  calculate mercier criterion
!  modified 9/30/97 to calculate peeling mode criterion
!-----------------------------------------------------------------------
      use elite_data, only: ns,outfile,runname,lrunname,rl,bpl,circumf, &
           pi,pprime_adj,qrefp,vpp,psiv_eq,f,vprime,shbishop,abishop, &
           shear,alpha,ffprime_adj,q_calc,analytic_peel,psigrid, &
           lhspeelstore,rhspeelstore,dmstore,nsurfdat,nxinterp,verbose
      implicit none
      real fa(ns),fb(ns),fc(ns)
      real fd(ns), fb2(ns)
      real qprim(1),fps2m(1),b2m(1),b2ps2m(1)
      real b2p(1)
      real f2b2m(1),alam(1),vpprim(1)
      real am1(1),am2(1),am3(1),am4(1),am5(1)
      real ar2(1),ar3(1),ar4(1),ar5(1)
      real am1n(1),am2n(1),am3n(1),am4n(1),am5n(1)
      real ar2n(1),ar3n(1),ar4n(1),ar5n(1)
      real hfactor(1),hcheck(1)
      real amercy(1),dmplus14(1)
      real gchi2(ns)
      real vddr(1)
      real fb3(ns),fb4(ns)
      real za(1),coeff(1)
      real rhspeel(1),lhspeel(1)
      integer njmg,nlocal,neqdata
      character*1 cideal,cres,cpeel
      real aq,bq,cq,det,sqr
      real ypn,yp1
      real suma,sumb,sumc,sumd,sume
      real sumb2,sumb3,sumb4
      real fact
      real ppn_dsk,pn_dsk
      real zthe,xthe,denom,xpsi,zpsi,psipsi
      real psix,psiz,bsquar,gpsith,gpsi2
      real thez,thex
      real pedge
      integer i,j
      real pprim
      real dum,vpmer
      real dchi
      real x(1,ns)
      real dresn(1),dres(1),didealn(1),dideal(1)
      real jacob(ns)
      real dmhat,test,dmhat2
      real rhs_c,rhs_pp,rhs_ffp
      integer isurf
!      real fluxfrac
      if (verbose .ge. 4) write(outfile,*)  "=============begin output from mercier===, surface=",isurf
!-----------------------------------------------------------------------
!  1.0 open file runname.jmg
!-----------------------------------------------------------------------
!      njmg=77
!!      open(unit=njmg,file="outjmg",status="unknown")
!     open(unit=njmg,file=runname(1:lrunname)//'.jmg',status="unknown")
!-----------------------------------------------------------------------
!  2.0  gchi2 commonly known as grad-psi-squared
!-----------------------------------------------------------------------
      i=1
      do  j=1,ns
            gchi2(j)=(rl(j)*bpl(j))**2
            jacob(j)=circumf/(2.*pi*bpl(j))
            x(i,j)=rl(j)
      end do
!-----------------------------------------------------------------------
! 3.0 need  pprime,qprim(1),vprime,vpp
!-----------------------------------------------------------------------
!      pprim=pprime/(4.*pi)
!      pprim=pprime_eq(isurf)/(4.*pi)
      pprim=pprime_adj(isurf)/(4.*pi)  ! use adjusted value 3/1/00
      qprim(1)=qrefp
      vpprim(i)=vpp
      if (verbose .ge. 4) write(outfile,*)' vpp=',vpp
!      fluxfrac=abs(psiv_eq(isurf)/psiv_eq(1))  
! fraction of outer surface flux for now, percenflux not used here
! 9/01 replace with actual psigrid read from .xtopsi

!-----------------------------------------------------------------------
! 4.0 calculate arrays for flux surface average
!     chi=(poloidal flux/2 pi) is the flux surface label
!-----------------------------------------------------------------------
      do j=1,ns
            bsquar=(gchi2(j)+f**2)/x(i,j)**2 ! note that sigl is f
            fa(j)=f/gchi2(j)*jacob(j)
            fb(j)=1./bsquar*jacob(j)
            fb2(j)=(bsquar/f)*jacob(j)
            fb3(j)=(gchi2(j)/x(i,j)**2)*jacob(j)
            fb4(j)=(f**2/x(i,j)**2)*jacob(j)
            fc(j)=bsquar/gchi2(j)*jacob(j)
            fd(j)=f**2/(bsquar*gchi2(j))*jacob(j)
      end do
!-----------------------------------------------------------------------
! 5.1 perform flux surface averages
!-----------------------------------------------------------------------
         suma=0.
         sumb=0.
         sumb2=0.
         sumb3=0.
         sumb4=0.
         sumc=0.
         sumd=0.
         sume=0.
         dchi=(2.*pi)/(ns-1.)
         do  j=2,ns
            suma=suma+0.5*(fa(j)+fa(j-1))*dchi
            sumb=sumb+0.5*(fb(j)+fb(j-1))*dchi
            sumb2=sumb2+0.5*(fb2(j)+fb2(j-1))*dchi
            sumb3=sumb3+0.5*(fb3(j)+fb3(j-1))*dchi
            sumb4=sumb4+0.5*(fb4(j)+fb4(j-1))*dchi
            sumc=sumc+0.5*(fc(j)+fc(j-1))*dchi
            sumd=sumd+0.5*(fd(j)+fd(j-1))*dchi
            sume=sume+0.5*(jacob(j)+jacob(j-1))*dchi
         end do
         fps2m(i)=suma/sume
         b2m(i)=sumb/sume
         b2p(i)=sume/sumb2      !b2p=f/<B^2>
         b2ps2m(i)=sumc/sume
         f2b2m(i)=sumd/sume
         vpmer=2.*pi*sume ! not used--just a debug check--should =vprime
!         write(6,*)"vpmer/vprime=",vpmer,vprime
!-----------------------------------------------------------------------
! 5.2 the Mercier criterion terms (ideal)
!-----------------------------------------------------------------------
         alam(i)=-(2.*pi)**2*qprim(i) ! note the minus sign
         am1(i)=alam(i)**2/4.
         am2(i)=4.*pi*pprim*vprime*alam(i)*fps2m(i)
         am3(i)=4.*pi*pprim*vprime*vpprim(i)*b2ps2m(i)
         am4(i)=-(4.*pi*pprim*vprime)**2*b2m(i)*b2ps2m(i)
         am5(i)=(4.*pi*pprim*vprime)**2* &
             (fps2m(i)**2-f2b2m(i)*b2ps2m(i))
!         write(6,*)' pp=',4.*pi*pprim,' vp=',vprime,' vpp=',vpprim(i)
!         write(6,*)' b2ps2m=',b2ps2m(i),' sumc=',sumc
!-----------------------------------------------------------------------
! there are 4 ideal criteria ( all equivalent ) defined below
! amercy > 0 is ideal stable
!-----------------------------------------------------------------------
         amercy(i)=am1(i)+am2(i)+am3(i)+am4(i)+am5(i)
!-----------------------------------------------------------------------
! dideal < 0 is ideal stable
! dmplus14=D(I)+1/4 < 1/4 is ideal stable
!-----------------------------------------------------------------------
         dmplus14(i)=-(am2(i)+am3(i)+am4(i)+am5(i))/alam(i)**2
         dideal(i)=dmplus14(i)-0.25
         dmhat=dmplus14(i)*shbishop**2/abishop
         dmhat2=dmplus14(i)*shear(isurf)**2/alpha(isurf)
!         write(6,*)' shbish=',shbishop,' abish=',abishop
         if (verbose .ge. 4) then
            write(outfile,*)' t1=',-am2(i)/alam(i)**2
            write(outfile,*)' t2=',-am3(i)/alam(i)**2
            write(outfile,*)' t3=',-am4(i)/alam(i)**2
            write(outfile,*)' t4=',-am5(i)/alam(i)**2
         endif
!         test=(k**3)*(12.+45.*k)*cos(gamma)/64.
!         test=dmhat/test
         test=0. ! no longer have k and gamma info
         if (verbose .ge. 4) then
            write(outfile,*)' Full Mercier D_m=',dmplus14(i)
            write(outfile,*) 'D_m*shear^2/alpha= ',dmhat2
!         write(93,*) dmplus14(i)
            write(outfile,'(1p,a,e12.3,a,e12.3)') "mercier: dmhat=", &
                 dmhat,"; ratio with wilson formula=",test
         endif
!         write(6,*)' qref=',qref
!     normalize so that am3n is the same as vpprim
         fact=4.*pi*pprim*vprime*b2ps2m(i)
         am1n(i)=am1(i)/fact
         am2n(i)=am2(i)/fact
         am3n(i)=am3(i)/fact
         am4n(i)=am4(i)/fact
         am5n(i)=am5(i)/fact
!-----------------------------------------------------------------------
!          apply the normalization to dideal which is v_dd/dres
! didealn < 0 is ideal stable assuming pprim < 0
!-----------------------------------------------------------------------
         didealn(i)=-dideal(i)*alam(i)**2/(fact*vprime**2)
!-----------------------------------------------------------------------
! 5.3      resistive interchange quantities
!            only ar2 and ar5 are changed
!-----------------------------------------------------------------------
         hfactor(i)=4*pi*pprim*vprime/alam(i)*(b2p(i) &
             *b2ps2m(i) - fps2m(i))
         ar2n(i)=alam(i)*b2p(i)
         ar3n(i)=vpprim(i)
         ar4n(i)=-4.*pi*pprim*vprime*b2m(i)
         ar5n(i)=4.*pi*pprim*vprime* &
           (-f2b2m(i)+2.*fps2m(i)*b2p(i)- &
           b2ps2m(i)*(b2p(i))**2)
!-----------------------------------------------------------------------
! three resistive criteria defined (all equivalent)
!   dres, dresn, vddr
!         dres=D(R)
!         dres(i)=-fact*(ar2n(i)+ar3n(i)+ar4n(i)+ar5n(i))/alam(i)**2
!-----------------------------------------------------------------------
         vddr(i)=(ar2n(i)+ar3n(i)+ar4n(i)+ar5n(i))/(vprime**2)
!     normalize so that we have the corresponding terms in D(I) and D(R)
         ar2(i)=ar2n(i)*fact
         ar3(i)=ar3n(i)*fact
         ar4(i)=ar4n(i)*fact
         ar5(i)=ar5n(i)*fact
         dres(i)=-(ar2(i)+ar3(i)+ar4(i)+ar5(i))/alam(i)**2
         dresn(i)=-dres(i)*alam(i)**2/(fact*vprime**2)
!-----------------------------------------------------------------------
! 5.4 hcheck will be zero if we did algebra correctly
!  dideal=E+F+H-1/4
!  dmplus14=E+F+H
!  dres=E+F+H^2=dideal+(H-1/2)^2=dmplus14+H^2-H
!-----------------------------------------------------------------------
         hcheck(i)=dmplus14(i)-hfactor(i)+hfactor(i)**2-dres(i)
!-----------------------------------------------------------------------
! 5.5 peeling mode
!-----------------------------------------------------------------------
         rhs_c=1.
         rhs_pp= -vprime/(2.*pi**2*qprim(i))*fps2m(i)
         rhs_ffp=-vprime/(2.*pi**2*qprim(i))/f*b2ps2m(i)
!         write(65,'(1p,3e16.8)') alam(i),am2(i),am3(i)
!         write(65,'(1p,3e16.8)') am4(i),am5(i)
!         write(65,'(1p,3e16.8)') rhs_c,rhs_pp,rhs_ffp
         rhspeel(i)=abs(1. -vprime/(2.*pi**2*qprim(i))* &
!             (pprime*fps2m(i)+ &  !check this carefully-rlm
             (pprime_adj(isurf)*fps2m(i)+ &
!             ffprime/f*b2ps2m(i)))
              ffprime_adj(isurf)/f*b2ps2m(i)))
         if(dmplus14(i).lt.0.25) then
            lhspeel(i)=sqrt(1.-4.*dmplus14(i))
         else
            lhspeel(i)=-1.
         endif
         if (verbose .ge. 2) write(outfile,*) &
              ' rhspeel=',rhspeel(i),'  lhspeel=',lhspeel(i)
!-----------------------------------------------------------------------
! 5.6 write out all stability results
!     assign s/u characters for stability
!-----------------------------------------------------------------------
!      nlocal=76
!      open(unit=nlocal,file="outlocal",status="unknown")
!      open(unit=nlocal,file=runname(1:lrunname)//'.local', &
!           status="unknown")
      if (verbose .ge. 2) write(outfile,'("  I  ","   % flux  ","   D(I)    ", &
          &"s/u","   D(R)  ","s/u"," (1-4D(I))^1/2"," peel term", &
          &" s/u")')
      if (isurf==1) then
         if (verbose .ge. 2) then
            write(6,*) ' '
            write(6,'("  I ","   % flux ","     q  ","      D(I)   ", &
                 &"s/u","    D(R) "," s/u"," (1-4D(I))^1/2"," peelterm", &
                 &" s/u")')
         endif
      endif
      if(dmplus14(i).lt.0.25) then
         cideal="s"
      else
         cideal="u"
      endif
      if(dres(i).lt.0.0) then
         cres="s"
      else
         cres="u"
      endif
      if(lhspeel(i).gt.rhspeel(i)) then
         cpeel="s"
      else
         cpeel="u"
      endif
      dum=1.0
      if (verbose .ge. 2) then 
         write(outfile,'(i4,1p,2e12.4,1x,a,e12.4,1x,a,2e12.4,1x,a)') &
              isurf,psigrid(isurf),dmplus14(i),cideal,dres(i),cres, &
              lhspeel(i),rhspeel(i),cpeel
         write(6,'(i4,1p,2e10.3,1e11.3,1x,a,e11.3,1x,a,2e12.3,2x,a)') &
              isurf,psigrid(isurf),q_calc(isurf),dmplus14(i),cideal, &
              dres(i),cres, &
              lhspeel(i),rhspeel(i),cpeel
      endif
      if (isurf==1) then
         lhspeelstore=lhspeel(isurf)
         rhspeelstore=rhspeel(isurf)
         dmstore=dmplus14(i)
      endif

! write line in *.surfdat file
      if (verbose .ge. 1) then
         write(nsurfdat,*) lhspeel,rhspeel,dmplus14(i),dres(i),dmhat2
         if (isurf == nxinterp) close(nsurfdat)
      endif

!      close(nlocal)
         
!-----------------------------------------------------------------------
!  6.0 print jmg results
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 6.1   the basic equilibrium quantities
!-----------------------------------------------------------------------
!      write(njmg,100)
 100  format(//,' *****  equilibrium quantities ***** ',//)
!      write(njmg,110)
 110  format('   i',"     psi    ","      q      ", &
         "     q'     " ,"      p'    " ,"      V'    " , &
         "      V''   "  ,/)
      dum=0.
!!      write(njmg,120) i,dum,qref,qprim(i),pprim,vprime,vpprim(i)
!      write(njmg,120) isurf,psiv_eq(isurf),q_calc(isurf),qprim(i), &
!           pprim,vprime,vpprim(i)
 120  format(i5,6e12.3)
!-----------------------------------------------------------------------
! 6.2   more equilibrium equnatities
!-----------------------------------------------------------------------
!      write(njmg,130)
 130  format(//,' ****  more equilibrium quantities *** ')     
!      write(njmg,140)
 140  format('   i',"            " ,"   fps2m    " , &
          "     b2m    " ,"   b2ps2m   "  ,"   f2b2m    " , &
          " pressure   "     ,/)
      dum=0.
!      write(njmg,150) isurf,vprime,fps2m(i),b2m(i),b2ps2m(i), &
!          f2b2m(i),dum
 150  format(i5,6e12.3)
!-----------------------------------------------------------------------
! 6.3   the D(I) and D(R) criteria and the check
!-----------------------------------------------------------------------
!      write(njmg,160)
 160  format(//,' ***** D(I)+1/4 , D(R), H , and Check *****',//)
!      write(njmg,170)
 170  format('   i',"  D(I)+1/4  " ,"   D(R)     " , &
        " V(dbl dag) " ,"  H factor  ","   hcheck   ",/)
!      write(njmg,180) isurf,dmplus14(i),dres(i),vddr(i), &
!              hfactor(i),hcheck(i)
 180  format(i5,5e12.3)
!-----------------------------------------------------------------------
! 6.4   first normalization for quantities in D(I)+1/4 
!-----------------------------------------------------------------------
!      write(njmg,190)
 190  format(//,' ***** the first normalization for D(I)****', &
         /,' ***** D(I)= -(am1+am2+am3+am4+am5)/4./am1   *****',//)
!      write(njmg,200)
 200  format('   i',"     am1    " ,"     am2    ", &
          "     am3    " ,"     am4    ","     am5    ",/)
!      write(njmg,180) isurf,am1(i),am2(i),am3(i),am4(i),am5(i)
!-----------------------------------------------------------------------
! 6.5  second normalization for quantities in D(I)+1/4
!-----------------------------------------------------------------------
!      write(njmg,210)
 210  format(//,' ***** the second normatlization for D(I) *****', &
         /, ' ***** am3n is normalized to be the same as V" *****',//)
!      write(njmg,220)
 220  format('   i',"    am1n    " ,"    am2n    ", &
         "  am3n(V')  " ,"    am4n    ","    am5n    ",/)
!      write(njmg,180) isurf,am1n(i),am2n(i),am3n(i),am4n(i),am5n(i)
!-----------------------------------------------------------------------
! 6.6  first normazliation for quantities in D(R)
!-----------------------------------------------------------------------
!      write(njmg,230)
 230  format(//,' ***** the first normalization for D(R)****', &
        /,' ***** D(R)= -(ar2+ar3+ar4+ar5)/4./am1     *****',//)
!      write(njmg,240)
 240  format('   i',"     am1    " ,"     ar2    ", &
          "     ar3    " ,"     ar4    ","     ar5    ",/)
!      write(njmg,180) isurf,am1(i),ar2(i),ar3(i),ar4(i),ar5(i)
!-----------------------------------------------------------------------
! 6.7  second normazliation for quantities in D(R)
!-----------------------------------------------------------------------
!      write(njmg,250)
 250  format(//,' ***** the second normatlization for D(R) *****', &
         /, ' ***** ar3n is normalized to be the same as V" *****',//)
!      write(njmg,260)
 260  format('   i',"    am1n    " ,"    ar2n    ", &
          "  ar3n(V')  " ,"    ar4n    ","    ar5n    ",/)
!      write(njmg,270) isurf,am1n(i),ar2n(i),ar3n(i),ar4n(i),ar5n(i)
 270  format(i5,5e12.3)
!      close(njmg)
      if (verbose .ge. 4) write(outfile,*)  &
           "====================end output from mercier========="

    if (analytic_peel) then
       write(6,*) 'analytic_peel flag has been set to output the'
       write(6,*) '  analytic peeling criterion lhspeel-rhspeel'
       write(6,*) '  rather than doing a stability calculation'
       write(6,*) '  (lhspeel-rhspeel>0 for stability in del->0 limit)'
       write(6,'(f20.14)') lhspeel(i)-rhspeel(i)
       stop
    endif

    return
      
end subroutine mercy



