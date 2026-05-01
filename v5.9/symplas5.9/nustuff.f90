
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ELITE (Edge Localized Instabilities in Tokamak Experiments)
!
! Authors:  P.B. Snyder, H.R. Wilson, R.L. Miller
!
! Copyright(c) 1999-2014, All rights reserved
!
! Code is under development, and the source code documentation and 
!  user's guide are presently in a preliminary form. The authors 
!  therefore request that users consult with the code authors, 
!  particularly before publishing ELITE results.   We also request 
!  that users report any errors in the documentation or code to the 
!  authors, and share any modifications made to the code or accompanying 
!  idl routines with the authors.  Please request permission from
!  P.B. Snyder or H.R. Wilson before distributing the code, as we wish 
!  to keep track of all users and changes to the code, and maintain a 
!  master version of the code to prevent unnecessary branching.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine nustuff(ixinterp)

! set up matrices for multiple flux surface equilibrium 2/8/00
!  translated to F90 1/31/00
!  including both 1/n terms and those 1/n^2 terms needed for Hermeticity 10/01
!  10/01 called only when n1term=2 (default value), otherwise call
!    oldnustuff below
!-----------------------------------------------------------------------
!     equilibrium quantities are constructed and fourier integrals
!     performed
!     this routine has been largely rewritten to take advantage of
!     Howard's integration by parts
!     changed various signs below to convert to sign convention with
!     psi=minimum on axis. B_p=Grad(phi) x Grad(psi)
!-----------------------------------------------------------------------
      use elite_data, only: ns,mmax,mmin,pprime_adj,ffprime_adj,f_eq, &
           gamscl,rl,dens,ne_eq,neprime_eq,outfile,bpl,uprimel,sinul, &
           d2psi1l,cosul,dpsi1l,ppp_eq,ffpp_eq,circumf,pi,q_eq,qp_eq, &
           qpp_eq,xinterp,m0,nn,q_calc,del,f_in,q_in,nxinterp,zl,n1term, &
           nowindow,nmwinhalf,maxdxinterp,mwindow,mindex,auarr, &
           updownsym,iuarr,bug,b2inv_edge,asurfu, &
           asurfup,isurfu,isurfup,runname,lrunname,ti_eq,tiprime_eq, &
           omegas,omegasn,omegaspi,gamscl2,ion_mass,kinkterms,funcal, &
           verbose
      implicit none
      integer ixinterp
      integer ir,ic,ibench
      real nu(ns),nup(ns),nupp(ns)
      real check1(ns),check2(ns)
      integer i,ik
      integer kk,jj,km
      integer maxrow
      parameter (kk=66)
      real test(ns),nutest(ns),testi2,fval
      real dchi
      real switchon,test1,test2
      real nu1,nu2
      real nuhat,dbpp,hfun,gfun
      real omeg(ns),omegp(ns),omegpp(ns)
      real sig(ns),sigp(ns)        ! parallel current
      real dfac(ns)           ! d/dpsi(R^2 B_p^2/(J B^2))
      real d2fac(ns)           ! (1/nu) d/dpsi(J R^2 B_p^2/B^2)
      real d3fac(ns)         ! d/dpsi (B_p/(R*B^2))^2
      real d4fac(ns),d5fac(ns),d6fac(ns),d7fac(ns),d8fac(ns)
      real dgdpsi(ns),dhdpsi(ns)
      real drbdchi(ns)            ! d/dchi (1/(R^2 B^2))
      real domprbdchi(ns)     ! d/dchi (omega^prime/(R^2 B^2))
      real dompbsqdchi(ns)     ! d/dchi (omega^prime/(B^2))
      real dsigpchi(ns)            ! (d/dchi) sigma-prime
      real dbfacdom(ns)        !(d/domega)(f^3Bp^2/(R^2B^4))
      real dnupdchi(ns)        !d/dchi (nu^prime/nu)
      real d2bsqdpsidchi(ns)   !d^2 B^2/(dchi dpsi)
      real qsurf,qsurfp,qsurfpp
      real bpval,rval,btval,zval
      real jval(ns),bsqval(ns)
      real dbsqp(ns),bpsqval(ns),rsqval(ns)
      integer im,ip,np,nk
      real bsqchi(ns)
      real phase,cosphase(2*(mmax-mmin)+1),sinphase(2*(mmax-mmin)+1)
      real kern(kk)
      real ksin(kk,2*(mmax-mmin)+1),kcos(kk,2*(mmax-mmin)+1)
      real dmnq,dmpnq
      real ucoefr1,ucoefr2,ucoefr3,ucoefr4,ucoefr7,ucoefr9,ucoefr11
      real ucoefi5,ucoefi6,ucoefi8
      real upcoefr3,upcoefr4
      real upcoefi6,upcoefi10
      real uppcoefr4
      integer m,mp,irow,icol,mtop,mlo
      real areal,aimag,apreal,apimag,appreal,appimag
      real arealdd,aimagdd,aprealdd,apimagdd,apprealdd,appimagdd
      real areald,aimagd,apreald,apimagd,appreald,appimagd
      real arealnd,aimagnd,aprealnd,apimagnd,apprealnd,appimagnd
      real arealmk,areal2mk,arealm2k,areal2m,areal2k,arealm
      real arealk,areal0
      real aimagmk,aimag2mk,aimagm2k,aimag2m,aimag2k,aimagm
      real aimagk,aimag0
      real aprealmk,apreal2mk,aprealm2k,apreal2m,apreal2k,aprealm
      real aprealk,apreal0
      real apimagmk,apimag2mk,apimagm2k,apimag2m,apimag2k,apimagm
      real apimagk,apimag0
      real apprealmk,appreal2mk,apprealm2k
      real appimagmk,appimag2mk,appimagm2k
      real irealm,irealk,ireal0
      real iimagm,iimagk,iimag0
      real iprealm,iprealk,ipreal0
      real ipimagm,ipimagk,ipimag0
      real ipprealm,ipprealk,ippreal0
      real ippimagm,ippimagk,ippimag0
      real asurfre,asurfim,apsurfre,apsurfim
      real isurfre,isurfim,ipsurfre,ipsurfim
      real domega
!      real check1,check2,check3
      real dbpdrho,dbtdrho,drhdpsi,drdrho,djbpdrho
      real pprime_loc,ffprime_loc,psi_1,psi_2,psi_3,f_loc,d_loc,dp_loc
      real Rc,psichk1,psichk2,nuppchk,qppfact,qpfact
! TEMP VARS FOR OMEGA* ESTIMATES
      real gamsqws,gamsqwsn


       ibench=0      ! 20/4/01  writes benchmarking file if set to 1 (HRW)
       pprime_loc=pprime_adj(ixinterp)   ! 3/1/00 use adjusted values
       ffprime_loc=ffprime_adj(ixinterp)  ! 3/1/00 use adjusted values
       f_loc=f_eq(ixinterp)
       if (ixinterp.eq.1) then
         gamscl=(rl(1)+rl(ns/2))**2/(f_loc/rl(1)+f_loc/rl(ns/2))**2
         gamscl2=(0.5*(rl(1)+rl(ns/2)))**4/f_loc**2
       end if

! local plasma density, and psi-derivative
       if (dens) then
         d_loc=ne_eq(ixinterp)/ne_eq(1)
         dp_loc=neprime_eq(ixinterp)/ne_eq(1)
       else
         d_loc=1.
         dp_loc=0.
       end if

       if (verbose .ge. 3) write(outfile,*) &
            'surface=',ixinterp,' pprime=',pprime_loc, &
            ' ffprime=',ffprime_loc,' f=',f_loc

!  TEMPORARY!!!!!  ESTIMATES OF OMEGA*
!    simple def of omega_*_i
       if (dens) then
          omegas(ixinterp)=2.9979e10*nn/ne_eq(ixinterp)/4.8032e-10* &
               (pprime_loc/(2.*4.*pi))  ! pprime_loc=4pi pprime
          omegasn(ixinterp)=2.9979e10*nn/ne_eq(ixinterp)/4.8032e-10* &
               (neprime_eq(ixinterp)*1.6022e-12*ti_eq(ixinterp))
          omegaspi(ixinterp)=omegasn(ixinterp)+ &
               2.9979e10*nn/ne_eq(ixinterp)/4.8032e-10* &
               (tiprime_eq(ixinterp)*1.6022e-12*ne_eq(ixinterp))
!  use ion_mass for relative ion mass here
!  4/06 correct to use edge density norm
          gamsqws=pi*ion_mass*1.6726e-24*ne_eq(1)* &
               omegas(ixinterp)**2
          gamsqwsn=pi*ion_mass*1.6726e-24*ne_eq(1)* &
               omegasn(ixinterp)**2
!          write(*,*) 'omegas=',omegas(ixinterp),'gamsqws=',gamsqws, &
!               'ion_mass=',ion_mass
          if (verbose .ge. 2) then
             write(outfile,*) 'omegas=',omegas(ixinterp),'gamsqws=',gamsqws, &
                  'ion_mass=',ion_mass
             write(outfile,*) 'omegasn=',omegasn(ixinterp),' omegaspi=', &
                  omegaspi(ixinterp)
          endif
!          write(*,*) 'omegas=',omegas(ixinterp),'gamsqws=',gamsqws
!          write(*,*) 'omegasn=',omegasn(ixinterp),'gamsqwsn=',gamsqwsn
!          write(outfile,*) 'omegas=',omegas(ixinterp),'gamsqws=',gamsqws
!          write(outfile,*) 'omegasn=',omegasn(ixinterp),'gamsqwsn=',gamsqwsn
!          write(*,*) 'pprime=',pprime_loc/(8.*pi), &
!               neprime_eq(ixinterp)*1.6022e-12*ti_eq(ixinterp), &
!               ne_eq(ixinterp)*1.6022e-12*tiprime_eq(ixinterp)
 !         write(124,*) ixinterp,pprime_loc,omegas,gamsqws
       endif

!-----------------------------------------------------------------------
!     1.0 define nu and first two psi derivatives
!     explicitly extract the psi dependence of nu to facilitate
!     taking psi derivatives
!     nu=f*L/(2*pi*R^2)*(1+nuhat*dpsi)/(bp+dbpp*dpsi)
!     these values have been checked 27/10/97
!     (10/20/00 nupp not being calculated suffieciently accurately...
!     Miller expansion improved to include O(rho^2) corrections for nupp
!     (see Miller et al PoP 5 (1998) 973)
!-----------------------------------------------------------------------
!      write(45,*) dpsi
      do i=1,ns
        psi_1=rl(i)*bpl(i)
        psi_2=0.5*((-rl(i)*uprimel(i)+sinul(i))*bpl(i)- &
          rl(i)**2*pprime_loc-ffprime_loc)
        psi_3=(1./6.)*((-2.*bpl(i)*sinul(i)+4.*psi_2+ffprime_loc)* &
          (sinul(i)/rl(i)-uprimel(i))-rl(i)**2*pprime_loc* &
          (sinul(i)/rl(i)+uprimel(i))-d2psi1l(i)+cosul(i)* &
          dpsi1l(i)/rl(i)-rl(i)*bpl(i)*(rl(i)**2*ppp_eq(ixinterp) &
          +ffpp_eq(ixinterp)))
        Rc=-1./uprimel(i)
        nu(i)= (f_loc*circumf)/(2.*pi*bpl(i)*rl(i)**2)
        nu1=-sinul(i)/rl(i)+uprimel(i)-2.*psi_2/psi_1
! 11/28/00 eliminate the dpsi1l term which cancels when the higher
!   order correction to dl_p is properly considered in the nu'' calculation
        nu2=(sinul(i)/rl(i)+2.*psi_2/psi_1)*(sinul(i)/rl(i)-uprimel(i)) &
            +4.*psi_2**2/psi_1**2-3.*psi_3/psi_1
            !-0.5*(dpsi1l(i)/psi_1)**2
        nup(i)=nu(i)*(ffprime_loc/f_loc**2+nu1/psi_1)
        nupp(i)=nu(i)*((ffpp_eq(ixinterp)-(ffprime_loc/f_loc)**2)/ &
          f_loc**2 +2.*ffprime_loc*nu1/(psi_1*f_loc**2)+ &
          2.*(nu2-nu1*psi_2/psi_1)/psi_1**2)
      end do
!-----------------------------------------------------------------------
!     1.1  can use nu to calculate q's and omega's
!     checked 27/10/97
!-----------------------------------------------------------------------
      omeg(1)=0.
      omegp(1)=0.
      omegpp(1)=0.
      dchi=2.*pi/(ns-1.)
      do i=2,ns
         omeg(i)  =omeg(i-1)  +dchi*(nu(i)  +nu(i-1)  )*0.5
         omegp(i) =omegp(i-1) +dchi*(nup(i) +nup(i-1) )*0.5
         omegpp(i)=omegpp(i-1)+dchi*(nupp(i)+nupp(i-1))*0.5
      end do
      qsurf=omeg(ns)/(2.*pi)
      qsurfp=omegp(ns)/(2.*pi)
      qsurfpp=omegpp(ns)/(2.*pi)
!      write(33,*)' surface=',ixinterp
!      write(33,*) "q/q'/q''=",qsurf,qsurfp,qsurfpp
!      write(33,*)" equilibrium q/q'/q''",q_eq(ixinterp), &
!        qp_eq(ixinterp),qpp_eq(ixinterp)

! temporary fix!!!!  multiply vupp by factor of qpp_eq/qsurfpp to 
!  get equilibrium qpp (which is approximately the derivative of
!  qsurf
      qppfact=qpp_eq(ixinterp)/qsurfpp
      qpfact=qp_eq(ixinterp)/qsurfp
!      do i=1,ns
!         nupp(i)=nupp(i)*qppfact
!         nup(i)=nup(i)*qpfact
!      enddo
!      do i=2,ns
!        omegpp(i)=omegpp(i-1)+dchi*(nupp(i)+nupp(i-1))*0.5
!        omegp(i) =omegp(i-1) +dchi*(nup(i) +nup(i-1) )*0.5
!     enddo
!     qsurfpp=omegpp(ns)/(2.*pi)
!     qsurfp=omegp(ns)/(2.*pi)
!     write(33,*) 'qppfact=',qppfact,'adjusted qsurfpp=',qsurfpp
!     write(33,*) 'qpfact=',qpfact,'adjusted qsurfp=',qsurfp

      xinterp(ixinterp)=m0-nn*qsurf
      if (verbose .ge. 3) write(outfile,*) 'qsurf=',qsurf,' q_calc=',q_calc(ixinterp), &
           'xinterp=',xinterp(ixinterp)
      if(ixinterp==1) then
!! 11/00 making change to make sure grid is consistent for all xinterps
!         del=xinterp(1)
!         xmin=del
!         if (min0(m0+nmwinhalf-int(xmin+0.5),mmax).ne.m1up) then
!            write(*,*) 'fixing roundoff in xmin has changed m1up'
!            stop
!         endif
!         write(*,*) 'adjusted del for roundoff=',del
         xinterp(ixinterp)=del
         if (verbose .ge. 2) write(outfile,*) 'reassigned xinterp=del to avoid roundoff problem'
!         write(75,*) 'm0-nn*qsurf=',m0-nn*qsurf, 'del=',del, &
!              (del-(m0-nn*qsurf))/del
!         write(75,*) 'f_loc=',f_loc,' f_in=',f_in, &
!              (f_in/f_loc)*qsurf, q_in,qsurf

!  TEMP!!!
!   change qsurf as well so it's consistent with xinterp=del
!         write to screen only if q changes significantly
         if ((qsurf-(m0-del)/nn)/qsurf > 0.0001) then
            write(*,*) 'modifying qsurf to get x=del on outer surface'
            write(*,*) 'q_in=',q_in,'qsurfold=',qsurf,'qsurfnew=',(m0-del)/nn
         endif
         if (verbose .ge. 2) then 
            write(outfile,*) &
              'modifying qsurf to get x=del on outer surface'
            write(outfile,*) 'q_in=',q_in,'qsurfold=',qsurf, &
                 'qsurfnew=',(m0-del)/nn
         endif
         qsurf=(m0-del)/nn
!         write(*,*) 'q_in=',q_in,'qsurf=',qsurf

      endif
      if (ixinterp.eq.1 .and. funcal) then
        open(73,file=runname(1:lrunname)//'.omegadat' )
      end if
      if (ixinterp.eq.1) then
        if (ibench.eq.1) then
          open(76,file='bench.dat')
          write(76,*)ns,nxinterp
        end if
      end if
      if (ibench.eq.1) then
        write(76,91)f_loc,ffprime_loc,pprime_loc
        write(76,91)qsurf,qsurfp,qsurfpp
      end if
      do i=1,ns
         omeg(i)=omeg(i)/qsurf
         if (funcal) write(73,91)rl(i),zl(i),bpl(i),omeg(i)
         omegp(i)=-qsurfp*omeg(i)/qsurf+omegp(i)/qsurf
         omegpp(i)=(-qsurfpp*omeg(i)-2*qsurfp*omegp(i)+omegpp(i))/qsurf
         if (ibench.eq.1) then
           write(76,91)rl(i),zl(i),bpl(i)
           write(76,91)omeg(i),omegp(i),omegpp(i)
           write(76,91)nu(i),nup(i),nupp(i)
         end if
      end do
 91   format(3e14.6)
      if (ixinterp.eq.nxinterp .and. funcal) then
        close(73)
        rewind(73)
      end if
!-----------------------------------------------------------------------
!     1.2 chi dependent equilibrium quantities appearing in equations
!     checked rsqval,jval,bsqval,bpsqval,dbsqp 27/10/97
!-----------------------------------------------------------------------

! 2/00 rho is now always zero, as is dpsi.  still need derivatives

      do i=1,ns

         psi_1=rl(i)*bpl(i)
         psi_2=0.5*((sinul(i)+rl(i)*(-uprimel(i)))*bpl(i) &
             -rl(i)**2*pprime_loc-ffprime_loc)
         drhdpsi=1./psi_1
         dbpdrho=-uprimel(i)*bpl(i)-rl(i)*pprime_loc-ffprime_loc/rl(i)
         bpval=bpl(i)
         rval=rl(i)
         rsqval(i)=rval*rval
         jval(i)=circumf/(2.*pi)/bpval
         dbtdrho=(ffprime_loc*bpl(i)/f_loc-f_loc*sinul(i)/rl(i)**2)
         btval=f_loc/rl(i)
         bsqval(i)=bpval**2+btval**2
         bpsqval(i)=bpval*bpval
         dbsqp(i)=2.*(bpl(i)*dbpdrho+f_loc*dbtdrho/rl(i))*drhdpsi
         drdrho=sinul(i)
         djbpdrho=jval(i)*bpl(i)*uprimel(i)
         fval=f_loc
! parallel current and its psi-derivative (10/23/00)
         sig(i)=-f_loc*pprime_loc/bsqval(i)-ffprime_loc/f_loc
         sigp(i)=(f_loc*pprime_loc*dbsqp(i)/bsqval(i)- &
                  ffprime_loc*pprime_loc/f_loc- &
                  f_loc*ppp_eq(ixinterp))/bsqval(i) &
           +ffprime_loc**2/f_loc**3-ffpp_eq(ixinterp)/f_loc
         if (ibench.eq.1) then
           write(76,91)sig(i),sigp(i),dbsqp(i)
         end if
! (d/dpsi) (R^2 B_p^2 /(J B^2))
         dfac(i)=(f_loc*Bpl(i)**2/(nu(i)*bsqval(i)))*(         &
                ffprime_loc/f_loc**2-nup(i)/nu(i)            &
                -dbsqp(i)/bsqval(i)+2.*dbpdrho/(rl(i)*bpl(i)**2))
! (1/nu) (d/dpsi) (R^2 B_p^2 J/B^2)
        d2fac(i)=(rl(i)**4*bpl(i)**2/(f_loc*bsqval(i)))*(            &
                -ffprime_loc/f_loc**2-dbsqp(i)/bsqval(i)+nup(i)/nu(i)   &
                +4.*drdrho/(rl(i)**2*bpl(i))+2.*dbpdrho/(rl(i)*bpl(i)**2))
! (d/dpsi) (B_p^2 /(R^2 B^4))
         d3fac(i)=(2.*bpl(i)/(rl(i)**3*bsqval(i)**2))*(           &
                dbpdrho/bpl(i)-drdrho/rl(i)-psi_1*dbsqp(i)/bsqval(i))
!  nu^2 (d/dpsi) (f^3 Bp^2/(nu^2 R^2 B^4))
         d4fac(i)=(f_loc**3*bpl(i)/(rl(i)**3*bsqval(i)**2))*(      &
                   3.*psi_1*ffprime_loc/f_loc**2                   &
                   +2.*dbpdrho/bpl(i)-2.*drdrho/rl(i)              &
                   -2.*psi_1*dbsqp(i)/bsqval(i)-2.*psi_1*nup(i)/nu(i))
!  nu (d/dpsi) (f p^prime Bp^2/(nu B^4))
         d5fac(i)=bpl(i)**2*(f_loc*ppp_eq(ixinterp)+               &
                             ffprime_loc*pprime_loc/f_loc)/bsqval(i)**2  &
                 +(f_loc*pprime_loc*bpl(i)/(rl(i)*bsqval(i)**2))*(       &
                  2.*dbpdrho/bpl(i)-2.*psi_1*dbsqp(i)/bsqval(i)          &
                  -psi_1*nup(i)/nu(i))
! (d/dpsi) (p^prime R^2 Bp^2/(J B^6))
        d6fac(i)=(bpl(i)**2/(nu(i)*bsqval(i)**3))*(                      &
                   ppp_eq(ixinterp)*f_loc+pprime_loc*ffprime_loc/f_loc)  &
                 +(f_loc*pprime_loc*bpl(i)/(rl(i)*nu(i)*bsqval(i)**3))*( &
                   2.*dbpdrho/bpl(i)-3.*psi_1*dbsqp(i)/bsqval(i)         &
                   -psi_1*nup(i)/nu(i))
! nu (d/dpsi) (f sig /(J B^2))
!  1/24/02 move sig(i) inside to avoid possible divide by zero
        d7fac(i)=(f_loc**2/(rl(i)**2*bsqval(i)))*(                &
                  2.*sig(i)*ffprime_loc/f_loc**2 &
                  +sigp(i)-sig(i)*nup(i)/nu(i)    &
                  -sig(i)*dbsqp(i)/bsqval(i) &
                  -2.*sig(i)*drdrho/(bpl(i)*rl(i)**2))
!        d7fac(i)=(f_loc**2*sig(i)/(rl(i)**2*bsqval(i)))*(                &
!                  2.*ffprime_loc/f_loc**2+sigp(i)/sig(i)-nup(i)/nu(i)    &
!                  -dbsqp(i)/bsqval(i)-2.*drdrho/(bpl(i)*rl(i)**2))
! (d/dpsi) (f^2 p^prime/(J B^6))
        d8fac(i)=(3.*f_loc*ffprime_loc*pprime_loc+                        &
                 ppp_eq(ixinterp)*f_loc**3)/(nu(i)*rl(i)**2*bsqval(i)**3) &
                -(f_loc**3*pprime_loc/(nu(i)*rl(i)**2*bsqval(i)**3))*(   &
                 nup(i)/nu(i)+3.*dbsqp(i)/bsqval(i)                      &
                 +2.*drdrho/(rl(i)**2*bpl(i)))
        dgdpsi(i)=(pprime_loc/bsqval(i)+                                  &
                  (nup(i)/nu(i))*(f_loc**2/(rl(i)**2*bsqval(i))))*        &
                  (rl(i)**4*bpl(i)**2*nu(i)/(f_loc*bsqval(i)))*(          &
                    4.*drdrho/(rl(i)**2*bpl(i))                           &
                   +2.*dbpdrho/(rl(i)*bpl(i)**2)                          &
                   +nup(i)/nu(i)-ffprime_loc/f_loc**2-dbsqp(i)/bsqval(i)) &
                +(rl(i)**4*bpl(i)**2*nu(i)/(f_loc*bsqval(i)))*(           &
                    ppp_eq(ixinterp)/bsqval(i)                            &
                   -pprime_loc*dbsqp(i)/bsqval(i)**2                      &
                   +f_loc**2*(nupp(i)/nu(i)-                                 &
                              (nup(i)/nu(i))**2)/(rl(i)**2*bsqval(i))     &
                   +(nup(i)*f_loc**2/(nu(i)*rl(i)**2*bsqval(i)))*(        &
                     2.*ffprime_loc/f_loc**2-dbsqp(i)/bsqval(i)           &
                     -2.*drdrho/(rl(i)**2*bpl(i))))

        dhdpsi(i)=(f_loc*rl(i)**2*bpl(i)**2/bsqval(i)**2)*(               &
                    ffprime_loc/f_loc**2                                  &
                   +2.*drdrho/(rl(i)**2*bpl(i))                           &
                   +2.*dbpdrho/(rl(i)*bpl(i)**2)-2.*dbsqp(i)/bsqval(i))
! calculate zval for writing to diagnostic output file surfaces.out
!    pbs 11/99
         zval=zl(i)
!         write(47,*) rval,zval,bpval
      end do
!-----------------------------------------------------------------------
!     1.3 need chi deriv of B^2
!         and 1/(R^2 B^2)  and omegp/(R^2 B^2)  (10/23/00)
!         and (nup/nu) and dbsqdpsi                         (9/13/01)
!-----------------------------------------------------------------------
      do i=1,ns
         ip=i+1
         im=i-1
         if(i.eq.ns) ip=2
         if(i.eq.1) im=ns-1
         bsqchi(i)=(bsqval(ip)-bsqval(im))/(2.*dchi)
         d2bsqdpsidchi(i)=(dbsqp(ip)-dbsqp(im))/(2.*dchi)
         drbdchi(i)=(1./(rsqval(ip)*bsqval(ip))- &
             1./(rsqval(im)*bsqval(im)))/(2.*dchi)
         dompbsqdchi(i)=(omegp(ip)/bsqval(ip)-omegp(im)/bsqval(im))/ &
             (2.*dchi)
!  debugging 3/20/14, may be sign error in 2nd part of this term
!         domprbdchi(i)=-((nu(i)*qsurfp/qsurf-nup(i))/ &
!             (qsurf*rsqval(i)*bsqval(i))+omegp(i)*drbdchi(i))
!         domprbdchi(i)=-((nu(i)*qsurfp/qsurf-nup(i)))/ &
!              (qsurf*rsqval(i)*bsqval(i))+omegp(i)*drbdchi(i)
         domprbdchi(i)=(omegp(ip)/(bsqval(ip)*rsqval(ip)) &
              -omegp(im)/(bsqval(im)*rsqval(im)))/ &
             (2.*dchi)


         dsigpchi(i)=(sigp(ip)-sigp(im))/(2.*dchi)
         dbfacdom(i)=(f_loc**3*qsurf/nu(i))*                         &
                     (bpl(ip)**2/(rsqval(ip)*bsqval(ip)**2)-  &
                  bpl(im)**2/(rsqval(im)*bsqval(im)**2))/(2.*dchi)
         dnupdchi(i)=(nup(ip)/nu(ip)-nup(im)/nu(im))/(2.*dchi)
      end do
!-----------------------------------------------------------------------
!     2.0 chi integrals--zero the integral summations
!     kcos(ik,ip) is Integral[T(ik)*Cos[(ip-1+mmin-mmax)*omeg]]
!     So ik is index of T matrix appearing in notes
!     and ip is mp-m +(1+mmax-mmin)
!     mmax i maximum m value
!     mmin is minimum m value
!     ipmax is 2*(mmax-mmin)+1
!-----------------------------------------------------------------------
      np=2*(mmax-mmin)+1
      nk=kk
      do ip=1,np
         do ik=1,nk
            kcos(ik,ip)=0.
            ksin(ik,ip)=0.
         end do
      end do
!-----------------------------------------------------------------------
!     2.1 perform the integrations in omega
!-----------------------------------------------------------------------
      do i=1,ns-1
         domega=nu(i)*dchi/qsurf
         do ip=1,np
            phase=(ip-1+mmin-mmax)*omeg(i) 
            cosphase(ip)=cos(phase)
            sinphase(ip)=sin(phase)
         end do
         kern(1)=f_loc/(bpsqval(i)*rsqval(i)**2)
         kern(2)=(qsurf*bpsqval(i)*nupp(i)*rsqval(i))/(bsqval(i)*jval(i))
         kern(3)=(bpsqval(i)*nup(i)*rsqval(i))/(bsqval(i)*jval(i))
         kern(4)=(f_loc*bpsqval(i))/bsqval(i)
         kern(5)=(bpsqval(i)*nup(i)*omegp(i)*rsqval(i))/(bsqval(i)*jval(i))
         kern(6)=(f_loc*bpsqval(i)*omegp(i))/bsqval(i)
         kern(7)=(f_loc*bpsqval(i)*omegp(i)**2)/bsqval(i)
         kern(8)=(f_loc*bpsqval(i)*omegpp(i))/bsqval(i)
         kern(9)=((pprime_loc + dbsqp(i)*0.5)*rsqval(i))/bsqval(i)
         kern(10)=(bsqchi(i)*rsqval(i))/(bsqval(i)**2*jval(i))
         kern(11)=(bsqchi(i)*omegp(i)*rsqval(i))/(bsqval(i)**2*jval(i))
         kern(12)=1.
         kern(13)=1./bsqval(i)
         kern(14)=nup(i)*dfac(i)
         kern(15)=nu(i)*omegp(i)*dfac(i)
         kern(16)=nu(i)*dfac(i)
         kern(17)=sig(i)*pprime_loc/bsqval(i)
         kern(18)=f_loc*pprime_loc*bpsqval(i)/bsqval(i)**2
         kern(19)=f_loc**2*sig(i)*omegp(i)/(rsqval(i)*bsqval(i))
         kern(20)=f_loc**2*sig(i)/(rsqval(i)*bsqval(i))
         kern(21)=f_loc**3*bpsqval(i)*omegp(i)/(rsqval(i)*bsqval(i)**2)
         kern(22)=f_loc**3*omegp(i)**2*bpsqval(i)/(rsqval(i)*bsqval(i)**2)
         kern(23)=f_loc**3*bpsqval(i)/(rsqval(i)*bsqval(i)**2)
         kern(24)=pprime_loc*rsqval(i)*bpsqval(i)*bsqchi(i)/ &
                  (jval(i)*bsqval(i)**3)
         kern(25)=pprime_loc*rsqval(i)*bpsqval(i)*omegp(i)*bsqchi(i)/ &
                  (jval(i)*bsqval(i)**3)
         kern(26)=f_loc**2*sig(i)*drbdchi(i)/nu(i)
         kern(27)=f_loc**2*sig(i)*domprbdchi(i)/nu(i)
         kern(28)=f_loc**3*bpsqval(i)*nup(i)/(rsqval(i)*nu(i) &
                  *bsqval(i)**2)
         kern(29)=f_loc**3*bpsqval(i)*omegp(i)*nup(i)/ &
                  (rsqval(i)*bsqval(i)**2*nu(i))
         kern(30)=dsigpchi(i)/nu(i)
         kern(31)=sigp(i)
         kern(32)=d_loc/bpsqval(i)
         kern(33)=d_loc*rsqval(i)**2*bpsqval(i)/bsqval(i)
         kern(34)=d_loc*rsqval(i)**2*bpsqval(i)*omegp(i)/bsqval(i)
         kern(35)=dp_loc*rsqval(i)**2*bpsqval(i)/(f_loc*bsqval(i)) &
                  +d_loc*d2fac(i)
         kern(36)=dp_loc*omegp(i)*rsqval(i)**2*bpsqval(i)/ &
                  (f_loc*bsqval(i))+d_loc*omegp(i)*d2fac(i)
         kern(37)=d_loc*rsqval(i)**2*bpsqval(i)*omegpp(i)/bsqval(i)
         kern(38)=d_loc*rsqval(i)**2*bpsqval(i)*omegp(i)**2/bsqval(i)
         hfun=f_loc*rsqval(i)*bpsqval(i)/bsqval(i)**2
         gfun=(nu(i)*rsqval(i)**2*bpsqval(i)/(f_loc*bsqval(i)))* &
                  (pprime_loc/bsqval(i)+ &
                  nup(i)*f_loc*f_loc/(nu(i)*rsqval(i)*bsqval(i)))
         kern(39)=d_loc*hfun*omegp(i)**2
         kern(40)=d_loc*hfun*omegp(i)
         kern(41)=d_loc*hfun
         kern(42)=d_loc*hfun*nupp(i)/nu(i)
         kern(43)=d_loc*gfun/nu(i)
         kern(44)=d_loc*hfun*nup(i)/nu(i)
         kern(45)=f_loc**3*bpsqval(i)*(nup(i)/(nu(i)*bsqval(i)))**2  &
                  /rsqval(i)
         kern(46)=f_loc**3*bpsqval(i)*omegpp(i)/(rsqval(i)*bsqval(i)**2)
         kern(47)=f_loc**3*bpsqval(i)*nupp(i)/                       &
                  (nu(i)*rsqval(i)*bsqval(i)**2)
         kern(48)=dbfacdom(i)*omegp(i)*nup(i)/nu(i)
         kern(49)=dbfacdom(i)*nupp(i)/nu(i)
         kern(50)=dbfacdom(i)*nup(i)/nu(i)
         kern(51)=d4fac(i)*omegp(i)
         kern(52)=d4fac(i)*nup(i)/nu(i)
         kern(53)=d4fac(i)
         kern(54)=d4fac(i)*dnupdchi(i)/nu(i)
         kern(55)=d5fac(i)
         kern(56)=pprime_loc*bpl(i)**2*f_loc*nup(i)*bsqchi(i)/        &
                  (bsqval(i)**3*nu(i)**2)
         kern(57)=d6fac(i)*bsqchi(i)                                 &
                 +(pprime_loc*rl(i)**2*bpl(i)**2/(jval(i)*bsqval(i)**3))  &
                     *d2bsqdpsidchi(i)
         kern(58)=f_loc**2*pprime_loc*nup(i)*bsqchi(i)/                   &
                  (jval(i)*nu(i)*bsqval(i)**3)
         kern(59)=d7fac(i)
         kern(60)=d8fac(i)*bsqchi(i)+                                     &
                  (f_loc**2*pprime_loc/(jval(i)*bsqval(i)**3))*           &
                   d2bsqdpsidchi(i)
         kern(61)=d_loc*hfun*omegpp(i)
         kern(62)=d_loc*hfun*omegp(i)*nup(i)/nu(i)
         kern(63)=(d_loc*dgdpsi(i)+dp_loc*gfun)/nu(i)
         kern(64)=omegp(i)*(d_loc*dhdpsi(i)+dp_loc*hfun)
         kern(65)=(nup(i)/nu(i))*(d_loc*dhdpsi(i)+dp_loc*hfun)
         kern(66)=d_loc*dhdpsi(i)+dp_loc*hfun
         do ik=1,nk
            do ip=1,np
               kcos(ik,ip)=kcos(ik,ip)+domega*cosphase(ip)*kern(ik)
               ksin(ik,ip)=ksin(ik,ip)+domega*sinphase(ip)*kern(ik)
            end do
         end do
      end do
!      write(86,*)' ixinterp=',ixinterp
!-----------------------------------------------------------------------
!     3.0 calculate full x-dependent matrix equation for Um and b.c
!
!     au(irow,icol)*um(irow,icol)+
!     aup(irow,icol)*d/dx[um(irow,icol)]+
!     aupp(irow,icol)*d^2/dx^2[um(irow,icol)]=0 
!
!     asurfu(irow,icol)*um(irow,icol)+
!     asurfup(irow,icol)*d/dx[um(irow,icol)]=
!     gam*agamma(irow,icol)*um(irow,icol)
!     note that um={Re[um[1]],Im[um[1]],Re[um[2],Im[um[2]],...}
!-----------------------------------------------------------------------
! PS 7/01 arrays are now lwindow x lwindow in size, redo indexing
!   to go from mtop to mlo rather then mmax to mmin, allow padding
!   on top and bottom (if not against mmax or mmin) to allow proper
!   offset in matgen
      mtop=mmax
      mlo=mmin
      if(.not.nowindow) then
         ! note that mtop is maxdxinterp above usual mup.  need to leave room
         !  for offset in matgen
         mtop=min0(m0+nmwinhalf-int(xinterp(ixinterp)+0.5)+maxdxinterp,mmax)
         if (mtop < (mmin+mwindow/mindex-1+maxdxinterp) ) &
              mtop=mmin+mwindow/mindex-1+maxdxinterp
         mtop=min0(mtop,mmax)
         mlo=mtop-mwindow/mindex+1-2*maxdxinterp  ! reduce by two to allow room
         if (mlo < mmin) mlo=mmin
      endif
!      do m=mmin,mmax
      do m=mlo,mtop
         dmnq=m-nn*qsurf
         do mp=mlo,mtop
            ip=mp-m+1+mmax-mmin
            dmpnq=mp-nn*qsurf
! label k denotes mp
! additional higher order terms added by HRW; 28 Sept 01
! terms to be multiplied by dmnq*dmpnq*u
            arealmk=kcos(1,ip)/qsurf-kcos(2,ip)/(nn*qsurf)**2 &
              -2.*m*ksin(5,ip)/(nn**2*qsurf)+m**2*kcos(7,ip)/(nn**2*qsurf) &
              -m*ksin(8,ip)/(nn**2*qsurf)                         &
              -kcos(14,ip)/(nn**2*qsurf)                          &
              -m*ksin(15,ip)/(nn**2*qsurf)                        &
              -4.*m*m*qsurfp*ksin(21,ip)/(nn*qsurf)**3            &
              -2.*(qsurfpp-2.*qsurfp**2/qsurf)*m*kcos(23,ip)/(nn*qsurf)**3  &
              -4.*m*qsurfp*kcos(28,ip)/(nn*qsurf)**3              &
              -2.*kcos(45,ip)/(nn**2*qsurf)                       &
              +m*kcos(47,ip)/(nn**3*qsurf**2)                    &
              -2.*m*kcos(48,ip)/(nn**3*qsurf**2)                  &
              +ksin(49,ip)/(nn**3*qsurf**2)                       &
              +mp*kcos(52,ip)/(nn**3*qsurf**2)                    &
              -2.*m*qsurfp*kcos(53,ip)/(nn**3*qsurf**3)          &
              -ksin(54,ip)/(nn**3*qsurf)                          &
              -kcos(55,ip)/(nn**2*qsurf)                          &
              -kcos(59,ip)/(nn**2*qsurf)
            aimagmk=ksin(1,ip)/qsurf-ksin(2,ip)/(nn*qsurf)**2 &
              +2.*m*kcos(5,ip)/(nn**2*qsurf)+m**2*ksin(7,ip)/(nn**2*qsurf) &
              +m*kcos(8,ip)/(nn**2*qsurf)                         &
              -ksin(14,ip)/(nn**2*qsurf)                          &
              +m*kcos(15,ip)/(nn**2*qsurf)                        &
              +4.*m*m*qsurfp*kcos(21,ip)/(nn*qsurf)**3            &
              -2.*(qsurfpp-2.*qsurfp**2/qsurf)*m*ksin(23,ip)/(nn*qsurf)**3 & 
              -4.*m*qsurfp*ksin(28,ip)/(nn*qsurf)**3              &
              -2.*ksin(45,ip)/(nn**2*qsurf)                       &
              +m*ksin(47,ip)/(nn**3*qsurf**2)                    &
              -2.*m*ksin(48,ip)/(nn**3*qsurf**2)                  &
              -kcos(49,ip)/(nn**3*qsurf**2)                       &
              +mp*ksin(52,ip)/(nn**3*qsurf**2)                    &
              -2.*m*qsurfp*ksin(53,ip)/(nn**3*qsurf**3)           &
              +kcos(54,ip)/(nn**3*qsurf)                          &
              -ksin(55,ip)/(nn**2*qsurf)                          &
              -ksin(59,ip)/(nn**2*qsurf)
! terms to be multiplied by dmnq**2*dmpnq*u
            areal2mk=-m*m*kcos(22,ip)/(nn**3*qsurf**2)            &
               +2.*m*ksin(29,ip)/(nn**3*qsurf**2)                 &
               +m*ksin(46,ip)/(nn**3*qsurf**2)                    &
               +m*ksin(51,ip)/(nn**3*qsurf**2)                    &
               +2.*kcos(52,ip)/(nn**3*qsurf**2)
            aimag2mk=-m*m*ksin(22,ip)/(nn**3*qsurf**2)            &
               -2.*m*kcos(29,ip)/(nn**3*qsurf**2)                 &
               -m*kcos(46,ip)/(nn**3*qsurf**2)                    &
               -m*kcos(51,ip)/(nn**3*qsurf**2)                    &
               +2.*ksin(52,ip)/(nn**3*qsurf**2)
! terms to be multiplied by dmnq*dmpnq**2*u
            arealm2k=-m*m*kcos(22,ip)/(nn**3*qsurf**2)            &
               +4.*m*ksin(29,ip)/(nn**3*qsurf**2)                 &
               +m*ksin(46,ip)/(nn**3*qsurf**2)                    &
               +2.*kcos(47,ip)/(nn**3*qsurf**2)                   &
               +m*ksin(51,ip)/(nn**3*qsurf**2)
            aimagm2k=-m*m*ksin(22,ip)/(nn**3*qsurf**2)            &
               -4.*m*kcos(29,ip)/(nn**3*qsurf**2)                 &
               -m*kcos(46,ip)/(nn**3*qsurf**2)                    &
               +2.*ksin(47,ip)/(nn**3*qsurf**2)                   &
               -m*kcos(51,ip)/(nn**3*qsurf**2)
! terms to be multiplied by dmnq**2*u
            areal2m=m*ksin(19,ip)/(nn**2*qsurf)
            aimag2m=-m*kcos(19,ip)/(nn**2*qsurf)
! terms to be multiplied by dmpnq**2*u
            areal2k=-m*ksin(19,ip)/(nn**2*qsurf)                  &
              -2.*m*m*qsurfp*ksin(21,ip)/(nn**3*qsurf**3)         &
              -(qsurfpp-2.*qsurfp**2/qsurf)*m*kcos(23,ip)/(nn*qsurf)**3  &
              -4.*m*qsurfp*kcos(28,ip)/(nn*qsurf)**3              &
              -m*qsurfp*kcos(53,ip)/(nn*qsurf)**3
            aimag2k=m*kcos(19,ip)/(nn**2*qsurf)                   &
              +2.*m*m*qsurfp*kcos(21,ip)/(nn*qsurf)**3            &
              -(qsurfpp-2.*qsurfp**2/qsurf)*m*ksin(23,ip)/(nn*qsurf)**3  &
              -4.*m*qsurfp*ksin(28,ip)/(nn*qsurf)**3              &
              -m*qsurfp*ksin(53,ip)/(nn*qsurf)**3
! terms to be mulitplied by dmnq*u
            arealm=-kcos(17,ip)/nn                                &
              -m*qsurfp*kcos(20,ip)/(nn*qsurf)**2                 &
              +m*kcos(25,ip)/nn**2                                &
              +m*kcos(27,ip)/nn**2                                &
              +kcos(31,ip)/nn                                     &
              -ksin(56,ip)/nn**2                                  &
              +ksin(58,ip)/nn**2
            if (.not. kinkterms) then
               arealm=arealm-kcos(31,ip)/nn
            endif
            aimagm=-ksin(17,ip)/nn                                &
              -m*qsurfp*ksin(20,ip)/(nn*qsurf)**2                 &
              +m*ksin(25,ip)/nn**2                                &
              +m*ksin(27,ip)/nn**2                                &
              +ksin(31,ip)/nn                                     &
              +kcos(56,ip)/nn**2                                  &
              -kcos(58,ip)/nn**2
            if (.not. kinkterms) then
               aimagm=aimagm-ksin(31,ip)/nn
            endif
! terms to be multiplied by dmpnq*u
            arealk=2.*m*qsurfp*kcos(3,ip)/(nn*qsurf)**2           &
              +m*qsurfpp*kcos(4,ip)/(nn*qsurf)**2                 &
              -2.*m*qsurfp**2/(nn**2*qsurf**3)*kcos(4,ip)         &
              +2.*m**2*qsurfp*ksin(6,ip)/(nn*qsurf)**2            &
              +m*qsurfp*kcos(16,ip)/(nn*qsurf)**2                 &
              -kcos(17,ip)/nn                                     &
              +m*qsurfp*kcos(20,ip)/(nn*qsurf)**2                 &
              +2.*m*m*qsurfp**2*kcos(23,ip)/(nn**3*qsurf**4)     &
              +m*kcos(25,ip)/nn**2                                &
              +m*kcos(27,ip)/nn**2                                &
              -2.*m*qsurfp*ksin(50,ip)/(nn*qsurf)**3              &
              -ksin(57,ip)/nn**2                                  &
              +ksin(60,ip)/nn**2
            aimagk=2.*m*qsurfp*ksin(3,ip)/(nn*qsurf)**2           &
              -2.*m*qsurfp**2*ksin(4,ip)/(nn**2*qsurf**3)         &
              +m*qsurfpp*ksin(4,ip)/(nn*qsurf)**2                 &
              -2.*m**2*qsurfp*kcos(6,ip)/(nn*qsurf)**2            &
              +m*qsurfp*ksin(16,ip)/(nn*qsurf)**2                 &
              -ksin(17,ip)/nn                                     &
              +m*qsurfp*ksin(20,ip)/(nn*qsurf)**2                 &
              +2.*m*m*qsurfp**2*ksin(23,ip)/(nn**3*qsurf**4)      &
              +m*ksin(25,ip)/nn**2                                &
              +m*ksin(27,ip)/nn**2                                &
              +2.*m*qsurfp*kcos(50,ip)/(nn*qsurf)**3              &
              +kcos(57,ip)/nn**2                                  &
              -kcos(60,ip)/nn**2
! terms to be multiplied by 1*u
            areal0=-2.*qsurf*pprime_loc*kcos(9,ip)/f_loc          &
              +pprime_loc*qsurf*m*kcos(11,ip)/nn                  &
              +m*qsurfp*ksin(24,ip)/(nn**2*qsurf)                 &
              +m*qsurfp*ksin(26,ip)/(nn**2*qsurf)                 &
              -qsurf*ksin(30,ip)/nn
            if (.not. kinkterms) then
               areal0=areal0+0.5*qsurf*ksin(30,ip)/nn
            endif
            aimag0=-2.*qsurf*pprime_loc*ksin(9,ip)/f_loc          &
              +pprime_loc*qsurf*m*ksin(11,ip)/nn                  &
              -m*qsurfp*kcos(24,ip)/(nn**2*qsurf)                 &
              -m*qsurfp*kcos(26,ip)/(nn**2*qsurf)                 &
              +qsurf*kcos(30,ip)/nn  
            if (.not. kinkterms) then
               aimag0=aimag0-0.5*qsurf*kcos(30,ip)/nn
            endif
! terms to be multiplied by dmnq*dmpnq*u'
            aprealmk=-2.*kcos(3,ip)/(nn**2*qsurf)                 &
              -2.*m*ksin(6,ip)/(nn**2*qsurf)                      &
              -kcos(16,ip)/(nn**2*qsurf)                          &
              -4.*m*qsurfp*kcos(23,ip)/(nn*qsurf)**3              &
              +2.*ksin(50,ip)/(nn**3*qsurf**2)
            apimagmk=-2.*ksin(3,ip)/(nn**2*qsurf)                 &
              +2.*m*kcos(6,ip)/(nn**2*qsurf)                      &
              -ksin(16,ip)/(nn**2*qsurf)                          &
              -4.*m*qsurfp*ksin(23,ip)/(nn*qsurf)**3              &
              -2.*kcos(50,ip)/(nn**3*qsurf**2)
! terms to be multiplied by dmnq**2*dmpnq*u'
            apreal2mk=2.*m*ksin(21,ip)/(nn**3*qsurf**2)           &
              +2.*kcos(28,ip)/(nn**3*qsurf**2)                    &
              +kcos(53,ip)/(nn**3*qsurf**2)
            apimag2mk=-2.*m*kcos(21,ip)/(nn**3*qsurf**2)          &
              +2.*ksin(28,ip)/(nn**3*qsurf**2)                    &
              +ksin(53,ip)/(nn**3*qsurf**2) 
! terms to be multiplied by dmnq*dmpnq**2*u'
            aprealm2k=2.*m*ksin(21,ip)/(nn**3*qsurf**2)           &
              +4.*kcos(28,ip)/(nn**3*qsurf**2)                    &
              +kcos(53,ip)/(nn**3*qsurf**2)
            apimagm2k=-2.*m*kcos(21,ip)/(nn**3*qsurf**2)          &
              +4.*ksin(28,ip)/(nn**3*qsurf**2)                    &
              +ksin(53,ip)/(nn**3*qsurf**2)
! terms to be multiplied by dmnq**2*u'
            apreal2m=kcos(20,ip)/(nn**2*qsurf)
            apimag2m=ksin(20,ip)/(nn**2*qsurf)
! terms to be multiplied by dmpnq**2*u'
            apreal2k=-2.*m*qsurfp*kcos(23,ip)/(nn*qsurf)**3       &
                     -kcos(20,ip)/(nn**2*qsurf)
            apimag2k=-2.*m*qsurfp*ksin(23,ip)/(nn*qsurf)**3       &
                     -ksin(20,ip)/(nn**2*qsurf)
! terms to be mulitplied by dmnq*u'
            aprealm=-ksin(24,ip)/nn**2-ksin(26,ip)/nn**2
            apimagm=kcos(24,ip)/nn**2+kcos(26,ip)/nn**2
! terms to be multiplied by dmpnq*u'
            aprealk=2.*m*qsurfp*kcos(4,ip)/(nn*qsurf)**2          &
              -ksin(24,ip)/nn**2                                  &
              -ksin(26,ip)/nn**2
            apimagk=2.*m*qsurfp*ksin(4,ip)/(nn*qsurf)**2          &
              +kcos(24,ip)/nn**2                                  &
              +kcos(26,ip)/nn**2
! terms to be multiplied by 1*u'
            apreal0=-qsurf*pprime_loc*ksin(10,ip)/nn
            apimag0=qsurf*pprime_loc*kcos(10,ip)/nn
! terms to be multiplied by dmnq*dmpnq*u"
            apprealmk=-kcos(4,ip)/(nn**2*qsurf)
            appimagmk=-ksin(4,ip)/(nn**2*qsurf)
! terms to be multiplied by dmnq**2*dmpnq*u"
            appreal2mk=kcos(23,ip)/(nn**3*qsurf**2)
            appimag2mk=ksin(23,ip)/(nn**3*qsurf**2)
! terms to be multiplied by dmnq*dmpnq**2*u"
            apprealm2k=kcos(23,ip)/(nn**3*qsurf**2)
            appimagm2k=ksin(23,ip)/(nn**3*qsurf**2)
!
!       write(97,*)' ixinterp=',ixinterp,' m=',m,' k=',mp
!       write(97,*)' nn=',nn,' qsurf=',qsurf
!       if (apprealmk.eq.0.) then
!       write(97,*)' term1=',apprealmk,' term2=',apprealm2k
!      else
!       write(97,*)' term1=',apprealmk,' ratio=', &
!           apprealm2k*(mp+m-2.*nn*qsurf) &
!                    /apprealmk
!       end if
!     write(86,*)' Ar:',arealmk,areal2mk,arealm2k
!     write(86,*)'    ',areal2m,areal2k,arealm
!     write(86,*)'    ',arealk,areal0
!     write(86,*)' Ai:',aimagmk,aimag2mk,aimagm2k
!     write(86,*)'    ',aimag2m,aimag2k,aimagm
!     write(86,*)'    ',aimagk,aimag0
!     write(86,*)' Apr:',aprealmk,apreal2mk,aprealm2k
!     write(86,*)'    ',apreal2m,apreal2k,aprealm
!     write(86,*)'    ',aprealk,apreal0
!     write(86,*)' Api:',apimagmk,apimag2mk,apimagm2k
!     write(86,*)'    ',apimag2m,apimag2k,apimagm
!     write(86,*)'    ',apimagk,apimag0
!     write(86,*)' Appr:',apprealmk,appreal2mk,apprealm2k
!     write(86,*)' Appi:',appimagmk,appimag2mk,appimagm2k

 !-----------------------------------------------------------------------
!            irow=mindex*(mp-mmin)+1 !mindex=1 for updownsym
!            icol=mindex*(m -mmin)+1 !mindex=2 for.not.updownsym
!     used to use above statements where irow=1 was mmin but 
!     now irow=1 is mmax
!-----------------------------------------------------------------------
!            irow=mindex*(mmax-mp)+1 !mindex=1 for updownsym
            irow=mindex*(mtop-mp)+1 !mindex=1 for updownsym
!            icol=mindex*(mmax -m)+1 !mindex=2 for.not.updownsym
            icol=mindex*(mtop -m)+1 !mindex=2 for.not.updownsym
            auarr(1,irow  ,icol  ,ixinterp)=arealmk
            auarr(2,irow  ,icol  ,ixinterp)=areal2mk
            auarr(3,irow  ,icol  ,ixinterp)=arealm2k
            auarr(4,irow  ,icol  ,ixinterp)=areal2m
            auarr(5,irow  ,icol  ,ixinterp)=areal2k
            auarr(6,irow  ,icol  ,ixinterp)=arealm
            auarr(7,irow  ,icol  ,ixinterp)=arealk
            auarr(8,irow  ,icol  ,ixinterp)=areal0
            auarr(9,irow  ,icol  ,ixinterp)= &
                 -nn*(qsurfp*aprealmk+qsurfpp*apprealmk)
            auarr(10,irow  ,icol  ,ixinterp)= &
                 -nn*(qsurfp*apreal2mk+qsurfpp*appreal2mk)
            auarr(11,irow  ,icol  ,ixinterp)= &
                 -nn*(qsurfp*aprealm2k+qsurfpp*apprealm2k)
            auarr(12,irow  ,icol  ,ixinterp)=-nn*qsurfp*apreal2m
            auarr(13,irow  ,icol  ,ixinterp)=-nn*qsurfp*apreal2k
            auarr(14,irow  ,icol  ,ixinterp)=-nn*qsurfp*aprealm
            auarr(15,irow  ,icol  ,ixinterp)=-nn*qsurfp*aprealk
            auarr(16,irow  ,icol  ,ixinterp)=-nn*qsurfp*apreal0
            auarr(17,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*apprealmk
            auarr(18,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*appreal2mk
            auarr(19,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*apprealm2k
! I want to test Hermitian relations for dW minimization
!            write(45,'(2i5,3e15.7)') m,mp &
!                ,(arealdd*  dmnq+areald  )*dmpnq+arealnd &
!                ,(aprealdd* dmnq+apreald )*dmpnq+aprealnd &
!                ,(apprealdd*dmnq+appreald)*dmpnq+apprealnd
            if(.not.updownsym) then
               auarr(1,irow+1,icol+1,ixinterp)=arealmk
               auarr(1,irow+1,icol  ,ixinterp)=aimagmk
               auarr(1,irow  ,icol+1,ixinterp)=-aimagmk
               auarr(2,irow+1,icol+1,ixinterp)=areal2mk
               auarr(2,irow+1,icol  ,ixinterp)=aimag2mk
               auarr(2,irow  ,icol+1,ixinterp)=-aimag2mk
               auarr(3,irow+1,icol+1,ixinterp)=arealm2k
               auarr(3,irow+1,icol  ,ixinterp)=aimagm2k
               auarr(3,irow  ,icol+1,ixinterp)=-aimagm2k
               auarr(4,irow+1,icol+1,ixinterp)=areal2m
               auarr(4,irow+1,icol  ,ixinterp)=aimag2m
               auarr(4,irow  ,icol+1,ixinterp)=-aimag2m
               auarr(5,irow+1,icol+1,ixinterp)=areal2k
               auarr(5,irow+1,icol  ,ixinterp)=aimag2k
               auarr(5,irow  ,icol+1,ixinterp)=-aimag2k
               auarr(6,irow+1,icol+1,ixinterp)=arealm
               auarr(6,irow+1,icol  ,ixinterp)=aimagm
               auarr(6,irow  ,icol+1,ixinterp)=-aimagm
               auarr(7,irow+1,icol+1,ixinterp)=arealk
               auarr(7,irow+1,icol  ,ixinterp)=aimagk
               auarr(7,irow  ,icol+1,ixinterp)=-aimagk
               auarr(8,irow+1,icol+1,ixinterp)=areal0
               auarr(8,irow+1,icol  ,ixinterp)=aimag0
               auarr(8,irow  ,icol+1,ixinterp)=-aimag0
               auarr(9,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*aprealmk+qsurfpp*apprealmk)
               auarr(9,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*apimagmk+qsurfpp*appimagmk)
               auarr(9,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*apimagmk+qsurfpp*appimagmk)
               auarr(10,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*apreal2mk+qsurfpp*appreal2mk)
               auarr(10,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*apimag2mk+qsurfpp*appimag2mk)
               auarr(10,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*apimag2mk+qsurfpp*appimag2mk)
               auarr(11,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*aprealm2k+qsurfpp*apprealm2k)
               auarr(11,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*apimagm2k+qsurfpp*appimagm2k)
               auarr(11,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*apimagm2k+qsurfpp*appimagm2k)
               auarr(12,irow+1,icol+1,ixinterp)=-nn*qsurfp*apreal2m
               auarr(12,irow+1,icol  ,ixinterp)=-nn*qsurfp*apimag2m
               auarr(12,irow  ,icol+1,ixinterp)= nn*qsurfp*apimag2m
               auarr(13,irow+1,icol+1,ixinterp)=-nn*qsurfp*apreal2k
               auarr(13,irow+1,icol  ,ixinterp)=-nn*qsurfp*apimag2k
               auarr(13,irow  ,icol+1,ixinterp)= nn*qsurfp*apimag2k
               auarr(14,irow+1,icol+1,ixinterp)=-nn*qsurfp*aprealm
               auarr(14,irow+1,icol  ,ixinterp)=-nn*qsurfp*apimagm
               auarr(14,irow  ,icol+1,ixinterp)= nn*qsurfp*apimagm
               auarr(15,irow+1,icol+1,ixinterp)=-nn*qsurfp*aprealk
               auarr(15,irow+1,icol  ,ixinterp)=-nn*qsurfp*apimagk
               auarr(15,irow  ,icol+1,ixinterp)= nn*qsurfp*apimagk
               auarr(16,irow+1,icol+1,ixinterp)=-nn*qsurfp*apreal0
               auarr(16,irow+1,icol  ,ixinterp)=-nn*qsurfp*apimag0
               auarr(16,irow  ,icol+1,ixinterp)=nn*qsurfp*apimag0
               auarr(17,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*apprealmk
               auarr(17,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*appimagmk
               auarr(17,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*appimagmk
               auarr(18,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*appreal2mk
               auarr(18,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*appimag2mk
               auarr(18,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*appimag2mk
               auarr(19,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*apprealm2k
               auarr(19,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*appimagm2k
               auarr(19,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*appimagm2k
            endif
!----------------------------------------------------------------------
! now do same thing for inertia terms....note these will be scaled by 
! gamsq (square of growth rate) in the shooting code
!---------------------------------------------------------------------
            irealm=m**2*kcos(39,ip)/nn**3                           &
                   -m*ksin(61,ip)/nn**3                             &
                   -2.*m*ksin(62,ip)/nn**3                           &
                   -m*ksin(64,ip)/nn**3
            iimagm=m**2*ksin(39,ip)/nn**3                              &
                   +m*kcos(61,ip)/nn**3                             &
                   +2.*m*kcos(62,ip)/nn**3                          &
                   +m*kcos(64,ip)/nn**3 
            irealk=m**2*kcos(39,ip)/nn**3                           &
                   -m*ksin(61,ip)/nn**3                             &
                   -m*ksin(64,ip)/nn**3
            iimagk=m**2*ksin(39,ip)/nn**3                           &
                   +m*kcos(61,ip)/nn**3                             &
                   +m*kcos(64,ip)/nn**3
            ireal0=-qsurf*kcos(32,ip)/f_loc                         &
                  +m*qsurf*ksin(36,ip)/nn**2                        &
                  +m*qsurf*ksin(37,ip)/(f_loc*nn**2)                &
                  -qsurf*m**2*kcos(38,ip)/(f_loc*nn**2)             &
                  +2.*m**2*qsurfp*ksin(40,ip)/(qsurf*nn**3)         &
                  +m*(qsurfpp-2.*qsurfp**2/qsurf)*kcos(41,ip)/(qsurf*nn**3) &
                  -m*kcos(42,ip)/nn**3                              &
                  +2.*m*qsurfp*kcos(44,ip)/(nn**3*qsurf)            &
                  +qsurf*kcos(63,ip)/nn**2                          &
                  -m*kcos(65,ip)/nn**3                              &
                  +m*qsurfp*kcos(66,ip)/(nn**3*qsurf)
            iimag0=-qsurf*ksin(32,ip)/f_loc                         &
                   -m*qsurf*kcos(36,ip)/nn**2                       &
                   -m*qsurf*kcos(37,ip)/(f_loc*nn**2)               &
                   -qsurf*m**2*ksin(38,ip)/(f_loc*nn**2)            &
                   -2.*m**2*qsurfp*kcos(40,ip)/(qsurf*nn**3)        &
                  +m*(qsurfpp-2.*qsurfp**2/qsurf)*ksin(41,ip)/(qsurf*nn**3) &
                  -m*ksin(42,ip)/nn**3                              &
                  +2.*m*qsurfp*ksin(44,ip)/(nn**3*qsurf)            &
                  +qsurf*ksin(63,ip)/nn**2                          &
                  -m*ksin(65,ip)/nn**3                              &
                  +m*qsurfp*ksin(66,ip)/(nn**3*qsurf)
            iprealk=-2.*m*ksin(40,ip)/nn**3                       &
                     -kcos(66,ip)/nn**3
            ipimagk=2.*m*kcos(40,ip)/nn**3                          &
                         -ksin(66,ip)/nn**3
            iprealm=-2.*m*ksin(40,ip)/nn**3                         &
                   -2.*kcos(44,ip)/nn**3                            &
                   -kcos(66,ip)/nn**3
            ipimagm=2.*m*kcos(40,ip)/nn**3                          &
                   -2.*ksin(44,ip)/nn**3                            &
                   -ksin(66,ip)/nn**3
            ipreal0=qsurf*kcos(35,ip)/nn**2                         &
                    +2.*m*qsurf*ksin(34,ip)/(f_loc*nn**2)           &
                    +2.*m*qsurfp*kcos(41,ip)/(nn**3*qsurf)
            ipimag0=qsurf*ksin(35,ip)/nn**2 &
                    -2.*m*qsurf*kcos(34,ip)/(f_loc*nn**2) &
                    +2.*m*qsurfp*ksin(41,ip)/(nn**3*qsurf)
            ipprealk=-kcos(41,ip)/nn**3
            ippimagk=-ksin(41,ip)/nn**3
            ipprealm=-kcos(41,ip)/nn**3
            ippimagm=-ksin(41,ip)/nn**3
            ippreal0=qsurf*kcos(33,ip)/(f_loc*nn**2)
            ippimag0=qsurf*ksin(33,ip)/(f_loc*nn**2)
            iuarr(2,irow  ,icol  ,ixinterp)=irealm
            iuarr(1,irow  ,icol  ,ixinterp)=irealk
            iuarr(3,irow  ,icol  ,ixinterp)=ireal0
            iuarr(5,irow  ,icol  ,ixinterp)=  &
                     -nn*(qsurfp*iprealm+qsurfpp*ipprealm)
            iuarr(4,irow  ,icol  ,ixinterp)=  &
                     -nn*(qsurfp*iprealk+qsurfpp*ipprealk)
            iuarr(6,irow  ,icol  ,ixinterp)=  &
                     -nn*(qsurfp*ipreal0+qsurfpp*ippreal0)
            iuarr(8,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*ipprealm
            iuarr(7,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*ipprealk
            iuarr(9,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*ippreal0
            if(.not.updownsym) then
               iuarr(2,irow+1,icol+1,ixinterp)=irealm
               iuarr(2,irow+1,icol  ,ixinterp)=iimagm
               iuarr(2,irow  ,icol+1,ixinterp)=-iimagm
               iuarr(1,irow+1,icol+1,ixinterp)=irealk
               iuarr(1,irow+1,icol  ,ixinterp)=iimagk
               iuarr(1,irow  ,icol+1,ixinterp)=-iimagk
               iuarr(3,irow+1,icol+1,ixinterp)=ireal0
               iuarr(3,irow+1,icol  ,ixinterp)=iimag0
               iuarr(3,irow  ,icol+1,ixinterp)=-iimag0
               iuarr(5,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*iprealm+qsurfpp*ipprealm)
               iuarr(5,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*ipimagm+qsurfpp*ippimagm)
               iuarr(5,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*ipimagm+qsurfpp*ippimagm)
               iuarr(4,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*iprealk+qsurfpp*ipprealk)
               iuarr(4,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*ipimagk+qsurfpp*ippimagk)
               iuarr(4,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*ipimagk+qsurfpp*ippimagk)
               iuarr(6,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*ipreal0+qsurfpp*ippreal0)
               iuarr(6,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*ipimag0+qsurfpp*ippimag0)
               iuarr(6,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*ipimag0+qsurfpp*ippimag0)
               iuarr(8,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*ipprealm
               iuarr(8,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*ippimagm
               iuarr(8,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*ippimagm
               iuarr(7,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*ipprealk
               iuarr(7,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*ippimagk
               iuarr(7,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*ippimagk
               iuarr(9,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*ippreal0
               iuarr(9,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*ippimag0
               iuarr(9,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*ippimag0
            endif



!-----------------------------------------------------------------------
! surface quantities for boundary condition
!-----------------------------------------------------------------------
!            if(dpsi.eq.0) then
            if((ixinterp==1).and.(bug(2)<10.)) then
!               if (m.eq.mmin .and. mp.eq.mmin) &
               if (m.eq.mlo .and. mp.eq.mlo) then
                  if (verbose .ge. 3) write(outfile,*) 'including surface terms'
               endif
               if(mp.eq.m) then ! yes this gets defined multiple times
                  b2inv_edge=kcos(13,ip)/(2.*pi) !for no good reason
!                  write(6,*) "b2inv_edge=",b2inv_edge
               endif
! 3/13/98 put fprime piece back in here--no longer eigenvalue
               asurfre=(-m*qsurfp*kcos(4,ip)/(nn*qsurf**2)         &
                              +f_loc*pprime_loc*kcos(13,ip)              &
                              +ffprime_loc*kcos(12,ip)/f_loc             &
                              +ksin(24,ip)/nn                            &
                              +ksin(26,ip)/nn   )                        &
                      +dmpnq*(m*qsurfp*kcos(23,ip)/(nn**2*qsurf**3)      &
                              +kcos(20,ip)/(nn*qsurf)   )                &
                      +dmnq*(kcos(3,ip)/(nn*qsurf)                       &
                              +m*ksin(6,ip)/(nn*qsurf)                   &
                              +kcos(18,ip)/(nn*qsurf)                    &
                              +2.*m*qsurfp*kcos(23,ip)/(nn**2*qsurf**3)  &
                              -mp*kcos(28,ip)/(nn**2*qsurf**2)           &
                              -ksin(50,ip)/(nn*qsurf)**2   )             &
                      +dmnq*dmpnq*(-m*ksin(21,ip)/(nn**2*qsurf**2)       &
                              -kcos(28,ip)/(nn*qsurf)**2   )             &
                      +dmnq*dmnq*(-m*ksin(21,ip)/(nn*qsurf)**2           &
                              -kcos(28,ip)/(nn*qsurf)**2   )
               asurfim=(-m*qsurfp*ksin(4,ip)/(nn*qsurf**2)         &
                              +f_loc*pprime_loc*ksin(13,ip)              &
                              +ffprime_loc*ksin(12,ip)/f_loc             &
                              -kcos(24,ip)/nn                            &
                              -kcos(26,ip)/nn)                           &
                      +dmpnq*(m*qsurfp*ksin(23,ip)/(nn**2*qsurf**3)      &
                              +ksin(20,ip)/(nn*qsurf)       )            &
                      +dmnq*(ksin(3,ip)/(nn*qsurf)                       &
                              -m*kcos(6,ip)/(nn*qsurf)                   &
                              +ksin(18,ip)/(nn*qsurf)                    &
                              +2.*m*qsurfp*ksin(23,ip)/(nn**2*qsurf**3)  &
                              -mp*ksin(28,ip)/(nn**2*qsurf**2)           &
                              +kcos(50,ip)/(nn*qsurf)**2)                &
                      +dmnq*dmpnq*(m*kcos(21,ip)/(nn**2*qsurf**2)       &
                              -ksin(28,ip)/(nn*qsurf)**2)                &
                      +dmnq*dmnq*(m*kcos(21,ip)/(nn*qsurf)**2           &
                              -ksin(28,ip)/(nn*qsurf)**2)
               asurfre=dmpnq*asurfre
               asurfim=dmpnq*asurfim
               apsurfre=dmpnq*(dmnq*kcos(4,ip)/(nn*qsurf) &
                    -dmnq*(dmnq+dmpnq)*kcos(23,ip)/(nn*qsurf)**2)
               apsurfim=dmpnq*(dmnq*ksin(4,ip)/(nn*qsurf) &
                    -dmnq*(dmnq+dmpnq)*ksin(23,ip)/(nn*qsurf)**2)
! inertia terms...
               isurfre=-qsurf*m*ksin(34,ip)/(f_loc*nn) &
                   -qsurf*kcos(43,ip)/nn &
                   +m*ksin(40,ip)*(dmpnq+dmnq)/nn**2 &
                   +m*kcos(44,ip)/nn**2 &
                   -m*qsurfp*kcos(41,ip)/(nn**2*qsurf)
               isurfim=qsurf*m*kcos(34,ip)/(f_loc*nn) &
                   -qsurf*ksin(43,ip)/nn &
                   -m*kcos(40,ip)*(dmpnq+dmnq)/nn**2 &
                   +m*ksin(44,ip)/nn**2 &
                   -m*qsurfp*ksin(41,ip)/(nn**2*qsurf)
               ipsurfre=-qsurf*kcos(33,ip)/(f_loc*nn) &
                   +kcos(41,ip)*(dmnq+dmpnq)/nn**2
               ipsurfim=-qsurf*ksin(33,ip)/(f_loc*nn) &
                   +ksin(41,ip)*(dmnq+dmpnq)/nn**2
               asurfu(irow  ,icol  )= asurfre
               asurfup(irow  ,icol  )=-nn*qsurfp*apsurfre
               isurfu(irow  ,icol  )= isurfre
               isurfup(irow  ,icol  )=-nn*qsurfp*ipsurfre
               if(.not.updownsym) then
                  asurfu(irow+1,icol+1)= asurfre
                  asurfu(irow+1,icol  )= asurfim
                  asurfu(irow  ,icol+1)=-asurfim
                  asurfup(irow+1,icol+1)=-nn*qsurfp*apsurfre
                  asurfup(irow+1,icol  )=-nn*qsurfp*apsurfim
                  asurfup(irow  ,icol+1)= nn*qsurfp*apsurfim
                  isurfu(irow+1,icol+1)= isurfre
                  isurfu(irow+1,icol  )= isurfim
                  isurfu(irow  ,icol+1)=-isurfim
                  isurfup(irow+1,icol+1)=-nn*qsurfp*ipsurfre
                  isurfup(irow+1,icol  )=-nn*qsurfp*ipsurfim
                  isurfup(irow  ,icol+1)= nn*qsurfp*ipsurfim
               endif
            endif
         end do
      end do
      if (ibench.eq.1) then
         if (ixinterp.eq.1) write(76,*)np
         do i=1,nk
            do ip=1,np
               write(76,92)kcos(i,ip),ksin(i,ip)
            end do
         end do
      end if
 92   format(3e14.6)
    return
      
end subroutine nustuff





subroutine oldnustuff(ixinterp)

! set up matrices for multiple flux surface equilibrium 2/8/00
!  translated to F90 1/31/00
!-----------------------------------------------------------------------
!     equilibrium quantities are constructed and fourier integrals
!     performed
!     this routine has been largely rewritten to take advantage of
!     Howard's integration by parts
!     changed various signs below to convert to sign convention with
!     psi=minimum on axis. B_p=Grad(phi) x Grad(psi)
!-----------------------------------------------------------------------
      use elite_data, only: ns,mmax,mmin,pprime_adj,ffprime_adj,f_eq, &
           gamscl,rl,dens,ne_eq,neprime_eq,outfile,bpl,uprimel,sinul, &
           d2psi1l,cosul,dpsi1l,ppp_eq,ffpp_eq,circumf,pi,q_eq,qp_eq, &
           qpp_eq,xinterp,m0,nn,q_calc,del,f_in,q_in,nxinterp,zl,n1term, &
           nowindow,nmwinhalf,maxdxinterp,mwindow,mindex,auarr, &
           updownsym,iuarr,bug,b2inv_edge,asurfu, &
           asurfup,isurfu,isurfup,runname,lrunname,gamscl2,funcal,verbose, &
           ti_eq,tiprime_eq,omegas,omegasn,omegaspi,ion_mass,kinkterms
      implicit none
      integer ixinterp
      integer ir,ic,ibench
      real nu(ns),nup(ns),nupp(ns)
      real check1(ns),check2(ns)
      integer i,ik
      integer kk,jj,km
      integer maxrow
      parameter (kk=56)
      real test(ns),nutest(ns),testi2,fval
      real dchi
      real switchon,test1,test2
      real nu1,nu2
      real nuhat,dbpp,hfun
      real omeg(ns),omegp(ns),omegpp(ns)
      real sig(ns),sigp(ns)        ! parallel current
      real dfac(ns)           ! d/dpsi(R^2 B_p^2/(J B^2))
      real d2fac(ns)           ! (1/nu) d/dpsi(J R^2 B_p^2/B^2)
      real d3fac(ns)         ! d/dpsi (B_p/(R*B^2))^2
      real drbdchi(ns)            ! d/dchi (1/(R^2 B^2))
      real domprbdchi(ns)     ! d/dchi (omega^prime/(R^2 B^2))
      real dompbsqdchi(ns)     ! d/dchi (omega^prime/(B^2))
      real dsigpchi(ns)            ! (d/dchi) sigma-prime
      real dbfacdom(ns)        !(d/domega)(f^3Bp^2/(R^2B^4))
      real qsurf,qsurfp,qsurfpp
      real bpval,rval,btval,zval
      real jval(ns),bsqval(ns)
      real dbsqp(ns),bpsqval(ns),rsqval(ns)
      integer im,ip,np,nk
      real bsqchi(ns)
      real phase,cosphase(2*(mmax-mmin)+1),sinphase(2*(mmax-mmin)+1)
      real kern(kk)
      real ksin(kk,2*(mmax-mmin)+1),kcos(kk,2*(mmax-mmin)+1)
      real dmnq,dmpnq
      real ucoefr1,ucoefr2,ucoefr3,ucoefr4,ucoefr7,ucoefr9,ucoefr11
      real ucoefi5,ucoefi6,ucoefi8
      real upcoefr3,upcoefr4
      real upcoefi6,upcoefi10
      real uppcoefr4
      integer m,mp,irow,icol,mtop,mlo
      real areal,aimag,apreal,apimag,appreal,appimag
      real arealdd,aimagdd,aprealdd,apimagdd,apprealdd,appimagdd
      real areald,aimagd,apreald,apimagd,appreald,appimagd
      real arealnd,aimagnd,aprealnd,apimagnd,apprealnd,appimagnd
      real arealmk,areal2mk,arealm2k,areal2m,areal2k,arealm
      real arealk,areal0
      real aimagmk,aimag2mk,aimagm2k,aimag2m,aimag2k,aimagm
      real aimagk,aimag0
      real aprealmk,apreal2mk,aprealm2k,apreal2m,apreal2k,aprealm
      real aprealk,apreal0
      real apimagmk,apimag2mk,apimagm2k,apimag2m,apimag2k,apimagm
      real apimagk,apimag0
      real apprealmk,appreal2mk,apprealm2k
      real appimagmk,appimag2mk,appimagm2k
      real irealm,irealk,ireal0
      real iimagm,iimagk,iimag0
      real iprealm,iprealk,ipreal0
      real ipimagm,ipimagk,ipimag0
      real ipprealm,ipprealk,ippreal0
      real ippimagm,ippimagk,ippimag0
      real asurfre,asurfim,apsurfre,apsurfim
      real isurfre,isurfim,ipsurfre,ipsurfim
      real domega
!      real check1,check2,check3
      real dbpdrho,dbtdrho,drhdpsi
      real pprime_loc,ffprime_loc,psi_1,psi_2,psi_3,f_loc,d_loc,dp_loc
      real Rc,psichk1,psichk2,nuppchk,qppfact,qpfact
      real gamsqws,gamsqwsn

       ibench=0      ! 20/4/01  writes benchmarking file if set to 1 (HRW)
       pprime_loc=pprime_adj(ixinterp)   ! 3/1/00 use adjusted values
       ffprime_loc=ffprime_adj(ixinterp)  ! 3/1/00 use adjusted values
       f_loc=f_eq(ixinterp)
       if (ixinterp.eq.1) then
         gamscl=(rl(1)+rl(ns/2))**2/(f_loc/rl(1)+f_loc/rl(ns/2))**2
         gamscl2=(0.5*(rl(1)+rl(ns/2)))**4/f_loc**2
       end if

! local plasma density, and psi-derivative
       if (dens) then
         d_loc=ne_eq(ixinterp)/ne_eq(1)
         dp_loc=neprime_eq(ixinterp)/ne_eq(1)
       else
         d_loc=1.
         dp_loc=0.
       end if

!    simple def of omega_*_i
       if (dens) then
          omegas(ixinterp)=2.9979e10*nn/ne_eq(ixinterp)/4.8032e-10* &
               (pprime_loc/(2.*4.*pi))  ! pprime_loc=4pi pprime
          omegasn(ixinterp)=2.9979e10*nn/ne_eq(ixinterp)/4.8032e-10* &
               (neprime_eq(ixinterp)*1.6022e-12*ti_eq(ixinterp))
          omegaspi(ixinterp)=omegasn(ixinterp)+ &
               2.9979e10*nn/ne_eq(ixinterp)/4.8032e-10* &
               (tiprime_eq(ixinterp)*1.6022e-12*ne_eq(ixinterp))
!  use ion_mass for relative ion mass here
!  4/06 correct to use edge density norm
          gamsqws=pi*ion_mass*1.6726e-24*ne_eq(1)* &
               omegas(ixinterp)**2
          gamsqwsn=pi*ion_mass*1.6726e-24*ne_eq(1)* &
               omegasn(ixinterp)**2
!          write(*,*) 'omegas=',omegas(ixinterp),'gamsqws=',gamsqws, &
!               'ion_mass=',ion_mass
          if (verbose .ge. 2) then
             write(outfile,*) 'omegas=',omegas(ixinterp),'gamsqws=',gamsqws, &
                  'ion_mass=',ion_mass
             write(outfile,*) 'omegasn=',omegasn(ixinterp),' omegaspi=', &
                  omegaspi(ixinterp)
          endif
       endif

       if (verbose .ge. 3) write(outfile,*) &
            'surface=',ixinterp,' pprime=',pprime_loc, &
            ' ffprime=',ffprime_loc,' f=',f_loc
!-----------------------------------------------------------------------
!     1.0 define nu and first two psi derivatives
!     explicitly extract the psi dependence of nu to facilitate
!     taking psi derivatives
!     nu=f*L/(2*pi*R^2)*(1+nuhat*dpsi)/(bp+dbpp*dpsi)
!     these values have been checked 27/10/97
!     (10/20/00 nupp not being calculated suffieciently accurately...
!     Miller expansion improved to include O(rho^2) corrections for nupp
!     (see Miller et al PoP 5 (1998) 973)
!-----------------------------------------------------------------------
!      write(45,*) dpsi
      do i=1,ns
        psi_1=rl(i)*bpl(i)
        psi_2=0.5*((-rl(i)*uprimel(i)+sinul(i))*bpl(i)- &
          rl(i)**2*pprime_loc-ffprime_loc)
        psi_3=(1./6.)*((-2.*bpl(i)*sinul(i)+4.*psi_2+ffprime_loc)* &
          (sinul(i)/rl(i)-uprimel(i))-rl(i)**2*pprime_loc* &
          (sinul(i)/rl(i)+uprimel(i))-d2psi1l(i)+cosul(i)* &
          dpsi1l(i)/rl(i)-rl(i)*bpl(i)*(rl(i)**2*ppp_eq(ixinterp) &
          +ffpp_eq(ixinterp)))
        Rc=-1./uprimel(i)
        nu(i)= (f_loc*circumf)/(2.*pi*bpl(i)*rl(i)**2)
        nu1=-sinul(i)/rl(i)+uprimel(i)-2.*psi_2/psi_1
! 11/28/00 eliminate the dpsi1l term which cancels when the higher
!   order correction to dl_p is properly considered in the nu'' calculation
        nu2=(sinul(i)/rl(i)+2.*psi_2/psi_1)*(sinul(i)/rl(i)-uprimel(i)) &
            +4.*psi_2**2/psi_1**2-3.*psi_3/psi_1
            !-0.5*(dpsi1l(i)/psi_1)**2
        nup(i)=nu(i)*(ffprime_loc/f_loc**2+nu1/psi_1)
        nupp(i)=nu(i)*((ffpp_eq(ixinterp)-(ffprime_loc/f_loc)**2)/ &
          f_loc**2 +2.*ffprime_loc*nu1/(psi_1*f_loc**2)+ &
          2.*(nu2-nu1*psi_2/psi_1)/psi_1**2)
      end do
!-----------------------------------------------------------------------
!     1.1  can use nu to calculate q's and omega's
!     checked 27/10/97
!-----------------------------------------------------------------------
      omeg(1)=0.
      omegp(1)=0.
      omegpp(1)=0.
      dchi=2.*pi/(ns-1.)
      do i=2,ns
         omeg(i)  =omeg(i-1)  +dchi*(nu(i)  +nu(i-1)  )*0.5
         omegp(i) =omegp(i-1) +dchi*(nup(i) +nup(i-1) )*0.5
         omegpp(i)=omegpp(i-1)+dchi*(nupp(i)+nupp(i-1))*0.5
      end do
      qsurf=omeg(ns)/(2.*pi)
      qsurfp=omegp(ns)/(2.*pi)
      qsurfpp=omegpp(ns)/(2.*pi)
!      write(33,*)' surface=',ixinterp
!      write(33,*) "q/q'/q''=",qsurf,qsurfp,qsurfpp
!      write(33,*)" equilibrium q/q'/q''",q_eq(ixinterp), &
!        qp_eq(ixinterp),qpp_eq(ixinterp)

! temporary fix!!!!  multiply vupp by factor of qpp_eq/qsurfpp to 
!  get equilibrium qpp (which is approximately the derivative of
!  qsurf
      qppfact=qpp_eq(ixinterp)/qsurfpp
      qpfact=qp_eq(ixinterp)/qsurfp
!      do i=1,ns
!         nupp(i)=nupp(i)*qppfact
!         nup(i)=nup(i)*qpfact
!      enddo
!      do i=2,ns
!        omegpp(i)=omegpp(i-1)+dchi*(nupp(i)+nupp(i-1))*0.5
!        omegp(i) =omegp(i-1) +dchi*(nup(i) +nup(i-1) )*0.5
!     enddo
!     qsurfpp=omegpp(ns)/(2.*pi)
!     qsurfp=omegp(ns)/(2.*pi)
!     write(33,*) 'qppfact=',qppfact,'adjusted qsurfpp=',qsurfpp
!     write(33,*) 'qpfact=',qpfact,'adjusted qsurfp=',qsurfp

      xinterp(ixinterp)=m0-nn*qsurf
      if (verbose .ge. 3) write(outfile,*) &
           'qsurf=',qsurf,' q_calc=',q_calc(ixinterp), &
           'xinterp=',xinterp(ixinterp)
      if(ixinterp==1) then
!! 11/00 making change to make sure grid is consistent for all xinterps
!         del=xinterp(1)
!         xmin=del
!         if (min0(m0+nmwinhalf-int(xmin+0.5),mmax).ne.m1up) then
!            write(*,*) 'fixing roundoff in xmin has changed m1up'
!            stop
!         endif
!         write(*,*) 'adjusted del for roundoff=',del
         xinterp(ixinterp)=del
         if (verbose .ge. 2) write(outfile,*) &
              'reassigned xinterp=del to avoid roundoff problem'
!         write(75,*) 'm0-nn*qsurf=',m0-nn*qsurf, 'del=',del, &
!              (del-(m0-nn*qsurf))/del
!         write(75,*) 'f_loc=',f_loc,' f_in=',f_in, &
!              (f_in/f_loc)*qsurf, q_in,qsurf

!  TEMP!!!
!   change qsurf as well so it's consistent with xinterp=del
!         write to screen only if q changes significantly
         if ((qsurf-(m0-del)/nn)/qsurf > 0.0001) then
            write(*,*) 'modifying qsurf to get x=del on outer surface'
            write(*,*) 'q_in=',q_in,'qsurfold=',qsurf,'qsurfnew=',(m0-del)/nn
         endif
         if (verbose .ge. 2) then
            write(outfile,*) 'modifying qsurf to get x=del on outer surface'
            write(outfile,*) 'q_in=',q_in,'qsurfold=',qsurf, &
                 'qsurfnew=',(m0-del)/nn
         endif
         qsurf=(m0-del)/nn
!         write(*,*) 'q_in=',q_in,'qsurf=',qsurf

      endif
      if (ixinterp.eq.1 .and. funcal) then
         open(73,file=runname(1:lrunname)//'.omegadat' )
      end if
      if (ixinterp.eq.1) then
        if (ibench.eq.1) then
          open(76,file='bench.dat')
          write(76,*)ns,nxinterp
        end if
      end if
      if (ibench.eq.1) then
        write(76,91)f_loc,ffprime_loc,pprime_loc
        write(76,91)qsurf,qsurfp,qsurfpp
      end if
      do i=1,ns
         omeg(i)=omeg(i)/qsurf
         if (funcal) write(73,91)rl(i),zl(i),bpl(i),omeg(i)
         omegp(i)=-qsurfp*omeg(i)/qsurf+omegp(i)/qsurf
         omegpp(i)=(-qsurfpp*omeg(i)-2*qsurfp*omegp(i)+omegpp(i))/qsurf
         if (ibench.eq.1) then
           write(76,91)rl(i),zl(i),bpl(i)
           write(76,91)omeg(i),omegp(i),omegpp(i)
           write(76,91)nu(i),nup(i),nupp(i)
         end if
      end do
 91   format(3e14.6)
      if (ixinterp.eq.nxinterp .and. funcal) then
        close(73)
        rewind(73)
      end if
!-----------------------------------------------------------------------
!     1.2 chi dependent equilibrium quantities appearing in equations
!     checked rsqval,jval,bsqval,bpsqval,dbsqp 27/10/97
!-----------------------------------------------------------------------

! 2/00 rho is now always zero, as is dpsi.  still need derivatives

      do i=1,ns

         psi_1=rl(i)*bpl(i)
         psi_2=0.5*((sinul(i)+rl(i)*(-uprimel(i)))*bpl(i) &
             -rl(i)**2*pprime_loc-ffprime_loc)
         drhdpsi=1./psi_1
         dbpdrho=-uprimel(i)*bpl(i)-rl(i)*pprime_loc-ffprime_loc/rl(i)
         bpval=bpl(i)
         rval=rl(i)
         rsqval(i)=rval*rval
         jval(i)=circumf/(2.*pi)/bpval
         dbtdrho=(ffprime_loc*bpl(i)/f_loc-f_loc*sinul(i)/rl(i)**2)
         btval=f_loc/rl(i)
         bsqval(i)=bpval**2+btval**2
         bpsqval(i)=bpval*bpval
         dbsqp(i)=2.*(bpl(i)*dbpdrho+f_loc*dbtdrho/rl(i))*drhdpsi
         fval=f_loc
! parallel current and its psi-derivative (10/23/00)
         sig(i)=-f_loc*pprime_loc/bsqval(i)-ffprime_loc/f_loc
         sigp(i)=(f_loc*pprime_loc*dbsqp(i)/bsqval(i)- &
                  ffprime_loc*pprime_loc/f_loc- &
                  f_loc*ppp_eq(ixinterp))/bsqval(i) &
           +ffprime_loc**2/f_loc**3-ffpp_eq(ixinterp)/f_loc
         if (ibench.eq.1) then
           write(76,91)sig(i),sigp(i),dbsqp(i)
         end if
! (d/dpsi) (R^2 B_p^2 /(J B^2))
         dfac(i)=(4.*psi_2-psi_1*(uprimel(i)+sinul(i)/rl(i)- &
                 2.*psi_2/psi_1)-psi_1**2*dbsqp(i)/bsqval(i))/ &
                 (jval(i)*bsqval(i))
! (1/nu) (d/dpsi) (R^2 B_p^2 J/B^2)
         d2fac(i)=rsqval(i)*psi_1*(sinul(i)/rl(i)+uprimel(i)+ &
                  2.*psi_2/psi_1)/(f_loc*bsqval(i))- &
                  rsqval(i)*psi_1**2*dbsqp(i)/(f_loc*bsqval(i)**2)
! (d/dpsi) (B_p^2 /(R^2 B^4))
         d3fac(i)=2.*dbpdrho/(rl(i)**3*bsqval(i)**2)- &
                  2.*bpl(i)*sinul(i)/(rl(i)**4*bsqval(i)**2)- &
                  2.*bpl(i)**2*dbsqp(i)/(rl(i)**2*bsqval(i)**3)
! calculate zval for writing to diagnostic output file surfaces.out
!    pbs 11/99
         zval=zl(i)
!         write(47,*) rval,zval,bpval
      end do
!-----------------------------------------------------------------------
!     1.3 need chi deriv of B^2
!         and 1/(R^2 B^2)  and omegp/(R^2 B^2)  (10/23/00)
!-----------------------------------------------------------------------
      do i=1,ns
         ip=i+1
         im=i-1
         if(i.eq.ns) ip=2
         if(i.eq.1) im=ns-1
         bsqchi(i)=(bsqval(ip)-bsqval(im))/(2.*dchi)
         drbdchi(i)=(1./(rsqval(ip)*bsqval(ip))- &
             1./(rsqval(im)*bsqval(im)))/(2.*dchi)
         dompbsqdchi(i)=(omegp(ip)/bsqval(ip)-omegp(im)/bsqval(im))/ &
             (2.*dchi)
!         domprbdchi(i)=-((nu(i)*qsurfp/qsurf-nup(i))/ &
!             (qsurf*rsqval(i)*bsqval(i))+omegp(i)*drbdchi(i))
!         fix to match fix in new nustuff
         domprbdchi(i)=(omegp(ip)/(bsqval(ip)*rsqval(ip)) &
              -omegp(im)/(bsqval(im)*rsqval(im)))/ &
              (2.*dchi)
         dsigpchi(i)=(sigp(ip)-sigp(im))/(2.*dchi)
         dbfacdom(i)=(f_loc**3*qsurf/nu(i))*                         &
                     (bpl(ip)**2/(rsqval(ip)*bsqval(ip)**2)-  &
                  bpl(im)**2/(rsqval(im)*bsqval(im)**2))/(2.*dchi)
      end do
!-----------------------------------------------------------------------
!     2.0 chi integrals--zero the integral summations
!     kcos(ik,ip) is Integral[T(ik)*Cos[(ip-1+mmin-mmax)*omeg]]
!     So ik is index of T matrix appearing in notes
!     and ip is mp-m +(1+mmax-mmin)
!     mmax i maximum m value
!     mmin is minimum m value
!     ipmax is 2*(mmax-mmin)+1
!-----------------------------------------------------------------------
      np=2*(mmax-mmin)+1
      nk=kk
      do ip=1,np
         do ik=1,nk
            kcos(ik,ip)=0.
            ksin(ik,ip)=0.
         end do
      end do
!-----------------------------------------------------------------------
!     2.1 perform the integrations in omega
!-----------------------------------------------------------------------
      do i=1,ns-1
         domega=nu(i)*dchi/qsurf
         do ip=1,np
            phase=(ip-1+mmin-mmax)*omeg(i) 
            cosphase(ip)=cos(phase)
            sinphase(ip)=sin(phase)
         end do
         kern(1)=f_loc/(bpsqval(i)*rsqval(i)**2)
         kern(2)=(qsurf*bpsqval(i)*nupp(i)*rsqval(i))/(bsqval(i)*jval(i))
         kern(3)=(bpsqval(i)*nup(i)*rsqval(i))/(bsqval(i)*jval(i))
         kern(4)=(f_loc*bpsqval(i))/bsqval(i)
         kern(5)=(bpsqval(i)*nup(i)*omegp(i)*rsqval(i))/(bsqval(i)*jval(i))
         kern(6)=(f_loc*bpsqval(i)*omegp(i))/bsqval(i)
         kern(7)=(f_loc*bpsqval(i)*omegp(i)**2)/bsqval(i)
         kern(8)=(f_loc*bpsqval(i)*omegpp(i))/bsqval(i)
         kern(9)=((pprime_loc + dbsqp(i)*0.5)*rsqval(i))/bsqval(i)
         kern(10)=(bsqchi(i)*rsqval(i))/(bsqval(i)**2*jval(i))
         kern(11)=(bsqchi(i)*omegp(i)*rsqval(i))/(bsqval(i)**2*jval(i))
         kern(12)=1.
         kern(13)=1./bsqval(i)
         kern(14)=nup(i)*dfac(i)
         kern(15)=nu(i)*omegp(i)*dfac(i)
         kern(16)=nu(i)*dfac(i)
         kern(17)=sig(i)*pprime_loc/bsqval(i)
         kern(18)=f_loc*pprime_loc*omegp(i)*bpsqval(i)/bsqval(i)**2
         kern(19)=f_loc*pprime_loc*bpsqval(i)/bsqval(i)**2
         kern(20)=f_loc**2*sig(i)*omegp(i)/(rsqval(i)*bsqval(i))
         kern(21)=f_loc**2*sig(i)/(rsqval(i)*bsqval(i))
         kern(22)=f_loc**3*bpsqval(i)*omegp(i)/(rsqval(i)*bsqval(i)**2)
         kern(23)=f_loc**3*omegp(i)**2*bpsqval(i)/(rsqval(i)*bsqval(i)**2)
         kern(24)=f_loc**3*bpsqval(i)/(rsqval(i)*bsqval(i)**2)
         kern(25)=f_loc*pprime_loc*nup(i)*bpsqval(i)/(nu(i)*bsqval(i)**2)
         kern(26)=f_loc**2*sig(i)*nup(i)/(nu(i)*rsqval(i)*bsqval(i))
         kern(27)=pprime_loc*rsqval(i)*bpsqval(i)*bsqchi(i)/ &
                  (jval(i)*bsqval(i)**3)
         kern(28)=pprime_loc*rsqval(i)*bpsqval(i)*omegp(i)*bsqchi(i)/ &
                  (jval(i)*bsqval(i)**3)
         kern(29)=f_loc**2*sig(i)*drbdchi(i)/nu(i)
         kern(30)=f_loc**2*sig(i)*domprbdchi(i)/nu(i)
         kern(31)=f_loc**2*sig(i)*omegp(i)*drbdchi(i)/nu(i)
         kern(32)=f_loc**3*bpsqval(i)*nup(i)/(rsqval(i)*nu(i) &
                  *bsqval(i)**2)
         kern(33)=f_loc**3*bpsqval(i)*omegp(i)*nup(i)/ &
                  (rsqval(i)*bsqval(i)**2*nu(i))
         kern(34)=f_loc*pprime_loc*bpsqval(i)* dompbsqdchi(i)/ &
                  (nu(i)*bsqval(i))
         kern(35)=f_loc*pprime_loc*bpsqval(i)*bsqchi(i)/ &
                  (nu(i)*bsqval(i)**3)
         kern(36)=dsigpchi(i)/nu(i)
         kern(37)=sigp(i)
         kern(38)=d_loc/bpsqval(i)
         kern(39)=d_loc*rsqval(i)**2*bpsqval(i)/bsqval(i)
         kern(40)=d_loc*rsqval(i)**2*bpsqval(i)*omegp(i)/bsqval(i)
         kern(41)=dp_loc*rsqval(i)**2*bpsqval(i)/(f_loc*bsqval(i)) &
                  +d_loc*d2fac(i)
         kern(42)=dp_loc*omegp(i)*rsqval(i)**2*bpsqval(i)/ &
                  (f_loc*bsqval(i))+d_loc*omegp(i)*d2fac(i)
         kern(43)=d_loc*rsqval(i)**2*bpsqval(i)*omegpp(i)/bsqval(i)
         kern(44)=d_loc*rsqval(i)**2*bpsqval(i)*omegp(i)**2/bsqval(i)
         hfun=f_loc*rsqval(i)*bpsqval(i)/bsqval(i)**2
         kern(45)=d_loc*hfun*omegp(i)**2
         kern(46)=d_loc*hfun*omegp(i)
         kern(47)=d_loc*hfun
         kern(48)=d_loc*hfun*nupp(i)/nu(i)
         kern(49)=(d_loc*rsqval(i)**2*bpsqval(i)/(f_loc*bsqval(i)))* &
                  (pprime_loc/bsqval(i)+ &
                  nup(i)*f_loc*f_loc/(nu(i)*rsqval(i)*bsqval(i)))
         kern(50)=d_loc*hfun*nup(i)/nu(i)
!         kern(51)=3.*ffprime_loc*kern(22)/f_loc**2+kern(24)*omegpp(i) &
!                  +kern(24)*nup(i)*omegp(i)/nu(i)+ &
!                  f_loc**3*d3fac(i)*omegp(i)
         kern(51)=f_loc**3*bpsqval(i)*(nup(i)/(nu(i)*bsqval(i)))**2  &
                  /rsqval(i)
         kern(52)=f_loc**3*bpsqval(i)*omegpp(i)/(rsqval(i)*bsqval(i)**2)
         kern(53)=f_loc**3*bpsqval(i)*nupp(i)/                       &
                  (nu(i)*rsqval(i)*bsqval(i)**2)
         kern(54)=dbfacdom(i)*omegp(i)*nup(i)/nu(i)
         kern(55)=dbfacdom(i)*nupp(i)/nu(i)
         kern(56)=dbfacdom(i)*nup(i)/nu(i)
!!$         kern(52)=3.*ffprime_loc*kern(24)/f_loc**2+kern(32) &
!!$                  +f_loc**3*d3fac(i)
!!$         kern(53)=3.*ffprime_loc*kern(32)/f_loc**2+ &
!!$                  kern(24)*nupp(i)/nu(i)+f_loc**3*d3fac(i)*nup(i)/nu(i)
         switchon=0.
!         do ik=14,37
         do ik=17,35
           kern(ik)=n1term*kern(ik)
         end do
!         do ik=22,24
!           kern(ik)=switchon*kern(ik)
!         end do
!         do ik=32,33
!           kern(ik)=switchon*kern(ik)
!         end do
!         kern(41)=n1term*kern(41)
!         kern(42)=n1term*kern(42)
!         kern(43)=n1term*kern(43)
         do ik=45,56
           kern(ik)=n1term*kern(ik)
         end do
!         do ikern=1,13
!            kern(ikern)=bsqval(i)/bpsqval(i)*kern(ikern)
!         end do
!-----------------------------------------------------------------------
! debug--zero out smallest terms in ordering
!-----------------------------------------------------------------------
!$$$         kern(2)=0.
!$$$         kern(3)=0.
!$$$         kern(5)=0.
!$$$         kern(8)=0.
!$$$         if(dpsi.eq.0.) then
!$$$            write(61,*) (kern(jj),jj=1,13)
!$$$            write(61,*) omegp(i)
!$$$            write(61,*) nup(i)
!$$$            write(61,*) check1(i)
!$$$            write(61,*) check2(i)
!$$$         end if
         do ik=1,nk
            do ip=1,np
               kcos(ik,ip)=kcos(ik,ip)+domega*cosphase(ip)*kern(ik)
!         if (ik.eq.1) then
!           write(86,*)'ip=',ip,' kern(1)=',kern(1),' kcos=',kcos(ik,ip)
!         end if
               ksin(ik,ip)=ksin(ik,ip)+domega*sinphase(ip)*kern(ik)
            end do
!            if (ik.gt.50) then
!              write(6,*)' ik=',ik,' kcos=',kcos(ik,ip),' ksin=',ksin(ik,ip)
!            end if
         end do
      end do
!      write(86,*)' ixinterp=',ixinterp
!      do ik=1,nk
!        write(86,*)' ******k=',ik
!        do ip=1,np
!          write(86,*)' ip=',ip,' kcos=',kcos(ik,ip),' ksin=',ksin(ik,ip)
!        end do
!      end do
!      write(6,*) "ksin(6,1)/ksin(6,3)",ksin(6,1),ksin(6,3)
!-----------------------------------------------------------------------
!     3.0 calculate full x-dependent matrix equation for Um and b.c
!
!     au(irow,icol)*um(irow,icol)+
!     aup(irow,icol)*d/dx[um(irow,icol)]+
!     aupp(irow,icol)*d^2/dx^2[um(irow,icol)]=0 
!
!     asurfu(irow,icol)*um(irow,icol)+
!     asurfup(irow,icol)*d/dx[um(irow,icol)]=
!     gam*agamma(irow,icol)*um(irow,icol)
!     note that um={Re[um[1]],Im[um[1]],Re[um[2],Im[um[2]],...}
!-----------------------------------------------------------------------
! PS 7/01 arrays are now lwindow x lwindow in size, redo indexing
!   to go from mtop to mlo rather then mmax to mmin, allow padding
!   on top and bottom (if not against mmax or mmin) to allow proper
!   offset in matgen
      mtop=mmax
      mlo=mmin
      if(.not.nowindow) then
         ! note that mtop is maxdxinterp above usual mup.  need to leave room
         !  for offset in matgen
         mtop=min0(m0+nmwinhalf-int(xinterp(ixinterp)+0.5)+maxdxinterp,mmax)
         if (mtop < (mmin+mwindow/mindex-1+maxdxinterp) ) &
              mtop=mmin+mwindow/mindex-1+maxdxinterp
         mtop=min0(mtop,mmax)
         mlo=mtop-mwindow/mindex+1-2*maxdxinterp  ! reduce by two to allow room
         if (mlo < mmin) mlo=mmin
      endif
!      do m=mmin,mmax
      do m=mlo,mtop
         dmnq=m-nn*qsurf
         do mp=mlo,mtop
            ip=mp-m+1+mmax-mmin
            dmpnq=mp-nn*qsurf
! label k denotes mp
! terms to be multiplied by dmnq*dmpnq*u
            arealmk=kcos(1,ip)/qsurf-kcos(2,ip)/(nn*qsurf)**2 &
              -2.*m*ksin(5,ip)/(nn**2*qsurf)+m**2*kcos(7,ip)/(nn**2*qsurf) &
              -m*ksin(8,ip)/(nn**2*qsurf)-m*ksin(15,ip)/(nn**2*qsurf) &
              -4.*m*m*qsurfp*ksin(22,ip)/(nn*qsurf)**3 &
              -4.*m*qsurfp*kcos(32,ip)/(nn*qsurf)**3              &
!    extra terms added by HRW 4 May 01
              +mp*kcos(53,ip)/(nn**3*qsurf**2)                    &
              -2.*m*kcos(54,ip)/(nn**3*qsurf**2)                  &
              +ksin(55,ip)/(nn**3*qsurf**2)                       &
              -2.*(qsurfpp-2.*qsurfp**2/qsurf)*m*kcos(24,ip)/(nn*qsurf)**3  
            aimagmk=ksin(1,ip)/qsurf-ksin(2,ip)/(nn*qsurf)**2 &
              +2.*m*kcos(5,ip)/(nn**2*qsurf)+m**2*ksin(7,ip)/(nn**2*qsurf) &
              +m*kcos(8,ip)/(nn**2*qsurf)+m*kcos(15,ip)/(nn**2*qsurf) &
              +4.*m*m*qsurfp*kcos(22,ip)/(nn*qsurf)**3 &
              -4.*m*qsurfp*ksin(32,ip)/(nn*qsurf)**3              &
!    extra terms added by HRW 4 May 01
              +mp*ksin(53,ip)/(nn**3*qsurf**2)                    &
              -2.*m*ksin(54,ip)/(nn**3*qsurf**2)                  &
              -kcos(55,ip)/(nn**3*qsurf**2)                       &
              -2.*(qsurfpp-2.*qsurfp**2/qsurf)*m*ksin(24,ip)/(nn*qsurf)**3  
! terms to be multiplied by dmnq**2*dmpnq*u
            areal2mk=-m*m*kcos(23,ip)/(nn**3*qsurf**2) &
               +2.*m*ksin(33,ip)/(nn**3*qsurf**2)                 &
!  extra terms (HRW, 4 May)
               +m*ksin(52,ip)/(nn**3*qsurf**2)                    &
               +kcos(53,ip)/(nn**3*qsurf**2)
            aimag2mk=-m*m*ksin(23,ip)/(nn**3*qsurf**2) &
               -2.*m*kcos(33,ip)/(nn**3*qsurf**2)                 &
!  extra terms (HRW, 4 May)
               -m*kcos(52,ip)/(nn**3*qsurf**2)                    &
               +ksin(53,ip)/(nn**3*qsurf**2)
! terms to be multiplied by dmnq*dmpnq**2*u
            arealm2k=-m*m*kcos(23,ip)/(nn**3*qsurf**2)            &
               +4.*m*ksin(33,ip)/(nn**3*qsurf**2)                 &
!  extra terms (HRW, 4 May)
               +m*ksin(52,ip)/(nn**3*qsurf**2)                    &
               +kcos(53,ip)/(nn**3*qsurf**2)
            aimagm2k=-m*m*ksin(23,ip)/(nn**3*qsurf**2)            &
               -4.*m*kcos(33,ip)/(nn**3*qsurf**2)                 &
!  extra terms (HRW, 4 May)
               -m*kcos(52,ip)/(nn**3*qsurf**2)                    &
               +ksin(53,ip)/(nn**3*qsurf**2)
! terms to be multiplied by dmnq**2*u
            areal2m=m*ksin(20,ip)/(nn**2*qsurf)
            aimag2m=-m*kcos(20,ip)/(nn**2*qsurf)
! terms to be multiplied by dmpnq**2*u
            areal2k=-m*ksin(20,ip)/(nn**2*qsurf) &
              -2.*m*m*qsurfp*ksin(22,ip)/(nn*qsurf)**3           &
!   extra terms (HRW, 4 May)
              -4.*m*qsurfp*kcos(32,ip)/(nn*qsurf)**3          &
              -(qsurfpp-2.*qsurfp**2/qsurf)*m*kcos(24,ip)/(nn*qsurf)**3
            aimag2k=m*kcos(20,ip)/(nn**2*qsurf) &
              +2.*m*m*qsurfp*kcos(22,ip)/(nn*qsurf)**3           &
!   extra terms (HRW, 4 May)
              -4.*m*qsurfp*ksin(32,ip)/(nn*qsurf)**3          &
              -(qsurfpp-2.*qsurfp**2/qsurf)*m*ksin(24,ip)/(nn*qsurf)**3
! terms to be mulitplied by dmnq*u
            arealm=-kcos(17,ip)/nn-m*qsurfp*kcos(21,ip)/(nn*qsurf)**2 &
              +m*kcos(28,ip)/nn**2+m*kcos(30,ip)/nn**2 +kcos(37,ip)/nn
            aimagm=-ksin(17,ip)/nn-m*qsurfp*ksin(21,ip)/(nn*qsurf)**2 &
              +m*ksin(28,ip)/nn**2+m*ksin(30,ip)/nn**2+ksin(37,ip)/nn
! terms to be multiplied by dmpnq*u
            arealk=2.*m*qsurfp*kcos(3,ip)/(nn*qsurf)**2 &
              -2.*m*qsurfp**2/(nn**2*qsurf**3)*kcos(4,ip) &
              +m*qsurfpp*kcos(4,ip)/(nn*qsurf)**2 &
              +2.*m**2*qsurfp*ksin(6,ip)/(nn*qsurf)**2 &
              +m*qsurfp*kcos(16,ip)/(nn*qsurf)**2 &
              -kcos(17,ip)/nn-m*qsurfp*kcos(19,ip)/(nn*qsurf)**2 &
              +(kcos(25,ip)+kcos(26,ip))/nn+m*kcos(31,ip)/nn**2 &
              -m*kcos(34,ip)/nn**2 &
              +2.*m*m*qsurfp**2*kcos(24,ip)/(nn**3*qsurf**4)    &
!  Extra terms added (HRW 4 May)
              -2.*m*qsurfp*ksin(56,ip)/(nn*qsurf)**3
            aimagk=2.*m*qsurfp*ksin(3,ip)/(nn*qsurf)**2 &
              -2.*m*qsurfp**2*ksin(4,ip)/(nn**2*qsurf**3) &
              +m*qsurfpp*ksin(4,ip)/(nn*qsurf)**2 &
              -2.*m**2*qsurfp*kcos(6,ip)/(nn*qsurf)**2 &
              +m*qsurfp*ksin(16,ip)/(nn*qsurf)**2 &
              -ksin(17,ip)/nn-m*qsurfp*ksin(19,ip)/(nn*qsurf)**2 &
              +(ksin(25,ip)+ksin(26,ip))/nn+m*ksin(31,ip)/nn**2 &
              -m*ksin(34,ip)/nn**2 &
              +2.*m*m*qsurfp**2*ksin(24,ip)/(nn**3*qsurf**4)    &
!  Extra terms added (HRW 4 May)
              +2.*m*qsurfp*kcos(56,ip)/(nn*qsurf)**3
! terms to be multiplied by 1*u
            areal0=-2.*qsurf*pprime_loc*kcos(9,ip)/f_loc &
              +pprime_loc*qsurf*m*kcos(11,ip)/nn &
              +m*qsurfp*ksin(27,ip)/(nn**2*qsurf) &
              +m*qsurfp*ksin(29,ip)/(nn**2*qsurf) &
              -qsurf*ksin(36,ip)/nn
            aimag0=-2.*qsurf*pprime_loc*ksin(9,ip)/f_loc &
              +pprime_loc*qsurf*m*ksin(11,ip)/nn &
              -m*qsurfp*kcos(27,ip)/(nn**2*qsurf) &
              -m*qsurfp*kcos(29,ip)/(nn**2*qsurf) &
              +qsurf*kcos(36,ip)/nn
! terms to be multiplied by dmnq*dmpnq*u'
            aprealmk=-2.*kcos(3,ip)/(nn**2*qsurf) &
              -2.*m*ksin(6,ip)/(nn**2*qsurf)-kcos(16,ip)/(nn**2*qsurf) &
              -4.*m*qsurfp*kcos(24,ip)/(nn*qsurf)**3                   &
!  Extra terms added (HRW 4 May)
              +2.*ksin(56,ip)/(nn**3*qsurf**2)
            apimagmk=-2.*ksin(3,ip)/(nn**2*qsurf) &
              +2.*m*kcos(6,ip)/(nn**2*qsurf)-ksin(16,ip)/(nn**2*qsurf)  &
              -4.*m*qsurfp*ksin(24,ip)/(nn*qsurf)**3                    &
!  Extra terms added (HRW 4 May)
              -2.*kcos(56,ip)/(nn**3*qsurf**2)
!              -6.*m*qsurfp*ksin(24,ip)/(nn*qsurf)**3 &
!              +(m-mp)*ksin(32,ip)/(nn**3*qsurf**2)
! terms to be multiplied by dmnq**2*dmpnq*u'
            apreal2mk=2.*m*ksin(22,ip)/(nn**3*qsurf**2) &
              +2.*kcos(32,ip)/(nn**3*qsurf**2)
            apimag2mk=-2.*m*kcos(22,ip)/(nn**3*qsurf**2) &
              +2.*ksin(32,ip)/(nn**3*qsurf**2) 
! terms to be multiplied by dmnq*dmpnq**2*u'
            aprealm2k=2.*m*ksin(22,ip)/(nn**3*qsurf**2)           &
              +4.*kcos(32,ip)/(nn**3*qsurf**2)
            apimagm2k=-2.*m*kcos(22,ip)/(nn**3*qsurf**2)          &
              +4.*ksin(32,ip)/(nn**3*qsurf**2) 
! terms to be multiplied by dmnq**2*u'
            apreal2m=kcos(21,ip)/(nn**2*qsurf)
            apimag2m=ksin(21,ip)/(nn**2*qsurf)
! terms to be multiplied by dmpnq**2*u'
            apreal2k=-kcos(21,ip)/(nn**2*qsurf) &
               -2.*m*qsurfp*kcos(24,ip)/(nn*qsurf)**3
            apimag2k=-ksin(21,ip)/(nn**2*qsurf) &
               -2.*m*qsurfp*ksin(24,ip)/(nn*qsurf)**3
! terms to be mulitplied by dmnq*u'
            aprealm=-ksin(27,ip)/nn**2-ksin(29,ip)/nn**2
            apimagm=kcos(27,ip)/nn**2+kcos(29,ip)/nn**2
! terms to be multiplied by dmpnq*u'
            aprealk=2.*m*qsurfp*kcos(4,ip)/(nn*qsurf)**2 &
              -ksin(29,ip)/nn**2-ksin(35,ip)/nn**2
            apimagk=2.*m*qsurfp*ksin(4,ip)/(nn*qsurf)**2 &
              +kcos(29,ip)/nn**2+kcos(35,ip)/nn**2
! terms to be multiplied by 1*u'
            apreal0=-qsurf*pprime_loc*ksin(10,ip)/nn
            apimag0=qsurf*pprime_loc*kcos(10,ip)/nn
! terms to be multiplied by dmnq*dmpnq*u"
            apprealmk=-kcos(4,ip)/(nn**2*qsurf)
            appimagmk=-ksin(4,ip)/(nn**2*qsurf)
! terms to be multiplied by dmnq**2*dmpnq*u"
            appreal2mk=kcos(24,ip)/(nn**3*qsurf**2)
            appimag2mk=ksin(24,ip)/(nn**3*qsurf**2)
! terms to be multiplied by dmnq*dmpnq**2*u"
            apprealm2k=kcos(24,ip)/(nn**3*qsurf**2)
            appimagm2k=ksin(24,ip)/(nn**3*qsurf**2)
!       write(97,*)' ixinterp=',ixinterp,' m=',m,' k=',mp
!       write(97,*)' nn=',nn,' qsurf=',qsurf
!       if (apprealmk.eq.0.) then
!       write(97,*)' term1=',apprealmk,' term2=',apprealm2k
!      else
!       write(97,*)' term1=',apprealmk,' ratio=', &
!           apprealm2k*(mp+m-2.*nn*qsurf) &
!                    /apprealmk
!       end if
!     write(86,*)' Ar:',arealmk,areal2mk,arealm2k
!     write(86,*)'    ',areal2m,areal2k,arealm
!     write(86,*)'    ',arealk,areal0
!     write(86,*)' Ai:',aimagmk,aimag2mk,aimagm2k
!     write(86,*)'    ',aimag2m,aimag2k,aimagm
!     write(86,*)'    ',aimagk,aimag0
!     write(86,*)' Apr:',aprealmk,apreal2mk,aprealm2k
!     write(86,*)'    ',apreal2m,apreal2k,aprealm
!     write(86,*)'    ',aprealk,apreal0
!     write(86,*)' Api:',apimagmk,apimag2mk,apimagm2k
!     write(86,*)'    ',apimag2m,apimag2k,apimagm
!     write(86,*)'    ',apimagk,apimag0
!     write(86,*)' Appr:',apprealmk,appreal2mk,apprealm2k
!     write(86,*)' Appi:',appimagmk,appimag2mk,appimagm2k

 !-----------------------------------------------------------------------
!            irow=mindex*(mp-mmin)+1 !mindex=1 for updownsym
!            icol=mindex*(m -mmin)+1 !mindex=2 for.not.updownsym
!     used to use above statements where irow=1 was mmin but 
!     now irow=1 is mmax
!-----------------------------------------------------------------------
!            irow=mindex*(mmax-mp)+1 !mindex=1 for updownsym
            irow=mindex*(mtop-mp)+1 !mindex=1 for updownsym
!            icol=mindex*(mmax -m)+1 !mindex=2 for.not.updownsym
            icol=mindex*(mtop -m)+1 !mindex=2 for.not.updownsym
            auarr(1,irow  ,icol  ,ixinterp)=arealmk
            auarr(2,irow  ,icol  ,ixinterp)=areal2mk
            auarr(3,irow  ,icol  ,ixinterp)=arealm2k
            auarr(4,irow  ,icol  ,ixinterp)=areal2m
            auarr(5,irow  ,icol  ,ixinterp)=areal2k
            auarr(6,irow  ,icol  ,ixinterp)=arealm
            auarr(7,irow  ,icol  ,ixinterp)=arealk
            auarr(8,irow  ,icol  ,ixinterp)=areal0
            auarr(9,irow  ,icol  ,ixinterp)= &
                 -nn*(qsurfp*aprealmk+qsurfpp*apprealmk)
            auarr(10,irow  ,icol  ,ixinterp)= &
                 -nn*(qsurfp*apreal2mk+qsurfpp*appreal2mk)
            auarr(11,irow  ,icol  ,ixinterp)= &
                 -nn*(qsurfp*aprealm2k+qsurfpp*apprealm2k)
            auarr(12,irow  ,icol  ,ixinterp)=-nn*qsurfp*apreal2m
            auarr(13,irow  ,icol  ,ixinterp)=-nn*qsurfp*apreal2k
            auarr(14,irow  ,icol  ,ixinterp)=-nn*qsurfp*aprealm
            auarr(15,irow  ,icol  ,ixinterp)=-nn*qsurfp*aprealk
            auarr(16,irow  ,icol  ,ixinterp)=-nn*qsurfp*apreal0
            auarr(17,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*apprealmk
            auarr(18,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*appreal2mk
            auarr(19,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*apprealm2k
! I want to test Hermitian relations for dW minimization
!            write(45,'(2i5,3e15.7)') m,mp &
!                ,(arealdd*  dmnq+areald  )*dmpnq+arealnd &
!                ,(aprealdd* dmnq+apreald )*dmpnq+aprealnd &
!                ,(apprealdd*dmnq+appreald)*dmpnq+apprealnd
            if(.not.updownsym) then
               auarr(1,irow+1,icol+1,ixinterp)=arealmk
               auarr(1,irow+1,icol  ,ixinterp)=aimagmk
               auarr(1,irow  ,icol+1,ixinterp)=-aimagmk
               auarr(2,irow+1,icol+1,ixinterp)=areal2mk
               auarr(2,irow+1,icol  ,ixinterp)=aimag2mk
               auarr(2,irow  ,icol+1,ixinterp)=-aimag2mk
               auarr(3,irow+1,icol+1,ixinterp)=arealm2k
               auarr(3,irow+1,icol  ,ixinterp)=aimagm2k
               auarr(3,irow  ,icol+1,ixinterp)=-aimagm2k
               auarr(4,irow+1,icol+1,ixinterp)=areal2m
               auarr(4,irow+1,icol  ,ixinterp)=aimag2m
               auarr(4,irow  ,icol+1,ixinterp)=-aimag2m
               auarr(5,irow+1,icol+1,ixinterp)=areal2k
               auarr(5,irow+1,icol  ,ixinterp)=aimag2k
               auarr(5,irow  ,icol+1,ixinterp)=-aimag2k
               auarr(6,irow+1,icol+1,ixinterp)=arealm
               auarr(6,irow+1,icol  ,ixinterp)=aimagm
               auarr(6,irow  ,icol+1,ixinterp)=-aimagm
               auarr(7,irow+1,icol+1,ixinterp)=arealk
               auarr(7,irow+1,icol  ,ixinterp)=aimagk
               auarr(7,irow  ,icol+1,ixinterp)=-aimagk
               auarr(8,irow+1,icol+1,ixinterp)=areal0
               auarr(8,irow+1,icol  ,ixinterp)=aimag0
               auarr(8,irow  ,icol+1,ixinterp)=-aimag0
               auarr(9,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*aprealmk+qsurfpp*apprealmk)
               auarr(9,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*apimagmk+qsurfpp*appimagmk)
               auarr(9,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*apimagmk+qsurfpp*appimagmk)
               auarr(10,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*apreal2mk+qsurfpp*appreal2mk)
               auarr(10,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*apimag2mk+qsurfpp*appimag2mk)
               auarr(10,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*apimag2mk+qsurfpp*appimag2mk)
               auarr(11,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*aprealm2k+qsurfpp*apprealm2k)
               auarr(11,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*apimagm2k+qsurfpp*appimagm2k)
               auarr(11,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*apimagm2k+qsurfpp*appimagm2k)
               auarr(12,irow+1,icol+1,ixinterp)=-nn*qsurfp*apreal2m
               auarr(12,irow+1,icol  ,ixinterp)=-nn*qsurfp*apimag2m
               auarr(12,irow  ,icol+1,ixinterp)= nn*qsurfp*apimag2m
               auarr(13,irow+1,icol+1,ixinterp)=-nn*qsurfp*apreal2k
               auarr(13,irow+1,icol  ,ixinterp)=-nn*qsurfp*apimag2k
               auarr(13,irow  ,icol+1,ixinterp)= nn*qsurfp*apimag2k
               auarr(14,irow+1,icol+1,ixinterp)=-nn*qsurfp*aprealm
               auarr(14,irow+1,icol  ,ixinterp)=-nn*qsurfp*apimagm
               auarr(14,irow  ,icol+1,ixinterp)= nn*qsurfp*apimagm
               auarr(15,irow+1,icol+1,ixinterp)=-nn*qsurfp*aprealk
               auarr(15,irow+1,icol  ,ixinterp)=-nn*qsurfp*apimagk
               auarr(15,irow  ,icol+1,ixinterp)= nn*qsurfp*apimagk
               auarr(16,irow+1,icol+1,ixinterp)=-nn*qsurfp*apreal0
               auarr(16,irow+1,icol  ,ixinterp)=-nn*qsurfp*apimag0
               auarr(16,irow  ,icol+1,ixinterp)=nn*qsurfp*apimag0
               auarr(17,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*apprealmk
               auarr(17,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*appimagmk
               auarr(17,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*appimagmk
               auarr(18,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*appreal2mk
               auarr(18,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*appimag2mk
               auarr(18,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*appimag2mk
               auarr(19,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*apprealm2k
               auarr(19,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*appimagm2k
               auarr(19,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*appimagm2k
            endif
!----------------------------------------------------------------------
! now do same thing for inertia terms....note these will be scaled by 
! gamsq (square of growth rate) in the shooting code
!---------------------------------------------------------------------
            irealk=m**2*kcos(45,ip)/nn**3
            iimagk=m**2*ksin(45,ip)/nn**3
            irealm=m**2*kcos(45,ip)/nn**3
            iimagm=m**2*ksin(45,ip)/nn**3
            ireal0=-qsurf*kcos(38,ip)/f_loc &
                   +m*qsurf*ksin(42,ip)/nn**2 &
                   +m*qsurf*ksin(43,ip)/(f_loc*nn**2) &
                   -qsurf*m**2*kcos(44,ip)/(f_loc*nn**2) &
                   +2.*m**2*qsurfp*ksin(46,ip)/(qsurf*nn**3)
            iimag0=-qsurf*ksin(38,ip)/f_loc &
                   -m*qsurf*kcos(42,ip)/nn**2 &
                   -m*qsurf*kcos(43,ip)/(f_loc*nn**2) &
                   -qsurf*m**2*ksin(44,ip)/(f_loc*nn**2) &
                   -2.*m**2*qsurfp*kcos(46,ip)/(qsurf*nn**3)
            iprealk=-2.*m*ksin(46,ip)/nn**3
            ipimagk=2.*m*kcos(46,ip)/nn**3
            iprealm=-2.*m*ksin(46,ip)/nn**3
            ipimagm=2.*m*kcos(46,ip)/nn**3
            ipreal0=qsurf*kcos(41,ip)/nn**2 &
                    +2.*m*qsurf*ksin(40,ip)/(f_loc*nn**2) &
                    +2.*m*qsurfp*kcos(47,ip)/(nn**3*qsurf)
            ipimag0=qsurf*ksin(41,ip)/nn**2 &
                    -2.*m*qsurf*kcos(40,ip)/(f_loc*nn**2) &
                    +2.*m*qsurfp*ksin(47,ip)/(nn**3*qsurf)
            ipprealk=-kcos(47,ip)/nn**3
            ippimagk=-ksin(47,ip)/nn**3
            ipprealm=-kcos(47,ip)/nn**3
            ippimagm=-ksin(47,ip)/nn**3
            ippreal0=qsurf*kcos(39,ip)/(f_loc*nn**2)
            ippimag0=qsurf*ksin(39,ip)/(f_loc*nn**2)
!            write(89,*)'iprealm,ipimagm:',iprealm,ipimagm
!            write(89,*)'iprealk,ipimagk:',iprealk,ipimagk
!            write(89,*)'ipreal0,ipimag0:',ipreal0,ipimag0
!!$     irealm=0.
!!$     iimagm=0.
!!$     irealk=0.
!!$     iimagk=0.
!!$     ireal0=0.
!!$     iimag0=0.
!!$     iprealm=0.
!!$     ipimagm=0.
!!$     iprealk=0.
!!$     ipimagk=0.
!!$     ipreal0=0.
!!$     ipimag0=0.
            iuarr(2,irow  ,icol  ,ixinterp)=irealm
            iuarr(1,irow  ,icol  ,ixinterp)=irealk
            iuarr(3,irow  ,icol  ,ixinterp)=ireal0
            iuarr(5,irow  ,icol  ,ixinterp)=  &
                     -nn*(qsurfp*iprealm+qsurfpp*ipprealm)
            iuarr(4,irow  ,icol  ,ixinterp)=  &
                     -nn*(qsurfp*iprealk+qsurfpp*ipprealk)
            iuarr(6,irow  ,icol  ,ixinterp)=  &
                     -nn*(qsurfp*ipreal0+qsurfpp*ippreal0)
            iuarr(8,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*ipprealm
            iuarr(7,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*ipprealk
            iuarr(9,irow  ,icol  ,ixinterp)=(nn*qsurfp)**2*ippreal0
            if(.not.updownsym) then
               iuarr(2,irow+1,icol+1,ixinterp)=irealm
               iuarr(2,irow+1,icol  ,ixinterp)=iimagm
               iuarr(2,irow  ,icol+1,ixinterp)=-iimagm
               iuarr(1,irow+1,icol+1,ixinterp)=irealk
               iuarr(1,irow+1,icol  ,ixinterp)=iimagk
               iuarr(1,irow  ,icol+1,ixinterp)=-iimagk
               iuarr(3,irow+1,icol+1,ixinterp)=ireal0
               iuarr(3,irow+1,icol  ,ixinterp)=iimag0
               iuarr(3,irow  ,icol+1,ixinterp)=-iimag0
               iuarr(5,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*iprealm+qsurfpp*ipprealm)
               iuarr(5,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*ipimagm+qsurfpp*ippimagm)
               iuarr(5,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*ipimagm+qsurfpp*ippimagm)
               iuarr(4,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*iprealk+qsurfpp*ipprealk)
               iuarr(4,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*ipimagk+qsurfpp*ippimagk)
               iuarr(4,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*ipimagk+qsurfpp*ippimagk)
               iuarr(6,irow+1,icol+1,ixinterp)= &
                   -nn*(qsurfp*ipreal0+qsurfpp*ippreal0)
               iuarr(6,irow+1,icol  ,ixinterp)= &
                   -nn*(qsurfp*ipimag0+qsurfpp*ippimag0)
               iuarr(6,irow  ,icol+1,ixinterp)= &
                   nn*(qsurfp*ipimag0+qsurfpp*ippimag0)
               iuarr(8,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*ipprealm
               iuarr(8,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*ippimagm
               iuarr(8,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*ippimagm
               iuarr(7,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*ipprealk
               iuarr(7,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*ippimagk
               iuarr(7,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*ippimagk
               iuarr(9,irow+1,icol+1,ixinterp)=(nn*qsurfp)**2*ippreal0
               iuarr(9,irow+1,icol  ,ixinterp)=(nn*qsurfp)**2*ippimag0
               iuarr(9,irow  ,icol+1,ixinterp)=-(nn*qsurfp)**2*ippimag0
            endif



!-----------------------------------------------------------------------
! surface quantities for boundary condition
!-----------------------------------------------------------------------
!            if(dpsi.eq.0) then
            if((ixinterp==1).and.(bug(2)<10.)) then
!               if (m.eq.mmin .and. mp.eq.mmin) &
               if (m.eq.mlo .and. mp.eq.mlo) then
                  if (verbose .ge. 3) write(outfile,*) 'including surface terms'
               endif
               if(mp.eq.m) then ! yes this gets defined multiple times
                  b2inv_edge=kcos(13,ip)/(2.*pi) !for no good reason
!                  write(6,*) "b2inv_edge=",b2inv_edge
               endif
! 3/13/98 put fprime piece back in here--no longer eigenvalue
               asurfre=(dmnq*kcos(3,ip))/(nn*qsurf) &
                   -(m*qsurfp*kcos(4,ip))/(nn*qsurf**2) &
                   +ffprime_loc/f_loc*kcos(12,ip) &
                   +f_loc*pprime_loc*kcos(13,ip) &
                   +(dmnq*m*ksin(6,ip))/(nn*qsurf) &
                   +dmpnq*kcos(21,ip)/(nn*qsurf) &
                   +ksin(29,ip)/nn &
                   +m*qsurfp*(2.*dmnq+dmpnq)*kcos(24,ip)/(nn**2*qsurf**3) &
                   -m*dmnq*kcos(32,ip)/(nn*qsurf)**2 &
                   -m*dmnq*(dmnq+dmpnq)*ksin(22,ip)/(nn*qsurf)**2 &
                   +dmnq*kcos(19,ip)/(nn*qsurf) &
                   +qsurf*ksin(35,ip)/(nn*qsurf)
               asurfre=dmpnq*asurfre
               asurfim=-dmnq*m*kcos(6,ip)/(nn*qsurf) &
                   +dmnq*ksin(3,ip)/(nn*qsurf) &
                   -m*qsurfp*ksin(4,ip)/(nn*qsurf**2) &
                   +ffprime_loc/f_loc*ksin(12,ip) &
                   +f_loc*pprime_loc*ksin(13,ip) &
                   +dmpnq*ksin(21,ip)/(nn*qsurf) &
                   -kcos(29,ip)/nn &
                   +m*qsurfp*(2.*dmnq+dmpnq)*ksin(24,ip)/(nn**2*qsurf**3) &
                   -m*dmnq*ksin(32,ip)/(nn*qsurf)**2 &
                   +m*dmnq*(dmnq+dmpnq)*kcos(22,ip)/(nn*qsurf)**2 &
                   +dmnq*ksin(19,ip)/(nn*qsurf) &
                   -qsurf*kcos(35,ip)/(nn*qsurf)
               asurfim=dmpnq*asurfim
! need to agree with shifted circle
!               if(ip.ne.2) then
!                  asurfre=asurfre-f_loc*pprime_loc*kcos(13,ip)
!               endif
!
               apsurfre=dmpnq*(dmnq*kcos(4,ip)/(nn*qsurf) &
                    -dmnq*(dmnq+dmpnq)*kcos(24,ip)/(nn*qsurf)**2)
               apsurfim=dmpnq*(dmnq*ksin(4,ip)/(nn*qsurf) &
                    -dmnq*(dmnq+dmpnq)*ksin(24,ip)/(nn*qsurf)**2)
! inertia terms...
               isurfre=-qsurf*m*ksin(40,ip)/(f_loc*nn) &
                   -qsurf*kcos(49,ip)/nn &
                   +m*ksin(46,ip)*(dmpnq+dmnq)/nn**2 &
                   +m*kcos(50,ip)/nn**2 &
                   -m*qsurfp*kcos(47,ip)/(nn**2*qsurf)
               isurfim=qsurf*m*kcos(40,ip)/(f_loc*nn) &
                   -qsurf*ksin(49,ip)/nn &
                   -m*kcos(46,ip)*(dmpnq+dmnq)/nn**2 &
                   +m*ksin(50,ip)/nn**2 &
                   -m*qsurfp*ksin(47,ip)/(nn**2*qsurf)
               ipsurfre=-qsurf*kcos(39,ip)/(f_loc*nn) &
                   +kcos(47,ip)*(dmnq+dmpnq)/nn**2
               ipsurfim=-qsurf*ksin(39,ip)/(f_loc*nn) &
                   +ksin(47,ip)*(dmnq+dmpnq)/nn**2
               asurfu(irow  ,icol  )= asurfre
               asurfup(irow  ,icol  )=-nn*qsurfp*apsurfre
               isurfu(irow  ,icol  )= isurfre
               isurfup(irow  ,icol  )=-nn*qsurfp*ipsurfre
               if(.not.updownsym) then
                  asurfu(irow+1,icol+1)= asurfre
                  asurfu(irow+1,icol  )= asurfim
                  asurfu(irow  ,icol+1)=-asurfim
                  asurfup(irow+1,icol+1)=-nn*qsurfp*apsurfre
                  asurfup(irow+1,icol  )=-nn*qsurfp*apsurfim
                  asurfup(irow  ,icol+1)= nn*qsurfp*apsurfim
                  isurfu(irow+1,icol+1)= isurfre
                  isurfu(irow+1,icol  )= isurfim
                  isurfu(irow  ,icol+1)=-isurfim
                  isurfup(irow+1,icol+1)=-nn*qsurfp*ipsurfre
                  isurfup(irow+1,icol  )=-nn*qsurfp*ipsurfim
                  isurfup(irow  ,icol+1)= nn*qsurfp*ipsurfim
               endif
            endif
         end do
      end do
      if (ibench.eq.1) then
         if (ixinterp.eq.1) write(76,*)np
         do i=1,nk
            do ip=1,np
               write(76,92)kcos(i,ip),ksin(i,ip)
            end do
         end do
      end if
 92   format(3e14.6)
    return
      
end subroutine oldnustuff








