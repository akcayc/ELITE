subroutine hirsh(isurf,pprime_loc,ffprime_loc)

!     ****************
!
!  Bootstrap current as calculated by Hirshman (Phys.Fluids 31
!  (1988) 3150). Modified to include collisional effects.
!
      use elite_data, only: q_eq,pi,circumf,ns,bpl,f,rl,te,tau,ne,etai, &
           zeff,bstrap,rmajor,aspect,indj,vloop
      implicit none
!      include 'i_plas.f'
      real t1,t2,t3
      real bk,test,mu0
      real pe,pep,pip,ti,rj0
      real a1e,a1i,a2e,a2i
      real zm,dox,x,rl310,rl320,alfi0
      real a13,b13,c13,a23,b23,c23,f1,f2,fa
      real rl31,rl32
      real eps0,eq,emas,epsfac,rt,coolog,botcol
      real topcol,colte,vthe,colti,vthi
      real tnue,cnue,tnui,cnui
      real bsqmax,bsqmin,bsq,dl,ftrap,bsqav,riav
      real nor,bdif,binv,dst,bigint,pitch,rlag
      real bfac,bmax,bmin,alfi,dchi
      real fac,pro,colt,sigsp,rnust,sigfac,coneo
      integer nst
      integer j,jp
      integer isurf ! equilibrium surface
      real pprime_loc,ffprime_loc  ! ffprime_loc will be returned to surfcalc
      real qsurf
!
!      f=f_eq(isurf)
      qsurf=q_eq(isurf)
      mu0=4.*pi*1.0d-7
!  Some flux-surface averages....
      dchi=circumf/(ns-1.)
      nor=0.
      riav=0.
      bsqav=0.
      do j=1,ns
        dl=dchi/bpl(j)
        bsq=(f/rl(j))**2+bpl(j)**2
        bsqav=bsqav+dl*bsq
        riav=riav+dl/rl(j)**2
        if (j.eq.1) then
          bsqmin=bsq
          bsqmax=bsq
        else
          if (bsqmin.gt.bsq) bsqmin=bsq
          if (bsqmax.lt.bsq) bsqmax=bsq
        end if
        nor=nor+dl
      end do
      bsqav=bsqav/nor
      riav=riav/nor
      bmax=sqrt(bsqmax)
      bmin=sqrt(bsqmin)
      bdif=0.
      do j=1,ns
        bsq=(f/rl(j))**2+bpl(j)**2
        test=1.-sqrt(bsq/bsqmax)
        if (test.le.0.) write(6,*)' test=',test
        bdif=bdif+dchi*sqrt(1.-sqrt(bsq/bsqmax))/bpl(j)
      end do
      binv=1./sqrt(bsqmax)
      nst=100
      dst=binv/nst
      bigint=0.
      do jp=1,nst
        pitch=(jp-1)*dst+0.5*dst
        rlag=0.
        do j=1,ns
          bsq=sqrt((f/rl(j))**2+bpl(j)**2)
        test=1.-pitch*bsq
        if (test.le.0.) write(6,*)' test=',test,' pitch=',pitch, &
            ' bsq=',bsq
          rlag=rlag+dchi*sqrt(1.-pitch*bsq)/bpl(j)
        end do
        rlag=rlag/nor
        bigint=bigint+pitch*dst/rlag
      end do
      ftrap=1.-3.*bsqav*bigint/4.
      x=ftrap/(1.-ftrap)
      bk=1.602
      ti=te/tau
      pe=bk*ne*te
      pep=(pprime_loc/mu0)*tau/(1.+tau)
      pip=(pprime_loc/mu0)-pep
!  pressure and temperature gradient bstrap drives
      a1e=pep/pe
      a1i=tau*pip/pe
      a2e=etai*a1e/(1.+etai)
      a2i=a2e
      rj0=f*pe
!  Assume ion charge=1
      zm=1.
      dox=1.414*zm+zm*zm+x*(0.754+2.657*zm+2.*zm*zm)+x*x* &
           (0.348+1.243*zm+zm*zm)
      rl310=rj0*x*(0.754+2.21*zm+zm*zm+x*(0.348+1.243*zm+zm*zm))/dox
      rl320=-rj0*x*(0.884+2.074*zm)/dox
      alfi0=-1.172/(1.+0.462*x)
!  Collisionality coeficients from Hinton,Hazeltine
      a13=0.027*zm*zm-0.211*zm+1.204
      b13=0.14*zm*zm-0.87*zm+1.8
      c13=0.097*zm*zm-0.67*zm+1.643
      a23=0.01*zm*zm-0.08*zm+0.64
      b23=0.088*zm*zm-0.535*zm+1.057
      c23=0.06*zm*zm-0.41*zm+0.97
!  Collisionalities:
      eps0=8.8542e-12
      eq=1.602e-19
      emas=9.1094e-31
      epsfac=(eps0/eq)**2
!  Electron collision time...
      rt=sqrt(emas)*sqrt(te*eq)
      coolog=24.-log(sqrt(ne*1.0e13)/te)
      botcol=zm*zm*ne*bk*coolog
      topcol=rt*epsfac*te
      colte=(3.*2.*pi*sqrt(2.*pi))*topcol/botcol
      vthe=sqrt(2.*eq*te/emas)
      bfac=bmax*bdif*((bmax+bmin)/(bmax-bmin))**2
      tnue=zeff*bfac/(4.*colte*vthe)
      cnue=zeff*sqrt(2.)*nor*bmax/(2.*pi*vthe*colte)
!  Ion collision time...
      rt=sqrt(ti*eq)
      botcol=zm*zm*zm*zm*ne*bk*coolog
      topcol=rt*epsfac*ti
      colti=(3.*4.*pi*sqrt(pi))*topcol/botcol
      vthi=sqrt(2.*eq*ti)
      tnui=zeff*bfac/(4.*colti*vthi)
      cnui=zeff*sqrt(2.)*nor*bmax/(2.*pi*vthi*colti)
      write(6,*)' tnue=',tnue,' cnue=',cnue,' tnui=',tnui, &
           ' cnui=',cnui
      f1=(1.+a13*sqrt(tnue)+b13*tnue)*(1.+c13*cnue)
      f2=(1.+a23*sqrt(tnue)+b23*tnue)*(1.+c23*cnue)
      rl31=rl310/f1
      rl32=(rl320+2.5*rl310)/f2-2.5*rl31
      fa=(1.+cnui**2)*(1.+cnue**2)
      alfi=((alfi0+0.35*sqrt(tnui))/(1.+0.7*sqrt(tnui)) &
           +2.1*cnui**2)/fa
      bstrap=-(rl31*(a1e+(a1i+alfi*a2i)/(tau*zm))+rl32*a2e)
!  Spitzer conductivity
      epsfac=(4.*pi*eps0/eq)**2
      pro=ne*bk
      rt=sqrt(emas)*sqrt(eq*te)
      botcol=sqrt(2.*pi)*pro*coolog
      topcol=rt*epsfac*te
      colt=(3./4.)*topcol/botcol
      fac=1.98*(ne*bk)*(eq/emas)
      sigsp=fac*colt
      rnust=sqrt(2.)*(rmajor*aspect**1.5)*qsurf/(colt*vthe)
!  Neoclassical enhancement
      call signeo(rnust,ftrap,zeff,sigfac)
      coneo=sigfac*sigsp/1.98
      indj=coneo*vloop*f*riav/(2.*pi)
      ffprime_loc=-mu0*f*(bstrap+indj)/bsqav-f*f*pprime_loc/bsqav
      t1=-mu0*f*(bstrap)/bsqav
      t2=-mu0*f*(indj)/bsqav
      t3=-f*f*pprime_loc/bsqav
      write(6,*) 't1= ',t1,' t2=',t2,' t3=',t3
      write(6,*)' b/strap=',bstrap/bsqav,' induc=',indj/bsqav
      return
      
end subroutine hirsh





