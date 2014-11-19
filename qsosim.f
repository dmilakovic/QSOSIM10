!=======================================================================
      subroutine qsosim(zqso,alpha,rmag,wstart,pw,
     +         sigblur,npix,nl,nhi4,z4,loglam,flux,
     +         noise,nnflux,flux_nc)
!=======================================================================
      !generate the artificial spectrum
      !input, real*8              : zqso, alpha, rmag, wstart, pw, sigblur
      !input, integer             : npix, nl
      !input, real*8 array(nl)    : nhi4, z4
      !output, real*8 array(npix) : loglam, flux, noise, nnflux, flux_nc
!=======================================================================
      implicit none
      integer :: i,j,npix,idum,numpix,nl
      integer,parameter :: nems=63
      real*8 :: wstart,wend,pw,zqso,zqsop1,const,rmag,f6182,alpha,nuplim
      real*8 :: sigblur,dark,gain,ronsq
      real*8 :: vlight,pi,h,g,c,d,t
      real*4 :: xmin,xmax,ymin,ymax,x,m,area
      real*4,dimension(npix)::lambda,en,qe,signal,pix,sky,ns,
     +                                 noabs,conv,snr
      real*8,dimension(npix)::loglam,flux,flux_nc,nnflux,noise
      real*8 :: nhi4(4000),z4(4000), throughx(26), throughy(26)
      real*8 :: lognhi,nhi,z,b,voigt
      real*8 :: wems(nems), relstr(nems), sigma(nems)
      character :: ion(nems)*2
      data pi/3.14159265/
      data vlight/299792.458/
      data h/6.62607e-27/       !ergs*s
      data throughx/3650.0,3850.0,4120.0,4400.0,4900.0,5230.0,5500.0,
     +     6000.0,6150.0,6500.0,7000.0,7130.0,7500.0,7550.0,
     +     7650.0,7900.0,8200.0,8500.0,8950.0,9050.0,9000.0,
     +     9300.0,9500.0,9750.0,10000.0,10250.0/
      data throughy/0.05,0.1,0.15,0.20,0.25,0.27,0.26,
     +     0.234,0.238,0.253,0.27,0.29,0.283,0.09,
     +     0.27,0.255,0.212,0.229,0.233,0.206,0.22,
     +     0.13,0.148,0.135,0.09,0.05/
      real*8 gasdev3,ran3
      external gasdev3, ran3, blur


c parameters
      wend=wstart+npix*pw
c underlying continuum
      m = 3631*1e-23*((vlight*1e3)/(0.6182e-6)**2)*1e-10 ! 10e-23 erg s-1 cm-2 Hz-1 --> erg s-1 cm-2 A-1
      m = alog10(m)/0.4
      f6182=10**(-(rmag-m)/2.5)
      const=f6182/(1.0/6182.0**(2.0+alpha))*1e17
      do i=1,npix
         loglam(i)=wstart+(i-1)*pw
         lambda(i)=10**loglam(i)
         flux(i)=const*(1.0/(lambda(i)**(alpha+2)))
      end do
c read emission lines
      open(unit=12,file='emission.dat',err=101)
      goto 102
 101  write (6,*) 'error opening file'
      stop
 102  do i=1,nems
         read(12,*) wems(i),relstr(i),sigma(i),ion(i)
      end do
      close(unit=12)
c add emission lines to flux
      zqsop1=zqso+1.
      do j=1,nems
         wems(j)=wems(j)*zqsop1
         relstr(j)=relstr(j)/100.
         do i=1,npix
            x=1.75*relstr(j)*const
     +           *(1.0/wems(j)**(2+alpha))
            g=x*exp(-.5*(((lambda(i)-wems(j))
     +           /(sigma(j)))**2))
            flux(i)=flux(i)+g
         end do
      end do
c keep unabsorbed spectrum
      do i=1,npix
         noabs(i)=real(flux(i))
      end do
c initialise random numbers
      idum=time()
c input absorption
      do i=1,nl
         lognhi = nhi4(i)
         nhi = 10**lognhi
         z = z4(i)
         b = 5*gasdev3(idum)+20
         call spvoigt(flux,lambda,npix,nhi,z,b,'H ','I   ')
      end do
c save uncolvolved flux real*8
      do i=1,npix
         flux_nc(i)=flux(i)
      end do
c instrumental blurring
c      call blur(flux,npix,sigblur)
c save convolved no-noise flux real*8
      do i=1,npix
         nnflux(i)=flux(i)
      end do
c add noise : 
      c=0
      numpix=npix*3                                                       !total number of pixels, Smee et al: spectrum width ~ 3 pix
      t=60*60                                                             !integration time [sec]
      area=pi*250**2                                                      !telescope collecting area [cm2]
      do i=1,npix
         pix(i)=(1./vlight)*69*lambda(i)                                  !pixel size [A], Bolton et al. 2012. : BOSS 69 km/s
         c = c + pix(i)
         d = gasdev3(idum)
         en(i)=vlight*h/((lambda(i))*1e-10)*1e17                          !incident photon energy [1e-17 ergs]
         qe(i)=0.20*exp(-.5*(((lambda(i)-5080.0)/8.1e2)**2))
     +        +0.17*exp(-.5*(((lambda(i)-7200.0)/7e2)**2))
     +        +0.15*exp(-.5*(((lambda(i)-9000.0)/8e2)**2))
     +        +0.10*exp(-.5*(((lambda(i)-7000.0)/20e2)**2))
     +        -0.19*exp(-.5*(((lambda(i)-7550.0)/20)**2))
         signal(i)=nnflux(i)*t*pix(i)*area/(en(i)*qe(i))                    !signal in incident photons
         sky(i)=abs(gasdev3(idum)*10**6.1)                                !sky signal in incident photons  (previously: 1e6)
!  ***   use different values for red and blue arms of the spectrographs
         if((lambda(i).ge.3600.0).and.(lambda(i).lt.6050.))then
            dark=(0.525+0.022*gasdev3(idum))/900                              !e-/pixel/s
            gain=1.02+0.01*gasdev3(idum)                                  !e-/ADU
         else
            dark=(1.065+0.161*gasdev3(idum))/900                              !e-/pixel/s
            gain=1.70+0.07*gasdev3(idum)                                  !e-/ADU
         end if
         ronsq=(2.25+0.2*gasdev3(idum))**2                                ! read-out**2 noise (e-) (per read-out!, interval = 55sec)
         ns(i)=d*sqrt(signal(i)+(sky(i)+dark*t)*numpix+ronsq*(t/55))        ! total noise in e-
         noise(i)=ns(i)*en(i)/qe(i)/(t*pix(i)*area)                       ! noise in erg s-1 cm-1 A-1
         flux(i) = nnflux(i)+noise(i)
         snr(i)=abs(nnflux(i)/noise(i))
c         write(6,'(3f12.3)') nnflux(i),noise(i),snr(i)
      end do
c plot      
      xmin=real(wstart)
      xmax=real(wend)
      ymin=0.0
      ymax=0.0
      do i=1,npix
c         write (6,*) i, flux(i)
         if (flux(i).gt.ymax) ymax=flux(i)
         if (flux(i).lt.ymin) ymin=flux(i)
      end do
      ymax=ymax+0.2*ymax
      ymin=ymin-0.2*ymin
c      call PGBEGIN(0,'/xserve',1,1)
c      call PGENV(xmin,xmax,ymin,ymax,0,1)
c      call PGLINE(npix,real(loglam),real(flux))
c      call PGEND
      
      end subroutine qsosim




















c Numercial Recipes subroutine
!======================================================================
      FUNCTION ran3(idum)
!======================================================================
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL*8 ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END

c Numercial Recipes subroutine
!======================================================================
      FUNCTION gasdev3(idum)
!======================================================================
      INTEGER idum
      REAL*8 gasdev3
c Uses ran3
      INTEGER iset
      REAL*8 ran3
      REAL*4 fac,gset,rsq,v1,v2
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran3(idum)-1.
        v2=2.*ran3(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev3=v2*fac
        iset=1
      else
        gasdev3=gset
        iset=0
      endif
      return
      END
!======================================================================
c BLUR does Gaussian blurring on array xbuf
	subroutine blur(xbuf,npt,sigma)
!======================================================================
	implicit none
	integer nfilt, npt, i, il, ilow, k
	real*8 xbuf(npt),ybuf(npt),work(512), sum
	real*8 xmns, xmnf, sigma, const, norm

	nfilt=int(10.0*sigma)+1
	if(nfilt.gt.511)stop ' too large a filter'
	if(npt.gt.262144)stop ' too many points in data array'
	if((nfilt/2)*2.eq.nfilt)nfilt=nfilt+1
c *** fill up blurring array
c	const=1.D0/(sqrt(8.D0*atan(1.D0))*sigma)
	const=1.0
	do i=1,nfilt
	   work(i)=const*exp(-0.5D0*(dble(i)-
     &	(dble(nfilt)+1.D0)/2.D0)**2.D0/sigma**2.D0)
	enddo
c *** set first and last edges equal
	il=nfilt/2
	ilow=max0(3,nfilt/4)
	ilow=(ilow/2)*2+1
	sum=0.D0
	do i=1,ilow
 	   sum=sum+xbuf(i)
	enddo
	xmns=sum/dble(ilow)
        sum=0.D0
	do i=1,ilow
 	   sum=sum+xbuf(npt+1-i)
	enddo
	xmnf=sum/dble(ilow)
c *** reflect edges before filtering
	do i=1,il
	   ybuf(i)=2.D0*xmns-xbuf(il+ilow+1-i)
 	   ybuf(npt+i+il)=2.D0*xmnf-xbuf(npt-i-ilow+1)
	enddo
	do i=1,npt
	   ybuf(i+il)=xbuf(i)
	enddo
c *** do blurring
	do k=1,npt
	   sum=0.D0
	   norm=0.D0
	   do i=1,nfilt
	      sum=sum+work(i)*ybuf(i+k-1)
	      norm=norm+work(i)
	   enddo
	   xbuf(k)=sum/norm
	enddo
	return
	end
