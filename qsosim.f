!=======================================================================
      subroutine qsosim(descriptor,zqso,alpha,rmag,wstart,pw,
     +         npix,nl,nhi4,z4,mask,novi4,loglam,flux,
     +         noise,nnflux,ncflux,noabs)
!=======================================================================
*      generate the artificial spectrum
      
*      input, real*8              : zqso         QSO redshift
*                                   alpha        spectral index
*                                   rmag         r magnitude
*                                   wstart       starting wavelength
*                                   pw           pixel width
*                                   sigblur      FWHM of the gaussian profile                 
*      input, integer             : npix         number of pixels
*                                   nl           number of H I lines
*      input, real*8 array(nl)    : nhi4         H I col densities
*                                   z4           H I redshifts
*      input, logical array(nl)   : mask         mask of H I lines with N_HI > 1e15 cm-2
*                                                (T if N_HI > 1e15 cm-2)

*      output, real*8 array(npix) : loglam       wavelengths in log
*                                   flux         flux in 10e-17 erg s-1 cm-2 A-1
*                                   noise        noise array
*                                   nnflux       no-noise flux
*                                   flux_nc      not convolved flux
*                                   noabs        flux with no absorption
!=======================================================================
      implicit none
      integer :: i,j,npix,idum,numpix,nl,s
      integer,parameter :: nems=63
      real*8 :: wstart,wend,pw,zqso,zqsop1,const,rmag,f6182,alpha,nuplim
      real*8 :: sigblur,dark,gain,ronsq
      real*8 :: vlight,pi,h,g,c,d,t
      real*4 :: xmin,xmax,ymin,ymax,x,m,area
      real*4,dimension(npix)::lambda,en,qe,signal,pix,sky,ns,
     +                                 conv,snr
      real*8,dimension(npix)::loglam,flux,ncflux,nnflux,noise,noabs
      real*8 ::  throughx(26), throughy(26)
      real*8,dimension(nl) :: nhi4,z4,novi4
      real*8 :: lognhi,nhi,novi,z,b,voigt
      real*8 :: wems(nems), relstr(nems), sigma(nems)
      logical,dimension(nl)::mask
      character :: ion(nems)*2, string*10, descriptor*15
      data pi/3.14159265/
      data vlight/299792.458/
      data h/6.62607e-27/             !units: ergs*s
*     throughput of SDSS BOSS, figure 38. in Smee et al. 2013.
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


* parameters
      wend=wstart+npix*pw
* underlying continuum
      m = 3631*1e-23*((vlight*1e3)/(0.6182e-6)**2)*1e-10 ! magnitude, conversion: 10e-23 erg s-1 cm-2 Hz-1 --> erg s-1 cm-2 A-1
      m = alog10(m)/0.4                                  
      f6182=10**(-(rmag-m)/2.5)
      const=f6182/(1.0/6182.0**(2.0+alpha))*1e17         !f_lambda = const * lambda**-(2+alpha) [1e17 erg s-1 cm-2 A-1]
      do i=1,npix
         loglam(i)=wstart+(i-1)*pw
         lambda(i)=10**loglam(i)
         flux(i)=const*(1.0/(lambda(i)**(alpha+2)))
      end do
* read emission lines from file
      open(unit=12,file='emission.dat',err=101)
      goto 102
 101  write (6,*) 'error opening file'
      stop
 102  do i=1,nems
         read(12,*) wems(i),relstr(i),sigma(i),ion(i)
      end do
      close(unit=12)
* add emission lines to flux
* modelled as gaussians with mean=wems(i) and sigma=sigma(i)
      zqsop1=zqso+1.
      do j=1,nems
         wems(j)=wems(j)*zqsop1
         relstr(j)=relstr(j)/100.
         do i=1,npix
            x=2.0*relstr(j)*const            !1.75*
     +           *(1.0/wems(j)**(2+alpha))
            g=x*exp(-.5*(((lambda(i)-wems(j))
     +           /(sigma(j)))**2))
            flux(i)=flux(i)+g
         end do
      end do
* keep unabsorbed spectrum
      do i=1,npix
         noabs(i)=flux(i)
      end do
* initialise random numbers
      idum=time()
* input H I absorption
      do i=1,nl
         lognhi = nhi4(i)
         nhi = 10**lognhi
         z = z4(i)
         b = 5*gasdev3(idum)+20
         call spvoigt(flux,lambda,npix,nhi,z,b,'H ','I   ')
      end do
* input O VI absorption
c      do j=1,nl
c            write (6,*) 'j,mask,nhi,novi',j,mask(j),nhi4(j),novi4(j)
c         end do
      do i=1,nl
         if (mask(i).eqv..true.) then
            novi = novi4(i)
            z = z4(i)
            b = 3*gasdev3(idum)+12
            call spvoigt(flux,lambda,npix,novi,z,b,'O ','VI  ')
         end if
      end do
* save uncolvolved flux real*8
      do i=1,npix
         ncflux(i)=flux(i)
      end do
* add instrumental blurring:
* BOSS resolution R~2000, sigma_v=vlight/(2.35*R)=64 km/s
* using log10(sigma_v) since the wavelength scale is log
      call blur(flux,npix,dble(alog10(64.0)))                             
                                                                          
* save convolved no-noise flux real*8
      do i=1,npix
         nnflux(i)=flux(i)
      end do
* add BOSS noise : 
      c=0
      numpix=npix*3                                                       !total number of pixels, Smee et al: spectrum width ~ 3 pix
      t=60*60                                                             !integration time, usually 4x15 min = 60x60 s [sec]
      area=pi*250**2                                                      !telescope collecting area [cm2], radius=250 cm
      do i=1,npix
         pix(i)=(1./vlight)*69*lambda(i)                                  !pixel size [A], Bolton et al. 2012. : BOSS pixels constant in velocity space: 69 km/s
         c = c + pix(i)
         d = gasdev3(idum)
         en(i)=vlight*h/((lambda(i))*1e-10)*1e17                          !incident photon energy [1e-17 ergs]
         qe(i)=0.20*exp(-.5*(((lambda(i)-5080.0)/8.1e2)**2))              !throughput model: 5 gaussians
     +        +0.17*exp(-.5*(((lambda(i)-7200.0)/7e2)**2))                !modelled after fig 38 in Smee et al. 2012
     +        +0.15*exp(-.5*(((lambda(i)-9000.0)/8e2)**2))
     +        +0.10*exp(-.5*(((lambda(i)-7000.0)/20e2)**2))
     +        -0.19*exp(-.5*(((lambda(i)-7550.0)/20)**2))
         signal(i)=nnflux(i)*t*pix(i)*area/(en(i)*qe(i))                  !signal in incident photons
         sky(i)=abs(gasdev3(idum)*10**6.2)                                !background signal in incident photons  (previously: 1e6)
!  ***   use different values for red and blue arms of the spectrographs
         if((lambda(i).ge.3600.0).and.(lambda(i).lt.6050.))then
            dark=(0.525+0.022*gasdev3(idum))/900                              !e-/pixel/s
            gain=1.02+0.01*gasdev3(idum)                                  !e-/ADU
         else
            dark=(1.065+0.161*gasdev3(idum))/900                              !e-/pixel/s
            gain=1.70+0.07*gasdev3(idum)                                  !e-/ADU
         end if
         ronsq=(2.25+0.2*gasdev3(idum))**2                                ! read-out**2 noise (e-) (per read-out!, interval = 55sec)
c         d = sign(1d0,d)
         ns(i)=d*sqrt(signal(i)+(sky(i)+dark*t)*numpix+ronsq*(t/55))        ! total noise in e-
         noise(i)=ns(i)*en(i)/qe(i)/(t*pix(i)*area)                        ! noise in erg s-1 cm-1 A-1
         flux(i) = nnflux(i)+noise(i)
         snr(i)=abs(nnflux(i)/noise(i))
c         write(6,'(3f12.3)') nnflux(i),noise(i),snr(i)
      end do
c plot      
      xmin=10**real(wstart)
      xmax=10**real(wend)
      ymin=0.0
      ymax=0.0
      do i=1,npix
c         write (6,*) i, flux(i)
         if (flux(i).gt.ymax) ymax=flux(i)
         if (flux(i).lt.ymin) ymin=flux(i)
      end do
      ymax=ymax+0.2*ymax
      ymin=ymin-0.2*ymin
      call PGBEGIN(0,'/null',1,1)
c      call PGBEGIN(0,'spectrum_sigma69.ps/cps',1,1)
      call PGSLW(3)
      call PGENV(xmin,xmax,ymin,ymax,0,0,1)
      call PGMTXT('T',1.0,0.5,0.5,descriptor)
      write(string,'(f5.2)') rmag
      call PGMTXT('T',1.0,1.0,1.0,'r='//trim(string))
      write(string,'(f8.6)') zqso
      call PGMTXT('T',1.0,0.0,0.0,'z='//trim(string))
      call PGSLW(1)
      call PGLINE(npix,real(lambda),real(flux))
      call PGSCI(3)
      call PGLINE(npix,real(lambda),real(noabs))
      call PGSCI(5)
      call PGLINE(npix,real(lambda),real(ncflux))
      call PGSCI(2); call PGSLW(3)
      call PGLINE(npix,real(lambda),real(nnflux))
      call PGSCI(1)
      call PGEND
      return
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
	real*8 xbuf(npt), sum
        real*8,dimension(:),allocatable :: work, ybuf
	real*8 xmns, xmnf, sigma, const, norm

	nfilt=int(10.0*sigma)+1
	if(nfilt.gt.811)stop ' too large a filter'
	if(npt.gt.262144)stop ' too many points in data array'
        allocate(work(nfilt))
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
        allocate(ybuf(npt+2*il))
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
c$$$
c$$$!======================================================================
c$$$      SUBROUTINE convlv(data,n,respns,m,isign,ans)
c$$$!======================================================================
c$$$      INTEGER isign,m,n,NMAX
c$$$      REAL data(n),respns(n)
c$$$      COMPLEX ans(n)
c$$$      PARAMETER (NMAX=4096) !Maximum anticipated size of FFT.
c$$$C     USES realft,twofft
c$$$C     Convolves or deconvolves a real data set data(1:n) (including any user-supplied zero
c$$$C     padding) with a response function respns, stored in wrap-around order in a real array of
c$$$C     length m < n. (m should be an odd integer.) Wrap-around order means that the rst half
c$$$C     of the array respns contains the impulse response function at positive times, while the
c$$$C     second half of the array contains the impulse response function at negative times, counting
c$$$C     down from the highest element respns(m). On input isign is +1 for convolution, −1
c$$$C     for deconvolution. The answer is returned in the rst n components of ans. However, ans
c$$$C     must be supplied in the calling program with length at least 2*n, for consistency with
c$$$C     twofft. n MUST be an integer power of two.
c$$$      INTEGER i,no2
c$$$      COMPLEX fft(NMAX)
c$$$      do i=1,(m-1)/2 !Put respns in array of length n.
c$$$         respns(n+1-i)=respns(m+1-i)
c$$$      end do
c$$$      do i=(m+3)/2,n-(m-1)/2 !Pad with zeros.
c$$$         respns(i)=0.0
c$$$      end do
c$$$      call twofft(data,respns,fft,ans,n)
c$$$      no2=n/2
c$$$      do i=1,no2+1
c$$$         if (isign.eq.1) then
c$$$            ans(i)=fft(i)*ans(i)/no2 !Multiply FFTs to convolve.
c$$$         else if (isign.eq.-1) then
c$$$            if (abs(ans(i)).eq.0.0) then 
c$$$               stop 'deconvolving at response zero in convlv'
c$$$            end if
c$$$            ans(i)=fft(i)/ans(i)/no2 !Divide FFTs to deconvolve.
c$$$         else
c$$$            stop 'no meaning for isign in convlv'
c$$$         endif
c$$$      end do
c$$$      ans(1)=cmplx(real(ans(1)),real(ans(no2+1))) !Pack last element with first for realft.
c$$$      call realft(ans,n,-1) !Inverse transform back to time domain.
c$$$      return
c$$$      END SUBROUTINE convlv
c$$$
c$$$!======================================================================
c$$$      SUBROUTINE twofft(data1,data2,fft1,fft2,n)
c$$$!======================================================================
c$$$      INTEGER n
c$$$      REAL data1(n),data2(n)
c$$$      COMPLEX fft1(n),fft2(n)
c$$$C     USES four1
c$$$C     Given two real input arrays data1(1:n) and data2(1:n), this routine calls four1 and
c$$$C     returns two complex output arrays, fft1(1:n) and fft2(1:n), each of complex length n
c$$$C     (i.e., real length 2*n), which contain the discrete Fourier transforms of the respective data
c$$$C     arrays. n MUST be an integer power of 2.
c$$$      INTEGER j,n2
c$$$      COMPLEX h1,h2,c1,c2
c$$$      c1=cmplx(0.5,0.0)
c$$$      c2=cmplx(0.0,-0.5)
c$$$      do j=1,n
c$$$         fft1(j)=cmplx(data1(j),data2(j)) !Pack the two real arrays into one complex array.
c$$$      end do                              
c$$$      call four1(fft1,n,1)                !Transform the complex array.
c$$$      fft2(1)=cmplx(aimag(fft1(1)),0.0)
c$$$      fft1(1)=cmplx(real(fft1(1)),0.0)
c$$$      n2=n+2
c$$$      do j=2,n/2+1
c$$$         h1=c1*(fft1(j)+conjg(fft1(n2-j))) !Use symmetries to separate the two trans
c$$$         h2=c2*(fft1(j)-conjg(fft1(n2-j)))    !forms.
c$$$         fft1(j)=h1                        !Ship them out in two complex arrays.
c$$$         fft1(n2-j)=conjg(h1)
c$$$         fft2(j)=h2
c$$$         fft2(n2-j)=conjg(h2)
c$$$      end do
c$$$      return
c$$$      END SUBROUTINE
c$$$
c$$$!======================================================================
c$$$      SUBROUTINE realft(data,n,isign)
c$$$!======================================================================
c$$$      INTEGER isign,n
c$$$      REAL data(n)
c$$$C     USES four1
c$$$C     Calculates the Fourier transform of a set of n real-valued data points. Replaces this data
c$$$C     (which is stored in array data(1:n)) by the positive frequency half of its complex Fourier
c$$$C     transform. The real-valued rst and last components of the complex transform are returned
c$$$C     as elements data(1) and data(2), respectively. n must be a power of 2. This routine
c$$$C     also calculates the inverse transform of a complex data array if it is the transform of real
c$$$C     data. (Result in this case must be multiplied by 2/n.)
c$$$      INTEGER i,i1,i2,i3,i4,n2p3
c$$$      REAL c1,c2,h1i,h1r,h2i,h2r,wis,wrs
c$$$      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
c$$$*     Double precision for the trigonometric recurrences.
c$$$      theta=3.141592653589793d0/dble(n/2)! Initialize the recurrence.
c$$$      c1=0.5
c$$$      if (isign.eq.1) then
c$$$         c2=-0.5
c$$$         call four1(data,n/2,+1) !The forward transform is here.
c$$$      else
c$$$         c2=0.5 !Otherwise set up for an inverse transform.
c$$$         theta=-theta
c$$$      endif
c$$$      wpr=-2.0d0*sin(0.5d0*theta)**2
c$$$      wpi=sin(theta)
c$$$      wr=1.0d0+wpr
c$$$      wi=wpi
c$$$      n2p3=n+3
c$$$      do i=2,n/4 !Case i=1 done separately below.
c$$$         i1=2*i-1
c$$$         i2=i1+1
c$$$         i3=n2p3-i2
c$$$         i4=i3+1
c$$$         wrs=sngl(wr)
c$$$         wis=sngl(wi)
c$$$         h1r=c1*(data(i1)+data(i3)) !The two separate transforms are separated out of data.
c$$$         h1i=c1*(data(i2)-data(i4)) 
c$$$         h2r=-c2*(data(i2)+data(i4))
c$$$         h2i=c2*(data(i1)-data(i3))
c$$$         data(i1)=h1r+wrs*h2r-wis*h2i !Here they are recombined to form the true trans
c$$$         data(i2)=h1i+wrs*h2i+wis*h2r !form of the original real data.
c$$$         data(i3)=h1r-wrs*h2r+wis*h2i
c$$$         data(i4)=-h1i+wrs*h2i+wis*h2r
c$$$         wtemp=wr !The recurrence.
c$$$         wr=wr*wpr-wi*wpi+wr
c$$$         wi=wi*wpr+wtemp*wpi+wi
c$$$      end do
c$$$      if (isign.eq.1) then
c$$$         h1r=data(1)
c$$$         data(1)=h1r+data(2)
c$$$         data(2)=h1r-data(2) !Squeeze the rst and last data together to get
c$$$                             ! them all else within the original array.
c$$$         h1r=data(1)
c$$$         data(1)=c1*(h1r+data(2))
c$$$         data(2)=c1*(h1r-data(2))
c$$$         call four1(data,n/2,-1) !This is the inverse transformfor the case isign=-1.
c$$$      endif
c$$$      return
c$$$      END SUBROUTINE 
c$$$!======================================================================
c$$$      SUBROUTINE four1(data,nn,isign)
c$$$!======================================================================
c$$$      INTEGER isign,nn
c$$$      REAL data(2*nn)
c$$$C     Replaces data(1:2*nn) by its discrete Fourier transform, if isign is input as 1; or replaces
c$$$C     data(1:2*nn) by nn times its inverse discrete Fourier transform, if isign is input as −1.
c$$$C     data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn
c$$$C     MUST be an integer power of 2 (this is not checked for!).
c$$$      INTEGER i,istep,j,m,mmax,n
c$$$      REAL tempi,tempr
c$$$      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp !Double precision for the trigonometn=
c$$$                                                 !2*nn ric recurrences.
c$$$      j=1
c$$$      do i=1,n,2 !This is the bit-reversal section of the routine.
c$$$         if(j.gt.i)then
c$$$            tempr=data(j) !Exchange the two complex numbers.
c$$$            tempi=data(j+1)
c$$$            data(j)=data(i)
c$$$            data(j+1)=data(i+1)
c$$$            data(i)=tempr
c$$$            data(i+1)=tempi
c$$$         endif
c$$$         m=n/2
c$$$ 17      if ((m.ge.2).and.(j.gt.m)) then
c$$$            j=j-m
c$$$            m=m/2
c$$$            goto 17
c$$$         endif
c$$$         j=j+m
c$$$      end do
c$$$      mmax=2 !Here begins the Danielson-Lanczos section of the routine.
c$$$ 18   if (n.gt.mmax) then       !Outer loop executed log2 nn times.
c$$$         istep=2*mmax
c$$$         theta=6.28318530717959d0/(isign*mmax) !Initialize for the trigonometric recurrence.
c$$$         wpr=-2.d0*sin(0.5d0*theta)**2
c$$$         wpi=sin(theta)
c$$$         wr=1.d0
c$$$         wi=0.d0
c$$$         do m=1,mmax,2 !Here are the two nested inner loops.
c$$$            do i=m,n,istep
c$$$               j=i+mmax !This is the Danielson-Lanczos formula:
c$$$               tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
c$$$               tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
c$$$               data(j)=data(i)-tempr
c$$$               data(j+1)=data(i+1)-tempi
c$$$               data(i)=data(i)+tempr
c$$$               data(i+1)=data(i+1)+tempi
c$$$            end do
c$$$            wtemp=wr !Trigonometric recurrence.
c$$$            wr=wr*wpr-wi*wpi+wr
c$$$            wi=wi*wpr+wtemp*wpi+wi
c$$$         end do
c$$$         mmax=istep
c$$$         goto 18 !Not yet done.
c$$$      endif !All done.
c$$$      return
c$$$      END SUBROUTINE
