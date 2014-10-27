c ----------------------------------------------------------------------------
c QSOSIM9. Prog. to make synthetic QSO spectrum - John Webb, UNSW, Dec 2013
c Uses VPFIT Voigt profile programs
c Compile using Makefile
c ----------------------------------------------------------------------------
	SUBROUTINE QSOSIM9(zqso,alpha,vmag,wstart,wend,dw,nc,nuplim,
     :          sigblur,s2n,inoise,dvavoid,npts,loglam,flux,
     :          flerr,nnflux,npoints,numlin,ni,nhi4,z4)
	  real*4 wda(262144),danoabs(262144),tau(262144)
	  real*4 da4(262144),da4conv(262144)
	  real*4 da_err4(262144),da4smno(262144)
	  real*8 da(262144),da_err(262144),da_err4mod(262144)
	  real*4 wems(62),relstr(62),sigma(62)
	  character*10 ion(62),med
	  character*20 cMed
	  real*8 nhi,b,z,z1p1,z2p1
	  real*4 w,ff,sum,zstart,zend
	  real*4 a12p5,rn,a13p75,a13p1,mbp1,zqsop1,pi
	  real*4 c,d,p,q,r,s,vlight,zleft,zright
	  real*8 lognhi, delta,mdelta
	  real*4 nhills(20),blls(20),zlls(20)
	  real*8, dimension(262144) :: en,ns,noise,signal,snr,qe,sky
	  real*8 :: zqso,alpha,vmag,wstart,wend,dw
	  real*8 :: s2n,dvavoid,nc,nuplim,sigblur
	  real*8, dimension(262144) :: loglam,lambda,flux
	  real*8, intent(out) :: flerr(262144), nnflux(262144)
          real*8 alm, fik,asm,ran3
	  real*8 npix,dark,gain,ron2,sig,t,pixsize,area
	  real*4 n1,n2,n3,median
	  real*4, dimension(3) :: bins
	  real*8, dimension(3) :: gamma
	  integer,dimension(3) :: ni
	  real*8, dimension(npoints) :: dummy, dummyy
	  real*8,dimension(numlin) :: nhi4,z4
	  integer npts,npoints,idum,numlls,iflag,inoise,index,numlin,mid

	  data pi/3.14159265/
	  data vlight/299792.458/
	  data h/6.62607e-27/   !ergs*s
	  data gamma/1.51,2.16,1.33/
	  data bins/12.0,14.0,17.0/
      external gasdev3
      common/vpc_ewllns/lbz(2500),lzz(2500),alm(2500),
     :                fik(2500),asm(2500),numel
c ----------------------------------------------------------------------------
c Emission line wavelengths
c	  data wems/1026.0,1216.0,1302.0,1335.0,1400.0,1549.0,1640.0,
c     +            1858.0,2000.0,2080.0,2140.0,2175.0,2225.0,2326.0,
c     +            2423.0,2798.0,2970.0,3130.0,3200.0,3346.0,3426.0,
c     +            3727.0,3869.0,3968.0,4068.0,4340.0,4861.0,4959.0,
c     +            5007.0,1240.0/
c Approx. relative emission line strengths - from eyeballing a few spectra
c	  data relstr/0.093,1.0,0.035,0.025,0.19,0.63,0.18,0.29,0.0049,
c     +              0.041,0.0034,0.0076,0.0047,0.06,0.022,0.34,
c     +              0.063,0.0073,0.0095,0.0052,0.01,0.0078,0.036,
c     +              0.013,0.028,0.13,0.22,0.0093,0.034,0.4/

c Max no. of pixels, hard-coded in array declarations
	  write (6,*) '--------------------------------------------------'
	  write (6,*) 'Started QSOSIM9!'
	  nptsmax=262144
	  nemlines=63
      open(unit=12,file='emission.dat',err=101)
	goto 102
 101	  write(6,*)' Error - no emission.dat, or format wrong'
	  stop
 102	  do i=1,nemlines
	     read(12,*) wems(i),relstr(i),sigma(i),ion(i)
	  end do
      close(unit=12)
c alpha is QSO spectral index, vmag is the V magnitude
c wstart, wend, dw are the start and end wavelengths and pixel size
c       nc is the column density cut-off (REAL, NOT log10)
	  zqsop1=zqso+1
	  f5550=10**(-(vmag+21.17)/2.5)
	  const=f5550/(1.0/5550.0**(2+alpha))
      sum=0.0
      npts=(wend-wstart)/dw
      if(npts.gt.nptsmax)then
      	write(6,*)' Max. no. of pixels = ',nptsmax
      	stop
      end if
	  do l=1,npts
c Make wavelength scale and underlying power-law QSO continuum
	     wda(l)=wstart+((l-1)*dw)
	     da(l)=const*(1.0/wda(l)**(2+alpha))
	  end do

c Put emission lines in. Guess/approximation for emission line width:
	  fwhm=60.0
	  sigma=fwhm/2.35

	  do m=1,nemlines
        wems(m)=wems(m)*zqsop1
	relstr(m)=relstr(m)/100.
		 do i=1,npts
           x=0.75*relstr(m)*const*(1.0/wems(m)**(2+alpha))
           g=x*exp(-.5*(((wda(i)-wems(m))/sigma(m))**2))
           da(i)=da(i)+g
         end do
      end do
	  write (6,*) 'Emission lines inputted!'
c Keep unabsorbed QSO spectrum
	  do i=1,npts
	    danoabs(i)=da(i)
	  end do

c Generate optical depth array in preparation for absorption input
	  do i=1,npts
	    tau(i)=0.0
	  end do
         write (*,*) 'Inputting absorption lines ...'
c Next section makes absorption lines
c g is from dn=A(1+z)^g dz, beta is from dn propto N^{-beta} dN.
c a13p75 is the value of A for lines above logN=13.75.
c gp1= gamma+1, mbp1=1-beta, nc is N cutoff, n is total no. of lines
      a13p75=10.0
      zstart=wstart/1215.67-1.
      zend=zqso
      z1p1=(wstart/1215.67)
      z2p1=zqsop1

c Calculate the total no. of lines
      n=numlin
	
c Initialise random numbers
      idum=time()
      do i=1,npoints
	 dummy(i)=dble(i)
      end do

c Call Voigt profile generator n times, once for each abs system
	  do i=1,numlin
	     
c Random selection of NHI (must have appropriate number of points in a defined
c range of NHI) and redshifts
	     nhi=10**nhi4(i)
	     z=z4(i)
! --------------------------------------------------------------------

c b-params.  Guess at sigma and mean of b distribution of 3 and 23 km/s.
        b = 5*gasdev3(idum)+20

c Put the forest lines in, but only if:
c 1. it is not if within the specified avoidance zone of each LLS, and
c 2. if its N(HI) is less than the specified upper limit for forest lines
      iflag=0
c      do j=1,numlls
c      zright=zlls(j) + dvavoid*(1.0+zlls(j))/vlight
c      zleft=zlls(j) - dvavoid*(1.0+zlls(j))/vlight
c      if(z.ge.zleft.and.z.le.zright)iflag=1
c      end do
      if(nhi.ge.nuplim)iflag=1
      if(iflag.eq.0)call spvoigt(da,wda,npts,nhi,z,b,'H ','I   ')

      end do
c End of loop for forest.  Now put the LLS's in
c      do j=1,numlls
c        call spvoigt (da,wda,npts,dble(nhills(j)),dble(zlls(j)),
c     :                dble(blls(j)),'H ','I   ')
c      end do

c Keep unconvolved real*4 spectrum
      do i=1,npts
        da4(i)=real(da(i))
      end do
Convolution with assumed Gaussian instrumental profile
	   write(6,*) "Convolving..."
	   call blur(da, npts, sigblur)
	   write (6,*) "Bluring complete..."
c Make real*4 array for pgplot
      do i=1,npts
        da4conv(i)=real(da(i))
      end do

c Make error array
      do i=1,npts
        da_err(i) = da(i)/dble(s2n) + 0.2*dble(danoabs(i))/dble(s2n)
        da_err4(i) = real(da_err(i))
      end do

c Add noise to real*4 convolved spectrum. 2 noise models.
c inoise=0 is constant.  inoise=1 gets worse towards blue. See qsosim9.pdf.
	npix=8000
	t=15*60
	pixsize=(wend-wstart)/npix
	area=pi*250**2
      if(inoise.eq.0)then
      do i=1,npts
        da_err4mod(i) = gasdev3(idum)*da_err4(i)
        da4smno(i) = da4conv(i) + da_err4mod(i)
      end do
      end if
      if(inoise.eq.1)then
      do i=1,npts
	 c=gasdev3(idum)
	 sigma=60.0
	 en(i)=vlight*h/(wda(i)*1e-10)          !incident photon energy [ergs]
	 qe(i)=0.20*exp(-.5*(((wda(i)-5000.0)/6e2)**2))
     +        +0.22*exp(-.5*(((wda(i)-7000.0)/7e2)**2))
     +        +0.22*exp(-.5*(((wda(i)-9000.0)/8e2)**2))
     +        +0.10*exp(-.5*(((wda(i)-6000.0)/9e2)**2))
     +        +0.02*alog10(wda(i)-3.1)
	 signal(i)=da4conv(i)*t*pixsize*area/en(i) !e-
	 sky(i)=abs(0.1*gasdev3(idum)*signal(i))   !e-
	 !use different values for red and blue arms of the spectrographs
	 if((wda(i).ge.3650).and.(wda(i).lt.6050))then
	    dark=(0.525+0.022*gasdev3(idum))        !e-/pixel/15min
	    gain=1.02+0.01*gasdev3(idum)            !e-/ADU
	 else
	    dark=(1.065+0.161*gasdev3(idum))        !e-/pixel/15min
	    gain=1.70+0.07*gasdev3(idum)            !e-/ADU
	 end if
	 !square of the read-out noise (e-)
	 ron2=(2.25+0.2*gasdev3(idum))**2
	 ! ns - noise in e-
        ns(i)=sign(1.0,c)*
     :         sqrt(signal(i)+(sky(i)+dark)*npix+ron2*(900/55))/qe(i)
	 !noise in ergs/s/cm^2/A
	 noise(i)=ns(i)*en(i)/(t*pixsize*area)
        da4smno(i) = da4conv(i)+noise(i)
	snr(i)=abs(da4conv(i)/noise(i))
	write(6,'(3d12.3)') da4smno(i),noise(i),snr(i)
	s=s+snr(i)
      end do
      end if
     
c Plot spectrum
      call PGBEGIN (0,'/cps',1,3)
      xmin=alog10(real(wstart))
      xmax=alog10(real(wend))
      ymin=0.0
      ymax=0.0
      do i=1,npts
        if(da4smno(i).gt.ymax)ymax=da4smno(i)
      end do
      ymin=ymax
      do i=1,npts
        if(da4smno(i).lt.ymin)ymin=da4smno(i)
      end do
      ymax=ymax+0.5*ymax
      ymin=ymin-0.5*ymin
	ypost=ymin+(ymax-ymin)*0.91
	ypos1=ymin+(ymax-ymin)*0.90
	ypos2=ymin+(ymax-ymin)*0.82
c Data to be returned into main program
	
 100	format(f12.5,2x,f12.5,2x,f12.5,2x,f12.5)
      do i=1,npts
	 lambda(i)=wda(i)
	 loglam(i)=alog10(wda(i))
	 flux(i)=da4smno(i)
	 flerr(i)=da_err4(i)
	 nnflux(i)=da4conv(i)
c         write(*,*)wda(i),flux(i),noise(i),nnflux(i)
      end do


      call PGENV (xmin,xmax,ymin,ymax,0,1)
      call PGLABEL ('Wavelength','f(lambda)','Log wavelengths')
      call pgline(npts,real(loglam),da4smno)
	do i=1,nemlines
	   call pgptxt(alog10(wems(i)),ypost,90.0,0.0,ion(i))
	   call pgmove(alog10(wems(i)),ypos1)
	   call pgdraw(alog10(wems(i)),ypos2)
	end do
      call pgslw(1)
      call pgsci(5)
      call pgline(npts,wda,da4conv)
      call pgsls(2)
      call pgsci(2)
      call pgline(npts,wda,danoabs)
      call pgsls(1)
      call pgsci(3)
      call pgline(npts,wda,da_err4)
	
	ymin=0.0
	ymax=0.0
	do i=1,npts
	   if(noise(i).gt.ymax)ymax=real(noise(i))
	end do
	ymin=ymax
	do i=1,npts
	   if(noise(i).lt.ymin)ymin=real(noise(i))
	end do
	ymax=ymax+0.5*ymax
	ymin=ymin+0.5*ymin
	call pgsci(1)
	call PGENV(real(wstart),real(wend),ymin,ymax,0,1)
	call PGLABEL('Wavelength','Noise','Signal-to-noise')
c	call PGPT(numlin,real(nhi4),real(z4),3)
	call pgline(npts,wda,real(noise))
c	call PGENV(0.0,real(npts),0.0,50.0,0,1)
c	call PGLABEL('log NHI','No of lines','Choice of NHI')
c	call PGBIN(3,bins,real(ni),.FALSE.)
	call PGHIST(npts,real(snr),0.0,400.0,100,0)
	call PGLABEL('S/N','#','Signal-to-noise histogram')
	if (mod(npts,2).eq.0) then
	mid=npts/2
	else 
	   mid=nint(npts/2+0.5)
	end if
	median=select(mid,npts,real(snr))
	write(med,'(f7.3)') median
	cMed='S/N median = '//med
	write (6,*) median
	call pgtext(60.0,800.0,cMed)
      call pgend

      return
      end

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
      REAL gasdev3
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
	1	(dble(nfilt)+1.D0)/2.D0)**2.D0/sigma**2.D0)
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
!======================================================================
	FUNCTION select(k,n,arr)
	INTEGER k,n
	REAL select,arr(n)
	!Returns the kth smallest value in the array arr(1:n). The input array will be rearranged to have this value in location arr(k), with all smaller elements moved to arr(1:k-1) (in arbitrary order) and all larger elements in arr[k+1..n] (also in arbitrary order).
!======================================================================
	INTEGER i,ir,j,l,mid
	REAL a,temp
	l=1
	ir=n
 1	if(ir-l.le.1)then 
	if(ir-l.eq.1)then
	if(arr(ir).lt.arr(l))then
	   temp=arr(l)
	   arr(l)=arr(ir)
	   arr(ir)=temp
	endif
	endif
	select=arr(k)
	return
	else
	   mid=(l+ir)/2 
	   temp=arr(mid)
	   arr(mid)=arr(l+1)
	   arr(l+1)=temp
	   if(arr(l).gt.arr(ir))then
	      temp=arr(l)
	      arr(l)=arr(ir)
	      arr(ir)=temp
	   endif
	   if(arr(l+1).gt.arr(ir))then
	      temp=arr(l+1)
	      arr(l+1)=arr(ir)
	      arr(ir)=temp
	   endif
	   if(arr(l).gt.arr(l+1))then
	      temp=arr(l)
	      arr(l)=arr(l+1)
	      arr(l+1)=temp
	   endif
	   i=l+1 
	   j=ir
	   a=arr(l+1)
 3	   continue 
	   i=i+1 
	   if(arr(i).lt.a)goto 3
 4	   continue
	   j=j-1
	   if(arr(j).gt.a)goto 4
	   if(j.lt.i)goto 5 
	   temp=arr(i)
	   arr(i)=arr(j)
	   arr(j)=temp
	   goto 3 
 5	   arr(l+1)=arr(j) 
	   arr(j)=a
	   if(j.ge.k)ir=j-1 
	   if(j.le.k)l=i
	endif
	goto 1
	END FUNCTION
