!=======================================================================
      subroutine qsosim(descriptor,zqso,alpha,rmag,wstart,pw,
     +         npix,nhil,nhi4,npnhi4,z4,znc4,loglam,flux,
     +         noise,ivar,nnflux,ncflux,noabs,nlin)
!=======================================================================
*      generate the artificial spectrum
      
*      input, real*8              : zqso         QSO redshift
*                                   alpha        spectral index
*                                   rmag         r magnitude
*                                   wstart       starting wavelength
*                                   pw           pixel width              
*      input, integer             : npix         number of pixels
*                                   nhil         number of H I lines
*      input, real*8 array(nl)    : nhi4         H I col densities
*                                   z4           H I redshifts

*      output, real*8 array(npix) : loglam       wavelengths in log
*                                   flux         flux in 10e-17 erg s-1 cm-2 A-1
*                                   noise        noise array
*                                   nnflux       no-noise flux
*                                   flux_nc      not convolved flux
*                                   noabs        flux with no absorption
*     output, integer             : nlin         number of lines in the spectrum
!=======================================================================
      implicit none
      integer :: i,j,npix,idum,numpix,nhil,ind,nlin
      integer :: view,log_wav,with_noise,with_conv
      integer,parameter :: nems=63
      real*8 :: wstart,wend,pw,zqso,zqsop1,const,rmag,f6182,alpha,nuplim
      real*8 :: sigblur,ronsq,sigma
      real*8 :: vlight,pi,planckh,g,d,t
      real*4 :: xmin,xmax,ymin,ymax,x,m0,area,xx,yy
      real*4,dimension(npix)::en,qe,signal,pix,sky,ns,
     +                                 conv,snr,deltaflux,plflux
      real*4,dimension(npix)::plotx,ploty,dark,gain
      real*8,dimension(npix)::loglam,lambda,flux,ncflux,nnflux,noise
      real*8,dimension(npix)::nmflux,noabs,mflux,ivar
      real*8,dimension(npix)::npflux,nccflux !no proximity effect flux, testing purposes
      real*8 ::  throughx(26), throughy(26)
      real*8,dimension(nhil) :: nhi4,z4,npnhi4,novi4,znc4,deltaz
      real*8 :: lognhi,nhi,novi,z,bpar,voigt,npnhi,znc,m
*     LINE LIST
*     define a new type of variable, containing the total line list: 
*             atom, ionisation, redshift, column density, b-parameter
      type linelist
         SEQUENCE
         character*2 atom
         character*4 ion
         real*8 colden
         real*8 rdf
         real*8 bpar
      end type
      type (linelist) :: llist(5000)
      type (linelist),allocatable :: line(:)
      common /linelist/llist
      real*8 colden
*     EMISSION LINES
      real*8 :: wems(nems), relstr(nems), sig(nems)
*     METALS
      LOGICAL,DIMENSION(:),allocatable::mask
      REAL*8,DIMENSION(30,20000) :: H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,
     :                      P,S,Cl,Ar,K,Ca,Sc,Ti,Va,Cr,Mn,Fe,Co,Ni,Cu,Zn
      common/metals/H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,
     :                      P,S,Cl,Ar,K,Ca,Sc,Ti,Va,Cr,Mn,Fe,Co,Ni,Cu,Zn

*     DATA
      real*8 throughput
      character :: spec(nems)*2, string*10, descriptor*15, name*35
      character :: atom*2,ion*4
      data pi/3.14159265/
      data vlight/299792.458/
      data planckh/6.62607e-27/             !units: ergs*s
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
      external gasdev3, ran3, blur,througput



* parameters
      wend=wstart+npix*pw
* ----------------------------------------------------------------------
*               P O W E R - L A W     C O N T I N U U M
* ----------------------------------------------------------------------
* underlying continuum
*     convert magnitude into physical flux
*     m0 - referent magnitude, conversion factor: 10e-23 erg s-1 cm-2 Hz-1 --> erg s-1 cm-2 A-1
      m0 = 3631*1e-23*((vlight*1e3)/(0.6182e-6)**2)*1e-10 
      m0 = alog10(m0)/0.4   
      f6182=10**(-(rmag-m0)/2.5)
      const=f6182/(1.0/6182.0**(2.0+alpha))*1e17         !f_lambda = const * lambda**-(2+alpha) [1e17 erg s-1 cm-2 A-1]
      do i=1,npix
         loglam(i)=wstart+(i-1)*pw
         lambda(i)=10**loglam(i)
         flux(i)=const*(1.0/(lambda(i)**(alpha+2)))
      end do
* save power-law
      plflux=flux

* ----------------------------------------------------------------------
*              Q S O    E M I S S I O N     L I N E     I N P U T
* ----------------------------------------------------------------------
* initialise random numbers
      idum=time()
* read emission lines from file
      open(unit=12,file='emission.dat',err=101)
      goto 102
 101  write (6,*) 'error opening file'
      stop
 102  do i=1,nems
         read(12,*) wems(i),relstr(i),sig(i),spec(i)
      end do
      close(unit=12)
* add emission lines to flux
* modelled as gaussians with mean=wems(i) and sigma=sigma(i)
      zqsop1=zqso+1.
      do j=1,nems
         wems(j)=wems(j)*zqsop1
         relstr(j)=relstr(j)/100.
         sigma=4*sig(j)+0.1*gasdev3(idum)*sig(j)
         do i=1,npix
            x=2.0*relstr(j)*const            !1.75*
     +           *(1.0/wems(j)**(2+alpha))
            g=x*exp(-.5*(((lambda(i)-wems(j))
     +           /(sigma))**2))
            flux(i)=flux(i)+g
         end do
      end do
*     keep unabsorbed spectrum
      noabs=flux
      npflux=flux
* ----------------------------------------------------------------------
*              A B S O R P T I O N     L I N E     I N P U T
* ----------------------------------------------------------------------
*     create a line list, using HI line list from 'assign'
*     'addsys' adds an HI line to the complete line list, and inputs 
*     metals if HI column density >= 10^17 cm^-2
*     nhi in linear, not log scale
      call rwllist(0)
      j=0
      do i=1,nhil
         j=j+1
         nhi = nhi4(i)
         z = z4(i)
         znc=znc4(i)
         npnhi = npnhi4(i)
         bpar = 5*gasdev3(idum)+22
         call addsys(j,nhi,z,bpar)
      end do

*     total number of lines = nlin
      nlin=j
*     save the line list into a new array that will be output and saved
      allocate(line(nlin))
      do i=1,nlin
         line(i)=llist(i)
      end do
*     prepare arrays for metal-only (mflux) and no-metal (nmflux) absorption 
      mflux = flux
      nmflux= flux
*     retrieve the complete line list and input lines using 'spvoigt'
      do i=1,nlin
         atom=llist(i)%atom
         ion=llist(i)%ion
         colden=llist(i)%colden
         z=llist(i)%rdf
         bpar=llist(i)%bpar
         call spvoigt(flux,lambda,npix,colden,z,bpar,atom,ion)
         if (atom.eq.'H '.and.ion.eq.'I   ') then
*           create spectrum with no metal absorption
            call spvoigt(nmflux,lambda,npix,colden,z,bpar,atom,ion)
         else if (atom.ne.'H ') then
*           create spectrum with metal only absorption
            call spvoigt(mflux,lambda,npix,colden,z,bpar,atom,ion)
         end if
      end do
*     no proximity effect
      j=0
      do i=1,nlin
         atom=llist(i)%atom
         ion=llist(i)%ion
         if (atom.eq.'H '.and.ion.eq.'I   ') then
            j=j+1
            z=llist(i)%rdf
            bpar=llist(i)%bpar
         end if
         colden=npnhi4(j)
         
c         call spvoigt(npflux,lambda,npix,colden,z,bpar,'H ','I   ')
         
      end do
* ----------------------------------------------------------------------
* save uncolvolved flux real*8
      ncflux=flux
      nccflux=npflux
* ----------------------------------------------------------------------
*             I N S T R U M E N T A L    B L U R R I N G
* ----------------------------------------------------------------------
* BOSS resolution R~2000, sigma_v=vlight/(2.35*R)=64 km/s
* using log10(sigma_v) since the wavelength scale is log                           
      call blur(flux,npix,dble(alog10(64.0)))                             
c      call blur(npflux,npix,dble(alog10(64.0)))  
c      call blur(mflux,npix,dble(alog10(64.0)))
*     offset the no prox flux to be better visible in the plot
* ----------------------------------------------------------------------
*               I N S T R U M E N T A L    N O I S E
* ----------------------------------------------------------------------
* save convolved no-noise flux real*8
      nnflux=flux
* add BOSS noise : 
      numpix=npix*3                                                       !total number of pixels, Smee et al: spectrum width ~ 3 pix
      t=60*90                                                            !integration time, usually 4x15 min = 60x60 [sec]
      area=pi*250**2                                                      !telescope collecting area [cm^2], radius=250 cm
      pix=(1./vlight)*69*lambda                                           !pixel size [A], Bolton et al. 2012. : BOSS pixels constant in velocity space: 69 km/s
      en=(vlight*1e3)*planckh/((lambda)*1e-10)*1e17                       !incident photon energy [1e-17 ergs]
c      qe=0.20*exp(-.5*(((lambda-5080.0)/8.1e2)**2))                       !throughput model: 5 gaussians
c     +     +0.17*exp(-.5*(((lambda-7200.0)/7e2)**2))                      !modelled after fig 38 in Smee et al. 2012
c     +     +0.15*exp(-.5*(((lambda-9000.0)/8e2)**2))
c     +     +0.10*exp(-.5*(((lambda-7000.0)/20e2)**2))
c     +     -0.19*exp(-.5*(((lambda-7550.0)/20)**2))
      do i=1,npix
         qe(i)=throughput(lambda(i))
      end do
      signal=nnflux*t*pix*area/(en*qe)                                    !signal [in incident photons]
      sky=abs(gasdev3(idum)*10**1.2)                                      !background signal in incident photons  (previously: 1e6.2)
      do i=1,npix
*     use different values for red and blue arms of the spectrographs
         if((lambda(i).ge.3600.0).and.(lambda(i).lt.6050.))then
            dark(i)=(0.525+0.022*gasdev3(idum))/900                              !e-/pixel/s
            gain(i)=1.02+0.01*gasdev3(idum)                                  !e-/ADU
         else
            dark(i)=(1.065+0.161*gasdev3(idum))/900                              !e-/pixel/s
            gain(i)=1.70+0.07*gasdev3(idum)                                  !e-/ADU
         end if
      end do
      ronsq=(2.25+0.2*gasdev3(idum))**2                                   ! read-out**2 noise (e-) (per read-out!, interval = 55sec)
      ns=sqrt(signal+sky+(dark*t)*numpix+ronsq*(t/55))                    ! standard deviation of noise in e-
      ns=ns*en/qe/(t*pix*area)                                            ! standard deviation of noise in units of flux
      ivar=1/ns**2                                                        ! inverse variance of noise (units flux)
      do i=1,npix
         d = gasdev3(idum)
         noise(i)=d*ns(i) ! noise in erg s-1 cm-1 A-1
      end do
      flux = nnflux+noise
      snr=abs(nnflux/noise)

* ----------------------------------------------------------------------
*               P L O T T I N G    T H E    S P E C T R U M
* ----------------------------------------------------------------------
c     PLOT
*     view       =1 : plot the spectrum in X-window 
*                =0 : plot the spectrum to a file 
*     log_wav    =1 : log wavelengths
*                =0 : lin wavelengths
*     with_noise =1 : flux with noise
*                =0 : flux with no noise
*     with_conv  =1 : convolved flux 
*                =0 : uncolvolved flux         
      view=0
      log_wav=0
      with_noise=1
      with_conv=0

c      call plotsp(view,log_wav,with_noise,descriptor,zqso,rmag,npix,
c     :loglam,lambda,noabs,plflux,flux,ncflux,nnflux,nmflux,mflux,npflux,
c     :noise)
c      write(6,*) 'Spectrum plotted in a PS file'


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

	nfilt=int(50.0*sigma)+1
	if(nfilt.gt.1011)stop ' too large a filter'
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
!======================================================================
      subroutine plotsp(view,log_wav,with_noise,descriptor,zqso,rmag,
     :npix,loglam,lambda,noabs,plflux,flux,ncflux,nnflux,nmflux,mflux,
     :npflux,noise)
*     plot the spectrum in a Xwindow or a PS file
!======================================================================
      integer log_wav,view,with_noise,npix
      real*8 zqso,rmag,z
      real*8,dimension(npix) :: loglam,lambda,flux,ncflux,nnflux,nmflux
      real*8,dimension(npix) :: mflux,noise,noabs,npflux
      real*4,dimension(npix) :: plotx,ploty,plflux
      real*4 xmin,xmax,ymin,ymax,loc
*     LINELIST
      real*8 colden,bpar,gasdev3
      type list
         SEQUENCE
         character*2 atom
         character*4 ion
         real*8 colden
         real*8 rdf
         real*8 bpar
      end type
      type (list) :: llist(5000)
      common /linelist/ llist
      character :: string*10, descriptor*15, name*35

      if (log_wav.eq.0) then
         plotx=real(lambda)
         xmin=minval(lambda)
         xmax=maxval(lambda)
      else if (log_wav.eq.1) then
         plotx=real(loglam)
         xmin=real(wstart)
         xmax=3.8               !real(wend)
      end if
      if (with_noise.eq.0) then
         ploty=real(nnflux)
      else if (with_noise.eq.1) then
         ploty=real(flux)
      end if 
      
      
      ymin=minval(ploty)
      ymax=maxval(ploty)
      ymax=ymax+0.1*abs(ymax)
      ymin=ymin-0.2*abs(ymin)
      
*     (1) plot flux
      name=trim('./plots/2_flux.ps')
      if (view.eq.0) then
         call PGBEGIN(0,name//'/cps',1,1)
      else if (view.eq.1) then
         call PGBEGIN(0,'/xserve',1,1)
      end if   
      call PGSLW(3)
      call PGENV(xmin,xmax,ymin,ymax,0,0,1)
      if (log_wav.eq.1) then
      call PGLABEL('log \gl','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','')
      else 
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      end if
      call PGMTXT('T',1.0,0.5,0.5,descriptor)
      write(string,'(f5.2)') rmag
      call PGMTXT('T',1.0,1.0,1.0,'r='//trim(string))
      write(string,'(f8.6)') zqso
      call PGMTXT('T',1.0,0.0,0.0,'z='//trim(string))
      call PGSLW(3)
      call PGSCI(2)
c      call PGLINE(npix,plotx,real(nmflux))
      call PGSCI(1)
      call PGLINE(npix,plotx,ploty)
      call PGEND

*     (2) plot no-noise flux
      ploty=real(nnflux)
      name=trim('./plots/2_nnflux.ps')
      call PGBEGIN(0,name//'/cps',1,1)
      call PGSLW(3)
      call PGENV(xmin,xmax,ymin,ymax,0,0,1)
      if (log_wav.eq.1) then
      call PGLABEL('log \gl','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','')
      else 
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      end if
      call PGMTXT('T',1.0,0.5,0.5,descriptor)
      write(string,'(f5.2)') rmag
      call PGMTXT('T',1.0,1.0,1.0,'r='//trim(string))
      write(string,'(f8.6)') zqso
      call PGMTXT('T',1.0,0.0,0.0,'z='//trim(string))
      call PGSLW(3)
      call PGSCI(2)
c      call PGLINE(npix,plotx,real(nmflux))
      call PGSCI(1)
      call PGLINE(npix,plotx,ploty)
      call PGEND

*     (3) plot uncolvolved flux
      ploty=real(ncflux)
      name=trim('./plots/2_ncflux.ps')
      call PGBEGIN(0,name//'/cps',1,1)
      call PGSLW(3)
      call PGENV(xmin,xmax,ymin,ymax,0,0,1)
      if (log_wav.eq.1) then
      call PGLABEL('log \gl','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','')
      else 
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      end if
      call PGMTXT('T',1.0,0.5,0.5,descriptor)
      write(string,'(f5.2)') rmag
      call PGMTXT('T',1.0,1.0,1.0,'r='//trim(string))
      write(string,'(f8.6)') zqso
      call PGMTXT('T',1.0,0.0,0.0,'z='//trim(string))
      call PGSLW(3)
      call PGSCI(2)
c      call PGLINE(npix,plotx,real(nmflux))
      call PGSCI(1)
      call PGLINE(npix,plotx,ploty)
      call PGEND

*     (4) plot power-law
      ploty=real(plflux)
      name=trim('./plots/2_plflux.ps')
      call PGBEGIN(0,name//'/cps',1,1)
      call PGSLW(3)
      call PGENV(xmin,xmax,ymin,ymax,0,0,1)
      if (log_wav.eq.1) then
      call PGLABEL('log \gl','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','')
      else 
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      end if
      call PGMTXT('T',1.0,0.5,0.5,descriptor)
      write(string,'(f5.2)') rmag
      call PGMTXT('T',1.0,1.0,1.0,'r='//trim(string))
      write(string,'(f8.6)') zqso
      call PGMTXT('T',1.0,0.0,0.0,'z='//trim(string))
      call PGSLW(3)
      call PGSCI(2)
c      call PGLINE(npix,plotx,real(nmflux))
      call PGSCI(1)
      call PGLINE(npix,plotx,ploty)
      call PGEND

*     (5) plot noabs
      ploty=real(noabs)
      name=trim('./plots/2_noabs.ps')
      call PGBEGIN(0,name//'/cps',1,1)
      call PGSLW(3)
      call PGENV(xmin,xmax,ymin,ymax,0,0,1)
      if (log_wav.eq.1) then
      call PGLABEL('log \gl','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','')
      else 
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      end if
      call PGMTXT('T',1.0,0.5,0.5,descriptor)
      write(string,'(f5.2)') rmag
      call PGMTXT('T',1.0,1.0,1.0,'r='//trim(string))
      write(string,'(f8.6)') zqso
      call PGMTXT('T',1.0,0.0,0.0,'z='//trim(string))
      call PGSLW(3)
      call PGSCI(2)
c      call PGLINE(npix,plotx,real(nmflux))
      call PGSCI(1)
      call PGLINE(npix,plotx,ploty)
      call PGEND

*     (6) plot metals
      ploty=real(mflux)
      name=trim('./plots/2_metal.ps')
      call PGBEGIN(0,name//'/cps',1,1)
      call PGSLW(3)
      call PGENV(xmin,xmax,ymin,ymax,0,0,1)
      if (log_wav.eq.1) then
      call PGLABEL('log \gl','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','')
      else 
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      end if
      call PGMTXT('T',1.0,0.5,0.5,descriptor)
      write(string,'(f5.2)') rmag
      call PGMTXT('T',1.0,1.0,1.0,'r='//trim(string))
      write(string,'(f8.6)') zqso
      call PGMTXT('T',1.0,0.0,0.0,'z='//trim(string))
      call PGSLW(3)
      call PGSCI(2)
c      call PGLINE(npix,plotx,real(nmflux))
      call PGSCI(1)
      call PGLINE(npix,plotx,ploty)
      call PGEND


*     (7) plot everything on a single page
      name=trim('./plots/1_combined.ps')
      call PGBEGIN(0,name//'/cps',1,3)
      call PGSLW(2)
      call PGENV(xmin,xmax,ymin,ymax,0,0)
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      call PGMTXT('t',1.5,0.5,0.5,'(a)')
      ploty=real(plflux)
      call PGLINE(npix,plotx,ploty)
      call PGPANL(1,2)
      call PGBOX('BCNTS',0.0,0,'BCNTS',0.0,0)
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      ploty=real(noabs)
      call PGMTXT('t',1.5,0.5,0.5,'(b)')
      call PGLINE(npix,plotx,ploty)
      call PGPANL(1,3)
      call PGBOX('BCNTS',0.0,0,'BCNTS',0.0,0)
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      call PGMTXT('t',1.5,0.5,0.5,'(c)')
      ploty=real(mflux)
      call PGLINE(npix,plotx,ploty)
      call PGEND


      name=trim('./plots/combined_6.ps')
      call PGBEGIN(0,name//'/cps',1,3)
      call PGSLW(2)
      call PGENV(xmin,xmax,ymin,ymax,0,0)
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      call PGMTXT('t',1.5,0.5,0.5,'(d)')
      ploty=real(ncflux)
      call PGLINE(npix,plotx,ploty)
      call PGPANL(1,2)
      call PGBOX('BCNTS',0.0,0,'BCNTS',0.0,0)
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      call PGMTXT('t',1.5,0.5,0.5,'(e)')
      ploty=real(nnflux)
      call PGLINE(npix,plotx,ploty)
      call PGPANL(1,3)
      call PGBOX('BCNTS',0.0,0,'BCNTS',0.0,0)
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      call PGMTXT('t',1.5,0.5,0.5,'(f)')
      ploty=real(flux)
      call PGLINE(npix,plotx,ploty)
      call PGEND

*     (8) plot no-noise flux in X-window
      ploty=real(flux)
      call PGBEGIN(0,'/xserve',1,1)
      call PGSLW(3)
      call PGENV(xmin,xmax,ymin,ymax,0,0,1)
      if (log_wav.eq.1) then
      call PGLABEL('log \gl','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','')
      else 
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      end if
      call PGMTXT('T',1.0,0.5,0.5,descriptor)
      write(string,'(f5.2)') rmag
      call PGMTXT('T',1.0,1.0,1.0,'r='//trim(string))
      write(string,'(f8.6)') zqso
      call PGMTXT('T',1.0,0.0,0.0,'z='//trim(string))
      call PGSLW(3)
      call PGSCI(2)
      do i=1,3000
         colden=llist(i)%colden
         if (colden.gt.1d18)then
            z=llist(i)%rdf
            loc=real((z+1))*912.0
            call pgsci(2)
            call pgmove(loc,0.95*ymax); call pgdraw(loc,0.0)
            loc=real((z+1))*1215.67
            call pgsci(3)
            call pgmove(loc,0.95*ymax); call pgdraw(loc,0.0)
         end if
      end do
      call PGSCI(1)
      call PGLINE(npix,plotx,ploty)
      call PGEND

*     (9) plot the zoom in of the proximity effect 
      ploty=real(nmflux)
      xmin = 1215.67*(1+0.90*real(zqso))
      xmax = 1215.67*(1+1.05*real(zqso))
      name=trim('./plots/zoom_pe.ps')
      call PGBEGIN(0,name//'/cps',1,1)
      call PGSLW(3)
      call PGENV(xmin,xmax,ymin,ymax,0,0,1)
      if (log_wav.eq.1) then
      call PGLABEL('log \gl','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','')
      else 
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      end if
      call PGMTXT('T',1.0,0.5,0.5,descriptor)
      write(string,'(f5.2)') rmag
      call PGMTXT('T',1.0,1.0,1.0,'r='//trim(string))
      write(string,'(f8.6)') zqso
      call PGMTXT('T',1.0,0.0,0.0,'z='//trim(string))
      call PGSLW(3)
      call PGSCI(1)
      call PGLINE(npix,plotx,ploty)
      call PGEND

*     (10) plot the zoom in with no proximity effect
      ploty=real(npflux)
      name=trim('./plots/zoom_npe.ps')
      call PGBEGIN(0,name//'/cps',1,1)
      call PGSLW(3)
      call PGENV(xmin,xmax,ymin,ymax,0,0,1)
      if (log_wav.eq.1) then
      call PGLABEL('log \gl','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','')
      else 
      call PGLAB('\gl [\A]','F\d\gl\u [10\u-17\d erg/s/cm\u2\d/\A]','') 
      end if
      call PGMTXT('T',1.0,0.5,0.5,descriptor)
      write(string,'(f5.2)') rmag
      call PGMTXT('T',1.0,1.0,1.0,'r='//trim(string))
      write(string,'(f8.6)') zqso
      call PGMTXT('T',1.0,0.0,0.0,'z='//trim(string))
      call PGSLW(3)
      call PGSCI(2)
      call PGSCI(1)
      call PGLINE(npix,plotx,ploty)
      call PGEND      



      return
      end subroutine plotsp
!======================================================================
      function throughput(wav)
!======================================================================
      real*8 wav,throughput,f
      integer,parameter :: m=4
      real*8 c(m),mu(m),sig(m)

      data c/0.247670970667,0.25436483658,0.176176984517,-0.14316777838/
      data mu/5002.72401967,7108.43270986,8897.07074493,7566.1928085/
      data sig/722.981333002,897.918049886,727.955532937,22.3997766177/

      f=0.0; throughput=0.0
      do i=1,10
         f=c(i)*exp(-(wav-mu(i))**2/(2*sig(i)**2))
         throughput=throughput+f
      end do
      return
      end function throughput
!======================================================================
      subroutine fmetals(x,y,z,ind)
*     locate a point in parametric space (Z,NHI,z) that is closest in 
*     value to the input and return the index of the line in the dtbase
!======================================================================
      integer i,j,k,xmin,ymin,zmin,ind
*     values: x = metallicity, y = log NHI, z = redshift
      real*8 :: x,y,z
      real*8,dimension(20000) :: xbuf,ybuf,zbuf,xdist,ydist,zdist,rad
      real*8 xval,yval,zval

      common/parspace/xbuf,ybuf,zbuf

      xdist=(x-xbuf)**2
      ydist=(y-ybuf)**2
      zdist=(z-zbuf)**2
      rad=xdist+ydist+zdist
      ind=minloc(rad,dim=1)
      xmin=minloc(xdist,dim=1)
      ymin=minloc(ydist,dim=1)
      zmin=minloc(zdist,dim=1)
      xval=xbuf(xmin)
      yval=ybuf(ymin)
      zval=zbuf(zmin)
      
      return
      end subroutine fmetals
      
!======================================================================
      subroutine addsys(j,nhi,z,bvel)
!======================================================================
      integer i,j,idum,ind
*     INPUTS
      real*8 :: nhi,z,bvel
*     METALS
      real*8,dimension(30,20000) :: H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,
     :                      P,S,Cl,Ar,K,Ca,Sc,Ti,Va,Cr,Mn,Fe,Co,Ni,Cu,Zn
      common/metals/H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,
     :                      P,S,Cl,Ar,K,Ca,Sc,Ti,Va,Cr,Mn,Fe,Co,Ni,Cu,Zn
*     LINELIST
      real*8 colden,bpar,gasdev3
      type list
         SEQUENCE
         character*2 atom
         character*4 ion
         real*8 colden
         real*8 rdf
         real*8 bpar
      end type
      type (list) :: llist(5000)
      common /linelist/ llist
*     LOCAL VARIABLES
      real*8 d,t,m,rdf,lognhi
      character atom*2,ion*4
      external gasdev3

      idum=time()
*     HI line parameters
      colden=nhi
      bpar=bvel
      rdf=z
      llist(j)=list('H ','I   ',colden,rdf,bpar)
*     if log NHI > 17, input metal lines as well
      if (nhi.ge.(1e17)) then
         s=gasdev3(idum); t=gasdev3(idum)
         m = -(0.22+0.03*d)*z - (0.65+0.09*t)
*        find the closest metal point in the database
         lognhi=alog10(real(nhi))
         call fmetals(m,lognhi,z,ind)
*        C:
         atom = 'C '
           j=j+1; ion='I   '; colden=C(1,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='II  '; colden=C(2,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='III '; colden=C(3,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='IV  '; colden=C(4,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
         atom = 'N '
           j=j+1; ion='I   '; colden=N(1,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='II  '; colden=N(2,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='III '; colden=N(3,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='V   '; colden=N(4,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
         atom = 'O'
           j=j+1; ion='I   '; colden=O(1,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='VI  '; colden=O(6,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
         atom = 'Mg'
           j=j+1; ion='I   ';colden=Mg(1,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='II  ';colden=Mg(2,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
         atom = 'Al'
           j=j+1; ion='II  ';colden=Al(2,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='III ';colden=Al(3,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
         atom = 'Si'
           j=j+1; ion='II  ';colden=Si(2,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='III ';colden=Si(3,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='IV  ';colden=Si(4,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
         atom = 'Cr'
           j=j+1; ion='II  ';colden=Cr(2,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
         atom = 'Mn'
           j=j+1; ion='II  ';colden=Mn(2,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
         atom = 'Fe'
           j=j+1; ion='II  ';colden=Fe(2,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
           j=j+1; ion='III ';colden=Fe(3,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
         atom = 'Ni'
           j=j+1; ion='II  ';colden=Ni(2,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
         atom = 'Zn'
           j=j+1; ion='II  ';colden=Zn(2,ind);bpar=1.5*gasdev3(idum)+6.3
           llist(j)=list(atom,ion,colden,rdf,bpar)
      end if
      return
      end subroutine addsys
!======================================================================
      subroutine rwllist(d)
*     delete the line list (remove artefacts)
!======================================================================
      integer d
      type list
         SEQUENCE
         character*2 atom
         character*4 ion
         real*8 colden
         real*8 rdf
         real*8 bpar
      end type
      type (list) :: llist(5000)
      common /linelist/ llist
      
      if (d.eq.0) then
         do i=1,3000
            llist%atom = '  '
            llist%ion = '    '
            llist%colden = 0.0
            llist%rdf = 0.0
            llist%bpar = 0.0
         end do
      else 
         write(6,*) 'WARNING: Line list was not reset!'
         write(6,*) 'WARNING: Possible artefacts'
      end if
      return
      end subroutine rwllist
!======================================================================
