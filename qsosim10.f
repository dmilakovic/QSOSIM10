!======================================================================
      PROGRAM qsosim10
c  PURPOSE: call other subroutines to generate artificial spectra!
c  OUTPUT:  artificial SDSS catalogue
!======================================================================
      IMPLICIT NONE
      CHARACTER :: infile*20, outfile*35, descriptor*15,command*20
      CHARACTER :: home*120,gohome*120,mockspec_folder*120,name*18
      CHARACTER :: val*6, valer*6
*     npoints - number of points for spline evaluation
*     nrows   - number of QSOs in the input file
*     nn      - number of QSOs to analyse the line distribution
*     dim     - number of points in the metal database
      INTEGER,PARAMETER :: npoints=5000,nrows=82701, nn=1000
      integer, parameter :: RegInt_K = selected_int_kind (12)
      INTEGER (kind=RegInt_K) :: i
      INTEGER :: nr,j,inoise,npts, nhil,nlin,option,flag,dim
      INTEGER :: z_w,id,iplate,imjd,ifiber,ipix,status,maxr
      REAL*8 :: wstart,wend,dw,mags(5)
      REAL*8 :: nc,nuplim,sigblur
      REAL*8 :: corr, zstart,zend
c      REAL*4 :: a,b,chi2,q,siga,sigb,xx(nn),yy(nn),xmin,xmax
      REAL*8,DIMENSION(3) :: bigA,gamma,numl
      REAL*4,DIMENSION(3,nn) :: np,sig
      CHARACTER,DIMENSION(nrows) :: SDSS_name*18
      INTEGER,DIMENSION(nrows) :: thing_id,plate,mjd,
     &                                   fiber,z_flag,npix
      REAL*8, DIMENSION(nrows) :: ra,dec,zqso,alpha,alpha_fit,
     &                                   begin_wave,rmag
      REAL*8,DIMENSION(nrows,5) :: psfmag
      REAL*8,DIMENSION(:),allocatable :: loglam,flux,ivar
      REAL*8,DIMENSION(:),allocatable :: noise, nnflux,ncflux, noabs
      REAL*8, DIMENSION(npoints) :: xs,ys,CDDF
      REAL*8,DIMENSION(:),allocatable :: nhi4,z4,npnhi4,znc4  !npz4 - testing purposes
      REAL*8,DIMENSION(:),allocatable :: nciv4, novi4
      
      EXTERNAL :: read_nrows, SDSS_readfits, qsosim, spline, readfits, 
     &            power_laws, writefits, assign, write_cloudy_input
     &            cloudy

!====================================================================== 
      home = '/Users/dm/Documents/GitHub/QSOSIM10'
! ---------------------------------------------------------------------
! DATA USED
! ---------------------------------------------------------------------
      ! from Kim et al 2014
      !     log NHI        log A         gamma
      !   [13.1,14.0]      1.52+-0.05   1.51+-0.09
      !   [14.0,17.0]      0.72+-0.08   2.16+-0.14
      ! from Unknown
      !   [17.0,22.0]      log(0.15133) 1.94+-0.18

      data bigA/33.11311,5.24807,0.15133/
      data gamma/1.51,2.16,1.94/

! ---------------------------------------------------------------------
! WRITE BASIC INFO ON THE PROGRAM
! ---------------------------------------------------------------------
      
      write (6,*)'===================================================='
      write (6,*)'                 Q S O S I M    10                  '
      write (6,*)'===================================================='
      write (6,*)'A  program  to  create  artificial  QSO spectra  and'
      write (6,*)"simulate the Sloan Digital Sky Survey's QSO data set"
      write (6,*)'Dated: 2015/2/6'
! ---------------------------------------------------------------------
! READ THE METAL DATABASE
! ---------------------------------------------------------------------
      write (6,*)'===================================================='
      write (6,*)'Reading the metal database'
      call metal_db(0)
      write (6,*)'===================================================='
! ---------------------------------------------------------------------
! READ INPUT PARAMETERS
! ---------------------------------------------------------------------
      home = '/Users/dm/Documents/GitHub/QSOSIM10'
      gohome = trim('cd '//home)
      mockspec_folder = trim(home)//'/mockspectra'
      write (6,*)'What would you like to create?'
      write (6,50)'(1)','SDSS'
c      write (6,50)'(2)','random (legacy)'
      write (6,*)'=================================================='
      write (6,'(1x,a13)',advance='no')'I choose... '
c      read (*,'(i1)') option
      option=1
      if (option.eq.1) then
         infile='DR10Q_r.fits'
         call read_nrows(home,infile,nr)
         write(6,*) 'Number of QSO spectra in the database: ',nr
         call SDSS_readfits(home,infile,nrows,nr,SDSS_name,RA,DEC,
     &         thing_id,plate,mjd,fiber,zqso,z_flag,alpha,alpha_fit,
     &         npix,begin_wave,psfmag,rmag)
         nc=1e12
         nuplim=1e22
         inoise=1
         sigblur=3.0
         dw=0.0001 !pixel size of SDSS in log(lambda), lambda in [A]
c      else if (option.eq.2) then        
c         infile='sin.fits'
c         call readfits(infile,wstart,wend,dw,nc,nuplim,inoise,dvavoid,
c     &         ra,dec,zqso,alpha,rmag,sigblur,s2n)
      end if
 50   format (3x,a3,2x,a15)
 100  format(f10.5,2x,f10.5,2x,f8.5,2x,f8.5)

! ---------------------------------------------------------------------
! GET THE SPLINE
! Column density distribution function of neutral hydrogen f(N_HI,X)
! is given by Prochaska et al. (2014) in the form of a Hermite cubic 
! spline. 'Spline' interpolates the spline to get the Column Density
! Distribution function - CDDF - at points xs
! ---------------------------------------------------------------------
      call spline(npoints,nc,nuplim,xs,ys,CDDF)
! ---------------------------------------------------------------------
! GET THE METAL DATABASE
! Saved in metals.csv
! ---------------------------------------------------------------------
      call metal_db(0)
! ---------------------------------------------------------------------
! GENERATE ARTIFICIAL SPECTRA
! use data read from the input file to generate spectra in a loop
! ---------------------------------------------------------------------
      flag=0
      do i=48555,48555 !1,nn
         !--------------------------------------------------------------
*      (1) Set the unique description of the QSO
         write (descriptor,"(i4.4,'-',i5.5,'-'i4.4)")
     &                     plate(i),mjd(i),fiber(i)
         write (6,*)'=================================================='
         write (6,*)'            Spectrum no. ',descriptor
         write (6,*)'=================================================='
         write (6,70) i,nr
 70      format ('       Iteration     ',i5,'/',i5)
         !--------------------------------------------------------------
*      (2) Set basic QSO data 
         !**** begin and end wavelengths for the spectrum, start and end
         !**** redshifts for the Lyman alpha forest
         wstart=begin_wave(i) 
         wend=begin_wave(i)+npix(i)*dw 
         zstart=(10**begin_wave(i)/1215.67)-1.
         zend=zqso(i)
         ipix=npix(i)   
         allocate(loglam(ipix));allocate(flux(ipix))
         allocate(noise(ipix));allocate(nnflux(ipix))
         allocate(ncflux(ipix));allocate(noabs(ipix))
         allocate(ivar(ipix))
         write (6,80) 'z_start =',zstart,'z_qso =',zend
         write (6,*)'u:',psfmag(i,1),'g:',psfmag(i,2),'r:',psfmag(i,3),
     :              'i:',psfmag(i,4),'z:',psfmag(i,5)
         write (6,*) 'SDSS name:',SDSS_name(i)
         write (6,*) 'alpha:',alpha(i), alpha_fit(i)
 80      format (1x,a9,f6.4,3x,a7,f6.4,3x)
*        Sometimes Lyman alpha forest can't be seen in the 
*        spectrum. When that happens, skip the spectrum generation.
         if (zstart.ge.zend.or.rmag(i).lt.0.0) then
            flag=flag+1
            write(6,*) 'zstart >= zend!'; goto 99
         end if
         !--------------------------------------------------------------               
*      (3) calculate the number of lines between zstart and zqso
         call power_laws(npoints,zstart,zend,xs,ys,CDDF,
     +                      bigA,gamma,nhil)
*        allocate memmory to: NHI , z , mask 
         allocate(nhi4(nhil)); allocate(z4(nhil)); allocate(znc4(nhil))
*        no proximity effect array (npnhi4)
         allocate(npnhi4(nhil))
*      (4) use the CDDF to sample the column density distribution
         call assign(npoints,zstart,zend,rmag(i),alpha_fit(i),xs,CDDF,
     +               nhil,nhi4,npnhi4,z4,znc4)
*      (5) generate artificial spectrum
         call qsosim(descriptor,zqso(i),alpha_fit(i),rmag(i),
     +         begin_wave(i),dw,npix(i),nhil,nhi4,npnhi4,z4,znc4,
     +         loglam,flux,noise,ivar,nnflux,ncflux,noabs,nlin)
c         read (*,'(i1)') option
*      (6) write qsosim output + qso general data into a fits file
         outfile='mockspec-'//descriptor//'.fits'
         
         do j=1,5
            mags(j)=psfmag(i,j)
         end do

         call writefits(home,mockspec_folder,outfile,ra(i),dec(i),
     +                  zqso(i),z_flag(i),alpha(i),alpha_fit(i),rmag(i),
     +                  SDSS_name(i),thing_id(i),plate(i),mjd(i),
     +                  fiber(i),npix(i),wstart,mags,loglam,flux,noise,
     +                  ivar,nnflux,ncflux,noabs,nlin)
         write (6,*)'--------------------------------------------------'

         !**** deallocate the arrays ****
         deallocate(nhi4); deallocate(z4) 
         deallocate(npnhi4); deallocate(znc4)
         deallocate(loglam);deallocate(flux)
         deallocate(noise);deallocate(nnflux)
         deallocate(ncflux);deallocate(noabs)
         deallocate(ivar)
         write(6,*) 'ARRAYS DEALLOCATED'
 99      continue
c         read (*,'(i1)') option
      end do
      
*     finished generating the spectra!
      write (6,*)'=================================================='
      write (6,*)'               DATA SET COMPLETED!                '
      write (6,*)'=================================================='
      write (6,*)'QSOs with flags :',flag

*     plot predicted values of the number of lines vs. actually assigned
c$$$      a=0.0;b=0.0;chi2=0.0;
c$$$      do i=1,3
c$$$         sig(i,:)=sqrt(nexp(i,:))
c$$$      end do
c$$$      do i=1,3
c$$$         chi2=0.0
c$$$         write(val,'(i1)') i
c$$$         name=trim('./plots/lines_'//val//'.ps')
c$$$         call fit(np(i,:),nexp(i,:),nn,sig(i,:),0,a,b,siga,sigb,chi2)
c$$$         chi2=chi2/(2*nn)
c$$$         write(6,*) a,b,chi2/(2*nn),siga,sigb
c$$$         xmax=max(maxval(np(i,:)),maxval(nexp(i,:)))
c$$$         do j=1,nn
c$$$            xx(j)=xmax*dble(j-1)/(nn-1)
c$$$            yy(j)=a+b*xx(j)
c$$$         end do
c$$$         call pgbegin(0,name//'/cps',1,1)
c$$$         call pgslw(3)
c$$$         
c$$$         call pgenv(0.0,xmax,0.0,xmax,1,0)
c$$$         call pglabel('Theoretical prediction','QSOSIM10','')
c$$$         if (i.eq.1) then
c$$$         call pgmtxt('T',0.5,0.0,0.0,'No. of lines - 12<log N\dHI\u<14')
c$$$         else if (i.eq.2) then
c$$$         call pgmtxt('T',0.5,0.0,0.0,'No. of lines - 14<log N\dHI\u<17')
c$$$         else if (i.eq.3) then
c$$$         call pgmtxt('T',0.5,0.0,0.0,'No. of lines - 17<log N\dHI\u<22')  
c$$$         end if
c$$$         write(val,'(f6.2)')a
c$$$         write(valer,'(f6.2)')siga
c$$$         call pgmtxt('T',0.5,1.0,1.0,'\(2175) = a + b\.\(2174)')
c$$$         call pgmtxt('T',-2.0,0.05,0.0,'a='//val//'\(2233)'//valer)
c$$$         write(val,'(f6.4)')b
c$$$         write(valer,'(f6.4)')sigb
c$$$         call pgmtxt('T',-3.5,0.05,0.0,'b='//val//'\(2233)'//valer)
c$$$         write(val,'(f6.2)')chi2
c$$$         call pgmtxt('T',-5.0,0.05,0.0,'\gx\u2\d='//val)
c$$$         call pgmtxt('T',-6.5,0.05,0.0,'N=1000')
c$$$         call pgline(nn,xx,yy)
c$$$         call pgsci(2)
c$$$         call pgpt(nn,np(i,:),nexp(i,:),17)
c$$$         call pgend
c$$$      end do
!======================================================================
      END PROGRAM qsosim10
!======================================================================
