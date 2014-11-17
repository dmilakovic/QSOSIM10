!======================================================================
      PROGRAM qsosim10
c  PURPOSE: call other subroutines to generate artificial spectra!
c  OUTPUT:  artificial SDSS catalogue
!======================================================================
      IMPLICIT NONE
      CHARACTER :: infile*20, outfile*30, descriptor*16,command*20
      CHARACTER :: home*120,name*18, splate*4,smjd*5,sfiber*4
      INTEGER,PARAMETER :: npoints=1000
      INTEGER :: i,j,inoise, numlin, npts, nl,s,option,nrows
      INTEGER :: z_w,id,iplate,imjd,ifiber,ipix
      REAL*8 :: wstart,wend,dw,mags(5)
      REAL*8 :: nc,nuplim,sigblur
      REAL*8 :: corr, zstart,zend
      REAL*8,DIMENSION(3) :: bigA,gamma,n
      INTEGER,DIMENSION(3) :: ni
      CHARACTER,DIMENSION(:),allocatable :: SDSS_name*18
      INTEGER,DIMENSION(:),allocatable :: thing_id,plate,mjd,
     &                                   fiber,z_flag,npix
      REAL*8, DIMENSION(:),allocatable :: ra,dec,zqso,alpha,alpha_fit,
     &                                   begin_wave,rmag,ovi,civ
c      REAL*8,DIMENSION(:,:),alloacatable :: psfmag
      REAL*8, DIMENSION(82701,5) :: psfmag
      REAL*8,DIMENSION(262144) :: lambda,flux, da4
      REAL*8,DIMENSION(262144) :: flerr, nnflux,flux_nc
      REAL*8, DIMENSION(:),ALLOCATABLE :: xs,ys,CDDF,H
      REAL*8,DIMENSION(:),allocatable :: nhi4,z4
      REAL*8,DIMENSION(:),allocatable :: nciv, novi
      INTEGER,DIMENSION(:),allocatable :: mask
      EXTERNAL :: qsosim9, spline, readfits, writefits, power_laws
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
      !   [17.0,22.0]      log(0.175)   1.33

      data bigA/33.11311,5.24807,0.1750/
      data gamma/1.51,2.16,1.33/

! ---------------------------------------------------------------------
! READ INPUT PARAMETERS
! ---------------------------------------------------------------------
      home = '/Users/dm/Documents/GitHub/QSOSIM10'
      write (6,*)'=================================================='
      write (6,*)'                    QSOSIM 10                     '
      write (6,*)'=================================================='
      write (6,*)'What would you like to create?'
      write (6,50)'(1)','SDSS'
      write (6,50)'(2)','random (legacy)'
      write (6,*)'=================================================='
      write (6,'(1x,a13)',advance='no')'I choose... '
      read (*,'(i1)') option
      if (option.eq.1) then
         infile='DR10Q_r.fits'
         call read_nrows(home,infile,nrows)
         write(6,*) nrows
         allocate(SDSS_name(nrows)); allocate(thing_id(nrows));
         allocate(plate(nrows));     allocate(mjd(nrows))
         allocate(fiber(nrows));     allocate(z_flag(nrows))
         allocate(npix(nrows));      allocate(ra(nrows))
         allocate(dec(nrows));       allocate(zqso(nrows))
         allocate(alpha(nrows));     allocate(alpha_fit(nrows))
         allocate(begin_wave(nrows));allocate(rmag(nrows))
c         allocate(psfmag(nrows,5)
         call SDSS_readfits(home,infile,nrows,SDSS_name,RA,DEC,thing_id,
     &         plate,mjd,fiber,zqso,z_flag,alpha,alpha_fit,npix,
     &         begin_wave,psfmag,rmag)
         nc=1e12
         nuplim=1e22
         inoise=1
         sigblur=3.0
         dw=0.0001 !pixel resolution of SDSS
c      else if (option.eq.2) then        
c         infile='sin.fits'
c         call readfits(infile,wstart,wend,dw,nc,nuplim,inoise,dvavoid,
c     &         ra,dec,zqso,alpha,rmag,sigblur,s2n)
      end if
 50   format (3x,a3,2x,a15)
 100  format(f10.5,2x,f10.5,2x,f8.5,2x,f8.5)

! ---------------------------------------------------------------------
! GET THE SPLINE
! ---------------------------------------------------------------------
      allocate(xs(npoints))
      allocate(ys(npoints))
      allocate(CDDF(npoints))
      call spline(npoints,nc,nuplim,xs,ys,CDDF)
! ---------------------------------------------------------------------
! GENERATE ARTIFICIAL SPECTRA
! ---------------------------------------------------------------------
      do i=1,1
         wstart=begin_wave(i)
         wend=begin_wave(i)+npix(i)*dw 
         zstart=(10**begin_wave(i)/1215.67)-1.
         zend=zqso(i)
         write (splate,'(I4.4)') plate(i)
         write (smjd,'(I5.5)') mjd(i)
         write (sfiber,'(I4.4)') fiber(i)         
         descriptor = splate//'-'//smjd//'-'//sfiber
         write (6,*)'=================================================='
         write (6,*)'               Spectrum no. ',descriptor
         write (6,*)'=================================================='

         call power_laws(npoints,zstart,zqso(i),xs,ys,CDDF,
     +                      bigA,gamma,corr,nl,ni)
         allocate(nhi4(nl)); allocate(z4(nl)); allocate(mask(nl))
         call assign(npoints,zstart,zqso(i),xs,CDDF,gamma,nl,ni,nhi4,z4)
         call write_cloudy_input(home,descriptor,zqso(i),nl,nhi4,mask,s)
         call cloudy(home,descriptor,mask,nl,s)
         allocate(ovi(s)); allocate(civ(s))
c         call read_cloudy_output()
         call qsosim9(zqso(i),alpha_fit(i),rmag(i),wstart,wend,dw,nuplim
     +         ,sigblur,inoise,npts,lambda,flux,
     +         flerr,nnflux,flux_nc,npoints,nl,ni,nhi4,z4)
         outfile='mockspec-'//descriptor//'.fits'
         do j=1,5
            mags(j)=psfmag(i,j)
         end do
         call writefits(outfile,ra(i),dec(i),zqso(i),z_flag(i),alpha(i),
     &                  alpha_fit(i),rmag(i),SDSS_name(i),thing_id(i),
     &                  plate(i),mjd(i),fiber(i),npix(i),wstart,
     &                  mags,npts,lambda,flux,flerr,nnflux,flux_nc)
         write (*,*)'--------------------------------------------------'
         deallocate(nhi4)
         deallocate(z4)
         deallocate(mask)
      end do
      
! ---------------------------------------------------------------------
! PLOT QSO SPECTRUM OF THE LAST SOURCE
! --------------------------------------------------------------------- 
      call PGBEGIN (0,'/null',1,1)
c      call PGSLW(1)
      call PGENV (12.0,22.0,0.0,1.0,0,1)
c      call PGLABEL ('lambda','flux','QSO spectrum')
      call pgline(npoints,real(xs),real(CDDF))
c      call pgsci(2)
c      call pgline(npts,real(lambda),real(nnflux))
c      call pgsci(1)
c      call PGENV(11.5,22.5,2.00,2.1,0,1)
c      call PGLABEL('log NHI','z','Random choice of NHI & redshift')
c      call PGPT(nl,real(nhi4),real(z4),3)
      call PGEND
! --------------------------------------------------------------------- 
! FREE ALLOCATED MEMORY
! --------------------------------------------------------------------- 
      deallocate(xs)
      deallocate(ys)
      deallocate(CDDF)
!======================================================================
      END PROGRAM qsosim10
!======================================================================
