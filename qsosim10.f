!======================================================================
      PROGRAM qsosim10
c  PURPOSE: call other subroutines to generate artificial spectra!
c  OUTPUT:  artificial SDSS catalogue
!======================================================================
      IMPLICIT NONE
      CHARACTER :: infile*20, outfile*20, descriptor*8,command*20
      CHARACTER :: home*120
      INTEGER,PARAMETER :: nrows=25
      INTEGER,PARAMETER :: npoints=1000
      INTEGER :: i,j,inoise, numlin, npts, nl,s
      REAL*8 :: wstart,wend,dw
      REAL*8 :: nc,nuplim,dvavoid
      REAL*8 :: corr, zstart,zend
      REAL*8,DIMENSION(3) :: bigA,gamma,n
      INTEGER,DIMENSION(3) :: ni
      REAL*8,DIMENSION(25) :: ra,dec,zqso,alpha,rmag,s2n,sigblur
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
      !     log NHI        log A        log gamma
      !   [13.1,14.0]      1.52+-0.05   1.51+-0.09
      !   [14.0,17.0]      0.72+-0.08   2.16+-0.14
      ! from Unknown
      !   [17.0,22.0]      log(0.175)   1.33

      data bigA/33.11311,5.24807,0.1750/
      data gamma/1.51,2.16,1.33/

! ---------------------------------------------------------------------
! READ INPUT PARAMETERS
! ---------------------------------------------------------------------
      infile='sin.fits'
      write (6,*)'=================================================='
      write (6,*)'                    QSOSIM 10                     '
      write (6,*)'=================================================='
      call readfits(infile,wstart,wend,dw,nc,nuplim,inoise,dvavoid,
     &             ra,dec,zqso,alpha,rmag,sigblur,s2n)
 100  format(f10.5,2x,f10.5,2x,f8.5,2x,f8.5)
      zstart=(wstart/1215.67)-1.
! ---------------------------------------------------------------------
! GET THE SPLINE
! ---------------------------------------------------------------------
      allocate(xs(npoints))
      allocate(ys(npoints))
      allocate(CDDF(npoints))
      write (6,*) nc, nuplim
      call spline(npoints,nc,nuplim,xs,ys,CDDF)
! ---------------------------------------------------------------------
! GENERATE ARTIFICIAL SPECTRA
! ---------------------------------------------------------------------
      do i=25,25
         zend=zqso(i)
         write (descriptor,'(I8.8)') i
         write (6,*)'=================================================='
         write (6,*)'               Spectrum no. ',descriptor
         write (6,*)'=================================================='

         call power_laws(npoints,zstart,zqso(i),xs,ys,CDDF,
     +                      bigA,gamma,corr,nl,ni)
         allocate(nhi4(nl))
         allocate(z4(nl))
         allocate(mask(nl))
         call assign(npoints,zstart,zqso(i),xs,CDDF,gamma,nl,ni,nhi4,z4)
         call write_cloudy_input(home,descriptor,zqso(i),nl,nhi4,mask,s)
         do j=1,nl
            write(6,*) j, mask(i)
         end do
         call cloudy(home,descriptor,mask,nl,s)
c         call read_cloudy_output()
         call qsosim9(home,zqso(i),alpha(i),rmag(i),wstart,wend,dw,nc,
     +         nuplim,sigblur(i),s2n(i),inoise,dvavoid,npts,lambda,flux,
     +         flerr,nnflux,flux_nc,npoints,nl,ni,nhi4,z4)
         outfile='spec-'//descriptor//'.fits'
         call writefits(outfile,ra(i),dec(i),zqso(i),alpha(i),rmag(i),
     &                     npts,lambda,flux,flerr,nnflux,flux_nc)
         write (*,*)'--------------------------------------------------'
         deallocate(nhi4)
         deallocate(z4)
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
