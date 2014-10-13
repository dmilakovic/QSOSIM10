!======================================================================
      PROGRAM qsosim10
c  PURPOSE: call other subroutines to generate artificial spectra!
c  OUTPUT:  artificial SDSS catalogue
!======================================================================
      IMPLICIT NONE
      CHARACTER :: infile*20, outfile*20, descriptor*6
      INTEGER,PARAMETER :: nrows=25
      INTEGER,PARAMETER :: npoints=10000
      INTEGER :: i,j,inoise, numlin, npts, nl
      INTEGER,DIMENSION(nrows) :: numlls
      REAL*8 :: wstart,wend,dw
      REAL*8 :: nc,nuplim,dvavoid
      REAL*8 :: corr, zstart,zend
      REAL*8,DIMENSION(3) :: bigA,gamma,n
      INTEGER,DIMENSION(3) :: ni
      REAL*8,DIMENSION(nrows) :: ra,dec,zqso,alpha,vmag,s2n,sigblur
      REAL*8 :: lambda(262144),flux(262144), da4(262144)
      REAL*8 :: flerr(262144), nnflux(262144)
      REAL*8, DIMENSION(npoints) :: xs,ys,CDDF,H
      REAL*8,DIMENSION(:),allocatable :: nhi4,z4
      EXTERNAL :: qsosim9, spline, readfits, writefits, power_laws
!====================================================================== 
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
     &             ra,dec,zqso,alpha,vmag,sigblur,s2n)
 100  format(f10.5,2x,f10.5,2x,f8.5,2x,f8.5)
      zstart=(wstart/1215.67)-1.
! ---------------------------------------------------------------------
! GET THE SPLINE
! ---------------------------------------------------------------------
      call spline(npoints,xs,ys,CDDF)
! ---------------------------------------------------------------------
! GENERATE ARTIFICIAL SPECTRA
! ---------------------------------------------------------------------
      do i=25,25
         zend=zqso(i)
         write (descriptor,'(I6.6)') i
         write (6,*)'=================================================='
         write (6,*)'               Spectrum no. ',descriptor
         write (6,*)'=================================================='

         call power_laws(npoints,zstart,zqso(i),xs,ys,CDDF,
     +                      bigA,gamma,corr,nl,ni)
         allocate(nhi4(nl))
         allocate(z4(nl))
         call assign(npoints,zstart,zqso(i),xs,CDDF,gamma,nl,ni,nhi4,z4)
         call qsosim9(zqso(i),alpha(i),vmag(i),wstart,wend,dw,nc,nuplim,
     +          sigblur(i),s2n(i),inoise,dvavoid,npts,lambda,flux,flerr,
     +          nnflux,npoints,nl,ni,nhi4,z4)
         outfile='spec-'//descriptor//'.fits'
         call writefits(outfile,ra(i),dec(i),zqso(i),alpha(i),npts,
     &                     lambda,flux,flerr,nnflux)
         write (*,*)'--------------------------------------------------'
      end do
! ---------------------------------------------------------------------
! PLOT QSO SPECTRUM OF THE LAST SOURCE
! --------------------------------------------------------------------- 
c      call PGBEGIN (0,'/xserve',1,2)
c      call PGSLW(1)
c      call PGENV (3550.,10500.,0.0,5e-15,0,1)
c      call PGLABEL ('lambda','flux','QSO spectrum')
c      call pgline(npts,real(lambda),real(flux))
c      call pgsci(2)
c      call pgline(npts,real(lambda),real(nnflux))
c      call pgsci(1)
c      call PGENV(11.5,22.5,2.00,2.1,0,1)
c      call PGLABEL('log NHI','z','Random choice of NHI & redshift')
c      call PGPT(nl,real(nhi4),real(z4),3)
c      call PGEND
!======================================================================
      END PROGRAM qsosim10
!======================================================================
