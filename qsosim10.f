!======================================================================
      PROGRAM qsosim10
c  PURPOSE: call other subroutines to generate artificial spectra!
c  OUTPUT:  artificial SDSS catalogue
!======================================================================
      IMPLICIT NONE
      CHARACTER :: infile*20, outfile*35, descriptor*15,command*20
      CHARACTER :: home*120,gohome*120,mockspec_folder*120,name*18
      INTEGER,PARAMETER :: npoints=10000,nrows=100000
      integer, parameter :: RegInt_K = selected_int_kind (12)
      integer (kind=RegInt_K) :: i
      INTEGER :: nr,j,inoise, numlin, npts, nl,s,option,zstart_flag
      INTEGER :: z_w,id,iplate,imjd,ifiber,ipix,status
      REAL*8 :: wstart,wend,dw,mags(5)
      REAL*8 :: nc,nuplim,sigblur
      REAL*8 :: corr, zstart,zend
      REAL*8,DIMENSION(3) :: bigA,gamma,n
      INTEGER,DIMENSION(3) :: ni
      CHARACTER,DIMENSION(nrows) :: SDSS_name*18
      INTEGER,DIMENSION(nrows) :: thing_id,plate,mjd,
     &                                   fiber,z_flag,npix
      REAL*8, DIMENSION(nrows) :: ra,dec,zqso,alpha,alpha_fit,
     &                                   begin_wave,rmag
      REAL*8,DIMENSION(nrows,5) :: psfmag
      REAL*8,DIMENSION(262144) :: loglam,flux
      REAL*4,DIMENSION(262144) :: pg_loglam, pg_flux
      REAL*8,DIMENSION(262144) :: noise, nnflux,flux_nc, noabs
      REAL*8, DIMENSION(npoints) :: xs,ys,CDDF,H
      REAL*8,DIMENSION(4000) :: nhi4,z4
      REAL*8,DIMENSION(:),allocatable :: nciv4, novi4
      LOGICAL,DIMENSION(4000)::mask
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
      !   [17.0,22.0]      log(0.175)   1.33

      data bigA/33.11311,5.24807,0.1750/
      data gamma/1.51,2.16,1.33/

! ---------------------------------------------------------------------
! READ INPUT PARAMETERS
! ---------------------------------------------------------------------
      home = '/Users/dm/Documents/GitHub/QSOSIM10'
      gohome = trim('cd '//home)
      mockspec_folder = trim(home)//'/mockspectra'
      write (6,*)'=================================================='
      write (6,*)'                    QSOSIM 10                     '
      write (6,*)'=================================================='
      write (6,*)'What would you like to create?'
      write (6,50)'(1)','SDSS'
      write (6,50)'(2)','random (legacy)'
      write (6,*)'=================================================='
      write (6,'(1x,a13)',advance='no')'I choose... '
c      read (*,'(i1)') option
      option=1
      if (option.eq.1) then
         infile='DR10Q_r.fits'
         call read_nrows(home,infile,nr)
         write(6,*) nrows
         call SDSS_readfits(home,infile,nrows,nr,SDSS_name,RA,DEC,
     &         thing_id,plate,mjd,fiber,zqso,z_flag,alpha,alpha_fit,
     &         npix,begin_wave,psfmag,rmag)
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
      call spline(npoints,nc,nuplim,xs,ys,CDDF)
! ---------------------------------------------------------------------
! GENERATE ARTIFICIAL SPECTRA
! ---------------------------------------------------------------------
      do i=1,1
         write (descriptor,"(i4.4,'-',i5.5,'-'i4.4)")
     &                     plate(i),mjd(i),fiber(i)
         write (6,*)'=================================================='
         write (6,*)'            Spectrum no. ',descriptor
         write (6,*)'=================================================='
         write (6,*) i,'th iteration'
         !**** QSO data
         wstart=begin_wave(i) 
         wend=begin_wave(i)+npix(i)*dw 
         zstart=(10**begin_wave(i)/1215.67)-1.
         zend=zqso(i)
         if (zstart.ge.zend) then
            zstart_flag=zstart_flag+1
            write(6,*) 'zstart >= zend!'; goto 99
         end if
         ipix=npix(i)         
         !**** calculate the number of lines between zstart and zqso ****
         call power_laws(npoints,zstart,zend,xs,ys,CDDF,
     +                      bigA,gamma,nl,ni)
         !**** assign zero values of nhi and z to lines ****
         do j=1,2000
            nhi4(j)=0.0
            z4(j)=0.0
         end do
         call assign(npoints,zstart,zend,xs,CDDF,nl,ni,nhi4,z4)
         !**** write files used by cloudy, call cloudy and read results ****
         call write_cloudy_input(home,descriptor,zend,nl,nhi4,s,mask)
c         call cloudy(home,descriptor,nl,s)
         allocate(novi4(nl)); 
c         allocate(nciv(s)); 
         call read_cloudy_output(home,descriptor,'O ','VI  ',nl,s,mask,
     &                           novi4)
         do j=1,nl
            write (6,*) 'j,mask,nhi,novi',j,mask(j),nhi4(j),novi4(j)
         end do
         !**** generate artificial spectrum ****
         call qsosim(zqso(i),alpha_fit(i),rmag(i),begin_wave(i),dw,
     +         sigblur,npix(i),nl,nhi4,z4,s,novi4,loglam,flux,
     +         noise,nnflux,flux_nc,noabs)
         !**** write qsosim output + qso general data into a fits file ****
         outfile='mockspec-'//descriptor//'.fits'
         do j=1,5
            mags(j)=psfmag(i,j)
         end do
         call writefits(home,mockspec_folder,outfile,ra(i),dec(i),
     +                  zqso(i),z_flag(i),alpha(i),alpha_fit(i),rmag(i),
     +                  SDSS_name(i),thing_id(i),plate(i),mjd(i),
     +                  fiber(i),npix(i),wstart,mags,loglam,flux,noise,
     +                  nnflux,flux_nc,noabs)
         write (6,*)'--------------------------------------------------'
         do j=1,262144
            pg_loglam(j)=real(loglam(j))
            pg_flux(j)=real(flux(j))
         end do
         !**** make sure that there are no artefacts in the next iteration ****
         do j=1,262144
            loglam(j)=0.0
            flux(j)=0.0
            noise(j)=0.0
            flux_nc(j)=0.0
            nnflux(j)=0.0
         end do
 99      continue
      end do
      write (6,*)'=================================================='
      write (6,*)'               DATA SET COMPLETED!                '
      write (6,*)'=================================================='
      write (6,*)'QSOs with zstart>zend :',zstart_flag
! ---------------------------------------------------------------------
! PLOT QSO SPECTRUM OF THE LAST SOURCE
! --------------------------------------------------------------------- 
c      call PGBEGIN (0,'/xserve',1,1)
c      call PGENV (3.55,4.02,0.0,20.0,0,1)
c      call PGLABEL('log NHI','z','Random selection of NHI-z')
c      call pgpt(nl,real(nhi4),real(z4),2)
c      call pgline(ipix,pg_loglam,pg_flux)
c      call PGEND
c      write (6,*) 'nothing else here!'
! --------------------------------------------------------------------- 
! FREE ALLOCATED MEMORY
! --------------------------------------------------------------------- 
c      deallocate(xs);        deallocate(ys);       deallocate(CDDF)
c      deallocate(SDSS_name); deallocate(thing_id); deallocate(plate)
c      deallocate(mjd);       deallocate(fiber);    deallocate(z_flag)
c      deallocate(npix);      deallocate(ra);       deallocate(dec)
c      deallocate(zqso);      deallocate(alpha);    deallocate(alpha_fit)
c      deallocate(begin_wave);deallocate(rmag);     deallocate(psfmag)
!======================================================================
      END PROGRAM qsosim10
!======================================================================
