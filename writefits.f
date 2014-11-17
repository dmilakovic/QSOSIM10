c-------------------------------------------------------------------------------
      SUBROUTINE writefits(outfile,ra,dec,zqso,z_w,alpha_nu,alpha,rmag,
     &                     name,id,iplate,imjd,ifiber,ipix,wstart,mags,
     &                     npts,lambda,flux,flerr,nnflux,flux_nc)
c     PURPOSE: create a fits file with parameters returned by qsosim9
c     OUTPUT:  spectrum.fits file with lambda,flux,sigma,no-noise-flux
c-------------------------------------------------------------------------------
c GENERAL DECLARATIONS
      IMPLICIT NONE
c Declare variables
      INTEGER :: unit,i,j,status,blocksize,bitpix,naxis,naxes,npts
      INTEGER :: fpixels(2),lpixels(2),nax(2),nullj
      INTEGER :: nrows, tfields, varidat, colnum, idum, inoise
      INTEGER :: z_w,id,iplate,imjd,ifiber,ipix
      CHARACTER :: name*18,nulls
      REAL*8,dimension(npts) :: lambda,flux,flerr,nnflux,flux_nc
      REAL*8 :: ra,dec,zqso,alpha,rmag,alpha_nu,wstart,mags(5),nulle
      CHARACTER*20 :: ttype1(15), tunit1(15), tform1(15)
      CHARACTER*20 :: ttype2(5), tunit2(5), tform2(5)
      CHARACTER*30 :: errtext,outfile, extname 
      LOGICAL :: simple, extend
c      ALLOCATE(lambda(npts),flux(npts),flerr(npts),nnflux(npts))
c Define parameters
      blocksize=2
      status=0
      simple=.true.
      bitpix=16
      naxis=0
      naxes=0
      extend=.true.
c-------------------------------------------------------------------------------
      call ftgiou(unit,status)
c Check if fits file exists and delete if does
      call DELETEFILE(outfile,status)
      call ftgerr(status,errtext)
      write (6,*) 'writefits'
      print *,status,' ',errtext
c Create a fits file
      call FTINIT(UNIT,outfile,blocksize,status)
      if (status.eq.0)then 
         print *,status,' Output file initialized'
      else 
         print *,status,' ',errtext
      end if
c Define primary array parameters
      call FTPHPR(UNIT,simple,bitpix,naxis,naxes,0,1,extend,status)
      call ftgerr(status,errtext)
      print *,status,' ',errtext
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c Define data to be inputted into fits file
      DATA ttype1/'RA','DEC','Z_VI','Z_WARNING','ALPHA_NU','ALPHA_FIT',
     &            'R_MAG','SDSS_name','THING_ID','PLATE','MJD','FIBER',
     &            'NPIX','BEGIN_WAVE','PSFMAG'/
      DATA ttype2/'LOGLAM','FLUX','IVAR','NNFLUX','FLUX_NC'/
      
      DATA tform1/'D','D','D','J','D','D','D','18A','J','J','J',
     &            'J','J','D','5E'/
      DATA tform2/'D','D','D','D','D'/
      DATA tunit1/'','','','','','','','','','','','','','',''/
      DATA tunit2/'log','1e-17 erg/s/cm^2/A','','',''/
c Create the first binary table HDU
      nrows=1
      tfields=15
      varidat=0
      nulle=0.0
      nullj=0
      nulls=''
      extname='General'
      call FTIBIN(unit,nrows,tfields,ttype1,tform1,tunit1,
     &            extname,varidat,status)
      write (6,*) 'inputting data'
c Rename header column names and assign values
      call FTPCLD(unit,1,1,1,nrows,ra,status)
      call FTPCLD(unit,2,1,1,nrows,dec,status)
      call FTPCLD(unit,3,1,1,nrows,zqso,status)
      call FTPCLJ(unit,4,1,1,nrows,z_w,status)
      call FTPCLD(unit,5,1,1,nrows,alpha_nu,status)
      call FTPCLD(unit,6,1,1,nrows,alpha,status)
      call FTPCLD(unit,7,1,1,nrows,rmag,status)
      call FTPCLS(unit,8,1,1,nrows,name,status)
      call FTPCLJ(unit,9,1,1,nrows,id,status)
      call FTPCLJ(unit,10,1,1,nrows,iplate,status)
      call FTPCLJ(unit,11,1,1,nrows,imjd,status)
      call FTPCLJ(unit,12,1,1,nrows,ifiber,status)
      call FTPCLJ(unit,13,1,1,nrows,ipix,status)
      call FTPCLD(unit,14,1,1,nrows,wstart,status)
      call FTPCLD(unit,15,1,1,5,mags,status)
      naxis=1
      nax=(/5,1/)
c      call FTPTDM(unit,15,naxis,nax,status)
c      call FTMKYD(unit,'psfmag',mags,5,'',status)
      fpixels=(/1,1/)
      lpixels=(/5,1/)
c      call FTPPRD(unit,15,1,5,mags,status)
c      call FTPSSD(unit,15,naxis,nax,fpixels,lpixels,mags,status)
      
c      call ftgerr(status,errtext)
c      print *,status,' ',errtext
c Create the second binary table HDU with data pertaining each QSO
      nrows=npts
 250  format(f10.5,2x,f10.5,2x,f10.5,2x,f10.5)
c      write (*,250)flerr(1:10)
      tfields=5
      extname='QSO'
      call FTIBIN(unit,nrows,tfields,ttype2,tform2,tunit2,
     &            extname,varidat,status)
      call FTPCLD(unit,1,1,1,nrows,lambda,status)
      call FTPCLD(unit,2,1,1,nrows,flux,status)
      call FTPCLD(unit,3,1,1,nrows,flerr,status)
      call FTPCLD(unit,4,1,1,nrows,nnflux,status)
      call FTPCLD(unit,5,1,1,nrows,flux_nc,status)
c Close fits file
      call FTCLOS(unit,status)
      call ftgerr(status,errtext)
      if (status.eq.0)then 
         print *,status,' File closed'
      else 
         print *,status,' ',errtext
      end if
c End program
      RETURN
      END SUBROUTINE writefits

c-------------------------------------------------------------------------------
c SUBROUTINES
      subroutine deletefile(outfile,status)
C  A simple little routine to delete a FITS file
      integer status,unit,blocksize
      character*30 outfile
      character*30 errtext
C  Simply return if status is greater than zero
      if (status .gt. 0)return
C  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)
C  Try to open the file, to see if it exists
      call ftopen(unit,outfile,1,blocksize,status)
      if (status .eq. 0)then
C         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq. 103)then
C         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
      else
C         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if
C  Free the unit number for later reuse
      call ftfiou(unit, status)
      end subroutine deletefile
c-------------------------------------------------------------------------------
