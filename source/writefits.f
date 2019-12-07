!=======================================================================
      SUBROUTINE writefits(home,folder,outfile,ra,dec,zqso,z_w,
     :                     alpha_nu,alpha,rmag,name,id,iplate,imjd,
     :                     ifiber,ipix,wstart,mags,lambda,flux,flerr,
     :                     ivar,nnflux,ncflux,noabs,nlin)
c     PURPOSE: create a fits file with parameters returned by qsosim9
c     OUTPUT:  spectrum.fits file with lambda,flux,sigma,no-noise-flux
!=======================================================================
*      input, character          :  home         path to home folder
*                                   folder       path to the folder where
*                                                the fits file will be saved
*                                   outfile      name of the output fits file
*                                   name         SDSS name of the object
*      input, real*8             :  ra           RA coordinate of QSO
*                                   dec          DEC coordinate of QSO
*                                   zqso         redshift of QSO
*                                   alpha_nu     spectral index, as per SDSS
*                                   alpha        spectral index, as per 
*                                                Sebastian Schmidt's program
*                                   rmag         SDSS r magnitude

*      input, integer            :  z_w          redshift warning flag
*                                   id           SDSS id
*                                   iplate       SDSS plate
*                                   imjd         SDSS mjd
*                                   ifiber       SDSS fiber
*                                   ipix         number of pixels in the spectrum
*                                   nline        number of lines in the spectrum
*      input, real*8 array(5)    :  mags         SDSS u,g,r,i,z magnitudes
*      input, real*8 array(npix) :  loglam       wavelengths in log
*                                   flux         flux in 10e-17 erg s-1 cm-2 A-1
*                                   noise        noise array
*                                   nnflux       no-noise flux
*                                   flux_nc      not convolved flux
*                                   noabs        flux with no absorption
!=======================================================================

c GENERAL DECLARATIONS
      IMPLICIT NONE
c Declare variables
      INTEGER :: unit,i,status,blocksize,bitpix,naxis,naxes
      INTEGER :: fpixels(2),lpixels(2),nax(2)
      INTEGER :: nrows, tfields, varidat, colnum, idum, inoise
      INTEGER :: z_w,id,iplate,imjd,ifiber,ipix,nlin
      CHARACTER :: name*18
      REAL*8,dimension(ipix) :: lambda,flux,flerr,ivar,
     &                                     nnflux,ncflux,noabs
      REAL*8 :: ra,dec,zqso,alpha,rmag,
     &                     alpha_nu,wstart,mags(5)
      REAL*8,dimension(3000)::cd,z,b
      CHARACTER::atom(3000)*2,ion(3000)*4
      CHARACTER*20 :: ttype1(15), tunit1(15), tform1(15)
      CHARACTER*20 :: ttype2(7), tunit2(7), tform2(7)
      CHARACTER*20 :: ttype3(5), tunit3(5), tform3(5)
      CHARACTER*30 :: errtext, extname 
      CHARACTER :: home*120,folder*120,outfile*35,mockfolder*220
      LOGICAL :: simple, extend

      type linelist
         SEQUENCE
         character*2 atom
         character*4 ion
         real*8 colden
         real*8 rdf
         real*8 bpar
      end type
      type (linelist) :: llist(3000)
      common /linelist/llist
c Define parameters, values same as in CFITSIO cookbook
      blocksize=2
      status=0
      simple=.true.
      bitpix=16
      naxis=0
      naxes=0
      extend=.true.
c-------------------------------------------------------------------------------
      call chdir(folder)
*     Check if FITS file exists and delete if does. Create a new FITS file and 
*     initialize it. 
      call deletefile(outfile,status)
      call FTGIOU(unit,status)
 555  call FTINIT(UNIT,outfile,blocksize,status)
      if (status.eq.0)then 
c         print *,status,' Output file initialized'
      else 
         print *,status,' ',errtext
      end if
c     Define primary array parameters
      call FTPHPR(UNIT,simple,bitpix,naxis,naxes,0,1,extend,status)
      call ftgerr(status,errtext)
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
*     Define data to be inputted into fits file: type, format & units
      DATA ttype1/'RA','DEC','Z_VI','Z_WARNING','ALPHA_NU','ALPHA_FIT',
     &            'R_MAG','SDSS_name','THING_ID','PLATE','MJD','FIBER',
     &            'NPIX','BEGIN_WAVE','PSFMAG'/
      DATA ttype2/'LOGLAM','FLUX','NOISE','IVAR','NNFLUX','NCFLUX',
     &            'NO_ABS'/
      DATA ttype3/'ATOM','ION','COLDEN','Z','B'/
      
      DATA tform1/'D','D','D','J','D','D','D','18A','J','J','J',
     &            'J','J','D','5E'/
      DATA tform2/'D','D','D','D','D','D','D'/
      DATA tform3/'2A','4A','D','D','D'/
      DATA tunit1/'','','','','','','','','','','','','','',''/
      DATA tunit2/'log10(A)','1e-17 erg/s/cm^2/A','1e-17 erg/s/cm^2/A',
     +  '1e-17 erg/s/cm^2/A','1e-17 erg/s/cm^2/A','1e-17 erg/s/cm^2/A',
     +  '1e-17 erg/s/cm^2/A'/
      DATA tunit3/'','','cm^-2','','km/s'/
*     Create the first binary table HDU
      nrows=1
      tfields=15
      varidat=0
      extname='QSO_INFO'
      call FTIBIN(unit,nrows,tfields,ttype1,tform1,tunit1,
     &            extname,varidat,status)
*     Put values in columns in the first binary table
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

*     Create the second binary table HDU with data pertaining the QSO
      nrows=ipix
      tfields=7
      extname='QSO_SPECTRUM'
      call FTIBIN(unit,nrows,tfields,ttype2,tform2,tunit2,
     &            extname,varidat,status)
      call FTPCLD(unit,1,1,1,nrows,lambda,status)
      call FTPCLD(unit,2,1,1,nrows,flux,status)
      call FTPCLD(unit,3,1,1,nrows,flerr,status)
      call FTPCLD(unit,4,1,1,nrows,ivar,status)
      call FTPCLD(unit,5,1,1,nrows,nnflux,status)
      call FTPCLD(unit,6,1,1,nrows,ncflux,status)
      call FTPCLD(unit,7,1,1,nrows,noabs,status)
      call FTPDAT(unit,status)

*     Create the third binary table HDU
      nrows=nlin
      tfields=5
      varidat=0
      extname='LINE_LIST'
      call FTIBIN(unit,nrows,tfields,ttype3,tform3,tunit3,
     &            extname,varidat,status)
*     Put values in columns in the first binary table
      atom=llist%atom
      ion=llist%ion
      cd=llist%colden
      z=llist%rdf
      b=llist%bpar
      call FTPCLS(unit,1,1,1,nrows,atom,status)
      call FTPCLS(unit,2,1,1,nrows,ion,status)
      call FTPCLD(unit,3,1,1,nrows,cd,status)
      call FTPCLD(unit,4,1,1,nrows,z,status)
      call FTPCLD(unit,5,1,1,nrows,b,status)

*     Close the FITS file
      call FTCLOS(unit,status)
      call FTFIOU(unit,status)
      call ftgerr(status,errtext)
      if (status.eq.0)then 
         write(6,*) 'QSO data & the spectrum saved in a FITS file.'
      else 
         print *,status,' ',errtext
      end if
*     End program
      call chdir(home)
      RETURN
      END SUBROUTINE writefits

c-------------------------------------------------------------------------------
c SUBROUTINES
      subroutine deletefile(filename,status)

C     A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*30 filename

C     simply return if status is greater than zero
      if (status .gt. 0)return
C     Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)
C     try to open the file, to see if it exists
      call ftopen(unit,filename,0,blocksize,status)
      if (status .eq. 0)then
C         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq.104)then
C         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      else
C         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if
C     free the unit number for later reuse
      call ftfiou(unit, status)
      return
      end subroutine
c-------------------------------------------------------------------------------
