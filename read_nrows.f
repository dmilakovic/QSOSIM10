c-------------------------------------------------------------------------------
      SUBROUTINE read_nrows(home,infile,nrs)!,SDSS_name,RA,DEC,
!     &         thing_id,plate,mjd,fiber,zqso,z_flag,alpha,alpha_fit,
!     &         npix,begin_wave,rmag)
c     PURPOSE: read fits file containing data for qsosim9
c     INPUT: filename  name of the fits file
c     OUTPUT: 
*            wstart    starting wavelength
*            wend      ending wavelength
*            dw        pixel size
*            nc        N(HI) lower cut-off
*            nuplim    N(HI) upper cut-off
*            inoise    inoise parameter
*            dvavoid   avoidance zone around each additional system in km/h
*            ra        RA coordinates (array)
*            dec       DEC coordinates (array)
*            zqso      Z (array)
*            alpha     spectral index (array)
*            vmag      V magnitude (array)
*            sigblur   spectral resolution (array)
c------------------------------------------------------------------------------
c GENERAL DECLARATIONS
c Declare variables
      INTEGER :: status,unit,readwrite,blocksize,hdutype,ntable,inoise
      INTEGER :: naxis,naxes(2)
      INTEGER :: felem,nelems,nullj,nfound,irow,colnum
      INTEGER :: nrs
      REAL :: nulle
      character :: infile*20,ttype(20)*10!, nhills,blls,zlls
      logical :: anynull
      character :: errtext*30,card*50,comment*30,home*120,nulls
c------------------------------------------------------------------------------
c Define parameters
      status=0
      call ftgiou(unit,status)
      readwrite=0
c------------------------------------------------------------------------------
c Open fits file to read
      call chdir(home)
      call ftopen(unit,infile,readwrite,blocksize,status)
      call ftgerr(status,errtext)
c      if (status.eq.0) then 
c         write (*,*)status,'File '//infile//' opened'
c      end if
c------------------------------------------------------------------------------
c Read contents of 'NoName' (ntable=2) binary table
c Read data from columns
      ntable=2
      call ftmahd(unit,ntable,hdutype,status)
      call FTGKYJ(unit,'NAXIS2',nrs,comment,status)
 100  format(2x,a10,d9.3,a4)
 125  format(2x,a10,f9.3,a4)
 150  format(2x,a10,i3)
c      write (6,*) 'Number of objects :',nrs

cc------------------------------------------------------------------------------
cc Close fits file
      call ftclos(unit,status)
      call ftfiou(unit,status)
      if (status.eq.0)then 
         return
      else 
         call FTGERR(status, errtext)
         print *,status,' ',errtext
      end if
c------------------------------------------------------------------------------
      stop
      END SUBROUTINE read_nrows
c------------------------------------------------------------------------------
