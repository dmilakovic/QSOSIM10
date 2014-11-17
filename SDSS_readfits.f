c-------------------------------------------------------------------------------
      SUBROUTINE SDSS_readfits(home,infile,nrows,SDSS_name,RA,DEC,
     &         thing_id,plate,mjd,fiber,zqso,z_flag,alpha,alpha_fit,
     &         npix,begin_wave,psfmag,rmag)
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
      INTEGER :: nrows,naxis,naxes(2),fpixels(2),lpixels(2),incs(2)
      INTEGER :: felem,nelems,nullj,nfound,irow,colnum
      REAL :: nulle
      REAL*8 :: wstart,wend,dw
      REAL*8 :: nc, nuplim,dvavoid,array(5)
      CHARACTER,DIMENSION(nrows) :: SDSS_name*18
      INTEGER,DIMENSION(nrows) :: thing_id,plate,mjd,
     &                                   fiber,z_flag,npix
      REAL*8, DIMENSION(nrows) :: ra,dec,zqso,alpha,alpha_fit,
     &                                   begin_wave,rmag
      REAL*8, DIMENSION(nrows,5) :: psfmag
      character :: infile*20,ttype(20)*10!, nhills,blls,zlls
      logical :: anynull
      character :: errtext*30,card*50,comment*30,home*120,nulls
c------------------------------------------------------------------------------
c Define parameters
      status=0
      unit=3
      readwrite=0
c------------------------------------------------------------------------------
c Open fits file to read
      call chdir(home)
      call ftopen(unit,infile,readwrite,blocksize,status)
      call ftgerr(status,errtext)
      if (status.eq.0) then 
         write (*,*)status,'File '//infile//' opened'
      end if
c------------------------------------------------------------------------------
c Read contents of 'NoName' (ntable=2) binary table
c Read data from columns
      ntable=2
      call ftmahd(unit,ntable,hdutype,status)
      call FTGKYJ(unit,'NAXIS2',nrows,comment,status)
c      call FTGKYD(unit,'wstart',wstart,comment,status)
c      call FTGKYD(unit,'wend',wend,comment,status)
c      call FTGKYD(unit,'dw',dw,comment,status)
c      call FTGKYD(unit,'nc',nc,comment,status)
c      call FTGKYD(unit,'nuplim',nuplim,comment,status)
c      call FTGKYJ(unit,'inoise',inoise,comment,status)
c      call FTGKYD(unit,'dvavoid',dvavoid,comment,status)
 100  format(2x,a10,d9.3,a4)
 125  format(2x,a10,f9.3,a4)
 150  format(2x,a10,i3)
c      write (6,*) 'Number of objects :',nrows 

c------------------------------------------------------------------------------
c      if (status.eq.0)then 
c         call FTGERR(status, errtext)
c         print *,status, errtext
c      else 
c         call FTGERR(status, errtext)
c         print *,status,' ',errtext
c      end if
c Read column data, one row at a time, and print them out
      felem=1
      nelems=1
      nulle=0.      
      nullj=0
      nulls=''
      naxis=1
      naxes=(/5*nrows,1/) 
      incs=(/1,1/)
c      write (*,300)'RA','DEC','Z','alpha'
      do irow=1,nrows
         call FTGCVS(unit,1,irow,felem,nelems,nulls,SDSS_name(irow),
     &        anynull,status)
         call FTGCVD(unit,2,irow,felem,nelems,nulle,ra(irow),
     &       anynull,status)
         call FTGCVD(unit,3,irow,felem,nelems,nulle,dec(irow),
     &       anynull,status)
         call FTGCVJ(unit,4,irow,felem,nelems,nullj,thing_id(irow),
     &       anynull,status)
         call FTGCVJ(unit,5,irow,felem,nelems,nullj,plate(irow),
     &       anynull,status)
         call FTGCVJ(unit,6,irow,felem,nelems,nullj,mjd(irow),
     &       anynull,status)
         call FTGCVJ(unit,7,irow,felem,nelems,nullj,fiber(irow),
     &       anynull,status)
         call FTGCVD(unit,8,irow,felem,nelems,nulle,zqso(irow),
     &       anynull,status)
         call FTGCVJ(unit,9,irow,felem,nelems,nullj,z_flag(irow),
     &       anynull,status)
         call FTGCVD(unit,10,irow,felem,nelems,nulle,alpha(irow),
     &       anynull,status)
         call FTGCVD(unit,11,irow,felem,nelems,nulle,alpha_fit(irow),
     &       anynull,status)
         call FTGCVJ(unit,12,irow,felem,nelems,nullj,npix(irow),
     &       anynull,status)
         call FTGCVD(unit,13,irow,felem,nelems,nulle,begin_wave(irow),
     &       anynull,status)
         fpixels=(/1,irow/)
         lpixels=(/5,irow/)
         call FTGSVD(unit,14,naxis,naxes,fpixels,lpixels,incs,nulle,
     &       array,anynull,status)
         do j=1,5
            psfmag(irow,j)=array(j)
         end do
         rmag(irow) = array(3)
c         write (*,200)irow,ra(irow),dec(irow),zqso(irow),alpha_fit(irow)
c     & ,rmag(irow),npix(irow)
c         write (6,*)irow,SDSS_name(irow),plate(irow),mjd(irow),
c     & fiber(irow),begin_wave(irow)
c         write (*,250) 'mags',array
      end do
 200  format (i6,4x,f10.5,2x,f10.3,2x,f10.5,2x,f10.5,2x,f10.5,2x,i5)
 250  format (2x,a4,4x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5)
 300  format (6x,a10,2x,a10,2x,a10,2x,a6,2x,a7)
cc------------------------------------------------------------------------------
cc Close fits file
      call ftclos(unit,status)
      if (status.eq.0)then 
         print *,status,' File closed'
      else 
         call FTGERR(status, errtext)
         print *,status,' ',errtext
      end if
c------------------------------------------------------------------------------
      RETURN
      END SUBROUTINE SDSS_readfits
c------------------------------------------------------------------------------
