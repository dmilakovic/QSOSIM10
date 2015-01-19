!=======================================================================
      subroutine write_cloudy_input(home,descriptor,z,nl,nhi4,z4,s,mask)
!=======================================================================
*      write files to be used by CLOUDY to calculate metal absorption

      
*      input, character           : home         path to the home folder
*                                   descriptor   unique QSO designation
*      input, real*8              : z            QSO's redshift
*      input, integer             : nl           number of H I lines
*      input, real*8 array(nl)    : nhi4         H I col densities
*                                 : z4           H I redshifts

*      output, integer            : s            number of H I lines with N_HI > 1e15 cm-2
*      output, logical array(nl)  : mask         mask of H I lines with N_HI > 1e15 cm-2
*                                                (T if N_HI > 1e15 cm-2)
!=======================================================================
      implicit none
      integer :: i,j,t,s,idum
      integer, intent(in) :: nl
      character :: filename*90, line*100, cwd*30, home*120, folder*90
      character ::descriptor*15,num*4,makedirectory*100,lognhi*10,zed*10
      character :: met*5
      real*8 :: m,c,d,gasdev3
      real*8,dimension(nl) :: nhi4, z4
      logical,dimension(nl) :: mask
      real*8,intent(in) :: z
      logical ex
      external gasdev3
! ----------------------------------------------------------------------
*     Create folder with object name and create input files for the
*     absorption systems
! -----------------------------------------------------------------------
      folder = trim(trim(home)//'/cloudy_in')
      inquire(file=folder,exist=ex)
      if (ex.eqv..true.) goto 102
      makedirectory = 'mkdir '//folder
      call system(makedirectory)
 102  folder = trim(trim(home)//'/cloudy_in/spec_'//descriptor)
      inquire(file=folder,exist=ex)
      if (ex.eqv..true.) goto 202
      makedirectory = 'mkdir '//trim(folder)
      call system(makedirectory)
 202  call chdir(folder)
      folder = trim(home)//'/cloudy_in/spec_'//descriptor
*     Initialize the counters and set the mask to be false for every 
*     H I absorption system. 's' counts the number of systems that satisfy
*     the masking condition. 't' counts the position from 1 to nl (num of 
*     H I systems).
      s=0; t=0
      do i=1,nl
         mask(i)=.false.
      end do
*     Using QSO's redshift, determine the cosmic metallicity m(z):
*     initialise random numbers
      idum=time()
      c = gasdev3(idum); d = gasdev3(idum)
*     determine metallicity (Rafelski 2012)
      m = -(0.22+0.03*c)*z - (0.65+0.09*d)
      write(6,*) 'METALLICITY == ',m,'c,d',c,d
*     For each H I absorption system, check whether the column density is
*     greater than 1e17 cm-2 and increase counter s by one & set mask=true.
*     Create a new file for that system and write CLOUDY input in it. 
      do i=1,nl
         if (nhi4(i).gt.17.0) then
            s=s+1
            mask(i)=.true.
            write(num,'(i4.4)') s
            write(zed,'(f5.3)') z4(i)
            write(lognhi,'(f5.2)') nhi4(i)
            write(met,'(f5.2)') m
            write(6,*) met
            filename = trim(folder)//'/system_'//num//'.in'
            open (unit=17,file=filename,status='replace')
 100        format (a70)
             line = 'c Lyman alpha system'; write (17,100) line
             line = 'cmb z='//trim(zed); write (17,100) line
             line = 'table hm05 z='//trim(zed); write (17,100) line
             line = 'hden -2';write(17,100) line
             line = 'metals '//trim(met); write(17,100) line
             line = 'stop neutral column density '//trim(lognhi)
                      write (17,100) line
             line = 'double'; write(17,100) line
             line = 'iterate to convergence'; write (17,100) line 
             line = 'print line faint -1'; write(17,100) line
             line = 'save overview last "system_'//num//'.ovr"'
                      write(17,100)line
             line = 'save results last "system_'//num//'.rlt"'
                      write(17,100)line
*           Close file
            close (unit=17)
         else
c            write (6,*) i, nhi4(i), 'x'
c            t=t+1
         end if
      end do
      write (6,'(5x,i4,3x,a20)') s,'clouds in the sky!'
      call chdir(home)
      return
      end subroutine
