      subroutine write_cloudy_input(home,descriptor,z,nl,nhi4,s)
! ----------------------------------------------------------------------
      implicit none
      integer :: i,j,t
      integer ::s
      integer, intent(in) :: nl
      character :: filename*90, line*100, cwd*30, home*120, folder*90
      character ::descriptor*15,num*4,makedirectory*100,lognhi*10,zed*10
      real*8,dimension(nl),intent(inout) :: nhi4
      integer,dimension(nl) :: mask
      real*8,intent(in) :: z
      logical ex
! ----------------------------------------------------------------------
! Create folder with object name and create input files for the absorption systems
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
      s=0; t=0
      do i=1,nl
         mask(i)=0
      end do
      do i=1,nl
c         write (6,*) i,nhi4(i)
         if (nhi4(i).gt.14.0) then
            s=s+1
            mask(i)=1
            write(num,'(i4.4)') s
            write(zed,'(f5.3)') z
            write(lognhi,'(f5.2)') nhi4(i)
            filename = trim(folder)//'/system_'//num//'.in'
c            write (6,*) i, nhi4(i), '_/'
            open (unit=17,file=filename,status='replace')
 100        format (a70)
             line = 'c Lyman alpha system'; write (17,100) line
             line = 'cmb z='//trim(zed); write (17,100) line
             line = 'table hm05 z='//trim(zed); write (17,100) line
             line = 'hden -2';write(17,100) line
             line = 'metals -1.5'; write(17,100) line
             line = 'stop neutral column density '//trim(lognhi)
                      write (17,100) line
             line = 'double'; write(17,100) line
             line = 'iterate to convergence'; write (17,100) line 
             line = 'print line faint -1'; write(17,100) line
c             line = 'save performance "igm_lalpha.per"'
c                      write (17,100) line
             line = 'save overview last "system_'//num//'.ovr"'
                      write(17,100)line
c             line = 'save dr last "igm_lalpha.dr"';write (17,100) line
             line = 'save results last "system_'//num//'.rlt"'
                      write(17,100)line
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
