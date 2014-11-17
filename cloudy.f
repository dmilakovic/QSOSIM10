      subroutine cloudy(home,descriptor,mask,nl,s)

      integer i,s
      character :: descriptor*16,num*4,ns*4
      character*120 :: path,home,folder,filename
      character command*200
      integer, dimension(nl) :: mask
      logical ex

      write (6,*) '--------------------------------------------------'
      write (6,*) '              The sky is CLOUDY                   '
      write (6,*) '--------------------------------------------------'
      path='/Users/dm/Programming/cloudy/c13.03/source/cloudy.exe'
      folder = trim(home)//'/cloudy_in/spec_'//descriptor//'/'
      write (ns,'(i4.4)') s
c      write (6,*) 'folder =', folder
      call chdir(folder)
      do i=1,1
         if (mask(i).eq.1) then
            write(num,'(i4.4)') i
            filename = trim('system_'//num//'.in')
c     write (6,*) filename
            inquire(file=filename, exist=ex)
c         write (6,*) ex
            if (ex.eqv..true.) then
               filename = trim('system_'//num)
               command = trim(path)//' -r '//filename
c            write (6,*) command
               write (6,*) 'Calling CLOUDY on system '//num//' of '//ns
               call system(command)
            end if
         end if
      end do
      call chdir(home)
      end subroutine
