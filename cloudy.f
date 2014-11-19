      subroutine cloudy(home,descriptor,nl,s)

      integer i,s,t,nl
      character :: descriptor*15,num*4,ns*4
      character :: path*120,home*120,folder*220,filename*150
      character command*200
      logical ex

      write (6,*) '--------------------------------------------------'
      write (6,*) '              The sky is CLOUDY                   '
      write (6,*) '--------------------------------------------------'
      path='/Users/dm/Programming/cloudy/c13.03/source/cloudy.exe'
      folder = trim(home)//'/cloudy_in/spec_'//descriptor//'/'
      write (ns,'(i4.4)') s

      call chdir(folder)
      call system('pwd')
      t=0
      do i=1,s
         t=t+1
         write(num,'(i4.4)') t
         filename = trim('system_'//num//'.in')
c      write (6,*) filename
         inquire(file=filename, exist=ex)
c     write (6,*) ex
         if (ex.eqv..true.) then
            filename = trim('system_'//num)
            command = trim(path)//' -r '//filename
c     write (6,*) command
            write (6,*) 'Calling CLOUDY on system '//num//' of '//ns
            call system(command)
         end if
         if (ex.eqv..false.) write (6,*) 'Cloudy file not found!!'
      
      end do
      call chdir(home)
      return
      end subroutine
