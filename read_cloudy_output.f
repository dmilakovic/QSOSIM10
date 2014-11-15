      subroutine read_cloudy_output(home,descriptor,atom,ion,nl,mask,
                   colden)
      
      character home*120, descriptor*8, atom*2, ion*4
      character folder*90, filename*90, line*120
      integer i,j,s
      real*8,dimension(nl) colden
      integer,dimension(nl) mask
      
      write (6,*) 'Touching clouds'
      folder = trim(trim(home)//'/cloudy_in/spec_'//descriptor)
      call chdir(folder)
      call system('pwd')
      s=0
      do i=1,nl
         if (mask(i).eq.1) then
            s=s+1
            colden(i)=1.0
         else
            colden(i)=0.0
         end if
      end do

      call chdir(home)
      end subroutine
