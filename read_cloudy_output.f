      subroutine read_cloudy_output(home,descriptor,atom,ion,s,mask,
     &              colden)
      
      character home*120, descriptor*15, atom*2, ion*4, num*4
      character folder*90, filename*120, line*15000
      integer i,j,s,nr,iostat
      integer, parameter :: maxrecs=100000
      real*8, dimension(s) :: colden
      logical :: ex
      
      write (6,*) 'Reaching the clouds'
      folder = trim(trim(home)//'/cloudy_in/spec_'//descriptor)
      call chdir(folder)
      call system('pwd')
      do i=1,s
         write(num,'(i4.4)') i
         filename = trim(folder)//'/system_'//num//'.rlt'
         filename = trim(filename)
         inquire(file=filename, exist=ex)
         write (6,*) ex
         if (ex.eqv..true.) then
            NR = 0
            write (6,*) filename
            OPEN(unit=15,file=filename,iostat=io)!,err=104)
 104        write (6,*) iostat
            DO J=1,maxrecs 
               READ(15,*,IOSTAT=ios) junk 
               IF (ios /= 0) EXIT
               IF (J == maxrecs) THEN
                  write(*,*) 'Error: Maximum number of records exceeded'
                  write(*,*) 'Exiting program now...'
                  STOP 
               ENDIF 
               NR = NR + 1
            ENDDO 
            write (6,*) 'NR',NR
            REWIND(15) 
!     Now read data into mydata 
            DO J=1,NR 
               READ(15,*) line
               WRITE(6,*) 'read line'
            ENDDO 
            CLOSE(1)
            colden(i)=1.0
            end if
         end do
         call chdir(home)
      end subroutine
