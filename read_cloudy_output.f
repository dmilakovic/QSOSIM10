      subroutine read_cloudy_output(home,descriptor,atom,ion,nl,s,mask,
     &              colden)
      
      character home*120, descriptor*15, atom*2, ion*4, num*4
      character folder*90, name*40,filename*120, srow*5
      character junk*25,text4*9
      integer i,j,s,nr,io,a,b,t,nl
      integer,dimension(495) :: atom_num,ion_num
      integer, parameter :: max=10000
      real*8, dimension(s) :: coldens
      real*8, dimension(495) :: cd
      logical,dimension(4000) :: mask
      logical :: ex
      
      write (6,*) 'Reaching the clouds'
      folder = trim(trim(home)//'/cloudy_in/spec_'//descriptor)
      call chdir(folder)
      call system('pwd')
      io=0
    
      t=0
      do i=1,nl
         write (6,*) 'mask =  ',mask(i)
         if (mask(i).eqv..true.) then
            t=t+1
            coldens(t)=0
            write(num,'(i4.4)') t
            name = 'system_'//num//'.rlt'
            filename = trim(folder)//'/system_'//num//'.rlt'
            filename = trim(filename)
            write (6,*) 'filename     ',filename
            inquire(file=filename, exist=ex)
            if (ex.eqv..true.) then
               NR = 0
               OPEN(unit=25,file=name,iostat=io)
               write(6,*) 'File opened'
               DO J=1,max 
                  READ(25,*,iostat=io) junk 
                  IF (io.ne.0) EXIT
                  IF (J == max) THEN
                     write(*,*) 'Error: Maximum exceeded'
                     write(*,*) 'Exiting program now...'
                     STOP 
                  ENDIF 
                  NR = NR + 1
               END DO 
               REWIND(25) 
!     Now read data into mydata 
               DO J=1,15
                  READ(25,*) junk
               END DO
               DO J=1,495
                  READ(25,'(2i2,f16.10)') atom_num(j),ion_num(j),cd(j)
c                  write (6,55) atom_num(j),ion_num(j),cd(j)
               end do
               CLOSE(25)
      write (6,*) 'closed file'
               call find_atom(atom,ion,a,b)
               do j=1,495
                  if (atom_num(j).eq.a.and.ion_num(j).eq.b) then
                     coldens(t)=cd(j)
                    write(6,*) 'atom_num,ion_num,coldens',
     &                   atom_num(j),ion_num(j),coldens(i)
                  end if
               end do
            end if
         else
            write (6,*) i,'not dense enough'
         end if
      end do
 55      format (i2,2x,i2,2x,d16.10)
         call chdir(home)
         end subroutine
!-----------------------------------------------------------------------
      subroutine find_atom(atom,ion,a,b)
      character atom*2,ion*4
      integer a,b
      if (atom.eq.'H ') then 
         a=1
      else if (atom.eq.'He') then 
         a=2
      else if (atom.eq.'Li') then 
         a=3
      else if (atom.eq.'Be') then 
         a=4
      else if (atom.eq.'B ') then 
         a=5
      else if (atom.eq.'C ') then 
         a=6
      else if (atom.eq.'N ') then 
         a=7
      else if (atom.eq.'O ') then 
         a=8
      else if (atom.eq.'F ') then 
         a=9
      else if (atom.eq.'Ne') then 
         a=10
      else if (atom.eq.'Na') then 
         a=11
      else if (atom.eq.'Mg') then 
         a=12
      else if (atom.eq.'Al') then 
         a=13
      else if (atom.eq.'Si') then 
         a=14
      else if (atom.eq.'P ') then 
         a=15
      else if (atom.eq.'S ') then 
         a=16
      else if (atom.eq.'Cl') then 
         a=17
      else if (atom.eq.'Ar') then 
         a=18
      else if (atom.eq.'K ') then 
         a=19
      else if (atom.eq.'Ca') then 
         a=20
      else if (atom.eq.'Sc') then 
         a=21
      else if (atom.eq.'Ti') then 
         a=22
      else if (atom.eq.'Va') then 
         a=23
      else if (atom.eq.'Cr') then 
         a=24
      else if (atom.eq.'Ma') then 
         a=25
      else if (atom.eq.'Fe') then 
         a=26
      else if (atom.eq.'Co') then 
         a=27
      else if (atom.eq.'Ni') then 
         a=28
      else if (atom.eq.'Cu') then 
         a=29
      else if (atom.eq.'Zn') then 
         a=30
      end if
      if (ion.eq.'I   ') then
         b=1
      else if(ion.eq.'II  ') then
         b=2
      else if(ion.eq.'III ') then
         b=3
      else if(ion.eq.'IV  ') then
         b=4
      else if(ion.eq.'V   ') then
         b=5
      else if(ion.eq.'VI  ') then
         b=6
      else if(ion.eq.'VII ') then
         b=7
      else if(ion.eq.'VIII') then
         b=8
      else if(ion.eq.'IX  ') then
         b=9
      else if(ion.eq.'X   ') then
         b=10
      else if(ion.eq.'XI  ') then
         b=11
      else if(ion.eq.'XII ') then
         b=12
      else if(ion.eq.'XIII') then
         b=13
      else if(ion.eq.'XIV ') then
         b=14
      else if(ion.eq.'XV  ') then
         b=15
      else if(ion.eq.'XVI ') then
         b=16
      else if(ion.eq.'XVII') then
         b=17
      end if
         
      end subroutine find_atom
