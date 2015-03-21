      subroutine metal_db(ind)

      implicit none
      integer,parameter :: d=20000
      integer :: i,j,l,m,lastchpos,ind,dim
      character*64 filnm
      character*200 cloudydb


      real*8,dimension(20000) :: met
      real*8,dimension(20000) :: rdf
      real*8,dimension(20000) :: lNHI
      real*8,dimension(30,20000) :: H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,
     :                      P,S,Cl,Ar,K,Ca,Sc,Ti,Va,Cr,Mn,Fe,Co,Ni,Cu,Zn
      common/metals/H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,
     :                      P,S,Cl,Ar,K,Ca,Sc,Ti,Va,Cr,Mn,Fe,Co,Ni,Cu,Zn
      common/parspace/met,lNHI,rdf
      external lastchpos

      dim=d
      if(ind.gt.0) goto 990
      cloudydb=' '
      call getenv ( 'CLOUDYDB', cloudydb )
      if ( cloudydb .eq. ' ' ) then
           cloudydb = './metals.csv'
      end if
      open(unit=18,file=cloudydb,status='old',err=99)
      goto 11
 99   i=lastchpos(cloudydb)
      if(cloudydb(i:i).eq.'/') then
          cloudydb=cloudydb(1:i)//'metals.csv'
         else
          cloudydb=cloudydb(1:i)//'/metals.csv'
      end if
      open(unit=18,file=cloudydb,status='old',err=990)
      goto 11
 990  write(6,*) 'CLOUDY database not found'
      stop
 11   continue
      inquire( 18, name=cloudydb )
      write ( *, * ) ' Using data from : '//cloudydb
      do i=1,dim
          read(18,*,end=97)met(i),lNHI(i),rdf(i),H(1,i),
     :He(1,i),He(2,i),Li(1,i),Li(2,i),Li(3,i),Be(1,i),
     :Be(2,i),Be(3,i),Be(4,i),B(1,i),B(2,i),B(3,i),B(4,i),
     :B(5,i),C(1,i),C(2,i),C(3,i),C(4,i),C(5,i),C(6,i),
     :N(1,i),N(2,i),N(3,i),N(4,i),N(5,i),N(6,i),N(7,i),O(1,i),
     :O(2,i),O(3,i),O(4,i),O(5,i),O(6,i),O(7,i),O(8,i),F(1,i),
     :F(2,i),F(3,i),F(4,i),F(5,i),F(6,i),F(7,i),F(8,i),F(9,i),
     :Ne(1,i),Ne(2,i),Ne(3,i),Ne(4,i),Ne(5,i),Ne(6,i),Ne(7,i),Ne(8,i),
     :Ne(9,i),Ne(10,i),Na(1,i),Na(2,i),Na(3,i),Na(4,i),Na(5,i),
     :Na(6,i),Na(7,i),Na(8,i),Na(9,i),Na(10,i),Na(11,i),
     :Mg(1,i),Mg(2,i),Mg(3,i),Mg(4,i),Mg(5,i),Mg(6,i),Mg(7,i),Mg(8,i),
     :Mg(9,i),Mg(10,i),Mg(11,i),Mg(12,i),Al(1,i),Al(2,i),
     :Al(3,i),Al(4,i),Al(5,i),Al(6,i),Al(7,i),Al(8,i),Al(9,i),Al(10,i),
     :Al(11,i),Al(12,i),Al(13,i),Si(1,i),Si(2,i),Si(3,i),
     :Si(4,i),Si(5,i),Si(6,i),Si(7,i),Si(8,i),Si(9,i),Si(10,i),Si(11,i),
     :Si(12,i),Si(13,i),Si(14,i),P(1,i),P(2,i),P(3,i),P(4,i),
     :P(5,i),P(6,i),P(7,i),P(8,i),P(9,i),P(10,i),P(11,i),P(12,i),
     :P(13,i),P(14,i),P(15,i),S(1,i),S(2,i),S(3,i),S(4,i),
     :S(5,i),S(6,i),S(7,i),S(8,i),S(9,i),S(10,i),S(11,i),S(12,i),
     :S(13,i),S(14,i),S(15,i),S(16,i),Cl(1,i),Cl(2,i),Cl(3,i),
     :Cl(4,i),Cl(5,i),Cl(6,i),Cl(7,i),Cl(8,i),Cl(9,i),Cl(10,i),Cl(11,i),
     :Cl(12,i),Cl(13,i),Cl(14,i),Cl(15,i),Cl(16,i),Cl(17,i),Ar(1,i),
     :Ar(2,i),Ar(3,i),Ar(4,i),Ar(5,i),Ar(6,i),
     :Ar(7,i),Ar(8,i),Ar(9,i),Ar(10,i),Ar(11,i),Ar(12,i),Ar(13,i),
     :Ar(14,i),Ar(15,i),Ar(16,i),Ar(17,i),Ar(18,i),K(1,i),
     :K(2,i),K(3,i),K(4,i),K(5,i),K(6,i),K(7,i),K(8,i),K(9,i),K(10,i),
     :K(11,i),K(12,i),K(13,i),K(14,i),K(15,i),K(16,i),K(17,i),K(18,i),
     :K(19,i),Ca(1,i),Ca(2,i),Ca(3,i),Ca(4,i),Ca(5,i),Ca(6,i),
     :Ca(7,i),Ca(8,i),Ca(9,i),Ca(10,i),Ca(11,i),Ca(12,i),Ca(13,i),
     :Ca(14,i),Ca(15,i),Ca(16,i),Ca(17,i),Ca(18,i),Ca(19,i),Ca(20,i),
     :Sc(1,i),Sc(2,i),Sc(3,i),Sc(4,i),Sc(5,i),Sc(6,i),Sc(7,i),
     :Sc(8,i),Sc(9,i),Sc(10,i),Sc(11,i),Sc(12,i),Sc(13,i),Sc(14,i),
     :Sc(15,i),Sc(16,i),Sc(17,i),Sc(18,i),Sc(19,i),Sc(20,i),Sc(21,i),
     :Ti(1,i),Ti(2,i),Ti(3,i),Ti(4,i),Ti(5,i),Ti(6,i),
     :Ti(7,i),Ti(8,i),Ti(9,i),Ti(10,i),Ti(11,i),Ti(12,i),Ti(13,i),
     :Ti(14,i),Ti(15,i),Ti(16,i),Ti(17,i),Ti(18,i),Ti(19,i),Ti(20,i),
     :Ti(21,i),Ti(22,i),Va(1,i),Va(2,i),Va(3,i),Va(4,i),
     :Va(5,i),Va(6,i),Va(7,i),Va(8,i),Va(9,i),Va(10,i),Va(11,i),
     :Va(12,i),Va(13,i),Va(14,i),Va(15,i),Va(16,i),Va(17,i),Va(18,i),
     :Va(19,i),VA(20,i),Va(21,i),Va(22,i),Va(23,i),Cr(1,i),
     :Cr(2,i),Cr(3,i),Cr(4,i),Cr(5,i),Cr(6,i),Cr(7,i),Cr(8,i),Cr(9,i),
     :Cr(10,i),Cr(11,i),Cr(12,i),Cr(13,i),Cr(14,i),Cr(15,i),Cr(16,i),
     :Cr(17,i),Cr(18,i),Cr(19,i),Cr(20,i),Cr(21,i),Cr(22,i),Cr(23,i),
     :Cr(24,i),Mn(1,i),Mn(2,i),Mn(3,i),Mn(4,i),Mn(5,i),Mn(6,i),
     :Mn(7,i),Mn(8,i),Mn(9,i),Mn(10,i),Mn(11,i),Mn(12,i),Mn(13,i),
     :Mn(14,i),Mn(15,i),Mn(16,i),Mn(17,i),Mn(18,i),Mn(19,i),Mn(20,i),
     :Mn(21,i),Mn(22,i),Mn(23,i),Mn(24,i),Mn(25,i),Fe(1,i),
     :Fe(2,i),Fe(3,i),Fe(4,i),Fe(5,i),Fe(6,i),Fe(7,i),Fe(8,i),Fe(9,i),
     :Fe(10,i),Fe(11,i),Fe(12,i),Fe(13,i),Fe(14,i),Fe(15,i),Fe(16,i),
     :Fe(17,i),Fe(18,i),Fe(19,i),Fe(20,i),Fe(21,i),Fe(22,i),Fe(23,i),
     :Fe(24,i),Fe(25,i),Fe(26,i),Co(1,i),Co(2,i),Co(3,i),
     :Co(4,i),Co(5,i),Co(6,i),Co(7,i),Co(8,i),Co(9,i),Co(10,i),Co(11,i),
     :Co(12,i),Co(13,i),Co(14,i),Co(15,i),Co(16,i),Co(17,i),Co(18,i),
     :Co(19,i),Co(20,i),Co(21,i),Co(22,i),Co(23,i),Co(24,i),Co(25,i),
     :Co(26,i),Co(27,i),Ni(1,i),Ni(2,i),Ni(3,i),Ni(4,i),
     :Ni(5,i),Ni(6,i),Ni(7,i),Ni(8,i),Ni(9,i),Ni(10,i),Ni(11,i),
     :Ni(12,i),Ni(13,i),Ni(14,i),Ni(15,i),Ni(16,i),Ni(17,i),Ni(18,i),
     :Ni(19,i),Ni(20,i),Ni(21,i),Ni(22,i),Ni(23,i),Ni(24,i),Ni(25,i),
     :Ni(26,i),Ni(27,i),Ni(28,i),Cu(1,i),Cu(2,i),Cu(3,i),
     :Cu(4,i),Cu(5,i),Cu(6,i),Cu(7,i),Cu(8,i),Cu(9,i),Cu(10,i),Cu(11,i),
     :Cu(12,i),Cu(13,i),Cu(14,i),Cu(15,i),Cu(16,i),Cu(17,i),Cu(18,i),
     :Cu(19,i),Cu(20,i),Cu(21,i),Cu(22,i),Cu(23,i),Cu(24,i),Cu(25,i),
     :Cu(26,i),Cu(27,i),Cu(28,i),Cu(29,i),Zn(1,i),Zn(2,i),
     :Zn(3,i),Zn(4,i)
*     DM (19/2/15):
*     works till here, but the rest isn't important anyway. keeping for bookkeeping purposes
c     :Zn(5,i),Zn(6,i),Zn(7,i),Zn(8,i),Zn(9,i),Zn(10,i),
c     :Zn(11,i),Zn(12,i),Zn(13,i),Zn(14,i),Zn(15,i),Zn(16,i),Zn(17,i),
c     :Zn(18,i),Zn(19,i),Zn(20,i),Zn(21,i),Zn(22,i),Zn(23,i),Zn(24,i),
c     :Zn(25,i),Zn(26,i),Zn(27,i),Zn(28,i),Zn(29,i),Zn(30,i)
      end do
 97   continue
c      write(6,*) 'METALS FOR QSOSIM10:'
c      do i=2000,dim
c         write(6,*)'---------------------------------------------------'
c         write(6,*) i,met(i),lNHI(i),rdf(i)
c         write(6,1) 'HI:',H(1,i)
c       write(6,1)'CI:',C(1,i),'CII:',C(2,i),'CIII:',C(3,i),'CIV:',C(4,i)
c         write(6,1) 'SiII:',Si(2,i),'SiIII:',Si(3,i),'SiIV:',Si(4,i)
c         write(6,1) 'MgI:',Mg(1,i),'MgII:',Mg(2,i)
c         write(6,1) 'FeI:',Fe(1,i),'FeII:',Fe(2,i),'FeIII:',Fe(3,i)
c         write(6,1) 'OI:',O(1,i),'OVI:',O(6,i)
c         write(6,1) 'AlII:',Al(2,i),'AlIII:',Al(3,i)
c         write(6,1) 'NiI:',Ni(1,i),'NiII:',Ni(2,i)
c       write(6,1) 'NI:',N(1,i),'NII:',N(2,i),'NIII:',N(3,i),'NV:',N(5,i)
c         write(6,1) 'MnII:',Mn(2,i),'ZnII:',Zn(2,i),'CrII:',Cr(2,i)
c      end do
 1    format (5(a6,e12.5,2x))
      rewind 18
      close(unit=18)
      return
      end subroutine metal_db
