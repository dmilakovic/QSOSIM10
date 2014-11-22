      subroutine assign(npoints,zstart,zqso,xs,CDDF,
     +                  nl,ni,nhi4,z4)

      integer npoints,nl,idum,index,i,j,count
      integer s1,s2,s3
      integer index12,index14,index17,index22
      integer,dimension(3) :: ni
      real*8,dimension(3) :: gamma,bins
      real*8,dimension(npoints) :: xs,CDDF
      real*8,dimension(4000) :: nhi4,z4,gamma_int
      real*8 :: c,d,p,q,ran3,dmin
      real*8 :: zstart,zqso,z1p1,z2p1
      real*8 :: nhi,lognhi,z,delta,x
      real*8 :: g_i, g_e
      external ran3
      data bins/12.0,14.0,17.0/
      data gamma/1.51,2.16,1.33/
! -------------------------------------------------------------------
! CALCULATE INDICES for array manipulation
! -------------------------------------------------------------------     
      dmin=0.1
      do i=1,npoints
         d12=abs(xs(i)-12.0)
         if(d12.lt.dmin)then
            dmin=d12
            index12=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npoints
         d14=abs(xs(i)-14.0)
         if(d14.lt.dmin)then
            dmin=d14
            index14=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npoints
         d17=abs(xs(i)-17.0)
         if(d17.lt.dmin)then
            dmin=d17
            index17=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npoints
         d22=abs(xs(i)-22.0)
         if(d22.lt.dmin)then
            dmin=d22
            index22=i
         end if
      end do
! -------------------------------------------------------------------
! CREATE ARRAYS CONTAINING DATA ON NHI AND Z         
! -------------------------------------------------------------------
      z1p1=zstart+1.
      z2p1=zqso+1.
      idum=time()
      index=0; i=0
      s1=0; s2=0; s3=0
      write (6,*) s1,ni(1)
      write (6,*) s2,ni(2)
      write (6,*) s3,ni(3)
      do while (s1.ne.ni(1))
         dmin=1.0
         delta=1.0
 1       c=ran3(idum)
         if (c.lt.1.0) then
! -------------------------------------------------------------------
	    do j=1,npoints
               delta=abs(c-CDDF(j))
               if (delta.lt.dmin) then
                  dmin=delta
                  index=j
               end if
            end do
         else
            goto 1
         end if
         lognhi=xs(index)
! --------------------------------------------------------------------
         if ((lognhi.ge.12.0).and.(lognhi.lt.14.0)) then
            i=i+1
            s1=s1+1
            nhi4(i)=xs(index)
            d=ran3(idum)
            call polint(bins,gamma,3,nhi4(i),g_i,g_e)
            g=g_i
            gp1=g+1.
            p=z2p1**gp1
            q=z1p1**gp1
            x=alog10(real(d*(p-q)+q))
            x=x/gp1
            z=10**x -1.0
            z4(i)=z
            gamma_int(i)=g_i
         end if
      end do
      do while (s2.ne.ni(2))
         dmin=1.0
         delta=1.0
 2       c=ran3(idum)
         if (c.lt.1.0) then
! -------------------------------------------------------------------
	    do j=1,npoints
               delta=abs(c-CDDF(j))
               if (delta.lt.dmin) then
                  dmin=delta
                  index=j
               end if
            end do
         else
            goto 2
         end if
         lognhi=xs(index)
! --------------------------------------------------------------------
         if ((lognhi.ge.14.0).and.(lognhi.lt.17.0)) then
            i=i+1
            s2=s2+1
            nhi4(i)=xs(index)
            d=ran3(idum)
            call polint(bins,gamma,3,nhi4(i),g_i,g_e)
            g=g_i
            gp1=g+1.
            p=z2p1**gp1
            q=z1p1**gp1
            x=alog10(real(d*(p-q)+q))
            x=x/gp1
            z=10**x -1.0
            z4(i)=z
            gamma_int(i)=g_i
         end if
      end do
      do while (s3.ne.ni(3))
         dmin=1.0
         delta=1.0
 3       c=ran3(idum)
         if (c.lt.1.0) then
! -------------------------------------------------------------------
	    do j=1,npoints
               delta=abs(c-CDDF(j))
               if (delta.lt.dmin) then
                  dmin=delta
                  index=j
               end if
            end do
         else
            goto 3
         end if
         lognhi=xs(index)
! --------------------------------------------------------------------
         if ((lognhi.ge.17.0).and.(lognhi.lt.22.0)) then
            i=i+1
            s3=s3+1
            nhi4(i)=xs(index)
            d=ran3(idum)
            call polint(bins,gamma,3,nhi4(i),g_i,g_e)
            g=g_i
            gp1=g+1.
            p=z2p1**gp1
            q=z1p1**gp1
            x=alog10(real(d*(p-q)+q))
            x=x/gp1
            z=10**x -1.0
            z4(i)=z
            gamma_int(i)=g_i
         end if
      end do
c      do j=1,nl
c         write(6,*) j, nhi4(j), z4(j), gamma_int(j)
c      end do
c      write (6,*) s1,s2,s3,i
      return
      end subroutine assign


      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=2) 
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do i=1,n
         dift=abs(x-xa(i))
         dift=abs(x-xa(i))
         if (dift.lt.dif) then
            ns=i
            dif=dift
         endif
         c(i)=ya(i) 
         d(i)=ya(i)
      end do 
      y=ya(ns) 
      ns=ns-1
      do m=1,n-1 
         do i=1,n-m 
            ho=xa(i)-x
            hp=xa(i+m)-x
            w=c(i+1)-d(i)
            den=ho-hp
            if(den.eq.0.) stop 'failure in polint'
            den=w/den
            d(i)=hp*den
            c(i)=ho*den
         end do
         if (2*ns.lt.n-m)then
            dy=c(ns+1)
         else
            dy=d(ns)
            ns=ns-1
         endif
         y=y+dy
      end do
      return
      END SUBROUTINE polint

