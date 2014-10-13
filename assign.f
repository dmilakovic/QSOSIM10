      subroutine assign(npoints,zstart,zqso,xs,CDDF,gamma,
     +                  numlin,ni,nhi4,z4)

      integer npoints,numlin,idum,index,i,j,count
      integer index12,index14,index17,index22
      integer,dimension(3) :: ni
      real*8,dimension(3) :: gamma
      real*8,dimension(npoints) :: xs,CDDF
      real*8,dimension(numlin) :: nhi4,z4
      real*8 :: c,d,p,q,ran3,dmin
      real*8 :: zstart,zqso,z1p1,z2p1
      real*8 :: nhi,lognhi,z,delta,x
      external ran3

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
      index=0
      do i=1,numlin
         dmin=1.0
         delta=1.0
         nhi4(i)=0
         z4(i)=0
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
         nhi4(i)=xs(index)
! --------------------------------------------------------------------
         d=ran3(idum)
         if ((lognhi.ge.12.0).and.(lognhi.lt.14.0)) then
            g=gamma(1)
            gp1=g+1.
            p=z2p1**gp1
            q=z1p1**gp1
            x=alog10(real(d*(p-q)+q))
            x=x/gp1
            z=10**x -1.0
         else if ((lognhi.ge.14.0).and.(lognhi.lt.17.0)) then
            g=gamma(2)
            gp1=g+1.
            p=z2p1**gp1
            q=z1p1**gp1
            x=alog10(real(d*(p-q)+q))
            x=x/gp1
            z=10**x -1.0
         else if ((lognhi.ge.17.0).and.(lognhi.lt.22.0)) then
            g=gamma(3)
            gp1=g+1.
            p=z2p1**gp1
            q=z1p1**gp1
            x=alog10(real(d*(p-q)+q))
            x=x/gp1
            z=10**x -1.0
         end if
         z4(i)=z
      end do 

      end subroutine assign
