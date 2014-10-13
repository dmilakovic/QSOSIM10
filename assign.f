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
      count=0
      do while (count.lt.ni(1))
 1       c=ran3(idum)  
         dmin=1.0
         if((c.ge.CDDF(index12)).and.(c.lt.CDDF(index14)))then 
            count=count+1
            do j=1,npoints
               delta=abs(c-CDDF(j))
               if (delta.lt.dmin) then
                  dmin=delta
                  index=j
               end if
            end do
            lognhi=xs(index)
            nhi=10**xs(index)
            nhi4(count)=xs(index)
            d=ran3(idum)	
            gp1=gamma(1)+1.
            p=z2p1**gp1
            q=z1p1**gp1
            x=dble(alog10(real(d*(p-q)+q)))
            x=x/gp1
            z=10**x -1.0
            z4(count)=z
         else
            goto 1
         end if
      end do 
      do while (count.lt.(ni(1)+ni(2)))
 2       c=ran3(idum)
         dmin=1.0
         if((c.ge.CDDF(index14)).and.(c.lt.CDDF(index17)))then
            count=count+1
            do j=1,npoints
               delta=abs(c-CDDF(j))
               if (delta.lt.dmin) then
                  dmin=delta
                  index=j
               end if
            end do
            lognhi=xs(index)
            nhi=10**xs(index)
            nhi4(count)=xs(index)
            d=ran3(idum)	
            gp1=gamma(2)+1.
            p=z2p1**gp1
            q=z1p1**gp1
            x=dble(alog10(real(d*(p-q)+q)))
            x=x/gp1
            z=10**x -1.0
            z4(count)=z
         else
            goto 2
         end if
      end do 
      do while (count.lt.(ni(1)+ni(2)+ni(3)))
 3       c=ran3(idum)
         dmin=1.0
         if((c.ge.CDDF(index17)).and.(c.lt.CDDF(index22)))then 
            count=count+1
            do j=1,npoints
               delta=abs(c-CDDF(j))
               if (delta.lt.dmin) then
                  dmin=delta
                  index=j
               end if
            end do
            lognhi=xs(index)
            nhi=10**xs(index)
            nhi4(count)=xs(index)
            d=ran3(idum)	
            gp1=gamma(3)+1.
            p=z2p1**gp1
            q=z1p1**gp1
            x=dble(alog10(real(d*(p-q)+q)))
            x=x/gp1
            z=10**x -1.0
            z4(count)=z
         else
            goto 3
         end if
      end do 

      end subroutine assign
