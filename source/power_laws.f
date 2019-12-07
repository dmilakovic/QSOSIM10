!=======================================================================
      subroutine power_laws(npts,zstart,zqso,xs,ys,CDDF,
     +                      bigA,gamma,nl)
!=======================================================================
* PURPOSE: use data from Kim et al. 2013 & Unknown to predict the number
*          of H I absorption systems in three log N(HI) bins. Return the
*          number of lines 'nl' into the main program.

*      input, real*8              : zstart       lowest redshift at which Lyman
*                                                alpha absorption can be seen
*                                 : zqso         QSO's redshift
*      input, real*8(3)           : bigA         values of A from A*(1+z)**gamma
*                                   gamma        values of gamma from A*(1+z)**gamma
*      input, integer             : npoints      number of points of CDDF
*      input, real*8 array(nl)    : xs           values of log N(H I) 
*                                   ys           values of f(N_HI,X)
*                                   CDDF         values of H I column density 
*                                                distribution function 
*      output, integer            : nl           number of H I lines
*      output, integer(3)         : ni           number of H I lines in bins:
*                                                log N(HI) = [12.0,14.0,17.0]
!=======================================================================
      implicit none
      integer :: index12,index13p1,index14,index17,index22,index12p75
      integer :: i,j,index,nl,npts,idum,count
      real*8,dimension(npts) :: xs,ys,CDDF
      real*8,dimension(3) :: bigA, gamma, n
      integer,dimension(3) :: ni
      real*8 :: dmin,d12,d12p75,d13p1,d14,d17,d22
      real*8 :: int12to14,int13p1to14,int12p75to14
      real*8 :: n1,n2,n3,total,gasdev3
      real*8 :: z1,z1p1,z2,z2p1,zqso, gp1
      real*8 :: zstart,corr,beta, dx
  
      external gasdev3
! -------------------------------------------------------------------
*     CALCULATE INDICES for data manipulation
*     Find elements of 'xs' that have the closest possible value to:
*     12.0, 12.75, 13.1, 14.0, 17.0, 22.0
! ------------------------------------------------------------------- 
      dmin=0.1
      do i=1,npts
         d12=abs(xs(i)-12.0)
         if(d12.lt.dmin)then
            dmin=d12
            index12=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npts
         d12p75=abs(xs(i)-12.75)
         if(d12p75.lt.dmin)then
            dmin=d12p75
            index12p75=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npts
         d13p1=abs(xs(i)-13.1)
         if(d13p1.lt.dmin)then
            dmin=d13p1
            index13p1=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npts
         d14=abs(xs(i)-14.0)
         if(d14.lt.dmin)then
            dmin=d14
            index14=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npts
         d17=abs(xs(i)-17.0)
         if(d17.lt.dmin)then
            dmin=d17
            index17=i
         end if
      end do
! -------------------------------------------------------------------
      dmin=0.1
      do i=1,npts
         d22=abs(xs(i)-22.0)
         if(d22.lt.dmin)then
            dmin=d22
            index22=i
         end if
      end do
! -------------------------------------------------------------------
! INTEGRATION to account for a gap in data from 12 to 13.1
! -------------------------------------------------------------------
*     Numerically integrate spline from 12 to 14 and from 13.1 to 14
*     in steps of dx to get the correction factor 'corr'. That way
*     we correct the value of A to include systems down to 1e12 cm-2.
      int12to14=0
      int12p75to14=0
      do i=index12,index14
         dx=10**xs(i+1)-10**xs(i)             
         int12to14=int12to14+10**ys(i)*dx
      end do
      do i=index12p75,index14
         dx=10**xs(i+1)-10**xs(i)             
         int12p75to14=int12p75to14+10**ys(i)*dx
      end do
      int13p1to14=0
      do i=index13p1,index14
         dx=10**xs(i+1)-10**xs(i)             
         int13p1to14=int13p1to14+10**ys(i)*dx
      end do
      corr=int12to14/int13p1to14

c      write (6,*) 'int[12.75->14] = ',int12p75to14
c      write (6,*) 'int[13.10->14] = ',int13p1to14
C     corr - correction factor

! -------------------------------------------------------------------
! DETERMINE THE NUMBER OF LINES IN EACH COLUMN DENSITY BIN         
! -------------------------------------------------------------------
      ! dn/dz = A(1+z)^gamma
      ! n = A/(1+gamma)*[(1+z2)^(1+gamma) - (1+z1)^(1+gamma)]
      ! n1(12->14) = corr*A1/(1+gamma1)*[(1+z2)^(1+gamma1) - (1+z1)^(1+gamma1)]
      ! n2(14->17) = A2/(1+gamma2)*[(1+z2)^(1+gamma2) - (1+z1)^(1+gamma2)]
      ! n3(17->22) = A3/(1+gamma3)*[(1+z2)^(1+gamma3) - (1+z1)^(1+gamma3)]
! -------------------------------------------------------------------
*     Initialize parameters
      z1=zstart
      z1p1=z1+1.
      z2=zqso
      z2p1=zqso+1.
c      beta=abs((ys(index14)-ys(index12))/(xs(index14)-xs(index12)))
c      write (6,*) 'Numerical correction factor=',corr
! -------------------------------------------------------------------
*     Predict the number of lines in three bins: 
*     n(bin) = integral( A(bin) * (1+z)*gamma(bin) ) in limits bin_down, bin_up
      do i=1,3
         gp1=gamma(i)+1.
         if(i.eq.1)then
            n(i) = corr*bigA(i)/(gp1)*((z2p1)**(gp1)-(z1p1)**(gp1))
         else
         n(i) = bigA(i)/(gp1)*((z2p1)**(gp1)-(z1p1)**(gp1))
         end if
         total=total+n(i)
      end do   
      ni(1)=nint(n(1))
      ni(2)=nint(n(2))
      ni(3)=nint(n(3))
      nl=ni(1)+ni(2)+ni(3)
c      write (6,*) 'A = corr*bigA/(1+gamma)'
c      write (6,*) 'A (numerical) =',corr*bigA(1)/(gamma(1)+1.)
c      write (6,*) 'gamma = ',gamma(1)
c      write (6,*) 'z1 =',z1
c      write (6,*) 'z2 =',z2
c      write (6,*) 'n [12.00->14.00] (predicted) =',nint(n(1))
c      write (6,*) 'n [14.00->17.00] (predicted) =',nint(n(2))
c      write (6,*) 'n [17.00->22.00] (predicted) =',nint(n(3))
c      write (6,*) 'Total number of lines (predicted) =', nl 
! -------------------------------------------------------------------
*     The predicted number of lines is normally distributed variable with
*     a mean=nl and sigma=sqrt(nl).
*     Therefore, we randomly sample a gaussian to determine the number 
*     of lines that will be returned into the main program.
      idum=time()
      nl=nl+sqrt(real(nl))*gasdev3(idum)
      write (6,*) 'Total number of H I lines =', nl 
      return
      end subroutine power_laws

!=======================================================================
      SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2)
      INTEGER mwt,ndata
      REAL*4 a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
C     USES gammq
C     Given a set of data points x(1:ndata),y(1:ndata) with individual standard deviations
C     sig(1:ndata),fit them to a straight line y = a + bx by minimizing chi2. Returned are
C     a,b and their respective probable uncertainties siga and sigb, the chi-square chi2, and
C     the goodness-of-fit probability q (that the fit would have chi2 this large or larger). If mwt=0
C     on input, then the standard deviations are assumed to be unavailable: q is returned as 1.0
C     and the normalization of chi2 is to unit standard deviation on all points.
      INTEGER i
      REAL*4 sigdat,ss,st2,sx,sxoss,sy,t,wt,gammq
      sx=0.
      sy=0.
      st2=0.
      b=0.
      if(mwt.ne.0) then
         ss=0.
         do i=1,ndata
            wt=1./(sig(i)**2)
            ss=ss+wt
            sx=sx+x(i)*wt
            sy=sy+y(i)*wt
         end do
      else
         do i=1,ndata 
            sx=sx+x(i)
            sy=sy+y(i)
         end do
         ss=float(ndata)
      endif
      sxoss=sx/ss
      if(mwt.ne.0) then
         do i=1,ndata
            t=(x(i)-sxoss)/sig(i)
            st2=st2+t*t
            b=b+t*y(i)/sig(i)
         end do
      else
         do i=1,ndata
            t=x(i)-sxoss
            st2=st2+t*t
            b=b+t*y(i)
         end do
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1./st2)
      chi2=0.
      q=1.
      if(mwt.eq.0) then
         do i=1,ndata
            chi2=chi2+(y(i)-a-b*x(i))**2
         end do
         sigdat=sqrt(chi2/(ndata-2)) 
         siga=siga*sigdat
         sigb=sigb*sigdat
      else
         do i=1,ndata
            chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
         end do
c         if(ndata.gt.2) q=gammq(0.5*(ndata-2),0.5*chi2)
      endif
      write(6,*) a,b,chi2
      return
      END SUBROUTINE fit
