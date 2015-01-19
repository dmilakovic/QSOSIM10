!=======================================================================
      subroutine assign(npts,zstart,zqso,xs,CDDF,
     +                  nl,nhi4,z4,na)
!=======================================================================
*      Using the CDDF interpolated from the spline, for a given redshift
*      range, fill the arrays of H I column densities and redshifts to 
*      be used in qsosim
                    
*      input, real*8              : zstart       lowest redshift at which Lyman
*                                                alpha absorption can be seen
*                                 : zqso         QSO's redshift
*      input, integer             : npts      number of points of CDDF
*                                   nl           number of H I lines
*      input, real*8 array(nl)    : xs           values of log N(H I) 
*                                   CDDF         values of H I column density 
*                                                distribution function 

*      output, real*8 array(nl)   : nhi4         H I column densities
*                                   z4           H I redshifts
!=======================================================================

      integer npts,nl,idum,index,i,j,count
      integer index12,index14,index17,index22
      integer,dimension(3) :: ni,na
      real*8,dimension(3) :: gamma
      real*8,dimension(npts) :: xs,ys,CDDF,zx,CDF1,CDF2,CDF3
      real*8,dimension(nl) :: nhi4,z4
      real*8 :: c,d,p,q,ran3,dmin
      real*8 :: zstart,zqso,z1p1,z2p1
      real*8 :: nhi,lognhi,z,delta,x

*     temporary auxiliary variables
      integer n1,n2,n3,t
      real*4 :: ymin, ymax
      character :: a

      external ran3,clustering
      data gamma/1.51,2.16,1.33/
      n1=0;n2=0;n3=0
! -------------------------------------------------------------------
! CREATE ARRAYS CONTAINING DATA ON NHI AND Z         
! -------------------------------------------------------------------
*     Initialize parameters
      z1p1=zstart+1.
      z2p1=zqso+1.
      idum=time()
      index=0
*     introduce clustering
      call clustering(npts,zstart,zqso,zx,CDF1,CDF2,CDF3)
*     for each of the nl lines, 
      do i=1,nl
*        reset parameters
         dmin=1.0
         delta=1.0
         nhi4(i)=0
         z4(i)=0
*        generate a random number and find it's position in the array
 1       c=ran3(idum)
         if (c.lt.1.0) then
! -------------------------------------------------------------------
	    call locate(npts,xs,CDDF,c,lognhi)
         else
            goto 1
         end if
*        determine its log NHI by reading the array
         nhi4(i)=lognhi
! --------------------------------------------------------------------
*        and randomly determine its redshift, depending on the column density
*        use CDFs that correspond to the column density        
         d=ran3(idum)
         if ((lognhi.ge.12.00).and.(lognhi.lt.14.0)) then
            n1=n1+1
            call locate(npts,zx,CDF1,d,z)
         else if ((lognhi.ge.14.0).and.(lognhi.lt.17.0)) then
            n2=n2+1
            call locate(npts,zx,CDF2,d,z)
         else if ((lognhi.ge.17.0).and.(lognhi.lt.22.0)) then
            n3=n3+1
            call locate(npts,zx,CDF3,d,z)
         end if
         z4(i)=z
      end do 
      t=n1+n2+n3
      write (6,*) 'n1,n2,n3,t',n1,n2,n3,t
      na(1)=n1
      na(2)=n2
      na(3)=n3
      ymin=0.97*zstart; ymax=1.03*zqso

*     plot the cumulative distributions of lines in redshift
*     (dm: 2015/01/20 - behaving as expected, commenting out)
c      call pgbegin(0,'/xserve',1,1)
c      call pgslw(3)
c      call pgenv(real(zstart),real(zqso),0.0,1.0,0,0)
c      call pgpt(nl,real(nhi4),real(z4),-2)
c      call pgline(npts,real(zx),real(CDF1)); call pgsci(2)
c      call pgline(npts,real(zx),real(CDF2)); call pgsci(3)
c      call pgline(npts,real(zx),real(CDF3))
c      call pgend
c      read(*,'(a)')a
      return
      end subroutine assign


!=======================================================================
      subroutine clustering(npt,z1,z2,xx,cCDF1,cCDF2,cCDF3)
!=======================================================================
      implicit none
      integer :: i,j,k,t,idum,nc,clin,npt
      integer,parameter :: ncl=20,p=100
      real :: xi0,v0,A,gamma,gammas(3),c,R
      real*8 :: z,v,zmin,zmax,nc_total,c_total
      real*8 :: gp1,z1p1,z2p1,ymin,ymax
      real*8 :: n,cl,h,nu,integrate,qsimp,tf
      real*8 :: ran3,z1,z2,dz,sum
      real*8,dimension(npt) :: xx,cCDF1,cCDF2,cCDF3
      real*8,dimension(npt) :: xy,yy,ncCDF
      real*8,dimension(p) :: zx,ncy,cy
      real*8,dimension(1000):: zcl
      character :: string*10, string2*10

      common /cls/ zcl
      common /vars/ xi0,v0,A,gamma,c,R

      external n,cl,h,nu,integrate,ran3,qsimp,tf

*     set parameters 
      data c/3e5/
      data gammas/1.51,2.16,1.33/
      A=33.4
      xi0=0.5
      v0=250.0

      idum=time()
      zmin=0.9*z1
      zmax=1.1*z2
*     set number of clusters and redshift limits
      t=nint(integrate(cl,zmin,zmax))
      write(6,*) 'Number of clusters:',t
      do i=1,t
         zcl(i)=zmin+ran3(idum)*(zmax-zmin)
c         write(6,*) i,zcl(i)!,h(zcl(i))*nu(zcl(i))
      end do
!-----------------------------------------------------------------------    
*     construct probability density functions by integrating
*     we only want a subset of the whole range
*     set new redshift limits and mark them on the plot
*     calculate the number of clusters in that region
c      clin=0
c      do i=1,t
c         if (zcl(i).ge.z1.and.zcl(i).le.z2)clin=clin+1
c      end do
*     the total number of lines between z1 and z2:
*     nc = no clustering
*     c = clustering
      nc_total=integrate(n,z1,z2)
      c_total=integrate(nu,z1,z2)
      
c      write(6,*) 'No of lines (no clustering): ',nc_total
c      write(6,*) 'No of lines (clustering): ',c_total
*     construct CDFs by evaluating at npt points
      do i=1,npt
         xx(i)=z1+dble(i-1)*(z2-z1)/(npt-1)
      end do
      sum=0.0
      gamma=gammas(1)
      c_total=integrate(nu,z1,z2)      
      do i=1,npt
         z=xx(i)
         if (i.eq.1) then
            dz=0
         else 
            dz=xx(i)-xx(i-1)
         end if
         sum=sum+nu(z)*dz
         cCDF1(i)=sum/c_total
      end do
      sum=0.0
      gamma=gammas(2)
      c_total=integrate(nu,z1,z2)      
      do i=1,npt
         z=xx(i)
         if (i.eq.1) then
            dz=0
         else 
            dz=xx(i)-xx(i-1)
         end if
         sum=sum+nu(z)*dz
         cCDF2(i)=sum/c_total
      end do
      sum=0.0
      gamma=gammas(3)
      c_total=integrate(nu,z1,z2)      
      do i=1,npt
         z=xx(i)
         if (i.eq.1) then
            dz=0
         else 
            dz=xx(i)-xx(i-1)
         end if
         sum=sum+nu(z)*dz
         cCDF3(i)=sum/c_total
      end do
      do i=1,npt
         write(6,*) i,xx(i),cCDF1(i),cCDF2(i),cCDF3(i)
      end do
      
      end subroutine clustering

!=======================================================================
      function n(z)
!=======================================================================
      common /vars/ xi0,v0,A,gamma,c,R
      real*8 n,z
      real A,gamma
      n=A*(1+z)**gamma
      return
      end function n
!=======================================================================
      function cl(z)
!=======================================================================
      common /vars/ xi0,v0,A,gamma,c,R
      real*8 cl,z
      real xi0,v0,c
      cl=1./(4*xi0*(1.+z)*(v0/c))
      return
      end function cl
!=======================================================================
      function nu(z)
!=======================================================================
      integer i,k,n,idum
      real*8 nu,z,aux,z_cl(1000),ran3,tmp
      external ran3
      common /vars/ xi0,v0,A,gamma,c,R
      common /cls/ z_cl
      
      R=v0/c*(1+z)
      aux=0.0
c      write(6,*) R
      do i=1,1000
         tmp=z_cl(i)
         aux=aux+exp(-abs(z-tmp)/R)*4*xi0*A*(1.+z)**gamma
      end do
      nu=aux
      return
      end function nu

!=======================================================================
      function integrate(f,a,b)
!=======================================================================
      implicit none
      integer,parameter :: npt=1000
      real*8 :: integrate, sumup,sumlow,sum, dx
      real*8 ::  f, a, b
      integer i
      real*8,dimension(npt):: x, y
      do i = 1,npt
         x(i)=a+(b-a)*dble(i-1)/(npt-1)
         y(i)=f(x(i))
      end do
      sumup=0.0
      sumlow=0.0
      do i=2,npt
         dx=x(i)-x(i-1)
         sumup=sumup+f(x(i))*dx
c         sumlow=sumlow+f(x(i+1))*dx
      end do
c      sum=0.5*(sumup+sumlow)
      integrate = sumup
      return
      end function integrate
!=======================================================================
      subroutine locate(npt,xbuf,ybuf,y, x)
!=======================================================================
*     Purpose: for a given value of the CDF 'y', search the CDF array 
*     'ybuf' to find the element that is closest in value. Read the 
*     corresponding value 'x' of the pair (x,y) from 'xbuf'.

*     INPUT,real*8 :        xbuf(npt), ybuf(npt), y
*     OUTPUT,real*8 :       x
!=======================================================================
      integer i,j,npt,index
      real*8,dimension(npt) :: xbuf, ybuf
      real*8 :: x,y,d,dmin

      dmin=1.0
      do i=1,npt
         d=abs(y-ybuf(i))
         if (d.lt.dmin) then
            dmin=d
            index=i
         end if
      end do
      x=xbuf(index)
      return
      end subroutine locate
