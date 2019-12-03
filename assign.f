!=======================================================================
      subroutine assign(npts,zstart,zqso,rmag,alpha,xs,CDDF,
     +                  nl,nhi4,npnhi4,z4,znc4)
!=======================================================================
*      Using the CDDF interpolated from the spline, for a given redshift
*      range, fill the arrays of H I column densities and redshifts to 
*      be used in qsosim
                    
*      input, real*8              : zstart       lowest redshift at which Lyman
*                                                alpha absorption can be seen
*                                 : zqso         QSO's redshift
*      input, integer             : npts         number of points of CDDF
*                                   nl           number of H I lines
*      input, real*8 array(nl)    : xs           values of log N(H I) 
*                                   CDDF         values of H I column density 
*                                                distribution function 

*      output, real*8 array(nl)   : nhi4         H I column densities
*                                   z4           H I redshifts
*                                   npz4         H I redshifts no proximity effect
!=======================================================================

      integer npts,nl,idum,index,i,j,count
      integer index12,index14,index17,index22
      integer,dimension(3) :: ni,numl
      real*8,dimension(3) :: gamma
      real*8,dimension(npts) :: xs,ys,CDDF,zx,CDF1,CDF2,CDF3
      real*8,dimension(npts) :: ncCDF1,ncCDF2,ncCDF3,npncCDF
      real*4,dimension(npts) :: plotx,ploty
      real,dimension(3*npts) :: omegas,dndxy
      real*8,dimension(nl) :: nhi4,z4,npnhi4,znc4
      real*8 :: c,d,p,q,ran3,dmin,omega,integrate,dndx
      real*8 :: zstart,zqso,z1p1,z2p1,rmag,alpha,pef
      real*8 :: nhi,lognhi,z,npz,delta,x,omegaz,znc

*     temporary auxiliary variables
      integer n1,n2,n3,t
      real*4 :: ymin, ymax
      character :: a, lbl*6,name*22

c      common /qso/ zq,rmag,alpha

      external ran3,clustering,omega,integrate,dndx
      data gamma/1.51,2.16,1.94/
      n1=0;n2=0;n3=0
! -------------------------------------------------------------------
! CREATE ARRAYS CONTAINING DATA ON NHI AND Z         
! -------------------------------------------------------------------
*     Initialize parameters
      z1p1=zstart+1.
      z2p1=zqso+1.
      idum=time()
      index=0
*     (dm: 2015/3/2 there's a bug that requires omega to be called before using it)
      pef=(1.+omega(dble(3.0)))**(-0.6)
*     introduce clustering, generate CDFs with clustering
      call clustering(npts,zx,zstart,zqso,rmag,alpha,CDF1,CDF2,CDF3)
*     and generate CDFs with no clustering
      call no_clustering(npts,zx,zstart,zqso,ncCDF1,ncCDF2,ncCDF3)
*     for each of the nl lines, 
      do i=1,nl
*        reset/initialise parameters
         dmin=1.0
         delta=1.0
         nhi4(i)=0
         z4(i)=0
         npnhi4(i)=0
*        generate a random number and find it's position in the array
 1       c=ran3(idum)
         if (c.lt.1.0) then
! -------------------------------------------------------------------
	    call locate(npts,xs,CDDF,c,lognhi)
         else
            goto 1
         end if
*        determine the log NHI by reading the array
         nhi=10**lognhi
! --------------------------------------------------------------------
*        randomly determine the line's redshift
*        depending on the column density use CDFs that corresponds to 
*        the column density bin      
         d=ran3(idum)
         if ((lognhi.ge.12.00).and.(lognhi.lt.14.0)) then
            n1=n1+1
            call locate(npts,zx,CDF1,d,z)
            call locate(npts,zx,ncCDF1,d,znc)
         else if ((lognhi.ge.14.0).and.(lognhi.lt.17.0)) then
            n2=n2+1
            call locate(npts,zx,CDF2,d,z)
            call locate(npts,zx,ncCDF2,d,znc)
         else if ((lognhi.ge.17.0).and.(lognhi.le.22.0)) then
            n3=n3+1
            call locate(npts,zx,CDF3,d,z)
            call locate(npts,zx,ncCDF3,d,znc)
         end if
         z4(i)=z
         znc4(i)=znc
         npnhi4(i)=nhi
         pef=(1.+omega(z))**(-1.0)
         nhi4(i)=nhi*pef
      end do 
      t=n1+n2+n3
c      write (6,*) 'n1,n2,n3,t',n1,n2,n3,t
      numl(1)=n1
      numl(2)=n2
      numl(3)=n3
      ymin=0.97*zstart; ymax=1.03*zqso

*     plot the cumulative distributions of lines in redshift
*     (dm: 2015/01/20 - behaving as expected, commenting out)
*     (dm: 2015/01/31 - included the no-proximity effect CDF (npCDF), for testing and comparison)

*     plot clustering of lines (presentation purposes)
c      write(lbl,'(i6)') nint(zqso*1e5)
c      write(6,*) lbl
c      name='./plots/pres_CDF_clus.ps'
c      call pgbegin(0,'/null',1,1)
c      call pgbegin(0,name//'/cps',1,1)
c      call pgslw(3)
c      call pgenv(real(zstart),real(zqso),0.0,1.15,0,0)
c      call pgenv(real(zstart),1.1*real(zstart),0.0,0.1,0,0)
c      call pglabel('\fiz','CDF','')
c      call pgpt(nl,real(nhi4),real(z4),-2)
c      do i=1,npts
c         gp1=gamma(1)+1.
c         p=(zstart+1.)**gp1
c         q=(zqso+1.)**gp1
c         npncCDF(i)=((zx(i)+1)**gp1-p)/(q-p)
c      end do
c      call pgslw(5); call pgsls(3)
c      call pgline(npts,real(zx),real(npncCDF))
c      call pgslw(3); call pgsls(1)
c      call pgline(npts,real(zx),real(CDF1)) 
c      call pgmove(3.8,0.26); call pgdraw(3.95,0.26); call pgsci(2)
c      call pgline(npts,real(zx),real(CDF2))
c      call pgmove(3.80,0.17); call pgdraw(3.95,0.17); call pgsci(3)
c      call pgline(npts,real(zx),real(CDF3));
c      call pgmove(3.80,0.08); call pgdraw(3.95,0.08); call pgsci(2)
c      call pgslw(1)
c      do i=1,nl
c         call PGMOVE(real(z4(i)),1.13)
c         call PGDRAW(real(z4(i)),1.08)
c      end do
c      call pgsci(4)
c      do i=1,nl
c         call PGMOVE(real(znc4(i)),1.06)
c         call PGDRAW(real(znc4(i)),1.01)
c      end do
c      call pgslw(3)
c      call pgsci(1); call pgsls(4)
c      call pgline(npts,real(zx),real(ncCDF1)); call pgsci(2)
c      call pgline(npts,real(zx),real(ncCDF2)); call pgsci(3)
c      call pgline(npts,real(zx),real(ncCDF3)); call pgsci(1)
c      
c      call pgslw(3); call pgsls(1)
c      call pgmtxt('b',-7.,.8,0.,'10\u12\d\(2243)N\dHI\u<10\u14\d')
c      call pgmtxt('b',-4.5,.8,0.,'10\u14\d\(2243)N\dHI\u<10\u17\d')
c      call pgmtxt('b',-2.,.8,0.,'10\u17\d\(2243)N\dHI\u\(2243)10\u22\d')
c      call pgmtxt('b',-12.4,0.6,0.0,'Clustering only')
c      call pgend
c
c      do i=1,3*npts
c         z=zstart+dble(i-1)*(zqso-zstart)/(3*npts-1)
c         o=omega(z)
c         omegas(i)=alog10(o)
c         dndxy(i)=(1.+o)**(-1.5)
c         write(6,*) i,z,omegas(i),dndxy(i)
c      end do
c      call pgenv(minval(omegas),maxval(omegas),0.0,2.0,0,0)
c      call pgsci(2)
c      call pgline(npts,omegas,dndxy)
c      call pgend
c      read(*,'(a)')a

*     plot the cumulative distribution function (presentation purposes)
c      name='./plots/pres_CDF_lines.ps'
c      call pgbegin(0,'/null',1,1)
c      call pgbegin(0,name//'/cps',1,1)
c      call pgslw(3)
c      call pgenv(real(zstart),real(zqso),0.0,1.15,0,0)
c      call pgenv(real(zstart),1.1*real(zstart),0.0,0.1,0,0)
c      call pglabel('\fiz','CDF','')
c      call pgpt(nl,real(nhi4),real(z4),-2)
c      do i=1,npts
c         gp1=gamma(1)+1.
c         p=(zstart+1.)**gp1
c         q=(zqso+1.)**gp1
c         npncCDF(i)=((zx(i)+1)**gp1-p)/(q-p)
c      end do
c      call pgslw(5); call pgsls(3)
c      call pgline(npts,real(zx),real(npncCDF))
c      call pgslw(3); call pgsls(1)
c      call pgline(npts,real(zx),real(CDF1)) 
c      call pgmove(3.8,0.26); call pgdraw(3.95,0.26); call pgsci(2)
c      call pgline(npts,real(zx),real(CDF2))
c      call pgmove(3.80,0.17); call pgdraw(3.95,0.17); call pgsci(3)
c      call pgline(npts,real(zx),real(CDF3));
c      call pgmove(3.80,0.08); call pgdraw(3.95,0.08); call pgsci(2)
c      call pgslw(1)
c      do i=1,nl
c         call PGMOVE(real(z4(i)),1.13)
c         call PGDRAW(real(z4(i)),1.08)
c      end do
c      call pgsci(4)
c      do i=1,nl
c         call PGMOVE(real(znc4(i)),1.06)
c         call PGDRAW(real(znc4(i)),1.01)
c      end do
c      call pgslw(3)
c      call pgsci(1); call pgsls(4)
c      call pgline(npts,real(zx),real(ncCDF1)); call pgsci(2)
c      call pgline(npts,real(zx),real(ncCDF2)); call pgsci(3)
c      call pgline(npts,real(zx),real(ncCDF3)); call pgsci(1)
c      
c      call pgslw(3); call pgsls(1)
c      call pgmtxt('b',-7.,.8,0.,'10\u12\d\(2243)N\dHI\u<10\u14\d')
c      call pgmtxt('b',-4.5,.8,0.,'10\u14\d\(2243)N\dHI\u<10\u17\d')
c      call pgmtxt('b',-2.,.8,0.,'10\u17\d\(2243)N\dHI\u\(2243)10\u22\d')
c      call pgmtxt('b',-12.4,0.6,0.0,'Clustering only')
c      call pgend
c
c      do i=1,3*npts
c         z=zstart+dble(i-1)*(zqso-zstart)/(3*npts-1)
c         o=omega(z)
c         omegas(i)=alog10(o)
c         dndxy(i)=(1.+o)**(-1.5)
c         write(6,*) i,z,omegas(i),dndxy(i)
c      end do
c      call pgenv(minval(omegas),maxval(omegas),0.0,2.0,0,0)
c      call pgsci(2)
c      call pgline(npts,omegas,dndxy)
c      call pgend
c
*     plot omega for the proximity effect
c      write(lbl,'(i6)') nint(zqso*1e5)
c      write(6,*) lbl
c      name='./plots/omega_z_'//lbl//'.ps'
c      call pgbegin(0,'/xserve',1,1)
c      call pgbegin(0,name//'/cps',1,1)
c      call pgslw(3)
c      call pgenv(0.90*real(zqso),1.05*real(zqso),0.0,1.1,0,0)
c      call pglabel('\fiz','\fn(1+\gw)\u-1\d','')
c      call pgsci(2); call pgsls(4); call pgslw(5)
c      call pgmove(real(zqso),0.0); call pgdraw(real(zqso),1.1)
c      call pgsci(1); call pgsls(1); call pgslw(3)
c      do i=1,npts
c         plotx(i)=real(zx(i))
c         ploty(i)=real((1.+omega(zx(i)))**(-1.0))
c      end do
c      call PGLINE(npts,plotx,ploty)
c      call PGTEXT(real(zqso-0.005),1.13,'\fiz\dQSO\u')
c      call PGEND
*     print out a line list
c      write(6,*) '---------- L I N E     L I S T ----------'
c      do i=1,nl
c         write(6,*) i, nhi4(i), z4(i), npz4(i)
c      end do
      return
      end subroutine assign


!=======================================================================
      subroutine clustering(npt,xx,z1,z2,rmag,alpha,cCDF1,cCDF2,cCDF3)
!=================        C L U S T E R I N G       ====================
!=================  P R O X I M I T Y   E F F E C T ====================
!=======================================================================
      implicit none
      integer :: i,j,k,t,idum,nc,clin,npt
      integer,parameter :: ncl=20,p=100
      real :: xi0,v0,A,gamma,gammas(3),c,R,beta
      real*8 :: z,v,zmin,zmax,nc_total,c_total,np_total,zq,rmag,rm,alpha
      real*8 :: gp1,z1p1,z2p1,ymin,ymax,J0,omega,alp
      real*8 :: n,cl,h,nu,nunp,integrate,qsimp,tf,gasdev3
      real*8 :: ran3,z1,z2,dz,sum,np_sum
      real*8,dimension(npt) :: xx,cCDF1,cCDF2,cCDF3
      real*8,dimension(npt) :: npCDF1,npCDF2,npCDF3
      real*8,dimension(npt) :: xy,yy,ncCDF
      real*8,dimension(p) :: zx,ncy,cy
      real*8,dimension(1000):: zcl
      character :: string*10, string2*10

*     common blocks: 'cls'  - cluster redshifts
*                    'vars' - clustering parameters & c
*                    'qso'  - quasar redshift, magnitude
      common /cls/ zcl
      common /vars/ xi0,v0,A,gamma,c,R
      common /qso/ zq,rm,alp
      
      external n,cl,h,nu,nunp,integrate,ran3,qsimp,tf,gasdev3

*     set parameters 
      data c/299792.458/ ! [km/s]
      data gammas/1.51,2.16,1.33/
      
      zq=1.0*z2
      rm=rmag
      alp=alpha
      
      A=33.4
      xi0=0.5
      v0=150.0

      idum=time()
      zmin=1.0*z1
      zmax=1.0*z2
*     calculate the number of clusters in redshift limits    
      t=nint(integrate(cl,zmin,zmax))
*     include counting error
      t=nint(t+sqrt(real(t))*gasdev3(idum))
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
*     NOTATION:
*     nc = no clustering + proximity effect
*     c = clustering + proximity effect
*     np = clustering + no proximity effect

      nc_total=integrate(n,z1,z2)
      c_total=integrate(nu,z1,z2)
      
c      write(6,*) 'No of lines (no clustering): ',nc_total
c      write(6,*) 'No of lines (clustering): ',c_total
*     construct CDFs by evaluating redshift z at npt points
      do i=1,npt
         xx(i)=z1+dble(i-1)*(z2-z1)/(npt-1)
      end do
      sum=0.0 ; np_sum=0.0
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
      sum=0.0 ; np_sum=0.0
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
      sum=0.0 ; np_sum=0.0
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

      call PGBEG(0,'/null',1,1)
      call PGENV(real(z1),real(z2),0.0,1.0,0,0)
      call PGSCI(2)
      call PGLINE(npt,real(xx),real(cCDF1))
      call PGSCI(5)
      call PGLINE(npt,real(xx),real(cCDF2))
      call PGSCI(7)
      call PGLINE(npt,real(xx),real(cCDF3))
      call PGEND
      end subroutine clustering
!=======================================================================
      function n(z)
!=======================================================================
*     differential distribution of HI absorbers per unit redshift, dn/dz
*     simple power-law model w/o clustering or proximity effect included
      common /vars/ xi0,v0,A,gamma,c,R
      real*8 n,z
      real A,gamma
      n=A*(1+z)**gamma
      return
      end function n
!=======================================================================
      function xi(z)
!=======================================================================
*     function to calculate the Two Point Correlation Function amplitude
*     at redshift z
      real*8 z,xi,beta,B,Nmin
      common /vars/ xi0,v0,A,gamma,c,R
      beta=1.6
      B=10**(9.57)
      Nmin=10**(12.75)
      xi=(1.+z)**gamma*(beta-1.0)/B*v0*(Nmin**(beta-1.0))
      !beta, B, Nmin from Kim et al. 2013: table 3, logNHI=12.75-18.0
      return 
      end function xi
!=======================================================================
      function cl(z)
!=======================================================================
*     differential distribution of clusters per unit redshift, dn(cl)/dz
*     see Barcons et Webb 1991
      common /vars/ xi0,v0,A,gamma,c,R
      real*8 cl,z,xi
      real xi0,v0,c
      external xi
      cl=1./(4*xi(z)*(1.+z)*(v0/c))
      return
      end function cl
!=======================================================================
      function nu(z)
!=======================================================================
*     differential distribution of HI absorbers per unit redshift dn/dz
*     model that includes clustering along the line of sight (1D) and
*     the proximity effect close to the QSO (1D)
      integer i,k,n,idum
      real*8 nu,z,int,z_cl(1000),ran3,tmp,omega,Fnu,J0,Jnu
      real*8 rm,alp,zq,pi,lumdist,m0,zap,om,xi
      external ran3,lumdist,uvbkg,xi
      data pi/3.14159265358979/
      common /vars/ xi0,v0,A,gamma,c,R
      common /cls/ z_cl
      common /qso/ zq,rm,alp
*     QSO redshift - zq
*     QSO r magnitude - rm
*     QSO spectral index alpha - alp

*     cloud clustering model (Barcons et al. 1991)
*     R = radius of the cluster (in redshift units)
*     xi0 = two-point correlation function, evaluated at v=0
      R=v0/c*(1+z)

*     int = auxilliary variable for calculating contribution from all
*     clusters
      int=0.0
      do i=1,1000
         tmp=z_cl(i)
         int=int+exp(-abs(z-tmp)/R)*4*xi0*(1.+z)**gamma
             !-0.5 is from 1-beta, beta=1.6 (Kim et al. 2013)
      end do
      nu=int
c      write(6,*) 'z=',z,'omega=',om,'dn/dX=',(1+om)**(-0.5) 
      return
      end function nu
!=======================================================================
      function nunp(z)      
!=======================================================================
*     the same as nu, but with no proximity effect, nunp = NU No Proximity
*     for testing purposes
      integer i,k,n,idum
      real*8 nunp,z,aux,z_cl(1000),ran3,tmp,omega,Fnu,J0,Jnu
      real*8 rm,alp,zq,pi,lumdist,m0,xi
      external ran3,lumdist,uvbkg,xi
      data pi/3.14159265358979/
      common /vars/ xi0,v0,A,gamma,c,R
      common /cls/ z_cl
      common /qso/ zq,rm,alp
      
      R=v0/c*(1+z)
      aux=0.0

      do i=1,1000
         tmp=z_cl(i)
         aux=aux+exp(-abs(z-tmp)/R)*4*xi0*A*(1.+z)**gamma
      end do
      nunp=aux
      return
      end function nunp
!=======================================================================
      function lumdist(z)
!=======================================================================
*     Calculate the luminosity distance of an object at a redshift z
      real*8 lumdist,z,E,integrate,OM,OR,OL,H0
      common /cosmology/ OM,OR,OL,H0
      common /vars/ xi0,v0,A,gamma,c,R
      external E,integrate

      lumdist=(1.+z)*c/H0*integrate(E,0d0,z)
      return
      end function lumdist
!=======================================================================
      function E(z)
!=======================================================================
      real*8 E,z,OM,OR,OL,H0
      common /cosmology/ OM,OR,OL,H0
      
*     from PYTHON : astropy.cosmology - Planck 2013
      OM=0.307
      OR=5.38412049426e-05
      OL=0.691391253393
      H0=67.8  ![km/s/Mpc]
      E=1./sqrt(OM*(1.+z)**3+OR*(1.+z)**2+OL)
      return
      end function E

!=======================================================================
      function integrate(f,a,b)
!=======================================================================
*     a __very__ simple integration routine
      implicit none
      integer,parameter :: npt=3000
      real*8 :: integrate, sumup,sumlow,sum,dx,dy
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
         dy=0.5*(y(i)+y(i-1))
         sumup=sumup+dy*dx
c         sumlow=sumlow+f(x(i+1))*dx
      end do
c      sum=0.5*(sumup+sumlow)
      integrate = sumup
      return
      end function integrate
!=======================================================================
      function uvbkg(z)
!=======================================================================
*     function to calculate the UV background at redshift z using data
*     from CUBA (Haardt & Madau 2012), publicly available on the webpage
*     http://www.ucolick.org/~pmadau/CUBA/Media/UVB.out
*     function is modelled by two Gaussians and a single lognormal distr

      real*8 z,gauss1,gauss2,mu1,mu2,mu3,sigma1,sigma2,sigma3
      real*8 C1,C2,C3,logn,uvbkg

      C1=1.90313321e-22
      mu1=2.17101606e+00
      sigma1=9.23796905e-01
      gauss1=C1*exp(-(z-mu1)**2/(2.*sigma1**2))
      
      C2=9.64175755e-23
      mu2=4.64514525e+00 
      sigma2=1.44949091e+00
      gauss2=C2*exp(-(z-mu2)**2/(2.*sigma2**2))

      C3=1.31483688e-22
      mu3=7.63122677e-01
      sigma3=5.96818518e-01
      logn=C3*exp(-(alog(real(z))-mu3)**2/(2.*sigma3**2))

      uvbkg=gauss1+gauss2+logn
      return
      end function uvbkg
!=======================================================================
      subroutine locate(npt,xbuf,ybuf,y, x)
!=======================================================================
*     For a given value of 'y', search the array 'ybuf' to find 
*     the element that is closest in value. Subroutine reads and returns
*     the corresponding value 'x' of the pair (x,y) from 'xbuf'.

*     INPUT,real*8 :        xbuf(npt), ybuf(npt), y
*     OUTPUT,real*8 :       x
!=======================================================================
      integer i,j,npt,index
      real*8,dimension(npt) :: xbuf, ybuf, d
      real*8 :: x,y,dmin

      d=abs(y-ybuf)
      index=minloc(d,dim=1)
      x=xbuf(index)
      return
      end subroutine locate

!=======================================================================
      function omega(z)
!=======================================================================
*     Calculate the proximity effect's omega, formula (6) Bajtlik et al. 1991
      real*8 z,omega,Jnu,Fnu,rmag,m0,f6182,f912,const,zap,uvbkg,pi
      real*8 alp,zq,rm,z_cl(1000),lumdist
      real xi0,v0,A,gamma,c,R
      external ran3,lumdist,uvbkg
      data pi/3.14159265358979/
      common /vars/ xi0,v0,A,gamma,c,R
      common /cls/ z_cl
      common /qso/ zq,rm,alp

*     conversion from SDSS r magnitude into physical flux
*     m0 is the referent magnitude
*     conversion: 1e-23 erg s-1 cm-2 Hz-1 --> 1e-23 erg s-1 cm-2 A-1
      m0 = 3631*1e-23!*((c*1e3)/(0.6182e-6)**2)*1e-10 ! 
      m0 = log10(m0)/0.4                                  
      f6182=10**(-(rm-m0)/2.5)
      const=f6182/(1.0/6182.0**(2.0+alp))       
*     continuum QSO flux is modelled by a power-law:
*     f_lambda = const * lambda**-(2+alpha)       [erg/s/cm^2/Hz]
*     physical flux from the QSO at the Lyman Limit (LL) i.e. 912 A
      f912=f6182*(6182.0/(912.0*(1.+z)))**(2.0+alp)
*     zap = redshift of the QSO as seen by an observer at redshift z
      zap = (zq+1.)/(z+1.)-1
*     Fnu = flux from the QSO at the LL           [erg/s/cm^2/Hz]
      Fnu=f912*(lumdist(zq)/lumdist(zap))**2/(1.+zq)   
c      write (6,*) 'Luminosity distance (Earth):',lumdist(zq), zq
c      write (6,*) 'Luminosity distance (cloud):',lumdist(zap), z
*     Jnu = UV background (at LL and redshift z) [erg/s/cm^2/Hz/sr]
      Jnu=uvbkg(z)!*(c*1e3)/(912e-10)**2*1e-10
*     omega = a measure of the proximity effect (see Bajtlik et al. 1988)
      omega=Fnu/(4*pi*Jnu)
c      write(6,*) 'm0 = ',m0
c      if (z.gt.0.99*zq) then
c         write(6,*) 'zq,zap,z',zq,zap,z
c         write(6,*) 'Fnu,Jnu =',fnu,jnu
c         write(6,*) 'f6182,f912 =',f6182,f912
c         write(6,*) 'omega =',omega
c      end if
      return
      end function omega

!=======================================================================
      function dndx(z)
      real*8 omega,z,dndx
      external omega
      
      dndx=(1.+omega(z))**(-0.5)

      return
      end function
!=======================================================================
      subroutine no_clustering(npt,xx,z1,z2,ncCDF1,ncCDF2,ncCDF3)
      
      integer npt
      real*8,dimension(npt) :: xx, ncCDF1, ncCDF2, ncCDF3
      real*8 :: z1,z2

*     case 12<logNHI<14
      gamma=1.51
      p=(1.+z1)**(gamma+1.)
      q=(1.+z2)**(gamma+1.)      
      ncCDF1=((1.+xx)**(gamma+1.)-p)/(q-p)
*     case 14<logNHI<17
      gamma=2.16
      p=(1.+z1)**(gamma+1.)
      q=(1.+z2)**(gamma+1.)      
      ncCDF2=((1.+xx)**(gamma+1.)-p)/(q-p)
*     case 17<logNHI<22
      gamma=1.33
      p=(1.+z1)**(gamma+1.)
      q=(1.+z2)**(gamma+1.)      
      ncCDF3=((1.+xx)**(gamma+1.)-p)/(q-p)
      call PGBEG(0,'/null',1,1)
      call PGENV(real(z1),real(z2),0.0,1.0,0,0)
      call PGSCI(2)
      call PGLINE(npt,real(xx),real(ncCDF1))
      call PGSCI(5)
      call PGLINE(npt,real(xx),real(ncCDF2))
      call PGSCI(7)
      call PGLINE(npt,real(xx),real(ncCDF3))
      call PGEND
      return
      end subroutine

      
