      subroutine vp_spread(flx,nlw,nw,ichunk, flb)
*
      implicit none
      include 'vp_sizes.f'
*
*     Convolve array with Gaussian (or other) profile. 
*     If local sigma is zero, just copies
*     
*     IN:
*     flx	r*4 array	data array
*     nlw	int		start channel for result
*     nw	int		end channel for result
*     ichunk	int		chunk number (for wavelengths)
*
*     OUT:
*     flb	r*4 array	flx convolved with instrument profile
*
*     COMMON:
*     Binned profile parameters, /vpc_profwts/, used if npinst.gt.0
*
*
      integer nlw,nw,ichunk
      double precision flx(*),flb(*)
*     double precision sighxx,res
*
*     LOCAL:
      integer i,j,k,nd,nlo,nup
      double precision bwhvard,df,fsigd,sigd,xdb
*     data value variables
      double precision s,sx,y,sm,smn
*
*     FUNCTIONS:
      double precision vp_wval,vpf_dvresn,dexpf
*     double precision erfcc
*     binned profile
      character*2 chprof
      double precision pfinst
      double precision wfinstd
      integer npinst,npin2
      common/vpc_profwts/wfinstd(maxipp,maxnch),pfinst(maxipp,maxnch),
     :        npinst(maxnch),npin2(maxnch),chprof(maxnch)
*
      double precision dnshft2(maxnch)
      common/vpc_shft2/dnshft2
*
*     smoothed data in array flx -> flb(nlw,nw)
*     
*     use binned profile if npinst.gt.0
      if(npinst(ichunk).gt.0.and.chprof(ichunk)(1:2).ne.'wt') then
*       pixel profile
        do i=nlw,nw
          sm=0.0d0
          smn=0.0d0
          do j=1,npinst(ichunk)
            k=i+j-npin2(ichunk)
            if(k.ge.1.and.k.le.nw) then
              sm=sm+flx(k)*pfinst(j,ichunk)
              smn=smn+pfinst(j,ichunk)
            end if
          end do
          if(smn.gt.0.0d0) then
            flb(i)=sm/smn
           else
*           raw value if weight is zero
            flb(i)=flx(i)
          end if
        end do
       else
*       no pixel instrument profile, so use continuous function
*       gaussian smoothing
        do i=nlw,nw
*         set the local sigma
          xdb=vp_wval(dble(i)+dnshft2(ichunk),ichunk)
          bwhvard=0.5d0*(vp_wval(dble(i+1)+dnshft2(ichunk),
     :              ichunk) - vp_wval(dble(i-1)+
     :                dnshft2(ichunk),ichunk))
          sigd=vpf_dvresn(xdb,ichunk)
*         for sigd zero, don't smooth
          if(sigd.gt.0.0d0) then
            fsigd=6.0d0*sigd
            df=fsigd*xdb/bwhvard
            nd=int(df)
            if(nd.le.0) nd=1
            nlo=i-nd
            nup=i+nd
            if(nlo.lt.1) nlo=1
            if(nup.gt.nw) nup=nw
            sx=0.0d0
            s=0.0d0
            do j=nlo,nup
              y=(vp_wval(dble(j)+dnshft2(ichunk),ichunk)/
     :                    xdb-1.0d0)/sigd
              y=dexpf(-y*y*0.5d0)
              sx=sx+y
              s=s+flx(j)*y
            end do
            if(sx.gt.0.0d0) then
              flb(i)=s/sx
             else
              flb(i)=flx(i)
            end if
           else
            flb(i)=flx(i)
          end if
        end do
      end if
      return
      end
