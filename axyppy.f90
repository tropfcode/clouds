subroutine AXYP_latlon_calc(im,jm,halfpole,axyp,lat,lon)

      real*8, parameter :: pi = 3.1415926535897932d0
      real*8, parameter :: twopi = 2d0*pi
      real*8, parameter :: radian = pi/180d0
      real*8, parameter :: radius = 6371000.
      integer, intent(in) :: im,jm, halfpole
      REAL*8, intent(out) :: axyp(im,jm)
      REAL*8, intent(out) :: lat(jm), lon(im)
      integer i,j
      REAL*8 dxyp(jm),fjeq,SINV,SINVm1      
      REAL*8 DLAT,DLAT_dg,DLON,DLON_dg


!     ModelE calculation
      FJEQ=0.5*(1+JM)
      DLON=TWOPI/(IM*1.)
      DLAT_DG=180./REAL(JM) 
      DLON_DG=360./REAL(IM)

      if (halfpole.gt.0) DLAT_DG=180./REAL(JM-1)   ! 1/2 box at pole

      DLAT=DLAT_DG*radian

      j=1
      SINV    = Sin (DLAT*(j+.5-FJEQ))
      DXYP(j) = RADIUS*RADIUS*DLON*(SINV+1)
      lat(j) = -90. + 0.5*DLAT_DG
      if (halfpole.gt.0) lat(j) = -90. + 0.25*DLAT_DG

      do j=2,jm-1
        SINVm1  = Sin (DLAT*(J-.5-FJEQ))
        SINV    = Sin (DLAT*(J+.5-FJEQ))
        DXYP(J) = RADIUS*RADIUS*DLON*(SINV-SINVm1)
        lat(j) = lat(j-1)+DLAT_DG
        if((j.eq.2).and.(halfpole.gt.0)) lat(j) = -90. + DLAT_DG
      enddo

      j=jm
      SINVm1  = Sin (DLAT*(J-.5-FJEQ))
      DXYP(J)= RADIUS*RADIUS*DLON*(1-SINVm1)
      lat(j) = lat(j-1)+DLAT_DG
      if (halfpole.gt.0) lat(j) = 90. - 0.25*DLAT_DG

      lon(1) = -180. + 0.5*DLON_DG
      do i=2,im
         lon(i)=lon(i-1)+DLON_DG
      enddo

      do i=1,im
         do j=1,jm
            axyp(i,j)=dxyp(j)
         enddo
      enddo
      
!      return
END subroutine AXYP_latlon_calc