subroutine DXYP_calc(im,jm,halfpole,dxyp)
! Calculates the surface area per gridbox, dxyp, given a resolution of lon(im) and lat(jm)
      real*8, parameter :: pi = 3.1415926535897932d0
      real*8 :: twopi = 2d0*pi
      real*8 :: radian = pi/180d0
      real*8 :: radius = 6371000.
      integer, intent(in) :: im,jm,halfpole
      real*8, intent(out) :: dxyp(jm)
      real*8 fjeq,SINV,SINVm1,DLAT,DLAT_dg,DLON
      integer j

!     ModelE calculation
      FJEQ=0.5*(1+JM)
      DLON=TWOPI/(IM*1.)
     
      DLAT_DG=180./REAL(JM) 

      if (halfpole.gt.0) then
        DLAT_DG=180./REAL(JM-1) ! 1/2 box at pole      
        print*,"Assuming half box at the pole"
      endif

      DLAT=DLAT_DG*radian

      j=1
      SINV    = Sin (DLAT*(j+.5-FJEQ))
      DXYP(j) = RADIUS*RADIUS*DLON*(SINV+1)

      do j=2,jm-1
        SINVm1  = Sin (DLAT*(J-.5-FJEQ))
        SINV    = Sin (DLAT*(J+.5-FJEQ))
        DXYP(J) = RADIUS*RADIUS*DLON*(SINV-SINVm1)
      enddo

      j=jm
      SINVm1  = Sin (DLAT*(J-.5-FJEQ))
      DXYP(J)= RADIUS*RADIUS*DLON*(1-SINVm1)
!      return
END subroutine DXYP_calc
