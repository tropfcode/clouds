SUBROUTINE work(E1ice,E1bed,axyp,dxyp,topo,bed,land,ice,zatmo)
      IMPLICIT none
      INTEGER, PARAMETER :: imE1 = 21600, jmE1 = 10800, nodata = -32768
      INTEGER, PARAMETER :: row1 = 288, col1 = 180, hyper = 0
      INTEGER, PARAMETER :: row2 = 144, col2 = 90, nx = 75, ny = 60
      INTEGER, INTENT(IN) :: E1ice(imE1,jmE1)
      INTEGER, INTENT(IN) :: E1bed(imE1,jmE1)
      REAL, INTENT(IN) :: axyp(row1,col1)
      REAL, INTENT(IN) :: dxyp(jmE1)
      REAL, INTENT(OUT) :: topo(row1,col1), bed(row1,col1)
      REAL, INTENT(OUT) :: land(row1,col1), ice(row1, col1)
      REAL, INTENT(OUT) :: zatmo(row1,col1)
      INTEGER :: ice_subgrid(nx,ny), bed_subgrid(nx,ny)
      REAL :: dxyp_slices(nx,ny,col1)
      REAL :: dxyp_slice(nx,ny)
      REAL :: ice_slice(nx,ny), bed_slice(nx,ny)
      REAL :: ice_fraction(nx,ny), bed_fraction(nx,ny)
      REAL :: ice_sumweight, bed_sumweight
      REAL :: icebed_diff(nx,ny)
      REAL :: landpoint, icepoint, tmparray(nx,ny)
      INTEGER :: nslices, start_col, end_col
      INTEGER :: start_row, end_row, i, j
      
      ! TEMP ARRAY SETUP FOR TESTING
      topo = 0
      bed = 0
      land = 0
      ice = 0
      
      ! Transform dxyp into 3D array copying values
      ! to produce slices, with each slice used for
      ! each nx,ny section of E1ice,E1bed
      dxyp_slices = 0
      nslices = jmE1/ny
      DO i=1,nslices
         DO j=1,nx
            dxyp_slices(j,:,i) = dxyp(i*ny:(i+1)*ny)
         END DO
      END DO
      
      DO i=1,SIZE(dxyp_slices,3)-1
         start_col = i*ny
         end_col = (i+1)*ny
         dxyp_slice = dxyp_slices(:,:,i)
         DO j=1,SIZE(E1ice,1)/nx
            start_row = j*nx
            end_row = (j+1)*nx

            ! Get subset of gridded data
            ice_subgrid = E1ice(start_row:end_row,start_col:end_col)
            bed_subgrid = E1bed(start_row:end_row,start_col:end_col)

            ! Make arrays filled with 0s and dxyp values for bed and ice
            WHERE(ice_subgrid.EQ.nodata)
               ice_slice = 0.0
            ELSEWHERE
               ice_slice=dxyp_slice
            END WHERE
            WHERE(bed_subgrid.EQ.nodata)
               bed_slice = 0.0
            ELSEWHERE
               bed_slice = dxyp_slice
            END WHERE

            ice_sumweight = SUM(ice_slice)
            bed_sumweight = SUM(bed_slice)
            ice_fraction = ice_slice/ice_sumweight
            bed_fraction = bed_slice/bed_sumweight

            ! calculate landpoint
            tmparray = 0.0
            WHERE(ice_subgrid.GE.0.0)
               tmparray = ice_fraction
            END WHERE
            landpoint = SUM(tmparray)
            land(j,i) = landpoint
            
            ! Calculate icepoint
            icebed_diff = (bed_subgrid-ice_subgrid)
            tmparray = 0.0
            WHERE(icebed_diff.LT.(0.000001) .AND. icebed_diff.GT.(-0.000001))
               tmparray = ice_fraction
            END WHERE
            icepoint = SUM(tmparray)
            ice(j,i) = icepoint
            
            ! Calculate topo
            topo(j,i) = SUM(ice_fraction*ice_subgrid)
            bed(j,i) = SUM(bed_fraction*bed_subgrid)
         END DO
      END DO
      
      ! For zatmo
      zatmo = 0
      WHERE(land.GT..5)
         zatmo = topo
      END WHERE
      PRINT*,'END OF WORK PROGRAM YAY'
END SUBROUTINE work
