MODULE cluster_subs
CONTAINS
   SUBROUTINE get_data(data_path, data_array)
      ! Purpose:
      !    Read centroids from netcdf file
      !
      ! Argument description:
      !    data_path - Full path to centroid data used for initalizing k-means clustering
      !    data_array - Loaded initialization centroids from full path

      USE netcdf
      CHARACTER(LEN = *), INTENT(IN) :: data_path
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: data_array(:,:)
      INTEGER :: status, ncid, dim1len, dim2len ! status for netcdf library, file id, first and second dimension length
      CHARACTER(LEN = 30) :: dimname ! name of dimension from netcdf file
      PRINT*,'FROM GET_DATA FUNCTION PRINTING DATA_PATH', data_path
 
      ! Open file
      status = NF90_OPEN(path = data_path, mode = NF90_NOWRITE, ncid = ncid)
  
      ! Obtain dimension lengths of data
      status = NF90_INQUIRE_DIMENSION(ncid, 1, dimname, dim1len)
      status = NF90_INQUIRE_DIMENSION(ncid, 2, dimname, dim2len)
      PRINT*,'dim1len and dim2len', dim1len, dim2len

      ! Allocate array, obtain data, and close file
      ! Note that dim2len and dim1len are flipped in allocation due to Fortran column major order
      ALLOCATE(data_array(dim2len, dim1len))
      status = NF90_GET_VAR(ncid, 1, data_array)
      PRINT*, 'OBTAINED DATA'
      status = NF90_CLOSE(ncid)
   END SUBROUTINE get_data
  
   SUBROUTINE save_data(file_path, n_iter, n_cent, init_cent, final_cent, sse, movement, counts)
      ! Purpose:
      !    Save results of k-means clustering into single netcdf file.
      !    File will contain initial centroids used in clustering, final centroids, 
      !    cluster counts for each iteration, sum of squared erros (SSE) for each iteration,
      !    and centroid movement for each iteration.
      !
      ! Argument description:
      !    file_path - Location of where to save cluster results
      !    n_iter - Number of iterations from cluster run
      !    n_cent - Number to centroids from cluster run
      !    init_cent - Initialization centroids for clustering
      !    final_cent - Final centroids resulting from clustering
      !    sse - Sum of Squared Erros for each iteration
      !    movement - Amount each centroid has moved for each iteration, calculated using Euclidean Distance
      !    counts - Membership count of each cluster for each iteration

      USE netcdf
      CHARACTER(LEN = *), INTENT(IN) :: file_path
      INTEGER, INTENT(IN) :: n_iter, n_cent
      REAL*8, INTENT(IN) :: init_cent(:,:), final_cent(:,:)
      REAL*8, INTENT(IN) :: sse(:,:), movement(:,:) 
      INTEGER, INTENT(IN) :: counts(:,:) 
      INTEGER :: status, ncid, init_id, final_id, sse_id, count_id, movement_id ! Variables used by netcdf library 
      INTEGER :: n_iterdim, n_centdim, pctau_binsdim, unlimiteddim ! Dimensions for netcdf
  
      ! Create file to save data to
      PRINT*, 'OPENING NETCDF FILE TO SAVE DATA'
      status = NF90_CREATE(file_path, NF90_NETCDF4, ncid)
  
      ! Create netcdf dimensions used by variables
      ! Number of pctau bins, hardcoded because constant
      PRINT*, 'CREATING DIMENSIONS'
      status = NF90_DEF_DIM(ncid, 'pctau_bins', 42, pctau_binsdim)
      status = NF90_DEF_DIM(ncid, 'n_centroids', n_cent, n_centdim)
      status = NF90_DEF_DIM(ncid, 'n_iterations', n_iter, n_iterdim)
  
      ! Create variables
      PRINT*, 'CREATING VARIABLES'
      status = NF90_DEF_VAR(ncid, 'init_centroids', NF90_DOUBLE, (/pctau_binsdim, n_centdim/), init_id)
      status = NF90_DEF_VAR(ncid, 'centroids', NF90_DOUBLE, (/pctau_binsdim, n_centdim/), final_id)
      status = NF90_DEF_VAR(ncid, 'counts', NF90_INT, (/n_centdim, n_iterdim/), count_id)
      status = NF90_DEF_VAR(ncid, 'sse', NF90_DOUBLE, (/n_centdim, n_iterdim/), sse_id)
      status = NF90_DEF_VAR(ncid, 'movement', NF90_DOUBLE, (/n_centdim, n_iterdim/), movement_id)
  
      ! Save data
      PRINT*, 'SAVING DATA'
      status = NF90_PUT_VAR(ncid, init_id, init_cent(:,:))
      status = NF90_PUT_VAR(ncid, final_id, final_cent(:,:))
      status = NF90_PUT_VAR(ncid, count_id, counts(:,:))
      status = NF90_PUT_VAR(ncid, sse_id, sse(:,:))
      status = NF90_PUT_VAR(ncid, movement_id, movement(:,:))
  
      ! Close file
      status = NF90_CLOSE(ncid)
   END SUBROUTINE save_data
  
  
   SUBROUTINE kmeans(nvec, ncent, vectors, centroids, cluster_totals, counts, sse)
      ! Purpose:
      !    Cluster ISCCP PC-TAU data using k-means method for a single iteration.
      !    Provide centroids, each of which defines a cluster centroid.
      !    For each PC-TAU histogram, the smallest Euclidean distance is found to determine cluster membership.
      !    The sum of all vectors in a cluster and cluster membership count are returned to calculate new centroids.
      !    Sum of Squared Distance (SSE) for each cluster is calculated and returned.
      !
      ! Argument description:
      !    nvec - Total number of vectors used in clustering.
      !    ncent - Total number of centroids.
      !    vectors - 2D matrix of flattend PC-TAU histograms. Each row is a unique histogram, each column a value in histogram.
      !    centroids - Centroids used to define number of clusters and for testing cluster membership of vectors.
      !    cluster_totals - Sum of member vectors in each cluster.
      !    counts - Total number of members in each cluster.
      !    sse - Sum of Squared Distance. SSE is summed distance of member centroids for a given cluster.

      INTEGER, INTENT(IN) :: nvec, ncent
      REAL*8, INTENT(IN) :: vectors(:,:)
      REAL*8, INTENT(IN) :: centroids(:,:)
      REAL*8, INTENT(OUT) :: cluster_totals(42,ncent)
      INTEGER, INTENT(OUT) :: counts(ncent)
      REAL*8, INTENT(OUT) :: sse(ncent)
      REAL*8 :: dist(ncent), start_time, end_time
      INTEGER :: i, j, ml
  
      PRINT*,'BEGIN CLUSTERING'
  
      ! Initialize variables and begin clock to time subroutine
      sse = 0.0
      counts = 0
      cluster_totals = 0.0
      CALL CPU_TIME(start_time)
  
      DO i=1,nvec
         dist = 0
         ! Find centroid that produces smallest euclidean distance
         DO j=1, ncent
            dist(j) = SQRT(SUM((vectors(:,i)-centroids(:,j))**2))
         END DO
         ! Get index of closest centroid
         ml = MINLOC(dist, DIM = 1)
         ! Total distance all members are from their centroid
         sse(ml) = sse(ml) + dist(ml)
         ! Keep a running sum of all members in each cluster
         cluster_totals(:,ml) = cluster_totals(:,ml)+vectors(:,i)
         ! Count the number of members in each cluster
         counts(ml) = counts(ml) + 1
      END DO
      CALL CPU_TIME(end_time)
      PRINT*,'TOTAL TIME FOR CLUSTER', end_time-start_time
   END SUBROUTINE kmeans
END MODULE cluster_subs
