MODULE test_module
      CONTAINS
         SUBROUTINE get_data(data_path, data_array)
            ! Purpose:
            !    Hide boilerplate code for data procurment from main program
            ! Record of revisions:
            !    DATE       PROGRAMMER     DESCRIPTION OF CHANGE
            !    ====       ==========     =====================
            !  2/19/2019    DEREK TROPF    DTYPE REAL -> REAL*8
            USE netcdf
            CHARACTER(LEN = *), INTENT(IN) :: data_path
            REAL*8, ALLOCATABLE, INTENT(INOUT) :: data_array(:,:)
            INTEGER :: dim1, dim2, ncid, status, dim1len, dim2len
            CHARACTER(LEN = 30) :: dimname

            ! Open file
            status = NF90_OPEN(path = data_path, mode = NF90_NOWRITE, ncid = ncid)

            ! Obtain dimension lengths of data
            status = NF90_INQUIRE_DIMENSION(ncid, 1, dimname, dim1len)
            status = NF90_INQUIRE_DIMENSION(ncid, 2, dimname, dim2len)
            PRINT*,'dim1len and dim2len', dim1len, dim2len

            ! Allocate array, obtain data, and close file
            ALLOCATE(data_array(dim2len, dim1len))
            status = NF90_GET_VAR(ncid, 1, data_array)
            PRINT*, 'Obtained Data'
            status = NF90_CLOSE(ncid)
         END SUBROUTINE get_data

         SUBROUTINE save_data(data_path, n_iter, n_cent, init_centroids, data_array, counts, sse, movement)
            ! Purpose
            ! -------
            ! Save data from k-means clustering.
            !
            ! Parameters
            ! ----------
            ! data_path: string
            !   Location and name to save data
            !
            ! n_iter: integer
            !   Number of iterations done for k-means clustering
            !
            ! n_cent: integer
            !   Number of centroids (k in k-means)
            !
            ! init_centroids: 2D array of Real*8
            !   Initialization centroids used to start k-means clustering
            !
            ! data_array: 2D array of Real*8
            !   Centroids produced from k-means clustering
            !
            ! counts: 2D array of integer
            !   Number of members in each cluster for each iteration
            !
            ! sse: 2D array of Real*8
            !   Sum of Squared Errors for each cluster for each iteration
            !
            ! movement: 2D array of Real*8
            !   Movement of each centroid after clustering each iteration
            USE netcdf
            CHARACTER(LEN = *), INTENT(IN) :: data_path
            INTEGER, INTENT(IN) :: n_iter, n_cent
            REAL*8, INTENT(IN) :: init_centroids(:,:), data_array(:,:), movement(:,:), sse(:,:)
            INTEGER, INTENT(IN) :: counts(:,:)
            INTEGER :: ncid, status, n_iterdim, n_centdim, pctau_binsdim, unlimiteddim
            INTEGER :: init_id, cent_id, sse_id, count_id, movement_id

            ! Create file to save data to
            PRINT*, 'OPENING NETCDF FILE TO SAVE DATA'
            status = NF90_CREATE(data_path, NF90_NETCDF4, ncid)

            ! Number of pctau bins, hardcoded because constant
            PRINT*, 'CREATING DIMENSIONS'
            status = NF90_DEF_DIM(ncid, 'pctau_bins', 42, pctau_binsdim)
            status = NF90_DEF_DIM(ncid, 'n_centroids', n_cent, n_centdim)
            status = NF90_DEF_DIM(ncid, 'n_iterations', n_iter, n_iterdim)

            ! Create variables
            PRINT*, 'CREATING VARIABLES'
            status = NF90_DEF_VAR(ncid, 'init_centroids', NF90_DOUBLE, (/pctau_binsdim, n_centdim/), init_id)
            status = NF90_DEF_VAR(ncid, 'centroids', NF90_DOUBLE, (/pctau_binsdim, n_centdim/), cent_id)
            status = NF90_DEF_VAR(ncid, 'counts', NF90_INT, (/n_centdim, n_iterdim/), count_id)
            status = NF90_DEF_VAR(ncid, 'sse', NF90_DOUBLE, (/n_centdim, n_iterdim/), sse_id)
            status = NF90_DEF_VAR(ncid, 'movement', NF90_DOUBLE, (/n_centdim, n_iterdim/), movement_id)

            ! Save data
            PRINT*, 'SAVING DATA'
            status = NF90_PUT_VAR(ncid, init_id, init_centroids(:,:))
            status = NF90_PUT_VAR(ncid, cent_id, data_array(:,:))
            status = NF90_PUT_VAR(ncid, count_id, counts(:,:))
            status = NF90_PUT_VAR(ncid, sse_id, sse(:,:))
            status = NF90_PUT_VAR(ncid, movement_id, movement(:,:))

            ! Close file
            status = NF90_CLOSE(ncid)
         END SUBROUTINE save_data


         SUBROUTINE kmeans(nvec, ncent, vectors, centroids, cluster_totals, counts, sse)
            ! Purpose
            ! -------
            ! Cluster HGG data based on smallest Euclidean Distance.
            ! Calculate centroid, number of members, and Sum of Squared Errors
            ! for each cluster
            !
            ! Parameters
            ! ----------
            ! nvec: integer
            !    Number of HGG gridboxes (vectors of shape (42))
            !
            ! ncent: integer
            !   Number of centroids used for clustering/number of clusters
            !
            ! vectors: 2D array of Real*8
            !   HGG data used for clustering.
            !
            ! centroids: 2D array of Real*8
            !   Initial centroids used to start clustering
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

            ! Clustering
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
END MODULE test_module
