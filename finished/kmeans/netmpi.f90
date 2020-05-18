! Program Description:
!    This program performs K-means clustering on ISCCP HGG PC-TAU histogram data in order to derive Weather States (WS).
!    The product made from this program is a netcdf file with initalization centroids, final centroids (or WS),
!    counts of membership for each cluster for each iteration, sum of squared errors (SSE) for each iteration,
!    and centroid movement for each cluster for each iteration.
!
! Program Implementation:
!    This program is meant to be ran on the NCCS NASA supercomputer and is submited as a batch job using SLURM via submit.sh file.
!    MPI is used to parallelize the program based on year -- there are as many processes as years of data.
!    Each process independetly runs the KMEANS subroutine from cluster_subs.f90. After completing the subroutine,
!    each process gathers cluster membership count (how many members are in each cluster), SSE (see cluster_subs.f90),
!    and cluster totals (sum of every member in each cluster). A new set of centroids is calculated for the next iteration
!    by dividing cluster totals by cluster membership count. Each iteration, centroid movement is calculated by
!    subtracting newly calculated centroids by the previous centroids.
!    This process repeats until niter number of iterations have been completed.
!
! Program Parameters:
!    niter : integer
!       Number of iterations to perform K-means clustering
!
!    ncent : integer
!       Number of centroids used in K-means clustering. Alternatively, the 'K' in K-means
!
!    startyear : integer
!       First year of data to use in clustering
!
!    init_path : character
!       Full path to file containing initialization centroids
!
!    data_path_base : character
!       Full path to directory containing data to be clustered
!
! Note:
!    This program assumes a shape of (N,42) for data to be clustered.
!    PC-TAU histograms in the ISCCP HGG dataset are of shape (7,6). However, this program requires that the
!    histograms are flattened to shape (42). If there are N histograms for a given year, then the final data
!    used in clustering is of shape (N,42) if all flattened histograms are stacked.

PROGRAM nettest
       USE cluster_subs
       USE MPI
       CHARACTER(LEN = 1024) :: fmt, datayear
       INTEGER :: numtasks, startyear, rank, len, ierr
       INTEGER :: i, j, k, ncent, niter, nvec
       CHARACTER(LEN = 100) :: cmd_arg(5), init_path, data_path, data_path_base
       REAL*8, ALLOCATABLE :: cluster_totals(:,:), movement(:,:), sse(:), sse_gather(:,:), sse_master(:,:)
       INTEGER, ALLOCATABLE :: counts(:), labels(:)
       REAL*8, ALLOCATABLE :: data_array(:,:), init_centroids(:,:), prev_centroids(:,:), centroids(:,:)
       REAL*8, ALLOCATABLE :: tmp_cluster(:,:), cluster_gather(:,:,:), centroid_bcast(:,:)
       INTEGER, ALLOCATABLE :: tmp_counts(:), counts_gather(:,:), counts_master(:,:)

       ! Setup MPI
       CALL MPI_INIT(ierr)
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

       ! Get comandline arguments
       CALL GET_COMMAND_ARGUMENT(1, cmd_arg(1))
       CALL GET_COMMAND_ARGUMENT(2, cmd_arg(2))
       CALL GET_COMMAND_ARGUMENT(3, cmd_arg(3))
       CALL GET_COMMAND_ARGUMENT(4, cmd_arg(4))
       CALL GET_COMMAND_ARGUMENT(5, cmd_arg(5))
       READ(cmd_arg(1),*) niter
       READ(cmd_arg(2),*) ncent
       READ(cmd_arg(3),*) startyear
       init_path = cmd_arg(4)
       data_path_base = cmd_arg(5)

       ! Format for acquiring correct data for each process
       fmt = '(I4)'
       WRITE(datayear,fmt) startyear+rank
       PRINT*, 'RANK YEAR NUMTASKS: ', rank, startyear+rank, numtasks

       ! Allocate Arrays
       ALLOCATE(cluster_totals(42,ncent))
       ALLOCATE(movement(ncent,niter))
       ALLOCATE(sse(ncent))
       ALLOCATE(counts(ncent))
       ALLOCATE(centroid_bcast(42,ncent))
       ALLOCATE(prev_centroids(42,ncent))
       ALLOCATE(centroids(42,ncent))
       ALLOCATE(tmp_cluster(42,ncent))
       ALLOCATE(tmp_counts(ncent))
       ALLOCATE(cluster_gather(42,ncent,numtasks))
       ALLOCATE(counts_gather(ncent,numtasks))
       ALLOCATE(sse_gather(ncent,numtasks))
       ALLOCATE(sse_master(ncent,niter))
       ALLOCATE(counts_master(ncent,niter))

       ! Get initial centroids
       PRINT*, 'OBTAINING INIT CENTROIDS', init_path
       CALL get_data(init_path, init_centroids)
       centroids(:,:) = init_centroids(:,:) 
       IF (rank .EQ. 0) THEN
          PRINT*,centroids(1,1:10)
       END IF
       IF (rank .EQ. 22) THEN
          PRINT*,centroids(1,1:10)
       END IF

       ! Obtain ISCCP data for single year
       data_path = trim(data_path_base)//trim(datayear)//'.nc'
       PRINT*, 'OBTAINING PCTAU DATA', data_path
       CALL get_data(data_path, data_array)
       nvec = SIZE(data_array, DIM = 2)

       ! Initalize variables
       sse_master = 0
       counts_master = 0

       DO i=1,niter
          !PRINT*,'ITERATION NUMBER AND RANK', i, rank
          !PRINT*,''
          sse = 0
          sse_gather = 0
          counts = 0
          cluster_gather = 0
          counts_gather = 0

          ! Cluster
          CALL kmeans(nvec, ncent, data_array, centroids, cluster_totals, counts, sse)

          ! Gather counts and cluster totals
          !PRINT*,'ABOUT TO DO MPI_ALLGATHER 1'
          CALL MPI_ALLGATHER(cluster_totals, 42*ncent, MPI_REAL8, cluster_gather, 42*ncent, MPI_REAL8, MPI_COMM_WORLD, ierr)
          !PRINT*,rank,'DID MPI_ALLGATHER'
          CALL MPI_ALLGATHER(counts, ncent, MPI_INTEGER, counts_gather, ncent, MPI_INTEGER, MPI_COMM_WORLD, ierr)
          !PRINT*,rank,'DID MPI_ALLGATHER 2'
          CALL MPI_ALLGATHER(sse, ncent, MPI_REAL8, sse_gather, ncent, MPI_REAL8, MPI_COMM_WORLD, ierr)
          !PRINT*,rank,'DID MPI_ALLGATHER 3'

          ! Keep track of previously used centroids
          prev_centroids(:,:) = centroids(:,:)
          ! Add up gathered values for cluster_totals
          tmp_cluster(:,:) = SUM(cluster_gather, DIM = 3)
          ! Add up gathered values
          tmp_counts = SUM(counts_gather, DIM = 2)
          ! Calculate new centroids and centroid movement
          DO j=1,ncent
             centroids(:,j) = tmp_cluster(:,j)/tmp_counts(j)
             movement(j,i) = SQRT(SUM((prev_centroids(:,j)-centroids(:,j))**2))
          END DO
          sse_master(:,i) = SUM(sse_gather, DIM = 2)
          counts_master(:,i) = tmp_counts(:)
          CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
       END DO

       ! Save all data
       IF (rank .EQ. 0) THEN
          CALL save_data('out'//trim(cmd_arg(2))//'.nc', niter, ncent, init_centroids, centroids, sse_master, movement, counts_master)
       END IF

       ! Deallocate Memory and finish program
       PRINT*, 'DEALLOCATING MEMORY'
       DEALLOCATE(cluster_gather)
       DEALLOCATE(counts_gather)
       DEALLOCATE(sse_gather)
       DEALLOCATE(cluster_totals)
       DEALLOCATE(sse)
       DEALLOCATE(counts)
       DEALLOCATE(init_centroids)
       DEALLOCATE(prev_centroids)
       DEALLOCATE(centroids)
       DEALLOCATE(data_array)
       DEALLOCATE(movement)
       DEALLOCATE(tmp_counts)
       DEALLOCATE(tmp_cluster) 
       DEALLOCATE(sse_master)
       DEALLOCATE(counts_master)

       ! Terminate MPI execution environment
       CALL MPI_FINALIZE(ierr)
END PROGRAM nettest
