PROGRAM for_kmeans
       USE test_module
       USE MPI
       CHARACTER(LEN = 1024) :: fmt, x1
       INTEGER :: numtasks, rank, len, ierr
       INTEGER :: i, j,k, ncent, niter, nvec
       CHARACTER(LEN = 30) :: cmd_arg(2), data_path
       REAL*8, ALLOCATABLE :: cluster_totals(:,:), movement(:,:), sse(:), sse_gather(:,:), sse_master(:,:)
       INTEGER, ALLOCATABLE :: counts(:), labels(:)
       REAL*8, ALLOCATABLE :: data_array(:,:), init_centroids(:,:), prev_centroids(:,:), centroids(:,:)
       REAL*8, ALLOCATABLE :: tmp_cluster(:,:), cluster_gather(:,:,:), centroid_bcast(:,:)
       INTEGER, ALLOCATABLE :: tmp_counts(:), counts_gather(:,:), counts_master(:,:)

       ! Setup MPI
       CALL MPI_INIT(ierr)
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

       ! Format for acquiring correct data for each process
       fmt = '(I4)'
       WRITE(x1,fmt) 1984+rank

       ! Get comandline arguments and alloate arrays
       CALL GET_COMMAND_ARGUMENT(1, cmd_arg(1))
       CALL GET_COMMAND_ARGUMENT(2, cmd_arg(2))
       READ(cmd_arg(1),*) niter
       READ(cmd_arg(2),*) ncent
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
       PRINT*, 'ARG2 THIS MUCH MATCH VALUE OF K', cmd_arg(2)
       CALL get_data('./newv3_ncdata/k'//trim(cmd_arg(2))//'.nc', init_centroids)
       centroids(:,:) = init_centroids(:,:) 

       ! Obtain ISCCP data for single year
       data_path = '../data/ncdata/'//trim(x1)//'.nc'
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
          CALL MPI_ALLGATHER(cluster_totals, 42*ncent, MPI_REAL8, cluster_gather, 42*ncent, MPI_REAL8, MPI_COMM_WORLD, ierr)
          CALL MPI_ALLGATHER(counts, ncent, MPI_INTEGER, counts_gather, ncent, MPI_INTEGER, MPI_COMM_WORLD, ierr)
          CALL MPI_ALLGATHER(sse, ncent, MPI_REAL8, sse_gather, ncent, MPI_REAL8, MPI_COMM_WORLD, ierr)

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
          CALL save_data('out'//trim(cmd_arg(2))//'.nc', niter, ncent, init_centroids, centroids, counts_master, sse_master, movement)
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
END PROGRAM for_kmeans
