SUBROUTINE cluster_label(vectors, centroids, nvec, ncent, ndim, clusters, label_count)
       integer, intent(in) :: nvec, ncent, ndim
       real, intent(in) :: vectors(nvec,ndim)
       real, intent(in) :: centroids(ncent,ndim)
       real, intent(out) :: clusters(ncent,ndim)
       integer, intent(out) :: label_count(ncent)
       integer :: labels(nvec)
       real :: dist(ncent)
       real :: vec(ndim)
       integer :: ml(1)

       DO i=1,nvec
          dist = 0
          vec = vectors(i,:)
          DO j=1, ncent
              dist(j) = SUM((vec - centroids(j,:))**2)
          END DO
          ml = MINLOC(dist)
          labels(i) = ml(1)
       END DO
       
       clusters = 0.0
       label_count = 0
       DO i=1,nvec
          vec = vectors(i,:)
          clusters(labels(i),:) = clusters(labels(i),:) + vec(:)
          label_count(labels(i)) = label_count(labels(i)) + 1
       END DO

END SUBROUTINE cluster_label
