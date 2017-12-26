subroutine compare(i,j,k, pctaudistnorm, histograms, masks, wsnums, wscount)
    integer, intent(in) :: i,j,k
    real, intent(in) :: pctaudistnorm(i,j,k)
    real, intent(in) :: histograms(12,i,j)
    logical, intent(in) :: masks(i,j,k)
    integer, intent(inout) :: wsnums(k,12)
    integer, intent(inout) :: wscount(k)
    integer :: m, n
    real, dimension(12) :: distances
    real, dimension(i,j) :: tmparr
    do m=1, k
        if (all(masks(:,:,m))) then
            cycle
        end if
        do n=1, 12
            where (.not. masks(:,:,m))
                tmparr = pctaudistnorm(:,:,m)-histograms(n,:,:)
            end where
            where (masks(:,:,m))
                tmparr = 0
            end where
            distances(n) = sqrt(sum(tmparr**2))
        end do
        wsnums(m,minloc(distances, 1)) = wsnums(m, minloc(distances,1)) + 1
        wscount(m) = wscount(m)+1
    end do
end subroutine compare


subroutine fdistance(array1, array2, mask, m, n, outnum)
    integer, intent(in) :: m, n
    real, dimension(m,n), intent(in) :: array1
    real, dimension(m,n), intent(in) :: array2
    logical, dimension(m,n), intent(in) :: mask
    real, intent(out) :: outnum
    real, dimension(m,n) :: array
    where (.not. mask)
         array = array1-array2
    end where
    where (mask)
         array = 0
    end where
    outnum = sqrt(sum(array**2))
end subroutine fdistance


subroutine fnorm(array1, val, mask, m, n, array)
    integer, intent(in) :: m, n
    real, dimension(m,n), intent(in) :: array1
    real, intent(in) :: val
    logical, dimension(m,n), intent(in) :: mask
    real, dimension(m,n), intent(out) :: array
    where (.not. mask)
         array = array1/val
    end where
end subroutine fnorm

subroutine fnormtotal(pctaudist, ntotal, masks, i, j, k, pctaudistnorm)
    integer, intent(in) :: i,j,k
    real, dimension(i,j,k), intent(in) :: pctaudist
    real, dimension(k), intent(in) :: ntotal
    logical, dimension(i,j,k), intent(in) :: masks
    real, dimension(i,j,k), intent(inout) :: pctaudistnorm
    integer :: m
    real, dimension(i,j) :: array
    do m = 1, k
        pctaudistnorm(:,:,m) = call subroutine fnorm(pctaudist(:,:,m),ntotal(m),masks(:,:,m), i, j)
    end do
    print *, 'hello'
end subroutine fnormtotal