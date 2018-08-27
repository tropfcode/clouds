subroutine reshape_data(pctaudist, masks, ntotal, sqlonbeg, sqlonend, eqlat, &
                        pctau_reshaped, masks_reshaped, ntotal_reshaped)
    !!!!!DESCRIPTION!!!!!
    ! Reshapes n_pctaudist,ntotal, and n_pctaudist mask from 
    ! raw HGG data in preparation to reduce resolution from
    ! 1x1 to 2x2 degrees. Conversion to Equal-Angle is also made for all.
    
    !!!!!RETURN VALUES!!!!!
    ! All return values are in Equal-Angle
    ! pctau_reshaped : n_pctaudist reshaped to (360,180,42)
    ! masks_reshaped : n_pctaudist mask reshaped to (360, 180, 42)
    ! ntotal_reshaped: ntotal reshaped to (360,180,42) for normalization
    
    real, dimension(6, 7, 41252), intent(in)      :: pctaudist
    logical, dimension(6, 7, 41252), intent(in)   :: masks
    real, dimension(41252), intent(in)            :: ntotal
    integer, dimension(41252), intent(in)         :: sqlonbeg
    integer, dimension(41252), intent(in)         :: sqlonend
    integer, dimension(41252), intent(in)         :: eqlat
    real, dimension(360, 180, 42), intent(out)    :: pctau_reshaped
    logical, dimension(360, 180, 42), intent(out) :: masks_reshaped
    real, dimension(360, 180, 42), intent(out)    :: ntotal_reshaped
    integer                                       :: i,j
    
    ! Reshape data and convert to Equal-Angle
    ! See algorithm for Equal-Angle conversion in HGG netCDF metadata.
    do i=1, 41252
        do j=sqlonbeg(i), sqlonend(i)
            pctau_reshaped(j, eqlat(i), :) = RESHAPE(pctaudist(:,:,i), (/42/))
            masks_reshaped(j, eqlat(i), :) = RESHAPE(masks(:,:,i), (/42/))
            ntotal_reshaped(j, eqlat(i),:) = ntotal(i)
        end do
    end do
end subroutine reshape_data


subroutine reduce_pctau(pctaudist_reshape, mask_reshape, pctaudist_reduced)
    ! ATTN: Input data must be reshaped before using this
    !!!!!DESCRIPTION!!!!!
    ! Reduces resolution of n_pctaudist HGG variable.
    ! Does this by turning quadrant of pixels into single pixel.
    
    !!!!!RETURN VALUES!!!!!
    ! pctaudist_reduced: 2x2 degree n_pctaudist with shape (180,90,42)
    
    real, dimension(360, 180, 42), intent(in)    :: pctaudist_reshape
    logical, dimension(360, 180, 42), intent(in) :: mask_reshape
    real, dimension(180, 90, 42), intent(out)    :: pctaudist_reduced
    real, dimension(360, 180, 42)                :: pctau_tmp
    real, dimension(180, 180, 42)                :: row_tmp
    integer :: i,j,row,col

    pctau_tmp = 0
    
    ! Reduce row dimension in half by adding adjacent row pairs 
    do i = 1, 42
        row = 1
        do j=1, 180
            ! Add together values that are not masked
            where(.NOT. mask_reshape(row:row+1,:,i))
                pctau_tmp(row:row+1,:,i) = pctaudist_reshape(row:row+1,:,i)
            end where
            row_tmp(j,:,i) = pctau_tmp(row,:,i)+pctau_tmp(row+1,:,i)
            row = row + 2
        end do
    end do

    ! Reduce column dimension in half by adding adjacent column pairs
    do i = 1, 42
        col = 1
        do j = 1, 90
            pctaudist_reduced(:,j,i) = row_tmp(:,col,i)+row_tmp(:,col+1,i)
            col = col + 2
        end do
    end do
        
end subroutine reduce_pctau

subroutine reduce_mask(mask_reshape, mask_reduced)
    ! ATTN: Input data must be reshaped before using this
    !!!!!DESCRIPTION!!!!!
    ! Reduces resolution of the mask of n_pctaudist HGG variable.
    ! If quadrant of pixels are all logical TRUE, then the reduced pixel
    ! is TRUE, otherwise if a single FALSE value exists the reduced pixel
    ! is FALSE (not masked).
    
    !!!!!RETURN VALUES!!!!!
    ! mask_reduced: 2x2 degree mask of n_pctaudist with shape (180,90,42)
    logical, dimension(360, 180, 42), intent(in) :: mask_reshape
    logical, dimension(180, 90, 42), intent(out) :: mask_reduced
    integer                                      :: i,j,k,row,col
    
    ! Make all values masked and unmask accordingly
    mask_reduced(:,:,:) = .TRUE.
    do k=1, 42
        row = 1
        do i=1,180
            col = 1
            do j=1,90
                ! Analyze adjacent quadrants of pixels to determine
                ! if a single FALSE value exists.
                if (.NOT. all(mask_reshape(row:row+1,col:col+1,k))) then
                    mask_reduced(i,j,k) = .FALSE.
                end if
                col = col + 2
            end do
            row = row + 2
        end do
    end do
end subroutine reduce_mask


subroutine reduce_ntotal(ntotal_reshape, mask_reshape, ntotal_reduced)
    ! ATTN: Input data must be reshaped before using this
    !!!!!DESCRIPTION!!!!!
    ! Reduces resolution of the mask of n_total HGG variable.
    ! Although n_total HGG variable is 1D, it is extended to match
    ! the shape of reduced n_pctaudist for easier normalization
    
    !!!!!RETURN VALUES!!!!!
    ! ntotal_reduced: 2x2 degree ntotal of shape (180, 90, 42)
    real, dimension(360, 180, 42), intent(in)    :: ntotal_reshape
    logical, dimension(360, 180, 42), intent(in) :: mask_reshape
    real, dimension(180, 90, 42), intent(out)    :: ntotal_reduced
    real, dimension(360, 180, 42)                :: ntotal_tmp
    real, dimension(180, 180, 42)                :: ntotal_row_tmp
    integer                                      :: i,j,row,col

    ! Make all values zero and then copy over unmasked values.
    ! This makes summing quadrant of pixels together easier 
    ! as addition by zero will not artificially inflate 
    ! pixel counts
    ntotal_tmp = 0
    ntotal_reduced = 0
    where (.NOT. mask_reshape)
        ntotal_tmp = ntotal_reshape
    end where
    
    ! Reduce row dimension in half by adding adjacent row pairs 
    do i = 1, 42
        row = 1
        do j=1, 180
            ! Add together values that are not masked
            ntotal_row_tmp(j,:,i) = ntotal_tmp(row,:,i)+ntotal_tmp(row+1,:,i)
            row = row + 2
        end do
    end do

    do i = 1, 42
        col = 1
        do j = 1, 90
            ! Reduce column dimension in half by adding adjacent column pairs
            ntotal_reduced(:,j,i) = ntotal_row_tmp(:,col,i)+ntotal_row_tmp(:,col+1,i)
            col = col + 2
        end do
    end do
end subroutine reduce_ntotal

subroutine fnorm(pctaudist, ntotal, masks, pctau_norm)
    ! ATTN: Must reshape and reduce input data first
    !!!!!DESCRIPTION!!!!!
    ! Converts n_pctaudist values to percentages

    !!!!!RETURN VALUES!!!!!
    ! pctau_norm: Original input array with decimal values less than 1 (percentages)

    real, dimension(180, 90, 42), intent(in)        :: pctaudist
    integer, dimension(180, 90, 42), intent(in)     :: ntotal
    logical, dimension(180, 90, 42), intent(in)     :: masks
    real, dimension(180, 90, 42), intent(out)       :: pctau_norm
    integer                                         :: i,j
    
    do i=1, 180
        do j=1, 90
            ! Convert only the unmasked values
            where (.not. masks(i,j,:))
                pctau_norm(i,j,:) = pctaudist(i,j,:)/ntotal(i,j,:)
            end where
        end do
    end do
end subroutine fnorm

subroutine findws(dist, wsdata, masks, wsnums)
    !!!!!DESCRIPTION!!!!!
    ! FILL IN LATER
    
    real, dimension(180, 90, 42), intent(in)     :: dist
    real, dimension(12, 6, 7), intent(in)        :: wsdata
    logical, dimension(180, 90, 42), intent(in)  :: masks
    integer, dimension(180, 90, 12), intent(out) :: wsnums
    real, dimension(12, 42)                      :: weather_states
    real, dimension(42)                          :: tmparr
    real, dimension(12)                          :: distances
    integer :: i, j, k, ws
    tmparr = 0
    wsnums = 0
    do m=1, 12
        weather_states(m,:) = RESHAPE(wsdata(m,:,:), (/42/))
    end do
    do i=1, 180
        do j= 1, 90
            distances = 0
            if (all(masks(i,j,:))) then
                cycle
            end if
            do k=1, 12
                where (masks(i,j,:))
                    tmparr = 0.0
                end where
                where (.not. masks(i,j,:))
                    tmparr = dist(i,j,:)-weather_states(k,:)
                end where
                distances(k) = sqrt(sum(tmparr**2))
            end do
            ws = minloc(distances,1)
            wsnums(i,j,ws) = wsnums(i,j,ws)+1
            print *, distances
        end do
    end do
end subroutine findws