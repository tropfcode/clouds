!Subroutines created for conversion and use in Python
!via f2py program

subroutine fnorm(pctaudist, ntotal, masks)
!!!!!!DESCRIPTION!!!!!!
! Convert pctaudist values to percentages by dividing
! values that correspond to false in masks by ntotal
!
!!!!!!RETURN VALUES!!!!!!
! pctaudist: Original input array with decimal values less than 1
    real, dimension(6, 7, 41252), intent(inout) :: pctaudist
!f2py intent(in,out) :: pctaudist
    integer, dimension(41252), intent(in) :: ntotal
    logical, dimension(6, 7, 41252), intent(in) :: masks
    do i=1, 41252
        where (.not. masks(:,:,i))
            pctaudist(:,:,i) = pctaudist(:,:,i)/ntotal(i)
        end where
    end do
end subroutine fnorm


subroutine findws(dist, histograms, masks, sqlonbeg, sqlonend, eqlat, k, wsnums)!, counts, sse)
!!!!!!DESCRIPTION!!!!!! 
!Find weather state of each gridbox in dist by finding
!the shortest euclidean distance of each gridbox histogram to each 
!weather state.
!
!!!!!INPUT VARIABLES!!!!!
!dist: n_pctaudist variable from ISCCP HGG data normalized by n_total variable
!histograms: Weather states to compute euclidean distanced against
!masks: Mask for dist variable. True values refer to no data
!sqlonbeg: sqlon_beg variable from ISCCP HGG data
!sqlonend: sqlon_end variable from ISCCP HGG data
!!eqlat: eqlat_index variable from ISCCP HGG data

!!!!!!RETURN VALUES!!!!!!
!wsnums: Count of each weather state for each gridbox
    
    integer, intent(in) :: k
    real, dimension(6, 7, 41252), intent(in) :: dist
    logical, dimension(6, 7, 41252), intent(in) :: masks
    integer, dimension(41252), intent(in) :: sqlonbeg
    integer, dimension(41252), intent(in) :: sqlonend
    integer, dimension(41252), intent(in) :: eqlat
    real, dimension(k, 6, 7), intent(in) :: histograms
    integer, dimension(41252), intent(out) :: wsnums
!f2py intent(out) :: wsnums
    real, dimension(6, 7) :: tmparr
    real, dimension(10) :: distances
    integer :: i, j, ws, wsplacement(k+2)

    tmparr = 0
    distances = 0
    wsnums = 0
    
    ! Create WS gridbox array
    DO i=1, k+2
       IF (i > k+1) THEN
          wsplacement(i) = -127
       ELSE
          wsplacement(i) = i
       END IF
    END DO

    do i=1, 41252
        if (all(masks(:,:,i))) then
            ws = k+2
        else if (all(dist(:,:,i)==0)) then
            ws = k+1
        else
            do j=1, k
                tmparr = dist(:,:,i)-histograms(j,:,:)
                distances(j) = sqrt(sum(tmparr**2))
            end do
            ws = minloc(distances,1)
        end if
        wsnums(i) = wsplacement(ws)
    end do
end subroutine findws
