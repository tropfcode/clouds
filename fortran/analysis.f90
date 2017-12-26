!Subroutines created for conversion and use in Python
!via f2py program

subroutine fnorm(pctaudist, ntotal, masks)
!!!!!!DESCRIPTION!!!!!!
!Convert pctaudist values to percentages by dividing
!values that correspond to false in masks by ntotal
!
!!!!!!RETURN VALUES!!!!!!
!pctaudist: Original input array with decimal values less than 1
    real, dimension(6, 7, 41252), intent(inout) :: pctaudist
!f2py intent(in,out) :: pctaudist
    integer, dimension(41252), intent(in) :: ntotal
    logical, dimension(6, 7, 41252), intent(in) :: masks
    pctaudist_norm = 0
    do i=1, 41252
        where (.not. masks(:,:,i))
            pctaudist(:,:,i) = pctaudist(:,:,i)/ntotal(i)
        end where
    end do
end subroutine fnorm

subroutine findws(dist, histograms, masks, wsnums, wscount, sqlonbeg, sqlonend, eqlat)
!!!!!!DESCRIPTION!!!!!! 
!Find weather state of each gridbox in dist by finding
!the shortest distance of each gridbox histogram to each 
!weather state in histograms.
!
!!!!!!RETURN VALUES!!!!!!
!wsnums: Count of each weather state for each gridbox
!wscount: Total count for each gridbox
    real, dimension(6, 7, 41252), intent(in) :: dist
    real, dimension(12, 6, 7), intent(in) :: histograms
    logical, dimension(6, 7, 41252), intent(in) :: masks
    integer, dimension(360, 180, 12), intent(inout) :: wsnums
!f2py intent(in, out) :: wsnums
    integer, dimension(360, 180), intent(inout) :: wscount
!f2py intent(in, out) :: wscount
    integer, dimension(41252), intent(in) :: sqlonbeg
    integer, dimension(41252), intent(in) :: sqlonend
    integer, dimension(41252), intent(in) :: eqlat
    real, dimension(6, 7) :: tmparr
    real, dimension(12) :: distances
    integer :: first, last, diff, i, j, ws
    tmparr = 0
    distances = 0
    do i=1, 41252
        diff = abs(sqlonbeg(i)-sqlonend(i))!+1
        last = last + diff
        if (all(masks(:,:,i))) then
            cycle
        end if
        do j=1, 12
            where (masks(:,:,i))
                tmparr = 0.0
            end where
            where (.not. masks(:,:,i))
                tmparr = dist(:,:,i)-histograms(j,:,:)
            end where
            distances(j) = sqrt(sum(tmparr**2))
        end do
        ws = minloc(distances,1)
        wsnums(sqlonbeg(i):sqlonend(i), eqlat(i), ws) = wsnums(sqlonbeg(i):sqlonend(i), eqlat(i), ws)+1
        wscount(sqlonbeg(i):sqlonend(i), eqlat(i)) = wscount(sqlonbeg(i):sqlonend(i), eqlat(i))+1 
    end do
end subroutine findws