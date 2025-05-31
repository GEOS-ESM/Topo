MODULE neighbor_search_mod
  IMPLICIT NONE

    type BlockType
      integer, allocatable :: indices(:)
      integer :: num_cells
     end type BlockType
  CONTAINS
   !====================================================================================
   ! Function: find_nearest_valid_neighbor
   !
   ! Purpose:
   !   Identify the nearest valid neighboring cell to a specified problematic cell
   !   within a structured block-based spatial search. Primarily used for correcting
   !   invalid or extreme values (fallback mechanism) by replacing them with
   !   geographically nearby, valid cell values.
   !
   ! Arguments:
   !   integer, intent(in) :: i
   !       Index of the problematic target cell that requires a valid neighbor.
   !
   !   real(8), intent(in) :: target_center_lon(:), target_center_lat(:)
   !       Arrays containing the longitudes (0–360 degrees) and latitudes (-90–90 degrees)
   !       of target grid cell centers.
   !
   !   logical, intent(in) :: valid_cells(:)
   !       Logical mask array indicating validity status of target grid cells:
   !       .true. if cell is valid, .false. otherwise.
   !
   !   integer, intent(in) :: num_lon_blocks, num_lat_blocks
   !       Total number of spatial blocks dividing the target grid longitudinally
   !       and latitudinally, respectively, used to accelerate the neighbor search.
   !
   !   real(8), intent(in) :: lon_block_size, lat_block_size
   !       Size (in degrees) of each spatial block along longitude and latitude.
   !
   !   type(BlockType), intent(in) :: blocks(:,:)
   !       2-dimensional array of BlockType structures, each containing indices
   !       of cells within its spatial block and the number of cells (`num_cells`).
   !
   !   integer, intent(in) :: max_search_radius
   !       Maximum number of blocks around the current block to search for a valid neighbor.
   !
   ! Returns:
   !   integer :: closest
   !       Index of the nearest valid neighboring cell found.
   !       Returns -1 if no valid neighbor is identified within the specified search radius.
   !
   ! Implementation details:
   !   - Spatial search is conducted incrementally by expanding radius around the problematic cell’s block.
   !   - Distances are computed based on simple squared Euclidean lat-lon differences.
   !   - The first valid cell found with minimal distance is selected as the neighbor.
   !   - Ensures robustness in stretched and irregular grids by systematically
   !     and efficiently locating fallback values.
   !
   ! Notes:
   !   - This function is critical for maintaining physical consistency in the final
   !     topography data, preventing anomalous elevations due to invalid cells seen with stretched grid runs.
   !====================================================================================

    function find_nearest_valid_neighbor(i, target_center_lon, target_center_lat, valid_cells, &
                                         num_lon_blocks, num_lat_blocks, lon_block_size, lat_block_size, &
                                         blocks, max_search_radius) result(closest)
        integer, intent(in) :: i, num_lon_blocks, num_lat_blocks
        real(8), intent(in) :: target_center_lon(:), target_center_lat(:), lon_block_size, lat_block_size
        logical, intent(in) :: valid_cells(:)
        type(BlockType), intent(in) :: blocks(:,:)
        integer, intent(in) :: max_search_radius
        integer :: closest
        integer :: ii, jj, k, idx, iblock, jblock, search_radius
        real(8) :: min_dist, dist, dlat, dlon
        logical :: found_valid

        min_dist = 1.0d30
        closest = -1
        found_valid = .false.
        search_radius = 1

        iblock = min(num_lon_blocks, max(1, int(target_center_lon(i) / lon_block_size) + 1))
        jblock = min(num_lat_blocks, max(1, int((target_center_lat(i) + 90.0d0) / lat_block_size) + 1))

        do while (.not. found_valid .and. search_radius <= max_search_radius)
          do ii = max(1, iblock - search_radius), min(num_lon_blocks, iblock + search_radius)
              do jj = max(1, jblock - search_radius), min(num_lat_blocks, jblock + search_radius)
                  do k = 1, blocks(ii,jj)%num_cells
                      idx = blocks(ii,jj)%indices(k)
                      if (valid_cells(idx) .and. idx /= i) then  ! exclude the problematic cell itself
                          dlat = target_center_lat(idx) - target_center_lat(i)
                          dlon = target_center_lon(idx) - target_center_lon(i)
                          dist = dlat**2 + dlon**2
                          if (dist < min_dist) then
                              min_dist = dist
                              closest = idx
                              found_valid = .true.
                          end if
                      endif
                  end do
              end do
          end do
            if (.not. found_valid) search_radius = search_radius + 1
        end do

       ! ! Simple confirmation print:
       ! if (found_valid) then
       !     write(*,*) "Nearest neighbor found for cell", i, "->", closest
       ! else
       !     write(*,*) "No nearest neighbor found for cell", i
       !endif

    end function find_nearest_valid_neighbor

END MODULE neighbor_search_mod

